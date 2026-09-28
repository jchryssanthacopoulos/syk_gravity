#!/usr/bin/env python3
r"""
memwatch.py -- run a command under a memory watchdog.  REQUIRED for every intensive process in this project
(see CLAUDE.md, "Numerical Experiments").  Standard library only.

    .venv/bin/python scripts/memwatch.py --limit-gb 12 [--log FILE] -- .venv/bin/python scripts/run_x.py --N 16
    .venv/bin/python scripts/memwatch.py --status          # jobs currently registered, reservations, usage

What it enforces
----------------
1. Per-job hard cap (--limit-gb).  The whole process tree of the command (it runs in its own process group) is
   measured every --poll seconds.  On macOS the metric is the *physical footprint* (proc_pid_rusage), which
   includes compressed pages; plain RSS undercounts badly on macOS once the compressor kicks in.  On Linux: RSS.
2. Project budget (--budget-gb, default $MEMWATCH_BUDGET_GB or 30).  Each job *reserves* its cap in a shared
   registry ($MEMWATCH_DIR, default $TMPDIR/memwatch_<uid>).  A launch whose cap would push the sum of reservations
   of live jobs over the budget is refused (exit 3), or waits for room with --wait.  So caps partition the budget:
   concurrent jobs can never collectively overcommit.  At runtime, if the *measured* total of all watched jobs
   exceeds the budget, the largest job is killed.
3. System floor (--min-free-gb, default 1.5).  If the OS reports less available memory than this, the job is
   killed regardless of its own usage (protects against memory used by things outside memwatch).

On a kill: SIGTERM to the process group, SIGKILL after 5 s; exit code 137.  Otherwise the command's exit code.
The peak footprint and the reason for any kill are printed to stderr and appended to --log (if given).
"""
from __future__ import annotations

import argparse
import ctypes
import fcntl
import json
import os
import signal
import subprocess
import sys
import tempfile
import time

GB = 1024 ** 3
KILLED = 137


# ----------------------------------------------------------------------------- memory measurement
class _RUsageV0(ctypes.Structure):
    _fields_ = [("uuid", ctypes.c_uint8 * 16)] + [(n, ctypes.c_uint64) for n in (
        "user_time", "system_time", "pkg_idle_wkups", "interrupt_wkups", "pageins", "wired_size",
        "resident_size", "phys_footprint", "proc_start_abstime", "proc_exit_abstime")]


_libproc = None
if sys.platform == "darwin":
    try:
        _libproc = ctypes.CDLL("/usr/lib/libproc.dylib")
    except OSError:
        _libproc = None


def proc_bytes(pid: int) -> int:
    """Physical footprint (macOS) or RSS (Linux / fallback) of one process, in bytes; 0 if it is gone."""
    if _libproc is not None:
        buf = _RUsageV0()
        if _libproc.proc_pid_rusage(pid, 0, ctypes.byref(buf)) == 0:
            return int(buf.phys_footprint)
    try:
        with open(f"/proc/{pid}/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) * 1024
    except OSError:
        pass
    try:
        out = subprocess.run(["ps", "-o", "rss=", "-p", str(pid)], capture_output=True, text=True).stdout.strip()
        return int(out) * 1024 if out else 0
    except (OSError, ValueError):
        return 0


def tree_pids(root: int) -> list[int]:
    out = subprocess.run(["ps", "-A", "-o", "pid=,ppid="], capture_output=True, text=True).stdout
    kids: dict[int, list[int]] = {}
    for line in out.splitlines():
        parts = line.split()
        if len(parts) == 2:
            kids.setdefault(int(parts[1]), []).append(int(parts[0]))
    pids, stack = [], [root]
    while stack:
        p = stack.pop()
        pids.append(p)
        stack.extend(kids.get(p, []))
    return pids


def tree_bytes(root: int) -> int:
    return sum(proc_bytes(p) for p in tree_pids(root))


def system_available_bytes() -> float:
    """Memory the OS considers available (macOS: kern.memorystatus_level % of RAM; Linux: MemAvailable)."""
    if sys.platform == "darwin":
        try:
            lvl = subprocess.run(["sysctl", "-n", "kern.memorystatus_level"], capture_output=True, text=True).stdout
            tot = subprocess.run(["sysctl", "-n", "hw.memsize"], capture_output=True, text=True).stdout
            return float(lvl) / 100.0 * float(tot)
        except (OSError, ValueError):
            return float("inf")
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) * 1024.0
    except OSError:
        pass
    return float("inf")


# ----------------------------------------------------------------------------- shared registry
def registry_dir() -> str:
    d = os.environ.get("MEMWATCH_DIR") or os.path.join(tempfile.gettempdir(), f"memwatch_{os.getuid()}")
    os.makedirs(d, exist_ok=True)
    return d


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True


def read_registry() -> list[dict]:
    """Entries of live watchers (stale files of dead watchers are removed)."""
    d, out = registry_dir(), []
    for name in os.listdir(d):
        if not name.endswith(".json"):
            continue
        path = os.path.join(d, name)
        try:
            with open(path) as f:
                e = json.load(f)
        except (OSError, ValueError):
            continue
        if _alive(int(e["watcher"])):
            out.append(e)
        else:
            try:
                os.remove(path)
            except OSError:
                pass
    return out


def write_entry(entry: dict) -> None:
    path = os.path.join(registry_dir(), f"{entry['watcher']}.json")
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(entry, f)
    os.replace(tmp, path)


def remove_entry(watcher: int) -> None:
    try:
        os.remove(os.path.join(registry_dir(), f"{watcher}.json"))
    except OSError:
        pass


class _Lock:
    def __enter__(self):
        self.f = open(os.path.join(registry_dir(), ".lock"), "w")
        fcntl.flock(self.f, fcntl.LOCK_EX)
        return self

    def __exit__(self, *exc):
        fcntl.flock(self.f, fcntl.LOCK_UN)
        self.f.close()


def status() -> None:
    reg = read_registry()
    print(f"registry: {registry_dir()}")
    print(f"system available: {system_available_bytes() / GB:.1f} GB")
    if not reg:
        print("no watched jobs")
        return
    tot_res = sum(e["limit_gb"] for e in reg)
    tot_use = sum(e.get("bytes", 0) for e in reg) / GB
    for e in sorted(reg, key=lambda e: e["started"]):
        print(f"  watcher {e['watcher']} child {e.get('child')}: {e.get('bytes', 0) / GB:6.2f} GB used "
              f"(peak {e.get('peak', 0) / GB:.2f}) / cap {e['limit_gb']:.1f} GB  budget {e['budget_gb']:.0f}  "
              f"| {e['cmd'][:100]}")
    print(f"total: {tot_use:.2f} GB used, {tot_res:.1f} GB reserved")


# ----------------------------------------------------------------------------- main loop
def kill_group(proc: subprocess.Popen) -> None:
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.wait()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--limit-gb", type=float, help="hard cap for this job's process tree (also its reservation)")
    ap.add_argument("--budget-gb", type=float, default=float(os.environ.get("MEMWATCH_BUDGET_GB", 30)),
                    help="total for all memwatch jobs of this user (default $MEMWATCH_BUDGET_GB or 30)")
    ap.add_argument("--min-free-gb", type=float, default=1.5, help="kill if the OS reports less available memory")
    ap.add_argument("--poll", type=float, default=1.0, help="seconds between measurements")
    ap.add_argument("--wait", action="store_true", help="wait for budget room instead of refusing to launch")
    ap.add_argument("--log", help="append a one-line summary (JSON) to this file")
    ap.add_argument("--status", action="store_true", help="print registered jobs and exit")
    ap.add_argument("cmd", nargs=argparse.REMAINDER, help="-- command ...")
    a = ap.parse_args()
    if a.status:
        status()
        return 0
    cmd = a.cmd[1:] if a.cmd[:1] == ["--"] else a.cmd
    if not cmd or a.limit_gb is None:
        ap.error("need --limit-gb and a command after --")
    if a.limit_gb > a.budget_gb:
        print(f"memwatch: cap {a.limit_gb} GB exceeds the budget {a.budget_gb} GB; refusing", file=sys.stderr)
        return 3
    me = os.getpid()
    cmd_str = " ".join(cmd)

    # admission: reserve the cap within the budget (atomic under the registry lock)
    while True:
        with _Lock():
            reserved = sum(e["limit_gb"] for e in read_registry())
            if reserved + a.limit_gb <= a.budget_gb + 1e-9:
                entry = dict(watcher=me, child=None, limit_gb=a.limit_gb, budget_gb=a.budget_gb, cmd=cmd_str,
                             started=time.time(), bytes=0, peak=0)
                write_entry(entry)
                break
        msg = (f"memwatch: {reserved:.1f} GB already reserved of {a.budget_gb:.0f} GB; "
               f"cannot reserve {a.limit_gb:.1f} GB more")
        if not a.wait:
            print(msg + " (use --wait to queue)", file=sys.stderr)
            return 3
        print(msg + "; waiting", file=sys.stderr, flush=True)
        time.sleep(30)

    t0, peak, reason = time.time(), 0, None
    proc = subprocess.Popen(cmd, start_new_session=True)
    entry["child"] = proc.pid
    fwd = lambda s, f: kill_group(proc)                                    # noqa: E731
    signal.signal(signal.SIGTERM, fwd)
    signal.signal(signal.SIGINT, fwd)
    print(f"memwatch: pid {proc.pid} cap {a.limit_gb:.1f} GB, budget {a.budget_gb:.0f} GB: {cmd_str}",
          file=sys.stderr, flush=True)
    try:
        while proc.poll() is None:
            used = tree_bytes(proc.pid)
            peak = max(peak, used)
            entry.update(bytes=used, peak=peak)
            write_entry(entry)
            if used > a.limit_gb * GB:
                reason = f"job footprint {used / GB:.2f} GB > cap {a.limit_gb:.1f} GB"
            else:
                reg = read_registry()
                total = sum(e.get("bytes", 0) for e in reg)
                if total > a.budget_gb * GB and max(reg, key=lambda e: e.get("bytes", 0))["watcher"] == me:
                    reason = f"all watched jobs {total / GB:.2f} GB > budget {a.budget_gb:.0f} GB (this is largest)"
                elif system_available_bytes() < a.min_free_gb * GB:
                    reason = f"system available memory < {a.min_free_gb} GB"
            if reason:
                print(f"memwatch: KILLING pid {proc.pid}: {reason}", file=sys.stderr, flush=True)
                kill_group(proc)
                break
            time.sleep(a.poll)
    finally:
        rc = proc.wait() if reason is None else KILLED
        remove_entry(me)
    summary = dict(cmd=cmd_str, rc=rc, killed=reason, peak_gb=round(peak / GB, 3), cap_gb=a.limit_gb,
                   secs=round(time.time() - t0, 1), ended=time.strftime("%Y-%m-%dT%H:%M:%S"))
    print(f"memwatch: done rc={rc} peak {peak / GB:.2f} GB (cap {a.limit_gb:.1f}) in {summary['secs']:.0f}s"
          + (f"  KILLED: {reason}" if reason else ""), file=sys.stderr, flush=True)
    if a.log:
        with open(a.log, "a") as f:
            f.write(json.dumps(summary) + "\n")
    return rc


if __name__ == "__main__":
    sys.exit(main())
