#!/usr/bin/env python3

import argparse
import itertools
import math

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import null_space

Q = 3


def basis(n, p):
    return list(itertools.combinations(range(n), p))


def wedge_sign(A, B):
    if set(A) & set(B):
        return 0
    inversions = sum(a > b for a in A for b in B)
    return -1 if inversions % 2 else 1


def random_couplings(n, rng):
    variance = math.factorial(Q - 1) / n ** (Q - 1)
    return {
        inds: np.sqrt(variance / 2) * (rng.normal() + 1j * rng.normal())
        for inds in itertools.combinations(range(n), Q)
    }


def restrict_couplings(C, n):
    return {inds: value for inds, value in C.items() if max(inds) < n}


def build_Q(n, p, C):
    src = basis(n, p)
    if p + Q > n:
        return np.zeros((0, len(src)), dtype=complex)

    dst = basis(n, p + Q)
    dst_index = {x: i for i, x in enumerate(dst)}
    M = np.zeros((len(dst), len(src)), dtype=complex)

    for j, S in enumerate(src):
        for T, c in C.items():
            if set(S) & set(T):
                continue
            U = tuple(sorted(T + S))
            M[dst_index[U], j] += c * wedge_sign(T, S)

    return M


def bps_basis(n, p, C):
    """Orthonormal basis for ker Q_p intersect ker(Q_{p-Q}^dagger)."""
    dim = len(basis(n, p))
    blocks = []

    Qp = build_Q(n, p, C)
    if Qp.shape[0] > 0:
        blocks.append(Qp)

    if p >= Q:
        Qprev = build_Q(n, p - Q, C)
        if Qprev.shape[0] > 0:
            blocks.append(Qprev.conj().T)

    if not blocks:
        return np.eye(dim, dtype=complex)

    return null_space(np.vstack(blocks), rcond=1e-10)


def embed_basis(n, p, B, channel):
    """Embed the old BPS basis into the appropriate N+1 sector."""
    old_basis = basis(n, p)
    P = p if channel == "down" else p + 1
    new_basis = basis(n + 1, P)
    new_index = {x: i for i, x in enumerate(new_basis)}

    E = np.zeros((len(new_basis), len(old_basis)), dtype=complex)

    for j, S in enumerate(old_basis):
        Snew = S if channel == "down" else tuple(sorted(S + (n,)))
        E[new_index[Snew], j] = 1.0

    return E @ B


def KT_eigenvalues(N, p, Cbig, channel):
    """Return the K_T eigenvalues for one disorder realization."""
    Cold = restrict_couplings(Cbig, N)
    B_old = bps_basis(N, p, Cold)

    P = p if channel == "down" else p + 1
    B_new = bps_basis(N + 1, P, Cbig)

    E = embed_basis(N, p, B_old, channel)
    K = E.conj().T @ B_new @ B_new.conj().T @ E
    K = (K + K.conj().T) / 2

    return np.clip(np.linalg.eigvalsh(K), 0.0, 1.0)


def wachter_parameters(N, p, channel, d):
    """Parameters in the convention used in the WIP."""
    D = math.comb(N, p)
    a = d / D

    if channel == "down":
        b = (N + 1 - p) / (N + 1)
    else:
        # The up channel has p -> p+1 in the enlarged theory.
        b = (N + 1 - (p + 1)) / (N + 1)

    return a, b, D


def wachter_density(x, a, b):
    x = np.asarray(x, dtype=float)

    lambda_minus = (
        np.sqrt(a * (1 - b)) - np.sqrt(b * (1 - a))
    ) ** 2
    lambda_plus = (
        np.sqrt(a * (1 - b)) + np.sqrt(b * (1 - a))
    ) ** 2

    rho = np.zeros_like(x)
    mask = (x > lambda_minus) & (x < lambda_plus)
    rho[mask] = np.sqrt(
        (lambda_plus - x[mask]) * (x[mask] - lambda_minus)
    ) / (2 * np.pi * a * x[mask] * (1 - x[mask]))

    return rho, lambda_minus, lambda_plus


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare the SYK K_T spectrum with the Wachter law."
    )

    parser.add_argument("N", type=int, help="Starting system size N.")
    parser.add_argument("p", type=int, help="Fermion number p.")
    parser.add_argument(
        "channel", choices=["up", "down"], help="Uplift channel."
    )
    parser.add_argument(
        "-r", "--realizations", type=int, default=1,
        help="Number of disorder realizations (default: 1)."
    )
    parser.add_argument(
        "-o", "--output", default=None,
        help="Output PDF filename."
    )
    parser.add_argument(
        "--seed", type=int, default=1234,
        help="Random seed (default: 1234)."
    )
    parser.add_argument(
        "--bins", type=int, default=40,
        help="Number of histogram bins (default: 40)."
    )
    parser.add_argument(
        "--keep-endpoints", action="store_true",
        help="Keep eigenvalues at lambda=0 and 1 in the histogram."
    )

    return parser.parse_args()


def main():
    args = parse_args()

    if args.N < 1:
        raise ValueError("N must be positive")
    if not (0 <= args.p <= args.N):
        raise ValueError("p must satisfy 0 <= p <= N")
    if args.channel == "up" and args.p >= args.N:
        raise ValueError("For the up channel require p < N")
    if args.realizations < 1:
        raise ValueError("--realizations must be >= 1")
    if args.bins < 1:
        raise ValueError("--bins must be >= 1")

    outfile = args.output or (
        f"KT_N{args.N}_p{args.p}_{args.channel}_R{args.realizations}.pdf"
    )
    if not outfile.lower().endswith(".pdf"):
        outfile += ".pdf"

    # First realization only used to determine the N-theory BPS dimension d.
    rng = np.random.default_rng(args.seed)
    Cbig = random_couplings(args.N + 1, rng)
    Cold = restrict_couplings(Cbig, args.N)
    d = bps_basis(args.N, args.p, Cold).shape[1]
    a, b, D = wachter_parameters(args.N, args.p, args.channel, d)

    print(f"N = {args.N}")
    print(f"p = {args.p}")
    print(f"channel = {args.channel}")
    print(f"realizations = {args.realizations}")
    print(f"d = {d}")
    print(f"D = {D}")
    print(f"a = {a:.6f}")
    print(f"b = {b:.6f}")

    all_vals = []

    # Re-use the first realization, then generate the rest.
    for r in range(args.realizations):
        if r == 0:
            C = Cbig
        else:
            C = random_couplings(args.N + 1, rng)

        vals = KT_eigenvalues(args.N, args.p, C, args.channel)
        all_vals.extend(vals)
        print(f"realization {r + 1}/{args.realizations}: {len(vals)} eigenvalues")

    all_vals = np.asarray(all_vals)

    # Remove protected endpoint atoms by default, matching the WIP bulk comparison.
    tol = 1e-10
    atom0 = np.count_nonzero(all_vals <= tol) / len(all_vals)
    atom1 = np.count_nonzero(all_vals >= 1 - tol) / len(all_vals)

    if args.keep_endpoints:
        plot_vals = all_vals
        plot_note = "all eigenvalues"
    else:
        plot_vals = all_vals[(all_vals > tol) & (all_vals < 1 - tol)]
        plot_note = "interior only"

    rho_x = np.linspace(1e-8, 1 - 1e-8, 10000)
    rho_w, lambda_minus, lambda_plus = wachter_density(rho_x, a, b)

    # Normalize the Wachter continuous piece to unit integral for a like-for-like
    # comparison with density=True on the interior ED histogram.
    support = (rho_x >= lambda_minus) & (rho_x <= lambda_plus)
    wachter_mass = np.trapezoid(rho_w[support], rho_x[support])
    rho_plot = rho_w / wachter_mass if wachter_mass > 0 else rho_w

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.hist(
        plot_vals,
        bins=args.bins,
        range=(0, 1),
        density=True,
        alpha=0.50,
        label="ED",
    )

    ax.plot(
        rho_x,
        rho_plot,
        color="black",
        linewidth=2.5,
        label="Wachter",
    )

    # Mark Wachter edges.
    ax.axvline(lambda_minus, color="black", linestyle=":", linewidth=1)
    ax.axvline(lambda_plus, color="black", linestyle=":", linewidth=1)

    ax.set_xlim(0, 1)
    ax.set_xlabel(r"$\lambda$")
    ax.set_ylabel(r"$\rho(\lambda)$")

    ax.set_title(
        rf"$K_T$: $N={args.N}\to{args.N+1}$, $p={args.p}$, {args.channel}"
        "\n"
        rf"$a={a:.4f}$, $b={b:.4f}$, $R={args.realizations}$"
    )

    ax.legend()

    # Small text box with precisely the quantities useful for comparison.
    ax.text(
        0.02, 0.97,
        rf"$\lambda_-={lambda_minus:.4f}$\n"
        rf"$\lambda_+={lambda_plus:.4f}$\n"
        rf"atom $\lambda=0$: {atom0:.3f}\n"
        rf"atom $\lambda=1$: {atom1:.3f}\n"
        rf"{plot_note}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=10,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    fig.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)

    print(f"Wachter lambda_- = {lambda_minus:.8f}")
    print(f"Wachter lambda_+ = {lambda_plus:.8f}")
    print(f"atom fraction at 0 = {atom0:.6f}")
    print(f"atom fraction at 1 = {atom1:.6f}")
    print(f"Saved: {outfile}")


if __name__ == "__main__":
    main()
