"""
kuramoto_clean.py
=================
Generates 4 clean figures for the Kuramoto model report:

  fig1_phase_snapshots.png  – phases on the unit circle (4 K values)
  fig2_r_timeseries.png     – r(t) time series
  fig3_r_vs_K.png           – synchronization transition curve
  fig4_finite_size.png      – finite-size comparison

Usage:
    python kuramoto_clean.py

Requires: numpy, matplotlib
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

# ── Style ─────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family"      : "serif",
    "mathtext.fontset" : "dejavuserif",
    "font.size"        : 12,
    "axes.labelsize"   : 13,
    "axes.titlesize"   : 12,
    "legend.fontsize"  : 11,
    "xtick.labelsize"  : 11,
    "ytick.labelsize"  : 11,
    "figure.dpi"       : 150,
    "axes.spines.top"  : False,
    "axes.spines.right": False,
})

# Colour palette
TEAL   = "#00827F"
RED    = "#C0392B"
GRAY   = "#7F8C8D"
DARK   = "#2C3E50"
COLORS = ["#2980B9", "#E67E22", "#27AE60", "#8E44AD"]   # for multi-line plots

RNG = np.random.default_rng(42)


# ══════════════════════════════════════════════════════════════════════════════
# Core simulation helpers
# ══════════════════════════════════════════════════════════════════════════════

def sample_lorentzian(N, gamma=1.0):
    return RNG.standard_cauchy(N) * gamma


def order_param(theta):
    z = np.mean(np.exp(1j * theta))
    return abs(z), np.angle(z)


def rhs(theta, omega, K):
    """Vectorised RHS using complex sum — O(N) not O(N²)."""
    S = np.sum(np.exp(1j * theta))
    return omega + (K / len(theta)) * np.imag(np.exp(-1j * theta) * S)


def simulate(N, K, T=50, dt=0.05, gamma=1.0):
    """
    Integrate the Kuramoto model with RK4.
    Returns (t, theta_history, r_array, psi_array).
    """
    omega  = sample_lorentzian(N, gamma)
    theta  = RNG.uniform(0, 2 * np.pi, N)
    steps  = int(T / dt)
    t      = np.arange(steps + 1) * dt

    theta_hist = np.empty((steps + 1, N))
    r_arr      = np.empty(steps + 1)
    psi_arr    = np.empty(steps + 1)

    theta_hist[0]    = theta
    r_arr[0], psi_arr[0] = order_param(theta)

    for s in range(steps):
        k1 = rhs(theta,               omega, K)
        k2 = rhs(theta + 0.5*dt*k1,   omega, K)
        k3 = rhs(theta + 0.5*dt*k2,   omega, K)
        k4 = rhs(theta +    dt*k3,     omega, K)
        theta = theta + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)

        theta_hist[s+1]          = theta
        r_arr[s+1], psi_arr[s+1] = order_param(theta)

    return t, theta_hist, r_arr, psi_arr


def steady_r(N, K, T=40, gamma=1.0, tail=0.2):
    """Return time-averaged r over the last `tail` fraction of the run."""
    _, _, r_arr, _ = simulate(N, K, T, gamma=gamma)
    cut = int((1 - tail) * len(r_arr))
    return float(np.mean(r_arr[cut:]))


# ══════════════════════════════════════════════════════════════════════════════
# Figure 1 – Phase portraits
# ══════════════════════════════════════════════════════════════════════════════

def make_fig1():
    print("[1/4] Phase portraits ...")
    K_vals = [0.5, 1.5, 2.5, 5.0]
    N, T   = 100, 50

    fig, axes = plt.subplots(1, 4, figsize=(14, 3.8))
    fig.suptitle(r"Oscillator phases on the unit circle   ($N=100$, $\gamma=1$, $K_c=2$)",
                 fontsize=13, fontweight="bold", y=1.02)

    for ax, K in zip(axes, K_vals):
        _, th_hist, r_arr, psi_arr = simulate(N, K, T)
        th    = th_hist[-1]
        r_val = r_arr[-1]
        psi   = psi_arr[-1]

        # unit circle
        ax.add_patch(Circle((0,0), 1, fill=False, color=GRAY, lw=1.0, ls="--"))

        # oscillator dots
        ax.scatter(np.cos(th), np.sin(th),
                   s=20, color=TEAL, alpha=0.75,
                   edgecolors="white", linewidths=0.3, zorder=3)

        # order parameter arrow
        ax.annotate(
            "", xy=(r_val*np.cos(psi), r_val*np.sin(psi)),
            xytext=(0, 0),
            arrowprops=dict(arrowstyle="-|>", color=RED, lw=2.2,
                            mutation_scale=14)
        )

        # dashed r-circle
        ax.add_patch(Circle((0,0), r_val, fill=False,
                             color=RED, lw=0.7, ls=":", alpha=0.5))

        ax.set_xlim(-1.35, 1.35); ax.set_ylim(-1.35, 1.35)
        ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)

        label = r"$K < K_c$" if K < 2.0 else (r"$K > K_c$" if K > 2.0 else r"$K \approx K_c$")
        ax.set_title(f"$K = {K}$   $r = {r_val:.2f}$\n{label}", fontsize=11)

    plt.tight_layout()
    plt.savefig("fig1_phase_snapshots.png", bbox_inches="tight", dpi=150)
    plt.close()
    print("   saved fig1_phase_snapshots.png")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 2 – r(t) time series
# ══════════════════════════════════════════════════════════════════════════════

def make_fig2():
    print("[2/4] Time series r(t) ...")
    K_vals = [0.5, 1.5, 2.0, 3.0, 5.0]
    N, T, Kc = 100, 60, 2.0

    fig, ax = plt.subplots(figsize=(9, 4.8))

    line_colors = ["#2980B9", "#7F8C8D", "#E67E22", "#27AE60", "#C0392B"]

    for K, col in zip(K_vals, line_colors):
        t, _, r_arr, _ = simulate(N, K, T)
        ax.plot(t, r_arr, color=col, lw=1.7, alpha=0.9, label=f"$K = {K}$")

        # analytical steady-state for K > Kc
        if K > Kc:
            r_th = np.sqrt(1 - Kc / K)
            ax.hlines(r_th, 0, T, colors=col, lw=1.1, ls=":", alpha=0.65)

    # simple annotations
    ax.text(1.5, 0.07, r"$K < K_c$: incoherent",
            color="#2980B9", fontsize=10)
    ax.text(36, 0.93, r"$K \gg K_c$: synchronized",
            color="#C0392B", fontsize=10, ha="right")
    ax.text(61, 0.60, "dotted = $r^*$ (theory)",
            color=GRAY, fontsize=9, style="italic", ha="right")

    ax.set_xlabel("Time  $t$")
    ax.set_ylabel("Order parameter  $r(t)$")
    ax.set_title("Time evolution of the order parameter", fontweight="bold")
    ax.set_ylim(-0.03, 1.05)
    ax.set_xlim(0, T)
    ax.legend(loc="center right", framealpha=0.9)

    plt.tight_layout()
    plt.savefig("fig2_r_timeseries.png", bbox_inches="tight", dpi=150)
    plt.close()
    print("   saved fig2_r_timeseries.png")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 3 – r vs K (synchronization transition)
# ══════════════════════════════════════════════════════════════════════════════

def make_fig3():
    print("[3/4] r vs K transition ...")
    N, gamma, Kc = 100, 1.0, 2.0

    K_pts  = np.arange(0.0, 6.2, 0.3)
    r_num  = np.array([steady_r(N, K, gamma=gamma) for K in K_pts])

    K_th   = np.linspace(Kc, 6.1, 300)
    r_th   = np.sqrt(1 - Kc / K_th)

    fig, ax = plt.subplots(figsize=(8, 5))

    # shaded regions
    ax.axvspan(0,   Kc,  alpha=0.06, color=TEAL)
    ax.axvspan(Kc, 6.2,  alpha=0.06, color=RED)

    ax.scatter(K_pts, r_num, s=48, color=TEAL, zorder=5,
               label=f"Simulation ($N = {N}$)", alpha=0.85)
    ax.plot(K_th, r_th, color=RED, lw=2.5,
            label=r"Theory: $r = \sqrt{1 - K_c/K}$")
    ax.axvline(Kc, color=GRAY, lw=1.4, ls="--", label=f"$K_c = {Kc}$")

    ax.text(0.25, 0.88, "Incoherent\nphase", color=TEAL, fontsize=10)
    ax.text(4.8,  0.08, "Synchronized\nphase", color=RED,  fontsize=10, ha="center")

    ax.set_xlabel("Coupling strength  $K$")
    ax.set_ylabel(r"Steady-state order parameter  $\langle r \rangle$")
    ax.set_title("Synchronization transition in the Kuramoto model",
                 fontweight="bold")
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlim(-0.1,  6.2)
    ax.legend(loc="upper left", framealpha=0.9)

    plt.tight_layout()
    plt.savefig("fig3_r_vs_K.png", bbox_inches="tight", dpi=150)
    plt.close()
    print("   saved fig3_r_vs_K.png")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 4 – Finite-size effects
# ══════════════════════════════════════════════════════════════════════════════

def make_fig4():
    print("[4/4] Finite-size effects ...")
    N_vals  = [20, 50, 100, 300]
    gamma   = 1.0
    Kc      = 2.0
    K_pts   = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0])

    K_th  = np.linspace(0, 6.1, 400)
    r_th  = np.where(K_th > Kc, np.sqrt(1 - Kc / np.maximum(K_th, 1e-9)), 0.0)

    fig, ax = plt.subplots(figsize=(8, 5))

    for N, col in zip(N_vals, COLORS):
        r_num = np.array([steady_r(N, K, gamma=gamma) for K in K_pts])
        ax.plot(K_pts, r_num, "o-", color=col, lw=1.8,
                ms=5.5, alpha=0.85, label=f"$N = {N}$")

    ax.plot(K_th, r_th, color=DARK, lw=2.5, ls="--",
            label=r"Theory ($N \to \infty$)", zorder=10)
    ax.axvline(Kc, color=GRAY, lw=1.2, ls=":", alpha=0.8)
    ax.text(Kc + 0.06, 1.04, r"$K_c = 2$", color=GRAY, fontsize=10)

    ax.set_xlabel("Coupling strength  $K$")
    ax.set_ylabel(r"$\langle r \rangle$")
    ax.set_title("Finite-size convergence to the thermodynamic limit",
                 fontweight="bold")
    ax.set_ylim(-0.05, 1.12)
    ax.set_xlim(-0.1,  6.2)
    ax.legend(loc="upper left", framealpha=0.9)

    plt.tight_layout()
    plt.savefig("fig4_finite_size.png", bbox_inches="tight", dpi=150)
    plt.close()
    print("   saved fig4_finite_size.png")


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("=" * 55)
    print("Kuramoto Model Simulation  (gamma=1, Kc=2)")
    print("=" * 55)
    make_fig1()
    make_fig2()
    make_fig3()
    make_fig4()
    print("=" * 55)
    print("Done. Place the 4 PNGs next to kuramoto_clean.tex")
    print("and compile:  pdflatex kuramoto_clean.tex  (run twice)")
    print("=" * 55)
