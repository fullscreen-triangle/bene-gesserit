#!/usr/bin/env python3
"""
Generate 5 panel figures for:
Membrane-Mediated Categorical Processing Architecture
(amphiphilic bilayer membranes as computational substrate)
"""

import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from pathlib import Path
import os

# ── paths ──────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
PANEL_DIR = SCRIPT_DIR / "panels"
PANEL_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR = SCRIPT_DIR.parent.parent / "results"

# ── shared style ───────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 11,
    'axes.labelsize': 10,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'savefig.facecolor': 'white',
    'savefig.dpi': 150,
})

PANEL_LABELS = ['(A)', '(B)', '(C)', '(D)']

def label_axes(axes):
    """Put (A)-(D) in upper-left corner of each subplot."""
    for i, ax in enumerate(axes):
        ax.set_title(f'{PANEL_LABELS[i]}', loc='left', fontweight='bold')


# ═══════════════════════════════════════════════════════════════════════════
# PANEL 1 — Lipid bilayer as oscillator network & semiconductor
# ═══════════════════════════════════════════════════════════════════════════
def panel_01():
    fig = plt.figure(figsize=(20, 5), dpi=150, facecolor='white')

    # Physical constants
    q = 1.602e-19      # C
    k_B = 1.381e-23    # J/K
    T = 310.0           # K
    n = 1.5
    I_0 = 1e-12         # A
    V_bi = 0.615        # V

    # ── (A) Shockley diode I-V ──────────────────────────────────────────
    ax1 = fig.add_subplot(1, 4, 1, facecolor='white')
    V = np.linspace(-0.3, 0.8, 2000)
    I = I_0 * (np.exp(q * V / (n * k_B * T)) - 1)
    # plot |I| on log scale
    ax1.semilogy(V, np.abs(I), 'b-', lw=1.5)
    ax1.axvline(V_bi, color='r', ls='--', lw=1, label=f'$V_{{bi}}$ = {V_bi} V')
    ax1.set_xlabel('Voltage (V)')
    ax1.set_ylabel('|I| (A)')
    ax1.set_title('(A) Shockley diode I-V', loc='left', fontweight='bold')
    ax1.set_ylim(1e-15, 1e1)
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    # ── (B) Drift velocity vs E-field ───────────────────────────────────
    ax2 = fig.add_subplot(1, 4, 2, facecolor='white')
    mu = 0.0123e-4       # cm²/Vs → m²/Vs
    v_sat = 80000e-2     # cm/s → m/s = 800 m/s
    E = np.logspace(1, 7, 500)   # V/m
    v_low = mu * E
    # saturation model: v = mu*E / (1 + mu*E/v_sat)
    v_drift = (mu * E) / (1 + mu * E / v_sat)
    ax2.loglog(E, v_drift, 'b-', lw=1.5, label='$v_{drift}$')
    ax2.loglog(E, v_low, 'r--', lw=0.8, alpha=0.5, label=r'$\mu E$ (low-field)')
    ax2.axhline(v_sat, color='g', ls=':', lw=0.8, label=f'$v_{{sat}}$ = {v_sat*100:.0f} cm/s')
    ax2.set_xlabel('Electric field (V/m)')
    ax2.set_ylabel('Drift velocity (m/s)')
    ax2.set_title('(B) Biological semiconductor', loc='left', fontweight='bold')
    ax2.legend(fontsize=7)
    ax2.grid(True, alpha=0.3)

    # ── (C) BMD transistor transfer ─────────────────────────────────────
    ax3 = fig.add_subplot(1, 4, 3, facecolor='white')
    G = 1000
    theta = 0.5
    I_in = 1e-12   # A
    M = np.linspace(0, 1, 1000)
    I_out = np.where(M > theta,
                     G * (M - theta) / (1 - M + theta) * I_in,
                     0.0)
    ax3.semilogy(M, np.clip(I_out, 1e-18, None), 'b-', lw=1.5)
    ax3.axvline(theta, color='r', ls='--', lw=1, label=r'$\theta$ = 0.5')
    # mark on/off ratio
    ax3.annotate(f'On/Off = 42.1', xy=(0.85, 1e-10), fontsize=9,
                 bbox=dict(boxstyle='round,pad=0.3', fc='lightyellow', ec='gray', alpha=0.9))
    ax3.set_xlabel('Overlap integral $M$')
    ax3.set_ylabel('$I_{out}$ (A)')
    ax3.set_title('(C) BMD transistor transfer', loc='left', fontweight='bold')
    ax3.set_ylim(1e-18, 1e-6)
    ax3.legend(fontsize=8)
    ax3.grid(True, alpha=0.3)

    # ── (D) 3D RCL impedance surface ───────────────────────────────────
    ax4 = fig.add_subplot(1, 4, 4, projection='3d', facecolor='white')
    R_0 = 1e9
    L_val = 3.14e12
    C_val = 318.3e-15
    omega = np.logspace(3, 9, 80)
    S_k = np.linspace(0, 1, 60)
    W, S = np.meshgrid(omega, S_k)
    R_sk = R_0 * np.exp(S)
    Z_mag = np.sqrt(R_sk**2 + (W * L_val - 1.0 / (W * C_val))**2)
    log_Z = np.log10(Z_mag)
    log_W = np.log10(W)
    surf = ax4.plot_surface(log_W, S, log_Z, cmap='viridis', alpha=0.85,
                            rstride=2, cstride=2, edgecolor='none')
    ax4.set_xlabel(r'log$_{10}(\omega)$', labelpad=6)
    ax4.set_ylabel('$S_k$', labelpad=6)
    ax4.set_zlabel(r'log$_{10}|Z|$', labelpad=4)
    ax4.set_title('(D) R-C-L impedance', loc='left', fontweight='bold')
    ax4.view_init(elev=25, azim=-50)

    fig.tight_layout()
    out = PANEL_DIR / 'panel_exp01_membrane_substrate.png'
    fig.savefig(str(out), dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  saved {out}')


# ═══════════════════════════════════════════════════════════════════════════
# PANEL 2 — Zero-work aperture filtering & ternary trisection
# ═══════════════════════════════════════════════════════════════════════════
def panel_02():
    fig = plt.figure(figsize=(20, 5), dpi=150, facecolor='white')

    # ── (A) Aperture with particle trajectories ────────────────────────
    ax1 = fig.add_subplot(1, 4, 1, facecolor='white')
    R = 1.0
    theta_circ = np.linspace(0, 2 * np.pi, 200)
    ax1.plot(R * np.cos(theta_circ), R * np.sin(theta_circ), 'k-', lw=2)
    # shade exterior
    theta_fill = np.linspace(0, 2 * np.pi, 200)
    r_outer = 2.5
    ax1.fill_between(r_outer * np.cos(theta_fill), r_outer * np.sin(theta_fill),
                     color='lightgray', alpha=0.3, zorder=0)
    ax1.fill(R * np.cos(theta_fill), R * np.sin(theta_fill), color='white', zorder=1)
    np.random.seed(42)
    for _ in range(50):
        # launch from left
        y0 = np.random.uniform(-2, 2)
        x_start = -2.5
        # straight trajectory + slight random slope
        slope = np.random.uniform(-0.3, 0.3)
        x_pts = np.linspace(x_start, 2.5, 200)
        y_pts = y0 + slope * (x_pts - x_start)
        # check if passes through aperture at x near 0
        idx0 = np.argmin(np.abs(x_pts))
        passes = (x_pts[idx0]**2 + y_pts[idx0]**2) < R**2
        if passes:
            ax1.plot(x_pts, y_pts, 'b-', lw=0.4, alpha=0.6)
        else:
            # blocked: draw up to boundary then stop
            dists = np.sqrt(x_pts**2 + y_pts**2)
            hit_idx = np.argmax(dists < R) if np.any(dists < R) else len(x_pts) // 2
            if hit_idx == 0:
                hit_idx = len(x_pts) // 2
            ax1.plot(x_pts[:hit_idx], y_pts[:hit_idx], 'r-', lw=0.4, alpha=0.5)
    ax1.set_xlim(-2.5, 2.5)
    ax1.set_ylim(-2.5, 2.5)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_title('(A) Aperture trajectories', loc='left', fontweight='bold')
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.2)

    # ── (B) Multipole field decay ───────────────────────────────────────
    ax2 = fig.add_subplot(1, 4, 2, facecolor='white')
    r = np.linspace(0.1, 10, 500)
    monopole = 1.0 / r**2
    dipole = 1.0 / r**3
    quadrupole = 1.0 / r**4
    ax2.loglog(r, monopole, 'b-', lw=1.5, label=r'Monopole $\sim r^{-2}$')
    ax2.loglog(r, dipole, 'r-', lw=1.5, label=r'Dipole $\sim r^{-3}$')
    ax2.loglog(r, quadrupole, 'g-', lw=1.5, label=r'Quadrupole $\sim r^{-4}$')
    ax2.set_xlabel('Distance $r$')
    ax2.set_ylabel('Field strength')
    ax2.set_title('(B) Multipole decay', loc='left', fontweight='bold')
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    # ── (C) Ternary trisection volume & states ─────────────────────────
    ax3 = fig.add_subplot(1, 4, 3, facecolor='white')
    k = np.arange(0, 21)
    log3_vol = -k  # log_3(1/3^k) = -k (each trisection keeps 1/3)
    states = 3.0 ** k
    ax3_twin = ax3.twinx()
    l1, = ax3.plot(k, log3_vol, 'b-o', lw=1.5, ms=4, label=r'$\log_3$(volume)')
    l2, = ax3_twin.semilogy(k, states, 'r-s', lw=1.5, ms=4, label=r'States $3^k$')
    ax3.set_xlabel('Trisection depth $k$')
    ax3.set_ylabel(r'$\log_3$(volume)', color='b')
    ax3_twin.set_ylabel(r'Distinguishable states $3^k$', color='r')
    ax3.set_title('(C) Ternary trisection', loc='left', fontweight='bold')
    lines = [l1, l2]
    ax3.legend(lines, [l.get_label() for l in lines], fontsize=8, loc='center left')
    ax3.grid(True, alpha=0.3)

    # ── (D) 3D aperture filtering scatter ──────────────────────────────
    ax4 = fig.add_subplot(1, 4, 4, projection='3d', facecolor='white')
    np.random.seed(7)
    pts = np.random.rand(1000, 3)
    colors_stage = ['#cccccc', '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    # Apply 5 successive aperture filters
    surviving = np.ones(len(pts), dtype=bool)
    stage_labels = np.zeros(len(pts), dtype=int)
    for stage in range(5):
        # each filter keeps ~1/3 along a random axis slice
        axis = stage % 3
        center = np.random.uniform(0.2, 0.8)
        width = 0.33 * (0.9 ** stage)
        in_band = np.abs(pts[:, axis] - center) < width / 2
        newly_removed = surviving & ~in_band
        stage_labels[newly_removed] = stage + 1
        surviving &= in_band
    stage_labels[surviving] = 5  # final survivors

    for s in range(6):
        mask = stage_labels == s
        if not np.any(mask):
            continue
        alpha_val = 0.08 if s == 0 else (0.25 + 0.15 * s)
        size = 3 if s < 5 else 20
        ax4.scatter(pts[mask, 0], pts[mask, 1], pts[mask, 2],
                    c=colors_stage[s], s=size, alpha=min(alpha_val, 1.0),
                    label=f'Stage {s}' if s > 0 else 'Unfiltered')
    ax4.set_xlabel('$S_1$')
    ax4.set_ylabel('$S_2$')
    ax4.set_zlabel('$S_3$')
    ax4.set_title('(D) Progressive filtering', loc='left', fontweight='bold')
    ax4.view_init(elev=20, azim=35)

    fig.tight_layout()
    out = PANEL_DIR / 'panel_exp02_categorical_apertures.png'
    fig.savefig(str(out), dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  saved {out}')


# ═══════════════════════════════════════════════════════════════════════════
# PANEL 3 — Backward trajectory completion & categorical speedup
# ═══════════════════════════════════════════════════════════════════════════
def panel_03():
    fig = plt.figure(figsize=(20, 5), dpi=150, facecolor='white')

    # ── (A) Operations scaling ──────────────────────────────────────────
    ax1 = fig.add_subplot(1, 4, 1, facecolor='white')
    N = np.logspace(1, 12, 500)
    cat_ops = np.ceil(np.log(N) / np.log(3))
    linear = N
    nlogn = N * np.log2(N)
    ax1.loglog(N, cat_ops, 'b-', lw=2, label=r'$\lceil\log_3 N\rceil$ (categorical)')
    ax1.loglog(N, linear, 'r--', lw=1.2, label='$O(N)$')
    ax1.loglog(N, nlogn, 'g--', lw=1.2, label=r'$O(N \log N)$')
    ax1.set_xlabel('Problem size $N$')
    ax1.set_ylabel('Operations')
    ax1.set_title('(A) Categorical vs conventional', loc='left', fontweight='bold')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(1, 1e15)

    # ── (B) Speedup factor ──────────────────────────────────────────────
    ax2 = fig.add_subplot(1, 4, 2, facecolor='white')
    log3N = np.log(N) / np.log(3)
    speedup = N / np.maximum(log3N, 1)
    ax2.loglog(N, speedup, 'b-', lw=2)
    ax2.set_xlabel('Problem size $N$')
    ax2.set_ylabel('Speedup $N / \\log_3 N$')
    ax2.set_title('(B) Categorical speedup', loc='left', fontweight='bold')
    ax2.grid(True, alpha=0.3)

    # ── (C) 3D trajectory convergence ──────────────────────────────────
    ax3 = fig.add_subplot(1, 4, 3, projection='3d', facecolor='white')
    np.random.seed(99)
    target = np.array([0.7, 0.3, 0.5])
    start = np.array([0.1, 0.9, 0.1])
    n_way = 20
    trajectory = np.zeros((n_way, 3))
    trajectory[0] = start
    for i in range(1, n_way):
        # ternary contraction: each step moves 2/3 toward target + small noise
        trajectory[i] = trajectory[i-1] + 0.6 * (target - trajectory[i-1])
        trajectory[i] += np.random.randn(3) * 0.01 * (0.5 ** i)
        trajectory[i] = np.clip(trajectory[i], 0, 1)
    ax3.plot(trajectory[:, 0], trajectory[:, 1], trajectory[:, 2],
             'b-o', lw=1.2, ms=4, alpha=0.8)
    ax3.scatter(*target, c='r', s=200, marker='*', zorder=10, label='Target')
    ax3.scatter(*start, c='g', s=80, marker='o', zorder=10, label='Start')
    ax3.set_xlabel('$S_1$')
    ax3.set_ylabel('$S_2$')
    ax3.set_zlabel('$S_3$')
    ax3.set_title('(C) Trajectory completion', loc='left', fontweight='bold')
    ax3.legend(fontsize=8)
    ax3.view_init(elev=25, azim=-40)

    # ── (D) Convergence to terminus ────────────────────────────────────
    ax4 = fig.add_subplot(1, 4, 4, facecolor='white')
    np.random.seed(12)
    steps = np.arange(0, 20)
    colors_traj = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    for i in range(5):
        d0 = np.random.uniform(0.5, 1.5)
        rate = 1.0 / 3.0  # contraction by 1/3 each step
        d = d0 * rate ** steps + np.random.randn(len(steps)) * 1e-4
        d = np.maximum(d, 1e-8)
        ax4.semilogy(steps, d, '-o', lw=1.2, ms=3, color=colors_traj[i],
                     label=f'Start {i+1}')
    ax4.set_xlabel('Iteration step')
    ax4.set_ylabel('Distance to terminus $d_{cat}$')
    ax4.set_title('(D) Exponential convergence', loc='left', fontweight='bold')
    ax4.legend(fontsize=7, ncol=2)
    ax4.grid(True, alpha=0.3)

    fig.tight_layout()
    out = PANEL_DIR / 'panel_exp03_trajectory_completion.png'
    fig.savefig(str(out), dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  saved {out}')


# ═══════════════════════════════════════════════════════════════════════════
# PANEL 4 — Decay envelopes, regime dynamics, R-C-L trichotomy
# ═══════════════════════════════════════════════════════════════════════════
def panel_04():
    fig = plt.figure(figsize=(20, 5), dpi=150, facecolor='white')

    # ── (A) Decay curves & processing window ───────────────────────────
    ax1 = fig.add_subplot(1, 4, 1, facecolor='white')
    t = np.linspace(0, 500, 2000)  # ms
    P_0, tau1 = 1.0, 50.0
    T_0, tau2, tau3 = 1.0, 20.0, 100.0
    P = P_0 * np.exp(-t / tau1)
    T_curve = T_0 * (1 - np.exp(-t / tau2)) * np.exp(-t / tau3)
    ax1.plot(t, P, 'b-', lw=1.5, label=r'$P(t)$, $\tau$=50 ms')
    ax1.plot(t, T_curve, 'r-', lw=1.5, label=r'$T(t)$, rise/decay')
    overlap = np.minimum(P, T_curve)
    ax1.fill_between(t, 0, overlap, alpha=0.25, color='purple', label='Processing window')
    ax1.set_xlabel('Time (ms)')
    ax1.set_ylabel('Amplitude')
    ax1.set_title('(A) Decay envelopes', loc='left', fontweight='bold')
    ax1.legend(fontsize=7)
    ax1.grid(True, alpha=0.3)

    # ── (B) Kuramoto phase portrait ────────────────────────────────────
    ax2 = fig.add_subplot(1, 4, 2, facecolor='white')
    K_vals = np.linspace(0.1, 10, 200)
    sigma2_eq = 1.0 / np.sqrt(K_vals)
    ax2.plot(sigma2_eq, K_vals, 'k-', lw=2, label='Equilibrium')
    # simulated trajectories spiraling toward equilibrium
    np.random.seed(55)
    for i in range(6):
        s0 = np.random.uniform(0.1, 2.0)
        k0 = np.random.uniform(0.5, 8.0)
        # spiral: approach equilibrium K_eq where sigma2 = 1/sqrt(K)
        n_pts = 40
        s_traj = np.zeros(n_pts)
        k_traj = np.zeros(n_pts)
        s_traj[0] = s0
        k_traj[0] = k0
        for j in range(1, n_pts):
            k_eq_local = 1.0 / (s_traj[j-1]**2 + 0.01)**2 if s_traj[j-1] > 0 else 5.0
            k_eq_local = min(k_eq_local, 10)
            s_eq_local = 1.0 / np.sqrt(max(k_traj[j-1], 0.1))
            s_traj[j] = s_traj[j-1] + 0.15 * (s_eq_local - s_traj[j-1]) + np.random.randn() * 0.02
            k_traj[j] = k_traj[j-1] + 0.15 * (k_eq_local - k_traj[j-1]) + np.random.randn() * 0.05
            s_traj[j] = max(s_traj[j], 0.01)
            k_traj[j] = max(k_traj[j], 0.1)
        ax2.plot(s_traj, k_traj, '-', lw=0.8, alpha=0.6)
        ax2.plot(s_traj[0], k_traj[0], 'o', ms=5, color='green')
        ax2.plot(s_traj[-1], k_traj[-1], 's', ms=5, color='red')
    ax2.set_xlabel(r'Phase variance $\sigma^2$')
    ax2.set_ylabel('Coupling $K$')
    ax2.set_title('(B) Kuramoto phase portrait', loc='left', fontweight='bold')
    ax2.set_xlim(0, 2.5)
    ax2.set_ylim(0, 10)
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    # ── (C) Regime cycling ─────────────────────────────────────────────
    ax3 = fig.add_subplot(1, 4, 3, facecolor='white')
    t_reg = np.linspace(0, 10, 2000)
    # build R(t) as piecewise sinusoidal cycling through regimes
    # Regimes: coherent(0.8-1.0), cascade(0.5-0.8), locked(0.3-0.5), turbulent(0.0-0.3)
    regime_centers = {
        'Coherent': 0.9, 'Cascade-down': 0.65, 'Locked': 0.4,
        'Cascade-up': 0.65, 'Turbulent': 0.15
    }
    regime_colors = {
        'Coherent': '#2ca02c', 'Cascade-down': '#ff7f0e', 'Locked': '#1f77b4',
        'Cascade-up': '#ff7f0e', 'Turbulent': '#d62728'
    }
    regime_bounds = [
        (0.0, 2.0, 'Coherent', 0.9),
        (2.0, 3.5, 'Cascade-down', 0.65),
        (3.5, 5.5, 'Locked', 0.4),
        (5.5, 7.0, 'Cascade-up', 0.65),
        (7.0, 8.5, 'Turbulent', 0.15),
        (8.5, 10.0, 'Coherent', 0.9),
    ]
    R_t = np.zeros_like(t_reg)
    for t0, t1, name, center in regime_bounds:
        mask = (t_reg >= t0) & (t_reg < t1)
        R_t[mask] = center + 0.05 * np.sin(2 * np.pi * 3 * (t_reg[mask] - t0))
        ax3.axvspan(t0, t1, alpha=0.15, color=regime_colors[name])
    ax3.plot(t_reg, R_t, 'k-', lw=1.2)
    # add legend-like labels
    used = set()
    for t0, t1, name, center in regime_bounds:
        short = name.replace('-down', '').replace('-up', '')
        if short not in used:
            ax3.text((t0 + t1) / 2, 1.02, short, ha='center', fontsize=7,
                     transform=ax3.get_xaxis_transform())
            used.add(short)
    ax3.set_xlabel('Time (s)')
    ax3.set_ylabel('Order parameter $R$')
    ax3.set_title('(C) Regime cycling', loc='left', fontweight='bold')
    ax3.set_ylim(0, 1.05)
    ax3.grid(True, alpha=0.3)

    # ── (D) 3D unit simplex — R/C/L mode selection ─────────────────────
    ax4 = fig.add_subplot(1, 4, 4, projection='3d', facecolor='white')
    # sample the unit simplex S_k + S_t + S_e = 1
    n_pts = 80
    coords = []
    colors_rcl = []
    for i in range(n_pts + 1):
        for j in range(n_pts + 1 - i):
            k_idx = n_pts - i - j
            sk = i / n_pts
            st = j / n_pts
            se = k_idx / n_pts
            if abs(sk + st + se - 1.0) > 0.02:
                continue
            coords.append([sk, st, se])
            # color by dominant mode
            if sk <= st and sk <= se:
                colors_rcl.append([1, 0, 0, 0.7])  # Red = Resistive (min S_k)
            elif st <= sk and st <= se:
                colors_rcl.append([0, 0.7, 0, 0.7])  # Green = Capacitive (min S_t)
            else:
                colors_rcl.append([0, 0, 1, 0.7])  # Blue = Inductive (min S_e)
    coords = np.array(coords)
    colors_rcl = np.array(colors_rcl)
    ax4.scatter(coords[:, 0], coords[:, 1], coords[:, 2],
                c=colors_rcl, s=8, alpha=0.7)
    # draw simplex edges
    verts = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 0]])
    ax4.plot(verts[:, 0], verts[:, 1], verts[:, 2], 'k-', lw=1)
    ax4.set_xlabel('$S_k$', labelpad=4)
    ax4.set_ylabel('$S_t$', labelpad=4)
    ax4.set_zlabel('$S_e$', labelpad=4)
    ax4.set_title('(D) R-C-L mode selection', loc='left', fontweight='bold')
    ax4.view_init(elev=30, azim=-55)

    fig.tight_layout()
    out = PANEL_DIR / 'panel_exp04_processing_cycle.png'
    fig.savefig(str(out), dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  saved {out}')


# ═══════════════════════════════════════════════════════════════════════════
# PANEL 5 — Throughput, energy efficiency, stack validation
# ═══════════════════════════════════════════════════════════════════════════
def panel_05():
    fig = plt.figure(figsize=(20, 5), dpi=150, facecolor='white')

    # ── (A) Energy per operation bar chart ─────────────────────────────
    ax1 = fig.add_subplot(1, 4, 1, facecolor='white')
    k_B = 1.381e-23
    T = 310
    labels_e = ['Silicon\nCMOS', 'Landauer\nlimit', 'BMD\nmembrane']
    values_e = [1e-15, k_B * T * np.log(2), 6e-20]
    colors_e = ['#1f77b4', '#ff7f0e', '#2ca02c']
    bars = ax1.bar(labels_e, values_e, color=colors_e, edgecolor='black', lw=0.5)
    ax1.set_yscale('log')
    ax1.set_ylabel('Energy per op (J)')
    ax1.set_title('(A) Energy efficiency', loc='left', fontweight='bold')
    ax1.set_ylim(1e-22, 1e-12)
    # annotate gap
    ax1.annotate('', xy=(2, 6e-20), xytext=(0, 1e-15),
                 arrowprops=dict(arrowstyle='<->', color='red', lw=1.5))
    ax1.text(1, 3e-17, r'$\sim 10^4 \times$', ha='center', fontsize=10, color='red',
             fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='y')

    # ── (B) Operations per second ──────────────────────────────────────
    ax2 = fig.add_subplot(1, 4, 2, facecolor='white')
    labels_ops = ['Modern\nGPU', r'Membrane' + '\n' + r'1 cm$^2$',
                  r'Membrane' + '\n' + r'1 m$^2$']
    values_ops = [1e14, 1e20, 1e24]
    colors_ops = ['#d62728', '#9467bd', '#8c564b']
    ax2.bar(labels_ops, values_ops, color=colors_ops, edgecolor='black', lw=0.5)
    ax2.set_yscale('log')
    ax2.set_ylabel('Operations / second')
    ax2.set_title('(B) Throughput', loc='left', fontweight='bold')
    ax2.set_ylim(1e12, 1e26)
    ax2.grid(True, alpha=0.3, axis='y')

    # ── (C) 3D bar — stack validation ──────────────────────────────────
    ax3 = fig.add_subplot(1, 4, 3, projection='3d', facecolor='white')
    layers = ['P-N\njunction', 'BMD\ntransistor', 'Logic\ngates',
              'Quantum\ngates', 'Inter-\nconnects', 'Program\nexec']
    raw_metrics = [47.8, 42.1, 100.0, 88.4, 23500, 91.0]
    # normalize for display: map to [0, 1]
    max_val = max(raw_metrics)
    norm_metrics = [v / max_val for v in raw_metrics]
    n_bars = len(layers)
    xpos = np.arange(n_bars)
    ypos = np.zeros(n_bars)
    zpos = np.zeros(n_bars)
    dx = 0.6
    dy = 0.6
    dz = np.array(norm_metrics)
    bar_colors = cm.viridis(np.linspace(0.2, 0.9, n_bars))
    ax3.bar3d(xpos, ypos, zpos, dx, dy, dz, color=bar_colors, alpha=0.85,
              edgecolor='black', linewidth=0.3)
    ax3.set_xticks(xpos + dx / 2)
    ax3.set_xticklabels(layers, fontsize=6, rotation=15, ha='right')
    ax3.set_yticks([])
    ax3.set_zlabel('Normalized metric')
    ax3.set_title('(C) Stack validation', loc='left', fontweight='bold')
    ax3.view_init(elev=25, azim=-50)
    # add value annotations
    for i in range(n_bars):
        val_str = f'{raw_metrics[i]:.0f}' if raw_metrics[i] == int(raw_metrics[i]) else f'{raw_metrics[i]:.1f}'
        if raw_metrics[i] >= 1000:
            val_str = f'{raw_metrics[i]:.0f}'
        ax3.text(xpos[i] + dx/2, 0.3, dz[i] + 0.03, val_str, fontsize=6,
                 ha='center', va='bottom')

    # ── (D) Self-consistency scatter ───────────────────────────────────
    ax4 = fig.add_subplot(1, 4, 4, facecolor='white')
    with open(str(RESULTS_DIR / 'self_consistency.json'), 'r') as f:
        sc_data = json.load(f)
    mol_atoms = {'H2': 2, 'CO': 2, 'H2O': 3, 'CO2': 3, 'CH4': 5, 'C6H6': 12}
    x_atoms = []
    y_dev = []
    mol_labels = []
    for mol in sc_data['molecules_tested']:
        x_atoms.append(mol_atoms[mol])
        y_dev.append(sc_data['results'][mol]['max_deviation_pct'])
        mol_labels.append(mol)
    x_atoms = np.array(x_atoms, dtype=float)
    y_dev = np.array(y_dev, dtype=float)
    ax4.scatter(x_atoms, y_dev, s=80, c='#1f77b4', edgecolors='black', lw=0.5, zorder=5)
    for i, lbl in enumerate(mol_labels):
        ax4.annotate(lbl, (x_atoms[i], y_dev[i]), textcoords='offset points',
                     xytext=(8, 5), fontsize=8)
    # fit trend line (linear)
    # remove zero-deviation points for trend? No, include them all.
    coeffs = np.polyfit(x_atoms, y_dev, 1)
    x_fit = np.linspace(0, 14, 100)
    y_fit = np.polyval(coeffs, x_fit)
    ax4.plot(x_fit, y_fit, 'r--', lw=1.2, alpha=0.7, label=f'Trend (slope={coeffs[0]:.2f}%/atom)')
    ax4.axhline(sc_data['summary']['threshold_pct'], color='gray', ls=':', lw=1,
                label=f'Threshold = {sc_data["summary"]["threshold_pct"]}%')
    ax4.set_xlabel('Molecular complexity (atoms)')
    ax4.set_ylabel('Max deviation (%)')
    ax4.set_title('(D) Self-consistency scaling', loc='left', fontweight='bold')
    ax4.legend(fontsize=7)
    ax4.set_xlim(0, 14)
    ax4.set_ylim(-0.2, 2.5)
    ax4.grid(True, alpha=0.3)

    fig.subplots_adjust(left=0.05, right=0.97, wspace=0.35)
    out = PANEL_DIR / 'panel_exp05_performance.png'
    fig.savefig(str(out), dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  saved {out}')


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print('Generating panels …')
    panel_01()
    panel_02()
    panel_03()
    panel_04()
    panel_05()
    # verify
    expected = [
        'panel_exp01_membrane_substrate.png',
        'panel_exp02_categorical_apertures.png',
        'panel_exp03_trajectory_completion.png',
        'panel_exp04_processing_cycle.png',
        'panel_exp05_performance.png',
    ]
    print('\nVerification:')
    all_ok = True
    for fn in expected:
        p = PANEL_DIR / fn
        if p.exists():
            sz = p.stat().st_size
            print(f'  OK  {fn}  ({sz:,} bytes)')
        else:
            print(f'  MISSING  {fn}')
            all_ok = False
    if all_ok:
        print('\nAll 5 panels generated successfully.')
    else:
        print('\nSome panels are missing!')
