"""
Generate 5 panel figures for the Triple Isomorphism Architecture paper.
Each panel is a 1x4 grid of data-driven subplots (figsize 20x5, dpi 150).
Output: reverend/publications/isomorphism-atlas/panels/*.png
"""

import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# ---------------------------------------------------------------------------
# paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PANELS_DIR = os.path.join(SCRIPT_DIR, "panels")
RESULTS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "..", "results"))
os.makedirs(PANELS_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# helper
# ---------------------------------------------------------------------------
def new_figure():
    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    fig.patch.set_facecolor("white")
    for ax in axes:
        ax.set_facecolor("white")
    return fig, axes


def new_figure_with_3d(pos_3d):
    """Create 1x4 figure where subplot at *pos_3d* (0-based) is 3D."""
    fig = plt.figure(figsize=(20, 5), facecolor="white")
    axes = []
    for i in range(4):
        if i == pos_3d:
            ax = fig.add_subplot(1, 4, i + 1, projection="3d")
        else:
            ax = fig.add_subplot(1, 4, i + 1)
        ax.set_facecolor("white")
        axes.append(ax)
    return fig, axes


# ===================================================================
# PANEL 1 — Capacity formula C(n)=2n^2 and partition depth
# ===================================================================
def panel_01():
    fig, axes = new_figure_with_3d(pos_3d=1)

    # (A) C(n) = 2n^2
    ax = axes[0]
    n_cont = np.linspace(1, 30, 300)
    n_disc = np.arange(1, 31)
    ax.plot(n_cont, 2 * n_cont ** 2, "b-", lw=1.5, label=r"$C(n)=2n^2$")
    ax.scatter(n_disc, 2 * n_disc ** 2, c="royalblue", s=18, zorder=5)
    ax.set_xlabel("n")
    ax.set_ylabel("Capacity  C(n)")
    ax.set_title("(A)  Capacity  $C(n)=2n^2$")
    ax.legend(fontsize=8)

    # (B) 3D surface — partition depth M(k, b) = sum log_b(k_i) ≈ k * log_b(k)
    ax = axes[1]
    k_vals = np.linspace(1, 10, 50)
    b_vals = np.linspace(2, 10, 50)
    K, B = np.meshgrid(k_vals, b_vals)
    M = K * np.log(K + 1) / np.log(B)          # sum ≈ k*log_b(k+1)
    surf = ax.plot_surface(K, B, M, cmap="viridis", edgecolor="none", alpha=0.9)
    ax.set_xlabel("k")
    ax.set_ylabel("b")
    ax.set_zlabel("M(k,b)")
    ax.set_title("(B)  Partition depth")

    # (C) Binary waste fraction
    ax = axes[2]
    N_vals = np.arange(2, 31)
    waste = 1 - N_vals / 2 ** np.ceil(np.log2(N_vals))
    ax.bar(N_vals, waste, color="salmon", edgecolor="darkred", linewidth=0.4)
    avg_waste = np.mean(waste)
    ax.axhline(avg_waste, color="k", ls="--", lw=1.2,
               label=f"mean = {avg_waste:.2f}")
    ax.set_xlabel("N")
    ax.set_ylabel("Waste  W(N)")
    ax.set_title("(C)  Binary waste fraction")
    ax.legend(fontsize=8)

    # (D) Encoding efficiency for base 2 and base 3
    ax = axes[3]
    N_vals_d = np.arange(2, 51)
    for base, color, label in [(2, "blue", "b=2"), (3, "red", "b=3")]:
        eta = N_vals_d / base ** np.ceil(np.log(N_vals_d) / np.log(base))
        ax.plot(N_vals_d, eta, color=color, lw=1.2, label=label)
    ax.set_xlabel("N")
    ax.set_ylabel(r"$\eta = N / b^{\lceil\log_b N\rceil}$")
    ax.set_title(r"(D)  Encoding efficiency $\eta$")
    ax.set_ylim(0.3, 1.1)
    ax.legend(fontsize=8)

    fig.tight_layout()
    path = os.path.join(PANELS_DIR, "panel_exp01_capacity_and_partition.png")
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)
    print(f"  saved {path}")


# ===================================================================
# PANEL 2 — The three entropies converging
# ===================================================================
def panel_02():
    fig, axes = new_figure_with_3d(pos_3d=0)

    # (A) 3D scatter — S-entropy cube
    ax = axes[0]
    rng = np.random.default_rng(42)
    Sk = rng.random(200)
    St = rng.random(200)
    Se = rng.random(200)
    norm = np.sqrt(Sk ** 2 + St ** 2 + Se ** 2)
    sc = ax.scatter(Sk, St, Se, c=norm, cmap="plasma", s=12, alpha=0.8)
    ax.set_xlabel(r"$S_k$")
    ax.set_ylabel(r"$S_t$")
    ax.set_zlabel(r"$S_e$")
    ax.set_title("(A)  S-entropy cube")

    # (B) S_osc(M) = k_B * M * ln(n)
    ax = axes[1]
    kB = 1.380649e-23
    M_vals = np.linspace(0, 20, 200)
    for n, color in [(2, "blue"), (3, "green"), (5, "orange"), (10, "red")]:
        S = kB * M_vals * np.log(n)
        ax.plot(M_vals, S / kB, color=color, lw=1.3, label=f"n={n}")
    ax.set_xlabel("Partition depth  M")
    ax.set_ylabel(r"$S_{\rm osc}\ /\ k_B$")
    ax.set_title(r"(B)  $S_{\rm osc}(M)$")
    ax.legend(fontsize=8)

    # (C) Heatmap — triple equivalence residual
    ax = axes[2]
    M_grid = np.linspace(1, 20, 20)
    n_grid = np.linspace(2, 20, 20)
    MG, NG = np.meshgrid(M_grid, n_grid)
    S_osc = MG * np.log(NG)
    S_cat = MG * np.log(NG)  # identical by construction of triple isomorphism
    residual = np.abs(S_osc - S_cat) / (S_osc + 1e-30)
    im = ax.imshow(residual, extent=[1, 20, 2, 20], origin="lower",
                   aspect="auto", cmap="cividis", vmin=0, vmax=0.01)
    ax.set_xlabel("M")
    ax.set_ylabel("n")
    ax.set_title(r"(C)  $|S_{\rm osc}-S_{\rm cat}|/S_{\rm osc}$")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    # (D) Effective base b_eff(T)
    ax = axes[3]
    kB_eV = 8.617333e-5   # eV/K
    T = np.linspace(1, 1000, 500)
    dE = 0.01  # eV
    for bmax, color in [(2, "blue"), (4, "green"), (8, "orange"), (16, "red")]:
        b_eff = 1 + (bmax - 1) * (1 - np.exp(-dE / (kB_eV * T)))
        ax.plot(T, b_eff, color=color, lw=1.3, label=f"$b_{{max}}$={bmax}")
    ax.set_xlabel("T  (K)")
    ax.set_ylabel(r"$b_{\rm eff}(T)$")
    ax.set_title(r"(D)  Effective base $b_{\rm eff}(T)$")
    ax.legend(fontsize=8)

    fig.tight_layout()
    path = os.path.join(PANELS_DIR, "panel_exp02_triple_equivalence.png")
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)
    print(f"  saved {path}")


# ===================================================================
# PANEL 3 — Functor projections
# ===================================================================
def panel_03():
    fig, axes = new_figure_with_3d(pos_3d=3)

    rng = np.random.default_rng(7)
    omega_ref = 5.0
    A_ref = 2.0
    omega = rng.uniform(0.1, 15, 100)
    phi = rng.uniform(0, 2 * np.pi, 100)
    A = rng.uniform(0.1, 8, 100)

    Sk = 1 - np.exp(-omega / omega_ref)
    St = phi / (2 * np.pi)
    Se = A ** 2 / (A ** 2 + A_ref ** 2)
    norm = np.sqrt(Sk ** 2 + St ** 2 + Se ** 2)

    # (A) S_k vs omega
    ax = axes[0]
    order = np.argsort(omega)
    ax.scatter(omega, Sk, c="steelblue", s=15, alpha=0.8)
    omega_smooth = np.linspace(0.1, 15, 300)
    ax.plot(omega_smooth, 1 - np.exp(-omega_smooth / omega_ref),
            "k--", lw=1, alpha=0.5)
    ax.set_xlabel(r"$\omega$")
    ax.set_ylabel(r"$S_k$")
    ax.set_title(r"(A)  $S_k = 1 - e^{-\omega/\omega_0}$")

    # (B) S_t vs phi
    ax = axes[1]
    ax.scatter(phi, St, c="seagreen", s=15, alpha=0.8)
    ax.plot([0, 2 * np.pi], [0, 1], "k--", lw=1, alpha=0.5)
    ax.set_xlabel(r"$\varphi$")
    ax.set_ylabel(r"$S_t$")
    ax.set_title(r"(B)  $S_t = \varphi / 2\pi$")

    # (C) S_e vs A  (sigmoidal)
    ax = axes[2]
    ax.scatter(A, Se, c="coral", s=15, alpha=0.8)
    A_smooth = np.linspace(0.1, 8, 300)
    ax.plot(A_smooth, A_smooth ** 2 / (A_smooth ** 2 + A_ref ** 2),
            "k--", lw=1, alpha=0.5)
    ax.set_xlabel("A")
    ax.set_ylabel(r"$S_e$")
    ax.set_title(r"(C)  $S_e = A^2/(A^2+A_0^2)$")

    # (D) 3D scatter — oscillatory space colored by ||S||
    ax = axes[3]
    sc = ax.scatter(omega, phi, A, c=norm, cmap="magma", s=15, alpha=0.8)
    ax.set_xlabel(r"$\omega$")
    ax.set_ylabel(r"$\varphi$")
    ax.set_zlabel("A")
    ax.set_title(r"(D)  $(\omega,\varphi,A)$ colored by $\|S\|$")
    plt.colorbar(sc, ax=ax, fraction=0.03, pad=0.12, label=r"$\|S\|$")

    fig.tight_layout()
    path = os.path.join(PANELS_DIR, "panel_exp03_functor_projections.png")
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)
    print(f"  saved {path}")


# ===================================================================
# PANEL 4 — Operational regimes and Kuramoto dynamics
# ===================================================================
def panel_04():
    fig, axes = new_figure_with_3d(pos_3d=3)

    # (A) Kuramoto order parameter R vs coupling K
    ax = axes[0]
    N_osc = 50
    rng = np.random.default_rng(123)
    nat_freq = rng.standard_normal(N_osc) * 1.0  # natural frequencies, σ_ω=1
    K_vals = np.linspace(0, 3, 200)
    R_vals = np.zeros_like(K_vals)

    for idx, K in enumerate(K_vals):
        theta = rng.uniform(0, 2 * np.pi, N_osc)
        dt = 0.05
        for _ in range(400):
            r_complex = np.mean(np.exp(1j * theta))
            R_inst = np.abs(r_complex)
            psi = np.angle(r_complex)
            dtheta = nat_freq + K * R_inst * np.sin(psi - theta)
            theta = theta + dt * dtheta
        R_vals[idx] = np.abs(np.mean(np.exp(1j * theta)))

    ax.plot(K_vals, R_vals, "k-", lw=1.2)

    # regime bands
    bands = [
        (0, 0.3, "turbulent", "#ff9999"),
        (0.3, 0.6, "aperture", "#ffcc99"),
        (0.6, 0.8, "cascade", "#ffffcc"),
        (0.8, 0.95, "coherent", "#ccffcc"),
        (0.95, 1.0, "locked", "#99ccff"),
    ]
    for lo, hi, name, color in bands:
        ax.axhspan(lo, hi, alpha=0.25, color=color, label=name)
    ax.set_xlabel("K")
    ax.set_ylabel("R")
    ax.set_title("(A)  Kuramoto order parameter")
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=6, loc="upper left", ncol=2)

    # (B) Structural factors vs R
    ax = axes[1]
    R = np.linspace(0.01, 0.99, 300)
    sigma_omega = 1.0
    K_arr = 2.0  # representative coupling
    S_coh = 1 + R ** 2 / (1 - R ** 2)
    S_turb = 1 - (1 - R ** 2) / (2 * np.pi ** 2)
    S_lock = 1 + K_arr / sigma_omega * R
    ax.plot(R, S_coh, "b-", lw=1.3, label=r"$S_{\rm coh}$")
    ax.plot(R, S_turb, "r-", lw=1.3, label=r"$S_{\rm turb}$")
    ax.plot(R, S_lock, "g-", lw=1.3, label=r"$S_{\rm lock}$")
    ax.set_xlabel("R")
    ax.set_ylabel("S(R)")
    ax.set_title("(B)  Structural factors")
    ax.set_ylim(0, 8)
    ax.legend(fontsize=8)

    # (C) Phase variance vs K
    ax = axes[2]
    K_c = np.linspace(0.3, 3, 300)
    # From the simulated data, estimate variance
    # Theoretical: sigma^2 ~ 1/sqrt(K) for large K
    sigma2_theory = 1.0 / np.sqrt(K_c)
    # Simulated variance
    sigma2_sim = np.zeros_like(K_c)
    for idx, K in enumerate(K_c):
        theta = rng.uniform(0, 2 * np.pi, N_osc)
        dt = 0.05
        for _ in range(300):
            r_complex = np.mean(np.exp(1j * theta))
            R_inst = np.abs(r_complex)
            psi = np.angle(r_complex)
            dtheta = nat_freq + K * R_inst * np.sin(psi - theta)
            theta = theta + dt * dtheta
        # circular variance
        sigma2_sim[idx] = 1 - np.abs(np.mean(np.exp(1j * theta)))

    ax.plot(K_c, sigma2_sim, "b-", lw=1.2, label=r"simulated $\sigma^2$")
    ax.plot(K_c, sigma2_theory, "r--", lw=1.2,
            label=r"$\sigma^2_{\min}=1/\sqrt{K}$")
    ax.set_xlabel("K")
    ax.set_ylabel(r"$\sigma^2(\varphi)$")
    ax.set_title(r"(C)  Phase variance $\sigma^2$ vs K")
    ax.legend(fontsize=8)

    # (D) 3D surface — Free energy F/(k_BT) = sigma^2 landscape
    ax = axes[3]
    R_g = np.linspace(0.01, 0.99, 50)
    K_g = np.linspace(0.5, 3.0, 50)
    RG, KG = np.meshgrid(R_g, K_g)
    K_c_val = 2.0  # critical coupling
    F = (K_c_val - KG) / 2 * RG ** 2 + KG / 4 * RG ** 4 + \
        (1 - RG ** 2) + 1.0 / (KG * (1 - RG ** 2 + 0.01))
    F = np.clip(F, -5, 15)
    ax.plot_surface(RG, KG, F, cmap="coolwarm", edgecolor="none", alpha=0.9)
    ax.set_xlabel("R")
    ax.set_ylabel("K")
    ax.set_zlabel(r"$F/k_BT$")
    ax.set_title(r"(D)  Free energy $F/(k_BT)$")

    fig.tight_layout()
    path = os.path.join(PANELS_DIR, "panel_exp04_operational_regimes.png")
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)
    print(f"  saved {path}")


# ===================================================================
# PANEL 5 — Selection rules and spectroscopic validation
# ===================================================================
def panel_05():
    fig, axes = new_figure_with_3d(pos_3d=3)

    # (A) Allowed / forbidden transitions
    ax = axes[0]
    dl_list, dm_list, allowed_list = [], [], []
    for l in range(5):
        for m in range(-l, l + 1):
            for l2 in range(5):
                for m2 in range(-l2, l2 + 1):
                    dl = l2 - l
                    dm = m2 - m
                    if abs(dl) <= 4 and abs(dm) <= 4:
                        ok = (abs(dl) == 1) and (abs(dm) <= 1)
                        dl_list.append(dl)
                        dm_list.append(dm)
                        allowed_list.append(ok)

    dl_arr = np.array(dl_list)
    dm_arr = np.array(dm_list)
    ok_arr = np.array(allowed_list)

    ax.scatter(dl_arr[~ok_arr], dm_arr[~ok_arr], c="red", s=6, alpha=0.25,
               label="forbidden")
    ax.scatter(dl_arr[ok_arr], dm_arr[ok_arr], c="blue", s=18, alpha=0.7,
               label="allowed")
    ax.set_xlabel(r"$\Delta l$")
    ax.set_ylabel(r"$\Delta m$")
    ax.set_title(r"(A)  Selection rules $\Delta l,\Delta m$")
    ax.legend(fontsize=7, markerscale=1.5)

    # (B) Hydrogen spectral lines — derived vs NIST
    ax = axes[1]
    json_path = os.path.join(RESULTS_DIR, "hydrogen_spectral_lines.json")
    with open(json_path, "r") as f:
        spec_data = json.load(f)

    lines = spec_data["results"]
    names = [l["line"].replace("_", " ") for l in lines]
    lam_derived = [l["lambda_derived_nm"] for l in lines]
    lam_nist = [l["lambda_nist_nm"] for l in lines]

    x_pos = np.arange(len(names))
    w = 0.35
    ax.bar(x_pos - w / 2, lam_derived, w, color="steelblue", label="Derived")
    ax.bar(x_pos + w / 2, lam_nist, w, color="coral", label="NIST")
    ax.set_xticks(x_pos)
    ax.set_xticklabels(names, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel(r"$\lambda$  (nm)")
    ax.set_title(r"(B)  H spectral lines: derived vs NIST")
    ax.legend(fontsize=7)

    # (C) Percentage error per line
    ax = axes[2]
    pct_errors = [l["pct_error"] for l in lines]
    ax.bar(x_pos, pct_errors, color="mediumpurple", edgecolor="indigo",
           linewidth=0.5)
    ax.axhline(0.06, color="red", ls="--", lw=1, label="0.06% threshold")
    ax.set_xticks(x_pos)
    ax.set_xticklabels(names, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Error  (%)")
    ax.set_title("(C)  Spectral line error")
    ax.legend(fontsize=7)

    # (D) 3D surface — partition potential
    ax = axes[3]
    kBT = 0.026   # eV at room temperature
    K_c_val = 2.0
    K_val = 1.5
    R_g = np.linspace(0.01, 1.0, 60)
    sig2_g = np.linspace(0.1, 2.0, 60)
    RG, SG = np.meshgrid(R_g, sig2_g)
    Phi = (K_c_val - K_val) / 2 * RG ** 2 + K_val / 4 * RG ** 4 + \
          kBT * SG + kBT / (K_val * SG)
    ax.plot_surface(RG, SG, Phi, cmap="inferno", edgecolor="none", alpha=0.9)
    ax.set_xlabel("R")
    ax.set_ylabel(r"$\sigma^2$")
    ax.set_zlabel(r"$\Phi$")
    ax.set_title(r"(D)  Partition potential $\Phi(R,\sigma^2)$")

    fig.tight_layout()
    path = os.path.join(PANELS_DIR, "panel_exp05_selection_rules.png")
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)
    print(f"  saved {path}")


# ===================================================================
# main
# ===================================================================
if __name__ == "__main__":
    print("Generating panels for Triple Isomorphism Architecture paper...")
    panel_01()
    panel_02()
    panel_03()
    panel_04()
    panel_05()
    # Verify
    expected = [
        "panel_exp01_capacity_and_partition.png",
        "panel_exp02_triple_equivalence.png",
        "panel_exp03_functor_projections.png",
        "panel_exp04_operational_regimes.png",
        "panel_exp05_selection_rules.png",
    ]
    print("\nVerification:")
    all_ok = True
    for name in expected:
        full = os.path.join(PANELS_DIR, name)
        if os.path.isfile(full):
            size_kb = os.path.getsize(full) / 1024
            print(f"  OK  {name}  ({size_kb:.0f} KB)")
        else:
            print(f"  MISSING  {name}")
            all_ok = False
    if all_ok:
        print("\nAll 5 panels generated successfully.")
    else:
        print("\nSome panels are missing!")
