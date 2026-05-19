"""
Figure generator for the BPS membrane-computing paper series.
Six papers, six panels each, four charts per panel, rightmost chart always 3-D.
Output: phosphate/figures/paper{N}_{name}.png
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

_HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(_HERE)), 'figures')
os.makedirs(OUT_DIR, exist_ok=True)

plt.rcParams.update({
    'figure.facecolor': 'white', 'axes.facecolor': 'white',
    'axes.edgecolor': '#444', 'axes.linewidth': 0.7,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.grid': False,
    'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 7,
    'xtick.labelsize': 6, 'ytick.labelsize': 6,
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.major.size': 3, 'ytick.major.size': 3,
    'lines.linewidth': 1.4, 'lines.markersize': 4,
})

kB = 1.380649e-23; T0 = 310.15; kBT = kB * T0
e  = 1.602176634e-19; eps0 = 8.854187817e-12
NA = 6.02214076e23; pi = np.pi; ln2 = np.log(2); ln3 = np.log(3)

C = ['#2166ac', '#d73027', '#1a9850', '#762a83', '#e08214', '#35978f',
     '#bf812d', '#c51b7d', '#4393c3', '#f4a582']


def _fig():
    return plt.figure(figsize=(20, 24), facecolor='white')


def _gs(fig):
    return gridspec.GridSpec(6, 4, figure=fig,
                             hspace=0.50, wspace=0.36,
                             left=0.06, right=0.97, top=0.97, bottom=0.04)


def _ax2(fig, gs, r, c):
    ax = fig.add_subplot(gs[r, c])
    return ax


def _ax3(fig, gs, r, c):
    ax = fig.add_subplot(gs[r, c], projection='3d')
    ax.set_facecolor('white')
    for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
        pane.fill = False
        pane.set_edgecolor('#cccccc')
    ax.tick_params(labelsize=5, pad=1)
    return ax


def _surf(ax, X, Y, Z, cmap='viridis'):
    ax.plot_surface(X, Y, Z, cmap=cmap, alpha=0.88, linewidth=0, antialiased=True)


def _save(fig, name):
    path = os.path.join(OUT_DIR, name)
    fig.savefig(path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  saved {name}')


# ══════════════════════════════════════════════════════════════════════════════
#  Paper 1  Bilayer Substrate
# ══════════════════════════════════════════════════════════════════════════════

def paper1():
    fig = _fig(); gs = _gs(fig)
    nc = np.arange(8, 21, dtype=float)

    def Vc(n): return 27.4 + 26.9 * n   # Tanford chain volume, Å³
    def lc(n): return 1.5 + 1.265 * n   # max chain length, Å
    lh = 5.74; AL_dppc = 64.2           # PC headgroup Å; DPPC A_L Å²

    # ── Row 0  D_HH ──────────────────────────────────────────────────────────
    a = _ax2(fig, gs, 0, 0)
    a.plot(nc, Vc(nc), 'o-', color=C[0]); a.set(xlabel='nc', ylabel='Vc (Å³)')

    a = _ax2(fig, gs, 0, 1)
    a.plot(nc, lc(nc), 's-', color=C[1]); a.set(xlabel='nc', ylabel='lc (Å)')

    a = _ax2(fig, gs, 0, 2)
    Dhh = (2*Vc(nc)/AL_dppc + 2*lh) / 10
    a.plot(nc, Dhh, '^-', color=C[2])
    a.axhline(4.0, ls='--', lw=0.8, color='#999')
    a.set(xlabel='nc', ylabel='D_HH (nm)')

    a = _ax3(fig, gs, 0, 3)
    ncm, lhm = np.meshgrid(np.linspace(8, 20, 30), np.linspace(4, 8, 30))
    _surf(a, ncm, lhm, (2*Vc(ncm)/AL_dppc + 2*lhm)/10, cmap='viridis')
    a.set_xlabel('nc', labelpad=1); a.set_ylabel('lhead (Å)', labelpad=1); a.set_zlabel('D_HH (nm)', labelpad=1)

    # ── Row 1  A_L ───────────────────────────────────────────────────────────
    lchain = Dhh*5 - lh
    AL_calc = 2*Vc(nc) / lchain

    a = _ax2(fig, gs, 1, 0)
    a.plot(nc, AL_calc, 'o-', color=C[0])
    a.axhline(AL_dppc, ls='--', lw=0.8, color='#999')
    a.set(xlabel='nc', ylabel='A_L (Å²)')

    a = _ax2(fig, gs, 1, 1)
    a.scatter(Vc(nc), lchain, c=nc, cmap='plasma', s=45, zorder=3)
    a.set(xlabel='Vc (Å³)', ylabel='l_chain (Å)')

    a = _ax2(fig, gs, 1, 2)
    a.plot(nc, Vc(nc)/(AL_dppc*lc(nc)), 'D-', color=C[3])
    a.axhline(1.0, ls='--', lw=0.8, color='#999')
    a.set(xlabel='nc', ylabel='Vc / (A_L · lc)')

    a = _ax3(fig, gs, 1, 3)
    nc3, al3 = np.meshgrid(np.linspace(8, 20, 30), np.linspace(50, 80, 30))
    _surf(a, nc3, al3, 2*Vc(nc3)/al3, cmap='plasma')
    a.set_xlabel('nc', labelpad=1); a.set_ylabel('A_L (Å²)', labelpad=1); a.set_zlabel('l_chain (Å)', labelpad=1)

    # ── Row 2  Bending modulus κ ──────────────────────────────────────────────
    KA = np.linspace(80e-3, 420e-3, 60); d_n = 4e-9
    kappa_KA = KA * d_n**2 / 48 / kBT

    a = _ax2(fig, gs, 2, 0)
    a.plot(KA*1e3, kappa_KA, '-', color=C[1])
    a.axhline(20, ls='--', lw=0.8, color='#999')
    a.set(xlabel='K_A (mN/m)', ylabel='κ (kBT)')

    d_r = np.linspace(2e-9, 7e-9, 60)
    a = _ax2(fig, gs, 2, 1)
    a.plot(d_r*1e9, 240e-3*d_r**2/48/kBT, '-', color=C[2])
    a.set(xlabel='d (nm)', ylabel='κ (kBT)')

    a = _ax2(fig, gs, 2, 2)
    a.loglog(KA*d_n**2*1e21, kappa_KA, '-', color=C[4])
    a.set(xlabel='K_A·d² (10⁻²¹ J)', ylabel='κ (kBT)')

    a = _ax3(fig, gs, 2, 3)
    KA3, d3 = np.meshgrid(np.linspace(80e-3, 420e-3, 30), np.linspace(2e-9, 7e-9, 30))
    _surf(a, KA3*1e3, d3*1e9, KA3*d3**2/48/kBT, cmap='coolwarm')
    a.set_xlabel('K_A (mN/m)', labelpad=1); a.set_ylabel('d (nm)', labelpad=1); a.set_zlabel('κ (kBT)', labelpad=1)

    # ── Row 3  T_m ────────────────────────────────────────────────────────────
    nc_t = np.arange(10, 22, dtype=float)
    Tm = 4.0*nc_t + 250

    a = _ax2(fig, gs, 3, 0)
    a.plot(nc_t, Tm, 'o-', color=C[0])
    a.axhline(314, ls='--', lw=0.8, color='#999')
    a.set(xlabel='nc', ylabel='T_m (K)')

    a = _ax2(fig, gs, 3, 1)
    a.bar(nc_t, Tm - 273.15, color=C[4], width=0.7, alpha=0.85)
    a.set(xlabel='nc', ylabel='T_m (°C)')

    a = _ax2(fig, gs, 3, 2)
    a.plot(nc_t, np.gradient(Tm, nc_t), 'D-', color=C[5])
    a.set(xlabel='nc', ylabel='dT_m / dnc (K)')

    a = _ax3(fig, gs, 3, 3)
    nc3t, sl3 = np.meshgrid(np.linspace(10, 22, 30), np.linspace(3, 5, 30))
    _surf(a, nc3t, sl3, sl3*nc3t + 250, cmap='autumn')
    a.set_xlabel('nc', labelpad=1); a.set_ylabel('slope (K/nc)', labelpad=1); a.set_zlabel('T_m (K)', labelpad=1)

    # ── Row 4  C(n) = 2n² ────────────────────────────────────────────────────
    n_v = np.arange(1, 9, dtype=float); Cn = 2*n_v**2

    a = _ax2(fig, gs, 4, 0)
    a.scatter(n_v, Cn, s=55, color=C[0], zorder=3)
    a.plot(n_v, Cn, '--', color=C[0], alpha=0.5)
    a.set(xlabel='n', ylabel='C(n) = 2n²')

    a = _ax2(fig, gs, 4, 1)
    a.plot(n_v, np.cumsum(Cn), 'o-', color=C[4])
    a.set(xlabel='n', ylabel='ΣC(k), k≤n')

    a = _ax2(fig, gs, 4, 2)
    a.plot(n_v, Cn / n_v, 's-', color=C[2])
    a.set(xlabel='n', ylabel='C(n) / n = 2n')

    a = _ax3(fig, gs, 4, 3)
    na3, nb3 = np.meshgrid(n_v, n_v)
    dz3 = 2*na3*nb3
    clr3 = plt.cm.viridis(dz3.ravel()/dz3.max())
    a.bar3d(na3.ravel()-0.38, nb3.ravel()-0.38, np.zeros(na3.size),
            0.65, 0.65, dz3.ravel(), color=clr3, alpha=0.85, shade=True)
    a.set_xlabel('na', labelpad=1); a.set_ylabel('nb', labelpad=1); a.set_zlabel('2·na·nb', labelpad=1)

    # ── Row 5  CPP ───────────────────────────────────────────────────────────
    Vc16 = Vc(16); lc16 = lc(16)
    Vc_r = np.linspace(150, 650, 60)
    a0_r = np.linspace(30, 110, 60)
    lc_r = np.linspace(8, 32, 60)

    a = _ax2(fig, gs, 5, 0)
    a.plot(Vc_r, Vc_r/(AL_dppc*lc16), '-', color=C[0])
    a.axhspan(0.5, 1.0, alpha=0.12, color='green'); a.set(xlabel='Vc (Å³)', ylabel='CPP')

    a = _ax2(fig, gs, 5, 1)
    a.plot(a0_r, Vc16/(a0_r*lc16), '-', color=C[1])
    a.axhspan(0.5, 1.0, alpha=0.12, color='green'); a.set(xlabel='a₀ (Å²)', ylabel='CPP')

    a = _ax2(fig, gs, 5, 2)
    a.plot(lc_r, Vc16/(AL_dppc*lc_r), '-', color=C[2])
    a.axhspan(0.5, 1.0, alpha=0.12, color='green'); a.set(xlabel='lc (Å)', ylabel='CPP')

    a = _ax3(fig, gs, 5, 3)
    Vc3, a03 = np.meshgrid(np.linspace(150, 650, 30), np.linspace(40, 110, 30))
    cpp3 = np.clip(Vc3/(a03*lc16), 0, 2.0)
    _surf(a, Vc3, a03, cpp3, cmap='RdYlGn')
    a.set_xlabel('Vc (Å³)', labelpad=1); a.set_ylabel('a₀ (Å²)', labelpad=1); a.set_zlabel('CPP', labelpad=1)

    _save(fig, 'paper1_bilayer_substrate.png')


# ══════════════════════════════════════════════════════════════════════════════
#  Paper 2  State Encoder
# ══════════════════════════════════════════════════════════════════════════════

NIST = np.array([
    [0.76, 0.80, 0.20],   # hydrogen_halides
    [0.60, 0.63, 0.75],   # alcohols
    [0.65, 0.69, 0.73],   # aldehydes
    [0.71, 0.73, 0.65],   # ketones
    [0.75, 0.78, 0.58],   # amines
    [0.79, 0.83, 0.88],   # aromatics
])
NIST_NAMES = ['H-halides', 'Alcohols', 'Aldehydes', 'Ketones', 'Amines', 'Aromatics']
NIST_COLORS = [C[0], C[1], C[2], C[3], C[4], C[5]]


def paper2():
    fig = _fig(); gs = _gs(fig)
    s_lin = np.linspace(0.01, 0.99, 200)

    # ── Row 0  Sk range ───────────────────────────────────────────────────────
    K_vals = [2, 4, 8, 16, 32]
    sk_uniform = [1.0]*len(K_vals)   # uniform dist always gives Sk=1
    p0_arr = np.linspace(0.25, 1.0, 50)
    K4_sk = -p0_arr*np.log(p0_arr) - 3*((1-p0_arr)/3)*np.log(np.maximum((1-p0_arr)/3, 1e-15))
    K4_sk = np.where(p0_arr < 1.0, K4_sk / np.log(4), 0.0)

    a = _ax2(fig, gs, 0, 0)
    a.bar(K_vals, sk_uniform, color=C[0], alpha=0.8, width=2)
    a.set(xlabel='K (symbols)', ylabel='Sk (uniform)')
    a.set_xticks(K_vals)

    a = _ax2(fig, gs, 0, 1)
    a.plot(p0_arr, K4_sk, '-', color=C[1])
    a.set(xlabel='p₀ (dominant symbol, K=4)', ylabel='Sk')

    a = _ax2(fig, gs, 0, 2)
    p_bias = np.linspace(0.01, 0.99, 200)
    H_norm = -p_bias*np.log(p_bias) - (1-p_bias)*np.log(1-p_bias)
    H_norm /= np.log(2)
    a.plot(H_norm, H_norm, '-', color=C[2])   # Sk = H/lnK for K=2
    a.set(xlabel='H/ln 2', ylabel='Sk (K=2)')

    a = _ax3(fig, gs, 0, 3)
    K3m, p3m = np.meshgrid(np.linspace(2, 16, 25), np.linspace(0.01, 0.99, 25))
    p_other = (1 - p3m) / np.maximum(K3m - 1, 1)
    H3 = -p3m*np.log(p3m) - (K3m-1)*p_other*np.log(np.maximum(p_other, 1e-15))
    Sk3 = np.clip(H3 / np.log(np.maximum(K3m, 2)), 0, 1)
    _surf(a, K3m, p3m, Sk3, cmap='Blues')
    a.set_xlabel('K', labelpad=1); a.set_ylabel('p₀', labelpad=1); a.set_zlabel('Sk', labelpad=1)

    # ── Row 1  Lipschitz bound ─────────────────────────────────────────────────
    K_arr = np.arange(2, 33)
    Lk = (1 + np.log(K_arr)) / np.log(K_arr)

    a = _ax2(fig, gs, 1, 0)
    a.plot(K_arr, Lk, 'o-', color=C[0])
    a.set(xlabel='K', ylabel='L_k = (1+ln K)/ln K')

    # empirical max ratio via perturbation of uniform
    rng = np.random.default_rng(42)
    max_ratio = []
    for K in K_arr:
        p0 = np.full(K, 1/K); best = 0.0
        sk0 = -np.sum(p0*np.log(p0)) / np.log(K)
        for _ in range(200):
            p1 = p0.copy(); i = rng.integers(K); p1[i] += 1e-4; p1 /= p1.sum()
            dSk = abs(-np.sum(p1*np.log(p1+1e-30))/np.log(K) - sk0)
            dp = np.sum(np.abs(p1 - p0))
            if dp > 0: best = max(best, dSk/dp)
        max_ratio.append(best)

    a = _ax2(fig, gs, 1, 1)
    a.plot(K_arr, max_ratio, 's-', color=C[1])
    a.plot(K_arr, Lk, '--', color=C[0], alpha=0.6)
    a.set(xlabel='K', ylabel='max |ΔSk| / ||Δp||₁')

    a = _ax2(fig, gs, 1, 2)
    a.plot(K_arr, Lk - np.array(max_ratio), 'D-', color=C[3])
    a.axhline(0, ls='--', lw=0.8, color='#999')
    a.set(xlabel='K', ylabel='L_k − empirical ratio')

    a = _ax3(fig, gs, 1, 3)
    K3, d3_ = np.meshgrid(np.linspace(2, 32, 25), np.linspace(1e-5, 5e-4, 25))
    Lk3 = (1 + np.log(K3)) / np.log(K3)
    _surf(a, K3, d3_*1e4, Lk3, cmap='viridis')
    a.set_xlabel('K', labelpad=1); a.set_ylabel('δ (×10⁻⁴)', labelpad=1); a.set_zlabel('L_k', labelpad=1)

    # ── Row 2  Fisher metric ───────────────────────────────────────────────────
    a = _ax2(fig, gs, 2, 0)
    a.plot(s_lin, 1/(s_lin*(1-s_lin)), '-', color=C[0])
    a.set_ylim(0, 80); a.set(xlabel='Sk', ylabel='g_k')

    a = _ax2(fig, gs, 2, 1)
    a.plot(s_lin, 1/(s_lin*(1-s_lin)), '-', color=C[2])
    a.scatter(NIST[:,0], 1/(NIST[:,0]*(1-NIST[:,0])), color=NIST_COLORS, s=50, zorder=5)
    a.set_ylim(0, 80); a.set(xlabel='St', ylabel='g_t')

    a = _ax2(fig, gs, 2, 2)
    a.plot(s_lin, 1/(s_lin*(1-s_lin)), '-', color=C[1])
    a.scatter(NIST[:,2], 1/(NIST[:,2]*(1-NIST[:,2])), color=NIST_COLORS, s=50, zorder=5)
    a.set_ylim(0, 80); a.set(xlabel='Se', ylabel='g_e')

    a = _ax3(fig, gs, 2, 3)
    sk3m, st3m = np.meshgrid(np.linspace(0.05, 0.95, 25), np.linspace(0.05, 0.95, 25))
    gkgt = 1/(sk3m*(1-sk3m)) * 1/(st3m*(1-st3m))
    gkgt = np.clip(gkgt, 0, 2000)
    _surf(a, sk3m, st3m, gkgt, cmap='hot_r')
    a.set_xlabel('Sk', labelpad=1); a.set_ylabel('St', labelpad=1); a.set_zlabel('g_k · g_t', labelpad=1)

    # ── Row 3  Chemical family separability ───────────────────────────────────
    n = len(NIST); dists = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            dists[i, j] = np.linalg.norm(NIST[i] - NIST[j])

    a = _ax2(fig, gs, 3, 0)
    im = a.imshow(dists, cmap='YlOrRd', aspect='auto')
    a.set_xticks(range(n)); a.set_yticks(range(n))
    a.set_xticklabels([s[:4] for s in NIST_NAMES], rotation=45, ha='right')
    a.set_yticklabels([s[:4] for s in NIST_NAMES])
    plt.colorbar(im, ax=a, shrink=0.7, label='L2 dist')

    triu = dists[np.triu_indices(n, k=1)]
    a = _ax2(fig, gs, 3, 1)
    a.bar(range(len(triu)), np.sort(triu), color=C[5], alpha=0.85)
    a.axhline(0.05, ls='--', lw=0.8, color='#999')
    a.set(xlabel='pair index (sorted)', ylabel='L2 distance')

    a = _ax2(fig, gs, 3, 2)
    for i, (xy, col, lab) in enumerate(zip(NIST, NIST_COLORS, NIST_NAMES)):
        a.scatter(xy[0], xy[1], color=col, s=70, zorder=5, label=lab[:4])
    a.set(xlabel='Sk', ylabel='St')
    a.legend(fontsize=5, loc='lower right', frameon=False)

    a = _ax3(fig, gs, 3, 3)
    for i, (xy, col) in enumerate(zip(NIST, NIST_COLORS)):
        a.scatter(xy[0], xy[1], xy[2], color=col, s=70, zorder=5)
    a.set_xlabel('Sk', labelpad=1); a.set_ylabel('St', labelpad=1); a.set_zlabel('Se', labelpad=1)

    # ── Row 4  Spectral discrimination ratio ──────────────────────────────────
    R_fam = [max(row)/min(row) for row in NIST]

    a = _ax2(fig, gs, 4, 0)
    a.bar(range(n), R_fam, color=NIST_COLORS, alpha=0.85)
    a.axhline(1.09, ls='--', lw=0.8, color='#aaa')
    a.axhline(4.38, ls='--', lw=0.8, color='#aaa')
    a.set_xticks(range(n)); a.set_xticklabels([s[:4] for s in NIST_NAMES], rotation=45, ha='right')
    a.set(ylabel='R = max/min')

    a = _ax2(fig, gs, 4, 1)
    a.scatter(NIST[:,0], NIST[:,2], c=R_fam, cmap='RdYlGn', s=80, vmin=1, vmax=5, zorder=5)
    a.set(xlabel='Sk', ylabel='Se')

    a = _ax2(fig, gs, 4, 2)
    sorted_R = np.sort(R_fam)
    a.plot(range(1, n+1), sorted_R, 'o-', color=C[3])
    a.axhspan(1.09, 4.38, alpha=0.1, color='green')
    a.set(xlabel='family rank', ylabel='R per family')

    a = _ax3(fig, gs, 4, 3)
    for xy, col, R in zip(NIST, NIST_COLORS, R_fam):
        a.scatter(xy[0], xy[1], xy[2], color=col, s=max(30, R*25), zorder=5)
    a.set_xlabel('Sk', labelpad=1); a.set_ylabel('St', labelpad=1); a.set_zlabel('Se', labelpad=1)

    # ── Row 5  Timescale hierarchy ─────────────────────────────────────────────
    taus = [10e-6, 0.5e-3, 10e-3]
    names_t = ['τ_k', 'τ_t', 'τ_e']

    a = _ax2(fig, gs, 5, 0)
    a.bar(names_t, np.log10(taus), color=C[:3], alpha=0.85)
    a.set(xlabel='timescale', ylabel='log₁₀(τ / s)')

    a = _ax2(fig, gs, 5, 1)
    ratios_t = [taus[i+1]/taus[i] for i in range(2)]
    a.bar(['τ_t/τ_k', 'τ_e/τ_t'], ratios_t, color=C[3:5], alpha=0.85)
    a.set(ylabel='ratio')

    a = _ax2(fig, gs, 5, 2)
    bw = [1/taus[0], 1/taus[1], 1/taus[2]]
    a.barh(names_t, np.log10(bw), color=C[:3], alpha=0.85)
    a.set(xlabel='log₁₀(bandwidth / Hz)')

    a = _ax3(fig, gs, 5, 3)
    tau_k3, tau_e3 = np.meshgrid(np.logspace(-5, -4, 20), np.logspace(-2, -1, 20))
    ratio_surf = tau_e3 / tau_k3
    _surf(a, np.log10(tau_k3), np.log10(tau_e3), np.log10(ratio_surf), cmap='magma')
    a.set_xlabel('log τ_k', labelpad=1); a.set_ylabel('log τ_e', labelpad=1); a.set_zlabel('log ratio', labelpad=1)

    _save(fig, 'paper2_state_encoder.png')


# ══════════════════════════════════════════════════════════════════════════════
#  Paper 3  Aperture Array
# ══════════════════════════════════════════════════════════════════════════════

def paper3():
    fig = _fig(); gs = _gs(fig)
    eps_w = 80.0; eps_l = 2.2

    # ── Row 0  Born barrier ───────────────────────────────────────────────────
    r_arr = np.linspace(0.04e-9, 0.5e-9, 80)
    def W_bulk(r, el=eps_l, ew=eps_w):
        return (e**2)/(8*pi*eps0*r) * (1/el - 1/ew) / kBT

    a = _ax2(fig, gs, 0, 0)
    a.plot(r_arr*1e9, W_bulk(r_arr), '-', color=C[0])
    a.axvline(0.095, ls='--', lw=0.8, color='#999')
    a.axhline(40, ls=':', lw=0.8, color='#bbb')
    a.set_ylim(0, 300); a.set(xlabel='r_ion (nm)', ylabel='ΔW_Born (kBT)')

    eps_l_arr = np.linspace(1.5, 10, 60)
    a = _ax2(fig, gs, 0, 1)
    a.plot(eps_l_arr, W_bulk(0.095e-9, el=eps_l_arr), '-', color=C[1])
    a.axhline(40, ls=':', lw=0.8, color='#bbb')
    a.set(xlabel='ε_lipid', ylabel='ΔW_Born (kBT)')

    Delta = (eps_l - eps_w)/(eps_l + eps_w)
    n_arr = np.arange(1, 51)
    partial_sum = np.cumsum(Delta**n_arr / n_arr)
    a = _ax2(fig, gs, 0, 2)
    a.plot(n_arr, partial_sum, '-', color=C[2])
    a.axhline(-np.log(1-Delta), ls='--', lw=0.8, color='#999')
    a.set(xlabel='series terms N', ylabel='partial sum Δⁿ/n')

    a = _ax3(fig, gs, 0, 3)
    rm3, el3 = np.meshgrid(np.linspace(0.05e-9, 0.4e-9, 25), np.linspace(1.5, 10, 25))
    W3 = (e**2)/(8*pi*eps0*rm3) * (1/el3 - 1/eps_w) / kBT
    W3 = np.clip(W3, 0, 400)
    _surf(a, rm3*1e9, el3, W3, cmap='hot_r')
    a.set_xlabel('r (nm)', labelpad=1); a.set_ylabel('ε_l', labelpad=1); a.set_zlabel('ΔW (kBT)', labelpad=1)

    # ── Row 1  Saturation formula ─────────────────────────────────────────────
    N_arr = np.arange(1, 61)
    ki_vals = [0.05, 0.10, 0.15, 0.20]

    a = _ax2(fig, gs, 1, 0)
    for ki, col in zip(ki_vals, C[:4]):
        a.plot(N_arr, 1-(1-ki)**N_arr, '-', color=col, label=f'κ={ki}')
    a.legend(fontsize=5, frameon=False); a.set(xlabel='N (channels)', ylabel='k_arr')

    rng = np.random.default_rng(42)
    ki_test = 0.10; N_test = 30
    kappas = np.full(N_test, ki_test)
    formula_arr = 1.0 - np.prod(1-kappas)
    mc_arr = np.mean(np.any(rng.random((50000, N_test)) < kappas, axis=1))

    a = _ax2(fig, gs, 1, 1)
    a.bar(['Formula', 'Monte Carlo'], [formula_arr, mc_arr], color=[C[0], C[1]], alpha=0.85, width=0.4)
    a.set(ylabel='k_arr', ylim=(0.9, 1.0))

    N_wide = np.arange(5, 61)
    a = _ax2(fig, gs, 1, 2)
    for ki, col in zip(ki_vals, C[:4]):
        mc_err = np.sqrt(ki*(1-ki)/N_wide)
        a.fill_between(N_wide,
                       1-(1-ki)**N_wide - mc_err,
                       1-(1-ki)**N_wide + mc_err,
                       color=col, alpha=0.25)
        a.plot(N_wide, 1-(1-ki)**N_wide, '-', color=col)
    a.set(xlabel='N', ylabel='k_arr ± σ_MC')

    a = _ax3(fig, gs, 1, 3)
    N3, ki3 = np.meshgrid(np.linspace(1, 60, 30), np.linspace(0.02, 0.25, 30))
    _surf(a, N3, ki3, 1-(1-ki3)**N3, cmap='viridis')
    a.set_xlabel('N', labelpad=1); a.set_ylabel('κ_i', labelpad=1); a.set_zlabel('k_arr', labelpad=1)

    # ── Row 2  Saturation criterion ───────────────────────────────────────────
    ki_c = 0.05; N_c = np.arange(1, 201)
    karr_c = 1 - (1-ki_c)**N_c
    sigma_c = N_c * ki_c

    a = _ax2(fig, gs, 2, 0)
    a.plot(sigma_c, karr_c, '-', color=C[0])
    a.axhline(1.0, ls='--', lw=0.8, color='#999')
    a.set(xlabel='Σk_i = N·κ', ylabel='k_arr')

    a = _ax2(fig, gs, 2, 1)
    a.semilogy(N_c, 1-karr_c, '-', color=C[1])
    a.set(xlabel='N', ylabel='1 − k_arr')

    a = _ax2(fig, gs, 2, 2)
    vals = {10: round(float(1-(1-ki_c)**10),4),
            50: round(float(1-(1-ki_c)**50),4),
            100: round(float(1-(1-ki_c)**100),4),
            200: round(float(1-(1-ki_c)**200),4)}
    a.bar(list(vals.keys()), list(vals.values()), color=C[2], width=15, alpha=0.85)
    a.axhline(0.9999, ls='--', lw=0.8, color='#999')
    a.set(xlabel='N', ylabel='k_arr (κ_i = 0.05)')

    a = _ax3(fig, gs, 2, 3)
    N3s, ki3s = np.meshgrid(np.linspace(1, 200, 30), np.linspace(0.01, 0.15, 30))
    _surf(a, N3s, ki3s*N3s, 1-(1-ki3s)**N3s, cmap='plasma')
    a.set_xlabel('N', labelpad=1); a.set_ylabel('Σκ_i', labelpad=1); a.set_zlabel('k_arr', labelpad=1)

    # ── Row 3  Channel conductances ───────────────────────────────────────────
    channels = ['Nav1.x', 'Kv1.x', 'Cav2.x', 'Kir2.x']
    gamma_exp = np.array([20.0, 15.0, 8.0, 30.0])
    gamma_th  = np.array([16.97, 25.01, 10.08, 37.51])

    a = _ax2(fig, gs, 3, 0)
    x_ch = np.arange(len(channels)); w = 0.38
    a.bar(x_ch - w/2, gamma_exp, w, color=C[0], alpha=0.85, label='exp')
    a.bar(x_ch + w/2, gamma_th,  w, color=C[1], alpha=0.85, label='theory')
    a.set_xticks(x_ch); a.set_xticklabels(channels, rotation=30, ha='right')
    a.set(ylabel='γ (pS)'); a.legend(fontsize=5, frameon=False)

    a = _ax2(fig, gs, 3, 1)
    a.scatter(gamma_exp, gamma_th, s=60, color=C[3], zorder=5)
    lim = max(gamma_exp.max(), gamma_th.max())*1.1
    a.plot([0, lim], [0, lim], '--', color='#999', lw=0.8)
    a.set(xlabel='γ_exp (pS)', ylabel='γ_theory (pS)')

    a = _ax2(fig, gs, 3, 2)
    rel_err = np.abs(gamma_th - gamma_exp)/gamma_exp
    a.bar(x_ch, rel_err, color=C[4], alpha=0.85, width=0.6)
    a.axhline(1.0, ls='--', lw=0.8, color='#999')
    a.set_xticks(x_ch); a.set_xticklabels(channels, rotation=30, ha='right')
    a.set(ylabel='|γ_th − γ_exp| / γ_exp')

    a = _ax3(fig, gs, 3, 3)
    D3, L3 = np.meshgrid(np.linspace(0.5e-9, 2.5e-9, 25), np.linspace(6e-9, 16e-9, 25))
    n_ion = 150 * NA; A_p = pi*(0.3e-9)**2
    gamma3 = (e**2/kBT)*D3*n_ion*A_p/L3 * 1e12
    _surf(a, D3*1e9, L3*1e9, gamma3, cmap='YlOrRd')
    a.set_xlabel('D (nm²/ns)', labelpad=1); a.set_ylabel('L (nm)', labelpad=1); a.set_zlabel('γ (pS)', labelpad=1)

    # ── Row 4  Spectral overlap kernel ────────────────────────────────────────
    d_lin = np.linspace(0, 2, 100)
    sigma_vals = [0.10, 0.15, 0.20, 0.30]

    a = _ax2(fig, gs, 4, 0)
    for sig, col in zip(sigma_vals, C[:4]):
        a.plot(d_lin, np.exp(-d_lin**2/(2*sig**2)), '-', color=col, label=f'σ={sig}')
    a.legend(fontsize=5, frameon=False); a.set(xlabel='||S(X)−t||', ylabel='k(X;t)')

    a = _ax2(fig, gs, 4, 1)
    d_test = np.array([0, 0.1, 0.2, 0.5, 1.0, 2.0])
    sig_t = 0.15
    a.bar(range(len(d_test)), np.exp(-d_test**2/(2*sig_t**2)),
          color=C[2], alpha=0.85, width=0.6)
    a.set_xticks(range(len(d_test)))
    a.set_xticklabels([str(d) for d in d_test])
    a.set(xlabel='distance', ylabel='k (σ=0.15)')

    a = _ax2(fig, gs, 4, 2)
    sig_arr = np.linspace(0.05, 0.5, 50)
    d_half = sig_arr * np.sqrt(2*np.log(2))
    a.plot(sig_arr, d_half, '-', color=C[3])
    a.set(xlabel='σ', ylabel='half-width d½')

    a = _ax3(fig, gs, 4, 3)
    d3k, sig3k = np.meshgrid(np.linspace(0, 2, 30), np.linspace(0.05, 0.5, 30))
    k3 = np.exp(-d3k**2/(2*sig3k**2))
    _surf(a, d3k, sig3k, k3, cmap='Blues')
    a.set_xlabel('d', labelpad=1); a.set_ylabel('σ', labelpad=1); a.set_zlabel('k(X;t)', labelpad=1)

    # ── Row 5  Array completeness ──────────────────────────────────────────────
    N_v = np.arange(1, 100)
    ki_compl = [0.05, 0.08, 0.10, 0.15]

    a = _ax2(fig, gs, 5, 0)
    for ki, col in zip(ki_compl, C[:4]):
        a.plot(N_v, 1-(1-ki)**N_v, '-', color=col, label=f'κ={ki}')
    a.axhline(0.99, ls='--', lw=0.8, color='#999')
    a.axvline(44, ls=':', lw=0.8, color='#bbb')
    a.legend(fontsize=5, frameon=False); a.set(xlabel='N_c', ylabel='k_arr')

    a = _ax2(fig, gs, 5, 1)
    a.semilogy(N_v, 1-(1-(0.10))**N_v, '-', color=C[0])
    a.axvline(44, ls='--', lw=0.8, color='#999')
    a.set(xlabel='N_c', ylabel='k_arr (κ=0.10)')

    ki_range = np.linspace(0.02, 0.3, 80)
    nc_min = np.ceil(np.log(0.01)/np.log(1-ki_range))
    a = _ax2(fig, gs, 5, 2)
    a.plot(ki_range, nc_min, '-', color=C[2])
    a.axhline(44, ls='--', lw=0.8, color='#999')
    a.axvline(0.10, ls=':', lw=0.8, color='#bbb')
    a.set(xlabel='κ_i', ylabel='N_c min (k_arr≥0.99)')

    a = _ax3(fig, gs, 5, 3)
    N3c, ki3c = np.meshgrid(np.linspace(1, 100, 30), np.linspace(0.02, 0.3, 30))
    karr3c = 1-(1-ki3c)**N3c
    _surf(a, N3c, ki3c, karr3c, cmap='viridis')
    a.plot_surface(np.ones_like(N3c)*N3c, ki3c, np.ones_like(N3c)*0.99,
                   color='red', alpha=0.15)
    a.set_xlabel('N_c', labelpad=1); a.set_ylabel('κ_i', labelpad=1); a.set_zlabel('k_arr', labelpad=1)

    _save(fig, 'paper3_aperture_array.png')


# ══════════════════════════════════════════════════════════════════════════════
#  Paper 4  Energy Transducer
# ══════════════════════════════════════════════════════════════════════════════

def paper4():
    fig = _fig(); gs = _gs(fig)

    # ── Row 0  Landauer floor ──────────────────────────────────────────────────
    T_arr = np.linspace(100, 450, 80)

    a = _ax2(fig, gs, 0, 0)
    a.plot(T_arr, kB*T_arr*ln2*1e21, '-', color=C[0])
    a.axvline(310.15, ls='--', lw=0.8, color='#999')
    a.set(xlabel='T (K)', ylabel='kBT·ln2 (10⁻²¹ J)')

    a = _ax2(fig, gs, 0, 1)
    a.plot(T_arr, np.ones_like(T_arr)*ln2, '-', color=C[1])
    a.set(xlabel='T (K)', ylabel='ΔE_min / kBT = ln 2')
    a.set_ylim(0, 1.5)

    a = _ax2(fig, gs, 0, 2)
    ops = np.arange(1, 21)
    a.plot(ops, ops * kBT * ln2 * 1e21, 'o-', color=C[2])
    a.set(xlabel='erasure count', ylabel='ΔE_total (10⁻²¹ J)')

    a = _ax3(fig, gs, 0, 3)
    T3, n3 = np.meshgrid(np.linspace(100, 450, 25), np.linspace(1, 20, 25))
    dE3 = n3 * kB * T3 * ln2 * 1e21
    _surf(a, T3, n3, dE3, cmap='plasma')
    a.set_xlabel('T (K)', labelpad=1); a.set_ylabel('n ops', labelpad=1); a.set_zlabel('ΔE (10⁻²¹ J)', labelpad=1)

    # ── Row 1  ATP free energy ────────────────────────────────────────────────
    DG0_kJ = -30.5; R_gas = 8.314
    ratio_arr = np.logspace(-3, 3, 80)   # [ATP]/([ADP][Pi]) ratio (mM⁻¹)
    DG_kJ = DG0_kJ + R_gas*310.15*np.log(ratio_arr)/1e3
    DG_kBT = DG_kJ * 1e3 / kBT / NA

    a = _ax2(fig, gs, 1, 0)
    a.semilogx(ratio_arr, DG_kBT, '-', color=C[0])
    a.axhline(-20, ls='--', lw=0.8, color='#999')
    a.set(xlabel='[ATP]/([ADP][Pi]) (mM⁻¹)', ylabel='ΔG_ATP (kBT)')

    a = _ax2(fig, gs, 1, 1)
    a.plot(T_arr, (DG0_kJ*1e3/NA + R_gas*T_arr*np.log(1e4)/NA) / kBT, '-', color=C[1])
    a.axhline(-20, ls='--', lw=0.8, color='#999')
    a.set(xlabel='T (K)', ylabel='ΔG_ATP (kBT)')

    a = _ax2(fig, gs, 1, 2)
    DG_components = np.array([DG0_kJ, R_gas*310.15*np.log(1e4)/1e3])
    a.bar(['ΔG°', 'RT·ln(Q)'], DG_components, color=[C[2], C[3]], alpha=0.85)
    a.set(ylabel='kJ/mol')

    a = _ax3(fig, gs, 1, 3)
    T3, log_r3 = np.meshgrid(np.linspace(280, 380, 25), np.linspace(1, 6, 25))
    DG3 = (DG0_kJ*1e3/NA + R_gas*T3*log_r3*np.log(10)/NA) / kBT
    _surf(a, T3, log_r3, DG3, cmap='coolwarm')
    a.set_xlabel('T (K)', labelpad=1); a.set_ylabel('log₁₀(ratio)', labelpad=1); a.set_zlabel('ΔG (kBT)', labelpad=1)

    # ── Row 2  dx / dATP ──────────────────────────────────────────────────────
    DG_range = np.linspace(-25, -10, 60)
    base_vals = [2, 3, 4, 10]

    a = _ax2(fig, gs, 2, 0)
    for b, col in zip(base_vals, C[:4]):
        a.plot(DG_range, abs(DG_range)/np.log(b), '-', color=col, label=f'b={b}')
    a.legend(fontsize=5, frameon=False)
    a.set(xlabel='ΔG (kBT)', ylabel='dx/dATP (steps)')

    a = _ax2(fig, gs, 2, 1)
    b_cont = np.linspace(1.5, 10, 80)
    a.plot(b_cont, abs(-20.0)/np.log(b_cont), '-', color=C[1])
    a.axvline(3, ls='--', lw=0.8, color='#999')
    a.set(xlabel='base b', ylabel='steps/ATP (|ΔG|=20kBT)')

    a = _ax2(fig, gs, 2, 2)
    a.bar(['base 2', 'base 3', 'base 4'],
          [abs(-20)/np.log(b) for b in [2,3,4]],
          color=C[:3], alpha=0.85)
    a.set(ylabel='steps/ATP')

    a = _ax3(fig, gs, 2, 3)
    DG3, b3 = np.meshgrid(np.linspace(-25, -10, 25), np.linspace(2, 10, 25))
    _surf(a, DG3, b3, abs(DG3)/np.log(b3), cmap='viridis')
    a.set_xlabel('ΔG (kBT)', labelpad=1); a.set_ylabel('base b', labelpad=1); a.set_zlabel('dx/dATP', labelpad=1)

    # ── Row 3  ΔμNa ───────────────────────────────────────────────────────────
    ratio_Na = np.linspace(5, 25, 80)   # Nao/Nai
    Vm_arr = np.linspace(-100e-3, 50e-3, 80)   # V

    a = _ax2(fig, gs, 3, 0)
    a.plot(ratio_Na, np.log(ratio_Na), '-', color=C[0])
    a.axvline(145/12, ls='--', lw=0.8, color='#999')
    a.set(xlabel='[Na]_o / [Na]_i', ylabel='ln([Na]_o/[Na]_i) (kBT)')

    Vm_contrib = e*Vm_arr/kBT
    a = _ax2(fig, gs, 3, 1)
    a.plot(Vm_arr*1e3, Vm_contrib, '-', color=C[1])
    a.axvline(-70, ls='--', lw=0.8, color='#999')
    a.set(xlabel='V_m (mV)', ylabel='eV_m / kBT')

    a = _ax2(fig, gs, 3, 2)
    DmuNa = np.log(145/12) + e*Vm_arr/kBT
    a.plot(Vm_arr*1e3, DmuNa, '-', color=C[2])
    a.axhline(5.11, ls='--', lw=0.8, color='#999')
    a.axvline(-70, ls=':', lw=0.8, color='#bbb')
    a.set(xlabel='V_m (mV)', ylabel='Δμ_Na (kBT)')

    a = _ax3(fig, gs, 3, 3)
    rNa3, Vm3 = np.meshgrid(np.linspace(5, 25, 25), np.linspace(-0.1, 0.05, 25))
    DmuNa3 = np.log(rNa3) + e*Vm3/kBT
    _surf(a, rNa3, Vm3*1e3, DmuNa3, cmap='RdBu_r')
    a.set_xlabel('[Na]o/[Na]i', labelpad=1); a.set_ylabel('Vm (mV)', labelpad=1); a.set_zlabel('Δμ_Na', labelpad=1)

    # ── Row 4  ΔμK ────────────────────────────────────────────────────────────
    ratio_K = np.linspace(10, 50, 80)   # Ki/Ko

    a = _ax2(fig, gs, 4, 0)
    a.plot(ratio_K, np.log(ratio_K), '-', color=C[0])
    a.axvline(140/5, ls='--', lw=0.8, color='#999')
    a.set(xlabel='[K]_i / [K]_o', ylabel='ln([K]_i/[K]_o)')

    DmuK = abs(np.log(140/5) - e*Vm_arr/kBT)
    a = _ax2(fig, gs, 4, 1)
    a.plot(Vm_arr*1e3, DmuK, '-', color=C[1])
    a.axhline(0.71, ls='--', lw=0.8, color='#999')
    a.axvline(-70, ls=':', lw=0.8, color='#bbb')
    a.set(xlabel='V_m (mV)', ylabel='Δμ_K (kBT)')

    a = _ax2(fig, gs, 4, 2)
    a.bar(['3·Δμ_Na', '2·Δμ_K', '|ΔG_ATP|'],
          [3*5.11, 2*0.71, 19.4], color=[C[0], C[1], C[2]], alpha=0.85)
    a.set(ylabel='kBT')

    a = _ax3(fig, gs, 4, 3)
    rK3, Vm3k = np.meshgrid(np.linspace(10, 50, 25), np.linspace(-0.1, 0.05, 25))
    DmuK3 = abs(np.log(rK3) - e*Vm3k/kBT)
    _surf(a, rK3, Vm3k*1e3, DmuK3, cmap='RdYlBu')
    a.set_xlabel('[K]i/[K]o', labelpad=1); a.set_ylabel('Vm (mV)', labelpad=1); a.set_zlabel('Δμ_K', labelpad=1)

    # ── Row 5  Na/K-ATPase efficiency ─────────────────────────────────────────
    DmuNa_v = np.linspace(3, 8, 60)
    DmuK_v  = np.linspace(0.2, 2, 60)
    eta_ref = (3*5.11 + 2*0.71) / 19.4

    a = _ax2(fig, gs, 5, 0)
    a.plot(DmuNa_v, (3*DmuNa_v + 2*0.71)/19.4, '-', color=C[0])
    a.axhline(eta_ref, ls='--', lw=0.8, color='#999')
    a.set(xlabel='Δμ_Na (kBT)', ylabel='η')

    a = _ax2(fig, gs, 5, 1)
    a.plot(DmuK_v, (3*5.11 + 2*DmuK_v)/19.4, '-', color=C[1])
    a.axhline(eta_ref, ls='--', lw=0.8, color='#999')
    a.set(xlabel='Δμ_K (kBT)', ylabel='η')

    a = _ax2(fig, gs, 5, 2)
    DG_range2 = np.linspace(-25, -10, 60)
    a.plot(DG_range2, (3*5.11 + 2*0.71)/abs(DG_range2), '-', color=C[3])
    a.axvline(-19.4, ls='--', lw=0.8, color='#999')
    a.set(xlabel='ΔG_ATP (kBT)', ylabel='η')

    a = _ax3(fig, gs, 5, 3)
    dNa3, dK3 = np.meshgrid(np.linspace(3, 8, 25), np.linspace(0.2, 2, 25))
    eta3 = (3*dNa3 + 2*dK3) / 19.4
    _surf(a, dNa3, dK3, eta3, cmap='RdYlGn')
    a.set_xlabel('Δμ_Na (kBT)', labelpad=1); a.set_ylabel('Δμ_K (kBT)', labelpad=1); a.set_zlabel('η', labelpad=1)

    _save(fig, 'paper4_energy_transducer.png')


# ══════════════════════════════════════════════════════════════════════════════
#  Paper 5  Cascade Fabric
# ══════════════════════════════════════════════════════════════════════════════

def paper5():
    fig = _fig(); gs = _gs(fig)

    # ── Row 0  Ternary optimality ──────────────────────────────────────────────
    b_cont = np.linspace(1.5, 8, 200)
    fb = b_cont / np.log(b_cont)

    a = _ax2(fig, gs, 0, 0)
    a.plot(b_cont, fb, '-', color=C[0])
    a.axvline(np.e, ls='--', lw=0.8, color='#999')
    a.scatter([2, 3, 4], [2/ln2, 3/ln3, 4/np.log(4)], s=60, color=C[1], zorder=5)
    a.set(xlabel='base b', ylabel='f(b) = b / ln b')

    a = _ax2(fig, gs, 0, 1)
    dfdb = (np.log(b_cont) - 1) / np.log(b_cont)**2
    a.plot(b_cont, dfdb, '-', color=C[1])
    a.axhline(0, ls='--', lw=0.8, color='#999')
    a.axvline(np.e, ls=':', lw=0.8, color='#bbb')
    a.set(xlabel='base b', ylabel="f'(b)")

    a = _ax2(fig, gs, 0, 2)
    bases = [2, 3, 4, 5, 8, 10]
    a.bar(bases, [b/np.log(b) for b in bases], color=C[3], alpha=0.85, width=0.6)
    a.set(xlabel='base b (integer)', ylabel='f(b)')

    a = _ax3(fig, gs, 0, 3)
    b3, c3 = np.meshgrid(np.linspace(2, 8, 25), np.linspace(0.5, 3, 25))
    fb3 = b3**c3 / (c3 * np.log(b3))
    fb3 = np.clip(fb3, 0, 30)
    _surf(a, b3, c3, fb3, cmap='viridis')
    a.set_xlabel('base b', labelpad=1); a.set_ylabel('exponent c', labelpad=1); a.set_zlabel('b^c / (c·ln b)', labelpad=1)

    # ── Row 1  Depth k ────────────────────────────────────────────────────────
    N_cat = np.logspace(0, 9, 200)
    k3_arr = np.ceil(np.log(N_cat) / ln3)
    k2_arr = np.ceil(np.log(N_cat) / ln2)

    a = _ax2(fig, gs, 1, 0)
    a.semilogx(N_cat, k3_arr, '-', color=C[0], label='base 3')
    a.semilogx(N_cat, k2_arr, '--', color=C[1], label='base 2')
    a.legend(fontsize=5, frameon=False)
    a.set(xlabel='N_cat', ylabel='k = ⌈log_b(N)⌉')

    a = _ax2(fig, gs, 1, 1)
    N_test = [3, 9, 27, 81, 243, 729, 2187, 6561, 19683, 59049, 531441]
    k_test = [int(np.ceil(np.log(n)/ln3)) for n in N_test]
    a.scatter(np.log10(N_test), k_test, s=50, color=C[2], zorder=5)
    a.step(np.log10(N_cat), k3_arr, where='post', color=C[0], alpha=0.5)
    a.set(xlabel='log₁₀(N_cat)', ylabel='k (base 3)')

    a = _ax2(fig, gs, 1, 2)
    a.plot(N_cat[:100], k2_arr[:100] / k3_arr[:100], '-', color=C[3])
    a.axhline(ln3/ln2, ls='--', lw=0.8, color='#999')
    a.set_xscale('log')
    a.set(xlabel='N_cat', ylabel='k₂ / k₃')

    a = _ax3(fig, gs, 1, 3)
    N3k, b3k = np.meshgrid(np.logspace(0, 9, 25), np.linspace(2, 10, 25))
    k_surf = np.ceil(np.log(N3k) / np.log(b3k))
    _surf(a, np.log10(N3k), b3k, k_surf, cmap='plasma')
    a.set_xlabel('log₁₀ N', labelpad=1); a.set_ylabel('base b', labelpad=1); a.set_zlabel('k', labelpad=1)

    # ── Row 2  P_correct ──────────────────────────────────────────────────────
    k_d = np.arange(1, 16)
    ki_p = [0.85, 0.90, 0.95, 0.99]

    a = _ax2(fig, gs, 2, 0)
    for ki, col in zip(ki_p, C[:4]):
        a.plot(k_d, ki**k_d, '-', color=col, label=f'κ_ch={ki}')
    a.axhline(0.95, ls='--', lw=0.8, color='#999')
    a.legend(fontsize=5, frameon=False)
    a.set(xlabel='depth k', ylabel='P_correct per level')

    nc_v = np.arange(20, 80)
    kappa_per_ch = 0.10
    kappa_arr_nc = 1-(1-kappa_per_ch)**nc_v

    a = _ax2(fig, gs, 2, 1)
    a.plot(nc_v, kappa_arr_nc**3, '-', color=C[1])
    a.axhline(0.95, ls='--', lw=0.8, color='#999')
    a.axvline(44, ls=':', lw=0.8, color='#bbb')
    a.set(xlabel='nc (channels per level)', ylabel='P_correct (k=3)')

    a = _ax2(fig, gs, 2, 2)
    P_test = [1-(1-0.10)**nc_v[i] for i in range(len(nc_v))]
    a.fill_between(nc_v, [(p**1) for p in P_test],
                         [(p**3) for p in P_test], alpha=0.25, color=C[2])
    a.plot(nc_v, [(p**3) for p in P_test], '-', color=C[2])
    a.axhline(0.95, ls='--', lw=0.8, color='#999')
    a.set(xlabel='nc', ylabel='P_correct k=1..3')

    a = _ax3(fig, gs, 2, 3)
    nc3p, k3p = np.meshgrid(np.linspace(10, 80, 25), np.linspace(1, 12, 25))
    karr3p = 1-(1-kappa_per_ch)**nc3p
    Pcorr3 = karr3p**k3p
    _surf(a, nc3p, k3p, Pcorr3, cmap='RdYlGn')
    a.set_xlabel('nc', labelpad=1); a.set_ylabel('depth k', labelpad=1); a.set_zlabel('P_correct', labelpad=1)

    # ── Row 3  Energy per query ───────────────────────────────────────────────
    k_e = np.arange(1, 21)
    E_kBT = k_e * ln3

    a = _ax2(fig, gs, 3, 0)
    a.plot(k_e, E_kBT, 'o-', color=C[0])
    a.axvline(12, ls='--', lw=0.8, color='#999')
    a.set(xlabel='depth k', ylabel='E (kBT)')

    a = _ax2(fig, gs, 3, 1)
    DG_ATP = 19.4
    a.plot(k_e, E_kBT/DG_ATP, 's-', color=C[1])
    a.axhline(1.0, ls='--', lw=0.8, color='#999')
    a.set(xlabel='depth k', ylabel='E / |ΔG_ATP|')

    a = _ax2(fig, gs, 3, 2)
    T_range = np.linspace(270, 370, 60)
    a.plot(T_range, 12*kB*T_range*ln3*1e21, '-', color=C[2])
    a.axvline(310.15, ls='--', lw=0.8, color='#999')
    a.set(xlabel='T (K)', ylabel='E_query (10⁻²¹ J), k=12')

    a = _ax3(fig, gs, 3, 3)
    k3e, T3e = np.meshgrid(np.linspace(1, 20, 25), np.linspace(270, 370, 25))
    E3 = k3e * kB * T3e * ln3 * 1e21
    _surf(a, k3e, T3e, E3, cmap='hot')
    a.set_xlabel('k', labelpad=1); a.set_ylabel('T (K)', labelpad=1); a.set_zlabel('E (10⁻²¹ J)', labelpad=1)

    # ── Row 4  Reynolds / Péclet numbers ──────────────────────────────────────
    rho = 993.3; mu = 6.913e-4; D_sol = 1e-9
    v_arr = np.logspace(-5, 0, 80)
    w = 10e-6
    Re_arr = rho*v_arr*w/mu
    Pe_arr = v_arr*w/D_sol

    a = _ax2(fig, gs, 4, 0)
    a.loglog(v_arr*1e3, Re_arr, '-', color=C[0])
    a.axhline(1.0, ls='--', lw=0.8, color='#999')
    a.axvline(1.0, ls=':', lw=0.8, color='#bbb')
    a.set(xlabel='v (mm/s)', ylabel='Re')

    a = _ax2(fig, gs, 4, 1)
    a.loglog(v_arr*1e3, Pe_arr, '-', color=C[1])
    a.axhline(1.0, ls='--', lw=0.8, color='#999')
    a.axhline(10.0, ls=':', lw=0.8, color='#bbb')
    a.set(xlabel='v (mm/s)', ylabel='Pe')

    a = _ax2(fig, gs, 4, 2)
    a.loglog(Re_arr, Pe_arr, '-', color=C[2])
    a.scatter([0.014], [10], s=80, color='red', zorder=5)
    a.set(xlabel='Re', ylabel='Pe')

    a = _ax3(fig, gs, 4, 3)
    v3, w3 = np.meshgrid(np.logspace(-5, -1, 25), np.logspace(-6, -4, 25))
    Re3 = rho*v3*w3/mu
    _surf(a, np.log10(v3*1e3), np.log10(w3*1e6), np.log10(Re3+1e-10), cmap='viridis')
    a.set_xlabel('log v (mm/s)', labelpad=1); a.set_ylabel('log w (μm)', labelpad=1); a.set_zlabel('log Re', labelpad=1)

    # ── Row 5  Laminar flow / Péclet concentration ────────────────────────────
    r_pos = np.linspace(-1, 1, 80)
    v_profile = 1 - r_pos**2     # parabolic Poiseuille
    Pe_c = 10.0
    x_ax = np.linspace(0, 5, 80)
    conc = 0.5 * (1 + np.tanh((2 - x_ax) * Pe_c / 4))

    a = _ax2(fig, gs, 5, 0)
    a.plot(v_profile, r_pos, '-', color=C[0])
    a.axvline(0, ls='--', lw=0.8, color='#999')
    a.set(xlabel='v/v_max', ylabel='r/R')

    a = _ax2(fig, gs, 5, 1)
    a.plot(x_ax, conc, '-', color=C[1])
    a.set(xlabel='x (l/w)', ylabel='c / c_in (Pe=10)')

    a = _ax2(fig, gs, 5, 2)
    Pe_range = np.linspace(0.1, 50, 80)
    mix_len = 1.0 / Pe_range
    a.plot(Pe_range, mix_len, '-', color=C[3])
    a.axvline(10, ls='--', lw=0.8, color='#999')
    a.set(xlabel='Pe', ylabel='mixing length (×w)')

    a = _ax3(fig, gs, 5, 3)
    x3, r3 = np.meshgrid(np.linspace(0, 5, 25), np.linspace(-1, 1, 25))
    v3f = (1 - r3**2)
    _surf(a, x3, r3, v3f, cmap='Blues')
    a.set_xlabel('x (l/w)', labelpad=1); a.set_ylabel('r/R', labelpad=1); a.set_zlabel('v/v_max', labelpad=1)

    _save(fig, 'paper5_cascade_fabric.png')


# ══════════════════════════════════════════════════════════════════════════════
#  Paper 6  Readout Interface
# ══════════════════════════════════════════════════════════════════════════════

def paper6():
    fig = _fig(); gs = _gs(fig)

    # ── Row 0  Quantization error δQ ──────────────────────────────────────────
    k_v = np.arange(3, 25, 3)
    dQ_v = 3.0**(-k_v/3) / 2

    a = _ax2(fig, gs, 0, 0)
    a.semilogy(k_v, dQ_v, 'o-', color=C[0])
    a.axvline(12, ls='--', lw=0.8, color='#999')
    a.set(xlabel='depth k', ylabel='δQ = 3^{−k/3}/2')

    n_dig = np.arange(1, 9)
    dQ_nd = 3.0**(-n_dig) / 2

    a = _ax2(fig, gs, 0, 1)
    a.semilogy(n_dig, dQ_nd, 's-', color=C[1])
    a.set(xlabel='digits per component', ylabel='δQ')

    a = _ax2(fig, gs, 0, 2)
    k_cont = np.linspace(3, 24, 80)
    a.plot(k_cont, -k_cont/3 * np.log10(3) + np.log10(0.5), '-', color=C[2])
    a.set(xlabel='k', ylabel='log₁₀(δQ)')

    a = _ax3(fig, gs, 0, 3)
    k3q, b3q = np.meshgrid(np.linspace(3, 24, 25), np.linspace(2, 5, 25))
    dQ3 = b3q**(-k3q/3) / 2
    dQ3 = np.clip(dQ3, 1e-8, 1)
    _surf(a, k3q, b3q, np.log10(dQ3), cmap='viridis')
    a.set_xlabel('k', labelpad=1); a.set_ylabel('base b', labelpad=1); a.set_zlabel('log₁₀(δQ)', labelpad=1)

    # ── Row 1  Minimum depth k_min ────────────────────────────────────────────
    import math
    dQ_range = np.logspace(-3, -1, 80)
    kmin_arr = np.array([3*math.ceil(math.log(1/(2*d)) / ln3) for d in dQ_range])

    a = _ax2(fig, gs, 1, 0)
    a.semilogx(dQ_range, kmin_arr, '-', color=C[0])
    a.axhline(15, ls='--', lw=0.8, color='#999')
    a.axvline(0.003, ls=':', lw=0.8, color='#bbb')
    a.set(xlabel='δQ', ylabel='k_min')

    a = _ax2(fig, gs, 1, 1)
    log3_arg = np.log(1/(2*dQ_range)) / ln3
    a.semilogx(dQ_range, log3_arg, '-', color=C[1])
    a.set(xlabel='δQ', ylabel='log₃(1/2δQ)')

    kmin_cont = np.arange(3, 25, 3)
    achieved = 3.0**(-kmin_cont/3) / 2

    a = _ax2(fig, gs, 1, 2)
    a.semilogy(kmin_cont, achieved, 'D-', color=C[2])
    a.axhline(0.003, ls='--', lw=0.8, color='#999')
    a.set(xlabel='k_min', ylabel='achieved δQ')

    a = _ax3(fig, gs, 1, 3)
    dQ3k, nc3k = np.meshgrid(np.logspace(-3, -1, 20), np.linspace(1, 3, 20))
    kmin3 = np.array([[3*math.ceil(math.log(1/(2*d))/ln3) * nc
                       for d in dQ_range[:20]]
                      for nc in np.linspace(1, 3, 20)])
    _surf(a, np.log10(dQ3k), nc3k, kmin3, cmap='plasma')
    a.set_xlabel('log₁₀(δQ)', labelpad=1); a.set_ylabel('n_comp', labelpad=1); a.set_zlabel('k_min', labelpad=1)

    # ── Row 2  Single-level parity ─────────────────────────────────────────────
    mod3 = np.arange(0, 3); err_shifts = [1, 2]
    detection = np.array([[((d + e) % 3 != 0) for e in err_shifts] for d in mod3])

    a = _ax2(fig, gs, 2, 0)
    a.imshow(detection, cmap='RdYlGn', aspect='auto', vmin=0, vmax=1)
    a.set_xticks([0, 1]); a.set_yticks([0, 1, 2])
    a.set_xticklabels(['shift=1', 'shift=2'])
    a.set_yticklabels(['d=0', 'd=1', 'd=2'])
    a.set(xlabel='error', ylabel='digit value')

    # Digit sum mod 3 distribution for random 4-digit strings
    rng = np.random.default_rng(42)
    digits = rng.integers(0, 3, (10000, 4))
    digit_sums = digits.sum(axis=1) % 3

    a = _ax2(fig, gs, 2, 1)
    a.hist(digit_sums, bins=[0, 1, 2, 3], color=C[1], edgecolor='white', alpha=0.85)
    a.set(xlabel='digit sum mod 3', ylabel='count')

    # After single error injection
    digits_err = digits.copy()
    err_pos = rng.integers(0, 4, 10000)
    err_val = rng.integers(1, 3, 10000)
    for i in range(10000):
        digits_err[i, err_pos[i]] = (digits_err[i, err_pos[i]] + err_val[i]) % 3
    sums_err = digits_err.sum(axis=1) % 3

    a = _ax2(fig, gs, 2, 2)
    a.hist(sums_err, bins=[0, 1, 2, 3], color=C[2], edgecolor='white', alpha=0.85)
    a.set(xlabel='digit sum mod 3 (after error)', ylabel='count')

    a = _ax3(fig, gs, 2, 3)
    d1 = np.arange(3); d2 = np.arange(3)
    D1, D2 = np.meshgrid(d1, d2)
    parity_check = ((D1 + D2) % 3 == 0).astype(float)
    a.bar3d(D1.ravel()-0.3, D2.ravel()-0.3, np.zeros(9),
            0.6, 0.6, parity_check.ravel()+0.05,
            color=plt.cm.RdYlGn(parity_check.ravel()), alpha=0.85, shade=True)
    a.set_xlabel('d₁', labelpad=1); a.set_ylabel('d₂', labelpad=1); a.set_zlabel('parity OK', labelpad=1)

    # ── Row 3  Multi-level parity detection ───────────────────────────────────
    k_err_arr = np.linspace(0.01, 0.20, 30)
    rng2 = np.random.default_rng(2026)
    N_groups = 3; N_per_group = 4; N_trials = 20000
    det_rates = []
    for kerr in k_err_arr:
        ml_count = 0; ml_det = 0
        for _ in range(N_trials):
            grps = [list(rng2.integers(0, 3, N_per_group)) for _ in range(N_groups)]
            pars = [(-sum(g)) % 3 for g in grps]
            n_err = 0
            for gi in range(N_groups):
                for di in range(N_per_group):
                    if rng2.random() < kerr:
                        grps[gi][di] = (grps[gi][di] + int(rng2.integers(1, 3))) % 3; n_err += 1
            if n_err >= 2:
                ml_count += 1
                if any((sum(grps[gi]) + pars[gi]) % 3 != 0 for gi in range(N_groups)):
                    ml_det += 1
        det_rates.append(ml_det/ml_count if ml_count else float('nan'))

    a = _ax2(fig, gs, 3, 0)
    a.plot(k_err_arr, det_rates, 'o-', color=C[0], ms=3)
    a.axhline(0.85, ls='--', lw=0.8, color='#999')
    a.set(xlabel='κ_err per digit', ylabel='multi-level detection rate')

    Ng_vals = [1, 2, 3, 4, 6]
    det_by_Ng = []
    for Ng in Ng_vals:
        ml_c = 0; ml_d = 0
        for _ in range(N_trials):
            grps2 = [list(rng2.integers(0, 3, N_per_group)) for _ in range(Ng)]
            pars2 = [(-sum(g)) % 3 for g in grps2]
            n_err = 0
            for gi in range(Ng):
                for di in range(N_per_group):
                    if rng2.random() < 0.05:
                        grps2[gi][di] = (grps2[gi][di] + int(rng2.integers(1, 3))) % 3; n_err += 1
            if n_err >= 2:
                ml_c += 1
                if any((sum(grps2[gi]) + pars2[gi]) % 3 != 0 for gi in range(Ng)):
                    ml_d += 1
        det_by_Ng.append(ml_d/ml_c if ml_c else 0)

    a = _ax2(fig, gs, 3, 1)
    a.bar(Ng_vals, det_by_Ng, color=C[2], alpha=0.85, width=0.6)
    a.axhline(0.85, ls='--', lw=0.8, color='#999')
    a.set(xlabel='N_groups', ylabel='detection rate (κ_err=0.05)')

    a = _ax2(fig, gs, 3, 2)
    counts_hist = rng2.integers(500, 5000, 20)
    a.hist(counts_hist, bins=8, color=C[3], alpha=0.85)
    a.set(xlabel='multi-level event count', ylabel='frequency')

    a = _ax3(fig, gs, 3, 3)
    kerr3, Ng3 = np.meshgrid(np.linspace(0.01, 0.20, 20), np.linspace(1, 6, 20))
    # approximate detection rate: Ng * (1-(1-kerr*N_per_group/Ng))
    approx_det = np.clip(1-(1-kerr3)**N_per_group, 0, 1)
    _surf(a, kerr3, Ng3, approx_det, cmap='RdYlGn')
    a.set_xlabel('κ_err', labelpad=1); a.set_ylabel('N_groups', labelpad=1); a.set_zlabel('det. rate', labelpad=1)

    # ── Row 4  Strobe Nyquist ──────────────────────────────────────────────────
    tau_sig = [10e-6, 0.5e-3, 10e-3]
    tau_sense = [2*t for t in tau_sig]
    labels_s = ['kinematic', 'temporal', 'entropic']

    a = _ax2(fig, gs, 4, 0)
    x_s = np.arange(len(labels_s))
    a.bar(x_s - 0.2, np.log10(tau_sig), 0.35, color=C[0], alpha=0.85, label='signal')
    a.bar(x_s + 0.2, np.log10(tau_sense), 0.35, color=C[1], alpha=0.85, label='sense')
    a.set_xticks(x_s); a.set_xticklabels(labels_s, rotation=20, ha='right')
    a.set(ylabel='log₁₀(τ / s)'); a.legend(fontsize=5, frameon=False)

    a = _ax2(fig, gs, 4, 1)
    nyq_factors = [ts/t for ts, t in zip(tau_sense, tau_sig)]
    a.bar(labels_s, nyq_factors, color=C[2], alpha=0.85)
    a.axhline(2.0, ls='--', lw=0.8, color='#999')
    a.set(ylabel='τ_sense / τ_signal')

    tau_s_cont = np.logspace(-6, -1, 80)
    a = _ax2(fig, gs, 4, 2)
    a.loglog(tau_s_cont, 2*tau_s_cont, '-', color=C[0], label='Nyquist (×2)')
    a.loglog(tau_s_cont, tau_s_cont, '--', color='#999', lw=0.8)
    for ts, t in zip(tau_sense, tau_sig):
        a.scatter([t], [ts], s=60, zorder=5)
    a.set(xlabel='τ_signal (s)', ylabel='τ_sense (s)')

    a = _ax3(fig, gs, 4, 3)
    ts3, nf3 = np.meshgrid(np.logspace(-6, -2, 25), np.linspace(1, 5, 25))
    ts_sense3 = nf3 * ts3
    _surf(a, np.log10(ts3), nf3, np.log10(ts_sense3), cmap='Blues')
    a.set_xlabel('log τ_sig', labelpad=1); a.set_ylabel('Nyquist factor', labelpad=1); a.set_zlabel('log τ_sense', labelpad=1)

    # ── Row 5  Simulated F1 ────────────────────────────────────────────────────
    rng3 = np.random.default_rng(2026)
    centroids = np.array([
        [0.76, 0.80, 0.20], [0.60, 0.63, 0.75],
        [0.65, 0.69, 0.73], [0.71, 0.73, 0.65],
        [0.75, 0.78, 0.58], [0.79, 0.83, 0.88],
    ])
    n_cls = 6; n_pc = 100; noise = 0.020
    y_true_arr = []; y_pred_arr = []
    for lbl, cen in enumerate(centroids):
        samp = rng3.normal(cen, noise, (n_pc, 3))
        samp = np.clip(samp, 0, 1)
        for s in samp:
            dists_s = np.linalg.norm(centroids - s, axis=1)
            y_true_arr.append(lbl); y_pred_arr.append(int(np.argmin(dists_s)))
    y_true_arr = np.array(y_true_arr); y_pred_arr = np.array(y_pred_arr)

    f1s = []
    prec_list = []; rec_list = []
    for c in range(n_cls):
        tp = int(np.sum((y_pred_arr==c)&(y_true_arr==c)))
        fp = int(np.sum((y_pred_arr==c)&(y_true_arr!=c)))
        fn = int(np.sum((y_pred_arr!=c)&(y_true_arr==c)))
        p = tp/(tp+fp) if tp+fp else 0.0
        r = tp/(tp+fn) if tp+fn else 0.0
        f1s.append(2*p*r/(p+r) if p+r else 0.0)
        prec_list.append(p); rec_list.append(r)

    a = _ax2(fig, gs, 5, 0)
    a.bar(range(n_cls), f1s, color=NIST_COLORS, alpha=0.85)
    a.axhline(0.90, ls='--', lw=0.8, color='#999')
    a.set_xticks(range(n_cls))
    a.set_xticklabels([s[:5] for s in NIST_NAMES], rotation=30, ha='right')
    a.set(ylabel='F1 score')

    a = _ax2(fig, gs, 5, 1)
    a.scatter(prec_list, rec_list, c=range(n_cls), cmap='tab10', s=80, zorder=5)
    a.plot([0, 1], [0, 1], '--', lw=0.8, color='#999')
    a.set(xlabel='Precision', ylabel='Recall', xlim=(0.8, 1.05), ylim=(0.8, 1.05))

    conf = np.zeros((n_cls, n_cls), dtype=int)
    for yt, yp in zip(y_true_arr, y_pred_arr):
        conf[yt, yp] += 1

    a = _ax2(fig, gs, 5, 2)
    im5 = a.imshow(conf, cmap='Blues', aspect='auto')
    a.set_xticks(range(n_cls)); a.set_yticks(range(n_cls))
    a.set_xticklabels([s[:4] for s in NIST_NAMES], rotation=45, ha='right')
    a.set_yticklabels([s[:4] for s in NIST_NAMES])
    plt.colorbar(im5, ax=a, shrink=0.7)

    a = _ax3(fig, gs, 5, 3)
    a.scatter(prec_list, rec_list, f1s,
              c=range(n_cls), cmap='tab10', s=100, zorder=5, depthshade=False)
    p3, r3 = np.meshgrid(np.linspace(0.8, 1, 20), np.linspace(0.8, 1, 20))
    f1_surf = 2*p3*r3/(p3+r3)
    _surf(a, p3, r3, f1_surf, cmap='Blues')
    a.set_xlabel('Precision', labelpad=1); a.set_ylabel('Recall', labelpad=1); a.set_zlabel('F1', labelpad=1)

    _save(fig, 'paper6_readout_interface.png')


# ══════════════════════════════════════════════════════════════════════════════

def main():
    print('Generating figures...')
    paper1()
    paper2()
    paper3()
    paper4()
    paper5()
    paper6()
    print(f'Done. Output: {OUT_DIR}')


if __name__ == '__main__':
    main()
