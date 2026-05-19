"""
Physical constants and system parameters used across all validation modules.
All values in SI units unless noted; derived convenience quantities labelled clearly.
"""

import math

# ── Fundamental constants ─────────────────────────────────────────────────────
kB          = 1.380649e-23      # J K⁻¹      Boltzmann
h_planck    = 6.62607015e-34    # J s         Planck
NA          = 6.02214076e23     # mol⁻¹       Avogadro
e_charge    = 1.60217663e-19    # C           elementary charge
eps0        = 8.854187817e-12   # F m⁻¹       permittivity of free space
R_gas       = 8.314462618       # J mol⁻¹ K⁻¹ molar gas constant
pi          = math.pi

# ── Physiological reference conditions ───────────────────────────────────────
T_physio    = 310.15            # K   (37 °C)
T_DPPC_gel  = 314.0             # K   DPPC gel↔fluid transition
T_DPPC_ref  = 323.15            # K   DPPC fluid-phase reference (50 °C, Nagle 2000)
kBT         = kB * T_physio     # J   (~4.283e-21 J)
kBT_kJ_mol  = R_gas * T_physio / 1000  # kJ mol⁻¹ (~2.578 kJ/mol)

# ── Derived information-theory constants ─────────────────────────────────────
ln2  = math.log(2)              # 0.6931
ln3  = math.log(3)              # 1.0986
ln_e = 1.0                      # ln(e) by definition

# ── DPPC lipid parameters (Tanford 1980; Nagle & Tristram-Nagle 2000) ────────
DPPC_NC         = 16            # palmitoyl chain carbons
DPPC_T_melt     = 314.0         # K   gel-to-fluid transition
DPPC_A_L        = 0.642         # nm² area per lipid (fluid phase, 323 K)
DPPC_d_c        = 4.0           # nm  hydrophobic core thickness
DPPC_kappa_kBT  = 20.0          # kBT bending modulus (Evans & Rawicz 1990)
DPPC_K_A        = 240e-3        # N m⁻¹  area-compressibility modulus (Evans 1990)
DPPC_S_param    = 0.55          # (-)  orientational order parameter, fluid phase
DPPC_l_head_A   = 5.74         # Ang  PC headgroup axial extension (Nagle 2000)

# ── Ion radii (Hille 2001) ────────────────────────────────────────────────────
r_Na_crystal    = 0.095e-9      # m   Na⁺ crystal radius
r_K_crystal     = 0.133e-9      # m   K⁺  crystal radius
r_Cl_crystal    = 0.181e-9      # m   Cl⁻ crystal radius
r_Na_hydrated   = 0.184e-9      # m   Na⁺ first-shell hydrated radius
r_K_hydrated    = 0.132e-9      # m   K⁺  first-shell hydrated radius

# ── Dielectric constants ──────────────────────────────────────────────────────
eps_water       = 80.0          # (−) water at 25 °C
eps_lipid       = 2.0           # (−) hydrocarbon core

# ── ATP energetics (physiological) ───────────────────────────────────────────
# [ATP]=5 mM, [ADP]=0.5 mM, [Pi]=10 mM, T=310 K
DG_ATP_kJ_mol   = -50.0        # kJ mol⁻¹  (Lehninger 2017; Kushmerick 1997)
DG_ATP_J        = abs(DG_ATP_kJ_mol) * 1e3 / NA  # J per molecule
DG_ATP_kBT      = DG_ATP_J / kBT  # ~20 kBT

# ── Na⁺/K⁺-ATPase ionic concentrations (physiological, human) ────────────────
Na_out_mM       = 145.0        # mM  extracellular Na⁺
Na_in_mM        = 12.0         # mM  intracellular Na⁺
K_out_mM        = 5.0          # mM  extracellular K⁺  (Attwell & Laughlin 2001)
K_in_mM         = 140.0        # mM  intracellular K⁺
V_rest_mV       = -70.0        # mV  resting membrane potential

# ── Microfluidic parameters ───────────────────────────────────────────────────
rho_water_37    = 993.3         # kg m⁻³  water at 37 °C
mu_water_37     = 6.913e-4      # Pa s    dynamic viscosity at 37 °C
D_ATP_m2s       = 6.0e-10       # m² s⁻¹  ATP free-solution diffusion (~600 µm²/s; Hubley 1996)
D_solute_m2s    = 1.0e-9        # m² s⁻¹  small organic solute (~MW 200 Da) in water at 37 C
v_flow          = 1.0e-3        # m s⁻¹   nominal flow velocity
w_channel_min   = 10.0e-6       # m       minimum viable channel width (10 µm)
L_junction      = 100.0e-6      # m       length of one cascade junction

# ── BPS partition-coordinate capacity ────────────────────────────────────────
def partition_capacity(n: int) -> int:
    """C(n) = 2n²  — number of states at partition depth n."""
    return 2 * n * n
