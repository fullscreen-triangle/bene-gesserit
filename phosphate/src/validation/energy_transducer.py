"""
Validation module — Paper 4: Energy Transducer
'Energy Transduction below the Landauer Floor:
 ATP Coupling and the dx/dATP Identity'
(phosphate/docs/energy-transducer/energy-transducer.tex)

Seven claims validated:
  et_01  Landauer floor ΔE_min = k_B T ln 2 at 310 K
  et_02  ATP free energy ΔG_ATP ≈ 20 k_B T at physiological conditions
  et_03  dx/dATP = ΔG_ATP / (k_B T ln 3) ≈ 18.2 ternary steps per ATP
  et_04  Electrochemical gradient ΔμNa per ion ≈ 5.11 k_B T (outward)
  et_05  Electrochemical gradient ΔμK  per ion ≈ 0.71 k_B T (inward)
  et_06  Na⁺/K⁺-ATPase efficiency η = (3ΔμNa + 2ΔμK) / ΔG_ATP ≈ 84%
  et_07  ATP diffusion: characteristic time < 100 ms across a 10-µm cell
"""

import math
from .constants import (
    kB, kBT, e_charge, R_gas, NA,
    ln2, ln3,
    T_physio,
    DG_ATP_kJ_mol, DG_ATP_kBT, DG_ATP_J,
    Na_out_mM, Na_in_mM, K_out_mM, K_in_mM, V_rest_mV,
    D_ATP_m2s,
)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _result(claim_id, description, predicted, expected, units,
            tolerance, tolerance_type, details=None):
    err = abs(predicted - expected)
    rel = err / abs(expected) if expected != 0 else float("inf")
    if tolerance_type == "relative":
        passed = rel <= tolerance
    elif tolerance_type == "absolute":
        passed = err <= tolerance
    elif tolerance_type == "upper_bound":
        passed = predicted <= expected
    else:
        passed = bool(predicted)
    return {
        "claim_id":       claim_id,
        "paper":          "energy_transducer",
        "description":    description,
        "predicted":      round(predicted, 6),
        "expected":       round(expected, 6),
        "units":          units,
        "absolute_error": round(err, 8),
        "relative_error": round(rel, 6),
        "tolerance":      tolerance,
        "tolerance_type": tolerance_type,
        "passed":         bool(passed),
        "details":        details or {},
    }


# ── et_01  Landauer floor ─────────────────────────────────────────────────────

def validate_landauer_floor() -> dict:
    """
    Minimum energy per irreversible binary operation (Landauer 1961):
        ΔE_min = k_B T ln 2
    Computed at T = 310 K (physiological).
    """
    T           = T_physio
    DeltaE_J    = kB * T * ln2
    DeltaE_kBT  = DeltaE_J / kBT     # should equal ln2 = 0.6931...

    return _result(
        "et_01",
        "Landauer floor ΔE_min = k_BT ln2 at 310 K",
        predicted=DeltaE_kBT,
        expected=ln2,
        units="kBT",
        tolerance=1e-10,
        tolerance_type="absolute",
        details={
            "T_K":          T,
            "kB_J_K":       kB,
            "DeltaE_J":     DeltaE_J,
            "DeltaE_kBT":   round(DeltaE_kBT, 8),
            "ln2":          round(ln2, 8),
            "formula":      "ΔE_min = k_B T ln 2  [Landauer 1961]",
        },
    )


# ── et_02  ATP free energy ────────────────────────────────────────────────────

def validate_atp_free_energy() -> dict:
    """
    Physiological free energy of ATP hydrolysis:
        ΔG_ATP ≈ −50 kJ mol⁻¹  →  ≈ 20 k_B T
    at [ATP]=5 mM, [ADP]=0.5 mM, [Pi]=10 mM, T=310 K.
    The claimed value is 20 k_B T; tolerance ±2 k_B T.
    """
    # Direct conversion from kJ/mol to kBT
    kBT_kJ_mol = R_gas * T_physio / 1000.0  # kJ mol⁻¹ per kBT
    DG_kBT     = abs(DG_ATP_kJ_mol) / kBT_kJ_mol

    return _result(
        "et_02",
        "ATP free energy |ΔG_ATP| ≈ 20 k_BT (physiological conditions)",
        predicted=DG_kBT,
        expected=20.0,
        units="kBT",
        tolerance=0.10,  # 10% relative tolerance
        tolerance_type="relative",
        details={
            "DG_kJ_mol":      DG_ATP_kJ_mol,
            "kBT_kJ_mol":     round(kBT_kJ_mol, 5),
            "DG_kBT":         round(DG_kBT, 4),
            "conditions":     "[ATP]=5mM, [ADP]=0.5mM, [Pi]=10mM, T=310K",
            "reference":      "Lehninger 2017; Kushmerick 1997",
        },
    )


# ── et_03  dx/dATP identity ───────────────────────────────────────────────────

def validate_dx_datp() -> dict:
    """
    Each ATP hydrolysis event drives dx/dATP ternary addressing steps:
        dx/dATP = |ΔG_ATP| / (k_B T ln 3)
    At |ΔG_ATP| = 20 k_B T:  dx/dATP = 20 / ln3 ≈ 18.21.
    """
    DG_kBT      = abs(DG_ATP_kJ_mol) / (R_gas * T_physio / 1000.0)
    dx_datp     = DG_kBT / ln3

    return _result(
        "et_03",
        "dx/dATP = |ΔG_ATP| / (k_BT ln 3) ≈ 18.2 ternary steps per ATP",
        predicted=dx_datp,
        expected=18.2,
        units="ternary steps / ATP",
        tolerance=0.05,
        tolerance_type="relative",
        details={
            "DG_kBT":       round(DG_kBT, 4),
            "ln3":          round(ln3, 6),
            "dx_datp":      round(dx_datp, 4),
            "bits_per_ATP": round(dx_datp * math.log2(3), 3),
            "formula":      "dx/dATP = |ΔG_ATP| / (k_BT ln 3)  [Paper 4 §5]",
        },
    )


# ── et_04  ΔμNa per ion ───────────────────────────────────────────────────────

def validate_delta_mu_na() -> dict:
    """
    Electrochemical potential driving force per Na⁺ ion expelled outward:
        Δμ_Na = k_B T ln([Na]_out / [Na]_in) + z·e·V_rest
    At physiological concentrations and V_rest = −70 mV.
    """
    z        = 1                              # Na⁺ valence
    ratio_Na = Na_out_mM / Na_in_mM          # 145/12 = 12.08
    V_rest_V = V_rest_mV * 1e-3              # −0.070 V

    # Chemical part: k_BT ln([Na]_out/[Na]_in)
    mu_chem_kBT  = math.log(ratio_Na)        # positive (outward concentration gradient)

    # Electrical part: z·e·V_rest / k_BT  (outward = against −70 mV)
    # For Na⁺ moving out (from −70 mV inside to 0 mV outside), ΔV = +70 mV
    mu_elec_kBT  = z * e_charge * abs(V_rest_V) / kBT

    delta_mu_Na_kBT = mu_chem_kBT + mu_elec_kBT

    return _result(
        "et_04",
        "Electrochemical gradient ΔμNa per ion ≈ 5.11 k_BT (outward)",
        predicted=delta_mu_Na_kBT,
        expected=5.11,
        units="kBT",
        tolerance=0.10,
        tolerance_type="absolute",
        details={
            "Na_out_mM":       Na_out_mM,
            "Na_in_mM":        Na_in_mM,
            "ratio_Na":        round(ratio_Na, 4),
            "mu_chem_kBT":     round(mu_chem_kBT, 4),
            "V_rest_mV":       V_rest_mV,
            "mu_elec_kBT":     round(mu_elec_kBT, 4),
            "delta_mu_Na_kBT": round(delta_mu_Na_kBT, 4),
            "formula":         "Δμ_Na = k_BT ln([Na]out/[Na]in) + z·e·|V_rest| / k_BT",
        },
    )


# ── et_05  ΔμK per ion ────────────────────────────────────────────────────────

def validate_delta_mu_k() -> dict:
    """
    Electrochemical driving force per K⁺ ion pumped inward:
        Δμ_K = k_B T ln([K]_in / [K]_out) − z·e·|V_rest|
    The K⁺ gradient partially opposes the electrical driving force.
    """
    z       = 1
    ratio_K = K_in_mM / K_out_mM           # 140/4 = 35
    V_rest_V = abs(V_rest_mV) * 1e-3       # 0.070 V

    # K⁺ moves inward: chemical gradient drives inward (Δμ_chem > 0)
    mu_chem_kBT  = math.log(ratio_K)       # positive

    # Electrical: K⁺ moving inward gains energy from −70 mV gradient
    mu_elec_kBT  = z * e_charge * V_rest_V / kBT  # positive

    # Net driving force (inward): chem − elec (electric opposes inward flow at −70 mV)
    # At Nernst equilibrium, Δμ_K = 0; resting potential is near K⁺ Nernst potential
    # True driving force per pump cycle = deficit from Nernst:
    V_K_nernst_V = (kBT / (z * e_charge)) * math.log(ratio_K)  # ~+93 mV (outward)
    # At V_rest = −70 mV, pump does work to maintain [K]_in > [K]_out against −70 mV
    delta_mu_K_kBT = abs(mu_chem_kBT - mu_elec_kBT)

    return _result(
        "et_05",
        "Electrochemical gradient ΔμK per ion ≈ 0.71 k_BT (inward, net)",
        predicted=delta_mu_K_kBT,
        expected=0.71,
        units="kBT",
        tolerance=0.15,
        tolerance_type="absolute",
        details={
            "K_in_mM":           K_in_mM,
            "K_out_mM":          K_out_mM,
            "ratio_K":           round(ratio_K, 2),
            "mu_chem_kBT":       round(mu_chem_kBT, 4),
            "V_rest_mV":         V_rest_mV,
            "mu_elec_kBT":       round(mu_elec_kBT, 4),
            "V_K_nernst_mV":     round(V_K_nernst_V * 1000, 1),
            "delta_mu_K_kBT":    round(delta_mu_K_kBT, 4),
            "formula":           "Δμ_K = |k_BT ln([K]in/[K]out) − z·e·|V_rest|/k_BT|",
        },
    )


# ── et_06  Pump efficiency ────────────────────────────────────────────────────

def validate_pump_efficiency() -> dict:
    """
    Na⁺/K⁺-ATPase pumps 3 Na⁺ out and 2 K⁺ in per ATP.
    Efficiency η = (3·Δμ_Na + 2·Δμ_K) / |ΔG_ATP|.
    Paper claims η ≈ 84%.
    """
    # Recompute gradients from first principles
    z        = 1
    V_rest_V = abs(V_rest_mV) * 1e-3

    mu_Na    = math.log(Na_out_mM / Na_in_mM) + z * e_charge * V_rest_V / kBT
    mu_K_raw = math.log(K_in_mM / K_out_mM)
    mu_elec  = z * e_charge * V_rest_V / kBT
    mu_K     = abs(mu_K_raw - mu_elec)

    DG_pump_kBT = 3.0 * mu_Na + 2.0 * mu_K
    DG_ATP_kBT_ = abs(DG_ATP_kJ_mol) / (R_gas * T_physio / 1000.0)
    eta         = DG_pump_kBT / DG_ATP_kBT_

    return _result(
        "et_06",
        "Na⁺/K⁺-ATPase efficiency η = (3ΔμNa + 2ΔμK) / |ΔG_ATP| ≈ 84%",
        predicted=eta,
        expected=0.84,
        units="dimensionless",
        tolerance=0.05,
        tolerance_type="absolute",
        details={
            "mu_Na_kBT":       round(mu_Na, 4),
            "mu_K_kBT":        round(mu_K, 4),
            "DG_pump_kBT":     round(DG_pump_kBT, 4),
            "DG_ATP_kBT":      round(DG_ATP_kBT_, 4),
            "eta":             round(eta, 4),
            "stoichiometry":   "3 Na⁺ out, 2 K⁺ in, 1 ATP [Skou 1957; Albers 1967]",
        },
    )


# ── et_07  ATP diffusion time ─────────────────────────────────────────────────

def validate_atp_diffusion() -> dict:
    """
    Characteristic diffusion time for ATP across a cell of diameter L = 10 µm:
        τ_diff = L² / (2D_ATP)
    with D_ATP ≈ 300 µm² s⁻¹ = 3×10⁻¹⁰ m² s⁻¹ (Hubley 1996).
    Must be < 100 ms (fast relative to signalling timescale).
    """
    L_m       = 10.0e-6         # 10 µm cell diameter
    D_ATP     = D_ATP_m2s       # 3×10⁻¹⁰ m² s⁻¹
    tau_s     = L_m ** 2 / (2 * D_ATP)
    tau_ms    = tau_s * 1000.0

    return _result(
        "et_07",
        "ATP diffusion time across 10-µm cell < 100 ms",
        predicted=tau_ms,
        expected=100.0,
        units="ms",
        tolerance=100.0,
        tolerance_type="upper_bound",
        details={
            "L_m":          L_m,
            "D_ATP_m2s":    D_ATP_m2s,
            "tau_s":        round(tau_s, 6),
            "tau_ms":       round(tau_ms, 4),
            "criterion":    "τ_diff < 100 ms  (fast supply to pump sites)",
            "reference":    "Hubley et al. 1996",
        },
    )


# ── Public entry point ────────────────────────────────────────────────────────

def run_all() -> list:
    return [
        validate_landauer_floor(),
        validate_atp_free_energy(),
        validate_dx_datp(),
        validate_delta_mu_na(),
        validate_delta_mu_k(),
        validate_pump_efficiency(),
        validate_atp_diffusion(),
    ]
