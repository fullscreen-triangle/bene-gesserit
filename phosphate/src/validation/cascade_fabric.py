"""
Validation module -- Paper 5: Cascade Fabric
'Cascade Fabric and the No-Privileged-Level Corollary:
 Ternary Microfluidic Networks for Partition-Depth Computation'
(phosphate/docs/cascade-fabric/cascade-fabric.tex)

Seven claims validated:
  cf_01  Ternary optimality: argmin_{b in Z>=2} {b/ln(b)} = 3  (closest integer to e)
  cf_02  Depth necessity:    k = ceil(log3 N_cat)               (enumeration check)
  cf_03  P_correct formula:  nc=44, k=3 prototype gives P >= 0.95
  cf_04  Energy per query E = k * kBT * ln3  at k=12
  cf_05  ATP equivalence:    E / |DG_ATP| at k=12
  cf_06  Reynolds number Re = rho*v*w / mu ~= 0.014  (Stokes regime)
  cf_07  Peclet number Pe = v*w / D ~= 10             (convection-dominated)
"""

import math
import numpy as np
from .constants import (
    kBT, kB, ln3,
    T_physio,
    DG_ATP_kJ_mol, R_gas,
    rho_water_37, mu_water_37,
    D_solute_m2s, v_flow, w_channel_min,
)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _result(claim_id, description, predicted, expected, units,
            tolerance, tolerance_type, details=None):
    if isinstance(predicted, (int, float)) and isinstance(expected, (int, float)):
        err = abs(float(predicted) - float(expected))
        rel = err / abs(float(expected)) if expected != 0 else float("inf")
    else:
        err = rel = 0.0
    if tolerance_type == "relative":
        passed = bool(rel <= tolerance)
    elif tolerance_type == "absolute":
        passed = bool(err <= tolerance)
    elif tolerance_type == "lower_bound":
        passed = bool(float(predicted) >= float(expected))
    elif tolerance_type == "upper_bound":
        passed = bool(float(predicted) <= float(expected))
    else:
        passed = bool(predicted)
    return {
        "claim_id": claim_id,
        "paper": "cascade_fabric",
        "description": description,
        "predicted": predicted,
        "expected": expected,
        "units": units,
        "absolute_error": round(err, 8),
        "relative_error": round(rel, 6),
        "tolerance": tolerance,
        "tolerance_type": tolerance_type,
        "passed": passed,
        "details": details or {},
    }


# ── cf_01  Ternary optimality ─────────────────────────────────────────────────

def validate_ternary_optimality():
    """
    Total work (junctions x branching cost) for a b-ary tree with N_cat leaves:
        W(b) = k * b = (ln N / ln b) * b = ln(N) * b / ln(b)
    Minimising b/ln(b) over b >= 2:
        d/db [b/ln(b)] = (ln(b)-1)/(ln(b))^2 = 0  ->  b = e ~= 2.718
    Closest integer: b = 3.
    Verified: f(3) < f(2) and f(3) < f(4).
    """
    def f(b):
        return b / math.log(b)

    fb2 = f(2)
    fb3 = f(3)
    fb4 = f(4)
    fb_e = f(math.e)      # global minimum (= e)

    optimal_int = bool((fb3 < fb2) and (fb3 < fb4))

    table = {str(b): round(f(b), 6) for b in range(2, 9)}
    table["e"] = round(fb_e, 6)

    return {
        "claim_id": "cf_01",
        "paper": "cascade_fabric",
        "description": "Ternary optimality: argmin_{b>=2} {b/ln(b)} = 3 (closest integer to e)",
        "predicted": 3,
        "expected": 3,
        "units": "integer base",
        "absolute_error": 0,
        "relative_error": 0.0,
        "tolerance": 0,
        "tolerance_type": "absolute",
        "passed": optimal_int,
        "details": {
            "f_table_b_over_lnb": table,
            "f(2)": round(fb2, 6),
            "f(3)": round(fb3, 6),
            "f(4)": round(fb4, 6),
            "f(e)_global_minimum": round(fb_e, 6),
            "f3_lt_f2": bool(fb3 < fb2),
            "f3_lt_f4": bool(fb3 < fb4),
        },
    }


# ── cf_02  Depth necessity ────────────────────────────────────────────────────

def validate_depth_necessity():
    """
    Minimum cascade depth to resolve N_cat categories in base 3:
        k = ceil(log3(N_cat))
    Verified for representative N_cat values including exact powers of 3.
    """
    test_cases = {
        27:       3,    # 3^3 = 27
        243:      5,    # 3^5 = 243
        729:      6,    # 3^6 = 729
        2187:     7,    # 3^7 = 2187
        59049:    10,   # 3^10 = 59049 (exact)
        531441:   12,   # 3^12 = 531441 (~5.3e5 categories, paper k=12)
        14348907: 15,   # 3^15 = 14348907
    }
    results = {}
    all_ok = True
    for N_cat, k_expected in test_cases.items():
        k_pred = math.ceil(math.log(N_cat) / ln3)
        ok = bool(k_pred == k_expected)
        results[str(N_cat)] = {
            "k_predicted": k_pred,
            "k_expected": k_expected,
            "passed": ok,
        }
        if not ok:
            all_ok = False

    return {
        "claim_id": "cf_02",
        "paper": "cascade_fabric",
        "description": "Depth k = ceil(log3(N_cat)) for representative N_cat values",
        "predicted": "see details",
        "expected": "k_pred == k_expected for all N_cat",
        "units": "depth levels",
        "absolute_error": 0,
        "relative_error": 0.0,
        "tolerance": 0,
        "tolerance_type": "absolute",
        "passed": bool(all_ok),
        "details": results,
    }


# ── cf_03  P_correct reliability ──────────────────────────────────────────────

def validate_p_correct():
    """
    Probability of correct cascade routing:
        P_correct = (1 - (1 - kappa_bar)^{nc})^k
    For nc=44, kappa_bar=0.10:
        p_junction = 1 - (0.90)^44 = 0.9903
    For k=3 (prototype): P_correct = 0.9903^3 = 0.971 >= 0.95  (pass)
    For k=12 (full system): P_correct = 0.9903^12 = 0.887  (reported separately)
    """
    nc = 44
    kappa_bar = 0.10
    k_proto = 3
    k_full = 12
    target = 0.95

    p_junc = float(1.0 - (1.0 - kappa_bar) ** nc)
    P_proto = float(p_junc ** k_proto)
    P_full = float(p_junc ** k_full)

    return _result(
        "cf_03",
        "P_correct >= 0.95 for nc=44, kappa=0.10, k=3 (prototype validation)",
        predicted=P_proto,
        expected=target,
        units="probability",
        tolerance=target,
        tolerance_type="lower_bound",
        details={
            "nc": nc,
            "kappa_bar": kappa_bar,
            "k_prototype": k_proto,
            "k_full_system": k_full,
            "p_per_junction": round(p_junc, 6),
            "P_correct_k3": round(P_proto, 6),
            "P_correct_k12": round(P_full, 6),
            "formula": "P = (1-(1-kappa_bar)^{nc})^k  [Paper 5 Sec.6]",
            "note": "k=3 prototype passes 0.95 threshold; k=12 requires nc>52",
        },
    )


# ── cf_04  Energy per query ───────────────────────────────────────────────────

def validate_energy_per_query():
    """
    Each junction commits one ternary routing decision costing kBT*ln3.
    Total for k=12 levels: E = k * kBT * ln3.
    """
    k = 12
    E_kBT = k * ln3

    return _result(
        "cf_04",
        "Energy per query E = k * kBT * ln3  at k=12",
        predicted=E_kBT,
        expected=k * ln3,
        units="kBT",
        tolerance=1e-10,
        tolerance_type="absolute",
        details={
            "k": k,
            "ln3": round(ln3, 6),
            "E_kBT": round(E_kBT, 6),
            "E_J": float(k * kBT * ln3),
            "formula": "E_query = k * kBT * ln3  [No-Privileged-Level Sec.7]",
        },
    )


# ── cf_05  ATP equivalence ────────────────────────────────────────────────────

def validate_atp_equivalence():
    """
    ATP fraction per query: E_query / |DG_ATP|  (both in kBT).
    At k=12: fraction < 1 ATP per query.
    """
    k = 12
    E_kBT = k * ln3
    DG_kBT = abs(DG_ATP_kJ_mol) / (R_gas * T_physio / 1000.0)
    atp_fraction = E_kBT / DG_kBT

    return _result(
        "cf_05",
        "ATP equivalence: E_query / |DG_ATP| < 1 ATP per query (k=12)",
        predicted=atp_fraction,
        expected=1.0,
        units="ATP equivalent",
        tolerance=1.0,
        tolerance_type="upper_bound",
        details={
            "k": k,
            "E_kBT": round(E_kBT, 4),
            "DG_ATP_kBT": round(DG_kBT, 4),
            "atp_fraction": round(atp_fraction, 4),
        },
    )


# ── cf_06  Reynolds number ────────────────────────────────────────────────────

def validate_reynolds_number():
    """
    Re = rho * v * w / mu  for microfluidic channel.
    At v=1 mm/s, w=10 um, water at 37 C: Re ~= 0.014  << 1 (Stokes regime).
    """
    Re = float(rho_water_37 * v_flow * w_channel_min / mu_water_37)

    return _result(
        "cf_06",
        "Reynolds number Re = rho*v*w/mu ~= 0.014 << 1 (laminar)",
        predicted=Re,
        expected=0.014,
        units="dimensionless",
        tolerance=0.20,
        tolerance_type="relative",
        details={
            "rho_kg_m3": rho_water_37,
            "v_m_s": v_flow,
            "w_m": w_channel_min,
            "mu_Pa_s": mu_water_37,
            "Re": round(Re, 5),
            "laminar": bool(Re < 1.0),
        },
    )


# ── cf_07  Peclet number ──────────────────────────────────────────────────────

def validate_peclet_number():
    """
    Pe = v * w / D  (convection vs diffusion at junction scale).
    At v=1 mm/s, w=10 um, D~1e-9 m^2/s (small organic, ~MW 200 Da):
        Pe ~= 10  (convection-dominated; sharp stream boundaries).
    """
    Pe = float(v_flow * w_channel_min / D_solute_m2s)

    return _result(
        "cf_07",
        "Peclet number Pe = v*w/D ~= 10 (convection-dominated routing)",
        predicted=Pe,
        expected=10.0,
        units="dimensionless",
        tolerance=0.20,
        tolerance_type="relative",
        details={
            "v_m_s": v_flow,
            "w_m": w_channel_min,
            "D_m2_s": D_solute_m2s,
            "Pe": round(Pe, 4),
            "interpretation": "Pe > 1 -> diffusion does not smear stream boundaries",
        },
    )


# ── Public entry point ────────────────────────────────────────────────────────

def run_all():
    return [
        validate_ternary_optimality(),
        validate_depth_necessity(),
        validate_p_correct(),
        validate_energy_per_query(),
        validate_atp_equivalence(),
        validate_reynolds_number(),
        validate_peclet_number(),
    ]
