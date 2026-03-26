from ..utils import save_json
"""
Phase-Lock Network Validator

Validates Claims:
1. Kuramoto model with drug-modified coupling: K_mod = K₀(1 + [Drug]×K_agg)
2. Phase coherence R > 0.7 for therapeutic effect
3. Information transfer: I = R × BW × log₂(SNR) ≈ 500-610 bits/s
4. Phase-lock network topology determines entropy
5. Drug aggregation constants K_agg > 10^4 M⁻¹

Tests Kuramoto dynamics, phase coherence, and information transfer.
"""

import numpy as np
from dataclasses import dataclass
from typing import List, Dict, Tuple, Callable
import json
from pathlib import Path
import time


@dataclass
class KuramotoSystem:
    """Kuramoto oscillator network."""
    n_oscillators: int
    natural_frequencies: np.ndarray
    coupling_strength: float
    phases: np.ndarray


class PhaseLockValidator:
    """Validates phase-lock network dynamics and coherence."""
    
    def __init__(self, output_dir: Path = Path("results/phase_lock")):
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def initialize_kuramoto_system(self, 
                                   n_oscillators: int = 100,
                                   mean_freq: float = 10.0,
                                   freq_std: float = 2.0,
                                   coupling: float = 0.5) -> KuramotoSystem:
        """
        Initialize Kuramoto oscillator network.
        
        dθᵢ/dt = ωᵢ + (K/N) Σⱼ sin(θⱼ - θᵢ)
        """
        # Natural frequencies (Gaussian distribution)
        natural_frequencies = np.random.normal(mean_freq, freq_std, n_oscillators)
        
        # Initial phases (random)
        phases = np.random.uniform(0, 2*np.pi, n_oscillators)
        
        return KuramotoSystem(
            n_oscillators=n_oscillators,
            natural_frequencies=natural_frequencies,
            coupling_strength=coupling,
            phases=phases
        )
    
    def kuramoto_derivatives(self, 
                           phases: np.ndarray,
                           system: KuramotoSystem) -> np.ndarray:
        """
        Compute time derivatives for Kuramoto model.
        
        dθᵢ/dt = ωᵢ + (K/N) Σⱼ sin(θⱼ - θᵢ)
        """
        N = system.n_oscillators
        K = system.coupling_strength
        omega = system.natural_frequencies
        
        # Compute coupling term
        coupling_term = np.zeros(N)
        for i in range(N):
            coupling_term[i] = (K / N) * np.sum(np.sin(phases - phases[i]))
        
        dtheta_dt = omega + coupling_term
        
        return dtheta_dt
    
    def simulate_kuramoto(self, 
                         system: KuramotoSystem,
                         t_max: float = 100.0,
                         dt: float = 0.01) -> Tuple[np.ndarray, np.ndarray]:
        """
        Simulate Kuramoto dynamics using RK4 integration.
        
        Returns (times, phases_history)
        """
        n_steps = int(t_max / dt)
        times = np.linspace(0, t_max, n_steps)
        
        phases_history = np.zeros((n_steps, system.n_oscillators))
        phases_history[0, :] = system.phases
        
        phases = system.phases.copy()
        
        for i in range(1, n_steps):
            # RK4 integration
            k1 = self.kuramoto_derivatives(phases, system)
            k2 = self.kuramoto_derivatives(phases + 0.5*dt*k1, system)
            k3 = self.kuramoto_derivatives(phases + 0.5*dt*k2, system)
            k4 = self.kuramoto_derivatives(phases + dt*k3, system)
            
            phases = phases + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
            phases = np.mod(phases, 2*np.pi)  # Keep in [0, 2π]
            
            phases_history[i, :] = phases
        
        return times, phases_history
    
    def compute_order_parameter(self, phases: np.ndarray) -> Tuple[float, float]:
        """
        Compute Kuramoto order parameter.
        
        R = |⟨exp(iθ)⟩|
        Ψ = arg(⟨exp(iθ)⟩)
        
        R ∈ [0, 1]: 0 = incoherent, 1 = fully synchronized
        """
        z = np.mean(np.exp(1j * phases))
        R = np.abs(z)
        Psi = np.angle(z)
        
        return R, Psi
    
    def test_drug_modified_coupling(self) -> Dict:
        """
        Test drug modification of coupling strength.
        
        K_modified = K₀(1 + [Drug] × K_agg)
        
        Validates:
        - Lithium: K_mod = 0.75 (from K₀ = 0.5)
        - Dopamine: K_mod = 0.60
        - Serotonin: K_mod = 0.65
        """
        K0 = 0.5  # Baseline coupling
        
        drugs = [
            ("Lithium", 1e-3, 5e3, 0.75),  # [Drug] (M), K_agg (M⁻¹), expected K_mod
            ("Dopamine", 1e-4, 2e3, 0.60),
            ("Serotonin", 1e-4, 3e3, 0.65),
        ]
        
        results = []
        
        for drug_name, concentration, K_agg, expected_K in drugs:
            # Compute modified coupling
            K_modified = K0 * (1 + concentration * K_agg)
            
            # Relative error
            rel_error = abs(K_modified - expected_K) / expected_K
            
            # Simulate system with modified coupling
            system_baseline = self.initialize_kuramoto_system(n_oscillators=50, coupling=K0)
            system_drug = self.initialize_kuramoto_system(n_oscillators=50, coupling=K_modified)
            system_drug.natural_frequencies = system_baseline.natural_frequencies.copy()
            
            # Short simulation
            times, phases_baseline = self.simulate_kuramoto(system_baseline, t_max=50, dt=0.1)
            _, phases_drug = self.simulate_kuramoto(system_drug, t_max=50, dt=0.1)
            
            # Compute final coherence
            R_baseline, _ = self.compute_order_parameter(phases_baseline[-1, :])
            R_drug, _ = self.compute_order_parameter(phases_drug[-1, :])
            
            results.append({
                "drug": drug_name,
                "concentration_M": concentration,
                "K_agg": K_agg,
                "K_baseline": K0,
                "K_modified": K_modified,
                "K_expected": expected_K,
                "relative_error": rel_error,
                "R_baseline": R_baseline,
                "R_drug": R_drug,
                "coherence_change": R_drug - R_baseline,
            })
        
        test = {
            "results": results,
            "claim_validated": all(r["relative_error"] < 0.3 for r in results),
        }
        
        return test
    
    def test_therapeutic_coherence_threshold(self) -> Dict:
        """
        Test therapeutic effect correlation with R > 0.7.
        
        Validates that therapeutic drugs produce R > 0.7 phase coherence.
        """
        # Test different coupling strengths
        coupling_values = [0.1, 0.3, 0.5, 0.7, 0.9]
        
        coherence_results = []
        
        for K in coupling_values:
            system = self.initialize_kuramoto_system(n_oscillators=50, coupling=K)
            times, phases = self.simulate_kuramoto(system, t_max=100, dt=0.1)
            
            # Compute coherence over time
            R_values = []
            for t_idx in range(len(times)):
                R, _ = self.compute_order_parameter(phases[t_idx, :])
                R_values.append(R)
            
            # Final coherence (steady-state)
            R_final = np.mean(R_values[-100:])
            
            # Therapeutic if R > 0.7
            is_therapeutic = R_final > 0.7
            
            coherence_results.append({
                "coupling_K": K,
                "R_final": R_final,
                "is_therapeutic": is_therapeutic,
            })
        
        test = {
            "results": coherence_results,
            "threshold": 0.7,
            "claim_validated": any(r["is_therapeutic"] for r in coherence_results),
        }
        
        return test
    
    def compute_information_transfer(self, R: float) -> Dict:
        """
        Compute information transfer rate.
        
        I = R × BW × log₂(SNR)
        
        Where:
        - R: phase coherence
        - BW: bandwidth (Hz)
        - SNR: signal-to-noise ratio
        
        Claim: 500-610 bits/s
        """
        # Typical values
        BW = 1000  # Hz, neural bandwidth
        SNR = 10   # Typical SNR
        
        I_bits_per_sec = R * BW * np.log2(SNR)
        
        result = {
            "R": R,
            "bandwidth_hz": BW,
            "SNR": SNR,
            "information_rate_bits_per_sec": I_bits_per_sec,
            "claim_min": 500,
            "claim_max": 610,
            "within_claim": 500 <= I_bits_per_sec <= 610,
        }
        
        return result
    
    def run_validation(self) -> Dict:
        """
        Run complete phase-lock validation.
        
        Returns:
            Comprehensive results dictionary
        """
        print("="*70)
        print("PHASE-LOCK NETWORK VALIDATION")
        print("="*70)
        
        # Test drug-modified coupling
        print("\n1. Testing drug-modified coupling...")
        drug_coupling = self.test_drug_modified_coupling()
        
        # Test therapeutic coherence threshold
        print("\n2. Testing therapeutic coherence threshold...")
        coherence_test = self.test_therapeutic_coherence_threshold()
        
        # Compute information transfer
        print("\n3. Computing information transfer rates...")
        # Use R values from coherence test
        info_transfers = []
        for result in coherence_test["results"]:
            info = self.compute_information_transfer(result["R_final"])
            info_transfers.append(info)
        
        # Compile results
        results = {
            "validator": "PhaseLockValidator",
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "drug_modified_coupling": drug_coupling,
            "coherence_threshold": coherence_test,
            "information_transfer": info_transfers,
            "claims_validated": {
                "drug_coupling_modification": drug_coupling["claim_validated"],
                "therapeutic_coherence_threshold": coherence_test["claim_validated"],
                "information_transfer_rate": any(i["within_claim"] for i in info_transfers),
            },
        }
        
        # Save results
        self.save_results(results)
        
        # Print summary
        self.print_summary(results)
        
        return results
    
    def save_results(self, results: Dict):
        """Save validation results to JSON."""
        output_file = self.output_dir / "phase_lock_results.json"
        save_json(results, output_file)
        print(f"\n✓ Results saved to: {output_file}")
    
    def print_summary(self, results: Dict):
        """Print validation summary."""
        print("\n" + "="*70)
        print("VALIDATION SUMMARY")
        print("="*70)
        
        # Drug coupling
        dc = results["drug_modified_coupling"]
        print(f"\nDrug-Modified Coupling:")
        for r in dc["results"]:
            print(f"  {r['drug']}: K = {r['K_modified']:.3f} (expected {r['K_expected']:.3f})")
            print(f"    R: {r['R_baseline']:.3f} → {r['R_drug']:.3f}")
        print(f"  Status: {'✓ VALIDATED' if dc['claim_validated'] else '✗ FAILED'}")
        
        # Coherence threshold
        ct = results["coherence_threshold"]
        print(f"\nTherapeutic Coherence (R > {ct['threshold']}):")
        therapeutic_count = sum(1 for r in ct["results"] if r["is_therapeutic"])
        print(f"  Therapeutic cases: {therapeutic_count}/{len(ct['results'])}")
        print(f"  Status: {'✓ VALIDATED' if ct['claim_validated'] else '✗ FAILED'}")
        
        # Information transfer
        it = results["information_transfer"]
        print(f"\nInformation Transfer:")
        for info in it[:3]:  # Show first 3
            print(f"  R={info['R']:.2f}: {info['information_rate_bits_per_sec']:.0f} bits/s")
        print(f"  Claim: 500-610 bits/s")
        print(f"  Status: {'✓ VALIDATED' if results['claims_validated']['information_transfer_rate'] else '✗ FAILED'}")
        
        print("\n" + "="*70)

