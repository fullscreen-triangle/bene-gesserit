//! Physical and biological constants for membrane dynamics simulation
//! 
//! This module contains fundamental constants used throughout the membrane
//! dynamics system, including physical constants, biological parameters,
//! and experimentally determined values.

/// Fundamental physical constants
pub mod physical {
    /// Gas constant (J/(mol·K))
    pub const GAS_CONSTANT: f64 = 8.314462618;
    
    /// Faraday constant (C/mol)
    pub const FARADAY_CONSTANT: f64 = 96485.33212;
    
    /// Avogadro's number (1/mol)
    pub const AVOGADRO_NUMBER: f64 = 6.02214076e23;
    
    /// Boltzmann constant (J/K)
    pub const BOLTZMANN_CONSTANT: f64 = 1.380649e-23;
    
    /// Elementary charge (C)
    pub const ELEMENTARY_CHARGE: f64 = 1.602176634e-19;
    
    /// Permittivity of free space (F/m)
    pub const EPSILON_0: f64 = 8.8541878128e-12;
    
    /// Relative permittivity of water
    pub const EPSILON_WATER: f64 = 80.1;
    
    /// Relative permittivity of lipid membrane
    pub const EPSILON_MEMBRANE: f64 = 2.1;
}

/// Biological membrane constants
pub mod membrane {
    /// Typical membrane thickness (m)
    pub const MEMBRANE_THICKNESS: f64 = 4e-9;  // 4 nm
    
    /// Typical membrane capacitance per unit area (F/m²)
    pub const SPECIFIC_CAPACITANCE: f64 = 1e-2;  // 1 μF/cm²
    
    /// Typical lipid area per molecule (m²)
    pub const LIPID_AREA: f64 = 0.7e-18;  // 0.7 nm²
    
    /// Water permeability coefficient (m/s)
    pub const WATER_PERMEABILITY: f64 = 2.4e-14;
    
    /// Typical membrane resistance (Ω·m²)
    pub const MEMBRANE_RESISTANCE: f64 = 1e6;
}

/// ATP and energy metabolism constants
pub mod atp {
    /// Standard free energy of ATP hydrolysis (J/mol)
    pub const ATP_HYDROLYSIS_ENERGY: f64 = -30500.0;  // -30.5 kJ/mol
    
    /// Typical cellular ATP concentration (M)
    pub const PHYSIOLOGICAL_ATP: f64 = 5e-3;  // 5 mM
    
    /// Typical cellular ADP concentration (M)
    pub const PHYSIOLOGICAL_ADP: f64 = 0.5e-3;  // 0.5 mM
    
    /// Typical cellular Pi concentration (M)
    pub const PHYSIOLOGICAL_PI: f64 = 1e-3;  // 1 mM
    
    /// ATP molecular weight (g/mol)
    pub const ATP_MOLECULAR_WEIGHT: f64 = 507.18;
    
    /// Energy per ATP molecule (J)
    pub const ATP_ENERGY_PER_MOLECULE: f64 = ATP_HYDROLYSIS_ENERGY / super::physical::AVOGADRO_NUMBER;
}

/// Ion channel and transport constants
pub mod transport {
    /// Na⁺/K⁺-ATPase stoichiometry (Na⁺ out per cycle)
    pub const NA_K_ATPASE_NA_STOICHIOMETRY: u8 = 3;
    
    /// Na⁺/K⁺-ATPase stoichiometry (K⁺ in per cycle)
    pub const NA_K_ATPASE_K_STOICHIOMETRY: u8 = 2;
    
    /// Ca²⁺-ATPase stoichiometry (Ca²⁺ out per cycle)
    pub const CA_ATPASE_STOICHIOMETRY: u8 = 1;
    
    /// Typical Na⁺/K⁺-ATPase turnover rate (cycles/s)
    pub const NA_K_ATPASE_TURNOVER: f64 = 100.0;
    
    /// Typical Ca²⁺-ATPase turnover rate (cycles/s)
    pub const CA_ATPASE_TURNOVER: f64 = 50.0;
    
    /// Voltage-gated channel activation energy (J/mol)
    pub const CHANNEL_ACTIVATION_ENERGY: f64 = 50000.0;  // 50 kJ/mol
}

/// Temperature-related constants
pub mod temperature {
    /// Physiological temperature (K)
    pub const PHYSIOLOGICAL: f64 = 310.15;  // 37°C
    
    /// Room temperature (K)
    pub const ROOM: f64 = 298.15;  // 25°C
    
    /// Minimum viable temperature for membrane function (K)
    pub const MIN_VIABLE: f64 = 273.15;  // 0°C
    
    /// Maximum viable temperature for membrane function (K)
    pub const MAX_VIABLE: f64 = 323.15;  // 50°C
    
    /// Q10 factor for biological processes
    pub const Q10_BIOLOGICAL: f64 = 2.0;
}

/// Simulation and numerical constants
pub mod simulation {
    /// Default simulation timestep (s)
    pub const DEFAULT_TIMESTEP: f64 = 1e-6;  // 1 μs
    
    /// Maximum stable timestep for membrane dynamics (s)
    pub const MAX_STABLE_TIMESTEP: f64 = 1e-4;  // 100 μs
    
    /// Minimum timestep for numerical stability (s)
    pub const MIN_TIMESTEP: f64 = 1e-9;  // 1 ns
    
    /// Convergence tolerance for iterative solvers
    pub const CONVERGENCE_TOLERANCE: f64 = 1e-12;
    
    /// Maximum iterations for convergence
    pub const MAX_ITERATIONS: usize = 1000;
    
    /// Default membrane patch area (m²)
    pub const DEFAULT_PATCH_AREA: f64 = 1e-9;  // 1 μm²
}

/// Physiological ion concentrations (M)
pub mod ion_concentrations {
    /// Intracellular Na⁺ concentration
    pub const NA_INSIDE: f64 = 15e-3;  // 15 mM
    
    /// Extracellular Na⁺ concentration
    pub const NA_OUTSIDE: f64 = 145e-3;  // 145 mM
    
    /// Intracellular K⁺ concentration
    pub const K_INSIDE: f64 = 140e-3;  // 140 mM
    
    /// Extracellular K⁺ concentration
    pub const K_OUTSIDE: f64 = 5e-3;  // 5 mM
    
    /// Intracellular Ca²⁺ concentration
    pub const CA_INSIDE: f64 = 100e-9;  // 100 nM
    
    /// Extracellular Ca²⁺ concentration
    pub const CA_OUTSIDE: f64 = 2e-3;  // 2 mM
    
    /// Intracellular Cl⁻ concentration
    pub const CL_INSIDE: f64 = 10e-3;  // 10 mM
    
    /// Extracellular Cl⁻ concentration
    pub const CL_OUTSIDE: f64 = 110e-3;  // 110 mM
}

/// Lipid properties and membrane composition
pub mod lipids {
    /// POPC phase transition temperature (K)
    pub const POPC_TRANSITION_TEMP: f64 = 271.15;  // -2°C
    
    /// POPE phase transition temperature (K)
    pub const POPE_TRANSITION_TEMP: f64 = 298.15;  // 25°C
    
    /// Cholesterol effect on membrane fluidity
    pub const CHOLESTEROL_FLUIDITY_FACTOR: f64 = 0.3;
    
    /// Typical lipid flip-flop rate (1/s)
    pub const LIPID_FLIP_RATE: f64 = 1e-5;  // Very slow
    
    /// Lateral diffusion coefficient in membrane (m²/s)
    pub const LATERAL_DIFFUSION: f64 = 1e-12;  // 1 μm²/s
}

/// Protein densities and properties
pub mod proteins {
    /// Typical Na⁺/K⁺-ATPase density (proteins/μm²)
    pub const NA_K_ATPASE_DENSITY: f64 = 1000.0;
    
    /// Typical Ca²⁺-ATPase density (proteins/μm²)
    pub const CA_ATPASE_DENSITY: f64 = 500.0;
    
    /// Typical ion channel density (channels/μm²)
    pub const ION_CHANNEL_DENSITY: f64 = 100.0;
    
    /// Protein molecular weight range (g/mol)
    pub const TYPICAL_PROTEIN_MW: f64 = 50000.0;  // 50 kDa
}

// Re-export commonly used constants at module level
pub use physical::{GAS_CONSTANT, FARADAY_CONSTANT, AVOGADRO_NUMBER};
pub use atp::{ATP_HYDROLYSIS_ENERGY, PHYSIOLOGICAL_ATP};
pub use temperature::PHYSIOLOGICAL as PHYSIOLOGICAL_TEMPERATURE;
pub use membrane::{MEMBRANE_THICKNESS, SPECIFIC_CAPACITANCE};

/// Calculate thermal voltage (RT/F) at given temperature
pub fn thermal_voltage(temperature: f64) -> f64 {
    (GAS_CONSTANT * temperature) / FARADAY_CONSTANT
}

/// Calculate Q10 temperature scaling factor
pub fn q10_factor(temp1: f64, temp2: f64, q10: f64) -> f64 {
    q10.powf((temp2 - temp1) / 10.0)
}

/// Validate temperature is within viable range
pub fn is_viable_temperature(temperature: f64) -> bool {
    temperature >= temperature::MIN_VIABLE && temperature <= temperature::MAX_VIABLE
}

/// Calculate membrane capacitance for given area
pub fn membrane_capacitance(area: f64) -> f64 {
    area * SPECIFIC_CAPACITANCE
}

/// Calculate number of molecules from concentration and volume
pub fn molecules_from_concentration(concentration: f64, volume: f64) -> u64 {
    (concentration * volume * AVOGADRO_NUMBER) as u64
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_thermal_voltage() {
        let vt = thermal_voltage(temperature::PHYSIOLOGICAL);
        // At 37°C, RT/F ≈ 26.7 mV
        assert!((vt - 0.0267).abs() < 0.001);
    }
    
    #[test]
    fn test_q10_factor() {
        let factor = q10_factor(298.15, 308.15, 2.0);  // 25°C to 35°C with Q10=2
        assert!((factor - 2.0).abs() < 0.1);
    }
    
    #[test]
    fn test_temperature_validation() {
        assert!(is_viable_temperature(temperature::PHYSIOLOGICAL));
        assert!(!is_viable_temperature(200.0));  // Too cold
        assert!(!is_viable_temperature(400.0));  // Too hot
    }
    
    #[test]
    fn test_membrane_capacitance() {
        let area = 1e-9;  // 1 μm²
        let capacitance = membrane_capacitance(area);
        assert!((capacitance - 1e-11).abs() < 1e-12);  // 10 pF
    }
    
    #[test]
    fn test_molecule_counting() {
        let conc = 1e-3;  // 1 mM
        let volume = 1e-15;  // 1 fL
        let molecules = molecules_from_concentration(conc, volume);
        assert!(molecules > 0);
    }
} 