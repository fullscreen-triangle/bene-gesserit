//! Core types for the Bene Gesserit membrane dynamics system
//! 
//! This module defines the fundamental data types used throughout the system,
//! including physical units, membrane states, and biological parameters.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Physical units and measurements
pub mod units {
    /// Membrane area in square meters
    pub type Area = f64;
    
    /// Temperature in Kelvin
    pub type Temperature = f64;
    
    /// Time in seconds
    pub type Time = f64;
    
    /// Voltage in volts
    pub type Voltage = f64;
    
    /// Current in amperes
    pub type Current = f64;
    
    /// Concentration in molar (mol/L)
    pub type Concentration = f64;
    
    /// Energy in joules
    pub type Energy = f64;
    
    /// Power in watts
    pub type Power = f64;
    
    /// Capacitance in farads
    pub type Capacitance = f64;
    
    /// Resistance in ohms
    pub type Resistance = f64;
    
    /// Conductance in siemens
    pub type Conductance = f64;
    
    /// Density (proteins per unit area, etc.)
    pub type Density = f64;
}

/// ATP-related types and measurements
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AtpMeasurement {
    /// ATP concentration in molar
    pub concentration: units::Concentration,
    /// Available ATP molecules (for discrete counting)
    pub molecules: u64,
    /// ATP consumption rate in mol/s
    pub consumption_rate: f64,
    /// ATP production rate in mol/s
    pub production_rate: f64,
}

/// Ion types commonly found in biological membranes
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum IonType {
    Sodium,      // Na⁺
    Potassium,   // K⁺
    Calcium,     // Ca²⁺
    Chloride,    // Cl⁻
    Magnesium,   // Mg²⁺
    Phosphate,   // PO₄³⁻
}

impl IonType {
    /// Get the charge of this ion type
    pub fn charge(&self) -> i8 {
        match self {
            Self::Sodium => 1,
            Self::Potassium => 1,
            Self::Calcium => 2,
            Self::Chloride => -1,
            Self::Magnesium => 2,
            Self::Phosphate => -3,
        }
    }
    
    /// Get the symbol for this ion
    pub fn symbol(&self) -> &'static str {
        match self {
            Self::Sodium => "Na⁺",
            Self::Potassium => "K⁺",
            Self::Calcium => "Ca²⁺",
            Self::Chloride => "Cl⁻",
            Self::Magnesium => "Mg²⁺",
            Self::Phosphate => "PO₄³⁻",
        }
    }
}

/// Ion concentrations on both sides of membrane
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct IonConcentrations {
    pub inside: HashMap<IonType, units::Concentration>,
    pub outside: HashMap<IonType, units::Concentration>,
}

impl IonConcentrations {
    /// Create physiological ion concentrations
    pub fn physiological() -> Self {
        let mut inside = HashMap::new();
        let mut outside = HashMap::new();
        
        // Typical mammalian cell concentrations (mM)
        inside.insert(IonType::Sodium, 0.015);      // 15 mM
        inside.insert(IonType::Potassium, 0.140);   // 140 mM
        inside.insert(IonType::Calcium, 0.0001);    // 0.1 μM
        inside.insert(IonType::Chloride, 0.010);    // 10 mM
        inside.insert(IonType::Magnesium, 0.001);   // 1 mM
        
        outside.insert(IonType::Sodium, 0.145);     // 145 mM
        outside.insert(IonType::Potassium, 0.005);  // 5 mM
        outside.insert(IonType::Calcium, 0.002);    // 2 mM
        outside.insert(IonType::Chloride, 0.110);   // 110 mM
        outside.insert(IonType::Magnesium, 0.002);  // 2 mM
        
        Self { inside, outside }
    }
    
    /// Calculate the Nernst potential for a given ion
    pub fn nernst_potential(&self, ion: IonType, temperature: units::Temperature) -> units::Voltage {
        use crate::constants::{GAS_CONSTANT, FARADAY_CONSTANT};
        
        let inside_conc = self.inside.get(&ion).copied().unwrap_or(0.0);
        let outside_conc = self.outside.get(&ion).copied().unwrap_or(0.0);
        
        if inside_conc <= 0.0 || outside_conc <= 0.0 {
            return 0.0;
        }
        
        let rt_over_zf = (GAS_CONSTANT * temperature) / (ion.charge() as f64 * FARADAY_CONSTANT);
        rt_over_zf * (outside_conc / inside_conc).ln()
    }
}

/// Lipid types in biological membranes
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum LipidType {
    /// Phosphatidylcholine
    POPC,
    /// Phosphatidylethanolamine
    POPE,
    /// Phosphatidylserine
    POPS,
    /// Sphingomyelin
    SM,
    /// Cholesterol
    Cholesterol,
    /// Phosphatidylinositol
    PI,
    /// Cardiolipin
    CL,
}

impl LipidType {
    /// Get the molecular weight in g/mol
    pub fn molecular_weight(&self) -> f64 {
        match self {
            Self::POPC => 760.076,
            Self::POPE => 717.992,
            Self::POPS => 761.967,
            Self::SM => 703.038,
            Self::Cholesterol => 386.654,
            Self::PI => 834.988,
            Self::CL => 1447.8,
        }
    }
    
    /// Get the headgroup charge at physiological pH
    pub fn charge(&self) -> i8 {
        match self {
            Self::POPC => 0,
            Self::POPE => 0,
            Self::POPS => -1,
            Self::SM => 0,
            Self::Cholesterol => 0,
            Self::PI => -1,
            Self::CL => -2,
        }
    }
}

/// Membrane lipid composition
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LipidComposition {
    pub fractions: HashMap<LipidType, f64>,
}

impl Default for LipidComposition {
    /// Default mammalian plasma membrane composition
    fn default() -> Self {
        let mut fractions = HashMap::new();
        fractions.insert(LipidType::POPC, 0.40);
        fractions.insert(LipidType::POPE, 0.25);
        fractions.insert(LipidType::POPS, 0.10);
        fractions.insert(LipidType::SM, 0.05);
        fractions.insert(LipidType::Cholesterol, 0.20);
        
        Self { fractions }
    }
}

impl LipidComposition {
    /// Validate that fractions sum to 1.0
    pub fn is_valid(&self) -> bool {
        let sum: f64 = self.fractions.values().sum();
        (sum - 1.0).abs() < 1e-6
    }
    
    /// Normalize fractions to sum to 1.0
    pub fn normalize(&mut self) {
        let sum: f64 = self.fractions.values().sum();
        if sum > 0.0 {
            for fraction in self.fractions.values_mut() {
                *fraction /= sum;
            }
        }
    }
}

/// Protein types in biological membranes
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ProteinType {
    /// Na⁺/K⁺-ATPase pump
    NaKATPase,
    /// Ca²⁺-ATPase pump
    CaATPase,
    /// Voltage-gated sodium channel
    VGSC,
    /// Voltage-gated potassium channel
    VGKC,
    /// Voltage-gated calcium channel
    VGCC,
    /// NMDA receptor
    NMDAR,
    /// AMPA receptor
    AMPAR,
    /// GABA receptor
    GABAR,
    /// Custom protein type
    Custom(String),
}

impl ProteinType {
    /// Check if this protein requires ATP
    pub fn requires_atp(&self) -> bool {
        matches!(self, Self::NaKATPase | Self::CaATPase)
    }
    
    /// Get typical ATP consumption rate (molecules/s per protein)
    pub fn atp_consumption_rate(&self) -> f64 {
        match self {
            Self::NaKATPase => 100.0,  // ~100 ATP/s
            Self::CaATPase => 50.0,    // ~50 ATP/s
            _ => 0.0,
        }
    }
    
    /// Get the protein's contribution to membrane conductance (S/protein)
    pub fn conductance(&self) -> units::Conductance {
        match self {
            Self::VGSC => 20e-12,      // 20 pS
            Self::VGKC => 10e-12,      // 10 pS
            Self::VGCC => 5e-12,       // 5 pS
            Self::NMDAR => 50e-12,     // 50 pS
            Self::AMPAR => 10e-12,     // 10 pS
            Self::GABAR => 30e-12,     // 30 pS
            _ => 0.0,
        }
    }
}

/// Current membrane state
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MembraneState {
    /// Membrane potential in volts
    pub voltage: units::Voltage,
    /// Current across membrane in amperes
    pub current: units::Current,
    /// Membrane capacitance in farads
    pub capacitance: units::Capacitance,
    /// Ion concentrations
    pub ion_concentrations: IonConcentrations,
    /// ATP measurements
    pub atp: AtpMeasurement,
    /// Temperature in Kelvin
    pub temperature: units::Temperature,
    /// Simulation time in seconds
    pub time: units::Time,
    /// Total energy consumed this step
    pub energy_consumed: units::Energy,
}

impl MembraneState {
    /// Create a new membrane state with physiological defaults
    pub fn physiological(area: units::Area, temperature: units::Temperature) -> Self {
        Self {
            voltage: -0.070,  // -70 mV resting potential
            current: 0.0,
            capacitance: area * 1e-2,  // 1 μF/cm² typical membrane capacitance
            ion_concentrations: IonConcentrations::physiological(),
            atp: AtpMeasurement {
                concentration: 5e-3,  // 5 mM ATP
                molecules: crate::constants::molecules_from_concentration(5e-3, area * 1e-6),
                consumption_rate: 0.0,
                production_rate: 0.0,
            },
            temperature,
            time: 0.0,
            energy_consumed: 0.0,
        }
    }
}

/// Circuit parameters derived from membrane state
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CircuitParameters {
    /// Equivalent circuit capacitance
    pub capacitance: units::Capacitance,
    /// Equivalent circuit resistance
    pub resistance: units::Resistance,
    /// Voltage sources (from ion pumps)
    pub voltage_sources: Vec<units::Voltage>,
    /// Current sources (from ion channels)
    pub current_sources: Vec<units::Current>,
    /// Circuit topology connections
    pub connections: Vec<(usize, usize)>,
    /// Update timestamp
    pub timestamp: units::Time,
}

/// Configuration for membrane simulation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MembraneConfig {
    /// Simulation timestep in seconds
    pub timestep: units::Time,
    /// Maximum simulation time
    pub max_time: units::Time,
    /// Temperature in Kelvin
    pub temperature: units::Temperature,
    /// Membrane area in m²
    pub area: units::Area,
    /// Initial ATP concentration
    pub initial_atp: units::Concentration,
    /// Enable ATP constraints
    pub atp_constrained: bool,
    /// Enable circuit interface
    pub circuit_interface: bool,
    /// Logging level
    pub log_level: String,
}

impl Default for MembraneConfig {
    fn default() -> Self {
        Self {
            timestep: 1e-6,           // 1 μs
            max_time: 1.0,            // 1 second
            temperature: 310.15,      // 37°C
            area: 1e-9,               // 1 μm²
            initial_atp: 5e-3,        // 5 mM
            atp_constrained: true,
            circuit_interface: true,
            log_level: "info".to_string(),
        }
    }
} 