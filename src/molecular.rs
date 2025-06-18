//! Molecular Layer - Lipid bilayer physics and protein interactions
//!
//! This module implements the lowest level of membrane dynamics, handling:
//! - Lipid bilayer mechanics and thermodynamics
//! - Membrane protein interactions and conformational changes
//! - Electrochemical gradients and ion transport
//! - ATP-driven processes with dx/dATP differential equations
//! - Membrane curvature and elastic properties

use crate::{
    constants::*,
    error::{MembraneError, Result},
    types::*,
};
use std::collections::HashMap;

/// Molecular-level membrane patch simulation
#[derive(Debug, Clone)]
pub struct MolecularMembrane {
    /// Current membrane state
    pub state: MembraneState,
    /// Lipid composition and properties
    pub lipid_composition: LipidComposition,
    /// Membrane proteins and their states
    pub proteins: HashMap<String, MembraneProtein>,
    /// Membrane geometry and curvature
    pub geometry: MembraneGeometry,
    /// ATP metabolism tracker
    pub atp_tracker: AtpTracker,
    /// Simulation configuration
    pub config: MembraneConfig,
}

/// Individual membrane protein with conformational states
#[derive(Debug, Clone)]
pub struct MembraneProtein {
    /// Protein type and properties
    pub protein_type: ProteinType,
    /// Current conformational state
    pub conformation: ProteinConformation,
    /// Position in membrane (x, y coordinates in nm)
    pub position: (f64, f64),
    /// Orientation angle (radians)
    pub orientation: f64,
    /// ATP binding sites and occupancy
    pub atp_sites: AtpBindingSites,
    /// Activity state (0.0 = inactive, 1.0 = fully active)
    pub activity: f64,
    /// Energy consumption history
    pub energy_history: Vec<f64>,
}

/// Protein conformational states
#[derive(Debug, Clone, PartialEq)]
pub enum ProteinConformation {
    /// Closed/inactive state
    Closed,
    /// Open/active state
    Open,
    /// Intermediate states for multi-step processes
    Intermediate(u8),
    /// ATP-bound states
    AtpBound,
    /// ADP-bound states
    AdpBound,
}

/// ATP binding sites on proteins
#[derive(Debug, Clone)]
pub struct AtpBindingSites {
    /// Number of ATP binding sites
    pub site_count: u8,
    /// Current ATP occupancy (0.0 to 1.0 per site)
    pub occupancy: Vec<f64>,
    /// Binding affinity constants (M⁻¹)
    pub binding_constants: Vec<f64>,
    /// ATP consumption rates per site (s⁻¹)
    pub consumption_rates: Vec<f64>,
}

/// Membrane geometry and curvature properties
#[derive(Debug, Clone)]
pub struct MembraneGeometry {
    /// Membrane area (m²)
    pub area: f64,
    /// Mean curvature (m⁻¹)
    pub mean_curvature: f64,
    /// Gaussian curvature (m⁻²)
    pub gaussian_curvature: f64,
    /// Local thickness variations (m)
    pub thickness_map: Vec<f64>,
    /// Elastic moduli
    pub elastic_properties: ElasticProperties,
}

/// Membrane elastic properties
#[derive(Debug, Clone)]
pub struct ElasticProperties {
    /// Bending modulus (J)
    pub bending_modulus: f64,
    /// Gaussian modulus (J)
    pub gaussian_modulus: f64,
    /// Area compressibility modulus (N/m)
    pub area_modulus: f64,
    /// Surface tension (N/m)
    pub surface_tension: f64,
}

/// ATP metabolism and energy tracking
#[derive(Debug, Clone)]
pub struct AtpTracker {
    /// Current ATP pool
    pub atp_pool: AtpMeasurement,
    /// ATP consumption by process type
    pub consumption_by_process: HashMap<String, f64>,
    /// Energy flux measurements
    pub energy_flux: EnergyFlux,
    /// Metabolic efficiency metrics
    pub efficiency: MetabolicEfficiency,
}

/// Energy flux measurements
#[derive(Debug, Clone)]
pub struct EnergyFlux {
    /// Total power consumption (W)
    pub total_power: f64,
    /// Power by protein type (W)
    pub power_by_protein: HashMap<ProteinType, f64>,
    /// Mechanical work done by membrane (J/s)
    pub mechanical_work: f64,
    /// Heat generation rate (W)
    pub heat_generation: f64,
}

/// Metabolic efficiency tracking
#[derive(Debug, Clone)]
pub struct MetabolicEfficiency {
    /// ATP utilization efficiency (0.0 to 1.0)
    pub atp_efficiency: f64,
    /// Energy coupling efficiency
    pub coupling_efficiency: f64,
    /// Waste heat fraction
    pub waste_heat_fraction: f64,
}

impl MolecularMembrane {
    /// Create a new molecular membrane with physiological defaults
    pub fn new(config: MembraneConfig) -> Result<Self> {
        let state = MembraneState::physiological(config.area, config.temperature);
        
        let geometry = MembraneGeometry {
            area: config.area,
            mean_curvature: 0.0,
            gaussian_curvature: 0.0,
            thickness_map: vec![MEMBRANE_THICKNESS; 100], // 10x10 grid
            elastic_properties: ElasticProperties::default(),
        };
        
        let atp_tracker = AtpTracker::new(config.initial_atp, config.area);
        
        Ok(Self {
            state,
            lipid_composition: LipidComposition::default(),
            proteins: HashMap::new(),
            geometry,
            atp_tracker,
            config,
        })
    }
    
    /// Add a membrane protein at specified position
    pub fn add_protein(&mut self, id: String, protein_type: ProteinType, position: (f64, f64)) -> Result<()> {
        let protein = MembraneProtein::new(protein_type, position)?;
        self.proteins.insert(id, protein);
        Ok(())
    }
    
    /// Simulate one timestep using ATP-based differential equations
    pub fn step(&mut self, dt: f64) -> Result<()> {
        // Validate timestep
        if dt <= 0.0 || dt > simulation::MAX_STABLE_TIMESTEP {
            return Err(MembraneError::InvalidParameter(
                format!("Invalid timestep: {}", dt)
            ));
        }
        
        // Update ATP availability
        self.update_atp_dynamics(dt)?;
        
        // Update protein conformations and activities
        self.update_protein_dynamics(dt)?;
        
        // Calculate ion transport and electrochemical gradients
        self.update_ion_transport(dt)?;
        
        // Update membrane geometry and mechanics
        self.update_membrane_mechanics(dt)?;
        
        // Calculate energy consumption and efficiency
        self.update_energy_accounting(dt)?;
        
        // Update simulation time
        self.state.time += dt;
        
        Ok(())
    }
    
    /// Update ATP dynamics using dx/dATP equations
    fn update_atp_dynamics(&mut self, dt: f64) -> Result<()> {
        let current_atp = self.state.atp.concentration;
        
        // Calculate ATP consumption by all proteins
        let mut total_consumption = 0.0;
        for protein in self.proteins.values() {
            if protein.protein_type.requires_atp() {
                let consumption_rate = protein.calculate_atp_consumption_rate(current_atp);
                total_consumption += consumption_rate * protein.activity;
            }
        }
        
        // ATP-constrained dynamics: dx/dATP instead of dx/dt
        let atp_availability_factor = if self.config.atp_constrained {
            // Michaelis-Menten kinetics for ATP availability
            current_atp / (current_atp + 1e-3) // Km ≈ 1 mM
        } else {
            1.0
        };
        
        // Update ATP concentration
        let datp_dt = -total_consumption * atp_availability_factor;
        let new_atp = (current_atp + datp_dt * dt).max(0.0);
        
        self.state.atp.concentration = new_atp;
        self.state.atp.consumption_rate = total_consumption;
        
        // Update molecule count
        self.state.atp.molecules = molecules_from_concentration(
            new_atp,
            self.config.area * MEMBRANE_THICKNESS
        );
        
        Ok(())
    }
    
    /// Update protein conformational dynamics
    fn update_protein_dynamics(&mut self, dt: f64) -> Result<()> {
        let atp_conc = self.state.atp.concentration;
        let temperature = self.state.temperature;
        let voltage = self.state.voltage;
        
        for protein in self.proteins.values_mut() {
            // Update ATP binding
            protein.update_atp_binding(atp_conc, dt)?;
            
            // Update conformational state
            protein.update_conformation(voltage, temperature, dt)?;
            
            // Update activity based on ATP availability and conformation
            protein.update_activity(atp_conc, dt)?;
        }
        
        Ok(())
    }
    
    /// Update ion transport and membrane potential
    fn update_ion_transport(&mut self, dt: f64) -> Result<()> {
        let mut total_current = 0.0;
        let temperature = self.state.temperature;
        
        // Calculate current from each protein
        for protein in self.proteins.values() {
            let protein_current = protein.calculate_current(&self.state.ion_concentrations, temperature);
            total_current += protein_current * protein.activity;
        }
        
        // Update membrane potential using capacitor equation
        let capacitance = self.state.capacitance;
        let dv_dt = total_current / capacitance;
        self.state.voltage += dv_dt * dt;
        self.state.current = total_current;
        
        // Update ion concentrations based on transport
        self.update_ion_concentrations(total_current, dt)?;
        
        Ok(())
    }
    
    /// Update ion concentrations from transport currents
    fn update_ion_concentrations(&mut self, current: f64, dt: f64) -> Result<()> {
        // Simplified ion concentration updates
        // In practice, this would track each ion species separately
        let volume_inside = self.config.area * MEMBRANE_THICKNESS;
        let volume_outside = volume_inside * 10.0; // Assume large extracellular volume
        
        // Update based on current and protein stoichiometry
        for protein in self.proteins.values() {
            match protein.protein_type {
                ProteinType::NaKATPase => {
                    // 3 Na+ out, 2 K+ in per ATP
                    let rate = protein.activity * protein.protein_type.atp_consumption_rate();
                    let na_flux = rate * 3.0 / AVOGADRO_NUMBER;
                    let k_flux = rate * 2.0 / AVOGADRO_NUMBER;
                    
                    // Update Na+ concentrations
                    if let Some(na_in) = self.state.ion_concentrations.inside.get_mut(&IonType::Sodium) {
                        *na_in -= na_flux * dt / volume_inside;
                    }
                    if let Some(na_out) = self.state.ion_concentrations.outside.get_mut(&IonType::Sodium) {
                        *na_out += na_flux * dt / volume_outside;
                    }
                    
                    // Update K+ concentrations
                    if let Some(k_in) = self.state.ion_concentrations.inside.get_mut(&IonType::Potassium) {
                        *k_in += k_flux * dt / volume_inside;
                    }
                    if let Some(k_out) = self.state.ion_concentrations.outside.get_mut(&IonType::Potassium) {
                        *k_out -= k_flux * dt / volume_outside;
                    }
                }
                _ => {} // Other protein types handled similarly
            }
        }
        
        Ok(())
    }
    
    /// Update membrane mechanics and geometry
    fn update_membrane_mechanics(&mut self, dt: f64) -> Result<()> {
        // Calculate forces from protein activity and ion gradients
        let osmotic_pressure = self.calculate_osmotic_pressure();
        let electrostatic_pressure = self.calculate_electrostatic_pressure();
        
        // Update membrane curvature based on protein distribution and forces
        self.update_membrane_curvature(osmotic_pressure, electrostatic_pressure, dt)?;
        
        // Update elastic properties based on lipid composition and temperature
        self.update_elastic_properties()?;
        
        Ok(())
    }
    
    /// Calculate osmotic pressure across membrane
    fn calculate_osmotic_pressure(&self) -> f64 {
        let mut osmolarity_inside = 0.0;
        let mut osmolarity_outside = 0.0;
        
        for (_, conc) in &self.state.ion_concentrations.inside {
            osmolarity_inside += conc;
        }
        for (_, conc) in &self.state.ion_concentrations.outside {
            osmolarity_outside += conc;
        }
        
        // van 't Hoff equation: π = RTΔc
        GAS_CONSTANT * self.state.temperature * (osmolarity_inside - osmolarity_outside)
    }
    
    /// Calculate electrostatic pressure from membrane potential
    fn calculate_electrostatic_pressure(&self) -> f64 {
        let voltage = self.state.voltage;
        let epsilon = physical::EPSILON_0 * physical::EPSILON_MEMBRANE;
        
        // Maxwell stress: P = ε₀εᵣE²/2
        0.5 * epsilon * (voltage / MEMBRANE_THICKNESS).powi(2)
    }
    
    /// Update membrane curvature from mechanical forces
    fn update_membrane_curvature(&mut self, osmotic_p: f64, electrostatic_p: f64, dt: f64) -> Result<()> {
        let total_pressure = osmotic_p + electrostatic_p;
        let bending_modulus = self.geometry.elastic_properties.bending_modulus;
        
        // Simplified curvature evolution
        let curvature_change = total_pressure / bending_modulus * dt;
        self.geometry.mean_curvature += curvature_change;
        
        // Limit curvature to physically reasonable values
        self.geometry.mean_curvature = self.geometry.mean_curvature.clamp(-1e6, 1e6);
        
        Ok(())
    }
    
    /// Update elastic properties based on composition and temperature
    fn update_elastic_properties(&mut self) -> Result<()> {
        let temp = self.state.temperature;
        let cholesterol_fraction = self.lipid_composition.fractions
            .get(&LipidType::Cholesterol)
            .copied()
            .unwrap_or(0.0);
        
        // Temperature and cholesterol effects on bending modulus
        let temp_factor = q10_factor(temperature::PHYSIOLOGICAL, temp, 1.5);
        let cholesterol_factor = 1.0 + cholesterol_fraction * 2.0;
        
        self.geometry.elastic_properties.bending_modulus = 20.0 * physical::BOLTZMANN_CONSTANT * temp * temp_factor * cholesterol_factor;
        
        Ok(())
    }
    
    /// Update energy accounting and efficiency metrics
    fn update_energy_accounting(&mut self, dt: f64) -> Result<()> {
        let mut total_power = 0.0;
        let mut power_by_protein = HashMap::new();
        
        // Calculate power consumption by each protein type
        for protein in self.proteins.values() {
            let protein_power = protein.calculate_power_consumption();
            total_power += protein_power;
            
            *power_by_protein.entry(protein.protein_type.clone()).or_insert(0.0) += protein_power;
        }
        
        // Update energy flux
        self.atp_tracker.energy_flux.total_power = total_power;
        self.atp_tracker.energy_flux.power_by_protein = power_by_protein;
        
        // Calculate mechanical work from membrane deformation
        let mechanical_work = self.calculate_mechanical_work(dt);
        self.atp_tracker.energy_flux.mechanical_work = mechanical_work;
        
        // Heat generation (energy not converted to work)
        self.atp_tracker.energy_flux.heat_generation = total_power - mechanical_work;
        
        // Update total energy consumed
        self.state.energy_consumed += total_power * dt;
        
        // Calculate efficiency metrics
        self.calculate_efficiency_metrics()?;
        
        Ok(())
    }
    
    /// Calculate mechanical work done by membrane
    fn calculate_mechanical_work(&self, dt: f64) -> f64 {
        // Work from membrane curvature changes
        let bending_work = self.geometry.elastic_properties.bending_modulus * 
            self.geometry.mean_curvature.powi(2) * self.geometry.area;
        
        // Work from area changes (osmotic work)
        let osmotic_work = self.calculate_osmotic_pressure() * self.geometry.area;
        
        (bending_work + osmotic_work) / dt
    }
    
    /// Calculate metabolic efficiency metrics
    fn calculate_efficiency_metrics(&mut self) -> Result<()> {
        let total_atp_consumption = self.state.atp.consumption_rate;
        let mechanical_work = self.atp_tracker.energy_flux.mechanical_work;
        let total_power = self.atp_tracker.energy_flux.total_power;
        
        // ATP utilization efficiency
        self.atp_tracker.efficiency.atp_efficiency = if total_atp_consumption > 0.0 {
            (total_atp_consumption * ATP_HYDROLYSIS_ENERGY.abs()) / (total_power + 1e-12)
        } else {
            0.0
        };
        
        // Energy coupling efficiency (work output / energy input)
        self.atp_tracker.efficiency.coupling_efficiency = if total_power > 0.0 {
            mechanical_work / total_power
        } else {
            0.0
        };
        
        // Waste heat fraction
        self.atp_tracker.efficiency.waste_heat_fraction = if total_power > 0.0 {
            self.atp_tracker.energy_flux.heat_generation / total_power
        } else {
            0.0
        };
        
        Ok(())
    }
    
    /// Get current circuit parameters for interface layer
    pub fn get_circuit_parameters(&self) -> CircuitParameters {
        let mut voltage_sources = Vec::new();
        let mut current_sources = Vec::new();
        let mut connections = Vec::new();
        
        // Add voltage sources from ion pumps
        for (i, protein) in self.proteins.values().enumerate() {
            match protein.protein_type {
                ProteinType::NaKATPase | ProteinType::CaATPase => {
                    voltage_sources.push(protein.calculate_pump_voltage(&self.state.ion_concentrations));
                }
                _ => {
                    current_sources.push(protein.calculate_current(&self.state.ion_concentrations, self.state.temperature));
                }
            }
            
            // Simple linear connection topology for now
            if i > 0 {
                connections.push((i - 1, i));
            }
        }
        
        CircuitParameters {
            capacitance: self.state.capacitance,
            resistance: 1.0 / self.calculate_total_conductance(),
            voltage_sources,
            current_sources,
            connections,
            timestamp: self.state.time,
        }
    }
    
    /// Calculate total membrane conductance
    fn calculate_total_conductance(&self) -> f64 {
        self.proteins.values()
            .map(|p| p.protein_type.conductance() * p.activity)
            .sum::<f64>()
            .max(1e-12) // Avoid division by zero
    }
}

impl MembraneProtein {
    /// Create a new membrane protein
    pub fn new(protein_type: ProteinType, position: (f64, f64)) -> Result<Self> {
        let atp_sites = AtpBindingSites::new(&protein_type);
        
        Ok(Self {
            protein_type,
            conformation: ProteinConformation::Closed,
            position,
            orientation: 0.0,
            atp_sites,
            activity: 0.0,
            energy_history: Vec::new(),
        })
    }
    
    /// Calculate ATP consumption rate based on current state
    pub fn calculate_atp_consumption_rate(&self, atp_concentration: f64) -> f64 {
        if !self.protein_type.requires_atp() {
            return 0.0;
        }
        
        // Michaelis-Menten kinetics
        let vmax = self.protein_type.atp_consumption_rate();
        let km = 1e-3; // 1 mM Km for ATP
        
        vmax * atp_concentration / (km + atp_concentration)
    }
    
    /// Update ATP binding to protein sites
    pub fn update_atp_binding(&mut self, atp_concentration: f64, dt: f64) -> Result<()> {
        for (i, occupancy) in self.atp_sites.occupancy.iter_mut().enumerate() {
            let binding_constant = self.atp_sites.binding_constants[i];
            let consumption_rate = self.atp_sites.consumption_rates[i];
            
            // Binding rate
            let binding_rate = binding_constant * atp_concentration * (1.0 - *occupancy);
            
            // Consumption rate (ATP hydrolysis)
            let hydrolysis_rate = consumption_rate * *occupancy;
            
            // Net change in occupancy
            let doccupancy_dt = binding_rate - hydrolysis_rate;
            *occupancy = (*occupancy + doccupancy_dt * dt).clamp(0.0, 1.0);
        }
        
        Ok(())
    }
    
    /// Update protein conformational state
    pub fn update_conformation(&mut self, voltage: f64, temperature: f64, dt: f64) -> Result<()> {
        let thermal_energy = physical::BOLTZMANN_CONSTANT * temperature;
        
        match self.protein_type {
            ProteinType::VGSC | ProteinType::VGKC | ProteinType::VGCC => {
                // Voltage-gated channel kinetics
                let activation_energy = transport::CHANNEL_ACTIVATION_ENERGY;
                let voltage_factor = (-voltage * FARADAY_CONSTANT / thermal_energy).exp();
                
                let open_probability = 1.0 / (1.0 + voltage_factor * (activation_energy / thermal_energy).exp());
                
                // Stochastic state transitions
                let transition_rate = 1000.0; // 1/ms typical
                if open_probability > 0.5 {
                    self.conformation = ProteinConformation::Open;
                } else {
                    self.conformation = ProteinConformation::Closed;
                }
            }
            ProteinType::NaKATPase | ProteinType::CaATPase => {
                // ATP-driven pump cycle
                let avg_occupancy: f64 = self.atp_sites.occupancy.iter().sum::<f64>() / self.atp_sites.site_count as f64;
                
                if avg_occupancy > 0.8 {
                    self.conformation = ProteinConformation::AtpBound;
                } else if avg_occupancy > 0.2 {
                    self.conformation = ProteinConformation::Intermediate(1);
                } else {
                    self.conformation = ProteinConformation::Closed;
                }
            }
            _ => {
                // Default behavior for other proteins
                self.conformation = ProteinConformation::Closed;
            }
        }
        
        Ok(())
    }
    
    /// Update protein activity based on conformation and ATP
    pub fn update_activity(&mut self, atp_concentration: f64, dt: f64) -> Result<()> {
        let base_activity = match self.conformation {
            ProteinConformation::Open => 1.0,
            ProteinConformation::AtpBound => 1.0,
            ProteinConformation::Intermediate(_) => 0.5,
            ProteinConformation::Closed => 0.0,
            ProteinConformation::AdpBound => 0.1,
        };
        
        // ATP availability factor
        let atp_factor = if self.protein_type.requires_atp() {
            atp_concentration / (atp_concentration + 1e-3)
        } else {
            1.0
        };
        
        // Smooth activity transitions
        let target_activity = base_activity * atp_factor;
        let tau = 1e-3; // 1 ms time constant
        let activity_change = (target_activity - self.activity) / tau * dt;
        
        self.activity = (self.activity + activity_change).clamp(0.0, 1.0);
        
        Ok(())
    }
    
    /// Calculate current generated by this protein
    pub fn calculate_current(&self, ion_concentrations: &IonConcentrations, temperature: f64) -> f64 {
        let conductance = self.protein_type.conductance();
        
        match self.protein_type {
            ProteinType::VGSC => {
                let nernst_v = ion_concentrations.nernst_potential(IonType::Sodium, temperature);
                conductance * nernst_v * self.activity
            }
            ProteinType::VGKC => {
                let nernst_v = ion_concentrations.nernst_potential(IonType::Potassium, temperature);
                conductance * nernst_v * self.activity
            }
            ProteinType::VGCC => {
                let nernst_v = ion_concentrations.nernst_potential(IonType::Calcium, temperature);
                conductance * nernst_v * self.activity
            }
            _ => 0.0, // Pumps don't contribute to passive current
        }
    }
    
    /// Calculate pump voltage for ATP-driven transporters
    pub fn calculate_pump_voltage(&self, ion_concentrations: &IonConcentrations) -> f64 {
        match self.protein_type {
            ProteinType::NaKATPase => {
                // Electrogenic pump: 3 Na+ out, 2 K+ in
                let na_nernst = ion_concentrations.nernst_potential(IonType::Sodium, PHYSIOLOGICAL_TEMPERATURE);
                let k_nernst = ion_concentrations.nernst_potential(IonType::Potassium, PHYSIOLOGICAL_TEMPERATURE);
                (3.0 * na_nernst - 2.0 * k_nernst) * self.activity
            }
            ProteinType::CaATPase => {
                let ca_nernst = ion_concentrations.nernst_potential(IonType::Calcium, PHYSIOLOGICAL_TEMPERATURE);
                ca_nernst * self.activity
            }
            _ => 0.0,
        }
    }
    
    /// Calculate power consumption of this protein
    pub fn calculate_power_consumption(&self) -> f64 {
        if !self.protein_type.requires_atp() {
            return 0.0;
        }
        
        let atp_rate = self.protein_type.atp_consumption_rate();
        let energy_per_atp = atp::ATP_ENERGY_PER_MOLECULE.abs();
        
        atp_rate * energy_per_atp * self.activity
    }
}

impl AtpBindingSites {
    /// Create ATP binding sites for a protein type
    pub fn new(protein_type: &ProteinType) -> Self {
        match protein_type {
            ProteinType::NaKATPase => Self {
                site_count: 1,
                occupancy: vec![0.0],
                binding_constants: vec![1e6], // 1/M
                consumption_rates: vec![100.0], // 1/s
            },
            ProteinType::CaATPase => Self {
                site_count: 1,
                occupancy: vec![0.0],
                binding_constants: vec![5e5], // 1/M
                consumption_rates: vec![50.0], // 1/s
            },
            _ => Self {
                site_count: 0,
                occupancy: vec![],
                binding_constants: vec![],
                consumption_rates: vec![],
            },
        }
    }
}

impl AtpTracker {
    /// Create a new ATP tracker
    pub fn new(initial_concentration: f64, area: f64) -> Self {
        let volume = area * MEMBRANE_THICKNESS;
        
        Self {
            atp_pool: AtpMeasurement {
                concentration: initial_concentration,
                molecules: molecules_from_concentration(initial_concentration, volume),
                consumption_rate: 0.0,
                production_rate: 0.0,
            },
            consumption_by_process: HashMap::new(),
            energy_flux: EnergyFlux {
                total_power: 0.0,
                power_by_protein: HashMap::new(),
                mechanical_work: 0.0,
                heat_generation: 0.0,
            },
            efficiency: MetabolicEfficiency {
                atp_efficiency: 0.0,
                coupling_efficiency: 0.0,
                waste_heat_fraction: 0.0,
            },
        }
    }
}

impl Default for ElasticProperties {
    fn default() -> Self {
        let kt = physical::BOLTZMANN_CONSTANT * PHYSIOLOGICAL_TEMPERATURE;
        
        Self {
            bending_modulus: 20.0 * kt,  // ~20 kT for typical membranes
            gaussian_modulus: -10.0 * kt, // Negative for saddle-splay
            area_modulus: 0.24,          // 240 mN/m typical
            surface_tension: 0.0,        // Relaxed membrane
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_molecular_membrane_creation() {
        let config = MembraneConfig::default();
        let membrane = MolecularMembrane::new(config).unwrap();
        
        assert_eq!(membrane.proteins.len(), 0);
        assert!(membrane.state.voltage < 0.0); // Negative resting potential
        assert!(membrane.state.atp.concentration > 0.0);
    }
    
    #[test]
    fn test_protein_addition() {
        let config = MembraneConfig::default();
        let mut membrane = MolecularMembrane::new(config).unwrap();
        
        membrane.add_protein("pump1".to_string(), ProteinType::NaKATPase, (0.0, 0.0)).unwrap();
        assert_eq!(membrane.proteins.len(), 1);
    }
    
    #[test]
    fn test_atp_consumption() {
        let protein = MembraneProtein::new(ProteinType::NaKATPase, (0.0, 0.0)).unwrap();
        let consumption = protein.calculate_atp_consumption_rate(5e-3); // 5 mM ATP
        assert!(consumption > 0.0);
        assert!(consumption <= 100.0); // Max rate
    }
    
    #[test]
    fn test_simulation_step() {
        let config = MembraneConfig::default();
        let mut membrane = MolecularMembrane::new(config).unwrap();
        
        membrane.add_protein("pump1".to_string(), ProteinType::NaKATPase, (0.0, 0.0)).unwrap();
        
        let initial_time = membrane.state.time;
        membrane.step(1e-6).unwrap(); // 1 μs step
        
        assert!(membrane.state.time > initial_time);
    }
} 