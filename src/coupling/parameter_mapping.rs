//! Parameter mapping between membrane biophysics and circuit elements
//!
//! This module provides the translation layer between biological membrane
//! parameters and electrical circuit equivalents, including:
//! - Biophysical to electrical parameter conversion
//! - Dynamic parameter updates based on membrane state
//! - Multi-scale parameter coupling
//! - Calibration and validation frameworks

use crate::{
    constants::*,
    error::{MembraneError, Result},
    types::*,
};
use std::collections::HashMap;

/// Parameter mapping configuration
#[derive(Debug, Clone)]
pub struct ParameterMapping {
    /// Mapping functions from biophysical to electrical parameters
    pub mappings: HashMap<String, MappingFunction>,
    /// Calibration parameters
    pub calibration: CalibrationParameters,
    /// Validation metrics
    pub validation: ValidationMetrics,
    /// Update frequency
    pub update_frequency: f64,
}

/// Mapping function definition
#[derive(Debug, Clone)]
pub struct MappingFunction {
    /// Function type
    pub function_type: MappingType,
    /// Input parameters (biophysical)
    pub inputs: Vec<String>,
    /// Output parameter (electrical)
    pub output: String,
    /// Function coefficients
    pub coefficients: Vec<f64>,
    /// Scaling factors
    pub scaling: ScalingParameters,
    /// Valid input ranges
    pub input_ranges: Vec<(f64, f64)>,
}

/// Types of mapping functions
#[derive(Debug, Clone, PartialEq)]
pub enum MappingType {
    /// Linear mapping: output = a * input + b
    Linear,
    /// Exponential mapping: output = a * exp(b * input)
    Exponential,
    /// Power law mapping: output = a * input^b
    PowerLaw,
    /// Sigmoidal mapping: output = a / (1 + exp(-b * (input - c)))
    Sigmoidal,
    /// Polynomial mapping: output = sum(a_i * input^i)
    Polynomial,
    /// Custom lookup table
    LookupTable(Vec<(f64, f64)>),
}

/// Scaling parameters for mapping functions
#[derive(Debug, Clone)]
pub struct ScalingParameters {
    /// Input scaling factor
    pub input_scale: f64,
    /// Output scaling factor
    pub output_scale: f64,
    /// Input offset
    pub input_offset: f64,
    /// Output offset
    pub output_offset: f64,
}

/// Calibration parameters for mapping validation
#[derive(Debug, Clone)]
pub struct CalibrationParameters {
    /// Reference experimental data
    pub reference_data: HashMap<String, Vec<(f64, f64)>>,
    /// Tolerance for validation
    pub tolerance: f64,
    /// R-squared threshold for goodness of fit
    pub r_squared_threshold: f64,
    /// Maximum allowed error
    pub max_error: f64,
}

/// Validation metrics for mapping quality
#[derive(Debug, Clone)]
pub struct ValidationMetrics {
    /// Root mean square error
    pub rmse: HashMap<String, f64>,
    /// R-squared values
    pub r_squared: HashMap<String, f64>,
    /// Maximum absolute error
    pub max_abs_error: HashMap<String, f64>,
    /// Mean absolute percentage error
    pub mape: HashMap<String, f64>,
}

/// Parameter mapper that handles all conversions
#[derive(Debug, Clone)]
pub struct ParameterMapper {
    /// Parameter mappings
    pub mapping: ParameterMapping,
    /// Current biophysical state
    pub biophysical_state: BiophysicalState,
    /// Current electrical parameters
    pub electrical_params: ElectricalParameters,
    /// Mapping history for trend analysis
    pub mapping_history: MappingHistory,
}

/// Current biophysical state
#[derive(Debug, Clone)]
pub struct BiophysicalState {
    /// Membrane properties
    pub membrane_area: f64,
    pub membrane_thickness: f64,
    pub lipid_composition: LipidComposition,
    pub protein_density: HashMap<ProteinType, f64>,
    
    /// Ion concentrations
    pub ion_concentrations: IonConcentrations,
    
    /// ATP state
    pub atp_concentration: f64,
    pub atp_consumption_rate: f64,
    
    /// Temperature
    pub temperature: f64,
    
    /// Protein activities
    pub protein_activities: HashMap<String, f64>,
}

/// Electrical circuit parameters
#[derive(Debug, Clone)]
pub struct ElectricalParameters {
    /// Membrane capacitance (F)
    pub capacitance: f64,
    /// Membrane resistance (Ω)
    pub resistance: f64,
    /// Voltage sources (V)
    pub voltage_sources: Vec<f64>,
    /// Current sources (A)
    pub current_sources: Vec<f64>,
    /// Conductances (S)
    pub conductances: HashMap<String, f64>,
}

/// History of parameter mappings
#[derive(Debug, Clone)]
pub struct MappingHistory {
    /// Time series of parameter values
    pub parameter_history: HashMap<String, Vec<(f64, f64)>>,
    /// Mapping quality over time
    pub mapping_quality: Vec<f64>,
    /// Update timestamps
    pub timestamps: Vec<f64>,
}

impl ParameterMapper {
    /// Create a new parameter mapper
    pub fn new() -> Self {
        Self {
            mapping: ParameterMapping::new(),
            biophysical_state: BiophysicalState::new(),
            electrical_params: ElectricalParameters::new(),
            mapping_history: MappingHistory::new(),
        }
    }
    
    /// Update biophysical state from membrane
    pub fn update_biophysical_state(&mut self, membrane_state: &MembraneState,
                                  lipid_composition: &LipidComposition,
                                  proteins: &HashMap<String, crate::molecular::MembraneProtein>) -> Result<()> {
        // Update basic membrane properties
        self.biophysical_state.temperature = membrane_state.temperature;
        self.biophysical_state.atp_concentration = membrane_state.atp.concentration;
        self.biophysical_state.ion_concentrations = membrane_state.ion_concentrations.clone();
        self.biophysical_state.lipid_composition = lipid_composition.clone();
        
        // Update protein activities
        self.biophysical_state.protein_activities.clear();
        for (protein_id, protein) in proteins {
            self.biophysical_state.protein_activities.insert(
                protein_id.clone(), 
                protein.activity
            );
        }
        
        // Calculate protein densities
        self.update_protein_densities(proteins)?;
        
        Ok(())
    }
    
    /// Update protein densities from protein map
    fn update_protein_densities(&mut self, proteins: &HashMap<String, crate::molecular::MembraneProtein>) -> Result<()> {
        self.biophysical_state.protein_density.clear();
        
        // Count proteins by type
        let mut protein_counts: HashMap<ProteinType, u32> = HashMap::new();
        for protein in proteins.values() {
            *protein_counts.entry(protein.protein_type).or_insert(0) += 1;
        }
        
        // Convert counts to densities (proteins per μm²)
        for (protein_type, count) in protein_counts {
            let density = count as f64 / (self.biophysical_state.membrane_area * 1e12); // Convert m² to μm²
            self.biophysical_state.protein_density.insert(protein_type, density);
        }
        
        Ok(())
    }
    
    /// Map all biophysical parameters to electrical parameters
    pub fn map_parameters(&mut self) -> Result<()> {
        // Map membrane capacitance
        self.electrical_params.capacitance = self.map_capacitance()?;
        
        // Map membrane resistance
        self.electrical_params.resistance = self.map_resistance()?;
        
        // Map voltage sources (from ion pumps)
        self.electrical_params.voltage_sources = self.map_voltage_sources()?;
        
        // Map current sources (from ion channels)
        self.electrical_params.current_sources = self.map_current_sources()?;
        
        // Map individual conductances
        self.electrical_params.conductances = self.map_conductances()?;
        
        Ok(())
    }
    
    /// Map membrane capacitance from biophysical properties
    fn map_capacitance(&self) -> Result<f64> {
        // Base membrane capacitance
        let base_capacitance = self.biophysical_state.membrane_area * SPECIFIC_CAPACITANCE;
        
        // Adjust for lipid composition
        let lipid_factor = self.calculate_lipid_capacitance_factor()?;
        
        // Adjust for temperature
        let temp_factor = self.calculate_temperature_factor()?;
        
        Ok(base_capacitance * lipid_factor * temp_factor)
    }
    
    /// Calculate lipid composition effect on capacitance
    fn calculate_lipid_capacitance_factor(&self) -> Result<f64> {
        let mut factor = 1.0;
        
        for (lipid_type, fraction) in &self.biophysical_state.lipid_composition.fractions {
            let lipid_factor = match lipid_type {
                LipidType::Cholesterol => 0.8, // Cholesterol reduces capacitance
                LipidType::POPC => 1.0,        // Reference lipid
                LipidType::POPE => 1.1,        // Slightly higher capacitance
                LipidType::POPS => 1.2,        // Charged lipid increases capacitance
                _ => 1.0,
            };
            
            factor += (lipid_factor - 1.0) * fraction;
        }
        
        Ok(factor)
    }
    
    /// Calculate temperature effect on electrical properties
    fn calculate_temperature_factor(&self) -> Result<f64> {
        // Q10 temperature dependence
        let temp_ref = PHYSIOLOGICAL_TEMPERATURE;
        let q10 = temperature::Q10_BIOLOGICAL;
        
        Ok(q10.powf((self.biophysical_state.temperature - temp_ref) / 10.0))
    }
    
    /// Map membrane resistance from biophysical properties
    fn map_resistance(&self) -> Result<f64> {
        // Base membrane resistance
        let base_resistance = membrane::MEMBRANE_RESISTANCE / self.biophysical_state.membrane_area;
        
        // Adjust for protein channels
        let channel_conductance = self.calculate_total_channel_conductance()?;
        let protein_resistance = if channel_conductance > 0.0 {
            1.0 / channel_conductance
        } else {
            f64::INFINITY
        };
        
        // Parallel resistance combination
        let total_resistance = 1.0 / (1.0 / base_resistance + 1.0 / protein_resistance);
        
        Ok(total_resistance)
    }
    
    /// Calculate total conductance from membrane proteins
    fn calculate_total_channel_conductance(&self) -> Result<f64> {
        let mut total_conductance = 0.0;
        
        for (protein_type, density) in &self.biophysical_state.protein_density {
            let single_channel_conductance = protein_type.conductance();
            let protein_conductance = single_channel_conductance * density * self.biophysical_state.membrane_area * 1e-12; // Convert μm² to m²
            total_conductance += protein_conductance;
        }
        
        Ok(total_conductance)
    }
    
    /// Map voltage sources from ion pumps
    fn map_voltage_sources(&self) -> Result<Vec<f64>> {
        let mut voltage_sources = Vec::new();
        
        // Na/K-ATPase voltage
        if let Some(nak_density) = self.biophysical_state.protein_density.get(&ProteinType::NaKATPase) {
            if *nak_density > 0.0 {
                let nak_voltage = self.calculate_nak_voltage()?;
                voltage_sources.push(nak_voltage);
            }
        }
        
        // Ca-ATPase voltage
        if let Some(ca_density) = self.biophysical_state.protein_density.get(&ProteinType::CaATPase) {
            if *ca_density > 0.0 {
                let ca_voltage = self.calculate_ca_voltage()?;
                voltage_sources.push(ca_voltage);
            }
        }
        
        Ok(voltage_sources)
    }
    
    /// Calculate Na/K-ATPase equivalent voltage
    fn calculate_nak_voltage(&self) -> Result<f64> {
        // Nernst potentials for Na+ and K+
        let na_nernst = self.biophysical_state.ion_concentrations.nernst_potential(
            IonType::Sodium, self.biophysical_state.temperature
        );
        let k_nernst = self.biophysical_state.ion_concentrations.nernst_potential(
            IonType::Potassium, self.biophysical_state.temperature
        );
        
        // Stoichiometry: 3 Na+ out, 2 K+ in
        let pump_voltage = 3.0 * na_nernst - 2.0 * k_nernst;
        
        // Scale by ATP availability
        let atp_factor = self.biophysical_state.atp_concentration / PHYSIOLOGICAL_ATP;
        
        Ok(pump_voltage * atp_factor)
    }
    
    /// Calculate Ca-ATPase equivalent voltage
    fn calculate_ca_voltage(&self) -> Result<f64> {
        let ca_nernst = self.biophysical_state.ion_concentrations.nernst_potential(
            IonType::Calcium, self.biophysical_state.temperature
        );
        
        // Scale by ATP availability
        let atp_factor = self.biophysical_state.atp_concentration / PHYSIOLOGICAL_ATP;
        
        Ok(ca_nernst * atp_factor)
    }
    
    /// Map current sources from ion channels
    fn map_current_sources(&self) -> Result<Vec<f64>> {
        let mut current_sources = Vec::new();
        
        // Voltage-gated channels
        for (protein_type, density) in &self.biophysical_state.protein_density {
            match protein_type {
                ProteinType::VGSC | ProteinType::VGKC | ProteinType::VGCC => {
                    let channel_current = self.calculate_channel_current(protein_type, *density)?;
                    if channel_current.abs() > 1e-15 { // Only add significant currents
                        current_sources.push(channel_current);
                    }
                },
                _ => {}, // Only channels contribute to current sources
            }
        }
        
        Ok(current_sources)
    }
    
    /// Calculate current from a specific channel type
    fn calculate_channel_current(&self, protein_type: &ProteinType, density: f64) -> Result<f64> {
        let single_channel_conductance = protein_type.conductance();
        let total_conductance = single_channel_conductance * density * self.biophysical_state.membrane_area * 1e-12;
        
        // Get driving force for the channel
        let driving_force = match protein_type {
            ProteinType::VGSC => {
                self.biophysical_state.ion_concentrations.nernst_potential(IonType::Sodium, self.biophysical_state.temperature)
                // Should subtract membrane voltage, but we don't have it in biophysical state
            },
            ProteinType::VGKC => {
                self.biophysical_state.ion_concentrations.nernst_potential(IonType::Potassium, self.biophysical_state.temperature)
            },
            ProteinType::VGCC => {
                self.biophysical_state.ion_concentrations.nernst_potential(IonType::Calcium, self.biophysical_state.temperature)
            },
            _ => 0.0,
        };
        
        // Get channel activity
        let activity = self.biophysical_state.protein_activities
            .get(&format!("{:?}", protein_type))
            .copied()
            .unwrap_or(0.1); // Default low activity
        
        Ok(total_conductance * driving_force * activity)
    }
    
    /// Map individual conductances
    fn map_conductances(&self) -> Result<HashMap<String, f64>> {
        let mut conductances = HashMap::new();
        
        for (protein_type, density) in &self.biophysical_state.protein_density {
            let conductance = protein_type.conductance() * density * self.biophysical_state.membrane_area * 1e-12;
            conductances.insert(format!("{:?}", protein_type), conductance);
        }
        
        Ok(conductances)
    }
    
    /// Get current electrical parameters
    pub fn get_electrical_parameters(&self) -> &ElectricalParameters {
        &self.electrical_params
    }
    
    /// Convert electrical parameters to circuit parameters
    pub fn to_circuit_parameters(&self) -> CircuitParameters {
        CircuitParameters {
            capacitance: self.electrical_params.capacitance,
            resistance: self.electrical_params.resistance,
            voltage_sources: self.electrical_params.voltage_sources.clone(),
            current_sources: self.electrical_params.current_sources.clone(),
            connections: Vec::new(), // Would need topology information
            timestamp: 0.0, // Would need current time
        }
    }
    
    /// Validate mapping against reference data
    pub fn validate_mapping(&mut self) -> Result<()> {
        // For each mapping, compare with reference data if available
        for (param_name, reference_data) in &self.mapping.calibration.reference_data {
            if let Some(predicted_values) = self.get_predicted_values(param_name) {
                let rmse = self.calculate_rmse(&predicted_values, reference_data);
                let r_squared = self.calculate_r_squared(&predicted_values, reference_data);
                
                self.mapping.validation.rmse.insert(param_name.clone(), rmse);
                self.mapping.validation.r_squared.insert(param_name.clone(), r_squared);
            }
        }
        
        Ok(())
    }
    
    /// Get predicted values for a parameter
    fn get_predicted_values(&self, param_name: &str) -> Option<Vec<f64>> {
        // This would depend on the specific parameter and available data
        // For now, return None (placeholder implementation)
        None
    }
    
    /// Calculate root mean square error
    fn calculate_rmse(&self, predicted: &[f64], actual: &[(f64, f64)]) -> f64 {
        if predicted.len() != actual.len() {
            return f64::INFINITY;
        }
        
        let mse: f64 = predicted.iter()
            .zip(actual.iter())
            .map(|(pred, (_, actual_val))| (pred - actual_val).powi(2))
            .sum::<f64>() / predicted.len() as f64;
        
        mse.sqrt()
    }
    
    /// Calculate R-squared
    fn calculate_r_squared(&self, predicted: &[f64], actual: &[(f64, f64)]) -> f64 {
        if predicted.len() != actual.len() || actual.is_empty() {
            return 0.0;
        }
        
        let actual_values: Vec<f64> = actual.iter().map(|(_, val)| *val).collect();
        let mean_actual: f64 = actual_values.iter().sum::<f64>() / actual_values.len() as f64;
        
        let ss_tot: f64 = actual_values.iter()
            .map(|val| (val - mean_actual).powi(2))
            .sum();
        
        let ss_res: f64 = predicted.iter()
            .zip(actual_values.iter())
            .map(|(pred, actual_val)| (actual_val - pred).powi(2))
            .sum();
        
        if ss_tot == 0.0 {
            1.0
        } else {
            1.0 - (ss_res / ss_tot)
        }
    }
}

impl ParameterMapping {
    pub fn new() -> Self {
        Self {
            mappings: HashMap::new(),
            calibration: CalibrationParameters::new(),
            validation: ValidationMetrics::new(),
            update_frequency: 1000.0, // 1 kHz
        }
    }
}

impl BiophysicalState {
    pub fn new() -> Self {
        Self {
            membrane_area: simulation::DEFAULT_PATCH_AREA,
            membrane_thickness: MEMBRANE_THICKNESS,
            lipid_composition: LipidComposition::default(),
            protein_density: HashMap::new(),
            ion_concentrations: IonConcentrations::physiological(),
            atp_concentration: PHYSIOLOGICAL_ATP,
            atp_consumption_rate: 0.0,
            temperature: PHYSIOLOGICAL_TEMPERATURE,
            protein_activities: HashMap::new(),
        }
    }
}

impl ElectricalParameters {
    pub fn new() -> Self {
        Self {
            capacitance: simulation::DEFAULT_PATCH_AREA * SPECIFIC_CAPACITANCE,
            resistance: membrane::MEMBRANE_RESISTANCE / simulation::DEFAULT_PATCH_AREA,
            voltage_sources: Vec::new(),
            current_sources: Vec::new(),
            conductances: HashMap::new(),
        }
    }
}

impl MappingHistory {
    pub fn new() -> Self {
        Self {
            parameter_history: HashMap::new(),
            mapping_quality: Vec::new(),
            timestamps: Vec::new(),
        }
    }
}

impl CalibrationParameters {
    pub fn new() -> Self {
        Self {
            reference_data: HashMap::new(),
            tolerance: 0.1, // 10% tolerance
            r_squared_threshold: 0.8,
            max_error: 0.2, // 20% max error
        }
    }
}

impl ValidationMetrics {
    pub fn new() -> Self {
        Self {
            rmse: HashMap::new(),
            r_squared: HashMap::new(),
            max_abs_error: HashMap::new(),
            mape: HashMap::new(),
        }
    }
}

impl Default for ScalingParameters {
    fn default() -> Self {
        Self {
            input_scale: 1.0,
            output_scale: 1.0,
            input_offset: 0.0,
            output_offset: 0.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_parameter_mapper_creation() {
        let mapper = ParameterMapper::new();
        assert_eq!(mapper.biophysical_state.temperature, PHYSIOLOGICAL_TEMPERATURE);
        assert!(mapper.electrical_params.capacitance > 0.0);
    }
    
    #[test]
    fn test_capacitance_mapping() {
        let mapper = ParameterMapper::new();
        let capacitance = mapper.map_capacitance().unwrap();
        assert!(capacitance > 0.0);
        
        // Should be close to theoretical value
        let expected = simulation::DEFAULT_PATCH_AREA * SPECIFIC_CAPACITANCE;
        assert!((capacitance / expected - 1.0).abs() < 0.5); // Within 50%
    }
    
    #[test]
    fn test_temperature_factor() {
        let mapper = ParameterMapper::new();
        let factor = mapper.calculate_temperature_factor().unwrap();
        assert!((factor - 1.0).abs() < 0.1); // Should be close to 1.0 at physiological temp
    }
}
