//! Membrane Phase Separation Module
//! 
//! This module simulates liquid-liquid phase separation in membranes,
//! including protein clustering, lipid domain formation, and phase transitions.

use crate::types::*;
use crate::constants::*;
use crate::error::BeneGesseritError;
use std::collections::HashMap;

/// Types of phase separation
#[derive(Debug, Clone, PartialEq)]
pub enum PhaseType {
    LiquidLiquid,    // Classical LLPS
    LiquidGel,       // Liquid-gel transition
    Crystalline,     // Ordered phase formation
    Amorphous,       // Disordered clustering
}

/// Phase separation mechanism
#[derive(Debug, Clone, PartialEq)]
pub enum SeparationMechanism {
    Spinodal,        // Spinodal decomposition
    Nucleation,      // Nucleation and growth
    Protein,         // Protein-driven clustering
    Lipid,           // Lipid-driven phase separation
    Electrostatic,   // Charge-based separation
}

/// Phase domain state
#[derive(Debug, Clone)]
pub struct PhaseDomain {
    pub domain_id: u64,
    pub phase_type: PhaseType,
    pub center_position: (f64, f64),
    pub radius: f64,
    pub density: f64,
    pub composition: HashMap<ComponentType, f64>,
    pub surface_tension: f64,
    pub stability: f64,
    pub age: f64,
    pub is_growing: bool,
}

/// Component types in phase domains
#[derive(Debug, Clone, PartialEq, Hash)]
pub enum ComponentType {
    Protein(ProteinType),
    Lipid(LipidType),
    Nucleic(String),    // RNA, DNA
    Metabolite(String), // Small molecules
}

/// Phase transition event
#[derive(Debug, Clone)]
pub struct PhaseTransition {
    pub transition_id: u64,
    pub from_phase: PhaseType,
    pub to_phase: PhaseType,
    pub trigger: TransitionTrigger,
    pub critical_parameter: f64,
    pub energy_barrier: f64,
    pub time: f64,
    pub affected_domains: Vec<u64>,
}

/// Triggers for phase transitions
#[derive(Debug, Clone, PartialEq)]
pub enum TransitionTrigger {
    Temperature,
    Concentration,
    pH,
    Pressure,
    ElectricField,
    ChemicalModification,
}

/// Interaction parameters for phase behavior
#[derive(Debug, Clone)]
pub struct InteractionMatrix {
    /// Interaction strengths between component types
    pub interactions: HashMap<(ComponentType, ComponentType), f64>,
    /// Temperature dependence of interactions
    pub temperature_coefficients: HashMap<(ComponentType, ComponentType), f64>,
    /// Concentration dependence
    pub concentration_coefficients: HashMap<ComponentType, f64>,
}

/// Phase separation simulator
pub struct PhaseSeparation {
    /// Active phase domains
    pub domains: HashMap<u64, PhaseDomain>,
    /// Interaction parameters
    pub interaction_matrix: InteractionMatrix,
    /// Component concentrations
    pub concentrations: HashMap<ComponentType, f64>,
    /// Phase transition history
    pub transition_history: Vec<PhaseTransition>,
    /// Nucleation sites
    pub nucleation_sites: Vec<(f64, f64)>,
    /// Temperature
    pub temperature: f64,
    /// pH
    pub ph: f64,
    /// Ionic strength
    pub ionic_strength: f64,
    /// Membrane area
    pub membrane_area: f64,
    /// Critical parameters for phase transitions
    pub critical_parameters: HashMap<PhaseType, f64>,
    /// Domain counter
    domain_counter: u64,
}

impl Default for PhaseSeparation {
    fn default() -> Self {
        let mut interactions = HashMap::new();
        let mut temp_coeffs = HashMap::new();
        let mut conc_coeffs = HashMap::new();
        let mut critical_params = HashMap::new();
        
        // Initialize some basic interaction parameters
        // Protein-protein interactions (attractive)
        interactions.insert(
            (ComponentType::Protein(ProteinType::Membrane), ComponentType::Protein(ProteinType::Membrane)),
            -2.0 // kJ/mol
        );
        
        // Lipid-lipid interactions
        interactions.insert(
            (ComponentType::Lipid(LipidType::Cholesterol), ComponentType::Lipid(LipidType::Sphingomyelin)),
            -1.5 // Raft formation
        );
        
        // Temperature coefficients
        temp_coeffs.insert(
            (ComponentType::Protein(ProteinType::Membrane), ComponentType::Protein(ProteinType::Membrane)),
            0.01 // kJ/mol/K
        );
        
        // Concentration coefficients
        conc_coeffs.insert(ComponentType::Protein(ProteinType::Membrane), 0.5);
        
        // Critical parameters
        critical_params.insert(PhaseType::LiquidLiquid, 0.1); // Critical concentration
        critical_params.insert(PhaseType::LiquidGel, 0.3);
        critical_params.insert(PhaseType::Crystalline, 0.5);
        
        Self {
            domains: HashMap::new(),
            interaction_matrix: InteractionMatrix {
                interactions,
                temperature_coefficients: temp_coeffs,
                concentration_coefficients: conc_coeffs,
            },
            concentrations: HashMap::new(),
            transition_history: Vec::new(),
            nucleation_sites: Vec::new(),
            temperature: TEMPERATURE,
            ph: 7.4,
            ionic_strength: 0.15, // M
            membrane_area: 1e-12, // m^2 (1 μm^2)
            critical_parameters: critical_params,
            domain_counter: 0,
        }
    }
}

impl PhaseSeparation {
    pub fn new() -> Self {
        Self::default()
    }
    
    /// Set component concentration
    pub fn set_concentration(&mut self, component: ComponentType, concentration: f64) {
        self.concentrations.insert(component, concentration);
    }
    
    /// Add interaction parameter
    pub fn add_interaction(&mut self, comp1: ComponentType, comp2: ComponentType, strength: f64) {
        self.interaction_matrix.interactions.insert((comp1.clone(), comp2.clone()), strength);
        self.interaction_matrix.interactions.insert((comp2, comp1), strength); // Symmetric
    }
    
    /// Add nucleation site
    pub fn add_nucleation_site(&mut self, x: f64, y: f64) {
        self.nucleation_sites.push((x, y));
    }
    
    /// Check for phase separation conditions
    pub fn check_phase_separation(&self) -> Result<Vec<PhaseType>, BeneGesseritError> {
        let mut possible_phases = Vec::new();
        
        for (phase_type, &critical_param) in &self.critical_parameters {
            if self.is_phase_favorable(phase_type, critical_param)? {
                possible_phases.push(phase_type.clone());
            }
        }
        
        Ok(possible_phases)
    }
    
    /// Check if a phase is thermodynamically favorable
    fn is_phase_favorable(&self, phase_type: &PhaseType, critical_param: f64) -> Result<bool, BeneGesseritError> {
        match phase_type {
            PhaseType::LiquidLiquid => {
                // Check if any component exceeds critical concentration
                for (component, &concentration) in &self.concentrations {
                    let effective_critical = critical_param * self.get_interaction_factor(component)?;
                    if concentration > effective_critical {
                        return Ok(true);
                    }
                }
            }
            PhaseType::LiquidGel => {
                // Requires higher concentration and specific interactions
                let total_protein_conc = self.get_total_protein_concentration();
                if total_protein_conc > critical_param && self.temperature < TEMPERATURE * 0.9 {
                    return Ok(true);
                }
            }
            PhaseType::Crystalline => {
                // Requires very high concentration and low temperature
                let max_conc = self.concentrations.values().cloned().fold(0.0, f64::max);
                if max_conc > critical_param && self.temperature < TEMPERATURE * 0.8 {
                    return Ok(true);
                }
            }
            PhaseType::Amorphous => {
                // Always possible at some level
                return Ok(true);
            }
        }
        
        Ok(false)
    }
    
    /// Get interaction factor for a component
    fn get_interaction_factor(&self, component: &ComponentType) -> Result<f64, BeneGesseritError> {
        let mut factor = 1.0;
        
        // Temperature dependence
        if let Some(&temp_coeff) = self.interaction_matrix.temperature_coefficients.get(&(component.clone(), component.clone())) {
            factor *= 1.0 + temp_coeff * (self.temperature - TEMPERATURE);
        }
        
        // Concentration dependence
        if let Some(&conc_coeff) = self.interaction_matrix.concentration_coefficients.get(component) {
            if let Some(&concentration) = self.concentrations.get(component) {
                factor *= 1.0 + conc_coeff * concentration;
            }
        }
        
        Ok(factor)
    }
    
    /// Get total protein concentration
    fn get_total_protein_concentration(&self) -> f64 {
        self.concentrations.iter()
            .filter_map(|(component, &concentration)| {
                match component {
                    ComponentType::Protein(_) => Some(concentration),
                    _ => None,
                }
            })
            .sum()
    }
    
    /// Initiate domain nucleation
    pub fn nucleate_domain(&mut self, phase_type: PhaseType, position: (f64, f64)) -> Result<u64, BeneGesseritError> {
        let domain_id = self.domain_counter;
        self.domain_counter += 1;
        
        let mut composition = HashMap::new();
        
        // Determine initial composition based on phase type
        match phase_type {
            PhaseType::LiquidLiquid => {
                // Enrich in components with favorable interactions
                for (component, &concentration) in &self.concentrations {
                    let enrichment = self.calculate_enrichment_factor(component, &phase_type)?;
                    composition.insert(component.clone(), concentration * enrichment);
                }
            }
            PhaseType::LiquidGel => {
                // Protein-rich phase
                for (component, &concentration) in &self.concentrations {
                    let enrichment = match component {
                        ComponentType::Protein(_) => 3.0,
                        _ => 0.5,
                    };
                    composition.insert(component.clone(), concentration * enrichment);
                }
            }
            _ => {
                // Default composition
                composition = self.concentrations.clone();
            }
        }
        
        let domain = PhaseDomain {
            domain_id,
            phase_type,
            center_position: position,
            radius: 10.0, // nm - initial radius
            density: self.calculate_phase_density(&composition),
            composition,
            surface_tension: self.calculate_surface_tension(&phase_type),
            stability: 1.0,
            age: 0.0,
            is_growing: true,
        };
        
        self.domains.insert(domain_id, domain);
        Ok(domain_id)
    }
    
    /// Calculate enrichment factor for a component in a phase
    fn calculate_enrichment_factor(&self, component: &ComponentType, phase_type: &PhaseType) -> Result<f64, BeneGesseritError> {
        let mut enrichment = 1.0;
        
        // Check interactions with other components
        for other_component in self.concentrations.keys() {
            if let Some(&interaction_strength) = self.interaction_matrix.interactions.get(&(component.clone(), other_component.clone())) {
                if interaction_strength < 0.0 { // Attractive interaction
                    enrichment += (-interaction_strength) * 0.1;
                }
            }
        }
        
        // Phase-specific modifiers
        match phase_type {
            PhaseType::LiquidLiquid => enrichment *= 2.0,
            PhaseType::LiquidGel => {
                if matches!(component, ComponentType::Protein(_)) {
                    enrichment *= 3.0;
                }
            }
            _ => {}
        }
        
        Ok(enrichment)
    }
    
    /// Calculate phase density
    fn calculate_phase_density(&self, composition: &HashMap<ComponentType, f64>) -> f64 {
        let total_concentration: f64 = composition.values().sum();
        total_concentration * 1000.0 // Convert to density units
    }
    
    /// Calculate surface tension for a phase type
    fn calculate_surface_tension(&self, phase_type: &PhaseType) -> f64 {
        match phase_type {
            PhaseType::LiquidLiquid => 0.001, // N/m
            PhaseType::LiquidGel => 0.005,
            PhaseType::Crystalline => 0.01,
            PhaseType::Amorphous => 0.0001,
        }
    }
    
    /// Simulate phase separation dynamics
    pub fn simulate_phase_dynamics(&mut self, dt: f64) -> Result<Vec<PhaseTransition>, BeneGesseritError> {
        let mut transitions = Vec::new();
        
        // Update existing domains
        self.update_domain_growth(dt)?;
        self.update_domain_stability(dt)?;
        
        // Check for new nucleation events
        let favorable_phases = self.check_phase_separation()?;
        for phase_type in favorable_phases {
            if self.should_nucleate(&phase_type)? {
                let position = self.select_nucleation_site();
                self.nucleate_domain(phase_type, position)?;
            }
        }
        
        // Check for phase transitions
        transitions.extend(self.check_phase_transitions(dt)?);
        
        // Domain coalescence
        self.check_domain_coalescence(dt)?;
        
        Ok(transitions)
    }
    
    /// Update domain growth
    fn update_domain_growth(&mut self, dt: f64) -> Result<(), BeneGesseritError> {
        for domain in self.domains.values_mut() {
            if domain.is_growing {
                let growth_rate = self.calculate_growth_rate(domain)?;
                domain.radius += growth_rate * dt;
                domain.age += dt;
                
                // Update composition due to growth
                self.update_domain_composition(domain)?;
                
                // Check if growth should stop
                if domain.radius > 100.0 || domain.stability < 0.5 {
                    domain.is_growing = false;
                }
            }
        }
        
        Ok(())
    }
    
    /// Calculate domain growth rate
    fn calculate_growth_rate(&self, domain: &PhaseDomain) -> Result<f64, BeneGesseritError> {
        let base_rate = 1.0; // nm/s
        
        // Driving force from supersaturation
        let supersaturation = self.calculate_supersaturation(&domain.phase_type)?;
        let driving_force = supersaturation.max(0.0);
        
        // Surface tension effects
        let surface_effect = 1.0 / (1.0 + domain.surface_tension * domain.radius);
        
        Ok(base_rate * driving_force * surface_effect)
    }
    
    /// Calculate supersaturation for a phase type
    fn calculate_supersaturation(&self, phase_type: &PhaseType) -> Result<f64, BeneGesseritError> {
        let critical_conc = self.critical_parameters.get(phase_type).unwrap_or(&0.1);
        let current_conc = self.get_effective_concentration(phase_type);
        
        Ok((current_conc - critical_conc) / critical_conc)
    }
    
    /// Get effective concentration for phase type
    fn get_effective_concentration(&self, phase_type: &PhaseType) -> f64 {
        match phase_type {
            PhaseType::LiquidLiquid => {
                self.concentrations.values().cloned().fold(0.0, f64::max)
            }
            PhaseType::LiquidGel => self.get_total_protein_concentration(),
            _ => self.concentrations.values().sum(),
        }
    }
    
    /// Update domain composition during growth
    fn update_domain_composition(&self, domain: &mut PhaseDomain) -> Result<(), BeneGesseritError> {
        for (component, concentration) in &mut domain.composition {
            let enrichment = self.calculate_enrichment_factor(component, &domain.phase_type)?;
            let bulk_conc = self.concentrations.get(component).unwrap_or(&0.0);
            *concentration = bulk_conc * enrichment;
        }
        
        domain.density = self.calculate_phase_density(&domain.composition);
        Ok(())
    }
    
    /// Update domain stability
    fn update_domain_stability(&mut self, dt: f64) -> Result<(), BeneGesseritError> {
        let mut domains_to_remove = Vec::new();
        
        for (domain_id, domain) in self.domains.iter_mut() {
            // Stability decreases with age and unfavorable conditions
            let stability_change = -0.01 * dt; // Base decay
            
            // Environmental effects
            let temp_effect = if self.temperature > TEMPERATURE * 1.1 { -0.05 * dt } else { 0.0 };
            let ph_effect = if (self.ph - 7.4).abs() > 1.0 { -0.02 * dt } else { 0.0 };
            
            domain.stability += stability_change + temp_effect + ph_effect;
            
            // Remove unstable domains
            if domain.stability < 0.1 {
                domains_to_remove.push(*domain_id);
            }
        }
        
        for domain_id in domains_to_remove {
            self.domains.remove(&domain_id);
        }
        
        Ok(())
    }
    
    /// Check for phase transitions
    fn check_phase_transitions(&mut self, dt: f64) -> Result<Vec<PhaseTransition>, BeneGesseritError> {
        let mut transitions = Vec::new();
        
        for domain in self.domains.values_mut() {
            // Check temperature-induced transitions
            if let Some(transition) = self.check_temperature_transition(domain)? {
                transitions.push(transition);
            }
            
            // Check concentration-induced transitions
            if let Some(transition) = self.check_concentration_transition(domain)? {
                transitions.push(transition);
            }
        }
        
        Ok(transitions)
    }
    
    /// Check for temperature-induced phase transitions
    fn check_temperature_transition(&self, domain: &PhaseDomain) -> Result<Option<PhaseTransition>, BeneGesseritError> {
        let critical_temp = match domain.phase_type {
            PhaseType::LiquidLiquid => TEMPERATURE * 0.95,
            PhaseType::LiquidGel => TEMPERATURE * 0.85,
            PhaseType::Crystalline => TEMPERATURE * 0.75,
            PhaseType::Amorphous => TEMPERATURE * 1.2,
        };
        
        if (self.temperature - critical_temp).abs() < 1.0 {
            let new_phase = if self.temperature < critical_temp {
                PhaseType::LiquidGel
            } else {
                PhaseType::LiquidLiquid
            };
            
            if new_phase != domain.phase_type {
                return Ok(Some(PhaseTransition {
                    transition_id: self.domain_counter,
                    from_phase: domain.phase_type.clone(),
                    to_phase: new_phase,
                    trigger: TransitionTrigger::Temperature,
                    critical_parameter: critical_temp,
                    energy_barrier: 10.0, // kJ/mol
                    time: 0.0,
                    affected_domains: vec![domain.domain_id],
                }));
            }
        }
        
        Ok(None)
    }
    
    /// Check for concentration-induced phase transitions
    fn check_concentration_transition(&self, domain: &PhaseDomain) -> Result<Option<PhaseTransition>, BeneGesseritError> {
        let total_conc = domain.composition.values().sum::<f64>();
        
        if total_conc > 0.5 && matches!(domain.phase_type, PhaseType::LiquidLiquid) {
            return Ok(Some(PhaseTransition {
                transition_id: self.domain_counter,
                from_phase: domain.phase_type.clone(),
                to_phase: PhaseType::LiquidGel,
                trigger: TransitionTrigger::Concentration,
                critical_parameter: 0.5,
                energy_barrier: 15.0,
                time: 0.0,
                affected_domains: vec![domain.domain_id],
            }));
        }
        
        Ok(None)
    }
    
    /// Check for domain coalescence
    fn check_domain_coalescence(&mut self, dt: f64) -> Result<(), BeneGesseritError> {
        let domain_ids: Vec<u64> = self.domains.keys().cloned().collect();
        
        for i in 0..domain_ids.len() {
            for j in (i + 1)..domain_ids.len() {
                let id1 = domain_ids[i];
                let id2 = domain_ids[j];
                
                if let (Some(domain1), Some(domain2)) = (self.domains.get(&id1), self.domains.get(&id2)) {
                    let distance = ((domain1.center_position.0 - domain2.center_position.0).powi(2) +
                                   (domain1.center_position.1 - domain2.center_position.1).powi(2)).sqrt();
                    
                    if distance < (domain1.radius + domain2.radius) && domain1.phase_type == domain2.phase_type {
                        // Coalescence occurs
                        self.coalesce_domains(id1, id2)?;
                        break;
                    }
                }
            }
        }
        
        Ok(())
    }
    
    /// Coalesce two domains
    fn coalesce_domains(&mut self, id1: u64, id2: u64) -> Result<(), BeneGesseritError> {
        if let (Some(domain1), Some(domain2)) = (self.domains.remove(&id1), self.domains.remove(&id2)) {
            let total_volume = domain1.radius.powi(3) + domain2.radius.powi(3);
            let new_radius = total_volume.powf(1.0/3.0);
            
            // Weighted average position
            let v1 = domain1.radius.powi(3);
            let v2 = domain2.radius.powi(3);
            let new_position = (
                (domain1.center_position.0 * v1 + domain2.center_position.0 * v2) / (v1 + v2),
                (domain1.center_position.1 * v1 + domain2.center_position.1 * v2) / (v1 + v2)
            );
            
            // Merge compositions
            let mut new_composition = HashMap::new();
            for (component, &conc1) in &domain1.composition {
                let conc2 = domain2.composition.get(component).unwrap_or(&0.0);
                new_composition.insert(component.clone(), (conc1 * v1 + conc2 * v2) / (v1 + v2));
            }
            
            let merged_domain = PhaseDomain {
                domain_id: self.domain_counter,
                phase_type: domain1.phase_type,
                center_position: new_position,
                radius: new_radius,
                density: self.calculate_phase_density(&new_composition),
                composition: new_composition,
                surface_tension: (domain1.surface_tension + domain2.surface_tension) / 2.0,
                stability: (domain1.stability + domain2.stability) / 2.0,
                age: domain1.age.min(domain2.age),
                is_growing: domain1.is_growing || domain2.is_growing,
            };
            
            self.domains.insert(self.domain_counter, merged_domain);
            self.domain_counter += 1;
        }
        
        Ok(())
    }
    
    /// Determine if nucleation should occur
    fn should_nucleate(&self, phase_type: &PhaseType) -> Result<bool, BeneGesseritError> {
        let supersaturation = self.calculate_supersaturation(phase_type)?;
        let nucleation_probability = (supersaturation - 1.0).max(0.0) * 0.01;
        
        Ok(self.random_f64() < nucleation_probability)
    }
    
    /// Select nucleation site
    fn select_nucleation_site(&self) -> (f64, f64) {
        if !self.nucleation_sites.is_empty() {
            let index = (self.random_f64() * self.nucleation_sites.len() as f64) as usize;
            self.nucleation_sites[index.min(self.nucleation_sites.len() - 1)]
        } else {
            // Random position
            (self.random_f64() * 1000.0, self.random_f64() * 1000.0) // nm
        }
    }
    
    /// Get phase separation statistics
    pub fn get_phase_statistics(&self) -> PhaseStatistics {
        let mut phase_counts = HashMap::new();
        let mut total_volume = 0.0;
        let mut average_stability = 0.0;
        
        for domain in self.domains.values() {
            *phase_counts.entry(domain.phase_type.clone()).or_insert(0) += 1;
            total_volume += domain.radius.powi(3) * 4.0 / 3.0 * std::f64::consts::PI;
            average_stability += domain.stability;
        }
        
        if !self.domains.is_empty() {
            average_stability /= self.domains.len() as f64;
        }
        
        PhaseStatistics {
            total_domains: self.domains.len(),
            phase_counts,
            total_volume,
            volume_fraction: total_volume / self.membrane_area,
            average_stability,
            transition_count: self.transition_history.len(),
        }
    }
    
    /// Set temperature
    pub fn set_temperature(&mut self, temperature: f64) {
        self.temperature = temperature;
    }
    
    /// Set pH
    pub fn set_ph(&mut self, ph: f64) {
        self.ph = ph;
    }
    
    /// Set ionic strength
    pub fn set_ionic_strength(&mut self, ionic_strength: f64) {
        self.ionic_strength = ionic_strength;
    }
    
    // Helper methods
    fn random_f64(&self) -> f64 {
        // Placeholder - would use proper RNG
        0.5
    }
}

/// Phase separation statistics
#[derive(Debug, Clone)]
pub struct PhaseStatistics {
    pub total_domains: usize,
    pub phase_counts: HashMap<PhaseType, usize>,
    pub total_volume: f64,
    pub volume_fraction: f64,
    pub average_stability: f64,
    pub transition_count: usize,
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_phase_separation_creation() {
        let phase_sep = PhaseSeparation::new();
        assert_eq!(phase_sep.domains.len(), 0);
        assert!(phase_sep.critical_parameters.len() > 0);
    }
    
    #[test]
    fn test_concentration_setting() {
        let mut phase_sep = PhaseSeparation::new();
        phase_sep.set_concentration(ComponentType::Protein(ProteinType::Membrane), 0.5);
        assert_eq!(phase_sep.concentrations.len(), 1);
    }
    
    #[test]
    fn test_domain_nucleation() {
        let mut phase_sep = PhaseSeparation::new();
        phase_sep.set_concentration(ComponentType::Protein(ProteinType::Membrane), 0.5);
        
        let result = phase_sep.nucleate_domain(PhaseType::LiquidLiquid, (100.0, 100.0));
        assert!(result.is_ok());
        assert_eq!(phase_sep.domains.len(), 1);
    }
}       