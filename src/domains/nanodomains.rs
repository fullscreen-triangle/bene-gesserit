//! Membrane Nanodomains Module
//! 
//! This module simulates the formation and dynamics of membrane nanodomains,
//! small-scale membrane organizations typically < 100 nm in size.

use crate::types::*;
use crate::constants::*;
use crate::error::BeneGesseritError;
use std::collections::HashMap;

/// Types of nanodomains
#[derive(Debug, Clone, PartialEq)]
pub enum NanodomainType {
    ProteinCluster,     // Protein-driven clustering
    LipidPatch,         // Lipid-driven organization
    Signaling,          // Signaling-specific domains
    Receptor,           // Receptor clustering domains
    Enzymatic,          // Enzyme activity domains
    Transient,          // Short-lived domains
}

/// Nanodomain state
#[derive(Debug, Clone, PartialEq)]
pub enum NanodomainState {
    Forming,            // Initial assembly
    Stable,             // Maintained structure
    Fluctuating,        // Dynamic assembly/disassembly
    Dissolving,         // Breaking apart
}

/// Individual nanodomain
#[derive(Debug, Clone)]
pub struct Nanodomain {
    pub domain_id: u64,
    pub domain_type: NanodomainType,
    pub state: NanodomainState,
    pub center_position: (f64, f64),
    pub radius: f64,
    pub protein_composition: HashMap<ProteinType, usize>,
    pub lipid_composition: HashMap<LipidType, f64>,
    pub stability: f64,
    pub lifetime: f64,
    pub formation_energy: f64,
    pub binding_affinity: f64,
    pub activity_level: f64,
    pub is_active: bool,
}

/// Nanodomain interaction
#[derive(Debug, Clone)]
pub struct NanodomainInteraction {
    pub interaction_id: u64,
    pub domain1_id: u64,
    pub domain2_id: u64,
    pub interaction_type: InteractionType,
    pub strength: f64,
    pub distance: f64,
    pub duration: f64,
    pub is_active: bool,
}

/// Types of nanodomain interactions
#[derive(Debug, Clone, PartialEq)]
pub enum InteractionType {
    Attractive,         // Domains attract each other
    Repulsive,          // Domains repel each other
    Cooperative,        // Mutual stabilization
    Competitive,        // Competition for resources
    Catalytic,          // One domain catalyzes the other
}

/// Nanodomain formation parameters
#[derive(Debug, Clone)]
pub struct NanodomainParams {
    pub min_size: f64,              // Minimum domain size (nm)
    pub max_size: f64,              // Maximum domain size (nm)
    pub formation_threshold: f64,    // Concentration threshold
    pub stability_threshold: f64,    // Minimum stability
    pub interaction_range: f64,      // Range for domain interactions (nm)
    pub fluctuation_rate: f64,       // Rate of size fluctuations
    pub dissolution_rate: f64,       // Rate of dissolution
}

/// Nanodomain system
pub struct Nanodomains {
    /// Active nanodomains
    pub domains: HashMap<u64, Nanodomain>,
    /// Domain interactions
    pub interactions: HashMap<u64, NanodomainInteraction>,
    /// Formation parameters
    pub params: NanodomainParams,
    /// Component concentrations
    pub concentrations: HashMap<ComponentId, f64>,
    /// Temperature
    pub temperature: f64,
    /// Membrane area
    pub membrane_area: f64,
    /// Formation events history
    pub formation_events: Vec<NanodomainEvent>,
    /// Domain counter
    domain_counter: u64,
    /// Interaction counter
    interaction_counter: u64,
}

/// Component identifier for nanodomains
#[derive(Debug, Clone, PartialEq, Hash)]
pub enum ComponentId {
    Protein(ProteinType),
    Lipid(LipidType),
    Complex(String),
}

/// Nanodomain formation event
#[derive(Debug, Clone)]
pub struct NanodomainEvent {
    pub event_id: u64,
    pub event_type: NanodomainEventType,
    pub domain_id: u64,
    pub time: f64,
    pub trigger: NanodomainTrigger,
    pub components_involved: Vec<ComponentId>,
}

/// Types of nanodomain events
#[derive(Debug, Clone, PartialEq)]
pub enum NanodomainEventType {
    Formation,
    Growth,
    Shrinkage,
    Dissolution,
    Interaction,
    Activation,
    Deactivation,
}

/// Triggers for nanodomain events
#[derive(Debug, Clone, PartialEq)]
pub enum NanodomainTrigger {
    ConcentrationChange,
    ProteinBinding,
    LipidModification,
    TemperatureChange,
    ActivityChange,
    ExternalSignal,
}

impl Default for NanodomainParams {
    fn default() -> Self {
        Self {
            min_size: 5.0,              // 5 nm minimum
            max_size: 100.0,            // 100 nm maximum
            formation_threshold: 0.05,   // 5% concentration
            stability_threshold: 0.3,    // 30% stability
            interaction_range: 50.0,     // 50 nm interaction range
            fluctuation_rate: 0.1,       // s^-1
            dissolution_rate: 0.01,      // s^-1
        }
    }
}

impl Default for Nanodomains {
    fn default() -> Self {
        Self {
            domains: HashMap::new(),
            interactions: HashMap::new(),
            params: NanodomainParams::default(),
            concentrations: HashMap::new(),
            temperature: TEMPERATURE,
            membrane_area: 1e-12, // m^2 (1 μm^2)
            formation_events: Vec::new(),
            domain_counter: 0,
            interaction_counter: 0,
        }
    }
}

impl Nanodomains {
    pub fn new() -> Self {
        Self::default()
    }
    
    /// Set component concentration
    pub fn set_concentration(&mut self, component: ComponentId, concentration: f64) {
        self.concentrations.insert(component, concentration);
    }
    
    /// Check if conditions favor nanodomain formation
    pub fn check_formation_conditions(&self, domain_type: &NanodomainType) -> bool {
        match domain_type {
            NanodomainType::ProteinCluster => {
                // Check protein concentrations
                let protein_conc: f64 = self.concentrations.iter()
                    .filter_map(|(id, &conc)| {
                        match id {
                            ComponentId::Protein(_) => Some(conc),
                            _ => None,
                        }
                    })
                    .sum();
                protein_conc > self.params.formation_threshold
            }
            NanodomainType::LipidPatch => {
                // Check lipid concentrations
                let lipid_conc: f64 = self.concentrations.iter()
                    .filter_map(|(id, &conc)| {
                        match id {
                            ComponentId::Lipid(_) => Some(conc),
                            _ => None,
                        }
                    })
                    .sum();
                lipid_conc > self.params.formation_threshold * 2.0
            }
            NanodomainType::Signaling => {
                // Check for signaling proteins
                self.concentrations.iter().any(|(id, &conc)| {
                    matches!(id, ComponentId::Protein(ProteinType::Kinase) | 
                                ComponentId::Protein(ProteinType::Phosphatase)) &&
                    conc > self.params.formation_threshold
                })
            }
            _ => true, // Other types can form under general conditions
        }
    }
    
    /// Nucleate new nanodomain
    pub fn nucleate_domain(
        &mut self,
        domain_type: NanodomainType,
        position: (f64, f64),
        initial_components: HashMap<ComponentId, f64>,
    ) -> Result<u64, BeneGesseritError> {
        let domain_id = self.domain_counter;
        self.domain_counter += 1;
        
        // Calculate initial properties
        let initial_radius = self.calculate_initial_size(&domain_type, &initial_components);
        let formation_energy = self.calculate_formation_energy(&domain_type, &initial_components);
        let binding_affinity = self.calculate_binding_affinity(&domain_type);
        
        // Separate proteins and lipids
        let mut protein_composition = HashMap::new();
        let mut lipid_composition = HashMap::new();
        
        for (component, amount) in initial_components {
            match component {
                ComponentId::Protein(protein_type) => {
                    protein_composition.insert(protein_type, amount as usize);
                }
                ComponentId::Lipid(lipid_type) => {
                    lipid_composition.insert(lipid_type, amount);
                }
                _ => {} // Handle complexes if needed
            }
        }
        
        let domain = Nanodomain {
            domain_id,
            domain_type: domain_type.clone(),
            state: NanodomainState::Forming,
            center_position: position,
            radius: initial_radius,
            protein_composition,
            lipid_composition,
            stability: 1.0,
            lifetime: 0.0,
            formation_energy,
            binding_affinity,
            activity_level: self.calculate_initial_activity(&domain_type),
            is_active: true,
        };
        
        self.domains.insert(domain_id, domain);
        
        // Record formation event
        let event = NanodomainEvent {
            event_id: self.domain_counter,
            event_type: NanodomainEventType::Formation,
            domain_id,
            time: 0.0,
            trigger: NanodomainTrigger::ConcentrationChange,
            components_involved: Vec::new(),
        };
        self.formation_events.push(event);
        
        Ok(domain_id)
    }
    
    /// Calculate initial domain size
    fn calculate_initial_size(&self, domain_type: &NanodomainType, components: &HashMap<ComponentId, f64>) -> f64 {
        let base_size = match domain_type {
            NanodomainType::ProteinCluster => 15.0,  // nm
            NanodomainType::LipidPatch => 25.0,
            NanodomainType::Signaling => 20.0,
            NanodomainType::Receptor => 18.0,
            NanodomainType::Enzymatic => 12.0,
            NanodomainType::Transient => 8.0,
        };
        
        // Scale by component amount
        let total_components: f64 = components.values().sum();
        let size_factor = (1.0 + total_components * 0.1).min(3.0);
        
        (base_size * size_factor).clamp(self.params.min_size, self.params.max_size)
    }
    
    /// Calculate formation energy
    fn calculate_formation_energy(&self, domain_type: &NanodomainType, components: &HashMap<ComponentId, f64>) -> f64 {
        let base_energy = match domain_type {
            NanodomainType::ProteinCluster => -5.0,  // kJ/mol (favorable)
            NanodomainType::LipidPatch => -3.0,
            NanodomainType::Signaling => -8.0,      // Highly favorable
            NanodomainType::Receptor => -6.0,
            NanodomainType::Enzymatic => -4.0,
            NanodomainType::Transient => -1.0,      // Weakly favorable
        };
        
        // Cooperative effects
        let component_count = components.len() as f64;
        let cooperative_bonus = -component_count * 0.5;
        
        base_energy + cooperative_bonus
    }
    
    /// Calculate binding affinity
    fn calculate_binding_affinity(&self, domain_type: &NanodomainType) -> f64 {
        match domain_type {
            NanodomainType::ProteinCluster => 2.0,
            NanodomainType::LipidPatch => 1.5,
            NanodomainType::Signaling => 3.0,
            NanodomainType::Receptor => 2.5,
            NanodomainType::Enzymatic => 2.0,
            NanodomainType::Transient => 0.5,
        }
    }
    
    /// Calculate initial activity level
    fn calculate_initial_activity(&self, domain_type: &NanodomainType) -> f64 {
        match domain_type {
            NanodomainType::Signaling => 2.0,      // Enhanced activity
            NanodomainType::Enzymatic => 1.8,
            NanodomainType::Receptor => 1.5,
            _ => 1.0,                               // Normal activity
        }
    }
    
    /// Simulate nanodomain dynamics
    pub fn simulate_nanodomain_dynamics(&mut self, dt: f64) -> Result<Vec<NanodomainEvent>, BeneGesseritError> {
        let mut events = Vec::new();
        
        // Update existing domains
        self.update_domain_states(dt)?;
        
        // Check for new domain formation
        for domain_type in [
            NanodomainType::ProteinCluster,
            NanodomainType::LipidPatch,
            NanodomainType::Signaling,
            NanodomainType::Receptor,
            NanodomainType::Enzymatic,
            NanodomainType::Transient,
        ] {
            if self.check_formation_conditions(&domain_type) && self.should_nucleate_domain(&domain_type)? {
                let position = self.select_nucleation_site();
                let components = self.select_initial_components(&domain_type);
                let domain_id = self.nucleate_domain(domain_type, position, components)?;
                events.push(self.formation_events.last().unwrap().clone());
            }
        }
        
        // Update domain interactions
        self.update_domain_interactions(dt)?;
        
        // Check for domain dissolution
        events.extend(self.check_domain_dissolution(dt)?);
        
        Ok(events)
    }
    
    /// Update domain states
    fn update_domain_states(&mut self, dt: f64) -> Result<(), BeneGesseritError> {
        for domain in self.domains.values_mut() {
            domain.lifetime += dt;
            
            match domain.state {
                NanodomainState::Forming => {
                    // Transition to stable or fluctuating
                    if domain.lifetime > 0.1 { // 100 ms formation time
                        domain.state = if domain.stability > 0.7 {
                            NanodomainState::Stable
                        } else {
                            NanodomainState::Fluctuating
                        };
                    }
                }
                NanodomainState::Stable => {
                    // Slowly decrease stability
                    domain.stability -= 0.001 * dt;
                    
                    if domain.stability < self.params.stability_threshold {
                        domain.state = NanodomainState::Fluctuating;
                    }
                }
                NanodomainState::Fluctuating => {
                    // Random size fluctuations
                    let fluctuation = (self.random_f64() - 0.5) * self.params.fluctuation_rate * dt;
                    domain.radius += fluctuation;
                    domain.radius = domain.radius.clamp(self.params.min_size, self.params.max_size);
                    
                    // Check for stabilization or dissolution
                    if domain.stability > 0.8 {
                        domain.state = NanodomainState::Stable;
                    } else if domain.stability < 0.2 {
                        domain.state = NanodomainState::Dissolving;
                    }
                }
                NanodomainState::Dissolving => {
                    // Shrink and decrease stability
                    domain.radius -= 2.0 * dt; // nm/s
                    domain.stability -= 0.1 * dt;
                    
                    if domain.radius < self.params.min_size || domain.stability < 0.1 {
                        domain.is_active = false;
                    }
                }
            }
            
            // Update activity based on state and composition
            domain.activity_level = self.calculate_current_activity(domain);
        }
        
        // Remove inactive domains
        self.domains.retain(|_, domain| domain.is_active);
        
        Ok(())
    }
    
    /// Calculate current activity level
    fn calculate_current_activity(&self, domain: &Nanodomain) -> f64 {
        let base_activity = match domain.domain_type {
            NanodomainType::Signaling => 2.0,
            NanodomainType::Enzymatic => 1.8,
            NanodomainType::Receptor => 1.5,
            _ => 1.0,
        };
        
        // Modulate by stability and size
        let stability_factor = domain.stability;
        let size_factor = (domain.radius / 20.0).min(2.0); // Optimal around 20 nm
        
        base_activity * stability_factor * size_factor
    }
    
    /// Update domain interactions
    fn update_domain_interactions(&mut self, dt: f64) -> Result<(), BeneGesseritError> {
        let domain_ids: Vec<u64> = self.domains.keys().cloned().collect();
        
        // Check for new interactions
        for i in 0..domain_ids.len() {
            for j in (i + 1)..domain_ids.len() {
                let id1 = domain_ids[i];
                let id2 = domain_ids[j];
                
                if let (Some(domain1), Some(domain2)) = (self.domains.get(&id1), self.domains.get(&id2)) {
                    let distance = self.calculate_distance(domain1, domain2);
                    
                    if distance < self.params.interaction_range {
                        // Check if interaction already exists
                        let existing_interaction = self.interactions.values()
                            .any(|interaction| {
                                (interaction.domain1_id == id1 && interaction.domain2_id == id2) ||
                                (interaction.domain1_id == id2 && interaction.domain2_id == id1)
                            });
                        
                        if !existing_interaction {
                            self.create_domain_interaction(id1, id2, distance)?;
                        }
                    }
                }
            }
        }
        
        // Update existing interactions
        let mut interactions_to_remove = Vec::new();
        for (interaction_id, interaction) in self.interactions.iter_mut() {
            interaction.duration += dt;
            
            // Check if domains still exist and are in range
            if let (Some(domain1), Some(domain2)) = (
                self.domains.get(&interaction.domain1_id),
                self.domains.get(&interaction.domain2_id)
            ) {
                let distance = self.calculate_distance(domain1, domain2);
                interaction.distance = distance;
                
                if distance > self.params.interaction_range {
                    interaction.is_active = false;
                    interactions_to_remove.push(*interaction_id);
                }
            } else {
                interactions_to_remove.push(*interaction_id);
            }
        }
        
        // Remove inactive interactions
        for interaction_id in interactions_to_remove {
            self.interactions.remove(&interaction_id);
        }
        
        Ok(())
    }
    
    /// Calculate distance between domains
    fn calculate_distance(&self, domain1: &Nanodomain, domain2: &Nanodomain) -> f64 {
        let dx = domain1.center_position.0 - domain2.center_position.0;
        let dy = domain1.center_position.1 - domain2.center_position.1;
        (dx * dx + dy * dy).sqrt()
    }
    
    /// Create domain interaction
    fn create_domain_interaction(&mut self, id1: u64, id2: u64, distance: f64) -> Result<(), BeneGesseritError> {
        let interaction_id = self.interaction_counter;
        self.interaction_counter += 1;
        
        // Determine interaction type based on domain types
        let interaction_type = if let (Some(domain1), Some(domain2)) = (self.domains.get(&id1), self.domains.get(&id2)) {
            match (&domain1.domain_type, &domain2.domain_type) {
                (NanodomainType::Signaling, NanodomainType::Receptor) => InteractionType::Cooperative,
                (NanodomainType::Enzymatic, NanodomainType::Signaling) => InteractionType::Catalytic,
                (NanodomainType::ProteinCluster, NanodomainType::ProteinCluster) => InteractionType::Attractive,
                (NanodomainType::Transient, _) => InteractionType::Repulsive,
                _ => InteractionType::Attractive,
            }
        } else {
            return Err(BeneGesseritError::InvalidState("Cannot create interaction for non-existent domains".to_string()));
        };
        
        let strength = self.calculate_interaction_strength(&interaction_type, distance);
        
        let interaction = NanodomainInteraction {
            interaction_id,
            domain1_id: id1,
            domain2_id: id2,
            interaction_type,
            strength,
            distance,
            duration: 0.0,
            is_active: true,
        };
        
        self.interactions.insert(interaction_id, interaction);
        Ok(())
    }
    
    /// Calculate interaction strength
    fn calculate_interaction_strength(&self, interaction_type: &InteractionType, distance: f64) -> f64 {
        let base_strength = match interaction_type {
            InteractionType::Attractive => 1.0,
            InteractionType::Repulsive => -0.5,
            InteractionType::Cooperative => 2.0,
            InteractionType::Competitive => -1.0,
            InteractionType::Catalytic => 1.5,
        };
        
        // Distance dependence (1/r^2)
        let distance_factor = 1.0 / (1.0 + (distance / 10.0).powi(2));
        
        base_strength * distance_factor
    }
    
    /// Check for domain dissolution
    fn check_domain_dissolution(&mut self, dt: f64) -> Result<Vec<NanodomainEvent>, BeneGesseritError> {
        let mut events = Vec::new();
        let mut domains_to_dissolve = Vec::new();
        
        for (domain_id, domain) in &self.domains {
            if matches!(domain.state, NanodomainState::Dissolving) && !domain.is_active {
                domains_to_dissolve.push(*domain_id);
            }
        }
        
        for domain_id in domains_to_dissolve {
            self.domains.remove(&domain_id);
            
            // Remove associated interactions
            self.interactions.retain(|_, interaction| {
                interaction.domain1_id != domain_id && interaction.domain2_id != domain_id
            });
            
            let event = NanodomainEvent {
                event_id: self.domain_counter,
                event_type: NanodomainEventType::Dissolution,
                domain_id,
                time: 0.0,
                trigger: NanodomainTrigger::ConcentrationChange,
                components_involved: Vec::new(),
            };
            events.push(event);
        }
        
        Ok(events)
    }
    
    /// Determine if domain nucleation should occur
    fn should_nucleate_domain(&self, domain_type: &NanodomainType) -> Result<bool, BeneGesseritError> {
        let base_probability = match domain_type {
            NanodomainType::ProteinCluster => 0.01,
            NanodomainType::LipidPatch => 0.005,
            NanodomainType::Signaling => 0.02,
            NanodomainType::Receptor => 0.015,
            NanodomainType::Enzymatic => 0.01,
            NanodomainType::Transient => 0.05,
        };
        
        // Modulate by current domain count (avoid overcrowding)
        let domain_count = self.domains.len();
        let crowding_factor = 1.0 / (1.0 + domain_count as f64 * 0.1);
        
        let nucleation_probability = base_probability * crowding_factor;
        Ok(self.random_f64() < nucleation_probability)
    }
    
    /// Select nucleation site
    fn select_nucleation_site(&self) -> (f64, f64) {
        // Random position - could be biased by existing domains
        (self.random_f64() * 1000.0, self.random_f64() * 1000.0) // nm
    }
    
    /// Select initial components for domain
    fn select_initial_components(&self, domain_type: &NanodomainType) -> HashMap<ComponentId, f64> {
        let mut components = HashMap::new();
        
        match domain_type {
            NanodomainType::ProteinCluster => {
                components.insert(ComponentId::Protein(ProteinType::Membrane), 5.0);
                components.insert(ComponentId::Lipid(LipidType::POPC), 0.3);
            }
            NanodomainType::LipidPatch => {
                components.insert(ComponentId::Lipid(LipidType::Cholesterol), 0.4);
                components.insert(ComponentId::Lipid(LipidType::Sphingomyelin), 0.3);
            }
            NanodomainType::Signaling => {
                components.insert(ComponentId::Protein(ProteinType::Kinase), 3.0);
                components.insert(ComponentId::Protein(ProteinType::Phosphatase), 2.0);
            }
            NanodomainType::Receptor => {
                components.insert(ComponentId::Protein(ProteinType::GPCR), 4.0);
            }
            NanodomainType::Enzymatic => {
                components.insert(ComponentId::Protein(ProteinType::Channel), 3.0);
            }
            NanodomainType::Transient => {
                components.insert(ComponentId::Protein(ProteinType::Membrane), 2.0);
            }
        }
        
        components
    }
    
    /// Get nanodomain statistics
    pub fn get_nanodomain_statistics(&self) -> NanodomainStatistics {
        let mut type_distribution = HashMap::new();
        let mut state_distribution = HashMap::new();
        let mut total_area = 0.0;
        let mut total_activity = 0.0;
        
        for domain in self.domains.values() {
            *type_distribution.entry(domain.domain_type.clone()).or_insert(0) += 1;
            *state_distribution.entry(domain.state.clone()).or_insert(0) += 1;
            total_area += std::f64::consts::PI * domain.radius.powi(2);
            total_activity += domain.activity_level;
        }
        
        NanodomainStatistics {
            total_domains: self.domains.len(),
            type_distribution,
            state_distribution,
            total_area,
            area_fraction: total_area / self.membrane_area,
            average_activity: if self.domains.is_empty() { 0.0 } else { total_activity / self.domains.len() as f64 },
            total_interactions: self.interactions.len(),
            formation_events: self.formation_events.len(),
        }
    }
    
    /// Set temperature
    pub fn set_temperature(&mut self, temperature: f64) {
        self.temperature = temperature;
    }
    
    // Helper methods
    fn random_f64(&self) -> f64 {
        // Placeholder - would use proper RNG
        0.5
    }
}

/// Nanodomain system statistics
#[derive(Debug, Clone)]
pub struct NanodomainStatistics {
    pub total_domains: usize,
    pub type_distribution: HashMap<NanodomainType, usize>,
    pub state_distribution: HashMap<NanodomainState, usize>,
    pub total_area: f64,
    pub area_fraction: f64,
    pub average_activity: f64,
    pub total_interactions: usize,
    pub formation_events: usize,
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_nanodomains_creation() {
        let nanodomains = Nanodomains::new();
        assert_eq!(nanodomains.domains.len(), 0);
        assert_eq!(nanodomains.interactions.len(), 0);
    }
    
    #[test]
    fn test_concentration_setting() {
        let mut nanodomains = Nanodomains::new();
        nanodomains.set_concentration(ComponentId::Protein(ProteinType::Membrane), 0.1);
        assert_eq!(nanodomains.concentrations.len(), 1);
    }
    
    #[test]
    fn test_formation_conditions() {
        let mut nanodomains = Nanodomains::new();
        nanodomains.set_concentration(ComponentId::Protein(ProteinType::Membrane), 0.1);
        
        assert!(nanodomains.check_formation_conditions(&NanodomainType::ProteinCluster));
    }
    
    #[test]
    fn test_domain_nucleation() {
        let mut nanodomains = Nanodomains::new();
        let mut components = HashMap::new();
        components.insert(ComponentId::Protein(ProteinType::Membrane), 5.0);
        
        let result = nanodomains.nucleate_domain(
            NanodomainType::ProteinCluster,
            (100.0, 100.0),
            components,
        );
        assert!(result.is_ok());
        assert_eq!(nanodomains.domains.len(), 1);
    }
}   