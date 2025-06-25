//! Lipid Raft Formation and Dynamics Module
//! 
//! This module simulates the formation and dynamics of lipid rafts,
//! specialized membrane domains enriched in cholesterol and sphingolipids.

use crate::types::*;
use crate::constants::*;
use crate::error::MembraneError;
use std::collections::HashMap;

/// Unique identifier for proteins
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ProteinId(pub u64);

/// Raft formation state
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum RaftState {
    Nucleating,    // Initial raft formation
    Growing,       // Raft expansion
    Mature,        // Stable raft
    Dissolving,    // Raft breakdown
    Coalescing,    // Raft-raft fusion
}

/// Raft size classification
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum RaftSize {
    Nano,          // < 50 nm
    Micro,         // 50-200 nm
    Macro,         // > 200 nm
}

/// Individual lipid raft
#[derive(Debug, Clone)]
pub struct LipidRaft {
    pub raft_id: u64,
    pub state: RaftState,
    pub size_class: RaftSize,
    pub center_position: (f64, f64),
    pub radius: f64,
    pub area: f64,
    pub cholesterol_content: f64,
    pub sphingolipid_content: f64,
    pub protein_content: HashMap<ProteinType, f64>,
    pub lipid_composition: HashMap<LipidType, f64>,
    pub order_parameter: f64,
    pub stability: f64,
    pub lifetime: f64,
    pub formation_time: f64,
    pub is_active: bool,
}

/// Raft-associated protein
#[derive(Debug, Clone)]
pub struct RaftProtein {
    pub protein_id: ProteinId,
    pub protein_type: ProteinType,
    pub raft_affinity: f64,
    pub current_raft: Option<u64>,
    pub membrane_insertion_depth: f64,
    pub lipid_interactions: HashMap<LipidType, f64>,
    pub clustering_tendency: f64,
    pub activity_in_raft: f64,
}

/// Raft formation parameters
#[derive(Debug, Clone)]
pub struct RaftFormationParams {
    pub cholesterol_threshold: f64,
    pub sphingolipid_threshold: f64,
    pub temperature_dependence: f64,
    pub nucleation_rate: f64,
    pub growth_rate: f64,
    pub dissolution_rate: f64,
    pub coalescence_probability: f64,
}

/// Lipid raft system
pub struct LipidRafts {
    /// Active rafts
    pub rafts: HashMap<u64, LipidRaft>,
    /// Raft-associated proteins
    pub raft_proteins: HashMap<ProteinId, RaftProtein>,
    /// Formation parameters
    pub formation_params: RaftFormationParams,
    /// Membrane lipid composition
    pub membrane_composition: HashMap<LipidType, f64>,
    /// Temperature
    pub temperature: f64,
    /// Membrane area
    pub membrane_area: f64,
    /// Raft formation history
    pub formation_events: Vec<RaftFormationEvent>,
    /// Raft counter
    raft_counter: u64,
}

/// Raft formation event
#[derive(Debug, Clone)]
pub struct RaftFormationEvent {
    pub event_id: u64,
    pub event_type: RaftEventType,
    pub raft_id: u64,
    pub time: f64,
    pub trigger: RaftTrigger,
    pub cholesterol_level: f64,
    pub protein_involvement: Vec<ProteinId>,
}

/// Types of raft events
#[derive(Debug, Clone, PartialEq)]
pub enum RaftEventType {
    Nucleation,
    Growth,
    Dissolution,
    Coalescence,
    ProteinRecruitment,
    ProteinRelease,
}

/// Triggers for raft formation/dissolution
#[derive(Debug, Clone, PartialEq)]
pub enum RaftTrigger {
    CholesterolIncrease,
    CholesterolDecrease,
    TemperatureChange,
    ProteinBinding,
    ProteinUnbinding,
    LipidModification,
}

impl Default for RaftFormationParams {
    fn default() -> Self {
        Self {
            cholesterol_threshold: 0.3,    // 30% cholesterol
            sphingolipid_threshold: 0.1,   // 10% sphingolipids
            temperature_dependence: 0.05,  // per degree K
            nucleation_rate: 0.1,          // s^-1
            growth_rate: 10.0,             // nm/s
            dissolution_rate: 0.01,        // s^-1
            coalescence_probability: 0.1,  // probability per collision
        }
    }
}

impl Default for LipidRafts {
    fn default() -> Self {
        let mut membrane_composition = HashMap::new();
        
        // Typical membrane composition
        membrane_composition.insert(LipidType::POPC, 0.4);
        membrane_composition.insert(LipidType::POPE, 0.2);
        membrane_composition.insert(LipidType::POPS, 0.1);
        membrane_composition.insert(LipidType::SM, 0.15);  // Changed from Sphingomyelin to SM
        membrane_composition.insert(LipidType::Cholesterol, 0.15);
        
        Self {
            rafts: HashMap::new(),
            raft_proteins: HashMap::new(),
            formation_params: RaftFormationParams::default(),
            membrane_composition,
            temperature: PHYSIOLOGICAL_TEMPERATURE,  // Changed from TEMPERATURE
            membrane_area: 1e-12, // m^2 (1 μm^2)
            formation_events: Vec::new(),
            raft_counter: 0,
        }
    }
}

impl LipidRafts {
    pub fn new() -> Self {
        Self::default()
    }
    
    /// Add raft-associated protein
    pub fn add_raft_protein(&mut self, protein_id: ProteinId, protein_type: ProteinType) -> Result<(), MembraneError> {
        let mut lipid_interactions = HashMap::new();
        let mut raft_affinity = 1.0;
        
        // Set protein-specific parameters
        match protein_type {
            ProteinType::Custom(ref name) if name == "GPI" => {
                // GPI-anchored proteins have high raft affinity
                raft_affinity = 5.0;
                lipid_interactions.insert(LipidType::Cholesterol, 2.0);
                lipid_interactions.insert(LipidType::SM, 1.5);  // Changed from Sphingomyelin to SM
            }
            ProteinType::Custom(ref name) if name == "Caveolin" => {
                // Caveolin proteins are raft markers
                raft_affinity = 8.0;
                lipid_interactions.insert(LipidType::Cholesterol, 3.0);
            }
            ProteinType::Custom(ref name) if name == "Flotillin" => {
                // Flotillin proteins associate with rafts
                raft_affinity = 3.0;
                lipid_interactions.insert(LipidType::SM, 2.0);  // Changed from Sphingomyelin to SM
            }
            _ => {
                // Default membrane protein
                raft_affinity = 0.5;
                lipid_interactions.insert(LipidType::POPC, 1.0);
            }
        }
        
        let raft_protein = RaftProtein {
            protein_id,
            protein_type,
            raft_affinity,
            current_raft: None,
            membrane_insertion_depth: 2.0, // nm
            lipid_interactions,
            clustering_tendency: 1.0,
            activity_in_raft: 1.0,
        };
        
        self.raft_proteins.insert(protein_id, raft_protein);
        Ok(())
    }
    
    /// Set membrane lipid composition
    pub fn set_lipid_composition(&mut self, composition: HashMap<LipidType, f64>) {
        self.membrane_composition = composition;
    }
    
    /// Check if conditions favor raft formation
    pub fn check_raft_formation_conditions(&self) -> bool {
        let cholesterol_level = self.membrane_composition.get(&LipidType::Cholesterol).unwrap_or(&0.0);
        let sphingolipid_level = self.membrane_composition.get(&LipidType::SM).unwrap_or(&0.0);  // Changed from Sphingomyelin to SM
        
        *cholesterol_level > self.formation_params.cholesterol_threshold &&
        *sphingolipid_level > self.formation_params.sphingolipid_threshold
    }
    
    /// Nucleate new raft
    pub fn nucleate_raft(&mut self, position: (f64, f64)) -> Result<u64, MembraneError> {
        let raft_id = self.raft_counter;
        self.raft_counter += 1;
        
        let cholesterol_content = *self.membrane_composition.get(&LipidType::Cholesterol).unwrap_or(&0.0);
        let sphingolipid_content = *self.membrane_composition.get(&LipidType::SM).unwrap_or(&0.0);  // Changed from Sphingomyelin to SM
        
        let mut lipid_composition = HashMap::new();
        // Rafts are enriched in cholesterol and sphingolipids
        lipid_composition.insert(LipidType::Cholesterol, cholesterol_content * 2.0);
        lipid_composition.insert(LipidType::SM, sphingolipid_content * 3.0);  // Changed from Sphingomyelin to SM
        lipid_composition.insert(LipidType::POPC, 0.1); // Depleted in unsaturated PC
        
        let initial_radius = 25.0_f64; // nm - typical nano-raft size
        let area = std::f64::consts::PI * initial_radius.powi(2);
        
        let raft = LipidRaft {
            raft_id,
            state: RaftState::Nucleating,
            size_class: RaftSize::Nano,
            center_position: position,
            radius: initial_radius,
            area,
            cholesterol_content: cholesterol_content * 2.0,
            sphingolipid_content: sphingolipid_content * 3.0,
            protein_content: HashMap::new(),
            lipid_composition,
            order_parameter: 0.8, // High order in rafts
            stability: 1.0,
            lifetime: 0.0,
            formation_time: 0.0, // Would be set to current time
            is_active: true,
        };
        
        self.rafts.insert(raft_id, raft);
        
        // Record formation event
        let event = RaftFormationEvent {
            event_id: self.raft_counter,
            event_type: RaftEventType::Nucleation,
            raft_id,
            time: 0.0,
            trigger: RaftTrigger::CholesterolIncrease,
            cholesterol_level: cholesterol_content,
            protein_involvement: Vec::new(),
        };
        self.formation_events.push(event);
        
        Ok(raft_id)
    }
    
    /// Simulate raft dynamics
    pub fn simulate_raft_dynamics(&mut self, dt: f64) -> Result<Vec<RaftFormationEvent>, MembraneError> {
        let mut events = Vec::new();
        
        // Update existing rafts
        self.update_raft_states(dt)?;
        
        // Check for new raft nucleation
        if self.check_raft_formation_conditions() && self.should_nucleate_raft()? {
            let position = self.select_nucleation_site();
            let _raft_id = self.nucleate_raft(position)?;
            events.push(self.formation_events.last().unwrap().clone());
        }
        
        // Update protein associations
        self.update_protein_raft_associations(dt)?;
        
        // Check for raft coalescence
        events.extend(self.check_raft_coalescence(dt)?);
        
        // Check for raft dissolution
        events.extend(self.check_raft_dissolution(dt)?);
        
        Ok(events)
    }
    
    /// Update raft states
    fn update_raft_states(&mut self, dt: f64) -> Result<(), MembraneError> {
        let raft_ids: Vec<u64> = self.rafts.keys().cloned().collect();
        
        for raft_id in raft_ids {
            if let Some(raft) = self.rafts.get_mut(&raft_id) {
                raft.lifetime += dt;
                
                match raft.state {
                    RaftState::Nucleating => {
                        // Transition to growing state
                        if raft.lifetime > 1.0 { // 1 second nucleation time
                            raft.state = RaftState::Growing;
                        }
                    }
                    RaftState::Growing => {
                        // Calculate growth rate first
                        let growth_rate = {
                            let base_rate = self.formation_params.growth_rate;
                            
                            // Temperature dependence
                            let temp_factor = 1.0 + self.formation_params.temperature_dependence * (self.temperature - PHYSIOLOGICAL_TEMPERATURE);
                            
                            // Cholesterol availability
                            let cholesterol_factor = self.membrane_composition.get(&LipidType::Cholesterol).unwrap_or(&0.0) / 0.3;
                            
                            // Size-dependent growth (smaller rafts grow faster)
                            let size_factor = 50.0 / raft.radius;
                            
                            base_rate * temp_factor * cholesterol_factor * size_factor
                        };
                        
                        // Grow the raft
                        raft.radius += growth_rate * dt;
                        raft.area = std::f64::consts::PI * raft.radius.powi(2);
                        
                        // Update size classification
                        raft.size_class = match raft.radius {
                            r if r < 50.0 => RaftSize::Nano,
                            r if r < 200.0 => RaftSize::Micro,
                            _ => RaftSize::Macro,
                        };
                        
                        // Transition to mature state
                        if raft.radius > 100.0 || raft.lifetime > 10.0 {
                            raft.state = RaftState::Mature;
                        }
                    }
                    RaftState::Mature => {
                        // Stable state - slowly decrease stability
                        raft.stability -= 0.001 * dt;
                        
                        if raft.stability < 0.5 {
                            raft.state = RaftState::Dissolving;
                        }
                    }
                    RaftState::Dissolving => {
                        // Shrink the raft
                        raft.radius -= 5.0 * dt; // nm/s dissolution rate
                        raft.area = std::f64::consts::PI * raft.radius.powi(2);
                        
                        if raft.radius < 10.0 {
                            raft.is_active = false;
                        }
                    }
                    RaftState::Coalescing => {
                        // Handled in coalescence function
                    }
                }
            }
        }
        
        // Remove inactive rafts
        self.rafts.retain(|_, raft| raft.is_active);
        
        Ok(())
    }
    
    /// Calculate raft growth rate
    fn calculate_growth_rate(&self, raft: &LipidRaft) -> Result<f64, MembraneError> {
        let base_rate = self.formation_params.growth_rate;
        
        // Temperature dependence
        let temp_factor = 1.0 + self.formation_params.temperature_dependence * (self.temperature - PHYSIOLOGICAL_TEMPERATURE);  // Changed from TEMPERATURE
        
        // Cholesterol availability
        let cholesterol_factor = self.membrane_composition.get(&LipidType::Cholesterol).unwrap_or(&0.0) / 0.3;
        
        // Size-dependent growth (smaller rafts grow faster)
        let size_factor = 50.0 / raft.radius;
        
        Ok(base_rate * temp_factor * cholesterol_factor * size_factor)
    }
    
    /// Update protein-raft associations
    fn update_protein_raft_associations(&mut self, dt: f64) -> Result<(), MembraneError> {
        // For each protein, check if it should associate with or dissociate from rafts
        let protein_ids: Vec<ProteinId> = self.raft_proteins.keys().cloned().collect();
        
        for protein_id in protein_ids {
            // First, collect the necessary data
            let (current_raft, protein_type) = {
                if let Some(protein) = self.raft_proteins.get(&protein_id) {
                    (protein.current_raft, protein.protein_type.clone())
                } else {
                    continue;
                }
            };
            
            if current_raft.is_none() {
                // Protein is not in a raft - check for association
                if let Some(raft_id) = self.find_nearby_raft(&protein_id)? {
                    let binding_probability = self.calculate_binding_probability_for_protein(&protein_id, raft_id)?;
                    if self.random_f64() < binding_probability * dt {
                        // Update protein
                        if let Some(protein) = self.raft_proteins.get_mut(&protein_id) {
                            protein.current_raft = Some(raft_id);
                        }
                        
                        // Add protein to raft
                        if let Some(raft) = self.rafts.get_mut(&raft_id) {
                            *raft.protein_content.entry(protein_type).or_insert(0.0) += 1.0;
                        }
                    }
                }
            } else {
                // Protein is in a raft - check for dissociation
                let dissociation_probability = self.calculate_dissociation_probability_for_protein(&protein_id)?;
                if self.random_f64() < dissociation_probability * dt {
                    if let Some(raft_id) = current_raft {
                        // Remove protein from raft
                        if let Some(raft) = self.rafts.get_mut(&raft_id) {
                            if let Some(count) = raft.protein_content.get_mut(&protein_type) {
                                *count -= 1.0;
                                if *count <= 0.0 {
                                    raft.protein_content.remove(&protein_type);
                                }
                            }
                        }
                    }
                    
                    // Update protein
                    if let Some(protein) = self.raft_proteins.get_mut(&protein_id) {
                        protein.current_raft = None;
                    }
                }
            }
        }
        
        Ok(())
    }
    
    /// Find nearby raft for a protein
    fn find_nearby_raft(&self, _protein_id: &ProteinId) -> Result<Option<u64>, MembraneError> {
        // Simplified - would need protein position in real implementation
        // For now, return any available raft
        for (raft_id, raft) in &self.rafts {
            if raft.is_active && matches!(raft.state, RaftState::Growing | RaftState::Mature) {
                return Ok(Some(*raft_id));
            }
        }
        Ok(None)
    }
    
    /// Calculate protein binding probability to raft
    fn calculate_binding_probability(&self, protein: &RaftProtein, raft_id: u64) -> Result<f64, MembraneError> {
        if let Some(raft) = self.rafts.get(&raft_id) {
            let base_probability = 0.1; // s^-1
            
            // Affinity factor
            let affinity_factor = protein.raft_affinity;
            
            // Lipid interaction factor
            let mut lipid_factor = 1.0;
            for (lipid_type, &interaction_strength) in &protein.lipid_interactions {
                if let Some(&lipid_fraction) = raft.lipid_composition.get(lipid_type) {
                    lipid_factor *= 1.0 + interaction_strength * lipid_fraction;
                }
            }
            
            Ok(base_probability * affinity_factor * lipid_factor)
        } else {
            Ok(0.0)
        }
    }
    
    /// Calculate protein binding probability to raft by protein ID
    fn calculate_binding_probability_for_protein(&self, protein_id: &ProteinId, raft_id: u64) -> Result<f64, MembraneError> {
        if let Some(protein) = self.raft_proteins.get(protein_id) {
            self.calculate_binding_probability(protein, raft_id)
        } else {
            Ok(0.0)
        }
    }
    
    /// Calculate protein dissociation probability from raft
    fn calculate_dissociation_probability(&self, protein: &RaftProtein) -> Result<f64, MembraneError> {
        let base_probability = 0.01; // s^-1
        
        // Lower dissociation for high-affinity proteins
        let affinity_factor = 1.0 / protein.raft_affinity;
        
        Ok(base_probability * affinity_factor)
    }
    
    /// Calculate protein dissociation probability from raft by protein ID
    fn calculate_dissociation_probability_for_protein(&self, protein_id: &ProteinId) -> Result<f64, MembraneError> {
        if let Some(protein) = self.raft_proteins.get(protein_id) {
            self.calculate_dissociation_probability(protein)
        } else {
            Ok(0.0)
        }
    }
    
    /// Check for raft coalescence
    fn check_raft_coalescence(&mut self, _dt: f64) -> Result<Vec<RaftFormationEvent>, MembraneError> {
        let mut events = Vec::new();
        let raft_ids: Vec<u64> = self.rafts.keys().cloned().collect();
        
        for i in 0..raft_ids.len() {
            for j in (i + 1)..raft_ids.len() {
                let id1 = raft_ids[i];
                let id2 = raft_ids[j];
                
                if let (Some(raft1), Some(raft2)) = (self.rafts.get(&id1), self.rafts.get(&id2)) {
                    let distance = self.calculate_raft_distance(raft1, raft2);
                    
                    if distance < (raft1.radius + raft2.radius) {
                        // Rafts are touching - check for coalescence
                        if self.random_f64() < self.formation_params.coalescence_probability {
                            let merged_raft_id = self.coalesce_rafts(id1, id2)?;
                            
                            let event = RaftFormationEvent {
                                event_id: self.raft_counter,
                                event_type: RaftEventType::Coalescence,
                                raft_id: merged_raft_id,
                                time: 0.0,
                                trigger: RaftTrigger::CholesterolIncrease,
                                cholesterol_level: *self.membrane_composition.get(&LipidType::Cholesterol).unwrap_or(&0.0),
                                protein_involvement: Vec::new(),
                            };
                            events.push(event);
                            break;
                        }
                    }
                }
            }
        }
        
        Ok(events)
    }
    
    /// Calculate distance between two rafts
    fn calculate_raft_distance(&self, raft1: &LipidRaft, raft2: &LipidRaft) -> f64 {
        let dx = raft1.center_position.0 - raft2.center_position.0;
        let dy = raft1.center_position.1 - raft2.center_position.1;
        (dx * dx + dy * dy).sqrt()
    }
    
    /// Coalesce two rafts
    fn coalesce_rafts(&mut self, id1: u64, id2: u64) -> Result<u64, MembraneError> {
        if let (Some(raft1), Some(raft2)) = (self.rafts.remove(&id1), self.rafts.remove(&id2)) {
            let total_area = raft1.area + raft2.area;
            let new_radius = (total_area / std::f64::consts::PI).sqrt();
            
            // Weighted average position
            let new_position = (
                (raft1.center_position.0 * raft1.area + raft2.center_position.0 * raft2.area) / total_area,
                (raft1.center_position.1 * raft1.area + raft2.center_position.1 * raft2.area) / total_area,
            );
            
            // Merge protein content
            let mut merged_proteins = raft1.protein_content.clone();
            for (protein_type, count) in raft2.protein_content {
                *merged_proteins.entry(protein_type).or_insert(0.0) += count;
            }
            
            // Merge lipid composition (weighted average)
            let mut merged_lipids = HashMap::new();
            for (lipid_type, frac1) in &raft1.lipid_composition {
                let frac2 = raft2.lipid_composition.get(lipid_type).unwrap_or(&0.0);
                merged_lipids.insert(*lipid_type, (frac1 * raft1.area + frac2 * raft2.area) / total_area);
            }
            
            let merged_raft = LipidRaft {
                raft_id: self.raft_counter,
                state: RaftState::Mature,
                size_class: match new_radius {
                    r if r < 50.0 => RaftSize::Nano,
                    r if r < 200.0 => RaftSize::Micro,
                    _ => RaftSize::Macro,
                },
                center_position: new_position,
                radius: new_radius,
                area: total_area,
                cholesterol_content: (raft1.cholesterol_content * raft1.area + raft2.cholesterol_content * raft2.area) / total_area,
                sphingolipid_content: (raft1.sphingolipid_content * raft1.area + raft2.sphingolipid_content * raft2.area) / total_area,
                protein_content: merged_proteins,
                lipid_composition: merged_lipids,
                order_parameter: (raft1.order_parameter + raft2.order_parameter) / 2.0,
                stability: (raft1.stability + raft2.stability) / 2.0,
                lifetime: raft1.lifetime.min(raft2.lifetime),
                formation_time: raft1.formation_time.min(raft2.formation_time),
                is_active: true,
            };
            
            let merged_id = self.raft_counter;
            self.raft_counter += 1;
            self.rafts.insert(merged_id, merged_raft);
            
            // Update protein associations
            for protein in self.raft_proteins.values_mut() {
                if protein.current_raft == Some(id1) || protein.current_raft == Some(id2) {
                    protein.current_raft = Some(merged_id);
                }
            }
            
            Ok(merged_id)
        } else {
            Err(MembraneError::ValidationError {
                field: "raft_coalescence".to_string(),
                reason: "Cannot coalesce non-existent rafts".to_string(),
            })
        }
    }
    
    /// Check for raft dissolution
    fn check_raft_dissolution(&mut self, _dt: f64) -> Result<Vec<RaftFormationEvent>, MembraneError> {
        let mut events = Vec::new();
        let mut rafts_to_dissolve = Vec::new();
        
        for (raft_id, raft) in &self.rafts {
            if matches!(raft.state, RaftState::Dissolving) && !raft.is_active {
                rafts_to_dissolve.push(*raft_id);
            }
        }
        
        for raft_id in rafts_to_dissolve {
            self.rafts.remove(&raft_id);
            
            // Release associated proteins
            for protein in self.raft_proteins.values_mut() {
                if protein.current_raft == Some(raft_id) {
                    protein.current_raft = None;
                }
            }
            
            let event = RaftFormationEvent {
                event_id: self.raft_counter,
                event_type: RaftEventType::Dissolution,
                raft_id,
                time: 0.0,
                trigger: RaftTrigger::CholesterolDecrease,
                cholesterol_level: *self.membrane_composition.get(&LipidType::Cholesterol).unwrap_or(&0.0),
                protein_involvement: Vec::new(),
            };
            events.push(event);
        }
        
        Ok(events)
    }
    
    /// Determine if raft nucleation should occur
    fn should_nucleate_raft(&self) -> Result<bool, MembraneError> {
        let cholesterol_level = self.membrane_composition.get(&LipidType::Cholesterol).unwrap_or(&0.0);
        let supersaturation = (cholesterol_level - self.formation_params.cholesterol_threshold) / self.formation_params.cholesterol_threshold;
        
        let nucleation_probability = self.formation_params.nucleation_rate * supersaturation.max(0.0);
        Ok(self.random_f64() < nucleation_probability)
    }
    
    /// Select nucleation site
    fn select_nucleation_site(&self) -> (f64, f64) {
        // Random position - in real implementation would consider local cholesterol concentration
        (self.random_f64() * 1000.0, self.random_f64() * 1000.0) // nm
    }
    
    /// Get raft statistics
    pub fn get_raft_statistics(&self) -> RaftStatistics {
        let mut size_distribution = HashMap::new();
        let mut state_distribution = HashMap::new();
        let mut total_area = 0.0;
        let mut total_proteins = 0;
        
        for raft in self.rafts.values() {
            *size_distribution.entry(raft.size_class.clone()).or_insert(0) += 1;
            *state_distribution.entry(raft.state.clone()).or_insert(0) += 1;
            total_area += raft.area;
            total_proteins += raft.protein_content.values().sum::<f64>() as usize;
        }
        
        RaftStatistics {
            total_rafts: self.rafts.len(),
            size_distribution,
            state_distribution,
            total_area,
            area_fraction: total_area / self.membrane_area,
            average_cholesterol: self.calculate_average_cholesterol(),
            total_raft_proteins: total_proteins,
            formation_events: self.formation_events.len(),
        }
    }
    
    /// Calculate average cholesterol content in rafts
    fn calculate_average_cholesterol(&self) -> f64 {
        if self.rafts.is_empty() {
            return 0.0;
        }
        
        let total_cholesterol: f64 = self.rafts.values().map(|raft| raft.cholesterol_content).sum();
        total_cholesterol / self.rafts.len() as f64
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

/// Raft system statistics
#[derive(Debug, Clone)]
pub struct RaftStatistics {
    pub total_rafts: usize,
    pub size_distribution: HashMap<RaftSize, usize>,
    pub state_distribution: HashMap<RaftState, usize>,
    pub total_area: f64,
    pub area_fraction: f64,
    pub average_cholesterol: f64,
    pub total_raft_proteins: usize,
    pub formation_events: usize,
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_lipid_rafts_creation() {
        let rafts = LipidRafts::new();
        assert_eq!(rafts.rafts.len(), 0);
        assert!(rafts.membrane_composition.len() > 0);
    }
    
    #[test]
    fn test_raft_protein_addition() {
        let mut rafts = LipidRafts::new();
        let result = rafts.add_raft_protein(ProteinId(1), ProteinType::Custom("GPI".to_string()));
        assert!(result.is_ok());
        assert_eq!(rafts.raft_proteins.len(), 1);
    }
    
    #[test]
    fn test_raft_formation_conditions() {
        let mut rafts = LipidRafts::new();
        
        // Set high cholesterol
        let mut composition = HashMap::new();
        composition.insert(LipidType::Cholesterol, 0.4);
        composition.insert(LipidType::SM, 0.2);  // Changed from Sphingomyelin to SM
        rafts.set_lipid_composition(composition);
        
        assert!(rafts.check_raft_formation_conditions());
    }
    
    #[test]
    fn test_raft_nucleation() {
        let mut rafts = LipidRafts::new();
        
        // Set favorable conditions
        let mut composition = HashMap::new();
        composition.insert(LipidType::Cholesterol, 0.4);
        composition.insert(LipidType::SM, 0.2);  // Changed from Sphingomyelin to SM
        rafts.set_lipid_composition(composition);
        
        let result = rafts.nucleate_raft((100.0, 100.0));
        assert!(result.is_ok());
        assert_eq!(rafts.rafts.len(), 1);
    }
}