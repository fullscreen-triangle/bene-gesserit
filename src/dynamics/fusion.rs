//! Membrane Fusion Dynamics Module
//! 
//! This module simulates membrane fusion processes including vesicle fusion,
//! SNARE-mediated fusion, and ATP-dependent fusion machinery.

use crate::types::*;
use crate::constants::*;
use crate::error::MembraneError;
use std::collections::HashMap;

/// Unique identifier for membranes
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct MembraneId(pub u64);

/// Unique identifier for proteins
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ProteinId(pub u64);

/// Types of membrane fusion
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum FusionType {
    Vesicle,        // Vesicle-to-membrane fusion
    Membrane,       // Membrane-to-membrane fusion
    Homotypic,      // Same compartment fusion
    Heterotypic,    // Different compartment fusion
}

/// Fusion mechanism
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum FusionMechanism {
    SNAREMediated,  // SNARE protein complex
    CalciumTriggered, // Calcium-dependent fusion
    pHTriggered,    // pH-dependent fusion
    Spontaneous,    // Spontaneous membrane fusion
}

/// Fusion stage
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum FusionStage {
    Docking,        // Initial membrane contact
    Priming,        // SNARE complex assembly
    Hemifusion,     // Outer leaflet fusion only
    FullFusion,     // Complete fusion pore formation
    PoreExpansion,  // Fusion pore expansion
}

/// Fusion event record
#[derive(Debug, Clone)]
pub struct FusionEvent {
    pub event_id: u64,
    pub fusion_type: FusionType,
    pub mechanism: FusionMechanism,
    pub stage: FusionStage,
    pub membrane1_id: MembraneId,
    pub membrane2_id: MembraneId,
    pub fusion_energy: f64,
    pub pore_size: f64,
    pub time: f64,
    pub atp_consumed: f64,
    pub calcium_required: f64,
}

/// SNARE protein complex
#[derive(Debug, Clone)]
pub struct SNAREComplex {
    pub complex_id: u64,
    pub v_snare: ProteinId,  // Vesicle SNARE
    pub t_snare: ProteinId,  // Target SNARE
    pub assembly_state: f64, // 0.0 to 1.0
    pub stability: f64,
    pub force_generated: f64,
    pub is_zipped: bool,
}

/// Fusion machinery state
#[derive(Debug, Clone)]
pub struct FusionMachinery {
    pub protein_id: ProteinId,
    pub protein_type: FusionProteinType,
    pub activity: f64,
    pub calcium_sensitivity: f64,
    pub atp_requirement: f64,
    pub membrane_affinity: f64,
    pub is_active: bool,
}

/// Types of fusion proteins
#[derive(Debug, Clone, PartialEq)]
pub enum FusionProteinType {
    Synaptotagmin,  // Calcium sensor
    Complexin,      // SNARE regulator
    Munc13,         // Priming factor
    Munc18,         // SNARE chaperone
    NSF,            // SNARE disassembly
    SNAP,           // NSF cofactor
}

/// Membrane fusion simulator
pub struct Fusion {
    /// Active fusion events
    pub active_fusions: HashMap<u64, FusionEvent>,
    /// SNARE complexes
    pub snare_complexes: HashMap<u64, SNAREComplex>,
    /// Fusion machinery proteins
    pub fusion_proteins: HashMap<ProteinId, FusionMachinery>,
    /// Fusion rates by type and mechanism
    pub fusion_rates: HashMap<(FusionType, FusionMechanism), f64>,
    /// Energy barriers for fusion stages
    pub stage_barriers: HashMap<FusionStage, f64>,
    /// Completed fusion events
    pub fusion_history: Vec<FusionEvent>,
    /// Current calcium concentration
    pub calcium_concentration: f64,
    /// Current ATP concentration
    pub atp_concentration: f64,
    /// pH value
    pub ph: f64,
    /// Temperature
    pub temperature: f64,
    /// Membrane curvature effects
    pub curvature_effects: HashMap<MembraneId, f64>,
    /// Event counter
    event_counter: u64,
}

impl Default for Fusion {
    fn default() -> Self {
        let mut fusion_rates = HashMap::new();
        let mut stage_barriers = HashMap::new();
        
        // Initialize fusion rates (s^-1)
        fusion_rates.insert((FusionType::Vesicle, FusionMechanism::SNAREMediated), 1.0);
        fusion_rates.insert((FusionType::Vesicle, FusionMechanism::CalciumTriggered), 10.0);
        fusion_rates.insert((FusionType::Membrane, FusionMechanism::SNAREMediated), 0.1);
        fusion_rates.insert((FusionType::Homotypic, FusionMechanism::Spontaneous), 1e-6);
        
        // Energy barriers (kJ/mol)
        stage_barriers.insert(FusionStage::Docking, 5.0);
        stage_barriers.insert(FusionStage::Priming, 15.0);
        stage_barriers.insert(FusionStage::Hemifusion, 25.0);
        stage_barriers.insert(FusionStage::FullFusion, 35.0);
        stage_barriers.insert(FusionStage::PoreExpansion, 10.0);
        
        Self {
            active_fusions: HashMap::new(),
            snare_complexes: HashMap::new(),
            fusion_proteins: HashMap::new(),
            fusion_rates,
            stage_barriers,
            fusion_history: Vec::new(),
            calcium_concentration: 1e-7, // 100 nM resting
            atp_concentration: PHYSIOLOGICAL_ATP,
            ph: 7.4,
            temperature: PHYSIOLOGICAL_TEMPERATURE,
            curvature_effects: HashMap::new(),
            event_counter: 0,
        }
    }
}

impl Fusion {
    pub fn new() -> Self {
        Self::default()
    }
    
    /// Add fusion protein
    pub fn add_fusion_protein(&mut self, protein_id: ProteinId, protein_type: FusionProteinType) -> Result<(), MembraneError> {
        let machinery = FusionMachinery {
            protein_id,
            protein_type: protein_type.clone(),
            activity: 1.0,
            calcium_sensitivity: match protein_type {
                FusionProteinType::Synaptotagmin => 1e-6, // 1 μM
                _ => 1e-5, // 10 μM
            },
            atp_requirement: match protein_type {
                FusionProteinType::NSF => ATP_HYDROLYSIS_ENERGY,
                FusionProteinType::Munc13 => ATP_HYDROLYSIS_ENERGY * 0.5,
                _ => 0.0,
            },
            membrane_affinity: 1.0,
            is_active: true,
        };
        
        self.fusion_proteins.insert(protein_id, machinery);
        Ok(())
    }
    
    /// Create SNARE complex
    pub fn create_snare_complex(&mut self, v_snare: ProteinId, t_snare: ProteinId) -> Result<u64, MembraneError> {
        let complex_id = self.event_counter;
        self.event_counter += 1;
        
        let complex = SNAREComplex {
            complex_id,
            v_snare,
            t_snare,
            assembly_state: 0.0,
            stability: 1.0,
            force_generated: 0.0,
            is_zipped: false,
        };
        
        self.snare_complexes.insert(complex_id, complex);
        Ok(complex_id)
    }
    
    /// Initiate fusion event
    pub fn initiate_fusion(
        &mut self,
        fusion_type: FusionType,
        mechanism: FusionMechanism,
        membrane1: MembraneId,
        membrane2: MembraneId,
    ) -> Result<u64, MembraneError> {
        let event_id = self.event_counter;
        self.event_counter += 1;
        
        let event = FusionEvent {
            event_id,
            fusion_type,
            mechanism,
            stage: FusionStage::Docking,
            membrane1_id: membrane1,
            membrane2_id: membrane2,
            fusion_energy: 0.0,
            pore_size: 0.0,
            time: 0.0,
            atp_consumed: 0.0,
            calcium_required: match mechanism {
                FusionMechanism::CalciumTriggered => 1e-6,
                _ => 0.0,
            },
        };
        
        self.active_fusions.insert(event_id, event);
        Ok(event_id)
    }
    
    /// Simulate fusion dynamics
    pub fn simulate_fusion(&mut self, dt: f64) -> Result<Vec<FusionEvent>, MembraneError> {
        let mut completed_events = Vec::new();
        
        // Update protein states
        self.update_protein_states();
        
        // Update SNARE complexes
        self.update_snare_complexes(dt)?;
        
        // Process active fusion events
        let mut events_to_remove = Vec::new();
        
        for (event_id, event) in self.active_fusions.iter_mut() {
            if self.advance_fusion_stage(event, dt)? {
                // Fusion completed
                completed_events.push(event.clone());
                events_to_remove.push(*event_id);
            }
        }
        
        // Remove completed events
        for event_id in events_to_remove {
            if let Some(event) = self.active_fusions.remove(&event_id) {
                self.fusion_history.push(event);
            }
        }
        
        Ok(completed_events)
    }
    
    /// Update protein activity states
    fn update_protein_states(&mut self) {
        for protein in self.fusion_proteins.values_mut() {
            match protein.protein_type {
                FusionProteinType::Synaptotagmin => {
                    protein.is_active = self.calcium_concentration > protein.calcium_sensitivity;
                    protein.activity = if protein.is_active {
                        self.calcium_concentration / (self.calcium_concentration + protein.calcium_sensitivity)
                    } else {
                        0.0
                    };
                }
                FusionProteinType::NSF | FusionProteinType::Munc13 => {
                    protein.is_active = self.atp_concentration > protein.atp_requirement;
                    protein.activity = if protein.is_active {
                        self.atp_concentration / (self.atp_concentration + protein.atp_requirement)
                    } else {
                        0.0
                    };
                }
                _ => {
                    protein.activity = 1.0;
                    protein.is_active = true;
                }
            }
        }
    }
    
    /// Update SNARE complex assembly
    fn update_snare_complexes(&mut self, dt: f64) -> Result<(), MembraneError> {
        for complex in self.snare_complexes.values_mut() {
            if !complex.is_zipped {
                // Assembly rate depends on protein regulators
                let base_rate = 1.0; // s^-1
                let regulatory_factor = self.calculate_snare_regulation();
                let assembly_rate = base_rate * regulatory_factor;
                
                complex.assembly_state += assembly_rate * dt;
                complex.assembly_state = complex.assembly_state.min(1.0);
                
                // Check if fully assembled
                if complex.assembly_state > 0.9 {
                    complex.is_zipped = true;
                    complex.force_generated = 35.0; // pN (piconewtons)
                }
            }
        }
        
        Ok(())
    }
    
    /// Calculate SNARE regulatory effects
    fn calculate_snare_regulation(&self) -> f64 {
        let mut regulation = 1.0;
        
        for protein in self.fusion_proteins.values() {
            if protein.is_active {
                match protein.protein_type {
                    FusionProteinType::Complexin => regulation *= 0.5, // Inhibitory
                    FusionProteinType::Munc13 => regulation *= 2.0,   // Stimulatory
                    FusionProteinType::Munc18 => regulation *= 1.5,   // Stimulatory
                    _ => {}
                }
            }
        }
        
        regulation
    }
    
    /// Advance fusion through stages
    fn advance_fusion_stage(&mut self, event: &mut FusionEvent, dt: f64) -> Result<bool, MembraneError> {
        let current_stage = event.stage.clone();
        let barrier = *self.stage_barriers.get(&current_stage).unwrap_or(&20.0);
        
        // Calculate transition probability
        let rate = self.calculate_stage_transition_rate(event, &current_stage)?;
        let probability = rate * dt;
        
        if self.random_f64() < probability {
            // Advance to next stage
            event.stage = self.get_next_stage(&current_stage);
            event.fusion_energy += barrier;
            
            // Consume ATP if required
            if let Some(atp_cost) = self.get_stage_atp_cost(&event.stage) {
                event.atp_consumed += atp_cost;
            }
            
            // Update pore size for full fusion stages
            match event.stage {
                FusionStage::FullFusion => event.pore_size = 1.0, // nm
                FusionStage::PoreExpansion => event.pore_size = 10.0, // nm
                _ => {}
            }
            
            // Check if fusion is complete
            if matches!(event.stage, FusionStage::PoreExpansion) {
                return Ok(true);
            }
        }
        
        Ok(false)
    }
    
    /// Calculate stage transition rate
    fn calculate_stage_transition_rate(&self, event: &FusionEvent, stage: &FusionStage) -> Result<f64, MembraneError> {
        let base_rate = match stage {
            FusionStage::Docking => 10.0,
            FusionStage::Priming => 5.0,
            FusionStage::Hemifusion => 2.0,
            FusionStage::FullFusion => 1.0,
            FusionStage::PoreExpansion => 0.5,
        };
        
        let mut rate = base_rate;
        
        // Mechanism-specific modifiers
        match event.mechanism {
            FusionMechanism::SNAREMediated => {
                let snare_effect = self.calculate_snare_effect();
                rate *= snare_effect;
            }
            FusionMechanism::CalciumTriggered => {
                if self.calcium_concentration > event.calcium_required {
                    rate *= 10.0; // Calcium acceleration
                } else {
                    rate *= 0.1; // Calcium requirement not met
                }
            }
            FusionMechanism::pHTriggered => {
                if self.ph < 6.0 {
                    rate *= 5.0; // Acidic pH promotes fusion
                }
            }
            FusionMechanism::Spontaneous => {
                rate *= 0.01; // Very slow
            }
        }
        
        // Curvature effects
        let curvature1 = self.curvature_effects.get(&event.membrane1_id).unwrap_or(&0.0);
        let curvature2 = self.curvature_effects.get(&event.membrane2_id).unwrap_or(&0.0);
        let curvature_factor = 1.0 + (curvature1 + curvature2) * 0.1;
        rate *= curvature_factor;
        
        Ok(rate)
    }
    
    /// Calculate SNARE complex contribution to fusion
    fn calculate_snare_effect(&self) -> f64 {
        let mut total_effect = 1.0;
        
        for complex in self.snare_complexes.values() {
            if complex.is_zipped {
                total_effect += complex.force_generated / 35.0; // Normalize by max force
            }
        }
        
        total_effect
    }
    
    /// Get next fusion stage
    fn get_next_stage(&self, current: &FusionStage) -> FusionStage {
        match current {
            FusionStage::Docking => FusionStage::Priming,
            FusionStage::Priming => FusionStage::Hemifusion,
            FusionStage::Hemifusion => FusionStage::FullFusion,
            FusionStage::FullFusion => FusionStage::PoreExpansion,
            FusionStage::PoreExpansion => FusionStage::PoreExpansion, // Terminal
        }
    }
    
    /// Get ATP cost for stage
    fn get_stage_atp_cost(&self, stage: &FusionStage) -> Option<f64> {
        match stage {
            FusionStage::Priming => Some(ATP_HYDROLYSIS_ENERGY * 2.0), // SNARE priming
            FusionStage::FullFusion => Some(ATP_HYDROLYSIS_ENERGY), // Pore formation
            _ => None,
        }
    }
    
    /// Get fusion statistics
    pub fn get_fusion_statistics(&self) -> FusionStatistics {
        let mut mechanism_counts = HashMap::new();
        let mut stage_counts = HashMap::new();
        let mut total_atp_consumed = 0.0;
        
        for event in &self.fusion_history {
            *mechanism_counts.entry(event.mechanism.clone()).or_insert(0) += 1;
            *stage_counts.entry(event.stage.clone()).or_insert(0) += 1;
            total_atp_consumed += event.atp_consumed;
        }
        
        FusionStatistics {
            total_fusions: self.fusion_history.len(),
            active_fusions: self.active_fusions.len(),
            mechanism_counts,
            stage_counts,
            total_atp_consumed,
            average_pore_size: self.calculate_average_pore_size(),
            snare_complexes: self.snare_complexes.len(),
        }
    }
    
    /// Calculate average pore size
    fn calculate_average_pore_size(&self) -> f64 {
        if self.fusion_history.is_empty() {
            return 0.0;
        }
        
        let total_size: f64 = self.fusion_history.iter()
            .map(|event| event.pore_size)
            .sum();
        
        total_size / self.fusion_history.len() as f64
    }
    
    /// Set calcium concentration
    pub fn set_calcium_concentration(&mut self, concentration: f64) {
        self.calcium_concentration = concentration;
    }
    
    /// Set ATP concentration
    pub fn set_atp_concentration(&mut self, concentration: f64) {
        self.atp_concentration = concentration;
    }
    
    /// Set pH
    pub fn set_ph(&mut self, ph: f64) {
        self.ph = ph;
    }
    
    /// Set membrane curvature
    pub fn set_membrane_curvature(&mut self, membrane_id: MembraneId, curvature: f64) {
        self.curvature_effects.insert(membrane_id, curvature);
    }
    
    // Helper methods
    fn random_f64(&self) -> f64 {
        // Placeholder - would use proper RNG
        0.5
    }
}

/// Fusion statistics
#[derive(Debug, Clone)]
pub struct FusionStatistics {
    pub total_fusions: usize,
    pub active_fusions: usize,
    pub mechanism_counts: HashMap<FusionMechanism, usize>,
    pub stage_counts: HashMap<FusionStage, usize>,
    pub total_atp_consumed: f64,
    pub average_pore_size: f64,
    pub snare_complexes: usize,
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_fusion_creation() {
        let fusion = Fusion::new();
        assert_eq!(fusion.active_fusions.len(), 0);
        assert!(fusion.fusion_rates.len() > 0);
    }
    
    #[test]
    fn test_snare_complex_creation() {
        let mut fusion = Fusion::new();
        let result = fusion.create_snare_complex(ProteinId(1), ProteinId(2));
        assert!(result.is_ok());
        assert_eq!(fusion.snare_complexes.len(), 1);
    }
    
    #[test]
    fn test_fusion_initiation() {
        let mut fusion = Fusion::new();
        let result = fusion.initiate_fusion(
            FusionType::Vesicle,
            FusionMechanism::SNAREMediated,
            MembraneId(1),
            MembraneId(2),
        );
        assert!(result.is_ok());
        assert_eq!(fusion.active_fusions.len(), 1);
    }
}       