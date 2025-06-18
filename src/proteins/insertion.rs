//! Membrane protein insertion mechanisms
//!
//! This module implements the detailed processes of protein insertion
//! into membranes, including:
//! - Co-translational and post-translational insertion
//! - Signal recognition particle (SRP) pathway
//! - Translocon-mediated insertion
//! - ATP-dependent insertion machinery

use crate::{
    constants::*,
    error::{MembraneError, Result},
    types::*,
};
use std::collections::HashMap;

/// Protein insertion pathway types
#[derive(Debug, Clone, PartialEq)]
pub enum InsertionPathway {
    /// Co-translational insertion (during translation)
    Cotranslational,
    /// Post-translational insertion (after translation)
    Posttranslational,
    /// Spontaneous insertion (for simple proteins)
    Spontaneous,
    /// Assisted insertion (with chaperones)
    Assisted,
}

/// Signal sequences for targeting proteins to membranes
#[derive(Debug, Clone)]
pub struct SignalSequence {
    /// Signal sequence type
    pub sequence_type: SignalType,
    /// Hydrophobicity index
    pub hydrophobicity: f64,
    /// Length of signal sequence
    pub length: u8,
    /// Recognition efficiency
    pub recognition_efficiency: f64,
}

/// Types of signal sequences
#[derive(Debug, Clone, PartialEq)]
pub enum SignalType {
    /// N-terminal signal sequence
    NTerminal,
    /// Internal signal sequence
    Internal,
    /// C-terminal signal sequence
    CTerminal,
    /// Multiple signal sequences
    Multiple,
}

/// Protein insertion state
#[derive(Debug, Clone)]
pub struct InsertionState {
    /// Current insertion stage
    pub stage: InsertionStage,
    /// Insertion pathway being used
    pub pathway: InsertionPathway,
    /// Signal sequence information
    pub signal_sequence: Option<SignalSequence>,
    /// Translocon being used
    pub translocon_id: Option<String>,
    /// Insertion progress (0.0 to 1.0)
    pub progress: f64,
    /// Energy consumed so far
    pub energy_consumed: f64,
    /// Time in current stage
    pub stage_time: f64,
}

/// Stages of protein insertion
#[derive(Debug, Clone, PartialEq)]
pub enum InsertionStage {
    /// Signal recognition
    SignalRecognition,
    /// Targeting to membrane
    Targeting,
    /// Translocon binding
    TransloconBinding,
    /// Membrane insertion
    Insertion,
    /// Signal sequence cleavage
    SignalCleavage,
    /// Folding and maturation
    Folding,
    /// Insertion complete
    Complete,
}

/// Insertion machinery complex
#[derive(Debug, Clone)]
pub struct InsertionMachinery {
    /// Signal recognition particle (SRP)
    pub srp: SignalRecognitionParticle,
    /// SRP receptor
    pub srp_receptor: SrpReceptor,
    /// Translocon complexes
    pub translocons: HashMap<String, TransloconComplex>,
    /// Insertion chaperones
    pub chaperones: HashMap<String, InsertionChaperone>,
    /// Signal peptidases
    pub peptidases: HashMap<String, SignalPeptidase>,
}

/// Signal recognition particle
#[derive(Debug, Clone)]
pub struct SignalRecognitionParticle {
    /// SRP binding state
    pub bound_ribosome: Option<String>,
    /// GTP binding state
    pub gtp_bound: bool,
    /// Activity state
    pub active: bool,
    /// ATP consumption rate
    pub atp_consumption_rate: f64,
}

/// SRP receptor
#[derive(Debug, Clone)]
pub struct SrpReceptor {
    /// Membrane-bound state
    pub membrane_bound: bool,
    /// SRP binding state
    pub srp_bound: bool,
    /// GTP binding state
    pub gtp_bound: bool,
}

/// Translocon complex for protein insertion
#[derive(Debug, Clone)]
pub struct TransloconComplex {
    /// Translocon type (Sec61, etc.)
    pub translocon_type: String,
    /// Current occupancy
    pub occupied: bool,
    /// Protein being inserted
    pub inserting_protein: Option<String>,
    /// Channel state (open/closed)
    pub channel_open: bool,
    /// Insertion rate
    pub insertion_rate: f64,
    /// ATP requirement per amino acid
    pub atp_per_residue: f64,
}

/// Chaperone assisting protein insertion
#[derive(Debug, Clone)]
pub struct InsertionChaperone {
    /// Chaperone type
    pub chaperone_type: String,
    /// Bound client proteins
    pub bound_clients: Vec<String>,
    /// ATP consumption rate
    pub atp_consumption_rate: f64,
    /// Folding assistance efficiency
    pub folding_efficiency: f64,
}

/// Signal peptidase for cleaving signal sequences
#[derive(Debug, Clone)]
pub struct SignalPeptidase {
    /// Peptidase type
    pub peptidase_type: String,
    /// Cleavage specificity
    pub cleavage_specificity: Vec<SignalType>,
    /// Cleavage rate
    pub cleavage_rate: f64,
}

/// Protein insertion manager
#[derive(Debug, Clone)]
pub struct InsertionManager {
    /// Proteins currently being inserted
    pub inserting_proteins: HashMap<String, InsertionState>,
    /// Insertion machinery
    pub machinery: InsertionMachinery,
    /// Insertion statistics
    pub statistics: InsertionStatistics,
}

/// Insertion statistics
#[derive(Debug, Clone)]
pub struct InsertionStatistics {
    /// Total proteins inserted
    pub total_inserted: u64,
    /// Insertion success rate
    pub success_rate: f64,
    /// Average insertion time
    pub avg_insertion_time: f64,
    /// Total ATP consumed
    pub total_atp_consumed: f64,
    /// Insertion failures
    pub failure_count: u64,
}

impl InsertionManager {
    /// Create a new insertion manager
    pub fn new() -> Self {
        Self {
            inserting_proteins: HashMap::new(),
            machinery: InsertionMachinery::new(),
            statistics: InsertionStatistics::new(),
        }
    }
    
    /// Initiate protein insertion
    pub fn initiate_insertion(&mut self, protein_id: String, protein_type: ProteinType,
                            signal_sequence: Option<SignalSequence>,
                            atp_availability: f64) -> Result<()> {
        // Determine insertion pathway
        let pathway = self.determine_insertion_pathway(&protein_type, &signal_sequence)?;
        
        // Check ATP requirements
        let required_atp = self.calculate_insertion_energy_requirement(&protein_type, &pathway);
        if atp_availability < required_atp {
            return Err(MembraneError::InsufficientAtp {
                required: required_atp,
                available: atp_availability,
            });
        }
        
        // Create insertion state
        let insertion_state = InsertionState {
            stage: InsertionStage::SignalRecognition,
            pathway,
            signal_sequence,
            translocon_id: None,
            progress: 0.0,
            energy_consumed: 0.0,
            stage_time: 0.0,
        };
        
        self.inserting_proteins.insert(protein_id, insertion_state);
        Ok(())
    }
    
    /// Determine appropriate insertion pathway
    fn determine_insertion_pathway(&self, protein_type: &ProteinType, 
                                 signal_sequence: &Option<SignalSequence>) -> Result<InsertionPathway> {
        match signal_sequence {
            Some(signal) => {
                match signal.sequence_type {
                    SignalType::NTerminal => Ok(InsertionPathway::Cotranslational),
                    SignalType::Internal => Ok(InsertionPathway::Posttranslational),
                    SignalType::CTerminal => Ok(InsertionPathway::Posttranslational),
                    SignalType::Multiple => Ok(InsertionPathway::Assisted),
                }
            },
            None => {
                // Simple proteins might insert spontaneously
                match protein_type {
                    ProteinType::VGSC | ProteinType::VGKC => Ok(InsertionPathway::Spontaneous),
                    _ => Ok(InsertionPathway::Assisted),
                }
            }
        }
    }
    
    /// Calculate energy requirement for insertion
    fn calculate_insertion_energy_requirement(&self, protein_type: &ProteinType, 
                                            pathway: &InsertionPathway) -> f64 {
        let base_energy = match protein_type {
            ProteinType::NaKATPase => 20.0 * ATP_ENERGY_PER_MOLECULE,
            ProteinType::CaATPase => 18.0 * ATP_ENERGY_PER_MOLECULE,
            ProteinType::VGSC => 12.0 * ATP_ENERGY_PER_MOLECULE,
            ProteinType::VGKC => 10.0 * ATP_ENERGY_PER_MOLECULE,
            ProteinType::VGCC => 14.0 * ATP_ENERGY_PER_MOLECULE,
            _ => 8.0 * ATP_ENERGY_PER_MOLECULE,
        };
        
        let pathway_modifier = match pathway {
            InsertionPathway::Cotranslational => 1.0,
            InsertionPathway::Posttranslational => 1.2,
            InsertionPathway::Spontaneous => 0.5,
            InsertionPathway::Assisted => 1.5,
        };
        
        base_energy * pathway_modifier
    }
    
    /// Update insertion processes
    pub fn update_insertions(&mut self, dt: f64, atp_pool: &mut AtpMeasurement) -> Result<()> {
        let mut completed_insertions = Vec::new();
        let mut failed_insertions = Vec::new();
        
        for (protein_id, insertion_state) in &mut self.inserting_proteins {
            insertion_state.stage_time += dt;
            
            match insertion_state.stage {
                InsertionStage::SignalRecognition => {
                    self.process_signal_recognition(protein_id, insertion_state, dt, atp_pool)?;
                },
                InsertionStage::Targeting => {
                    self.process_targeting(protein_id, insertion_state, dt, atp_pool)?;
                },
                InsertionStage::TransloconBinding => {
                    self.process_translocon_binding(protein_id, insertion_state, dt, atp_pool)?;
                },
                InsertionStage::Insertion => {
                    self.process_insertion(protein_id, insertion_state, dt, atp_pool)?;
                },
                InsertionStage::SignalCleavage => {
                    self.process_signal_cleavage(protein_id, insertion_state, dt, atp_pool)?;
                },
                InsertionStage::Folding => {
                    self.process_folding(protein_id, insertion_state, dt, atp_pool)?;
                },
                InsertionStage::Complete => {
                    completed_insertions.push(protein_id.clone());
                },
            }
            
            // Check for insertion failure
            if self.has_insertion_failed(insertion_state) {
                failed_insertions.push(protein_id.clone());
            }
        }
        
        // Remove completed and failed insertions
        for protein_id in completed_insertions {
            self.inserting_proteins.remove(&protein_id);
            self.statistics.total_inserted += 1;
        }
        
        for protein_id in failed_insertions {
            self.inserting_proteins.remove(&protein_id);
            self.statistics.failure_count += 1;
        }
        
        // Update statistics
        self.update_statistics();
        
        Ok(())
    }
    
    /// Process signal recognition stage
    fn process_signal_recognition(&mut self, protein_id: &str, state: &mut InsertionState,
                                dt: f64, atp_pool: &mut AtpMeasurement) -> Result<()> {
        let recognition_time = 5.0; // 5 seconds
        let atp_cost = 2.0 * ATP_ENERGY_PER_MOLECULE;
        
        if state.stage_time > recognition_time {
            if atp_pool.molecules as f64 * ATP_ENERGY_PER_MOLECULE >= atp_cost {
                // Consume ATP for signal recognition
                let atp_molecules = (atp_cost / ATP_ENERGY_PER_MOLECULE) as u64;
                atp_pool.molecules -= atp_molecules;
                state.energy_consumed += atp_cost;
                
                // Advance to targeting stage
                state.stage = InsertionStage::Targeting;
                state.stage_time = 0.0;
                state.progress = 0.1;
                
                // Bind SRP if using cotranslational pathway
                if state.pathway == InsertionPathway::Cotranslational {
                    self.machinery.srp.bound_ribosome = Some(protein_id.to_string());
                    self.machinery.srp.active = true;
                }
            }
        }
        
        Ok(())
    }
    
    /// Process targeting stage
    fn process_targeting(&mut self, protein_id: &str, state: &mut InsertionState,
                       dt: f64, atp_pool: &mut AtpMeasurement) -> Result<()> {
        let targeting_time = 10.0; // 10 seconds
        let atp_cost = 3.0 * ATP_ENERGY_PER_MOLECULE;
        
        if state.stage_time > targeting_time {
            if atp_pool.molecules as f64 * ATP_ENERGY_PER_MOLECULE >= atp_cost {
                let atp_molecules = (atp_cost / ATP_ENERGY_PER_MOLECULE) as u64;
                atp_pool.molecules -= atp_molecules;
                state.energy_consumed += atp_cost;
                
                state.stage = InsertionStage::TransloconBinding;
                state.stage_time = 0.0;
                state.progress = 0.3;
            }
        }
        
        Ok(())
    }
    
    /// Process translocon binding stage
    fn process_translocon_binding(&mut self, protein_id: &str, state: &mut InsertionState,
                                dt: f64, atp_pool: &mut AtpMeasurement) -> Result<()> {
        // Find available translocon
        if state.translocon_id.is_none() {
            if let Some(translocon_id) = self.find_available_translocon()? {
                state.translocon_id = Some(translocon_id.clone());
                
                // Bind to translocon
                if let Some(translocon) = self.machinery.translocons.get_mut(&translocon_id) {
                    translocon.occupied = true;
                    translocon.inserting_protein = Some(protein_id.to_string());
                    translocon.channel_open = true;
                }
            }
        }
        
        let binding_time = 3.0; // 3 seconds
        let atp_cost = 1.0 * ATP_ENERGY_PER_MOLECULE;
        
        if state.stage_time > binding_time && state.translocon_id.is_some() {
            if atp_pool.molecules as f64 * ATP_ENERGY_PER_MOLECULE >= atp_cost {
                let atp_molecules = (atp_cost / ATP_ENERGY_PER_MOLECULE) as u64;
                atp_pool.molecules -= atp_molecules;
                state.energy_consumed += atp_cost;
                
                state.stage = InsertionStage::Insertion;
                state.stage_time = 0.0;
                state.progress = 0.5;
            }
        }
        
        Ok(())
    }
    
    /// Process insertion stage
    fn process_insertion(&mut self, protein_id: &str, state: &mut InsertionState,
                       dt: f64, atp_pool: &mut AtpMeasurement) -> Result<()> {
        let insertion_time = 20.0; // 20 seconds for full insertion
        let atp_cost_per_second = 0.5 * ATP_ENERGY_PER_MOLECULE;
        
        let atp_cost = atp_cost_per_second * dt;
        let atp_molecules = (atp_cost / ATP_ENERGY_PER_MOLECULE) as u64;
        
        if atp_pool.molecules >= atp_molecules {
            atp_pool.molecules -= atp_molecules;
            state.energy_consumed += atp_cost;
            
            // Update insertion progress
            let progress_increment = dt / insertion_time;
            state.progress += progress_increment * 0.3; // 30% of total progress
            
            if state.stage_time > insertion_time {
                state.stage = InsertionStage::SignalCleavage;
                state.stage_time = 0.0;
                state.progress = 0.8;
            }
        }
        
        Ok(())
    }
    
    /// Process signal cleavage stage
    fn process_signal_cleavage(&mut self, protein_id: &str, state: &mut InsertionState,
                             dt: f64, atp_pool: &mut AtpMeasurement) -> Result<()> {
        // Only cleave if signal sequence is present
        if state.signal_sequence.is_some() {
            let cleavage_time = 2.0; // 2 seconds
            let atp_cost = 1.0 * ATP_ENERGY_PER_MOLECULE;
            
            if state.stage_time > cleavage_time {
                if atp_pool.molecules as f64 * ATP_ENERGY_PER_MOLECULE >= atp_cost {
                    let atp_molecules = (atp_cost / ATP_ENERGY_PER_MOLECULE) as u64;
                    atp_pool.molecules -= atp_molecules;
                    state.energy_consumed += atp_cost;
                    
                    state.stage = InsertionStage::Folding;
                    state.stage_time = 0.0;
                    state.progress = 0.9;
                }
            }
        } else {
            // Skip cleavage if no signal sequence
            state.stage = InsertionStage::Folding;
            state.stage_time = 0.0;
            state.progress = 0.9;
        }
        
        Ok(())
    }
    
    /// Process folding stage
    fn process_folding(&mut self, protein_id: &str, state: &mut InsertionState,
                     dt: f64, atp_pool: &mut AtpMeasurement) -> Result<()> {
        let folding_time = 15.0; // 15 seconds
        let atp_cost = 2.0 * ATP_ENERGY_PER_MOLECULE;
        
        if state.stage_time > folding_time {
            if atp_pool.molecules as f64 * ATP_ENERGY_PER_MOLECULE >= atp_cost {
                let atp_molecules = (atp_cost / ATP_ENERGY_PER_MOLECULE) as u64;
                atp_pool.molecules -= atp_molecules;
                state.energy_consumed += atp_cost;
                
                // Release translocon
                if let Some(translocon_id) = &state.translocon_id {
                    if let Some(translocon) = self.machinery.translocons.get_mut(translocon_id) {
                        translocon.occupied = false;
                        translocon.inserting_protein = None;
                        translocon.channel_open = false;
                    }
                }
                
                state.stage = InsertionStage::Complete;
                state.progress = 1.0;
            }
        }
        
        Ok(())
    }
    
    /// Find available translocon
    fn find_available_translocon(&self) -> Result<Option<String>> {
        for (id, translocon) in &self.machinery.translocons {
            if !translocon.occupied {
                return Ok(Some(id.clone()));
            }
        }
        Ok(None)
    }
    
    /// Check if insertion has failed
    fn has_insertion_failed(&self, state: &InsertionState) -> bool {
        // Insertion fails if stuck in a stage too long
        match state.stage {
            InsertionStage::SignalRecognition => state.stage_time > 30.0,
            InsertionStage::Targeting => state.stage_time > 60.0,
            InsertionStage::TransloconBinding => state.stage_time > 30.0,
            InsertionStage::Insertion => state.stage_time > 120.0,
            InsertionStage::SignalCleavage => state.stage_time > 15.0,
            InsertionStage::Folding => state.stage_time > 90.0,
            InsertionStage::Complete => false,
        }
    }
    
    /// Update insertion statistics
    fn update_statistics(&mut self) {
        let total_attempts = self.statistics.total_inserted + self.statistics.failure_count;
        if total_attempts > 0 {
            self.statistics.success_rate = self.statistics.total_inserted as f64 / total_attempts as f64;
        }
        
        // Calculate average insertion time from completed insertions
        if !self.inserting_proteins.is_empty() {
            let total_time: f64 = self.inserting_proteins.values()
                .map(|state| state.stage_time)
                .sum();
            self.statistics.avg_insertion_time = total_time / self.inserting_proteins.len() as f64;
        }
    }
}

impl InsertionMachinery {
    /// Create new insertion machinery
    pub fn new() -> Self {
        let mut translocons = HashMap::new();
        translocons.insert("Sec61".to_string(), TransloconComplex::new("Sec61".to_string()));
        translocons.insert("Oxa1".to_string(), TransloconComplex::new("Oxa1".to_string()));
        
        let mut chaperones = HashMap::new();
        chaperones.insert("BiP".to_string(), InsertionChaperone::new("BiP".to_string()));
        chaperones.insert("Hsp70".to_string(), InsertionChaperone::new("Hsp70".to_string()));
        
        let mut peptidases = HashMap::new();
        peptidases.insert("SPase".to_string(), SignalPeptidase::new("SPase".to_string()));
        
        Self {
            srp: SignalRecognitionParticle::new(),
            srp_receptor: SrpReceptor::new(),
            translocons,
            chaperones,
            peptidases,
        }
    }
}

impl SignalRecognitionParticle {
    pub fn new() -> Self {
        Self {
            bound_ribosome: None,
            gtp_bound: false,
            active: false,
            atp_consumption_rate: 1.0,
        }
    }
}

impl SrpReceptor {
    pub fn new() -> Self {
        Self {
            membrane_bound: true,
            srp_bound: false,
            gtp_bound: false,
        }
    }
}

impl TransloconComplex {
    pub fn new(translocon_type: String) -> Self {
        Self {
            translocon_type,
            occupied: false,
            inserting_protein: None,
            channel_open: false,
            insertion_rate: 0.05, // proteins/second
            atp_per_residue: 0.1,
        }
    }
}

impl InsertionChaperone {
    pub fn new(chaperone_type: String) -> Self {
        Self {
            chaperone_type,
            bound_clients: Vec::new(),
            atp_consumption_rate: 2.0,
            folding_efficiency: 0.85,
        }
    }
}

impl SignalPeptidase {
    pub fn new(peptidase_type: String) -> Self {
        Self {
            peptidase_type,
            cleavage_specificity: vec![SignalType::NTerminal],
            cleavage_rate: 1.0,
        }
    }
}

impl InsertionStatistics {
    pub fn new() -> Self {
        Self {
            total_inserted: 0,
            success_rate: 1.0,
            avg_insertion_time: 0.0,
            total_atp_consumed: 0.0,
            failure_count: 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_insertion_manager_creation() {
        let manager = InsertionManager::new();
        assert!(manager.inserting_proteins.is_empty());
        assert!(!manager.machinery.translocons.is_empty());
    }
    
    #[test]
    fn test_pathway_determination() {
        let manager = InsertionManager::new();
        let signal = SignalSequence {
            sequence_type: SignalType::NTerminal,
            hydrophobicity: 0.8,
            length: 20,
            recognition_efficiency: 0.9,
        };
        
        let pathway = manager.determine_insertion_pathway(
            &ProteinType::NaKATPase, 
            &Some(signal)
        ).unwrap();
        
        assert_eq!(pathway, InsertionPathway::Cotranslational);
    }
}
