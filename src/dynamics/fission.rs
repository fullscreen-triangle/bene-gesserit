//! Membrane Fission Dynamics Module
//! 
//! This module simulates membrane fission processes including vesicle budding,
//! dynamin-mediated scission, and ATP-dependent fission machinery.

use crate::types::*;
use crate::constants::*;
use crate::error::BeneGesseritError;
use std::collections::HashMap;

/// Types of membrane fission
#[derive(Debug, Clone, PartialEq)]
pub enum FissionType {
    Endocytosis,    // Clathrin-mediated endocytosis
    Exocytosis,     // Vesicle budding from organelles
    Mitochondrial,  // Mitochondrial division
    Peroxisomal,    // Peroxisome division
    Spontaneous,    // Spontaneous membrane budding
}

/// Fission mechanism
#[derive(Debug, Clone, PartialEq)]
pub enum FissionMechanism {
    DynaminMediated,  // Dynamin GTPase activity
    ESCRTMediated,    // ESCRT complex machinery
    ClathrinMediated, // Clathrin coat-mediated
    CurvatureInduced, // Curvature-driven fission
    Spontaneous,      // Spontaneous scission
}

/// Fission stage
#[derive(Debug, Clone, PartialEq)]
pub enum FissionStage {
    Initiation,     // Membrane deformation begins
    Budding,        // Bud formation
    Constriction,   // Neck constriction
    Scission,       // Final membrane severing
    Release,        // Vesicle release
}

/// Fission event record
#[derive(Debug, Clone)]
pub struct FissionEvent {
    pub event_id: u64,
    pub fission_type: FissionType,
    pub mechanism: FissionMechanism,
    pub stage: FissionStage,
    pub parent_membrane: MembraneId,
    pub vesicle_id: Option<VesicleId>,
    pub neck_diameter: f64,
    pub vesicle_size: f64,
    pub energy_cost: f64,
    pub time: f64,
    pub atp_consumed: f64,
    pub gtp_consumed: f64,
}

/// Dynamin assembly state
#[derive(Debug, Clone)]
pub struct DynaminAssembly {
    pub assembly_id: u64,
    pub dynamin_proteins: Vec<ProteinId>,
    pub helix_pitch: f64,
    pub constriction_force: f64,
    pub gtp_binding_state: f64,
    pub assembly_completeness: f64,
    pub is_active: bool,
}

/// ESCRT complex state
#[derive(Debug, Clone)]
pub struct ESCRTComplex {
    pub complex_id: u64,
    pub escrt_i: Vec<ProteinId>,
    pub escrt_ii: Vec<ProteinId>,
    pub escrt_iii: Vec<ProteinId>,
    pub vps4_atpase: Option<ProteinId>,
    pub assembly_state: f64,
    pub scission_activity: f64,
    pub is_assembled: bool,
}

/// Clathrin coat state
#[derive(Debug, Clone)]
pub struct ClathrinCoat {
    pub coat_id: u64,
    pub clathrin_proteins: Vec<ProteinId>,
    pub adaptors: Vec<ProteinId>,
    pub coat_curvature: f64,
    pub polymerization_state: f64,
    pub is_complete: bool,
}

/// Fission machinery state
#[derive(Debug, Clone)]
pub struct FissionMachinery {
    pub protein_id: ProteinId,
    pub protein_type: FissionProteinType,
    pub activity: f64,
    pub nucleotide_binding: f64,
    pub membrane_affinity: f64,
    pub oligomerization_state: f64,
    pub is_active: bool,
}

/// Types of fission proteins
#[derive(Debug, Clone, PartialEq)]
pub enum FissionProteinType {
    Dynamin,        // GTPase for scission
    Amphiphysin,    // BAR domain protein
    Endophilin,     // N-BAR protein
    SNX9,           // Sorting nexin
    Clathrin,       // Coat protein
    AP2,            // Adaptor protein
    ESCRT_I,        // ESCRT-I components
    ESCRT_II,       // ESCRT-II components
    ESCRT_III,      // ESCRT-III components
    VPS4,           // AAA+ ATPase
}

/// Membrane fission simulator
pub struct Fission {
    /// Active fission events
    pub active_fissions: HashMap<u64, FissionEvent>,
    /// Dynamin assemblies
    pub dynamin_assemblies: HashMap<u64, DynaminAssembly>,
    /// ESCRT complexes
    pub escrt_complexes: HashMap<u64, ESCRTComplex>,
    /// Clathrin coats
    pub clathrin_coats: HashMap<u64, ClathrinCoat>,
    /// Fission machinery proteins
    pub fission_proteins: HashMap<ProteinId, FissionMachinery>,
    /// Fission rates by type and mechanism
    pub fission_rates: HashMap<(FissionType, FissionMechanism), f64>,
    /// Energy barriers for fission stages
    pub stage_barriers: HashMap<FissionStage, f64>,
    /// Completed fission events
    pub fission_history: Vec<FissionEvent>,
    /// Current ATP concentration
    pub atp_concentration: f64,
    /// Current GTP concentration
    pub gtp_concentration: f64,
    /// Membrane tension effects
    pub membrane_tensions: HashMap<MembraneId, f64>,
    /// Temperature
    pub temperature: f64,
    /// Event counter
    event_counter: u64,
}

impl Default for Fission {
    fn default() -> Self {
        let mut fission_rates = HashMap::new();
        let mut stage_barriers = HashMap::new();
        
        // Initialize fission rates (s^-1)
        fission_rates.insert((FissionType::Endocytosis, FissionMechanism::DynaminMediated), 0.1);
        fission_rates.insert((FissionType::Endocytosis, FissionMechanism::ClathrinMediated), 0.05);
        fission_rates.insert((FissionType::Exocytosis, FissionMechanism::ESCRTMediated), 0.2);
        fission_rates.insert((FissionType::Mitochondrial, FissionMechanism::DynaminMediated), 0.01);
        fission_rates.insert((FissionType::Spontaneous, FissionMechanism::CurvatureInduced), 1e-6);
        
        // Energy barriers (kJ/mol)
        stage_barriers.insert(FissionStage::Initiation, 10.0);
        stage_barriers.insert(FissionStage::Budding, 20.0);
        stage_barriers.insert(FissionStage::Constriction, 30.0);
        stage_barriers.insert(FissionStage::Scission, 40.0);
        stage_barriers.insert(FissionStage::Release, 5.0);
        
        Self {
            active_fissions: HashMap::new(),
            dynamin_assemblies: HashMap::new(),
            escrt_complexes: HashMap::new(),
            clathrin_coats: HashMap::new(),
            fission_proteins: HashMap::new(),
            fission_rates,
            stage_barriers,
            fission_history: Vec::new(),
            atp_concentration: ATP_CONCENTRATION,
            gtp_concentration: 1e-4, // 100 μM GTP
            membrane_tensions: HashMap::new(),
            temperature: TEMPERATURE,
            event_counter: 0,
        }
    }
}

impl Fission {
    pub fn new() -> Self {
        Self::default()
    }
    
    /// Add fission protein
    pub fn add_fission_protein(&mut self, protein_id: ProteinId, protein_type: FissionProteinType) -> Result<(), BeneGesseritError> {
        let machinery = FissionMachinery {
            protein_id,
            protein_type: protein_type.clone(),
            activity: 1.0,
            nucleotide_binding: match protein_type {
                FissionProteinType::Dynamin => self.gtp_concentration,
                FissionProteinType::VPS4 => self.atp_concentration,
                _ => 0.0,
            },
            membrane_affinity: match protein_type {
                FissionProteinType::Amphiphysin | FissionProteinType::Endophilin => 2.0,
                _ => 1.0,
            },
            oligomerization_state: 0.0,
            is_active: true,
        };
        
        self.fission_proteins.insert(protein_id, machinery);
        Ok(())
    }
    
    /// Create dynamin assembly
    pub fn create_dynamin_assembly(&mut self, dynamin_proteins: Vec<ProteinId>) -> Result<u64, BeneGesseritError> {
        let assembly_id = self.event_counter;
        self.event_counter += 1;
        
        let assembly = DynaminAssembly {
            assembly_id,
            dynamin_proteins,
            helix_pitch: 13.0, // nm
            constriction_force: 0.0,
            gtp_binding_state: 0.0,
            assembly_completeness: 0.0,
            is_active: false,
        };
        
        self.dynamin_assemblies.insert(assembly_id, assembly);
        Ok(assembly_id)
    }
    
    /// Create ESCRT complex
    pub fn create_escrt_complex(
        &mut self,
        escrt_i: Vec<ProteinId>,
        escrt_ii: Vec<ProteinId>,
        escrt_iii: Vec<ProteinId>,
        vps4: Option<ProteinId>,
    ) -> Result<u64, BeneGesseritError> {
        let complex_id = self.event_counter;
        self.event_counter += 1;
        
        let complex = ESCRTComplex {
            complex_id,
            escrt_i,
            escrt_ii,
            escrt_iii,
            vps4_atpase: vps4,
            assembly_state: 0.0,
            scission_activity: 0.0,
            is_assembled: false,
        };
        
        self.escrt_complexes.insert(complex_id, complex);
        Ok(complex_id)
    }
    
    /// Create clathrin coat
    pub fn create_clathrin_coat(&mut self, clathrin_proteins: Vec<ProteinId>, adaptors: Vec<ProteinId>) -> Result<u64, BeneGesseritError> {
        let coat_id = self.event_counter;
        self.event_counter += 1;
        
        let coat = ClathrinCoat {
            coat_id,
            clathrin_proteins,
            adaptors,
            coat_curvature: 0.0,
            polymerization_state: 0.0,
            is_complete: false,
        };
        
        self.clathrin_coats.insert(coat_id, coat);
        Ok(coat_id)
    }
    
    /// Initiate fission event
    pub fn initiate_fission(
        &mut self,
        fission_type: FissionType,
        mechanism: FissionMechanism,
        parent_membrane: MembraneId,
    ) -> Result<u64, BeneGesseritError> {
        let event_id = self.event_counter;
        self.event_counter += 1;
        
        let event = FissionEvent {
            event_id,
            fission_type,
            mechanism,
            stage: FissionStage::Initiation,
            parent_membrane,
            vesicle_id: None,
            neck_diameter: 50.0, // nm - initial neck diameter
            vesicle_size: 100.0, // nm - target vesicle diameter
            energy_cost: 0.0,
            time: 0.0,
            atp_consumed: 0.0,
            gtp_consumed: 0.0,
        };
        
        self.active_fissions.insert(event_id, event);
        Ok(event_id)
    }
    
    /// Simulate fission dynamics
    pub fn simulate_fission(&mut self, dt: f64) -> Result<Vec<FissionEvent>, BeneGesseritError> {
        let mut completed_events = Vec::new();
        
        // Update protein states
        self.update_protein_states();
        
        // Update assemblies and complexes
        self.update_dynamin_assemblies(dt)?;
        self.update_escrt_complexes(dt)?;
        self.update_clathrin_coats(dt)?;
        
        // Process active fission events
        let mut events_to_remove = Vec::new();
        
        for (event_id, event) in self.active_fissions.iter_mut() {
            if self.advance_fission_stage(event, dt)? {
                // Fission completed
                completed_events.push(event.clone());
                events_to_remove.push(*event_id);
            }
        }
        
        // Remove completed events
        for event_id in events_to_remove {
            if let Some(event) = self.active_fissions.remove(&event_id) {
                self.fission_history.push(event);
            }
        }
        
        Ok(completed_events)
    }
    
    /// Update protein activity states
    fn update_protein_states(&mut self) {
        for protein in self.fission_proteins.values_mut() {
            match protein.protein_type {
                FissionProteinType::Dynamin => {
                    protein.nucleotide_binding = self.gtp_concentration / (self.gtp_concentration + 1e-5);
                    protein.is_active = self.gtp_concentration > 1e-6;
                }
                FissionProteinType::VPS4 => {
                    protein.nucleotide_binding = self.atp_concentration / (self.atp_concentration + 1e-5);
                    protein.is_active = self.atp_concentration > ATP_CONCENTRATION * 0.1;
                }
                _ => {
                    protein.is_active = true;
                }
            }
            
            if protein.is_active {
                protein.activity = protein.nucleotide_binding.max(0.1);
            } else {
                protein.activity = 0.0;
            }
        }
    }
    
    /// Update dynamin assembly states
    fn update_dynamin_assemblies(&mut self, dt: f64) -> Result<(), BeneGesseritError> {
        for assembly in self.dynamin_assemblies.values_mut() {
            if !assembly.is_active {
                // Check if enough GTP-bound dynamin is available
                let gtp_binding = self.gtp_concentration / (self.gtp_concentration + 1e-5);
                assembly.gtp_binding_state = gtp_binding;
                
                if gtp_binding > 0.5 {
                    assembly.assembly_completeness += 2.0 * dt; // Assembly rate
                    assembly.assembly_completeness = assembly.assembly_completeness.min(1.0);
                    
                    if assembly.assembly_completeness > 0.8 {
                        assembly.is_active = true;
                        assembly.constriction_force = 15.0; // pN per dynamin
                    }
                }
            } else {
                // Active constriction
                assembly.helix_pitch -= 1.0 * dt; // Constriction rate
                assembly.helix_pitch = assembly.helix_pitch.max(2.0); // Minimum pitch
                
                // Increase constriction force
                assembly.constriction_force = (13.0 - assembly.helix_pitch) * 2.0;
            }
        }
        
        Ok(())
    }
    
    /// Update ESCRT complex states
    fn update_escrt_complexes(&mut self, dt: f64) -> Result<(), BeneGesseritError> {
        for complex in self.escrt_complexes.values_mut() {
            if !complex.is_assembled {
                // Sequential assembly: ESCRT-I -> ESCRT-II -> ESCRT-III
                let assembly_rate = 1.0 * dt;
                complex.assembly_state += assembly_rate;
                complex.assembly_state = complex.assembly_state.min(1.0);
                
                if complex.assembly_state > 0.9 {
                    complex.is_assembled = true;
                }
            } else {
                // Active scission by ESCRT-III filaments
                if complex.vps4_atpase.is_some() {
                    let atp_factor = self.atp_concentration / (self.atp_concentration + ATP_CONCENTRATION * 0.1);
                    complex.scission_activity = 5.0 * atp_factor; // ATP-dependent activity
                } else {
                    complex.scission_activity = 1.0; // Basal activity
                }
            }
        }
        
        Ok(())
    }
    
    /// Update clathrin coat states
    fn update_clathrin_coats(&mut self, dt: f64) -> Result<(), BeneGesseritError> {
        for coat in self.clathrin_coats.values_mut() {
            if !coat.is_complete {
                // Clathrin polymerization
                let polymerization_rate = 2.0 * dt;
                coat.polymerization_state += polymerization_rate;
                coat.polymerization_state = coat.polymerization_state.min(1.0);
                
                // Increase curvature as coat forms
                coat.coat_curvature = coat.polymerization_state * 0.02; // 1/nm
                
                if coat.polymerization_state > 0.95 {
                    coat.is_complete = true;
                }
            }
        }
        
        Ok(())
    }
    
    /// Advance fission through stages
    fn advance_fission_stage(&mut self, event: &mut FissionEvent, dt: f64) -> Result<bool, BeneGesseritError> {
        let current_stage = event.stage.clone();
        let barrier = *self.stage_barriers.get(&current_stage).unwrap_or(&25.0);
        
        // Calculate transition probability
        let rate = self.calculate_stage_transition_rate(event, &current_stage)?;
        let probability = rate * dt;
        
        if self.random_f64() < probability {
            // Advance to next stage
            event.stage = self.get_next_stage(&current_stage);
            event.energy_cost += barrier;
            
            // Update neck diameter and consume nucleotides
            match event.stage {
                FissionStage::Budding => {
                    event.neck_diameter = 30.0; // nm
                }
                FissionStage::Constriction => {
                    event.neck_diameter = 10.0; // nm
                    if matches!(event.mechanism, FissionMechanism::DynaminMediated) {
                        event.gtp_consumed += 10.0; // GTP molecules per dynamin
                    }
                }
                FissionStage::Scission => {
                    event.neck_diameter = 0.0; // Complete scission
                    event.vesicle_id = Some(VesicleId(self.event_counter));
                    self.event_counter += 1;
                    
                    if matches!(event.mechanism, FissionMechanism::ESCRTMediated) {
                        event.atp_consumed += ATP_HYDROLYSIS_ENERGY * 5.0; // VPS4 activity
                    }
                }
                FissionStage::Release => {
                    // Vesicle released
                }
                _ => {}
            }
            
            // Check if fission is complete
            if matches!(event.stage, FissionStage::Release) {
                return Ok(true);
            }
        }
        
        Ok(false)
    }
    
    /// Calculate stage transition rate
    fn calculate_stage_transition_rate(&self, event: &FissionEvent, stage: &FissionStage) -> Result<f64, BeneGesseritError> {
        let base_rate = match stage {
            FissionStage::Initiation => 5.0,
            FissionStage::Budding => 2.0,
            FissionStage::Constriction => 1.0,
            FissionStage::Scission => 0.5,
            FissionStage::Release => 10.0,
        };
        
        let mut rate = base_rate;
        
        // Mechanism-specific modifiers
        match event.mechanism {
            FissionMechanism::DynaminMediated => {
                let dynamin_effect = self.calculate_dynamin_effect();
                rate *= dynamin_effect;
            }
            FissionMechanism::ESCRTMediated => {
                let escrt_effect = self.calculate_escrt_effect();
                rate *= escrt_effect;
            }
            FissionMechanism::ClathrinMediated => {
                let clathrin_effect = self.calculate_clathrin_effect();
                rate *= clathrin_effect;
            }
            FissionMechanism::CurvatureInduced => {
                // Depends on membrane tension
                let tension = self.membrane_tensions.get(&event.parent_membrane).unwrap_or(&0.0);
                rate *= 1.0 + tension * 0.1;
            }
            FissionMechanism::Spontaneous => {
                rate *= 0.001; // Very slow
            }
        }
        
        Ok(rate)
    }
    
    /// Calculate dynamin assembly contribution
    fn calculate_dynamin_effect(&self) -> f64 {
        let mut total_effect = 1.0;
        
        for assembly in self.dynamin_assemblies.values() {
            if assembly.is_active {
                total_effect += assembly.constriction_force / 15.0; // Normalize by max force
            }
        }
        
        total_effect
    }
    
    /// Calculate ESCRT complex contribution
    fn calculate_escrt_effect(&self) -> f64 {
        let mut total_effect = 1.0;
        
        for complex in self.escrt_complexes.values() {
            if complex.is_assembled {
                total_effect += complex.scission_activity / 5.0; // Normalize by max activity
            }
        }
        
        total_effect
    }
    
    /// Calculate clathrin coat contribution
    fn calculate_clathrin_effect(&self) -> f64 {
        let mut total_effect = 1.0;
        
        for coat in self.clathrin_coats.values() {
            if coat.is_complete {
                total_effect += coat.coat_curvature * 50.0; // Curvature contribution
            }
        }
        
        total_effect
    }
    
    /// Get next fission stage
    fn get_next_stage(&self, current: &FissionStage) -> FissionStage {
        match current {
            FissionStage::Initiation => FissionStage::Budding,
            FissionStage::Budding => FissionStage::Constriction,
            FissionStage::Constriction => FissionStage::Scission,
            FissionStage::Scission => FissionStage::Release,
            FissionStage::Release => FissionStage::Release, // Terminal
        }
    }
    
    /// Get fission statistics
    pub fn get_fission_statistics(&self) -> FissionStatistics {
        let mut mechanism_counts = HashMap::new();
        let mut type_counts = HashMap::new();
        let mut total_atp_consumed = 0.0;
        let mut total_gtp_consumed = 0.0;
        
        for event in &self.fission_history {
            *mechanism_counts.entry(event.mechanism.clone()).or_insert(0) += 1;
            *type_counts.entry(event.fission_type.clone()).or_insert(0) += 1;
            total_atp_consumed += event.atp_consumed;
            total_gtp_consumed += event.gtp_consumed;
        }
        
        FissionStatistics {
            total_fissions: self.fission_history.len(),
            active_fissions: self.active_fissions.len(),
            mechanism_counts,
            type_counts,
            total_atp_consumed,
            total_gtp_consumed,
            average_vesicle_size: self.calculate_average_vesicle_size(),
            dynamin_assemblies: self.dynamin_assemblies.len(),
            escrt_complexes: self.escrt_complexes.len(),
            clathrin_coats: self.clathrin_coats.len(),
        }
    }
    
    /// Calculate average vesicle size
    fn calculate_average_vesicle_size(&self) -> f64 {
        if self.fission_history.is_empty() {
            return 0.0;
        }
        
        let total_size: f64 = self.fission_history.iter()
            .map(|event| event.vesicle_size)
            .sum();
        
        total_size / self.fission_history.len() as f64
    }
    
    /// Set ATP concentration
    pub fn set_atp_concentration(&mut self, concentration: f64) {
        self.atp_concentration = concentration;
    }
    
    /// Set GTP concentration
    pub fn set_gtp_concentration(&mut self, concentration: f64) {
        self.gtp_concentration = concentration;
    }
    
    /// Set membrane tension
    pub fn set_membrane_tension(&mut self, membrane_id: MembraneId, tension: f64) {
        self.membrane_tensions.insert(membrane_id, tension);
    }
    
    // Helper methods
    fn random_f64(&self) -> f64 {
        // Placeholder - would use proper RNG
        0.5
    }
}

/// Fission statistics
#[derive(Debug, Clone)]
pub struct FissionStatistics {
    pub total_fissions: usize,
    pub active_fissions: usize,
    pub mechanism_counts: HashMap<FissionMechanism, usize>,
    pub type_counts: HashMap<FissionType, usize>,
    pub total_atp_consumed: f64,
    pub total_gtp_consumed: f64,
    pub average_vesicle_size: f64,
    pub dynamin_assemblies: usize,
    pub escrt_complexes: usize,
    pub clathrin_coats: usize,
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_fission_creation() {
        let fission = Fission::new();
        assert_eq!(fission.active_fissions.len(), 0);
        assert!(fission.fission_rates.len() > 0);
    }
    
    #[test]
    fn test_dynamin_assembly_creation() {
        let mut fission = Fission::new();
        let dynamins = vec![ProteinId(1), ProteinId(2), ProteinId(3)];
        let result = fission.create_dynamin_assembly(dynamins);
        assert!(result.is_ok());
        assert_eq!(fission.dynamin_assemblies.len(), 1);
    }
    
    #[test]
    fn test_fission_initiation() {
        let mut fission = Fission::new();
        let result = fission.initiate_fission(
            FissionType::Endocytosis,
            FissionMechanism::DynaminMediated,
            MembraneId(1),
        );
        assert!(result.is_ok());
        assert_eq!(fission.active_fissions.len(), 1);
    }
}   