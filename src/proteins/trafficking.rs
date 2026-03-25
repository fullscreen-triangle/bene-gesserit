//! Protein trafficking and membrane insertion dynamics
//!
//! This module simulates the complex processes of protein trafficking
//! to and from the membrane, including:
//! - Protein insertion and removal from membranes
//! - Vesicular trafficking pathways
//! - Endocytosis and exocytosis of membrane proteins
//! - ATP-dependent trafficking processes

use crate::{
    constants::*,
    error::{MembraneError, Result},
    types::*,
};
use std::collections::{HashMap, VecDeque};

/// Trafficking pathway types
#[derive(Debug, Clone, PartialEq)]
pub enum TraffickingPathway {
    /// Direct insertion from cytoplasm
    DirectInsertion,
    /// ER-Golgi-plasma membrane pathway
    SecretoryPathway,
    /// Endocytic recycling
    EndocyticRecycling,
    /// Lysosomal degradation
    LysosomalDegradation,
    /// Autophagy pathway
    Autophagy,
}

/// Vesicle types for protein trafficking
#[derive(Debug, Clone, PartialEq)]
pub enum VesicleType {
    /// COPII vesicles (ER to Golgi)
    COPII,
    /// COPI vesicles (intra-Golgi and Golgi to ER)
    COPI,
    /// Clathrin-coated vesicles
    Clathrin,
    /// Caveolin-coated vesicles
    Caveolin,
    /// Uncoated vesicles
    Uncoated,
}

/// Protein trafficking state
#[derive(Debug, Clone)]
pub struct TraffickingState {
    /// Current cellular location
    pub location: CellularLocation,
    /// Trafficking pathway being followed
    pub pathway: TraffickingPathway,
    /// Associated vesicle (if any)
    pub vesicle: Option<TraffickingVesicle>,
    /// Time since trafficking began
    pub trafficking_time: f64,
    /// Energy consumed in trafficking
    pub energy_consumed: f64,
}

/// Cellular locations for protein trafficking
#[derive(Debug, Clone, PartialEq)]
pub enum CellularLocation {
    /// Endoplasmic reticulum
    ER,
    /// Golgi apparatus
    Golgi,
    /// Plasma membrane
    PlasmaMembrane,
    /// Early endosome
    EarlyEndosome,
    /// Late endosome
    LateEndosome,
    /// Lysosome
    Lysosome,
    /// Cytoplasm
    Cytoplasm,
    /// In transit (vesicle-bound)
    InTransit,
}

/// Trafficking vesicle
#[derive(Debug, Clone)]
pub struct TraffickingVesicle {
    /// Type of vesicle
    pub vesicle_type: VesicleType,
    /// Proteins being transported
    pub cargo_proteins: Vec<String>,
    /// Vesicle diameter (m)
    pub diameter: f64,
    /// Transport velocity (m/s)
    pub velocity: f64,
    /// ATP cost per unit distance
    pub atp_cost_per_distance: f64,
    /// Current position in trafficking pathway
    pub pathway_position: f64, // 0.0 to 1.0
}

/// Membrane insertion machinery
#[derive(Debug, Clone)]
pub struct InsertionMachinery {
    /// Translocon complexes for insertion
    pub translocons: HashMap<String, Translocon>,
    /// Chaperone proteins assisting insertion
    pub chaperones: HashMap<String, Chaperone>,
    /// ATP pools for insertion energy
    pub atp_pools: HashMap<CellularLocation, f64>,
}

/// Translocon complex for protein insertion
#[derive(Debug, Clone)]
pub struct Translocon {
    /// Translocon type (Sec61, Oxa1, etc.)
    pub translocon_type: String,
    /// Current occupancy state
    pub occupied: bool,
    /// Protein being inserted
    pub inserting_protein: Option<String>,
    /// Insertion rate (proteins/second)
    pub insertion_rate: f64,
    /// ATP consumption per insertion
    pub atp_per_insertion: f64,
}

/// Chaperone protein assisting trafficking
#[derive(Debug, Clone)]
pub struct Chaperone {
    /// Chaperone type (BiP, calnexin, etc.)
    pub chaperone_type: String,
    /// Client proteins bound
    pub bound_clients: Vec<String>,
    /// ATP consumption rate
    pub atp_consumption_rate: f64,
    /// Folding assistance efficiency
    pub folding_efficiency: f64,
}

/// Protein trafficking manager
#[derive(Debug, Clone)]
pub struct TraffickingManager {
    /// Proteins currently being trafficked
    pub trafficking_proteins: HashMap<String, TraffickingState>,
    /// Active vesicles
    pub active_vesicles: HashMap<String, TraffickingVesicle>,
    /// Insertion machinery
    pub insertion_machinery: InsertionMachinery,
    /// Trafficking statistics
    pub statistics: TraffickingStatistics,
}

/// Trafficking statistics and metrics
#[derive(Debug, Clone)]
pub struct TraffickingStatistics {
    /// Total proteins trafficked
    pub total_trafficked: u64,
    /// Trafficking success rate
    pub success_rate: f64,
    /// Average trafficking time by pathway
    pub avg_trafficking_times: HashMap<TraffickingPathway, f64>,
    /// Total ATP consumed in trafficking
    pub total_atp_consumed: f64,
    /// Trafficking errors and failures
    pub error_count: u64,
}

impl TraffickingManager {
    /// Create a new trafficking manager
    pub fn new() -> Self {
        Self {
            trafficking_proteins: HashMap::new(),
            active_vesicles: HashMap::new(),
            insertion_machinery: InsertionMachinery::new(),
            statistics: TraffickingStatistics::new(),
        }
    }
    
    /// Initiate protein trafficking
    pub fn initiate_trafficking(&mut self, protein_id: String, 
                               target_location: CellularLocation,
                               atp_availability: f64) -> Result<()> {
        // Determine appropriate trafficking pathway
        let pathway = self.determine_trafficking_pathway(&target_location)?;
        
        // Check ATP availability for trafficking
        let required_atp = self.calculate_trafficking_energy(&pathway, &target_location);
        if atp_availability < required_atp {
            return Err(MembraneError::InsufficientAtp {
                required: required_atp,
                available: atp_availability,
            });
        }
        
        // Create trafficking state
        let trafficking_state = TraffickingState {
            location: CellularLocation::Cytoplasm, // Assume starting point
            pathway,
            vesicle: None,
            trafficking_time: 0.0,
            energy_consumed: 0.0,
        };
        
        self.trafficking_proteins.insert(protein_id, trafficking_state);
        Ok(())
    }
    
    /// Determine appropriate trafficking pathway
    fn determine_trafficking_pathway(&self, target: &CellularLocation) -> Result<TraffickingPathway> {
        match target {
            CellularLocation::PlasmaMembrane => Ok(TraffickingPathway::SecretoryPathway),
            CellularLocation::ER => Ok(TraffickingPathway::DirectInsertion),
            CellularLocation::Golgi => Ok(TraffickingPathway::SecretoryPathway),
            CellularLocation::Lysosome => Ok(TraffickingPathway::LysosomalDegradation),
            CellularLocation::EarlyEndosome => Ok(TraffickingPathway::EndocyticRecycling),
            CellularLocation::LateEndosome => Ok(TraffickingPathway::EndocyticRecycling),
            _ => Err(MembraneError::InvalidParameter(format!("Invalid target location: {:?}", target))),
        }
    }
    
    /// Calculate energy required for trafficking
    fn calculate_trafficking_energy(&self, pathway: &TraffickingPathway, 
                                  target: &CellularLocation) -> f64 {
        // ATP costs based on pathway complexity and distance
        match pathway {
            TraffickingPathway::DirectInsertion => 5.0 * ATP_ENERGY_PER_MOLECULE,
            TraffickingPathway::SecretoryPathway => 20.0 * ATP_ENERGY_PER_MOLECULE,
            TraffickingPathway::EndocyticRecycling => 15.0 * ATP_ENERGY_PER_MOLECULE,
            TraffickingPathway::LysosomalDegradation => 10.0 * ATP_ENERGY_PER_MOLECULE,
            TraffickingPathway::Autophagy => 25.0 * ATP_ENERGY_PER_MOLECULE,
        }
    }
    
    /// Update trafficking processes
    pub fn update_trafficking(&mut self, dt: f64, atp_pool: &mut AtpMeasurement) -> Result<()> {
        let mut completed_trafficking = Vec::new();
        
        for (protein_id, trafficking_state) in &mut self.trafficking_proteins {
            trafficking_state.trafficking_time += dt;
            
            // Update trafficking progress
            match &trafficking_state.pathway {
                TraffickingPathway::SecretoryPathway => {
                    self.update_secretory_trafficking(protein_id, trafficking_state, dt, atp_pool)?;
                },
                TraffickingPathway::EndocyticRecycling => {
                    self.update_endocytic_trafficking(protein_id, trafficking_state, dt, atp_pool)?;
                },
                TraffickingPathway::DirectInsertion => {
                    self.update_direct_insertion(protein_id, trafficking_state, dt, atp_pool)?;
                },
                _ => {
                    // Default simple trafficking model
                    self.update_simple_trafficking(protein_id, trafficking_state, dt, atp_pool)?;
                },
            }
            
            // Check if trafficking is complete
            if self.is_trafficking_complete(trafficking_state) {
                completed_trafficking.push(protein_id.clone());
            }
        }
        
        // Remove completed trafficking
        for protein_id in completed_trafficking {
            if let Some(state) = self.trafficking_proteins.remove(&protein_id) {
                self.statistics.total_trafficked += 1;
                self.statistics.total_atp_consumed += state.energy_consumed;
            }
        }
        
        Ok(())
    }
    
    /// Update secretory pathway trafficking
    fn update_secretory_trafficking(&mut self, protein_id: &str, 
                                  state: &mut TraffickingState, 
                                  dt: f64, atp_pool: &mut AtpMeasurement) -> Result<()> {
        let stage_duration = 300.0; // 5 minutes per stage (simplified)
        let atp_cost_per_stage = 5.0 * ATP_ENERGY_PER_MOLECULE;
        
        // Progress through secretory pathway stages
        match state.location {
            CellularLocation::Cytoplasm => {
                if state.trafficking_time > stage_duration {
                    if atp_pool.molecules as f64 >= 5.0 {
                        state.location = CellularLocation::ER;
                        state.energy_consumed += atp_cost_per_stage;
                        atp_pool.molecules -= 5;
                    }
                }
            },
            CellularLocation::ER => {
                if state.trafficking_time > 2.0 * stage_duration {
                    if atp_pool.molecules as f64 >= 5.0 {
                        state.location = CellularLocation::Golgi;
                        state.energy_consumed += atp_cost_per_stage;
                        atp_pool.molecules -= 5;
                    }
                }
            },
            CellularLocation::Golgi => {
                if state.trafficking_time > 3.0 * stage_duration {
                    if atp_pool.molecules as f64 >= 5.0 {
                        state.location = CellularLocation::PlasmaMembrane;
                        state.energy_consumed += atp_cost_per_stage;
                        atp_pool.molecules -= 5;
                    }
                }
            },
            _ => {},
        }
        
        Ok(())
    }
    
    /// Update endocytic trafficking
    fn update_endocytic_trafficking(&mut self, protein_id: &str,
                                  state: &mut TraffickingState,
                                  dt: f64, atp_pool: &mut AtpMeasurement) -> Result<()> {
        let stage_duration = 120.0; // 2 minutes per stage
        let atp_cost = 3.0 * ATP_ENERGY_PER_MOLECULE;
        
        match state.location {
            CellularLocation::PlasmaMembrane => {
                if state.trafficking_time > stage_duration {
                    if atp_pool.molecules as f64 >= 3.0 {
                        state.location = CellularLocation::EarlyEndosome;
                        state.energy_consumed += atp_cost;
                        atp_pool.molecules -= 3;
                    }
                }
            },
            CellularLocation::EarlyEndosome => {
                if state.trafficking_time > 2.0 * stage_duration {
                    // Decision point: recycling vs degradation
                    let recycling_probability = 0.7; // 70% recycling rate
                    if rand::random::<f64>() < recycling_probability {
                        state.location = CellularLocation::PlasmaMembrane; // Recycle
                    } else {
                        state.location = CellularLocation::LateEndosome; // Degrade
                    }
                    state.energy_consumed += atp_cost;
                    atp_pool.molecules = atp_pool.molecules.saturating_sub(3);
                }
            },
            _ => {},
        }
        
        Ok(())
    }
    
    /// Update direct insertion
    fn update_direct_insertion(&mut self, protein_id: &str,
                             state: &mut TraffickingState,
                             dt: f64, atp_pool: &mut AtpMeasurement) -> Result<()> {
        let insertion_time = 60.0; // 1 minute for direct insertion
        let atp_cost = 2.0 * ATP_ENERGY_PER_MOLECULE;
        
        if state.trafficking_time > insertion_time {
            if atp_pool.molecules as f64 >= 2.0 {
                state.location = CellularLocation::PlasmaMembrane;
                state.energy_consumed += atp_cost;
                atp_pool.molecules -= 2;
            }
        }
        
        Ok(())
    }
    
    /// Update simple trafficking model
    fn update_simple_trafficking(&mut self, protein_id: &str,
                               state: &mut TraffickingState,
                               dt: f64, atp_pool: &mut AtpMeasurement) -> Result<()> {
        let default_time = 180.0; // 3 minutes
        let atp_cost = 4.0 * ATP_ENERGY_PER_MOLECULE;
        
        if state.trafficking_time > default_time {
            if atp_pool.molecules as f64 >= 4.0 {
                state.location = CellularLocation::PlasmaMembrane;
                state.energy_consumed += atp_cost;
                atp_pool.molecules -= 4;
            }
        }
        
        Ok(())
    }
    
    /// Check if trafficking is complete
    fn is_trafficking_complete(&self, state: &TraffickingState) -> bool {
        match state.pathway {
            TraffickingPathway::SecretoryPathway => state.location == CellularLocation::PlasmaMembrane,
            TraffickingPathway::DirectInsertion => state.location == CellularLocation::PlasmaMembrane,
            TraffickingPathway::EndocyticRecycling => {
                state.location == CellularLocation::PlasmaMembrane || 
                state.location == CellularLocation::Lysosome
            },
            TraffickingPathway::LysosomalDegradation => state.location == CellularLocation::Lysosome,
            TraffickingPathway::Autophagy => state.location == CellularLocation::Lysosome,
        }
    }
    
    /// Insert protein into membrane
    pub fn insert_protein_into_membrane(&mut self, protein_id: String,
                                      protein_type: ProteinType,
                                      position: (f64, f64),
                                      atp_pool: &mut AtpMeasurement) -> Result<()> {
        // Check if translocon is available
        let translocon_id = self.find_available_translocon()?;
        
        // Calculate insertion energy requirement
        let insertion_energy = self.calculate_insertion_energy(&protein_type);
        let required_atp = (insertion_energy / ATP_ENERGY_PER_MOLECULE) as u64;
        
        if atp_pool.molecules < required_atp {
            return Err(MembraneError::InsufficientAtp {
                required: insertion_energy,
                available: atp_pool.molecules as f64 * ATP_ENERGY_PER_MOLECULE,
            });
        }
        
        // Consume ATP for insertion
        atp_pool.molecules -= required_atp;
        
        // Mark translocon as occupied
        if let Some(translocon) = self.insertion_machinery.translocons.get_mut(&translocon_id) {
            translocon.occupied = true;
            translocon.inserting_protein = Some(protein_id.clone());
        }
        
        // Update statistics
        self.statistics.total_trafficked += 1;
        self.statistics.total_atp_consumed += insertion_energy;
        
        Ok(())
    }
    
    /// Find available translocon for insertion
    fn find_available_translocon(&self) -> Result<String> {
        for (id, translocon) in &self.insertion_machinery.translocons {
            if !translocon.occupied {
                return Ok(id.clone());
            }
        }
        Err(MembraneError::ResourceUnavailable("No available translocons".to_string()))
    }
    
    /// Calculate energy required for protein insertion
    fn calculate_insertion_energy(&self, protein_type: &ProteinType) -> f64 {
        // Energy based on protein size and complexity
        match protein_type {
            ProteinType::NaKATPase => 15.0 * ATP_ENERGY_PER_MOLECULE, // Large, complex
            ProteinType::CaATPase => 12.0 * ATP_ENERGY_PER_MOLECULE,
            ProteinType::VGSC => 8.0 * ATP_ENERGY_PER_MOLECULE,
            ProteinType::VGKC => 6.0 * ATP_ENERGY_PER_MOLECULE,
            ProteinType::VGCC => 7.0 * ATP_ENERGY_PER_MOLECULE,
            _ => 5.0 * ATP_ENERGY_PER_MOLECULE, // Default
        }
    }
}

impl InsertionMachinery {
    /// Create new insertion machinery
    pub fn new() -> Self {
        let mut translocons = HashMap::new();
        translocons.insert("Sec61".to_string(), Translocon::new("Sec61".to_string()));
        translocons.insert("Oxa1".to_string(), Translocon::new("Oxa1".to_string()));
        
        let mut chaperones = HashMap::new();
        chaperones.insert("BiP".to_string(), Chaperone::new("BiP".to_string()));
        chaperones.insert("Calnexin".to_string(), Chaperone::new("Calnexin".to_string()));
        
        let mut atp_pools = HashMap::new();
        atp_pools.insert(CellularLocation::ER, PHYSIOLOGICAL_ATP);
        atp_pools.insert(CellularLocation::PlasmaMembrane, PHYSIOLOGICAL_ATP);
        
        Self {
            translocons,
            chaperones,
            atp_pools,
        }
    }
}

impl Translocon {
    /// Create new translocon
    pub fn new(translocon_type: String) -> Self {
        Self {
            translocon_type,
            occupied: false,
            inserting_protein: None,
            insertion_rate: 0.1, // proteins/second
            atp_per_insertion: 10.0,
        }
    }
}

impl Chaperone {
    /// Create new chaperone
    pub fn new(chaperone_type: String) -> Self {
        Self {
            chaperone_type,
            bound_clients: Vec::new(),
            atp_consumption_rate: 2.0, // ATP/second
            folding_efficiency: 0.8,
        }
    }
}

impl TraffickingStatistics {
    /// Create new statistics tracker
    pub fn new() -> Self {
        Self {
            total_trafficked: 0,
            success_rate: 1.0,
            avg_trafficking_times: HashMap::new(),
            total_atp_consumed: 0.0,
            error_count: 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_trafficking_manager_creation() {
        let manager = TraffickingManager::new();
        assert!(manager.trafficking_proteins.is_empty());
        assert!(!manager.insertion_machinery.translocons.is_empty());
    }
    
    #[test]
    fn test_trafficking_pathway_determination() {
        let manager = TraffickingManager::new();
        let pathway = manager.determine_trafficking_pathway(&CellularLocation::PlasmaMembrane).unwrap();
        assert_eq!(pathway, TraffickingPathway::SecretoryPathway);
    }
    
    #[test]
    fn test_energy_calculation() {
        let manager = TraffickingManager::new();
        let energy = manager.calculate_trafficking_energy(
            &TraffickingPathway::SecretoryPathway, 
            &CellularLocation::PlasmaMembrane
        );
        assert!(energy > 0.0);
    }
}
