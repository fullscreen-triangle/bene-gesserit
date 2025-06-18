//! Protein-protein and protein-lipid interactions in membrane dynamics
//!
//! This module implements the complex interactions between membrane proteins
//! and their lipid environment, including:
//! - Protein-protein binding and complex formation
//! - Protein-lipid interactions and lipid shell effects
//! - Cooperative binding and allosteric effects
//! - Membrane-mediated protein interactions

use crate::{
    constants::*,
    error::{MembraneError, Result},
    types::*,
};
use std::collections::HashMap;

/// Interaction types between membrane components
#[derive(Debug, Clone, PartialEq)]
pub enum InteractionType {
    /// Direct protein-protein binding
    ProteinProtein {
        binding_energy: f64,  // kJ/mol
        cooperativity: f64,   // Hill coefficient
    },
    /// Protein-lipid interactions
    ProteinLipid {
        lipid_preference: HashMap<LipidType, f64>,
        binding_sites: u8,
    },
    /// Membrane-mediated interactions
    MembraneMediated {
        curvature_sensitivity: f64,
        thickness_sensitivity: f64,
    },
    /// Electrostatic interactions
    Electrostatic {
        charge_product: f64,
        screening_length: f64, // Debye length
    },
}

/// Protein interaction network
#[derive(Debug, Clone)]
pub struct InteractionNetwork {
    /// Map of protein pairs to their interaction types
    pub interactions: HashMap<(String, String), InteractionType>,
    /// Current binding states
    pub binding_states: HashMap<String, BindingState>,
    /// Interaction energies and forces
    pub interaction_energies: HashMap<String, f64>,
    /// Cooperative effects
    pub cooperativity_factors: HashMap<String, f64>,
}

/// Binding state of a protein
#[derive(Debug, Clone)]
pub struct BindingState {
    /// Proteins bound to this protein
    pub bound_proteins: Vec<String>,
    /// Lipids in the first shell
    pub lipid_shell: LipidShell,
    /// Current binding energy
    pub total_binding_energy: f64,
    /// Allosteric state modifications
    pub allosteric_effects: AllostericEffects,
}

/// Lipid shell around a protein
#[derive(Debug, Clone)]
pub struct LipidShell {
    /// Lipid composition in first shell
    pub first_shell: HashMap<LipidType, u32>,
    /// Lipid composition in second shell
    pub second_shell: HashMap<LipidType, u32>,
    /// Total lipid count in shells
    pub total_lipids: u32,
    /// Shell reorganization energy
    pub reorganization_energy: f64,
}

/// Allosteric effects on protein function
#[derive(Debug, Clone)]
pub struct AllostericEffects {
    /// Activity modulation factor (0.0 to 2.0, 1.0 = no effect)
    pub activity_modulation: f64,
    /// ATP binding affinity change
    pub atp_affinity_change: f64,
    /// Conformational preference changes
    pub conformation_bias: HashMap<String, f64>,
}

impl InteractionNetwork {
    /// Create a new interaction network
    pub fn new() -> Self {
        Self {
            interactions: HashMap::new(),
            binding_states: HashMap::new(),
            interaction_energies: HashMap::new(),
            cooperativity_factors: HashMap::new(),
        }
    }
    
    /// Add an interaction between two proteins
    pub fn add_interaction(&mut self, protein1: String, protein2: String, interaction: InteractionType) -> Result<()> {
        let key = if protein1 <= protein2 {
            (protein1.clone(), protein2.clone())
        } else {
            (protein2.clone(), protein1.clone())
        };
        
        self.interactions.insert(key, interaction);
        
        // Initialize binding states if not present
        self.binding_states.entry(protein1).or_insert_with(BindingState::new);
        self.binding_states.entry(protein2).or_insert_with(BindingState::new);
        
        Ok(())
    }
    
    /// Calculate interaction energies for all protein pairs
    pub fn calculate_interaction_energies(&mut self, proteins: &HashMap<String, crate::molecular::MembraneProtein>, 
                                         lipid_composition: &LipidComposition,
                                         temperature: f64) -> Result<()> {
        self.interaction_energies.clear();
        
        for ((protein1_id, protein2_id), interaction_type) in &self.interactions {
            let protein1 = proteins.get(protein1_id)
                .ok_or_else(|| MembraneError::InvalidParameter(format!("Protein {} not found", protein1_id)))?;
            let protein2 = proteins.get(protein2_id)
                .ok_or_else(|| MembraneError::InvalidParameter(format!("Protein {} not found", protein2_id)))?;
            
            let distance = self.calculate_distance(protein1.position, protein2.position);
            let energy = self.calculate_pair_interaction_energy(
                interaction_type, protein1, protein2, distance, temperature
            )?;
            
            let key = format!("{}-{}", protein1_id, protein2_id);
            self.interaction_energies.insert(key, energy);
        }
        
        // Update lipid shell interactions
        self.update_lipid_shell_interactions(proteins, lipid_composition, temperature)?;
        
        Ok(())
    }
    
    /// Calculate distance between two proteins
    fn calculate_distance(&self, pos1: (f64, f64), pos2: (f64, f64)) -> f64 {
        let dx = pos1.0 - pos2.0;
        let dy = pos1.1 - pos2.1;
        (dx * dx + dy * dy).sqrt()
    }
    
    /// Calculate interaction energy between a protein pair
    fn calculate_pair_interaction_energy(&self, 
                                       interaction_type: &InteractionType,
                                       protein1: &crate::molecular::MembraneProtein,
                                       protein2: &crate::molecular::MembraneProtein,
                                       distance: f64,
                                       temperature: f64) -> Result<f64> {
        match interaction_type {
            InteractionType::ProteinProtein { binding_energy, cooperativity } => {
                // Simple exponential decay with distance
                let r0 = 5e-9; // 5 nm characteristic distance
                let energy = binding_energy * (-distance / r0).exp();
                
                // Apply cooperativity effects
                let cooperative_factor = self.calculate_cooperativity_factor(
                    &protein1.protein_type, &protein2.protein_type, *cooperativity
                );
                
                Ok(energy * cooperative_factor)
            },
            
            InteractionType::Electrostatic { charge_product, screening_length } => {
                // Screened Coulomb interaction
                let energy = charge_product * ELEMENTARY_CHARGE.powi(2) / 
                    (4.0 * std::f64::consts::PI * physical::EPSILON_0 * physical::EPSILON_WATER * distance) *
                    (-distance / screening_length).exp();
                
                Ok(energy)
            },
            
            InteractionType::MembraneMediated { curvature_sensitivity, thickness_sensitivity } => {
                // Membrane-mediated interactions depend on local membrane properties
                // For now, return a distance-dependent interaction
                let energy = -curvature_sensitivity * (-distance / 10e-9).exp();
                Ok(energy)
            },
            
            InteractionType::ProteinLipid { .. } => {
                // Protein-lipid interactions handled separately
                Ok(0.0)
            },
        }
    }
    
    /// Calculate cooperativity factor based on protein types
    fn calculate_cooperativity_factor(&self, protein1: &ProteinType, protein2: &ProteinType, cooperativity: f64) -> f64 {
        // Simplified cooperativity model
        // In reality, this would depend on the specific protein complex being formed
        match (protein1, protein2) {
            (ProteinType::NaKATPase, ProteinType::NaKATPase) => cooperativity.powf(2.0),
            (ProteinType::VGSC, ProteinType::VGKC) => 0.8, // Slight negative cooperativity
            _ => 1.0, // No cooperative effects
        }
    }
    
    /// Update lipid shell interactions around proteins
    fn update_lipid_shell_interactions(&mut self, 
                                     proteins: &HashMap<String, crate::molecular::MembraneProtein>,
                                     lipid_composition: &LipidComposition,
                                     temperature: f64) -> Result<()> {
        for (protein_id, protein) in proteins {
            let binding_state = self.binding_states.entry(protein_id.clone())
                .or_insert_with(BindingState::new);
            
            // Calculate lipid shell composition based on protein preferences
            binding_state.lipid_shell = self.calculate_lipid_shell(protein, lipid_composition)?;
            
            // Calculate reorganization energy
            binding_state.lipid_shell.reorganization_energy = 
                self.calculate_lipid_reorganization_energy(&binding_state.lipid_shell, temperature);
        }
        
        Ok(())
    }
    
    /// Calculate lipid shell composition around a protein
    fn calculate_lipid_shell(&self, protein: &crate::molecular::MembraneProtein, 
                           bulk_composition: &LipidComposition) -> Result<LipidShell> {
        let mut first_shell = HashMap::new();
        let mut second_shell = HashMap::new();
        
        // Estimate lipid counts based on protein size and preferences
        let first_shell_radius = 2e-9; // 2 nm
        let second_shell_radius = 4e-9; // 4 nm
        
        let first_shell_area = std::f64::consts::PI * first_shell_radius.powi(2);
        let second_shell_area = std::f64::consts::PI * (second_shell_radius.powi(2) - first_shell_radius.powi(2));
        
        let lipids_per_nm2 = 1.4; // Approximate lipid density
        let first_shell_lipids = (first_shell_area * 1e18 * lipids_per_nm2) as u32;
        let second_shell_lipids = (second_shell_area * 1e18 * lipids_per_nm2) as u32;
        
        // Distribute lipids based on bulk composition (simplified)
        for (lipid_type, fraction) in &bulk_composition.fractions {
            first_shell.insert(*lipid_type, (first_shell_lipids as f64 * fraction) as u32);
            second_shell.insert(*lipid_type, (second_shell_lipids as f64 * fraction) as u32);
        }
        
        Ok(LipidShell {
            first_shell,
            second_shell,
            total_lipids: first_shell_lipids + second_shell_lipids,
            reorganization_energy: 0.0, // Will be calculated separately
        })
    }
    
    /// Calculate energy cost of lipid shell reorganization
    fn calculate_lipid_reorganization_energy(&self, lipid_shell: &LipidShell, temperature: f64) -> f64 {
        // Simplified model based on entropy loss and enthalpic contributions
        let entropy_cost = BOLTZMANN_CONSTANT * temperature * (lipid_shell.total_lipids as f64).ln();
        let enthalpic_contribution = 2e-21; // ~2 kT per lipid
        
        entropy_cost + enthalpic_contribution * (lipid_shell.total_lipids as f64)
    }
    
    /// Update binding states based on current interactions
    pub fn update_binding_states(&mut self, proteins: &HashMap<String, crate::molecular::MembraneProtein>, 
                                temperature: f64, dt: f64) -> Result<()> {
        for (protein_id, binding_state) in &mut self.binding_states {
            // Calculate allosteric effects based on current binding
            binding_state.allosteric_effects = self.calculate_allosteric_effects(
                protein_id, &binding_state.bound_proteins, proteins
            )?;
            
            // Update total binding energy
            binding_state.total_binding_energy = self.calculate_total_binding_energy(protein_id)?;
        }
        
        Ok(())
    }
    
    /// Calculate allosteric effects from protein binding
    fn calculate_allosteric_effects(&self, protein_id: &str, 
                                  bound_proteins: &[String],
                                  proteins: &HashMap<String, crate::molecular::MembraneProtein>) -> Result<AllostericEffects> {
        let mut activity_modulation = 1.0;
        let mut atp_affinity_change = 1.0;
        let mut conformation_bias = HashMap::new();
        
        // Example allosteric effects (simplified)
        for bound_protein_id in bound_proteins {
            if let Some(bound_protein) = proteins.get(bound_protein_id) {
                match bound_protein.protein_type {
                    ProteinType::NaKATPase => {
                        activity_modulation *= 1.2; // Positive cooperativity
                        atp_affinity_change *= 0.8; // Reduced ATP affinity
                    },
                    ProteinType::VGSC => {
                        activity_modulation *= 0.9; // Slight inhibition
                    },
                    _ => {},
                }
            }
        }
        
        Ok(AllostericEffects {
            activity_modulation,
            atp_affinity_change,
            conformation_bias,
        })
    }
    
    /// Calculate total binding energy for a protein
    fn calculate_total_binding_energy(&self, protein_id: &str) -> Result<f64> {
        let mut total_energy = 0.0;
        
        // Sum up all interactions involving this protein
        for (key, energy) in &self.interaction_energies {
            if key.contains(protein_id) {
                total_energy += energy;
            }
        }
        
        // Add lipid shell reorganization energy
        if let Some(binding_state) = self.binding_states.get(protein_id) {
            total_energy += binding_state.lipid_shell.reorganization_energy;
        }
        
        Ok(total_energy)
    }
}

impl BindingState {
    /// Create a new binding state
    pub fn new() -> Self {
        Self {
            bound_proteins: Vec::new(),
            lipid_shell: LipidShell::new(),
            total_binding_energy: 0.0,
            allosteric_effects: AllostericEffects::default(),
        }
    }
}

impl LipidShell {
    /// Create a new empty lipid shell
    pub fn new() -> Self {
        Self {
            first_shell: HashMap::new(),
            second_shell: HashMap::new(),
            total_lipids: 0,
            reorganization_energy: 0.0,
        }
    }
}

impl Default for AllostericEffects {
    fn default() -> Self {
        Self {
            activity_modulation: 1.0,
            atp_affinity_change: 1.0,
            conformation_bias: HashMap::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_interaction_network_creation() {
        let network = InteractionNetwork::new();
        assert!(network.interactions.is_empty());
        assert!(network.binding_states.is_empty());
    }
    
    #[test]
    fn test_add_interaction() {
        let mut network = InteractionNetwork::new();
        let interaction = InteractionType::ProteinProtein {
            binding_energy: -50.0,
            cooperativity: 2.0,
        };
        
        network.add_interaction("protein1".to_string(), "protein2".to_string(), interaction).unwrap();
        assert_eq!(network.interactions.len(), 1);
        assert_eq!(network.binding_states.len(), 2);
    }
    
    #[test]
    fn test_distance_calculation() {
        let network = InteractionNetwork::new();
        let distance = network.calculate_distance((0.0, 0.0), (3.0, 4.0));
        assert!((distance - 5.0).abs() < 1e-10);
    }
}
