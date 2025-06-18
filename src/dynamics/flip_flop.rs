//! Membrane Flip-Flop Dynamics Module
//! 
//! This module simulates lipid flip-flop (transbilayer translocation) dynamics,
//! including spontaneous and protein-assisted flip-flop events.

use crate::types::*;
use crate::constants::*;
use crate::error::BeneGesseritError;
use std::collections::HashMap;

/// Represents the mechanism of flip-flop
#[derive(Debug, Clone, PartialEq)]
pub enum FlipFlopMechanism {
    Spontaneous,
    Flippase,      // ATP-dependent inward translocation
    Floppase,      // ATP-dependent outward translocation
    Scramblase,    // Bidirectional, calcium-activated
}

/// Flip-flop event record
#[derive(Debug, Clone)]
pub struct FlipFlopEvent {
    pub lipid_id: LipidId,
    pub lipid_type: LipidType,
    pub mechanism: FlipFlopMechanism,
    pub from_leaflet: Leaflet,
    pub to_leaflet: Leaflet,
    pub energy_barrier: f64,
    pub time: f64,
    pub atp_consumed: f64,
}

/// Leaflet specification
#[derive(Debug, Clone, PartialEq)]
pub enum Leaflet {
    Outer,
    Inner,
}

/// Flip-flop dynamics simulator
pub struct FlipFlop {
    /// Lipid distribution by leaflet
    pub leaflet_distribution: HashMap<Leaflet, HashMap<LipidType, usize>>,
    /// Flip-flop rates by lipid type and mechanism
    pub flip_rates: HashMap<(LipidType, FlipFlopMechanism), f64>,
    /// Energy barriers for spontaneous flip-flop
    pub energy_barriers: HashMap<LipidType, f64>,
    /// Flippase proteins (ATP-dependent)
    pub flippases: HashMap<ProteinId, FlippaseState>,
    /// Floppase proteins (ATP-dependent)
    pub floppases: HashMap<ProteinId, FloppaseState>,
    /// Scramblase proteins (calcium-activated)
    pub scramblases: HashMap<ProteinId, ScramblaseState>,
    /// Flip-flop event history
    pub event_history: Vec<FlipFlopEvent>,
    /// Current ATP concentration
    pub atp_concentration: f64,
    /// Current calcium concentration
    pub calcium_concentration: f64,
    /// Membrane asymmetry metrics
    pub asymmetry_index: f64,
    /// Temperature
    pub temperature: f64,
}

/// Flippase protein state
#[derive(Debug, Clone)]
pub struct FlippaseState {
    pub protein_id: ProteinId,
    pub activity: f64,
    pub atp_binding_affinity: f64,
    pub substrate_specificity: HashMap<LipidType, f64>,
    pub translocation_rate: f64,
    pub is_active: bool,
}

/// Floppase protein state
#[derive(Debug, Clone)]
pub struct FloppaseState {
    pub protein_id: ProteinId,
    pub activity: f64,
    pub atp_binding_affinity: f64,
    pub substrate_specificity: HashMap<LipidType, f64>,
    pub translocation_rate: f64,
    pub is_active: bool,
}

/// Scramblase protein state
#[derive(Debug, Clone)]
pub struct ScramblaseState {
    pub protein_id: ProteinId,
    pub activity: f64,
    pub calcium_sensitivity: f64,
    pub substrate_specificity: HashMap<LipidType, f64>,
    pub bidirectional_rate: f64,
    pub is_activated: bool,
}

impl Default for FlipFlop {
    fn default() -> Self {
        let mut flip_rates = HashMap::new();
        let mut energy_barriers = HashMap::new();
        
        // Initialize spontaneous flip-flop rates (very slow for most lipids)
        for lipid_type in [LipidType::POPC, LipidType::POPE, LipidType::POPS, LipidType::Sphingomyelin, LipidType::Cholesterol] {
            flip_rates.insert((lipid_type, FlipFlopMechanism::Spontaneous), match lipid_type {
                LipidType::Cholesterol => 1e-3, // Relatively fast
                LipidType::POPC => 1e-6,        // Very slow
                LipidType::POPE => 1e-7,        // Extremely slow (charged)
                LipidType::POPS => 1e-8,        // Extremely slow (charged)
                LipidType::Sphingomyelin => 1e-7, // Very slow
                _ => 1e-6,
            });
            
            // Energy barriers (kJ/mol)
            energy_barriers.insert(lipid_type, match lipid_type {
                LipidType::Cholesterol => 15.0,
                LipidType::POPC => 80.0,
                LipidType::POPE => 120.0,
                LipidType::POPS => 150.0,
                LipidType::Sphingomyelin => 100.0,
                _ => 80.0,
            });
        }
        
        Self {
            leaflet_distribution: HashMap::new(),
            flip_rates,
            energy_barriers,
            flippases: HashMap::new(),
            floppases: HashMap::new(),
            scramblases: HashMap::new(),
            event_history: Vec::new(),
            atp_concentration: ATP_CONCENTRATION,
            calcium_concentration: 1e-7, // 100 nM resting Ca2+
            asymmetry_index: 0.0,
            temperature: TEMPERATURE,
        }
    }
}

impl FlipFlop {
    pub fn new() -> Self {
        Self::default()
    }
    
    /// Add a flippase protein
    pub fn add_flippase(&mut self, protein_id: ProteinId, specificity: HashMap<LipidType, f64>) -> Result<(), BeneGesseritError> {
        let flippase = FlippaseState {
            protein_id,
            activity: 1.0,
            atp_binding_affinity: 1e-6, // 1 μM
            substrate_specificity: specificity,
            translocation_rate: 10.0, // s^-1
            is_active: true,
        };
        
        self.flippases.insert(protein_id, flippase);
        Ok(())
    }
    
    /// Add a floppase protein
    pub fn add_floppase(&mut self, protein_id: ProteinId, specificity: HashMap<LipidType, f64>) -> Result<(), BeneGesseritError> {
        let floppase = FloppaseState {
            protein_id,
            activity: 1.0,
            atp_binding_affinity: 1e-6, // 1 μM
            substrate_specificity: specificity,
            translocation_rate: 8.0, // s^-1
            is_active: true,
        };
        
        self.floppases.insert(protein_id, floppase);
        Ok(())
    }
    
    /// Add a scramblase protein
    pub fn add_scramblase(&mut self, protein_id: ProteinId, specificity: HashMap<LipidType, f64>) -> Result<(), BeneGesseritError> {
        let scramblase = ScramblaseState {
            protein_id,
            activity: 1.0,
            calcium_sensitivity: 1e-6, // 1 μM
            substrate_specificity: specificity,
            bidirectional_rate: 50.0, // s^-1
            is_activated: self.calcium_concentration > 1e-6,
        };
        
        self.scramblases.insert(protein_id, scramblase);
        Ok(())
    }
    
    /// Initialize leaflet distribution
    pub fn initialize_leaflets(&mut self, outer_lipids: HashMap<LipidType, usize>, inner_lipids: HashMap<LipidType, usize>) {
        self.leaflet_distribution.insert(Leaflet::Outer, outer_lipids);
        self.leaflet_distribution.insert(Leaflet::Inner, inner_lipids);
        self.calculate_asymmetry();
    }
    
    /// Calculate membrane asymmetry index
    pub fn calculate_asymmetry(&mut self) {
        let outer = self.leaflet_distribution.get(&Leaflet::Outer).unwrap_or(&HashMap::new());
        let inner = self.leaflet_distribution.get(&Leaflet::Inner).unwrap_or(&HashMap::new());
        
        let mut total_asymmetry = 0.0;
        let mut total_lipids = 0;
        
        for lipid_type in [LipidType::POPC, LipidType::POPE, LipidType::POPS, LipidType::Sphingomyelin, LipidType::Cholesterol] {
            let outer_count = *outer.get(&lipid_type).unwrap_or(&0) as f64;
            let inner_count = *inner.get(&lipid_type).unwrap_or(&0) as f64;
            let total = outer_count + inner_count;
            
            if total > 0.0 {
                let asymmetry = (outer_count - inner_count) / total;
                total_asymmetry += asymmetry.abs();
                total_lipids += 1;
            }
        }
        
        self.asymmetry_index = if total_lipids > 0 {
            total_asymmetry / total_lipids as f64
        } else {
            0.0
        };
    }
    
    /// Simulate flip-flop events over time
    pub fn simulate_flip_flop(&mut self, dt: f64) -> Result<Vec<FlipFlopEvent>, BeneGesseritError> {
        let mut events = Vec::new();
        
        // Update protein states based on current conditions
        self.update_protein_states();
        
        // Simulate spontaneous flip-flop
        events.extend(self.simulate_spontaneous_flip_flop(dt)?);
        
        // Simulate protein-mediated flip-flop
        events.extend(self.simulate_flippase_activity(dt)?);
        events.extend(self.simulate_floppase_activity(dt)?);
        events.extend(self.simulate_scramblase_activity(dt)?);
        
        // Update leaflet distributions
        for event in &events {
            self.apply_flip_flop_event(event)?;
        }
        
        // Update asymmetry
        self.calculate_asymmetry();
        
        // Store events
        self.event_history.extend(events.clone());
        
        Ok(events)
    }
    
    /// Update protein activity states
    fn update_protein_states(&mut self) {
        // Update flippases
        for flippase in self.flippases.values_mut() {
            flippase.is_active = self.atp_concentration > flippase.atp_binding_affinity;
            if flippase.is_active {
                flippase.activity = self.atp_concentration / (self.atp_concentration + flippase.atp_binding_affinity);
            } else {
                flippase.activity = 0.0;
            }
        }
        
        // Update floppases
        for floppase in self.floppases.values_mut() {
            floppase.is_active = self.atp_concentration > floppase.atp_binding_affinity;
            if floppase.is_active {
                floppase.activity = self.atp_concentration / (self.atp_concentration + floppase.atp_binding_affinity);
            } else {
                floppase.activity = 0.0;
            }
        }
        
        // Update scramblases
        for scramblase in self.scramblases.values_mut() {
            scramblase.is_activated = self.calcium_concentration > scramblase.calcium_sensitivity;
            if scramblase.is_activated {
                scramblase.activity = self.calcium_concentration / (self.calcium_concentration + scramblase.calcium_sensitivity);
            } else {
                scramblase.activity = 0.0;
            }
        }
    }
    
    /// Simulate spontaneous flip-flop events
    fn simulate_spontaneous_flip_flop(&self, dt: f64) -> Result<Vec<FlipFlopEvent>, BeneGesseritError> {
        let mut events = Vec::new();
        
        for (leaflet, lipids) in &self.leaflet_distribution {
            for (lipid_type, &count) in lipids {
                if let Some(&rate) = self.flip_rates.get(&(*lipid_type, FlipFlopMechanism::Spontaneous)) {
                    let probability = rate * dt;
                    let expected_events = count as f64 * probability;
                    
                    // Poisson sampling (simplified)
                    let num_events = if expected_events < 0.1 {
                        if self.random_f64() < expected_events { 1 } else { 0 }
                    } else {
                        (expected_events + 0.5) as usize
                    };
                    
                    for _ in 0..num_events {
                        let to_leaflet = match leaflet {
                            Leaflet::Outer => Leaflet::Inner,
                            Leaflet::Inner => Leaflet::Outer,
                        };
                        
                        let event = FlipFlopEvent {
                            lipid_id: self.generate_lipid_id(),
                            lipid_type: *lipid_type,
                            mechanism: FlipFlopMechanism::Spontaneous,
                            from_leaflet: leaflet.clone(),
                            to_leaflet,
                            energy_barrier: *self.energy_barriers.get(lipid_type).unwrap_or(&80.0),
                            time: 0.0, // Current time would be tracked externally
                            atp_consumed: 0.0,
                        };
                        
                        events.push(event);
                    }
                }
            }
        }
        
        Ok(events)
    }
    
    /// Simulate flippase-mediated flip-flop (outer to inner)
    fn simulate_flippase_activity(&self, dt: f64) -> Result<Vec<FlipFlopEvent>, BeneGesseritError> {
        let mut events = Vec::new();
        
        if let Some(outer_lipids) = self.leaflet_distribution.get(&Leaflet::Outer) {
            for flippase in self.flippases.values() {
                if !flippase.is_active {
                    continue;
                }
                
                for (lipid_type, &count) in outer_lipids {
                    if let Some(&specificity) = flippase.substrate_specificity.get(lipid_type) {
                        let rate = flippase.translocation_rate * flippase.activity * specificity;
                        let probability = rate * dt;
                        let expected_events = count as f64 * probability;
                        
                        let num_events = if expected_events < 0.1 {
                            if self.random_f64() < expected_events { 1 } else { 0 }
                        } else {
                            (expected_events + 0.5) as usize
                        };
                        
                        for _ in 0..num_events {
                            let event = FlipFlopEvent {
                                lipid_id: self.generate_lipid_id(),
                                lipid_type: *lipid_type,
                                mechanism: FlipFlopMechanism::Flippase,
                                from_leaflet: Leaflet::Outer,
                                to_leaflet: Leaflet::Inner,
                                energy_barrier: 0.0, // Protein-assisted
                                time: 0.0,
                                atp_consumed: ATP_HYDROLYSIS_ENERGY,
                            };
                            
                            events.push(event);
                        }
                    }
                }
            }
        }
        
        Ok(events)
    }
    
    /// Simulate floppase-mediated flip-flop (inner to outer)
    fn simulate_floppase_activity(&self, dt: f64) -> Result<Vec<FlipFlopEvent>, BeneGesseritError> {
        let mut events = Vec::new();
        
        if let Some(inner_lipids) = self.leaflet_distribution.get(&Leaflet::Inner) {
            for floppase in self.floppases.values() {
                if !floppase.is_active {
                    continue;
                }
                
                for (lipid_type, &count) in inner_lipids {
                    if let Some(&specificity) = floppase.substrate_specificity.get(lipid_type) {
                        let rate = floppase.translocation_rate * floppase.activity * specificity;
                        let probability = rate * dt;
                        let expected_events = count as f64 * probability;
                        
                        let num_events = if expected_events < 0.1 {
                            if self.random_f64() < expected_events { 1 } else { 0 }
                        } else {
                            (expected_events + 0.5) as usize
                        };
                        
                        for _ in 0..num_events {
                            let event = FlipFlopEvent {
                                lipid_id: self.generate_lipid_id(),
                                lipid_type: *lipid_type,
                                mechanism: FlipFlopMechanism::Floppase,
                                from_leaflet: Leaflet::Inner,
                                to_leaflet: Leaflet::Outer,
                                energy_barrier: 0.0, // Protein-assisted
                                time: 0.0,
                                atp_consumed: ATP_HYDROLYSIS_ENERGY,
                            };
                            
                            events.push(event);
                        }
                    }
                }
            }
        }
        
        Ok(events)
    }
    
    /// Simulate scramblase-mediated bidirectional flip-flop
    fn simulate_scramblase_activity(&self, dt: f64) -> Result<Vec<FlipFlopEvent>, BeneGesseritError> {
        let mut events = Vec::new();
        
        for scramblase in self.scramblases.values() {
            if !scramblase.is_activated {
                continue;
            }
            
            // Bidirectional translocation
            for (from_leaflet, lipids) in &self.leaflet_distribution {
                for (lipid_type, &count) in lipids {
                    if let Some(&specificity) = scramblase.substrate_specificity.get(lipid_type) {
                        let rate = scramblase.bidirectional_rate * scramblase.activity * specificity;
                        let probability = rate * dt;
                        let expected_events = count as f64 * probability;
                        
                        let num_events = if expected_events < 0.1 {
                            if self.random_f64() < expected_events { 1 } else { 0 }
                        } else {
                            (expected_events + 0.5) as usize
                        };
                        
                        for _ in 0..num_events {
                            let to_leaflet = match from_leaflet {
                                Leaflet::Outer => Leaflet::Inner,
                                Leaflet::Inner => Leaflet::Outer,
                            };
                            
                            let event = FlipFlopEvent {
                                lipid_id: self.generate_lipid_id(),
                                lipid_type: *lipid_type,
                                mechanism: FlipFlopMechanism::Scramblase,
                                from_leaflet: from_leaflet.clone(),
                                to_leaflet,
                                energy_barrier: 0.0, // Protein-assisted
                                time: 0.0,
                                atp_consumed: 0.0, // Calcium-activated, not ATP-dependent
                            };
                            
                            events.push(event);
                        }
                    }
                }
            }
        }
        
        Ok(events)
    }
    
    /// Apply flip-flop event to update leaflet distributions
    fn apply_flip_flop_event(&mut self, event: &FlipFlopEvent) -> Result<(), BeneGesseritError> {
        // Remove from source leaflet
        if let Some(from_leaflet) = self.leaflet_distribution.get_mut(&event.from_leaflet) {
            if let Some(count) = from_leaflet.get_mut(&event.lipid_type) {
                if *count > 0 {
                    *count -= 1;
                }
            }
        }
        
        // Add to destination leaflet
        self.leaflet_distribution
            .entry(event.to_leaflet.clone())
            .or_insert_with(HashMap::new)
            .entry(event.lipid_type)
            .and_modify(|count| *count += 1)
            .or_insert(1);
        
        Ok(())
    }
    
    /// Get flip-flop statistics
    pub fn get_flip_flop_statistics(&self) -> FlipFlopStatistics {
        let mut mechanism_counts = HashMap::new();
        let mut lipid_type_counts = HashMap::new();
        let mut total_atp_consumed = 0.0;
        
        for event in &self.event_history {
            *mechanism_counts.entry(event.mechanism.clone()).or_insert(0) += 1;
            *lipid_type_counts.entry(event.lipid_type).or_insert(0) += 1;
            total_atp_consumed += event.atp_consumed;
        }
        
        FlipFlopStatistics {
            total_events: self.event_history.len(),
            mechanism_counts,
            lipid_type_counts,
            total_atp_consumed,
            asymmetry_index: self.asymmetry_index,
            leaflet_distribution: self.leaflet_distribution.clone(),
        }
    }
    
    /// Set ATP concentration
    pub fn set_atp_concentration(&mut self, concentration: f64) {
        self.atp_concentration = concentration;
    }
    
    /// Set calcium concentration
    pub fn set_calcium_concentration(&mut self, concentration: f64) {
        self.calcium_concentration = concentration;
    }
    
    // Helper methods
    fn random_f64(&self) -> f64 {
        // Placeholder - would use proper RNG
        0.5
    }
    
    fn generate_lipid_id(&self) -> LipidId {
        // Placeholder - would generate unique ID
        LipidId(42)
    }
}

/// Flip-flop statistics
#[derive(Debug, Clone)]
pub struct FlipFlopStatistics {
    pub total_events: usize,
    pub mechanism_counts: HashMap<FlipFlopMechanism, usize>,
    pub lipid_type_counts: HashMap<LipidType, usize>,
    pub total_atp_consumed: f64,
    pub asymmetry_index: f64,
    pub leaflet_distribution: HashMap<Leaflet, HashMap<LipidType, usize>>,
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_flip_flop_creation() {
        let flip_flop = FlipFlop::new();
        assert_eq!(flip_flop.asymmetry_index, 0.0);
        assert!(flip_flop.flip_rates.len() > 0);
    }
    
    #[test]
    fn test_flippase_addition() {
        let mut flip_flop = FlipFlop::new();
        let mut specificity = HashMap::new();
        specificity.insert(LipidType::POPC, 1.0);
        
        let result = flip_flop.add_flippase(ProteinId(1), specificity);
        assert!(result.is_ok());
        assert_eq!(flip_flop.flippases.len(), 1);
    }
    
    #[test]
    fn test_asymmetry_calculation() {
        let mut flip_flop = FlipFlop::new();
        
        let mut outer = HashMap::new();
        outer.insert(LipidType::POPC, 100);
        
        let mut inner = HashMap::new();
        inner.insert(LipidType::POPC, 50);
        
        flip_flop.initialize_leaflets(outer, inner);
        assert!(flip_flop.asymmetry_index > 0.0);
    }
}           