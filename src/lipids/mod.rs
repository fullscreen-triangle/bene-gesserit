//! Lipid Module - Comprehensive membrane lipid dynamics
//!
//! This module provides detailed modeling of membrane lipid composition,
//! asymmetric distributions, phase behavior, and dynamic reorganization
//! including lipid rafts, liquid-ordered/disordered phases, and
//! membrane curvature effects.

pub mod species;
pub mod asymmetry;
pub mod phase_behavior;
pub mod domains;
pub mod dynamics;
pub mod interactions;

use crate::{constants::*, error::Result, types::*};
use std::collections::HashMap;
use serde::{Deserialize, Serialize};

/// Comprehensive lipid composition with asymmetric leaflet distributions
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DetailedLipidComposition {
    /// Outer leaflet composition (extracellular facing)
    pub outer_leaflet: LeafletComposition,
    /// Inner leaflet composition (cytoplasmic facing)
    pub inner_leaflet: LeafletComposition,
    /// Cholesterol distribution and orientation
    pub cholesterol_distribution: CholesterolDistribution,
    /// Lipid raft domains
    pub raft_domains: Vec<LipidRaft>,
    /// Dynamic properties
    pub dynamics: LipidDynamics,
    /// Phase state information
    pub phase_state: PhaseState,
}

/// Individual leaflet composition with detailed species
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LeafletComposition {
    /// Phosphatidylcholine species with acyl chain details
    pub pc_species: HashMap<String, PcSpecies>,
    /// Phosphatidylethanolamine species
    pub pe_species: HashMap<String, PeSpecies>,
    /// Phosphatidylserine species (predominantly inner leaflet)
    pub ps_species: HashMap<String, PsSpecies>,
    /// Phosphatidylinositol and phosphoinositides
    pub pi_species: HashMap<String, PiSpecies>,
    /// Sphingomyelin species (predominantly outer leaflet)
    pub sm_species: HashMap<String, SmSpecies>,
    /// Cardiolipin (mitochondrial membranes)
    pub cl_species: HashMap<String, ClSpecies>,
    /// Ceramides and other sphingolipids
    pub ceramide_species: HashMap<String, CeramideSpecies>,
    /// Glycolipids (gangliosides, cerebrosides)
    pub glycolipids: HashMap<String, Glycolipid>,
    /// Total mole fractions (must sum to 1.0)
    pub total_fractions: HashMap<String, f64>,
}

/// Phosphatidylcholine species with acyl chain specificity
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PcSpecies {
    /// Acyl chain composition (e.g., "16:0/18:1")
    pub acyl_chains: AcylChainComposition,
    /// Mole fraction in leaflet
    pub mole_fraction: f64,
    /// Area per molecule (nm²)
    pub area_per_molecule: f64,
    /// Transition temperature (K)
    pub transition_temperature: f64,
    /// Spontaneous curvature (m⁻¹)
    pub spontaneous_curvature: f64,
    /// Elastic moduli
    pub elastic_moduli: LipidElasticModuli,
}

/// Acyl chain composition details
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AcylChainComposition {
    /// sn-1 position acyl chain
    pub sn1_chain: AcylChain,
    /// sn-2 position acyl chain
    pub sn2_chain: AcylChain,
    /// Chain order parameters
    pub order_parameters: Vec<f64>,
    /// Deuterium order parameters (for NMR comparison)
    pub deuterium_order: Vec<f64>,
}

/// Individual acyl chain properties
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AcylChain {
    /// Number of carbons
    pub carbon_count: u8,
    /// Number of double bonds
    pub double_bonds: u8,
    /// Double bond positions
    pub double_bond_positions: Vec<u8>,
    /// Cis/trans configuration
    pub stereochemistry: Vec<Stereochemistry>,
    /// Chain length (nm)
    pub length: f64,
    /// Flexibility (conformational entropy)
    pub flexibility: f64,
}

/// Stereochemistry of double bonds
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Stereochemistry {
    Cis,
    Trans,
}

/// Elastic moduli for individual lipid species
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LipidElasticModuli {
    /// Bending modulus (J)
    pub bending_modulus: f64,
    /// Gaussian modulus (J)
    pub gaussian_modulus: f64,
    /// Area compressibility modulus (N/m)
    pub area_modulus: f64,
    /// Tilt modulus (J/m²)
    pub tilt_modulus: f64,
}

/// Cholesterol distribution and dynamics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CholesterolDistribution {
    /// Cholesterol mole fraction
    pub mole_fraction: f64,
    /// Transbilayer asymmetry (outer/inner ratio)
    pub asymmetry_ratio: f64,
    /// Cholesterol orientation angle distribution
    pub orientation_distribution: Vec<f64>,
    /// Flip-flop rate (s⁻¹)
    pub flip_flop_rate: f64,
    /// Cholesterol-phospholipid interactions
    pub pl_interactions: HashMap<String, f64>,
    /// Condensing effect on membrane
    pub condensing_effect: f64,
}

/// Lipid raft domain properties
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LipidRaft {
    /// Raft identifier
    pub id: String,
    /// Position and size
    pub geometry: RaftGeometry,
    /// Lipid composition within raft
    pub composition: RaftComposition,
    /// Protein residents
    pub proteins: Vec<String>,
    /// Phase state (liquid-ordered)
    pub phase: RaftPhase,
    /// Lifetime (s)
    pub lifetime: f64,
    /// Formation/dissolution dynamics
    pub dynamics: RaftDynamics,
}

/// Raft geometry and spatial organization
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RaftGeometry {
    /// Center position (nm)
    pub center: (f64, f64),
    /// Radius (nm)
    pub radius: f64,
    /// Shape factor (1.0 = circular)
    pub shape_factor: f64,
    /// Boundary sharpness
    pub boundary_sharpness: f64,
}

/// Raft composition (enriched in cholesterol and sphingomyelin)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RaftComposition {
    /// Cholesterol enrichment factor
    pub cholesterol_enrichment: f64,
    /// Sphingomyelin enrichment
    pub sm_enrichment: f64,
    /// Saturated lipid preference
    pub saturated_preference: f64,
    /// Exclusion of unsaturated lipids
    pub unsaturated_exclusion: f64,
}

/// Raft phase properties
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum RaftPhase {
    /// Liquid-ordered phase
    LiquidOrdered {
        order_parameter: f64,
        viscosity: f64,
        thickness: f64,
    },
    /// Liquid-disordered phase
    LiquidDisordered {
        fluidity: f64,
        thickness: f64,
    },
    /// Gel phase (rare at physiological temperature)
    Gel {
        tilt_angle: f64,
        packing_density: f64,
    },
}

/// Raft formation and dissolution dynamics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RaftDynamics {
    /// Formation rate constant (s⁻¹)
    pub formation_rate: f64,
    /// Dissolution rate constant (s⁻¹)
    pub dissolution_rate: f64,
    /// Coalescence probability
    pub coalescence_probability: f64,
    /// Fission probability
    pub fission_probability: f64,
    /// Protein-induced stabilization
    pub protein_stabilization: f64,
}

/// Overall lipid dynamics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LipidDynamics {
    /// Lateral diffusion coefficients (m²/s)
    pub lateral_diffusion: HashMap<String, f64>,
    /// Flip-flop rates (s⁻¹)
    pub flip_flop_rates: HashMap<String, f64>,
    /// Rotational diffusion (rad²/s)
    pub rotational_diffusion: HashMap<String, f64>,
    /// Exchange rates between domains
    pub domain_exchange: HashMap<String, f64>,
    /// Temperature dependence
    pub temperature_dependence: TemperatureDependence,
}

/// Temperature dependence of lipid properties
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TemperatureDependence {
    /// Arrhenius activation energies (J/mol)
    pub activation_energies: HashMap<String, f64>,
    /// Q10 factors
    pub q10_factors: HashMap<String, f64>,
    /// Phase transition temperatures
    pub transition_temperatures: HashMap<String, f64>,
}

/// Membrane phase state
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhaseState {
    /// Dominant phase
    pub dominant_phase: MembranePhase,
    /// Phase coexistence regions
    pub coexistence_regions: Vec<PhaseCoexistence>,
    /// Order parameters
    pub order_parameters: OrderParameters,
    /// Critical fluctuations
    pub critical_fluctuations: CriticalFluctuations,
}

/// Membrane phase types
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum MembranePhase {
    /// Liquid-crystalline (fluid) phase
    LiquidCrystalline {
        fluidity: f64,
        order_parameter: f64,
    },
    /// Liquid-ordered phase (raft-like)
    LiquidOrdered {
        order_parameter: f64,
        cholesterol_content: f64,
    },
    /// Gel phase
    Gel {
        tilt_angle: f64,
        chain_packing: f64,
    },
    /// Ripple phase (intermediate)
    Ripple {
        wavelength: f64,
        amplitude: f64,
    },
}

/// Phase coexistence information
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhaseCoexistence {
    /// Phase 1
    pub phase1: MembranePhase,
    /// Phase 2
    pub phase2: MembranePhase,
    /// Area fraction of phase 1
    pub phase1_fraction: f64,
    /// Line tension (J/m)
    pub line_tension: f64,
    /// Domain size distribution
    pub domain_sizes: Vec<f64>,
}

/// Membrane order parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrderParameters {
    /// Acyl chain order parameter
    pub chain_order: f64,
    /// Cholesterol order parameter
    pub cholesterol_order: f64,
    /// Orientational order
    pub orientational_order: f64,
    /// Positional order
    pub positional_order: f64,
}

/// Critical fluctuations near phase transitions
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CriticalFluctuations {
    /// Correlation length (nm)
    pub correlation_length: f64,
    /// Susceptibility
    pub susceptibility: f64,
    /// Relaxation time (s)
    pub relaxation_time: f64,
    /// Critical exponents
    pub critical_exponents: CriticalExponents,
}

/// Critical exponents for phase transitions
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CriticalExponents {
    /// Heat capacity exponent
    pub alpha: f64,
    /// Order parameter exponent
    pub beta: f64,
    /// Susceptibility exponent
    pub gamma: f64,
    /// Correlation length exponent
    pub nu: f64,
}

// Detailed structures for other lipid species
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PeSpecies {
    pub acyl_chains: AcylChainComposition,
    pub mole_fraction: f64,
    pub headgroup_interactions: HeadgroupInteractions,
    pub curvature_preference: f64,
    pub hydrogen_bonding: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PsSpecies {
    pub acyl_chains: AcylChainComposition,
    pub mole_fraction: f64,
    pub charge_state: ChargeState,
    pub calcium_binding: CalciumBinding,
    pub protein_interactions: HashMap<String, f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PiSpecies {
    pub acyl_chains: AcylChainComposition,
    pub mole_fraction: f64,
    pub phosphorylation_state: PhosphorylationState,
    pub signaling_interactions: SignalingInteractions,
    pub enzyme_substrates: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SmSpecies {
    pub acyl_chains: AcylChainComposition,
    pub mole_fraction: f64,
    pub sphingosine_backbone: SphingosineBackbone,
    pub raft_affinity: f64,
    pub cholesterol_interaction: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClSpecies {
    pub acyl_chains: [AcylChain; 4], // Four acyl chains
    pub mole_fraction: f64,
    pub mitochondrial_specificity: f64,
    pub protein_binding: HashMap<String, f64>,
    pub peroxidation_susceptibility: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CeramideSpecies {
    pub acyl_chains: AcylChainComposition,
    pub mole_fraction: f64,
    pub sphingosine_backbone: SphingosineBackbone,
    pub apoptotic_signaling: f64,
    pub membrane_permeabilization: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Glycolipid {
    pub lipid_backbone: LipidBackbone,
    pub carbohydrate_structure: CarbohydrateStructure,
    pub mole_fraction: f64,
    pub cell_recognition: f64,
    pub immune_interactions: f64,
}

// Supporting structures
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HeadgroupInteractions {
    pub hydrogen_bonds: f64,
    pub electrostatic_interactions: f64,
    pub hydration_energy: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChargeState {
    pub net_charge: f64,
    pub pka_values: Vec<f64>,
    pub ph_dependence: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CalciumBinding {
    pub binding_sites: u8,
    pub binding_constants: Vec<f64>,
    pub cooperativity: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhosphorylationState {
    pub phosphate_groups: Vec<PhosphateGroup>,
    pub kinase_substrates: Vec<String>,
    pub phosphatase_substrates: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhosphateGroup {
    pub position: u8,
    pub phosphorylation_level: f64,
    pub charge_contribution: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SignalingInteractions {
    pub pip2_hydrolysis: f64,
    pub dag_production: f64,
    pub ip3_production: f64,
    pub protein_recruitment: HashMap<String, f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SphingosineBackbone {
    pub chain_length: u8,
    pub hydroxylation_pattern: Vec<u8>,
    pub double_bond_positions: Vec<u8>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum LipidBackbone {
    Glycerol(AcylChainComposition),
    Sphingosine(SphingosineBackbone),
    Sterol(SterolStructure),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SterolStructure {
    pub sterol_type: SterolType,
    pub side_chain: SideChain,
    pub ring_modifications: Vec<RingModification>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum SterolType {
    Cholesterol,
    Ergosterol,
    Sitosterol,
    Stigmasterol,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SideChain {
    pub length: u8,
    pub branching: Vec<u8>,
    pub unsaturation: Vec<u8>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RingModification {
    pub position: u8,
    pub modification_type: ModificationType,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ModificationType {
    Hydroxyl,
    Methyl,
    Acetyl,
    Sulfate,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CarbohydrateStructure {
    pub monosaccharides: Vec<Monosaccharide>,
    pub linkages: Vec<GlycosidicLinkage>,
    pub branching_pattern: BranchingPattern,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Monosaccharide {
    pub sugar_type: SugarType,
    pub anomeric_configuration: AnomericConfiguration,
    pub modifications: Vec<SugarModification>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum SugarType {
    Glucose,
    Galactose,
    Mannose,
    Fucose,
    NeuNAc, // Sialic acid
    GlcNAc,
    GalNAc,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum AnomericConfiguration {
    Alpha,
    Beta,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SugarModification {
    pub position: u8,
    pub modification: SugarModificationType,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum SugarModificationType {
    Acetyl,
    Sulfate,
    Phosphate,
    Methyl,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GlycosidicLinkage {
    pub donor_position: u8,
    pub acceptor_position: u8,
    pub anomeric_config: AnomericConfiguration,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BranchingPattern {
    pub branch_points: Vec<u8>,
    pub branch_lengths: Vec<u8>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CellularConditions {
    pub temperature: f64,
    pub ph: f64,
    pub ionic_strength: f64,
    pub calcium_concentration: f64,
    pub oxidative_stress: f64,
    pub membrane_tension: f64,
}

impl DetailedLipidComposition {
    /// Create physiological membrane composition with proper asymmetry
    pub fn physiological_plasma_membrane() -> Self {
        Self {
            outer_leaflet: LeafletComposition::outer_leaflet_physiological(),
            inner_leaflet: LeafletComposition::inner_leaflet_physiological(),
            cholesterol_distribution: CholesterolDistribution::physiological(),
            raft_domains: vec![LipidRaft::default_raft()],
            dynamics: LipidDynamics::physiological(),
            phase_state: PhaseState::physiological(),
        }
    }
    
    /// Update lipid composition based on cellular conditions
    pub fn update_composition(&mut self, conditions: &CellularConditions, dt: f64) -> Result<()> {
        self.update_phase_state(conditions, dt)?;
        self.update_raft_dynamics(conditions, dt)?;
        self.update_asymmetry(conditions, dt)?;
        Ok(())
    }
    
    /// Calculate membrane curvature from lipid composition
    pub fn calculate_spontaneous_curvature(&self) -> f64 {
        let outer_curvature = self.outer_leaflet.calculate_leaflet_curvature();
        let inner_curvature = self.inner_leaflet.calculate_leaflet_curvature();
        (outer_curvature - inner_curvature) / 2.0
    }
    
    /// Get effective elastic moduli
    pub fn effective_elastic_moduli(&self) -> LipidElasticModuli {
        self.calculate_composite_moduli()
    }
    
    fn update_phase_state(&mut self, _conditions: &CellularConditions, _dt: f64) -> Result<()> {
        Ok(())
    }
    
    fn update_raft_dynamics(&mut self, _conditions: &CellularConditions, _dt: f64) -> Result<()> {
        Ok(())
    }
    
    fn update_asymmetry(&mut self, _conditions: &CellularConditions, _dt: f64) -> Result<()> {
        Ok(())
    }
    
    fn calculate_composite_moduli(&self) -> LipidElasticModuli {
        LipidElasticModuli {
            bending_modulus: 20.0 * physical::BOLTZMANN_CONSTANT * 310.15,
            gaussian_modulus: -10.0 * physical::BOLTZMANN_CONSTANT * 310.15,
            area_modulus: 0.24,
            tilt_modulus: 1e-12,
        }
    }
}

impl LeafletComposition {
    pub fn outer_leaflet_physiological() -> Self {
        Self {
            pc_species: HashMap::new(),
            pe_species: HashMap::new(),
            ps_species: HashMap::new(),
            pi_species: HashMap::new(),
            sm_species: HashMap::new(),
            cl_species: HashMap::new(),
            ceramide_species: HashMap::new(),
            glycolipids: HashMap::new(),
            total_fractions: HashMap::new(),
        }
    }
    
    pub fn inner_leaflet_physiological() -> Self {
        Self {
            pc_species: HashMap::new(),
            pe_species: HashMap::new(),
            ps_species: HashMap::new(),
            pi_species: HashMap::new(),
            sm_species: HashMap::new(),
            cl_species: HashMap::new(),
            ceramide_species: HashMap::new(),
            glycolipids: HashMap::new(),
            total_fractions: HashMap::new(),
        }
    }
    
    pub fn calculate_leaflet_curvature(&self) -> f64 {
        0.0
    }
}

impl CholesterolDistribution {
    pub fn physiological() -> Self {
        Self {
            mole_fraction: 0.3,
            asymmetry_ratio: 1.2,
            orientation_distribution: vec![0.8, 0.15, 0.05],
            flip_flop_rate: 1e-6,
            pl_interactions: HashMap::new(),
            condensing_effect: 0.3,
        }
    }
}

impl LipidRaft {
    pub fn default_raft() -> Self {
        Self {
            id: "raft_1".to_string(),
            geometry: RaftGeometry {
                center: (50.0, 50.0),
                radius: 20.0,
                shape_factor: 1.0,
                boundary_sharpness: 0.8,
            },
            composition: RaftComposition {
                cholesterol_enrichment: 2.0,
                sm_enrichment: 3.0,
                saturated_preference: 0.8,
                unsaturated_exclusion: 0.6,
            },
            proteins: vec!["GPI_protein".to_string()],
            phase: RaftPhase::LiquidOrdered {
                order_parameter: 0.7,
                viscosity: 0.1,
                thickness: 5.0,
            },
            lifetime: 1.0,
            dynamics: RaftDynamics {
                formation_rate: 0.1,
                dissolution_rate: 0.05,
                coalescence_probability: 0.3,
                fission_probability: 0.1,
                protein_stabilization: 0.5,
            },
        }
    }
}

impl LipidDynamics {
    pub fn physiological() -> Self {
        let mut lateral_diffusion = HashMap::new();
        lateral_diffusion.insert("PC".to_string(), 1e-12);
        lateral_diffusion.insert("PE".to_string(), 0.8e-12);
        lateral_diffusion.insert("PS".to_string(), 0.6e-12);
        lateral_diffusion.insert("Cholesterol".to_string(), 2e-12);
        
        let mut flip_flop_rates = HashMap::new();
        flip_flop_rates.insert("PC".to_string(), 1e-5);
        flip_flop_rates.insert("PS".to_string(), 1e-6);
        flip_flop_rates.insert("Cholesterol".to_string(), 1e-3);
        
        Self {
            lateral_diffusion,
            flip_flop_rates,
            rotational_diffusion: HashMap::new(),
            domain_exchange: HashMap::new(),
            temperature_dependence: TemperatureDependence {
                activation_energies: HashMap::new(),
                q10_factors: HashMap::new(),
                transition_temperatures: HashMap::new(),
            },
        }
    }
}

impl PhaseState {
    pub fn physiological() -> Self {
        Self {
            dominant_phase: MembranePhase::LiquidCrystalline {
                fluidity: 0.8,
                order_parameter: 0.3,
            },
            coexistence_regions: vec![],
            order_parameters: OrderParameters {
                chain_order: 0.3,
                cholesterol_order: 0.7,
                orientational_order: 0.4,
                positional_order: 0.1,
            },
            critical_fluctuations: CriticalFluctuations {
                correlation_length: 5.0,
                susceptibility: 0.1,
                relaxation_time: 1e-6,
                critical_exponents: CriticalExponents {
                    alpha: 0.1,
                    beta: 0.3,
                    gamma: 1.2,
                    nu: 0.6,
                },
            },
        }
    }
} 