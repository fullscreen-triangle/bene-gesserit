//! Circuit Interface Layer - Membrane → Hierarchical Probabilistic Electric Circuits
//!
//! This module provides the critical interface between biological membrane dynamics
//! and the Nebuchadnezzar hierarchical probabilistic electric circuit system.
//! It converts membrane states to circuit parameters using ATP-based differential
//! equations (dx/dATP instead of dx/dt).

use crate::{
    constants::*,
    error::{MembraneError, Result},
    molecular::MolecularMembrane,
    types::*,
};
use std::collections::HashMap;

/// Circuit interface for Nebuchadnezzar integration
#[derive(Debug, Clone)]
pub struct CircuitInterface {
    /// Membrane system being interfaced
    pub membrane: MolecularMembrane,
    /// Current circuit representation
    pub circuit_state: HierarchicalCircuit,
    /// ATP-based circuit dynamics
    pub atp_dynamics: AtpCircuitDynamics,
    /// Probabilistic circuit parameters
    pub probabilistic_params: ProbabilisticCircuitParams,
    /// Interface configuration
    pub config: CircuitInterfaceConfig,
}

/// Hierarchical probabilistic electric circuit representation
#[derive(Debug, Clone)]
pub struct HierarchicalCircuit {
    /// Circuit hierarchy levels (molecular → mesoscale → cellular)
    pub levels: Vec<CircuitLevel>,
    /// Inter-level connections and coupling
    pub level_connections: Vec<LevelConnection>,
    /// Global circuit properties
    pub global_properties: GlobalCircuitProperties,
    /// ATP-driven circuit evolution
    pub atp_evolution: AtpCircuitEvolution,
}

/// Individual circuit level in hierarchy
#[derive(Debug, Clone)]
pub struct CircuitLevel {
    /// Level identifier (0=molecular, 1=mesoscale, 2=cellular)
    pub level_id: u8,
    /// Circuit nodes at this level
    pub nodes: Vec<CircuitNode>,
    /// Circuit edges/connections
    pub edges: Vec<CircuitEdge>,
    /// Level-specific properties
    pub properties: LevelProperties,
    /// ATP constraints at this level
    pub atp_constraints: AtpConstraints,
}

/// Circuit node with probabilistic properties
#[derive(Debug, Clone)]
pub struct CircuitNode {
    /// Node unique identifier
    pub id: String,
    /// Node position in circuit space
    pub position: (f64, f64, f64),
    /// Electrical properties
    pub electrical: ElectricalProperties,
    /// Probabilistic properties
    pub probabilistic: ProbabilisticProperties,
    /// ATP dependence
    pub atp_dependence: AtpDependence,
    /// Biological source (protein, lipid domain, etc.)
    pub biological_source: BiologicalSource,
}

/// Circuit edge/connection
#[derive(Debug, Clone)]
pub struct CircuitEdge {
    /// Source node ID
    pub source: String,
    /// Target node ID
    pub target: String,
    /// Connection strength/conductance
    pub conductance: f64,
    /// Probabilistic connection weight
    pub probability: f64,
    /// ATP-dependent connection strength
    pub atp_modulation: f64,
    /// Connection type
    pub connection_type: ConnectionType,
}

/// Electrical properties of circuit nodes
#[derive(Debug, Clone)]
pub struct ElectricalProperties {
    /// Voltage (V)
    pub voltage: f64,
    /// Current (A)
    pub current: f64,
    /// Resistance (Ω)
    pub resistance: f64,
    /// Capacitance (F)
    pub capacitance: f64,
    /// Inductance (H) - for dynamic effects
    pub inductance: f64,
}

/// Probabilistic properties for Nebuchadnezzar integration
#[derive(Debug, Clone)]
pub struct ProbabilisticProperties {
    /// Probability distribution parameters
    pub distribution: ProbabilityDistribution,
    /// Uncertainty bounds
    pub uncertainty: UncertaintyBounds,
    /// Correlation with other nodes
    pub correlations: HashMap<String, f64>,
    /// Temporal evolution probability
    pub evolution_probability: f64,
}

/// Probability distribution for circuit parameters
#[derive(Debug, Clone)]
pub enum ProbabilityDistribution {
    /// Gaussian distribution (mean, std_dev)
    Gaussian(f64, f64),
    /// Beta distribution (alpha, beta)
    Beta(f64, f64),
    /// Exponential distribution (lambda)
    Exponential(f64),
    /// Custom distribution with samples
    Custom(Vec<f64>),
}

/// Uncertainty bounds for probabilistic parameters
#[derive(Debug, Clone)]
pub struct UncertaintyBounds {
    /// Lower bound (confidence interval)
    pub lower: f64,
    /// Upper bound (confidence interval)
    pub upper: f64,
    /// Confidence level (0.0 to 1.0)
    pub confidence: f64,
}

/// ATP dependence of circuit parameters
#[derive(Debug, Clone)]
pub struct AtpDependence {
    /// Base value without ATP constraint
    pub base_value: f64,
    /// ATP sensitivity coefficient
    pub atp_sensitivity: f64,
    /// Michaelis-Menten parameters for ATP dependence
    pub km_atp: f64,
    /// Hill coefficient for cooperativity
    pub hill_coefficient: f64,
    /// Current ATP-modulated value
    pub current_value: f64,
}

/// Biological source of circuit element
#[derive(Debug, Clone)]
pub enum BiologicalSource {
    /// Ion channel protein
    IonChannel(ProteinType),
    /// ATP-driven pump
    AtpPump(ProteinType),
    /// Lipid domain/raft
    LipidDomain(LipidComposition),
    /// Membrane patch
    MembranePatch(f64), // area
    /// Protein complex
    ProteinComplex(Vec<ProteinType>),
}

/// Connection types in circuit
#[derive(Debug, Clone)]
pub enum ConnectionType {
    /// Direct electrical connection
    Electrical,
    /// Capacitive coupling
    Capacitive,
    /// ATP-mediated coupling
    AtpMediated,
    /// Probabilistic correlation
    Probabilistic,
    /// Hierarchical level coupling
    Hierarchical,
}

/// ATP-based circuit dynamics (dx/dATP equations)
#[derive(Debug, Clone)]
pub struct AtpCircuitDynamics {
    /// Current ATP concentration affecting circuit
    pub atp_concentration: f64,
    /// ATP consumption rate by circuit elements
    pub atp_consumption_rate: f64,
    /// Circuit parameter derivatives w.r.t. ATP
    pub datp_derivatives: HashMap<String, f64>,
    /// ATP-constrained evolution equations
    pub evolution_equations: Vec<AtpEvolutionEquation>,
}

/// Individual ATP evolution equation: dx/dATP = f(x, ATP, t)
#[derive(Debug, Clone)]
pub struct AtpEvolutionEquation {
    /// Variable name being evolved
    pub variable: String,
    /// Current variable value
    pub current_value: f64,
    /// Derivative function parameters
    pub derivative_params: DerivativeParams,
    /// ATP dependence type
    pub atp_dependence_type: AtpDependenceType,
}

/// Parameters for derivative calculation
#[derive(Debug, Clone)]
pub struct DerivativeParams {
    /// Linear coefficient
    pub linear_coeff: f64,
    /// Nonlinear coefficient
    pub nonlinear_coeff: f64,
    /// Coupling coefficients with other variables
    pub coupling_coeffs: HashMap<String, f64>,
    /// Time constant
    pub time_constant: f64,
}

/// Types of ATP dependence
#[derive(Debug, Clone)]
pub enum AtpDependenceType {
    /// Linear dependence: dx/dATP = k * ATP
    Linear(f64),
    /// Michaelis-Menten: dx/dATP = Vmax * ATP / (Km + ATP)
    MichaelisMenten(f64, f64), // Vmax, Km
    /// Hill equation: dx/dATP = Vmax * ATP^n / (Km^n + ATP^n)
    Hill(f64, f64, f64), // Vmax, Km, n
    /// Exponential: dx/dATP = k * exp(-E/RT) * ATP
    Exponential(f64, f64), // k, E/RT
}

/// Probabilistic circuit parameters for Nebuchadnezzar
#[derive(Debug, Clone)]
pub struct ProbabilisticCircuitParams {
    /// Node probability distributions
    pub node_distributions: HashMap<String, ProbabilityDistribution>,
    /// Edge probability distributions
    pub edge_distributions: HashMap<String, ProbabilityDistribution>,
    /// Correlation matrix between circuit elements
    pub correlation_matrix: Vec<Vec<f64>>,
    /// Temporal evolution probabilities
    pub temporal_evolution: TemporalEvolution,
}

/// Temporal evolution of probabilistic parameters
#[derive(Debug, Clone)]
pub struct TemporalEvolution {
    /// Time-dependent probability functions
    pub time_functions: Vec<TimeFunction>,
    /// Markov chain transition matrices
    pub transition_matrices: Vec<Vec<Vec<f64>>>,
    /// Stochastic differential equation parameters
    pub sde_params: SdeParams,
}

/// Time-dependent probability function
#[derive(Debug, Clone)]
pub struct TimeFunction {
    /// Function type
    pub function_type: TimeFunctionType,
    /// Parameters
    pub parameters: Vec<f64>,
    /// Time range
    pub time_range: (f64, f64),
}

/// Types of time functions
#[derive(Debug, Clone)]
pub enum TimeFunctionType {
    /// Exponential decay
    Exponential,
    /// Oscillatory
    Sinusoidal,
    /// Power law
    PowerLaw,
    /// Custom function
    Custom,
}

/// Stochastic differential equation parameters
#[derive(Debug, Clone)]
pub struct SdeParams {
    /// Drift coefficients
    pub drift_coeffs: Vec<f64>,
    /// Diffusion coefficients
    pub diffusion_coeffs: Vec<f64>,
    /// Noise correlation matrix
    pub noise_correlation: Vec<Vec<f64>>,
}

/// Configuration for circuit interface
#[derive(Debug, Clone)]
pub struct CircuitInterfaceConfig {
    /// Number of hierarchy levels
    pub hierarchy_levels: u8,
    /// ATP constraint strength (0.0 = no constraint, 1.0 = full constraint)
    pub atp_constraint_strength: f64,
    /// Probabilistic update frequency
    pub probabilistic_update_freq: f64,
    /// Circuit resolution (nodes per unit area)
    pub circuit_resolution: f64,
    /// Enable Nebuchadnezzar integration
    pub nebuchadnezzar_integration: bool,
}

/// Level-specific properties
#[derive(Debug, Clone)]
pub struct LevelProperties {
    /// Spatial scale (m)
    pub spatial_scale: f64,
    /// Temporal scale (s)
    pub temporal_scale: f64,
    /// Number of nodes at this level
    pub node_count: usize,
    /// Characteristic impedance
    pub characteristic_impedance: f64,
}

/// ATP constraints at circuit level
#[derive(Debug, Clone)]
pub struct AtpConstraints {
    /// Minimum ATP required for level function
    pub min_atp_required: f64,
    /// ATP allocation to this level
    pub atp_allocation: f64,
    /// ATP efficiency at this level
    pub atp_efficiency: f64,
}

/// Connection between hierarchy levels
#[derive(Debug, Clone)]
pub struct LevelConnection {
    /// Source level
    pub source_level: u8,
    /// Target level
    pub target_level: u8,
    /// Coupling strength
    pub coupling_strength: f64,
    /// ATP-dependent coupling
    pub atp_coupling: f64,
}

/// Global circuit properties
#[derive(Debug, Clone)]
pub struct GlobalCircuitProperties {
    /// Total circuit impedance
    pub total_impedance: f64,
    /// Global ATP consumption
    pub global_atp_consumption: f64,
    /// Circuit stability measure
    pub stability_measure: f64,
    /// Information flow capacity
    pub information_capacity: f64,
}

/// ATP-driven circuit evolution
#[derive(Debug, Clone)]
pub struct AtpCircuitEvolution {
    /// Evolution time constant
    pub time_constant: f64,
    /// ATP depletion rate
    pub atp_depletion_rate: f64,
    /// Circuit adaptation rate
    pub adaptation_rate: f64,
}

impl CircuitInterface {
    /// Create new circuit interface from membrane
    pub fn new(membrane: MolecularMembrane, config: CircuitInterfaceConfig) -> Result<Self> {
        let circuit_state = Self::build_hierarchical_circuit(&membrane, &config)?;
        let atp_dynamics = Self::initialize_atp_dynamics(&membrane)?;
        let probabilistic_params = Self::initialize_probabilistic_params(&circuit_state)?;
        
        Ok(Self {
            membrane,
            circuit_state,
            atp_dynamics,
            probabilistic_params,
            config,
        })
    }
    
    /// Build hierarchical circuit from membrane state
    fn build_hierarchical_circuit(membrane: &MolecularMembrane, config: &CircuitInterfaceConfig) -> Result<HierarchicalCircuit> {
        let mut levels = Vec::new();
        
        // Level 0: Molecular level (individual proteins and lipids)
        let molecular_level = Self::build_molecular_level(membrane)?;
        levels.push(molecular_level);
        
        // Level 1: Mesoscale level (protein clusters, lipid domains)
        if config.hierarchy_levels > 1 {
            let mesoscale_level = Self::build_mesoscale_level(membrane)?;
            levels.push(mesoscale_level);
        }
        
        // Level 2: Cellular level (whole membrane patches)
        if config.hierarchy_levels > 2 {
            let cellular_level = Self::build_cellular_level(membrane)?;
            levels.push(cellular_level);
        }
        
        // Build inter-level connections
        let level_connections = Self::build_level_connections(&levels)?;
        
        let global_properties = GlobalCircuitProperties {
            total_impedance: Self::calculate_total_impedance(&levels),
            global_atp_consumption: membrane.state.atp.consumption_rate,
            stability_measure: Self::calculate_stability(&levels),
            information_capacity: Self::calculate_information_capacity(&levels),
        };
        
        let atp_evolution = AtpCircuitEvolution {
            time_constant: 1e-3, // 1 ms
            atp_depletion_rate: membrane.state.atp.consumption_rate,
            adaptation_rate: 0.1, // 10% per update
        };
        
        Ok(HierarchicalCircuit {
            levels,
            level_connections,
            global_properties,
            atp_evolution,
        })
    }
    
    /// Build molecular level circuit (Level 0)
    fn build_molecular_level(membrane: &MolecularMembrane) -> Result<CircuitLevel> {
        let mut nodes = Vec::new();
        let mut edges = Vec::new();
        
        // Create nodes for each protein
        for (protein_id, protein) in &membrane.proteins {
            let node = CircuitNode {
                id: protein_id.clone(),
                position: (protein.position.0, protein.position.1, 0.0),
                electrical: ElectricalProperties {
                    voltage: membrane.state.voltage,
                    current: protein.calculate_current(&membrane.state.ion_concentrations, membrane.state.temperature),
                    resistance: 1.0 / protein.protein_type.conductance().max(1e-12),
                    capacitance: membrane.state.capacitance / membrane.proteins.len() as f64,
                    inductance: 1e-9, // Small inductance for dynamic effects
                },
                probabilistic: ProbabilisticProperties {
                    distribution: ProbabilityDistribution::Gaussian(protein.activity, 0.1),
                    uncertainty: UncertaintyBounds {
                        lower: protein.activity - 0.2,
                        upper: protein.activity + 0.2,
                        confidence: 0.95,
                    },
                    correlations: HashMap::new(),
                    evolution_probability: protein.activity,
                },
                atp_dependence: AtpDependence {
                    base_value: protein.protein_type.conductance(),
                    atp_sensitivity: if protein.protein_type.requires_atp() { 1.0 } else { 0.0 },
                    km_atp: 1e-3, // 1 mM
                    hill_coefficient: 1.0,
                    current_value: protein.protein_type.conductance() * protein.activity,
                },
                biological_source: BiologicalSource::IonChannel(protein.protein_type.clone()),
            };
            nodes.push(node);
        }
        
        // Create edges between nearby proteins
        for (i, node1) in nodes.iter().enumerate() {
            for (j, node2) in nodes.iter().enumerate() {
                if i != j {
                    let distance = ((node1.position.0 - node2.position.0).powi(2) + 
                                  (node1.position.1 - node2.position.1).powi(2)).sqrt();
                    
                    if distance < 10.0 { // 10 nm coupling distance
                        let conductance = 1e-12 / distance; // Distance-dependent coupling
                        let edge = CircuitEdge {
                            source: node1.id.clone(),
                            target: node2.id.clone(),
                            conductance,
                            probability: (-distance / 5.0).exp(), // Exponential decay
                            atp_modulation: 1.0,
                            connection_type: ConnectionType::Electrical,
                        };
                        edges.push(edge);
                    }
                }
            }
        }
        
        let properties = LevelProperties {
            spatial_scale: 1e-9, // nm scale
            temporal_scale: 1e-6, // μs scale
            node_count: nodes.len(),
            characteristic_impedance: 1e6, // MΩ typical for single proteins
        };
        
        let atp_constraints = AtpConstraints {
            min_atp_required: 1e-4, // 0.1 mM
            atp_allocation: 0.8, // 80% of ATP to molecular level
            atp_efficiency: 0.3, // 30% efficiency
        };
        
        Ok(CircuitLevel {
            level_id: 0,
            nodes,
            edges,
            properties,
            atp_constraints,
        })
    }
    
    /// Build mesoscale level circuit (Level 1)
    fn build_mesoscale_level(membrane: &MolecularMembrane) -> Result<CircuitLevel> {
        // Aggregate molecular nodes into mesoscale domains
        let mut nodes = Vec::new();
        
        // Create lipid domain nodes
        let domain_node = CircuitNode {
            id: "lipid_domain_1".to_string(),
            position: (0.0, 0.0, 0.0),
            electrical: ElectricalProperties {
                voltage: membrane.state.voltage,
                current: membrane.state.current,
                resistance: 1e3, // kΩ for domain
                capacitance: membrane.state.capacitance,
                inductance: 1e-6, // μH
            },
            probabilistic: ProbabilisticProperties {
                distribution: ProbabilityDistribution::Beta(2.0, 2.0),
                uncertainty: UncertaintyBounds {
                    lower: 0.1,
                    upper: 0.9,
                    confidence: 0.90,
                },
                correlations: HashMap::new(),
                evolution_probability: 0.5,
            },
            atp_dependence: AtpDependence {
                base_value: membrane.state.capacitance,
                atp_sensitivity: 0.5,
                km_atp: 2e-3, // 2 mM
                hill_coefficient: 2.0,
                current_value: membrane.state.capacitance,
            },
            biological_source: BiologicalSource::LipidDomain(membrane.lipid_composition.clone()),
        };
        nodes.push(domain_node);
        
        let properties = LevelProperties {
            spatial_scale: 1e-6, // μm scale
            temporal_scale: 1e-3, // ms scale
            node_count: nodes.len(),
            characteristic_impedance: 1e3, // kΩ for domains
        };
        
        let atp_constraints = AtpConstraints {
            min_atp_required: 1e-3, // 1 mM
            atp_allocation: 0.15, // 15% of ATP to mesoscale
            atp_efficiency: 0.5, // 50% efficiency
        };
        
        Ok(CircuitLevel {
            level_id: 1,
            nodes,
            edges: Vec::new(), // Simplified for now
            properties,
            atp_constraints,
        })
    }
    
    /// Build cellular level circuit (Level 2)
    fn build_cellular_level(membrane: &MolecularMembrane) -> Result<CircuitLevel> {
        let mut nodes = Vec::new();
        
        // Single node representing entire membrane patch
        let patch_node = CircuitNode {
            id: "membrane_patch".to_string(),
            position: (0.0, 0.0, 0.0),
            electrical: ElectricalProperties {
                voltage: membrane.state.voltage,
                current: membrane.state.current,
                resistance: 1.0 / membrane.calculate_total_conductance(),
                capacitance: membrane.state.capacitance,
                inductance: 1e-3, // mH
            },
            probabilistic: ProbabilisticProperties {
                distribution: ProbabilityDistribution::Gaussian(membrane.state.voltage, 0.01),
                uncertainty: UncertaintyBounds {
                    lower: membrane.state.voltage - 0.02,
                    upper: membrane.state.voltage + 0.02,
                    confidence: 0.99,
                },
                correlations: HashMap::new(),
                evolution_probability: 0.8,
            },
            atp_dependence: AtpDependence {
                base_value: membrane.state.voltage,
                atp_sensitivity: 0.1,
                km_atp: 5e-3, // 5 mM
                hill_coefficient: 1.0,
                current_value: membrane.state.voltage,
            },
            biological_source: BiologicalSource::MembranePatch(membrane.config.area),
        };
        nodes.push(patch_node);
        
        let properties = LevelProperties {
            spatial_scale: 1e-3, // mm scale
            temporal_scale: 1e-1, // 100 ms scale
            node_count: nodes.len(),
            characteristic_impedance: 1.0, // Ω for whole patch
        };
        
        let atp_constraints = AtpConstraints {
            min_atp_required: 5e-3, // 5 mM
            atp_allocation: 0.05, // 5% of ATP to cellular level
            atp_efficiency: 0.8, // 80% efficiency
        };
        
        Ok(CircuitLevel {
            level_id: 2,
            nodes,
            edges: Vec::new(),
            properties,
            atp_constraints,
        })
    }
    
    /// Build connections between hierarchy levels
    fn build_level_connections(levels: &[CircuitLevel]) -> Result<Vec<LevelConnection>> {
        let mut connections = Vec::new();
        
        for i in 0..levels.len() - 1 {
            let connection = LevelConnection {
                source_level: levels[i].level_id,
                target_level: levels[i + 1].level_id,
                coupling_strength: 0.1, // 10% coupling
                atp_coupling: 0.5, // 50% ATP-dependent
            };
            connections.push(connection);
        }
        
        Ok(connections)
    }
    
    /// Initialize ATP dynamics for circuit
    fn initialize_atp_dynamics(membrane: &MolecularMembrane) -> Result<AtpCircuitDynamics> {
        let mut evolution_equations = Vec::new();
        
        // Voltage evolution equation: dV/dATP
        let voltage_eq = AtpEvolutionEquation {
            variable: "voltage".to_string(),
            current_value: membrane.state.voltage,
            derivative_params: DerivativeParams {
                linear_coeff: -0.01, // Voltage decreases with ATP depletion
                nonlinear_coeff: 0.0,
                coupling_coeffs: HashMap::new(),
                time_constant: 1e-3,
            },
            atp_dependence_type: AtpDependenceType::MichaelisMenten(0.1, 1e-3),
        };
        evolution_equations.push(voltage_eq);
        
        // Current evolution equation: dI/dATP
        let current_eq = AtpEvolutionEquation {
            variable: "current".to_string(),
            current_value: membrane.state.current,
            derivative_params: DerivativeParams {
                linear_coeff: 1e-9, // Current increases with ATP availability
                nonlinear_coeff: 0.0,
                coupling_coeffs: HashMap::new(),
                time_constant: 1e-6,
            },
            atp_dependence_type: AtpDependenceType::Linear(1e-9),
        };
        evolution_equations.push(current_eq);
        
        Ok(AtpCircuitDynamics {
            atp_concentration: membrane.state.atp.concentration,
            atp_consumption_rate: membrane.state.atp.consumption_rate,
            datp_derivatives: HashMap::new(),
            evolution_equations,
        })
    }
    
    /// Initialize probabilistic parameters
    fn initialize_probabilistic_params(circuit: &HierarchicalCircuit) -> Result<ProbabilisticCircuitParams> {
        let mut node_distributions = HashMap::new();
        let mut edge_distributions = HashMap::new();
        
        for level in &circuit.levels {
            for node in &level.nodes {
                node_distributions.insert(node.id.clone(), node.probabilistic.distribution.clone());
            }
            for edge in &level.edges {
                let edge_id = format!("{}_{}", edge.source, edge.target);
                edge_distributions.insert(edge_id, ProbabilityDistribution::Exponential(edge.probability));
            }
        }
        
        // Simple correlation matrix (identity for now)
        let n_nodes = node_distributions.len();
        let mut correlation_matrix = vec![vec![0.0; n_nodes]; n_nodes];
        for i in 0..n_nodes {
            correlation_matrix[i][i] = 1.0;
        }
        
        let temporal_evolution = TemporalEvolution {
            time_functions: vec![TimeFunction {
                function_type: TimeFunctionType::Exponential,
                parameters: vec![1.0, -0.1], // exp(-0.1*t)
                time_range: (0.0, 10.0),
            }],
            transition_matrices: vec![correlation_matrix.clone()],
            sde_params: SdeParams {
                drift_coeffs: vec![0.0; n_nodes],
                diffusion_coeffs: vec![0.1; n_nodes],
                noise_correlation: correlation_matrix.clone(),
            },
        };
        
        Ok(ProbabilisticCircuitParams {
            node_distributions,
            edge_distributions,
            correlation_matrix,
            temporal_evolution,
        })
    }
    
    /// Update circuit interface with new membrane state
    pub fn update(&mut self, dt: f64) -> Result<()> {
        // Update membrane first
        self.membrane.step(dt)?;
        
        // Update ATP dynamics using dx/dATP equations
        self.update_atp_dynamics(dt)?;
        
        // Update circuit parameters based on new membrane state
        self.update_circuit_parameters()?;
        
        // Update probabilistic parameters
        self.update_probabilistic_parameters(dt)?;
        
        Ok(())
    }
    
    /// Update ATP-based circuit dynamics
    fn update_atp_dynamics(&mut self, dt: f64) -> Result<()> {
        let current_atp = self.membrane.state.atp.concentration;
        let atp_change = current_atp - self.atp_dynamics.atp_concentration;
        
        // Update each evolution equation
        for eq in &mut self.atp_dynamics.evolution_equations {
            let derivative = self.calculate_atp_derivative(eq, current_atp)?;
            
            // dx/dATP * dATP gives dx
            let variable_change = derivative * atp_change;
            eq.current_value += variable_change;
            
            // Store derivative for analysis
            self.atp_dynamics.datp_derivatives.insert(eq.variable.clone(), derivative);
        }
        
        self.atp_dynamics.atp_concentration = current_atp;
        self.atp_dynamics.atp_consumption_rate = self.membrane.state.atp.consumption_rate;
        
        Ok(())
    }
    
    /// Calculate dx/dATP derivative for evolution equation
    fn calculate_atp_derivative(&self, eq: &AtpEvolutionEquation, atp_conc: f64) -> Result<f64> {
        match &eq.atp_dependence_type {
            AtpDependenceType::Linear(k) => Ok(*k),
            AtpDependenceType::MichaelisMenten(vmax, km) => {
                Ok(vmax * km / (km + atp_conc).powi(2))
            }
            AtpDependenceType::Hill(vmax, km, n) => {
                let atp_n = atp_conc.powf(*n);
                let km_n = km.powf(*n);
                Ok(vmax * n * atp_n * km_n / (km_n + atp_n).powi(2) / atp_conc)
            }
            AtpDependenceType::Exponential(k, e_rt) => {
                Ok(k * (-e_rt).exp())
            }
        }
    }
    
    /// Update circuit parameters from membrane state
    fn update_circuit_parameters(&mut self) -> Result<()> {
        // Update molecular level (Level 0)
        if let Some(molecular_level) = self.circuit_state.levels.get_mut(0) {
            for (node_id, protein) in &self.membrane.proteins {
                if let Some(node) = molecular_level.nodes.iter_mut().find(|n| n.id == *node_id) {
                    // Update electrical properties
                    node.electrical.current = protein.calculate_current(
                        &self.membrane.state.ion_concentrations,
                        self.membrane.state.temperature
                    );
                    node.electrical.resistance = 1.0 / protein.protein_type.conductance().max(1e-12);
                    
                    // Update ATP dependence
                    node.atp_dependence.current_value = 
                        node.atp_dependence.base_value * protein.activity;
                    
                    // Update probabilistic properties
                    node.probabilistic.distribution = ProbabilityDistribution::Gaussian(
                        protein.activity,
                        0.1 * (1.0 - protein.activity) // Higher uncertainty when inactive
                    );
                }
            }
        }
        
        // Update global properties
        self.circuit_state.global_properties.global_atp_consumption = 
            self.membrane.state.atp.consumption_rate;
        self.circuit_state.global_properties.total_impedance = 
            Self::calculate_total_impedance(&self.circuit_state.levels);
        
        Ok(())
    }
    
    /// Update probabilistic parameters
    fn update_probabilistic_parameters(&mut self, dt: f64) -> Result<()> {
        // Update temporal evolution
        for time_func in &mut self.probabilistic_params.temporal_evolution.time_functions {
            match time_func.function_type {
                TimeFunctionType::Exponential => {
                    // Update exponential decay parameters based on ATP
                    let atp_factor = self.atp_dynamics.atp_concentration / PHYSIOLOGICAL_ATP;
                    time_func.parameters[1] = -0.1 * atp_factor; // Decay rate proportional to ATP
                }
                _ => {} // Other function types handled similarly
            }
        }
        
        // Update SDE parameters
        let atp_noise_factor = (self.atp_dynamics.atp_concentration / PHYSIOLOGICAL_ATP).sqrt();
        for coeff in &mut self.probabilistic_params.temporal_evolution.sde_params.diffusion_coeffs {
            *coeff = 0.1 * atp_noise_factor; // Noise scales with ATP availability
        }
        
        Ok(())
    }
    
    /// Get circuit parameters for Nebuchadnezzar integration
    pub fn get_nebuchadnezzar_parameters(&self) -> NebuchadnezzarParams {
        NebuchadnezzarParams {
            hierarchy_levels: self.circuit_state.levels.len() as u8,
            circuit_nodes: self.extract_circuit_nodes(),
            circuit_edges: self.extract_circuit_edges(),
            atp_dynamics: self.atp_dynamics.clone(),
            probabilistic_params: self.probabilistic_params.clone(),
            global_impedance: self.circuit_state.global_properties.total_impedance,
            stability_measure: self.circuit_state.global_properties.stability_measure,
        }
    }
    
    /// Extract circuit nodes for external use
    fn extract_circuit_nodes(&self) -> Vec<NebuchadnezzarNode> {
        let mut nodes = Vec::new();
        
        for level in &self.circuit_state.levels {
            for node in &level.nodes {
                let neb_node = NebuchadnezzarNode {
                    id: node.id.clone(),
                    level: level.level_id,
                    position: node.position,
                    voltage: node.electrical.voltage,
                    current: node.electrical.current,
                    impedance: complex_impedance(
                        node.electrical.resistance,
                        node.electrical.capacitance,
                        node.electrical.inductance,
                        1000.0 // 1 kHz frequency
                    ),
                    atp_dependence: node.atp_dependence.atp_sensitivity,
                    probability_distribution: node.probabilistic.distribution.clone(),
                };
                nodes.push(neb_node);
            }
        }
        
        nodes
    }
    
    /// Extract circuit edges for external use
    fn extract_circuit_edges(&self) -> Vec<NebuchadnezzarEdge> {
        let mut edges = Vec::new();
        
        for level in &self.circuit_state.levels {
            for edge in &level.edges {
                let neb_edge = NebuchadnezzarEdge {
                    source: edge.source.clone(),
                    target: edge.target.clone(),
                    conductance: edge.conductance,
                    probability: edge.probability,
                    atp_modulation: edge.atp_modulation,
                    connection_type: edge.connection_type.clone(),
                };
                edges.push(neb_edge);
            }
        }
        
        edges
    }
    
    /// Calculate total circuit impedance
    fn calculate_total_impedance(levels: &[CircuitLevel]) -> f64 {
        let mut total_conductance = 0.0;
        
        for level in levels {
            for node in &level.nodes {
                total_conductance += 1.0 / node.electrical.resistance;
            }
        }
        
        1.0 / total_conductance.max(1e-12)
    }
    
    /// Calculate circuit stability measure
    fn calculate_stability(levels: &[CircuitLevel]) -> f64 {
        // Simplified stability measure based on node activity variance
        let mut activity_sum = 0.0;
        let mut activity_sq_sum = 0.0;
        let mut count = 0;
        
        for level in levels {
            for node in &level.nodes {
                let activity = node.atp_dependence.current_value / node.atp_dependence.base_value.max(1e-12);
                activity_sum += activity;
                activity_sq_sum += activity * activity;
                count += 1;
            }
        }
        
        if count > 1 {
            let mean = activity_sum / count as f64;
            let variance = activity_sq_sum / count as f64 - mean * mean;
            1.0 / (1.0 + variance) // Higher stability = lower variance
        } else {
            1.0
        }
    }
    
    /// Calculate information capacity
    fn calculate_information_capacity(levels: &[CircuitLevel]) -> f64 {
        // Shannon capacity based on signal-to-noise ratio
        let mut total_capacity = 0.0;
        
        for level in levels {
            let bandwidth = 1.0 / level.properties.temporal_scale; // Hz
            let signal_power = level.nodes.iter()
                .map(|n| n.electrical.voltage.powi(2) / n.electrical.resistance)
                .sum::<f64>();
            let noise_power = physical::BOLTZMANN_CONSTANT * PHYSIOLOGICAL_TEMPERATURE * bandwidth;
            
            if noise_power > 0.0 {
                let snr = signal_power / noise_power;
                total_capacity += bandwidth * (1.0 + snr).log2();
            }
        }
        
        total_capacity
    }
}

/// Parameters for Nebuchadnezzar system integration
#[derive(Debug, Clone)]
pub struct NebuchadnezzarParams {
    /// Number of hierarchy levels
    pub hierarchy_levels: u8,
    /// Circuit nodes
    pub circuit_nodes: Vec<NebuchadnezzarNode>,
    /// Circuit edges
    pub circuit_edges: Vec<NebuchadnezzarEdge>,
    /// ATP dynamics
    pub atp_dynamics: AtpCircuitDynamics,
    /// Probabilistic parameters
    pub probabilistic_params: ProbabilisticCircuitParams,
    /// Global impedance
    pub global_impedance: f64,
    /// Stability measure
    pub stability_measure: f64,
}

/// Circuit node for Nebuchadnezzar integration
#[derive(Debug, Clone)]
pub struct NebuchadnezzarNode {
    /// Node ID
    pub id: String,
    /// Hierarchy level
    pub level: u8,
    /// 3D position
    pub position: (f64, f64, f64),
    /// Voltage
    pub voltage: f64,
    /// Current
    pub current: f64,
    /// Complex impedance
    pub impedance: (f64, f64), // (real, imaginary)
    /// ATP dependence strength
    pub atp_dependence: f64,
    /// Probability distribution
    pub probability_distribution: ProbabilityDistribution,
}

/// Circuit edge for Nebuchadnezzar integration
#[derive(Debug, Clone)]
pub struct NebuchadnezzarEdge {
    /// Source node ID
    pub source: String,
    /// Target node ID
    pub target: String,
    /// Conductance
    pub conductance: f64,
    /// Connection probability
    pub probability: f64,
    /// ATP modulation strength
    pub atp_modulation: f64,
    /// Connection type
    pub connection_type: ConnectionType,
}

impl Default for CircuitInterfaceConfig {
    fn default() -> Self {
        Self {
            hierarchy_levels: 3,
            atp_constraint_strength: 1.0,
            probabilistic_update_freq: 100.0, // Hz
            circuit_resolution: 1e6, // nodes/m²
            nebuchadnezzar_integration: true,
        }
    }
}

/// Calculate complex impedance Z = R + j(ωL - 1/ωC)
fn complex_impedance(r: f64, c: f64, l: f64, frequency: f64) -> (f64, f64) {
    let omega = 2.0 * std::f64::consts::PI * frequency;
    let reactance = omega * l - 1.0 / (omega * c);
    (r, reactance)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecular::MolecularMembrane;
    
    #[test]
    fn test_circuit_interface_creation() {
        let config = MembraneConfig::default();
        let membrane = MolecularMembrane::new(config).unwrap();
        let interface_config = CircuitInterfaceConfig::default();
        
        let interface = CircuitInterface::new(membrane, interface_config).unwrap();
        assert_eq!(interface.circuit_state.levels.len(), 3);
    }
    
    #[test]
    fn test_atp_dynamics_update() {
        let config = MembraneConfig::default();
        let membrane = MolecularMembrane::new(config).unwrap();
        let interface_config = CircuitInterfaceConfig::default();
        
        let mut interface = CircuitInterface::new(membrane, interface_config).unwrap();
        let initial_atp = interface.atp_dynamics.atp_concentration;
        
        interface.update(1e-6).unwrap(); // 1 μs step
        
        // ATP should change due to consumption
        assert_ne!(interface.atp_dynamics.atp_concentration, initial_atp);
    }
    
    #[test]
    fn test_nebuchadnezzar_parameters() {
        let config = MembraneConfig::default();
        let membrane = MolecularMembrane::new(config).unwrap();
        let interface_config = CircuitInterfaceConfig::default();
        
        let interface = CircuitInterface::new(membrane, interface_config).unwrap();
        let params = interface.get_nebuchadnezzar_parameters();
        
        assert_eq!(params.hierarchy_levels, 3);
        assert!(params.global_impedance > 0.0);
    }
    
    #[test]
    fn test_complex_impedance() {
        let z = complex_impedance(100.0, 1e-6, 1e-3, 1000.0);
        assert!(z.0 > 0.0); // Resistance
        assert!(z.1.abs() > 0.0); // Reactance
    }
} 