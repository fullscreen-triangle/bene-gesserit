# Circuit Interface Layer - Nebuchadnezzar Integration

## Core Philosophy

**"Membranes define the circuit topology; ATP consumption drives the dynamics"**

The **Circuit Interface Layer** serves as the bridge between biophysical membrane dynamics and Nebuchadnezzar's hierarchical probabilistic electric circuit system, operating under **metacognitive orchestrator supervision**. This layer translates membrane properties into dynamic circuit parameters, updates circuit topology based on membrane changes, and enables bidirectional coupling between ATP consumption and membrane dynamics using `dx/dATP` equations.

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────────────────┐
│                    CIRCUIT INTERFACE LAYER                             │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  ┌─────────────────────────────────────────────────────────────────┐   │
│  │               MEMBRANE → CIRCUIT MAPPING                       │   │
│  │                                                                 │   │
│  │  ┌──────────────────┐    ┌──────────────────┐    ┌───────────┐ │   │
│  │  │ Lipid Bilayer    │───►│ RC Circuits      │───►│ Circuit   │ │   │
│  │  │ Properties       │    │ (Dynamic)        │    │ Topology  │ │   │
│  │  └──────────────────┘    └──────────────────┘    └───────────┘ │   │
│  │                                                                 │   │
│  │  ┌──────────────────┐    ┌──────────────────┐    ┌───────────┐ │   │
│  │  │ Membrane         │───►│ Active Circuit   │───►│ ATP       │ │   │
│  │  │ Proteins         │    │ Elements         │    │ Coupling  │ │   │
│  │  └──────────────────┘    └──────────────────┘    └───────────┘ │   │
│  └─────────────────────────────────────────────────────────────────┘   │
│                                   ↕                                     │
│  ┌─────────────────────────────────────────────────────────────────┐   │
│  │               NEBUCHADNEZZAR INTERFACE                         │   │
│  │                                                                 │   │
│  │  ┌──────────────────┐    ┌──────────────────┐    ┌───────────┐ │   │
│  │  │ dx/dATP          │◄──►│ Circuit State    │◄──►│ ATP Pool  │ │   │
│  │  │ Equations        │    │ Updates          │    │ Dynamics  │ │   │
│  │  └──────────────────┘    └──────────────────┘    └───────────┘ │   │
│  │                                                                 │   │
│  │  ┌──────────────────┐    ┌──────────────────┐    ┌───────────┐ │   │
│  │  │ Probabilistic    │◄──►│ Hierarchical     │◄──►│ Dynamic   │ │   │
│  │  │ Circuit Priors   │    │ Abstraction      │    │ Expansion │ │   │
│  │  └──────────────────┘    └──────────────────┘    └───────────┘ │   │
│  └─────────────────────────────────────────────────────────────────┘   │
│                                   ↕                                     │
│  ┌─────────────────────────────────────────────────────────────────┐   │
│  │               BIDIRECTIONAL COUPLING                           │   │
│  │                                                                 │   │
│  │  ┌──────────────────┐    ┌──────────────────┐    ┌───────────┐ │   │
│  │  │ ATP Consumption  │───►│ Membrane         │───►│ Circuit   │ │   │
│  │  │ → Membrane       │    │ Property         │    │ Parameter │ │   │
│  │  │ Changes          │    │ Updates          │    │ Updates   │ │   │
│  │  └──────────────────┘    └──────────────────┘    └───────────┘ │   │
│  │                                                                 │   │
│  │  ┌──────────────────┐    ┌──────────────────┐    ┌───────────┐ │   │
│  │  │ Circuit State    │───►│ Membrane         │───►│ ATP       │ │   │
│  │  │ → ATP            │    │ Process          │    │ Demand    │ │   │
│  │  │ Demand           │    │ Activation       │    │ Updates   │ │   │
│  │  └──────────────────┘    └──────────────────┘    └───────────┘ │   │
│  └─────────────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────────────┘
```

## Membrane → Circuit Parameter Mapping

### Core Mapping Engine

```rust
pub struct MembraneCircuitMapper {
    // Membrane state inputs
    molecular_layer: MolecularLayerInterface,
    mesoscale_layer: MesoscaleLayerInterface,
    cellular_layer: CellularLayerInterface,
    
    // Circuit parameter calculators
    capacitance_calculator: CapacitanceCalculator,
    resistance_calculator: ResistanceCalculator,
    topology_generator: TopologyGenerator,
    
    // ATP coupling
    atp_membrane_coupling: AtpMembraneCoupling,
    
    // Nebuchadnezzar interface
    nebuchadnezzar_interface: NebuchadnezzarCircuitInterface,
    
    // Mapping history for optimization
    parameter_history: CircuitParameterHistory,
    mapping_performance: MappingPerformanceTracker,
}

impl MembraneCircuitMapper {
    pub fn map_membrane_to_circuit(&self, membrane_state: &MembraneState, atp_state: &AtpState) -> CircuitTopology {
        // 1. Extract membrane electrical properties from all layers
        let molecular_properties = self.molecular_layer.get_electrical_properties();
        let mesoscale_properties = self.mesoscale_layer.get_domain_properties();
        let cellular_properties = self.cellular_layer.get_network_properties();
        
        // 2. Calculate circuit elements from membrane patches
        let circuit_elements = self.generate_circuit_elements(
            &molecular_properties,
            &mesoscale_properties,
            &cellular_properties,
            atp_state,
        );
        
        // 3. Determine circuit topology from membrane connectivity
        let topology = self.topology_generator.generate_from_membrane_connectivity(
            &membrane_state.connectivity_graph,
            &circuit_elements,
        );
        
        // 4. Create ATP coupling matrix
        let atp_coupling = self.atp_membrane_coupling.create_coupling_matrix(
            membrane_state,
            atp_state,
            &circuit_elements,
        );
        
        CircuitTopology {
            elements: circuit_elements,
            connections: topology.connections,
            boundary_conditions: topology.boundary_conditions,
            atp_coupling,
            dynamic_update_rules: self.create_dynamic_update_rules(membrane_state),
        }
    }
    
    fn generate_circuit_elements(
        &self,
        molecular: &MolecularElectricalProperties,
        mesoscale: &MesoscaleDomainProperties,
        cellular: &CellularNetworkProperties,
        atp_state: &AtpState,
    ) -> Vec<CircuitElement> {
        let mut elements = Vec::new();
        
        // Membrane patches as RC circuits
        for patch in &molecular.membrane_patches {
            let capacitance = self.capacitance_calculator.calculate_patch_capacitance(
                patch,
                &mesoscale.local_properties[&patch.id],
                atp_state,
            );
            
            let resistance = self.resistance_calculator.calculate_patch_resistance(
                patch,
                &mesoscale.local_properties[&patch.id],
                atp_state,
            );
            
            elements.push(CircuitElement::MembraneRC {
                id: patch.id.clone(),
                capacitance: DynamicCapacitance::new(capacitance, patch.atp_dependencies.clone()),
                resistance: DynamicResistance::new(resistance, patch.atp_dependencies.clone()),
                area: patch.area,
                location: patch.location.clone(),
            });
        }
        
        // Proteins as active elements
        for protein in &molecular.membrane_proteins {
            let protein_circuit = self.convert_protein_to_circuit_element(protein, atp_state);
            elements.push(protein_circuit);
        }
        
        // Membrane contact sites as coupling elements
        for contact in &cellular.contact_sites {
            let coupling_element = self.create_contact_coupling_element(contact, atp_state);
            elements.push(coupling_element);
        }
        
        elements
    }
}
```

### Dynamic Circuit Parameter Updates

```rust
pub struct DynamicCircuitParameters {
    // Parameter evolution equations
    capacitance_evolution: CapacitanceEvolution,
    resistance_evolution: ResistanceEvolution,
    topology_evolution: TopologyEvolution,
    
    // ATP-dependent update rules
    atp_update_rules: HashMap<CircuitElementId, AtpUpdateRule>,
    
    // Membrane feedback mechanisms
    membrane_feedback: MembraneFeedbackMechanism,
}

impl DynamicCircuitParameters {
    pub fn update_circuit_for_atp_consumption(
        &mut self,
        circuit: &mut CircuitTopology,
        atp_consumption: &AtpConsumption,
        dt: f64,
    ) -> CircuitUpdateResult {
        let mut update_summary = CircuitUpdateResult::new();
        
        // Update each circuit element based on ATP consumption
        for element in &mut circuit.elements {
            match element {
                CircuitElement::MembraneRC { capacitance, resistance, .. } => {
                    let cap_update = self.capacitance_evolution.update_for_atp(
                        capacitance,
                        &atp_consumption.spatial_distribution,
                        dt,
                    );
                    let res_update = self.resistance_evolution.update_for_atp(
                        resistance,
                        &atp_consumption.spatial_distribution,
                        dt,
                    );
                    
                    update_summary.add_parameter_changes(cap_update, res_update);
                }
                
                CircuitElement::AtpPump { pump_parameters, .. } => {
                    let pump_update = pump_parameters.update_for_atp_availability(
                        atp_consumption.local_atp_depletion,
                        dt,
                    );
                    update_summary.add_pump_update(pump_update);
                }
                
                CircuitElement::AtpSynthase { synthase_parameters, .. } => {
                    let synthase_update = synthase_parameters.update_for_atp_production(
                        atp_consumption.local_atp_demand,
                        dt,
                    );
                    update_summary.add_synthase_update(synthase_update);
                }
                
                _ => {}
            }
        }
        
        // Update circuit topology if membrane connectivity changed
        if atp_consumption.topology_affecting_processes.is_empty() == false {
            let topology_update = self.topology_evolution.update_for_membrane_changes(
                circuit,
                &atp_consumption.topology_affecting_processes,
                dt,
            );
            update_summary.add_topology_changes(topology_update);
        }
        
        update_summary
    }
    
    pub fn predict_parameter_evolution(
        &self,
        current_circuit: &CircuitTopology,
        atp_trajectory: &AtpTrajectory,
        time_horizon: f64,
    ) -> PredictedCircuitEvolution {
        // Predict how circuit parameters will evolve over time
        let mut predicted_evolution = PredictedCircuitEvolution::new(time_horizon);
        
        let dt = 0.001; // 1 ms time steps
        let steps = (time_horizon / dt) as usize;
        
        let mut current_state = current_circuit.clone();
        
        for step in 0..steps {
            let time = step as f64 * dt;
            let atp_consumption = atp_trajectory.get_consumption_at_time(time);
            
            let update_result = self.update_circuit_for_atp_consumption(
                &mut current_state,
                &atp_consumption,
                dt,
            );
            
            predicted_evolution.add_time_point(time, current_state.clone(), update_result);
        }
        
        predicted_evolution
    }
}
```

## Nebuchadnezzar Interface Implementation

### ATP-Based Differential Equations for Membrane Dynamics

```rust
pub struct MembraneAtpDifferentialEquations {
    // Membrane state variables
    membrane_variables: MembraneStateVariables,
    
    // ATP coupling coefficients
    atp_coupling_matrix: AtpCouplingMatrix,
    
    // Rate equations
    lipid_synthesis_rates: LipidSynthesisRates,
    protein_insertion_rates: ProteinInsertionRates,
    membrane_remodeling_rates: MembraneRemodelingRates,
    
    // Circuit parameter dependencies
    circuit_parameter_derivatives: CircuitParameterDerivatives,
}

impl MembraneAtpDifferentialEquations {
    pub fn calculate_membrane_derivatives(&self, atp_consumption_rate: f64) -> MembraneDerivatives {
        let mut derivatives = MembraneDerivatives::new();
        
        // Lipid concentration changes: dLipid/dATP
        derivatives.lipid_concentrations = self.lipid_synthesis_rates.calculate_datp_rates(
            atp_consumption_rate,
            &self.membrane_variables.current_lipid_concentrations,
        );
        
        // Protein density changes: dProtein/dATP
        derivatives.protein_densities = self.protein_insertion_rates.calculate_datp_rates(
            atp_consumption_rate,
            &self.membrane_variables.current_protein_densities,
        );
        
        // Membrane curvature changes: dCurvature/dATP
        derivatives.membrane_curvature = self.membrane_remodeling_rates.calculate_curvature_datp(
            atp_consumption_rate,
            self.membrane_variables.current_curvature,
        );
        
        // Membrane potential changes: dV/dATP
        derivatives.membrane_potential = self.calculate_potential_datp_dependence(
            atp_consumption_rate,
            &self.membrane_variables,
        );
        
        // Circuit parameter derivatives: dR/dATP, dC/dATP
        derivatives.circuit_parameters = self.circuit_parameter_derivatives.calculate_datp_derivatives(
            atp_consumption_rate,
            &derivatives,
        );
        
        derivatives
    }
    
    fn calculate_potential_datp_dependence(
        &self,
        atp_consumption_rate: f64,
        membrane_vars: &MembraneStateVariables,
    ) -> f64 {
        // Membrane potential changes due to ATP-dependent ion pump activity
        let pump_contribution = self.calculate_pump_potential_contribution(atp_consumption_rate);
        let conductance_changes = self.calculate_conductance_potential_effects(
            atp_consumption_rate,
            membrane_vars,
        );
        
        pump_contribution + conductance_changes
    }
    
    pub fn integrate_with_nebuchadnezzar(
        &mut self,
        nebuchadnezzar_state: &NebuchadnezzarState,
    ) -> IntegrationResult {
        // Bidirectional coupling between membrane and Nebuchadnezzar
        
        // 1. Get ATP consumption from Nebuchadnezzar circuit simulation
        let atp_consumption = nebuchadnezzar_state.get_current_atp_consumption_rate();
        
        // 2. Calculate membrane changes due to ATP consumption
        let membrane_derivatives = self.calculate_membrane_derivatives(atp_consumption);
        
        // 3. Update membrane state
        self.membrane_variables.integrate_derivatives(&membrane_derivatives);
        
        // 4. Calculate new circuit parameters from updated membrane state
        let new_circuit_params = self.calculate_updated_circuit_parameters();
        
        // 5. Send updated circuit parameters back to Nebuchadnezzar
        let circuit_update = CircuitUpdate {
            new_parameters: new_circuit_params,
            topology_changes: self.calculate_topology_changes(),
            boundary_condition_updates: self.calculate_boundary_updates(),
        };
        
        IntegrationResult {
            membrane_state_change: membrane_derivatives,
            circuit_parameter_updates: circuit_update,
            atp_demand_feedback: self.calculate_atp_demand_changes(),
            integration_stability: self.assess_integration_stability(),
        }
    }
}
```

### Hierarchical Circuit Abstraction Interface

```rust
pub struct HierarchicalMembraneCircuitInterface {
    // Abstraction levels
    molecular_level_circuits: MolecularLevelCircuits,
    mesoscale_level_circuits: MesoscaleLevelCircuits,
    cellular_level_circuits: CellularLevelCircuits,
    
    // Abstraction control
    abstraction_controller: AbstractionController,
    expansion_criteria: ExpansionCriteria,
    
    // Probabilistic node management
    probabilistic_nodes: HashMap<NodeId, ProbabilisticMembraneNode>,
    detailed_expansions: HashMap<NodeId, DetailedCircuitExpansion>,
}

impl HierarchicalMembraneCircuitInterface {
    pub fn create_hierarchical_representation(
        &mut self,
        full_membrane_circuit: &CircuitTopology,
        computational_budget: ComputationalBudget,
    ) -> HierarchicalCircuitRepresentation {
        // Start with most abstract representation
        let mut hierarchical_repr = HierarchicalCircuitRepresentation::new();
        
        // Create probabilistic nodes for membrane regions
        for region in self.identify_membrane_regions(full_membrane_circuit) {
            let probabilistic_node = self.create_probabilistic_membrane_node(region);
            hierarchical_repr.add_node(probabilistic_node);
        }
        
        // Determine which nodes to expand based on importance and budget
        let expansion_candidates = self.abstraction_controller.select_expansion_candidates(
            &hierarchical_repr,
            &self.expansion_criteria,
            computational_budget,
        );
        
        // Expand high-priority nodes to detailed circuits
        for candidate in expansion_candidates {
            let detailed_circuit = self.expand_node_to_detailed_circuit(
                &candidate,
                full_membrane_circuit,
            );
            hierarchical_repr.replace_node_with_detailed_circuit(candidate.node_id, detailed_circuit);
        }
        
        hierarchical_repr
    }
    
    fn create_probabilistic_membrane_node(&self, region: &MembraneRegion) -> ProbabilisticMembraneNode {
        // Create uncertain representation of membrane region
        let avg_capacitance = region.calculate_average_capacitance();
        let avg_resistance = region.calculate_average_resistance();
        
        // Estimate uncertainty based on spatial variation
        let capacitance_uncertainty = region.calculate_capacitance_variation();
        let resistance_uncertainty = region.calculate_resistance_variation();
        
        ProbabilisticMembraneNode {
            node_id: region.id.clone(),
            location: region.centroid(),
            area: region.total_area(),
            
            // Probabilistic electrical properties
            capacitance_distribution: ProbabilityDistribution::normal(
                avg_capacitance,
                capacitance_uncertainty,
            ),
            resistance_distribution: ProbabilityDistribution::normal(
                avg_resistance,
                resistance_uncertainty,
            ),
            
            // ATP dependencies
            atp_sensitivity: region.calculate_atp_sensitivity(),
            atp_consumption_rate: region.estimate_atp_consumption(),
            
            // Expansion criteria
            importance_score: self.calculate_region_importance(region),
            uncertainty_score: capacitance_uncertainty + resistance_uncertainty,
            computational_cost: self.estimate_expansion_cost(region),
        }
    }
    
    pub fn adaptive_expansion(
        &mut self,
        hierarchical_repr: &mut HierarchicalCircuitRepresentation,
        optimization_results: &OptimizationResults,
    ) -> ExpansionResult {
        // Adaptively expand nodes based on optimization sensitivity
        let sensitive_nodes = optimization_results.identify_sensitive_parameters();
        
        let mut expansion_result = ExpansionResult::new();
        
        for node_id in sensitive_nodes {
            if let Some(probabilistic_node) = self.probabilistic_nodes.get(&node_id) {
                // Check if expansion would improve optimization accuracy
                let expansion_benefit = self.estimate_expansion_benefit(
                    probabilistic_node,
                    &optimization_results,
                );
                
                if expansion_benefit > self.expansion_criteria.benefit_threshold {
                    let detailed_circuit = self.expand_node_to_detailed_circuit(
                        probabilistic_node,
                        &hierarchical_repr.get_full_circuit(),
                    );
                    
                    hierarchical_repr.replace_node_with_detailed_circuit(
                        node_id,
                        detailed_circuit,
                    );
                    
                    expansion_result.add_expansion(node_id, expansion_benefit);
                }
            }
        }
        
        expansion_result
    }
}
```

## Probabilistic Circuit Priors

### Membrane-Based Circuit Priors

```rust
pub struct MembraneCircuitPriors {
    // Prior knowledge from membrane biophysics
    lipid_bilayer_priors: LipidBilayerPriors,
    protein_function_priors: ProteinFunctionPriors,
    membrane_organization_priors: MembraneOrganizationPriors,
    
    // Experimental parameter distributions
    experimental_parameter_db: ExperimentalParameterDatabase,
    
    // Bayesian prior update mechanisms
    prior_update_engine: BayesianPriorUpdateEngine,
}

impl MembraneCircuitPriors {
    pub fn generate_circuit_parameter_priors(
        &self,
        membrane_composition: &MembraneComposition,
        temperature: f64,
        atp_concentration: f64,
    ) -> CircuitParameterPriors {
        let mut priors = CircuitParameterPriors::new();
        
        // Capacitance priors from lipid bilayer physics
        let capacitance_prior = self.lipid_bilayer_priors.generate_capacitance_prior(
            membrane_composition,
            temperature,
        );
        priors.add_capacitance_prior(capacitance_prior);
        
        // Resistance priors from membrane permeability
        let resistance_prior = self.lipid_bilayer_priors.generate_resistance_prior(
            membrane_composition,
            temperature,
        );
        priors.add_resistance_prior(resistance_prior);
        
        // Protein circuit element priors
        for protein_type in membrane_composition.get_protein_types() {
            let protein_prior = self.protein_function_priors.generate_protein_circuit_prior(
                protein_type,
                atp_concentration,
            );
            priors.add_protein_prior(protein_type, protein_prior);
        }
        
        // Membrane organization priors (domains, contacts, etc.)
        let organization_priors = self.membrane_organization_priors.generate_organization_priors(
            membrane_composition,
        );
        priors.add_organization_priors(organization_priors);
        
        priors
    }
    
    pub fn update_priors_from_experimental_data(
        &mut self,
        experimental_data: &ExperimentalData,
        membrane_conditions: &MembraneConditions,
    ) -> PriorUpdateResult {
        // Update Bayesian priors based on new experimental observations
        let likelihood = self.calculate_experimental_likelihood(
            experimental_data,
            membrane_conditions,
        );
        
        let prior_updates = self.prior_update_engine.bayesian_update(
            &self.current_priors(),
            &likelihood,
        );
        
        // Apply updates to all prior distributions
        self.lipid_bilayer_priors.update_from_bayesian_results(&prior_updates);
        self.protein_function_priors.update_from_bayesian_results(&prior_updates);
        self.membrane_organization_priors.update_from_bayesian_results(&prior_updates);
        
        // Store experimental data for future reference
        self.experimental_parameter_db.add_experimental_data(
            experimental_data.clone(),
            membrane_conditions.clone(),
        );
        
        PriorUpdateResult {
            updated_parameters: prior_updates.get_updated_parameters(),
            confidence_changes: prior_updates.get_confidence_changes(),
            prediction_improvements: self.assess_prediction_improvements(&prior_updates),
        }
    }
}

#[derive(Debug, Clone)]
pub struct LipidBilayerPriors {
    // Capacitance prior distributions
    specific_capacitance_prior: ProbabilityDistribution<f64>, // F/m²
    thickness_dependence: ThicknessDependencePrior,
    temperature_dependence: TemperatureDependencePrior,
    
    // Resistance prior distributions  
    ion_permeability_priors: HashMap<IonType, ProbabilityDistribution<f64>>,
    lipid_composition_effects: LipidCompositionEffectPriors,
    
    // Uncertainty quantification
    measurement_uncertainties: MeasurementUncertainties,
    model_uncertainties: ModelUncertainties,
}

impl LipidBilayerPriors {
    pub fn generate_capacitance_prior(
        &self,
        composition: &MembraneComposition,
        temperature: f64,
    ) -> CapacitancePrior {
        // Start with base specific capacitance distribution
        let mut capacitance_dist = self.specific_capacitance_prior.clone();
        
        // Adjust for membrane composition
        let composition_factor = self.calculate_composition_capacitance_factor(composition);
        capacitance_dist = capacitance_dist.scale(composition_factor);
        
        // Adjust for temperature
        let temperature_factor = self.temperature_dependence.calculate_capacitance_factor(temperature);
        capacitance_dist = capacitance_dist.scale(temperature_factor);
        
        // Include thickness uncertainty
        let thickness_uncertainty = self.thickness_dependence.estimate_thickness_uncertainty(composition);
        capacitance_dist = capacitance_dist.convolve_with_uncertainty(thickness_uncertainty);
        
        CapacitancePrior {
            distribution: capacitance_dist,
            temperature: temperature,
            composition: composition.clone(),
            confidence: self.calculate_prior_confidence(composition, temperature),
        }
    }
}
```

## Example: Complete Membrane-Circuit Integration

```rust
// Example: Complete workflow for membrane-circuit integration
pub struct MembraneCircuitIntegrationExample;

impl MembraneCircuitIntegrationExample {
    pub fn run_complete_integration_example() -> Result<IntegrationResults, IntegrationError> {
        // 1. Initialize membrane system
        let mut membrane_system = MembraneSystem::new()
            .with_lipid_composition(LipidComposition::mammalian_plasma_membrane())
            .with_proteins(vec![
                ProteinType::SodiumPotassiumPump,
                ProteinType::VoltageGatedSodiumChannel,
                ProteinType::VoltageGatedPotassiumChannel,
                ProteinType::AtpSynthase,
            ])
            .with_temperature(310.0) // 37°C
            .with_membrane_area(1e-9) // 1 μm²
            .build()?;
        
        // 2. Initialize ATP pool
        let mut atp_pool = AtpPool::physiological()
            .with_concentration(5.0) // 5 mM ATP
            .with_energy_charge(0.85)
            .build();
        
        // 3. Create circuit interface
        let mut circuit_interface = CircuitInterfaceLayer::new()
            .with_membrane_system(&membrane_system)
            .with_atp_coupling(AtpCouplingMode::Full)
            .with_hierarchy_enabled(true)
            .build()?;
        
        // 4. Generate initial circuit topology
        let initial_circuit = circuit_interface.map_membrane_to_circuit(
            &membrane_system.state(),
            &atp_pool.state(),
        );
        
        // 5. Create Nebuchadnezzar system with membrane-derived circuit
        let mut nebuchadnezzar = NebuchadnezzarSystem::new()
            .with_circuit_topology(initial_circuit)
            .with_atp_pool(atp_pool.clone())
            .with_solver(SolverType::AtpBasedRungeKutta4)
            .with_objective(OptimizationObjective::MaximizeATPEfficiency)
            .build()?;
        
        // 6. Coupled simulation loop
        let mut integration_results = IntegrationResults::new();
        let simulation_time = 1.0; // 1 second
        let dt_atp = 0.001; // 1 ms ATP time steps
        
        for step in 0..(simulation_time / dt_atp) as usize {
            let current_time = step as f64 * dt_atp;
            
            // Nebuchadnezzar ATP-based integration
            let nebuchadnezzar_result = nebuchadnezzar.integrate_step_atp(dt_atp)?;
            
            // Update membrane system based on ATP consumption
            let membrane_update = membrane_system.update_for_atp_consumption(
                &nebuchadnezzar_result.atp_consumption,
                dt_atp,
            )?;
            
            // Update circuit parameters from membrane changes
            let circuit_update = circuit_interface.update_circuit_from_membrane_changes(
                &mut nebuchadnezzar.circuit,
                &membrane_update.membrane_state_changes,
            )?;
            
            // Bidirectional feedback
            let feedback_result = circuit_interface.apply_bidirectional_feedback(
                &mut membrane_system,
                &mut nebuchadnezzar,
            )?;
            
            // Record results
            integration_results.add_time_point(IntegrationTimePoint {
                time: current_time,
                membrane_state: membrane_system.state().clone(),
                circuit_state: nebuchadnezzar.circuit.state().clone(),
                atp_state: nebuchadnezzar.atp_pool.state().clone(),
                integration_metrics: feedback_result.metrics,
            });
            
            // Adaptive hierarchy adjustment
            if step % 100 == 0 { // Every 100 ms
                circuit_interface.adaptive_hierarchy_adjustment(
                    &nebuchadnezzar.optimization_state(),
                )?;
            }
        }
        
        // 7. Analysis and optimization
        let optimization_results = nebuchadnezzar.finalize_optimization()?;
        integration_results.add_optimization_results(optimization_results);
        
        // 8. Generate final circuit parameters optimized for ATP efficiency
        let optimized_circuit = circuit_interface.generate_optimized_circuit(
            &integration_results.final_membrane_state(),
            &optimization_results,
        )?;
        
        integration_results.set_final_optimized_circuit(optimized_circuit);
        
        Ok(integration_results)
    }
    
    pub fn analyze_integration_results(results: &IntegrationResults) -> IntegrationAnalysis {
        IntegrationAnalysis {
            atp_efficiency_over_time: results.calculate_atp_efficiency_trajectory(),
            membrane_parameter_evolution: results.calculate_membrane_parameter_changes(),
            circuit_topology_changes: results.identify_circuit_topology_evolution(),
            optimization_convergence: results.analyze_optimization_convergence(),
            stability_metrics: results.calculate_stability_metrics(),
            biological_realism_assessment: results.assess_biological_realism(),
        }
    }
}
```

---

The **Circuit Interface Layer** provides the essential bridge between biophysical membrane dynamics and Nebuchadnezzar's hierarchical probabilistic electric circuits, enabling authentic biological simulation through ATP-based differential equations while maintaining computational efficiency through adaptive hierarchical abstraction. 

## Orchestrator Communication Protocol

### Incoming Orchestrator Commands

The circuit interface receives orchestration signals that influence membrane-circuit mapping:

```python
class OrchestratorInterface:
    def receive_context_update(self, context_layer_state):
        """Context layer provides system-wide membrane awareness"""
        self.global_membrane_context = context_layer_state
        self.adjust_circuit_priorities()
    
    def receive_reasoning_directive(self, reasoning_decision):
        """Reasoning layer coordinates parameter mapping strategies"""
        self.parameter_mapping_strategy = reasoning_decision
        self.update_mapping_algorithms()
    
    def receive_intuition_prediction(self, predicted_changes):
        """Intuition layer predicts membrane behavior changes"""
        self.preload_circuit_adaptations(predicted_changes)
    
    def receive_atp_allocation(self, atp_budget):
        """V8 metabolism pipeline sets ATP constraints"""
        self.membrane_atp_budget = atp_budget
        self.prioritize_essential_processes()
```

### Outgoing Status Reports

The interface continuously reports membrane state to the orchestrator:

```python
def report_to_orchestrator(self):
    return {
        'membrane_voltage_state': self.current_membrane_potentials,
        'atp_consumption_rate': self.calculate_atp_usage(),
        'circuit_parameter_changes': self.recent_parameter_updates,
        'biological_constraint_violations': self.validate_realism(),
        'predicted_membrane_transitions': self.forecast_changes()
    }
``` 