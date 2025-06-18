# Molecular Layer - Fundamental Membrane Biophysics

## Core Philosophy

**"Every membrane property emerges from molecular interactions - and every molecular interaction can be expressed as circuit parameters"**

The **Molecular Layer** implements the fundamental biophysical processes that determine membrane electrical properties. This layer translates lipid bilayer physics, protein-membrane interactions, membrane curvature dynamics, and electrochemical gradients into the circuit parameters that Nebuchadnezzar requires for ATP-based differential equations.

## Lipid Bilayer Physics

### Membrane Capacitance Calculation

The bilayer capacitance is the most fundamental electrical property of membranes:

```rust
pub struct LipidBilayerCapacitance {
    // Physical constants
    vacuum_permittivity: f64,       // ε₀ = 8.854 × 10⁻¹² F/m
    
    // Composition-dependent properties
    lipid_composition: LipidComposition,
    cholesterol_fraction: f64,      // 0.0 - 1.0
    temperature: f64,               // Kelvin
    
    // Dynamic properties
    membrane_thickness: f64,        // nm - varies with composition
    dielectric_profile: DielectricProfile,
    
    // ATP-dependent modifications
    atp_dependent_lipid_synthesis: AtpLipidSynthesis,
    membrane_asymmetry: MembraneAsymmetry,
}

impl LipidBilayerCapacitance {
    pub fn calculate_specific_capacitance(&self) -> f64 {
        // Multi-layer dielectric model for biological membranes
        let effective_permittivity = self.calculate_effective_permittivity();
        let effective_thickness = self.calculate_effective_thickness();
        
        self.vacuum_permittivity * effective_permittivity / (effective_thickness * 1e-9)
    }
    
    fn calculate_effective_permittivity(&self) -> f64 {
        // Weighted average based on lipid composition
        let mut effective_permittivity = 0.0;
        
        for (lipid_type, fraction) in &self.lipid_composition.lipid_fractions {
            let lipid_permittivity = match lipid_type {
                LipidType::PhosphatidylCholine => 2.3,  // PC headgroup
                LipidType::PhosphatidylSerine => 3.1,   // PS headgroup (charged)
                LipidType::Cholesterol => 2.1,          // Cholesterol
                LipidType::Sphingomyelin => 2.5,        // SM headgroup
                _ => 2.2, // Default value
            };
            effective_permittivity += lipid_permittivity * fraction;
        }
        
        // Temperature correction
        let temperature_factor = 1.0 + 0.002 * (self.temperature - 298.0); // ~0.2% per K
        effective_permittivity * temperature_factor
    }
    
    fn calculate_effective_thickness(&self) -> f64 {
        // Bilayer thickness depends on acyl chain length and saturation
        let base_thickness = self.calculate_base_thickness();
        let cholesterol_correction = 1.0 + 0.15 * self.cholesterol_fraction; // Cholesterol increases thickness
        let temperature_correction = 1.0 - 0.001 * (self.temperature - 298.0); // Higher T decreases thickness
        
        base_thickness * cholesterol_correction * temperature_correction
    }
    
    pub fn update_for_atp_consumption(&mut self, atp_consumed: f64, dt: f64) {
        // ATP consumption affects membrane through lipid synthesis
        let lipid_synthesis_rate = self.atp_dependent_lipid_synthesis.calculate_synthesis_rate(atp_consumed);
        
        // Update membrane composition
        self.lipid_composition.update_from_synthesis(lipid_synthesis_rate, dt);
        
        // Recalculate thickness and permittivity
        self.membrane_thickness = self.calculate_effective_thickness();
    }
}
```

### Membrane Resistance and Permeability

```rust
pub struct MembranePermeability {
    // Ion-specific permeabilities (m/s)
    sodium_permeability: f64,
    potassium_permeability: f64,
    chloride_permeability: f64,
    calcium_permeability: f64,
    
    // Composition dependencies
    lipid_composition: LipidComposition,
    cholesterol_effect: CholesterolEffect,
    temperature: f64,
    
    // ATP-dependent modifications
    atp_dependent_flippases: Vec<AtpFlippase>,
    asymmetric_composition: AsymmetricComposition,
}

impl MembranePermeability {
    pub fn calculate_ion_resistance(&self, ion_type: IonType, membrane_area: f64) -> f64 {
        let permeability = match ion_type {
            IonType::Sodium => self.sodium_permeability,
            IonType::Potassium => self.potassium_permeability,
            IonType::Chloride => self.chloride_permeability,
            IonType::Calcium => self.calcium_permeability,
        };
        
        // Apply composition corrections
        let composition_factor = self.calculate_composition_factor(ion_type);
        let temperature_factor = self.calculate_temperature_factor();
        
        let effective_permeability = permeability * composition_factor * temperature_factor;
        
        // Convert permeability to resistance: R = 1/(P × A × F²/RT × z²)
        let faraday_constant = 96485.0; // C/mol
        let gas_constant = 8.314; // J/(mol⋅K)
        let ion_charge = ion_type.charge();
        
        1.0 / (effective_permeability * membrane_area * 
               faraday_constant.powi(2) / (gas_constant * self.temperature) * 
               ion_charge.powi(2))
    }
    
    fn calculate_composition_factor(&self, ion_type: IonType) -> f64 {
        // Different lipids have different effects on ion permeability
        let mut factor = 1.0;
        
        // Cholesterol generally decreases permeability
        factor *= 1.0 - 0.5 * self.cholesterol_effect.permeability_reduction * 
                  self.lipid_composition.cholesterol_fraction();
        
        // Charged lipids affect charged ion permeability
        if ion_type.is_cation() {
            let negative_lipid_fraction = self.lipid_composition.get_fraction(LipidType::PhosphatidylSerine);
            factor *= 1.0 + 0.3 * negative_lipid_fraction; // PS increases cation permeability
        }
        
        factor
    }
    
    pub fn update_for_atp_flippase_activity(&mut self, atp_consumed: f64) {
        // ATP-dependent flippases change membrane asymmetry
        for flippase in &mut self.atp_dependent_flippases {
            let flippase_activity = flippase.calculate_activity(atp_consumed);
            self.asymmetric_composition.update_from_flippase_activity(flippase_activity);
        }
        
        // Asymmetry affects permeability
        self.update_permeabilities_from_asymmetry();
    }
}
```

## Protein-Membrane Interactions

### Transmembrane Protein Electrical Models

```rust
pub struct TransmembraneProteinElectrics {
    protein_id: String,
    protein_type: ProteinType,
    
    // Structural properties
    transmembrane_domains: u8,
    hydrophobic_mismatch: f64,    // nm - difference between protein and membrane thickness
    protein_radius: f64,          // nm - effective radius
    
    // Electrical properties
    intrinsic_conductance: f64,   // S - protein conductance
    voltage_dependence: VoltageGating,
    ion_selectivity: IonSelectivity,
    
    // ATP-dependent properties
    atp_binding_sites: Vec<AtpBindingSite>,
    conformational_states: HashMap<AtpConcentration, ProteinConformation>,
    
    // Membrane perturbation
    lipid_annulus: LipidAnnulus,  // Perturbed lipids around protein
    membrane_curvature_induction: f64,
}

impl TransmembraneProteinElectrics {
    pub fn calculate_equivalent_circuit(&self, membrane_voltage: f64, atp_concentration: f64) -> ProteinEquivalentCircuit {
        match &self.protein_type {
            ProteinType::IonChannel { gating_mechanism, .. } => {
                self.calculate_ion_channel_circuit(membrane_voltage, atp_concentration)
            }
            ProteinType::AtpPump { .. } => {
                self.calculate_atp_pump_circuit(membrane_voltage, atp_concentration)
            }
            ProteinType::AtpSynthase { .. } => {
                self.calculate_atp_synthase_circuit(membrane_voltage, atp_concentration)
            }
            _ => ProteinEquivalentCircuit::default(),
        }
    }
    
    fn calculate_ion_channel_circuit(&self, voltage: f64, atp_conc: f64) -> ProteinEquivalentCircuit {
        // Calculate open probability from voltage and ATP dependence
        let open_probability = self.calculate_open_probability(voltage, atp_conc);
        
        // Effective conductance
        let effective_conductance = self.intrinsic_conductance * open_probability;
        
        // Channel as voltage-dependent resistor
        ProteinEquivalentCircuit::VoltageGatedResistor {
            max_conductance: self.intrinsic_conductance,
            current_conductance: effective_conductance,
            reversal_potential: self.calculate_reversal_potential(),
            gating_dynamics: self.calculate_gating_dynamics(voltage, atp_conc),
        }
    }
    
    fn calculate_atp_pump_circuit(&self, voltage: f64, atp_conc: f64) -> ProteinEquivalentCircuit {
        // ATP pump as current source with internal resistance
        let pump_current = self.calculate_pump_current(voltage, atp_conc);
        let internal_resistance = self.calculate_pump_resistance(atp_conc);
        let atp_consumption_rate = self.calculate_atp_consumption_rate(pump_current);
        
        ProteinEquivalentCircuit::CurrentSource {
            current: pump_current,
            internal_resistance,
            voltage_dependence: self.calculate_pump_voltage_dependence(),
            atp_consumption_rate,
            power_consumption: pump_current * voltage,
        }
    }
    
    fn calculate_atp_synthase_circuit(&self, voltage: f64, atp_conc: f64) -> ProteinEquivalentCircuit {
        // ATP synthase as reversible current source/sink
        let synthase_current = self.calculate_synthase_current(voltage, atp_conc);
        let internal_resistance = self.calculate_synthase_resistance(atp_conc);
        let atp_production_rate = self.calculate_atp_production_rate(synthase_current);
        
        ProteinEquivalentCircuit::ReversibleCurrentSource {
            current: synthase_current,
            internal_resistance,
            forward_bias_voltage: self.calculate_forward_bias_voltage(atp_conc),
            reverse_bias_voltage: self.calculate_reverse_bias_voltage(atp_conc),
            atp_production_rate,
            power_generation: -synthase_current * voltage, // Negative for power generation
        }
    }
    
    pub fn calculate_membrane_perturbation(&self) -> MembranePertu rbation {
        // How protein affects surrounding membrane
        let thickness_perturbation = self.calculate_thickness_perturbation();
        let curvature_induction = self.calculate_curvature_induction();
        let lipid_ordering = self.calculate_lipid_ordering_effect();
        
        MembraneP erturbation {
            affected_radius: self.protein_radius + 2.0, // ~2 nm perturbation radius
            thickness_change: thickness_perturbation,
            curvature_change: curvature_induction,
            lipid_ordering_change: lipid_ordering,
            capacitance_change: self.calculate_capacitance_perturbation(),
            resistance_change: self.calculate_resistance_perturbation(),
        }
    }
}
```

### ATP-Dependent Protein Conformational Changes

```rust
pub struct AtpDependentConformations {
    protein_id: String,
    
    // ATP binding sites and their effects
    atp_binding_sites: Vec<AtpBindingSite>,
    cooperative_binding: CooperativeBinding,
    
    // Conformational states
    conformational_states: HashMap<ConformationId, ProteinConformation>,
    transition_rates: HashMap<(ConformationId, ConformationId), f64>,
    
    // Electrical consequences
    conformation_conductances: HashMap<ConformationId, f64>,
    conformation_capacitances: HashMap<ConformationId, f64>,
}

#[derive(Debug, Clone)]
pub struct AtpBindingSite {
    site_id: String,
    binding_affinity: f64,        // M⁻¹ - association constant
    cooperativity_factor: f64,    // Hill coefficient
    conformational_effect: ConformationalEffect,
}

#[derive(Debug, Clone)]
pub enum ConformationalEffect {
    ConductanceChange(f64),       // Factor by which conductance changes
    GatingEffect(GatingEffect),   // Changes in voltage/time dependence
    SelectivityChange(SelectivityChange), // Changes in ion selectivity
    ATPConsumptionChange(f64),    // Changes in ATP consumption rate
}

impl AtpDependentConformations {
    pub fn calculate_current_conformation(&self, atp_concentration: f64) -> ProteinConformation {
        // Calculate occupancy of each ATP binding site
        let site_occupancies = self.calculate_site_occupancies(atp_concentration);
        
        // Determine most probable conformation
        let conformation_probabilities = self.calculate_conformation_probabilities(&site_occupancies);
        
        // Weighted average of conformational properties
        self.average_conformational_properties(&conformation_probabilities)
    }
    
    fn calculate_site_occupancies(&self, atp_conc: f64) -> Vec<f64> {
        let mut occupancies = Vec::new();
        
        for site in &self.atp_binding_sites {
            // Hill equation for cooperative binding
            let occupancy = (site.binding_affinity * atp_conc).powf(site.cooperativity_factor) /
                          (1.0 + (site.binding_affinity * atp_conc).powf(site.cooperativity_factor));
            occupancies.push(occupancy);
        }
        
        occupancies
    }
    
    fn calculate_conformation_probabilities(&self, site_occupancies: &[f64]) -> HashMap<ConformationId, f64> {
        let mut probabilities = HashMap::new();
        
        // Each combination of site occupancies corresponds to a conformation
        for (conformation_id, conformation) in &self.conformational_states {
            let mut probability = 1.0;
            
            for (site_idx, required_occupancy) in conformation.required_site_states.iter().enumerate() {
                if *required_occupancy {
                    probability *= site_occupancies[site_idx];
                } else {
                    probability *= 1.0 - site_occupancies[site_idx];
                }
            }
            
            probabilities.insert(*conformation_id, probability);
        }
        
        // Normalize probabilities
        let total_probability: f64 = probabilities.values().sum();
        for probability in probabilities.values_mut() {
            *probability /= total_probability;
        }
        
        probabilities
    }
    
    pub fn calculate_atp_dependent_conductance(&self, atp_concentration: f64) -> f64 {
        let conformation_probabilities = self.calculate_conformation_probabilities(
            &self.calculate_site_occupancies(atp_concentration)
        );
        
        let mut weighted_conductance = 0.0;
        for (conformation_id, probability) in conformation_probabilities {
            let conductance = self.conformation_conductances[&conformation_id];
            weighted_conductance += conductance * probability;
        }
        
        weighted_conductance
    }
}
```

## Membrane Curvature Dynamics

### Curvature-Electrical Property Coupling

```rust
pub struct MembraneCurvatureDynamics {
    // Curvature properties
    spontaneous_curvature: f64,   // m⁻¹ - intrinsic curvature preference
    bending_modulus: f64,         // J - resistance to bending
    gaussian_modulus: f64,        // J - resistance to topology changes
    
    // Current membrane shape
    mean_curvature: f64,          // m⁻¹ - average curvature
    gaussian_curvature: f64,      // m⁻² - topology measure
    
    // Protein-curvature coupling
    curvature_sensing_proteins: Vec<CurvatureSensingProtein>,
    curvature_inducing_proteins: Vec<CurvatureInducingProtein>,
    
    // ATP-dependent curvature generation
    atp_dependent_curvature_generators: Vec<AtpCurvatureGenerator>,
}

impl MembraneCurvatureDynamics {
    pub fn calculate_curvature_electrical_effects(&self) -> CurvatureElectricalEffects {
        // Curvature affects membrane electrical properties
        let capacitance_effect = self.calculate_curvature_capacitance_effect();
        let resistance_effect = self.calculate_curvature_resistance_effect();
        let protein_activity_effect = self.calculate_curvature_protein_effects();
        
        CurvatureElectricalEffects {
            capacitance_modification: capacitance_effect,
            resistance_modification: resistance_effect,
            protein_activity_modifications: protein_activity_effect,
        }
    }
    
    fn calculate_curvature_capacitance_effect(&self) -> f64 {
        // Curvature changes effective membrane area and thickness
        let area_factor = 1.0 + self.mean_curvature.abs() * 0.1; // Curved membranes have larger effective area
        let thickness_factor = 1.0 - self.mean_curvature.abs() * 0.05; // Curvature slightly reduces effective thickness
        
        area_factor / thickness_factor
    }
    
    fn calculate_curvature_resistance_effect(&self) -> f64 {
        // Curvature affects lipid packing and permeability
        let packing_factor = 1.0 + self.mean_curvature.abs() * 0.2; // Curvature loosens packing
        let permeability_factor = 1.0 + packing_factor * 0.3; // Looser packing increases permeability
        
        1.0 / permeability_factor // Resistance is inverse of permeability
    }
    
    fn calculate_curvature_protein_effects(&self) -> HashMap<String, f64> {
        let mut effects = HashMap::new();
        
        // Curvature-sensing proteins change activity with curvature
        for protein in &self.curvature_sensing_proteins {
            let activity_change = protein.calculate_curvature_response(self.mean_curvature);
            effects.insert(protein.protein_id.clone(), activity_change);
        }
        
        effects
    }
    
    pub fn update_curvature_for_atp_consumption(&mut self, atp_consumed: f64, dt: f64) {
        // ATP-dependent processes can generate membrane curvature
        for generator in &mut self.atp_dependent_curvature_generators {
            let curvature_change = generator.calculate_curvature_generation(atp_consumed, dt);
            self.apply_curvature_change(curvature_change);
        }
        
        // Update protein activities based on new curvature
        self.update_protein_curvature_responses();
    }
    
    pub fn calculate_curvature_energy(&self) -> f64 {
        // Helfrich bending energy
        let bending_energy = 0.5 * self.bending_modulus * 
                           (self.mean_curvature - self.spontaneous_curvature).powi(2);
        let gaussian_energy = self.gaussian_modulus * self.gaussian_curvature;
        
        bending_energy + gaussian_energy
    }
}

#[derive(Debug, Clone)]
pub struct CurvatureSensingProtein {
    protein_id: String,
    curvature_sensitivity: f64,   // Response magnitude per unit curvature
    optimal_curvature: f64,       // m⁻¹ - preferred curvature
    activity_modulation: ActivityModulation,
}

impl CurvatureSensingProtein {
    pub fn calculate_curvature_response(&self, membrane_curvature: f64) -> f64 {
        // Activity changes based on deviation from optimal curvature
        let curvature_deviation = (membrane_curvature - self.optimal_curvature).abs();
        let activity_factor = 1.0 + self.curvature_sensitivity * curvature_deviation;
        
        match self.activity_modulation {
            ActivityModulation::Increase => activity_factor,
            ActivityModulation::Decrease => 1.0 / activity_factor,
            ActivityModulation::Optimal => {
                // Bell-shaped response - maximum activity at optimal curvature
                (-0.5 * curvature_deviation.powi(2)).exp()
            }
        }
    }
}
```

## Electrochemical Gradients

### Ion Distribution and Membrane Potential

```rust
pub struct ElectrochemicalGradients {
    // Ion concentrations (mM)
    intracellular_concentrations: HashMap<IonType, f64>,
    extracellular_concentrations: HashMap<IonType, f64>,
    
    // Membrane properties
    membrane_potential: f64,      // mV
    temperature: f64,             // K
    
    // ATP-dependent ion pumps
    atp_pumps: Vec<AtpIonPump>,
    
    // Passive leak channels
    leak_channels: HashMap<IonType, f64>, // Permeabilities in m/s
}

impl ElectrochemicalGradients {
    pub fn calculate_nernst_potential(&self, ion_type: IonType) -> f64 {
        // Nernst equation: E = (RT/zF) * ln([ion]_out / [ion]_in)
        let gas_constant = 8.314; // J/(mol⋅K)
        let faraday_constant = 96485.0; // C/mol
        let ion_charge = ion_type.charge();
        
        let conc_out = self.extracellular_concentrations[&ion_type];
        let conc_in = self.intracellular_concentrations[&ion_type];
        
        (gas_constant * self.temperature) / (ion_charge as f64 * faraday_constant) * 
        (conc_out / conc_in).ln() * 1000.0 // Convert to mV
    }
    
    pub fn calculate_goldman_potential(&self) -> f64 {
        // Goldman-Hodgkin-Katz equation for membrane potential
        let pk = self.leak_channels[&IonType::Potassium];
        let pna = self.leak_channels[&IonType::Sodium];
        let pcl = self.leak_channels[&IonType::Chloride];
        
        let k_out = self.extracellular_concentrations[&IonType::Potassium];
        let k_in = self.intracellular_concentrations[&IonType::Potassium];
        let na_out = self.extracellular_concentrations[&IonType::Sodium];
        let na_in = self.intracellular_concentrations[&IonType::Sodium];
        let cl_out = self.extracellular_concentrations[&IonType::Chloride];
        let cl_in = self.intracellular_concentrations[&IonType::Chloride];
        
        let numerator = pk * k_out + pna * na_out + pcl * cl_in; // Note: Cl inside for anions
        let denominator = pk * k_in + pna * na_in + pcl * cl_out;
        
        let rt_f = 8.314 * self.temperature / 96485.0;
        rt_f * (numerator / denominator).ln() * 1000.0 // Convert to mV
    }
    
    pub fn update_for_atp_pump_activity(&mut self, atp_consumed: f64, dt: f64) {
        // Update ion concentrations based on ATP pump activity
        for pump in &mut self.atp_pumps {
            let pump_activity = pump.calculate_activity(atp_consumed, self.membrane_potential);
            let ion_flux = pump.calculate_ion_flux(pump_activity, dt);
            
            // Update concentrations
            for (ion_type, flux) in ion_flux {
                self.intracellular_concentrations.entry(ion_type)
                    .and_modify(|conc| *conc += flux.intracellular_change);
                self.extracellular_concentrations.entry(ion_type)
                    .and_modify(|conc| *conc += flux.extracellular_change);
            }
        }
        
        // Recalculate membrane potential
        self.membrane_potential = self.calculate_goldman_potential();
    }
    
    pub fn calculate_driving_force(&self, ion_type: IonType) -> f64 {
        // Electrochemical driving force: V_m - E_ion
        let nernst_potential = self.calculate_nernst_potential(ion_type);
        self.membrane_potential - nernst_potential
    }
    
    pub fn calculate_ion_current(&self, ion_type: IonType, conductance: f64) -> f64 {
        // Ohmic current: I = g * (V_m - E_ion)
        let driving_force = self.calculate_driving_force(ion_type);
        conductance * driving_force / 1000.0 // Convert mV to V
    }
}

#[derive(Debug, Clone)]
pub struct AtpIonPump {
    pump_type: PumpType,
    max_turnover_rate: f64,       // s⁻¹
    atp_stoichiometry: f64,       // ATP consumed per cycle
    ion_stoichiometry: HashMap<IonType, i8>, // Ions moved per cycle (+ out, - in)
    voltage_dependence: f64,      // Voltage dependence factor
    atp_km: f64,                  // mM - Michaelis constant for ATP
}

impl AtpIonPump {
    pub fn calculate_activity(&self, atp_concentration: f64, membrane_voltage: f64) -> f64 {
        // Michaelis-Menten kinetics for ATP dependence
        let atp_factor = atp_concentration / (self.atp_km + atp_concentration);
        
        // Voltage dependence (some pumps are electrogenic)
        let voltage_factor = if self.voltage_dependence != 0.0 {
            (-self.voltage_dependence * membrane_voltage / 1000.0).exp() // Convert mV to V
        } else {
            1.0
        };
        
        self.max_turnover_rate * atp_factor * voltage_factor
    }
    
    pub fn calculate_atp_consumption_rate(&self, activity: f64) -> f64 {
        activity * self.atp_stoichiometry
    }
}
```

## Circuit Parameter Integration

### Dynamic Circuit Parameter Calculation

```rust
pub struct MolecularToCircuitMapper {
    // Layer components
    bilayer_physics: LipidBilayerCapacitance,
    protein_electrics: HashMap<String, TransmembraneProteinElectrics>,
    curvature_dynamics: MembraneCurvatureDynamics,
    gradients: ElectrochemicalGradients,
    
    // Mapping parameters
    membrane_area: f64,           // m² - total membrane area
    temperature: f64,             // K
    
    // ATP coupling
    current_atp_concentration: f64, // mM
}

impl MolecularToCircuitMapper {
    pub fn calculate_circuit_parameters(&self) -> CircuitParameters {
        // Membrane capacitance
        let specific_capacitance = self.bilayer_physics.calculate_specific_capacitance();
        let curvature_effects = self.curvature_dynamics.calculate_curvature_electrical_effects();
        let total_capacitance = specific_capacitance * self.membrane_area * 
                               curvature_effects.capacitance_modification;
        
        // Membrane resistance (parallel combination of all ion pathways)
        let total_conductance = self.calculate_total_membrane_conductance();
        let total_resistance = 1.0 / total_conductance;
        
        // Protein circuit elements
        let protein_circuits = self.calculate_protein_circuit_elements();
        
        // Electrochemical potentials
        let equilibrium_potentials = self.calculate_equilibrium_potentials();
        
        CircuitParameters {
            membrane_capacitance: total_capacitance,
            membrane_resistance: total_resistance,
            protein_elements: protein_circuits,
            equilibrium_potentials,
            membrane_potential: self.gradients.membrane_potential,
            atp_dependencies: self.calculate_atp_dependencies(),
        }
    }
    
    fn calculate_total_membrane_conductance(&self) -> f64 {
        let mut total_conductance = 0.0;
        
        // Passive leak conductances
        for (ion_type, permeability) in &self.gradients.leak_channels {
            let ion_conductance = self.permeability_to_conductance(*permeability, *ion_type);
            total_conductance += ion_conductance;
        }
        
        // Protein conductances
        for (protein_id, protein) in &self.protein_electrics {
            let protein_conductance = protein.calculate_atp_dependent_conductance(
                self.current_atp_concentration
            );
            total_conductance += protein_conductance;
        }
        
        // Curvature effects
        let curvature_effects = self.curvature_dynamics.calculate_curvature_electrical_effects();
        total_conductance *= curvature_effects.resistance_modification;
        
        total_conductance
    }
    
    fn permeability_to_conductance(&self, permeability: f64, ion_type: IonType) -> f64 {
        // Convert permeability (m/s) to conductance (S)
        let faraday_constant = 96485.0; // C/mol
        let gas_constant = 8.314; // J/(mol⋅K)
        let ion_charge = ion_type.charge();
        
        let ion_concentration = (self.gradients.intracellular_concentrations[&ion_type] + 
                               self.gradients.extracellular_concentrations[&ion_type]) / 2.0; // Average
        
        permeability * self.membrane_area * ion_concentration * 
        faraday_constant.powi(2) / (gas_constant * self.temperature) * 
        ion_charge.powi(2)
    }
    
    pub fn update_for_atp_changes(&mut self, atp_consumption: f64, new_atp_concentration: f64, dt: f64) {
        // Update all molecular components for ATP changes
        self.bilayer_physics.update_for_atp_consumption(atp_consumption, dt);
        
        for protein in self.protein_electrics.values_mut() {
            protein.update_for_atp_changes(atp_consumption, new_atp_concentration);
        }
        
        self.curvature_dynamics.update_curvature_for_atp_consumption(atp_consumption, dt);
        self.gradients.update_for_atp_pump_activity(atp_consumption, dt);
        
        self.current_atp_concentration = new_atp_concentration;
    }
}
```

## Example: ATP Synthase Circuit Modeling

```rust
// Example: Detailed ATP synthase circuit model
pub fn create_atp_synthase_circuit() -> TransmembraneProteinElectrics {
    let mut atp_synthase = TransmembraneProteinElectrics::new("atp_synthase_f0f1");
    
    // F0 motor properties
    atp_synthase.set_f0_properties(F0Properties {
        c_ring_stoichiometry: 8,        // 8 c-subunits
        proton_binding_sites: 8,        // One per c-subunit
        rotational_resistance: 1e-21,   // Nm⋅s - rotational friction
        torque_per_proton: 40e-21,      // Nm - torque generated per proton
    });
    
    // F1 motor properties  
    atp_synthase.set_f1_properties(F1Properties {
        catalytic_sites: 3,             // 3 β subunits
        atp_binding_cooperativity: 2.0, // Hill coefficient
        atp_synthesis_rate: 100.0,      // s⁻¹ maximum rate
        atp_km: 0.5,                    // mM
    });
    
    // Coupling between F0 and F1
    atp_synthase.set_coupling(F0F1Coupling {
        gear_ratio: 8.0 / 3.0,          // 8 protons per 3 ATP
        coupling_efficiency: 0.95,       // 95% mechanical coupling
        slip_rate: 0.01,                // 1% slip under load
    });
    
    atp_synthase
}

// Circuit parameter calculation for ATP synthase
impl TransmembraneProteinElectrics {
    pub fn calculate_atp_synthase_circuit(&self, proton_gradient: f64, atp_concentration: f64) -> ProteinEquivalentCircuit {
        let f0_properties = self.f0_properties.as_ref().unwrap();
        let f1_properties = self.f1_properties.as_ref().unwrap();
        let coupling = self.coupling.as_ref().unwrap();
        
        // Calculate driving torque from proton gradient
        let proton_chemical_potential = 2.3 * 8.314 * 310.0 * proton_gradient.log10(); // J/mol
        let driving_torque = f0_properties.torque_per_proton * f0_properties.proton_binding_sites;
        
        // Calculate load torque from ATP synthesis
        let atp_load = f1_properties.calculate_atp_synthesis_load(atp_concentration);
        
        // Net torque and rotation rate
        let net_torque = driving_torque - atp_load / coupling.gear_ratio;
        let rotation_rate = net_torque / f0_properties.rotational_resistance;
        
        // ATP synthesis rate
        let atp_synthesis_rate = rotation_rate * coupling.coupling_efficiency / coupling.gear_ratio;
        
        // Equivalent circuit
        ProteinEquivalentCircuit::AtpSynthase {
            proton_current: rotation_rate * f0_properties.c_ring_stoichiometry,
            atp_production_rate: atp_synthesis_rate,
            internal_resistance: self.calculate_synthase_internal_resistance(),
            torque_constant: driving_torque / proton_gradient,
            back_emf: atp_load,
        }
    }
}
```

---

The **Molecular Layer** provides the fundamental biophysical foundation for all membrane electrical properties, translating lipid-protein interactions, membrane curvature, and electrochemical gradients into the circuit parameters that enable Nebuchadnezzar's ATP-based differential equation modeling of biological systems. 