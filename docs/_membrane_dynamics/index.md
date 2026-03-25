---
layout: default
title: "Membrane Dynamics"
permalink: /membrane-dynamics/
description: "Biological quantum computing through membrane dynamics, Maxwell's demons, and circuit interfaces"
---

# Membrane Dynamics: Biological Quantum Computing

The **Membrane Dynamics** module provides the biophysical foundation for biological quantum computing in the Bene Gesserit framework. By treating cellular membranes as dynamic quantum circuits, we enable room-temperature quantum computation through Environment-Assisted Quantum Transport (ENAQT) while maintaining ATP-constrained authenticity.

---

## Revolutionary Membrane Computing

### Membranes as Quantum Circuits

Traditional biology views membranes as passive barriers. **Bene Gesserit treats membranes as active quantum computing elements:**

```rust
// Traditional membrane model
struct PassiveMembrane {
    lipid_composition: Vec<Lipid>,
    proteins: Vec<Protein>,
    permeability: f64,
}

// Bene Gesserit quantum membrane
struct QuantumMembrane {
    // Quantum computing elements
    quantum_states: Vec<MembraneQuantumState>,
    coherence_time: f64,
    entanglement_network: EntanglementGraph,
    
    // ENAQT properties
    environment_coupling: f64,
    transport_enhancement: f64,
    optimal_frequencies: Vec<f64>,
    
    // ATP integration
    atp_dependent_gates: Vec<AtpQuantumGate>,
    energy_landscapes: HashMap<String, EnergyLandscape>,
    
    // Circuit representation
    equivalent_circuit: QuantumCircuitTopology,
    dynamic_parameters: AtpDependentParameters,
}
```

### Environment-Assisted Quantum Transport (ENAQT)

The key insight is that **biological environments enhance rather than destroy quantum coherence:**

<div class="enaqt-visualization">
  <div class="enaqt-panel">
    <h4>🌊 Traditional Decoherence</h4>
    <div class="decoherence-graph">
      <p><strong>Coherence ∝ exp(-γt)</strong></p>
      <p>Environmental coupling destroys quantum effects</p>
      <ul>
        <li>Exponential decay</li>
        <li>Short coherence times</li>
        <li>No quantum advantage</li>
      </ul>
    </div>
  </div>
  
  <div class="enaqt-panel">
    <h4>⚡ ENAQT Enhancement</h4>
    <div class="enhancement-graph">
      <p><strong>Transport ∝ (1 + γ_optimal)</strong></p>
      <p>Optimal coupling enhances quantum transport</p>
      <ul>
        <li>Transport enhancement</li>
        <li>Robust quantum effects</li>
        <li>Biological temperatures</li>
      </ul>
    </div>
  </div>
</div>

#### ENAQT Implementation

```rust
pub struct EnaqtTransportSystem {
    pub system_hamiltonian: Array2<Complex<f64>>,
    pub environment_spectral_density: SpectralDensity,
    pub coupling_strength: f64,
    pub temperature: f64,
    pub transport_sites: Vec<TransportSite>,
}

impl EnaqtTransportSystem {
    pub fn calculate_transport_efficiency(&self) -> f64 {
        // Calculate bare (uncoupled) transport
        let bare_efficiency = self.calculate_bare_transport();
        
        // Calculate environmental enhancement
        let enhancement_factor = self.calculate_enhancement_factor();
        
        // ENAQT gives multiplicative enhancement
        bare_efficiency * (1.0 + enhancement_factor)
    }
    
    pub fn calculate_enhancement_factor(&self) -> f64 {
        let optimal_coupling = self.find_optimal_coupling_strength();
        let coupling_ratio = self.coupling_strength / optimal_coupling;
        
        // Enhancement peaks at optimal coupling
        match coupling_ratio {
            r if r <= 1.0 => r, // Linear increase to optimum
            r => 1.0 / r,       // Decrease beyond optimum
        }
    }
    
    pub fn evolve_with_environment(&mut self, dt: f64) -> QuantumState {
        // Construct total Hamiltonian including environment
        let h_total = self.construct_total_hamiltonian();
        
        // Time evolution with environmental assistance
        let evolution_operator = (-Complex::i() * h_total * dt).mapv(|x| x.exp());
        
        // Apply evolution to quantum state
        let new_state = evolution_operator.dot(&self.current_state);
        self.current_state = new_state.clone();
        
        // Return enhanced quantum state
        QuantumState {
            amplitudes: new_state,
            coherence_measure: self.calculate_coherence(),
            transport_fidelity: self.calculate_transport_fidelity(),
            enhancement_factor: self.calculate_enhancement_factor(),
        }
    }
}
```

---

## Biological Maxwell's Demons

**Information catalysts** that create order through pattern selection while maintaining thermodynamic consistency.

### The iCat Framework

```
iCat = ℑ_input ∘ ℑ_output
```

Where:
- `ℑ_input`: Input pattern recognition filter
- `ℑ_output`: Output channel targeting system  
- `∘`: Composition operator (information flow)

#### BMD Implementation

```rust
pub trait BiologicalMaxwellsDemon {
    fn recognize_patterns(&mut self, input: &[f64]) -> PatternRecognitionResult;
    fn select_output_channel(&self, pattern: &Pattern) -> OutputChannel;
    fn catalyze_information(&mut self, input: &[f64]) -> CatalysisResult;
    fn measure_efficiency(&self) -> f64;
    fn track_degradation(&mut self, time_step: f64);
    fn ensure_thermodynamic_consistency(&self) -> bool;
}

pub struct MembraneMaxwellsDemon {
    // Pattern recognition system
    pub pattern_memory: PatternRecognitionMemory<MembranePattern>,
    pub recognition_threshold: f64,
    pub learning_rate: f64,
    
    // Information catalysis
    pub catalysis_efficiency: f64,
    pub energy_cost_per_bit: f64,
    pub information_throughput: f64,
    
    // Thermodynamic compliance
    pub entropy_production_rate: f64,
    pub heat_dissipation: f64,
    pub haldane_relation_checker: HaldaneRelationChecker,
    
    // ATP coupling
    pub atp_consumption_rate: f64,
    pub atp_efficiency: f64,
}

impl BiologicalMaxwellsDemon for MembraneMaxwellsDemon {
    fn catalyze_information(&mut self, membrane_state: &[f64]) -> CatalysisResult {
        // 1. Recognize patterns in membrane state
        let pattern_result = self.recognize_patterns(membrane_state);
        
        // 2. Select appropriate output channel
        let output_channel = self.select_output_channel(&pattern_result.pattern());
        
        // 3. Ensure thermodynamic consistency
        if !self.ensure_thermodynamic_consistency() {
            return CatalysisResult::ThermodynamicViolation;
        }
        
        // 4. Calculate ATP cost
        let atp_cost = self.calculate_atp_cost(&output_channel);
        
        // 5. Perform information catalysis
        let catalysis_efficiency = self.measure_efficiency();
        
        // 6. Update system state
        self.update_catalysis_metrics(catalysis_efficiency);
        
        CatalysisResult::Success {
            input_pattern: pattern_result.pattern(),
            output_channel,
            efficiency: catalysis_efficiency,
            atp_cost,
            entropy_produced: self.entropy_production_rate,
            information_processed: output_channel.information_content(),
        }
    }
    
    fn ensure_thermodynamic_consistency(&self) -> bool {
        // Check Haldane relation for detailed balance
        self.haldane_relation_checker.verify_detailed_balance() &&
        // Ensure positive entropy production
        self.entropy_production_rate > 0.0 &&
        // Verify energy conservation
        self.verify_energy_conservation()
    }
}
```

### Pattern Recognition in Membranes

```rust
#[derive(Debug, Clone)]
pub struct MembranePattern {
    pub lipid_composition_signature: Vec<f64>,
    pub protein_distribution_pattern: Vec<f64>,
    pub curvature_profile: Vec<f64>,
    pub electrical_potential_landscape: Vec<f64>,
    pub oscillation_frequencies: Vec<f64>,
    pub atp_concentration_gradients: Vec<f64>,
}

impl MembranePattern {
    pub fn extract_from_membrane(membrane: &QuantumMembrane) -> Self {
        MembranePattern {
            lipid_composition_signature: membrane.analyze_lipid_signature(),
            protein_distribution_pattern: membrane.analyze_protein_clustering(),
            curvature_profile: membrane.calculate_curvature_profile(),
            electrical_potential_landscape: membrane.calculate_potential_landscape(),
            oscillation_frequencies: membrane.extract_oscillation_spectrum(),
            atp_concentration_gradients: membrane.measure_atp_gradients(),
        }
    }
    
    pub fn similarity(&self, other: &MembranePattern) -> f64 {
        // Multi-dimensional pattern similarity
        let lipid_sim = cosine_similarity(&self.lipid_composition_signature, 
                                         &other.lipid_composition_signature);
        let protein_sim = cosine_similarity(&self.protein_distribution_pattern,
                                           &other.protein_distribution_pattern);
        let curvature_sim = cosine_similarity(&self.curvature_profile,
                                             &other.curvature_profile);
        let potential_sim = cosine_similarity(&self.electrical_potential_landscape,
                                             &other.electrical_potential_landscape);
        let oscillation_sim = cosine_similarity(&self.oscillation_frequencies,
                                               &other.oscillation_frequencies);
        let atp_sim = cosine_similarity(&self.atp_concentration_gradients,
                                       &other.atp_concentration_gradients);
        
        // Weighted average of similarities
        (lipid_sim * 0.2 + protein_sim * 0.2 + curvature_sim * 0.15 + 
         potential_sim * 0.15 + oscillation_sim * 0.15 + atp_sim * 0.15)
    }
}
```

---

## ATP-Constrained Membrane Dynamics

### dx/dATP Formulation

Moving from time-based to energy-based dynamics:

```rust
pub struct AtpConstrainedMembraneDynamics {
    pub membrane_state: MembraneState,
    pub atp_pool: AtpPool,
    pub atp_consumption_rates: HashMap<ProcessType, f64>,
    pub energy_landscapes: HashMap<ProcessType, EnergyLandscape>,
}

impl AtpConstrainedMembraneDynamics {
    pub fn calculate_membrane_derivatives(&self, atp_consumption: f64) -> MembraneDerivatives {
        let mut derivatives = MembraneDerivatives::new();
        
        // Lipid synthesis: dLipid/dATP
        derivatives.lipid_concentration = 
            self.atp_consumption_rates[&ProcessType::LipidSynthesis] * atp_consumption;
        
        // Protein insertion: dProtein/dATP
        derivatives.protein_density = 
            self.atp_consumption_rates[&ProcessType::ProteinInsertion] * atp_consumption;
        
        // Membrane curvature: dCurvature/dATP
        derivatives.membrane_curvature = 
            self.calculate_curvature_atp_dependence(atp_consumption);
        
        // Quantum state evolution: dQuantumState/dATP
        derivatives.quantum_state_amplitudes = 
            self.calculate_quantum_atp_coupling(atp_consumption);
        
        // Maxwell's demon activity: dDemonEfficiency/dATP
        derivatives.demon_efficiency = 
            self.calculate_demon_atp_dependence(atp_consumption);
        
        derivatives
    }
    
    pub fn integrate_atp_step(&mut self, atp_budget: f64) -> Result<IntegrationResult, AtpError> {
        // Check ATP availability
        if !self.atp_pool.can_consume(atp_budget) {
            return Err(AtpError::InsufficientAtp);
        }
        
        // Calculate derivatives
        let derivatives = self.calculate_membrane_derivatives(atp_budget);
        
        // Update membrane state
        self.membrane_state.apply_derivatives(&derivatives, atp_budget);
        
        // Consume ATP
        self.atp_pool.consume(atp_budget);
        
        // Update quantum coherence based on new membrane state
        let coherence_update = self.update_quantum_coherence();
        
        // Update Maxwell's demon efficiency
        let demon_update = self.update_demon_efficiency();
        
        Ok(IntegrationResult {
            membrane_state: self.membrane_state.clone(),
            atp_consumed: atp_budget,
            quantum_coherence: coherence_update.coherence,
            demon_efficiency: demon_update.efficiency,
            entropy_production: self.calculate_entropy_production(),
        })
    }
}
```

### Oscillatory Membrane Dynamics

```rust
pub struct OscillatoryMembrane {
    pub oscillators: Vec<MembraneOscillator>,
    pub coupling_matrix: Array2<f64>,
    pub endpoint_distributions: HashMap<String, EndpointDistribution>,
    pub atp_modulated_frequencies: HashMap<String, f64>,
}

impl OscillatoryMembrane {
    pub fn evolve_oscillations(&mut self, atp_consumption: f64) -> OscillationResult {
        let mut endpoint_updates = HashMap::new();
        
        for oscillator in &mut self.oscillators {
            // ATP modulates oscillation frequency
            let atp_modulated_freq = self.atp_modulated_frequencies[&oscillator.name] * 
                                   (1.0 + 0.1 * atp_consumption);
            
            // Evolve oscillator
            let oscillation_result = oscillator.evolve_with_frequency(atp_modulated_freq);
            
            // Check for endpoint reached
            if oscillation_result.endpoint_reached {
                let endpoint = OscillationEndpoint {
                    oscillator_name: oscillator.name.clone(),
                    position: oscillation_result.final_position,
                    velocity: oscillation_result.final_velocity,
                    energy: oscillation_result.final_energy,
                    atp_consumed: atp_consumption,
                    probability: oscillation_result.endpoint_probability,
                };
                
                // Update endpoint distribution
                self.endpoint_distributions
                    .get_mut(&oscillator.name)
                    .unwrap()
                    .add_endpoint(endpoint);
            }
        }
        
        // Calculate total entropy from endpoint distributions
        let total_entropy = self.calculate_oscillatory_entropy();
        
        OscillationResult {
            oscillation_endpoints: endpoint_updates,
            total_entropy,
            atp_efficiency: self.calculate_atp_efficiency(),
        }
    }
    
    pub fn calculate_oscillatory_entropy(&self) -> f64 {
        self.endpoint_distributions
            .values()
            .map(|dist| dist.calculate_entropy())
            .sum()
    }
}
```

---

## Circuit Interface Layer

### Membrane-to-Circuit Mapping

```rust
pub struct MembraneCircuitMapper {
    pub capacitance_mapper: CapacitanceMapper,
    pub resistance_mapper: ResistanceMapper,
    pub inductance_mapper: InductanceMapper,
    pub quantum_gate_mapper: QuantumGateMapper,
}

impl MembraneCircuitMapper {
    pub fn map_membrane_to_quantum_circuit(&self, membrane: &QuantumMembrane) -> QuantumCircuit {
        // Map membrane patches to circuit elements
        let mut circuit_elements = Vec::new();
        
        // Lipid bilayer as quantum capacitor
        for patch in &membrane.membrane_patches {
            let capacitance = self.capacitance_mapper.calculate_quantum_capacitance(patch);
            let quantum_capacitor = QuantumCircuitElement::QuantumCapacitor {
                capacitance,
                quantum_state: patch.quantum_state.clone(),
                coherence_time: patch.coherence_time,
            };
            circuit_elements.push(quantum_capacitor);
        }
        
        // Membrane proteins as quantum gates
        for protein in &membrane.membrane_proteins {
            let quantum_gate = self.quantum_gate_mapper.protein_to_quantum_gate(protein);
            circuit_elements.push(QuantumCircuitElement::QuantumGate(quantum_gate));
        }
        
        // Ion channels as variable resistors
        for channel in &membrane.ion_channels {
            let resistance = self.resistance_mapper.calculate_quantum_resistance(channel);
            let quantum_resistor = QuantumCircuitElement::QuantumResistor {
                resistance,
                quantum_conductance: channel.quantum_conductance(),
                atp_dependence: channel.atp_dependence.clone(),
            };
            circuit_elements.push(quantum_resistor);
        }
        
        // ATP pumps as quantum batteries
        for pump in &membrane.atp_pumps {
            let quantum_battery = QuantumCircuitElement::QuantumBattery {
                voltage: pump.electrochemical_potential(),
                quantum_efficiency: pump.quantum_efficiency(),
                atp_consumption_rate: pump.atp_consumption_rate,
            };
            circuit_elements.push(quantum_battery);
        }
        
        QuantumCircuit {
            elements: circuit_elements,
            topology: self.determine_circuit_topology(membrane),
            quantum_entanglement: membrane.entanglement_network.clone(),
            atp_coupling: self.create_atp_coupling_matrix(membrane),
        }
    }
}
```

---

## Interactive Membrane Visualization

<div class="membrane-demo">
  <div class="demo-controls">
    <button class="control-button active" onclick="showMembraneDemo('enaqt')">ENAQT Transport</button>
    <button class="control-button" onclick="showMembraneDemo('maxwell')">Maxwell's Demons</button>
    <button class="control-button" onclick="showMembraneDemo('oscillations')">Oscillatory Dynamics</button>
    <button class="control-button" onclick="showMembraneDemo('circuits')">Quantum Circuits</button>
  </div>

  <div id="enaqt-demo" class="membrane-content active">
    <h4>⚡ Environment-Assisted Quantum Transport</h4>
    <div class="demo-visualization">
      <div class="transport-pathway">
        <div class="quantum-site">Site A</div>
        <div class="transport-arrow">→</div>
        <div class="environment-coupling">Environment<br/>Coupling</div>
        <div class="transport-arrow">→</div>
        <div class="quantum-site">Site B</div>
      </div>
      <div class="transport-equations">
        <p><strong>Enhancement Factor:</strong></p>
        <code>η = 1 + (γ/γ_optimal) for γ ≤ γ_optimal</code>
        <p><strong>Transport Efficiency:</strong></p>
        <code>T_eff = T_bare × η</code>
      </div>
    </div>
    <pre><code class="language-rust">// ENAQT implementation example
let mut enaqt_system = EnaqtTransportSystem::new()
    .with_coupling_strength(0.5)  // Optimal coupling
    .with_temperature(310.0)      // Biological temperature
    .with_transport_sites(vec![
        TransportSite::new("donor", energy: 0.0),
        TransportSite::new("acceptor", energy: -0.2),
    ]);

let transport_result = enaqt_system.calculate_transport_efficiency();
println!("Transport enhancement: {:.2}x", transport_result.enhancement_factor);</code></pre>
  </div>

  <div id="maxwell-demo" class="membrane-content">
    <h4>🧠 Biological Maxwell's Demons</h4>
    <div class="demo-visualization">
      <div class="demon-process">
        <div class="demon-input">Input<br/>Patterns</div>
        <div class="demon-arrow">→</div>
        <div class="demon-brain">Pattern<br/>Recognition</div>
        <div class="demon-arrow">→</div>
        <div class="demon-output">Output<br/>Channels</div>
      </div>
      <div class="demon-equation">
        <p><strong>Information Catalysis:</strong></p>
        <code>iCat = ℑ_input ∘ ℑ_output</code>
        <p><strong>Thermodynamic Cost:</strong></p>
        <code>ΔS_universe ≥ k ln(2) per bit</code>
      </div>
    </div>
    <pre><code class="language-rust">// Maxwell's demon implementation
let mut membrane_demon = MembraneMaxwellsDemon::new()
    .with_pattern_memory_size(1000)
    .with_recognition_threshold(0.8)
    .with_atp_coupling(true);

let catalysis_result = membrane_demon.catalyze_information(&membrane_state);
match catalysis_result {
    CatalysisResult::Success { efficiency, atp_cost, .. } => {
        println!("Information catalyzed with {:.1}% efficiency", efficiency * 100.0);
        println!("ATP cost: {:.3} mM", atp_cost);
    }
    _ => println!("Catalysis failed"),
}</code></pre>
  </div>

  <div id="oscillations-demo" class="membrane-content">
    <h4>🌊 Oscillatory Dynamics</h4>
    <div class="demo-visualization">
      <div class="oscillation-timeline">
        <div class="oscillation-phase">Initiation</div>
        <div class="oscillation-arrow">→</div>
        <div class="oscillation-phase">Energy<br/>Dissipation</div>
        <div class="oscillation-arrow">→</div>
        <div class="oscillation-phase">Endpoint<br/>Selection</div>
      </div>
      <div class="entropy-calculation">
        <p><strong>Oscillatory Entropy:</strong></p>
        <code>S = k ln Ω_endpoints</code>
        <p><strong>ATP Coupling:</strong></p>
        <code>f(ATP) = f_0 × (1 + α × [ATP])</code>
      </div>
    </div>
    <pre><code class="language-rust">// Oscillatory membrane dynamics
let mut oscillatory_membrane = OscillatoryMembrane::new()
    .with_oscillators(vec![
        MembraneOscillator::new("lipid_waves", frequency: 1.0),
        MembraneOscillator::new("protein_fluctuations", frequency: 10.0),
        MembraneOscillator::new("atp_oscillations", frequency: 0.1),
    ])
    .with_atp_modulation(true);

let oscillation_result = oscillatory_membrane.evolve_oscillations(atp_budget);
println!("Oscillatory entropy: {:.3} k_B", oscillation_result.total_entropy);</code></pre>
  </div>

  <div id="circuits-demo" class="membrane-content">
    <h4>🔌 Quantum Circuit Representation</h4>
    <div class="demo-visualization">
      <div class="circuit-diagram">
        <div class="circuit-element">Lipid<br/>Capacitor</div>
        <div class="circuit-connection">—</div>
        <div class="circuit-element">Protein<br/>Gate</div>
        <div class="circuit-connection">—</div>
        <div class="circuit-element">Ion Channel<br/>Resistor</div>
      </div>
      <div class="circuit-equations">
        <p><strong>Quantum Capacitance:</strong></p>
        <code>C_q = e²(∂n/∂μ)</code>
        <p><strong>Quantum Conductance:</strong></p>
        <code>G_q = (2e²/h) × T</code>
      </div>
    </div>
    <pre><code class="language-rust">// Membrane-to-circuit mapping
let circuit_mapper = MembraneCircuitMapper::new();
let quantum_circuit = circuit_mapper.map_membrane_to_quantum_circuit(&membrane);

// Execute quantum computation
let computation_result = quantum_circuit.execute_computation(
    &quantum_algorithm,
    atp_budget: 100.0,
    coherence_threshold: 0.8,
);

println!("Quantum fidelity: {:.3}", computation_result.fidelity);</code></pre>
  </div>
</div>

---

## Applications

### 🧬 Protein Folding Prediction
Use membrane quantum computation to predict protein folding pathways with unprecedented accuracy.

### 💊 Drug Discovery
Simulate drug-membrane interactions using quantum-enhanced molecular dynamics.

### 🧠 Neural Computation
Model synaptic transmission as quantum information processing in membrane circuits.

### 🔋 Bioenergetics
Optimize ATP synthesis and consumption through quantum-enhanced metabolic modeling.

---

## Next Steps

<div class="membrane-nav">
  <a href="/membrane-dynamics/biological-maxwell-demons/" class="nav-card">
    <h4>🧠 Biological Maxwell's Demons</h4>
    <p>Information catalysts, pattern recognition, and thermodynamic consistency</p>
  </a>
  
  <a href="/membrane-dynamics/circuit-interface-layer/" class="nav-card">
    <h4>🔌 Circuit Interface Layer</h4>
    <p>Membrane-to-circuit mapping, quantum gates, and ATP coupling</p>
  </a>
  
  <a href="/membrane-dynamics/quickstart-example/" class="nav-card">
    <h4>🚀 Quickstart Example</h4>
    <p>Complete membrane quantum computation tutorial</p>
  </a>
  
  <a href="/examples/" class="nav-card">
    <h4>🔬 Practical Examples</h4>
    <p>Real-world membrane computing applications</p>
  </a>
</div>

---

<style>
.enaqt-visualization {
  display: grid;
  grid-template-columns: 1fr 1fr;
  gap: 2rem;
  margin: 2rem 0;
}

.enaqt-panel {
  background: white;
  padding: 1.5rem;
  border-radius: 0.5rem;
  box-shadow: 0 2px 10px rgba(0,0,0,0.1);
}

.enaqt-panel h4 {
  margin-bottom: 1rem;
  text-align: center;
}

.decoherence-graph, .enhancement-graph {
  text-align: center;
}

.decoherence-graph {
  border-left: 3px solid #dc3545;
}

.enhancement-graph {
  border-left: 3px solid #28a745;
}

.membrane-demo {
  background: #f8f9fa;
  border-radius: 0.5rem;
  padding: 1.5rem;
  margin: 2rem 0;
}

.demo-controls {
  display: flex;
  gap: 0.5rem;
  margin-bottom: 1rem;
  flex-wrap: wrap;
}

.control-button {
  padding: 0.5rem 1rem;
  background: #e9ecef;
  border: none;
  border-radius: 0.25rem;
  cursor: pointer;
  transition: all 0.3s ease;
}

.control-button.active {
  background: #007bff;
  color: white;
}

.membrane-content {
  display: none;
}

.membrane-content.active {
  display: block;
}

.demo-visualization {
  background: white;
  padding: 1rem;
  border-radius: 0.5rem;
  margin: 1rem 0;
  display: flex;
  flex-direction: column;
  align-items: center;
  gap: 1rem;
}

.transport-pathway, .demon-process, .oscillation-timeline, .circuit-diagram {
  display: flex;
  align-items: center;
  gap: 1rem;
  flex-wrap: wrap;
  justify-content: center;
}

.quantum-site, .demon-input, .demon-brain, .demon-output, 
.oscillation-phase, .circuit-element {
  background: #007bff;
  color: white;
  padding: 0.5rem 1rem;
  border-radius: 0.25rem;
  text-align: center;
  min-width: 80px;
}

.environment-coupling {
  background: #28a745;
  color: white;
  padding: 0.5rem 1rem;
  border-radius: 0.25rem;
  text-align: center;
}

.transport-arrow, .demon-arrow, .oscillation-arrow {
  font-size: 1.5rem;
  color: #666;
}

.circuit-connection {
  font-size: 1.5rem;
  color: #666;
  font-weight: bold;
}

.transport-equations, .demon-equation, .entropy-calculation, .circuit-equations {
  background: #f8f9fa;
  padding: 1rem;
  border-radius: 0.25rem;
  text-align: center;
}

.membrane-nav {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
  gap: 1.5rem;
  margin: 2rem 0;
}

.nav-card {
  background: white;
  padding: 1.5rem;
  border-radius: 0.5rem;
  box-shadow: 0 2px 8px rgba(0,0,0,0.1);
  text-decoration: none;
  color: inherit;
  transition: all 0.3s ease;
  border-left: 4px solid #6f42c1;
}

.nav-card:hover {
  transform: translateY(-3px);
  box-shadow: 0 4px 15px rgba(0,0,0,0.15);
}

.nav-card h4 {
  margin-bottom: 0.5rem;
  color: #6f42c1;
}

.nav-card p {
  color: #666;
  line-height: 1.6;
}

pre {
  background: #2d3748;
  color: #e2e8f0;
  padding: 1rem;
  border-radius: 0.5rem;
  overflow-x: auto;
  margin: 1rem 0;
}

code {
  font-family: 'Fira Code', 'Monaco', 'Consolas', monospace;
  font-size: 0.9em;
}
</style>

<script>
function showMembraneDemo(demoId) {
  // Hide all demos
  const demos = document.querySelectorAll('.membrane-content');
  demos.forEach(demo => demo.classList.remove('active'));
  
  // Hide all control buttons
  const buttons = document.querySelectorAll('.control-button');
  buttons.forEach(button => button.classList.remove('active'));
  
  // Show selected demo
  document.getElementById(demoId + '-demo').classList.add('active');
  
  // Activate corresponding button
  event.target.classList.add('active');
}
</script>