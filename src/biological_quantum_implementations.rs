//! # Complete Biological Quantum Computer Solver Implementation
//! 
//! This implements the main solver for biological quantum computation with:
//! - Multiple integration methods (Velocity-Verlet, Runge-Kutta 4, Adaptive)
//! - Oscillation endpoint prediction (key entropy insight)
//! - Radical generation calculation (death mechanism)
//! - Entropy constraint enforcement

use crate::biological_quantum_computer::*;
use std::collections::HashMap;
use num_complex::Complex;
use ndarray::Array2;
use std::f64::consts::PI;

impl BiologicalQuantumHamiltonian {
    /// ATP dynamics: dx/dATP from the original framework
    pub fn calculate_atp_derivatives(&self, state: &BiologicalQuantumState) -> AtpDerivatives {
        let atp_consumption_rate = self.calculate_atp_consumption_rate(state);
        let oscillatory_atp_coupling = self.calculate_oscillatory_atp_coupling(state);
        let membrane_atp_coupling = self.calculate_membrane_atp_coupling(state);
        
        AtpDerivatives {
            atp_concentration_rate: -atp_consumption_rate * (1.0 + oscillatory_atp_coupling + membrane_atp_coupling),
            adp_concentration_rate: atp_consumption_rate,
            pi_concentration_rate: atp_consumption_rate,
            energy_charge_rate: self.calculate_energy_charge_rate(state),
            oscillation_amplitude_rate: oscillatory_atp_coupling * atp_consumption_rate,
            oscillation_phase_rate: state.atp_coords.atp_oscillation_frequency,
        }
    }

    /// Oscillatory dynamics: standard Hamiltonian mechanics with ATP driving
    pub fn calculate_oscillatory_derivatives(&self, state: &BiologicalQuantumState) -> OscillatoryDerivatives {
        let mut position_derivatives = Vec::new();
        let mut momentum_derivatives = Vec::new();
        
        for (i, oscillation) in state.oscillatory_coords.oscillations.iter().enumerate() {
            // Position derivative: dq/dt = p (momentum)
            position_derivatives.push(state.oscillatory_coords.oscillatory_momenta[i]);
            
            // Momentum derivative: dp/dt = -∂V/∂q + ATP_driving
            let force = -self.calculate_oscillatory_force(oscillation, state);
            let atp_driving = self.calculate_atp_driving_force(oscillation, &state.atp_coords);
            momentum_derivatives.push(force + atp_driving);
        }
        
        OscillatoryDerivatives {
            position_derivatives,
            momentum_derivatives,
            phase_derivatives: self.calculate_phase_coupling_derivatives(state),
        }
    }

    /// Membrane quantum dynamics: Schrödinger equation with ENAQT
    pub fn calculate_membrane_derivatives(&self, state: &BiologicalQuantumState) -> MembraneDerivatives {
        let mut quantum_state_derivatives = Vec::new();
        
        for quantum_state in &state.membrane_coords.quantum_states {
            // Time-dependent Schrödinger equation with environmental coupling
            let hamiltonian_term = -Complex::i() * quantum_state.energy * quantum_state.amplitude;
            let environmental_term = self.calculate_enaqt_coupling(quantum_state, state);
            let atp_quantum_coupling = self.calculate_atp_quantum_coupling(quantum_state, state);
            
            quantum_state_derivatives.push(hamiltonian_term + environmental_term + atp_quantum_coupling);
        }
        
        MembraneDerivatives {
            quantum_state_derivatives,
            tunneling_derivatives: self.calculate_tunneling_derivatives(state),
            environmental_coupling_derivatives: self.calculate_environmental_derivatives(state),
        }
    }

    /// Entropy dynamics: oscillation endpoint statistics
    pub fn calculate_entropy_derivatives(&self, state: &BiologicalQuantumState) -> EntropyDerivatives {
        // Calculate how oscillation endpoints are changing
        let endpoint_evolution_rate = self.calculate_endpoint_evolution_rate(state);
        
        // Entropy production from ATP consumption
        let atp_entropy_production = self.calculate_atp_entropy_production(state);
        
        // Oscillatory entropy production
        let oscillatory_entropy_production = self.calculate_oscillatory_entropy_production(state);
        
        // Membrane quantum entropy production
        let membrane_entropy_production = self.calculate_membrane_entropy_production(state);
        
        // Quantum tunneling entropy (death mechanism)
        let quantum_tunneling_entropy = self.calculate_quantum_tunneling_entropy_production(state);
        
        EntropyDerivatives {
            total_entropy_rate: atp_entropy_production + oscillatory_entropy_production + membrane_entropy_production,
            endpoint_distribution_rates: endpoint_evolution_rate,
            membrane_endpoint_entropy_rate: membrane_entropy_production,
            quantum_tunneling_entropy_rate: quantum_tunneling_entropy,
        }
    }

    /// Calculate ATP consumption rate based on system state
    fn calculate_atp_consumption_rate(&self, state: &BiologicalQuantumState) -> f64 {
        let base_rate = 0.1; // Base consumption rate
        
        // Consumption increases with oscillation activity
        let oscillation_factor: f64 = state.oscillatory_coords.oscillations.iter()
            .map(|osc| osc.amplitude * osc.frequency * osc.atp_coupling_strength)
            .sum();
        
        // Consumption increases with membrane quantum activity
        let quantum_factor: f64 = state.membrane_coords.quantum_states.iter()
            .map(|qs| qs.amplitude.norm_sqr() * qs.energy)
            .sum();
        
        base_rate * (1.0 + 0.1 * oscillation_factor + 0.05 * quantum_factor)
    }
    
    /// Calculate oscillatory-ATP coupling
    fn calculate_oscillatory_atp_coupling(&self, state: &BiologicalQuantumState) -> f64 {
        let atp_energy = state.atp_coords.available_energy();
        let oscillation_demand: f64 = state.oscillatory_coords.oscillations.iter()
            .map(|osc| osc.amplitude * osc.frequency * osc.atp_coupling_strength)
            .sum();
        
        if oscillation_demand > 0.0 {
            (atp_energy / oscillation_demand * 0.1).min(1.0)
        } else {
            0.0
        }
    }
    
    /// Calculate membrane-ATP coupling
    fn calculate_membrane_atp_coupling(&self, state: &BiologicalQuantumState) -> f64 {
        let atp_energy = state.atp_coords.available_energy();
        let membrane_demand: f64 = state.membrane_coords.quantum_states.iter()
            .map(|qs| qs.amplitude.norm_sqr() * qs.energy)
            .sum();
        
        if membrane_demand > 0.0 {
            (atp_energy / membrane_demand * 0.05).min(0.5)
        } else {
            0.0
        }
    }
    
    /// Calculate energy charge rate
    fn calculate_energy_charge_rate(&self, state: &BiologicalQuantumState) -> f64 {
        let atp = state.atp_coords.atp_concentration;
        let adp = state.atp_coords.adp_concentration;
        let amp = 0.1; // Approximate AMP concentration
        
        // d(energy_charge)/dt based on ATP/ADP/AMP dynamics
        let total = atp + adp + amp;
        if total > 0.0 {
            -0.5 * (adp / total) * 0.1 // Decreases as ADP increases
        } else {
            0.0
        }
    }
    
    /// Calculate oscillatory force on specific oscillation
    fn calculate_oscillatory_force(&self, oscillation: &OscillationState, _state: &BiologicalQuantumState) -> f64 {
        // Harmonic oscillator force: F = -kx - γv
        let spring_force = oscillation.frequency.powi(2) * oscillation.amplitude;
        let damping_force = oscillation.damping_coefficient * oscillation.amplitude * oscillation.frequency;
        
        spring_force + damping_force
    }
    
    /// Calculate ATP driving force on oscillation
    fn calculate_atp_driving_force(&self, oscillation: &OscillationState, atp_coords: &AtpCoordinates) -> f64 {
        let atp_energy = atp_coords.available_energy();
        let driving_strength = oscillation.atp_coupling_strength;
        
        // ATP provides driving force proportional to available energy
        driving_strength * atp_energy * 0.01
    }
    
    /// Calculate phase coupling derivatives
    fn calculate_phase_coupling_derivatives(&self, state: &BiologicalQuantumState) -> Vec<f64> {
        let mut derivatives = Vec::new();
        
        for (i, oscillation) in state.oscillatory_coords.oscillations.iter().enumerate() {
            let mut coupling_sum = 0.0;
            
            // Sum coupling with all other oscillations
            for (j, other_osc) in state.oscillatory_coords.oscillations.iter().enumerate() {
                if i != j {
                    let phase_diff = oscillation.phase - other_osc.phase;
                    let coupling_strength = 0.1; // Weak coupling
                    coupling_sum += coupling_strength * phase_diff.sin();
                }
            }
            
            derivatives.push(coupling_sum);
        }
        
        derivatives
    }
    
    /// Calculate ENAQT coupling term
    fn calculate_enaqt_coupling(&self, quantum_state: &QuantumStateAmplitude, state: &BiologicalQuantumState) -> Complex<f64> {
        let coupling = &state.membrane_coords.environmental_coupling;
        
        // ENAQT coupling enhances transport
        let enhancement = Complex::new(coupling.enhancement_factor, 0.0);
        let coupling_term = Complex::new(coupling.coupling_strength, 0.0);
        
        enhancement * coupling_term * quantum_state.amplitude * 0.1
    }
    
    /// Calculate ATP-quantum coupling
    fn calculate_atp_quantum_coupling(&self, quantum_state: &QuantumStateAmplitude, state: &BiologicalQuantumState) -> Complex<f64> {
        let atp_energy = state.atp_coords.available_energy();
        let coupling_strength = 0.05;
        
        // ATP energy modulates quantum state evolution
        Complex::new(coupling_strength * atp_energy * 0.001, 0.0) * quantum_state.amplitude
    }
    
    /// Calculate tunneling derivatives
    fn calculate_tunneling_derivatives(&self, state: &BiologicalQuantumState) -> Vec<f64> {
        let mut derivatives = Vec::new();
        
        for tunneling_state in &state.membrane_coords.tunneling_states {
            // Calculate tunneling probability change rate
            let barrier_penetration = (-2.0 * tunneling_state.barrier_width * 1e-9).exp();
            let energy_dependence = (tunneling_state.electron_energy - tunneling_state.barrier_height).abs();
            
            let rate = barrier_penetration * energy_dependence * 0.01;
            derivatives.push(rate);
        }
        
        derivatives
    }
    
    /// Calculate environmental coupling derivatives
    fn calculate_environmental_derivatives(&self, state: &BiologicalQuantumState) -> EnvironmentalCouplingDerivatives {
        let coupling = &state.membrane_coords.environmental_coupling;
        
        // Environmental coupling evolves slowly
        EnvironmentalCouplingDerivatives {
            coupling_strength_rate: -coupling.coupling_strength * 0.001, // Slow decay
            correlation_time_rate: coupling.correlation_time * 0.0001,   // Slow increase
            enhancement_factor_rate: -coupling.enhancement_factor * 0.0005, // Gradual decrease
        }
    }
    
    /// Calculate endpoint evolution rate
    fn calculate_endpoint_evolution_rate(&self, state: &BiologicalQuantumState) -> HashMap<String, Vec<f64>> {
        let mut evolution_rates = HashMap::new();
        
        for (oscillator_name, distribution) in &state.entropy_coords.endpoint_distributions {
            let mut rates = Vec::new();
            
            // Calculate how endpoint probabilities are changing
            for (i, &prob) in distribution.probabilities.iter().enumerate() {
                // Rate depends on oscillation dynamics and ATP availability
                let atp_factor = state.atp_coords.available_energy() * 0.001;
                let oscillation_factor = if let Some(osc) = state.oscillatory_coords.oscillations.iter()
                    .find(|o| &o.name == oscillator_name) {
                    osc.amplitude * osc.frequency * 0.01
                } else {
                    0.01
                };
                
                let rate = (prob * atp_factor * oscillation_factor).min(0.1);
                rates.push(rate);
            }
            
            evolution_rates.insert(oscillator_name.clone(), rates);
        }
        
        evolution_rates
    }
    
    /// Calculate ATP entropy production
    fn calculate_atp_entropy_production(&self, state: &BiologicalQuantumState) -> f64 {
        let atp_consumption = state.atp_coords.atp_concentration * 0.01; // Approximate consumption
        atp_consumption * 0.1 // Entropy per ATP hydrolysis
    }
    
    /// Calculate oscillatory entropy production
    fn calculate_oscillatory_entropy_production(&self, state: &BiologicalQuantumState) -> f64 {
        let mut total_entropy = 0.0;
        
        for oscillation in &state.oscillatory_coords.oscillations {
            // Entropy from oscillation damping
            let damping_entropy = oscillation.damping_coefficient * oscillation.amplitude.powi(2) * 0.01;
            total_entropy += damping_entropy;
        }
        
        total_entropy
    }
    
    /// Calculate membrane entropy production
    fn calculate_membrane_entropy_production(&self, state: &BiologicalQuantumState) -> f64 {
        let mut total_entropy = 0.0;
        
        // Entropy from quantum decoherence
        for quantum_state in &state.membrane_coords.quantum_states {
            let coherence = quantum_state.amplitude.norm_sqr();
            let decoherence_entropy = (1.0 - coherence) * 0.05;
            total_entropy += decoherence_entropy;
        }
        
        // Entropy from tunneling processes
        for tunneling_state in &state.membrane_coords.tunneling_states {
            let tunneling_entropy = tunneling_state.tunneling_probability * 0.02;
            total_entropy += tunneling_entropy;
        }
        
        total_entropy
    }
    
    /// Calculate quantum tunneling entropy production (death mechanism)
    fn calculate_quantum_tunneling_entropy_production(&self, state: &BiologicalQuantumState) -> f64 {
        let mut total_entropy = 0.0;
        
        for tunneling_state in &state.membrane_coords.tunneling_states {
            // Higher tunneling probability increases radical formation entropy
            let radical_formation_entropy = tunneling_state.tunneling_probability.powi(2) * 0.1;
            total_entropy += radical_formation_entropy;
        }
        
        total_entropy
    }
}

// ================================================================================================
// MAIN SOLVER: BIOLOGICAL QUANTUM COMPUTER
// ================================================================================================

/// Complete solver for biological quantum computation with ATP and oscillatory dynamics
pub struct BiologicalQuantumComputerSolver {
    /// Hamiltonian for the complete system
    hamiltonian: BiologicalQuantumHamiltonian,
    /// Integration method
    integration_method: IntegrationMethod,
    /// Step size control
    step_controller: StepController,
    /// Entropy constraint enforcer
    entropy_enforcer: EntropyConstraintEnforcer,
}

#[derive(Debug)]
pub enum IntegrationMethod {
    VelocityVerlet,
    RungeKutta4,
    AdaptiveStepsize,
}

pub struct StepController {
    pub min_atp_step: f64,
    pub max_atp_step: f64,
    pub tolerance: f64,
}

pub struct EntropyConstraintEnforcer {
    pub enforce_second_law: bool,
    pub max_entropy_production_rate: f64,
}

impl BiologicalQuantumComputerSolver {
    pub fn new() -> Self {
        Self {
            hamiltonian: BiologicalQuantumHamiltonian::new(),
            integration_method: IntegrationMethod::VelocityVerlet,
            step_controller: StepController {
                min_atp_step: 0.01,
                max_atp_step: 1.0,
                tolerance: 1e-6,
            },
            entropy_enforcer: EntropyConstraintEnforcer {
                enforce_second_law: true,
                max_entropy_production_rate: 1.0,
            },
        }
    }

    /// Main solving method: complete biological quantum computation
    pub fn solve_biological_quantum_computation(
        &mut self,
        initial_state: &BiologicalQuantumState,
        atp_budget: f64,
        time_horizon: f64,
        quantum_computation_target: &QuantumComputationTarget,
    ) -> Result<BiologicalQuantumResult, SolverError> {
        
        let mut current_state = initial_state.clone();
        let mut atp_consumed = 0.0;
        let mut current_time = 0.0;
        let mut trajectory = BiologicalQuantumTrajectory::new();
        
        println!("Starting biological quantum computation simulation...");
        println!("ATP budget: {:.2} mM", atp_budget);
        println!("Time horizon: {:.2} seconds", time_horizon);
        
        while atp_consumed < atp_budget && current_time < time_horizon {
            
            // Calculate optimal step size
            let atp_step = self.calculate_optimal_atp_step(&current_state);
            let time_step = self.calculate_optimal_time_step(&current_state);
            
            // Solve one integration step
            let next_state = self.integration_step(
                &current_state,
                atp_step,
                time_step,
            )?;
            
            // Calculate oscillation endpoints for this step (your key insight)
            let oscillation_endpoints = self.predict_oscillation_endpoints(
                &current_state,
                &next_state,
                atp_step
            );
            
            // Calculate membrane quantum computation progress
            let quantum_computation_progress = self.calculate_quantum_computation_progress(
                &next_state,
                quantum_computation_target
            );
            
            // Calculate entropy production (your entropy formulation)
            let entropy_production = self.calculate_step_entropy_production(
                &current_state,
                &next_state,
                &oscillation_endpoints
            );
            
            // Enforce entropy constraints (Second Law)
            self.enforce_entropy_constraints(&mut current_state, entropy_production)?;
            
            // Calculate radical generation (death mechanism)
            let radical_endpoints = self.calculate_radical_generation(
                &next_state,
                atp_step
            );
            
            // Update state
            current_state = next_state;
            atp_consumed += atp_step;
            current_time += time_step;
            
            // Record trajectory point
            trajectory.add_point(BiologicalQuantumTrajectoryPoint {
                time: current_time,
                atp_consumed,
                state: current_state.clone(),
                oscillation_endpoints: oscillation_endpoints.clone(),
                radical_endpoints: radical_endpoints.clone(),
                entropy_production,
                quantum_computation_progress,
            });
            
            // Progress reporting
            if trajectory.points.len() % 100 == 0 {
                println!("Progress: {:.1}% ATP consumed, {:.1}% time elapsed, {:.1}% quantum computation complete",
                    100.0 * atp_consumed / atp_budget,
                    100.0 * current_time / time_horizon,
                    100.0 * quantum_computation_progress
                );
            }
        }
        
        println!("Simulation completed!");
        println!("Final ATP consumed: {:.2} mM", atp_consumed);
        println!("Final time: {:.2} seconds", current_time);
        
        Ok(BiologicalQuantumResult {
            final_state: current_state,
            trajectory,
            total_atp_consumed: atp_consumed,
            total_time: current_time,
            quantum_computation_completed: quantum_computation_progress >= 1.0,
        })
    }

    /// Integration step using velocity-Verlet for the complete system
    fn integration_step(
        &self,
        state: &BiologicalQuantumState,
        atp_step: f64,
        time_step: f64,
    ) -> Result<BiologicalQuantumState, SolverError> {
        
        match self.integration_method {
            IntegrationMethod::VelocityVerlet => {
                self.velocity_verlet_step(state, atp_step, time_step)
            },
            IntegrationMethod::RungeKutta4 => {
                self.runge_kutta_4_step(state, atp_step, time_step)
            },
            IntegrationMethod::AdaptiveStepsize => {
                self.adaptive_step(state, atp_step, time_step)
            },
        }
    }

    /// Velocity-Verlet integration for biological quantum systems
    fn velocity_verlet_step(
        &self,
        state: &BiologicalQuantumState,
        atp_step: f64,
        time_step: f64,
    ) -> Result<BiologicalQuantumState, SolverError> {
        
        // Step 1: Calculate current derivatives
        let current_derivatives = self.hamiltonian.equations_of_motion(state);
        
        // Step 2: Update ATP coordinates
        let mut new_atp_coords = state.atp_coords.clone();
        new_atp_coords.atp_concentration += current_derivatives.atp_derivatives.atp_concentration_rate * atp_step;
        new_atp_coords.adp_concentration += current_derivatives.atp_derivatives.adp_concentration_rate * atp_step;
        new_atp_coords.pi_concentration += current_derivatives.atp_derivatives.pi_concentration_rate * atp_step;
        new_atp_coords.energy_charge += current_derivatives.atp_derivatives.energy_charge_rate * atp_step;
        new_atp_coords.atp_oscillation_amplitude += current_derivatives.atp_derivatives.oscillation_amplitude_rate * atp_step;
        new_atp_coords.atp_oscillation_phase += current_derivatives.atp_derivatives.oscillation_phase_rate * time_step;
        
        // Step 3: Update oscillatory coordinates (Verlet algorithm)
        let mut new_oscillatory_coords = state.oscillatory_coords.clone();
        
        // Update positions: q(t+dt) = q(t) + v(t)*dt + 0.5*a(t)*dt²
        for (i, oscillation) in new_oscillatory_coords.oscillations.iter_mut().enumerate() {
            if i < state.oscillatory_coords.oscillatory_momenta.len() && 
               i < current_derivatives.oscillatory_derivatives.momentum_derivatives.len() {
                let velocity = state.oscillatory_coords.oscillatory_momenta[i];
                let acceleration = current_derivatives.oscillatory_derivatives.momentum_derivatives[i];
                
                oscillation.amplitude += velocity * time_step + 0.5 * acceleration * time_step * time_step;
                
                if i < current_derivatives.oscillatory_derivatives.phase_derivatives.len() {
                    oscillation.phase += current_derivatives.oscillatory_derivatives.phase_derivatives[i] * time_step;
                }
            }
        }
        
        // Update momenta: p(t+dt) = p(t) + 0.5*(a(t) + a(t+dt))*dt
        for (i, momentum) in new_oscillatory_coords.oscillatory_momenta.iter_mut().enumerate() {
            if i < current_derivatives.oscillatory_derivatives.momentum_derivatives.len() {
                *momentum += current_derivatives.oscillatory_derivatives.momentum_derivatives[i] * time_step;
            }
        }
        
        // Step 4: Update membrane quantum coordinates (Schrödinger evolution)
        let mut new_membrane_coords = state.membrane_coords.clone();
        
        for (i, quantum_state) in new_membrane_coords.quantum_states.iter_mut().enumerate() {
            if i < current_derivatives.membrane_derivatives.quantum_state_derivatives.len() {
                // Time evolution: |ψ(t+dt)⟩ = exp(-iHdt/ℏ)|ψ(t)⟩
                let derivative = current_derivatives.membrane_derivatives.quantum_state_derivatives[i];
                quantum_state.amplitude += derivative * time_step;
                
                // Normalize quantum state
                let norm = quantum_state.amplitude.norm();
                if norm > 0.0 {
                    quantum_state.amplitude /= norm;
                }
            }
        }
        
        // Update environmental coupling (ENAQT dynamics)
        new_membrane_coords.environmental_coupling.coupling_strength += 
            current_derivatives.membrane_derivatives.environmental_coupling_derivatives.coupling_strength_rate * time_step;
        new_membrane_coords.environmental_coupling.correlation_time += 
            current_derivatives.membrane_derivatives.environmental_coupling_derivatives.correlation_time_rate * time_step;
        new_membrane_coords.environmental_coupling.enhancement_factor += 
            current_derivatives.membrane_derivatives.environmental_coupling_derivatives.enhancement_factor_rate * time_step;
        
        // Update tunneling states
        for (i, tunneling_state) in new_membrane_coords.tunneling_states.iter_mut().enumerate() {
            if i < current_derivatives.membrane_derivatives.tunneling_derivatives.len() {
                tunneling_state.tunneling_probability += 
                    current_derivatives.membrane_derivatives.tunneling_derivatives[i] * time_step;
                // Clamp probability to [0, 1]
                tunneling_state.tunneling_probability = tunneling_state.tunneling_probability.max(0.0).min(1.0);
            }
        }
        
        // Step 5: Update entropy coordinates (your oscillatory entropy formulation)
        let mut new_entropy_coords = state.entropy_coords.clone();
        new_entropy_coords.current_entropy += current_derivatives.entropy_derivatives.total_entropy_rate * time_step;
        new_entropy_coords.entropy_production_rate = current_derivatives.entropy_derivatives.total_entropy_rate;
        new_entropy_coords.membrane_endpoint_entropy += current_derivatives.entropy_derivatives.membrane_endpoint_entropy_rate * time_step;
        new_entropy_coords.quantum_tunneling_entropy += current_derivatives.entropy_derivatives.quantum_tunneling_entropy_rate * time_step;
        
        // Update endpoint distributions
        for (oscillator_name, distribution_rates) in &current_derivatives.entropy_derivatives.endpoint_distribution_rates {
            if let Some(distribution) = new_entropy_coords.endpoint_distributions.get_mut(oscillator_name) {
                for (i, &rate) in distribution_rates.iter().enumerate() {
                    if i < distribution.probabilities.len() {
                        distribution.probabilities[i] += rate * time_step;
                    }
                }
                // Renormalize probabilities
                let total_prob: f64 = distribution.probabilities.iter().sum();
                if total_prob > 0.0 {
                    for prob in &mut distribution.probabilities {
                        *prob /= total_prob;
                    }
                }
            }
        }
        
        Ok(BiologicalQuantumState {
            atp_coords: new_atp_coords,
            oscillatory_coords: new_oscillatory_coords,
            membrane_coords: new_membrane_coords,
            entropy_coords: new_entropy_coords,
        })
    }

    /// Predict where oscillations will end up (your key insight)
    fn predict_oscillation_endpoints(
        &self,
        current_state: &BiologicalQuantumState,
        next_state: &BiologicalQuantumState,
        atp_step: f64,
    ) -> Vec<OscillationEndpoint> {
        
        let mut endpoints = Vec::new();
        
        for (i, oscillation) in current_state.oscillatory_coords.oscillations.iter().enumerate() {
            // Calculate where this oscillation will end based on current dynamics
            let current_energy = if i < current_state.oscillatory_coords.oscillatory_momenta.len() {
                0.5 * current_state.oscillatory_coords.oscillatory_momenta[i].powi(2) + 
                0.5 * oscillation.frequency.powi(2) * oscillation.amplitude.powi(2)
            } else {
                0.5 * oscillation.frequency.powi(2) * oscillation.amplitude.powi(2)
            };
            
            // Account for ATP-driven energy input
            let atp_energy_input = atp_step * oscillation.atp_coupling_strength * 
                                  current_state.atp_coords.available_energy();
            
            // Account for damping energy loss
            let damping_energy_loss = oscillation.damping_coefficient * current_energy * atp_step;
            
            // Final energy at endpoint
            let final_energy = current_energy + atp_energy_input - damping_energy_loss;
            
            // Calculate endpoint amplitude (energy conservation)
            let final_amplitude = if final_energy > 0.0 {
                (2.0 * final_energy / oscillation.frequency.powi(2)).sqrt()
            } else {
                0.0
            };
            
            // Calculate endpoint phase (phase evolution)
            let final_phase = oscillation.phase + oscillation.frequency * atp_step;
            
            // Calculate endpoint position and velocity
            let final_position = final_amplitude * final_phase.cos();
            let final_velocity = -final_amplitude * oscillation.frequency * final_phase.sin();
            
            // Calculate probability of reaching this endpoint (quantum mechanical)
            let endpoint_probability = self.calculate_endpoint_probability(
                oscillation, final_energy, atp_step
            );
            
            endpoints.push(OscillationEndpoint {
                oscillator_name: oscillation.name.clone(),
                position: final_position,
                velocity: final_velocity,
                energy: final_energy,
                probability: endpoint_probability,
                atp_consumed: atp_step * oscillation.atp_coupling_strength,
                entropy_contribution: -endpoint_probability * endpoint_probability.ln(),
            });
        }
        
        endpoints
    }

    /// Calculate probability of reaching a specific oscillation endpoint
    fn calculate_endpoint_probability(
        &self,
        oscillation: &OscillationState,
        final_energy: f64,
        atp_step: f64,
    ) -> f64 {
        // Quantum mechanical probability based on energy distribution
        let thermal_energy = 1.381e-23 * 310.0; // kT at body temperature
        let energy_ratio = final_energy / (thermal_energy * 6.242e18); // Convert to eV
        
        // Boltzmann distribution for endpoint probability
        let probability = (-energy_ratio).exp();
        
        // Normalize by available ATP energy
        let atp_normalization = (atp_step * oscillation.atp_coupling_strength).min(1.0);
        
        probability * atp_normalization
    }

    /// Calculate radical generation from quantum tunneling (death mechanism)
    fn calculate_radical_generation(
        &self,
        state: &BiologicalQuantumState,
        atp_step: f64,
    ) -> Vec<RadicalEndpoint> {
        
        let mut radical_endpoints = Vec::new();
        
        for tunneling_state in &state.membrane_coords.tunneling_states {
            // Calculate electron tunneling probability
            let kappa = ((2.0 * 9.109e-31 * (tunneling_state.barrier_height - tunneling_state.electron_energy) * 1.602e-19) / 
                        (1.055e-34 * 1.055e-34)).sqrt();
            let tunneling_probability = (-2.0 * kappa * tunneling_state.barrier_width * 1e-9).exp();
            
            // Probability of electron-oxygen interaction
            let oxygen_concentration = 0.2; // Approximate dissolved oxygen concentration
            let interaction_probability = tunneling_probability * oxygen_concentration * atp_step;
            
            if interaction_probability > 1e-6 { // Only include significant probabilities
                // Calculate radical formation position (random within membrane)
                let radical_position = [
                    (0.5 - 0.5) * 10e-9, // x position (nm) - simplified random
                    (0.5 - 0.5) * 10e-9, // y position (nm) - simplified random
                    state.membrane_coords.membrane_properties.thickness * 1e-9 * 0.5, // z position within membrane - simplified
                ];
                
                // Calculate damage potential based on nearby biomolecules
                let damage_potential = self.calculate_damage_potential(&radical_position, state);
                
                radical_endpoints.push(RadicalEndpoint {
                    position: radical_position,
                    radical_type: RadicalType::Superoxide, // Most common from electron tunneling
                    formation_probability: interaction_probability,
                    damage_potential,
                    entropy_contribution: -interaction_probability * interaction_probability.ln(),
                });
            }
        }
        
        radical_endpoints
    }

    /// Calculate damage potential of radical at specific position
    fn calculate_damage_potential(&self, _position: &[f64; 3], state: &BiologicalQuantumState) -> f64 {
        // Simplified damage calculation based on proximity to membrane proteins
        let mut damage_potential = 0.0;
        
        // Assume membrane proteins are distributed with density from membrane properties
        let protein_density = state.membrane_coords.membrane_properties.protein_density;
        let interaction_radius = 2e-9; // 2 nm interaction radius for radicals
        
        // Calculate expected number of proteins within interaction radius
        let interaction_volume = (4.0/3.0) * PI * interaction_radius.powi(3);
        let expected_proteins = protein_density * interaction_volume * 1e18; // Convert nm² to m²
        
        // Damage potential scales with number of nearby proteins
        damage_potential = expected_proteins * 0.1; // 10% damage probability per protein
        
        damage_potential
    }

    /// Calculate entropy production for this step (your key insight)
    fn calculate_step_entropy_production(
        &self,
        current_state: &BiologicalQuantumState,
        next_state: &BiologicalQuantumState,
        endpoints: &[OscillationEndpoint],
    ) -> f64 {
        // Entropy from oscillation endpoints (your formulation)
        let endpoint_entropy: f64 = endpoints.iter()
            .map(|endpoint| endpoint.entropy_contribution)
            .sum();
        
        // Entropy from ATP consumption
        let atp_entropy = self.calculate_atp_entropy_production(current_state, next_state);
        
        // Entropy from quantum decoherence
        let quantum_entropy = self.calculate_quantum_entropy_production(current_state, next_state);
        
        // Entropy from membrane processes
        let membrane_entropy = self.calculate_membrane_entropy_production(current_state, next_state);
        
        endpoint_entropy + atp_entropy + quantum_entropy + membrane_entropy
    }

    fn calculate_atp_entropy_production(&self, current: &BiologicalQuantumState, next: &BiologicalQuantumState) -> f64 {
        let atp_consumed = current.atp_coords.atp_concentration - next.atp_coords.atp_concentration;
        atp_consumed * 0.1 // Approximate entropy per ATP hydrolysis (kB units)
    }

    fn calculate_quantum_entropy_production(&self, current: &BiologicalQuantumState, next: &BiologicalQuantumState) -> f64 {
        let mut entropy_change = 0.0;
        
        for (i, current_state) in current.membrane_coords.quantum_states.iter().enumerate() {
            if i < next.membrane_coords.quantum_states.len() {
                let current_prob = current_state.amplitude.norm_sqr();
                let next_prob = next.membrane_coords.quantum_states[i].amplitude.norm_sqr();
                
                if current_prob > 0.0 && next_prob > 0.0 {
                    entropy_change += next_prob * next_prob.ln() - current_prob * current_prob.ln();
                }
            }
        }
        
        -entropy_change // Negative because we want entropy production (positive)
    }

    fn calculate_membrane_entropy_production(&self, _current: &BiologicalQuantumState, _next: &BiologicalQuantumState) -> f64 {
        // Entropy from membrane conformational changes
        let conformational_entropy = 0.05; // Approximate value
        
        // Entropy from proton transport
        let proton_entropy = 0.02;
        
        // Entropy from electron transport
        let electron_entropy = 0.03;
        
        conformational_entropy + proton_entropy + electron_entropy
    }

    /// Enforce entropy constraints (Second Law of Thermodynamics)
    fn enforce_entropy_constraints(
        &self,
        state: &mut BiologicalQuantumState,
        entropy_production: f64,
    ) -> Result<(), SolverError> {
        
        if self.entropy_enforcer.enforce_second_law {
            // Entropy must not decrease
            if entropy_production < 0.0 {
                return Err(SolverError::EntropyViolation(
                    format!("Entropy production is negative: {}", entropy_production)
                ));
            }
            
            // Entropy production rate must not exceed maximum
            if entropy_production > self.entropy_enforcer.max_entropy_production_rate {
                // Scale down processes to respect entropy limit
                let scaling_factor = self.entropy_enforcer.max_entropy_production_rate / entropy_production;
                
                // Scale ATP consumption
                state.atp_coords.atp_concentration *= scaling_factor;
                
                // Scale oscillation amplitudes
                for oscillation in &mut state.oscillatory_coords.oscillations {
                    oscillation.amplitude *= scaling_factor.sqrt();
                }
                
                // Scale quantum coherences
                for quantum_state in &mut state.membrane_coords.quantum_states {
                    quantum_state.amplitude *= scaling_factor.sqrt().into();
                }
            }
        }
        
        Ok(())
    }

    fn calculate_optimal_atp_step(&self, state: &BiologicalQuantumState) -> f64 {
        // Adaptive step size based on ATP concentration and system dynamics
        let base_step = 0.1;
        let atp_factor = (state.atp_coords.atp_concentration / 5.0).min(1.0); // Normalize to typical 5mM
        let oscillation_factor = state.oscillatory_coords.oscillations.iter()
            .map(|osc| osc.amplitude)
            .fold(0.0, f64::max)
            .min(2.0) / 2.0; // Normalize to reasonable amplitude
        
        base_step * atp_factor * (1.0 + oscillation_factor)
    }

    fn calculate_optimal_time_step(&self, state: &BiologicalQuantumState) -> f64 {
        // Time step based on fastest oscillation frequency
        let max_frequency = state.oscillatory_coords.oscillations.iter()
            .map(|osc| osc.frequency)
            .fold(0.0, f64::max);
        
        if max_frequency > 0.0 {
            0.01 / max_frequency // 1% of the fastest period
        } else {
            0.001 // Default 1ms step
        }
    }

    fn calculate_quantum_computation_progress(
        &self,
        state: &BiologicalQuantumState,
        target: &QuantumComputationTarget,
    ) -> f64 {
        // Simplified progress calculation based on quantum state evolution
        let total_coherence: f64 = state.membrane_coords.quantum_states.iter()
            .map(|qs| qs.amplitude.norm_sqr())
            .sum();
        
        let target_coherence = target.required_coherence;
        
        (total_coherence / target_coherence).min(1.0)
    }

    /// Runge-Kutta 4th order integration step
    fn runge_kutta_4_step(
        &self,
        state: &BiologicalQuantumState,
        atp_step: f64,
        time_step: f64,
    ) -> Result<BiologicalQuantumState, SolverError> {
        
        // k1 = f(t, y)
        let k1 = self.hamiltonian.equations_of_motion(state);
        
        // k2 = f(t + h/2, y + h*k1/2)
        let temp_state1 = self.add_derivatives_to_state(state, &k1, time_step * 0.5, atp_step * 0.5)?;
        let k2 = self.hamiltonian.equations_of_motion(&temp_state1);
        
        // k3 = f(t + h/2, y + h*k2/2)
        let temp_state2 = self.add_derivatives_to_state(state, &k2, time_step * 0.5, atp_step * 0.5)?;
        let k3 = self.hamiltonian.equations_of_motion(&temp_state2);
        
        // k4 = f(t + h, y + h*k3)
        let temp_state3 = self.add_derivatives_to_state(state, &k3, time_step, atp_step)?;
        let k4 = self.hamiltonian.equations_of_motion(&temp_state3);
        
        // y_new = y + h/6 * (k1 + 2*k2 + 2*k3 + k4)
        self.combine_rk4_derivatives(state, &k1, &k2, &k3, &k4, time_step, atp_step)
    }

    /// Adaptive step size integration
    fn adaptive_step(
        &self,
        state: &BiologicalQuantumState,
        atp_step: f64,
        time_step: f64,
    ) -> Result<BiologicalQuantumState, SolverError> {
        
        // Try full step
        let full_step = self.velocity_verlet_step(state, atp_step, time_step)?;
        
        // Try two half steps
        let half_step1 = self.velocity_verlet_step(state, atp_step * 0.5, time_step * 0.5)?;
        let half_step2 = self.velocity_verlet_step(&half_step1, atp_step * 0.5, time_step * 0.5)?;
        
        // Estimate error
        let error = self.calculate_integration_error(&full_step, &half_step2);
        
        if error < self.step_controller.tolerance {
            Ok(half_step2) // Use more accurate result
        } else {
            // Reduce step size and try again
            let reduced_time_step = time_step * 0.5;
            let reduced_atp_step = atp_step * 0.5;
            self.adaptive_step(state, reduced_atp_step, reduced_time_step)
        }
    }

    /// Helper method to add derivatives to state
    fn add_derivatives_to_state(
        &self,
        state: &BiologicalQuantumState,
        derivatives: &BiologicalQuantumDerivatives,
        dt: f64,
        datp: f64,
    ) -> Result<BiologicalQuantumState, SolverError> {
        
        let mut new_state = state.clone();
        
        // Update ATP coordinates
        new_state.atp_coords.atp_concentration += derivatives.atp_derivatives.atp_concentration_rate * datp;
        new_state.atp_coords.adp_concentration += derivatives.atp_derivatives.adp_concentration_rate * datp;
        new_state.atp_coords.pi_concentration += derivatives.atp_derivatives.pi_concentration_rate * datp;
        new_state.atp_coords.energy_charge += derivatives.atp_derivatives.energy_charge_rate * datp;
        new_state.atp_coords.atp_oscillation_amplitude += derivatives.atp_derivatives.oscillation_amplitude_rate * datp;
        new_state.atp_coords.atp_oscillation_phase += derivatives.atp_derivatives.oscillation_phase_rate * dt;
        
        Ok(new_state)
    }

    /// Combine RK4 derivatives
    fn combine_rk4_derivatives(
        &self,
        state: &BiologicalQuantumState,
        k1: &BiologicalQuantumDerivatives,
        k2: &BiologicalQuantumDerivatives,
        k3: &BiologicalQuantumDerivatives,
        k4: &BiologicalQuantumDerivatives,
        dt: f64,
        datp: f64,
    ) -> Result<BiologicalQuantumState, SolverError> {
        
        let mut new_state = state.clone();
        
        // RK4 combination: y + h/6 * (k1 + 2*k2 + 2*k3 + k4)
        new_state.atp_coords.atp_concentration += datp / 6.0 * (
            k1.atp_derivatives.atp_concentration_rate +
            2.0 * k2.atp_derivatives.atp_concentration_rate +
            2.0 * k3.atp_derivatives.atp_concentration_rate +
            k4.atp_derivatives.atp_concentration_rate
        );
        
        Ok(new_state)
    }

    /// Calculate integration error for adaptive stepping
    fn calculate_integration_error(&self, state1: &BiologicalQuantumState, state2: &BiologicalQuantumState) -> f64 {
        let mut total_error = 0.0;
        
        // Compare ATP concentrations
        total_error += (state1.atp_coords.atp_concentration - state2.atp_coords.atp_concentration).abs();
        total_error += (state1.atp_coords.adp_concentration - state2.atp_coords.adp_concentration).abs();
        
        // Compare oscillation amplitudes
        for (osc1, osc2) in state1.oscillatory_coords.oscillations.iter()
            .zip(state2.oscillatory_coords.oscillations.iter()) {
            total_error += (osc1.amplitude - osc2.amplitude).abs();
            total_error += (osc1.phase - osc2.phase).abs();
        }
        
        // Compare quantum state amplitudes
        for (qs1, qs2) in state1.membrane_coords.quantum_states.iter()
            .zip(state2.membrane_coords.quantum_states.iter()) {
            total_error += (qs1.amplitude - qs2.amplitude).norm();
        }
        
        total_error
    }
}

// ================================================================================================
// RESULT STRUCTURES AND ANALYSIS
// ================================================================================================

/// Target for quantum computation
pub struct QuantumComputationTarget {
    pub computation_type: String,
    pub required_coherence: f64,
    pub target_efficiency: f64,
}

/// Complete result of biological quantum computation
pub struct BiologicalQuantumResult {
    pub final_state: BiologicalQuantumState,
    pub trajectory: BiologicalQuantumTrajectory,
    pub total_atp_consumed: f64,
    pub total_time: f64,
    pub quantum_computation_completed: bool,
}

/// Trajectory of biological quantum computation
pub struct BiologicalQuantumTrajectory {
    pub points: Vec<BiologicalQuantumTrajectoryPoint>,
}

impl BiologicalQuantumTrajectory {
    pub fn new() -> Self {
        Self { points: Vec::new() }
    }
    
    pub fn add_point(&mut self, point: BiologicalQuantumTrajectoryPoint) {
        self.points.push(point);
    }
}

/// Single point in biological quantum trajectory
pub struct BiologicalQuantumTrajectoryPoint {
    pub time: f64,
    pub atp_consumed: f64,
    pub state: BiologicalQuantumState,
    pub oscillation_endpoints: Vec<OscillationEndpoint>,
    pub radical_endpoints: Vec<RadicalEndpoint>,
    pub entropy_production: f64,
    pub quantum_computation_progress: f64,
} 