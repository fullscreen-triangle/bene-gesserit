use crate::types::*;
use crate::error::BeneGesseritError;
use crate::biological_quantum_computer::BiologicalQuantumComputerSolver;
use crate::biological_maxwell_demons::*;
use ndarray::{Array1, Array2};
use num_complex::Complex;

/// Enhanced derivatives that include information flow from BMDs
#[derive(Debug, Clone)]
pub struct EnhancedDerivatives {
    pub atp_derivatives: AtpDerivatives,
    pub oscillatory_derivatives: OscillatoryDerivatives,
    pub membrane_quantum_derivatives: MembraneQuantumDerivatives,
    pub entropy_derivatives: OscillatoryEntropyDerivatives,
    pub information_flow: InformationFlow,
}

/// Information flow metrics between BMDs
#[derive(Debug, Clone)]
pub struct InformationFlow {
    pub atp_to_oscillatory: f64,
    pub oscillatory_to_quantum: f64,
    pub quantum_to_atp: f64,
    pub total_information_processing: f64,
    pub catalytic_efficiency: f64,
}

/// Enhanced trajectory point with BMD information
#[derive(Debug, Clone)]
pub struct EnhancedTrajectoryPoint {
    pub time: f64,
    pub state: BiologicalQuantumState,
    pub atp_remaining: f64,
    pub bmd_metrics: CombinedBmdMetrics,
    pub information_flow: InformationFlow,
    pub energy_allocation: EnergyAllocation,
    pub oscillation_control: OscillationControl,
    pub quantum_operation: QuantumOperation,
}

/// Combined metrics from all BMDs
#[derive(Debug, Clone)]
pub struct CombinedBmdMetrics {
    pub atp_metrics: InformationCatalysisMetrics,
    pub oscillatory_metrics: InformationCatalysisMetrics,
    pub quantum_metrics: InformationCatalysisMetrics,
    pub overall_efficiency: f64,
    pub system_coherence: f64,
}

impl CombinedBmdMetrics {
    pub fn new(
        atp_metrics: InformationCatalysisMetrics,
        oscillatory_metrics: InformationCatalysisMetrics,
        quantum_metrics: InformationCatalysisMetrics,
    ) -> Self {
        let overall_efficiency = (
            atp_metrics.calculate_overall_efficiency() *
            oscillatory_metrics.calculate_overall_efficiency() *
            quantum_metrics.calculate_overall_efficiency()
        ).powf(1.0/3.0); // Geometric mean
        
        let system_coherence = (
            atp_metrics.targeting_accuracy +
            oscillatory_metrics.targeting_accuracy +
            quantum_metrics.targeting_accuracy
        ) / 3.0;
        
        Self {
            atp_metrics,
            oscillatory_metrics,
            quantum_metrics,
            overall_efficiency,
            system_coherence,
        }
    }
}

/// Enhanced trajectory with BMD information
#[derive(Debug, Clone)]
pub struct EnhancedQuantumTrajectory {
    pub points: Vec<EnhancedTrajectoryPoint>,
    pub total_computation_time: f64,
    pub total_atp_consumed: f64,
    pub quantum_targets_achieved: Vec<bool>,
    pub information_processing_efficiency: f64,
}

impl EnhancedQuantumTrajectory {
    pub fn new() -> Self {
        Self {
            points: Vec::new(),
            total_computation_time: 0.0,
            total_atp_consumed: 0.0,
            quantum_targets_achieved: Vec::new(),
            information_processing_efficiency: 0.0,
        }
    }
    
    pub fn add_point(&mut self, point: EnhancedTrajectoryPoint) {
        self.total_computation_time = point.time;
        self.total_atp_consumed += point.energy_allocation.energy_cost;
        self.points.push(point);
    }
    
    pub fn finalize(&mut self, targets_achieved: Vec<bool>) {
        self.quantum_targets_achieved = targets_achieved;
        
        // Calculate overall information processing efficiency
        if !self.points.is_empty() {
            let total_efficiency: f64 = self.points.iter()
                .map(|p| p.bmd_metrics.overall_efficiency)
                .sum();
            self.information_processing_efficiency = total_efficiency / self.points.len() as f64;
        }
    }
}

/// BMD-Enhanced Solver that integrates all Maxwell's Demons
pub struct BmdEnhancedSolver {
    /// Core biological quantum computer
    pub core_solver: BiologicalQuantumComputerSolver,
    /// ATP Maxwell's demon
    pub atp_demon: AtpMaxwellsDemon,
    /// Oscillatory Maxwell's demon
    pub oscillatory_demon: OscillatoryMaxwellsDemon,
    /// Membrane quantum Maxwell's demon
    pub quantum_demon: MembraneQuantumMaxwellsDemon,
    /// Combined catalysis metrics
    pub combined_metrics: CombinedBmdMetrics,
    /// BMD renewal thresholds
    pub renewal_threshold: f64,
}

impl BmdEnhancedSolver {
    pub fn new(num_oscillators: usize, num_qubits: usize) -> Self {
        let core_solver = BiologicalQuantumComputerSolver::new();
        let atp_demon = AtpMaxwellsDemon::new();
        let oscillatory_demon = OscillatoryMaxwellsDemon::new(num_oscillators);
        let quantum_demon = MembraneQuantumMaxwellsDemon::new(num_qubits);
        
        let combined_metrics = CombinedBmdMetrics::new(
            atp_demon.metrics.clone(),
            oscillatory_demon.metrics.clone(),
            quantum_demon.metrics.clone(),
        );
        
        Self {
            core_solver,
            atp_demon,
            oscillatory_demon,
            quantum_demon,
            combined_metrics,
            renewal_threshold: 0.8,
        }
    }
    
    /// Main solving method with information catalysis
    pub fn solve_with_information_catalysis(
        &mut self,
        initial_state: BiologicalQuantumState,
        atp_budget: f64,
        time_horizon: f64,
        quantum_targets: &[Complex<f64>],
    ) -> Result<EnhancedQuantumTrajectory, BeneGesseritError> {
        
        let mut trajectory = EnhancedQuantumTrajectory::new();
        let mut current_state = initial_state;
        let mut remaining_atp = atp_budget;
        let dt = time_horizon / 1000.0; // Adaptive step size
        
        for step in 0..1000 {
            // Check if any BMD needs renewal
            self.check_and_renew_bmds()?;
            
            // 1. ATP Maxwell's demon processing
            let atp_allocation = self.atp_demon.catalytic_cycle(
                current_state.atp_coordinates.clone()
            )?;
            
            // 2. Oscillatory Maxwell's demon processing
            let oscillation_control = self.oscillatory_demon.catalytic_cycle(
                current_state.oscillatory_coordinates.clone()
            )?;
            
            // 3. Quantum Maxwell's demon processing
            let quantum_operation = self.quantum_demon.catalytic_cycle(
                current_state.membrane_quantum_coordinates.clone()
            )?;
            
            // 4. Apply information-guided evolution
            let enhanced_derivatives = self.calculate_bmd_enhanced_derivatives(
                &current_state,
                &atp_allocation,
                &oscillation_control,
                &quantum_operation
            )?;
            
            // 5. Evolve state using enhanced derivatives
            current_state = self.evolve_state_with_bmd(
                current_state,
                enhanced_derivatives.clone(),
                dt
            )?;
            
            // 6. Update ATP budget based on BMD decisions
            remaining_atp -= atp_allocation.energy_cost;
            
            // 7. Update combined metrics
            self.update_combined_metrics();
            
            // 8. Record trajectory point
            let trajectory_point = EnhancedTrajectoryPoint {
                time: step as f64 * dt,
                state: current_state.clone(),
                atp_remaining: remaining_atp,
                bmd_metrics: self.combined_metrics.clone(),
                information_flow: enhanced_derivatives.information_flow,
                energy_allocation: atp_allocation,
                oscillation_control,
                quantum_operation,
            };
            trajectory.add_point(trajectory_point);
            
            // 9. Check termination conditions
            if remaining_atp <= 0.0 {
                break;
            }
            
            if self.quantum_targets_achieved(&current_state, quantum_targets) {
                break;
            }
        }
        
        // Finalize trajectory
        let targets_achieved = self.evaluate_quantum_targets(&current_state, quantum_targets);
        trajectory.finalize(targets_achieved);
        
        Ok(trajectory)
    }
    
    /// Calculate BMD-enhanced derivatives
    fn calculate_bmd_enhanced_derivatives(
        &self,
        state: &BiologicalQuantumState,
        atp_allocation: &EnergyAllocation,
        oscillation_control: &OscillationControl,
        quantum_operation: &QuantumOperation,
    ) -> Result<EnhancedDerivatives, BeneGesseritError> {
        
        // Base derivatives from core solver
        let base_derivatives = self.core_solver.calculate_derivatives(state)?;
        
        // ATP enhancement based on BMD allocation decisions
        let enhanced_atp_derivatives = self.enhance_atp_derivatives(
            &base_derivatives.atp_derivatives,
            atp_allocation
        );
        
        // Oscillatory enhancement based on pattern recognition
        let enhanced_oscillatory_derivatives = self.enhance_oscillatory_derivatives(
            &base_derivatives.oscillatory_derivatives,
            oscillation_control
        );
        
        // Quantum enhancement based on information processing
        let enhanced_quantum_derivatives = self.enhance_quantum_derivatives(
            &base_derivatives.membrane_quantum_derivatives,
            quantum_operation
        );
        
        // Calculate information flow
        let information_flow = self.calculate_information_flow(state);
        
        Ok(EnhancedDerivatives {
            atp_derivatives: enhanced_atp_derivatives,
            oscillatory_derivatives: enhanced_oscillatory_derivatives,
            membrane_quantum_derivatives: enhanced_quantum_derivatives,
            entropy_derivatives: base_derivatives.entropy_derivatives,
            information_flow,
        })
    }
    
    /// Enhance ATP derivatives based on BMD allocation decisions
    fn enhance_atp_derivatives(
        &self,
        base_derivatives: &AtpDerivatives,
        allocation: &EnergyAllocation,
    ) -> AtpDerivatives {
        let enhancement_factor = allocation.allocation_confidence;
        
        AtpDerivatives {
            concentration_rate: base_derivatives.concentration_rate * enhancement_factor,
            energy_charge_rate: base_derivatives.energy_charge_rate * enhancement_factor,
            oscillation_coupling_rate: base_derivatives.oscillation_coupling_rate * enhancement_factor,
        }
    }
    
    /// Enhance oscillatory derivatives based on BMD control decisions
    fn enhance_oscillatory_derivatives(
        &self,
        base_derivatives: &OscillatoryDerivatives,
        control: &OscillationControl,
    ) -> OscillatoryDerivatives {
        let mut enhanced = base_derivatives.clone();
        
        // Apply frequency adjustments
        for (i, freq_rate) in enhanced.frequency_rates.iter_mut().enumerate() {
            *freq_rate += control.frequency_adjustment * control.control_confidence;
        }
        
        // Apply amplitude scaling
        for (i, amp_rate) in enhanced.amplitude_rates.iter_mut().enumerate() {
            *amp_rate *= control.amplitude_scaling;
        }
        
        // Apply phase shifts
        for (i, phase_rate) in enhanced.phase_rates.iter_mut().enumerate() {
            *phase_rate += control.phase_shift * control.control_confidence;
        }
        
        enhanced
    }
    
    /// Enhance quantum derivatives based on BMD operation decisions
    fn enhance_quantum_derivatives(
        &self,
        base_derivatives: &MembraneQuantumDerivatives,
        operation: &QuantumOperation,
    ) -> MembraneQuantumDerivatives {
        let mut enhanced = base_derivatives.clone();
        
        // Apply quantum operation enhancements
        let enhancement = operation.success_probability;
        
        for (i, amp_rate) in enhanced.amplitude_rates.iter_mut().enumerate() {
            if operation.target_qubits.contains(&i) {
                *amp_rate *= Complex::new(enhancement, 0.0);
            }
        }
        
        // Enhance tunneling rates
        for (i, tunnel_rate) in enhanced.tunneling_rates.iter_mut().enumerate() {
            if operation.target_qubits.contains(&i) {
                *tunnel_rate *= enhancement;
            }
        }
        
        enhanced
    }
    
    /// Calculate information flow between BMDs
    fn calculate_information_flow(&self, state: &BiologicalQuantumState) -> InformationFlow {
        let atp_efficiency = self.atp_demon.information_efficiency();
        let osc_efficiency = self.oscillatory_demon.information_efficiency();
        let quantum_efficiency = self.quantum_demon.information_efficiency();
        
        let atp_to_oscillatory = atp_efficiency * state.atp_coordinates.energy_charge;
        let oscillatory_to_quantum = osc_efficiency * state.oscillatory_coordinates.frequencies.iter().sum::<f64>();
        let quantum_to_atp = quantum_efficiency * state.membrane_quantum_coordinates.quantum_amplitudes.iter()
            .map(|amp| amp.norm_sqr()).sum::<f64>();
        
        let total_information_processing = atp_to_oscillatory + oscillatory_to_quantum + quantum_to_atp;
        let catalytic_efficiency = (atp_efficiency * osc_efficiency * quantum_efficiency).powf(1.0/3.0);
        
        InformationFlow {
            atp_to_oscillatory,
            oscillatory_to_quantum,
            quantum_to_atp,
            total_information_processing,
            catalytic_efficiency,
        }
    }
    
    /// Evolve state using BMD-enhanced derivatives
    fn evolve_state_with_bmd(
        &self,
        mut state: BiologicalQuantumState,
        derivatives: EnhancedDerivatives,
        dt: f64,
    ) -> Result<BiologicalQuantumState, BeneGesseritError> {
        
        // Evolve ATP coordinates
        state.atp_coordinates.concentration += derivatives.atp_derivatives.concentration_rate * dt;
        state.atp_coordinates.energy_charge += derivatives.atp_derivatives.energy_charge_rate * dt;
        state.atp_coordinates.oscillation_coupling += derivatives.atp_derivatives.oscillation_coupling_rate * dt;
        
        // Evolve oscillatory coordinates
        for (i, freq) in state.oscillatory_coordinates.frequencies.iter_mut().enumerate() {
            if i < derivatives.oscillatory_derivatives.frequency_rates.len() {
                *freq += derivatives.oscillatory_derivatives.frequency_rates[i] * dt;
            }
        }
        
        for (i, amp) in state.oscillatory_coordinates.amplitudes.iter_mut().enumerate() {
            if i < derivatives.oscillatory_derivatives.amplitude_rates.len() {
                *amp += derivatives.oscillatory_derivatives.amplitude_rates[i] * dt;
            }
        }
        
        for (i, phase) in state.oscillatory_coordinates.phases.iter_mut().enumerate() {
            if i < derivatives.oscillatory_derivatives.phase_rates.len() {
                *phase += derivatives.oscillatory_derivatives.phase_rates[i] * dt;
            }
        }
        
        // Evolve quantum coordinates
        for (i, amp) in state.membrane_quantum_coordinates.quantum_amplitudes.iter_mut().enumerate() {
            if i < derivatives.membrane_quantum_derivatives.amplitude_rates.len() {
                *amp += derivatives.membrane_quantum_derivatives.amplitude_rates[i] * dt;
            }
        }
        
        // Ensure physical constraints
        state.atp_coordinates.concentration = state.atp_coordinates.concentration.max(0.0);
        state.atp_coordinates.energy_charge = state.atp_coordinates.energy_charge.clamp(0.0, 1.0);
        
        Ok(state)
    }
    
    /// Check if BMDs need renewal and renew them
    fn check_and_renew_bmds(&mut self) -> Result<(), BeneGesseritError> {
        if self.atp_demon.needs_renewal() {
            self.atp_demon = AtpMaxwellsDemon::new();
        }
        
        if self.oscillatory_demon.needs_renewal() {
            let num_oscillators = self.oscillatory_demon.frequency_filters.len();
            self.oscillatory_demon = OscillatoryMaxwellsDemon::new(num_oscillators);
        }
        
        if self.quantum_demon.needs_renewal() {
            let num_qubits = self.quantum_demon.quantum_recognition_operators[0].dim().0;
            self.quantum_demon = MembraneQuantumMaxwellsDemon::new(num_qubits);
        }
        
        Ok(())
    }
    
    /// Update combined metrics from all BMDs
    fn update_combined_metrics(&mut self) {
        self.combined_metrics = CombinedBmdMetrics::new(
            self.atp_demon.metrics.clone(),
            self.oscillatory_demon.metrics.clone(),
            self.quantum_demon.metrics.clone(),
        );
    }
    
    /// Check if quantum targets are achieved
    fn quantum_targets_achieved(&self, state: &BiologicalQuantumState, targets: &[Complex<f64>]) -> bool {
        if targets.is_empty() {
            return false;
        }
        
        for (i, target) in targets.iter().enumerate() {
            if i < state.membrane_quantum_coordinates.quantum_amplitudes.len() {
                let current_amp = state.membrane_quantum_coordinates.quantum_amplitudes[i];
                let distance = (current_amp - target).norm();
                if distance < 0.1 { // Threshold for achievement
                    return true;
                }
            }
        }
        
        false
    }
    
    /// Evaluate how well quantum targets were achieved
    fn evaluate_quantum_targets(&self, state: &BiologicalQuantumState, targets: &[Complex<f64>]) -> Vec<bool> {
        let mut achieved = Vec::new();
        
        for (i, target) in targets.iter().enumerate() {
            if i < state.membrane_quantum_coordinates.quantum_amplitudes.len() {
                let current_amp = state.membrane_quantum_coordinates.quantum_amplitudes[i];
                let distance = (current_amp - target).norm();
                achieved.push(distance < 0.1);
            } else {
                achieved.push(false);
            }
        }
        
        achieved
    }
    
    /// Get current system efficiency
    pub fn system_efficiency(&self) -> f64 {
        self.combined_metrics.overall_efficiency
    }
    
    /// Get current system coherence
    pub fn system_coherence(&self) -> f64 {
        self.combined_metrics.system_coherence
    }
    
    /// Get total catalytic cycles performed
    pub fn total_catalytic_cycles(&self) -> u64 {
        self.atp_demon.metrics.cycle_count +
        self.oscillatory_demon.metrics.cycle_count +
        self.quantum_demon.metrics.cycle_count
    }
}

/// Factory function for creating BMD-enhanced solver with default parameters
pub fn create_bmd_enhanced_solver() -> BmdEnhancedSolver {
    BmdEnhancedSolver::new(10, 5) // Default: 10 oscillators, 5 qubits
}

/// Factory function for creating BMD-enhanced solver with custom parameters
pub fn create_custom_bmd_enhanced_solver(num_oscillators: usize, num_qubits: usize) -> BmdEnhancedSolver {
    BmdEnhancedSolver::new(num_oscillators, num_qubits)
} 