//! Feedback loops between membrane dynamics and circuit parameters
//!
//! This module implements the dynamic feedback mechanisms that couple
//! membrane state to circuit behavior, including:
//! - Voltage-dependent membrane properties
//! - Activity-dependent ATP consumption
//! - Adaptive circuit parameter tuning
//! - Homeostatic regulation mechanisms

use crate::{
    constants::*,
    error::{MembraneError, Result},
    types::*,
};
use std::collections::HashMap;

/// Types of feedback loops in the system
#[derive(Debug, Clone, PartialEq)]
pub enum FeedbackType {
    /// Negative feedback (stabilizing)
    Negative,
    /// Positive feedback (amplifying)
    Positive,
    /// Adaptive feedback (learning)
    Adaptive,
    /// Homeostatic feedback (maintaining set point)
    Homeostatic,
}

/// Feedback loop structure
#[derive(Debug, Clone)]
pub struct FeedbackLoop {
    /// Unique identifier for the feedback loop
    pub id: String,
    /// Type of feedback
    pub feedback_type: FeedbackType,
    /// Input parameter (what triggers the feedback)
    pub input_parameter: String,
    /// Output parameter (what gets modified)
    pub output_parameter: String,
    /// Feedback strength (gain)
    pub gain: f64,
    /// Time constant for feedback response
    pub time_constant: f64,
    /// Current feedback value
    pub current_value: f64,
    /// Setpoint for homeostatic feedback
    pub setpoint: Option<f64>,
    /// Integration history for adaptive feedback
    pub integration_history: Vec<f64>,
}

/// Feedback system manager
#[derive(Debug, Clone)]
pub struct FeedbackSystem {
    /// Active feedback loops
    pub feedback_loops: HashMap<String, FeedbackLoop>,
    /// System state tracking
    pub state_tracker: StateTracker,
    /// Feedback statistics
    pub statistics: FeedbackStatistics,
    /// Control parameters
    pub control_params: ControlParameters,
}

/// State tracker for monitoring system variables
#[derive(Debug, Clone)]
pub struct StateTracker {
    /// Current membrane voltage
    pub membrane_voltage: f64,
    /// Current ATP concentration
    pub atp_concentration: f64,
    /// Ion channel activities
    pub channel_activities: HashMap<String, f64>,
    /// Protein conformational states
    pub protein_states: HashMap<String, f64>,
    /// Circuit parameters
    pub circuit_params: HashMap<String, f64>,
    /// Historical values for trend analysis
    pub history: HashMap<String, Vec<f64>>,
}

/// Control parameters for feedback system
#[derive(Debug, Clone)]
pub struct ControlParameters {
    /// Maximum allowed feedback gain
    pub max_gain: f64,
    /// Minimum time constant
    pub min_time_constant: f64,
    /// Maximum time constant
    pub max_time_constant: f64,
    /// History length for integration
    pub history_length: usize,
    /// Feedback update rate
    pub update_rate: f64,
}

/// Feedback system statistics
#[derive(Debug, Clone)]
pub struct FeedbackStatistics {
    /// Total feedback loops active
    pub active_loops: usize,
    /// Average feedback strength
    pub avg_feedback_strength: f64,
    /// System stability measure
    pub stability_index: f64,
    /// Adaptation rate
    pub adaptation_rate: f64,
    /// Total energy consumed by feedback
    pub total_energy_consumed: f64,
}

impl FeedbackSystem {
    /// Create a new feedback system
    pub fn new() -> Self {
        Self {
            feedback_loops: HashMap::new(),
            state_tracker: StateTracker::new(),
            statistics: FeedbackStatistics::new(),
            control_params: ControlParameters::default(),
        }
    }
    
    /// Add a feedback loop to the system
    pub fn add_feedback_loop(&mut self, loop_def: FeedbackLoop) -> Result<()> {
        // Validate feedback loop parameters
        self.validate_feedback_loop(&loop_def)?;
        
        // Add to system
        self.feedback_loops.insert(loop_def.id.clone(), loop_def);
        
        Ok(())
    }
    
    /// Validate feedback loop parameters
    fn validate_feedback_loop(&self, feedback_loop: &FeedbackLoop) -> Result<()> {
        // Check gain limits
        if feedback_loop.gain.abs() > self.control_params.max_gain {
            return Err(MembraneError::InvalidParameter(
                format!("Feedback gain {} exceeds maximum {}", 
                       feedback_loop.gain, self.control_params.max_gain)
            ));
        }
        
        // Check time constant limits
        if feedback_loop.time_constant < self.control_params.min_time_constant ||
           feedback_loop.time_constant > self.control_params.max_time_constant {
            return Err(MembraneError::InvalidParameter(
                format!("Time constant {} outside allowed range [{}, {}]",
                       feedback_loop.time_constant,
                       self.control_params.min_time_constant,
                       self.control_params.max_time_constant)
            ));
        }
        
        Ok(())
    }
    
    /// Update all feedback loops
    pub fn update_feedback(&mut self, membrane_state: &MembraneState,
                         circuit_params: &mut CircuitParameters,
                         dt: f64) -> Result<()> {
        // Update state tracker
        self.state_tracker.update_state(membrane_state, circuit_params);
        
        // Process each feedback loop
        for (loop_id, feedback_loop) in &mut self.feedback_loops {
            self.process_feedback_loop(loop_id, feedback_loop, circuit_params, dt)?;
        }
        
        // Update system statistics
        self.update_statistics();
        
        Ok(())
    }
    
    /// Process a single feedback loop
    fn process_feedback_loop(&mut self, loop_id: &str, feedback_loop: &mut FeedbackLoop,
                           circuit_params: &mut CircuitParameters, dt: f64) -> Result<()> {
        // Get input value
        let input_value = self.get_parameter_value(&feedback_loop.input_parameter)?;
        
        // Calculate feedback signal
        let feedback_signal = self.calculate_feedback_signal(feedback_loop, input_value, dt)?;
        
        // Apply feedback to output parameter
        self.apply_feedback_to_parameter(&feedback_loop.output_parameter, 
                                       feedback_signal, circuit_params)?;
        
        // Update feedback loop state
        feedback_loop.current_value = feedback_signal;
        
        // Update integration history for adaptive feedback
        if feedback_loop.feedback_type == FeedbackType::Adaptive {
            feedback_loop.integration_history.push(feedback_signal);
            if feedback_loop.integration_history.len() > self.control_params.history_length {
                feedback_loop.integration_history.remove(0);
            }
        }
        
        Ok(())
    }
    
    /// Calculate feedback signal based on feedback type
    fn calculate_feedback_signal(&self, feedback_loop: &FeedbackLoop, 
                               input_value: f64, dt: f64) -> Result<f64> {
        match feedback_loop.feedback_type {
            FeedbackType::Negative => {
                // Simple proportional negative feedback
                Ok(-feedback_loop.gain * input_value)
            },
            
            FeedbackType::Positive => {
                // Simple proportional positive feedback
                Ok(feedback_loop.gain * input_value)
            },
            
            FeedbackType::Homeostatic => {
                // PID-style homeostatic control
                if let Some(setpoint) = feedback_loop.setpoint {
                    let error = input_value - setpoint;
                    let proportional = feedback_loop.gain * error;
                    
                    // Add integral term (simplified)
                    let integral = if !feedback_loop.integration_history.is_empty() {
                        let sum: f64 = feedback_loop.integration_history.iter().sum();
                        0.1 * feedback_loop.gain * sum * dt
                    } else {
                        0.0
                    };
                    
                    // Add derivative term (simplified)
                    let derivative = if feedback_loop.integration_history.len() >= 2 {
                        let last_idx = feedback_loop.integration_history.len() - 1;
                        let d_error = feedback_loop.integration_history[last_idx] - 
                                     feedback_loop.integration_history[last_idx - 1];
                        0.05 * feedback_loop.gain * d_error / dt
                    } else {
                        0.0
                    };
                    
                    Ok(-(proportional + integral + derivative))
                } else {
                    Err(MembraneError::InvalidParameter(
                        "Homeostatic feedback requires setpoint".to_string()
                    ))
                }
            },
            
            FeedbackType::Adaptive => {
                // Adaptive feedback with learning
                let base_signal = feedback_loop.gain * input_value;
                
                // Adaptation based on history
                let adaptation_factor = if !feedback_loop.integration_history.is_empty() {
                    let variance = self.calculate_variance(&feedback_loop.integration_history);
                    1.0 + 0.1 * variance // Increase gain if high variance
                } else {
                    1.0
                };
                
                Ok(base_signal * adaptation_factor)
            },
        }
    }
    
    /// Calculate variance of a vector
    fn calculate_variance(&self, values: &[f64]) -> f64 {
        if values.len() < 2 {
            return 0.0;
        }
        
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        let variance = values.iter()
            .map(|x| (x - mean).powi(2))
            .sum::<f64>() / values.len() as f64;
        
        variance
    }
    
    /// Get parameter value from system state
    fn get_parameter_value(&self, parameter_name: &str) -> Result<f64> {
        match parameter_name {
            "membrane_voltage" => Ok(self.state_tracker.membrane_voltage),
            "atp_concentration" => Ok(self.state_tracker.atp_concentration),
            param if param.starts_with("channel_") => {
                let channel_name = &param[8..]; // Remove "channel_" prefix
                self.state_tracker.channel_activities.get(channel_name)
                    .copied()
                    .ok_or_else(|| MembraneError::InvalidParameter(
                        format!("Channel activity {} not found", channel_name)
                    ))
            },
            param if param.starts_with("protein_") => {
                let protein_name = &param[8..]; // Remove "protein_" prefix
                self.state_tracker.protein_states.get(protein_name)
                    .copied()
                    .ok_or_else(|| MembraneError::InvalidParameter(
                        format!("Protein state {} not found", protein_name)
                    ))
            },
            param if param.starts_with("circuit_") => {
                let circuit_param = &param[8..]; // Remove "circuit_" prefix
                self.state_tracker.circuit_params.get(circuit_param)
                    .copied()
                    .ok_or_else(|| MembraneError::InvalidParameter(
                        format!("Circuit parameter {} not found", circuit_param)
                    ))
            },
            _ => Err(MembraneError::InvalidParameter(
                format!("Unknown parameter: {}", parameter_name)
            )),
        }
    }
    
    /// Apply feedback signal to output parameter
    fn apply_feedback_to_parameter(&mut self, parameter_name: &str, 
                                 feedback_signal: f64, 
                                 circuit_params: &mut CircuitParameters) -> Result<()> {
        match parameter_name {
            "circuit_capacitance" => {
                circuit_params.capacitance += feedback_signal;
                circuit_params.capacitance = circuit_params.capacitance.max(1e-15); // Minimum capacitance
            },
            "circuit_resistance" => {
                circuit_params.resistance += feedback_signal;
                circuit_params.resistance = circuit_params.resistance.max(1e3); // Minimum resistance
            },
            param if param.starts_with("voltage_source_") => {
                if let Ok(index) = param[15..].parse::<usize>() {
                    if index < circuit_params.voltage_sources.len() {
                        circuit_params.voltage_sources[index] += feedback_signal;
                    }
                }
            },
            param if param.starts_with("current_source_") => {
                if let Ok(index) = param[15..].parse::<usize>() {
                    if index < circuit_params.current_sources.len() {
                        circuit_params.current_sources[index] += feedback_signal;
                    }
                }
            },
            _ => {
                return Err(MembraneError::InvalidParameter(
                    format!("Cannot apply feedback to parameter: {}", parameter_name)
                ));
            },
        }
        
        Ok(())
    }
    
    /// Create common feedback loops for membrane-circuit coupling
    pub fn create_standard_feedback_loops(&mut self) -> Result<()> {
        // Voltage-dependent capacitance feedback
        let voltage_capacitance_loop = FeedbackLoop {
            id: "voltage_capacitance".to_string(),
            feedback_type: FeedbackType::Negative,
            input_parameter: "membrane_voltage".to_string(),
            output_parameter: "circuit_capacitance".to_string(),
            gain: 1e-12, // Small gain to avoid instability
            time_constant: 0.001, // 1 ms response time
            current_value: 0.0,
            setpoint: None,
            integration_history: Vec::new(),
        };
        
        self.add_feedback_loop(voltage_capacitance_loop)?;
        
        // ATP-dependent resistance feedback
        let atp_resistance_loop = FeedbackLoop {
            id: "atp_resistance".to_string(),
            feedback_type: FeedbackType::Homeostatic,
            input_parameter: "atp_concentration".to_string(),
            output_parameter: "circuit_resistance".to_string(),
            gain: 1e6, // Resistance change per ATP unit
            time_constant: 0.1, // 100 ms response time
            current_value: 0.0,
            setpoint: Some(PHYSIOLOGICAL_ATP), // Maintain physiological ATP
            integration_history: Vec::new(),
        };
        
        self.add_feedback_loop(atp_resistance_loop)?;
        
        // Activity-dependent adaptation
        let activity_adaptation_loop = FeedbackLoop {
            id: "activity_adaptation".to_string(),
            feedback_type: FeedbackType::Adaptive,
            input_parameter: "membrane_voltage".to_string(),
            output_parameter: "voltage_source_0".to_string(),
            gain: 0.1,
            time_constant: 10.0, // 10 s adaptation time
            current_value: 0.0,
            setpoint: None,
            integration_history: Vec::new(),
        };
        
        self.add_feedback_loop(activity_adaptation_loop)?;
        
        Ok(())
    }
    
    /// Update system statistics
    fn update_statistics(&mut self) {
        self.statistics.active_loops = self.feedback_loops.len();
        
        if !self.feedback_loops.is_empty() {
            // Calculate average feedback strength
            let total_strength: f64 = self.feedback_loops.values()
                .map(|loop_| loop_.current_value.abs())
                .sum();
            self.statistics.avg_feedback_strength = total_strength / self.feedback_loops.len() as f64;
            
            // Estimate stability index (simplified)
            let gains: Vec<f64> = self.feedback_loops.values()
                .map(|loop_| loop_.gain.abs())
                .collect();
            let max_gain = gains.iter().copied().fold(0.0, f64::max);
            self.statistics.stability_index = (1.0 / (1.0 + max_gain)).min(1.0);
        }
    }
    
    /// Check system stability
    pub fn check_stability(&self) -> bool {
        // Simple stability check based on feedback gains
        let total_positive_gain: f64 = self.feedback_loops.values()
            .filter(|loop_| loop_.feedback_type == FeedbackType::Positive)
            .map(|loop_| loop_.gain.abs())
            .sum();
        
        let total_negative_gain: f64 = self.feedback_loops.values()
            .filter(|loop_| loop_.feedback_type == FeedbackType::Negative || 
                           loop_.feedback_type == FeedbackType::Homeostatic)
            .map(|loop_| loop_.gain.abs())
            .sum();
        
        // System is stable if negative feedback dominates
        total_negative_gain > total_positive_gain
    }
    
    /// Get feedback loop by ID
    pub fn get_feedback_loop(&self, loop_id: &str) -> Option<&FeedbackLoop> {
        self.feedback_loops.get(loop_id)
    }
    
    /// Remove feedback loop
    pub fn remove_feedback_loop(&mut self, loop_id: &str) -> Option<FeedbackLoop> {
        self.feedback_loops.remove(loop_id)
    }
}

impl StateTracker {
    /// Create a new state tracker
    pub fn new() -> Self {
        Self {
            membrane_voltage: 0.0,
            atp_concentration: PHYSIOLOGICAL_ATP,
            channel_activities: HashMap::new(),
            protein_states: HashMap::new(),
            circuit_params: HashMap::new(),
            history: HashMap::new(),
        }
    }
    
    /// Update state from membrane and circuit parameters
    pub fn update_state(&mut self, membrane_state: &MembraneState, 
                       circuit_params: &CircuitParameters) {
        // Update membrane state
        self.membrane_voltage = membrane_state.voltage;
        self.atp_concentration = membrane_state.atp.concentration;
        
        // Update circuit parameters
        self.circuit_params.insert("capacitance".to_string(), circuit_params.capacitance);
        self.circuit_params.insert("resistance".to_string(), circuit_params.resistance);
        
        // Store history for trend analysis
        self.store_history("membrane_voltage", self.membrane_voltage);
        self.store_history("atp_concentration", self.atp_concentration);
    }
    
    /// Store parameter history
    fn store_history(&mut self, parameter: &str, value: f64) {
        let history = self.history.entry(parameter.to_string()).or_insert_with(Vec::new);
        history.push(value);
        
        // Keep only recent history
        if history.len() > 1000 {
            history.remove(0);
        }
    }
}

impl Default for ControlParameters {
    fn default() -> Self {
        Self {
            max_gain: 100.0,
            min_time_constant: 0.001, // 1 ms
            max_time_constant: 100.0, // 100 s
            history_length: 100,
            update_rate: 1000.0, // 1 kHz
        }
    }
}

impl FeedbackStatistics {
    pub fn new() -> Self {
        Self {
            active_loops: 0,
            avg_feedback_strength: 0.0,
            stability_index: 1.0,
            adaptation_rate: 0.0,
            total_energy_consumed: 0.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_feedback_system_creation() {
        let system = FeedbackSystem::new();
        assert!(system.feedback_loops.is_empty());
        assert_eq!(system.statistics.active_loops, 0);
    }
    
    #[test]
    fn test_feedback_loop_validation() {
        let mut system = FeedbackSystem::new();
        
        let valid_loop = FeedbackLoop {
            id: "test_loop".to_string(),
            feedback_type: FeedbackType::Negative,
            input_parameter: "membrane_voltage".to_string(),
            output_parameter: "circuit_capacitance".to_string(),
            gain: 1.0,
            time_constant: 1.0,
            current_value: 0.0,
            setpoint: None,
            integration_history: Vec::new(),
        };
        
        assert!(system.add_feedback_loop(valid_loop).is_ok());
    }
    
    #[test]
    fn test_stability_check() {
        let mut system = FeedbackSystem::new();
        
        // Add a stabilizing negative feedback loop
        let negative_loop = FeedbackLoop {
            id: "stabilizing".to_string(),
            feedback_type: FeedbackType::Negative,
            input_parameter: "membrane_voltage".to_string(),
            output_parameter: "circuit_capacitance".to_string(),
            gain: 2.0,
            time_constant: 1.0,
            current_value: 0.0,
            setpoint: None,
            integration_history: Vec::new(),
        };
        
        system.add_feedback_loop(negative_loop).unwrap();
        assert!(system.check_stability());
    }
}
