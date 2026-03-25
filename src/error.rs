//! Error handling for the Bene Gesserit membrane dynamics system
//! 
//! This module provides comprehensive error types for all aspects of membrane
//! simulation, from basic validation to complex biophysical constraint violations.

use std::fmt;
use thiserror::Error;

/// Result type used throughout the membrane dynamics system
pub type Result<T> = std::result::Result<T, MembraneError>;

/// Comprehensive error types for membrane dynamics simulation
#[derive(Error, Debug)]
pub enum MembraneError {
    /// System does not meet minimum requirements for membrane simulation
    #[error("Insufficient memory: required {required}GB, available {available}GB")]
    InsufficientMemory { required: f64, available: f64 },
    
    /// Generic system information error
    #[error("System information error: {0}")]
    SystemInfo(String),
    
    /// Configuration validation errors
    #[error("Invalid configuration: {parameter} = {value}, expected {expected}")]
    InvalidConfiguration {
        parameter: String,
        value: String,
        expected: String,
    },
    
    /// Membrane patch creation and validation errors
    #[error("Invalid membrane patch: {reason}")]
    InvalidMembranePatch { reason: String },
    
    /// ATP pool and energy constraint errors
    #[error("ATP constraint violation: {operation} requires {required} ATP, only {available} available")]
    AtpConstraintViolation {
        operation: String,
        required: f64,
        available: f64,
    },
    
    /// Protein-related errors (invalid types, densities, conformations)
    #[error("Protein error: {protein_type} - {reason}")]
    ProteinError {
        protein_type: String,
        reason: String,
    },
    
    /// Lipid bilayer and membrane composition errors
    #[error("Lipid bilayer error: {reason}")]
    LipidBilayerError { reason: String },
    
    /// Electrochemical gradient and ion distribution errors
    #[error("Electrochemical error: {ion} gradient violation - {details}")]
    ElectrochemicalError { ion: String, details: String },
    
    /// Temperature-dependent process errors
    #[error("Temperature constraint violation: {temperature}K outside viable range [{min}K, {max}K]")]
    TemperatureConstraint {
        temperature: f64,
        min: f64,
        max: f64,
    },
    
    /// Simulation timestep and numerical stability errors
    #[error("Simulation error: {reason} (timestep: {timestep}s)")]
    SimulationError { reason: String, timestep: f64 },
    
    /// Circuit interface and parameter mapping errors
    #[error("Circuit interface error: {component} mapping failed - {reason}")]
    CircuitInterfaceError { component: String, reason: String },
    
    /// External system integration errors
    #[error("External system error: {system} - {error}")]
    ExternalSystemError { system: String, error: String },
    
    /// Orchestrator communication errors
    #[error("Orchestrator communication error: {operation} failed - {reason}")]
    OrchestratorError { operation: String, reason: String },
    
    /// Nebuchadnezzar circuit system errors
    #[error("Nebuchadnezzar error: {reason}")]
    NebuchadnezzarError { reason: String },
    
    /// Biological constraint violations (non-physical states)
    #[error("Biological constraint violation: {constraint} - {details}")]
    BiologicalConstraint { constraint: String, details: String },
    
    /// Numerical computation errors (convergence, precision, overflow)
    #[error("Numerical error: {computation} failed - {reason}")]
    NumericalError { computation: String, reason: String },
    
    /// Parallel processing and concurrency errors
    #[error("Parallel processing error: {reason}")]
    ParallelProcessingError { reason: String },
    
    /// File I/O and serialization errors
    #[error("I/O error: {operation} - {source}")]
    IoError {
        operation: String,
        #[source]
        source: std::io::Error,
    },
    
    /// JSON serialization/deserialization errors
    #[error("Serialization error: {context} - {source}")]
    SerializationError {
        context: String,
        #[source]
        source: serde_json::Error,
    },
    
    /// HTTP/Network communication errors
    #[error("Network error: {operation} - {source}")]
    NetworkError {
        operation: String,
        #[source]
        source: reqwest::Error,
    },
    
    /// WebSocket communication errors
    #[error("WebSocket error: {operation} - {reason}")]
    WebSocketError { operation: String, reason: String },
    
    /// Parse errors for numeric values
    #[error("Parse error: cannot parse '{value}' as {expected_type}")]
    ParseError {
        value: String,
        expected_type: String,
        #[source]
        source: std::num::ParseIntError,
    },
    
    /// Generic validation errors
    #[error("Validation error: {field} - {reason}")]
    ValidationError { field: String, reason: String },
    
    /// Environment container errors
    #[error("Environment error: {reason}")]
    Environment { reason: String },
}

/// Environment container specific errors
#[derive(Error, Debug)]
pub enum EnvironmentError {
    /// Energy charge constraint violation
    #[error("Energy charge violation: outside acceptable range")]
    EnergyChargeViolation,
    
    /// Temperature constraint violation
    #[error("Temperature constraint violation: {temperature}K outside viable range")]
    TemperatureViolation { temperature: f64 },
    
    /// pH constraint violation
    #[error("pH constraint violation: {ph} outside physiological range")]
    PhViolation { ph: f64 },
    
    /// Osmolarity constraint violation
    #[error("Osmolarity constraint violation: {osmolarity} mOsm/L outside viable range")]
    OsmolarityViolation { osmolarity: f64 },
    
    /// General environment constraint violation
    #[error("Environment constraint violation: {constraint}")]
    ConstraintViolation { constraint: String },
    
    /// Environment interaction processing error
    #[error("Environment interaction error: {reason}")]
    InteractionError { reason: String },
}

impl MembraneError {
    /// Create a new ATP constraint violation error
    pub fn atp_insufficient(operation: &str, required: f64, available: f64) -> Self {
        Self::AtpConstraintViolation {
            operation: operation.to_string(),
            required,
            available,
        }
    }
    
    /// Create a new protein error
    pub fn protein_error(protein_type: &str, reason: &str) -> Self {
        Self::ProteinError {
            protein_type: protein_type.to_string(),
            reason: reason.to_string(),
        }
    }
    
    /// Create a new biological constraint violation
    pub fn biological_violation(constraint: &str, details: &str) -> Self {
        Self::BiologicalConstraint {
            constraint: constraint.to_string(),
            details: details.to_string(),
        }
    }
    
    /// Create a new simulation error
    pub fn simulation_error(reason: &str, timestep: f64) -> Self {
        Self::SimulationError {
            reason: reason.to_string(),
            timestep,
        }
    }
    
    /// Create a new circuit interface error
    pub fn circuit_interface_error(component: &str, reason: &str) -> Self {
        Self::CircuitInterfaceError {
            component: component.to_string(),
            reason: reason.to_string(),
        }
    }
    
    /// Check if this error is recoverable (can retry the operation)
    pub fn is_recoverable(&self) -> bool {
        match self {
            // Recoverable errors - temporary conditions
            Self::AtpConstraintViolation { .. } => true,
            Self::NetworkError { .. } => true,
            Self::WebSocketError { .. } => true,
            Self::OrchestratorError { .. } => true,
            Self::ParallelProcessingError { .. } => true,
            
            // Non-recoverable errors - fundamental problems
            Self::InsufficientMemory { .. } => false,
            Self::InvalidConfiguration { .. } => false,
            Self::BiologicalConstraint { .. } => false,
            Self::TemperatureConstraint { .. } => false,
            Self::ValidationError { .. } => false,
            
            // Context-dependent
            _ => false,
        }
    }
    
    /// Get the error category for logging and metrics
    pub fn category(&self) -> ErrorCategory {
        match self {
            Self::InsufficientMemory { .. } | Self::SystemInfo(_) => ErrorCategory::System,
            Self::InvalidConfiguration { .. } | Self::ValidationError { .. } => ErrorCategory::Configuration,
            Self::AtpConstraintViolation { .. } => ErrorCategory::Energy,
            Self::ProteinError { .. } | Self::LipidBilayerError { .. } => ErrorCategory::Biochemical,
            Self::ElectrochemicalError { .. } => ErrorCategory::Electrochemical,
            Self::TemperatureConstraint { .. } => ErrorCategory::Physical,
            Self::SimulationError { .. } | Self::NumericalError { .. } => ErrorCategory::Numerical,
            Self::CircuitInterfaceError { .. } => ErrorCategory::Circuit,
            Self::OrchestratorError { .. } | Self::ExternalSystemError { .. } => ErrorCategory::External,
            Self::BiologicalConstraint { .. } => ErrorCategory::Biological,
            Self::IoError { .. } | Self::NetworkError { .. } | Self::WebSocketError { .. } => ErrorCategory::Communication,
            _ => ErrorCategory::Other,
        }
    }
}

/// Error categories for classification and handling
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ErrorCategory {
    System,
    Configuration,
    Energy,
    Biochemical,
    Electrochemical,
    Physical,
    Numerical,
    Circuit,
    External,
    Biological,
    Communication,
    Other,
}

impl fmt::Display for ErrorCategory {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::System => write!(f, "System"),
            Self::Configuration => write!(f, "Configuration"),
            Self::Energy => write!(f, "Energy"),
            Self::Biochemical => write!(f, "Biochemical"),
            Self::Electrochemical => write!(f, "Electrochemical"),
            Self::Physical => write!(f, "Physical"),
            Self::Numerical => write!(f, "Numerical"),
            Self::Circuit => write!(f, "Circuit"),
            Self::External => write!(f, "External"),
            Self::Biological => write!(f, "Biological"),
            Self::Communication => write!(f, "Communication"),
            Self::Other => write!(f, "Other"),
        }
    }
}

// Implement conversions from common error types
impl From<std::io::Error> for MembraneError {
    fn from(err: std::io::Error) -> Self {
        Self::IoError {
            operation: "file operation".to_string(),
            source: err,
        }
    }
}

impl From<serde_json::Error> for MembraneError {
    fn from(err: serde_json::Error) -> Self {
        Self::SerializationError {
            context: "JSON processing".to_string(),
            source: err,
        }
    }
}

impl From<reqwest::Error> for MembraneError {
    fn from(err: reqwest::Error) -> Self {
        Self::NetworkError {
            operation: "HTTP request".to_string(),
            source: err,
        }
    }
}

impl From<std::num::ParseIntError> for MembraneError {
    fn from(err: std::num::ParseIntError) -> Self {
        Self::ParseError {
            value: "unknown".to_string(),
            expected_type: "integer".to_string(),
            source: err,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_error_categories() {
        let atp_error = MembraneError::atp_insufficient("pump", 100.0, 50.0);
        assert_eq!(atp_error.category(), ErrorCategory::Energy);
        
        let protein_error = MembraneError::protein_error("NaKATPase", "invalid conformation");
        assert_eq!(protein_error.category(), ErrorCategory::Biochemical);
        
        let bio_error = MembraneError::biological_violation("membrane potential", "exceeds physiological range");
        assert_eq!(bio_error.category(), ErrorCategory::Biological);
    }
    
    #[test]
    fn test_error_recoverability() {
        let atp_error = MembraneError::atp_insufficient("pump", 100.0, 50.0);
        assert!(atp_error.is_recoverable());
        
        let config_error = MembraneError::InvalidConfiguration {
            parameter: "temperature".to_string(),
            value: "-100".to_string(),
            expected: "> 0".to_string(),
        };
        assert!(!config_error.is_recoverable());
    }
    
    #[test]
    fn test_error_display() {
        let error = MembraneError::atp_insufficient("Na/K pump", 150.0, 100.0);
        let error_str = format!("{}", error);
        assert!(error_str.contains("Na/K pump"));
        assert!(error_str.contains("150"));
        assert!(error_str.contains("100"));
    }
} 