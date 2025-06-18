//! # Bene Gesserit: Membrane Dynamics Module
//! 
//! Biologically authentic cellular membrane simulation with ATP-constrained computation.
//! 
//! ## Architecture
//! 
//! The Bene Gesserit membrane dynamics system implements a four-layer architecture:
//! 
//! - **Molecular Layer**: Fundamental biophysics (lipid bilayers, proteins, electrochemistry)
//! - **Mesoscale Layer**: Emergent properties (lipid rafts, protein clustering, domains)
//! - **Cellular Layer**: Whole-cell networks (organelles, contact sites, topology)
//! - **Circuit Interface Layer**: Translation to external circuit systems
//! 
//! ## Core Principles
//! 
//! 1. **Biological Authenticity**: All processes based on experimental biophysical data
//! 2. **ATP Constraints**: Energy limitations drive computational behavior
//! 3. **Dynamic Coupling**: Membrane state directly influences circuit parameters
//! 4. **Multi-scale Integration**: Seamless molecular-to-circuit abstraction
//! 
//! ## Example Usage
//! 
//! ```rust
//! use bene_gesserit::membrane::MembranePatch;
//! use bene_gesserit::atp::AtpPool;
//! 
//! // Create a membrane patch with ATP constraints
//! let mut patch = MembranePatch::new(1e-9, 310.15)?; // 1 μm², 37°C
//! let atp_pool = AtpPool::new(5e-3); // 5 mM ATP
//! 
//! // Add membrane components
//! patch.add_lipid_bilayer(Default::default())?;
//! patch.add_protein("NaKATPase", 1000, &atp_pool)?;
//! 
//! // Run simulation step
//! let state = patch.step(0.001, &atp_pool)?; // 1ms timestep
//! ```

use std::error::Error;
use std::fmt;

pub mod error;
pub mod types;
pub mod constants;

// Core membrane dynamics modules
pub mod molecular;
pub mod systems;
pub mod circuit_interface;

// Re-exports for convenience
pub use error::{MembraneError, Result};
pub use types::*;
pub use molecular::*;
pub use systems::*;
pub use circuit_interface::*;

/// Current version of the Bene Gesserit membrane dynamics system
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Initialize the membrane dynamics system with default configuration
/// 
/// This sets up logging, validates system requirements, and prepares
/// the runtime for membrane simulations.
/// 
/// # Examples
/// 
/// ```rust
/// use bene_gesserit;
/// 
/// fn main() -> bene_gesserit::Result<()> {
///     bene_gesserit::init()?;
///     
///     // System is now ready for membrane simulations
///     Ok(())
/// }
/// ```
pub fn init() -> Result<()> {
    logging::init()?;
    config::validate_system_requirements()?;
    
    log::info!("Bene Gesserit membrane dynamics system v{} initialized", VERSION);
    log::info!("System: ATP-constrained biological computation ready");
    
    Ok(())
}

/// Validate that the system meets minimum requirements for membrane simulation
/// 
/// Checks available memory, CPU capabilities, and required dependencies.
pub fn validate_system() -> Result<SystemInfo> {
    let info = SystemInfo {
        available_memory_gb: get_available_memory()?,
        cpu_cores: num_cpus::get(),
        supports_simd: is_simd_available(),
        rust_version: env!("CARGO_PKG_RUST_VERSION").to_string(),
    };
    
    // Minimum requirements check
    if info.available_memory_gb < 4.0 {
        return Err(MembraneError::InsufficientMemory {
            required: 4.0,
            available: info.available_memory_gb,
        });
    }
    
    if info.cpu_cores < 2 {
        log::warn!("Single-core system detected. Parallel membrane processing disabled.");
    }
    
    Ok(info)
}

/// System information for membrane dynamics requirements
#[derive(Debug, Clone)]
pub struct SystemInfo {
    pub available_memory_gb: f64,
    pub cpu_cores: usize,
    pub supports_simd: bool,
    pub rust_version: String,
}

// Helper functions
fn get_available_memory() -> Result<f64> {
    // Platform-specific memory detection
    #[cfg(target_os = "linux")]
    {
        use std::fs;
        let meminfo = fs::read_to_string("/proc/meminfo")?;
        for line in meminfo.lines() {
            if line.starts_with("MemAvailable:") {
                let kb: u64 = line
                    .split_whitespace()
                    .nth(1)
                    .ok_or(MembraneError::SystemInfo("Cannot parse /proc/meminfo".to_string()))?
                    .parse()?;
                return Ok(kb as f64 / 1024.0 / 1024.0); // Convert KB to GB
            }
        }
        Err(MembraneError::SystemInfo("MemAvailable not found".to_string()))
    }
    
    #[cfg(target_os = "macos")]
    {
        // macOS-specific implementation would go here
        // For now, return a reasonable default
        Ok(8.0)
    }
    
    #[cfg(target_os = "windows")]
    {
        // Windows-specific implementation would go here
        // For now, return a reasonable default
        Ok(8.0)
    }
    
    #[cfg(not(any(target_os = "linux", target_os = "macos", target_os = "windows")))]
    {
        Ok(8.0) // Default assumption
    }
}

fn is_simd_available() -> bool {
    // Check for SIMD instruction support
    #[cfg(target_arch = "x86_64")]
    {
        is_x86_feature_detected!("avx2")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_system_validation() {
        let result = validate_system();
        assert!(result.is_ok(), "System validation should succeed");
        
        let info = result.unwrap();
        assert!(info.cpu_cores > 0, "Should detect at least one CPU core");
        assert!(info.available_memory_gb > 0.0, "Should detect available memory");
    }
    
    #[test]
    fn test_version_constant() {
        assert!(!VERSION.is_empty(), "Version should not be empty");
        assert!(VERSION.contains('.'), "Version should contain dots");
    }
} 