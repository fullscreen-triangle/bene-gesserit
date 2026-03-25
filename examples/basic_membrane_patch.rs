//! Basic Membrane Patch Example
//! 
//! Demonstrates the core functionality of the Bene Gesserit membrane dynamics system:
//! - Creating a molecular membrane patch
//! - Adding ATP-driven proteins
//! - Running ATP-constrained simulation
//! - Extracting circuit parameters for Nebuchadnezzar integration

use bene_gesserit::{
    molecular::MolecularMembrane,
    circuit_interface::{CircuitInterface, CircuitInterfaceConfig},
    types::{MembraneConfig, ProteinType},
    constants::PHYSIOLOGICAL_TEMPERATURE,
    Result,
};

fn main() -> Result<()> {
    // Initialize logging
    env_logger::init();
    
    println!("🧬 Bene Gesserit Membrane Dynamics Example");
    println!("=========================================");
    
    // Create membrane configuration
    let mut config = MembraneConfig::default();
    config.area = 1e-9;  // 1 μm² patch
    config.temperature = PHYSIOLOGICAL_TEMPERATURE;  // 37°C
    config.initial_atp = 5e-3;  // 5 mM ATP
    config.atp_constrained = true;  // Enable ATP constraints
    config.max_time = 0.01;  // 10 ms simulation
    config.timestep = 1e-6;  // 1 μs steps
    
    println!("📐 Membrane Configuration:");
    println!("   Area: {:.1e} m² ({:.1} μm²)", config.area, config.area * 1e12);
    println!("   Temperature: {:.1} K ({:.1}°C)", config.temperature, config.temperature - 273.15);
    println!("   Initial ATP: {:.1} mM", config.initial_atp * 1000.0);
    println!("   ATP Constrained: {}", config.atp_constrained);
    println!();
    
    // Create molecular membrane
    let mut membrane = MolecularMembrane::new(config.clone())?;
    println!("🔬 Created molecular membrane patch");
    
    // Add membrane proteins
    membrane.add_protein("pump1".to_string(), ProteinType::NaKATPase, (0.0, 0.0))?;
    membrane.add_protein("pump2".to_string(), ProteinType::CaATPase, (5.0, 0.0))?;
    membrane.add_protein("channel1".to_string(), ProteinType::VGSC, (0.0, 5.0))?;
    membrane.add_protein("channel2".to_string(), ProteinType::VGKC, (5.0, 5.0))?;
    
    println!("🧪 Added membrane proteins:");
    println!("   - Na⁺/K⁺-ATPase pump (ATP-driven)");
    println!("   - Ca²⁺-ATPase pump (ATP-driven)");
    println!("   - Voltage-gated Na⁺ channel");
    println!("   - Voltage-gated K⁺ channel");
    println!();
    
    // Create circuit interface for Nebuchadnezzar integration
    let interface_config = CircuitInterfaceConfig::default();
    let mut circuit_interface = CircuitInterface::new(membrane, interface_config)?;
    
    println!("⚡ Created circuit interface with {} hierarchy levels", 
             circuit_interface.circuit_state.levels.len());
    println!();
    
    // Run simulation
    println!("🚀 Running ATP-constrained membrane simulation...");
    println!("Time (ms) | ATP (mM) | Voltage (mV) | Current (pA) | Power (pW)");
    println!("----------|----------|--------------|--------------|----------");
    
    let mut time = 0.0;
    let dt = config.timestep;
    let max_steps = (config.max_time / dt) as usize;
    let print_interval = max_steps / 20;  // Print 20 updates
    
    for step in 0..max_steps {
        // Update circuit interface (includes membrane step)
        circuit_interface.update(dt)?;
        
        time += dt;
        
        // Print progress
        if step % print_interval == 0 {
            let membrane_state = &circuit_interface.membrane.state;
            let atp_mm = membrane_state.atp.concentration * 1000.0;  // Convert to mM
            let voltage_mv = membrane_state.voltage * 1000.0;  // Convert to mV
            let current_pa = membrane_state.current * 1e12;  // Convert to pA
            let power_pw = circuit_interface.atp_dynamics.atp_consumption_rate * 
                          bene_gesserit::constants::atp::ATP_ENERGY_PER_MOLECULE.abs() * 1e12;  // Convert to pW
            
            println!("{:8.2} | {:8.2} | {:11.1} | {:11.1} | {:9.1}",
                     time * 1000.0, atp_mm, voltage_mv, current_pa, power_pw);
        }
    }
    
    println!();
    println!("✅ Simulation completed!");
    println!();
    
    // Extract final results
    let final_state = &circuit_interface.membrane.state;
    let atp_tracker = &circuit_interface.membrane.atp_tracker;
    
    println!("📊 Final Membrane State:");
    println!("   Voltage: {:.1} mV", final_state.voltage * 1000.0);
    println!("   Current: {:.1} pA", final_state.current * 1e12);
    println!("   ATP: {:.2} mM", final_state.atp.concentration * 1000.0);
    println!("   Energy consumed: {:.2} fJ", final_state.energy_consumed * 1e15);
    println!();
    
    println!("⚡ Energy Efficiency:");
    println!("   ATP efficiency: {:.1}%", atp_tracker.efficiency.atp_efficiency * 100.0);
    println!("   Coupling efficiency: {:.1}%", atp_tracker.efficiency.coupling_efficiency * 100.0);
    println!("   Waste heat: {:.1}%", atp_tracker.efficiency.waste_heat_fraction * 100.0);
    println!();
    
    // Get Nebuchadnezzar integration parameters
    let neb_params = circuit_interface.get_nebuchadnezzar_parameters();
    
    println!("🔗 Nebuchadnezzar Integration Parameters:");
    println!("   Hierarchy levels: {}", neb_params.hierarchy_levels);
    println!("   Circuit nodes: {}", neb_params.circuit_nodes.len());
    println!("   Circuit edges: {}", neb_params.circuit_edges.len());
    println!("   Global impedance: {:.1e} Ω", neb_params.global_impedance);
    println!("   Stability measure: {:.3}", neb_params.stability_measure);
    println!();
    
    // Show circuit hierarchy
    println!("🏗️  Circuit Hierarchy:");
    for (level_idx, level) in circuit_interface.circuit_state.levels.iter().enumerate() {
        let level_name = match level_idx {
            0 => "Molecular",
            1 => "Mesoscale", 
            2 => "Cellular",
            _ => "Unknown",
        };
        
        println!("   Level {}: {} ({} nodes, {:.1e} m scale)",
                 level_idx, level_name, level.nodes.len(), level.properties.spatial_scale);
        
        // Show ATP constraints
        println!("     ATP allocation: {:.1}%", level.atp_constraints.atp_allocation * 100.0);
        println!("     ATP efficiency: {:.1}%", level.atp_constraints.atp_efficiency * 100.0);
    }
    println!();
    
    // Show protein activities
    println!("🧪 Protein Activities:");
    for (protein_id, protein) in &circuit_interface.membrane.proteins {
        println!("   {}: {:.1}% active ({:?})", 
                 protein_id, 
                 protein.activity * 100.0,
                 protein.conformation);
    }
    println!();
    
    println!("🎯 Example completed successfully!");
    println!("The membrane dynamics system is ready for integration with your Nebuchadnezzar");
    println!("hierarchical probabilistic electric circuit system using ATP-based differential");
    println!("equations (dx/dATP instead of dx/dt).");
    
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_basic_membrane_patch() {
        let result = main();
        assert!(result.is_ok(), "Basic membrane patch example should run successfully");
    }
} 