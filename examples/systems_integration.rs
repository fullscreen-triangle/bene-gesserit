//! Systems Integration Example
//! 
//! Demonstrates the complete 4-level hierarchy of the Bene Gesserit membrane dynamics:
//! - Level 0: Molecular (proteins, lipids)
//! - Level 1: Mesoscale (domains, clusters) 
//! - Level 2: Cellular (membrane patches)
//! - Level 3: Systems (O₂ uptake, temperature, ionic gradients, fluid flow)

use bene_gesserit::{
    molecular::MolecularMembrane,
    systems::{SystemsLevel, EnvironmentalConditions},
    circuit_interface::{CircuitInterface, CircuitInterfaceConfig},
    types::{MembraneConfig, ProteinType},
    constants::PHYSIOLOGICAL_TEMPERATURE,
    Result,
};

fn main() -> Result<()> {
    env_logger::init();
    
    println!("🌟 Bene Gesserit Systems Integration Example");
    println!("===========================================");
    println!();
    
    // Environmental conditions
    let environment = EnvironmentalConditions {
        ambient_temperature: 298.15,     // 25°C room temp
        atmospheric_pressure: 101325.0,  // 1 atm
        oxygen_partial_pressure: 21000.0, // 21% O₂
        co2_partial_pressure: 40.0,      // 0.04% CO₂
        humidity: 0.6,                   // 60% RH
        external_ion_concentrations: std::collections::HashMap::new(),
    };
    
    println!("🌍 Environmental Conditions:");
    println!("   Temperature: {:.1}°C", environment.ambient_temperature - 273.15);
    println!("   O₂ pressure: {:.0} Pa ({:.1}%)", environment.oxygen_partial_pressure, 
             environment.oxygen_partial_pressure / 1013.25);
    println!("   CO₂ pressure: {:.1} Pa", environment.co2_partial_pressure);
    println!("   Humidity: {:.0}%", environment.humidity * 100.0);
    println!();
    
    // Create membrane configuration
    let mut config = MembraneConfig::default();
    config.area = 1e-9;  // 1 μm² patch
    config.temperature = PHYSIOLOGICAL_TEMPERATURE;
    config.initial_atp = 5e-3;  // 5 mM ATP
    config.atp_constrained = true;
    config.max_time = 1.0;  // 1 second simulation
    config.timestep = 1e-3;  // 1 ms steps
    
    // Create molecular membrane with proteins
    let mut membrane = MolecularMembrane::new(config.clone())?;
    membrane.add_protein("na_k_pump_1".to_string(), ProteinType::NaKATPase, (0.0, 0.0))?;
    membrane.add_protein("na_k_pump_2".to_string(), ProteinType::NaKATPase, (10.0, 0.0))?;
    membrane.add_protein("ca_pump".to_string(), ProteinType::CaATPase, (5.0, 5.0))?;
    membrane.add_protein("na_channel".to_string(), ProteinType::VGSC, (0.0, 10.0))?;
    membrane.add_protein("k_channel".to_string(), ProteinType::VGKC, (10.0, 10.0))?;
    
    println!("🔬 Molecular Level Setup:");
    println!("   Membrane area: {:.1} μm²", config.area * 1e12);
    println!("   Proteins: {} total", membrane.proteins.len());
    println!("   - 2x Na⁺/K⁺-ATPase pumps (ATP-driven)");
    println!("   - 1x Ca²⁺-ATPase pump (ATP-driven)");
    println!("   - 1x Voltage-gated Na⁺ channel");
    println!("   - 1x Voltage-gated K⁺ channel");
    println!();
    
    // Create circuit interface with 4 hierarchy levels
    let interface_config = CircuitInterfaceConfig {
        hierarchy_levels: 4,  // Include systems level
        atp_constraint_strength: 1.0,
        probabilistic_update_freq: 10.0,  // 10 Hz for systems level
        circuit_resolution: 1e6,
        nebuchadnezzar_integration: true,
    };
    
    let mut circuit_interface = CircuitInterface::new(membrane, interface_config)?;
    
    println!("⚡ Circuit Interface Created:");
    println!("   Hierarchy levels: {}", circuit_interface.circuit_state.levels.len());
    for (i, level) in circuit_interface.circuit_state.levels.iter().enumerate() {
        let level_name = match i {
            0 => "Molecular",
            1 => "Mesoscale",
            2 => "Cellular", 
            3 => "Systems",
            _ => "Unknown",
        };
        println!("   Level {}: {} ({} nodes, {:.1e} m scale, {:.1e} s time)",
                 i, level_name, level.nodes.len(), 
                 level.properties.spatial_scale, level.properties.temporal_scale);
    }
    println!();
    
    // Create systems level
    let mut systems = SystemsLevel::new(environment);
    
    println!("🏥 Systems Level Initialized:");
    println!("   Metabolic rate: {:.1} W", systems.metabolism.metabolic_rate);
    println!("   O₂ uptake: {:.2e} mol/s", systems.metabolism.oxygen_uptake);
    println!("   CO₂ production: {:.2e} mol/s", systems.metabolism.co2_production);
    println!("   Core temperature: {:.1}°C", systems.thermal.core_temperature - 273.15);
    println!("   Blood flow: {:.2e} m³/s", systems.fluid_dynamics.blood_flow);
    println!();
    
    // Run integrated simulation
    println!("🚀 Running integrated 4-level simulation...");
    println!();
    println!("Time | ATP  | Voltage | O₂ Uptake | Temp | Blood Flow | Efficiency");
    println!("(s)  | (mM) | (mV)    | (μmol/s)  | (°C) | (mL/s)     | (%)");
    println!("-----|------|---------|-----------|------|------------|----------");
    
    let dt = config.timestep;
    let max_steps = (config.max_time / dt) as usize;
    let print_interval = max_steps / 20;
    
    for step in 0..max_steps {
        let time = step as f64 * dt;
        
        // Update circuit interface (molecular → cellular levels)
        circuit_interface.update(dt)?;
        
        // Update systems level based on circuit interfaces
        systems.update(&[circuit_interface.clone()], dt)?;
        
        // Print progress
        if step % print_interval == 0 {
            let membrane_state = &circuit_interface.membrane.state;
            let atp_mm = membrane_state.atp.concentration * 1000.0;
            let voltage_mv = membrane_state.voltage * 1000.0;
            let o2_uptake_umol = systems.metabolism.oxygen_uptake * 1e6;
            let temp_c = systems.thermal.core_temperature - 273.15;
            let blood_flow_ml = systems.fluid_dynamics.blood_flow * 1e6;
            let efficiency = systems.energy_distribution.system_efficiency * 100.0;
            
            println!("{:4.1} | {:4.1} | {:7.1} | {:9.2} | {:4.1} | {:10.2} | {:8.1}",
                     time, atp_mm, voltage_mv, o2_uptake_umol, temp_c, blood_flow_ml, efficiency);
        }
    }
    
    println!();
    println!("✅ Simulation completed!");
    println!();
    
    // Final systems analysis
    println!("📊 Final Systems State:");
    println!();
    
    println!("🫁 Respiratory System:");
    println!("   O₂ uptake: {:.2e} mol/s ({:.1} mL/min)", 
             systems.metabolism.oxygen_uptake,
             systems.metabolism.oxygen_uptake * 22.4 * 60.0);
    println!("   CO₂ production: {:.2e} mol/s ({:.1} mL/min)",
             systems.metabolism.co2_production,
             systems.metabolism.co2_production * 22.4 * 60.0);
    println!("   Respiratory quotient: {:.2}", systems.metabolism.respiratory_quotient);
    if systems.metabolism.lactate_production > 0.0 {
        println!("   Lactate production: {:.2e} mol/s (anaerobic)", systems.metabolism.lactate_production);
    }
    println!();
    
    println!("🌡️  Thermal System:");
    println!("   Core temperature: {:.2}°C", systems.thermal.core_temperature - 273.15);
    println!("   Local temperature: {:.2}°C", systems.thermal.local_temperature - 273.15);
    println!("   Heat production: {:.1} W", systems.thermal.heat_production);
    println!("   Heat loss: {:.1} W", systems.thermal.heat_loss);
    println!("   Temperature gradient: {:.0} K/m", systems.thermal.temperature_gradient);
    println!();
    
    println!("⚡ Ionic Transport:");
    for (ion, flux) in &systems.ionic_transport.bulk_fluxes {
        println!("   {}: {:.2e} mol/s", ion.symbol(), flux);
    }
    println!("   Total osmolarity: {:.3} osmol/L", systems.ionic_transport.total_osmolarity);
    println!("   Water flux: {:.2e} m³/s", systems.ionic_transport.water_flux);
    println!();
    
    println!("🩸 Fluid Dynamics:");
    println!("   Blood flow: {:.2e} m³/s ({:.1} mL/s)", 
             systems.fluid_dynamics.blood_flow,
             systems.fluid_dynamics.blood_flow * 1e6);
    println!("   Lymphatic flow: {:.2e} m³/s", systems.fluid_dynamics.lymphatic_flow);
    println!("   Capillary pressure: {:.0} Pa ({:.1} mmHg)", 
             systems.fluid_dynamics.capillary_pressure,
             systems.fluid_dynamics.capillary_pressure / 133.322);
    println!("   Interstitial pressure: {:.0} Pa ({:.1} mmHg)",
             systems.fluid_dynamics.interstitial_pressure,
             systems.fluid_dynamics.interstitial_pressure / 133.322);
    println!();
    
    println!("🏠 Homeostasis:");
    println!("   pH: {:.2} (target: {:.2})", 
             systems.homeostasis.ph_control.current_ph,
             systems.homeostasis.ph_control.target_ph);
    println!("   Free Ca²⁺: {:.0} nM", systems.homeostasis.calcium_control.free_calcium * 1e9);
    println!("   Energy charge: {:.3}", systems.homeostasis.energy_control.energy_charge);
    println!("   ATP/ADP ratio: {:.1}", systems.homeostasis.energy_control.atp_adp_ratio);
    println!();
    
    println!("⚡ Energy Distribution:");
    println!("   System efficiency: {:.1}%", systems.energy_distribution.system_efficiency * 100.0);
    println!("   Total ATP pool: {:.2e} mol", systems.energy_distribution.total_atp_pool);
    println!();
    
    // Get complete Nebuchadnezzar parameters
    let neb_params = circuit_interface.get_nebuchadnezzar_parameters();
    let systems_params = systems.get_systems_circuit_parameters();
    
    println!("🔗 Complete Nebuchadnezzar Integration:");
    println!("   Hierarchy levels: {}", neb_params.hierarchy_levels);
    println!("   Total circuit nodes: {}", neb_params.circuit_nodes.len());
    println!("   Total circuit edges: {}", neb_params.circuit_edges.len());
    println!("   Global impedance: {:.2e} Ω", neb_params.global_impedance);
    println!("   Stability measure: {:.3}", neb_params.stability_measure);
    println!();
    
    println!("🏗️  Systems Circuit Nodes:");
    println!("   Metabolic impedance: {:.2e} Ω", systems_params.metabolic_node.metabolic_impedance);
    println!("   Thermal resistance: {:.2e} K/W", systems_params.thermal_node.thermal_resistance);
    println!("   Ionic conductance: {:.2e} S", systems_params.ionic_node.ionic_conductance);
    println!("   Hydraulic resistance: {:.2e} Pa·s/m³", systems_params.fluid_node.hydraulic_resistance);
    println!();
    
    println!("🎯 Integration Summary:");
    println!("   ✅ Molecular level: {} proteins with ATP constraints", 
             circuit_interface.membrane.proteins.len());
    println!("   ✅ Mesoscale level: Lipid domains and protein clusters");
    println!("   ✅ Cellular level: Whole membrane patch dynamics");
    println!("   ✅ Systems level: O₂ uptake, thermal regulation, fluid flow");
    println!("   ✅ Circuit interface: Ready for Nebuchadnezzar integration");
    println!();
    
    println!("🌟 The complete 4-level hierarchical system is now operational!");
    println!("Your Nebuchadnezzar system can interface with all levels from molecular");
    println!("proteins to systems-level physiological flows using ATP-based differential");
    println!("equations (dx/dATP) that capture the true energy constraints of biology.");
    
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_systems_integration() {
        let result = main();
        assert!(result.is_ok(), "Systems integration example should run successfully");
    }
} 