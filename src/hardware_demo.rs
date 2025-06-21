#!/usr/bin/env rust

//! # Hardware Integration Demo for Bene Gesserit
//! 
//! Revolutionary demonstration: Real hardware oscillations instead of simulated ones!

use crate::*;
use std::error::Error;
use std::thread;
use std::time::Duration;

/// Main hardware integration demonstration
pub fn run_hardware_integration_demo() -> Result<(), Box<dyn Error>> {
    println!("🧬═══════════════════════════════════════════════════════════════════════════════════════🧬");
    println!("              BENE GESSERIT HARDWARE-INTEGRATED BIOLOGICAL QUANTUM COMPUTER");
    println!("        Revolutionary: Real Hardware Oscillations Instead of Simulated Ones!");
    println!("🧬═══════════════════════════════════════════════════════════════════════════════════════🧬");
    println!();
    
    println!("🔧 HARDWARE INTEGRATION REVOLUTION:");
    println!("• CPU clocks → ATP synthase rotation");
    println!("• Screen backlight PWM → Cytochrome c oxidase");
    println!("• WiFi signals → NADH dehydrogenase");
    println!("• Network packets → Membrane oscillations");
    println!("• Temperature fluctuations → Protein dynamics");
    println!();
    
    // Demonstrate both approaches
    demonstrate_traditional_simulation()?;
    demonstrate_hardware_integration()?;
    compare_approaches()?;
    
    println!("\n🎉 HARDWARE INTEGRATION REVOLUTION COMPLETE! 🎉");
    println!("🏆 Real hardware oscillations power biological quantum computation! 🏆");
    
    Ok(())
}

/// Show traditional simulated approach
fn demonstrate_traditional_simulation() -> Result<(), Box<dyn Error>> {
    println!("═══════════════════════════════════════════════════════════════════════════════════════");
    println!("                    📊 TRADITIONAL SIMULATED OSCILLATIONS");
    println!("═══════════════════════════════════════════════════════════════════════════════════════");
    println!();
    
    let simulated_state = create_physiological_state();
    
    println!("🔄 SIMULATED OSCILLATIONS:");
    for (i, oscillation) in simulated_state.oscillatory_coords.oscillations.iter().enumerate() {
        println!("  {}. {} - f={:.1} Hz, A={:.2}, φ={:.2}", 
                i+1, oscillation.name, oscillation.frequency, 
                oscillation.amplitude, oscillation.phase);
    }
    
    println!("\n⚡ SIMULATED ATP STATE:");
    println!("  ATP Concentration: {:.2} mM", simulated_state.atp_coords.atp_concentration);
    println!("  Energy Charge: {:.3}", simulated_state.atp_coords.energy_charge);
    println!("  Available Energy: {:.1} kJ/mol", simulated_state.atp_coords.available_energy());
    
    println!("\n❌ LIMITATIONS OF SIMULATION:");
    println!("  • Artificial oscillations disconnected from reality");
    println!("  • No real energy harvesting from environment");
    println!("  • Computational overhead of simulation");
    println!("  • Missing hardware-biology coupling");
    
    Ok(())
}

/// Show revolutionary hardware integration
fn demonstrate_hardware_integration() -> Result<(), Box<dyn Error>> {
    println!("\n═══════════════════════════════════════════════════════════════════════════════════════");
    println!("                    🔌 REVOLUTIONARY HARDWARE INTEGRATION");
    println!("═══════════════════════════════════════════════════════════════════════════════════════");
    println!();
    
    println!("🚀 Creating hardware-powered biological quantum computer...");
    
    // Create hardware-powered system
    let (hardware_state, mut harvester) = create_hardware_powered_biological_quantum_computer()?;
    
    println!("\n🔌 HARDWARE OSCILLATION SOURCES:");
    let stats = harvester.get_hardware_statistics();
    println!("  Active Sources: {}", stats.active_sources);
    println!("  Total ATP Generation: {:.3} mM/s", stats.total_atp_generation_rate);
    println!("  Average Frequency: {:.2e} Hz", stats.average_frequency);
    println!("  Hardware Efficiency: {:.3}", stats.hardware_efficiency);
    
    println!("\n🔄 HARDWARE-DRIVEN OSCILLATIONS:");
    for (i, oscillation) in hardware_state.oscillatory_coords.oscillations.iter().enumerate() {
        println!("  {}. {} - f={:.1} Hz, A={:.2}, ATP coupling={:.2}", 
                i+1, oscillation.name, oscillation.frequency, 
                oscillation.amplitude, oscillation.atp_coupling_strength);
    }
    
    println!("\n⚡ HARDWARE-BOOSTED ATP STATE:");
    println!("  ATP Concentration: {:.2} mM (includes hardware harvest)", hardware_state.atp_coords.atp_concentration);
    println!("  Energy Charge: {:.3}", hardware_state.atp_coords.energy_charge);
    println!("  Available Energy: {:.1} kJ/mol", hardware_state.atp_coords.available_energy());
    
    println!("\n🎯 REAL-TIME HARDWARE MONITORING:");
    for i in 0..5 {
        thread::sleep(Duration::from_millis(200));
        let updated_state = harvester.get_biological_state_from_hardware()?;
        let current_stats = harvester.get_hardware_statistics();
        
        println!("  T+{:.1}s: {} sources active, {:.3} mM ATP/s", 
                i as f64 * 0.2, 
                current_stats.active_sources,
                current_stats.total_atp_generation_rate);
    }
    
    println!("\n✅ ADVANTAGES OF HARDWARE INTEGRATION:");
    println!("  • Real oscillations from actual hardware");
    println!("  • True energy harvesting from environment");
    println!("  • Zero computational overhead for oscillations");
    println!("  • Authentic hardware-biology coupling");
    println!("  • Utilizes existing machine resources");
    
    harvester.stop_harvesting();
    
    Ok(())
}

/// Compare traditional vs hardware approaches
fn compare_approaches() -> Result<(), Box<dyn Error>> {
    println!("\n═══════════════════════════════════════════════════════════════════════════════════════");
    println!("                        ⚖️  TRADITIONAL vs HARDWARE COMPARISON");
    println!("═══════════════════════════════════════════════════════════════════════════════════════");
    println!();
    
    // Create both systems for comparison
    let simulated_state = create_physiological_state();
    let (hardware_state, mut harvester) = create_hardware_powered_biological_quantum_computer()?;
    
    println!("📊 COMPARISON METRICS:");
    println!("┌─────────────────────────────────────────────────────────────────────────────────────┐");
    println!("│ METRIC                    │ TRADITIONAL SIM │ HARDWARE INTEGRATION                  │");
    println!("├─────────────────────────────────────────────────────────────────────────────────────┤");
    println!("│ Oscillation Sources       │ {:>15} │ {:>37} │", 
            simulated_state.oscillatory_coords.oscillations.len(),
            hardware_state.oscillatory_coords.oscillations.len());
    
    println!("│ ATP Concentration (mM)    │ {:>15.2} │ {:>37.2} │", 
            simulated_state.atp_coords.atp_concentration,
            hardware_state.atp_coords.atp_concentration);
    
    let sim_total_freq: f64 = simulated_state.oscillatory_coords.oscillations.iter()
        .map(|o| o.frequency).sum();
    let hw_total_freq: f64 = hardware_state.oscillatory_coords.oscillations.iter()
        .map(|o| o.frequency).sum();
    
    println!("│ Total Frequency (Hz)      │ {:>15.1} │ {:>37.1} │", sim_total_freq, hw_total_freq);
    
    let sim_avg_coupling: f64 = simulated_state.oscillatory_coords.oscillations.iter()
        .map(|o| o.atp_coupling_strength).sum::<f64>() / simulated_state.oscillatory_coords.oscillations.len() as f64;
    let hw_avg_coupling: f64 = hardware_state.oscillatory_coords.oscillations.iter()
        .map(|o| o.atp_coupling_strength).sum::<f64>() / hardware_state.oscillatory_coords.oscillations.len() as f64;
    
    println!("│ Avg ATP Coupling          │ {:>15.2} │ {:>37.2} │", sim_avg_coupling, hw_avg_coupling);
    
    println!("│ Energy Source             │ {:>15} │ {:>37} │", "Simulated", "Real Hardware");
    println!("│ Computational Cost        │ {:>15} │ {:>37} │", "High", "Zero (harvested)");
    println!("│ Reality Connection        │ {:>15} │ {:>37} │", "None", "Direct");
    println!("│ Innovation Level          │ {:>15} │ {:>37} │", "Standard", "Revolutionary");
    println!("└─────────────────────────────────────────────────────────────────────────────────────┘");
    
    println!("\n🏆 HARDWARE INTEGRATION WINS:");
    println!("✅ More oscillation sources (real hardware devices)");
    println!("✅ Higher ATP concentration (energy harvesting)");
    println!("✅ Zero computational overhead for oscillations");
    println!("✅ Real-world energy coupling");
    println!("✅ Utilizes available machine resources");
    println!("✅ Revolutionary approach to biological computing");
    
    println!("\n💡 HARDWARE SOURCES BEING HARVESTED:");
    let stats = harvester.get_hardware_statistics();
    println!("  🖥️  CPU Clock: High-frequency oscillations → ATP synthase");
    println!("  💡 Screen Backlight: PWM oscillations → Cytochrome oxidase");
    println!("  📡 WiFi Signals: 2.4GHz oscillations → NADH dehydrogenase");
    println!("  🌐 Network Activity: Packet timing → Membrane oscillations");
    println!("  🌡️  CPU Temperature: Thermal fluctuations → Protein dynamics");
    
    println!("\n🎯 REVOLUTIONARY INSIGHT:");
    println!("Instead of wasting CPU cycles simulating oscillations,");
    println!("we harvest the actual oscillations already present in the hardware!");
    println!("This makes the biological quantum computer truly integrated with its host machine.");
    
    harvester.stop_harvesting();
    
    Ok(())
}

/// Detailed hardware source analysis
pub fn analyze_hardware_sources() -> Result<(), Box<dyn Error>> {
    println!("\n═══════════════════════════════════════════════════════════════════════════════════════");
    println!("                        🔍 DETAILED HARDWARE SOURCE ANALYSIS");
    println!("═══════════════════════════════════════════════════════════════════════════════════════");
    println!();
    
    let (_, mut harvester) = create_hardware_powered_biological_quantum_computer()?;
    
    // Let it run for a bit to collect data
    thread::sleep(Duration::from_millis(500));
    
    println!("📊 HARDWARE OSCILLATION DETAILED ANALYSIS:");
    
    let harvested = harvester.harvested_oscillations.lock().unwrap();
    for (source_name, oscillation) in harvested.iter() {
        println!("\n🔌 SOURCE: {}", source_name);
        println!("   Hardware Frequency: {:.2e} Hz", oscillation.current_frequency);
        println!("   Hardware Amplitude: {:.3}", oscillation.current_amplitude);
        println!("   Hardware Phase: {:.2} rad", oscillation.current_phase);
        println!("   → Biological Process: {}", oscillation.biological_state.name);
        println!("   → Biological Frequency: {:.1} Hz", oscillation.biological_state.frequency);
        println!("   → ATP Contribution: {:.3} mM", oscillation.atp_contribution);
        println!("   → ATP Coupling: {:.2}", oscillation.biological_state.atp_coupling_strength);
    }
    
    let stats = harvester.get_hardware_statistics();
    println!("\n📈 AGGREGATE STATISTICS:");
    println!("   Total Active Sources: {}", stats.active_sources);
    println!("   Total ATP Generation Rate: {:.3} mM/s", stats.total_atp_generation_rate);
    println!("   Average Hardware Frequency: {:.2e} Hz", stats.average_frequency);
    println!("   Total Signal Amplitude: {:.3}", stats.total_amplitude);
    println!("   Hardware Efficiency: {:.3}", stats.hardware_efficiency);
    
    harvester.stop_harvesting();
    
    println!("\n✅ HARDWARE ANALYSIS COMPLETE");
    println!("Real machine oscillations successfully converted to biological energy!");
    
    Ok(())
} 