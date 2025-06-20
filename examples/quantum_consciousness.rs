use bene_gesserit::*;
use bene_gesserit::quantum_extensions::*;
use ndarray::Array2;
use num_complex::Complex;
use std::collections::HashMap;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🧠 QUANTUM CONSCIOUSNESS EMERGENCE EXAMPLE");
    println!("Demonstrating consciousness as integrated oscillatory entropy");
    println!();

    // Create biological quantum state
    let state = create_consciousness_state();
    
    // Calculate extended quantum measures
    let measures = ExtendedQuantumMeasures::calculate_all(&state);
    
    println!("📊 QUANTUM BIOLOGY MEASURES:");
    println!("  Consciousness Φ: {:.6}", measures.consciousness_phi);
    println!("  Evolutionary Fitness: {:.6}", measures.evolutionary_fitness);
    println!("  Temporal Correlations: {:.6}", measures.temporal_correlations);
    println!("  Cosmic Entanglement: {:.6e}", measures.cosmic_entanglement);
    println!("  Total Extended Energy: {:.6}", measures.total_extended_energy());
    println!();

    // Demonstrate consciousness enhancement
    let enhanced_state = enhance_consciousness(&state);
    let enhanced_measures = ExtendedQuantumMeasures::calculate_all(&enhanced_state);
    
    println!("🌟 ENHANCED CONSCIOUSNESS STATE:");
    println!("  Enhanced Consciousness Φ: {:.6}", enhanced_measures.consciousness_phi);
    println!("  Enhancement Factor: {:.2}x", enhanced_measures.consciousness_phi / measures.consciousness_phi);
    
    if enhanced_measures.consciousness_phi > 0.1 {
        println!("  ✨ CONSCIOUSNESS EMERGENCE ACHIEVED!");
    }

    Ok(())
}

fn create_consciousness_state() -> BiologicalQuantumState {
    let atp_coords = AtpCoordinates {
        atp_concentration: 5.0,
        adp_concentration: 1.0,
        pi_concentration: 5.0,
        energy_charge: 0.85,
        atp_oscillation_amplitude: 0.5,
        atp_oscillation_phase: 0.0,
        atp_oscillation_frequency: 10.0,
    };

    let oscillations = vec![
        OscillationState::new("Neural_Oscillation", 1.0, 0.0, 40.0),  // Gamma waves
        OscillationState::new("Metabolic_Rhythm", 0.8, std::f64::consts::PI/3.0, 1.0),
        OscillationState::new("Quantum_Coherence", 1.2, std::f64::consts::PI/2.0, 100.0),
    ];

    let oscillatory_coords = OscillatoryCoordinates {
        oscillations,
        oscillatory_momenta: vec![0.1, -0.05, 0.15],
        phase_coupling_matrix: Array2::eye(3) * 0.2,
        membrane_oscillations: vec![],
    };

    let quantum_states = vec![
        QuantumStateAmplitude::new("Microtubule_A", Complex::new(0.8, 0.2)),
        QuantumStateAmplitude::new("Microtubule_B", Complex::new(0.6, 0.4)),
        QuantumStateAmplitude::new("Neural_Quantum", Complex::new(0.7, 0.3)),
    ];

    let environmental_coupling = EnvironmentalCoupling {
        coupling_strength: 0.1,
        correlation_time: 1e-12,
        temperature: 310.0,
        enhancement_factor: 2.0,
    };

    let membrane_coords = MembraneQuantumCoordinates {
        quantum_states,
        environmental_coupling,
        tunneling_states: vec![],
        membrane_properties: MembraneProperties {
            thickness: 4e-9,
            dielectric_constant: 2.1,
            protein_density: 1000.0,
            lipid_composition: LipidComposition {
                phospholipid_fraction: 0.6,
                cholesterol_fraction: 0.3,
                other_lipids_fraction: 0.1,
            },
        },
    };

    let mut endpoint_distributions = HashMap::new();
    endpoint_distributions.insert("Neural_Oscillation".to_string(), EndpointDistribution {
        positions: vec![0.0, 0.5, 1.0],
        probabilities: vec![0.3, 0.4, 0.3],
        velocities: vec![0.0, 0.1, 0.0],
        energies: vec![1.0, 2.0, 1.0],
    });

    let entropy_coords = OscillatoryEntropyCoordinates {
        endpoint_distributions,
        current_entropy: 0.5,
        entropy_production_rate: 0.05,
        membrane_endpoint_entropy: 0.3,
        quantum_tunneling_entropy: 0.1,
    };

    BiologicalQuantumState {
        atp_coords,
        oscillatory_coords,
        membrane_coords,
        entropy_coords,
    }
}

fn enhance_consciousness(state: &BiologicalQuantumState) -> BiologicalQuantumState {
    let mut enhanced_state = state.clone();
    
    // Enhance oscillation amplitudes for greater integration
    for osc in &mut enhanced_state.oscillatory_coords.oscillations {
        osc.amplitude *= 1.5;
        osc.atp_coupling_strength *= 1.2;
    }
    
    // Enhance quantum coherence
    for qs in &mut enhanced_state.membrane_coords.quantum_states {
        qs.amplitude *= 1.3;
    }
    
    // Increase ATP availability
    enhanced_state.atp_coords.energy_charge *= 1.1;
    
    enhanced_state
} 