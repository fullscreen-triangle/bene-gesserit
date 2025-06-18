<p align="center">
  <h1 align="center">Bene Gesserit</h1>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/License-MIT-blue.svg" alt="License">
  <img src="https://img.shields.io/badge/Rust-1.70+-orange.svg?logo=rust" alt="Rust Version">
  <img src="https://img.shields.io/badge/Python-3.9+-blue.svg?logo=python" alt="Python Version">
  <img src="https://img.shields.io/badge/Status-Development-yellow.svg" alt="Development Status">
  <img src="https://img.shields.io/badge/Version-0.1.0-green.svg" alt="Version">
  <img src="https://img.shields.io/badge/Platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg" alt="Platform">
  <img src="https://img.shields.io/badge/Membrane-Biophysics-purple.svg" alt="Membrane Biophysics">
  <img src="https://img.shields.io/badge/ATP-Constrained-red.svg" alt="ATP Constrained">
  <img src="https://img.shields.io/badge/Circuit-Integration-teal.svg" alt="Circuit Integration">
  <img src="https://img.shields.io/badge/Biological-Authentic-brightgreen.svg" alt="Biological Authentic">
</p>

<p align="center">
  <img src="assets/img/bene-gesserit.png" alt="Bene Gesserit Logo" width="400">
</p>

<p align="center">
  <em>"When a man dips his hand in a river, upon withdrawing, he is no longer the same man, and the river is no longer the same river"</em>
</p>

**"Membranes define the circuit topology; ATP consumption drives the dynamics"**

The **Membrane Dynamics Module** provides biologically authentic cellular membrane simulation that integrates with external systems. This module translates real membrane biophysics into dynamic circuit parameters for ATP-based differential equations, enabling biologically realistic artificial intelligence through authentic cellular constraints.

## What This Module Does

The Membrane Dynamics Module simulates authentic cellular membrane behavior and translates it into computational parameters:

- **Real Membrane Physics**: Lipid bilayers, ion channels, ATP pumps with authentic biophysical properties
- **ATP-Based Dynamics**: Uses `dx/dATP` equations instead of traditional `dx/dt` for energy-constrained modeling
- **Circuit Parameter Generation**: Converts membrane states into dynamic circuit parameters for external systems
- **Biological Constraints**: Imposes authentic metabolic limitations on computational processes

## Architecture Overview

```
Membrane Dynamics Module
├── Molecular Layer (Rust Core)
│   ├── Lipid Bilayer Physics
│   ├── Protein-Membrane Interactions  
│   ├── Membrane Curvature Dynamics
│   └── Electrochemical Gradients
├── Mesoscale Layer (Rust/Python)
│   ├── Lipid Raft Formation
│   ├── Protein Clustering
│   ├── Membrane Domain Organization
│   └── Local Membrane Properties
├── Cellular Layer (Python Extensions)
│   ├── Organelle Membrane Networks
│   ├── Membrane Contact Sites
│   ├── Whole-Cell Membrane Topology
│   └── Membrane Remodeling
└── Circuit Interface Layer
    ├── Membrane → Circuit Parameter Mapping
    ├── Dynamic Circuit Topology Updates
    ├── ATP-Based Differential Equations
    └── External System Integration
```

## Key Features

### 🔬 Authentic Membrane Biophysics
- **Lipid Bilayer Modeling**: Real phospholipid behavior with temperature-dependent fluidity
- **Ion Channel Dynamics**: Voltage-gated and ligand-gated channels with realistic kinetics
- **ATP Pump Mechanics**: Na⁺/K⁺-ATPase, Ca²⁺-ATPase with authentic energy consumption
- **Membrane Curvature**: Protein-induced membrane bending and tubulation

### ⚡ ATP-Based Computation
- **Energy Constraints**: All processes limited by ATP availability
- **Metabolic Realism**: ATP consumption rates match biological values
- **Energy Allocation**: Dynamic ATP budgeting between membrane processes
- **Efficiency Optimization**: Minimizes ATP waste while maintaining function

### 🔌 Circuit Integration
- **Dynamic Parameter Mapping**: Membrane properties → circuit parameters in real-time
- **Topology Updates**: Membrane changes alter circuit connectivity
- **Hierarchical Abstraction**: Multi-scale circuit representation
- **Bidirectional Coupling**: Circuit state influences membrane behavior

### 🧮 External System Integration
- **Nebuchadnezzar Circuits**: Seamless integration with hierarchical probabilistic circuits
- **Orchestrator Interface**: Managed operation under external cognitive systems
- **Modular Design**: Clean APIs for integration with other biological AI components

## Quick Start

### Prerequisites
- Rust 1.70+ (for high-performance molecular layer)
- Python 3.9+ (for mesoscale and cellular layers)
- External system connections (optional for standalone use)

### Basic Usage

```python
from membrane_dynamics import MembranePatch, CircuitInterface

# Create a membrane patch
membrane = MembranePatch(
    area=1e-9,  # 1 μm² patch
    temperature=310.15,  # 37°C
    atp_concentration=5e-3  # 5 mM ATP
)

# Add membrane components
membrane.add_lipid_bilayer(
    composition={'POPC': 0.4, 'POPE': 0.3, 'Cholesterol': 0.3}
)

membrane.add_protein(
    protein_type='NaKATPase',
    density=1000,  # proteins/μm²
    atp_consumption_rate=100  # ATP/s per protein
)

# Initialize circuit interface
circuit = CircuitInterface(
    membrane_patch=membrane,
    update_frequency=1000  # Hz
)

# Run simulation
for t in range(1000):  # 1 second at 1ms steps
    # Update membrane state
    membrane_state = membrane.step(dt=0.001)
    
    # Generate circuit parameters
    circuit_params = circuit.update_from_membrane(membrane_state)
    
    print(f"Time: {t}ms, Voltage: {membrane_state.voltage:.2f}mV, "
          f"ATP used: {membrane_state.atp_consumed:.1f}")
```

## Documentation

### Core Membrane Dynamics
- [Architecture Overview](docs/membrane-dynamics/index.md)
- [Molecular Layer Implementation](docs/membrane-dynamics/molecular-layer.md)
- [Circuit Interface](docs/membrane-dynamics/circuit-interface-layer.md)
- [Quickstart Example](docs/membrane-dynamics/quickstart-example.md)

### External Integration
- [Orchestrator Integration](docs/membrane-dynamics/orchestrator-integration.md) - For managed operation
- [Nebuchadnezzar Circuits](docs/membrane-dynamics/circuit-interface-layer.md) - Circuit system integration

## Key Innovations

### 1. Authentic Membrane Biophysics
Unlike abstract neural networks, this module implements real cellular membrane physics with authentic:
- Lipid phase transitions and fluidity changes
- Protein conformational dynamics
- Ion electrochemical gradients
- ATP-dependent processes

### 2. ATP-Constrained Computation
Computation is limited by ATP availability, creating realistic metabolic constraints:
- Processes compete for limited ATP resources
- Energy efficiency becomes a computational objective
- System behavior emerges from bioenergetic limitations

### 3. Membrane-Circuit Translation
Direct translation of membrane biophysics into circuit parameters:
- Membrane capacitance → circuit capacitance
- Ion channel conductance → circuit resistance
- ATP pump activity → dynamic voltage sources
- Membrane topology → circuit connectivity

### 4. Multi-Scale Integration
Seamless integration across biological scales:
- Molecular: Individual protein dynamics
- Mesoscale: Protein clustering and lipid rafts
- Cellular: Whole-cell membrane networks
- Circuit: Abstract electrical representation

## System Requirements

### Computational
- **Memory**: 8GB+ RAM (16GB+ for large membrane patches)
- **CPU**: Multi-core processor (parallel patch processing)
- **Storage**: 5GB+ for molecular dynamics data
- **GPU**: Optional, for accelerated molecular simulations

### Dependencies
- **Rust toolchain**: For high-performance molecular layer
- **Python scientific stack**: NumPy, SciPy, Matplotlib
- **Optional External Systems**:
  - Nebuchadnezzar circuit system
  - Metacognitive orchestrator
  - ATP budget management systems

## Use Cases

### Standalone Membrane Simulation
- Research into membrane biophysics
- Drug-membrane interaction studies
- Membrane protein function analysis
- Lipid raft dynamics investigation

### Biological AI Integration
- Authentic constraints for artificial neural networks
- Metabolic limitations in AI systems
- Biologically realistic circuit modeling
- Energy-efficient computation research

### Educational Applications
- Teaching membrane biophysics
- Demonstrating ATP-dependent processes
- Visualizing membrane dynamics
- Understanding biological constraints

## Contributing

This module focuses on **biological authenticity** above all else. Contributions should:

### Maintain Biological Realism
- Use experimentally validated parameters
- Implement authentic biophysical mechanisms
- Preserve energy conservation principles
- Respect thermodynamic constraints

### Performance Optimization
- Optimize without sacrificing accuracy
- Implement efficient algorithms for large-scale simulations
- Maintain real-time capability for circuit integration
- Balance computational cost with biological detail

### Documentation Standards
- Document biological basis for all implementations
- Provide experimental validation where possible
- Include performance benchmarks
- Maintain clear API documentation

## Philosophy

*"The membrane is the fundamental unit of biological computation"*

This module recognizes that cellular membranes are not just barriers but active computational elements. By implementing authentic membrane biophysics, we create AI systems that operate within genuine biological constraints rather than abstract mathematical frameworks.

Membranes define what's possible - they set the energy costs, the time scales, and the fundamental limits of biological computation. This module brings those constraints into artificial systems, creating more realistic and potentially more capable AI.

## License

MIT License - See [LICENSE](LICENSE) for details.

---

**Membrane Dynamics**: Where cellular biophysics meets computational intelligence.