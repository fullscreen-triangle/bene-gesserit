<p align="center">
  <h1 align="center">Bene Gesserit</h1>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/License-MIT%20%7C%20Apache--2.0-blue.svg" alt="License">
  <img src="https://img.shields.io/badge/Rust-1.70+-orange.svg?logo=rust" alt="Rust Version">
  <img src="https://img.shields.io/badge/Status-Development-yellow.svg" alt="Development Status">
  <img src="https://img.shields.io/badge/Platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg" alt="Platform">
</p>

<p align="center">
  <img src="assets/img/bene-gesserit.png" alt="Bene Gesserit Logo" width="400">
</p>

**Bene Gesserit** is a categorical processing architecture in which amphiphilic bilayer membranes serve as the computational substrate. The framework derives from three axioms — bounded phase space, no null state, and finite observational resolution — and proves that oscillatory dynamics, categorical partition structure, and biological membrane processes are the same mathematical object viewed through three equivalent projections.

## Theoretical Foundation

The architecture rests on the **Bounded Phase Space (BPS)** framework. From a single geometric axiom — *all physical systems occupy bounded regions of phase space admitting partition and nesting* — the following are derived with zero free parameters:

### Triple Equivalence Theorem

Oscillatory entropy, categorical entropy, and partition entropy are algebraically identical:

$$S_{\text{osc}} = S_{\text{cat}} = S_{\text{part}} = k_B \mathcal{M} \ln n$$

where $\mathcal{M}$ is partition depth and $n$ is the state count at that depth. This is not analogy — it is mathematical identity proved by three independent derivations.

### Triple Isomorphism

Three categories — **O** (oscillatory), **C** (categorical), **B** (biological) — are connected by explicit functors $F_{OC}$, $F_{CB}$, $F_{BO}$ forming a triangle of equivalences:

$$F_{BO} \circ F_{CB} \circ F_{OC} \cong \text{Id}_O$$

Every component in the framework is **one object with three equivalent faces**:

| Oscillatory | Categorical | Biological |
|---|---|---|
| Coupled oscillator state | Trajectory-terminus-memory triple | Membrane conformational trajectory |
| Phase-lock synchronization | Path-independent convergence | Synaptic junction coupling |
| Frequency bandpass | Categorical aperture ($W = 0$) | Ion channel selectivity filter |
| First resonance completion | Shortest categorical path wins | Signal priority by membrane proximity |
| R-C-L impedance | $S_k / S_t / S_e$ minimization | Knowledge / temporal / evolution processing |
| Decay envelope overlap | Integration window | Cortical processing window |
| Anharmonic drift | Poincaré deviation ($\delta > 0$ a.s.) | Non-recurrent conformational cycling |
| Resonant mode selector | BMD filtering operator ($\eta \sim 10^{9}$) | Pattern-recognition catalysis |
| Carrier density transition | Depletion region ($V_{bi} = 615$ mV) | Lipid raft boundary P-N junction |
| Phase accumulation | $\int (dH/dt)^+ dt$ | Partition depth integration |

## Architecture

The processing paradigm is **trajectory completion** — not forward simulation. Given a target terminus in S-entropy space, the system navigates backward to the penultimate state and completes:

$$O(x) = C(x) = P(x) = \text{Resolve}(\mathcal{A}_x)$$

Observation, computation, and processing are identical operations: categorical address resolution.

### Categorical Speedup

Binary processors traverse partition space at fixed depth $\Delta\mathcal{M} = 1$ per operation, wasting ~31% of partition capacity. The membrane architecture operates at variable effective base $b_{\text{eff}}(T)$ and achieves:

- **Binary**: $O(N)$ sequential steps
- **Categorical**: $O(\log_3 N)$ ternary navigations

For $N = 10^6$: speedup $\approx 79{,}000\times$. For $N = 10^{12}$: speedup $\approx 4 \times 10^{10}\times$.

### Processing Stack

```
Layer 7: Interface
         Problem encoding into S-entropy coordinates, resonant coupling I/O,
         gear-reduced clock hierarchy

Layer 6: Runtime
         Penultimate state scheduler (priority = 1/d_cat), trajectory completion,
         triple equivalence runtime verification: dM/dt = ω/(2π/M) = 1/⟨τ_p⟩

Layer 5: Memory
         Categorical memory with S-entropy addressing, 3^k ternary hierarchy,
         precision-by-difference navigation, tier placement by categorical distance

Layer 4: Hybrid Microfluidic Circuit Dynamics
         Trajectory-terminus-memory triples, decay envelope intersection,
         R-C-L mode selection by S-entropy minimization, Poincaré deviation,
         five operational regimes (coherent/turbulent/cascade/aperture/phase-locked)

Layer 3: Circuits
         Tri-dimensional logic gates (AND/OR/XOR on S-coordinates simultaneously),
         ALU via frequency superposition (add) and modulation (multiply),
         virtual foundry with femtosecond processor lifecycle

Layer 2: Membrane Substrate
         Lipid bilayer (d = 4.0 nm, A_L = 0.64 nm², κ = 20 k_BT),
         biological semiconductor P-N junctions (V_bi = 615 mV),
         BMD transistor (on/off = 42.1, switching < 1 μs, 10⁶× efficiency over silicon),
         categorical apertures (zero-work filtering, W = 0),
         Kuramoto phase-lock networks

Layer 1: S-Entropy Space
         Coordinates (S_k, S_t, S_e) ∈ [0,1]³, ternary trit addressing,
         five operational regimes, neural partition Lagrangian L_NPL,
         Euler-Lagrange dynamics, Noether conservation laws

Layer 0: Axioms
         Bounded phase space, no null state, finite observational resolution,
         partition coordinates (n, l, m, s), capacity C(n) = 2n²,
         partition depth M = Σ log_b(k_i), triple equivalence
```

### Key Properties

- **No von Neumann bottleneck**: the membrane is simultaneously processor, memory, and interconnect
- **Zero-work filtering**: categorical apertures select by geometry, not energy — Landauer's erasure principle does not apply
- **Sensing = computing**: the Operation Equivalence theorem proves that coupling to an external signal IS computation
- **Five operational regimes**: coherent ($R > 0.8$), phase-locked ($R > 0.95$), cascade ($0.3 < R < 0.8$), aperture-dominated, and turbulent ($R < 0.3$), classified by the Kuramoto order parameter
- **Computational completeness**: tri-dimensional gates provide Boolean completeness, ALU provides arithmetic completeness, trajectory completion provides algorithmic completeness

## Quantitative Predictions

All derived from the three axioms with zero free parameters fitted to membrane data:

| Prediction | Derived Value | Experimental Value |
|---|---|---|
| Hydrogen bond energy | 27 kJ/mol | 20–25 kJ/mol |
| Bilayer thickness | 4.0 nm | 3.7–4.3 nm |
| Area per lipid | 0.64 nm² | 0.60–0.72 nm² |
| Bending modulus | 20 $k_BT$ | 15–25 $k_BT$ |
| Membrane potential | −70 mV | −70 mV |
| P-N junction $V_{bi}$ | 615 mV | — |
| Proton transfer frequency | 4.06 × 10¹³ Hz | 4.00 × 10¹³ Hz |
| BMD energy per operation | 6 × 10⁻²⁰ J | — |
| Processing throughput | ~10²⁰ ops/s/cm² | — |

## Publications

- **Triple Isomorphism Architecture**: Formal equivalence of oscillatory, categorical, and biological descriptions in bounded phase space systems — [`reverend/publications/isomorphism-atlas/`](reverend/publications/isomorphism-atlas/)

- **Membrane-Mediated Categorical Processing Architecture**: Trajectory completion, aperture-based filtering, and partition gradient optimization on amphiphilic bilayer substrates — [`reverend/publications/membrane-cognitive-architecture/`](reverend/publications/membrane-cognitive-architecture/)

## Project Structure

```
bene-gesserit/
├── src/                        # Rust implementation
├── reverend/
│   ├── publications/           # Formal papers (LaTeX)
│   │   ├── isomorphism-atlas/
│   │   └── membrane-cognitive-architecture/
│   └── sources/                # Source derivations
│       ├── lipid-derivation.tex
│       └── categorical-converter.tex
├── examples/                   # Usage examples
├── benches/                    # Performance benchmarks
├── docs/                       # Documentation
├── Cargo.toml
└── README.md
```

## Building

```bash
# Prerequisites: Rust 1.70+
cargo build --release
cargo test
cargo bench
```

## License

Dual-licensed under [MIT](LICENSE-MIT) and [Apache 2.0](LICENSE-APACHE).
