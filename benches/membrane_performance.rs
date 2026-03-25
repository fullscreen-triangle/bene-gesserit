use criterion::{black_box, criterion_group, criterion_main, Criterion};
use bene_gesserit::*;

fn benchmark_biological_quantum_state_creation(c: &mut Criterion) {
    c.bench_function("create_physiological_state", |b| {
        b.iter(|| {
            black_box(BiologicalQuantumState::new_physiological())
        })
    });
}

fn benchmark_total_energy_calculation(c: &mut Criterion) {
    let state = BiologicalQuantumState::new_physiological();
    
    c.bench_function("total_energy", |b| {
        b.iter(|| {
            black_box(state.total_energy())
        })
    });
}

fn benchmark_equations_of_motion(c: &mut Criterion) {
    let state = BiologicalQuantumState::new_physiological();
    
    c.bench_function("equations_of_motion", |b| {
        b.iter(|| {
            black_box(state.equations_of_motion())
        })
    });
}

fn benchmark_biological_quantum_solver(c: &mut Criterion) {
    let mut solver = create_biological_quantum_solver();
    let initial_state = create_physiological_state();
    let target = biological_quantum_implementations::QuantumComputationTarget::new();
    
    c.bench_function("solve_biological_quantum_computation", |b| {
        b.iter(|| {
            black_box(solver.solve_biological_quantum_computation(
                &initial_state,
                1000.0, // ATP budget
                1.0,    // time horizon
                &target
            ))
        })
    });
}

fn benchmark_glycolysis_quantum_computer(c: &mut Criterion) {
    c.bench_function("create_glycolysis_quantum_computer", |b| {
        b.iter(|| {
            black_box(glycolysis_quantum_computer::create_glycolysis_quantum_computer())
        })
    });
}

fn benchmark_extended_quantum_biology(c: &mut Criterion) {
    let mut solver = create_biological_quantum_solver();
    let initial_state = create_physiological_state();
    
    c.bench_function("solve_with_extensions", |b| {
        b.iter(|| {
            black_box(solver.solve_with_extensions(
                &initial_state,
                1000.0, // ATP budget
                1.0,    // time horizon
            ))
        })
    });
}

criterion_group!(
    benches,
    benchmark_biological_quantum_state_creation,
    benchmark_total_energy_calculation,
    benchmark_equations_of_motion,
    benchmark_biological_quantum_solver,
    benchmark_glycolysis_quantum_computer,
    benchmark_extended_quantum_biology
);

criterion_main!(benches); 