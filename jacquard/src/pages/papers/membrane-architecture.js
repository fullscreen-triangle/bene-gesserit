import Head from "next/head";
import Layout from "@/components/Layout";
import TransitionEffect from "@/components/TransitionEffect";
import dynamic from "next/dynamic";

const PharmaKuramotoChart = dynamic(() => import("@/components/charts/PharmaKuramotoChart"), { ssr: false });
const ComputingStackChart = dynamic(() => import("@/components/charts/ComputingStackChart"), { ssr: false });
const SpeedupChart = dynamic(() => import("@/components/charts/SpeedupChart"), { ssr: false });

const MembraneArchitecture = () => {
  return (
    <>
      <Head>
        <title>Membrane-Mediated Categorical Processing Architecture &mdash; Bene Gesserit</title>
        <meta name="description" content="On the Geometric Consequences of Aperture Traversal Mechanisms on Amphiphilic Lipid Bilayer Substrates" />
      </Head>
      <TransitionEffect />
      <main className="w-full mb-16 flex flex-col items-center justify-center dark:text-light">
        <Layout className="pt-16">
          <article className="max-w-4xl mx-auto">
            <h1 className="text-4xl font-bold mb-2 leading-tight">
              On the Geometric Consequences of Aperture Traversal Mechanisms on Amphiphilic Lipid Bilayer Substrates: Membrane-Mediated Categorical Processing Architecture
            </h1>
            <p className="text-gray-400 mb-8">Kundai Farai Sachikonye &bull; AIMe Registry for Artificial Intelligence</p>

            {/* Abstract */}
            <div className="bg-gray-800/50 border border-gray-700 rounded-lg p-6 mb-12">
              <h2 className="text-lg font-semibold mb-3 text-gray-300">Abstract</h2>
              <p className="text-gray-300 leading-relaxed">
                We derive a complete categorical processing architecture in which amphiphilic bilayer membranes serve as the computational substrate, with hybrid microfluidic circuit dynamics replacing binary logic. The architecture rests on three axioms&mdash;bounded phase space, no null state, and finite observational resolution&mdash;from which we prove the Triple Equivalence Theorem <span className="font-mono text-cyan-400">S_osc = S_cat = S_part = k_B M ln n</span>, establishing that oscillatory, categorical, and partition descriptions of membrane dynamics are faces of a single mathematical object. We formalize categorical processing as trajectory completion in S-entropy space <span className="font-mono text-cyan-400">(S_k, S_t, S_e) &isin; [0,1]&sup3;</span>, prove the Operation Equivalence <span className="font-mono text-cyan-400">O(x) = C(x) = P(x)</span>&mdash;observation, computation, and processing are identical categorical address resolutions&mdash;and derive a categorical speedup: <span className="font-mono text-cyan-400">O(log_3 N)</span> navigations for problems requiring O(N) sequential binary steps. The architecture comprises: (i) the amphiphilic bilayer as oscillator network with each lipid executing ~10&sup;11; conformational operations per second; (ii) biological semiconductor physics at lipid raft boundaries yielding a BMD transistor with on/off ratio 42.1; (iii) categorical apertures performing zero-work filtering; (iv) backward trajectory completion for O(log_3 N) categorical navigation; (v) tri-dimensional logic on S-entropy coordinates; (vi) an R-C-L processing trichotomy; (vii) five operational regimes classified by the Kuramoto order parameter; and (viii) partition gradient optimization. The architecture eliminates the von Neumann bottleneck&mdash;the membrane is simultaneously processor, memory, and interconnect. Throughput: ~10&sup;20; operations per second per cm&sup2;, energy efficiency 10&sup6;&times; silicon per bit.
              </p>
            </div>

            {/* Introduction */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">1. Introduction</h2>
              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">The Binary Processing Bottleneck</h3>
              <p className="text-gray-300 leading-relaxed mb-4">All conventional processors are base-2 partition traversal systems. Each elementary operation advances partition depth by exactly &Delta;M = 1 at base b = 2, producing entropy &Delta;S = k_B ln 2 per bit transition. The von Neumann architecture separates memory from processing, requiring data to traverse a shared bus&mdash;the von Neumann bottleneck. Binary quantization imposes a systematic inefficiency: encoding N states in base b = 2 requires &lceil;log_2 N&rceil; bits, with an average waste fraction of approximately 31%.</p>
              <div className="bg-gray-900 border border-gray-700 rounded p-4 my-6 text-center font-mono text-gray-200">
                <p>&langle;W&rangle; = 1 &minus; 1/ln 2 &asymp; 0.31</p>
              </div>

              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">Biological Membrane Systems as Categorical Processors</h3>
              <p className="text-gray-300 leading-relaxed mb-4">The amphiphilic bilayer is not a passive boundary but an active computational substrate. Each lipid executes ~10&sup;11; conformational oscillations per second. The effective computational base is not fixed at b = 2 but varies continuously with temperature:</p>
              <div className="bg-gray-900 border border-gray-700 rounded p-4 my-6 text-center font-mono text-gray-200">b_eff(T) = exp(S_config / (k_B M))</div>
              <p className="text-gray-300 leading-relaxed mb-4">At physiological temperature, b_eff &asymp; 3, yielding ternary partition structure. The key identity underlying the entire architecture is:</p>
              <div className="bg-gray-900 border border-gray-700 rounded p-4 my-6 text-center font-mono text-gray-200 text-lg">O(x) = C(x) = P(x) = Resolve(A_x)</div>
              <p className="text-gray-300 leading-relaxed mb-4">Observation, computation, and processing are identical operations: categorical address resolution.</p>
            </section>

            {/* Foundations */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">2. Foundations</h2>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Axiom 1 (Bounded Phase Space)</p>
                <p className="text-gray-300">All physical systems occupy bounded regions of phase space admitting partition and nesting.</p>
              </div>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Axiom 2 (No Null State)</p>
                <p className="text-gray-300">No physical system occupies the null partition state M = 0. Every system has M &ge; M_min &gt; 0.</p>
              </div>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Axiom 3 (Finite Observational Resolution)</p>
                <p className="text-gray-300">Every observation resolves partition structure to finite depth M_obs &lt; &infin;.</p>
              </div>

              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Triple Equivalence)</p>
                <p className="font-mono text-center text-gray-200 text-lg">S_osc = S_cat = S_part = k_B M ln n</p>
                <p className="text-gray-300 mt-2">Every component of the architecture has three faces&mdash;oscillatory, categorical, and biological&mdash;that are the same mathematical object. A membrane oscillation IS a categorical morphism IS a partition refinement.</p>
              </div>

              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Variance&ndash;Free Energy Identity)</p>
                <p className="font-mono text-center text-gray-200">F = k_BT &middot; &sigma;&sup2;(&phi;)</p>
                <p className="text-gray-300 mt-2">The Helmholtz free energy equals k_BT times the phase variance. Hybrid microfluidic circuit dynamics IS variance minimization.</p>
              </div>
            </section>

            {/* Membrane Substrate */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">3. The Membrane as Computational Substrate</h2>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Membrane Necessity)</p>
                <p className="text-gray-300">Any bounded aqueous system at finite temperature generates amphiphilic bilayer boundaries as minimum-cost partition boundaries, with d = 4.0 nm, A_L = 0.64 nm&sup2;, &kappa; = 20 k_BT. All three quantities are derived from partition geometry with zero free parameters.</p>
              </div>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Membrane Oscillator Network)</p>
                <p className="text-gray-300 mb-2">A lipid membrane of area A constitutes a coupled oscillator network with total throughput:</p>
                <p className="font-mono text-center text-gray-200">&Phi;(A) = (A/A_L) &middot; &nu;_conf &asymp; 1.56 &times; 10&sup;25; A [m&sup;&minus;2; s&sup;&minus;1;]</p>
                <p className="text-gray-300 mt-2">For A = 1 cm&sup2;: &Phi; &asymp; 10&sup;20; ops/s (conservative, error-corrected estimate).</p>
              </div>

              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">Biological Semiconductor Physics</h3>
              <p className="text-gray-300 leading-relaxed mb-4">In the membrane oscillator network, carriers are classified as: <strong>P-type carriers</strong> (oscillatory holes&mdash;vacancies where local partition depth M is reduced) and <strong>N-type carriers</strong> (molecular oscillators that increase local partition depth). Lipid raft boundaries form P-N junctions with built-in potential V_bi = 615 mV and rectification ratio 42.1, yielding the BMD transistor.</p>
            </section>

            {/* Categorical Apertures */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">4. Categorical Aperture Architecture</h2>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Zero-Work Aperture Filtering)</p>
                <p className="text-gray-300">The work performed by a categorical aperture on any trajectory passing through it is identically zero: W_A = 0. Landauer&apos;s erasure principle does not apply because the aperture performs no erasure&mdash;it is a fixed geometric constraint, not a measurement apparatus.</p>
              </div>

              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">Aperture Taxonomy by Multipole Order</h3>
              <div className="overflow-x-auto my-6">
                <table className="w-full text-left border-collapse">
                  <thead><tr className="border-b border-gray-700"><th className="p-3 text-gray-400">Type</th><th className="p-3 text-gray-400">l</th><th className="p-3 text-gray-400">Description</th><th className="p-3 text-gray-400">Hill coefficient</th></tr></thead>
                  <tbody>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Monopole</td><td className="p-3 text-gray-300">0</td><td className="p-3 text-gray-300">Radially symmetric, couples to single target</td><td className="p-3 text-gray-300">n_H = 1</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Dipole</td><td className="p-3 text-gray-300">1</td><td className="p-3 text-gray-300">Two-lobed, couples to two targets</td><td className="p-3 text-gray-300">n_H = 2</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Quadrupole</td><td className="p-3 text-gray-300">2</td><td className="p-3 text-gray-300">Four-lobed, tetrameric structure</td><td className="p-3 text-gray-300">n_H = 4</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">General multipole</td><td className="p-3 text-gray-300">l</td><td className="p-3 text-gray-300">(2l+1)-lobed, Y_l^m angular dependence</td><td className="p-3 text-gray-300">n_H = 2^l</td></tr>
                  </tbody>
                </table>
              </div>

              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Ternary Trisection)</p>
                <p className="text-gray-300 mb-2">Each categorical aperture performs a ternary trisection of S-entropy space, refining along one of three dimensions (S_k, S_t, or S_e). A cascade of k apertures navigates to a volume element of measure 3&sup;&minus;k;.</p>
                <p className="font-mono text-center text-gray-200 mt-2">V_k = 3&sup;&minus;k;</p>
              </div>

              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Aperture Composition)</p>
                <p className="text-gray-300">Multiple apertures compose with multiplicative selectivity: &eta;_total = &Pi; &eta;_i. For n = 10 stages each with &eta;_i = 10&sup;10;: &eta;_total = 10&sup;100;. The total work remains zero: W_total = 0.</p>
              </div>
            </section>

            {/* Trajectory Completion */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">5. Trajectory Completion Architecture</h2>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Operation Equivalence)</p>
                <p className="font-mono text-center text-gray-200">O(x) = C(x) = P(x) = Resolve(A_x)</p>
                <p className="text-gray-300 mt-2">Observation, computation, and processing are identical operations: categorical address resolution. A membrane that couples to an external signal has already performed categorical processing on it. The coupling process IS the computation. This eliminates the sensor-processor separation that is a bottleneck in conventional architectures.</p>
              </div>

              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">Backward Trajectory Completion</h3>
              <p className="text-gray-300 leading-relaxed mb-4">Given a desired terminus &Gamma;_f, backward trajectory completion constructs a trajectory by: (1) specifying &Gamma;_f as an S-entropy coordinate, (2) identifying the penultimate state P_pen within one partition step, (3) applying the completion morphism. The penultimate state always exists because the ternary tree has branching factor 3, with at least 2 siblings per parent cell.</p>

              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Categorical Speedup)</p>
                <p className="text-gray-300 mb-2">Navigating to a target coordinate in S-entropy space with precision &epsilon; requires:</p>
                <p className="font-mono text-center text-gray-200">k = O(log_3(1/&epsilon;))</p>
                <p className="text-gray-300 mt-2">trit operations. For a task requiring N binary steps: O(log_3 N) categorical navigations. Speedup factor:</p>
                <p className="font-mono text-center text-gray-200 mt-1">S(N) = N / log_3 N</p>
                <div className="mt-2 text-gray-300">
                  <p>N = 10&sup6;: speedup &asymp; 79,500</p>
                  <p>N = 10&sup;12;: speedup &asymp; 3.98 &times; 10&sup;10;</p>
                </div>
                <p className="text-gray-300 mt-2">This is not a constant-factor improvement but a complexity class change from O(N) to O(log N).</p>
              </div>

              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Corollary (P = NP for Categorical Problems)</p>
                <p className="text-gray-300">For any decision problem whose instances are representable as categorical partition structures in bounded phase space, the complexity of finding a solution equals the complexity of verifying a solution. The operation compute(p_0, C) = verify(p_k, C) because both are identical categorical address resolutions.</p>
              </div>

              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (&epsilon;-Boundary Existence)</p>
                <p className="text-gray-300">The exact final state at d_cat = 0 is unreachable by any finite sequence of categorical operations. The physical solution exists at the &epsilon;-boundary: 0 &lt; d_cat &le; &epsilon;. This is the partition-theoretic analogue of the incompleteness phenomenon: the system can approach any target with arbitrary precision but cannot exactly reach it in finite steps. The &epsilon;-boundary itself is the physical answer.</p>
              </div>
            </section>

            {/* R-C-L Trichotomy */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">6. The R-C-L Processing Trichotomy</h2>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (R-C-L Trichotomy)</p>
                <p className="text-gray-300">Every processing element simultaneously functions as:</p>
                <ul className="list-disc list-inside text-gray-300 mt-2 space-y-2">
                  <li><strong>Resistor in S_k</strong>: R(S_k) = R_0 exp(S_k / S_k_ref)&mdash;dissipates energy during knowledge-dimension partition resolution. High S_k (coarse) is easy; low S_k (fine) faces exponentially increasing barriers.</li>
                  <li><strong>Capacitor in S_t</strong>: C(S_t) = C_0(1 &minus; S_t / S_t_max)&mdash;accumulates temporal information. At S_t = 0 (maximum temporal resolution), full capacitance. At S_t = S_t_max, no further storage.</li>
                  <li><strong>Inductor in S_e</strong>: L(S_e) = L_0 exp(S_e / S_e_ref)&mdash;resists changes in evolutionary direction. High evolution entropy creates trajectory inertia.</li>
                </ul>
                <p className="text-gray-300 mt-2">All three behaviors operate simultaneously because S_k, S_t, S_e are independent coordinates.</p>
              </div>

              <div className="bg-gray-900 border border-gray-700 rounded p-4 my-6 text-center font-mono text-gray-200">
                Z(&omega;, S) = R(S_k) + j&omega;L(S_e) + 1/(j&omega;C(S_t))
              </div>
              <p className="text-gray-300 leading-relaxed mb-4">The minimum impedance (resonance) occurs at &omega;_0 = 1/&radic;(LC), where only the resistive (knowledge) dimension contributes to the partition cost and the temporal and evolutionary dimensions are perfectly balanced.</p>
            </section>

            {/* Operational Regimes */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">7. Five Operational Regimes</h2>
              <p className="text-gray-300 leading-relaxed mb-4">The operational regime is classified by the Kuramoto order parameter:</p>
              <div className="bg-gray-900 border border-gray-700 rounded p-4 my-6 text-center font-mono text-gray-200">
                R = |N&sup;&minus;1; &Sigma;_j exp(i&phi;_j)|
              </div>
              <div className="overflow-x-auto my-6">
                <table className="w-full text-left border-collapse">
                  <thead><tr className="border-b border-gray-700"><th className="p-3 text-gray-400">Regime</th><th className="p-3 text-gray-400">R</th><th className="p-3 text-gray-400">Description</th><th className="p-3 text-gray-400">Structural factor S</th></tr></thead>
                  <tbody>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Phase-locked</td><td className="p-3 text-gray-300">&gt;0.95</td><td className="p-3 text-gray-300">Maximum synchronization, max throughput</td><td className="p-3 text-gray-300">M ln n</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Coherent</td><td className="p-3 text-gray-300">0.8&ndash;0.95</td><td className="p-3 text-gray-300">High synchronization, coordinated processing</td><td className="p-3 text-gray-300">R&sup2; M ln n</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Cascade</td><td className="p-3 text-gray-300">0.3&ndash;0.8</td><td className="p-3 text-gray-300">Hierarchical multi-scale processing</td><td className="p-3 text-gray-300">R M ln n / ln(1/R)</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Aperture-dominated</td><td className="p-3 text-gray-300">varies</td><td className="p-3 text-gray-300">Geometry controls flow via topology</td><td className="p-3 text-gray-300">&mdash;</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Turbulent</td><td className="p-3 text-gray-300">&le;0.3</td><td className="p-3 text-gray-300">Low sync, max entropy, exploratory</td><td className="p-3 text-gray-300">ln(1/R)</td></tr>
                  </tbody>
                </table>
              </div>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Critical Coupling)</p>
                <p className="text-gray-300 mb-2">The transition between incoherent and coherent regimes occurs at:</p>
                <p className="font-mono text-center text-gray-200">K_c = 2&sigma;_&omega;</p>
                <p className="text-gray-300 mt-2">This is a second-order (Landau) phase transition with R ~ &radic;((K &minus; K_c)/K_c) near the critical point.</p>
              </div>
              <p className="text-gray-300 leading-relaxed mb-4">Each regime satisfies an equation of state:</p>
              <div className="bg-gray-900 border border-gray-700 rounded p-4 my-6 text-center font-mono text-gray-200">PV = Nk_BT &middot; S(M, n, R)</div>
              <p className="text-gray-300 leading-relaxed mb-4">A complete processing cycle traverses the regimes: Coherent &rarr; Cascade &rarr; Phase-locked &rarr; Cascade &rarr; Turbulent &rarr; Coherent, providing both exploitation (coherent/phase-locked) and exploration (turbulent).</p>
            </section>

            <div className="my-12 bg-gray-900/50 rounded-lg p-4 border border-gray-700">
              <PharmaKuramotoChart />
            </div>

            {/* Performance */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">8. Performance and Completeness</h2>
              <div className="overflow-x-auto my-6">
                <table className="w-full text-left border-collapse">
                  <thead><tr className="border-b border-gray-700"><th className="p-3 text-gray-400">Metric</th><th className="p-3 text-gray-400">Membrane Architecture</th><th className="p-3 text-gray-400">Silicon (comparison)</th></tr></thead>
                  <tbody>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Throughput per cm&sup2;</td><td className="p-3 text-cyan-400 font-mono">~10&sup;20; ops/s</td><td className="p-3 text-gray-300">~10&sup;14; ops/s</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Energy per operation</td><td className="p-3 text-cyan-400 font-mono">~6&times;10&sup;&minus;20; J</td><td className="p-3 text-gray-300">~10&sup;&minus;15; J</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Energy efficiency ratio</td><td className="p-3 text-cyan-400 font-mono">10&sup6;&times; better per bit</td><td className="p-3 text-gray-300">baseline</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Speedup for N=10&sup6;</td><td className="p-3 text-cyan-400 font-mono">~10&sup5;&times;</td><td className="p-3 text-gray-300">baseline</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">BMD transistor on/off</td><td className="p-3 text-cyan-400 font-mono">42.1</td><td className="p-3 text-gray-300">~10&sup5; (MOSFET)</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">BMD switching time</td><td className="p-3 text-cyan-400 font-mono">&lt; 1 &mu;s</td><td className="p-3 text-gray-300">~1 ns</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Von Neumann bottleneck</td><td className="p-3 text-cyan-400">Eliminated</td><td className="p-3 text-gray-300">Present</td></tr>
                  </tbody>
                </table>
              </div>

              <p className="text-gray-300 leading-relaxed mb-4">The architecture achieves Boolean completeness (ternary NAND/NOR), arithmetic completeness (addition, multiplication on S-entropy coordinates), algorithmic completeness (any Turing-computable function via categorical navigation), and resource completeness (the membrane provides processor, memory, and interconnect as a single unified substrate).</p>
            </section>

            <div className="my-12 bg-gray-900/50 rounded-lg p-4 border border-gray-700">
              <ComputingStackChart />
            </div>

            <div className="my-12 bg-gray-900/50 rounded-lg p-4 border border-gray-700">
              <SpeedupChart />
            </div>

            {/* Path-Independent Convergence */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">9. Partition Gradient Optimization</h2>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Path-Independent Convergence)</p>
                <p className="text-gray-300 mb-2">When the sufficiency condition ||&nabla;L(&gamma;(t))|| &lt; &epsilon; is satisfied, all trajectories satisfying the gradient flow converge to the same terminus &Gamma;*, regardless of initial conditions. Convergence is exponential:</p>
                <p className="font-mono text-center text-gray-200">||&delta;(t)||&sup2; &le; ||&delta;(0)||&sup2; exp(&minus;2&alpha;&lambda;_min t)</p>
              </div>
              <p className="text-gray-300 leading-relaxed mb-4">When multiple trajectories are initiated in parallel toward the same terminus, the trajectory with shortest path length completes first. This provides automatic priority: shorter paths complete first without any explicit scheduling logic&mdash;eliminating the need for a scheduler.</p>
            </section>

            {/* Discussion */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">10. Discussion</h2>
              <p className="text-gray-300 leading-relaxed mb-4">The membrane-mediated categorical processing architecture eliminates the von Neumann bottleneck by unifying processor, memory, and interconnect in a single physical substrate. The membrane is simultaneously the computing element (lipid conformational oscillations execute partition operations), the memory (phase patterns and domain structures store information), and the interconnect (lateral diffusion and collective modes propagate signals).</p>
              <p className="text-gray-300 leading-relaxed mb-4">All results follow from three axioms with zero free parameters. The architecture provides a physical foundation for understanding why biological systems compute so efficiently: they exploit the categorical regime where computation, observation, and processing are identical operations. The categorical speedup from O(N) to O(log_3 N) is not a theoretical curiosity but an explanation for the observed computational power of biological systems&mdash;from immune recognition to neural pattern completion to protein folding.</p>
              <p className="text-gray-300 leading-relaxed mb-4">The five operational regimes provide a complete dynamics vocabulary: phase-locked for consolidation, coherent for coordinated processing, cascade for multi-scale integration, aperture-dominated for filtering, and turbulent for exploration. The regime cycling protocol provides both exploitation and exploration without external control&mdash;the system self-organizes through the Kuramoto synchronization dynamics of its constituent oscillators.</p>
            </section>
          </article>
        </Layout>
      </main>
    </>
  );
};

export default MembraneArchitecture;
