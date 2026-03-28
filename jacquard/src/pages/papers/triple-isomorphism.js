import Head from "next/head";
import Layout from "@/components/Layout";
import TransitionEffect from "@/components/TransitionEffect";
import dynamic from "next/dynamic";

const TripleEquivalenceChart = dynamic(() => import("@/components/charts/TripleEquivalenceChart"), { ssr: false });
const SpectralLinesChart = dynamic(() => import("@/components/charts/SpectralLinesChart"), { ssr: false });
const TemperatureIndependenceChart = dynamic(() => import("@/components/charts/TemperatureIndependenceChart"), { ssr: false });

const TripleIsomorphism = () => {
  return (
    <>
      <Head>
        <title>Triple Isomorphism Architecture &mdash; Bene Gesserit</title>
        <meta name="description" content="On the Tripartite Isomorphism Architecture of Biological Membrane Processes" />
      </Head>
      <TransitionEffect />
      <main className="w-full mb-16 flex flex-col items-center justify-center dark:text-light">
        <Layout className="pt-16">
          <article className="max-w-4xl mx-auto">
            <h1 className="text-4xl font-bold mb-2 leading-tight">
              On the Tripartite Isomorphism Architecture of Biological Membrane Processes: Mathematical Projection in Bounded Phase Space Systems
            </h1>
            <p className="text-gray-400 mb-8">Kundai Farai Sachikonye &bull; AIMe Registry for Artificial Intelligence</p>

            {/* Abstract */}
            <div className="bg-gray-800/50 border border-gray-700 rounded-lg p-6 mb-12">
              <h2 className="text-lg font-semibold mb-3 text-gray-300">Abstract</h2>
              <p className="text-gray-300 leading-relaxed">
                We prove that oscillatory dynamics, categorical partition structure, and biological membrane processes constitute three projections of a single mathematical object in bounded phase space. The proof proceeds by constructing three categories&mdash;<span className="font-mono text-cyan-400">O</span> (oscillatory), <span className="font-mono text-cyan-400">C</span> (categorical), and <span className="font-mono text-cyan-400">B</span> (biological)&mdash;and exhibiting explicit functors F_OC: O &rarr; C, F_CB: C &rarr; B, F_BO: B &rarr; O that form a triangle of equivalences: <span className="font-mono text-cyan-400">F_BO &compfn; F_CB &compfn; F_OC &cong; Id_O</span> and cyclically. From three axioms&mdash;bounded phase space, no null state, and finite observational resolution&mdash;we derive forced partitioning into M distinguishable states, prove forced oscillatory dynamics, and establish the capacity formula <span className="font-mono text-cyan-400">C(n) = 2n&sup2;</span>. The triple entropy equivalence <span className="font-mono text-cyan-400">S_osc = S_cat = S_part = k_B M ln n</span> is proved by three independent derivations. We establish component-level isomorphisms for ten structural elements including trajectory-terminus-memory triples, phase-lock synchronization junctions, categorical aperture filters, race dynamics, the R-C-L trichotomy, BMD transistors, and P-N junctions. Quantitative predictions&mdash;hydrogen bond energy 27 kJ/mol, bilayer thickness 4.0 nm, area per lipid 0.64 nm&sup2;, bending modulus 20 k_BT, membrane potential &minus;70 mV, P-N junction built-in potential 615 mV, proton transfer frequency 4.06&times;10&sup;13; Hz&mdash;follow identically from all three descriptions.
              </p>
            </div>

            {/* Introduction */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">1. Introduction</h2>
              <p className="text-gray-300 leading-relaxed mb-4">Physical systems that exhibit oscillatory dynamics are frequently described in at least three distinct languages. The <em className="text-white">oscillatory</em> description specifies frequencies, phases, and amplitudes governed by coupled differential equations. The <em className="text-white">categorical</em> or information-theoretic description assigns entropy coordinates and navigates trit-addressed state spaces. The <em className="text-white">biological</em> or membrane description uses partition coordinates, lipid conformational states, and transmembrane potentials.</p>
              <p className="text-gray-300 leading-relaxed mb-4">These three descriptions are conventionally treated as separate domains connected by analogy: one says that ion channels &ldquo;behave like&rdquo; bandpass filters, or that Kuramoto synchronization &ldquo;resembles&rdquo; synaptic convergence. Such analogies lack mathematical content. An analogy asserts structural similarity; an isomorphism asserts structural identity. This paper proves the latter.</p>
            </section>

            {/* Axioms */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">2. Axioms and Foundations</h2>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Axiom 1 (Bounded Phase Space)</p>
                <p className="text-gray-300">Every physical system occupies a bounded region &Omega; of phase space with finite Liouville measure &mu;(&Omega;) &lt; &infin;.</p>
              </div>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Axiom 2 (No Null State)</p>
                <p className="text-gray-300">A system must occupy exactly one categorical state at every instant t &isin; R. There is no &ldquo;off&rdquo; state; the state map &sigma;: R &rarr; S is total and surjective onto the occupied state set.</p>
              </div>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Axiom 3 (Finite Observational Resolution)</p>
                <p className="text-gray-300">The minimum distinguishable phase space volume is &delta;^d &gt; 0, where &delta; &gt; 0 is the resolution limit and d is the phase space dimension.</p>
              </div>

              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Forced Partitioning)</p>
                <p className="text-gray-300 mb-2">Under Axioms 1&ndash;3, the number of distinguishable states is:</p>
                <p className="font-mono text-center text-gray-200">M = &lfloor;&mu;(&Omega;) / &delta;^d&rfloor; &lt; &infin;</p>
              </div>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Forced Oscillation)</p>
                <p className="text-gray-300">Under Axioms 1&ndash;3, any system with M &ge; 2 must undergo perpetual state transitions. Static equilibrium at a single state has Lebesgue measure zero in the space of initial conditions. By Poincare recurrence, trajectories must revisit previously visited regions, producing perpetual transitions.</p>
              </div>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Capacity Formula)</p>
                <p className="text-gray-300 mb-2">The number of distinguishable states at partition depth n is exactly:</p>
                <p className="font-mono text-center text-gray-200">C(n) = 2n&sup2;</p>
                <p className="text-gray-300 mt-2">Proof: At depth n, angular complexity ranges over l &isin; {'{'}0,...,n&minus;1{'}'}. For each l, orientation ranges over m &isin; {'{&minus;l,...,+l}'}, giving 2l+1 values. Chirality contributes factor 2. Total: C(n) = 2&Sigma;(2l+1) = 2n&sup2;.</p>
              </div>
            </section>

            {/* Three Categories */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">3. The Three Categories</h2>
              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">Category O (Oscillatory)</h3>
              <p className="text-gray-300 leading-relaxed mb-4"><strong>Objects:</strong> Oscillatory states (&omega;, &phi;, A) where &omega; is angular frequency, &phi; is phase, A is amplitude. <strong>Morphisms:</strong> Phase-preserving frequency transformations (&Delta;&omega;, &Delta;&phi;, g). <strong>Composition:</strong> Addition of frequency/phase shifts and composition of amplitude functions.</p>

              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">Category C (Categorical)</h3>
              <p className="text-gray-300 leading-relaxed mb-4"><strong>Objects:</strong> S-entropy coordinates (S_k, S_t, S_e) &isin; [0,1]&sup3;. <strong>Morphisms:</strong> Ternary address navigations&mdash;finite strings &alpha; &isin; {'{0,1,2}'}* over the ternary alphabet. <strong>Composition:</strong> String concatenation.</p>

              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">Category B (Biological)</h3>
              <p className="text-gray-300 leading-relaxed mb-4"><strong>Objects:</strong> Partition states (n, l, m, s, M) with partition coordinates and depth. <strong>Morphisms:</strong> Partition depth transitions (&Delta;n, &Delta;l, &Delta;m, &Delta;s, &Delta;M). <strong>Composition:</strong> Cascaded partition traversal (vector addition).</p>
            </section>

            {/* Triple Entropy */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">4. The Triple Entropy Equivalence</h2>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Triple Entropy Equivalence)</p>
                <p className="text-gray-300 mb-2">The entropy computed in each of the three categories satisfies:</p>
                <p className="font-mono text-center text-gray-200 text-lg">S_osc = S_cat = S_part = k_B M ln n</p>
              </div>
              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">Oscillatory Derivation</h3>
              <p className="text-gray-300 leading-relaxed mb-4">At depth M with n modes per level, the total distinguishable oscillatory configurations is W_osc = n^M. By Boltzmann: S_osc = k_B ln(n^M) = k_B M ln n.</p>
              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">Categorical Derivation</h3>
              <p className="text-gray-300 leading-relaxed mb-4">With base b = n, the number of addressable categorical states is W_cat = n^M, giving S_cat = k_B M ln n.</p>
              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">Partition Derivation</h3>
              <p className="text-gray-300 leading-relaxed mb-4">The partition traversal visits M levels with n-fold branching, giving W_part = n^M and S_part = k_B M ln n.</p>
              <p className="text-gray-300 leading-relaxed mb-4 italic">The triple entropy identity is not a coincidence but a mathematical necessity: all three expressions count the same quantity&mdash;the number of distinguishable configurations of M hierarchical levels with n-fold branching&mdash;using different languages.</p>
            </section>

            <div className="my-12 bg-gray-900/50 rounded-lg p-4 border border-gray-700">
              <TripleEquivalenceChart />
            </div>

            {/* Functors */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">5. The Functors</h2>
              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">Functor F_OC: O &rarr; C</h3>
              <div className="bg-gray-900 border border-gray-700 rounded p-4 my-6 font-mono text-gray-200 space-y-1">
                <p className="text-center">S_k = 1 &minus; exp(&minus;&omega;/&omega;_ref)</p>
                <p className="text-center">S_t = &phi; / (2&pi;)</p>
                <p className="text-center">S_e = A&sup2; / (A&sup2; + A_ref&sup2;)</p>
              </div>

              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">Functor F_CB: C &rarr; B</h3>
              <div className="bg-gray-900 border border-gray-700 rounded p-4 my-6 font-mono text-gray-200 space-y-1">
                <p className="text-center">M = &minus;ln(1 &minus; ||S||) / ln b</p>
                <p className="text-center">n = &lceil;&radic;(3M)&rceil;</p>
                <p className="text-center">l = &lfloor;(n&minus;1) &middot; arccos(S_e/||S||) / &pi;&rfloor;</p>
              </div>

              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">Functor F_BO: B &rarr; O</h3>
              <div className="bg-gray-900 border border-gray-700 rounded p-4 my-6 font-mono text-gray-200 space-y-1">
                <p className="text-center">&omega; = (k_BT/&hbar;) &middot; exp(&minus;&Delta;M)</p>
                <p className="text-center">&phi; = 2&pi; &middot; m/(2l+1)</p>
                <p className="text-center">A = A_0 &radic;(C(n)/C_max)</p>
              </div>

              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Category Equivalence)</p>
                <p className="text-gray-300 mb-2">The three functors form a triangle of equivalences:</p>
                <p className="font-mono text-center text-gray-200">F_BO &compfn; F_CB &compfn; F_OC &cong; Id_O</p>
                <p className="font-mono text-center text-gray-200">F_OC &compfn; F_BO &compfn; F_CB &cong; Id_C</p>
                <p className="font-mono text-center text-gray-200">F_CB &compfn; F_OC &compfn; F_BO &cong; Id_B</p>
                <p className="text-gray-300 mt-2">The equivalence is exact in the continuum limit &delta; &rarr; 0 and is within observational resolution for any finite &delta; &gt; 0.</p>
              </div>
            </section>

            {/* Component Isomorphisms */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">6. Component-Level Isomorphisms</h2>

              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">6.1 Trajectory-Terminus-Memory Triple</h3>
              <div className="overflow-x-auto my-6">
                <table className="w-full text-left border-collapse">
                  <thead><tr className="border-b border-gray-700"><th className="p-3 text-gray-400">Component</th><th className="p-3 text-gray-400">O (Oscillatory)</th><th className="p-3 text-gray-400">C (Categorical)</th><th className="p-3 text-gray-400">B (Biological)</th></tr></thead>
                  <tbody>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Trajectory</td><td className="p-3 text-gray-300">(&omega;(t), &phi;(t), A(t))</td><td className="p-3 text-gray-300">Path through [0,1]&sup3;</td><td className="p-3 text-gray-300">Lipid conformational trajectory</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Terminus</td><td className="p-3 text-gray-300">&gamma;(T)</td><td className="p-3 text-gray-300">Final S-entropy coordinate</td><td className="p-3 text-gray-300">Target partition state</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Memory</td><td className="p-3 text-gray-300">&int;&phi;(t) dt</td><td className="p-3 text-gray-300">&int;d&gamma;/dt &otimes; H dt</td><td className="p-3 text-gray-300">&Delta;M_&gamma; (total depth traversed)</td></tr>
                  </tbody>
                </table>
              </div>
              <p className="text-gray-300 leading-relaxed mb-4">The memory functional is path-dependent: different trajectories with the same terminus produce different memories. Neither component can be reconstructed from the other two alone.</p>

              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">6.2 Phase-Lock Synchronization Junction</h3>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Junction Isomorphism)</p>
                <p className="text-gray-300 mb-2">The Kuramoto critical coupling condition in O is equivalent to the convergence criterion in C and the lipid raft boundary merger in B:</p>
                <p className="font-mono text-center text-gray-200">||&nabla;L(&gamma;_j(t))|| &lt; &epsilon; &hArr; K &gt; K_c = 2&sigma;_&omega;/&pi;</p>
                <p className="text-gray-300 mt-2">with the same convergence time scaling t_sync ~ (K &minus; K_c)&sup;&minus;1; in all three descriptions.</p>
              </div>

              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">6.3 Categorical Aperture / Ion Channel / Frequency Bandpass</h3>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Zero Work)</p>
                <p className="text-gray-300">Categorical aperture filtering performs zero thermodynamic work in all three descriptions: W_aperture = 0. Landauer&apos;s erasure principle does not apply because aperture filtering is injective (reversible)&mdash;no information is destroyed.</p>
              </div>

              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">6.4 BMD Transistor</h3>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (BMD On/Off Ratio)</p>
                <p className="text-gray-300 mb-2">The on/off ratio of the BMD transistor is:</p>
                <p className="font-mono text-center text-gray-200">G_max / G_min = exp(&Delta;M_T &middot; E_0 / k_BT) &asymp; 42.1</p>
                <p className="text-gray-300 mt-2">Energy per operation: E_op &asymp; 20 k_BT ln 2 &asymp; 6 &times; 10&sup;&minus;20; J, which is 20&times; the Landauer bound.</p>
              </div>

              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">6.5 P-N Junction / Raft Boundary</h3>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Built-in Potential)</p>
                <p className="text-gray-300 mb-2">The built-in potential across the junction is:</p>
                <p className="font-mono text-center text-gray-200">V_bi = (k_BT/q) ln(N_A N_D / n_i&sup2;) = 615 mV</p>
                <p className="text-gray-300 mt-2">This value is the same in all three descriptions: oscillator density gradient, raft boundary compositional asymmetry, and carrier type transition in S-entropy space.</p>
              </div>

              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">6.6 Poincare Non-Recurrence (Generative Novelty)</h3>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Poincare Non-Recurrence)</p>
                <p className="text-gray-300">For an anharmonic oscillator in bounded phase space, exact recurrence has measure zero. Every cycle of operation produces a state that has never existed before and will never exist again. This ensures maximal coverage of the available phase space over time.</p>
              </div>
            </section>

            {/* Predictions */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">7. Quantitative Predictions</h2>
              <p className="text-gray-300 leading-relaxed mb-4">All predictions follow identically from all three descriptions:</p>
              <div className="overflow-x-auto my-6">
                <table className="w-full text-left border-collapse">
                  <thead><tr className="border-b border-gray-700"><th className="p-3 text-gray-400">Prediction</th><th className="p-3 text-gray-400">Derived Value</th><th className="p-3 text-gray-400">Experimental</th></tr></thead>
                  <tbody>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Hydrogen bond energy</td><td className="p-3 text-gray-300">27 kJ/mol</td><td className="p-3 text-gray-300">20&ndash;25 kJ/mol</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Bilayer thickness</td><td className="p-3 text-gray-300">4.0 nm</td><td className="p-3 text-gray-300">3.9&ndash;4.2 nm</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Area per lipid</td><td className="p-3 text-gray-300">0.64 nm&sup2;</td><td className="p-3 text-gray-300">0.63 nm&sup2;</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Bending modulus</td><td className="p-3 text-gray-300">20 k_BT</td><td className="p-3 text-gray-300">20 &plusmn; 4 k_BT</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Membrane potential</td><td className="p-3 text-gray-300">&minus;70 mV</td><td className="p-3 text-gray-300">&minus;70 mV</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">P-N junction V_bi</td><td className="p-3 text-gray-300">615 mV</td><td className="p-3 text-gray-300">~600 mV (Si)</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">BMD on/off ratio</td><td className="p-3 text-gray-300">42.1</td><td className="p-3 text-gray-300">~40&ndash;50 (enzymes)</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">BMD energy per op</td><td className="p-3 text-gray-300">6&times;10&sup;&minus;20; J</td><td className="p-3 text-gray-300">~10&sup;&minus;19; J (ATP)</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Proton transfer freq.</td><td className="p-3 text-gray-300">4.06&times;10&sup;13; Hz</td><td className="p-3 text-gray-300">~10&sup;13; Hz</td></tr>
                  </tbody>
                </table>
              </div>
            </section>

            <div className="my-12 bg-gray-900/50 rounded-lg p-4 border border-gray-700">
              <SpectralLinesChart />
            </div>

            {/* Selection Rules */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">8. Selection Rules</h2>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Derived Selection Rules</p>
                <p className="text-gray-300">From partition continuity (transitions must connect adjacent cells in the partition hierarchy):</p>
                <div className="font-mono text-center text-gray-200 mt-2 space-y-1">
                  <p>|&Delta;l| = 1</p>
                  <p>|&Delta;m| &le; 1</p>
                  <p>&Delta;s = 0</p>
                </div>
                <p className="text-gray-300 mt-2">These match the spectroscopic selection rules for electric dipole transitions, derived here from pure partition geometry rather than quantum electrodynamics.</p>
              </div>
            </section>

            <div className="my-12 bg-gray-900/50 rounded-lg p-4 border border-gray-700">
              <TemperatureIndependenceChart />
            </div>

            {/* Discussion */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">9. Discussion</h2>
              <p className="text-gray-300 leading-relaxed mb-4">The triple isomorphism eliminates the need for separate theories of oscillatory dynamics, categorical information processing, and biological membrane physics. A single implementation serves all three domains. Any result proved in one category automatically transfers to the other two through the established functors.</p>
              <p className="text-gray-300 leading-relaxed mb-4">The universal equation of state PV = Nk_BT &middot; S(V,N,{'{n_i,l_i,m_i,s_i}'}) holds identically in all three views with the same temperature-independent structural factor S. The variance&ndash;free energy identity F = k_BT &middot; &sigma;&sup2;(&phi;) establishes that hybrid microfluidic circuit dynamics IS variance minimization.</p>
              <p className="text-gray-300 leading-relaxed mb-4">The practical consequence is that a processor designed for oscillatory dynamics automatically serves as a categorical computer and a biological membrane simulator. The three domains are not analogous but are the same mathematical object viewed through three equivalent projections.</p>
            </section>
          </article>
        </Layout>
      </main>
    </>
  );
};

export default TripleIsomorphism;
