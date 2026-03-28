import Head from "next/head";
import Layout from "@/components/Layout";
import TransitionEffect from "@/components/TransitionEffect";
import dynamic from "next/dynamic";

const SpeedupChart = dynamic(() => import("@/components/charts/SpeedupChart"), { ssr: false });
const CapacityChart = dynamic(() => import("@/components/charts/CapacityChart"), { ssr: false });

const CategoricalConverter = () => {
  return (
    <>
      <Head>
        <title>Membrane-Mediated Categorical Conversion &mdash; Bene Gesserit</title>
        <meta name="description" content="Non-Quantizable Processors via Partition Extinction on Binary Substrates" />
      </Head>
      <TransitionEffect />
      <main className="w-full mb-16 flex flex-col items-center justify-center dark:text-light">
        <Layout className="pt-16">
          <article className="max-w-4xl mx-auto">
            <h1 className="text-4xl font-bold mb-2 leading-tight">
              Membrane-Mediated Categorical Conversion: Non-Quantizable Processors via Partition Extinction on Binary Substrates
            </h1>
            <p className="text-gray-400 mb-8">Kundai Farai Sachikonye &bull; AIMe Registry for Artificial Intelligence</p>

            {/* Abstract */}
            <div className="bg-gray-800/50 border border-gray-700 rounded-lg p-6 mb-12">
              <h2 className="text-lg font-semibold mb-3 text-gray-300">Abstract</h2>
              <p className="text-gray-300 leading-relaxed">
                We introduce the categorical converter&mdash;a processor architecture that uses virtual membranes to convert binary digital computation into continuous categorical computation without intermediate analog stages. The architecture rests on three pillars derived from first principles. First, we prove that any binary processor is a state counter operating at base b = 2, generating entropy <span className="font-mono text-cyan-400">&Delta;S = k_B ln 2</span> per clock transition, and that this state counting is isomorphic to partition traversal in bounded phase space. Second, we establish the Partition Extinction Theorem: when temperature drops below a critical threshold T_c for a group of binary states, those states become categorically indistinguishable and merge into a single continuous categorical state&mdash;the binary quantization dissolves. Third, we prove that virtual membranes&mdash;boundary conditions in partition space implemented through oscillatory dynamics&mdash;serve as transducers between the binary and categorical regimes, with the capacity formula <span className="font-mono text-cyan-400">C(n) = 2n&sup2;</span> providing the natural bijection between binary state counts and ternary partition depth. The resulting processor has no fixed quantizable unit: its computational resolution is a continuous function of temperature. We derive the complete thermodynamics of categorical conversion, prove that the conversion is reversible with efficiency bounded by Carnot, and establish that a membrane of area A implements a categorical processor with <span className="font-mono text-cyan-400">~10&sup;20;</span> partition operations per second. The framework unifies digital and analog computation as limiting cases of a single partition-depth continuum and provides a physical foundation for biological computation. We present the categorical speedup: problems requiring O(n) binary steps require only O(log n) categorical navigations.
              </p>
            </div>

            {/* Introduction */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">1. Introduction</h2>
              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">The Quantization Assumption</h3>
              <p className="text-gray-300 leading-relaxed mb-4">All conventional processors rest on a single assumption: computation requires quantization. The bit&mdash;a system in one of exactly two distinguishable states&mdash;is the atomic unit of digital computation. Logic gates combine bits through Boolean operations. Clock cycles synchronize state transitions. The entire edifice of computer science is built on the premise that information must be discretized into countable, combinable units before it can be processed.</p>
              <p className="text-gray-300 leading-relaxed mb-4">This assumption has been spectacularly productive, yet it is not a law of nature&mdash;it is an engineering choice. Nature computes without bits. Protein folding, neural signaling, gene regulation, and cellular decision-making all proceed through continuous dynamics, not discrete logic. We prove that binary digital computation and continuous categorical computation are limiting cases of a single framework governed by partition depth, and that membranes serve as the natural transducers between these regimes.</p>
            </section>

            {/* Binary Processor */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">2. The Binary Processor as State Counter</h2>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Binary Processor Identity)</p>
                <p className="text-gray-300">Any binary digital processor with clock frequency f and word length w is a state counter operating at base b = 2 with:</p>
                <ul className="list-disc list-inside text-gray-300 mt-2 space-y-1">
                  <li>Partition depth per clock cycle: <span className="font-mono text-cyan-400">M_cycle = w</span></li>
                  <li>Entropy production per cycle: <span className="font-mono text-cyan-400">&Delta;S_cycle = w &middot; k_B ln 2</span></li>
                  <li>State counting rate: <span className="font-mono text-cyan-400">R = f &middot; w</span> partition transitions per second</li>
                  <li>Total partition depth traversed in time t: <span className="font-mono text-cyan-400">M(t) = f &middot; w &middot; t</span></li>
                </ul>
              </div>
              <p className="text-gray-300 leading-relaxed mb-4">A single bit is a partition at depth n = 1 with capacity C(1) = 2. For a 64-bit processor at f = 3 GHz: R = 1.92 &times; 10&sup;11; partition transitions/s, and entropy production rate ~1.84 &times; 10&sup;&minus;12; W/K.</p>

              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Binary Quantization Constraint)</p>
                <p className="text-gray-300">A binary processor is constrained to partition depth increments of exactly &Delta;M = 1 per elementary operation. Encoding N states requires &lceil;log_2 N&rceil; bits, wasting 2^&lceil;log_2 N&rceil; &minus; N states. The average waste fraction:</p>
                <p className="font-mono text-center text-gray-200 mt-2">&langle;W&rangle; = 1 &minus; 1/ln 2 &asymp; 0.31</p>
                <p className="text-gray-300 mt-2">Binary processors forfeit approximately 31% of their partition capacity when encoding non-power-of-2 structures.</p>
              </div>

              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">The Landauer Limit and Binary Entropy</h3>
              <div className="bg-gray-900 border border-gray-700 rounded p-4 my-6 text-center font-mono text-gray-200">
                <p>E_min = k_BT ln 2 &asymp; 2.87 &times; 10&sup;&minus;21; J at T = 300 K</p>
                <p className="text-sm text-gray-500 mt-1">Current processors dissipate ~10&sup;&minus;15; J per transition&mdash;10&sup6; times the Landauer limit.</p>
              </div>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Corollary (Base-Dependent Landauer Limit)</p>
                <p className="text-gray-300">For a processor operating at base b: E_min(b) = k_BT ln b. The information per unit energy is I/E = 1/(k_BT ln 2) = const. The Landauer efficiency is base-independent. The advantage of ternary encoding lies in structural matching, not energy efficiency.</p>
              </div>
            </section>

            {/* Categorical Conversion */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">3. The Categorical Conversion Theorem</h2>
              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">The Effective Partition Base</h3>
              <div className="bg-gray-900 border border-gray-700 rounded p-4 my-6 text-center font-mono text-gray-200">
                b_eff(T) = 1 + (b_max &minus; 1)(1 &minus; e^(&minus;&Delta;E / k_BT))
              </div>
              <p className="text-gray-300 leading-relaxed mb-4">This is the central quantity of the paper. It defines a continuous variable that interpolates between: b_eff = 1 (partition extinction), b_eff = 2 (binary), b_eff = 3 (ternary, the natural base of S-entropy coordinates), and b_eff &rarr; &infin; (fully continuous categorical computation).</p>

              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Categorical Conversion)</p>
                <p className="text-gray-300 mb-2">A binary processor with state energy spacing &Delta;E can be converted to a categorical processor operating at effective base b_eff by controlling the temperature. The conversion passes through three regimes:</p>
                <ul className="list-disc list-inside text-gray-300 mt-2 space-y-1">
                  <li><strong>Binary regime</strong> (T &gt;&gt; &Delta;E/k_B): b_eff &asymp; b_max. Standard digital computation.</li>
                  <li><strong>Ternary regime</strong> (T &asymp; T_3): b_eff = 3. Optimal categorical computation.</li>
                  <li><strong>Extinction regime</strong> (T &lt;&lt; &Delta;E/k_B): b_eff &rarr; 1. Partition operations become undefined.</li>
                </ul>
              </div>
              <p className="text-gray-300 leading-relaxed mb-4">For a 64-bit silicon processor with &Delta;E = 0.1 eV: T_3 &asymp; 10&sup;36; K. This absurd temperature confirms that conventional silicon processors cannot be thermally converted to categorical mode. Categorical conversion requires a substrate where &Delta;E is comparable to k_BT.</p>

              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">The Biological Operating Point</h3>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Proposition (Biological Categorical Conversion)</p>
                <p className="text-gray-300">Lipid membranes operate naturally at the ternary crossover temperature. The gauche energy &Delta;E_g &asymp; 2.1 kJ/mol gives a collective ternary crossover T_3 that brackets 310&ndash;500 K when including interchain cooperative corrections. Biological membranes operate at or slightly below T_3, in the regime where conformational states are partially distinguishable&mdash;the categorical computing regime.</p>
              </div>
              <p className="text-gray-300 leading-relaxed mb-4 italic">The biological operating temperature is not arbitrary. It is tuned to place the lipid membrane at the ternary crossover&mdash;the regime of maximum categorical computational power. Evolution did not &ldquo;choose&rdquo; 37&deg;C for protein stability alone; it chose it because this temperature places the membrane processor at optimal partition resolution.</p>
            </section>

            {/* Membrane Transducer */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">4. The Membrane as Transducer</h2>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Membrane Transduction)</p>
                <p className="text-gray-300 mb-2">A virtual membrane separating a binary region (b = 2) from a categorical region (b = b_c) performs base conversion with information conservation:</p>
                <p className="font-mono text-center text-gray-200 mt-2">I_binary = I_categorical</p>
                <p className="font-mono text-center text-gray-200">M_cat = n_bits / log_2(b_c)</p>
                <p className="text-gray-300 mt-2">For binary-to-ternary conversion: M_ternary = n_bits / 1.585 &asymp; 0.631 &middot; n_bits. A 64-bit binary word converts to ternary partition depth M &asymp; 40.4. The information content is identical; the representation is more compact.</p>
              </div>

              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">The Capacity Formula as Base Converter</h3>
              <div className="overflow-x-auto my-6">
                <table className="w-full text-left border-collapse">
                  <thead><tr className="border-b border-gray-700"><th className="p-3 text-gray-400">n</th><th className="p-3 text-gray-400">C(n)</th><th className="p-3 text-gray-400">C_tot(n)</th><th className="p-3 text-gray-400">log_2 C_tot</th><th className="p-3 text-gray-400">log_3 C_tot</th></tr></thead>
                  <tbody>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">1</td><td className="p-3 text-gray-300">2</td><td className="p-3 text-gray-300">2</td><td className="p-3 text-gray-300">1.00</td><td className="p-3 text-gray-300">0.63</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">2</td><td className="p-3 text-gray-300">8</td><td className="p-3 text-gray-300">10</td><td className="p-3 text-gray-300">3.32</td><td className="p-3 text-gray-300">2.10</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">5</td><td className="p-3 text-gray-300">50</td><td className="p-3 text-gray-300">110</td><td className="p-3 text-gray-300">6.78</td><td className="p-3 text-gray-300">4.28</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">10</td><td className="p-3 text-gray-300">200</td><td className="p-3 text-gray-300">770</td><td className="p-3 text-gray-300">9.59</td><td className="p-3 text-gray-300">6.05</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">20</td><td className="p-3 text-gray-300">800</td><td className="p-3 text-gray-300">5740</td><td className="p-3 text-gray-300">12.49</td><td className="p-3 text-gray-300">7.88</td></tr>
                  </tbody>
                </table>
              </div>
              <p className="text-gray-300 leading-relaxed mb-4">The ratio log_3 / log_2 = 1/log_2(3) = 0.631 is constant for all n. The capacity formula provides a perfect base conversion because it encodes the same partition structure in both bases.</p>
            </section>

            <div className="my-12 bg-gray-900/50 rounded-lg p-4 border border-gray-700">
              <CapacityChart />
            </div>

            {/* Non-Quantizable Processor */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">5. The Non-Quantizable Processor</h2>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Non-Quantizability)</p>
                <p className="text-gray-300 mb-2">The categorical processor has no fixed computational unit. Its partition depth increment per operation is:</p>
                <p className="font-mono text-center text-gray-200">&delta;M(T) = 1 / ln b_eff(T)</p>
                <p className="text-gray-300 mt-2">which is a continuous function of temperature. At no temperature is there a &ldquo;natural unit&rdquo; of computation.</p>
              </div>
              <div className="overflow-x-auto my-6">
                <table className="w-full text-left border-collapse">
                  <thead><tr className="border-b border-gray-700"><th className="p-3 text-gray-400">Regime</th><th className="p-3 text-gray-400">T/T_c</th><th className="p-3 text-gray-400">b_eff</th><th className="p-3 text-gray-400">&delta;M</th></tr></thead>
                  <tbody>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Binary</td><td className="p-3 text-gray-300">~10</td><td className="p-3 text-gray-300">2</td><td className="p-3 text-gray-300">1.443</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Ternary</td><td className="p-3 text-gray-300">~3</td><td className="p-3 text-gray-300">3</td><td className="p-3 text-gray-300">0.910</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Pentary</td><td className="p-3 text-gray-300">~2</td><td className="p-3 text-gray-300">5</td><td className="p-3 text-gray-300">0.621</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Near extinction</td><td className="p-3 text-gray-300">~1.01</td><td className="p-3 text-gray-300">100</td><td className="p-3 text-gray-300">0.217</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Continuous</td><td className="p-3 text-gray-300">&rarr;1</td><td className="p-3 text-gray-300">&rarr;&infin;</td><td className="p-3 text-gray-300">&rarr;0</td></tr>
                  </tbody>
                </table>
              </div>

              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">S-Entropy Navigation</h3>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (S-Entropy Navigation)</p>
                <p className="text-gray-300 mb-2">Computation in the categorical processor is navigation through S-entropy space. The computational state S(t) evolves according to:</p>
                <p className="font-mono text-center text-gray-200">dS/dt = &minus;&nabla;_S F(S) + &xi;(t)</p>
                <p className="text-gray-300 mt-2">where F(S) is the categorical free energy landscape encoding the problem, and &xi;(t) is thermal noise. The solution is found when |S &minus; S*| &lt; &epsilon;.</p>
              </div>
            </section>

            {/* Thermodynamics */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">6. Categorical Thermodynamics</h2>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Conversion Entropy)</p>
                <p className="font-mono text-center text-gray-200">&Delta;S_conversion = N k_B (ln b_c &minus; ln 2)</p>
                <p className="text-gray-300 mt-2">For ternary (b_c = 3): &Delta;S = N k_B ln(3/2) &gt; 0. The conversion to higher base increases entropy and is thermodynamically spontaneous at all temperatures.</p>
              </div>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Categorical Carnot Bound)</p>
                <p className="font-mono text-center text-gray-200">&eta;_cat &le; 1 &minus; T_extinction / T_binary</p>
                <p className="text-gray-300 mt-2">For a biological membrane operating between T_binary &asymp; 500 K (fully fluid) and T_extinction &asymp; 270 K (gel phase): &eta;_max = 0.46. At the biological operating point T = 310 K: &eta;_bio &asymp; 0.38. This is remarkably close to the actual efficiency of ATP synthesis (~40%), suggesting that biological energy conversion operates near the categorical Carnot limit.</p>
              </div>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Computational Arrow)</p>
                <p className="text-gray-300">Categorical computation has an intrinsic direction: from higher to lower partition depth. This direction is the arrow of computation, which is equivalent to the thermodynamic arrow of time. The passage of time is the accumulation of partition transitions; computation is the directed accumulation toward a solution.</p>
              </div>
            </section>

            {/* Categorical Speedup */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">7. The Categorical Speedup</h2>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Categorical Speedup)</p>
                <p className="text-gray-300">For a search problem in a space of N states:</p>
                <ul className="list-disc list-inside text-gray-300 mt-2 space-y-1">
                  <li>Binary sequential search: <span className="font-mono text-cyan-400">O(N)</span> partition transitions</li>
                  <li>Binary parallel search: <span className="font-mono text-cyan-400">O(N/p)</span> transitions with p processors</li>
                  <li>Categorical navigation: <span className="font-mono text-cyan-400">O(log N)</span> partition transitions</li>
                </ul>
                <p className="text-gray-300 mt-2">The categorical speedup is exponential: N / log N.</p>
              </div>
              <p className="text-gray-300 leading-relaxed mb-4">Starting from the full volume V_0 = 1 (unit cube in S-entropy space), after k steps the remaining volume is V_k = V_0 / 2^k. The solution is found when V_k &lt; &epsilon;&sup3;. For N states uniformly distributed in [0,1]&sup3;:</p>
              <div className="bg-gray-900 border border-gray-700 rounded p-4 my-6 text-center font-mono text-gray-200">
                <p>k = log_2 N</p>
                <p className="mt-2">Speedup for N = 10&sup6;: 10&sup6; / 20 = 50,000&times;</p>
                <p>Speedup for N = 10&sup;12;: 10&sup;12; / 40 = 2.5 &times; 10&sup;10;&times;</p>
              </div>

              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">The &epsilon;-Boundary and Precision</h3>
              <div className="bg-gray-900 border border-gray-700 rounded p-4 my-6 text-center font-mono text-gray-200">
                <p>|S &minus; S*| &ge; &epsilon; = G / N^(1/3)</p>
                <p className="mt-2">G (Godelian residue) &asymp; 1.2 &times; 10&sup;&minus;3; at T = 310 K</p>
                <p className="mt-2">&epsilon; &asymp; 1.2 &times; 10&sup;&minus;5; for N = 10&sup6; (~17 bits of precision)</p>
              </div>
            </section>

            <div className="my-12 bg-gray-900/50 rounded-lg p-4 border border-gray-700">
              <SpeedupChart />
            </div>

            {/* Biological Implementation */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">8. Biological Implementation</h2>
              <div className="border-l-4 border-cyan-500 bg-gray-800/30 p-4 my-6 rounded-r">
                <p className="font-semibold text-cyan-400 mb-2">Theorem (Biological Categorical Processing)</p>
                <p className="text-gray-300">A eukaryotic cell membrane of area A &asymp; 1000 &mu;m&sup2; implements a categorical processor with:</p>
                <ul className="list-disc list-inside text-gray-300 mt-2 space-y-1">
                  <li>Processing rate: <span className="font-mono text-cyan-400">R &asymp; 10&sup;20; ops/s</span></li>
                  <li>Effective base: <span className="font-mono text-cyan-400">b_eff &asymp; 2.5&ndash;3 at T = 310 K</span></li>
                  <li>Partition depth per cycle: <span className="font-mono text-cyan-400">M &asymp; 1.3&ndash;1.6</span></li>
                  <li>Computational precision: <span className="font-mono text-cyan-400">~17 bits equivalent</span></li>
                </ul>
                <p className="text-gray-300 mt-2">This exceeds the world&apos;s fastest supercomputer (~10&sup;18; FLOPS) by three orders of magnitude.</p>
              </div>

              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">Lipid Rafts as Computational Domains</h3>
              <p className="text-gray-300 leading-relaxed mb-4">Lipid rafts function as categorical coprocessors with enhanced precision at the cost of reduced processing rate. In the raft domain: conformational partition depth is reduced, b_eff is lower, processing rate decreases, but precision increases. The raft trades speed for accuracy. This is the physical basis of signal transduction: signaling proteins concentrate in rafts because the higher precision enables reliable detection of weak signals.</p>
              <div className="bg-gray-900 border border-gray-700 rounded p-4 my-6 text-center font-mono text-gray-200">R_raft &lt; R_bulk &nbsp;&nbsp;but&nbsp;&nbsp; &epsilon;_raft &lt; &epsilon;_bulk</div>

              <h3 className="text-xl font-semibold mb-3 mt-6 text-gray-200">Action Potentials as Categorical Computation Cycles</h3>
              <p className="text-gray-300 leading-relaxed mb-4">The neuronal action potential is a complete categorical computation cycle: (1) Input via synaptic signals, (2) Threshold crossing as a bifurcation in F(S), (3) Navigation from resting to peak in S-entropy space, (4) Output propagation down the axon, (5) Reset during refractory period. The total categorical depth per action potential: M_AP &asymp; 1.5 &times; 10&sup8; (~5 &times; 10&sup7; bits equivalent).</p>
            </section>

            {/* Discussion */}
            <section className="mb-12">
              <h2 className="text-2xl font-bold mb-4 border-b border-gray-700 pb-2">9. Discussion: Unification of Digital and Analog</h2>
              <div className="overflow-x-auto my-6">
                <table className="w-full text-left border-collapse">
                  <thead><tr className="border-b border-gray-700"><th className="p-3 text-gray-400">Property</th><th className="p-3 text-gray-400">Digital</th><th className="p-3 text-gray-400">Categorical</th><th className="p-3 text-gray-400">Analog</th></tr></thead>
                  <tbody>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Base b</td><td className="p-3 text-gray-300">2</td><td className="p-3 text-gray-300">b_eff(T)</td><td className="p-3 text-gray-300">&infin;</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Resolution</td><td className="p-3 text-gray-300">&Delta;M = 1</td><td className="p-3 text-gray-300">&delta;M(T)</td><td className="p-3 text-gray-300">&delta;M &rarr; 0</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Noise sensitivity</td><td className="p-3 text-gray-300">Low</td><td className="p-3 text-gray-300">Moderate</td><td className="p-3 text-gray-300">High</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Precision per step</td><td className="p-3 text-gray-300">1 bit</td><td className="p-3 text-gray-300">log_2 b_eff bits</td><td className="p-3 text-gray-300">&infin; (ideal)</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Scaling</td><td className="p-3 text-gray-300">O(N)</td><td className="p-3 text-gray-300">O(log N)</td><td className="p-3 text-gray-300">Undefined</td></tr>
                    <tr className="border-b border-gray-800"><td className="p-3 text-gray-300">Physical substrate</td><td className="p-3 text-gray-300">Silicon</td><td className="p-3 text-gray-300">Membranes</td><td className="p-3 text-gray-300">Continuous fields</td></tr>
                  </tbody>
                </table>
              </div>
              <p className="text-gray-300 leading-relaxed mb-4">The categorical regime occupies the sweet spot: sufficient resolution for precise computation (~1.5 bits/step), low enough noise for reliable operation, and logarithmic scaling through S-entropy navigation. Biology uses membranes for computation because they naturally operate at the ternary crossover, provide massive parallelism (~10&sup9; lipid oscillators per cell), enable adaptive resolution through temperature fluctuations, and offer domain specialization via lipid rafts.</p>
            </section>
          </article>
        </Layout>
      </main>
    </>
  );
};

export default CategoricalConverter;
