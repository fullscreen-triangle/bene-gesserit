import Head from "next/head";
import { motion } from "framer-motion";
import Layout from "@/components/Layout";
import AnimatedText from "@/components/AnimatedText";

const container = {
  hidden: { opacity: 0 },
  visible: {
    opacity: 1,
    transition: { staggerChildren: 0.12 },
  },
};

const item = {
  hidden: { opacity: 0, y: 20 },
  visible: { opacity: 1, y: 0, transition: { duration: 0.5 } },
};

const sections = [
  {
    title: "Overview",
    content:
      "The Bene Gesserit framework provides a rigorous mathematical foundation for categorical processing on amphiphilic bilayer substrates. It unifies lipid bilayer thermodynamics, categorical computation, and bounded phase space dynamics under a single algebraic framework, enabling novel approaches to computation grounded in membrane physics.",
  },
  {
    title: "Core Concepts",
    content:
      "At its heart, the framework rests on three pillars: the derivation of computational primitives from lipid bilayer dynamics, a formal categorical converter mapping between physical and mathematical domains, and the triple isomorphism theorem establishing equivalence across all three foundational domains. These pillars support a membrane cognitive architecture capable of robust, adaptive, parallel processing.",
  },
  {
    title: "Architecture",
    content:
      "The framework is structured as a layered system. The physical layer models amphiphilic bilayer behavior using thermodynamic potentials and statistical mechanics. The categorical layer provides the abstract mathematical machinery for reasoning about computations. The bridge layer implements the functorial converter that maintains coherence between physical and categorical representations.",
  },
  {
    title: "API Reference",
    content:
      "The framework exposes interfaces for defining membrane configurations, executing categorical transformations, querying phase space trajectories, and composing complex computations from primitive operations. Full API documentation with usage examples will be published alongside the reference implementation.",
  },
  {
    title: "Getting Started",
    content:
      "To begin working with the Bene Gesserit framework, clone the repository from GitHub and follow the setup instructions in the README. The framework requires a standard scientific computing environment and provides both programmatic and interactive interfaces for exploring membrane-based categorical computation.",
  },
];

export default function Framework() {
  return (
    <>
      <Head>
        <title>Framework — Bene Gesserit</title>
        <meta
          name="description"
          content="Bene Gesserit framework architecture, API reference, and usage guide."
        />
      </Head>
      <main className="flex w-full flex-col items-center justify-center dark:text-light">
        <Layout className="pt-16">
          <AnimatedText
            text="Framework"
            className="!text-6xl mb-16 lg:!text-5xl sm:!text-4xl"
          />
          <motion.div
            className="space-y-8"
            variants={container}
            initial="hidden"
            animate="visible"
          >
            {sections.map((section, index) => (
              <motion.div
                key={index}
                variants={item}
                className="rounded-xl border border-solid border-dark/20 dark:border-light/20
                bg-light dark:bg-dark p-8"
              >
                <h2 className="text-2xl font-bold mb-4 text-dark dark:text-light">
                  {section.title}
                </h2>
                <p className="text-base text-dark/70 dark:text-light/70 leading-relaxed">
                  {section.content}
                </p>
              </motion.div>
            ))}
          </motion.div>
        </Layout>
      </main>
    </>
  );
}
