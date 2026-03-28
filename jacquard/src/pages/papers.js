import Head from "next/head";
import Link from "next/link";
import { motion } from "framer-motion";
import Layout from "@/components/Layout";
import AnimatedText from "@/components/AnimatedText";

const papers = [
  {
    slug: "lipid-derivation",
    title: "Lipid Derivation",
    abstract:
      "This paper establishes the mathematical foundations for deriving computational primitives from lipid bilayer dynamics. By formalizing the thermodynamic properties of amphiphilic molecules in bounded phase spaces, we demonstrate that lipid membrane configurations naturally encode categorical structures suitable for information processing. The derivation proceeds from first principles of statistical mechanics through to constructive proofs of computational universality.",
  },
  {
    slug: "categorical-converter",
    title: "Categorical Converter",
    abstract:
      "We present a formal categorical converter that maps between lipid bilayer phase states and computational categories. The converter establishes functorial relationships between the physical dynamics of membrane systems and abstract categorical constructions, enabling rigorous translation of biological membrane computations into mathematical frameworks. This bridge allows tools from category theory to be applied directly to the analysis of amphiphilic substrate processing.",
  },
  {
    slug: "triple-isomorphism",
    title: "Triple Isomorphism",
    abstract:
      "This paper proves a triple isomorphism theorem connecting three fundamental domains: lipid bilayer thermodynamics, categorical computation, and bounded phase space dynamics. We show that the algebraic structure of membrane phase transitions is isomorphic to both a specific class of categorical transformations and a family of bounded dynamical systems. This three-way correspondence provides a unified mathematical foundation for the Bene Gesserit framework.",
  },
  {
    slug: "membrane-architecture",
    title: "Membrane Cognitive Architecture",
    abstract:
      "We propose a cognitive architecture grounded in the physical principles of biological membranes. The architecture leverages the self-organizing properties of amphiphilic bilayers to implement categorical processing in bounded phase spaces. We describe the architectural components, their interactions, and demonstrate how membrane-based computation achieves properties including robustness, adaptability, and efficient parallel processing that emerge naturally from the underlying physics.",
  },
];

const container = {
  hidden: { opacity: 0 },
  visible: {
    opacity: 1,
    transition: {
      staggerChildren: 0.15,
    },
  },
};

const item = {
  hidden: { opacity: 0, y: 30 },
  visible: {
    opacity: 1,
    y: 0,
    transition: { duration: 0.6 },
  },
};

export default function Papers() {
  return (
    <>
      <Head>
        <title>Papers — Bene Gesserit</title>
        <meta
          name="description"
          content="Research papers underpinning the Bene Gesserit bounded phase space framework."
        />
      </Head>
      <main className="flex w-full flex-col items-center justify-center dark:text-light">
        <Layout className="pt-16">
          <AnimatedText
            text="Research Papers"
            className="!text-6xl mb-16 lg:!text-5xl sm:!text-4xl"
          />
          <motion.div
            className="grid grid-cols-2 gap-8 lg:grid-cols-1"
            variants={container}
            initial="hidden"
            animate="visible"
          >
            {papers.map((paper) => (
              <motion.div key={paper.slug} variants={item}>
                <Link href={`/papers/${paper.slug}`}>
                  <div
                    className="group rounded-2xl border border-solid border-dark/20 dark:border-light/20
                    bg-light dark:bg-dark p-8 hover:border-dark/50 dark:hover:border-light/50
                    transition-all duration-300 cursor-pointer h-full"
                  >
                    <h2
                      className="text-2xl font-bold mb-4 text-dark dark:text-light
                      group-hover:text-primary dark:group-hover:text-primaryDark transition-colors duration-300"
                    >
                      {paper.title}
                    </h2>
                    <p className="text-base text-dark/70 dark:text-light/70 leading-relaxed">
                      {paper.abstract}
                    </p>
                    <span
                      className="inline-block mt-4 text-sm font-semibold text-primary dark:text-primaryDark
                      group-hover:translate-x-2 transition-transform duration-300"
                    >
                      Read paper &rarr;
                    </span>
                  </div>
                </Link>
              </motion.div>
            ))}
          </motion.div>
        </Layout>
      </main>
    </>
  );
}
