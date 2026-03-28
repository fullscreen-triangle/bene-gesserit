import Head from "next/head";
import { motion } from "framer-motion";
import Layout from "@/components/Layout";
import AnimatedText from "@/components/AnimatedText";
import Link from "next/link";

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
    title: "Licensing",
    content:
      "The Bene Gesserit framework is available under a dual licensing model. Academic and research use is supported under an open license, while commercial applications require a separate commercial license. This model ensures broad accessibility for the scientific community while sustaining continued development and support for production deployments.",
  },
  {
    title: "Research Partnerships",
    content:
      "We actively seek research partnerships with institutions working at the intersection of membrane biophysics, category theory, and computational architecture. Partners gain early access to framework developments, dedicated support, and opportunities for co-authorship on publications arising from collaborative work.",
  },
  {
    title: "Investment Opportunities",
    content:
      "The Bene Gesserit framework represents a fundamentally new approach to computation grounded in rigorous mathematics and physical science. Investment supports the continued development of the core framework, expansion of the research program, and translation of theoretical advances into practical applications across domains including scientific computing, cognitive systems, and novel hardware architectures.",
  },
  {
    title: "Contact",
    content:
      "For licensing inquiries, partnership proposals, or investment discussions, please reach out through our GitHub repository or contact the development team directly. We welcome conversations with anyone interested in the future of membrane-based categorical computation.",
  },
];

export default function Invest() {
  return (
    <>
      <Head>
        <title>Invest — Bene Gesserit</title>
        <meta
          name="description"
          content="Licensing, partnerships, and investment opportunities for the Bene Gesserit framework."
        />
      </Head>
      <main className="flex w-full flex-col items-center justify-center dark:text-light">
        <Layout className="pt-16">
          <AnimatedText
            text="Invest"
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

            <motion.div variants={item} className="text-center pt-8">
              <Link
                href="https://github.com/fullscreen-triangle/bene-gesserit"
                target="_blank"
                className="inline-block rounded-lg border-2 border-solid border-dark dark:border-light
                bg-dark dark:bg-light text-light dark:text-dark px-8 py-3 text-lg font-semibold
                hover:bg-transparent hover:text-dark dark:hover:bg-transparent dark:hover:text-light
                transition-all duration-300"
              >
                View on GitHub
              </Link>
            </motion.div>
          </motion.div>
        </Layout>
      </main>
    </>
  );
}
