import Head from "next/head";
import { motion } from "framer-motion";
import Layout from "@/components/Layout";
import Link from "next/link";

const sectionVariant = {
  hidden: { opacity: 0, y: 20 },
  visible: (i) => ({
    opacity: 1,
    y: 0,
    transition: { delay: 0.2 + i * 0.1, duration: 0.5 },
  }),
};

const PaperPage = ({ title, abstract, sections = [], charts = [] }) => {
  return (
    <>
      <Head>
        <title>{title} — Bene Gesserit</title>
        <meta name="description" content={abstract} />
      </Head>
      <main className="flex w-full flex-col items-center justify-center dark:text-light">
        <Layout className="pt-16">
          {/* Back link */}
          <Link
            href="/papers"
            className="text-sm text-dark/60 dark:text-light/60 hover:text-primary dark:hover:text-primaryDark
            transition-colors duration-300 mb-8 inline-block"
          >
            &larr; Back to Papers
          </Link>

          {/* Title */}
          <motion.h1
            className="text-5xl font-bold mb-8 text-dark dark:text-light lg:text-4xl md:text-3xl"
            initial={{ opacity: 0, y: -20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.6 }}
          >
            {title}
          </motion.h1>

          {/* Abstract */}
          {abstract && (
            <motion.div
              className="rounded-xl border border-solid border-dark/20 dark:border-light/20
              bg-light dark:bg-dark/50 p-8 mb-12"
              initial={{ opacity: 0, y: 20 }}
              animate={{ opacity: 1, y: 0 }}
              transition={{ delay: 0.2, duration: 0.5 }}
            >
              <h2 className="text-lg font-semibold uppercase tracking-wider text-primary dark:text-primaryDark mb-4">
                Abstract
              </h2>
              <p className="text-base leading-relaxed text-dark/80 dark:text-light/80">
                {abstract}
              </p>
            </motion.div>
          )}

          {/* Sections */}
          {sections.map((section, index) => (
            <motion.div
              key={index}
              className="mb-10"
              custom={index}
              variants={sectionVariant}
              initial="hidden"
              animate="visible"
            >
              <h2 className="text-3xl font-bold mb-4 text-dark dark:text-light lg:text-2xl">
                {section.title}
              </h2>
              <div className="prose prose-lg dark:prose-invert max-w-none text-dark/80 dark:text-light/80 leading-relaxed">
                {section.content}
              </div>

              {/* Render chart if one exists at this index */}
              {charts[index] && (
                <div className="my-8 rounded-xl border border-solid border-dark/10 dark:border-light/10 p-6 bg-light dark:bg-dark/30">
                  {charts[index]}
                </div>
              )}
            </motion.div>
          ))}

          {/* Render any remaining charts not matched to sections */}
          {charts.length > sections.length && (
            <div className="mt-8 space-y-8">
              {charts.slice(sections.length).map((chart, index) => (
                <div
                  key={`extra-chart-${index}`}
                  className="rounded-xl border border-solid border-dark/10 dark:border-light/10 p-6 bg-light dark:bg-dark/30"
                >
                  {chart}
                </div>
              ))}
            </div>
          )}
        </Layout>
      </main>
    </>
  );
};

export default PaperPage;
