import Head from "next/head";
import Image from "next/image";
import { motion } from "framer-motion";

export default function Home() {
  return (
    <>
      <Head>
        <title>Bene Gesserit — Bounded Phase Space Framework</title>
        <meta
          name="description"
          content="Bene Gesserit: A bounded phase space framework for categorical processing on amphiphilic bilayer substrates."
        />
      </Head>
      <div className="relative w-full h-screen overflow-hidden">
        {/* Fullscreen background image */}
        <Image
          src="/bene-gesserit.png"
          alt="Bene Gesserit"
          fill
          className="object-cover"
          priority
        />
        {/* Dark overlay */}
        <div className="absolute inset-0 bg-black/60" />
        {/* Centered text */}
        <motion.div
          className="absolute inset-0 flex flex-col items-center justify-center text-center z-10"
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          transition={{ duration: 2 }}
        >
          <h1 className="text-7xl font-bold text-white tracking-widest mb-4 md:text-5xl sm:text-4xl">
            BENE GESSERIT
          </h1>
          <p className="text-xl text-gray-300 tracking-wide md:text-lg sm:text-base">
            Bounded Phase Space Framework
          </p>
          <p className="text-sm text-gray-400 mt-2 sm:text-xs sm:px-4">
            Categorical Processing on Amphiphilic Bilayer Substrates
          </p>
        </motion.div>
      </div>
    </>
  );
}
