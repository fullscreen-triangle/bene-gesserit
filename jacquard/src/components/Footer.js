import Link from "next/link";
import React from "react";
import Layout from "./Layout";

const Footer = () => {
  return (
    <footer
      className="w-full border-t-2 border-solid border-dark
      font-medium text-lg dark:text-light dark:border-light sm:text-base"
    >
      <Layout className="py-8 flex items-center justify-between lg:flex-col lg:py-6">
        <span>{new Date().getFullYear()} &copy; Bene Gesserit Framework</span>

        <Link
          href="https://github.com/fullscreen-triangle/bene-gesserit"
          target="_blank"
          className="underline underline-offset-2"
        >
          GitHub
        </Link>
      </Layout>
    </footer>
  );
};

export default Footer;
