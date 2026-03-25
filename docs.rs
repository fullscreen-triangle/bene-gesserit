//! # Bene Gesserit Documentation Build Configuration
//! 
//! This file configures documentation generation for docs.rs

#![cfg_attr(docsrs, feature(doc_cfg))]

// Enable all features for documentation builds on docs.rs
#[cfg(docsrs)]
compile_error!("This file is only for documentation configuration") 