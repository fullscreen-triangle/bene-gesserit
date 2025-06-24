use std::env;

fn main() {
    // Enable documentation features for docs.rs builds
    if env::var("DOCS_RS").is_ok() {
        println!("cargo:rustc-cfg=docsrs");
    }

    // Check for native math library support
    if cfg!(target_os = "linux") || cfg!(target_os = "macos") {
        println!("cargo:rustc-link-lib=m");
    }

    // Enable SIMD optimizations if available
    if is_x86_feature_detected("avx2") {
        println!("cargo:rustc-cfg=has_avx2");
    }
    
    if is_x86_feature_detected("fma") {
        println!("cargo:rustc-cfg=has_fma");
    }

    // Set optimization flags for scientific computing
    println!("cargo:rustc-env=CARGO_CFG_TARGET_FEATURE=+fma,+avx2");
    
    // Rerun if any of these files change
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=Cargo.toml");
    println!("cargo:rerun-if-env-changed=DOCS_RS");
}

#[cfg(target_arch = "x86_64")]
fn is_x86_feature_detected(feature: &str) -> bool {
    match feature {
        "avx2" => std::arch::is_x86_feature_detected!("avx2"),
        "fma" => std::arch::is_x86_feature_detected!("fma"),
        _ => false,
    }
}

#[cfg(not(target_arch = "x86_64"))]
fn is_x86_feature_detected(_feature: &str) -> bool {
    false
} 