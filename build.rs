use std::env;
use std::fs;
use std::path::PathBuf;

fn main() {
    println!("cargo:rerun-if-env-changed=PHASE_MNV_REQUIRE_BINDGEN");

    let target_arch = env::var("CARGO_CFG_TARGET_ARCH").unwrap_or_default();
    let fermi_supported_arch = matches!(target_arch.as_str(), "x86_64");
    if !fermi_supported_arch {
        println!(
            "cargo:warning=fermi-lite FFI disabled on target_arch={target_arch}; vendored ksw.c requires x86_64/SSE2"
        );
        return;
    }

    let fermi_sources = [
        "vendor/fermi-lite/bfc.c",
        "vendor/fermi-lite/bubble.c",
        "src/fermi_lite_shim.c",
        "vendor/fermi-lite/htab.c",
        "vendor/fermi-lite/ksw.c",
        "vendor/fermi-lite/kthread.c",
        "vendor/fermi-lite/mag.c",
        "vendor/fermi-lite/misc.c",
        "vendor/fermi-lite/mrope.c",
        "vendor/fermi-lite/rld0.c",
        "vendor/fermi-lite/rle.c",
        "vendor/fermi-lite/rope.c",
        "vendor/fermi-lite/unitig.c",
    ];

    for path in &fermi_sources {
        println!("cargo:rerun-if-changed={path}");
    }
    println!("cargo:rerun-if-changed=vendor/fermi-lite/fml.h");
    println!("cargo:rerun-if-changed=vendor/fermi-lite/internal.h");
    println!("cargo:rerun-if-changed=src/fermi_lite_shim.c");
    println!("cargo:rerun-if-changed=src/fermi_lite_bindings.rs");

    let out_path = PathBuf::from(env::var("OUT_DIR").expect("OUT_DIR is set by cargo"));
    let bindings_path = out_path.join("fermi_lite_bindings.rs");
    let old_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let generated = std::panic::catch_unwind(|| {
        bindgen::Builder::default()
            .header("vendor/fermi-lite/fml.h")
            .allowlist_type("bseq1_t")
            .allowlist_type("magopt_t")
            .allowlist_type("fml_opt_t")
            .allowlist_type("fml_ovlp_t")
            .allowlist_type("fml_utg_t")
            .allowlist_function("fml_opt_init")
            .allowlist_function("fml_assemble")
            .allowlist_function("fml_utg_destroy")
            .allowlist_var("fm_verbose")
            .generate()
    })
    .ok()
    .and_then(Result::ok);
    std::panic::set_hook(old_hook);
    if let Some(bindings) = generated {
        bindings
            .write_to_file(&bindings_path)
            .expect("failed to write fermi-lite bindings");
    } else if env::var_os("PHASE_MNV_REQUIRE_BINDGEN").is_some() {
        panic!("bindgen failed for vendor/fermi-lite/fml.h and PHASE_MNV_REQUIRE_BINDGEN is set");
    } else {
        println!(
            "cargo:warning=bindgen unavailable or failed; using checked-in fermi-lite fallback bindings"
        );
        fs::copy("src/fermi_lite_bindings.rs", &bindings_path)
            .expect("failed to copy fallback fermi-lite bindings");
    }

    let mut build = cc::Build::new();
    build
        .include("vendor/fermi-lite")
        .define("_GNU_SOURCE", None)
        .flag_if_supported("-Wno-unused-function")
        .flag_if_supported("-Wno-sign-compare")
        .flag_if_supported("-Wno-missing-field-initializers")
        .flag_if_supported("-Wno-unused-parameter")
        .flag_if_supported("-Wno-unused-but-set-variable")
        .flag_if_supported("-Wno-unused-result");
    for path in &fermi_sources {
        build.file(path);
    }
    build.compile("fermi_lite");

    if std::env::var("CARGO_CFG_TARGET_OS").as_deref() != Ok("windows") {
        println!("cargo:rustc-link-lib=m");
        println!("cargo:rustc-link-lib=pthread");
    }
}
