use std::env;
use std::path::PathBuf;

fn main() {
    // Enable recompilation when files change
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=zbessel.cc");
    println!("cargo:rerun-if-changed=zbessel.h");
    println!("cargo:rerun-if-changed=zbessel.hh");

    // Compile C++ source files
    let mut build = cc::Build::new();
    build
        .cpp(true)
        .file("zbessel.cc")
        .include(".")
        .include("zbessel");

    // Set C++17 standard and compiler-specific flags
    if build.get_compiler().is_like_msvc() {
        // MSVC specific flags
        build.flag("/std:c++17");
        build.flag("/wd4996"); // Disable deprecated function warnings
        build.flag("/wd4244"); // Disable conversion warnings
        build.flag("/wd4267"); // Disable size_t conversion warnings
        build.flag("/wd4305"); // Disable truncation warnings
        // Define __restrict for MSVC compatibility
        build.define("__restrict__", "__restrict");
    } else {
        // GCC/Clang specific flags
        build
            .flag("-std=c++17")
            .flag("-Wno-defaulted-function-deleted")
            .flag("-Wno-duplicate-decl-specifier")
            .flag("-Wno-error")
            .flag("-w");
    }

    build.compile("zbessel");

    // Generate Rust bindings using bindgen
    let bindings = bindgen::Builder::default()
        .header("zbessel.h")
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .generate()
        .expect("Failed to generate bindings");

    // Save the generated bindings to OUT_DIR
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Failed to write bindings to file!");
}
