use std::env;
use std::path::PathBuf;

fn main() {
    // Compile C++ source files
    cc::Build::new()
        .cpp(true)
        .file("zbessel.cc")
        .include(".")
        .include("zbessel")
        .flag("-std=c++17")
        .flag("-Wno-defaulted-function-deleted")
        .flag("-Wno-duplicate-decl-specifier")
        .flag("-Wno-error")
        .flag("-w")
        .compile("zbessel");

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
