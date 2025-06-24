use std::env;
use std::path::PathBuf;

fn main() {
    // C++ 소스 파일들을 컴파일합니다
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

    // bindgen을 사용하여 Rust 바인딩을 생성합니다
    let bindings = bindgen::Builder::default()
        .header("zbessel.h")
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .generate()
        .expect("바인딩 생성에 실패했습니다");

    // 생성된 바인딩을 OUT_DIR에 저장합니다
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("바인딩을 파일에 쓰는데 실패했습니다!");
}
