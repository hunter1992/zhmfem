fn main() {
    // Setting path to pardiso's .so file
    let lib_path = "/opt/panua-pardiso-20240229-linux/lib/";
    println!("cargo:rustc-link-search=native={}", lib_path);

    // Setting rpath(runtime path)
    println!("cargo:rustc-link-arg=-Wl,-rpath,{}", lib_path);

    // link type
    println!("cargo:rustc-link-lib=dylib=pardiso");
}
