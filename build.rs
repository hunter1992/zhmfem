fn main() {
    // Setting path to pardiso's .so file
    let lib_path = "/home/zhm/Documents/Scripts/Rust/Learn/Calc/PARDISO/panua_ffi_dynamic/lib/";
    println!("cargo:rustc-link-search=native={}", lib_path);

    // Setting rpath(runtime path)
    println!("cargo:rustc-link-arg=-Wl,-rpath,{}", lib_path);

    // link type
    println!("cargo:rustc-link-lib=dylib=pardiso");
}
