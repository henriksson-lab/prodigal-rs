fn main() {
    // All C code has been replaced by Rust.
    // Only zlib is still needed for gzip-compressed input support.
    println!("cargo:rustc-link-lib=z");
    println!("cargo:rustc-link-lib=m");
}
