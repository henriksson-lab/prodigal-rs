fn main() {
    // All C code has been replaced by Rust.
    // zlib replaced by flate2 crate.
    // Only libm is still linked for math functions used via libc::fprintf etc.
    println!("cargo:rustc-link-lib=m");
}
