use std::os::raw::{c_char, c_int};

extern "C" {
    /// The original Prodigal main() function, renamed via -Dmain=prodigal_c_main
    /// at compile time. Accepts standard C argc/argv.
    pub fn prodigal_c_main(argc: c_int, argv: *const *const c_char) -> c_int;
}
