pub mod ffi;

use std::ffi::CString;
use std::os::raw::c_char;

/// Run Prodigal with the given command-line arguments.
///
/// `args` should be the full argv including the program name as the first element.
///
/// # Panics
///
/// Panics if any argument contains an interior NUL byte.
///
/// # Safety
///
/// The underlying C code calls `exit()` on errors and on `-h`/`-v` flags.
/// Use [`run_prodigal_safe`] to catch those exits in a subprocess if needed.
pub fn run_prodigal(args: &[&str]) -> i32 {
    let c_strings: Vec<CString> = args
        .iter()
        .map(|s| CString::new(*s).expect("argument contains NUL byte"))
        .collect();

    let c_ptrs: Vec<*const c_char> = c_strings.iter().map(|s| s.as_ptr()).collect();

    unsafe { ffi::prodigal_c_main(c_ptrs.len() as i32, c_ptrs.as_ptr()) }
}
