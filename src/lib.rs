pub mod types;
pub mod reader;
pub mod bitmap;
pub mod training;
pub mod training_data;
pub mod sequence;
pub mod node;
pub mod dprog;
pub mod gene;
pub mod metagenomic;
pub mod pipeline;

/// Run Prodigal with the given command-line arguments.
///
/// `args` should be the full argv including the program name as the first element.
///
/// # Safety
///
/// The underlying code calls `exit()` on errors and on `-h`/`-v` flags.
pub fn run_prodigal(args: &[&str]) -> i32 {
    let owned: Vec<String> = args.iter().map(|s| s.to_string()).collect();
    unsafe { pipeline::run_pipeline(&owned) }
}
