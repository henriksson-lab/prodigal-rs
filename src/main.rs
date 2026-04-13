use clap::Parser;
use std::process;

/// Prokaryotic gene prediction using dynamic programming.
///
/// Prodigal v2.6.3 — fully rewritten in Rust.
#[derive(Parser, Debug)]
#[command(name = "prodigal-rs", version = "2.6.3")]
struct Cli {
    /// Write protein translations to the selected file
    #[arg(short = 'a')]
    trans_file: Option<String>,

    /// Closed ends — do not allow genes to run off edges
    #[arg(short = 'c')]
    closed: bool,

    /// Write nucleotide sequences of genes to the selected file
    #[arg(short = 'd')]
    nuc_file: Option<String>,

    /// Output format: gbk, gff, sco, or gca (default gbk)
    #[arg(short = 'f')]
    output_format: Option<String>,

    /// Translation table (default 11)
    #[arg(short = 'g')]
    trans_table: Option<i32>,

    /// Input FASTA/Genbank file (default stdin)
    #[arg(short = 'i')]
    input_file: Option<String>,

    /// Treat runs of N as masked sequence
    #[arg(short = 'm')]
    mask: bool,

    /// Bypass Shine-Dalgarno trainer and force full motif scan
    #[arg(short = 'n')]
    force_nonsd: bool,

    /// Output file (default stdout)
    #[arg(short = 'o')]
    output_file: Option<String>,

    /// Procedure: single or meta (default single)
    #[arg(short = 'p')]
    mode: Option<String>,

    /// Run quietly (suppress normal stderr output)
    #[arg(short = 'q')]
    quiet: bool,

    /// Write all potential genes (with scores) to the selected file
    #[arg(short = 's')]
    start_file: Option<String>,

    /// Training file (read if exists, write if not)
    #[arg(short = 't')]
    train_file: Option<String>,
}

fn parse_output_format(s: &str) -> i32 {
    match s.to_lowercase().as_str() {
        "0" | "gbk" => 0,
        "1" | "gca" => 1,
        "2" | "sco" => 2,
        "3" | "gff" => 3,
        _ => {
            eprintln!("\nInvalid output format specified.");
            process::exit(15);
        }
    }
}

fn parse_mode(s: &str) -> bool {
    match s.as_bytes().first() {
        Some(b'0') | Some(b's') | Some(b'S') => false,
        Some(b'1') | Some(b'm') | Some(b'M') => true,
        _ => {
            eprintln!("\nInvalid meta/single genome type specified.");
            process::exit(15);
        }
    }
}

fn validate_trans_table(tt: i32) -> i32 {
    if tt < 1 || tt > 25 || tt == 7 || tt == 8 || (tt >= 17 && tt <= 20) {
        eprintln!("\nInvalid translation table specified.");
        process::exit(15);
    }
    tt
}

fn main() {
    let cli = Cli::parse();

    let config = prodigal_rs::pipeline::PipelineConfig {
        input_file: cli.input_file,
        output_file: cli.output_file,
        trans_file: cli.trans_file,
        nuc_file: cli.nuc_file,
        start_file: cli.start_file,
        train_file: cli.train_file,
        output_format: cli.output_format.as_deref().map_or(0, parse_output_format),
        trans_table: cli.trans_table.map_or(0, validate_trans_table),
        closed: cli.closed,
        do_mask: cli.mask,
        force_nonsd: cli.force_nonsd,
        is_meta: cli.mode.as_deref().map_or(false, parse_mode),
        quiet: cli.quiet,
    };

    let rc = unsafe { prodigal_rs::pipeline::run_pipeline(&config) };

    process::exit(rc);
}
