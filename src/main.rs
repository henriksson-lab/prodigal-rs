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
    trans_table: Option<String>,

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

fn main() {
    let cli = Cli::parse();

    let mut args: Vec<String> = vec!["prodigal".to_string()];

    if let Some(ref f) = cli.trans_file {
        args.push("-a".into());
        args.push(f.clone());
    }
    if cli.closed {
        args.push("-c".into());
    }
    if let Some(ref f) = cli.nuc_file {
        args.push("-d".into());
        args.push(f.clone());
    }
    if let Some(ref f) = cli.output_format {
        args.push("-f".into());
        args.push(f.clone());
    }
    if let Some(ref g) = cli.trans_table {
        args.push("-g".into());
        args.push(g.clone());
    }
    if let Some(ref f) = cli.input_file {
        args.push("-i".into());
        args.push(f.clone());
    }
    if cli.mask {
        args.push("-m".into());
    }
    if cli.force_nonsd {
        args.push("-n".into());
    }
    if let Some(ref f) = cli.output_file {
        args.push("-o".into());
        args.push(f.clone());
    }
    if let Some(ref p) = cli.mode {
        args.push("-p".into());
        args.push(p.clone());
    }
    if cli.quiet {
        args.push("-q".into());
    }
    if let Some(ref f) = cli.start_file {
        args.push("-s".into());
        args.push(f.clone());
    }
    if let Some(ref f) = cli.train_file {
        args.push("-t".into());
        args.push(f.clone());
    }

    let rc = unsafe { prodigal_rs::pipeline::run_pipeline(&args) };

    process::exit(rc);
}
