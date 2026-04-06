//! Integration tests that run both the original C prodigal binary and the Rust
//! wrapper, then assert that all output files are byte-identical.
//!
//! Prerequisites: the C binary must be built at `Prodigal/prodigal` before
//! running these tests (run `make` in the `Prodigal/` directory).

use std::path::{Path, PathBuf};
use std::process::Command;
use tempfile::TempDir;

/// Path to the original C binary (built via `make` in the Prodigal directory).
fn c_binary() -> PathBuf {
    let p = Path::new(env!("CARGO_MANIFEST_DIR")).join("Prodigal/prodigal");
    assert!(
        p.exists(),
        "C prodigal binary not found at {}. Run `make` in Prodigal/ first.",
        p.display()
    );
    p
}

/// Path to the Rust wrapper binary (built by cargo).
fn rust_binary() -> PathBuf {
    let p = Path::new(env!("CARGO_BIN_EXE_prodigal-rs"));
    assert!(p.exists(), "Rust binary not found at {}", p.display());
    p.to_path_buf()
}

/// Path to the sample FASTA input.
fn sample_input() -> PathBuf {
    let p = Path::new(env!("CARGO_MANIFEST_DIR")).join("Prodigal/anthus_aco.fas");
    assert!(p.exists(), "Sample input not found at {}", p.display());
    p
}

/// Run prodigal (C or Rust) with the given extra arguments.
/// Returns the exit status code.
fn run_prodigal(binary: &Path, extra_args: &[&str]) -> i32 {
    let output = Command::new(binary)
        .args(extra_args)
        .output()
        .expect("failed to execute prodigal");

    output.status.code().unwrap_or(-1)
}

/// Assert two files are byte-identical, with a helpful diff message on failure.
fn assert_files_identical(path_a: &Path, path_b: &Path) {
    let a = std::fs::read(path_a).unwrap_or_else(|e| panic!("cannot read {}: {e}", path_a.display()));
    let b = std::fs::read(path_b).unwrap_or_else(|e| panic!("cannot read {}: {e}", path_b.display()));

    if a != b {
        // Show first differing line for debugging
        let lines_a: Vec<&str> = std::str::from_utf8(&a)
            .unwrap_or("<binary>")
            .lines()
            .collect();
        let lines_b: Vec<&str> = std::str::from_utf8(&b)
            .unwrap_or("<binary>")
            .lines()
            .collect();

        for (i, (la, lb)) in lines_a.iter().zip(lines_b.iter()).enumerate() {
            if la != lb {
                panic!(
                    "Files differ at line {} (0-indexed):\n  C:    {}\n  Rust: {}\n  C file:    {}\n  Rust file: {}",
                    i, la, lb, path_a.display(), path_b.display()
                );
            }
        }

        if lines_a.len() != lines_b.len() {
            panic!(
                "Files have different line counts: C={} Rust={}\n  C file:    {}\n  Rust file: {}",
                lines_a.len(),
                lines_b.len(),
                path_a.display(),
                path_b.display()
            );
        }

        panic!(
            "Files differ (binary content) but no text diff found.\n  C file:    {}\n  Rust file: {}",
            path_a.display(),
            path_b.display()
        );
    }
}

/// Helper: run both binaries with the same args and compare all output files.
fn compare_run(label: &str, common_args: &[&str], output_files: &[&str]) {
    let tmp_c = TempDir::new().unwrap();
    let tmp_rs = TempDir::new().unwrap();

    // Build full arg lists, substituting OUTDIR placeholder
    let build_args = |tmp: &TempDir| -> Vec<String> {
        common_args
            .iter()
            .map(|a| a.replace("OUTDIR", tmp.path().to_str().unwrap()))
            .collect()
    };

    let c_args = build_args(&tmp_c);
    let rs_args = build_args(&tmp_rs);

    let c_str_args: Vec<&str> = c_args.iter().map(|s| s.as_str()).collect();
    let rs_str_args: Vec<&str> = rs_args.iter().map(|s| s.as_str()).collect();

    let rc_c = run_prodigal(&c_binary(), &c_str_args);
    let rc_rs = run_prodigal(&rust_binary(), &rs_str_args);

    assert_eq!(
        rc_c, rc_rs,
        "[{label}] Exit codes differ: C={rc_c}, Rust={rc_rs}"
    );

    for filename in output_files {
        let c_file = tmp_c.path().join(filename);
        let rs_file = tmp_rs.path().join(filename);

        assert!(
            c_file.exists(),
            "[{label}] C output file missing: {}",
            c_file.display()
        );
        assert!(
            rs_file.exists(),
            "[{label}] Rust output file missing: {}",
            rs_file.display()
        );

        assert_files_identical(&c_file, &rs_file);
    }

    eprintln!("[{label}] PASSED - all {} output files identical", output_files.len());
}

#[test]
fn test_single_genome_gbk_output() {
    let input = sample_input();
    compare_run(
        "single-gbk",
        &[
            "-i", input.to_str().unwrap(),
            "-o", "OUTDIR/output.gbk",
            "-q",
        ],
        &["output.gbk"],
    );
}

#[test]
fn test_single_genome_gff_output() {
    let input = sample_input();
    compare_run(
        "single-gff",
        &[
            "-i", input.to_str().unwrap(),
            "-o", "OUTDIR/output.gff",
            "-f", "gff",
            "-q",
        ],
        &["output.gff"],
    );
}

#[test]
fn test_single_genome_sco_output() {
    let input = sample_input();
    compare_run(
        "single-sco",
        &[
            "-i", input.to_str().unwrap(),
            "-o", "OUTDIR/output.sco",
            "-f", "sco",
            "-q",
        ],
        &["output.sco"],
    );
}

#[test]
fn test_single_genome_with_translations() {
    let input = sample_input();
    compare_run(
        "single-translations",
        &[
            "-i", input.to_str().unwrap(),
            "-o", "OUTDIR/output.gbk",
            "-a", "OUTDIR/proteins.faa",
            "-q",
        ],
        &["output.gbk", "proteins.faa"],
    );
}

#[test]
fn test_single_genome_with_nucleotides() {
    let input = sample_input();
    compare_run(
        "single-nucleotides",
        &[
            "-i", input.to_str().unwrap(),
            "-o", "OUTDIR/output.gbk",
            "-d", "OUTDIR/genes.fna",
            "-q",
        ],
        &["output.gbk", "genes.fna"],
    );
}

#[test]
fn test_single_genome_with_starts() {
    let input = sample_input();
    compare_run(
        "single-starts",
        &[
            "-i", input.to_str().unwrap(),
            "-o", "OUTDIR/output.gbk",
            "-s", "OUTDIR/starts.txt",
            "-q",
        ],
        &["output.gbk", "starts.txt"],
    );
}

#[test]
fn test_single_genome_all_outputs() {
    let input = sample_input();
    compare_run(
        "single-all",
        &[
            "-i", input.to_str().unwrap(),
            "-o", "OUTDIR/output.gff",
            "-f", "gff",
            "-a", "OUTDIR/proteins.faa",
            "-d", "OUTDIR/genes.fna",
            "-s", "OUTDIR/starts.txt",
            "-q",
        ],
        &["output.gff", "proteins.faa", "genes.fna", "starts.txt"],
    );
}

#[test]
fn test_single_genome_closed_ends() {
    let input = sample_input();
    compare_run(
        "single-closed",
        &[
            "-i", input.to_str().unwrap(),
            "-o", "OUTDIR/output.gbk",
            "-c",
            "-q",
        ],
        &["output.gbk"],
    );
}

#[test]
fn test_single_genome_nonsd() {
    let input = sample_input();
    compare_run(
        "single-nonsd",
        &[
            "-i", input.to_str().unwrap(),
            "-o", "OUTDIR/output.gbk",
            "-n",
            "-q",
        ],
        &["output.gbk"],
    );
}

#[test]
fn test_single_genome_mask() {
    let input = sample_input();
    compare_run(
        "single-mask",
        &[
            "-i", input.to_str().unwrap(),
            "-o", "OUTDIR/output.gbk",
            "-m",
            "-q",
        ],
        &["output.gbk"],
    );
}

#[test]
fn test_single_genome_trans_table_4() {
    let input = sample_input();
    compare_run(
        "single-tt4",
        &[
            "-i", input.to_str().unwrap(),
            "-o", "OUTDIR/output.gbk",
            "-g", "4",
            "-q",
        ],
        &["output.gbk"],
    );
}

#[test]
fn test_metagenomic_gbk() {
    let input = sample_input();
    compare_run(
        "meta-gbk",
        &[
            "-i", input.to_str().unwrap(),
            "-o", "OUTDIR/output.gbk",
            "-p", "meta",
            "-q",
        ],
        &["output.gbk"],
    );
}

#[test]
fn test_metagenomic_gff() {
    let input = sample_input();
    compare_run(
        "meta-gff",
        &[
            "-i", input.to_str().unwrap(),
            "-o", "OUTDIR/output.gff",
            "-f", "gff",
            "-p", "meta",
            "-q",
        ],
        &["output.gff"],
    );
}

#[test]
fn test_metagenomic_all_outputs() {
    let input = sample_input();
    compare_run(
        "meta-all",
        &[
            "-i", input.to_str().unwrap(),
            "-o", "OUTDIR/output.gff",
            "-f", "gff",
            "-a", "OUTDIR/proteins.faa",
            "-d", "OUTDIR/genes.fna",
            "-s", "OUTDIR/starts.txt",
            "-p", "meta",
            "-q",
        ],
        &["output.gff", "proteins.faa", "genes.fna", "starts.txt"],
    );
}

#[test]
fn test_training_file_roundtrip() {
    let input = sample_input();
    let tmp_c = TempDir::new().unwrap();
    let tmp_rs = TempDir::new().unwrap();

    // Step 1: Generate training file with C binary
    let rc = run_prodigal(
        &c_binary(),
        &[
            "-i", input.to_str().unwrap(),
            "-t", tmp_c.path().join("training.trn").to_str().unwrap(),
            "-q",
        ],
    );
    assert_eq!(rc, 0, "C training file generation failed");

    // Step 2: Generate training file with Rust binary
    let rc = run_prodigal(
        &rust_binary(),
        &[
            "-i", input.to_str().unwrap(),
            "-t", tmp_rs.path().join("training.trn").to_str().unwrap(),
            "-q",
        ],
    );
    assert_eq!(rc, 0, "Rust training file generation failed");

    // Training files should be identical
    assert_files_identical(
        &tmp_c.path().join("training.trn"),
        &tmp_rs.path().join("training.trn"),
    );

    // Step 3: Use C training file with Rust binary and vice versa — outputs should match
    let rc_c = run_prodigal(
        &c_binary(),
        &[
            "-i", input.to_str().unwrap(),
            "-t", tmp_c.path().join("training.trn").to_str().unwrap(),
            "-o", tmp_c.path().join("output.gbk").to_str().unwrap(),
            "-q",
        ],
    );
    let rc_rs = run_prodigal(
        &rust_binary(),
        &[
            "-i", input.to_str().unwrap(),
            "-t", tmp_rs.path().join("training.trn").to_str().unwrap(),
            "-o", tmp_rs.path().join("output.gbk").to_str().unwrap(),
            "-q",
        ],
    );
    assert_eq!(rc_c, 0);
    assert_eq!(rc_rs, 0);

    assert_files_identical(
        &tmp_c.path().join("output.gbk"),
        &tmp_rs.path().join("output.gbk"),
    );

    eprintln!("[training-roundtrip] PASSED");
}
