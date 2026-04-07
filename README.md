# Prodigal (Rust)

Pure Rust rewrite of [Prodigal](https://github.com/hyattpd/Prodigal) v2.6.3 — a prokaryotic gene prediction tool. Produces byte-identical output to the original C implementation, with no C dependencies.

The use of Rust enables the integration of Prodigal as a library in your project, and in applications such as webassembly

## Installation

```bash
cargo install prodigal-rs
```

Or build from source:

```bash
cargo build --release
```

No C compiler, zlib, or other system libraries required.

## Library usage

```toml
[dependencies]
prodigal-rs = "0.1"
```

### Metagenomic mode (simplest)

Predict genes using pre-trained models — no training step needed:

```rust
use prodigal_rs::predict_meta;

fn main() -> Result<(), prodigal_rs::ProdigalError> {
    let sequence = b"ATGCGATCGATCGATCG...";
    let genes = predict_meta(sequence)?;

    for gene in &genes {
        println!(
            "{}-{} ({}) score={:.1} conf={:.1}%",
            gene.begin, gene.end, gene.strand, gene.score, gene.confidence
        );
    }
    Ok(())
}
```

### Single genome mode

Train on the genome first, then predict on individual contigs:

```rust
use prodigal_rs::{train, predict};

fn main() -> Result<(), prodigal_rs::ProdigalError> {
    // Train on the full genome (>= 20,000 bp required)
    let genome = std::fs::read("genome.fasta")?;
    let training = train(&genome)?;

    // Predict genes on each contig using the trained model
    let contig = b"ATGCGATCGATCG...";
    let genes = predict(contig, &training)?;

    for gene in &genes {
        println!(
            "{}-{} {} {} gc={:.3}",
            gene.begin, gene.end, gene.strand, gene.start_codon, gene.gc_content
        );
    }

    // Save/load training for later reuse
    training.save("genome.trn")?;
    let training = prodigal_rs::TrainingData::load("genome.trn")?;

    Ok(())
}
```

### Batch processing (parallel)

For processing many sequences, `MetaPredictor` caches the 50 models and evaluates them in parallel:

```rust
use prodigal_rs::MetaPredictor;

fn main() -> Result<(), prodigal_rs::ProdigalError> {
    let predictor = MetaPredictor::new()?;

    for record in my_fasta_records {
        let genes = predictor.predict(&record.sequence)?;
        println!("{}: {} genes", record.name, genes.len());
    }
    Ok(())
}
```

### Custom configuration

```rust
use prodigal_rs::{predict_meta_with, ProdigalConfig};

fn main() -> Result<(), prodigal_rs::ProdigalError> {
    let config = ProdigalConfig {
        translation_table: 4,    // Mycoplasma genetic code
        closed_ends: true,       // No genes running off edges
        ..Default::default()
    };

    let genes = predict_meta_with(b"ATGCGATCG...", &config)?;
    Ok(())
}
```

### Gene prediction results

Each `PredictedGene` contains:

| Field | Type | Description |
|-------|------|-------------|
| `begin` | `usize` | 1-indexed start position |
| `end` | `usize` | 1-indexed end position (inclusive) |
| `strand` | `Strand` | `Forward` or `Reverse` |
| `start_codon` | `StartCodon` | `ATG`, `GTG`, `TTG`, or `Edge` |
| `partial` | `(bool, bool)` | Left/right partial (runs off edge) |
| `rbs_motif` | `String` | RBS motif (e.g. "AGGAG" or "None") |
| `rbs_spacer` | `String` | RBS spacer distance |
| `gc_content` | `f64` | GC content of the gene |
| `confidence` | `f64` | Confidence score (50-100) |
| `score` | `f64` | Total score |
| `cscore` | `f64` | Coding potential score |
| `sscore` | `f64` | Start score |
| `rscore` | `f64` | RBS score |
| `uscore` | `f64` | Upstream composition score |
| `tscore` | `f64` | Start codon type score |

## CLI usage

`prodigal-rs` accepts the same flags as the original `prodigal`:

```bash
# Single genome, GenBank output
prodigal-rs -i genome.fasta -o genes.gbk

# Metagenomic mode, GFF output with protein translations
prodigal-rs -i contigs.fasta -p meta -f gff -o genes.gff -a proteins.faa

# All outputs at once
prodigal-rs -i genome.fasta -o genes.gff -f gff -a proteins.faa -d genes.fna -s starts.txt

# Generate and reuse a training file
prodigal-rs -i genome.fasta -t genome.trn
prodigal-rs -i genome.fasta -t genome.trn -o genes.gbk
```

Supports gzip-compressed input files transparently.

## Testing

```bash
# Build the C binary first (needed only for comparison tests)
cd Prodigal && make && cd ..

cargo test
```

28 tests cover: the high-level API (metagenomic prediction, single-genome training + prediction, training save/load, error handling, custom config), byte-identical CLI output vs the original C binary across all output formats and flag combinations, and struct layout verification.

## Performance

~2x faster than the original C implementation (gcc -O3) on typical workloads.

## License

GPL-3.0 (same as the original Prodigal).

## Credits

Original Prodigal by Doug Hyatt, University of Tennessee / Oak Ridge National Lab.
Based on Prodigal commit [`c1e2d36`](https://github.com/hyattpd/Prodigal/commit/c1e2d361479cc1b18175ea79ebd8ff10411c46cb) (v2.6.3).
