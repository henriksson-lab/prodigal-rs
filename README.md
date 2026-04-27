# Prodigal-rs (Rust)

Pure Rust translation of [Prodigal](https://github.com/hyattpd/Prodigal) — a prokaryotic gene prediction tool.

Based on Prodigal commit [`c1e2d36`](https://github.com/hyattpd/Prodigal/commit/c1e2d361479cc1b18175ea79ebd8ff10411c46cb) (v2.6.3).

* 2026-04-23 - Code is ready for testing on real data, but stay vigilant for remaining errors. Compare with regular Prodigal on your data before considering this as a replacement. 
* This translation is only marginally faster than the original

## This is an LLM-mediated faithful (hopefully) translation, not the original code!

Most users should probably first see if the existing original code works for them, unless they have reason otherwise. The original source
may have newer features and it has had more love in terms of fixing bugs. In fact, we aim to replicate bugs if they are present, for the
sake of reproducibility! (but then we might have added a few more in the process)

There are however cases when you might prefer this Rust version. We generally agree with [this page](https://rewrites.bio/)
but more specifically:
* We have had many issues with ensuring that our software works using existing containers (Docker, PodMan, Singularity). One size does not fit all and it eats our resources trying to keep up with every way of delivering software
* Common package managers do not work well. It was great when we had a few Linux distributions with stable procedures, but now there are just too many ecosystems (Homebrew, Conda). Conda has an NP-complete resolver which does not scale. Homebrew is only so-stable. And our dependencies in Python still break. These can no longer be considered professional serious options. Meanwhile, Cargo enables multiple versions of packages to be available, even within the same program(!)
* The future is the web. We deploy software in the web browser, and until now that has meant Javascript. This is a language where even the == operator is broken. Typescript is one step up, but a game changer is the ability to compile Rust code into webassembly, enabling performance and sharing of code with the backend. Translating code to Rust enables new ways of deployment and running code in the browser has especial benefits for science - researchers do not have deep pockets to run servers, so pushing compute to the user enables deployment that otherwise would be impossible
* Old CLI-based utilities are bad for the environment(!). A large amount of compute resources are spent creating and communicating via small files, which we can bypass by using code as libraries. Even better, we can avoid frequent reloading of databases by hoisting this stage, with up to 100x speedups in some cases. Less compute means faster compute and less electricity wasted
* LLM-mediated translations may actually be safer to use than the original code. This article shows that [running the same code on different operating systems can give somewhat different answers](https://doi.org/10.1038/nbt.3820). This is a gap that Rust+Cargo can reduce. Typesafe interfaces also reduce coding mistakes and error handling, as opposed to typical command-line scripting

But:

* **This approach should still be considered experimental**. The LLM technology is immature and has sharp corners. But there are opportunities to reap, and the genie is not going back to the bottle. This translation is as much aimed to learn how to improve the technology and get feedback on the results.
* Translations are not endorsed by the original authors unless otherwise noted. **Do not send bug reports to the original developers**. Use our Github issues page instead.
* **Do not trust the benchmarks on this page**. They are used to help evaluate the translation. If you want improved performance, you generally have to use this code as a library, and use the additional tricks it offers. We generally accept performance losses in order to reduce our dependency issues
* **Check the original Github pages for information about the package**. This README is kept sparse on purpose. It is not meant to be the primary source of information

## Translation audit status

The original C source is kept in `Prodigal/`. The root `ccc_mapping.toml` maps original C functions to their Rust counterparts for use with [`code-complexity-comparator`](../code-complexity-comparator). Most translated functions intentionally keep the original function name; the mapping adds path pins and maps C `main` to Rust `run_pipeline`.

Known audit notes:

* The C metagenome initializers `initialize_metagenome_0` through `initialize_metagenome_49` are represented in Rust as binary fixtures loaded by `load_metagenome`, so they are not one-function-per-function mappings in `ccc_mapping.toml`.
* `dprog` omits one defensive guard from the original simple-overlap untangling pass. The C code skips when no matching stop node is found; the Rust loop can walk before the node array in that unexpected case.
* CLI compatibility is not exact: uppercase aliases such as `-A`, `-C`, etc. and lowercase `-v` are accepted by the original but rejected by the current Clap-based CLI.

The last local audit run used `ccc-rs analyze/compare/missing/constants-diff` against `Prodigal/` and `src/`, plus `cargo test`. Comparator output still contains expected noise from Rust formatting syntax and byte literals, so direct source inspection is still required for flagged functions.



## Installation

```bash
cargo install prodigal-rs
```

Or build from source:

```bash
cargo build --release
```

## Library usage

```toml
[dependencies]
prodigal-rs = "0.3"
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

To reuse an application-level Rayon pool, pass an `Arc<rayon::ThreadPool>` to the predictor. The pool should be built with the recommended worker stack size:

```rust
use std::sync::Arc;
use prodigal_rs::{MetaPredictor, META_PREDICTOR_STACK_SIZE};

fn main() -> Result<(), prodigal_rs::ProdigalError> {
    let pool = Arc::new(
        rayon::ThreadPoolBuilder::new()
            .num_threads(8)
            .stack_size(META_PREDICTOR_STACK_SIZE)
            .build()
            .expect("failed to build Rayon pool"),
    );

    let predictor = MetaPredictor::with_thread_pool(pool)?;
    let genes = predictor.predict(b"ATGCGATCG...")?;
    println!("{} genes", genes.len());
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

`prodigal-rs` accepts the main lowercase flags from the original `prodigal`:

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

Supports gzip-compressed input files transparently. Exact CLI compatibility is still pending; see the audit notes above.

## Testing

```bash
# Build the C binary first (needed only for comparison tests)
cd Prodigal && make && cd ..

cargo test
```

49 default tests currently cover: the high-level API (metagenomic prediction, single-genome training + prediction, training save/load, error handling, custom config), byte-identical CLI output vs the original C binary across selected output formats and flag combinations, and struct layout verification.

The comparison test target also includes larger real-data parity checks marked `#[ignore]` so the default suite stays fast. These run full-size sibling-repo fixtures without truncation and compare C and Rust outputs byte-for-byte:

```bash
cargo test --test compare_outputs -- --ignored full_single
cargo test --test compare_outputs -- --ignored full_meta
```

The latest local release check also ran:

```bash
cargo clippy --all-targets -- -D warnings
cargo package --allow-dirty
```

## Performance

Local release-build comparison against the bundled original C `Prodigal/prodigal`, using `/usr/bin/time` on real FASTA inputs already present in the workspace. These are evaluation benchmarks only; they are not meant as a general performance claim across environments.

| Dataset | Mode | Input | C wall | Rust wall | Rust/C | C max RSS | Rust max RSS | Output |
|-------|------|------:|------:|---------:|------:|----------:|-------------:|--------|
| `prokka genome` | `single` | 6.7 MB | 15.69 s | 13.03 s | 0.83x | 172344 KB | 281600 KB | identical GFF |
| `prokka genome` | `meta` | 6.7 MB | 74.82 s | 62.35 s | 0.83x | 127708 KB | 233600 KB | identical GFF |
| `priestia` | `single` | 5.4 MB | 8.56 s | 6.87 s | 0.80x | 120136 KB | 234560 KB | identical GFF |
| `priestia` | `meta` | 5.4 MB | 30.95 s | 27.04 s | 0.87x | 96968 KB | 206080 KB | identical GFF |

Current memory note: the Rust port is faster on these measured runs, but it currently uses about 1.6-2.1x the peak RSS. A likely major reason is that the Rust node buffer grows with `Vec::resize(..., zeroed())`, which materializes the entire new tail in memory, whereas the original C code uses `realloc` and only touches the subset of nodes that are actually populated.

## License

GPL-3.0 (same as the original Prodigal).

## Citing

Hyatt, D., Chen, GL., LoCascio, P.F. et al. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics 11, 119 (2010). doi:10.1186/1471-2105-11-119.
