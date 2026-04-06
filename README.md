# Prodigal (Rust)

Rust wrapper around [Prodigal](https://github.com/hyattpd/Prodigal) v2.6.3 — a prokaryotic gene prediction tool. The C implementation is compiled via FFI and exposed both as a CLI (`prodigal-rs`) and as a Rust library crate.

## Building

Requires: `gcc`, `zlib` (development headers), and a Rust toolchain.

```bash
# Build the original C binary (needed for tests)
cd Prodigal && make && cd ..

# Build the Rust CLI
cargo build --release
```

## CLI usage

`prodigal-rs` accepts the same flags as the original `prodigal`:

```bash
# Single genome, GenBank output
prodigal-rs -i genome.fasta -o genes.gbk

# Metagenomic mode, GFF output with protein translations
prodigal-rs -i contigs.fasta -p meta -f gff -o genes.gff -a proteins.faa

# Generate and reuse a training file
prodigal-rs -i genome.fasta -t genome.trn
prodigal-rs -i genome.fasta -t genome.trn -o genes.gbk
```

## Library usage

Add to your `Cargo.toml`:

```toml
[dependencies]
prodigal = { path = "." }
```

### Example: run Prodigal on a FASTA file

```rust
use prodigal::run_prodigal;

fn main() {
    let rc = run_prodigal(&[
        "prodigal",
        "-i", "genome.fasta",
        "-o", "output.gff",
        "-f", "gff",
        "-a", "proteins.faa",
        "-q",
    ]);
    assert_eq!(rc, 0, "Prodigal exited with code {rc}");
}
```

### Example: metagenomic mode

```rust
use prodigal::run_prodigal;

fn main() {
    let rc = run_prodigal(&[
        "prodigal",
        "-i", "contigs.fasta",
        "-p", "meta",
        "-f", "gff",
        "-o", "genes.gff",
        "-a", "proteins.faa",
        "-d", "genes.fna",
        "-q",
    ]);
    assert_eq!(rc, 0);
}
```

### Example: using the C FFI directly

```rust
use std::ffi::CString;
use prodigal::ffi::prodigal_c_main;

fn main() {
    let args = ["prodigal", "-i", "input.fasta", "-o", "output.gbk", "-q"];
    let c_strings: Vec<CString> = args.iter()
        .map(|s| CString::new(*s).unwrap())
        .collect();
    let c_ptrs: Vec<*const i8> = c_strings.iter()
        .map(|s| s.as_ptr())
        .collect();

    let rc = unsafe {
        prodigal_c_main(c_ptrs.len() as i32, c_ptrs.as_ptr())
    };
    assert_eq!(rc, 0);
}
```

## Testing

Integration tests verify that the Rust wrapper produces byte-identical output to the original C binary across all output formats and flag combinations:

```bash
cargo test
```

15 tests cover: all 4 output formats (gbk/gff/sco/gca), protein translations (`-a`), nucleotide sequences (`-d`), start scores (`-s`), closed ends (`-c`), masked sequences (`-m`), non-SD mode (`-n`), alternate translation tables (`-g`), metagenomic mode (`-p meta`), and training file round-trips (`-t`).

## License

GPL-3.0 (same as the original Prodigal).
