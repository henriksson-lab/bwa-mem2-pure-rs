# bwa-mem2-pure-rs

A faithful Rust translation of `bwa-mem2`

* 2026-08-01: Further unix-isms removed. CI added
* 2026-07-02: Safer API for programs integrating bwa-mem2 algorithm, and bug fix to unsafe memory handling. Benchmark updated
* 2026-07-01: Large improvements in paired-end RSS. The shipped CLI now constrains glibc allocator arenas before startup, and the translated seed storage follows upstream's arena-backed layout more closely
* 2026-05-31: new audit; many edge cases now handled better
* 2026-05-15: the `bwa-mem2-rs` binary can be built with mimalloc for about 10% better wall time
* 2026-04-30: Index generation now as fast as original
* 2026-04-26: Passing all tests so far, speed on par with original. More testing needed though - use on your own risk!

Anyone interested in BWA-MEM2 might also be interested in [BWA-MEM3](https://github.com/fg-labs/bwa-mem3). BWA-MEM3 adds features and has correctness fixes but is in C++; while our BWA-MEM2-RS crate aims to replicate the original BWA-MEM2 C++ codebase behavior faithfully, but in Rust.

## This is an LLM-mediated faithful (hopefully) translation, not the original code! 

Most users should probably first see if the existing original code works for them, unless they have reason otherwise. The original source
may have newer features and it has had more love in terms of fixing bugs. In fact, we aim to replicate bugs if they are present, for the
sake of reproducibility! (but then we might have added a few more in the process)

There are however cases when you might prefer this Rust version. We generally agree with [this manifesto](https://rewrites.bio/) but more specifically:
* We have had many issues with ensuring that our software works using existing containers (Docker, PodMan, Singularity). One size does not fit all and it eats our resources trying to keep up with every way of delivering software
* Common package managers do not work well. It was great when we had a few Linux distributions with stable procedures, but now there are just too many ecosystems (Homebrew, Conda). Conda has an NP-complete resolver which does not scale. Homebrew is only so-stable. And our dependencies in Python still break. These can no longer be considered professional serious options. Meanwhile, Cargo enables multiple versions of packages to be available, even within the same program(!)
* The future is the web. We deploy software in the web browser, and until now that has meant Javascript. This is a language where even the == operator is broken. Typescript is one step up, but a game changer is the ability to compile Rust code into webassembly, enabling performance and sharing of code with the backend. Translating code to Rust enables new ways of deployment and running code in the browser has especial benefits for science - researchers do not have deep pockets to run servers, so pushing compute to the user enables deployment that otherwise would be impossible
* Old CLI-based utilities are bad for the environment(!). A large amount of compute resources are spent creating and communicating via small files, which we can bypass by using code as libraries. Even better, we can avoid frequent reloading of databases by hoisting this stage, with up to 100x speedups in some cases. Less compute means faster compute and less electricity wasted
* LLM-mediated translations may actually be safer to use than the original code. This article shows that [running the same code on different operating systems can give somewhat different answers](https://doi.org/10.1038/nbt.3820). This is a gap that Rust+Cargo can reduce. Typesafe interfaces also reduce coding mistakes and error handling, as opposed to typical command-line scripting

But:

* **This approach should still be considered experimental**. The LLM technology is immature and has sharp corners. But there are opportunities to reap, and the genie is not going back into the bottle. This translation is as much aimed to learn how to improve the technology and get feedback on the results.
* Translations are not endorsed by the original authors unless otherwise noted. **Do not send bug reports to the original developers**. Use our Github issues page instead.
* **Do not trust the benchmarks on this page**. They are used to help evaluate the translation. If you want improved performance, you generally have to use this code as a library, and use the additional tricks it offers. We generally accept performance losses in order to reduce our dependency issues
* **Check the original Github pages for information about the package**. This README is kept sparse on purpose. It is not meant to be the primary source of information
* **If you are the author of the original code and wish to move to Rust, you can obtain ownership of this repository and crate**. Until then, our commitment is to offer an as-faithful-as-possible translation of a snapshot of your code. If we find serious bugs, we will report them to you. Otherwise we will just replicate them, to ensure comparability across studies that claim to use package XYZ v.666. Think of this like a fancy Ubuntu .deb-package of your software - that is how we treat it

This blurb might be out of date. Go to [this page](https://github.com/henriksson-lab/rustification) for the latest information and further information about how we approach translation

## Building the Binary

A release build is required — the workload is CPU-bound and a debug build is roughly an order of magnitude slower.

```sh
cargo build --release --bin bwa-mem2-rs --features mimalloc
```

The binary is written to `target/release/bwa-mem2-rs` and mirrors the upstream subcommands:

```sh
./target/release/bwa-mem2-rs index [-p prefix] <ref.fa>
./target/release/bwa-mem2-rs mem -t <threads> -o out.sam <prefix> <r1.fq> [r2.fq]
```

### Input Pipes and Remote Inputs

This Rust port intentionally differs from upstream `bwa-mem2` for command-style
inputs. Upstream accepts filenames beginning with `<` and runs the rest through
the system shell, for example `'< gzip -dc reads.fq.gz'`. This port does not
support that form, because it depends on Unix shell behavior and is awkward to
support correctly on Windows.

Use a normal shell pipeline into stdin instead:

```sh
gzip -dc reads.fq.gz | ./target/release/bwa-mem2-rs mem -t 4 -o out.sam ref -
```

The `-` input still means stdin, matching upstream behavior. Local `.gz` FASTQ
and FASTA inputs are detected and decompressed directly. Remote `http://` and
`https://` inputs are downloaded with Rust code rather than external `curl` or
`wget`. `ftp://` inputs are not supported.

The `mimalloc` feature is required for the shipped CLI binary and is not enabled by default for library users. For profiling builds that keep release optimizations plus debug symbols, use the `profiling` profile: `cargo build --profile profiling --bin bwa-mem2-rs --features mimalloc` (output at `target/profiling/bwa-mem2-rs`).

## Library Usage

Add the crate to your `Cargo.toml`:

```toml
[dependencies]
bwa-mem2-pure-rs = "0.2"
```

Load an existing `bwa-mem2` index once, then align paired reads in batches:

```rust
use std::path::Path;

use bwa_mem2_pure_rs::mem_api::{MemAligner, MemReadPair};
use bwa_mem2_pure_rs::output::SharedWriterOutput;

fn main() -> Result<(), String> {
    let index_prefix = Path::new("ref/ecoli_rel606");
    let mut aligner = MemAligner::builder(index_prefix)
        .threads(2)
        .build()?;

    let pairs = vec![MemReadPair {
        name: "read-1",
        r1: b"ACGTACGTACGT",
        q1: b"FFFFFFFFFFFF",
        r2: b"TGCATGCATGCA",
        q2: b"FFFFFFFFFFFF",
    }];

    print!("{}", aligner.sam_header()?);
    for sam_record in aligner.align_pairs(&pairs)? {
        print!("{sam_record}");
    }

    // Or capture output instead of writing to process stdout/stderr.
    let captured = SharedWriterOutput::with_stream_labels(Vec::new());
    aligner.write_sam_for_pairs(&pairs, &captured)?;
    let captured_text = String::from_utf8(captured.into_inner().unwrap()).unwrap();
    assert!(captured_text.contains("[stdout] @SQ"));

    Ok(())
}
```

The repository also includes a compile-checked example:

```sh
cargo run --example mem_api -- ref/ecoli_rel606
```

The index files must already exist for `index_prefix`, for example from `bwa-mem2-rs index -p ref/ecoli_rel606 ref/ecoli_rel606.fasta`. For server applications that already use Rayon, pass an existing `Arc<rayon::ThreadPool>` with `.thread_pool(pool)` to share it instead of creating an internal pool. The `output` module also provides `StdioOutput` and `SharedWriterOutput` for stdout/stderr-style library capture.

Library consumers do not get [`mimalloc`](https://docs.rs/mimalloc/latest/mimalloc/) through default features. We recommend registering mimalloc, or another high-performance allocator such as `jemallocator`, as the global allocator in your own binary if alignment throughput matters. The default glibc allocator scales poorly under the per-thread allocation pressure of `mem`; on the 700k-read E. coli fixture, mimalloc reduces wall time by roughly 11% at `-t 1` and 17% at `-t 4`. The shipped `bwa-mem2-rs` binary uses mimalloc when built with `--features mimalloc`.

The published crate only ships the `bwa-mem2-rs` binary plus the library. Additional `src/bin/dump-*.rs` tools in the git repository (e.g. `dump-pe-batch-read`) are development-only and intentionally excluded from the crates.io package — build them from a git checkout if needed.

## Benchmark Setup

Original benchmark baseline: vendored upstream BWA-MEM2 from
`https://github.com/bwa-mem2/bwa-mem2`, commit `97978f950c3a`
(`v2.3-1-g97978f9`).

### 2026-07-14 Overnight Rerun

The current benchmark audit reran the available paired-end mapping cases with
three repeats each, using `scripts/benchmark-median-real-data.sh` and comparing
SAM bodies against the vendored C++ binary. The HOVD index-construction
benchmark could not be rerun because its documented source FASTA
`/husky/henriksson/atrandi/kraken_ref/hovd/rep2/HOVD-geneseqences.fasta` was
not present locally.

```bash
OUT_DIR=/tmp/pres_rustification_bwa_median \
TMPDIR=/tmp/pres_rustification_bwa_tmp \
RUNS=3 \
RUN_CPP=1 \
BENCH_K=20000000 \
OPTIONAL_READ_LIMIT=5000 \
scripts/benchmark-median-real-data.sh
```

| Benchmark | Runs | Rust median | C++ median | Original/Rust speedup | Parity |
|---|---:|---:|---:|---:|---|
| `tutorial_bwa_small/default_pe_t1` | 3 | 10.919 s | 12.624 s | 1.156x | SAM body identical |
| `ecoli_miseq_2x250/default_pe_t1` | 3 | 1.309 s | 2.439 s | 1.863x | SAM body identical |
| `ecoli_miseq_2x250/seed_len_k5_pe` | 3 | 30.906 s | 27.779 s | 0.899x | SAM body identical |

Across the completed paired-end mapping rows, the arithmetic mean
original/Rust wall-time speedup is 1.31x. RSS was not captured separately by
this median script, so the 2026-07-14 roll-up leaves the memory ratio as `NA`.

### Historical Index Benchmark

The historical README index-construction benchmark used the first approximately
50 MB of HOVD FASTA records from
`/husky/henriksson/atrandi/kraken_ref/hovd/rep2/HOVD-geneseqences.fasta`,
written to `.tmp/index_bench_hovd50/hovd_first_records_50mb.fa`. That source
file was missing during the 2026-07-14 audit, so these numbers are not part of
the fresh rerun.

| Command | Wall time | Max RSS |
|---|---:|---:|
| `bwa-mem2-rs index` | 33.90s | 1,215,016 KB |
| upstream `bwa-mem2.avx512bw index` | 33.09s | 1,278,304 KB |
| Rust/C++ ratio | 1.02x | 0.95x |

The generated `.amb`, `.ann`, `.pac`, `.0123`, and `.bwt.2bit.64` files were
md5-identical between Rust and upstream C++ for that historical fixture.

## License

MIT (derived from original code)

## Citation

Vasimuddin Md, Sanchit Misra, Heng Li, Srinivas Aluru. Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. IEEE Parallel and Distributed Processing Symposium (IPDPS), 2019. <https://doi.org/10.1109/IPDPS.2019.00041>
