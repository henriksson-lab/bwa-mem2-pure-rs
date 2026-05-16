# bwa-mem2-pure-rs

A faithful Rust translation of `bwa-mem2`

* 2026-05-15: the `bwa-mem2-rs` binary can be built with mimalloc for about 10% better wall time
* 2026-04-30: Index generation now as fast as original
* 2026-04-26: Passing all tests so far, speed on par with original. More testing needed though - use on your own risk!

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
        name: "read-1".to_string(),
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

The index files must already exist for `index_prefix`, for example from `bwa-mem2-rs index -p ref/ecoli_rel606 ref/ecoli_rel606.fasta`. For server applications that already use Rayon, pass an existing `Arc<rayon::ThreadPool>` with `.thread_pool(pool)` to share it instead of creating an internal pool. The `output` module also provides `StdioOutput` and `SharedWriterOutput` for stdout/stderr-style library capture.

Library consumers do not get [`mimalloc`](https://docs.rs/mimalloc/latest/mimalloc/) through default features. We recommend registering mimalloc, or another high-performance allocator such as `jemallocator`, as the global allocator in your own binary if alignment throughput matters. The default glibc allocator scales poorly under the per-thread allocation pressure of `mem`; on the 700k-read E. coli fixture, mimalloc reduces wall time by roughly 11% at `-t 1` and 17% at `-t 4`. The shipped `bwa-mem2-rs` binary uses mimalloc when built with `--features mimalloc`.

The published crate only ships the `bwa-mem2-rs` binary plus the library. Additional `src/bin/dump-*.rs` tools in the git repository (e.g. `dump-pe-batch-read`) are development-only and intentionally excluded from the crates.io package — build them from a git checkout if needed.

## Benchmark Setup

These numbers were measured on the real paired-end tutorial fixture in `external/tutorial_bwa_small`:

- Reads:
  - `external/tutorial_bwa_small/SRR2584863_1.trim.sub.fastq`
  - `external/tutorial_bwa_small/SRR2584863_2.trim.sub.fastq`
- Reference:
  - `external/tutorial_bwa_small/ecoli_rel606.fasta.gz`
- Indexed reference prefix:
  - `.tmp/real_bench/ecoli_rel606`
- Command shape:

```bash
target/release/bwa-mem2-rs mem -t 2 .tmp/real_bench/ecoli_rel606 <reads_1.fq> <reads_2.fq>
```

The full-fixture comparisons below were validated by exact SAM-body comparison against upstream output.

Additional real-data conformance coverage can be prepared with:

```bash
scripts/download-real-world-fixtures.py --extra miseq-2x250 --extra-read-limit 5000
STRICT_KNOWN_GAPS=1 READ_LIMIT=2000 scripts/conformance-real-data-matrix.sh
```

The optional `ecoli_miseq_2x250` fixture streams a 5,000-pair subset from ENA run `SRR13321180` (Illumina MiSeq 2x250 E. coli) and compares Rust against upstream C++ using the same REL606 reference index. The matrix also keeps a single-end pass over the tutorial fixture to exercise SE behavior.

The focused SIMD16 regression guard for the non-default scoring case is:

```bash
READ_LIMIT=2000 scripts/regression-simd16-a5.sh
```

## Full-Dataset Speed Comparison

Full local real paired-end fixture, 175,000 paired reads, `-K 20000000`.

| Command | `bwa-mem2-rs` | upstream `bwa-mem2.avx512bw` | Rust vs upstream |
|---|---:|---:|---:|
| `mem -t 1` | 15.49s | 18.50s | 1.19x faster, about 16% less wall time |
| `mem -t 2` | 10.51s | 10.35s | 1.02x slower, about 2% more wall time |

Interpretation:

- These timings exclude index construction.
- Rust and upstream SAM bodies are byte-identical for both `-t 1` and `-t 2`, ignoring header command-line differences.
- Older benchmark notes used a smaller/different setup and should not be compared directly to this full-fixture result.
- Re-run locally with `scripts/benchmark-real-data.sh`; it emits `.tmp/real_bench/report.tsv` when the upstream C++ binary is available.

## Index Generation Speed Comparison

Index construction benchmark on the first approximately 50 MB of HOVD FASTA records from `/husky/henriksson/atrandi/kraken_ref/hovd/rep2/HOVD-geneseqences.fasta`, written to `.tmp/index_bench_hovd50/hovd_first_records_50mb.fa`.

| Command | Wall time | Max RSS |
|---|---:|---:|
| `bwa-mem2-rs index` | 33.90s | 1,215,016 KB |
| upstream `bwa-mem2.avx512bw index` | 33.09s | 1,278,304 KB |

The generated `.amb`, `.ann`, `.pac`, `.0123`, and `.bwt.2bit.64` files were md5-identical between Rust and upstream C++ for this fixture.

## License

MIT (derived from original code)

## Citation

Vasimuddin Md, Sanchit Misra, Heng Li, Srinivas Aluru. Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. IEEE Parallel and Distributed Processing Symposium (IPDPS), 2019. <https://doi.org/10.1109/IPDPS.2019.00041>
