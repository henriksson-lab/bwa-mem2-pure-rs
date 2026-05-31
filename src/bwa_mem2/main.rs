#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Port of `bwa-mem2/src/main.h` + `bwa-mem2/src/main.cpp`.

use crate::bwa_mem2::bwa::BWA_PG;
use crate::bwa_mem2::bwtindex::bwa_index;
use crate::bwa_mem2::fastmap::main_mem;

// --- main.cpp ---

const PACKAGE_VERSION: &str = "2.2.1";

fn usage_text() -> &'static str {
    "Usage: bwa-mem2 <command> <arguments>\nCommands:\n  index         create index\n  mem           alignment\n  version       print version number\n"
}

/// Print the top-level `bwa-mem2 <command> <arguments>` usage banner to
/// stderr and return `1`. Matches `usage` in `bwa-mem2/src/main.cpp`.
#[doc = "Original function: usage:43"]
pub fn usage() -> i32 {
    eprint!("{}", usage_text());
    1
}

/// Top-level CLI dispatcher mirroring `main` in `bwa-mem2/src/main.cpp`.
///
/// The C++ version additionally:
///   * calibrates `proc_freq` by reading `__rdtsc()` around a 1-second sleep
///     and uses it for runtime profiling (`tprof[MEM][0]`),
///   * prints the active SIMD ISA (`AVX512`/`AVX2`/`AVX`/`SSE4`) at startup,
///   * prints `BATCH_SIZE`, `MAX_SEQ_LEN_*`, `SEEDS_PER_READ`, `SIMD_WIDTH*`
///     after a successful run.
/// None of these are required for SAM-output parity, so the Rust port omits
/// them.
#[doc = "Original function: main:53"]
pub fn main(argv: &[String]) -> i32 {
    if argv.len() < 2 {
        return usage();
    }

    match argv[1].as_str() {
        "index" => bwa_index(&argv[1..]),
        "mem" => {
            // Build the SAM @PG line from the command line so downstream
            // record-emission paths can stamp the SAM header with the exact
            // invocation (matches the C++ `ksprintf(&pg, "@PG\tID:bwa-mem2..."`
            // sequence).
            let cmdline = argv.join(" ");
            *BWA_PG.lock().expect("bwa_pg lock") = Some(format!(
                "@PG\tID:bwa-mem2\tPN:bwa-mem2\tVN:{}\tCL:{}\n",
                PACKAGE_VERSION, cmdline
            ));
            let ret = main_mem(&argv[1..]);
            *BWA_PG.lock().expect("bwa_pg lock") = None;
            ret
        }
        "version" => {
            println!("{PACKAGE_VERSION}");
            0
        }
        other => {
            eprintln!("ERROR: unknown command '{other}'");
            1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{main, usage, usage_text, PACKAGE_VERSION};
    use std::fs;
    use std::io::Cursor;
    use std::time::{SystemTime, UNIX_EPOCH};

    use crate::bwa_mem2::bntseq::bns_fasta2bntseq;
    use crate::bwa_mem2::fmi_search::FMI_search;

    #[test]
    fn usage_text_matches_cpp_commands() {
        let text = usage_text();
        assert!(text.starts_with("Usage: bwa-mem2 <command> <arguments>\n"));
        assert!(text.contains("  index         create index\n"));
        assert!(text.contains("  mem           alignment\n"));
        assert!(text.contains("  version       print version number\n"));
        assert_eq!(usage(), 1);
    }

    #[test]
    fn main_rejects_unknown_command() {
        let argv = vec!["bwa-mem2".to_string(), "wat".to_string()];
        assert_eq!(main(&argv), 1);
    }

    #[test]
    fn main_dispatches_index() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_main_index_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        let fasta = dir.join("ref.fa");
        let prefix = dir.join("idx");
        fs::write(&fasta, b">chr1\nACGT\n").expect("write fasta");

        let argv = vec![
            "bwa-mem2".to_string(),
            "index".to_string(),
            "-p".to_string(),
            prefix.to_str().expect("utf8").to_string(),
            fasta.to_str().expect("utf8").to_string(),
        ];
        assert_eq!(main(&argv), 0);
        assert!(std::path::Path::new(&format!("{}.pac", prefix.display())).is_file());
        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn main_dispatches_mem() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_main_mem_dispatch_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("ref");
        let fasta = b">chr1\nACGTACGTACGT\n";
        bns_fasta2bntseq(
            Cursor::new(fasta.as_slice()),
            prefix.to_str().expect("utf8"),
            1,
        );

        let mut fmi = FMI_search::ctor(prefix.to_str().expect("utf8"));
        assert_eq!(fmi.build_index(), 0);
        fmi.load_index();
        fmi.dtor();

        let reads = dir.join("reads.fq");
        fs::write(&reads, b"@mm1\nACGTAC\n+\nIIIIII\n").expect("write reads");
        let out = dir.join("out.sam");
        let argv = vec![
            "bwa-mem2".to_string(),
            "mem".to_string(),
            "-o".to_string(),
            out.to_str().expect("utf8").to_string(),
            "-t".to_string(),
            "1".to_string(),
            "-k".to_string(),
            "2".to_string(),
            "-T".to_string(),
            "1".to_string(),
            prefix.to_str().expect("utf8").to_string(),
            reads.to_str().expect("utf8").to_string(),
        ];
        let rc = {
            let _guard = crate::bwa_mem2::fastmap::MAIN_MEM_TEST_LOCK
                .lock()
                .expect("main_mem test lock");
            main(&argv)
        };
        assert_eq!(rc, 0);
        let sam = fs::read_to_string(&out).expect("sam output");
        assert!(sam.contains("@PG\tID:bwa-mem2\tPN:bwa-mem2\tVN:"), "{sam}");
        assert!(sam.contains("mm1\t"), "{sam}");
        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn version_constant_is_expected() {
        assert_eq!(PACKAGE_VERSION, "2.2.1");
    }
}
