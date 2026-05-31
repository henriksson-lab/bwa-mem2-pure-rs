#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Port of `bwa-mem2/src/bwtindex.cpp`.

use crate::bwa_mem2::bntseq::{bns_fasta2bntseq, gzip_or_plain_reader_from_read};
use crate::bwa_mem2::fmi_search::FMI_search;
use std::fs::File;
use std::io::{BufReader, Read};
use std::time::Instant;

// --- bwtindex.cpp ---

/// The "index" subcommand entry point (parses `-p prefix` then calls
/// `bwa_idx_build`). Mirrors `bwa_index` in `bwa-mem2/src/bwtindex.cpp`.
#[doc = "Original function: bwa_index:43"]
pub fn bwa_index(argv: &[String]) -> i32 {
    let mut prefix: Option<String> = None;
    let mut fa: Option<String> = None;
    let mut i = 1;
    let mut end_options = false;
    while i < argv.len() {
        if !end_options && argv[i] == "--" {
            end_options = true;
        } else if !end_options && argv[i].starts_with("-p") && argv[i] != "-p" {
            prefix = Some(argv[i][2..].to_string());
        } else if !end_options && argv[i] == "-p" {
            i += 1;
            if i >= argv.len() {
                return 1;
            }
            prefix = Some(argv[i].clone());
        } else if !end_options && argv[i].starts_with('-') && argv[i] != "-" {
            return 1;
        } else if fa.is_none() {
            fa = Some(argv[i].clone());
        } else {
            // Upstream only requires at least one FASTA operand; extra positional
            // arguments are left unused after getopt has found options.
        }
        i += 1;
    }

    let Some(fa) = fa else {
        eprintln!("Usage: bwa-mem2 index [-p prefix] <in.fasta>");
        return 1;
    };

    let prefix = prefix.unwrap_or_else(|| fa.clone());
    bwa_idx_build(&fa, &prefix)
}

/// Pack the FASTA at `fa` into the `<prefix>.pac`/`.ann`/`.amb` triplet and
/// then build the FM-index over it (writing `<prefix>.bwt.2bit.64` and the
/// other index files). Mirrors `bwa_idx_build` in
/// `bwa-mem2/src/bwtindex.cpp`.
///
/// # Arguments
/// * `fa` - path to the input FASTA (gzip auto-detected by stream magic).
/// * `prefix` - output index prefix.
#[doc = "Original function: bwa_idx_build:61"]
pub fn bwa_idx_build(fa: &str, prefix: &str) -> i32 {
    // Nucleotide indexing: pack the FASTA into a .pac/.ann/.amb triplet, then
    // build the FM-index over it.
    let input: Box<dyn Read> = if fa == "-" {
        Box::new(std::io::stdin())
    } else {
        Box::new(File::open(fa).unwrap_or_else(|e| panic!("fail to open file '{fa}' : {e}")))
    };
    let reader: Box<dyn Read> = gzip_or_plain_reader_from_read(input)
        .unwrap_or_else(|e| panic!("fail to open file '{fa}' : {e}"));

    let t = Instant::now();
    eprint!("[bwa_index] Pack FASTA... ");
    let _l_pac = bns_fasta2bntseq(BufReader::new(reader), prefix, 1);
    eprintln!("{:.2} sec", t.elapsed().as_secs_f32());

    let mut fmi = FMI_search::ctor(prefix);
    let rc = fmi.build_index();
    fmi.dtor();
    rc
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::io::Write;
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::{bwa_idx_build, bwa_index};
    use flate2::write::GzEncoder;
    use flate2::Compression;

    fn temp_dir() -> std::path::PathBuf {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_bwtindex_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        dir
    }

    #[test]
    fn bwa_idx_build_writes_pack_and_marker_files() {
        let dir = temp_dir();
        let fasta = dir.join("ref.fa");
        let prefix = dir.join("idx");
        fs::write(&fasta, b">chr1\nACGT\n").expect("write fasta");

        assert_eq!(
            bwa_idx_build(
                fasta.to_str().expect("utf8"),
                prefix.to_str().expect("utf8"),
            ),
            0
        );
        assert!(std::path::Path::new(&format!("{}.pac", prefix.display())).is_file());
        assert!(std::path::Path::new(&format!("{}.bwt.2bit.64", prefix.display())).is_file());
        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn bwa_idx_build_detects_gzip_without_suffix() {
        let dir = temp_dir();
        let fasta = dir.join("ref-compressed");
        let prefix = dir.join("idx_gz_magic");
        let mut gz = GzEncoder::new(Vec::new(), Compression::default());
        gz.write_all(b">chr1\nACGT\n").expect("write gz");
        fs::write(&fasta, gz.finish().expect("finish gz")).expect("write fasta");

        assert_eq!(
            bwa_idx_build(
                fasta.to_str().expect("utf8"),
                prefix.to_str().expect("utf8"),
            ),
            0
        );
        assert!(std::path::Path::new(&format!("{}.pac", prefix.display())).is_file());
        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn bwa_index_parses_prefix_flag() {
        let dir = temp_dir();
        let fasta = dir.join("ref.fa");
        let prefix = dir.join("pref");
        fs::write(&fasta, b">chr1\nACGT\n").expect("write fasta");
        let argv = vec![
            "bwa-mem2".to_string(),
            "-p".to_string(),
            prefix.to_str().expect("utf8").to_string(),
            fasta.to_str().expect("utf8").to_string(),
        ];
        assert_eq!(bwa_index(&argv), 0);
        assert!(std::path::Path::new(&format!("{}.pac", prefix.display())).is_file());
        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn bwa_index_accepts_prefix_after_fasta_like_getopt_permutation() {
        let dir = temp_dir();
        let fasta = dir.join("ref.fa");
        let prefix = dir.join("pref_late");
        fs::write(&fasta, b">chr1\nACGT\n").expect("write fasta");
        let argv = vec![
            "bwa-mem2".to_string(),
            fasta.to_str().expect("utf8").to_string(),
            "-p".to_string(),
            prefix.to_str().expect("utf8").to_string(),
        ];
        assert_eq!(bwa_index(&argv), 0);
        assert!(std::path::Path::new(&format!("{}.pac", prefix.display())).is_file());
        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn bwa_index_honors_getopt_option_terminator() {
        let dir = temp_dir();
        let fasta = dir.join("ref.fa");
        let prefix = dir.join("pref_dashdash");
        fs::write(&fasta, b">chr1\nACGT\n").expect("write fasta");
        let argv = vec![
            "bwa-mem2".to_string(),
            "-p".to_string(),
            prefix.to_str().expect("utf8").to_string(),
            "--".to_string(),
            fasta.to_str().expect("utf8").to_string(),
        ];
        assert_eq!(bwa_index(&argv), 0);
        assert!(std::path::Path::new(&format!("{}.pac", prefix.display())).is_file());
        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn bwa_index_accepts_attached_prefix_value_like_getopt() {
        let dir = temp_dir();
        let fasta = dir.join("ref.fa");
        let prefix = dir.join("pref_attached");
        fs::write(&fasta, b">chr1\nACGT\n").expect("write fasta");
        let argv = vec![
            "bwa-mem2".to_string(),
            format!("-p{}", prefix.to_str().expect("utf8")),
            fasta.to_str().expect("utf8").to_string(),
        ];
        assert_eq!(bwa_index(&argv), 0);
        assert!(std::path::Path::new(&format!("{}.pac", prefix.display())).is_file());
        fs::remove_dir_all(&dir).expect("cleanup");
    }
}
