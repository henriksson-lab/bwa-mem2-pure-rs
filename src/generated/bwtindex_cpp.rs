#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Generated scaffold for `bwa-mem2/src/bwtindex.cpp`.

use std::fs::File;
use std::io::{BufReader, Read};
use std::time::Instant;

use flate2::read::MultiGzDecoder;

use crate::generated::bntseq_cpp::bns_fasta2bntseq;
use crate::generated::fmi_search_cpp::FMI_search;

#[doc = "Original function: bwa_index:43"]
pub fn bwa_index(argv: &[String]) -> i32 {
    let mut prefix: Option<String> = None;
    let mut i = 1;
    while i < argv.len() {
        if argv[i] == "-p" {
            i += 1;
            if i >= argv.len() {
                return 1;
            }
            prefix = Some(argv[i].clone());
        } else {
            break;
        }
        i += 1;
    }

    if i >= argv.len() {
        eprintln!("Usage: bwa-mem2 index [-p prefix] <in.fasta>");
        return 1;
    }

    let fa = &argv[i];
    let prefix = prefix.unwrap_or_else(|| fa.clone());
    bwa_idx_build(fa, &prefix)
}

#[doc = "Original function: bwa_idx_build:61"]
pub fn bwa_idx_build(fa: &str, prefix: &str) -> i32 {
    let file = File::open(fa).unwrap_or_else(|e| panic!("fail to open file '{fa}' : {e}"));
    let reader: Box<dyn Read> = if fa.ends_with(".gz") {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };

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
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::{bwa_idx_build, bwa_index};

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
}
