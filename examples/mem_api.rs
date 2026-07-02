use std::env;
use std::path::PathBuf;

use bwa_mem2_pure_rs::mem_api::{MemAligner, MemReadPair};

fn main() -> Result<(), String> {
    let index_prefix = env::args_os()
        .nth(1)
        .map(PathBuf::from)
        .ok_or_else(|| "usage: cargo run --example mem_api -- <index-prefix>".to_string())?;

    let mut aligner = MemAligner::builder(&index_prefix).threads(2).build()?;

    let pairs = [MemReadPair {
        name: "example-read-1",
        r1: b"ACGTACGTACGTACGTACGTACGTACGTACGT",
        q1: b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        r2: b"TGCATGCATGCATGCATGCATGCATGCATGCA",
        q2: b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
    }];

    print!("{}", aligner.sam_header()?);
    aligner.align_pairs_into_indexed(&pairs, |pair_index, sam_record| {
        let _input_pair = &pairs[pair_index];
        print!("{sam_record}");
        Ok(())
    })?;

    Ok(())
}
