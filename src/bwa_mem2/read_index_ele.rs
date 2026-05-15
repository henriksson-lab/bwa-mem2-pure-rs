#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Port of `bwa-mem2/src/read_index_ele.h` + `bwa-mem2/src/read_index_ele.cpp`.

use crate::bwa_mem2::bntseq::{bns_destroy, bns_restore, bntseq_t};
use std::io::Read;
use std::path::Path;

// --- read_index_ele.h ---

pub const BWA_IDX_BWT: i32 = 0x1;
pub const BWA_IDX_BNS: i32 = 0x2;
pub const BWA_IDX_PAC: i32 = 0x4;
pub const BWA_IDX_ALL: i32 = 0x7;

#[doc = "Original struct: bwaidx_fm_t (bwa-mem2/src/read_index_ele.h)"]
#[derive(Debug, Default)]
pub struct bwaidx_fm_t {
    pub bns: Option<bntseq_t>,
    pub pac: Vec<u8>,
    pub is_shm: i32,
    pub l_mem: i64,
    pub mem: Vec<u8>,
}

#[doc = "Original class: indexEle (bwa-mem2/src/read_index_ele.h)"]
#[derive(Debug, Default)]
pub struct indexEle {
    pub idx: bwaidx_fm_t,
}

// --- read_index_ele.cpp ---

impl indexEle {
    #[doc = "Original function: indexEle::indexEle:40"]
    pub fn ctor() -> Self {
        Self {
            idx: bwaidx_fm_t::default(),
        }
    }

    #[doc = "Original function: indexEle::~indexEle:46"]
    pub fn dtor(&mut self) {
        if self.idx.mem.is_empty() {
            if let Some(bns) = self.idx.bns.as_mut() {
                bns_destroy(bns);
            }
            self.idx.bns = None;
            self.idx.pac.clear();
        } else {
            if let Some(bns) = self.idx.bns.as_mut() {
                bns.anns.clear();
                bns.n_seqs = 0;
            }
            if self.idx.is_shm == 0 {
                self.idx.mem.clear();
            }
            self.idx.bns = None;
        }
    }

    #[doc = "Original function: indexEle::bwa_idx_load_ele:60"]
    pub fn bwa_idx_load_ele(&mut self, hint: &str, which: i32) {
        let prefix = hint.to_string();
        eprintln!("* Index prefix: {}", prefix);

        if which & BWA_IDX_BNS != 0 {
            let bns = bns_restore(&prefix);
            let alt_count = bns.anns.iter().filter(|ann| ann.is_alt != 0).count();
            eprintln!("* Read {} ALT contigs", alt_count);

            if which & BWA_IDX_PAC != 0 {
                let pac_len = usize::try_from(bns.l_pac / 4 + 1).expect("pac length overflow");
                let mut pac = vec![0; pac_len];
                let mut file = bns
                    .fp_pac
                    .as_ref()
                    .expect("missing pac file")
                    .try_clone()
                    .expect("clone pac file");
                // concatenated 2-bit encoded sequence
                file.read_exact(&mut pac).expect("read pac");
                self.idx.pac = pac;
            }

            self.idx.bns = Some(bns);
        }
    }

    #[doc = "Original function: indexEle::bwa_idx_infer_prefix:96"]
    pub fn bwa_idx_infer_prefix(&self, hint: &str) -> Option<String> {
        let suffix64 = format!("{hint}.64.bwt");
        if Path::new(&suffix64).is_file() {
            return Some(format!("{hint}.64"));
        }

        let suffix = format!("{hint}.bwt");
        if Path::new(&suffix).is_file() {
            return Some(hint.to_string());
        }

        None
    }
}

#[cfg(test)]
mod tests {
    use std::fs::{self, File};
    use std::io::Write;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    use crate::bwa_mem2::bntseq::{bntann1_t, bntseq_t};

    use super::indexEle;

    fn temp_prefix(name: &str) -> PathBuf {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_{name}_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        dir.join("ref")
    }

    #[test]
    fn infer_prefix_prefers_64_bwt() {
        let prefix = temp_prefix("infer64");
        File::create(prefix.with_extension("64.bwt")).expect("create .64.bwt");
        File::create(prefix.with_extension("bwt")).expect("create .bwt");
        let ele = indexEle::ctor();
        assert_eq!(
            ele.bwa_idx_infer_prefix(prefix.to_str().expect("utf8 path")),
            Some(format!("{}.64", prefix.display()))
        );
        fs::remove_dir_all(prefix.parent().expect("temp dir")).expect("cleanup");
    }

    #[test]
    fn infer_prefix_falls_back_to_plain_bwt() {
        let prefix = temp_prefix("infer");
        File::create(prefix.with_extension("bwt")).expect("create .bwt");
        let ele = indexEle::ctor();
        assert_eq!(
            ele.bwa_idx_infer_prefix(prefix.to_str().expect("utf8 path")),
            Some(prefix.display().to_string())
        );
        fs::remove_dir_all(prefix.parent().expect("temp dir")).expect("cleanup");
    }

    #[test]
    fn infer_prefix_returns_none_when_missing() {
        let prefix = temp_prefix("missing");
        let ele = indexEle::ctor();
        assert_eq!(
            ele.bwa_idx_infer_prefix(prefix.to_str().expect("utf8 path")),
            None
        );
        fs::remove_dir_all(prefix.parent().expect("temp dir")).expect("cleanup");
    }

    #[test]
    fn destructor_clears_owned_bns_and_pac() {
        let mut ele = indexEle::ctor();
        ele.idx.pac = vec![1, 2, 3];
        ele.idx.bns = Some(bntseq_t {
            n_seqs: 1,
            anns: vec![bntann1_t::default()],
            ..Default::default()
        });
        ele.dtor();
        assert!(ele.idx.bns.is_none());
        assert!(ele.idx.pac.is_empty());
    }

    #[test]
    fn destructor_keeps_shm_mem_but_drops_bns() {
        let mut ele = indexEle::ctor();
        ele.idx.is_shm = 1;
        ele.idx.mem = vec![1, 2, 3];
        ele.idx.bns = Some(bntseq_t {
            n_seqs: 1,
            anns: vec![bntann1_t::default()],
            ..Default::default()
        });
        ele.dtor();
        assert!(ele.idx.bns.is_none());
        assert_eq!(ele.idx.mem, vec![1, 2, 3]);
    }

    #[test]
    fn load_ele_reads_bns_and_pac() {
        let prefix = temp_prefix("load_ele");
        fs::write(
            prefix.with_extension("ann"),
            "12 2 7\n1 chr1 comment one\n0 5 0\n2 chr2 comment two\n5 7 1\n",
        )
        .expect("write ann");
        fs::write(prefix.with_extension("amb"), "12 2 0\n").expect("write amb");
        let mut pac = File::create(prefix.with_extension("pac")).expect("write pac");
        pac.write_all(&[1, 2, 3, 4]).expect("pac bytes");
        fs::write(prefix.with_extension("alt"), "chr2\n").expect("write alt");

        let mut ele = indexEle::ctor();
        ele.bwa_idx_load_ele(
            prefix.to_str().expect("utf8"),
            super::BWA_IDX_BNS | super::BWA_IDX_PAC,
        );

        assert_eq!(ele.idx.pac, vec![1, 2, 3, 4]);
        let bns = ele.idx.bns.as_ref().expect("loaded bns");
        assert_eq!(bns.anns.len(), 2);
        assert_eq!(bns.anns[1].is_alt, 1);
        fs::remove_dir_all(prefix.parent().expect("temp dir")).expect("cleanup");
    }
}
