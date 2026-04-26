#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Generated scaffold for `bwa-mem2/src/bntseq.h`.

use std::fs::File;

#[doc = "Original struct: bntann1_t (bwa-mem2/src/bntseq.h)"]
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct bntann1_t {
    pub offset: i64,
    pub len: i32,
    pub n_ambs: i32,
    pub gi: u32,
    pub is_alt: i32,
    pub name: String,
    pub anno: String,
}

#[doc = "Original struct: bntamb1_t (bwa-mem2/src/bntseq.h)"]
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct bntamb1_t {
    pub offset: i64,
    pub len: i32,
    pub amb: u8,
}

#[doc = "Original struct: bntseq_t (bwa-mem2/src/bntseq.h)"]
#[derive(Debug, Default)]
pub struct bntseq_t {
    pub l_pac: i64,
    pub n_seqs: i32,
    pub seed: u32,
    pub anns: Vec<bntann1_t>,
    pub n_holes: i32,
    pub ambs: Vec<bntamb1_t>,
    pub fp_pac: Option<File>,
}

#[doc = "Original function: bns_depos:87"]
#[inline]
pub fn bns_depos(bns: &bntseq_t, pos: i64, is_rev: &mut i32) -> i64 {
    if pos >= bns.l_pac {
        *is_rev = 1;
        (bns.l_pac << 1) - 1 - pos
    } else {
        *is_rev = 0;
        pos
    }
}

#[cfg(test)]
mod tests {
    use super::{bns_depos, bntseq_t};

    #[test]
    fn bns_depos_tracks_forward_and_reverse_coordinates() {
        let bns = bntseq_t {
            l_pac: 10,
            ..Default::default()
        };
        let mut is_rev = -1;
        assert_eq!(bns_depos(&bns, 3, &mut is_rev), 3);
        assert_eq!(is_rev, 0);
        assert_eq!(bns_depos(&bns, 12, &mut is_rev), 7);
        assert_eq!(is_rev, 1);
    }
}
