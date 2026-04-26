#![allow(dead_code, non_snake_case, non_camel_case_types, non_upper_case_globals)]

//! Generated scaffold for `bwa-mem2/src/read_index_ele.h`.

use crate::generated::bntseq_h::bntseq_t;

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
