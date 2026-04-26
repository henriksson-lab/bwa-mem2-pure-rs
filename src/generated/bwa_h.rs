#![allow(dead_code, non_snake_case, non_camel_case_types, non_upper_case_globals)]

//! Generated scaffold for `bwa-mem2/src/bwa.h`.

#[doc = "Original struct: bseq1_t (bwa-mem2/src/bwa.h)"]
#[derive(Debug, Default, Clone)]
pub struct bseq1_t {
    pub l_seq: i32,
    pub id: i32,
    pub name: Option<String>,
    pub comment: Option<String>,
    pub seq: Option<String>,
    pub qual: Option<String>,
    pub sam: Option<String>,
    pub seq_nt4: Vec<u8>,
}

#[doc = "Original struct: bwaidx_t (bwa-mem2/src/bwa.h)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct bwaidx_t {
    pub _opaque: (),
}
