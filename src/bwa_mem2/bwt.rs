#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Port of `bwa-mem2/src/bwt.h`.

// --- bwt.h ---

#[doc = "Original struct: bwt2_t (bwa-mem2/src/bwt.h)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct bwt2_t {
    pub _opaque: (),
}

#[doc = "Original struct: bwt_t (bwa-mem2/src/bwt.h)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct bwt_t {
    pub _opaque: (),
}

#[doc = "Original struct: bwtintv_t (bwa-mem2/src/bwt.h)"]
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct bwtintv_t {
    pub x: [u64; 3],
    pub info: u64,
}

#[doc = "Original struct: bwtintv_v (bwa-mem2/src/bwt.h)"]
#[derive(Debug, Default, Clone)]
pub struct bwtintv_v {
    pub n: usize,
    pub m: usize,
    pub a: Vec<bwtintv_t>,
}
