#![allow(dead_code, non_snake_case, non_camel_case_types, non_upper_case_globals)]

//! Generated scaffold for `bwa-mem2/src/kswv.h`.

pub const MAX_SEQ_LEN_REF_SAM: i32 = 2048;
pub const MAX_SEQ_LEN_QER_SAM: i32 = 512;
pub const SIMD_WIDTH8: usize = 64;
pub const SIMD_WIDTH16: usize = 32;
pub const DEFAULT_AMBIG: i8 = -1;

pub use crate::generated::bandedswa_h::SeqPair;

#[doc = "Original struct: dnaSeqPair (bwa-mem2/src/kswv.h)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct dnaSeqPair {
    pub _opaque: (),
}

#[doc = "Original struct: kswq_t (bwa-mem2/src/kswv.h)"]
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct kswq_t {
    pub qlen: i32,
    pub slen: i32,
    pub shift: u8,
    pub mdiff: u8,
    pub max: u8,
    pub size: u8,
    pub qp_u8: Vec<i8>,
    pub qp_i16: Vec<i16>,
}

pub use crate::generated::ksw_h::kswr_t;

#[doc = "Original class: kswv (bwa-mem2/src/kswv.h)"]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct kswv {
    pub m: i32,
    pub o_del: i32,
    pub o_ins: i32,
    pub e_del: i32,
    pub e_ins: i32,
    pub w_match: i8,
    pub w_mismatch: i8,
    pub w_open: i8,
    pub w_extend: i8,
    pub w_ambig: i8,
    pub F8: Vec<u8>,
    pub H8_0: Vec<u8>,
    pub H8_max: Vec<u8>,
    pub H8_1: Vec<u8>,
    pub rowMax8: Vec<u8>,
    pub F16: Vec<i16>,
    pub H16_0: Vec<i16>,
    pub H16_max: Vec<i16>,
    pub H16_1: Vec<i16>,
    pub rowMax16: Vec<i16>,
    pub maxRefLen: i32,
    pub maxQerLen: i32,
    pub g_qmax: i32,
    pub sort1Ticks: i64,
    pub setupTicks: i64,
    pub swTicks: i64,
    pub sort2Ticks: i64,
}

impl Default for kswv {
    fn default() -> Self {
        Self {
            m: 0,
            o_del: 0,
            o_ins: 0,
            e_del: 0,
            e_ins: 0,
            w_match: 0,
            w_mismatch: 0,
            w_open: 0,
            w_extend: 0,
            w_ambig: 0,
            F8: Vec::new(),
            H8_0: Vec::new(),
            H8_max: Vec::new(),
            H8_1: Vec::new(),
            rowMax8: Vec::new(),
            F16: Vec::new(),
            H16_0: Vec::new(),
            H16_max: Vec::new(),
            H16_1: Vec::new(),
            rowMax16: Vec::new(),
            maxRefLen: 0,
            maxQerLen: 0,
            g_qmax: 0,
            sort1Ticks: 0,
            setupTicks: 0,
            swTicks: 0,
            sort2Ticks: 0,
        }
    }
}
