#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Generated scaffold for `bwa-mem2/src/bwamem.h`.

use crate::generated::bandedswa_h::SeqPair;
use crate::generated::bwa_h::bseq1_t;
use crate::generated::bwt_h::bwtintv_v;
use crate::generated::fmi_search_cpp::FMI_search;
use crate::generated::fmi_search_h::SMEM;

#[doc = "Original struct: abc (bwa-mem2/src/bwamem.h)"]
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct abc {
    pub rbeg: i64,
    pub qbeg: i32,
    pub len: i32,
    pub score: i32,
    pub done: i8,
    pub aln: i32,
}

#[doc = "Original struct: mem_aln_t (bwa-mem2/src/bwamem.h)"]
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct mem_aln_t {
    pub pos: i64,
    pub rid: i32,
    pub flag: i32,
    pub is_rev: u32,
    pub is_alt: u32,
    pub mapq: u32,
    pub NM: u32,
    pub n_cigar: i32,
    pub cigar: Vec<u32>,
    pub md: Vec<u8>,
    pub XA: Option<String>,
    pub score: i32,
    pub sub: i32,
    pub alt_sc: i32,
}

#[doc = "Original struct: mem_alnreg_t (bwa-mem2/src/bwamem.h)"]
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct mem_alnreg_t {
    pub rb: i64,
    pub re: i64,
    pub qb: i32,
    pub qe: i32,
    pub rid: i32,
    pub c: usize,
    pub score: i32,
    pub truesc: i32,
    pub sub: i32,
    pub alt_sc: i32,
    pub csub: i32,
    pub sub_n: i32,
    pub w: i32,
    pub seedcov: i32,
    pub secondary: i32,
    pub secondary_all: i32,
    pub seedlen0: i32,
    pub n_comp: u32,
    pub is_alt: u32,
    pub frac_rep: f32,
    pub hash: u64,
    pub flg: i32,
}

#[doc = "Original struct: mem_alnreg_v (bwa-mem2/src/bwamem.h)"]
#[derive(Debug, Default, Clone, PartialEq)]
pub struct mem_alnreg_v {
    pub n: usize,
    pub m: usize,
    pub a: Vec<mem_alnreg_t>,
}

#[doc = "Original struct: mem_cache (bwa-mem2/src/bwamem.h)"]
#[derive(Debug, Default, Clone)]
pub struct mem_cache {
    pub seqPairArrayAux: Vec<Vec<SeqPair>>,
    pub seqPairArrayLeft128: Vec<Vec<SeqPair>>,
    pub seqPairArrayRight128: Vec<Vec<SeqPair>>,
    pub wsize: Vec<i64>,
    pub wsize_buf_ref: Vec<i64>,
    pub wsize_buf_qer: Vec<i64>,
    pub seqBufLeftRef: Vec<Vec<u8>>,
    pub seqBufRightRef: Vec<Vec<u8>>,
    pub seqBufLeftQer: Vec<Vec<u8>>,
    pub seqBufRightQer: Vec<Vec<u8>>,
    pub matchArray: Vec<Vec<SMEM>>,
    pub min_intv_ar: Vec<Vec<i32>>,
    pub rid: Vec<Vec<i32>>,
    pub lim: Vec<Vec<i32>>,
    pub query_pos_ar: Vec<Vec<i16>>,
    pub enc_qdb: Vec<Vec<u8>>,
    pub wsize_mem: Vec<i64>,
    pub wsize_mem_s: Vec<i64>,
    pub wsize_mem_r: Vec<i64>,
}

#[doc = "Original struct: mem_chain_t (bwa-mem2/src/bwamem.h)"]
#[derive(Debug, Default, Clone, PartialEq)]
pub struct mem_chain_t {
    pub seqid: i32,
    pub cseed: i32,
    pub n: i32,
    pub m: i32,
    pub first: i32,
    pub rid: i32,
    pub w: u32,
    pub kept: u32,
    pub is_alt: u32,
    pub frac_rep: f32,
    pub pos: i64,
    pub seeds: Vec<mem_seed_t>,
}

#[doc = "Original struct: mem_chain_v (bwa-mem2/src/bwamem.h)"]
#[derive(Debug, Default, Clone, PartialEq)]
pub struct mem_chain_v {
    pub n: usize,
    pub m: usize,
    pub cc: usize,
    pub a: Vec<mem_chain_t>,
}

#[doc = "Original struct: mem_opt_t (bwa-mem2/src/bwamem.h)"]
#[derive(Debug, Clone)]
pub struct mem_opt_t {
    pub a: i32,
    pub b: i32,
    pub o_del: i32,
    pub e_del: i32,
    pub o_ins: i32,
    pub e_ins: i32,
    pub pen_unpaired: i32,
    pub pen_clip5: i32,
    pub pen_clip3: i32,
    pub w: i32,
    pub zdrop: i32,
    pub max_mem_intv: u64,
    pub T: i32,
    pub flag: i32,
    pub min_seed_len: i32,
    pub min_chain_weight: i32,
    pub max_chain_extend: i32,
    pub split_factor: f32,
    pub split_width: i32,
    pub max_occ: i32,
    pub max_chain_gap: i32,
    pub n_threads: i32,
    pub chunk_size: i64,
    pub mask_level: f32,
    pub drop_ratio: f32,
    pub XA_drop_ratio: f32,
    pub mask_level_redun: f32,
    pub mapQ_coef_len: f32,
    pub mapQ_coef_fac: i32,
    pub max_ins: i32,
    pub max_matesw: i32,
    pub max_XA_hits: i32,
    pub max_XA_hits_alt: i32,
    pub mat: [i8; 25],
}

impl Default for mem_opt_t {
    fn default() -> Self {
        Self {
            a: 0,
            b: 0,
            o_del: 0,
            e_del: 0,
            o_ins: 0,
            e_ins: 0,
            pen_unpaired: 0,
            pen_clip5: 0,
            pen_clip3: 0,
            w: 0,
            zdrop: 0,
            max_mem_intv: 0,
            T: 0,
            flag: 0,
            min_seed_len: 0,
            min_chain_weight: 0,
            max_chain_extend: 0,
            split_factor: 0.0,
            split_width: 0,
            max_occ: 0,
            max_chain_gap: 0,
            n_threads: 0,
            chunk_size: 0,
            mask_level: 0.0,
            drop_ratio: 0.0,
            XA_drop_ratio: 0.0,
            mask_level_redun: 0.0,
            mapQ_coef_len: 0.0,
            mapQ_coef_fac: 0,
            max_ins: 0,
            max_matesw: 0,
            max_XA_hits: 0,
            max_XA_hits_alt: 0,
            mat: [0; 25],
        }
    }
}

#[doc = "Original struct: mem_pestat_t (bwa-mem2/src/bwamem.h)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct mem_pestat_t {
    pub low: i32,
    pub high: i32,
    pub failed: i32,
    pub avg: f64,
    pub std: f64,
}

#[doc = "Original struct: mem_seed_t (bwa-mem2/src/bwamem.h)"]
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct mem_seed_t {
    pub rbeg: i64,
    pub qbeg: i32,
    pub len: i32,
    pub score: i32,
    pub done: i8,
    pub aln: i32,
}

#[doc = "Original struct: smem_aux_t (bwa-mem2/src/bwamem.h)"]
#[derive(Debug, Default, Clone)]
pub struct smem_aux_t {
    pub mem: bwtintv_v,
    pub mem1: bwtintv_v,
    pub tmpv: [Box<bwtintv_v>; 2],
}

#[doc = "Original struct: worker_t (bwa-mem2/src/bwamem.h)"]
#[derive(Debug, Default)]
pub struct worker_t {
    pub opt: Option<Box<mem_opt_t>>,
    pub pes: [mem_pestat_t; 4],
    pub aux: Vec<smem_aux_t>,
    pub seqs: Vec<bseq1_t>,
    pub regs: Vec<mem_alnreg_v>,
    pub n_processed: i64,
    pub chain_ar: Vec<mem_chain_v>,
    pub mmc: mem_cache,
    pub seedBuf: Vec<mem_seed_t>,
    pub seedBufSize: i64,
    pub auxSeedBuf: Vec<mem_seed_t>,
    pub auxSeedBufSize: i64,
    pub ref_string: Vec<u8>,
    pub nthreads: i16,
    pub nreads: i32,
    pub fmi: Option<FMI_search>,
}

#[doc = "Original function: abc:114"]
pub fn abc() {
    crate::support::stub::<()>("abc")
}
