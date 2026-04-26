#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Generated scaffold for `bwa-mem2/src/fastmap.h`.

use crate::generated::bwa_h::bseq1_t;
use crate::generated::bwamem_h::{mem_opt_t, mem_pestat_t};
use crate::generated::kseq_h::kseq_t;
use crate::generated::utils_cpp::ErrFile;

#[doc = "Original struct: ktp_aux_t (bwa-mem2/src/fastmap.h)"]
#[derive(Debug, Default)]
pub struct ktp_aux_t {
    pub ks: Option<kseq_t>,
    pub ks2: Option<kseq_t>,
    pub opt: Option<Box<mem_opt_t>>,
    pub pes0: Option<Box<[mem_pestat_t; 4]>>,
    pub n_processed: i64,
    pub copy_comment: i32,
    pub my_ntasks: i64,
    pub ntasks: i64,
    pub task_size: i64,
    pub actual_chunk_size: i64,
    pub fp: Option<ErrFile>,
    pub ref_string: Vec<u8>,
}

#[doc = "Original struct: ktp_data_t (bwa-mem2/src/fastmap.h)"]
#[derive(Debug, Default, Clone)]
pub struct ktp_data_t {
    pub n_seqs: i32,
    pub seqs: Vec<bseq1_t>,
    pub sam_lines: Vec<Box<str>>,
}
