#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Generated scaffold for `bwa-mem2/src/bandedSWA.h`.

#[doc = "Original class: BandedPairWiseSW (bwa-mem2/src/bandedSWA.h)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct BandedPairWiseSW {
    pub _opaque: (),
}

#[doc = "Original struct: OutScore (bwa-mem2/src/bandedSWA.h)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct OutScore {
    pub _opaque: (),
}

#[doc = "Original struct: SeqPair (bwa-mem2/src/bandedSWA.h)"]
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct SeqPair {
    pub idr: i32,
    pub idq: i32,
    pub id: i32,
    pub len1: i32,
    pub len2: i32,
    pub h0: i32,
    pub seqid: i32,
    pub regid: i32,
    pub score: i32,
    pub tle: i32,
    pub gtle: i32,
    pub qle: i32,
    pub gscore: i32,
    pub max_off: i32,
}

#[doc = "Original struct: dnaOutScore (bwa-mem2/src/bandedSWA.h)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct dnaOutScore {
    pub _opaque: (),
}

#[doc = "Original struct: dnaSeqPair (bwa-mem2/src/bandedSWA.h)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct dnaSeqPair {
    pub _opaque: (),
}

#[doc = "Original struct: eh_t (bwa-mem2/src/bandedSWA.h)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct eh_t {
    pub _opaque: (),
}
