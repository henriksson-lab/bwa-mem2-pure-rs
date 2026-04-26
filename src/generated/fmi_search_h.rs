#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Generated scaffold for `bwa-mem2/src/FMI_search.h`.

#[doc = "Original struct: CP_OCC (bwa-mem2/src/FMI_search.h)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct CP_OCC {
    pub _opaque: (),
}

#[doc = "Original class: FMI_search (bwa-mem2/src/FMI_search.h)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct FMI_search {
    pub _opaque: (),
}

#[doc = "Original struct: SMEM (bwa-mem2/src/FMI_search.h)"]
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct SMEM {
    pub rid: u32,
    pub m: u32,
    pub n: u32,
    pub k: i64,
    pub l: i64,
    pub s: i64,
}

#[doc = "Original struct: checkpoint_occ_scalar (bwa-mem2/src/FMI_search.h)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct checkpoint_occ_scalar {
    pub _opaque: (),
}

#[doc = "Original struct: smem_struct (bwa-mem2/src/FMI_search.h)"]
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct smem_struct {
    pub rid: u32,
    pub m: u32,
    pub n: u32,
    pub k: i64,
    pub l: i64,
    pub s: i64,
}

#[doc = "Original function: _mm_countbits_64:61"]
pub fn mm_countbits_64(_arg0: crate::support::Opaque) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("_mm_countbits_64")
}
