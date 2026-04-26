#![allow(dead_code, non_snake_case, non_camel_case_types, non_upper_case_globals)]

//! Generated scaffold for `bwa-mem2/src/ksw.h`.

pub const KSW_XBYTE: i32 = 0x10000;
pub const KSW_XSTOP: i32 = 0x20000;
pub const KSW_XSUBO: i32 = 0x40000;
pub const KSW_XSTART: i32 = 0x80000;

#[doc = "Original struct: _kswq_t (bwa-mem2/src/ksw.h)"]
#[derive(Debug, Default, Clone)]
pub struct _kswq_t {
    pub qlen: i32,
    pub slen: i32,
    pub shift: u8,
    pub mdiff: u8,
    pub max: u8,
    pub size: u8,
    pub query: Vec<u8>,
    pub m: i32,
    pub mat: Vec<i8>,
}

#[doc = "Original struct: kswr_t (bwa-mem2/src/ksw.h)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct kswr_t {
    pub score: i32,
    pub te: i32,
    pub qe: i32,
    pub score2: i32,
    pub te2: i32,
    pub tb: i32,
    pub qb: i32,
}
