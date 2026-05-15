#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Port of `bwa-mem2/src/kvec.h`.

// --- kvec.h ---

#[doc = "Original struct: kvecint (bwa-mem2/src/kvec.h)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct kvecint {
    pub _opaque: (),
}
