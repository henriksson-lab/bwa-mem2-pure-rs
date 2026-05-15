#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Port of `bwa-mem2/src/ksort.h`.

// --- ksort.h ---

#[doc = "Original struct: ks_isort_stack_t (bwa-mem2/src/ksort.h)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct ks_isort_stack_t {
    pub _opaque: (),
}
