#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Port of `bwa-mem2/src/memcpy_bwamem.h` + `bwa-mem2/src/memcpy_bwamem.cpp`.


// --- memcpy_bwamem.cpp ---

#[doc = "Original function: memcpy_bwamem:32"]
pub fn memcpy_bwamem(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
    _arg2: crate::support::Opaque,
    _arg3: crate::support::Opaque,
    _arg4: crate::support::Opaque,
    _arg5: crate::support::Opaque,
) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("memcpy_bwamem")
}
