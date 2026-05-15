#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Port of `bwa-mem2/src/kthread.h` + `bwa-mem2/src/kthread.cpp`.

// --- kthread.h ---

#[doc = "Original struct: kt_for_t (bwa-mem2/src/kthread.h)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct kt_for_t {
    pub _opaque: (),
}

#[doc = "Original struct: ktf_worker_t (bwa-mem2/src/kthread.h)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct ktf_worker_t {
    pub _opaque: (),
}

#[doc = "Original struct: ktp_t (bwa-mem2/src/kthread.h)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct ktp_t {
    pub _opaque: (),
}

#[doc = "Original struct: ktp_worker_t (bwa-mem2/src/kthread.h)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct ktp_worker_t {
    pub _opaque: (),
}

// --- kthread.cpp ---

#[doc = "Original function: steal_work:41"]
pub fn steal_work(_arg0: crate::support::Opaque) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("steal_work")
}

#[doc = "Original function: ktf_worker:53"]
pub fn ktf_worker(_arg0: crate::support::Opaque) {
    crate::support::stub::<()>("ktf_worker")
}

#[doc = "Original function: kt_for:80"]
pub fn kt_for(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
    _arg2: crate::support::Opaque,
) {
    crate::support::stub::<()>("kt_for")
}
