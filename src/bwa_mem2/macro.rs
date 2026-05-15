#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Port of `bwa-mem2/src/macro.h`.

// --- macro.h ---

pub const BATCH_SIZE: usize = 512;
pub const CACHE_LINE: usize = 16;
pub const MAX_THREADS: usize = 256;
pub const SEEDS_PER_CHAIN: usize = 1;
