#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Generated scaffold for `bwa-mem2/src/utils.h`.

#[doc = "Original struct: pair64_t (bwa-mem2/src/utils.h)"]
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct pair64_t {
    pub x: u64,
    pub y: u64,
}

#[doc = "Original struct: uint64_v (bwa-mem2/src/utils.h)"]
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct uint64_v {
    pub n: usize,
    pub m: usize,
    pub a: Vec<u64>,
}

#[doc = "Original struct: pair64_v (bwa-mem2/src/utils.h)"]
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct pair64_v {
    pub n: usize,
    pub m: usize,
    pub a: Vec<pair64_t>,
}

#[doc = "Original function: __rdtsc:53"]
pub fn rdtsc__L53(_arg0: crate::support::Opaque) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("__rdtsc")
}

#[doc = "Original function: __rdtsc:60"]
pub fn rdtsc__L60(_arg0: crate::support::Opaque) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("__rdtsc")
}

#[doc = "Original function: hash_64:117"]
#[inline]
pub fn hash_64(mut key: u64) -> u64 {
    key = key.wrapping_add(!(key << 32));
    key ^= key >> 22;
    key = key.wrapping_add(!(key << 13));
    key ^= key >> 8;
    key = key.wrapping_add(key << 3);
    key ^= key >> 15;
    key = key.wrapping_add(!(key << 27));
    key ^= key >> 31;
    key
}

#[cfg(test)]
mod tests {
    use super::hash_64;

    #[test]
    fn hash_64_matches_reference_values() {
        assert_eq!(hash_64(0), 7_654_268_697_807_496_793);
        assert_eq!(hash_64(1), 2_320_827_452_992_767_577);
        assert_eq!(hash_64(0x0123_4567_89ab_cdef), 2_602_480_231_338_512_580);
    }
}
