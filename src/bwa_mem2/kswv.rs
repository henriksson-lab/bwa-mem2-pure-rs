#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Port of `bwa-mem2/src/kswv.h` + `bwa-mem2/src/kswv.cpp`.

use crate::bwa_mem2::ksw::{KSW_XBYTE, KSW_XSTART, ksw_align2, ksw_i16_slices, ksw_qmax, ksw_u8_slices};

// --- kswv.h ---

pub const MAX_SEQ_LEN_REF_SAM: i32 = 2048;
pub const MAX_SEQ_LEN_QER_SAM: i32 = 512;
pub const SIMD_WIDTH8: usize = 64;
pub const SIMD_WIDTH16: usize = 32;
pub const DEFAULT_AMBIG: i8 = -1;

pub use crate::bwa_mem2::bandedswa::SeqPair;

#[doc = "Original struct: dnaSeqPair (bwa-mem2/src/kswv.h)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct dnaSeqPair {
    pub _opaque: (),
}

#[doc = "Original struct: kswq_t (bwa-mem2/src/kswv.h)"]
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct kswq_t {
    pub qlen: i32,
    pub slen: i32,
    pub shift: u8,
    pub mdiff: u8,
    pub max: u8,
    pub size: u8,
    pub qp_u8: Vec<i8>,
    pub qp_i16: Vec<i16>,
}

pub use crate::bwa_mem2::ksw::kswr_t;

#[doc = "Original class: kswv (bwa-mem2/src/kswv.h)"]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct kswv {
    pub m: i32,
    pub o_del: i32,
    pub o_ins: i32,
    pub e_del: i32,
    pub e_ins: i32,
    pub w_match: i8,
    pub w_mismatch: i8,
    pub w_open: i8,
    pub w_extend: i8,
    pub w_ambig: i8,
    pub F8: Vec<u8>,
    pub H8_0: Vec<u8>,
    pub H8_max: Vec<u8>,
    pub H8_1: Vec<u8>,
    pub rowMax8: Vec<u8>,
    pub F16: Vec<i16>,
    pub H16_0: Vec<i16>,
    pub H16_max: Vec<i16>,
    pub H16_1: Vec<i16>,
    pub rowMax16: Vec<i16>,
    pub maxRefLen: i32,
    pub maxQerLen: i32,
    pub g_qmax: i32,
    pub sort1Ticks: i64,
    pub setupTicks: i64,
    pub swTicks: i64,
    pub sort2Ticks: i64,
}

impl Default for kswv {
    fn default() -> Self {
        Self {
            m: 0,
            o_del: 0,
            o_ins: 0,
            e_del: 0,
            e_ins: 0,
            w_match: 0,
            w_mismatch: 0,
            w_open: 0,
            w_extend: 0,
            w_ambig: 0,
            F8: Vec::new(),
            H8_0: Vec::new(),
            H8_max: Vec::new(),
            H8_1: Vec::new(),
            rowMax8: Vec::new(),
            F16: Vec::new(),
            H16_0: Vec::new(),
            H16_max: Vec::new(),
            H16_1: Vec::new(),
            rowMax16: Vec::new(),
            maxRefLen: 0,
            maxQerLen: 0,
            g_qmax: 0,
            sort1Ticks: 0,
            setupTicks: 0,
            swTicks: 0,
            sort2Ticks: 0,
        }
    }
}

// --- kswv.cpp ---

#[cfg(target_arch = "x86_64")]
use core::arch::x86_64::*;

const G_DEFR: kswr_t = kswr_t {
    score: 0,
    te: -1,
    qe: -1,
    score2: -1,
    te2: -1,
    tb: -1,
    qb: -1,
};

// Thread-local scratch buffers for kswvBatchWrapper8/16. These were allocated per call
// (padded_pairs, seq1_soa, seq2_soa); per-thread reuse amortizes the heap traffic since
// kswv is invoked once per worker_sam call, ~50-200 times per program run.
thread_local! {
    static KSWV_BATCH_SCRATCH: std::cell::RefCell<(
        Vec<SeqPair>,
        Vec<u8>,
        Vec<u8>,
        Vec<i16>,
        Vec<i16>,
    )> = const {
        std::cell::RefCell::new((Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new()))
    };
}

#[cfg(target_arch = "x86_64")]
thread_local! {
    static KSWV16_AVX_SCRATCH: std::cell::RefCell<(
        Vec<i16>,
        Vec<i16>,
        Vec<i16>,
        Vec<i16>,
    )> = const {
        std::cell::RefCell::new((Vec::new(), Vec::new(), Vec::new(), Vec::new()))
    };
}

#[cfg(target_arch = "x86_64")]
thread_local! {
    static KSWV8_AVX_SCRATCH: std::cell::RefCell<(
        Vec<u8>,
        Vec<u8>,
        Vec<u8>,
        Vec<u8>,
    )> = const {
        std::cell::RefCell::new((Vec::new(), Vec::new(), Vec::new(), Vec::new()))
    };
}

fn disable_kswv8_avx512() -> bool {
    static FLAG: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *FLAG.get_or_init(|| std::env::var_os("BWA_DISABLE_KSWV8").is_some())
}

fn disable_kswv16_avx512() -> bool {
    static FLAG: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *FLAG.get_or_init(|| std::env::var_os("BWA_DISABLE_KSWV16").is_some())
}

#[doc = "Original function: parseCmdLine:1636"]
pub fn parseCmdLine(_arg0: crate::support::Opaque, _arg1: crate::support::Opaque) {
    crate::support::stub::<()>("parseCmdLine")
}

#[doc = "Original function: loadPairs:1686"]
pub fn loadPairs(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
    _arg2: crate::support::Opaque,
    _arg3: crate::support::Opaque,
) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("loadPairs")
}

#[doc = "Original function: find_stats:1738"]
pub fn find_stats(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
    _arg2: crate::support::Opaque,
    _arg3: crate::support::Opaque,
    _arg4: crate::support::Opaque,
) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("find_stats")
}

#[doc = "Original function: main:1751"]
pub fn main(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("main")
}

fn max_i8(a: i8, b: i8) -> i8 {
    if a > b {
        a
    } else {
        b
    }
}

impl kswv {
    // constructor
    #[doc = "Original function: kswv::kswv:114"]
    pub fn ctor(
        o_del: i32,
        e_del: i32,
        o_ins: i32,
        e_ins: i32,
        w_match: i8,
        w_mismatch: i8,
        numThreads: i32,
        maxRefLen: Option<i32>,
        maxQerLen: Option<i32>,
    ) -> Self {
        let num_threads = numThreads.max(1) as usize;
        let max_ref_len = maxRefLen.unwrap_or(MAX_SEQ_LEN_REF_SAM) + 16;
        let max_qer_len = maxQerLen.unwrap_or(MAX_SEQ_LEN_QER_SAM) + 16;
        let max_ref_us = max_ref_len.max(0) as usize;
        let max_qer_us = max_qer_len.max(0) as usize;
        let ref_cap16 = max_ref_us * SIMD_WIDTH16 * num_threads;
        let qer_cap16 = max_qer_us * SIMD_WIDTH16 * num_threads;
        let ref_cap8 = max_ref_us * SIMD_WIDTH8 * num_threads;
        let qer_cap8 = max_qer_us * SIMD_WIDTH8 * num_threads;
        let g_qmax = i32::from(max_i8(max_i8(w_match, w_mismatch), DEFAULT_AMBIG));
        Self {
            m: 5,
            o_del,
            o_ins,
            e_del,
            e_ins,
            w_match,
            w_mismatch,
            // w_open/w_extend are redundant copies of o_del/e_del; used in vector code.
            w_open: i8::try_from(o_del).unwrap_or(i8::MAX),
            w_extend: i8::try_from(e_del).unwrap_or(i8::MAX),
            w_ambig: DEFAULT_AMBIG,
            F8: vec![0; qer_cap8],
            H8_0: vec![0; qer_cap8],
            H8_max: vec![0; qer_cap8],
            H8_1: vec![0; qer_cap8],
            rowMax8: vec![0; ref_cap8],
            F16: vec![0; qer_cap16],
            H16_0: vec![0; qer_cap16],
            H16_max: vec![0; qer_cap16],
            H16_1: vec![0; qer_cap16],
            rowMax16: vec![0; ref_cap16],
            maxRefLen: max_ref_len,
            maxQerLen: max_qer_len,
            g_qmax,
            sort1Ticks: 0,
            setupTicks: 0,
            swTicks: 0,
            sort2Ticks: 0,
        }
    }

    // destructor
    #[doc = "Original function: kswv::~kswv:157"]
    pub fn dtor(&mut self) {
        self.F8.clear();
        self.H8_0.clear();
        self.H8_max.clear();
        self.H8_1.clear();
        self.rowMax8.clear();
        self.F16.clear();
        self.H16_0.clear();
        self.H16_max.clear();
        self.H16_1.clear();
        self.rowMax16.clear();
    }

    // Vector u8 Smith-Waterman entry point (AVX-512 path).
    #[doc = "Original function: kswv::getScores8:164"]
    pub fn getScores8(
        &self,
        pairArray: &mut [SeqPair],
        seqBufRef: &[u8],
        seqBufQer: &[u8],
        aln: &mut [kswr_t],
        numPairs: i32,
        nthreads: u16,
        phase: i32,
    ) {
        self.kswvBatchWrapper8(
            pairArray, seqBufRef, seqBufQer, aln, numPairs, nthreads, phase,
        );
    }

    // Pack SIMD_WIDTH8 SeqPairs into SoA (struct-of-arrays) buffers and run the AVX-512 u8
    // SW kernel. Rounds the batch up to SIMD_WIDTH8 with padded "id=ii, len1=len2=0" lanes
    // so the kernel can run a full vector per chunk.
    #[doc = "Original function: kswv::kswvBatchWrapper8:177"]
    pub fn kswvBatchWrapper8(
        &self,
        pairArray: &mut [SeqPair],
        seqBufRef: &[u8],
        seqBufQer: &[u8],
        aln: &mut [kswr_t],
        numPairs: i32,
        _numThreads: u16,
        phase: i32,
    ) {
        let total = numPairs as usize;
        KSWV_BATCH_SCRATCH.with(|cell| {
            let mut scratch = cell.borrow_mut();
            let (padded, seq1_soa, seq2_soa, _, _) = &mut *scratch;
            let ref_cap = self.maxRefLen as usize * SIMD_WIDTH8;
            let qer_cap = self.maxQerLen as usize * SIMD_WIDTH8;
            if seq1_soa.len() < ref_cap {
                seq1_soa.resize(ref_cap, 0xff);
            }
            if seq2_soa.len() < qer_cap {
                seq2_soa.resize(qer_cap, 0xff);
            }

            let mut idx = 0_usize;
            while idx + SIMD_WIDTH8 <= total {
                let lanes = &pairArray[idx..idx + SIMD_WIDTH8];
                self.kswv_batch8_chunk(
                    lanes, seqBufRef, seqBufQer, seq1_soa, seq2_soa, aln, idx, numPairs, phase,
                );
                idx += SIMD_WIDTH8;
            }
            if idx < total {
                padded.clear();
                padded.resize(SIMD_WIDTH8, SeqPair::default());
                let tail = &pairArray[idx..total];
                padded[..tail.len()].copy_from_slice(tail);
                for (lane, p) in padded[tail.len()..].iter_mut().enumerate() {
                    let id = idx + tail.len() + lane;
                    p.id = id as i32;
                    p.regid = p.id;
                }
                self.kswv_batch8_chunk(
                    &padded[..SIMD_WIDTH8],
                    seqBufRef,
                    seqBufQer,
                    seq1_soa,
                    seq2_soa,
                    aln,
                    idx,
                    numPairs,
                    phase,
                );
            }
        });
    }

    #[inline]
    fn kswv_batch8_chunk(
        &self,
        lanes: &[SeqPair],
        seqBufRef: &[u8],
        seqBufQer: &[u8],
        seq1_soa: &mut [u8],
        seq2_soa: &mut [u8],
        aln: &mut [kswr_t],
        offset: usize,
        numPairs: i32,
        phase: i32,
    ) {
        debug_assert_eq!(lanes.len(), SIMD_WIDTH8);
        let mut max_len1 = 0_i32;
        let mut max_len2 = 0_i32;

        // Compute max_len1/max_len2 first so we can size-validate the SoA writes once.
        for sp in lanes.iter().copied() {
            max_len1 = max_len1.max(sp.len1);
            let quanta = (sp.len2 + 16 - 1) / 16 * 16;
            max_len2 = max_len2.max(quanta);
        }
        let max_len1_usize = max_len1 as usize;
        let max_len2_usize = max_len2 as usize;
        // Caller (kswvBatchWrapper8) sized seq1_soa to maxRefLen*SIMD_WIDTH8 and seq2_soa to
        // maxQerLen*SIMD_WIDTH8. Per-row writes use indices k*SIMD_WIDTH8 + j with k <= max_len{1,2}
        // and j < SIMD_WIDTH8, so the maximum write offset is max_len*SIMD_WIDTH8 + (SIMD_WIDTH8-1)
        // = (max_len+1)*SIMD_WIDTH8 - 1. Validate the SoA buffers cover that range up front.
        let s1_used = (max_len1_usize + 1) * SIMD_WIDTH8;
        let s2_used = (max_len2_usize + 1) * SIMD_WIDTH8;
        assert!(seq1_soa.len() >= s1_used);
        assert!(seq2_soa.len() >= s2_used);
        // Pre-fill the chunk's used range to 0xff (the AMBIG sentinel). This replaces the per-lane
        // strided pad loops at indices k*SIMD_WIDTH8 + j for k in [len, max_len], which would issue
        // 64-stride byte writes; a contiguous fill is far more cache-friendly.
        seq1_soa[..s1_used].fill(0xff);
        seq2_soa[..s2_used].fill(0xff);
        let seq1_soa_ptr = seq1_soa.as_mut_ptr();
        let seq2_soa_ptr = seq2_soa.as_mut_ptr();

        // C++ kswv.cpp line 274: AMBR=AMBIG_=4 (identity for seq1). Skip the conditional.
        for (j, sp) in lanes.iter().copied().enumerate() {
            let idr = sp.idr as usize;
            let len1 = sp.len1 as usize;
            let seq1 = unsafe { seqBufRef.get_unchecked(idr..idr + len1) };
            for k in 0..len1 {
                let b = unsafe { *seq1.get_unchecked(k) };
                unsafe { *seq1_soa_ptr.add(k * SIMD_WIDTH8 + j) = b };
            }
        }

        // C++ kswv.cpp line 300: AMBQ=8 (so the seq2 conditional 4→8 is real, not identity).
        // Branchless `b + (b & 4)` yields b for b in {0,1,2,3} and 8 for b=4.
        for (j, sp) in lanes.iter().copied().enumerate() {
            let idq = sp.idq as usize;
            let len2 = sp.len2 as usize;
            let quanta = ((sp.len2 + 16 - 1) / 16 * 16) as usize;
            let seq2 = unsafe { seqBufQer.get_unchecked(idq..idq + len2) };
            for k in 0..len2 {
                let b = unsafe { *seq2.get_unchecked(k) };
                let v = b + (b & 4);
                unsafe { *seq2_soa_ptr.add(k * SIMD_WIDTH8 + j) = v };
            }
            for k in len2..quanta {
                unsafe { *seq2_soa_ptr.add(k * SIMD_WIDTH8 + j) = 5 };
            }
        }

        self.kswv512_u8(
            seq1_soa,
            seq2_soa,
            max_len1 as i16,
            max_len2 as i16,
            lanes,
            aln,
            offset as i32,
            0,
            numPairs,
            phase,
        );
    }

    // Vectorized 16-bit (i16 lane) SW entry point.
    #[doc = "Original function: kswv::getScores16:713"]
    pub fn getScores16(
        &self,
        pairArray: &mut [SeqPair],
        seqBufRef: &[u8],
        seqBufQer: &[u8],
        aln: &mut [kswr_t],
        numPairs: i32,
        nthreads: u16,
        phase: i32,
    ) {
        self.kswvBatchWrapper16(
            pairArray, seqBufRef, seqBufQer, aln, numPairs, nthreads, phase,
        );
    }

    // Pack SIMD_WIDTH16 SeqPairs into SoA buffers and run the AVX-512 i16 SW kernel; same
    // padding / rounding scheme as the u8 wrapper.
    #[doc = "Original function: kswv::kswvBatchWrapper16:725"]
    pub fn kswvBatchWrapper16(
        &self,
        pairArray: &mut [SeqPair],
        seqBufRef: &[u8],
        seqBufQer: &[u8],
        aln: &mut [kswr_t],
        numPairs: i32,
        _numThreads: u16,
        phase: i32,
    ) {
        let total = numPairs as usize;
        #[cfg(target_arch = "x86_64")]
        {
            if std::is_x86_feature_detected!("avx512bw") && !disable_kswv16_avx512() {
                KSWV_BATCH_SCRATCH.with(|cell| {
                    let mut scratch = cell.borrow_mut();
                    let (padded, _, _, seq1_soa, seq2_soa) = &mut *scratch;
                    let mut idx = 0_usize;
                    while idx + SIMD_WIDTH16 <= total {
                        let lanes = &pairArray[idx..idx + SIMD_WIDTH16];
                        self.kswv_batch16_chunk_avx512(
                            lanes, seqBufRef, seqBufQer, seq1_soa, seq2_soa, aln, idx, numPairs,
                            phase,
                        );
                        idx += SIMD_WIDTH16;
                    }
                    if idx < total {
                        padded.clear();
                        padded.resize(SIMD_WIDTH16, SeqPair::default());
                        let tail = &pairArray[idx..total];
                        padded[..tail.len()].copy_from_slice(tail);
                        for (lane, p) in padded[tail.len()..].iter_mut().enumerate() {
                            let id = idx + tail.len() + lane;
                            p.id = id as i32;
                            p.regid = p.id;
                        }
                        self.kswv_batch16_chunk_avx512(
                            &padded[..SIMD_WIDTH16],
                            seqBufRef,
                            seqBufQer,
                            seq1_soa,
                            seq2_soa,
                            aln,
                            idx,
                            numPairs,
                            phase,
                        );
                    }
                });
                return;
            }
        }
        let mut mat = [0_i8; 25];
        self.bwa_fill_scmat(&mut mat);
        let q_max = ksw_qmax(self.m, &mat);
        for chunk in pairArray[..total].chunks(SIMD_WIDTH16) {
            self.kswv_scalar_lanes_i16(chunk, seqBufRef, seqBufQer, aln, &mat, q_max, phase);
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[inline]
    fn kswv_batch16_chunk_avx512(
        &self,
        lanes: &[SeqPair],
        seqBufRef: &[u8],
        seqBufQer: &[u8],
        seq1_soa: &mut Vec<i16>,
        seq2_soa: &mut Vec<i16>,
        aln: &mut [kswr_t],
        offset: usize,
        numPairs: i32,
        phase: i32,
    ) {
        debug_assert_eq!(lanes.len(), SIMD_WIDTH16);
        let mut max_len1 = 0_i32;
        let mut max_len2 = 0_i32;

        for &sp in lanes {
            max_len1 = max_len1.max(sp.len1);
            let quanta = (sp.len2 + 8 - 1) / 8 * 8;
            max_len2 = max_len2.max(quanta);
        }
        let rows1 = (max_len1 + 1) as usize;
        let rows2 = (max_len2 + 1) as usize;
        seq1_soa.clear();
        seq1_soa.resize(rows1 * SIMD_WIDTH16, -1);
        seq2_soa.clear();
        seq2_soa.resize(rows2 * SIMD_WIDTH16, -1);

        // seq1_soa and seq2_soa were already filled with -1 by the resize() above. The trailing
        // pad loops only need to write the non-(-1) sentinels: 26 in the [len2, quanta) range
        // for seq2. The k in [len1, max_len1] and k in [quanta, max_len2] loops are redundant.
        // Branchless transforms for the AMBIG (b=4) case:
        //   seq1: 4→15. Decompose: 15-4=11; add 11 when bit-4 of b is set.
        //   seq2: 4→16. Decompose: 16-4=12; add 12 when bit-4 of b is set.
        for (lane, sp) in lanes.iter().copied().enumerate() {
            let seq1 = &seqBufRef[sp.idr as usize..];
            let seq1_ptr = seq1_soa.as_mut_ptr();
            let len1 = sp.len1 as usize;
            for k in 0..len1 {
                let base = unsafe { *seq1.get_unchecked(k) };
                let bit4 = (base & 4) >> 2; // 0 or 1
                let v = i16::from(base) + i16::from(bit4) * 11;
                unsafe { *seq1_ptr.add(k * SIMD_WIDTH16 + lane) = v };
            }

            let seq2 = &seqBufQer[sp.idq as usize..];
            let seq2_ptr = seq2_soa.as_mut_ptr();
            let len2 = sp.len2 as usize;
            let quanta = ((sp.len2 + 8 - 1) / 8 * 8) as usize;
            for k in 0..len2 {
                let base = unsafe { *seq2.get_unchecked(k) };
                let bit4 = (base & 4) >> 2; // 0 or 1
                let v = i16::from(base) + i16::from(bit4) * 12;
                unsafe { *seq2_ptr.add(k * SIMD_WIDTH16 + lane) = v };
            }
            for k in len2..quanta {
                unsafe { *seq2_ptr.add(k * SIMD_WIDTH16 + lane) = 26 };
            }
        }

        unsafe {
            self.kswv512_16_avx512(
                seq1_soa,
                seq2_soa,
                max_len1 as i16,
                max_len2 as i16,
                lanes,
                aln,
                offset as i32,
                numPairs,
                phase,
            );
        }
    }

    #[inline]
    fn kswv_scalar_lanes_i16(
        &self,
        lanes: &[SeqPair],
        seqBufRef: &[u8],
        seqBufQer: &[u8],
        aln: &mut [kswr_t],
        mat: &[i8; 25],
        q_max: u8,
        phase: i32,
    ) {
        unsafe {
            for &sp in lanes {
                let ind = sp.regid as usize;
                let target_start = sp.idr as usize;
                let query_start = sp.idq as usize;
                let target_len = sp.len1 as usize;
                let query_len = sp.len2 as usize;
                let target =
                    std::slice::from_raw_parts(seqBufRef.as_ptr().add(target_start), target_len);
                let query =
                    std::slice::from_raw_parts(seqBufQer.as_ptr().add(query_start), query_len);
                let ks = ksw_i16_slices(
                    query, self.m, mat, q_max, sp.len1, target, self.o_del, self.e_del, self.o_ins,
                    self.e_ins, sp.h0,
                );
                let a = aln.get_unchecked_mut(ind);
                if phase != 0 {
                    if a.score == ks.score {
                        a.tb = a.te - ks.te;
                        a.qb = a.qe - ks.qe;
                    }
                } else {
                    a.score = ks.score;
                    a.te = ks.te;
                    a.qe = ks.qe;
                    a.score2 = ks.score2;
                    a.te2 = ks.te2;
                }
            }
        }
    }

    #[inline]
    fn kswv_scalar_lanes_u8(
        &self,
        lanes: &[SeqPair],
        seqBufRef: &[u8],
        seqBufQer: &[u8],
        aln: &mut [kswr_t],
        mat: &[i8; 25],
        q_max: u8,
        phase: i32,
    ) {
        for &sp in lanes {
            let ind = sp.regid as usize;
            let target_start = sp.idr as usize;
            let query_start = sp.idq as usize;
            let target_len = sp.len1 as usize;
            let query_len = sp.len2 as usize;
            let target = &seqBufRef[target_start..target_start + target_len];
            let query = &seqBufQer[query_start..query_start + query_len];
            let ks = ksw_u8_slices(
                query, self.m, mat, q_max, sp.len1, target, self.o_del, self.e_del, self.o_ins,
                self.e_ins, sp.h0,
            );
            if phase != 0 {
                if aln[ind].score == ks.score {
                    aln[ind].tb = aln[ind].te - ks.te;
                    aln[ind].qb = aln[ind].qe - ks.qe;
                }
            } else {
                aln[ind].score = ks.score.min(255);
                aln[ind].te = ks.te;
                aln[ind].qe = ks.qe;
                aln[ind].score2 = ks.score2;
                aln[ind].te2 = ks.te2;
            }
        }
    }

    // *********************************** Vectorized Code 16 bit *****************************
    // 16 bit lanes. Scalar fallback path used when AVX-512BW is unavailable; the AVX-512
    // implementation lives in `kswv512_16_avx512_impl` below.
    #[doc = "Original function: kswv::kswv512_16:933"]
    #[inline]
    pub fn kswv512_16(
        &self,
        seq1SoA: &[i16],
        seq2SoA: &[i16],
        _nrow: i16,
        _ncol: i16,
        p: &[SeqPair],
        aln: &mut [kswr_t],
        po_ind: i32,
        _tid: u16,
        numPairs: i32,
        phase: i32,
    ) -> i32 {
        let mut mat = [0_i8; 25];
        self.bwa_fill_scmat(&mut mat);
        // ksw_qmax depends only on `mat`, not the per-lane query — hoist it.
        let q_max = ksw_qmax(self.m, &mat);
        // Hoist the per-lane decode buffers outside the loop so capacity reuses.
        let mut target_buf: Vec<u8> = Vec::new();
        let mut query_buf: Vec<u8> = Vec::new();
        let po_ind_us = po_ind as usize;
        let num_pairs_us = numPairs as usize;
        for lane in 0..SIMD_WIDTH16 {
            if po_ind_us + lane >= num_pairs_us {
                break;
            }
            let sp = p[lane];
            let ind = sp.regid as usize;
            decode_soa_lane_into(seq1SoA, lane, sp.len1 as usize, &mut target_buf);
            decode_soa_lane_into(seq2SoA, lane, sp.len2 as usize, &mut query_buf);
            let target = &target_buf[..];
            let query = &query_buf[..];
            // Forward-only SW (no internal reverse pass) — matches SIMD kernel semantics.
            // Borrow query/mat directly to avoid the two clones inside ksw_qinit + ksw_i16.
            let ks = ksw_i16_slices(
                query, self.m, &mat, q_max, sp.len1, target, self.o_del, self.e_del, self.o_ins,
                self.e_ins, sp.h0,
            );
            if phase != 0 {
                if aln[ind].score == ks.score {
                    aln[ind].tb = aln[ind].te - ks.te;
                    aln[ind].qb = aln[ind].qe - ks.qe;
                }
            } else {
                aln[ind].score = ks.score;
                aln[ind].te = ks.te;
                aln[ind].qe = ks.qe;
                aln[ind].score2 = ks.score2;
                aln[ind].te2 = ks.te2;
            }
        }
        1
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx512bw")]
    #[allow(clippy::too_many_arguments)]
    unsafe fn kswv512_16_avx512(
        &self,
        seq1_soa: &[i16],
        seq2_soa: &[i16],
        nrow: i16,
        ncol: i16,
        p: &[SeqPair],
        aln: &mut [kswr_t],
        po_ind: i32,
        numPairs: i32,
        phase: i32,
    ) -> i32 {
        KSWV16_AVX_SCRATCH.with(|cell| {
            let mut scratch = cell.borrow_mut();
            let (h0, h1, f_buf, row_max) = &mut *scratch;
            unsafe {
                self.kswv512_16_avx512_impl(
                    seq1_soa, seq2_soa, nrow, ncol, p, aln, po_ind, numPairs, phase, h0, h1, f_buf,
                    row_max,
                )
            }
        })
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx512bw")]
    #[allow(clippy::too_many_arguments)]
    unsafe fn kswv512_16_avx512_impl(
        &self,
        seq1_soa: &[i16],
        seq2_soa: &[i16],
        nrow: i16,
        ncol: i16,
        p: &[SeqPair],
        aln: &mut [kswr_t],
        po_ind: i32,
        numPairs: i32,
        phase: i32,
        h0: &mut Vec<i16>,
        h1: &mut Vec<i16>,
        f_buf: &mut Vec<i16>,
        row_max: &mut Vec<i16>,
    ) -> i32 {
        #[inline]
        unsafe fn load_i16x32(ptr: *const i16) -> __m512i {
            unsafe { _mm512_loadu_si512(ptr as *const __m512i) }
        }
        #[inline]
        unsafe fn store_i16x32(ptr: *mut i16, v: __m512i) {
            unsafe { _mm512_storeu_si512(ptr as *mut __m512i, v) }
        }
        #[inline]
        unsafe fn cmpge_epi16_mask(a: __m512i, b: __m512i) -> u32 {
            unsafe { _mm512_cmpgt_epi16_mask(a, b) | _mm512_cmpeq_epi16_mask(a, b) }
        }

        let nrow_us = nrow.max(0) as usize;
        let ncol_us = ncol.max(0) as usize;
        let dp_len = (ncol_us + 1) * SIMD_WIDTH16;
        h0.clear();
        h0.resize(dp_len, 0);
        if h1.len() < dp_len {
            h1.resize(dp_len, 0);
        } else {
            h1.truncate(dp_len);
        }
        f_buf.clear();
        f_buf.resize(dp_len, 0);
        row_max.clear();
        row_max.resize(nrow_us.max(1) * SIMD_WIDTH16, -1);

        let zero = unsafe { _mm512_setzero_si512() };
        let one = unsafe { _mm512_set1_epi16(1) };
        let minus1 = unsafe { _mm512_set1_epi16(-1) };
        let e_del = unsafe { _mm512_set1_epi16(self.e_del as i16) };
        let oe_del = unsafe { _mm512_set1_epi16((self.o_del + self.e_del) as i16) };
        let e_ins = unsafe { _mm512_set1_epi16(self.e_ins as i16) };
        let oe_ins = unsafe { _mm512_set1_epi16((self.o_ins + self.e_ins) as i16) };

        let mut perm = [0_i16; SIMD_WIDTH16];
        perm[0] = i16::from(self.w_match);
        perm[1] = i16::from(self.w_mismatch);
        perm[2] = i16::from(self.w_mismatch);
        perm[3] = i16::from(self.w_mismatch);
        for item in &mut perm[12..16] {
            *item = i16::from(self.w_ambig);
        }
        for item in &mut perm[16..20] {
            *item = i16::from(self.w_ambig);
        }
        perm[31] = i16::from(self.w_ambig);
        let perm_v = unsafe { load_i16x32(perm.as_ptr()) };

        let mut minsc = [0_i16; SIMD_WIDTH16];
        let mut endsc = [0_i16; SIMD_WIDTH16];
        let mut minsc_mask_a = 0_u32;
        let mut endsc_mask_a = 0_u32;
        for lane in 0..SIMD_WIDTH16 {
            let xtra = p[lane].h0;
            let val = if xtra & crate::bwa_mem2::ksw::KSW_XSUBO != 0 {
                xtra & 0xffff
            } else {
                0x10000
            };
            if val <= i32::from(i16::MAX) {
                minsc[lane] = val as i16;
                minsc_mask_a |= 1_u32 << lane;
            }
            let val = if xtra & crate::bwa_mem2::ksw::KSW_XSTOP != 0 {
                xtra & 0xffff
            } else {
                0x10000
            };
            if val <= i32::from(i16::MAX) {
                endsc[lane] = val as i16;
                endsc_mask_a |= 1_u32 << lane;
            }
        }
        let minsc_v = unsafe { load_i16x32(minsc.as_ptr()) };
        let endsc_v = unsafe { load_i16x32(endsc.as_ptr()) };

        let mut gmax = zero;
        let mut te_v = unsafe { _mm512_set1_epi16(-1) };
        let mut qe_v = zero;
        let mut exit0 = u32::MAX;
        let mut max_v: __m512i;
        let mut imax_v;
        let mut pimax_v = zero;
        let mut mask_v = 0_u32;
        let mut minsc_mask = 0_u32;
        let mut i_v = zero;
        let mut limit = nrow_us;
        let mut rows_done = 0_usize;

        // DP loop (corresponds to MAIN_SAM_CODE16_OPT in kswv.cpp):
        //   m11 = h00 + sbt    (sbt = shuffled match/mismatch from xor(s1,s2))
        //   h11 = max(m11, e11, f11, 0)
        //   e11 = max(h11 - oe_ins, e11 - e_ins)
        //   f21 = max(h11 - oe_del, f11 - e_del)
        // imax / iqe track the row-max and its column.
        for i in 0..nrow_us {
            let s1 = unsafe { load_i16x32(seq1_soa.as_ptr().add(i * SIMD_WIDTH16)) };
            let mut e11 = zero;
            imax_v = zero;
            let mut iqe_v = unsafe { _mm512_set1_epi16(-1) };
            let mut l_v = zero;
            for j in 0..ncol_us {
                let h00 = unsafe { load_i16x32(h0.as_ptr().add(j * SIMD_WIDTH16)) };
                let s2 = unsafe { load_i16x32(seq2_soa.as_ptr().add(j * SIMD_WIDTH16)) };
                let f11 = unsafe { load_i16x32(f_buf.as_ptr().add((j + 1) * SIMD_WIDTH16)) };
                let xor_v = unsafe { _mm512_xor_si512(s1, s2) };
                let sbt = unsafe { _mm512_permutexvar_epi16(xor_v, perm_v) };
                // Detect AMBIG (i16 = -1) at i16 granularity then fuse the zero-blend with the
                // add: m11 = !invalid ? add(h00, sbt) : 0. Saves 1 SIMD op vs the prior
                // add → movepi8 → mask_blend_epi8 chain.
                let or_v = unsafe { _mm512_or_si512(s1, s2) };
                let invalid_i16: u32 = unsafe { _mm512_movepi16_mask(or_v) };
                let m11 = unsafe { _mm512_mask_add_epi16(zero, !invalid_i16, h00, sbt) };
                // Fold the max chain into a tree to halve the critical path: original was
                // `((m11 max e11) max f11) max 0` (chain 3 deep); now `(m11 max e11) max
                // (f11 max 0)` is two parallel maxes followed by one combiner (depth 2).
                let me = unsafe { _mm512_max_epi16(m11, e11) };
                let f0 = unsafe { _mm512_max_epi16(f11, zero) };
                let h11 = unsafe { _mm512_max_epi16(me, f0) };
                let cmp0 = unsafe { _mm512_cmpgt_epi16_mask(h11, imax_v) };
                imax_v = unsafe { _mm512_max_epi16(imax_v, h11) };
                iqe_v = unsafe { _mm512_mask_blend_epi16(cmp0, iqe_v, l_v) };
                let gap_e = unsafe { _mm512_sub_epi16(h11, oe_ins) };
                e11 = unsafe { _mm512_sub_epi16(e11, e_ins) };
                e11 = unsafe { _mm512_max_epi16(gap_e, e11) };
                let gap_d = unsafe { _mm512_sub_epi16(h11, oe_del) };
                let mut f21 = unsafe { _mm512_sub_epi16(f11, e_del) };
                f21 = unsafe { _mm512_max_epi16(gap_d, f21) };
                unsafe { store_i16x32(h1.as_mut_ptr().add((j + 1) * SIMD_WIDTH16), h11) };
                unsafe { store_i16x32(f_buf.as_mut_ptr().add((j + 1) * SIMD_WIDTH16), f21) };
                l_v = unsafe { _mm512_add_epi16(l_v, one) };
            }

            // Block I: store the previous row's max for the score2 search below.
            if i > 0 {
                let mut msk = unsafe { _mm512_cmpgt_epi16_mask(imax_v, pimax_v) } | mask_v;
                // Combined: stored = (!msk && minsc_mask && exit0) ? pimax : minus1. Saves 2
                // blends vs the original 3-blend chain. Note: minus1 is the "fill" sentinel here
                // (vs zero in the u8 wrapper).
                let keep = !msk & minsc_mask & exit0;
                let stored = unsafe { _mm512_mask_blend_epi16(keep, minus1, pimax_v) };
                unsafe { store_i16x32(row_max.as_mut_ptr().add((i - 1) * SIMD_WIDTH16), stored) };
                msk &= u32::MAX;
                mask_v = !msk;
            }
            pimax_v = imax_v;
            minsc_mask = unsafe { cmpge_epi16_mask(imax_v, minsc_v) } & minsc_mask_a;

            // Block II: update per-lane gmax/te/qe; lanes where gmax >= endsc are marked dead
            // (cleared from exit0) so they no longer contribute to subsequent rows.
            let cmp0 = unsafe { _mm512_cmpgt_epi16_mask(imax_v, gmax) } & exit0;
            gmax = unsafe { _mm512_mask_blend_epi16(cmp0, gmax, imax_v) };
            te_v = unsafe { _mm512_mask_blend_epi16(cmp0, te_v, i_v) };
            qe_v = unsafe { _mm512_mask_blend_epi16(cmp0, qe_v, iqe_v) };
            let stop_mask = unsafe { cmpge_epi16_mask(gmax, endsc_v) } & endsc_mask_a;
            exit0 &= !stop_mask;
            rows_done = i + 1;
            if exit0 == 0 {
                limit = i;
                break;
            }
            std::mem::swap(h0, h1);
            i_v = unsafe { _mm512_add_epi16(i_v, one) };
        }

        if rows_done > 0 {
            // Combined: stored = (!mask_v && minsc_mask && exit0) ? pimax : minus1. Saves 2 blends.
            let keep = !mask_v & minsc_mask & exit0;
            let stored = unsafe { _mm512_mask_blend_epi16(keep, minus1, pimax_v) };
            unsafe {
                store_i16x32(
                    row_max.as_mut_ptr().add((rows_done - 1) * SIMD_WIDTH16),
                    stored,
                )
            };
        }

        let mut score = [0_i16; SIMD_WIDTH16];
        let mut te = [0_i16; SIMD_WIDTH16];
        let mut qe = [0_i16; SIMD_WIDTH16];
        unsafe {
            store_i16x32(score.as_mut_ptr(), gmax);
            store_i16x32(te.as_mut_ptr(), te_v);
            store_i16x32(qe.as_mut_ptr(), qe_v);
        }

        for lane in 0..SIMD_WIDTH16 {
            if po_ind + lane as i32 >= numPairs {
                break;
            }
            let ind = p[lane].regid as usize;
            if phase != 0 {
                if aln[ind].score == i32::from(score[lane]) {
                    aln[ind].tb = aln[ind].te - i32::from(te[lane]);
                    aln[ind].qb = aln[ind].qe - i32::from(qe[lane]);
                }
            } else {
                aln[ind].score = i32::from(score[lane]);
                aln[ind].te = i32::from(te[lane]);
                aln[ind].qe = i32::from(qe[lane]);
            }
        }
        if phase != 0 {
            return 1;
        }

        // Score2 and te2: scan rowMax for the best score outside the [low, high] band around
        // te. The band radius is (score + qmax - 1) / qmax, mirroring C++ kswv.cpp:613.
        let qmax = self.g_qmax.max(1);
        let mut low = [0_i16; SIMD_WIDTH16];
        let mut high = [0_i16; SIMD_WIDTH16];
        let mut maxl = 0_i32;
        let mut minh = nrow_us as i32;
        for lane in 0..SIMD_WIDTH16 {
            let val = (i32::from(score[lane]) + qmax - 1) / qmax;
            let lo = i32::from(te[lane]) - val;
            let hi = i32::from(te[lane]) + val;
            low[lane] = lo as i16;
            high[lane] = hi as i16;
            maxl = maxl.max(lo);
            minh = minh.min(hi);
        }
        max_v = unsafe { _mm512_set1_epi16(-1) };
        te_v = unsafe { _mm512_set1_epi16(-1) };
        let low_v = unsafe { load_i16x32(low.as_ptr()) };
        let high_v = unsafe { load_i16x32(high.as_ptr()) };

        for row in 0..(maxl.max(0) as usize) {
            let row_i = unsafe { _mm512_set1_epi16(row as i16) };
            let rmax = unsafe { load_i16x32(row_max.as_ptr().add(row * SIMD_WIDTH16)) };
            let mask1 = unsafe { _mm512_cmpgt_epi16_mask(low_v, row_i) };
            let mask2 = unsafe { _mm512_cmpgt_epi16_mask(rmax, max_v) } & mask1;
            max_v = unsafe { _mm512_mask_blend_epi16(mask2, max_v, rmax) };
            te_v = unsafe { _mm512_mask_blend_epi16(mask2, te_v, row_i) };
        }

        // Bounded scan above the band: require row < rlen so we don't read past valid rowMax.
        // The `rlen` mask was added in upstream to plug a bug where padding rows could become
        // the score2 candidate.
        let mut rlen = [0_i16; SIMD_WIDTH16];
        for lane in 0..SIMD_WIDTH16 {
            rlen[lane] = p[lane].len1 as i16;
        }
        let rlen_v = unsafe { load_i16x32(rlen.as_ptr()) };
        let start = (minh + 1).max(0) as usize;
        for row in start..limit {
            let row_i = unsafe { _mm512_set1_epi16(row as i16) };
            let rmax = unsafe { load_i16x32(row_max.as_ptr().add(row * SIMD_WIDTH16)) };
            let mask1 = unsafe { _mm512_cmpgt_epi16_mask(row_i, high_v) };
            let mut mask2 = unsafe { _mm512_cmpgt_epi16_mask(rmax, max_v) } & mask1;
            mask2 &= unsafe { _mm512_cmpgt_epi16_mask(rlen_v, row_i) };
            max_v = unsafe { _mm512_mask_blend_epi16(mask2, max_v, rmax) };
            te_v = unsafe { _mm512_mask_blend_epi16(mask2, te_v, row_i) };
        }

        let mut score2 = [0_i16; SIMD_WIDTH16];
        let mut te2 = [0_i16; SIMD_WIDTH16];
        unsafe {
            store_i16x32(score2.as_mut_ptr(), max_v);
            store_i16x32(te2.as_mut_ptr(), te_v);
        }
        for lane in 0..SIMD_WIDTH16 {
            if po_ind + lane as i32 >= numPairs {
                break;
            }
            let ind = p[lane].regid as usize;
            aln[ind].score2 = i32::from(score2[lane]);
            aln[ind].te2 = i32::from(te2[lane]);
        }
        1
    }

    // Vectorized u8 SW kernel — entry point. Dispatches to the AVX-512BW path when available
    // (`kswv512_u8_avx512_impl`) and otherwise falls back to per-lane scalar Farrar SW via
    // `ksw_u8_slices`. C++ uses an unsigned byte to store scores with a 256-shift so the
    // SW recurrence works in unsigned saturated arithmetic; on saturation (score >= 255 - shift)
    // the lane is marked dead and the alignment result will be retried at i16.
    #[doc = "Original function: kswv::kswv512_u8:371"]
    #[inline]
    pub fn kswv512_u8(
        &self,
        seq1SoA: &[u8],
        seq2SoA: &[u8],
        _nrow: i16,
        _ncol: i16,
        p: &[SeqPair],
        aln: &mut [kswr_t],
        po_ind: i32,
        _tid: u16,
        numPairs: i32,
        phase: i32,
    ) -> i32 {
        #[cfg(target_arch = "x86_64")]
        {
            if std::is_x86_feature_detected!("avx512bw") && !disable_kswv8_avx512() {
                return unsafe {
                    self.kswv512_u8_avx512(
                        seq1SoA, seq2SoA, _nrow, _ncol, p, aln, po_ind, numPairs, phase,
                    )
                };
            }
        }
        let mut mat = [0_i8; 25];
        self.bwa_fill_scmat(&mut mat);
        let q_max = ksw_qmax(self.m, &mat);
        let mut target_buf: Vec<u8> = Vec::new();
        let mut query_buf: Vec<u8> = Vec::new();
        let po_ind_us = po_ind as usize;
        let num_pairs_us = numPairs as usize;
        for lane in 0..SIMD_WIDTH8 {
            if po_ind_us + lane >= num_pairs_us {
                break;
            }
            let sp = p[lane];
            let ind = sp.regid as usize;
            decode_soa_lane_u8_into(seq1SoA, lane, sp.len1 as usize, &mut target_buf);
            decode_soa_lane_u8_query_into(seq2SoA, lane, sp.len2 as usize, &mut query_buf);
            let target = &target_buf[..];
            let query = &query_buf[..];
            // Forward-only SW (no internal reverse pass) — matches SIMD kernel semantics.
            // Borrow query/mat directly to avoid the two clones inside ksw_qinit + ksw_u8.
            let ks = ksw_u8_slices(
                query, self.m, &mat, q_max, sp.len1, target, self.o_del, self.e_del, self.o_ins,
                self.e_ins, sp.h0,
            );
            if phase != 0 {
                if aln[ind].score == ks.score {
                    aln[ind].tb = aln[ind].te - ks.te;
                    aln[ind].qb = aln[ind].qe - ks.qe;
                }
            } else {
                aln[ind].score = ks.score.min(255);
                aln[ind].te = ks.te;
                aln[ind].qe = ks.qe;
                aln[ind].score2 = ks.score2;
                aln[ind].te2 = ks.te2;
            }
        }
        1
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx512bw")]
    #[allow(clippy::too_many_arguments)]
    unsafe fn kswv512_u8_avx512(
        &self,
        seq1_soa: &[u8],
        seq2_soa: &[u8],
        nrow: i16,
        ncol: i16,
        p: &[SeqPair],
        aln: &mut [kswr_t],
        po_ind: i32,
        numPairs: i32,
        phase: i32,
    ) -> i32 {
        KSWV8_AVX_SCRATCH.with(|cell| {
            let mut scratch = cell.borrow_mut();
            let (h0, h1, f_buf, row_max) = &mut *scratch;
            unsafe {
                self.kswv512_u8_avx512_impl(
                    seq1_soa, seq2_soa, nrow, ncol, p, aln, po_ind, numPairs, phase, h0, h1, f_buf,
                    row_max,
                )
            }
        })
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx512bw")]
    #[allow(clippy::too_many_arguments)]
    unsafe fn kswv512_u8_avx512_impl(
        &self,
        seq1_soa: &[u8],
        seq2_soa: &[u8],
        nrow: i16,
        ncol: i16,
        p: &[SeqPair],
        aln: &mut [kswr_t],
        po_ind: i32,
        numPairs: i32,
        phase: i32,
        h0: &mut Vec<u8>,
        h1: &mut Vec<u8>,
        f_buf: &mut Vec<u8>,
        row_max: &mut Vec<u8>,
    ) -> i32 {
        #[inline]
        unsafe fn load_u8x64(ptr: *const u8) -> __m512i {
            unsafe { _mm512_loadu_si512(ptr as *const __m512i) }
        }
        #[inline]
        unsafe fn store_u8x64(ptr: *mut u8, v: __m512i) {
            unsafe { _mm512_storeu_si512(ptr as *mut __m512i, v) }
        }
        #[inline]
        unsafe fn load_i16x32(ptr: *const i16) -> __m512i {
            unsafe { _mm512_loadu_si512(ptr as *const __m512i) }
        }
        #[inline]
        unsafe fn cmpge_epu8_mask(a: __m512i, b: __m512i) -> u64 {
            unsafe { _mm512_cmpgt_epu8_mask(a, b) | _mm512_cmpeq_epu8_mask(a, b) }
        }
        #[inline]
        unsafe fn cmpgt_epi16_mask(a: __m512i, b: __m512i) -> u32 {
            unsafe { _mm512_cmpgt_epi16_mask(a, b) }
        }

        let nrow_us = nrow.max(0) as usize;
        let ncol_us = ncol.max(0) as usize;
        let dp_len = (ncol_us + 1) * SIMD_WIDTH8;
        h0.clear();
        h0.resize(dp_len, 0);
        if h1.len() < dp_len {
            h1.resize(dp_len, 0);
        } else {
            h1.truncate(dp_len);
            h1.fill(0);
        }
        f_buf.clear();
        f_buf.resize(dp_len, 0);
        row_max.clear();
        row_max.resize(nrow_us.max(1) * SIMD_WIDTH8, 0);

        let zero = unsafe { _mm512_setzero_si512() };
        let one = unsafe { _mm512_set1_epi8(1) };
        let mut shift = self.w_match.min(self.w_mismatch).min(self.w_ambig) as u8;
        let mdiff = self.w_match.max(self.w_mismatch).max(self.w_ambig) as u8;
        let qmax = mdiff;
        shift = 0_u8.wrapping_sub(shift);
        let sft = unsafe { _mm512_set1_epi8(shift as i8) };
        let cmax = unsafe { _mm512_set1_epi8(-1) };

        let mut perm = [0_i8; SIMD_WIDTH8];
        perm[0] = self.w_match;
        perm[1] = self.w_mismatch;
        perm[2] = self.w_mismatch;
        perm[3] = self.w_mismatch;
        for item in &mut perm[4..13] {
            *item = self.w_ambig;
        }
        for item in perm.iter_mut().take(16) {
            *item = item.wrapping_add(shift as i8);
        }
        for idx in 16..SIMD_WIDTH8 {
            perm[idx] = perm[(idx - 16) & 15];
        }
        let perm_v = unsafe { load_u8x64(perm.as_ptr() as *const u8) };

        let mut minsc = [0_u8; SIMD_WIDTH8];
        let mut endsc = [0_u8; SIMD_WIDTH8];
        let mut minsc_mask_a = 0_u64;
        let mut endsc_mask_a = 0_u64;
        for lane in 0..SIMD_WIDTH8 {
            let xtra = p[lane].h0;
            let val = if xtra & crate::bwa_mem2::ksw::KSW_XSUBO != 0 {
                xtra & 0xffff
            } else {
                0x10000
            };
            if val <= u8::MAX as i32 {
                minsc[lane] = val as u8;
                minsc_mask_a |= 1_u64 << lane;
            }
            let val = if xtra & crate::bwa_mem2::ksw::KSW_XSTOP != 0 {
                xtra & 0xffff
            } else {
                0x10000
            };
            if val <= u8::MAX as i32 {
                endsc[lane] = val as u8;
                endsc_mask_a |= 1_u64 << lane;
            }
        }
        let minsc_v = unsafe { load_u8x64(minsc.as_ptr()) };
        let endsc_v = unsafe { load_u8x64(endsc.as_ptr()) };
        let e_del = unsafe { _mm512_set1_epi8(self.e_del as i8) };
        let oe_del = unsafe { _mm512_set1_epi8((self.o_del + self.e_del) as i8) };
        let e_ins = unsafe { _mm512_set1_epi8(self.e_ins as i8) };
        let oe_ins = unsafe { _mm512_set1_epi8((self.o_ins + self.e_ins) as i8) };
        let five = unsafe { _mm512_set1_epi8(5) };

        let mut gmax = zero;
        let mut te_lo = unsafe { _mm512_set1_epi16(-1) };
        let mut te_hi = unsafe { _mm512_set1_epi16(-1) };
        let mut qe_v = zero;
        let mut exit0 = u64::MAX;
        let mut pimax = zero;
        let mut mask = 0_u64;
        let mut minsc_mask = 0_u64;
        let mut limit = nrow_us;
        let mut rows_done = 0_usize;

        // DP loop (corresponds to MAIN_SAM_CODE8_OPT in kswv.cpp). Like the i16 variant,
        // but in saturated unsigned u8 arithmetic with a `shift` bias so subtraction
        // doesn't underflow; gmax + shift >= 255 triggers per-lane saturation exit.
        for i in 0..nrow_us {
            let s1 = unsafe { load_u8x64(seq1_soa.as_ptr().add(i * SIMD_WIDTH8)) };
            let mut e11 = zero;
            let mut imax = zero;
            let mut iqe = unsafe { _mm512_set1_epi8(-1) };
            let mut l_v = zero;
            for j in 0..ncol_us {
                let h00 = unsafe { load_u8x64(h0.as_ptr().add(j * SIMD_WIDTH8)) };
                let s2 = unsafe { load_u8x64(seq2_soa.as_ptr().add(j * SIMD_WIDTH8)) };
                let f11 = unsafe { load_u8x64(f_buf.as_ptr().add((j + 1) * SIMD_WIDTH8)) };
                let xor_v = unsafe { _mm512_xor_si512(s1, s2) };
                let mut sbt = unsafe { _mm512_shuffle_epi8(perm_v, xor_v) };
                let cmpq = unsafe { _mm512_cmpeq_epu8_mask(s2, five) };
                sbt = unsafe { _mm512_mask_blend_epi8(cmpq, sbt, sft) };
                let or_v = unsafe { _mm512_or_si512(s1, s2) };
                let invalid = unsafe { _mm512_movepi8_mask(or_v) };
                let m11_raw = unsafe { _mm512_adds_epu8(h00, sbt) };
                // Fused: m11 = (!invalid) ? subs_epu8(m11_raw, sft) : 0. Saves the separate
                // blend(invalid, m11, zero) before the subs.
                let m11 = unsafe { _mm512_mask_subs_epu8(zero, !invalid, m11_raw, sft) };
                let mut h11 = unsafe { _mm512_max_epu8(m11, e11) };
                h11 = unsafe { _mm512_max_epu8(h11, f11) };
                let cmp0 = unsafe { _mm512_cmpgt_epu8_mask(h11, imax) };
                imax = unsafe { _mm512_max_epu8(imax, h11) };
                iqe = unsafe { _mm512_mask_blend_epi8(cmp0, iqe, l_v) };
                let gap_e = unsafe { _mm512_subs_epu8(h11, oe_ins) };
                e11 = unsafe { _mm512_subs_epu8(e11, e_ins) };
                e11 = unsafe { _mm512_max_epu8(gap_e, e11) };
                let gap_d = unsafe { _mm512_subs_epu8(h11, oe_del) };
                let mut f21 = unsafe { _mm512_subs_epu8(f11, e_del) };
                f21 = unsafe { _mm512_max_epu8(gap_d, f21) };
                unsafe { store_u8x64(h1.as_mut_ptr().add((j + 1) * SIMD_WIDTH8), h11) };
                unsafe { store_u8x64(f_buf.as_mut_ptr().add((j + 1) * SIMD_WIDTH8), f21) };
                l_v = unsafe { _mm512_add_epi8(l_v, one) };
            }

            if i > 0 {
                let msk = unsafe { _mm512_cmpgt_epu8_mask(imax, pimax) } | mask;
                // Combined: pimax_out = (!msk && minsc_mask && exit0) ? pimax : 0. Saves 2 blends
                // vs the original 3-blend chain.
                let keep = !msk & minsc_mask & exit0;
                let stored = unsafe { _mm512_mask_blend_epi8(keep, zero, pimax) };
                unsafe { store_u8x64(row_max.as_mut_ptr().add((i - 1) * SIMD_WIDTH8), stored) };
                mask = !msk;
            }
            pimax = imax;
            minsc_mask = unsafe { cmpge_epu8_mask(imax, minsc_v) } & minsc_mask_a;

            let cmp0 = unsafe { _mm512_cmpgt_epu8_mask(imax, gmax) } & exit0;
            gmax = unsafe { _mm512_mask_blend_epi8(cmp0, gmax, imax) };
            let i_v = unsafe { _mm512_set1_epi16(i as i16) };
            te_lo = unsafe { _mm512_mask_blend_epi16(cmp0 as u32, te_lo, i_v) };
            te_hi = unsafe { _mm512_mask_blend_epi16((cmp0 >> SIMD_WIDTH16) as u32, te_hi, i_v) };
            qe_v = unsafe { _mm512_mask_blend_epi8(cmp0, qe_v, iqe) };

            let stop_mask = unsafe { cmpge_epu8_mask(gmax, endsc_v) } & endsc_mask_a;
            let left = unsafe { _mm512_adds_epu8(gmax, sft) };
            let sat_mask = unsafe { cmpge_epu8_mask(left, cmax) };
            exit0 &= !(stop_mask | sat_mask);
            rows_done = i + 1;
            if exit0 == 0 {
                limit = i;
                break;
            }
            std::mem::swap(h0, h1);
        }

        if rows_done > 0 {
            // Combined: stored = (!mask && minsc_mask && exit0) ? pimax : 0. Saves 2 blends.
            let keep = !mask & minsc_mask & exit0;
            let stored = unsafe { _mm512_mask_blend_epi8(keep, zero, pimax) };
            unsafe {
                store_u8x64(
                    row_max.as_mut_ptr().add((rows_done - 1) * SIMD_WIDTH8),
                    stored,
                )
            };
        }

        let mut score = [0_u8; SIMD_WIDTH8];
        let mut te = [0_i16; SIMD_WIDTH8];
        let mut qe = [0_u8; SIMD_WIDTH8];
        unsafe {
            store_u8x64(score.as_mut_ptr(), gmax);
            _mm512_storeu_si512(te.as_mut_ptr() as *mut __m512i, te_lo);
            _mm512_storeu_si512(te.as_mut_ptr().add(SIMD_WIDTH16) as *mut __m512i, te_hi);
            store_u8x64(qe.as_mut_ptr(), qe_v);
        }

        let mut live = 0_i32;
        for lane in 0..SIMD_WIDTH8 {
            if po_ind + lane as i32 >= numPairs {
                break;
            }
            let ind = p[lane].regid as usize;
            if phase != 0 {
                if aln[ind].score == i32::from(score[lane]) {
                    aln[ind].tb = aln[ind].te - i32::from(te[lane]);
                    aln[ind].qb = aln[ind].qe - i32::from(qe[lane]);
                }
            } else {
                aln[ind].score = i32::from(score[lane]);
                aln[ind].te = i32::from(te[lane]);
                aln[ind].qe = i32::from(qe[lane]);
                if score[lane] != u8::MAX {
                    qe[lane] = 1;
                    live += 1;
                } else {
                    qe[lane] = 0;
                }
            }
        }
        if phase != 0 || live == 0 {
            return 1;
        }

        let qmax = i32::from(qmax).max(1);
        let mut low = [0_i16; SIMD_WIDTH8];
        let mut high = [0_i16; SIMD_WIDTH8];
        let mut maxl = 0_i32;
        let mut minh = nrow_us as i32;
        for lane in 0..SIMD_WIDTH8 {
            let val = (i32::from(score[lane]) + qmax - 1) / qmax;
            low[lane] = (i32::from(te[lane]) - val) as i16;
            high[lane] = (i32::from(te[lane]) + val) as i16;
            if qe[lane] != 0 {
                maxl = maxl.max(i32::from(low[lane]));
                minh = minh.min(i32::from(high[lane]));
            }
        }

        let mut max_v = zero;
        te_lo = unsafe { _mm512_set1_epi16(-1) };
        te_hi = unsafe { _mm512_set1_epi16(-1) };
        let low_lo = unsafe { load_i16x32(low.as_ptr()) };
        let high_lo = unsafe { load_i16x32(high.as_ptr()) };
        let low_hi = unsafe { load_i16x32(low.as_ptr().add(SIMD_WIDTH16)) };
        let high_hi = unsafe { load_i16x32(high.as_ptr().add(SIMD_WIDTH16)) };

        for row in 0..(maxl.max(0) as usize) {
            let row_i = unsafe { _mm512_set1_epi16(row as i16) };
            let rmax = unsafe { load_u8x64(row_max.as_ptr().add(row * SIMD_WIDTH8)) };
            let mask_lo = unsafe { cmpgt_epi16_mask(low_lo, row_i) };
            let mask_hi = unsafe { cmpgt_epi16_mask(low_hi, row_i) };
            let mut mask2 = unsafe { _mm512_cmpgt_epu8_mask(rmax, max_v) };
            let mask1 = u64::from(mask_lo) | (u64::from(mask_hi) << SIMD_WIDTH16);
            mask2 &= mask1;
            max_v = unsafe { _mm512_mask_blend_epi8(mask2, max_v, rmax) };
            te_lo = unsafe { _mm512_mask_blend_epi16(mask2 as u32, te_lo, row_i) };
            te_hi =
                unsafe { _mm512_mask_blend_epi16((mask2 >> SIMD_WIDTH16) as u32, te_hi, row_i) };
        }

        let mut rlen = [0_i16; SIMD_WIDTH8];
        for lane in 0..SIMD_WIDTH8 {
            rlen[lane] = p[lane].len1 as i16;
        }
        let rlen_lo = unsafe { load_i16x32(rlen.as_ptr()) };
        let rlen_hi = unsafe { load_i16x32(rlen.as_ptr().add(SIMD_WIDTH16)) };
        let start = (minh + 1).max(0) as usize;
        for row in start..limit {
            let row_i = unsafe { _mm512_set1_epi16(row as i16) };
            let rmax = unsafe { load_u8x64(row_max.as_ptr().add(row * SIMD_WIDTH8)) };
            let mask_lo = unsafe { cmpgt_epi16_mask(row_i, high_lo) };
            let mask_hi = unsafe { cmpgt_epi16_mask(row_i, high_hi) };
            let mut mask2 = unsafe { _mm512_cmpgt_epu8_mask(rmax, max_v) };
            let mask1 = u64::from(mask_lo) | (u64::from(mask_hi) << SIMD_WIDTH16);
            let rmask_lo = unsafe { cmpgt_epi16_mask(rlen_lo, row_i) };
            let rmask_hi = unsafe { cmpgt_epi16_mask(rlen_hi, row_i) };
            let rmask = u64::from(rmask_lo) | (u64::from(rmask_hi) << SIMD_WIDTH16);
            mask2 &= mask1 & rmask;
            max_v = unsafe { _mm512_mask_blend_epi8(mask2, max_v, rmax) };
            te_lo = unsafe { _mm512_mask_blend_epi16(mask2 as u32, te_lo, row_i) };
            te_hi =
                unsafe { _mm512_mask_blend_epi16((mask2 >> SIMD_WIDTH16) as u32, te_hi, row_i) };
        }

        let mut score2 = [0_u8; SIMD_WIDTH8];
        let mut te2 = [0_i16; SIMD_WIDTH8];
        unsafe {
            store_u8x64(score2.as_mut_ptr(), max_v);
            _mm512_storeu_si512(te2.as_mut_ptr() as *mut __m512i, te_lo);
            _mm512_storeu_si512(te2.as_mut_ptr().add(SIMD_WIDTH16) as *mut __m512i, te_hi);
        }
        for lane in 0..SIMD_WIDTH8 {
            if po_ind + lane as i32 >= numPairs {
                break;
            }
            let ind = p[lane].regid as usize;
            if qe[lane] != 0 {
                aln[ind].score2 = if score2[lane] == 0 {
                    -1
                } else {
                    i32::from(score2[lane])
                };
                aln[ind].te2 = i32::from(te2[lane]);
            } else {
                aln[ind].score2 = -1;
                aln[ind].te2 = -1;
            }
        }
        1
    }

    // **************************************Scalar code***************************************
    // This is the original SW code from bwa-mem. We are keeping both, 8-bit and 16-bit
    // implementations, here for benchmarking purposes. The interface to the code is very
    // simple and similar to the one we used above. By default the C++ build disables this.
    #[doc = "Original function: kswv::bwa_fill_scmat:1231"]
    #[inline]
    pub fn bwa_fill_scmat(&self, mat: &mut [i8; 25]) {
        let mut k = 0_usize;
        for i in 0..4 {
            for j in 0..4 {
                mat[k] = if i == j {
                    self.w_match
                } else {
                    self.w_mismatch
                };
                k += 1;
            }
            mat[k] = self.w_ambig; // ambiguous base
            k += 1;
        }
        for item in mat.iter_mut().skip(k).take(5) {
            *item = self.w_ambig;
        }
    }

    /// Initialize the query data structure.
    ///
    /// * `size`  — Number of bytes used to store a score; valid values are 1 or 2.
    /// * `qlen`  — Length of the query sequence.
    /// * `query` — Query sequence.
    /// * `m`     — Size of the alphabet.
    /// * `mat`   — Scoring matrix in a one-dimensional array.
    ///
    /// Returns the query data structure.
    ///
    /// An example: p=8, qlen=19, slen=3 and segmentation:
    ///   {{0,3,6,9,12,15,18,-1},{1,4,7,10,13,16,-1,-1},{2,5,8,11,14,17,-1,-1}}
    #[doc = "Original function: kswv::ksw_qinit:1256"]
    pub fn ksw_qinit(&self, size: i32, qlen: i32, query: &[u8], m: i32, mat: &[i8]) -> kswq_t {
        let size = if size > 1 { 2 } else { 1 };
        let p = 8 * (3 - size); // # values per __m128i
        let slen = (qlen + p - 1) / p; // segmented length
        let qlen_usize = usize::try_from(qlen).expect("qlen");
        let slen_usize = usize::try_from(slen).expect("slen");
        let p_usize = usize::try_from(p).expect("p");
        let m_usize = usize::try_from(m).expect("m");

        // compute shift: find the minimum and maximum score
        let mut shift = i8::MAX;
        let mut mdiff = 0_i8;
        for &score in mat.iter().take(m_usize * m_usize) {
            if score < shift {
                shift = score;
            }
            if score > mdiff {
                mdiff = score;
            }
        }
        let max = mdiff;
        // NB: shift is uint8_t (wraparound bias = 256 - min)
        let shift_u8 = (0_i16).wrapping_sub(i16::from(shift)) as u8;
        // difference between the min and max scores after biasing
        let mdiff_u8 = mdiff.wrapping_add(shift_u8 as i8) as u8;

        let mut qp_u8 = Vec::new();
        let mut qp_i16 = Vec::new();
        if size == 1 {
            qp_u8.reserve(slen_usize * p_usize * m_usize);
            for a in 0..m_usize {
                let ma = &mat[a * m_usize..(a + 1) * m_usize];
                let nlen = slen_usize * p_usize;
                for i in 0..slen_usize {
                    let mut k = i;
                    while k < nlen {
                        let v = if k >= qlen_usize {
                            0
                        } else {
                            i16::from(ma[usize::from(query[k])]) + i16::from(shift_u8)
                        };
                        qp_u8.push(i8::try_from(v).expect("qp_u8"));
                        k += slen_usize;
                    }
                }
            }
        } else {
            qp_i16.reserve(slen_usize * p_usize * m_usize);
            for a in 0..m_usize {
                let ma = &mat[a * m_usize..(a + 1) * m_usize];
                let nlen = slen_usize * p_usize;
                for i in 0..slen_usize {
                    let mut k = i;
                    while k < nlen {
                        let v = if k >= qlen_usize {
                            0
                        } else {
                            i16::from(ma[usize::from(query[k])])
                        };
                        qp_i16.push(v);
                        k += slen_usize;
                    }
                }
            }
        }

        kswq_t {
            qlen,
            slen,
            shift: shift_u8,
            mdiff: mdiff_u8,
            max: u8::try_from(max).expect("max"),
            size: u8::try_from(size).expect("size"),
            qp_u8,
            qp_i16,
        }
    }

    // Scalar u8 SW (the first gap costs -(_o+_e)). Original SSE2 SW kernel from bwa-mem,
    // here adapted to call into ksw_align2 for the scalar/scalar-replacement path.
    #[doc = "Original function: kswv::kswvScalar_u8:1306"]
    pub fn kswvScalar_u8(
        &self,
        q: &kswq_t,
        tlen: i32,
        target: &[u8],
        _o_del: i32,
        _e_del: i32,
        _o_ins: i32,
        _e_ins: i32,
        xtra: i32,
    ) -> kswr_t {
        let mut mat = [0_i8; 25];
        self.bwa_fill_scmat(&mut mat);
        let mut query = reconstruct_query_u8(q);
        let mut target = target[..usize::try_from(tlen).expect("tlen")].to_vec();
        let ks = ksw_align2(
            q.qlen,
            &mut query,
            tlen,
            &mut target,
            self.m,
            &mat,
            self.o_del,
            self.e_del,
            self.o_ins,
            self.e_ins,
            xtra | KSW_XSTART,
            None,
        );
        kswr_t {
            score: ks.score,
            te: ks.te,
            qe: ks.qe,
            score2: ks.score2,
            te2: ks.te2,
            tb: ks.tb,
            qb: ks.qb,
        }
    }

    // Scalar i16 SW (the first gap costs -(_o+_e)). Same role as kswvScalar_u8 but i16 lanes.
    #[doc = "Original function: kswv::kswvScalar_i16:1434"]
    pub fn kswvScalar_i16(
        &self,
        q: &kswq_t,
        tlen: i32,
        target: &[u8],
        _o_del: i32,
        _e_del: i32,
        _o_ins: i32,
        _e_ins: i32,
        xtra: i32,
    ) -> kswr_t {
        let mut mat = [0_i8; 25];
        self.bwa_fill_scmat(&mut mat);
        let mut query = reconstruct_query_i16(q, &mat, self.m);
        let mut target = target[..usize::try_from(tlen).expect("tlen")].to_vec();
        let ks = ksw_align2(
            q.qlen,
            &mut query,
            tlen,
            &mut target,
            self.m,
            &mat,
            self.o_del,
            self.e_del,
            self.o_ins,
            self.e_ins,
            xtra | KSW_XSTART,
            None,
        );
        kswr_t {
            score: ks.score,
            te: ks.te,
            qe: ks.qe,
            score2: ks.score2,
            te2: ks.te2,
            tb: ks.tb,
            qb: ks.qb,
        }
    }

    // -------------------------------------------------------------
    // kswc scalar, wrapper function, the interface. Iterates over the SeqPair batch and
    // dispatches each pair to kswvScalar_u8 / kswvScalar_i16 (selected by `sw` and the
    // KSW_XBYTE flag in p->h0).
    // -------------------------------------------------------------
    #[doc = "Original function: kswv::kswvScalarWrapper:1550"]
    pub fn kswvScalarWrapper(
        &self,
        seqPairArray: &mut [SeqPair],
        seqBufRef: &[u8],
        seqBufQer: &[u8],
        aln: &mut [kswr_t],
        numPairs: i32,
        _nthreads: i32,
        sw: bool,
        _tid: i32,
    ) {
        let mut mat = [0_i8; 25];
        self.bwa_fill_scmat(&mut mat);
        for p in seqPairArray
            .iter()
            .take(usize::try_from(numPairs).expect("numPairs"))
        {
            let target_start = usize::try_from(p.idr).expect("idr");
            let query_start = usize::try_from(p.idq).expect("idq");
            let tlen = usize::try_from(p.len1).expect("len1");
            let qlen = usize::try_from(p.len2).expect("len2");
            let target = &seqBufRef[target_start..target_start + tlen];
            let query = &seqBufQer[query_start..query_start + qlen];
            let q = self.ksw_qinit(
                if sw {
                    2
                } else if (p.h0 & KSW_XBYTE) != 0 {
                    1
                } else {
                    2
                },
                p.len2,
                query,
                self.m,
                &mat,
            );
            let ks = if sw {
                self.kswvScalar_i16(
                    &q, p.len1, target, self.o_del, self.e_del, self.o_ins, self.e_ins, p.h0,
                )
            } else {
                self.kswvScalar_u8(
                    &q, p.len1, target, self.o_del, self.e_del, self.o_ins, self.e_ins, p.h0,
                )
            };
            let slot = &mut aln[usize::try_from(p.regid).expect("regid")];
            *slot = ks;
        }
    }
}

fn reconstruct_query_u8(q: &kswq_t) -> Vec<u8> {
    let qlen = usize::try_from(q.qlen).expect("qlen");
    let slen = usize::try_from(q.slen).expect("slen");
    let p = 16_usize;
    let mut query = vec![4_u8; qlen];
    for pos in 0..qlen {
        let idx = linearized_index(pos, slen);
        let mut best_base = 4_u8;
        let mut best_score = i8::MIN;
        for base in 0..4_usize {
            let v = q.qp_u8[base * slen * p + idx];
            if v > best_score {
                best_score = v;
                best_base = u8::try_from(base).expect("base");
            }
        }
        query[pos] = best_base;
    }
    query
}

fn reconstruct_query_i16(q: &kswq_t, mat: &[i8; 25], m: i32) -> Vec<u8> {
    let qlen = usize::try_from(q.qlen).expect("qlen");
    let slen = usize::try_from(q.slen).expect("slen");
    let p = 8_usize;
    let m_usize = usize::try_from(m).expect("m");
    let mut query = vec![4_u8; qlen];
    for pos in 0..qlen {
        for base in 0..4_usize {
            let stored = q.qp_i16[base * slen * p + linearized_index_i16(pos, slen)];
            let expect = i16::from(mat[base * m_usize + base]);
            if stored == expect {
                query[pos] = u8::try_from(base).expect("base");
                break;
            }
        }
    }
    query
}

fn linearized_index(pos: usize, slen: usize) -> usize {
    let lane = pos / slen;
    let i = pos % slen;
    i + lane * slen
}

fn linearized_index_i16(pos: usize, slen: usize) -> usize {
    let lane = pos / slen;
    let i = pos % slen;
    i + lane * slen
}

fn decode_soa_lane(soa: &[i16], lane: usize, len: usize) -> Vec<u8> {
    let mut seq = Vec::with_capacity(len);
    decode_soa_lane_into(soa, lane, len, &mut seq);
    seq
}

#[inline]
fn decode_soa_lane_into(soa: &[i16], lane: usize, len: usize, out: &mut Vec<u8>) {
    out.clear();
    out.reserve(len);
    for row in 0..len {
        let v = soa[row * SIMD_WIDTH16 + lane];
        let base = if (0..=3).contains(&v) { v as u8 } else { 4 };
        out.push(base);
    }
}

fn decode_soa_lane_u8(soa: &[u8], lane: usize, len: usize) -> Vec<u8> {
    let mut seq = Vec::with_capacity(len);
    decode_soa_lane_u8_into(soa, lane, len, &mut seq);
    seq
}

#[inline]
fn decode_soa_lane_u8_into(soa: &[u8], lane: usize, len: usize, out: &mut Vec<u8>) {
    out.clear();
    out.reserve(len);
    for row in 0..len {
        let v = soa[row * SIMD_WIDTH8 + lane];
        out.push(if v <= 3 { v } else { 4 });
    }
}

fn decode_soa_lane_u8_query(soa: &[u8], lane: usize, len: usize) -> Vec<u8> {
    let mut seq = Vec::with_capacity(len);
    decode_soa_lane_u8_query_into(soa, lane, len, &mut seq);
    seq
}

#[inline]
fn decode_soa_lane_u8_query_into(soa: &[u8], lane: usize, len: usize, out: &mut Vec<u8>) {
    out.clear();
    out.reserve(len);
    for row in 0..len {
        let v = soa[row * SIMD_WIDTH8 + lane];
        out.push(if v <= 3 { v } else { 4 });
    }
}

#[cfg(test)]
mod tests {
    use super::{
        decode_soa_lane, decode_soa_lane_u8, decode_soa_lane_u8_query, kswr_t, kswv, SeqPair,
        SIMD_WIDTH16, SIMD_WIDTH8,
    };

    fn score_mat() -> [i8; 25] {
        let mut mat = [0_i8; 25];
        for i in 0..4 {
            for j in 0..4 {
                mat[i * 5 + j] = if i == j { 1 } else { -4 };
            }
            mat[i * 5 + 4] = -1;
        }
        for j in 0..5 {
            mat[20 + j] = -1;
        }
        mat
    }

    #[test]
    fn ctor_allocates_expected_buffers() {
        let k = kswv::ctor(6, 1, 6, 1, 1, -4, 2, Some(32), Some(16));
        assert_eq!(k.m, 5);
        assert_eq!(k.maxRefLen, 48);
        assert_eq!(k.maxQerLen, 32);
        assert!(!k.F16.is_empty());
        assert!(!k.H16_0.is_empty());
        assert!(!k.rowMax16.is_empty());
    }

    #[test]
    fn ksw_qinit_builds_segmented_profile_shapes() {
        let k = kswv::ctor(6, 1, 6, 1, 1, -4, 1, None, None);
        let mat = score_mat();
        let q = k.ksw_qinit(1, 5, &[0, 1, 2, 3, 0], 5, &mat);
        assert_eq!(q.qlen, 5);
        assert_eq!(q.size, 1);
        assert!(!q.qp_u8.is_empty());
    }

    #[test]
    fn getScores8_matches_scalar_wrapper_outputs() {
        let k = kswv::ctor(6, 1, 6, 1, 1, -4, 1, None, None);
        let mut pairs = vec![SeqPair {
            idr: 0,
            idq: 0,
            len1: 4,
            len2: 4,
            h0: 0x10000 | 0x80000,
            regid: 0,
            ..Default::default()
        }];
        let seq_ref = vec![0_u8, 1, 2, 3];
        let seq_qer = vec![0_u8, 1, 2, 3];
        let mut aln = vec![kswr_t::default(); 1];
        k.getScores8(&mut pairs, &seq_ref, &seq_qer, &mut aln, 1, 1, 0);
        assert!(aln[0].score > 0);
        assert_eq!(aln[0].te, 3);
        assert_eq!(aln[0].qe, 3);
    }

    #[test]
    fn getScores8_batch_wrapper_matches_scalar_lanes_on_full_and_tail_chunks() {
        let k = kswv::ctor(6, 1, 6, 1, 1, -4, 1, Some(64), Some(64));
        let mut pairs = Vec::new();
        let mut seq_ref = Vec::new();
        let mut seq_qer = Vec::new();
        for lane in 0..(SIMD_WIDTH8 + 7) {
            let len1 = 9 + (lane % 11) as i32;
            let len2 = 8 + (lane % 9) as i32;
            let idr = i32::try_from(seq_ref.len()).expect("idr");
            let idq = i32::try_from(seq_qer.len()).expect("idq");
            for i in 0..usize::try_from(len1).expect("len1") {
                seq_ref.push(((i + lane) & 3) as u8);
            }
            for i in 0..usize::try_from(len2).expect("len2") {
                seq_qer.push(((i + lane + usize::from(lane % 5 == 0)) & 3) as u8);
            }
            pairs.push(SeqPair {
                idr,
                idq,
                len1,
                len2,
                h0: if lane % 3 == 0 { 0x10000 | 3 } else { 0x80000 },
                regid: lane as i32,
                ..Default::default()
            });
        }

        let mut mat = [0_i8; 25];
        k.bwa_fill_scmat(&mut mat);
        let q_max = crate::bwa_mem2::ksw::ksw_qmax(k.m, &mat);
        let mut scalar_aln = vec![kswr_t::default(); pairs.len()];
        for chunk in pairs.chunks(SIMD_WIDTH8) {
            k.kswv_scalar_lanes_u8(chunk, &seq_ref, &seq_qer, &mut scalar_aln, &mat, q_max, 0);
        }

        let mut batch_pairs = pairs.clone();
        let mut batch_aln = vec![kswr_t::default(); pairs.len()];
        k.getScores8(
            &mut batch_pairs,
            &seq_ref,
            &seq_qer,
            &mut batch_aln,
            i32::try_from(pairs.len()).expect("numPairs"),
            1,
            0,
        );
        for (lane, (batch, scalar)) in batch_aln.iter().zip(&scalar_aln).enumerate() {
            assert_eq!(batch.score, scalar.score, "score lane {lane}");
            assert_eq!(batch.te, scalar.te, "te lane {lane}");
            assert_eq!(batch.qe, scalar.qe, "qe lane {lane}");
            assert_eq!(batch.score2, scalar.score2, "score2 lane {lane}");
            assert_eq!(batch.te2, scalar.te2, "te2 lane {lane}");
            assert_eq!(batch.tb, scalar.tb, "tb lane {lane}");
            assert_eq!(batch.qb, scalar.qb, "qb lane {lane}");
        }
    }

    #[test]
    fn getScores16_matches_scalar_wrapper_outputs() {
        let k = kswv::ctor(6, 1, 6, 1, 1, -4, 1, None, None);
        let mut pairs = vec![SeqPair {
            idr: 0,
            idq: 0,
            len1: 4,
            len2: 4,
            h0: 0x80000,
            regid: 0,
            ..Default::default()
        }];
        let seq_ref = vec![0_u8, 1, 2, 3];
        let seq_qer = vec![0_u8, 1, 2, 3];
        let mut aln = vec![kswr_t::default(); 1];
        k.getScores16(&mut pairs, &seq_ref, &seq_qer, &mut aln, 1, 1, 0);
        assert!(aln[0].score > 0);
        assert_eq!(aln[0].te, 3);
        assert_eq!(aln[0].qe, 3);
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn getScores16_avx512_matches_scalar_lanes_on_mixed_batch() {
        if !std::is_x86_feature_detected!("avx512bw") {
            return;
        }
        let k = kswv::ctor(6, 1, 6, 1, 1, -4, 1, Some(64), Some(64));
        let mut pairs = Vec::new();
        let mut seq_ref = Vec::new();
        let mut seq_qer = Vec::new();
        for lane in 0..SIMD_WIDTH16 {
            let len1 = 12 + (lane % 7) as i32;
            let len2 = 10 + (lane % 5) as i32;
            let idr = i32::try_from(seq_ref.len()).expect("idr");
            let idq = i32::try_from(seq_qer.len()).expect("idq");
            for i in 0..usize::try_from(len1).expect("len1") {
                seq_ref.push(((i + lane) & 3) as u8);
            }
            for i in 0..usize::try_from(len2).expect("len2") {
                seq_qer.push(((i + lane + usize::from(lane % 3 == 0)) & 3) as u8);
            }
            pairs.push(SeqPair {
                idr,
                idq,
                len1,
                len2,
                h0: if lane % 4 == 0 { 0x40000 | 4 } else { 0x80000 },
                regid: lane as i32,
                ..Default::default()
            });
        }

        let mut mat = [0_i8; 25];
        k.bwa_fill_scmat(&mut mat);
        let q_max = crate::bwa_mem2::ksw::ksw_qmax(k.m, &mat);
        let mut scalar_aln = vec![kswr_t::default(); SIMD_WIDTH16];
        k.kswv_scalar_lanes_i16(&pairs, &seq_ref, &seq_qer, &mut scalar_aln, &mat, q_max, 0);

        let mut avx_aln = vec![kswr_t::default(); SIMD_WIDTH16];
        let mut seq1_soa = Vec::new();
        let mut seq2_soa = Vec::new();
        k.kswv_batch16_chunk_avx512(
            &pairs,
            &seq_ref,
            &seq_qer,
            &mut seq1_soa,
            &mut seq2_soa,
            &mut avx_aln,
            0,
            SIMD_WIDTH16 as i32,
            0,
        );
        for (lane, (avx, scalar)) in avx_aln.iter().zip(&scalar_aln).enumerate() {
            assert_eq!(avx.score, scalar.score, "score lane {lane}");
            assert_eq!(avx.te, scalar.te, "te lane {lane}");
            assert_eq!(avx.qe, scalar.qe, "qe lane {lane}");
            assert_eq!(avx.score2, scalar.score2, "score2 lane {lane}");
            assert_eq!(avx.te2, scalar.te2, "te2 lane {lane}");
        }
    }

    #[test]
    fn decode_soa_lane_restores_lane_sequence() {
        let mut soa = vec![-1_i16; 4 * SIMD_WIDTH16];
        soa[0] = 0;
        soa[SIMD_WIDTH16] = 1;
        soa[2 * SIMD_WIDTH16] = 15;
        soa[3 * SIMD_WIDTH16] = 3;
        assert_eq!(decode_soa_lane(&soa, 0, 4), vec![0, 1, 4, 3]);
    }

    #[test]
    fn kswv512_16_updates_phase_zero_and_phase_one_fields() {
        let k = kswv::ctor(6, 1, 6, 1, 1, -4, 1, None, None);
        let p = vec![
            SeqPair {
                len1: 4,
                len2: 4,
                h0: 0x80000,
                regid: 0,
                ..Default::default()
            };
            SIMD_WIDTH16
        ];
        let mut seq1 = vec![-1_i16; 4 * SIMD_WIDTH16];
        let mut seq2 = vec![-1_i16; 8 * SIMD_WIDTH16];
        for lane in 0..SIMD_WIDTH16 {
            seq1[lane] = 0;
            seq1[SIMD_WIDTH16 + lane] = 1;
            seq1[2 * SIMD_WIDTH16 + lane] = 2;
            seq1[3 * SIMD_WIDTH16 + lane] = 3;
            seq2[lane] = 0;
            seq2[SIMD_WIDTH16 + lane] = 1;
            seq2[2 * SIMD_WIDTH16 + lane] = 2;
            seq2[3 * SIMD_WIDTH16 + lane] = 3;
        }
        let mut aln = vec![kswr_t::default(); SIMD_WIDTH16];
        k.kswv512_16(&seq1, &seq2, 4, 8, &p, &mut aln, 0, 0, 1, 0);
        assert!(aln[0].score > 0);
        assert_eq!(aln[0].te, 3);
        assert_eq!(aln[0].qe, 3);

        let mut rev1 = seq1.clone();
        let mut rev2 = seq2.clone();
        rev1[0] = 3;
        rev1[SIMD_WIDTH16] = 2;
        rev1[2 * SIMD_WIDTH16] = 1;
        rev1[3 * SIMD_WIDTH16] = 0;
        rev2[0] = 3;
        rev2[SIMD_WIDTH16] = 2;
        rev2[2 * SIMD_WIDTH16] = 1;
        rev2[3 * SIMD_WIDTH16] = 0;
        k.kswv512_16(&rev1, &rev2, 4, 8, &p, &mut aln, 0, 0, 1, 1);
        assert_eq!(aln[0].tb, 0);
        assert_eq!(aln[0].qb, 0);
    }

    #[test]
    fn decode_soa_lane_u8_restores_lane_sequence() {
        let mut soa = vec![0xff_u8; 4 * SIMD_WIDTH8];
        soa[0] = 0;
        soa[SIMD_WIDTH8] = 1;
        soa[2 * SIMD_WIDTH8] = 4;
        soa[3 * SIMD_WIDTH8] = 3;
        assert_eq!(decode_soa_lane_u8(&soa, 0, 4), vec![0, 1, 4, 3]);
        assert_eq!(decode_soa_lane_u8_query(&soa, 0, 4), vec![0, 1, 4, 3]);
    }

    #[test]
    fn kswv512_u8_updates_phase_zero_and_phase_one_fields() {
        let k = kswv::ctor(6, 1, 6, 1, 1, -4, 1, None, None);
        let p = vec![
            SeqPair {
                len1: 4,
                len2: 4,
                h0: 0x10000 | 0x80000,
                regid: 0,
                ..Default::default()
            };
            SIMD_WIDTH8
        ];
        let mut seq1 = vec![0xff_u8; 4 * SIMD_WIDTH8];
        let mut seq2 = vec![0xff_u8; 16 * SIMD_WIDTH8];
        for lane in 0..SIMD_WIDTH8 {
            seq1[lane] = 0;
            seq1[SIMD_WIDTH8 + lane] = 1;
            seq1[2 * SIMD_WIDTH8 + lane] = 2;
            seq1[3 * SIMD_WIDTH8 + lane] = 3;
            seq2[lane] = 0;
            seq2[SIMD_WIDTH8 + lane] = 1;
            seq2[2 * SIMD_WIDTH8 + lane] = 2;
            seq2[3 * SIMD_WIDTH8 + lane] = 3;
        }
        let mut aln = vec![kswr_t::default(); SIMD_WIDTH8];
        k.kswv512_u8(&seq1, &seq2, 4, 16, &p, &mut aln, 0, 0, 1, 0);
        assert!(aln[0].score > 0);
        assert_eq!(aln[0].te, 3);
        assert_eq!(aln[0].qe, 3);

        let mut rev1 = seq1.clone();
        let mut rev2 = seq2.clone();
        rev1[0] = 3;
        rev1[SIMD_WIDTH8] = 2;
        rev1[2 * SIMD_WIDTH8] = 1;
        rev1[3 * SIMD_WIDTH8] = 0;
        rev2[0] = 3;
        rev2[SIMD_WIDTH8] = 2;
        rev2[2 * SIMD_WIDTH8] = 1;
        rev2[3 * SIMD_WIDTH8] = 0;
        k.kswv512_u8(&rev1, &rev2, 4, 16, &p, &mut aln, 0, 0, 1, 1);
        assert_eq!(aln[0].tb, 0);
        assert_eq!(aln[0].qb, 0);
    }
}
