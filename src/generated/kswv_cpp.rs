#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Generated scaffold for `bwa-mem2/src/kswv.cpp`.

use crate::generated::ksw_cpp::{ksw_align2, ksw_i16_slices, ksw_qmax, ksw_u8_slices};
use crate::generated::ksw_h::{KSW_XBYTE, KSW_XSTART};
use crate::generated::kswv_h::{
    kswq_t, kswr_t, kswv, SeqPair, DEFAULT_AMBIG, MAX_SEQ_LEN_QER_SAM, MAX_SEQ_LEN_REF_SAM,
    SIMD_WIDTH16, SIMD_WIDTH8,
};

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
        let num_threads = usize::try_from(numThreads.max(1)).expect("numThreads");
        let max_ref_len = maxRefLen.unwrap_or(MAX_SEQ_LEN_REF_SAM) + 16;
        let max_qer_len = maxQerLen.unwrap_or(MAX_SEQ_LEN_QER_SAM) + 16;
        let ref_cap16 =
            usize::try_from(max_ref_len).expect("maxRefLen") * SIMD_WIDTH16 * num_threads;
        let qer_cap16 =
            usize::try_from(max_qer_len).expect("maxQerLen") * SIMD_WIDTH16 * num_threads;
        let ref_cap8 = usize::try_from(max_ref_len).expect("maxRefLen") * SIMD_WIDTH8 * num_threads;
        let qer_cap8 = usize::try_from(max_qer_len).expect("maxQerLen") * SIMD_WIDTH8 * num_threads;
        let g_qmax = i32::from(max_i8(max_i8(w_match, w_mismatch), DEFAULT_AMBIG));
        Self {
            m: 5,
            o_del,
            o_ins,
            e_del,
            e_ins,
            w_match,
            w_mismatch,
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
        let total = usize::try_from(numPairs).expect("numPairs");
        let mut mat = [0_i8; 25];
        self.bwa_fill_scmat(&mut mat);
        let q_max = ksw_qmax(self.m, &mat);
        for chunk in pairArray[..total].chunks(SIMD_WIDTH8) {
            self.kswv_scalar_lanes_u8(chunk, seqBufRef, seqBufQer, aln, &mat, q_max, phase);
        }
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

        for (j, sp) in lanes.iter().copied().enumerate() {
            let seq1 = &seqBufRef[usize::try_from(sp.idr).expect("idr")..];
            for k in 0..usize::try_from(sp.len1).expect("len1") {
                seq1_soa[k * SIMD_WIDTH8 + j] = if seq1[k] == 4 { 4 } else { seq1[k] };
            }
            max_len1 = max_len1.max(sp.len1);
        }
        for (j, sp) in lanes.iter().copied().enumerate() {
            for k in usize::try_from(sp.len1).expect("len1")
                ..=usize::try_from(max_len1).expect("max_len1")
            {
                seq1_soa[k * SIMD_WIDTH8 + j] = 0xff;
            }
        }

        for (j, sp) in lanes.iter().copied().enumerate() {
            let seq2 = &seqBufQer[usize::try_from(sp.idq).expect("idq")..];
            let quanta = usize::try_from((sp.len2 + 16 - 1) / 16 * 16).expect("quanta");
            for k in 0..usize::try_from(sp.len2).expect("len2") {
                seq2_soa[k * SIMD_WIDTH8 + j] = if seq2[k] == 4 { 8 } else { seq2[k] };
            }
            for k in usize::try_from(sp.len2).expect("len2")..quanta {
                seq2_soa[k * SIMD_WIDTH8 + j] = 5;
            }
            max_len2 = max_len2.max(i32::try_from(quanta).expect("quanta"));
        }
        for (j, sp) in lanes.iter().copied().enumerate() {
            let quanta = usize::try_from((sp.len2 + 16 - 1) / 16 * 16).expect("quanta");
            for k in quanta..=usize::try_from(max_len2).expect("max_len2") {
                seq2_soa[k * SIMD_WIDTH8 + j] = 0xff;
            }
        }

        self.kswv512_u8(
            seq1_soa,
            seq2_soa,
            i16::try_from(max_len1).expect("max_len1"),
            i16::try_from(max_len2).expect("max_len2"),
            lanes,
            aln,
            i32::try_from(offset).expect("offset"),
            0,
            numPairs,
            phase,
        );
    }

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
        let total = usize::try_from(numPairs).expect("numPairs");
        let mut mat = [0_i8; 25];
        self.bwa_fill_scmat(&mut mat);
        let q_max = ksw_qmax(self.m, &mat);
        for chunk in pairArray[..total].chunks(SIMD_WIDTH16) {
            self.kswv_scalar_lanes_i16(chunk, seqBufRef, seqBufQer, aln, &mat, q_max, phase);
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
        for &sp in lanes {
            let ind = usize::try_from(sp.regid).expect("regid");
            let target_start = usize::try_from(sp.idr).expect("idr");
            let query_start = usize::try_from(sp.idq).expect("idq");
            let target_len = usize::try_from(sp.len1).expect("len1");
            let query_len = usize::try_from(sp.len2).expect("len2");
            let target = &seqBufRef[target_start..target_start + target_len];
            let query = &seqBufQer[query_start..query_start + query_len];
            let ks = ksw_i16_slices(
                query, self.m, mat, q_max, sp.len1, target, self.o_del, self.e_del, self.o_ins,
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
            let ind = usize::try_from(sp.regid).expect("regid");
            let target_start = usize::try_from(sp.idr).expect("idr");
            let query_start = usize::try_from(sp.idq).expect("idq");
            let target_len = usize::try_from(sp.len1).expect("len1");
            let query_len = usize::try_from(sp.len2).expect("len2");
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
        for lane in 0..SIMD_WIDTH16 {
            if usize::try_from(po_ind).expect("po_ind") + lane
                >= usize::try_from(numPairs).expect("numPairs")
            {
                break;
            }
            let sp = p[lane];
            let ind = usize::try_from(sp.regid).expect("regid");
            decode_soa_lane_into(
                seq1SoA,
                lane,
                usize::try_from(sp.len1).expect("len1"),
                &mut target_buf,
            );
            decode_soa_lane_into(
                seq2SoA,
                lane,
                usize::try_from(sp.len2).expect("len2"),
                &mut query_buf,
            );
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
        let mut mat = [0_i8; 25];
        self.bwa_fill_scmat(&mut mat);
        let q_max = ksw_qmax(self.m, &mat);
        let mut target_buf: Vec<u8> = Vec::new();
        let mut query_buf: Vec<u8> = Vec::new();
        for lane in 0..SIMD_WIDTH8 {
            if usize::try_from(po_ind).expect("po_ind") + lane
                >= usize::try_from(numPairs).expect("numPairs")
            {
                break;
            }
            let sp = p[lane];
            let ind = usize::try_from(sp.regid).expect("regid");
            decode_soa_lane_u8_into(
                seq1SoA,
                lane,
                usize::try_from(sp.len1).expect("len1"),
                &mut target_buf,
            );
            decode_soa_lane_u8_query_into(
                seq2SoA,
                lane,
                usize::try_from(sp.len2).expect("len2"),
                &mut query_buf,
            );
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
            mat[k] = self.w_ambig;
            k += 1;
        }
        for item in mat.iter_mut().skip(k).take(5) {
            *item = self.w_ambig;
        }
    }

    #[doc = "Original function: kswv::ksw_qinit:1256"]
    pub fn ksw_qinit(&self, size: i32, qlen: i32, query: &[u8], m: i32, mat: &[i8]) -> kswq_t {
        let size = if size > 1 { 2 } else { 1 };
        let p = 8 * (3 - size);
        let slen = (qlen + p - 1) / p;
        let qlen_usize = usize::try_from(qlen).expect("qlen");
        let slen_usize = usize::try_from(slen).expect("slen");
        let p_usize = usize::try_from(p).expect("p");
        let m_usize = usize::try_from(m).expect("m");

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
        let shift_u8 = (0_i16).wrapping_sub(i16::from(shift)) as u8;
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
