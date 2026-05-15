#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Port of `bwa-mem2/src/ksw.h` + `bwa-mem2/src/ksw.cpp`.

// --- ksw.h ---

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

// --- ksw.cpp ---

#[doc = "Original struct: eh_t (bwa-mem2/src/ksw.cpp)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct eh_t {
    pub h: i32,
    pub e: i32,
}

// Thread-local scratch for ksw_global2 (z, qp, eh). The z traceback matrix is the
// largest reuser — typically n_col * tlen ≈ 22KB per call, allocated 100K-150K times
// over a 50K-read run (per-alignment cigar gen via bwa_gen_cigar2).
thread_local! {
    static KSW_GLOBAL_SCRATCH: std::cell::RefCell<(Vec<u8>, Vec<i8>, Vec<eh_t>)> =
        const { std::cell::RefCell::new((Vec::new(), Vec::new(), Vec::new())) };
    // Reused across ksw_extend2 calls. eh/qp grow to ~qlen*m bytes per call. With chain
    // extension scalar fallback hitting this path, pooling avoids per-pair Vec allocations.
    static KSW_EXTEND_SCRATCH: std::cell::RefCell<(Vec<eh_t>, Vec<i8>)> =
        const { std::cell::RefCell::new((Vec::new(), Vec::new())) };
}

const MINUS_INF: i32 = -0x4000_0000;
const G_DEFR: kswr_t = kswr_t {
    score: 0,
    te: -1,
    qe: -1,
    score2: -1,
    te2: -1,
    tb: -1,
    qb: -1,
};

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
#[doc = "Original function: ksw_qinit:62"]
#[inline]
pub fn ksw_qinit(size: i32, qlen: i32, query: &[u8], m: i32, mat: &[i8]) -> _kswq_t {
    let size = if size > 1 { 2 } else { 1 };
    let p = 8 * (3 - size); // # values per __m128i
    let slen = (qlen + p - 1) / p; // segmented length
    // find the minimum and maximum score (mdiff = max-min is the dynamic range)
    let mut shift = i8::MAX;
    let mut mdiff = i8::MIN;
    for &score in mat.iter().take((m * m) as usize) {
        if score < shift {
            shift = score;
        }
        if score > mdiff {
            mdiff = score;
        }
    }
    let max = mdiff as u8;
    _kswq_t {
        qlen,
        slen,
        shift: (0_i16).wrapping_sub(i16::from(shift)) as u8,
        mdiff: mdiff.wrapping_add((0_i16).wrapping_sub(i16::from(shift)) as i8) as u8,
        max,
        size: size as u8,
        query: query[..qlen as usize].to_vec(),
        m,
        mat: mat.to_vec(),
    }
}

#[doc = "Original function: ksw_u8:111"]
pub fn ksw_u8(
    q: &_kswq_t,
    tlen: i32,
    target: &[u8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    xtra: i32,
) -> kswr_t {
    let tlen_usize = tlen as usize;
    let endsc_in = if (xtra & crate::bwa_mem2::ksw::KSW_XSTOP) != 0 {
        xtra & 0xffff
    } else {
        i32::MAX
    };
    let (mut r, row_maxes) = local_align_affine_with_rows_endsc(
        &q.query,
        &target[..tlen_usize],
        q.m,
        &q.mat,
        o_del,
        e_del,
        o_ins,
        e_ins,
        endsc_in,
    );
    if r.score > 255 {
        r.score = 255;
    }
    let minsc = if (xtra & crate::bwa_mem2::ksw::KSW_XSUBO) != 0 {
        xtra & 0xffff
    } else {
        0x10000
    };
    let endsc = if (xtra & crate::bwa_mem2::ksw::KSW_XSTOP) != 0 {
        xtra & 0xffff
    } else {
        0x10000
    };
    // C++ ksw.cpp:213 wraps the score2/te2 search in `if (r.score != 255)` — when SW saturates
    // at 255, both stay at -1 (G_DEFR defaults).
    if r.score != 255 && r.score >= minsc {
        let radius = (r.score + i32::from(q.max).saturating_sub(1)) / i32::from(q.max.max(1));
        let low = r.te - radius;
        let high = r.te + radius;
        let b = build_score_peaks(&row_maxes, minsc);
        let mut alt = G_DEFR;
        for &(imax, te) in &b {
            if (te < low || te > high) && imax > alt.score2 {
                alt.score2 = imax.min(255);
                alt.te2 = te;
            }
        }
        r.score2 = alt.score2;
        r.te2 = alt.te2;
    }
    if r.score >= endsc && endsc != 0x10000 {
        return r;
    }
    r
}

#[doc = "Original function: ksw_i16:234"]
pub fn ksw_i16(
    q: &_kswq_t,
    tlen: i32,
    target: &[u8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    xtra: i32,
) -> kswr_t {
    let tlen_usize = tlen as usize;
    let endsc_in = if (xtra & crate::bwa_mem2::ksw::KSW_XSTOP) != 0 {
        xtra & 0xffff
    } else {
        i32::MAX
    };
    let (mut r, row_maxes) = local_align_affine_with_rows_endsc(
        &q.query,
        &target[..tlen_usize],
        q.m,
        &q.mat,
        o_del,
        e_del,
        o_ins,
        e_ins,
        endsc_in,
    );
    let minsc = if (xtra & crate::bwa_mem2::ksw::KSW_XSUBO) != 0 {
        xtra & 0xffff
    } else {
        0x10000
    };
    let endsc = if (xtra & crate::bwa_mem2::ksw::KSW_XSTOP) != 0 {
        xtra & 0xffff
    } else {
        0x10000
    };
    if r.score >= minsc {
        let radius = (r.score + i32::from(q.max).saturating_sub(1)) / i32::from(q.max.max(1));
        let low = r.te - radius;
        let high = r.te + radius;
        let b = build_score_peaks(&row_maxes, minsc);
        let mut best2 = -1;
        let mut te2 = -1;
        for &(imax, te) in &b {
            if (te < low || te > high) && imax > best2 {
                best2 = imax;
                te2 = te;
            }
        }
        r.score2 = best2;
        r.te2 = te2;
    }
    if r.score >= endsc && endsc != 0x10000 {
        return r;
    }
    r
}

// ============================================================================================
// SSE2 striped Farrar SW kernel — Rust port of bwa-mem2/src/ksw.cpp:111 (`ksw_u8`).
//
// The C++ kernel processes 16 query positions per __m128i in a striped layout:
//   stripe i, lane k -> query[i + k*slen]   (i ∈ 0..slen, k ∈ 0..16)
// where slen = ceil(qlen/16). Per outer iteration over the target, the inner loop walks
// the slen stripes performing the SW recurrence in SIMD, then a "lazy F" propagation
// loop fixes up any cells where F dominates.
//
// We mirror the C++ exactly to preserve byte-identical output (max position tie-breaking
// included). Produces a `kswr_t` matching the scalar port; falls back to scalar
// `local_align_affine_with_rows_endsc` if SSE2 is unavailable (always available on x86_64).
// ============================================================================================
// SSE2 striped Farrar SW kernel for i16 — Rust port of bwa-mem2/src/ksw.cpp:234 (`ksw_i16`).
// Same striped layout as ksw_u8 but 8 i16 lanes per __m128i (instead of 16 u8). No `shift`
// since i16 has plenty of headroom.
#[cfg(target_arch = "x86_64")]
mod ksw_i16_sse2 {
    use super::{kswr_t, G_DEFR};
    use crate::bwa_mem2::ksw::{KSW_XSTOP, KSW_XSUBO};
    use core::arch::x86_64::*;
    use std::cell::RefCell;

    #[derive(Default)]
    struct Scratch {
        qp: Vec<i16>,
        h0: Vec<__m128i>,
        h1: Vec<__m128i>,
        e_buf: Vec<__m128i>,
        hmax_buf: Vec<__m128i>,
        peaks: Vec<(i32, i32)>,
    }

    thread_local! {
        static SCRATCH: RefCell<Scratch> = RefCell::new(Scratch::default());
    }

    fn build_qp_into(query: &[u8], m: i32, mat: &[i8], qp: &mut Vec<i16>) -> i32 {
        let qlen = query.len() as i32;
        let slen = (qlen + 7) / 8;
        let m_usize = m as usize;
        let nlen = slen * 8;
        let total = m_usize * (slen as usize) * 8;
        qp.clear();
        if qp.capacity() < total {
            qp.reserve(total - qp.capacity());
        }
        // Every profile slot is assigned below, including padding cells, so avoid
        // zero-initializing the buffer before immediately overwriting it.
        unsafe {
            qp.set_len(total);
        }
        let mut t = 0_usize;
        for a in 0..m_usize {
            let ma = &mat[a * m_usize..a * m_usize + m_usize];
            for i in 0..slen {
                let mut k = i;
                while k < nlen {
                    qp[t] = if k >= qlen {
                        0
                    } else {
                        i16::from(ma[usize::from(query[k as usize])])
                    };
                    t += 1;
                    k += slen;
                }
            }
        }
        slen
    }

    #[target_feature(enable = "sse2")]
    unsafe fn hmax_i16(mut x: __m128i) -> i32 {
        x = _mm_max_epi16(x, _mm_srli_si128::<8>(x));
        x = _mm_max_epi16(x, _mm_srli_si128::<4>(x));
        x = _mm_max_epi16(x, _mm_srli_si128::<2>(x));
        _mm_extract_epi16::<0>(x)
    }

    pub(super) fn run(
        query: &[u8],
        m: i32,
        mat: &[i8],
        q_max: u8,
        tlen: i32,
        target: &[u8],
        o_del: i32,
        e_del: i32,
        o_ins: i32,
        e_ins: i32,
        xtra: i32,
    ) -> kswr_t {
        SCRATCH.with(|cell| {
            let mut scratch = cell.borrow_mut();
            run_with_scratch(
                query,
                m,
                mat,
                q_max,
                tlen,
                target,
                o_del,
                e_del,
                o_ins,
                e_ins,
                xtra,
                &mut scratch,
            )
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn run_with_scratch(
        query: &[u8],
        m: i32,
        mat: &[i8],
        q_max: u8,
        tlen: i32,
        target: &[u8],
        o_del: i32,
        e_del: i32,
        o_ins: i32,
        e_ins: i32,
        xtra: i32,
        scratch: &mut Scratch,
    ) -> kswr_t {
        let slen = build_qp_into(query, m, mat, &mut scratch.qp);
        let slen_us = slen as usize;
        let tlen_us = tlen as usize;
        let minsc = if xtra & KSW_XSUBO != 0 {
            xtra & 0xffff
        } else {
            0x10000
        };
        let endsc = if xtra & KSW_XSTOP != 0 {
            xtra & 0xffff
        } else {
            0x10000
        };

        let zero = unsafe { _mm_setzero_si128() };
        scratch.h0.clear();
        scratch.h0.resize(slen_us, zero);
        scratch.h1.clear();
        if scratch.h1.capacity() < slen_us {
            scratch.h1.reserve(slen_us - scratch.h1.capacity());
        }
        // h1 is assigned for every stripe before it is read by the lazy-F pass.
        unsafe {
            scratch.h1.set_len(slen_us);
        }
        scratch.e_buf.clear();
        scratch.e_buf.resize(slen_us, zero);
        scratch.hmax_buf.clear();
        scratch.hmax_buf.resize(slen_us, zero);
        scratch.peaks.clear();

        let mut gmax: i32 = 0;
        let mut te: i32 = -1;

        unsafe {
            let oe_del_v = _mm_set1_epi16((o_del + e_del) as i16);
            let e_del_v = _mm_set1_epi16(e_del as i16);
            let oe_ins_v = _mm_set1_epi16((o_ins + e_ins) as i16);
            let e_ins_v = _mm_set1_epi16(e_ins as i16);

            'outer: for i in 0..tlen_us {
                let mut f = zero;
                let mut max_v = zero;
                let s_base = (target[i] as usize) * slen_us;
                let mut h = _mm_slli_si128::<2>(*scratch.h0.get_unchecked(slen_us - 1));
                for j in 0..slen_us {
                    let s_v = _mm_loadu_si128(
                        scratch.qp.as_ptr().add(s_base * 8 + j * 8) as *const __m128i
                    );
                    h = _mm_adds_epi16(h, s_v);
                    let mut e = *scratch.e_buf.get_unchecked(j);
                    h = _mm_max_epi16(h, e);
                    h = _mm_max_epi16(h, f);
                    max_v = _mm_max_epi16(max_v, h);
                    *scratch.h1.get_unchecked_mut(j) = h;
                    e = _mm_subs_epu16(e, e_del_v);
                    let t = _mm_subs_epu16(h, oe_del_v);
                    e = _mm_max_epi16(e, t);
                    *scratch.e_buf.get_unchecked_mut(j) = e;
                    f = _mm_subs_epu16(f, e_ins_v);
                    let t = _mm_subs_epu16(h, oe_ins_v);
                    f = _mm_max_epi16(f, t);
                    h = *scratch.h0.get_unchecked(j);
                }
                'lazyf: for _k in 0..16 {
                    f = _mm_slli_si128::<2>(f);
                    for j in 0..slen_us {
                        let mut h = *scratch.h1.get_unchecked(j);
                        h = _mm_max_epi16(h, f);
                        *scratch.h1.get_unchecked_mut(j) = h;
                        let h2 = _mm_subs_epu16(h, oe_ins_v);
                        f = _mm_subs_epu16(f, e_ins_v);
                        // C++: !movemask(cmpgt(f, h)) → exit when no F lane exceeds H.
                        if _mm_movemask_epi8(_mm_cmpgt_epi16(f, h2)) == 0 {
                            break 'lazyf;
                        }
                    }
                }
                let imax = hmax_i16(max_v);
                if imax >= minsc {
                    let cur_i = i as i32;
                    if let Some(last) = scratch.peaks.last_mut() {
                        if last.1 + 1 == cur_i {
                            if last.0 < imax {
                                *last = (imax, cur_i);
                            }
                        } else {
                            scratch.peaks.push((imax, cur_i));
                        }
                    } else {
                        scratch.peaks.push((imax, cur_i));
                    }
                }
                if imax > gmax {
                    gmax = imax;
                    te = i as i32;
                    for j in 0..slen_us {
                        *scratch.hmax_buf.get_unchecked_mut(j) = *scratch.h1.get_unchecked(j);
                    }
                    if gmax >= endsc {
                        break 'outer;
                    }
                }
                std::mem::swap(&mut scratch.h0, &mut scratch.h1);
            }
        }

        let mut r = G_DEFR;
        r.score = gmax;
        r.te = te;
        // qe via Hmax scan (ksw.cpp:320-324).
        let total = slen_us * 8;
        let mut max_v = -1_i32;
        let mut qe = -1_i32;
        unsafe {
            let hmax_words = scratch.hmax_buf.as_ptr() as *const u16;
            for i in 0..total {
                let v = *hmax_words.add(i) as i32;
                let pos = (i / 8) + (i % 8) * slen_us;
                if v > max_v {
                    max_v = v;
                    qe = pos as i32;
                } else if v == max_v && (pos as i32) < qe {
                    qe = pos as i32;
                }
            }
        }
        r.qe = qe;

        if r.score >= minsc {
            let radius = (r.score + i32::from(q_max).saturating_sub(1)) / i32::from(q_max.max(1));
            let low = r.te - radius;
            let high = r.te + radius;
            let mut best2 = -1;
            let mut te2 = -1;
            for &(imax, te2_cand) in &scratch.peaks {
                if (te2_cand < low || te2_cand > high) && imax > best2 {
                    best2 = imax;
                    te2 = te2_cand;
                }
            }
            r.score2 = best2;
            r.te2 = te2;
        }
        r
    }
}

#[cfg(target_arch = "x86_64")]
mod ksw_u8_sse2 {
    use super::{kswr_t, G_DEFR};
    use crate::bwa_mem2::ksw::{KSW_XSTOP, KSW_XSUBO};
    use core::arch::x86_64::*;
    use std::cell::RefCell;

    #[derive(Default)]
    struct Scratch {
        qp: Vec<u8>,
        h0: Vec<__m128i>,
        h1: Vec<__m128i>,
        e_buf: Vec<__m128i>,
        hmax_buf: Vec<__m128i>,
        peaks: Vec<(i32, i32)>,
    }

    thread_local! {
        static SCRATCH: RefCell<Scratch> = RefCell::new(Scratch::default());
    }

    /// Build the C++-style SoA query profile: per-base block of slen __m128i, each holding
    /// `query[stripe + lane*slen]` (or 0 padding) plus a `q_shift` bias so that subtraction
    /// in the SW recurrence works in u8 arithmetic.
    ///
    /// Returns `(slen, q_shift)`. `qp` length is `m * slen * 16` u8.
    fn build_qp_into(query: &[u8], m: i32, mat: &[i8], qp: &mut Vec<u8>) -> (i32, u8) {
        let qlen = query.len() as i32;
        let slen = (qlen + 15) / 16; // ceil
        let m_usize = m as usize;
        // Find min mat value (becomes q_shift = 256 - min, in u8 wraparound).
        let mut shift_min: i8 = i8::MAX;
        for &s in mat.iter().take(m_usize * m_usize) {
            if s < shift_min {
                shift_min = s;
            }
        }
        let q_shift = (0_i16).wrapping_sub(i16::from(shift_min)) as u8;
        let nlen = slen * 16;
        qp.clear();
        qp.resize(m_usize * (slen as usize) * 16, 0);
        let mut t = 0_usize;
        for a in 0..m_usize {
            let ma = &mat[a * m_usize..a * m_usize + m_usize];
            for i in 0..slen {
                let mut k = i;
                while k < nlen {
                    let v: u8 = if k >= qlen {
                        0
                    } else {
                        ma[usize::from(query[k as usize])] as u8
                    };
                    qp[t] = v.wrapping_add(q_shift);
                    t += 1;
                    k += slen;
                }
            }
        }
        (slen, q_shift)
    }

    /// Horizontal max of 16 u8 lanes — returns the max as u32.
    #[target_feature(enable = "sse2")]
    unsafe fn hmax_u8(mut x: __m128i) -> u32 {
        x = _mm_max_epu8(x, _mm_srli_si128::<8>(x));
        x = _mm_max_epu8(x, _mm_srli_si128::<4>(x));
        x = _mm_max_epu8(x, _mm_srli_si128::<2>(x));
        x = _mm_max_epu8(x, _mm_srli_si128::<1>(x));
        (_mm_extract_epi16::<0>(x) & 0x00ff) as u32
    }

    /// SSE2 striped Farrar SW. Direct port of ksw.cpp:111 ksw_u8.
    pub(super) fn run(
        query: &[u8],
        m: i32,
        mat: &[i8],
        q_max: u8,
        tlen: i32,
        target: &[u8],
        o_del: i32,
        e_del: i32,
        o_ins: i32,
        e_ins: i32,
        xtra: i32,
    ) -> kswr_t {
        SCRATCH.with(|cell| {
            let mut scratch = cell.borrow_mut();
            run_with_scratch(
                query,
                m,
                mat,
                q_max,
                tlen,
                target,
                o_del,
                e_del,
                o_ins,
                e_ins,
                xtra,
                &mut scratch,
            )
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn run_with_scratch(
        query: &[u8],
        m: i32,
        mat: &[i8],
        q_max: u8,
        tlen: i32,
        target: &[u8],
        o_del: i32,
        e_del: i32,
        o_ins: i32,
        e_ins: i32,
        xtra: i32,
        scratch: &mut Scratch,
    ) -> kswr_t {
        let (slen, q_shift) = build_qp_into(query, m, mat, &mut scratch.qp);
        let slen_us = slen as usize;
        let tlen_us = tlen as usize;

        let minsc = if xtra & KSW_XSUBO != 0 {
            xtra & 0xffff
        } else {
            0x10000
        };
        let endsc = if xtra & KSW_XSTOP != 0 {
            xtra & 0xffff
        } else {
            0x10000
        };

        // Aligned scratch — Vec<__m128i> guarantees 16-byte alignment.
        let zero_m = unsafe { _mm_setzero_si128() };
        scratch.h0.clear();
        scratch.h0.resize(slen_us, zero_m);
        scratch.h1.clear();
        scratch.h1.resize(slen_us, zero_m);
        scratch.e_buf.clear();
        scratch.e_buf.resize(slen_us, zero_m);
        scratch.hmax_buf.clear();
        scratch.hmax_buf.resize(slen_us, zero_m);

        // Track score peaks (b[] in C++): list of (imax, i) where each entry holds the
        // running max within a contiguous run of i values that have imax >= minsc.
        scratch.peaks.clear();
        let mut gmax: u32 = 0;
        let mut te: i32 = -1;

        let qp = scratch.qp.as_slice();
        let mut h0 = &mut scratch.h0;
        let mut h1 = &mut scratch.h1;
        let e_buf = &mut scratch.e_buf;
        let hmax_buf = &mut scratch.hmax_buf;
        let peaks = &mut scratch.peaks;

        unsafe {
            let oe_del_v = _mm_set1_epi8((o_del + e_del) as i8);
            let e_del_v = _mm_set1_epi8(e_del as i8);
            let oe_ins_v = _mm_set1_epi8((o_ins + e_ins) as i8);
            let e_ins_v = _mm_set1_epi8(e_ins as i8);
            let shift_v = _mm_set1_epi8(q_shift as i8);
            let zero = _mm_setzero_si128();

            // the core loop
            'outer: for i in 0..tlen_us {
                let mut f = zero;
                let mut max_v = zero;
                // s is the 1st score vector for target[i]
                let s_base = (target[i] as usize) * slen_us;
                // h = H0[slen-1] shifted left by 1 byte (logical "<<" = SIMD right-shift in
                // little-endian byte order; matches `_mm_slli_si128(h, 1)` in C++).
                // h={2,5,8,11,14,17,-1,-1} in the segmentation example; this is H(i-1,-1).
                let mut h = _mm_slli_si128::<1>(*h0.get_unchecked(slen_us - 1));
                let mut j = 0;
                while j < slen_us {
                    // SW cells are computed in the following order:
                    //   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
                    //   E(i+1,j) = max{H(i,j)-q, E(i,j)-r}
                    //   F(i,j+1) = max{H(i,j)-q, F(i,j)-r}
                    let s_v_bytes = qp.as_ptr().add(s_base * 16 + j * 16) as *const __m128i;
                    let s_v = _mm_loadu_si128(s_v_bytes);
                    // compute H'(i,j); note that at the beginning, h=H'(i-1,j-1)
                    h = _mm_adds_epu8(h, s_v);
                    h = _mm_subs_epu8(h, shift_v); // h = H'(i-1,j-1) + S(i,j)
                    let mut e = *e_buf.get_unchecked(j); // e = E'(i,j)
                    h = _mm_max_epu8(h, e);
                    h = _mm_max_epu8(h, f); // h = H'(i,j)
                    max_v = _mm_max_epu8(max_v, h); // set max
                    *h1.get_unchecked_mut(j) = h; // save to H'(i,j)
                    // now compute E'(i+1,j) = max(E - e_del, H - oe_del)
                    e = _mm_subs_epu8(e, e_del_v); // e = E'(i,j) - e_del
                    let t = _mm_subs_epu8(h, oe_del_v); // h = H'(i,j) - o_del - e_del
                    e = _mm_max_epu8(e, t); // e = E'(i+1,j)
                    *e_buf.get_unchecked_mut(j) = e; // save to E'(i+1,j)
                    // now compute F'(i, j+1) = max(F - e_ins, H - oe_ins)
                    f = _mm_subs_epu8(f, e_ins_v);
                    let t = _mm_subs_epu8(h, oe_ins_v); // h = H'(i,j) - o_ins - e_ins
                    f = _mm_max_epu8(f, t);
                    // get H'(i-1,j) and prepare for the next stripe iteration.
                    h = *h0.get_unchecked(j); // h = H'(i-1,j)
                    j += 1;
                }
                // NB: we do not need to set E(i,j) as we disallow adjacent insertion and then deletion.
                // Lazy F propagation (SWPS3-style); H(i,j) updated in the lazy-F loop cannot exceed max.
                // 16 max passes are sufficient for u8.
                'lazyf: for _k in 0..16 {
                    f = _mm_slli_si128::<1>(f);
                    for j in 0..slen_us {
                        let mut h = *h1.get_unchecked(j);
                        h = _mm_max_epu8(h, f);
                        *h1.get_unchecked_mut(j) = h;
                        let h2 = _mm_subs_epu8(h, oe_ins_v);
                        f = _mm_subs_epu8(f, e_ins_v);
                        let cmp = _mm_movemask_epi8(_mm_cmpeq_epi8(_mm_subs_epu8(f, h2), zero));
                        if cmp == 0xffff {
                            break 'lazyf;
                        }
                    }
                }
                let imax = hmax_u8(max_v) as i32;

                // b[] peak update (matches ksw.cpp:194-202).
                if imax >= minsc {
                    let cur_i = i as i32;
                    if let Some(last) = peaks.last_mut() {
                        if last.1 + 1 == cur_i {
                            if last.0 < imax {
                                *last = (imax, cur_i);
                            }
                        } else {
                            peaks.push((imax, cur_i));
                        }
                    } else {
                        peaks.push((imax, cur_i));
                    }
                }
                if imax as u32 > gmax {
                    gmax = imax as u32;
                    te = i as i32;
                    // Snapshot Hmax = H1.
                    for j in 0..slen_us {
                        *hmax_buf.get_unchecked_mut(j) = *h1.get_unchecked(j);
                    }
                    // Saturation / endsc early exit.
                    if gmax + u32::from(q_shift) >= 255 || (gmax as i32) >= endsc {
                        break 'outer;
                    }
                }
                // Swap H0 ↔ H1.
                std::mem::swap(&mut h0, &mut h1);
            }
        }

        let mut r = G_DEFR;
        // r.score saturates at 255 if the SW result is >= 255 - shift.
        r.score = if (gmax + u32::from(q_shift)) < 255 {
            gmax as i32
        } else {
            255
        };
        r.te = te;

        if r.score != 255 {
            // Find r.qe by scanning Hmax (un-stripe to query coordinates and find the rightmost
            // max position). C++ ksw.cpp:217-218: ties prefer the smaller i index.
            let total = slen_us * 16;
            let hmax_bytes = hmax_buf.as_ptr() as *const u8;
            let mut max_v = -1_i32;
            let mut qe = -1_i32;
            unsafe {
                for i in 0..total {
                    let v = *hmax_bytes.add(i) as i32;
                    let pos = (i / 16) + (i % 16) * slen_us;
                    if v > max_v {
                        max_v = v;
                        qe = pos as i32;
                    } else if v == max_v && (pos as i32) < qe {
                        qe = pos as i32;
                    }
                }
            }
            r.qe = qe;

            // score2/te2 search using the b[] peaks (matches C++ ksw.cpp:222-228).
            if r.score >= minsc {
                let radius =
                    (r.score + i32::from(q_max).saturating_sub(1)) / i32::from(q_max.max(1));
                let low = r.te - radius;
                let high = r.te + radius;
                let mut alt = G_DEFR;
                for &(imax, te2) in peaks.iter() {
                    if (te2 < low || te2 > high) && imax > alt.score2 {
                        alt.score2 = imax.min(255);
                        alt.te2 = te2;
                    }
                }
                r.score2 = alt.score2;
                r.te2 = alt.te2;
            }
        }
        r
    }
}

// Borrowing fast-paths for ksw_u8/ksw_i16 used by kswv hot path. The default ksw_u8 builds a
// `_kswq_t` by cloning `query` and `mat` into Vecs, but the SW kernel only needs slices —
// these variants accept the inputs by reference and avoid two heap allocations per pair lane.
// `q_max` is the per-pair max score from `mat` (i.e. ksw_qinit's `q->max`), needed for the
// score2 radius computation. `m` is the alphabet size (5 for DNA + N).
#[inline]
pub(crate) fn ksw_u8_slices(
    query: &[u8],
    m: i32,
    mat: &[i8],
    q_max: u8,
    tlen: i32,
    target: &[u8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    xtra: i32,
) -> kswr_t {
    #[cfg(target_arch = "x86_64")]
    {
        // SSE2 striped Farrar kernel — matches C++ ksw_u8 byte-for-byte.
        return ksw_u8_sse2::run(
            query, m, mat, q_max, tlen, target, o_del, e_del, o_ins, e_ins, xtra,
        );
    }
    #[cfg(not(target_arch = "x86_64"))]
    ksw_u8_slices_scalar(
        query, m, mat, q_max, tlen, target, o_del, e_del, o_ins, e_ins, xtra,
    )
}

#[cfg(not(target_arch = "x86_64"))]
fn ksw_u8_slices_scalar(
    query: &[u8],
    m: i32,
    mat: &[i8],
    q_max: u8,
    tlen: i32,
    target: &[u8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    xtra: i32,
) -> kswr_t {
    let tlen_usize = tlen as usize;
    let endsc_in = if (xtra & crate::bwa_mem2::ksw::KSW_XSTOP) != 0 {
        xtra & 0xffff
    } else {
        i32::MAX
    };
    let (mut r, row_maxes) = local_align_affine_with_rows_endsc(
        query,
        &target[..tlen_usize],
        m,
        mat,
        o_del,
        e_del,
        o_ins,
        e_ins,
        endsc_in,
    );
    if r.score > 255 {
        r.score = 255;
    }
    let minsc = if (xtra & crate::bwa_mem2::ksw::KSW_XSUBO) != 0 {
        xtra & 0xffff
    } else {
        0x10000
    };
    let endsc = if (xtra & crate::bwa_mem2::ksw::KSW_XSTOP) != 0 {
        xtra & 0xffff
    } else {
        0x10000
    };
    if r.score != 255 && r.score >= minsc {
        let radius = (r.score + i32::from(q_max).saturating_sub(1)) / i32::from(q_max.max(1));
        let low = r.te - radius;
        let high = r.te + radius;
        let b = build_score_peaks(&row_maxes, minsc);
        let mut alt = G_DEFR;
        for &(imax, te) in &b {
            if (te < low || te > high) && imax > alt.score2 {
                alt.score2 = imax.min(255);
                alt.te2 = te;
            }
        }
        r.score2 = alt.score2;
        r.te2 = alt.te2;
    }
    if r.score >= endsc && endsc != 0x10000 {
        return r;
    }
    r
}

#[inline]
pub(crate) fn ksw_i16_slices(
    query: &[u8],
    m: i32,
    mat: &[i8],
    q_max: u8,
    tlen: i32,
    target: &[u8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    xtra: i32,
) -> kswr_t {
    #[cfg(target_arch = "x86_64")]
    {
        return ksw_i16_sse2::run(
            query, m, mat, q_max, tlen, target, o_del, e_del, o_ins, e_ins, xtra,
        );
    }
    #[cfg(not(target_arch = "x86_64"))]
    ksw_i16_slices_scalar(
        query, m, mat, q_max, tlen, target, o_del, e_del, o_ins, e_ins, xtra,
    )
}

#[cfg(not(target_arch = "x86_64"))]
fn ksw_i16_slices_scalar(
    query: &[u8],
    m: i32,
    mat: &[i8],
    q_max: u8,
    tlen: i32,
    target: &[u8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    xtra: i32,
) -> kswr_t {
    let tlen_usize = tlen as usize;
    let endsc_in = if (xtra & crate::bwa_mem2::ksw::KSW_XSTOP) != 0 {
        xtra & 0xffff
    } else {
        i32::MAX
    };
    let (mut r, row_maxes) = local_align_affine_with_rows_endsc(
        query,
        &target[..tlen_usize],
        m,
        mat,
        o_del,
        e_del,
        o_ins,
        e_ins,
        endsc_in,
    );
    let minsc = if (xtra & crate::bwa_mem2::ksw::KSW_XSUBO) != 0 {
        xtra & 0xffff
    } else {
        0x10000
    };
    let endsc = if (xtra & crate::bwa_mem2::ksw::KSW_XSTOP) != 0 {
        xtra & 0xffff
    } else {
        0x10000
    };
    if r.score >= minsc {
        let radius = (r.score + i32::from(q_max).saturating_sub(1)) / i32::from(q_max.max(1));
        let low = r.te - radius;
        let high = r.te + radius;
        let b = build_score_peaks(&row_maxes, minsc);
        let mut best2 = -1;
        let mut te2 = -1;
        for &(imax, te) in &b {
            if (te < low || te > high) && imax > best2 {
                best2 = imax;
                te2 = te;
            }
        }
        r.score2 = best2;
        r.te2 = te2;
    }
    if r.score >= endsc && endsc != 0x10000 {
        return r;
    }
    r
}

// Compute `q->max` from a scoring matrix (mirrors ksw_qinit's mdiff calculation).
pub(crate) fn ksw_qmax(m: i32, mat: &[i8]) -> u8 {
    let mut mdiff = i8::MIN;
    for &score in mat.iter().take((m * m) as usize) {
        if score > mdiff {
            mdiff = score;
        }
    }
    mdiff as u8
}

// Build the `b[]` peak list of C++ ksw_u8/i16 (ksw.cpp:198-205).
// For each row i with imax >= minsc:
//   - if b is empty or the previous entry's stored i is NOT exactly i-1: append (imax, i).
//   - else if imax is strictly greater than the previous entry's stored imax: update last to (imax, i).
//   - else: do nothing (leaves the entry's i frozen at the position of the running max).
// This collapses contiguous runs of qualifying rows into single peaks. Critical for matching
// C++'s score2/te2 behavior — Rust must NOT scan all `row_maxes` (that finds extra candidates).
fn build_score_peaks(row_maxes: &[i32], minsc: i32) -> Vec<(i32, i32)> {
    let mut b: Vec<(i32, i32)> = Vec::new();
    for (i, &imax) in row_maxes.iter().enumerate() {
        if imax < minsc {
            continue;
        }
        let i = i as i32;
        if let Some(last) = b.last_mut() {
            if last.1 + 1 == i {
                if last.0 < imax {
                    *last = (imax, i);
                }
            } else {
                b.push((imax, i));
            }
        } else {
            b.push((imax, i));
        }
    }
    b
}

#[doc = "Original function: revseq:340"]
#[inline]
pub fn revseq(l: i32, s: &mut [u8]) {
    s[..l as usize].reverse();
}

fn local_align_affine(
    query: &[u8],
    target: &[u8],
    m: i32,
    mat: &[i8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
) -> kswr_t {
    local_align_affine_with_rows(query, target, m, mat, o_del, e_del, o_ins, e_ins).0
}

fn local_align_affine_with_rows(
    query: &[u8],
    target: &[u8],
    m: i32,
    mat: &[i8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
) -> (kswr_t, Vec<i32>) {
    local_align_affine_with_rows_endsc(query, target, m, mat, o_del, e_del, o_ins, e_ins, i32::MAX)
}

fn local_align_affine_with_rows_endsc(
    query: &[u8],
    target: &[u8],
    m: i32,
    mat: &[i8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    endsc: i32,
) -> (kswr_t, Vec<i32>) {
    let qlen = query.len();
    let tlen = target.len();
    let m_usize = m as usize;
    let mut h_prev = vec![0_i32; qlen + 1];
    let mut h_curr = vec![0_i32; qlen + 1];
    let mut e_prev = vec![MINUS_INF; qlen + 1];
    let mut e_curr = vec![MINUS_INF; qlen + 1];
    let mut row_maxes = vec![0_i32; tlen];
    let mut best = G_DEFR;

    let oe_del = o_del + e_del;
    let oe_ins = o_ins + e_ins;
    // Inner-loop indexing uses get_unchecked since j∈1..=qlen and all DP buffers are
    // exactly qlen+1 — bounds checks dominated this hot path (~25% perf win).
    for i in 1..=tlen {
        let mut row_max = 0_i32;
        let target_base = usize::from(target[i - 1]);
        let mat_row = &mat[target_base * m_usize..target_base * m_usize + m_usize];
        let mut f_left = MINUS_INF;
        h_curr[0] = 0;
        e_curr[0] = MINUS_INF;
        // Snapshot the relevant subslices so the optimizer sees stable bounds.
        let qslice = &query[..qlen];
        let hp = &h_prev[..=qlen];
        let ep = &e_prev[..=qlen];
        // SAFETY: qlen+1 elements, j∈1..=qlen, j-1∈0..qlen.
        for j in 1..=qlen {
            unsafe {
                let score =
                    i32::from(*mat_row.get_unchecked(usize::from(*qslice.get_unchecked(j - 1))));
                let e_ij = (*hp.get_unchecked(j) - oe_del).max(*ep.get_unchecked(j) - e_del);
                let f_ij = (*h_curr.get_unchecked(j - 1) - oe_ins).max(f_left - e_ins);
                let diag = *hp.get_unchecked(j - 1) + score;
                let h_ij = 0.max(diag.max(e_ij.max(f_ij)));
                *e_curr.get_unchecked_mut(j) = e_ij;
                f_left = f_ij;
                *h_curr.get_unchecked_mut(j) = h_ij;
                if h_ij > row_max {
                    row_max = h_ij;
                }
                if h_ij > best.score {
                    best.score = h_ij;
                    best.te = (i - 1) as i32;
                    best.qe = (j - 1) as i32;
                }
            }
        }
        row_maxes[i - 1] = row_max;
        if best.score >= endsc {
            // Match SIMD per-row KSW_XSTOP early termination.
            for k in i..tlen {
                row_maxes[k] = 0;
            }
            break;
        }
        std::mem::swap(&mut h_prev, &mut h_curr);
        std::mem::swap(&mut e_prev, &mut e_curr);
    }
    (best, row_maxes)
}

/// Aligning two sequences.
///
/// When `xtra==0`, `ksw_align()` uses a signed two-byte integer to store a score and only
/// finds the best score and the end positions. The 2nd best score or the start positions
/// are not attempted. The default behavior can be tuned by setting `KSW_X*` flags:
///
///   `KSW_XBYTE`:  use an unsigned byte to store a score. If overflow occurs,
///                 `kswr_t::score` will be set to 255.
///   `KSW_XSUBO`:  track the 2nd best score and the ending position on the target
///                 if the 2nd best is higher than `(xtra&0xffff)`.
///   `KSW_XSTOP`:  stop if the maximum score is above `(xtra&0xffff)`.
///   `KSW_XSTART`: find the start positions.
///
/// When `*qry==NULL`, `ksw_align()` will compute and allocate the query profile and on
/// return `*qry` will point to the profile. If one query is aligned against multiple
/// target sequences, `*qry` should be set to NULL during the first call and freed after
/// the last call. Note that `qry` can equal 0; in that case the query profile is
/// deallocated in `ksw_align()`.
#[doc = "Original function: ksw_align2:347"]
pub fn ksw_align2(
    qlen: i32,
    query: &mut [u8],
    tlen: i32,
    target: &mut [u8],
    m: i32,
    mat: &[i8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    xtra: i32,
    qry: Option<&mut Option<_kswq_t>>,
) -> kswr_t {
    let size = if (xtra & KSW_XBYTE) != 0 { 1 } else { 2 };
    let qlen_us = qlen as usize;
    let tlen_us = tlen as usize;
    let owned_q = if qry.is_none() {
        Some(ksw_qinit(size, qlen, &query[..qlen_us], m, mat))
    } else {
        None
    };
    let q_ref: &_kswq_t = if let Some(qry_slot) = qry {
        if qry_slot.is_none() {
            *qry_slot = Some(ksw_qinit(size, qlen, &query[..qlen_us], m, mat));
        }
        qry_slot.as_ref().expect("qry slot initialized")
    } else {
        owned_q.as_ref().expect("owned query profile")
    };

    let func_result = if q_ref.size == 2 {
        ksw_i16(q_ref, tlen, &target[..tlen_us], o_del, e_del, o_ins, e_ins, xtra)
    } else {
        ksw_u8(q_ref, tlen, &target[..tlen_us], o_del, e_del, o_ins, e_ins, xtra)
    };
    let mut r = func_result;
    if (xtra & KSW_XSTART) == 0 || ((xtra & KSW_XSUBO) != 0 && r.score < (xtra & 0xffff)) {
        return r;
    }

    // +1 because qe/te points to the exact end, not the position after the end
    revseq(r.qe + 1, query);
    revseq(r.te + 1, target);
    let qe1_us = (r.qe + 1) as usize;
    let te1_us = (r.te + 1) as usize;
    let q2 = ksw_qinit(size, r.qe + 1, &query[..qe1_us], m, mat);
    let rr = if q2.size == 2 {
        ksw_i16(
            &q2,
            r.te + 1,
            &target[..te1_us],
            o_del,
            e_del,
            o_ins,
            e_ins,
            KSW_XSTOP | r.score,
        )
    } else {
        ksw_u8(
            &q2,
            r.te + 1,
            &target[..te1_us],
            o_del,
            e_del,
            o_ins,
            e_ins,
            KSW_XSTOP | r.score,
        )
    };
    revseq(r.qe + 1, query);
    revseq(r.te + 1, target);
    if r.score == rr.score {
        r.tb = r.te - rr.te;
        r.qb = r.qe - rr.qe;
    }
    r
}

#[doc = "Original function: ksw_align2_orig_bak:383"]
pub fn ksw_align2_orig_bak(
    qlen: i32,
    query: &mut [u8],
    tlen: i32,
    target: &mut [u8],
    m: i32,
    mat: &[i8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    xtra: i32,
    qry: Option<&mut Option<_kswq_t>>,
) -> kswr_t {
    ksw_align2(
        qlen, query, tlen, target, m, mat, o_del, e_del, o_ins, e_ins, xtra, qry,
    )
}

#[doc = "Original function: ksw_align:419"]
pub fn ksw_align(
    qlen: i32,
    query: &mut [u8],
    tlen: i32,
    target: &mut [u8],
    m: i32,
    mat: &[i8],
    gapo: i32,
    gape: i32,
    xtra: i32,
    qry: Option<&mut Option<_kswq_t>>,
) -> kswr_t {
    ksw_align2(
        qlen, query, tlen, target, m, mat, gapo, gape, gapo, gape, xtra, qry,
    )
}

// ********************
// *** SW extension ***
// ********************
//
// Extend alignment from the seed (h0). Aligns query and target assuming their upstream
// sequences (not provided) have been aligned with score h0. On return [0,*qle) on the
// query and [0,*tle) on the target are aligned together. If *gscore>=0, *gscore keeps the
// best score for an end-to-end alignment of the query, and *gtle is its target position;
// these help the caller decide between end-to-end and partial hits.
#[doc = "Original function: ksw_extend2:432"]
pub fn ksw_extend2(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i32,
    mat: &[i8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    mut w: i32,
    end_bonus: i32,
    zdrop: i32,
    h0: i32,
    qle: Option<&mut i32>,
    tle: Option<&mut i32>,
    gtle: Option<&mut i32>,
    gscore: Option<&mut i32>,
    max_off: Option<&mut i32>,
) -> i32 {
    assert!(h0 > 0);
    let qlen_usize = qlen as usize;
    let tlen_usize = tlen as usize;
    let m_usize = m as usize;
    let oe_del = o_del + e_del;
    let oe_ins = o_ins + e_ins;
    // allocate memory (eh = score array, qp = query profile)
    let (mut eh, mut qp) = KSW_EXTEND_SCRATCH.with(|c| std::mem::take(&mut *c.borrow_mut()));
    eh.clear();
    eh.resize(qlen_usize + 1, eh_t::default());
    qp.clear();
    qp.resize(qlen_usize * m_usize, 0);

    // generate the query profile
    let mut idx = 0_usize;
    for k in 0..m_usize {
        let p = &mat[k * m_usize..(k + 1) * m_usize];
        for &q in query.iter().take(qlen_usize) {
            qp[idx] = p[usize::from(q)];
            idx += 1;
        }
    }

    // fill the first row
    eh[0].h = h0;
    if qlen_usize >= 1 {
        eh[1].h = if h0 > oe_ins { h0 - oe_ins } else { 0 };
    }
    let mut j = 2_i32;
    while j <= qlen && eh[(j - 1) as usize].h > e_ins {
        let ju = j as usize;
        eh[ju].h = eh[(j - 1) as usize].h - e_ins;
        j += 1;
    }

    // adjust w if it is too large
    let mut max_sc = 0_i32;
    for &score in mat.iter().take(m_usize * m_usize) {
        max_sc = max_sc.max(i32::from(score));
    }
    let max_ins = (((qlen * max_sc + end_bonus - o_ins) as f64) / f64::from(e_ins) + 1.0) as i32;
    w = w.min(max_ins.max(1));
    let max_del = (((qlen * max_sc + end_bonus - o_del) as f64) / f64::from(e_del) + 1.0) as i32;
    w = w.min(max_del.max(1)); // TODO: is this necessary?

    // DP loop
    let mut best = h0;
    let mut best_i = -1;
    let mut best_j = -1;
    let mut best_ie = -1;
    let mut best_gscore = -1;
    let mut best_off = 0;
    let mut beg = 0_i32;
    let mut end = qlen;

    for i in 0..tlen_usize {
        let mut f = 0_i32;
        let mut h1;
        let mut row_best = 0_i32;
        let mut row_best_j = -1_i32;
        let qprof = &qp
            [usize::from(target[i]) * qlen_usize..usize::from(target[i]) * qlen_usize + qlen_usize];

        let i_i32 = i as i32;
        // apply the band and the constraint (if provided)
        if beg < i_i32 - w {
            beg = i_i32 - w;
        }
        if end > i_i32 + w + 1 {
            end = i_i32 + w + 1;
        }
        if end > qlen {
            end = qlen;
        }

        // compute the first column
        if beg == 0 {
            h1 = (h0 - (o_del + e_del * (i_i32 + 1))).max(0);
        } else {
            h1 = 0;
        }

        let mut j_i32 = beg;
        while j_i32 < end {
            // At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1).
            // Similar to SSE2-SW, cells are computed in the following order:
            //   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
            //   E(i+1,j) = max{H(i,j)-gapo, E(i,j)} - gape
            //   F(i,j+1) = max{H(i,j)-gapo, F(i,j)} - gape
            let j_usize = j_i32 as usize;
            let p = &mut eh[j_usize];
            let mut mm = p.h; // get H(i-1,j-1)
            let mut e = p.e; // get E(i-1,j)
            p.h = h1; // set H(i,j-1) for the next row
            // separating H and M to disallow a cigar like "100M3I3D20M"
            mm = if mm != 0 {
                mm + i32::from(qprof[j_usize])
            } else {
                0
            };
            // e and f are guaranteed to be non-negative, so h>=0 even if M<0
            let h = mm.max(e).max(f);
            h1 = h; // save H(i,j) to h1 for the next column
            // record the position where max score is achieved (m is stored at eh[mj+1])
            if h >= row_best {
                row_best_j = j_i32;
                row_best = h;
            }
            let mut t = (mm - oe_del).max(0);
            e = (e - e_del).max(t); // computed E(i+1,j)
            p.e = e; // save E(i+1,j) for the next row
            t = (mm - oe_ins).max(0);
            f = (f - e_ins).max(t); // computed F(i,j+1)
            j_i32 += 1;
        }

        let end_us = end as usize;
        eh[end_us].h = h1;
        eh[end_us].e = 0;
        if j_i32 == qlen {
            // C++ ksw.cpp:504-506 uses `gscore > h1 ? max_ie : i` ternary which updates on TIE.
            if h1 >= best_gscore {
                best_ie = i_i32;
                best_gscore = h1;
            }
        }
        if row_best == 0 {
            break;
        }
        if row_best > best {
            best = row_best;
            best_i = i_i32;
            best_j = row_best_j;
            best_off = best_off.max((row_best_j - i_i32).abs());
        } else if zdrop > 0 {
            let lhs = i_i32 - best_i;
            let rhs = row_best_j - best_j;
            if lhs > rhs {
                if best - row_best - (lhs - rhs) * e_del > zdrop {
                    break;
                }
            } else if best - row_best - (rhs - lhs) * e_ins > zdrop {
                break;
            }
        }

        // update beg and end for the next round
        let mut new_beg = beg;
        while new_beg < end {
            let p = &eh[new_beg as usize];
            if p.h != 0 || p.e != 0 {
                break;
            }
            new_beg += 1;
        }
        beg = new_beg;
        let mut new_end = end;
        while new_end >= beg {
            let p = &eh[new_end as usize];
            if p.h != 0 || p.e != 0 {
                break;
            }
            new_end -= 1;
        }
        end = (new_end + 2).min(qlen);
        // beg = 0; end = qlen; // uncomment this line for debugging
    }

    if let Some(qle) = qle {
        *qle = best_j + 1;
    }
    if let Some(tle) = tle {
        *tle = best_i + 1;
    }
    if let Some(gtle) = gtle {
        *gtle = best_ie + 1;
    }
    if let Some(gscore) = gscore {
        *gscore = best_gscore;
    }
    if let Some(max_off) = max_off {
        *max_off = best_off;
    }
    KSW_EXTEND_SCRATCH.with(|c| *c.borrow_mut() = (eh, qp));
    best
}

#[doc = "Original function: ksw_extend:535"]
pub fn ksw_extend(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i32,
    mat: &[i8],
    gapo: i32,
    gape: i32,
    w: i32,
    end_bonus: i32,
    zdrop: i32,
    h0: i32,
    qle: Option<&mut i32>,
    tle: Option<&mut i32>,
    gtle: Option<&mut i32>,
    gscore: Option<&mut i32>,
    max_off: Option<&mut i32>,
) -> i32 {
    ksw_extend2(
        qlen, query, tlen, target, m, mat, gapo, gape, gapo, gape, w, end_bonus, zdrop, h0, qle,
        tle, gtle, gscore, max_off,
    )
}

#[doc = "Original function: push_cigar:546"]
#[inline]
pub fn push_cigar(n_cigar: &mut i32, _m_cigar: &mut i32, cigar: &mut Vec<u32>, op: i32, len: i32) {
    let idx = (*n_cigar - 1) as usize;
    if *n_cigar == 0 || op != (cigar[idx] & 0xf) as i32 {
        cigar.push(((len as u32) << 4) | (op as u32));
        *n_cigar += 1;
    } else {
        cigar[idx] += (len as u32) << 4;
    }
}

// ********************
// * Global alignment *
// ********************
//
// Banded global alignment with optional CIGAR traceback. `w` is the band width; when
// `n_cigar_` and `cigar_` are both provided, the backtrack matrix `z` is allocated
// (z[i,j] keeps h for the current cell and e/f for the next cell, packed as f<<4|e<<2|h
// — in principle we could halve the memory but the backtrack would be slightly more
// complex). Returns the alignment score.
#[doc = "Original function: ksw_global2:558"]
pub fn ksw_global2(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i32,
    mat: &[i8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    w: i32,
    n_cigar_: Option<&mut i32>,
    cigar_: Option<&mut Vec<u32>>,
) -> i32 {
    let qlen_usize = qlen as usize;
    let tlen_usize = tlen as usize;
    let m_usize = m as usize;
    let oe_del = o_del + e_del;
    let oe_ins = o_ins + e_ins;
    let mut n_cigar_opt = n_cigar_;
    if let Some(n_cigar) = n_cigar_opt.as_deref_mut() {
        *n_cigar = 0;
    }

    // allocate memory; n_col is the maximum #columns of the backtrack matrix
    let n_col = qlen.min(2 * w + 1);
    let n_col_usize = n_col as usize;
    KSW_GLOBAL_SCRATCH.with(|cell| {
        let mut buf = cell.borrow_mut();
        let (z, qp, eh) = &mut *buf;
        let need_z = n_cigar_opt.is_some() && cigar_.is_some();
        if need_z {
            let z_len = n_col_usize * tlen_usize;
            z.clear();
            z.resize(z_len, 0);
        } else {
            z.clear();
        }
        // Skip the resize-with-zero step: every qp byte is overwritten by extend below.
        qp.clear();
        qp.reserve(qlen_usize * m_usize);
        eh.clear();
        eh.resize(qlen_usize + 1, eh_t::default());

        // generate the query profile
        for k in 0..m_usize {
            let p = &mat[k * m_usize..(k + 1) * m_usize];
            qp.extend(query.iter().take(qlen_usize).map(|&q| p[usize::from(q)]));
        }

        // fill the first row
        eh[0].h = 0;
        eh[0].e = MINUS_INF;
        let mut j = 1_i32;
        while j <= qlen && j <= w {
            let ju = j as usize;
            eh[ju].h = -(o_ins + e_ins * j);
            eh[ju].e = MINUS_INF;
            j += 1;
        }
        // everything is -inf outside the band
        while j <= qlen {
            let ju = j as usize;
            eh[ju].h = MINUS_INF;
            eh[ju].e = MINUS_INF;
            j += 1;
        }

        // DP loop — target sequence is in the outer loop
        for i in 0..tlen {
            let mut f = MINUS_INF;
            let beg = if i > w { i - w } else { 0 };
            // only loop through [beg,end) of the query sequence
            let end = if i + w + 1 < qlen { i + w + 1 } else { qlen };
            let mut h1 = if beg == 0 {
                -(o_del + e_del * (i + 1))
            } else {
                MINUS_INF
            };
            let i_us = i as usize;
            let q = &qp[usize::from(target[i_us]) * qlen_usize..];
            // Inner-DP bounds-check elimination: j ∈ [beg, end) ⊆ [0, qlen), eh has qlen+1 entries
            // and q has at least qlen entries, so j and j as usize are in range.
            let beg_us = beg as usize;
            let end_us = end as usize;
            if !z.is_empty() {
                let zi_base = i_us * n_col_usize;
                for j in beg_us..end_us {
                    // At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1).
                    // Cells are computed in the following order:
                    //   M(i,j)   = H(i-1,j-1) + S(i,j)
                    //   H(i,j)   = max{M(i,j), E(i,j), F(i,j)}
                    //   E(i+1,j) = max{M(i,j)-gapo, E(i,j)} - gape
                    //   F(i,j+1) = max{M(i,j)-gapo, F(i,j)} - gape
                    // We have to separate M(i,j); otherwise the direction may not be recorded correctly.
                    // However, a CIGAR like "10M3I3D10M" allowed by local() is disallowed by global().
                    // Such a CIGAR may occur, in theory, if mismatch_penalty > 2*gap_ext_penalty + 2*gap_open_penalty/k.
                    // In practice, this should happen very rarely given a reasonable scoring system.
                    unsafe {
                        let p = eh.get_unchecked_mut(j);
                        let mut m_score = p.h;
                        let mut e = p.e;
                        p.h = h1;
                        m_score += i32::from(*q.get_unchecked(j));
                        // d = direction
                        let mut d = if m_score >= e { 0_u8 } else { 1_u8 };
                        let mut h = if m_score >= e { m_score } else { e };
                        d = if h >= f { d } else { 2_u8 };
                        h = if h >= f { h } else { f };
                        h1 = h;
                        let t = m_score - oe_del;
                        e -= e_del;
                        d |= if e > t { 1_u8 << 2 } else { 0 };
                        e = if e > t { e } else { t };
                        p.e = e;
                        let t = m_score - oe_ins;
                        f -= e_ins;
                        // if we want to halve the memory, use one bit only, instead of two
                        d |= if f > t { 2_u8 << 4 } else { 0 };
                        f = if f > t { f } else { t };
                        // z[i,j] keeps h for the current cell and e/f for the next cell
                        *z.get_unchecked_mut(zi_base + j - beg_us) = d;
                    }
                }
            } else {
                for j in beg_us..end_us {
                    unsafe {
                        let p = eh.get_unchecked_mut(j);
                        let mut m_score = p.h;
                        let mut e = p.e;
                        p.h = h1;
                        m_score += i32::from(*q.get_unchecked(j));
                        let mut h = if m_score >= e { m_score } else { e };
                        h = if h >= f { h } else { f };
                        h1 = h;
                        let t = m_score - oe_del;
                        e -= e_del;
                        e = if e > t { e } else { t };
                        p.e = e;
                        let t = m_score - oe_ins;
                        f -= e_ins;
                        f = if f > t { f } else { t };
                    }
                }
            }
            let endu = end as usize;
            eh[endu].h = h1;
            eh[endu].e = MINUS_INF;
        }

        let score = eh[qlen_usize].h;
        if !z.is_empty() {
            // backtrack
            let mut n_cigar = 0_i32;
            let mut m_cigar = 0_i32;
            let mut which = 0_u8;
            // Push directly into the caller's Vec (assumed empty) to avoid a local cigar Vec.
            let cigar = cigar_.expect("cigar slot required when traceback enabled");
            cigar.clear();
            // (i,k) points to the last cell
            let mut i = tlen - 1;
            let mut k = (if i + w + 1 < qlen { i + w + 1 } else { qlen }) - 1;
            while i >= 0 && k >= 0 {
                let zi = (i as usize) * n_col_usize
                    + (k - if i > w { i - w } else { 0 }) as usize;
                which = (z[zi] >> (which << 1)) & 3;
                if which == 0 {
                    push_cigar(&mut n_cigar, &mut m_cigar, cigar, 0, 1);
                    i -= 1;
                    k -= 1;
                } else if which == 1 {
                    push_cigar(&mut n_cigar, &mut m_cigar, cigar, 2, 1);
                    i -= 1;
                } else {
                    push_cigar(&mut n_cigar, &mut m_cigar, cigar, 1, 1);
                    k -= 1;
                }
            }
            if i >= 0 {
                push_cigar(&mut n_cigar, &mut m_cigar, cigar, 2, i + 1);
            }
            if k >= 0 {
                push_cigar(&mut n_cigar, &mut m_cigar, cigar, 1, k + 1);
            }
            // reverse CIGAR
            cigar.reverse();
            if let Some(n_cigar_slot) = n_cigar_opt {
                *n_cigar_slot = n_cigar;
            }
        }
        score
    })
}

#[doc = "Original function: ksw_global:670"]
pub fn ksw_global(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i32,
    mat: &[i8],
    gapo: i32,
    gape: i32,
    w: i32,
    n_cigar_: Option<&mut i32>,
    cigar_: Option<&mut Vec<u32>>,
) -> i32 {
    ksw_global2(
        qlen, query, tlen, target, m, mat, gapo, gape, gapo, gape, w, n_cigar_, cigar_,
    )
}

#[doc = "Original function: main:706"]
pub fn main(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("main")
}

#[cfg(test)]
mod tests {
    use super::{
        ksw_align, ksw_align2, ksw_extend, ksw_extend2, ksw_global, ksw_global2, ksw_i16,
        ksw_qinit, ksw_u8, push_cigar,
    };
    use crate::bwa_mem2::ksw::{KSW_XSTART, KSW_XSUBO};

    fn score_mat(match_score: i8, mismatch_penalty: i8) -> [i8; 25] {
        let mut mat = [-1_i8; 25];
        let mut k = 0;
        for i in 0..4 {
            for j in 0..4 {
                mat[k] = if i == j {
                    match_score
                } else {
                    -mismatch_penalty
                };
                k += 1;
            }
            mat[k] = -1;
            k += 1;
        }
        mat
    }

    #[test]
    fn build_score_peaks_collapses_consecutive_runs() {
        use super::build_score_peaks;
        // C++ ksw.cpp:198-205 b[] semantics: consecutive runs collapse to peaks.
        // Run [50, 60, 55, 70] (i=0..3) collapses; gap at i=4 starts new run.
        // Expected: when running, we update entry's i only when imax strictly increases.
        // i=0: append (50, 0) → b=[(50,0)]
        // i=1: 0+1==1, 50<60 → update to (60, 1) → b=[(60,1)]
        // i=2: 1+1==2, 60<55 false → no update → b=[(60,1)]
        // i=3: 1+1==2 != 3 → append (70, 3) → b=[(60,1),(70,3)]
        // i=4: imax=0 < minsc=1 → skip
        // i=5: append (40, 5) → b=[(60,1),(70,3),(40,5)]
        let row_maxes = [50, 60, 55, 70, 0, 40];
        let b = build_score_peaks(&row_maxes, 1);
        assert_eq!(b, vec![(60, 1), (70, 3), (40, 5)]);
    }

    #[test]
    fn build_score_peaks_differs_from_full_row_scan() {
        use super::build_score_peaks;
        // The whole point: Rust used to scan all `row_maxes` directly. C++ b[] keeps
        // only run-peaks. They diverge when a run has an interior dip (lower than the
        // running max but still >= minsc), or when scores rise then fall.
        // row_maxes [10, 20, 15, 25] with minsc=8:
        //   Full scan: yields {10, 20, 15, 25} (4 candidates)
        //   b[] peaks: i=0 append (10,0); i=1 update (20,1); i=2 20<15 false no update;
        //              i=3 last.i=1 != 2, append (25,3). Result: [(20,1), (25,3)] (2 candidates)
        let row_maxes = [10, 20, 15, 25];
        let b = build_score_peaks(&row_maxes, 8);
        assert_eq!(b, vec![(20, 1), (25, 3)]);
    }

    #[test]
    fn build_score_peaks_filters_below_minsc() {
        use super::build_score_peaks;
        // Below-minsc rows are filtered entirely; runs span only qualifying rows.
        let row_maxes = [30, 40, 5, 50];
        let b = build_score_peaks(&row_maxes, 10);
        // i=0: append (30, 0)
        // i=1: 0+1==1, 30<40 → update to (40, 1)
        // i=2: 5 < 10, skip
        // i=3: 1+1==2 != 3 → append (50, 3)
        assert_eq!(b, vec![(40, 1), (50, 3)]);
    }

    #[test]
    fn push_cigar_merges_adjacent_same_ops() {
        let mut n_cigar = 0;
        let mut m_cigar = 0;
        let mut cigar = Vec::new();
        push_cigar(&mut n_cigar, &mut m_cigar, &mut cigar, 0, 2);
        push_cigar(&mut n_cigar, &mut m_cigar, &mut cigar, 0, 3);
        push_cigar(&mut n_cigar, &mut m_cigar, &mut cigar, 1, 1);
        assert_eq!(n_cigar, 2);
        assert_eq!(cigar, vec![(5 << 4), (1 << 4) | 1]);
    }

    #[test]
    fn ksw_global2_aligns_exact_match_with_single_m_cigar() {
        let mat = score_mat(2, 3);
        let query = [0_u8, 1, 2, 3];
        let target = [0_u8, 1, 2, 3];
        let mut n_cigar = 0;
        let mut cigar = Vec::new();
        let score = ksw_global2(
            4,
            &query,
            4,
            &target,
            5,
            &mat,
            5,
            1,
            5,
            1,
            0,
            Some(&mut n_cigar),
            Some(&mut cigar),
        );
        assert_eq!(score, 8);
        assert_eq!(n_cigar, 1);
        assert_eq!(cigar, vec![4 << 4]);
    }

    #[test]
    fn ksw_global_wrapper_uses_symmetric_gap_penalties() {
        let mat = score_mat(1, 4);
        let query = [0_u8, 1, 2];
        let target = [0_u8, 1, 2];
        let score = ksw_global(3, &query, 3, &target, 5, &mat, 6, 1, 0, None, None);
        assert_eq!(score, 3);
    }

    #[test]
    fn ksw_align2_finds_local_alignment_endpoints_and_starts() {
        let mat = score_mat(2, 3);
        let mut query = [0_u8, 1, 2, 3];
        let mut target = [3_u8, 3, 0, 1, 2, 3, 0];
        let r = ksw_align2(
            4,
            &mut query,
            7,
            &mut target,
            5,
            &mat,
            5,
            1,
            5,
            1,
            KSW_XSTART,
            None,
        );
        assert_eq!(r.score, 8);
        assert_eq!(r.qe, 3);
        assert_eq!(r.te, 5);
        assert_eq!(r.qb, 0);
        assert_eq!(r.tb, 2);
    }

    #[test]
    fn ksw_align_wrapper_matches_symmetric_local_alignment() {
        let mat = score_mat(2, 3);
        let mut query = [0_u8, 1, 2];
        let mut target = [3_u8, 0, 1, 2, 3];
        let r = ksw_align(3, &mut query, 5, &mut target, 5, &mat, 5, 1, 0, None);
        assert_eq!(r.score, 6);
        assert_eq!(r.qe, 2);
        assert_eq!(r.te, 3);
    }

    #[test]
    fn ksw_qinit_tracks_query_profile_metadata() {
        let mat = score_mat(2, 3);
        let query = [0_u8, 1, 2, 3, 0];
        let q = ksw_qinit(1, 5, &query, 5, &mat);
        assert_eq!(q.qlen, 5);
        assert_eq!(q.slen, 1);
        assert_eq!(q.size, 1);
        assert_eq!(q.shift, 3);
        assert_eq!(q.mdiff, 5);
        assert_eq!(q.max, 2);
        assert_eq!(q.query, query);
        assert_eq!(q.m, 5);
        assert_eq!(&q.mat[..25], &mat);
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn sse2_ksw_u8_matches_scalar_on_simple_alignment() {
        let mat = score_mat(2, 3);
        let query = [0_u8, 1, 2, 3];
        let target = [3_u8, 0, 1, 2, 3, 3];
        let q8 = ksw_qinit(1, 4, &query, 5, &mat);
        let scalar = ksw_u8(&q8, 6, &target, 5, 1, 5, 1, 0);
        let simd = super::ksw_u8_sse2::run(&query, 5, &mat, q8.max, 6, &target, 5, 1, 5, 1, 0);
        assert_eq!(simd.score, scalar.score, "score");
        assert_eq!(simd.te, scalar.te, "te");
        assert_eq!(simd.qe, scalar.qe, "qe");
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn sse2_ksw_u8_matches_scalar_with_subo_search() {
        // Use a long enough query that there are no padding lanes (qlen=16, slen=1).
        // With qlen<16 the striped layout has phantom lanes that propagate scores via
        // slli into out-of-range positions; this is C++ ksw_u8's actual behavior, but it
        // makes the score2/te2 peak list diverge from the scalar (which only tracks real
        // cells). Test on aligned-length input to keep SIMD and scalar comparable.
        let mat = score_mat(2, 3);
        let query: [u8; 16] = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let target: Vec<u8> = (0..32).map(|i| (i & 3) as u8).collect();
        let q8 = ksw_qinit(1, 16, &query, 5, &mat);
        let xtra = KSW_XSUBO | 1;
        let scalar = ksw_u8(&q8, 32, &target, 5, 1, 5, 1, xtra);
        let simd = super::ksw_u8_sse2::run(&query, 5, &mat, q8.max, 32, &target, 5, 1, 5, 1, xtra);
        assert_eq!(simd.score, scalar.score, "score");
        assert_eq!(simd.te, scalar.te, "te");
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn sse2_ksw_i16_slices_reuses_scratch_and_matches_scalar() {
        let mat = score_mat(2, 3);
        let cases: [(&[u8], &[u8]); 2] = [
            (&[0, 1, 2, 3, 0, 1, 2, 3], &[3, 0, 1, 2, 3, 0, 1, 2, 3]),
            (&[0, 1, 2, 3], &[3, 3, 0, 1, 2, 3]),
        ];
        for (query, target) in cases {
            let q16 = ksw_qinit(2, query.len() as i32, query, 5, &mat);
            let scalar = ksw_i16(&q16, target.len() as i32, target, 5, 1, 5, 1, KSW_XSUBO | 1);
            let simd = super::ksw_i16_slices(
                query,
                5,
                &mat,
                q16.max,
                target.len() as i32,
                target,
                5,
                1,
                5,
                1,
                KSW_XSUBO | 1,
            );
            assert_eq!(simd.score, scalar.score, "score");
            assert_eq!(simd.te, scalar.te, "te");
            assert_eq!(simd.qe, scalar.qe, "qe");
            assert_eq!(simd.score2, scalar.score2, "score2");
            assert_eq!(simd.te2, scalar.te2, "te2");
        }
    }

    #[test]
    fn ksw_u8_and_i16_find_same_exact_local_alignment() {
        let mat = score_mat(2, 3);
        let query = [0_u8, 1, 2, 3];
        let target = [3_u8, 0, 1, 2, 3, 3];
        let q8 = ksw_qinit(1, 4, &query, 5, &mat);
        let q16 = ksw_qinit(2, 4, &query, 5, &mat);
        let r8 = ksw_u8(&q8, 6, &target, 5, 1, 5, 1, 0);
        let r16 = ksw_i16(&q16, 6, &target, 5, 1, 5, 1, 0);
        assert_eq!(r8.score, 8);
        assert_eq!(r8.qe, 3);
        assert_eq!(r8.te, 4);
        assert_eq!(r16.score, r8.score);
        assert_eq!(r16.qe, r8.qe);
        assert_eq!(r16.te, r8.te);
    }

    #[test]
    fn ksw_u8_and_i16_report_second_best_outside_primary_window() {
        let mat = score_mat(2, 3);
        let query = [0_u8, 1];
        let target = [0_u8, 1, 3, 3, 3, 3, 0, 1];
        let q8 = ksw_qinit(1, 2, &query, 5, &mat);
        let q16 = ksw_qinit(2, 2, &query, 5, &mat);
        let xtra = KSW_XSUBO | 1;
        let r8 = ksw_u8(&q8, 8, &target, 5, 1, 5, 1, xtra);
        let r16 = ksw_i16(&q16, 8, &target, 5, 1, 5, 1, xtra);
        assert_eq!(r8.score, 4);
        assert_eq!(r8.te, 1);
        assert_eq!(r8.score2, 4);
        assert!(r8.te2 >= 4);
        assert_eq!(r16.score, r8.score);
        assert_eq!(r16.te, r8.te);
        assert_eq!(r16.score2, r8.score2);
        assert_eq!(r16.te2, r8.te2);
    }

    #[test]
    fn ksw_extend2_reports_best_and_global_end_positions() {
        let mat = score_mat(2, 3);
        let query = [0_u8, 1, 2, 3];
        let target = [0_u8, 1, 2, 3, 3];
        let mut qle = -1;
        let mut tle = -1;
        let mut gtle = -1;
        let mut gscore = -1;
        let mut max_off = -1;
        let score = ksw_extend2(
            4,
            &query,
            5,
            &target,
            5,
            &mat,
            5,
            1,
            5,
            1,
            10,
            0,
            20,
            1,
            Some(&mut qle),
            Some(&mut tle),
            Some(&mut gtle),
            Some(&mut gscore),
            Some(&mut max_off),
        );
        assert_eq!(score, 9);
        assert_eq!(qle, 4);
        assert_eq!(tle, 4);
        assert_eq!(gtle, 4);
        assert_eq!(gscore, 9);
        assert_eq!(max_off, 0);
    }

    #[test]
    fn ksw_extend_wrapper_uses_symmetric_gap_penalties() {
        let mat = score_mat(2, 3);
        let query = [0_u8, 1, 2];
        let target = [0_u8, 1, 2, 3];
        let mut qle2 = -1;
        let mut tle2 = -1;
        let mut gtle2 = -1;
        let mut gscore2 = -1;
        let mut max_off2 = -1;
        let score2 = ksw_extend2(
            3,
            &query,
            4,
            &target,
            5,
            &mat,
            6,
            1,
            6,
            1,
            8,
            0,
            15,
            1,
            Some(&mut qle2),
            Some(&mut tle2),
            Some(&mut gtle2),
            Some(&mut gscore2),
            Some(&mut max_off2),
        );
        let mut qle = -1;
        let mut tle = -1;
        let mut gtle = -1;
        let mut gscore = -1;
        let mut max_off = -1;
        let score = ksw_extend(
            3,
            &query,
            4,
            &target,
            5,
            &mat,
            6,
            1,
            8,
            0,
            15,
            1,
            Some(&mut qle),
            Some(&mut tle),
            Some(&mut gtle),
            Some(&mut gscore),
            Some(&mut max_off),
        );
        assert_eq!(score, score2);
        assert_eq!(qle, qle2);
        assert_eq!(tle, tle2);
        assert_eq!(gtle, gtle2);
        assert_eq!(gscore, gscore2);
        assert_eq!(max_off, max_off2);
    }
}
