#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Port of `bwa-mem2/src/bwamem_pair.cpp`.

use crate::bwa_mem2::bandedswa::SeqPair;
use crate::bwa_mem2::bntseq::bntseq_t;
use crate::bwa_mem2::bwa::bseq1_t;
use crate::bwa_mem2::bwamem::{
    mem_aln2sam, mem_aln_t, mem_alnreg_t, mem_alnreg_v, mem_approx_mapq_se, mem_cache,
    mem_mark_primary_se, mem_opt_t, mem_pestat_t, mem_reg2aln, mem_reg2sam, mem_reorder_primary5,
    mem_sort_dedup_patch, nt4_cow, query_string_for_aln,
};
use crate::bwa_mem2::bwamem_extra::mem_gen_alt;
use crate::bwa_mem2::kstring::kstring_t;
use crate::bwa_mem2::ksw::{ksw_align2, kswr_t, KSW_XBYTE, KSW_XSTART, KSW_XSTOP, KSW_XSUBO};
use crate::bwa_mem2::kswv::kswv;
use crate::bwa_mem2::utils::{hash_64, pair64_t};

// --- bwamem_pair.cpp ---

const MIN_RATIO: f64 = 0.8;
const MIN_DIR_CNT: usize = 10;
const MIN_DIR_RATIO: f64 = 0.05;
const OUTLIER_BOUND: f64 = 2.0;
const MAPPING_BOUND: f64 = 3.0;
const MAX_STDDEV: f64 = 4.0;
const MEM_F_ALL: i32 = 0x8;
const MEM_F_NOPAIRING: i32 = 0x4;
const MEM_F_PRIMARY5: i32 = 0x800;

#[inline]
fn max_matesw_rescue_limit(opt: &mem_opt_t) -> usize {
    opt.max_matesw.max(0) as usize
}

#[inline]
fn take_kstring_boxed(str_: &mut kstring_t) -> String {
    let mut bytes = std::mem::take(&mut str_.s);
    bytes.truncate(str_.l);
    str_.l = 0;
    str_.m = 0;
    // SAFETY: every byte in `str_` was written by kput*/ksprintf helpers (ASCII output) or by
    // FASTQ sequence/quality bytes which the formats define as printable ASCII. Skipping the
    // O(n) UTF-8 validation pass eliminates ~hundreds of MB of byte scans on a 700K-read run.
    // Skip into_boxed_str() too — keeps Vec<u8> capacity (no shrink-to-fit alloc/memcpy per
    // SAM record). The .as_deref()/.as_str() readers downstream work uniformly with String.
    unsafe { String::from_utf8_unchecked(bytes) }
}

fn debug_trace_read(name: Option<&str>) -> bool {
    static TARGET: std::sync::OnceLock<Option<String>> = std::sync::OnceLock::new();
    let Some(target) = TARGET
        .get_or_init(|| std::env::var("BWA_MEM2_RS_TRACE_READ").ok())
        .as_deref()
    else {
        return false;
    };
    name == Some(target)
}

#[inline]
fn erfc_approx(x: f64) -> f64 {
    // Match libm's erfc precisely (used by upstream C++ via <math.h>); a low-precision
    // Chebyshev approximation here was producing ~1-unit differences in mem_pair's
    // pair-score quantization, which propagated to MAPQ divergences vs upstream.
    extern "C" {
        fn erfc(x: f64) -> f64;
    }
    unsafe { erfc(x) }
}

#[inline]
fn raw_mapq(diff: i32, a: i32) -> i32 {
    (6.02 * (diff as f64) / (a as f64) + 0.499) as i32
}

/// Infer the relative orientation and distance of two mapping positions.
///
/// # Arguments
/// * `l_pac` - packed-reference length (single-strand)
/// * `b1`, `b2` - reference coordinates of the two mates (in concatenated fwd|rev space)
/// * `dist` - written with the (always non-negative) distance between the mates
///
/// # Returns
/// Orientation code in `{0 = FF, 1 = FR, 2 = RF, 3 = RR}`.
#[doc = "Original function: mem_infer_dir:58"]
#[inline]
pub fn mem_infer_dir(l_pac: i64, b1: i64, b2: i64, dist: &mut i64) -> i32 {
    let r1 = b1 >= l_pac;
    let r2 = b2 >= l_pac;
    // p2 is the coordinate of read 2 on the read 1 strand
    let p2 = if r1 == r2 { b2 } else { (l_pac << 1) - 1 - b2 };
    *dist = if p2 > b1 { p2 - b1 } else { b1 - p2 };
    (if r1 == r2 { 0 } else { 1 }) ^ (if p2 > b1 { 0 } else { 3 })
}

/// Compute the suboptimal-alignment score used for MAPQ calibration.
///
/// Walks the secondary hits and returns the score of the first one whose query-side overlap
/// with the primary is significant (overlap length >= `opt.mask_level * min(overlap length)`).
/// Falls back to `opt.min_seed_len * opt.a` when no such secondary exists.
#[doc = "Original function: cal_sub:67"]
#[inline]
pub fn cal_sub(opt: &mem_opt_t, r: &mem_alnreg_v) -> i32 {
    // choose unique alignment
    for j in 1..r.n {
        let b_max = r.a[j].qb.max(r.a[0].qb);
        let e_min = r.a[j].qe.min(r.a[0].qe);
        // have overlap on the query
        if e_min > b_max {
            let min_l = (r.a[j].qe - r.a[j].qb).min(r.a[0].qe - r.a[0].qb);
            // significant overlap
            if (e_min - b_max) as f32 >= (min_l as f32) * opt.mask_level {
                return r.a[j].score;
            }
        }
    }
    opt.min_seed_len * opt.a
}

/// Estimate the paired-end insert-size distribution from a chunk of SE candidates.
///
/// For each of the four orientation classes (FF, FR, RF, RR) collect insert sizes of pairs
/// whose primary hits are confidently unique (no significant secondary overlap; see `cal_sub`),
/// then derive low/high cutoffs and mean+stddev from the 25/50/75 percentiles after trimming
/// with `OUTLIER_BOUND * IQR`. Orientations with too few supporting pairs or with counts below
/// `MIN_DIR_RATIO` of the dominant orientation are marked failed.
#[doc = "Original function: mem_pestat:81"]
pub fn mem_pestat(
    opt: &mem_opt_t,
    l_pac: i64,
    n: i32,
    regs: &[mem_alnreg_v],
    pes: &mut [mem_pestat_t; 4],
) {
    *pes = [mem_pestat_t::default(); 4];
    let mut isize: [Vec<u64>; 4] = std::array::from_fn(|_| Vec::new());

    for i in 0..((n as usize) >> 1) {
        let r0 = &regs[i << 1];
        let r1 = &regs[i << 1 | 1];
        if r0.n == 0 || r1.n == 0 {
            continue;
        }
        if (cal_sub(opt, r0) as f64) > MIN_RATIO * (r0.a[0].score as f64) {
            continue;
        }
        if (cal_sub(opt, r1) as f64) > MIN_RATIO * (r1.a[0].score as f64) {
            continue;
        }
        // not on the same chr
        if r0.a[0].rid != r1.a[0].rid {
            continue;
        }
        let mut is = 0_i64;
        let dir = mem_infer_dir(l_pac, r0.a[0].rb, r1.a[0].rb, &mut is) as usize;
        if is != 0 && is <= i64::from(opt.max_ins) {
            isize[dir].push(is as u64);
        }
    }

    // TODO: this block is nearly identical to the one in bwtsw2_pair.c. It would be better to
    // merge these two (upstream note).
    for d in 0..4 {
        let q = &mut isize[d];
        let r = &mut pes[d];
        // skip orientation as there are not enough pairs
        if q.len() < MIN_DIR_CNT {
            r.failed = 1;
            continue;
        }
        q.sort_unstable();
        let p25 = q[(0.25 * q.len() as f64 + 0.499) as usize] as f64;
        let p50 = q[(0.50 * q.len() as f64 + 0.499) as usize] as f64;
        let p75 = q[(0.75 * q.len() as f64 + 0.499) as usize] as f64;
        let mut low = (p25 - OUTLIER_BOUND * (p75 - p25) + 0.499) as i32;
        if low < 1 {
            low = 1;
        }
        let high = (p75 + OUTLIER_BOUND * (p75 - p25) + 0.499) as i32;
        r.low = low;
        r.high = high;

        let inside: Vec<f64> = q
            .iter()
            .map(|&x| x as f64)
            .filter(|&x| x >= r.low as f64 && x <= r.high as f64)
            .collect();
        assert!(!inside.is_empty(), "mem_pestat empty filtered insert set");
        let _ = p50;
        r.avg = inside.iter().sum::<f64>() / inside.len() as f64;
        r.std = (inside
            .iter()
            .map(|x| (x - r.avg) * (x - r.avg))
            .sum::<f64>()
            / inside.len() as f64)
            .sqrt();

        r.low = (p25 - MAPPING_BOUND * (p75 - p25) + 0.499) as i32;
        r.high = (p75 + MAPPING_BOUND * (p75 - p25) + 0.499) as i32;
        let low_bound = (r.avg - MAX_STDDEV * r.std + 0.499) as i32;
        let high_bound = (r.avg + MAX_STDDEV * r.std + 0.499) as i32;
        if r.low > low_bound {
            r.low = low_bound;
        }
        if r.high < high_bound {
            r.high = high_bound;
        }
        if r.low < 1 {
            r.low = 1;
        }
    }

    let max_pairs = isize.iter().map(Vec::len).max().unwrap_or(0);
    for d in 0..4 {
        if pes[d].failed == 0 && (isize[d].len() as f64) < (max_pairs as f64) * MIN_DIR_RATIO {
            pes[d].failed = 1;
        }
    }
}

/// Mate-rescue Smith-Waterman.
///
/// Given one mate's best alignment `a`, search for the other mate by extracting the genomic
/// window around `a` for each orientation that has not already produced a consistent pair in
/// `ma`, then running `ksw_align2`. Any high-scoring hits are inserted into `ma` (sorted by
/// descending score) and the deduplication pass is run.
///
/// # Returns
/// The number of orientations for which SW was actually executed.
#[doc = "Original function: mem_matesw:150"]
pub fn mem_matesw(
    opt: &mem_opt_t,
    bns: &bntseq_t,
    pac: &[u8],
    pes: &[mem_pestat_t; 4],
    a: &mem_alnreg_t,
    l_ms: i32,
    ms: &[u8],
    ma: &mut mem_alnreg_v,
) -> i32 {
    let l_pac = bns.l_pac;
    let mut skip = [0_i32; 4];
    for r in 0..4 {
        skip[r] = if pes[r].failed != 0 { 1 } else { 0 };
    }

    // check which orientation has been found (already produces a consistent pair distance)
    for i in 0..ma.n {
        let mut dist = 0_i64;
        let r = mem_infer_dir(l_pac, a.rb, ma.a[i].rb, &mut dist) as usize;
        if dist >= i64::from(pes[r].low) && dist <= i64::from(pes[r].high) {
            skip[r] = 1;
        }
    }
    // consistent pair exist; no need to perform SW
    if skip.iter().copied().sum::<i32>() == 4 {
        return 0;
    }

    let mut n = 0_i32;
    let l_ms_usize = l_ms as usize;
    // Reuse thread-local buffers across mem_matesw calls — was per-call alloc × ~25K-100K
    // per PE batch.
    let (mut rev_buf, mut query_buf, mut ref_buf) =
        MEM_MATESW_SCRATCH.with(|c| std::mem::take(&mut *c.borrow_mut()));
    for r in 0..4 {
        if skip[r] != 0 {
            continue;
        }
        // is_rev: whether to reverse complement the mate
        // is_larger: whether the mate has the larger coordinate
        let is_rev = ((r >> 1) != (r & 1)) as i32;
        let is_larger = ((r >> 1) == 0) as i32;
        let seq: &[u8] = if is_rev != 0 {
            // Reverse-complement: walk ms backward, complement each base via LUT.
            // LUT maps 0,1,2,3,4,5+ -> 3,2,1,0,4,4 (any value >= 4 maps to AMBIG=4).
            const RC_LUT: [u8; 8] = [3, 2, 1, 0, 4, 4, 4, 4];
            rev_buf.clear();
            rev_buf.extend(
                ms[..l_ms_usize]
                    .iter()
                    .rev()
                    .map(|&c| unsafe { *RC_LUT.get_unchecked((c as usize).min(7)) }),
            );
            &rev_buf
        } else {
            &ms[..l_ms_usize]
        };

        let mut rb;
        let mut re;
        if is_rev == 0 {
            // if on the same strand, end position should be larger to make room for the seq length
            rb = if is_larger != 0 {
                a.rb + i64::from(pes[r].low)
            } else {
                a.rb - i64::from(pes[r].high)
            };
            re = if is_larger != 0 {
                a.rb + i64::from(pes[r].high)
            } else {
                a.rb - i64::from(pes[r].low)
            };
            re += i64::from(l_ms);
        } else {
            // similarly on opposite strands: the start position is biased by l_ms instead
            rb = if is_larger != 0 {
                a.rb + i64::from(pes[r].low)
            } else {
                a.rb - i64::from(pes[r].high)
            };
            rb -= i64::from(l_ms);
            re = if is_larger != 0 {
                a.rb + i64::from(pes[r].high)
            } else {
                a.rb - i64::from(pes[r].low)
            };
        }
        if rb < 0 {
            rb = 0;
        }
        if re > (l_pac << 1) {
            re = l_pac << 1;
        }
        let mut rid = -1_i32;
        let ref_ok = if rb < re {
            let mid = (rb + re) >> 1;
            crate::bwa_mem2::bntseq::bns_fetch_seq_into(
                bns,
                pac,
                &mut rb,
                mid,
                &mut re,
                &mut rid,
                &mut ref_buf,
            );
            true
        } else {
            ref_buf.clear();
            false
        };

        // no funny things happening: rid matched and the window is at least seed-long
        if ref_ok && a.rid == rid && re - rb >= i64::from(opt.min_seed_len) {
            let xtra = KSW_XSUBO
                | KSW_XSTART
                | if l_ms * opt.a < 250 { KSW_XBYTE } else { 0 }
                | (opt.min_seed_len * opt.a);
            query_buf.clear();
            query_buf.extend_from_slice(seq);
            let span = (re - rb) as i32;
            let aln = ksw_align2(
                l_ms,
                &mut query_buf,
                span,
                &mut ref_buf,
                5,
                &opt.mat,
                opt.o_del,
                opt.e_del,
                opt.o_ins,
                opt.e_ins,
                xtra,
                None,
            );
            // something goes wrong if aln.qb < 0
            if aln.score >= opt.min_seed_len && aln.qb >= 0 {
                let mut b = mem_alnreg_t::default();
                b.rid = a.rid;
                b.is_alt = a.is_alt;
                b.qb = if is_rev != 0 {
                    l_ms - (aln.qe + 1)
                } else {
                    aln.qb
                };
                b.qe = if is_rev != 0 {
                    l_ms - aln.qb
                } else {
                    aln.qe + 1
                };
                b.rb = if is_rev != 0 {
                    (l_pac << 1) - (rb + i64::from(aln.te) + 1)
                } else {
                    rb + i64::from(aln.tb)
                };
                b.re = if is_rev != 0 {
                    (l_pac << 1) - (rb + i64::from(aln.tb))
                } else {
                    rb + i64::from(aln.te) + 1
                };
                b.score = aln.score;
                b.csub = aln.score2;
                b.secondary = -1;
                b.seedcov = ((b.re - b.rb).min(i64::from(b.qe - b.qb)) >> 1) as i32;

                // make room for a new element, then move b s.t. ma stays sorted by score
                ma.a.push(b);
                ma.n = ma.a.len();
                ma.m = ma.a.len();

                // find the insertion point
                let mut insert_at = ma.n - 1;
                for i in 0..ma.n - 1 {
                    if ma.a[i].score < b.score {
                        insert_at = i;
                        break;
                    }
                }
                if insert_at < ma.n - 1 {
                    ma.a[insert_at..].rotate_right(1);
                    ma.a[insert_at] = b;
                }
            }
            n += 1;
        }
        if n != 0 {
            ma.n = mem_sort_dedup_patch(opt, None, None, None, ma.n as i32, &mut ma.a) as usize;
            ma.m = ma.a.len();
        }
    }
    MEM_MATESW_SCRATCH.with(|c| *c.borrow_mut() = (rev_buf, query_buf, ref_buf));
    n
}

thread_local! {
    // Reused across mem_pair calls (~25K per 50K PE batch). v holds per-alignment pairs sorted
    // by ref position; u holds candidate pair scores. Per-call allocations would otherwise show
    // up in PE-rescue churn.
    static MEM_PAIR_SCRATCH: std::cell::RefCell<(Vec<pair64_t>, Vec<pair64_t>)> =
        const { std::cell::RefCell::new((Vec::new(), Vec::new())) };

    // Reused across mem_sam_pe_batch calls (per batch). Holds the rev-pass pair list built
    // from phase-0 results.
    static MEM_PE_BATCH_PHASE1_PAIRS: std::cell::RefCell<Vec<SeqPair>> =
        const { std::cell::RefCell::new(Vec::new()) };

    // Reused across mem_matesw* calls (per pair × 4 directions). Hold the rev-complement
    // query, copied query for ksw_align2, and reference sequence buffer.
    static MEM_MATESW_SCRATCH: std::cell::RefCell<(Vec<u8>, Vec<u8>, Vec<u8>)> =
        const { std::cell::RefCell::new((Vec::new(), Vec::new(), Vec::new())) };

    // Reused across mem_sam_pe / mem_sam_pe_batch_post calls (per pair). Each call previously
    // allocated `b: [Vec<mem_alnreg_t>; 2]` via collect; ~700K Vec allocations per 700K-read run.
    static MEM_SAM_PE_B_SCRATCH: std::cell::RefCell<(Vec<mem_alnreg_t>, Vec<mem_alnreg_t>)> =
        const { std::cell::RefCell::new((Vec::new(), Vec::new())) };

    // Reused across mem_sam_pe / mem_sam_pe_batch_post calls (per pair). Holds aa (the per-end
    // mem_aln_t lists). Each call previously allocated 2 fresh Vec<mem_aln_t>; ~700K * 2 Vec
    // allocations per 700K-read run.
    static MEM_SAM_PE_AA_SCRATCH: std::cell::RefCell<(Vec<mem_aln_t>, Vec<mem_aln_t>)> =
        const { std::cell::RefCell::new((Vec::new(), Vec::new())) };
}

// Drain h/g/aa mem_aln_t cigar/md buffers into the global pools so they're reused on next call
// instead of being freed when these stack-locals go out of scope. Also returns the aa Vec
// capacity to the per-pair pool.
#[inline]
fn drain_pe_aln_pools(
    h: &mut [mem_aln_t; 2],
    g: &mut [mem_aln_t; 2],
    aa: &mut [Vec<mem_aln_t>; 2],
) {
    for i in 0..2 {
        crate::bwa_mem2::bwa::return_cigar_buf(std::mem::take(&mut h[i].cigar));
        crate::bwa_mem2::bwa::return_md_buf(std::mem::take(&mut h[i].md));
        crate::bwa_mem2::bwa::return_cigar_buf(std::mem::take(&mut g[i].cigar));
        crate::bwa_mem2::bwa::return_md_buf(std::mem::take(&mut g[i].md));
        for aln in aa[i].drain(..) {
            crate::bwa_mem2::bwa::return_cigar_buf(aln.cigar);
            crate::bwa_mem2::bwa::return_md_buf(aln.md);
        }
    }
    // Return aa Vec capacity to the thread-local pool. drain(..) above already cleared them.
    let aa0 = std::mem::take(&mut aa[0]);
    let aa1 = std::mem::take(&mut aa[1]);
    MEM_SAM_PE_AA_SCRATCH.with(|c| *c.borrow_mut() = (aa0, aa1));
}

#[inline]
fn take_pe_aa_pools() -> [Vec<mem_aln_t>; 2] {
    let (mut a0, mut a1) = MEM_SAM_PE_AA_SCRATCH.with(|c| std::mem::take(&mut *c.borrow_mut()));
    a0.clear();
    a1.clear();
    [a0, a1]
}

/// Pair two reads' single-end candidate alignment sets into the best proper pair.
///
/// Builds a flat list keyed by `(rid, forward-strand position)` over both ends' primary hits,
/// sorted so a single sweep finds candidate pairs whose distance falls in
/// `pes[dir].low..high`. For each candidate, the pair score is the sum of SE scores plus a
/// Gaussian-fit insert-size term
/// `.721 * log(2 * erfc(|ns| / sqrt(2))) * opt.a`, where `ns` is the insert size's z-score,
/// quantized with the C++ `(int)(x + 0.499)` bias-and-truncate rule and clamped to `>= 0`.
///
/// # Returns
/// The best pair score. Writes `z[0]`/`z[1]` with the chosen indices and reports `subo`/`n_sub`
/// for MAPQ calibration.
#[doc = "Original function: mem_pair:285"]
pub fn mem_pair(
    opt: &mem_opt_t,
    bns: &bntseq_t,
    _pac: &[u8],
    pes: &[mem_pestat_t; 4],
    a: &mut [mem_alnreg_v; 2],
    id: i32,
    sub: &mut i32,
    n_sub: &mut i32,
    z: &mut [i32; 2],
    n_pri: &[i32; 2],
) -> i32 {
    MEM_PAIR_SCRATCH.with(|cell| {
        let mut buf = cell.borrow_mut();
        let (v, u) = &mut *buf;
        v.clear();
        u.clear();
        mem_pair_inner(opt, bns, pes, a, id, sub, n_sub, z, n_pri, v, u)
    })
}

#[allow(clippy::too_many_arguments)]
fn mem_pair_inner(
    opt: &mem_opt_t,
    bns: &bntseq_t,
    pes: &[mem_pestat_t; 4],
    a: &mut [mem_alnreg_v; 2],
    id: i32,
    sub: &mut i32,
    n_sub: &mut i32,
    z: &mut [i32; 2],
    n_pri: &[i32; 2],
    v: &mut Vec<pair64_t>,
    u: &mut Vec<pair64_t>,
) -> i32 {
    let l_pac = bns.l_pac;

    // loop through read number: pack each candidate's (rid, fwd-strand pos) into x and
    // (score, idx, strand-of-rb, read-side r) into y
    for r in 0..2_usize {
        for i in 0..(n_pri[r] as usize) {
            let e = &a[r].a[i];
            // forward position
            let mut x = if e.rb < l_pac {
                e.rb
            } else {
                (l_pac << 1) - 1 - e.rb
            };
            x = (i64::from(e.rid) << 32) | (x - bns.anns[e.rid as usize].offset);
            let y = ((e.score as u64) << 32)
                | ((i as u64) << 2)
                | (u64::from((e.rb >= l_pac) as u8) << 1)
                | (r as u64);
            v.push(pair64_t { x: x as u64, y });
        }
    }

    v.sort_by(|lhs, rhs| lhs.x.cmp(&rhs.x).then(lhs.y.cmp(&rhs.y)));
    // y_last keeps the last hit for each of 4 (read-side, strand-of-rb) classes
    let mut y_last = [-1_i32; 4];
    for i in 0..v.len() {
        // loop through direction: orientation = (read-side r) << 1 | strand-of-rb-of-v[i]
        for r in 0..2_u64 {
            let dir = ((r << 1) | ((v[i].y >> 1) & 1)) as usize;
            // invalid orientation
            if pes[dir].failed != 0 {
                continue;
            }
            let which = ((r << 1) | ((v[i].y & 1) ^ 1)) as usize;
            // no previous hits in the complementary class
            if y_last[which] < 0 {
                continue;
            }
            // TODO: this is a O(n^2) solution in the worst case; remember to check if this loop
            // takes a lot of time (upstream: "I doubt").
            for k in (0..=(y_last[which] as usize)).rev() {
                if (v[k].y & 3) != (which as u64) {
                    continue;
                }
                let dist = (v[i].x as i64) - (v[k].x as i64);
                if dist > i64::from(pes[dir].high) {
                    break;
                }
                if dist < i64::from(pes[dir].low) {
                    continue;
                }
                // ns is the dist's z-score against the orientation's fitted Normal
                let ns = (dist as f64 - pes[dir].avg) / pes[dir].std;
                let mut q = ((v[i].y >> 32) + (v[k].y >> 32)) as f64;
                // .721 = 1/log(4) — converts ln(2*erfc(...)) into log4 units so it scales like a
                // mapping-quality bit count.
                let erfc_term = erfc_approx(ns.abs() * std::f64::consts::FRAC_1_SQRT_2);
                q += 0.721 * (2.0_f64 * erfc_term).ln() * (opt.a as f64);
                // Match C++: (int)(... + 0.499) truncates toward zero, then clamp >= 0.
                let q_int = (q + 0.499) as i32;
                let q = q_int.max(0) as u64;
                let pair_y = ((k as u64) << 32) | (i as u64);
                let pair_x = (q << 32) | (hash_64(pair_y ^ ((id as u64) << 8)) & 0xffff_ffff);
                u.push(pair64_t {
                    x: pair_x,
                    y: pair_y,
                });
            }
        }
        y_last[(v[i].y & 3) as usize] = i as i32;
    }

    // found at least one proper pair
    if !u.is_empty() {
        let mut tmp = opt.a + opt.b;
        tmp = tmp.max(opt.o_del + opt.e_del);
        tmp = tmp.max(opt.o_ins + opt.e_ins);
        u.sort_by(|lhs, rhs| lhs.x.cmp(&rhs.x).then(lhs.y.cmp(&rhs.y)));
        let best = *u.last().expect("best pair");
        let i = (best.y >> 32) as usize;
        let k = (best.y as u32) as usize;
        // index of the best pair (per read side)
        z[(v[i].y & 1) as usize] = ((v[i].y >> 2) & 0x3fff_ffff) as i32;
        z[(v[k].y & 1) as usize] = ((v[k].y >> 2) & 0x3fff_ffff) as i32;
        let ret = (best.x >> 32) as i32;
        *sub = if u.len() > 1 {
            (u[u.len() - 2].x >> 32) as i32
        } else {
            0
        };
        *n_sub = 0;
        for pair in u.iter().rev().skip(1) {
            if *sub - (pair.x >> 32) as i32 <= tmp {
                *n_sub += 1;
            }
        }
        ret
    } else {
        *sub = 0;
        *n_sub = 0;
        0
    }
}

/// Paired-end SAM-output entry point.
///
/// Runs mate-rescue SW (`mem_matesw`) for each candidate above the unpaired-penalty floor on
/// each end (unless `MEM_F_NO_RESCUE`), then primary marking, then pairing via `mem_pair`. When
/// a proper pair survives the multi-pair gate, MAPQs are computed (`raw_mapq` from `o - subo`,
/// capped at 60, capped at the tandem-repeat score, plus a `frac_rep` haircut), and per-record
/// SAM lines are emitted via `mem_aln2sam`. The "no_pairing" fallback emits two independent
/// SE-style records and may still flag them as a proper pair if their best hits land within
/// `pes[d].low..high`.
#[doc = "Original function: mem_sam_pe:353"]
pub fn mem_sam_pe(
    opt: &mem_opt_t,
    bns: &bntseq_t,
    pac: &[u8],
    pes: &[mem_pestat_t; 4],
    id: u64,
    s: &mut [bseq1_t; 2],
    a: &mut [mem_alnreg_v; 2],
) -> i32 {
    let mut n = 0_i32;
    let mut z = [0_i32; 2];
    let mut subo = 0_i32;
    let mut n_sub = 0_i32;
    let mut extra_flag = 1_i32;
    let mut n_pri = [0_i32; 2];
    // Pre-size for two SAM lines (~1KB) — typical paired-end SAM record pair fits without growth.
    let mut str_ = kstring_t::default();
    crate::bwa_mem2::kstring::ks_resize(&mut str_, 1024);
    let mut h: [mem_aln_t; 2] = std::array::from_fn(|_| mem_aln_t::default());
    let mut g: [mem_aln_t; 2] = std::array::from_fn(|_| mem_aln_t::default());
    let mut aa: [Vec<mem_aln_t>; 2] = take_pe_aa_pools();

    // flag 0x20 == MEM_F_NO_RESCUE; without it, perform SW for the best alignment(s) on each end
    if (opt.flag & 0x20) == 0 {
        let (mut b0, mut b1) = MEM_SAM_PE_B_SCRATCH.with(|c| std::mem::take(&mut *c.borrow_mut()));
        // b[i] = candidates within pen_unpaired of the best — these are the seeds for mate-SW
        b0.clear();
        b0.extend(
            a[0].a
                .iter()
                .take(a[0].n)
                .filter(|reg| reg.score >= a[0].a[0].score - opt.pen_unpaired)
                .copied(),
        );
        b1.clear();
        b1.extend(
            a[1].a
                .iter()
                .take(a[1].n)
                .filter(|reg| reg.score >= a[1].a[0].score - opt.pen_unpaired)
                .copied(),
        );
        let b = [&b0, &b1];
        let max_matesw = max_matesw_rescue_limit(opt);
        for i in 0..2_usize {
            // Hoist nt4_cow outside the j loop — same mate seq each iter.
            let mate_nt4 = nt4_cow(&s[1 - i]);
            for j in 0..b[i].len().min(max_matesw) {
                let val = mem_matesw(
                    opt,
                    bns,
                    pac,
                    pes,
                    &b[i][j],
                    s[1 - i].l_seq,
                    &mate_nt4,
                    &mut a[1 - i],
                );
                n += val;
            }
        }
        MEM_SAM_PE_B_SCRATCH.with(|c| *c.borrow_mut() = (b0, b1));
    }

    n_pri[0] = mem_mark_primary_se(opt, a[0].n as i32, &mut a[0].a, (id << 1) as i64);
    n_pri[1] = mem_mark_primary_se(opt, a[1].n as i32, &mut a[1].a, ((id << 1) | 1) as i64);
    if (opt.flag & MEM_F_PRIMARY5) != 0 {
        mem_reorder_primary5(opt.T, &mut a[0]);
        mem_reorder_primary5(opt.T, &mut a[1]);
    }
    if debug_trace_read(s[0].name.as_deref()) {
        eprintln!(
            "[trace::pe_batch_post:after_primary] read={} id={} n_pri={:?} side0={:?} side1={:?}",
            s[0].name.as_deref().unwrap_or(""),
            id,
            n_pri,
            a[0].a
                .iter()
                .take(a[0].n.min(12))
                .map(|r| (
                    r.score,
                    r.sub,
                    r.csub,
                    r.qb,
                    r.qe,
                    r.rb,
                    r.re,
                    r.secondary,
                    r.secondary_all
                ))
                .collect::<Vec<_>>(),
            a[1].a
                .iter()
                .take(a[1].n.min(12))
                .map(|r| (
                    r.score,
                    r.sub,
                    r.csub,
                    r.qb,
                    r.qe,
                    r.rb,
                    r.re,
                    r.secondary,
                    r.secondary_all
                ))
                .collect::<Vec<_>>(),
        );
    }

    if (opt.flag & MEM_F_NOPAIRING) == 0 {
        let o = if n_pri[0] > 0 && n_pri[1] > 0 {
            mem_pair(
                opt, bns, pac, pes, a, id as i32, &mut subo, &mut n_sub, &mut z, &n_pri,
            )
        } else {
            0
        };
        if debug_trace_read(s[0].name.as_deref()) {
            eprintln!(
                "[trace::pe_batch_post:after_pair] read={} id={} n_pri={:?} o={} subo={} n_sub={} z={:?} a0={:?} a1={:?}",
                s[0].name.as_deref().unwrap_or(""),
                id,
                n_pri,
                o,
                subo,
                n_sub,
                z,
                a[0]
                    .a
                    .iter()
                    .take(a[0].n.min(16))
                    .map(|r| (r.score, r.rb, r.qb, r.qe, r.sub, r.csub, r.secondary, r.secondary_all))
                    .collect::<Vec<_>>(),
                a[1]
                    .a
                    .iter()
                    .take(a[1].n.min(16))
                    .map(|r| (r.score, r.rb, r.qb, r.qe, r.sub, r.csub, r.secondary, r.secondary_all))
                    .collect::<Vec<_>>(),
            );
        }
        if o > 0 {
            // check if an end has multiple hits even after mate-SW
            let mut is_multi = [0_i32; 2];
            for i in 0..2_usize {
                let mut j = 1_usize;
                while j < n_pri[i] as usize {
                    if a[i].a[j].secondary < 0 && a[i].a[j].score >= opt.T {
                        break;
                    }
                    j += 1;
                }
                is_multi[i] = if j < n_pri[i] as usize { 1 } else { 0 };
            }
            // TODO (upstream): in rare cases, the true hit may be long but with low score
            if is_multi == [0, 0] {
                // compute mapQ for the best SE hit
                let score_un = a[0].a[0].score + a[1].a[0].score - opt.pen_unpaired;
                subo = subo.max(score_un);
                let mut q_pe = raw_mapq(o - subo, opt.a);
                if n_sub > 0 {
                    q_pe -= (4.343 * ((n_sub + 1) as f64).ln() + 0.499) as i32;
                }
                q_pe = q_pe.clamp(0, 60);
                q_pe = (q_pe as f64
                    * (1.0 - 0.5 * (a[0].a[0].frac_rep as f64 + a[1].a[0].frac_rep as f64))
                    + 0.499) as i32;

                // the following assumes no split hits
                let mut q_se = [0_i32; 2];
                if o > score_un {
                    // paired alignment is preferred
                    for i in 0..2_usize {
                        let zi = z[i] as usize;
                        let secondary = a[i].a[zi].secondary;
                        let secondary_score = if secondary >= 0 {
                            Some(a[i].a[secondary as usize].score)
                        } else {
                            None
                        };
                        let c = &mut a[i].a[zi];
                        if let Some(score) = secondary_score {
                            c.sub = score;
                            c.secondary = -2;
                        }
                        q_se[i] = mem_approx_mapq_se(opt, c);
                    }
                    q_se[0] = if q_se[0] > q_pe {
                        q_se[0]
                    } else if q_pe < q_se[0] + 40 {
                        q_pe
                    } else {
                        q_se[0] + 40
                    };
                    q_se[1] = if q_se[1] > q_pe {
                        q_se[1]
                    } else if q_pe < q_se[1] + 40 {
                        q_pe
                    } else {
                        q_se[1] + 40
                    };
                    extra_flag |= 2;
                    // cap at the tandem repeat score
                    for i in 0..2_usize {
                        let c = &a[i].a[z[i] as usize];
                        q_se[i] = q_se[i].min(raw_mapq(c.score - c.csub, opt.a));
                    }
                } else {
                    // the unpaired alignment is preferred
                    z = [0, 0];
                    q_se[0] = mem_approx_mapq_se(opt, &a[0].a[0]);
                    q_se[1] = mem_approx_mapq_se(opt, &a[1].a[0]);
                }
                if debug_trace_read(s[0].name.as_deref()) {
                    eprintln!(
                        "[trace::pe_batch_post:chosen] read={} id={} z={:?} q_se={:?} extra_flag={} chosen0={:?} chosen1={:?}",
                        s[0].name.as_deref().unwrap_or(""),
                        id,
                        z,
                        q_se,
                        extra_flag,
                        a[0]
                            .a
                            .get(usize::try_from(z[0].max(0)).expect("z0"))
                            .map(|r| (r.score, r.sub, r.csub, r.secondary, r.secondary_all)),
                        a[1]
                            .a
                            .get(usize::try_from(z[1].max(0)).expect("z1"))
                            .map(|r| (r.score, r.sub, r.csub, r.secondary, r.secondary_all)),
                    );
                }

                // switch secondary and primary if both of them are non-ALT
                for i in 0..2_usize {
                    let zi = z[i] as usize;
                    let k = a[i].a[zi].secondary_all;
                    if k >= 0 && k < n_pri[i] {
                        for j in 0..a[i].n {
                            if a[i].a[j].secondary_all == k || j == k as usize {
                                a[i].a[j].secondary_all = z[i];
                            }
                        }
                        a[i].a[zi].secondary_all = -1;
                    }
                }

                let mut xa = if (opt.flag & MEM_F_ALL) == 0 {
                    [
                        mem_gen_alt(
                            opt,
                            bns,
                            pac,
                            &a[0],
                            s[0].l_seq,
                            query_string_for_aln(&s[0]).as_ref(),
                        ),
                        mem_gen_alt(
                            opt,
                            bns,
                            pac,
                            &a[1],
                            s[1].l_seq,
                            query_string_for_aln(&s[1]).as_ref(),
                        ),
                    ]
                } else {
                    [vec![None; a[0].n], vec![None; a[1].n]]
                };

                for i in 0..2_usize {
                    let zi = z[i] as usize;
                    h[i] = mem_reg2aln(
                        opt,
                        bns,
                        pac,
                        s[i].l_seq,
                        query_string_for_aln(&s[i]).as_ref(),
                        Some(&a[i].a[zi]),
                    );
                    h[i].mapq = q_se[i].max(0) as u32;
                    h[i].flag |= (0x40 << i) | extra_flag;
                    // Take XA from xa rather than clone — xa goes out of scope at the end of
                    // this block so the take is safe.
                    h[i].XA = std::mem::take(&mut xa[i][zi]);
                    aa[i].push(h[i].clone());
                    // the read has ALT hits
                    if (n_pri[i] as usize) < a[i].n {
                        let p = &a[i].a[n_pri[i] as usize];
                        if p.score >= opt.T && p.secondary < 0 && p.is_alt != 0 {
                            g[i] = mem_reg2aln(
                                opt,
                                bns,
                                pac,
                                s[i].l_seq,
                                query_string_for_aln(&s[i]).as_ref(),
                                Some(p),
                            );
                            g[i].flag |= 0x800 | (0x40 << i) | extra_flag;
                            g[i].XA = std::mem::take(&mut xa[i][n_pri[i] as usize]);
                            aa[i].push(g[i].clone());
                        }
                    }
                }

                // write SAM: read1 hits, then read2 hits — each record carries the mate's primary
                // aln in `m` so flag/insert-size/MRNM fields can be filled
                for i in 0..aa[0].len() {
                    mem_aln2sam(
                        opt,
                        bns,
                        &mut str_,
                        &s[0],
                        aa[0].len() as i32,
                        &aa[0],
                        i as i32,
                        Some(&h[1]),
                    );
                }
                s[0].sam = Some(take_kstring_boxed(&mut str_));
                for i in 0..aa[1].len() {
                    mem_aln2sam(
                        opt,
                        bns,
                        &mut str_,
                        &s[1],
                        aa[1].len() as i32,
                        &aa[1],
                        i as i32,
                        Some(&h[0]),
                    );
                }
                s[1].sam = Some(take_kstring_boxed(&mut str_));
                assert_eq!(s[0].name, s[1].name, "paired reads have different names");
                drain_pe_aln_pools(&mut h, &mut g, &mut aa);
                return n;
            }
        }
    }

    for i in 0..2_usize {
        let mut which = -1_i32;
        if a[i].n > 0 {
            if a[i].a[0].score >= opt.T {
                which = 0;
            } else if (n_pri[i] as usize) < a[i].n && a[i].a[n_pri[i] as usize].score >= opt.T {
                which = n_pri[i];
            }
        }
        h[i] = if which >= 0 {
            mem_reg2aln(
                opt,
                bns,
                pac,
                s[i].l_seq,
                query_string_for_aln(&s[i]).as_ref(),
                Some(&a[i].a[which as usize]),
            )
        } else {
            mem_reg2aln(
                opt,
                bns,
                pac,
                s[i].l_seq,
                query_string_for_aln(&s[i]).as_ref(),
                None,
            )
        };
    }
    // if the top hits from the two ends constitute a proper pair, flag it.
    if (opt.flag & MEM_F_NOPAIRING) == 0 && h[0].rid == h[1].rid && h[0].rid >= 0 {
        let mut dist = 0_i64;
        let d = mem_infer_dir(bns.l_pac, a[0].a[0].rb, a[1].a[0].rb, &mut dist) as usize;
        if pes[d].failed == 0 && dist >= i64::from(pes[d].low) && dist <= i64::from(pes[d].high) {
            extra_flag |= 2;
        }
    }
    mem_reg2sam(
        opt,
        bns,
        pac,
        &mut s[0],
        &mut a[0],
        0x41 | extra_flag,
        Some(&h[1]),
    );
    mem_reg2sam(
        opt,
        bns,
        pac,
        &mut s[1],
        &mut a[1],
        0x81 | extra_flag,
        Some(&h[0]),
    );
    assert_eq!(s[0].name, s[1].name, "paired reads have different names");
    drain_pe_aln_pools(&mut h, &mut g, &mut aa);
    n
}

/// Pre-pass for the batched mate-rescue path.
///
/// Enqueues every `(pair, direction)` SW job into the thread-local `SeqPair` arrays in `mmc` so
/// a later single SIMD pass (`mem_sam_pe_batch`) processes them all at once. NEW, batching —
/// replaces per-pair `ksw_align2` calls.
#[doc = "Original function: mem_sam_pe_batch_pre:553"]
pub fn mem_sam_pe_batch_pre(
    opt: &mem_opt_t,
    bns: &bntseq_t,
    pac: &[u8],
    pes: &[mem_pestat_t; 4],
    _id: u64,
    s: &[bseq1_t; 2],
    a: &[mem_alnreg_v; 2],
    mmc: &mut mem_cache,
    pcnt: &mut i32,
    gcnt: &mut i32,
    maxRefLen: &mut i32,
    maxQerLen: &mut i32,
    tid: usize,
) -> i32 {
    if (opt.flag & 0x20) != 0 {
        return 1;
    }
    for i in 0..2_usize {
        if a[i].n == 0 {
            continue;
        }
        let min_score = a[i].a[0].score - opt.pen_unpaired;
        let mut rescue_count = 0_usize;
        let max_matesw = max_matesw_rescue_limit(opt);
        // Hoist nt4_cow outside the per-reg loop — same mate seq each iter.
        let mate_nt4 = nt4_cow(&s[1 - i]);
        for reg in a[i].a.iter().take(a[i].n) {
            if reg.score < min_score {
                continue;
            }
            if rescue_count >= max_matesw {
                break;
            }
            *pcnt = mem_matesw_batch_pre(
                opt,
                bns,
                pac,
                pes,
                reg,
                s[1 - i].l_seq,
                &mate_nt4,
                &a[1 - i],
                mmc,
                *pcnt,
                *gcnt,
                maxRefLen,
                maxQerLen,
                tid,
            );
            *gcnt += 4;
            rescue_count += 1;
        }
    }
    1
}

/// In-place reverse of the first `l` bytes of `s`.
///
/// Used to flip query / reference between Smith-Waterman phases in the batched mate-rescue
/// path.
#[doc = "Original function: revseq:604"]
pub fn revseq(l: i32, s: &mut [u8]) {
    let l = l.max(0) as usize;
    let mut i = 0_usize;
    while i < (l >> 1) {
        s.swap(i, l - 1 - i);
        i += 1;
    }
}

/// Run the batched mate-rescue Smith-Waterman.
///
/// Equivalent to `ksw_align2` but vectorized over many pairs at once. Two-phase: phase-0
/// forward extension via `get_scores8/16` by size class, then reverse the hit windows in place
/// and run phase-1 to recover the start coordinates. Output goes into `aln` indexed by
/// `SeqPair.regid`.
#[doc = "Original function: mem_sam_pe_batch:612"]
pub fn mem_sam_pe_batch(
    opt: &mem_opt_t,
    mmc: &mut mem_cache,
    pcnt: i32,
    pcnt8: i32,
    aln: &mut [kswr_t],
    maxRefLen: i32,
    maxQerLen: i32,
    tid: usize,
) -> i32 {
    let total = pcnt.max(0) as usize;
    if total == 0 {
        return 1;
    }
    let pcnt8_usize = pcnt8.max(0) as usize;
    let pcnt16_usize = total.saturating_sub(pcnt8_usize);

    for r in aln.iter_mut().take(total) {
        r.tb = -1;
        r.qb = -1;
    }

    let pwsw = kswv::ctor(
        opt.o_del,
        opt.e_del,
        opt.o_ins,
        opt.e_ins,
        i8::try_from(opt.a).unwrap_or(i8::MAX),
        i8::try_from(-opt.b).unwrap_or(i8::MIN),
        1,
        Some(maxRefLen),
        Some(maxQerLen),
    );

    // Phase 0: forward pass on each size class
    if pcnt8_usize > 0 {
        let pairs = &mut mmc.seqPairArrayLeft128[tid][..pcnt8_usize];
        pwsw.get_scores8(
            pairs,
            &mmc.seqBufLeftRef[tid],
            &mmc.seqBufLeftQer[tid],
            aln,
            pcnt8_usize as i32,
            1,
            0,
        );
    }
    if pcnt16_usize > 0 {
        let pairs = &mut mmc.seqPairArrayLeft128[tid][pcnt8_usize..total];
        pwsw.get_scores16(
            pairs,
            &mmc.seqBufLeftRef[tid],
            &mmc.seqBufLeftQer[tid],
            aln,
            pcnt16_usize as i32,
            1,
            0,
        );
    }

    // Phase 1 prep: reverse buffers in place, build the rev-pass pair list (u8 first, then i16)
    MEM_PE_BATCH_PHASE1_PAIRS.with(|cell| {
        let mut buf = cell.borrow_mut();
        let phase1_pairs: &mut Vec<SeqPair> = &mut buf;
        phase1_pairs.clear();
        phase1_pairs.reserve(total);
        let mut pos8 = 0_usize;
        let mut pos16 = 0_usize;
        for class in 0..2_usize {
            let (start, count) = if class == 0 {
                (0, pcnt8_usize)
            } else {
                (pcnt8_usize, pcnt16_usize)
            };
            for i in 0..count {
                let mut sp = mmc.seqPairArrayLeft128[tid][start + i];
                let ind = sp.regid as usize;
                let r = aln[ind];
                let xtra = sp.h0;
                if (xtra & KSW_XSTART) == 0
                    || ((xtra & KSW_XSUBO) != 0 && r.score < (xtra & 0xffff))
                {
                    continue;
                }
                sp.h0 = KSW_XSTOP | r.score;
                sp.len1 = r.te + 1;
                sp.len2 = r.qe + 1;
                let qs_start = sp.idq as usize;
                let rs_start = sp.idr as usize;
                revseq(r.qe + 1, &mut mmc.seqBufLeftQer[tid][qs_start..]);
                revseq(r.te + 1, &mut mmc.seqBufLeftRef[tid][rs_start..]);
                phase1_pairs.push(sp);
                if class == 0 {
                    pos8 += 1;
                } else {
                    pos16 += 1;
                }
            }
        }

        // Phase 1: i16 first then u8 (matches upstream order)
        if pos16 > 0 {
            let pairs = &mut phase1_pairs[pos8..pos8 + pos16];
            pwsw.get_scores16(
                pairs,
                &mmc.seqBufLeftRef[tid],
                &mmc.seqBufLeftQer[tid],
                aln,
                pos16 as i32,
                1,
                1,
            );
        }
        if pos8 > 0 {
            let pairs = &mut phase1_pairs[..pos8];
            pwsw.get_scores8(
                pairs,
                &mmc.seqBufLeftRef[tid],
                &mmc.seqBufLeftQer[tid],
                aln,
                pos8 as i32,
                1,
                1,
            );
        }
    });

    let _ = ksw_align2;
    1
}

/// Post-pass for the batched mate-rescue path.
///
/// Consumes the SIMD SW results in `myaln` produced by `mem_sam_pe_batch`, then mirrors
/// `mem_sam_pe`'s pair scoring / MAPQ / SAM emission. The `mem_matesw_batch_post` call replays
/// each `(pair, direction)` using the cached alignment (or falls back to `ksw_align2` when the
/// pre-pass marked the slot as `-1`).
#[doc = "Original function: mem_sam_pe_batch_post:713"]
pub fn mem_sam_pe_batch_post(
    opt: &mem_opt_t,
    bns: &bntseq_t,
    pac: &[u8],
    pes: &[mem_pestat_t; 4],
    id: u64,
    s: &mut [bseq1_t; 2],
    a: &mut [mem_alnreg_v; 2],
    myaln: &[kswr_t],
    mmc: &mut mem_cache,
    gcnt: &mut i32,
    tid: usize,
) -> i32 {
    let mut n = 0_i32;
    let mut z = [0_i32; 2];
    let mut subo = 0_i32;
    let mut n_sub = 0_i32;
    let mut extra_flag = 1_i32;
    let mut n_pri = [0_i32; 2];
    // Pre-size for two SAM lines (~1KB) — typical paired-end SAM record pair fits without growth.
    let mut str_ = kstring_t::default();
    crate::bwa_mem2::kstring::ks_resize(&mut str_, 1024);
    let mut h: [mem_aln_t; 2] = std::array::from_fn(|_| mem_aln_t::default());
    let mut g: [mem_aln_t; 2] = std::array::from_fn(|_| mem_aln_t::default());
    let mut aa: [Vec<mem_aln_t>; 2] = take_pe_aa_pools();
    // flag 0x20 == MEM_F_NO_RESCUE; without it, replay the batched SW results into mem_alnreg_v.
    if (opt.flag & 0x20) == 0 {
        let (mut b0, mut b1) = MEM_SAM_PE_B_SCRATCH.with(|c| std::mem::take(&mut *c.borrow_mut()));
        // b[i] = candidates within pen_unpaired of the best — seeds for mate-SW
        b0.clear();
        b0.extend(
            a[0].a
                .iter()
                .take(a[0].n)
                .filter(|reg| reg.score >= a[0].a[0].score - opt.pen_unpaired)
                .copied(),
        );
        b1.clear();
        b1.extend(
            a[1].a
                .iter()
                .take(a[1].n)
                .filter(|reg| reg.score >= a[1].a[0].score - opt.pen_unpaired)
                .copied(),
        );
        let b = [&b0, &b1];
        // Reset per-reg .flg before the batched post pass (used by dedup downstream).
        for reg in &mut a[0].a[..a[0].n] {
            reg.flg = 0;
        }
        for reg in &mut a[1].a[..a[1].n] {
            reg.flg = 0;
        }
        let max_matesw = max_matesw_rescue_limit(opt);
        for i in 0..2_usize {
            for j in 0..b[i].len().min(max_matesw) {
                let mate_nt4 = nt4_cow(&s[1 - i]);
                let val = mem_matesw_batch_post(
                    opt,
                    bns,
                    pac,
                    pes,
                    &b[i][j],
                    s[1 - i].l_seq,
                    &mate_nt4,
                    &mut a[1 - i],
                    myaln,
                    *gcnt,
                    &mmc.seqPairArrayAux[tid],
                );
                n += val;
                *gcnt += 4;
            }
        }
        MEM_SAM_PE_B_SCRATCH.with(|c| *c.borrow_mut() = (b0, b1));
    }

    n_pri[0] = mem_mark_primary_se(opt, a[0].n as i32, &mut a[0].a, (id << 1) as i64);
    n_pri[1] = mem_mark_primary_se(opt, a[1].n as i32, &mut a[1].a, ((id << 1) | 1) as i64);
    if (opt.flag & MEM_F_PRIMARY5) != 0 {
        mem_reorder_primary5(opt.T, &mut a[0]);
        mem_reorder_primary5(opt.T, &mut a[1]);
    }
    if debug_trace_read(s[0].name.as_deref()) {
        eprintln!(
            "[trace::pe_batch_post2:after_primary] read={} id={} n_pri={:?} side0={:?} side1={:?}",
            s[0].name.as_deref().unwrap_or(""),
            id,
            n_pri,
            a[0].a
                .iter()
                .take(a[0].n.min(12))
                .map(|r| (
                    r.score,
                    r.sub,
                    r.csub,
                    r.qb,
                    r.qe,
                    r.rb,
                    r.re,
                    r.secondary,
                    r.secondary_all
                ))
                .collect::<Vec<_>>(),
            a[1].a
                .iter()
                .take(a[1].n.min(12))
                .map(|r| (
                    r.score,
                    r.sub,
                    r.csub,
                    r.qb,
                    r.qe,
                    r.rb,
                    r.re,
                    r.secondary,
                    r.secondary_all
                ))
                .collect::<Vec<_>>(),
        );
    }

    if (opt.flag & MEM_F_NOPAIRING) == 0 {
        let o = if n_pri[0] > 0 && n_pri[1] > 0 {
            mem_pair(
                opt, bns, pac, pes, a, id as i32, &mut subo, &mut n_sub, &mut z, &n_pri,
            )
        } else {
            0
        };
        if debug_trace_read(s[0].name.as_deref()) {
            eprintln!(
                "[trace::pe_batch_post2:after_pair] read={} id={} n_pri={:?} o={} subo={} n_sub={} z={:?} a0={:?} a1={:?}",
                s[0].name.as_deref().unwrap_or(""),
                id,
                n_pri,
                o,
                subo,
                n_sub,
                z,
                a[0]
                    .a
                    .iter()
                    .take(a[0].n.min(16))
                    .map(|r| (r.score, r.rb, r.qb, r.qe, r.sub, r.csub, r.secondary, r.secondary_all))
                    .collect::<Vec<_>>(),
                a[1]
                    .a
                    .iter()
                    .take(a[1].n.min(16))
                    .map(|r| (r.score, r.rb, r.qb, r.qe, r.sub, r.csub, r.secondary, r.secondary_all))
                    .collect::<Vec<_>>(),
            );
        }
        if o > 0 {
            // check if an end has multiple hits even after mate-SW
            let mut is_multi = [0_i32; 2];
            for i in 0..2_usize {
                let mut j = 1_usize;
                while j < n_pri[i] as usize {
                    if a[i].a[j].secondary < 0 && a[i].a[j].score >= opt.T {
                        break;
                    }
                    j += 1;
                }
                is_multi[i] = if j < n_pri[i] as usize { 1 } else { 0 };
            }
            // TODO (upstream): in rare cases, the true hit may be long but with low score
            if is_multi == [0, 0] {
                // compute mapQ for the best SE hit
                let score_un = a[0].a[0].score + a[1].a[0].score - opt.pen_unpaired;
                subo = subo.max(score_un);
                let mut q_pe = raw_mapq(o - subo, opt.a);
                if n_sub > 0 {
                    q_pe -= (4.343 * ((n_sub + 1) as f64).ln() + 0.499) as i32;
                }
                q_pe = q_pe.clamp(0, 60);
                q_pe = (q_pe as f64
                    * (1.0 - 0.5 * (a[0].a[0].frac_rep as f64 + a[1].a[0].frac_rep as f64))
                    + 0.499) as i32;

                // the following assumes no split hits
                let mut q_se = [0_i32; 2];
                if o > score_un {
                    // paired alignment is preferred
                    for i in 0..2_usize {
                        let zi = z[i] as usize;
                        let secondary = a[i].a[zi].secondary;
                        let secondary_score = if secondary >= 0 {
                            Some(a[i].a[secondary as usize].score)
                        } else {
                            None
                        };
                        let c = &mut a[i].a[zi];
                        if let Some(score) = secondary_score {
                            c.sub = score;
                            c.secondary = -2;
                        }
                        q_se[i] = mem_approx_mapq_se(opt, c);
                    }
                    q_se[0] = if q_se[0] > q_pe {
                        q_se[0]
                    } else if q_pe < q_se[0] + 40 {
                        q_pe
                    } else {
                        q_se[0] + 40
                    };
                    q_se[1] = if q_se[1] > q_pe {
                        q_se[1]
                    } else if q_pe < q_se[1] + 40 {
                        q_pe
                    } else {
                        q_se[1] + 40
                    };
                    extra_flag |= 2;
                    // cap at the tandem repeat score
                    for i in 0..2_usize {
                        let c = &a[i].a[z[i] as usize];
                        q_se[i] = q_se[i].min(raw_mapq(c.score - c.csub, opt.a));
                    }
                } else {
                    // the unpaired alignment is preferred
                    z = [0, 0];
                    q_se[0] = mem_approx_mapq_se(opt, &a[0].a[0]);
                    q_se[1] = mem_approx_mapq_se(opt, &a[1].a[0]);
                }
                if debug_trace_read(s[0].name.as_deref()) {
                    eprintln!(
                        "[trace::pe_batch_post2:chosen] read={} id={} z={:?} q_se={:?} extra_flag={} chosen0={:?} chosen1={:?}",
                        s[0].name.as_deref().unwrap_or(""),
                        id,
                        z,
                        q_se,
                        extra_flag,
                        a[0]
                            .a
                            .get(usize::try_from(z[0].max(0)).expect("z0"))
                            .map(|r| (r.score, r.sub, r.csub, r.secondary, r.secondary_all)),
                        a[1]
                            .a
                            .get(usize::try_from(z[1].max(0)).expect("z1"))
                            .map(|r| (r.score, r.sub, r.csub, r.secondary, r.secondary_all)),
                    );
                }

                // switch secondary and primary if both of them are non-ALT
                for i in 0..2_usize {
                    let zi = z[i] as usize;
                    let k = a[i].a[zi].secondary_all;
                    if k >= 0 && k < n_pri[i] {
                        for j in 0..a[i].n {
                            if a[i].a[j].secondary_all == k || j == k as usize {
                                a[i].a[j].secondary_all = z[i];
                            }
                        }
                        a[i].a[zi].secondary_all = -1;
                    }
                }

                let mut xa = if (opt.flag & MEM_F_ALL) == 0 {
                    [
                        mem_gen_alt(
                            opt,
                            bns,
                            pac,
                            &a[0],
                            s[0].l_seq,
                            query_string_for_aln(&s[0]).as_ref(),
                        ),
                        mem_gen_alt(
                            opt,
                            bns,
                            pac,
                            &a[1],
                            s[1].l_seq,
                            query_string_for_aln(&s[1]).as_ref(),
                        ),
                    ]
                } else {
                    [vec![None; a[0].n], vec![None; a[1].n]]
                };

                for i in 0..2_usize {
                    let zi = z[i] as usize;
                    h[i] = mem_reg2aln(
                        opt,
                        bns,
                        pac,
                        s[i].l_seq,
                        query_string_for_aln(&s[i]).as_ref(),
                        Some(&a[i].a[zi]),
                    );
                    h[i].mapq = q_se[i].max(0) as u32;
                    h[i].flag |= (0x40 << i) | extra_flag;
                    // Take XA from xa rather than clone — xa goes out of scope at the end of
                    // this block so the take is safe.
                    h[i].XA = std::mem::take(&mut xa[i][zi]);
                    aa[i].push(h[i].clone());
                    // the read has ALT hits
                    if (n_pri[i] as usize) < a[i].n {
                        let p = &a[i].a[n_pri[i] as usize];
                        if p.score >= opt.T && p.secondary < 0 && p.is_alt != 0 {
                            g[i] = mem_reg2aln(
                                opt,
                                bns,
                                pac,
                                s[i].l_seq,
                                query_string_for_aln(&s[i]).as_ref(),
                                Some(p),
                            );
                            g[i].flag |= 0x800 | (0x40 << i) | extra_flag;
                            g[i].XA = std::mem::take(&mut xa[i][n_pri[i] as usize]);
                            aa[i].push(g[i].clone());
                        }
                    }
                }

                // write SAM: read1 hits, then read2 hits — each record carries the mate's primary
                // aln in `m` so flag/insert-size/MRNM fields can be filled
                for i in 0..aa[0].len() {
                    mem_aln2sam(
                        opt,
                        bns,
                        &mut str_,
                        &s[0],
                        aa[0].len() as i32,
                        &aa[0],
                        i as i32,
                        Some(&h[1]),
                    );
                }
                s[0].sam = Some(take_kstring_boxed(&mut str_));
                for i in 0..aa[1].len() {
                    mem_aln2sam(
                        opt,
                        bns,
                        &mut str_,
                        &s[1],
                        aa[1].len() as i32,
                        &aa[1],
                        i as i32,
                        Some(&h[0]),
                    );
                }
                s[1].sam = Some(take_kstring_boxed(&mut str_));
                assert_eq!(s[0].name, s[1].name, "paired reads have different names");
                drain_pe_aln_pools(&mut h, &mut g, &mut aa);
                return n;
            }
        }
    }

    for i in 0..2_usize {
        let mut which = -1_i32;
        if a[i].n > 0 {
            if a[i].a[0].score >= opt.T {
                which = 0;
            } else if (n_pri[i] as usize) < a[i].n && a[i].a[n_pri[i] as usize].score >= opt.T {
                which = n_pri[i];
            }
        }
        h[i] = if which >= 0 {
            mem_reg2aln(
                opt,
                bns,
                pac,
                s[i].l_seq,
                query_string_for_aln(&s[i]).as_ref(),
                Some(&a[i].a[which as usize]),
            )
        } else {
            mem_reg2aln(
                opt,
                bns,
                pac,
                s[i].l_seq,
                query_string_for_aln(&s[i]).as_ref(),
                None,
            )
        };
    }
    // if the top hits from the two ends constitute a proper pair, flag it.
    if (opt.flag & MEM_F_NOPAIRING) == 0 && h[0].rid == h[1].rid && h[0].rid >= 0 {
        let mut dist = 0_i64;
        let d = mem_infer_dir(bns.l_pac, a[0].a[0].rb, a[1].a[0].rb, &mut dist) as usize;
        if pes[d].failed == 0 && dist >= i64::from(pes[d].low) && dist <= i64::from(pes[d].high) {
            extra_flag |= 2;
        }
    }
    mem_reg2sam(
        opt,
        bns,
        pac,
        &mut s[0],
        &mut a[0],
        0x41 | extra_flag,
        Some(&h[1]),
    );
    mem_reg2sam(
        opt,
        bns,
        pac,
        &mut s[1],
        &mut a[1],
        0x81 | extra_flag,
        Some(&h[0]),
    );
    assert_eq!(s[0].name, s[1].name, "paired reads have different names");
    drain_pe_aln_pools(&mut h, &mut g, &mut aa);
    n
}

/// Batched variant of `mem_matesw`.
///
/// Instead of running `ksw_align2` immediately, emits a `SeqPair` describing the
/// `(mate seq, ref window)` job into the per-thread `seqPairArrayLeft128` / `seqBufLeft` arrays
/// and records the resulting `pcnt` index in `seqPairArrayAux[gcnt + r]` (or `-1` if skipped).
/// `mem_sam_pe_batch` then runs SIMD SW over all enqueued pairs and `mem_matesw_batch_post`
/// consumes the results.
#[doc = "Original function: mem_matesw_batch_pre:930"]
pub fn mem_matesw_batch_pre(
    opt: &mem_opt_t,
    bns: &bntseq_t,
    pac: &[u8],
    pes: &[mem_pestat_t; 4],
    a: &mem_alnreg_t,
    l_ms: i32,
    ms: &[u8],
    ma: &mem_alnreg_v,
    mmc: &mut mem_cache,
    mut pcnt: i32,
    gcnt: i32,
    maxRefLen: &mut i32,
    maxQerLen: &mut i32,
    tid: usize,
) -> i32 {
    let l_pac = bns.l_pac;
    let mut skip = [0_i32; 4];
    for r in 0..4 {
        skip[r] = if pes[r].failed != 0 { 1 } else { 0 };
    }
    // check which orientation has been found
    for i in 0..ma.n {
        let mut dist = 0_i64;
        let r = mem_infer_dir(l_pac, a.rb, ma.a[i].rb, &mut dist) as usize;
        if dist >= i64::from(pes[r].low) && dist <= i64::from(pes[r].high) {
            skip[r] = 1;
        }
    }
    // consistent pair exist; no need to perform SW — write -1 sentinels into the aux slots so
    // mem_matesw_batch_post knows to skip this group.
    if skip.iter().copied().sum::<i32>() == 4 {
        let gcnt_us = gcnt as usize;
        while mmc.seqPairArrayAux[tid].len() < gcnt_us + 4 {
            mmc.seqPairArrayAux[tid].push(SeqPair {
                id: -1,
                ..Default::default()
            });
        }
        for r in 0..4_usize {
            mmc.seqPairArrayAux[tid][gcnt_us + r].id = -1;
        }
        return pcnt;
    }

    let l_ms_usize = l_ms as usize;
    // Reuse thread-local scratch across the 4-direction r loop AND across mem_matesw_batch_pre
    // calls (one per pair × per-batch). query_buf is unused here but we still take/restore the
    // tuple to keep the API consistent across mem_matesw* sites.
    let (mut rev_buf, _query_buf_unused, mut ref_seq_buf) =
        MEM_MATESW_SCRATCH.with(|c| std::mem::take(&mut *c.borrow_mut()));
    for r in 0..4_usize {
        if skip[r] != 0 {
            while mmc.seqPairArrayAux[tid].len() < (gcnt as usize + r + 1) {
                mmc.seqPairArrayAux[tid].push(SeqPair {
                    id: -1,
                    ..Default::default()
                });
            }
            mmc.seqPairArrayAux[tid][gcnt as usize + r].id = -1;
            continue;
        }
        // is_rev: whether to reverse complement the mate
        // is_larger: whether the mate has the larger coordinate
        let is_rev = ((r >> 1) != (r & 1)) as i32;
        let is_larger = ((r >> 1) == 0) as i32;
        let seq: &[u8] = if is_rev != 0 {
            // Reverse-complement: walk ms backward, complement each base via LUT.
            // LUT maps 0,1,2,3,4,5+ -> 3,2,1,0,4,4 (any value >= 4 maps to AMBIG=4).
            const RC_LUT: [u8; 8] = [3, 2, 1, 0, 4, 4, 4, 4];
            rev_buf.clear();
            rev_buf.extend(
                ms[..l_ms_usize]
                    .iter()
                    .rev()
                    .map(|&c| unsafe { *RC_LUT.get_unchecked((c as usize).min(7)) }),
            );
            &rev_buf
        } else {
            &ms[..l_ms_usize]
        };

        let mut rb;
        let mut re;
        if is_rev == 0 {
            // if on the same strand, end position should be larger to make room for the seq length
            rb = if is_larger != 0 {
                a.rb + i64::from(pes[r].low)
            } else {
                a.rb - i64::from(pes[r].high)
            };
            re = if is_larger != 0 {
                a.rb + i64::from(pes[r].high)
            } else {
                a.rb - i64::from(pes[r].low)
            };
            re += i64::from(l_ms);
        } else {
            // similarly on opposite strands: the start position is biased by l_ms instead
            rb = if is_larger != 0 {
                a.rb + i64::from(pes[r].low)
            } else {
                a.rb - i64::from(pes[r].high)
            };
            rb -= i64::from(l_ms);
            re = if is_larger != 0 {
                a.rb + i64::from(pes[r].high)
            } else {
                a.rb - i64::from(pes[r].low)
            };
        }
        if rb < 0 {
            rb = 0;
        }
        if re > (l_pac << 1) {
            re = l_pac << 1;
        }
        let mut rid = -1_i32;
        let ref_ok = if rb < re {
            let mid = (rb + re) >> 1;
            crate::bwa_mem2::bntseq::bns_fetch_seq_into(
                bns,
                pac,
                &mut rb,
                mid,
                &mut re,
                &mut rid,
                &mut ref_seq_buf,
            );
            true
        } else {
            ref_seq_buf.clear();
            false
        };
        let ref_seq: Option<&[u8]> = if ref_ok { Some(&ref_seq_buf) } else { None };
        // no funny things happening: rid matched and the window is at least seed-long
        if a.rid == rid && re - rb >= i64::from(opt.min_seed_len) {
            let xtra = KSW_XSUBO
                | KSW_XSTART
                | if l_ms * opt.a < 250 { KSW_XBYTE } else { 0 }
                | (opt.min_seed_len * opt.a);
            let (ref_offset, qer_offset) = if pcnt != 0 {
                let prev = mmc.seqPairArrayLeft128[tid][(pcnt - 1) as usize];
                (prev.idr + prev.len1, prev.idq + prev.len2)
            } else {
                (0_i32, 0_i32)
            };
            let len1 = (re - rb) as i32;
            let len2 = l_ms;
            *maxRefLen = (*maxRefLen).max(len1);
            *maxQerLen = (*maxQerLen).max(len2);
            let sp = SeqPair {
                h0: xtra,
                idq: qer_offset,
                idr: ref_offset,
                len1,
                len2,
                regid: pcnt,
                ..Default::default()
            };
            if mmc.seqPairArrayLeft128.len() <= tid {
                mmc.seqPairArrayLeft128.resize_with(tid + 1, Vec::new);
            }
            if mmc.seqPairArrayRight128.len() <= tid {
                mmc.seqPairArrayRight128.resize_with(tid + 1, Vec::new);
            }
            if mmc.seqPairArrayAux.len() <= tid {
                mmc.seqPairArrayAux.resize_with(tid + 1, Vec::new);
            }
            let gcnt_us = gcnt as usize;
            while mmc.seqPairArrayAux[tid].len() < gcnt_us + r + 1 {
                mmc.seqPairArrayAux[tid].push(SeqPair {
                    id: -1,
                    ..Default::default()
                });
            }
            mmc.seqPairArrayAux[tid][gcnt_us + r].id = pcnt;
            let ref_start = ref_offset as usize;
            let qer_start = qer_offset as usize;
            let ref_end = ref_start + len1 as usize;
            let qer_end = qer_start + len2 as usize;
            if mmc.seqBufLeftRef[tid].len() < ref_end {
                mmc.seqBufLeftRef[tid].resize(ref_end, 0);
            }
            if mmc.seqBufLeftQer[tid].len() < qer_end {
                mmc.seqBufLeftQer[tid].resize(qer_end, 0);
            }
            mmc.seqBufLeftRef[tid][ref_start..ref_end]
                .copy_from_slice(ref_seq.expect("reference seq"));
            mmc.seqBufLeftQer[tid][qer_start..qer_end].copy_from_slice(seq);
            let sp_index = pcnt as usize;
            if mmc.seqPairArrayLeft128[tid].len() <= sp_index {
                mmc.seqPairArrayLeft128[tid].resize(sp_index + 1, SeqPair::default());
            }
            mmc.seqPairArrayLeft128[tid][sp_index] = sp;
            pcnt += 1;
        } else {
            let gcnt_us = gcnt as usize;
            while mmc.seqPairArrayAux[tid].len() < gcnt_us + r + 1 {
                mmc.seqPairArrayAux[tid].push(SeqPair {
                    id: -1,
                    ..Default::default()
                });
            }
            mmc.seqPairArrayAux[tid][gcnt_us + r].id = -1;
        }
    }
    MEM_MATESW_SCRATCH.with(|c| *c.borrow_mut() = (rev_buf, _query_buf_unused, ref_seq_buf));
    pcnt
}

/// Consume the batched SW results for one `(pair, direction-quadruple)` group.
///
/// For each direction `r`, look up the `pcnt` index recorded by `mem_matesw_batch_pre` in
/// `gar[gcnt + r]`: if non-negative, reuse `myaln[index]`; if `-1` (rare re-routing case), fall
/// back to a synchronous `ksw_align2` here. Successful alignments are inserted into `ma` sorted
/// by descending score.
#[doc = "Original function: mem_matesw_batch_post:1095"]
pub fn mem_matesw_batch_post(
    opt: &mem_opt_t,
    bns: &bntseq_t,
    pac: &[u8],
    pes: &[mem_pestat_t; 4],
    a: &mem_alnreg_t,
    l_ms: i32,
    ms: &[u8],
    ma: &mut mem_alnreg_v,
    myaln: &[kswr_t],
    gcnt: i32,
    gar: &[SeqPair],
) -> i32 {
    let l_pac = bns.l_pac;
    let mut skip = [0_i32; 4];
    for r in 0..4 {
        skip[r] = if pes[r].failed != 0 { 1 } else { 0 };
    }
    // check which orientation has been found
    for i in 0..ma.n {
        let mut dist = 0_i64;
        let r = mem_infer_dir(l_pac, a.rb, ma.a[i].rb, &mut dist) as usize;
        if dist >= i64::from(pes[r].low) && dist <= i64::from(pes[r].high) {
            skip[r] = 1;
        }
    }
    // consistent pair exist; no need to perform SW
    if skip.iter().copied().sum::<i32>() == 4 {
        return 0;
    }

    let mut n = 0_i32;
    let l_ms_usize = l_ms as usize;
    // Reserve up to 4 extra mem_alnreg_t slots in ma.a before potential per-direction push —
    // amortizes the per-rescue realloc that perf showed at 0.62% in mem_matesw_batch_post.
    if ma.a.capacity() < ma.a.len() + 4 {
        ma.a.reserve(ma.a.len() + 4 - ma.a.capacity());
    }
    // Reuse thread-local scratch across the 4-direction r loop AND across mem_matesw_batch_post
    // calls (one per pair × per-batch).
    let (mut rev_buf, mut query_buf, mut ref_seq_buf) =
        MEM_MATESW_SCRATCH.with(|c| std::mem::take(&mut *c.borrow_mut()));
    for r in 0..4_usize {
        if skip[r] != 0 {
            continue;
        }
        // is_rev: whether to reverse complement the mate
        // is_larger: whether the mate has the larger coordinate
        let is_rev = ((r >> 1) != (r & 1)) as i32;
        let is_larger = ((r >> 1) == 0) as i32;
        let seq: &[u8] = if is_rev != 0 {
            // Reverse-complement: walk ms backward, complement each base via LUT.
            // LUT maps 0,1,2,3,4,5+ -> 3,2,1,0,4,4 (any value >= 4 maps to AMBIG=4).
            const RC_LUT: [u8; 8] = [3, 2, 1, 0, 4, 4, 4, 4];
            rev_buf.clear();
            rev_buf.extend(
                ms[..l_ms_usize]
                    .iter()
                    .rev()
                    .map(|&c| unsafe { *RC_LUT.get_unchecked((c as usize).min(7)) }),
            );
            &rev_buf
        } else {
            &ms[..l_ms_usize]
        };
        let mut rb;
        let mut re;
        if is_rev == 0 {
            // if on the same strand, end position should be larger to make room for the seq length
            rb = if is_larger != 0 {
                a.rb + i64::from(pes[r].low)
            } else {
                a.rb - i64::from(pes[r].high)
            };
            re = if is_larger != 0 {
                a.rb + i64::from(pes[r].high)
            } else {
                a.rb - i64::from(pes[r].low)
            };
            re += i64::from(l_ms);
        } else {
            // similarly on opposite strands
            rb = if is_larger != 0 {
                a.rb + i64::from(pes[r].low)
            } else {
                a.rb - i64::from(pes[r].high)
            };
            rb -= i64::from(l_ms);
            re = if is_larger != 0 {
                a.rb + i64::from(pes[r].high)
            } else {
                a.rb - i64::from(pes[r].low)
            };
        }
        if rb < 0 {
            rb = 0;
        }
        if re > (l_pac << 1) {
            re = l_pac << 1;
        }
        let mut rid = -1_i32;
        let ref_ok = if rb < re {
            let mid = (rb + re) >> 1;
            crate::bwa_mem2::bntseq::bns_fetch_seq_into(
                bns,
                pac,
                &mut rb,
                mid,
                &mut re,
                &mut rid,
                &mut ref_seq_buf,
            );
            true
        } else {
            ref_seq_buf.clear();
            false
        };
        // no funny things happening: rid matched and the window is at least seed-long
        if a.rid == rid && re - rb >= i64::from(opt.min_seed_len) {
            let xtra = KSW_XSUBO
                | KSW_XSTART
                | if l_ms * opt.a < 250 { KSW_XBYTE } else { 0 }
                | (opt.min_seed_len * opt.a);
            // Re-routing: encountered -ve index for gcnt+r — fall back to ksw_align2 here.
            let index = gar.get(gcnt as usize + r).map(|sp| sp.id).unwrap_or(-1);
            let aln = if index < 0 {
                debug_assert!(ref_ok, "ksw fallback path requires ref_seq_buf populated");
                query_buf.clear();
                query_buf.extend_from_slice(seq);
                ksw_align2(
                    l_ms,
                    &mut query_buf,
                    (re - rb) as i32,
                    &mut ref_seq_buf,
                    5,
                    &opt.mat,
                    opt.o_del,
                    opt.e_del,
                    opt.o_ins,
                    opt.e_ins,
                    xtra,
                    None,
                )
            } else {
                myaln[index as usize]
            };
            // something goes wrong if aln.qb < 0
            if aln.score >= opt.min_seed_len && aln.qb >= 0 {
                let b = mem_alnreg_t {
                    rid: a.rid,
                    is_alt: a.is_alt,
                    qb: if is_rev != 0 {
                        l_ms - (aln.qe + 1)
                    } else {
                        aln.qb
                    },
                    qe: if is_rev != 0 {
                        l_ms - aln.qb
                    } else {
                        aln.qe + 1
                    },
                    rb: if is_rev != 0 {
                        (l_pac << 1) - (rb + i64::from(aln.te) + 1)
                    } else {
                        rb + i64::from(aln.tb)
                    },
                    re: if is_rev != 0 {
                        (l_pac << 1) - (rb + i64::from(aln.tb))
                    } else {
                        rb + i64::from(aln.te) + 1
                    },
                    score: aln.score,
                    csub: aln.score2,
                    secondary: -1,
                    seedcov: (((if is_rev != 0 {
                        (l_pac << 1) - (rb + i64::from(aln.tb))
                    } else {
                        rb + i64::from(aln.te) + 1
                    }) - (if is_rev != 0 {
                        (l_pac << 1) - (rb + i64::from(aln.te) + 1)
                    } else {
                        rb + i64::from(aln.tb)
                    }))
                    .min(i64::from(
                        if is_rev != 0 {
                            l_ms - aln.qb
                        } else {
                            aln.qe + 1
                        } - if is_rev != 0 {
                            l_ms - (aln.qe + 1)
                        } else {
                            aln.qb
                        },
                    )) >> 1) as i32,
                    ..Default::default()
                };
                // make room for a new element, then move b s.t. ma stays sorted by score
                ma.a.push(b);
                ma.n = ma.a.len();
                ma.m = ma.a.len();

                // find the insertion point
                let mut insert_at = ma.n - 1;
                for i in 0..ma.n - 1 {
                    if ma.a[i].score < b.score {
                        insert_at = i;
                        break;
                    }
                }
                if insert_at < ma.n - 1 {
                    ma.a[insert_at..].rotate_right(1);
                    ma.a[insert_at] = b;
                }
            }
            n += 1;
        }
        if n != 0 {
            ma.n = mem_sort_dedup_patch(opt, None, None, None, ma.n as i32, &mut ma.a) as usize;
            ma.m = ma.a.len();
        }
    }
    MEM_MATESW_SCRATCH.with(|c| *c.borrow_mut() = (rev_buf, query_buf, ref_seq_buf));
    n
}

#[cfg(test)]
mod tests {
    use super::{
        cal_sub, max_matesw_rescue_limit, mem_infer_dir, mem_matesw, mem_pair, mem_pestat,
        mem_sam_pe, mem_sam_pe_batch, mem_sam_pe_batch_post, mem_sam_pe_batch_pre, revseq,
    };
    use crate::bwa_mem2::bandedswa::SeqPair;
    use crate::bwa_mem2::bntseq::{bntann1_t, bntseq_t};
    use crate::bwa_mem2::bwa::bseq1_t;
    use crate::bwa_mem2::bwamem::mem_opt_init;
    use crate::bwa_mem2::bwamem::{mem_alnreg_t, mem_alnreg_v, mem_cache, mem_pestat_t};
    use crate::bwa_mem2::ksw::{kswr_t, KSW_XSTART};

    fn pack_seq(seq: &[u8]) -> Vec<u8> {
        let mut pac = vec![0_u8; (seq.len() + 3) / 4];
        for (i, &base) in seq.iter().enumerate() {
            let shift = (((!i64::try_from(i).expect("i")) & 3) << 1) as u8;
            pac[i >> 2] |= base << shift;
        }
        pac
    }

    #[test]
    fn max_matesw_rescue_limit_clamps_before_usize_conversion() {
        let mut opt = (*mem_opt_init()).clone();

        opt.max_matesw = -3;
        assert_eq!(max_matesw_rescue_limit(&opt), 0);

        opt.max_matesw = 0;
        assert_eq!(max_matesw_rescue_limit(&opt), 0);

        opt.max_matesw = 7;
        assert_eq!(max_matesw_rescue_limit(&opt), 7);
    }

    #[test]
    fn mem_infer_dir_matches_forward_reverse_orientation_cases() {
        let mut dist = 0_i64;
        assert_eq!(mem_infer_dir(100, 10, 30, &mut dist), 0);
        assert_eq!(dist, 20);
        assert_eq!(mem_infer_dir(100, 10, 170, &mut dist), 1);
        assert_eq!(dist, 19);
    }

    #[test]
    fn cal_sub_returns_first_significant_overlap_or_seed_floor() {
        let opt = mem_opt_init();
        let r = mem_alnreg_v {
            n: 2,
            m: 2,
            a: vec![
                mem_alnreg_t {
                    qb: 0,
                    qe: 10,
                    score: 20,
                    ..Default::default()
                },
                mem_alnreg_t {
                    qb: 2,
                    qe: 9,
                    score: 15,
                    ..Default::default()
                },
            ],
        };
        assert_eq!(cal_sub(&opt, &r), 15);

        let r = mem_alnreg_v {
            n: 1,
            m: 1,
            a: vec![mem_alnreg_t {
                qb: 0,
                qe: 10,
                score: 20,
                ..Default::default()
            }],
        };
        assert_eq!(cal_sub(&opt, &r), opt.min_seed_len * opt.a);
    }

    #[test]
    fn mem_pestat_infers_supported_orientation_distribution() {
        let mut opt = (*mem_opt_init()).clone();
        opt.max_ins = 1000;
        let mut regs = Vec::new();
        for i in 0..12 {
            let left = mem_alnreg_v {
                n: 1,
                m: 1,
                a: vec![mem_alnreg_t {
                    rid: 0,
                    rb: 100 + i * 10,
                    re: 150 + i * 10,
                    qb: 0,
                    qe: 50,
                    score: 100,
                    ..Default::default()
                }],
            };
            let right = mem_alnreg_v {
                n: 1,
                m: 1,
                a: vec![mem_alnreg_t {
                    rid: 0,
                    rb: 160 + i * 10,
                    re: 210 + i * 10,
                    qb: 0,
                    qe: 50,
                    score: 100,
                    ..Default::default()
                }],
            };
            regs.push(left);
            regs.push(right);
        }
        let mut pes = [mem_pestat_t::default(); 4];
        mem_pestat(&opt, 1_000, regs.len() as i32, &regs, &mut pes);
        assert_eq!(pes[0].failed, 0);
        assert!(pes[0].low >= 1);
        assert!(pes[0].high >= pes[0].low);
        assert!(pes[0].avg > 0.0);
    }

    #[test]
    fn mem_matesw_skips_when_consistent_pair_already_exists() {
        let opt = mem_opt_init();
        let bns = bntseq_t {
            l_pac: 1000,
            n_seqs: 1,
            anns: vec![bntann1_t {
                offset: 0,
                len: 1000,
                name: "chr1".into(),
                ..Default::default()
            }],
            ..Default::default()
        };
        let pac = pack_seq(&vec![0_u8; 1000]);
        let pes = [
            mem_pestat_t {
                low: 10,
                high: 30,
                failed: 0,
                avg: 20.0,
                std: 1.0,
            },
            mem_pestat_t {
                failed: 1,
                ..Default::default()
            },
            mem_pestat_t {
                failed: 1,
                ..Default::default()
            },
            mem_pestat_t {
                failed: 1,
                ..Default::default()
            },
        ];
        let a = mem_alnreg_t {
            rid: 0,
            rb: 100,
            re: 110,
            ..Default::default()
        };
        let mut ma = mem_alnreg_v {
            n: 1,
            m: 1,
            a: vec![mem_alnreg_t {
                rid: 0,
                rb: 120,
                re: 130,
                score: 10,
                ..Default::default()
            }],
        };
        let ms = vec![0_u8; 20];
        assert_eq!(mem_matesw(&opt, &bns, &pac, &pes, &a, 20, &ms, &mut ma), 0);
        assert_eq!(ma.n, 1);
    }

    #[test]
    fn mem_matesw_adds_rescued_alignment_for_matching_window() {
        let mut opt = (*mem_opt_init()).clone();
        opt.min_seed_len = 2;
        let ref_seq = vec![0_u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let bns = bntseq_t {
            l_pac: i64::try_from(ref_seq.len()).expect("len"),
            n_seqs: 1,
            anns: vec![bntann1_t {
                offset: 0,
                len: i32::try_from(ref_seq.len()).expect("len"),
                name: "chr1".into(),
                ..Default::default()
            }],
            ..Default::default()
        };
        let pac = pack_seq(&ref_seq);
        let pes = [
            mem_pestat_t {
                low: 2,
                high: 6,
                failed: 0,
                avg: 4.0,
                std: 1.0,
            },
            mem_pestat_t {
                failed: 1,
                ..Default::default()
            },
            mem_pestat_t {
                failed: 1,
                ..Default::default()
            },
            mem_pestat_t {
                failed: 1,
                ..Default::default()
            },
        ];
        let a = mem_alnreg_t {
            rid: 0,
            rb: 0,
            re: 4,
            is_alt: 0,
            ..Default::default()
        };
        let mut ma = mem_alnreg_v::default();
        let ms = vec![2_u8, 3, 0, 1];
        let rescued = mem_matesw(&opt, &bns, &pac, &pes, &a, 4, &ms, &mut ma);
        assert!(rescued > 0);
        assert!(ma.n > 0);
        assert_eq!(ma.a[0].rid, 0);
        assert!(ma.a[0].score >= opt.min_seed_len);
    }

    #[test]
    fn mem_pair_finds_best_proper_pair_and_reports_suboptimal_count() {
        let opt = mem_opt_init();
        let bns = bntseq_t {
            l_pac: 1000,
            n_seqs: 1,
            anns: vec![bntann1_t {
                offset: 0,
                len: 1000,
                name: "chr1".into(),
                ..Default::default()
            }],
            ..Default::default()
        };
        let pes = [
            mem_pestat_t {
                failed: 1,
                ..Default::default()
            },
            mem_pestat_t {
                low: 15,
                high: 30,
                failed: 0,
                avg: 20.0,
                std: 3.0,
            },
            mem_pestat_t {
                failed: 1,
                ..Default::default()
            },
            mem_pestat_t {
                failed: 1,
                ..Default::default()
            },
        ];
        let _s = [bseq1_t::default(), bseq1_t::default()];
        let mut regs = [
            mem_alnreg_v {
                n: 2,
                m: 2,
                a: vec![
                    mem_alnreg_t {
                        rid: 0,
                        rb: 1000 + 100,
                        score: 40,
                        ..Default::default()
                    },
                    mem_alnreg_t {
                        rid: 0,
                        rb: 1000 + 150,
                        score: 20,
                        ..Default::default()
                    },
                ],
            },
            mem_alnreg_v {
                n: 2,
                m: 2,
                a: vec![
                    mem_alnreg_t {
                        rid: 0,
                        rb: 880,
                        score: 39,
                        ..Default::default()
                    },
                    mem_alnreg_t {
                        rid: 0,
                        rb: 840,
                        score: 15,
                        ..Default::default()
                    },
                ],
            },
        ];
        let n_pri = [2, 2];
        let mut sub = -1;
        let mut n_sub = -1;
        let mut z = [-1, -1];
        let ret = mem_pair(
            &opt,
            &bns,
            &[],
            &pes,
            &mut regs,
            7,
            &mut sub,
            &mut n_sub,
            &mut z,
            &n_pri,
        );
        assert!(ret > 0);
        assert_eq!(z, [0, 0]);
        assert!(sub >= 0);
        assert!(n_sub >= 0);
    }

    #[test]
    fn mem_pair_returns_zero_when_no_supported_pair_exists() {
        let opt = mem_opt_init();
        let bns = bntseq_t {
            l_pac: 1000,
            n_seqs: 1,
            anns: vec![bntann1_t {
                offset: 0,
                len: 1000,
                name: "chr1".into(),
                ..Default::default()
            }],
            ..Default::default()
        };
        let pes = [mem_pestat_t {
            failed: 1,
            ..Default::default()
        }; 4];
        let _s = [bseq1_t::default(), bseq1_t::default()];
        let mut regs = [
            mem_alnreg_v {
                n: 1,
                m: 1,
                a: vec![mem_alnreg_t {
                    rid: 0,
                    rb: 100,
                    score: 20,
                    ..Default::default()
                }],
            },
            mem_alnreg_v {
                n: 1,
                m: 1,
                a: vec![mem_alnreg_t {
                    rid: 0,
                    rb: 200,
                    score: 20,
                    ..Default::default()
                }],
            },
        ];
        let mut sub = -1;
        let mut n_sub = -1;
        let mut z = [-1, -1];
        let ret = mem_pair(
            &opt,
            &bns,
            &[],
            &pes,
            &mut regs,
            1,
            &mut sub,
            &mut n_sub,
            &mut z,
            &[1, 1],
        );
        assert_eq!(ret, 0);
        assert_eq!(sub, 0);
        assert_eq!(n_sub, 0);
    }

    #[test]
    fn mem_sam_pe_no_pairing_falls_back_to_single_end_records() {
        let mut opt = (*mem_opt_init()).clone();
        opt.flag |= super::MEM_F_NOPAIRING;
        opt.T = 1;
        let ref_seq = vec![0_u8, 1, 2, 3, 0, 1, 2, 3];
        let bns = bntseq_t {
            l_pac: 8,
            n_seqs: 1,
            anns: vec![bntann1_t {
                offset: 0,
                len: 8,
                name: "chr1".into(),
                ..Default::default()
            }],
            ..Default::default()
        };
        let pac = pack_seq(&ref_seq);
        let pes = [mem_pestat_t::default(); 4];
        let mut s = [
            bseq1_t {
                l_seq: 4,
                name: Some("pair1".into()),
                seq: Some("ACGT".into()),
                qual: Some("IIII".into()),
                ..Default::default()
            },
            bseq1_t {
                l_seq: 4,
                name: Some("pair1".into()),
                seq: Some("ACGT".into()),
                qual: Some("JJJJ".into()),
                ..Default::default()
            },
        ];
        let mut a = [
            mem_alnreg_v {
                n: 1,
                m: 1,
                a: vec![mem_alnreg_t {
                    rid: 0,
                    rb: 0,
                    re: 4,
                    qb: 0,
                    qe: 4,
                    score: 4,
                    truesc: 4,
                    w: 100,
                    ..Default::default()
                }],
            },
            mem_alnreg_v {
                n: 1,
                m: 1,
                a: vec![mem_alnreg_t {
                    rid: 0,
                    rb: 4,
                    re: 8,
                    qb: 0,
                    qe: 4,
                    score: 4,
                    truesc: 4,
                    w: 100,
                    ..Default::default()
                }],
            },
        ];
        let rescued = mem_sam_pe(&opt, &bns, &pac, &pes, 1, &mut s, &mut a);
        assert_eq!(rescued, 0);
        assert!(
            s[0].sam.as_deref().expect("sam").starts_with("pair1\t65\t"),
            "{}",
            s[0].sam.as_deref().unwrap()
        );
        assert!(
            s[1].sam
                .as_deref()
                .expect("sam")
                .starts_with("pair1\t129\t"),
            "{}",
            s[1].sam.as_deref().unwrap()
        );
    }

    #[test]
    fn mem_sam_pe_writes_preferred_proper_pair_records() {
        let mut opt = (*mem_opt_init()).clone();
        opt.T = 1;
        let ref_seq = vec![0_u8, 1, 2, 3, 0, 1, 2, 3];
        let bns = bntseq_t {
            l_pac: 8,
            n_seqs: 1,
            anns: vec![bntann1_t {
                offset: 0,
                len: 8,
                name: "chr1".into(),
                ..Default::default()
            }],
            ..Default::default()
        };
        let pac = pack_seq(&ref_seq);
        let pes = [
            mem_pestat_t {
                failed: 1,
                ..Default::default()
            },
            mem_pestat_t {
                low: 1,
                high: 10,
                failed: 0,
                avg: 3.0,
                std: 1.0,
            },
            mem_pestat_t {
                failed: 1,
                ..Default::default()
            },
            mem_pestat_t {
                failed: 1,
                ..Default::default()
            },
        ];
        let mut s = [
            bseq1_t {
                l_seq: 4,
                name: Some("pair2".into()),
                seq: Some("ACGT".into()),
                qual: Some("IIII".into()),
                ..Default::default()
            },
            bseq1_t {
                l_seq: 4,
                name: Some("pair2".into()),
                seq: Some("ACGT".into()),
                qual: Some("JJJJ".into()),
                ..Default::default()
            },
        ];
        let mut a = [
            mem_alnreg_v {
                n: 1,
                m: 1,
                a: vec![mem_alnreg_t {
                    rid: 0,
                    rb: 8 + 4,
                    re: 8 + 8,
                    qb: 0,
                    qe: 4,
                    score: 40,
                    truesc: 4,
                    w: 100,
                    ..Default::default()
                }],
            },
            mem_alnreg_v {
                n: 1,
                m: 1,
                a: vec![mem_alnreg_t {
                    rid: 0,
                    rb: 0,
                    re: 4,
                    qb: 0,
                    qe: 4,
                    score: 39,
                    truesc: 4,
                    w: 100,
                    ..Default::default()
                }],
            },
        ];
        let rescued = mem_sam_pe(&opt, &bns, &pac, &pes, 2, &mut s, &mut a);
        assert_eq!(rescued, 0);
        let sam0 = s[0].sam.as_deref().expect("sam0");
        let sam1 = s[1].sam.as_deref().expect("sam1");
        assert!(
            sam0.starts_with("pair2\t99\t") || sam0.starts_with("pair2\t83\t"),
            "{sam0}"
        );
        assert!(
            sam1.starts_with("pair2\t147\t") || sam1.starts_with("pair2\t163\t"),
            "{sam1}"
        );
    }

    #[test]
    fn revseq_reverses_prefix_in_place() {
        let mut seq = vec![1_u8, 2, 3, 4, 5];
        revseq(4, &mut seq);
        assert_eq!(seq, vec![4, 3, 2, 1, 5]);
    }

    #[test]
    fn mem_sam_pe_batch_path_matches_scalar_mem_sam_pe() {
        let mut opt = (*mem_opt_init()).clone();
        opt.min_seed_len = 2;
        let ref_seq = vec![0_u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let pac = pack_seq(&ref_seq);
        let bns = bntseq_t {
            l_pac: i64::try_from(ref_seq.len()).expect("len"),
            n_seqs: 1,
            anns: vec![bntann1_t {
                offset: 0,
                len: i32::try_from(ref_seq.len()).expect("len"),
                name: "chr1".into(),
                ..Default::default()
            }],
            ..Default::default()
        };
        let pes = [
            mem_pestat_t {
                low: 3,
                high: 8,
                failed: 0,
                avg: 5.0,
                std: 1.0,
            },
            mem_pestat_t {
                failed: 1,
                ..Default::default()
            },
            mem_pestat_t {
                failed: 1,
                ..Default::default()
            },
            mem_pestat_t {
                failed: 1,
                ..Default::default()
            },
        ];
        let s0 = bseq1_t {
            name: Some("pairb".into()),
            seq: Some("ACGT".into()),
            qual: Some("IIII".into()),
            l_seq: 4,
            ..Default::default()
        };
        let s1 = bseq1_t {
            name: Some("pairb".into()),
            seq: Some("CGTA".into()),
            qual: Some("JJJJ".into()),
            l_seq: 4,
            ..Default::default()
        };
        let a0 = mem_alnreg_v {
            n: 1,
            m: 1,
            a: vec![mem_alnreg_t {
                rid: 0,
                rb: 0,
                re: 4,
                qb: 0,
                qe: 4,
                score: 20,
                seedcov: 4,
                secondary: -1,
                secondary_all: -1,
                ..Default::default()
            }],
        };
        let a1 = mem_alnreg_v {
            n: 1,
            m: 1,
            a: vec![mem_alnreg_t {
                rid: 0,
                rb: 5,
                re: 9,
                qb: 0,
                qe: 4,
                score: 18,
                seedcov: 4,
                secondary: -1,
                secondary_all: -1,
                ..Default::default()
            }],
        };

        let mut scalar_s = [s0.clone(), s1.clone()];
        let mut scalar_a = [a0.clone(), a1.clone()];
        mem_sam_pe(&opt, &bns, &pac, &pes, 7, &mut scalar_s, &mut scalar_a);

        let mut batch_s = [s0, s1];
        let mut batch_a = [a0, a1];
        let mut mmc = mem_cache {
            seqPairArrayAux: vec![Vec::new()],
            seqPairArrayLeft128: vec![Vec::new()],
            seqPairArrayRight128: vec![Vec::new()],
            wsize: vec![0],
            wsize_buf_ref: vec![0],
            wsize_buf_qer: vec![0],
            seqBufLeftRef: vec![Vec::new()],
            seqBufRightRef: vec![Vec::new()],
            seqBufLeftQer: vec![Vec::new()],
            seqBufRightQer: vec![Vec::new()],
            ..Default::default()
        };
        let mut pcnt = 0_i32;
        let mut gcnt = 0_i32;
        let mut max_ref_len = 0_i32;
        let mut max_qer_len = 0_i32;
        mem_sam_pe_batch_pre(
            &opt,
            &bns,
            &pac,
            &pes,
            7,
            &batch_s,
            &batch_a,
            &mut mmc,
            &mut pcnt,
            &mut gcnt,
            &mut max_ref_len,
            &mut max_qer_len,
            0,
        );
        let mut aln = vec![kswr_t::default(); usize::try_from(pcnt.max(1)).expect("pcnt")];
        mem_sam_pe_batch(
            &opt,
            &mut mmc,
            pcnt,
            0,
            &mut aln,
            max_ref_len,
            max_qer_len,
            0,
        );
        gcnt = 0;
        mem_sam_pe_batch_post(
            &opt,
            &bns,
            &pac,
            &pes,
            7,
            &mut batch_s,
            &mut batch_a,
            &aln,
            &mut mmc,
            &mut gcnt,
            0,
        );

        assert_eq!(batch_s[0].sam, scalar_s[0].sam);
        assert_eq!(batch_s[1].sam, scalar_s[1].sam);
    }

    #[test]
    fn mem_sam_pe_batch_phase1_trims_reference_to_forward_te() {
        let opt = (*mem_opt_init()).clone();
        let mut mmc = mem_cache {
            seqPairArrayLeft128: vec![vec![SeqPair {
                idr: 0,
                idq: 0,
                id: 0,
                len1: 8,
                len2: 4,
                h0: KSW_XSTART,
                regid: 0,
                ..Default::default()
            }]],
            seqBufLeftRef: vec![vec![0, 1, 2, 3, 3, 2, 1, 0]],
            seqBufLeftQer: vec![vec![0, 1, 2, 3]],
            ..Default::default()
        };
        let mut aln = vec![kswr_t::default(); 1];

        mem_sam_pe_batch(&opt, &mut mmc, 1, 1, &mut aln, 8, 4, 0);

        assert_eq!(aln[0].score, 4);
        assert_eq!(aln[0].te, 3);
        assert_eq!(aln[0].qe, 3);
        assert_eq!(aln[0].tb, 0);
        assert_eq!(aln[0].qb, 0);
    }

    #[test]
    fn mem_sam_pe_batch_path_matches_scalar_with_preallocated_worker_buffers() {
        let mut opt = (*mem_opt_init()).clone();
        opt.flag = 0;
        opt.T = 0;
        opt.pen_unpaired = 100;
        opt.max_matesw = 8;
        opt.min_seed_len = 2;
        let ref_seq = vec![0_u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let bns = bntseq_t {
            l_pac: i64::try_from(ref_seq.len()).expect("len"),
            n_seqs: 1,
            anns: vec![bntann1_t {
                offset: 0,
                len: i32::try_from(ref_seq.len()).expect("len"),
                name: "chr1".into(),
                ..Default::default()
            }],
            ..Default::default()
        };
        let pac = pack_seq(&ref_seq);
        let pes = [
            mem_pestat_t {
                low: 1,
                high: 8,
                failed: 0,
                avg: 4.0,
                std: 1.0,
            },
            mem_pestat_t {
                failed: 1,
                ..Default::default()
            },
            mem_pestat_t {
                failed: 1,
                ..Default::default()
            },
            mem_pestat_t {
                failed: 1,
                ..Default::default()
            },
        ];
        let s0 = bseq1_t {
            name: Some("pairc".into()),
            seq: Some("ACGT".into()),
            qual: Some("IIII".into()),
            l_seq: 4,
            ..Default::default()
        };
        let s1 = bseq1_t {
            name: Some("pairc".into()),
            seq: Some("CGTA".into()),
            qual: Some("JJJJ".into()),
            l_seq: 4,
            ..Default::default()
        };
        let a0 = mem_alnreg_v {
            n: 1,
            m: 1,
            a: vec![mem_alnreg_t {
                rid: 0,
                rb: 0,
                re: 4,
                qb: 0,
                qe: 4,
                score: 20,
                seedcov: 4,
                secondary: -1,
                secondary_all: -1,
                ..Default::default()
            }],
        };
        let a1 = mem_alnreg_v {
            n: 1,
            m: 1,
            a: vec![mem_alnreg_t {
                rid: 0,
                rb: 5,
                re: 9,
                qb: 0,
                qe: 4,
                score: 18,
                seedcov: 4,
                secondary: -1,
                secondary_all: -1,
                ..Default::default()
            }],
        };

        let mut scalar_s = [s0.clone(), s1.clone()];
        let mut scalar_a = [a0.clone(), a1.clone()];
        mem_sam_pe(&opt, &bns, &pac, &pes, 9, &mut scalar_s, &mut scalar_a);

        let mut batch_s = [s0, s1];
        let mut batch_a = [a0, a1];
        let mut mmc = mem_cache {
            seqPairArrayAux: vec![vec![SeqPair::default(); 128]],
            seqPairArrayLeft128: vec![vec![SeqPair::default(); 128]],
            seqPairArrayRight128: vec![vec![SeqPair::default(); 128]],
            wsize: vec![128],
            wsize_buf_ref: vec![1024],
            wsize_buf_qer: vec![1024],
            seqBufLeftRef: vec![vec![0; 1024]],
            seqBufRightRef: vec![vec![0; 1024]],
            seqBufLeftQer: vec![vec![0; 1024]],
            seqBufRightQer: vec![vec![0; 1024]],
            ..Default::default()
        };
        let mut pcnt = 0_i32;
        let mut gcnt = 0_i32;
        let mut max_ref_len = 0_i32;
        let mut max_qer_len = 0_i32;
        mem_sam_pe_batch_pre(
            &opt,
            &bns,
            &pac,
            &pes,
            9,
            &batch_s,
            &batch_a,
            &mut mmc,
            &mut pcnt,
            &mut gcnt,
            &mut max_ref_len,
            &mut max_qer_len,
            0,
        );
        assert!(pcnt > 0);
        assert!(mmc.seqPairArrayLeft128[0][0].len1 > 0);
        let mut aln = vec![kswr_t::default(); usize::try_from(pcnt.max(1)).expect("pcnt")];
        mem_sam_pe_batch(
            &opt,
            &mut mmc,
            pcnt,
            0,
            &mut aln,
            max_ref_len,
            max_qer_len,
            0,
        );
        assert!(aln
            .iter()
            .take(usize::try_from(pcnt).expect("pcnt"))
            .any(|x| x.score > 0));
        gcnt = 0;
        mem_sam_pe_batch_post(
            &opt,
            &bns,
            &pac,
            &pes,
            9,
            &mut batch_s,
            &mut batch_a,
            &aln,
            &mut mmc,
            &mut gcnt,
            0,
        );

        assert_eq!(batch_s[0].sam, scalar_s[0].sam);
        assert_eq!(batch_s[1].sam, scalar_s[1].sam);
    }
}
