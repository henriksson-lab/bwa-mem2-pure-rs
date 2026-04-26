#![allow(dead_code, non_snake_case, non_camel_case_types, non_upper_case_globals)]

//! Generated scaffold for `bwa-mem2/src/bwamem_pair.cpp`.

use crate::generated::bandedswa_h::SeqPair;
use crate::generated::bntseq_cpp::bns_fetch_seq;
use crate::generated::bntseq_h::bntseq_t;
use crate::generated::bwa_h::bseq1_t;
use crate::generated::bwamem_cpp::{
    mem_aln2sam, mem_approx_mapq_se, mem_mark_primary_se, mem_reg2aln, mem_reg2sam,
    mem_reorder_primary5, mem_sort_dedup_patch, nt4_cow,
};
use crate::generated::bwamem_extra_cpp::mem_gen_alt;
use crate::generated::bwamem_h::{
    mem_aln_t, mem_alnreg_t, mem_alnreg_v, mem_cache, mem_opt_t, mem_pestat_t,
};
use crate::generated::kstring_h::kstring_t;
use crate::generated::ksw_cpp::ksw_align2;
use crate::generated::ksw_h::{kswr_t, KSW_XBYTE, KSW_XSTART, KSW_XSTOP, KSW_XSUBO};
use crate::generated::kswv_h::kswv;
use crate::generated::utils_h::{hash_64, pair64_t};

const MIN_RATIO: f64 = 0.8;
const MIN_DIR_CNT: usize = 10;
const MIN_DIR_RATIO: f64 = 0.05;
const OUTLIER_BOUND: f64 = 2.0;
const MAPPING_BOUND: f64 = 3.0;
const MAX_STDDEV: f64 = 4.0;
const MEM_F_ALL: i32 = 0x8;
const MEM_F_NOPAIRING: i32 = 0x4;
const MEM_F_PRIMARY5: i32 = 0x800;

fn debug_trace_read(name: Option<&str>) -> bool {
    let Ok(target) = std::env::var("BWA_MEM2_RS_TRACE_READ") else {
        return false;
    };
    name == Some(target.as_str())
}

fn erfc_approx(x: f64) -> f64 {
    // Match libm's erfc precisely (used by upstream C++ via <math.h>); a low-precision
    // Chebyshev approximation here was producing ~1-unit differences in mem_pair's
    // pair-score quantization, which propagated to MAPQ divergences vs upstream.
    extern "C" {
        fn erfc(x: f64) -> f64;
    }
    unsafe { erfc(x) }
}

fn raw_mapq(diff: i32, a: i32) -> i32 {
    (6.02 * (diff as f64) / (a as f64) + 0.499) as i32
}

#[doc = "Original function: mem_infer_dir:58"]
pub fn mem_infer_dir(l_pac: i64, b1: i64, b2: i64, dist: &mut i64) -> i32 {
    let r1 = b1 >= l_pac;
    let r2 = b2 >= l_pac;
    let p2 = if r1 == r2 { b2 } else { (l_pac << 1) - 1 - b2 };
    *dist = if p2 > b1 { p2 - b1 } else { b1 - p2 };
    (if r1 == r2 { 0 } else { 1 }) ^ (if p2 > b1 { 0 } else { 3 })
}

#[doc = "Original function: cal_sub:67"]
pub fn cal_sub(opt: &mem_opt_t, r: &mem_alnreg_v) -> i32 {
    for j in 1..r.n {
        let b_max = r.a[j].qb.max(r.a[0].qb);
        let e_min = r.a[j].qe.min(r.a[0].qe);
        if e_min > b_max {
            let min_l = (r.a[j].qe - r.a[j].qb).min(r.a[0].qe - r.a[0].qb);
            if (e_min - b_max) as f32 >= (min_l as f32) * opt.mask_level {
                return r.a[j].score;
            }
        }
    }
    opt.min_seed_len * opt.a
}

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

    for i in 0..(usize::try_from(n).expect("n") >> 1) {
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
        if r0.a[0].rid != r1.a[0].rid {
            continue;
        }
        let mut is = 0_i64;
        let dir = usize::try_from(mem_infer_dir(l_pac, r0.a[0].rb, r1.a[0].rb, &mut is)).expect("dir");
        if is != 0 && is <= i64::from(opt.max_ins) {
            isize[dir].push(u64::try_from(is).expect("insert size"));
        }
    }

    for d in 0..4 {
        let q = &mut isize[d];
        let r = &mut pes[d];
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
        r.std = (inside.iter().map(|x| (x - r.avg) * (x - r.avg)).sum::<f64>() / inside.len() as f64).sqrt();

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

    for i in 0..ma.n {
        let mut dist = 0_i64;
        let r = usize::try_from(mem_infer_dir(l_pac, a.rb, ma.a[i].rb, &mut dist)).expect("dir");
        if dist >= i64::from(pes[r].low) && dist <= i64::from(pes[r].high) {
            skip[r] = 1;
        }
    }
    if skip.iter().copied().sum::<i32>() == 4 {
        return 0;
    }

    let mut n = 0_i32;
    let l_ms_usize = usize::try_from(l_ms).expect("l_ms");
    for r in 0..4 {
        if skip[r] != 0 {
            continue;
        }
        let is_rev = ((r >> 1) != (r & 1)) as i32;
        let is_larger = ((r >> 1) == 0) as i32;
        let seq = if is_rev != 0 {
            let mut rev = vec![0_u8; l_ms_usize];
            for i in 0..l_ms_usize {
                rev[l_ms_usize - 1 - i] = if ms[i] < 4 { 3 - ms[i] } else { 4 };
            }
            rev
        } else {
            ms[..l_ms_usize].to_vec()
        };

        let mut rb;
        let mut re;
        if is_rev == 0 {
            rb = if is_larger != 0 { a.rb + i64::from(pes[r].low) } else { a.rb - i64::from(pes[r].high) };
            re = if is_larger != 0 { a.rb + i64::from(pes[r].high) } else { a.rb - i64::from(pes[r].low) };
            re += i64::from(l_ms);
        } else {
            rb = if is_larger != 0 { a.rb + i64::from(pes[r].low) } else { a.rb - i64::from(pes[r].high) };
            rb -= i64::from(l_ms);
            re = if is_larger != 0 { a.rb + i64::from(pes[r].high) } else { a.rb - i64::from(pes[r].low) };
        }
        if rb < 0 {
            rb = 0;
        }
        if re > (l_pac << 1) {
            re = l_pac << 1;
        }
        let mut rid = -1_i32;
        let ref_seq = if rb < re {
            let mid = (rb + re) >> 1;
            Some(bns_fetch_seq(bns, pac, &mut rb, mid, &mut re, &mut rid))
        } else {
            None
        };

        if a.rid == rid && re - rb >= i64::from(opt.min_seed_len) {
            let xtra = KSW_XSUBO | KSW_XSTART | if l_ms * opt.a < 250 { KSW_XBYTE } else { 0 } | (opt.min_seed_len * opt.a);
            let mut query = seq.clone();
            let mut target = ref_seq.expect("reference seq");
            let aln = ksw_align2(
                l_ms,
                &mut query,
                i32::try_from(re - rb).expect("ref span"),
                &mut target,
                5,
                &opt.mat,
                opt.o_del,
                opt.e_del,
                opt.o_ins,
                opt.e_ins,
                xtra,
                None,
            );
            if aln.score >= opt.min_seed_len && aln.qb >= 0 {
                let mut b = mem_alnreg_t::default();
                b.rid = a.rid;
                b.is_alt = a.is_alt;
                b.qb = if is_rev != 0 { l_ms - (aln.qe + 1) } else { aln.qb };
                b.qe = if is_rev != 0 { l_ms - aln.qb } else { aln.qe + 1 };
                b.rb = if is_rev != 0 { (l_pac << 1) - (rb + i64::from(aln.te) + 1) } else { rb + i64::from(aln.tb) };
                b.re = if is_rev != 0 { (l_pac << 1) - (rb + i64::from(aln.tb)) } else { rb + i64::from(aln.te) + 1 };
                b.score = aln.score;
                b.csub = aln.score2;
                b.secondary = -1;
                b.seedcov = i32::try_from((b.re - b.rb).min(i64::from(b.qe - b.qb)) >> 1).expect("seedcov");

                ma.a.push(b);
                ma.n = ma.a.len();
                ma.m = ma.a.len();

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
            ma.n = usize::try_from(mem_sort_dedup_patch(opt, None, None, None, i32::try_from(ma.n).expect("n"), &mut ma.a)).expect("n");
            ma.m = ma.a.len();
        }
    }
    n
}

#[doc = "Original function: mem_pair:285"]
pub fn mem_pair(
    opt: &mem_opt_t,
    bns: &bntseq_t,
    _pac: &[u8],
    pes: &[mem_pestat_t; 4],
    _s: &[bseq1_t; 2],
    a: &mut [mem_alnreg_v; 2],
    id: i32,
    sub: &mut i32,
    n_sub: &mut i32,
    z: &mut [i32; 2],
    n_pri: &[i32; 2],
) -> i32 {
    let mut v = Vec::<pair64_t>::new();
    let mut u = Vec::<pair64_t>::new();
    let l_pac = bns.l_pac;

    for r in 0..2_usize {
        for i in 0..usize::try_from(n_pri[r]).expect("n_pri") {
            let e = &a[r].a[i];
            let mut x = if e.rb < l_pac { e.rb } else { (l_pac << 1) - 1 - e.rb };
            x = (i64::from(e.rid) << 32) | (x - bns.anns[usize::try_from(e.rid).expect("rid")].offset);
            let y = (u64::try_from(e.score).expect("score") << 32)
                | (u64::try_from(i).expect("i") << 2)
                | (u64::from((e.rb >= l_pac) as u8) << 1)
                | u64::try_from(r).expect("r");
            v.push(pair64_t { x: u64::try_from(x).expect("x"), y });
        }
    }

    v.sort_by(|lhs, rhs| lhs.x.cmp(&rhs.x).then(lhs.y.cmp(&rhs.y)));
    let mut y_last = [-1_i32; 4];
    for i in 0..v.len() {
        for r in 0..2_u64 {
            let dir = usize::try_from((r << 1) | ((v[i].y >> 1) & 1)).expect("dir");
            if pes[dir].failed != 0 {
                continue;
            }
            let which = usize::try_from((r << 1) | ((v[i].y & 1) ^ 1)).expect("which");
            if y_last[which] < 0 {
                continue;
            }
            for k in (0..=usize::try_from(y_last[which]).expect("y_last")).rev() {
                if (v[k].y & 3) != u64::try_from(which).expect("which") {
                    continue;
                }
                let dist = i64::try_from(v[i].x).expect("x") - i64::try_from(v[k].x).expect("x");
                if dist > i64::from(pes[dir].high) {
                    break;
                }
                if dist < i64::from(pes[dir].low) {
                    continue;
                }
                let ns = (dist as f64 - pes[dir].avg) / pes[dir].std;
                let mut q = ((v[i].y >> 32) + (v[k].y >> 32)) as f64;
                let erfc_term = erfc_approx(ns.abs() * std::f64::consts::FRAC_1_SQRT_2);
                q += 0.721 * (2.0_f64 * erfc_term).ln() * (opt.a as f64);
                // Match C++: (int)(... + 0.499) truncates toward zero, then clamp >= 0.
                let q_int = (q + 0.499) as i32;
                let q = u64::try_from(q_int.max(0)).expect("q");
                let pair_y = (u64::try_from(k).expect("k") << 32) | u64::try_from(i).expect("i");
                let pair_x = (q << 32) | (hash_64(pair_y ^ (u64::try_from(id).expect("id") << 8)) & 0xffff_ffff);
                u.push(pair64_t { x: pair_x, y: pair_y });
            }
        }
        y_last[usize::try_from(v[i].y & 3).expect("which")] = i32::try_from(i).expect("i");
    }

    if !u.is_empty() {
        let mut tmp = opt.a + opt.b;
        tmp = tmp.max(opt.o_del + opt.e_del);
        tmp = tmp.max(opt.o_ins + opt.e_ins);
        u.sort_by(|lhs, rhs| lhs.x.cmp(&rhs.x).then(lhs.y.cmp(&rhs.y)));
        let best = *u.last().expect("best pair");
        let i = usize::try_from(best.y >> 32).expect("i");
        let k = usize::try_from(best.y as u32).expect("k");
        z[usize::try_from(v[i].y & 1).expect("z")] = i32::try_from((v[i].y >> 2) & 0x3fff_ffff).expect("z idx");
        z[usize::try_from(v[k].y & 1).expect("z")] = i32::try_from((v[k].y >> 2) & 0x3fff_ffff).expect("z idx");
        let ret = i32::try_from(best.x >> 32).expect("ret");
        *sub = if u.len() > 1 { i32::try_from(u[u.len() - 2].x >> 32).expect("sub") } else { 0 };
        *n_sub = 0;
        for pair in u.iter().rev().skip(1) {
            if *sub - i32::try_from(pair.x >> 32).expect("pair score") <= tmp {
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
    let mut str_ = kstring_t::default();
    let mut h: [mem_aln_t; 2] = std::array::from_fn(|_| mem_aln_t::default());
    let mut g: [mem_aln_t; 2] = std::array::from_fn(|_| mem_aln_t::default());
    let mut aa: [Vec<mem_aln_t>; 2] = std::array::from_fn(|_| Vec::new());

    if (opt.flag & 0x20) == 0 {
        let b: [Vec<mem_alnreg_t>; 2] = std::array::from_fn(|i| {
            a[i].a
                .iter()
                .take(a[i].n)
                .filter(|reg| reg.score >= a[i].a[0].score - opt.pen_unpaired)
                .copied()
                .collect()
        });
        for i in 0..2_usize {
            for j in 0..b[i].len().min(usize::try_from(opt.max_matesw).expect("max_matesw")) {
                let mate_nt4 = nt4_cow(&s[1 - i]);
                let val = mem_matesw(opt, bns, pac, pes, &b[i][j], s[1 - i].l_seq, &mate_nt4, &mut a[1 - i]);
                n += val;
            }
        }
    }

    n_pri[0] = mem_mark_primary_se(opt, i32::try_from(a[0].n).expect("n"), &mut a[0].a, i64::try_from(id << 1).expect("id"));
    n_pri[1] = mem_mark_primary_se(
        opt,
        i32::try_from(a[1].n).expect("n"),
        &mut a[1].a,
        i64::try_from((id << 1) | 1).expect("id"),
    );
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
            a[0]
                .a
                .iter()
                .take(a[0].n.min(12))
                .map(|r| (r.score, r.sub, r.csub, r.qb, r.qe, r.rb, r.re, r.secondary, r.secondary_all))
                .collect::<Vec<_>>(),
            a[1]
                .a
                .iter()
                .take(a[1].n.min(12))
                .map(|r| (r.score, r.sub, r.csub, r.qb, r.qe, r.rb, r.re, r.secondary, r.secondary_all))
                .collect::<Vec<_>>(),
        );
    }

    if (opt.flag & MEM_F_NOPAIRING) == 0 {
        let s_view = [s[0].clone(), s[1].clone()];
        let o = if n_pri[0] > 0 && n_pri[1] > 0 {
            mem_pair(opt, bns, pac, pes, &s_view, a, i32::try_from(id).expect("id"), &mut subo, &mut n_sub, &mut z, &n_pri)
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
            let mut is_multi = [0_i32; 2];
            for i in 0..2_usize {
                let mut j = 1_usize;
                while j < usize::try_from(n_pri[i]).expect("n_pri") {
                    if a[i].a[j].secondary < 0 && a[i].a[j].score >= opt.T {
                        break;
                    }
                    j += 1;
                }
                is_multi[i] = if j < usize::try_from(n_pri[i]).expect("n_pri") { 1 } else { 0 };
            }
            if is_multi == [0, 0] {
                let score_un = a[0].a[0].score + a[1].a[0].score - opt.pen_unpaired;
                subo = subo.max(score_un);
                let mut q_pe = raw_mapq(o - subo, opt.a);
                if n_sub > 0 {
                    q_pe -= (4.343 * ((n_sub + 1) as f64).ln() + 0.499) as i32;
                }
                q_pe = q_pe.clamp(0, 60);
                q_pe = (q_pe as f64 * (1.0 - 0.5 * (a[0].a[0].frac_rep as f64 + a[1].a[0].frac_rep as f64)) + 0.499) as i32;

                let mut q_se = [0_i32; 2];
                if o > score_un {
                    for i in 0..2_usize {
                        let zi = usize::try_from(z[i]).expect("z");
                        let secondary = a[i].a[zi].secondary;
                        let secondary_score = if secondary >= 0 {
                            Some(a[i].a[usize::try_from(secondary).expect("secondary")].score)
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
                    q_se[0] = if q_se[0] > q_pe { q_se[0] } else if q_pe < q_se[0] + 40 { q_pe } else { q_se[0] + 40 };
                    q_se[1] = if q_se[1] > q_pe { q_se[1] } else if q_pe < q_se[1] + 40 { q_pe } else { q_se[1] + 40 };
                    extra_flag |= 2;
                    for i in 0..2_usize {
                        let c = &a[i].a[usize::try_from(z[i]).expect("z")];
                        q_se[i] = q_se[i].min(raw_mapq(c.score - c.csub, opt.a));
                    }
                } else {
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

                for i in 0..2_usize {
                    let zi = usize::try_from(z[i]).expect("z");
                    let k = a[i].a[zi].secondary_all;
                    if k >= 0 && k < n_pri[i] {
                        for j in 0..a[i].n {
                            if a[i].a[j].secondary_all == k || j == usize::try_from(k).expect("k") {
                                a[i].a[j].secondary_all = z[i];
                            }
                        }
                        a[i].a[zi].secondary_all = -1;
                    }
                }

                let xa = if (opt.flag & MEM_F_ALL) == 0 {
                    [
                        mem_gen_alt(opt, bns, pac, &a[0], s[0].l_seq, s[0].seq.as_deref().unwrap_or("")),
                        mem_gen_alt(opt, bns, pac, &a[1], s[1].l_seq, s[1].seq.as_deref().unwrap_or("")),
                    ]
                } else {
                    [vec![None; a[0].n], vec![None; a[1].n]]
                };

                for i in 0..2_usize {
                    let zi = usize::try_from(z[i]).expect("z");
                    h[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq.as_deref().unwrap_or(""), Some(&a[i].a[zi]));
                    h[i].mapq = u32::try_from(q_se[i]).expect("mapq");
                    h[i].flag |= (0x40 << i) | extra_flag;
                    h[i].XA = xa[i][zi].clone();
                    aa[i].push(h[i].clone());
                    if usize::try_from(n_pri[i]).expect("n_pri") < a[i].n {
                        let p = &a[i].a[usize::try_from(n_pri[i]).expect("n_pri")];
                        if p.score >= opt.T && p.secondary < 0 && p.is_alt != 0 {
                            g[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq.as_deref().unwrap_or(""), Some(p));
                            g[i].flag |= 0x800 | (0x40 << i) | extra_flag;
                            g[i].XA = xa[i][usize::try_from(n_pri[i]).expect("n_pri")].clone();
                            aa[i].push(g[i].clone());
                        }
                    }
                }

                for i in 0..aa[0].len() {
                    mem_aln2sam(opt, bns, &mut str_, &s[0], i32::try_from(aa[0].len()).expect("len"), &aa[0], i32::try_from(i).expect("i"), Some(&h[1]));
                }
                s[0].sam = Some(str_.as_str().to_owned());
                str_.l = 0;
                if !str_.s.is_empty() {
                    str_.s[0] = 0;
                }
                for i in 0..aa[1].len() {
                    mem_aln2sam(opt, bns, &mut str_, &s[1], i32::try_from(aa[1].len()).expect("len"), &aa[1], i32::try_from(i).expect("i"), Some(&h[0]));
                }
                s[1].sam = Some(str_.as_str().to_owned());
                assert_eq!(s[0].name, s[1].name, "paired reads have different names");
                return n;
            }
        }
    }

    for i in 0..2_usize {
        let mut which = -1_i32;
        if a[i].n > 0 {
            if a[i].a[0].score >= opt.T {
                which = 0;
            } else if usize::try_from(n_pri[i]).expect("n_pri") < a[i].n
                && a[i].a[usize::try_from(n_pri[i]).expect("n_pri")].score >= opt.T
            {
                which = n_pri[i];
            }
        }
        h[i] = if which >= 0 {
            mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq.as_deref().unwrap_or(""), Some(&a[i].a[usize::try_from(which).expect("which")]))
        } else {
            mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq.as_deref().unwrap_or(""), None)
        };
    }
    if (opt.flag & MEM_F_NOPAIRING) == 0 && h[0].rid == h[1].rid && h[0].rid >= 0 {
        let mut dist = 0_i64;
        let d = usize::try_from(mem_infer_dir(bns.l_pac, a[0].a[0].rb, a[1].a[0].rb, &mut dist)).expect("dir");
        if pes[d].failed == 0 && dist >= i64::from(pes[d].low) && dist <= i64::from(pes[d].high) {
            extra_flag |= 2;
        }
    }
    mem_reg2sam(opt, bns, pac, &mut s[0], &mut a[0], 0x41 | extra_flag, Some(&h[1]));
    mem_reg2sam(opt, bns, pac, &mut s[1], &mut a[1], 0x81 | extra_flag, Some(&h[0]));
    assert_eq!(s[0].name, s[1].name, "paired reads have different names");
    n
}

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
        for reg in a[i]
            .a
            .iter()
            .take(a[i].n)
            .filter(|reg| reg.score >= a[i].a[0].score - opt.pen_unpaired)
            .take(usize::try_from(opt.max_matesw).expect("max_matesw"))
        {
            let mate_nt4 = nt4_cow(&s[1 - i]);
            *pcnt = mem_matesw_batch_pre(
                opt, bns, pac, pes, reg, s[1 - i].l_seq, &mate_nt4, &a[1 - i], mmc, *pcnt, *gcnt,
                maxRefLen, maxQerLen, tid,
            );
            *gcnt += 4;
        }
    }
    1
}

#[doc = "Original function: revseq:604"]
pub fn revseq(l: i32, s: &mut [u8]) {
    let l = usize::try_from(l).expect("l");
    let mut i = 0_usize;
    while i < (l >> 1) {
        s.swap(i, l - 1 - i);
        i += 1;
    }
}

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
    let total = usize::try_from(pcnt).expect("pcnt");
    if total == 0 {
        return 1;
    }
    let pcnt8_usize = usize::try_from(pcnt8.max(0)).expect("pcnt8");
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
        pwsw.getScores8(
            pairs,
            &mmc.seqBufLeftRef[tid],
            &mmc.seqBufLeftQer[tid],
            aln,
            i32::try_from(pcnt8_usize).expect("pcnt8"),
            1,
            0,
        );
    }
    if pcnt16_usize > 0 {
        let pairs = &mut mmc.seqPairArrayLeft128[tid][pcnt8_usize..total];
        pwsw.getScores16(
            pairs,
            &mmc.seqBufLeftRef[tid],
            &mmc.seqBufLeftQer[tid],
            aln,
            i32::try_from(pcnt16_usize).expect("pcnt16"),
            1,
            0,
        );
    }

    // Phase 1 prep: reverse buffers in place, build the rev-pass pair list (u8 first, then i16)
    let mut phase1_pairs: Vec<SeqPair> = Vec::with_capacity(total);
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
            let ind = usize::try_from(sp.regid).expect("regid");
            let r = aln[ind];
            let xtra = sp.h0;
            if (xtra & KSW_XSTART) == 0
                || ((xtra & KSW_XSUBO) != 0 && r.score < (xtra & 0xffff))
            {
                continue;
            }
            sp.h0 = KSW_XSTOP | r.score;
            sp.len2 = r.qe + 1;
            let qs_start = usize::try_from(sp.idq).expect("idq");
            let rs_start = usize::try_from(sp.idr).expect("idr");
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
        pwsw.getScores16(
            pairs,
            &mmc.seqBufLeftRef[tid],
            &mmc.seqBufLeftQer[tid],
            aln,
            i32::try_from(pos16).expect("pos16"),
            1,
            1,
        );
    }
    if pos8 > 0 {
        let pairs = &mut phase1_pairs[..pos8];
        pwsw.getScores8(
            pairs,
            &mmc.seqBufLeftRef[tid],
            &mmc.seqBufLeftQer[tid],
            aln,
            i32::try_from(pos8).expect("pos8"),
            1,
            1,
        );
    }

    let _ = ksw_align2;
    1
}

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
    let gar: Vec<i32> = mmc.seqPairArrayAux[tid].iter().map(|sp| sp.id).collect();
    let mut n = 0_i32;
    let mut z = [0_i32; 2];
    let mut subo = 0_i32;
    let mut n_sub = 0_i32;
    let mut extra_flag = 1_i32;
    let mut n_pri = [0_i32; 2];
    let mut str_ = kstring_t::default();
    let mut h: [mem_aln_t; 2] = std::array::from_fn(|_| mem_aln_t::default());
    let mut g: [mem_aln_t; 2] = std::array::from_fn(|_| mem_aln_t::default());
    let mut aa: [Vec<mem_aln_t>; 2] = std::array::from_fn(|_| Vec::new());
    if (opt.flag & 0x20) == 0 {
        let b: [Vec<mem_alnreg_t>; 2] = std::array::from_fn(|i| {
            a[i].a
                .iter()
                .take(a[i].n)
                .filter(|reg| reg.score >= a[i].a[0].score - opt.pen_unpaired)
                .copied()
                .collect()
        });
        for reg in &mut a[0].a[..a[0].n] {
            reg.flg = 0;
        }
        for reg in &mut a[1].a[..a[1].n] {
            reg.flg = 0;
        }
        for i in 0..2_usize {
            for j in 0..b[i].len().min(usize::try_from(opt.max_matesw).expect("max_matesw")) {
                let mate_nt4 = nt4_cow(&s[1 - i]);
                let val = mem_matesw_batch_post(
                    opt, bns, pac, pes, &b[i][j], s[1 - i].l_seq, &mate_nt4, &mut a[1 - i], myaln,
                    *gcnt, &gar, mmc,
                );
                n += val;
                *gcnt += 4;
            }
        }
    }

    n_pri[0] = mem_mark_primary_se(opt, i32::try_from(a[0].n).expect("n"), &mut a[0].a, i64::try_from(id << 1).expect("id"));
    n_pri[1] = mem_mark_primary_se(
        opt,
        i32::try_from(a[1].n).expect("n"),
        &mut a[1].a,
        i64::try_from((id << 1) | 1).expect("id"),
    );
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
            a[0]
                .a
                .iter()
                .take(a[0].n.min(12))
                .map(|r| (r.score, r.sub, r.csub, r.qb, r.qe, r.rb, r.re, r.secondary, r.secondary_all))
                .collect::<Vec<_>>(),
            a[1]
                .a
                .iter()
                .take(a[1].n.min(12))
                .map(|r| (r.score, r.sub, r.csub, r.qb, r.qe, r.rb, r.re, r.secondary, r.secondary_all))
                .collect::<Vec<_>>(),
        );
    }

    if (opt.flag & MEM_F_NOPAIRING) == 0 {
        let s_view = [s[0].clone(), s[1].clone()];
        let o = if n_pri[0] > 0 && n_pri[1] > 0 {
            mem_pair(opt, bns, pac, pes, &s_view, a, i32::try_from(id).expect("id"), &mut subo, &mut n_sub, &mut z, &n_pri)
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
            let mut is_multi = [0_i32; 2];
            for i in 0..2_usize {
                let mut j = 1_usize;
                while j < usize::try_from(n_pri[i]).expect("n_pri") {
                    if a[i].a[j].secondary < 0 && a[i].a[j].score >= opt.T {
                        break;
                    }
                    j += 1;
                }
                is_multi[i] = if j < usize::try_from(n_pri[i]).expect("n_pri") { 1 } else { 0 };
            }
            if is_multi == [0, 0] {
                let score_un = a[0].a[0].score + a[1].a[0].score - opt.pen_unpaired;
                subo = subo.max(score_un);
                let mut q_pe = raw_mapq(o - subo, opt.a);
                if n_sub > 0 {
                    q_pe -= (4.343 * ((n_sub + 1) as f64).ln() + 0.499) as i32;
                }
                q_pe = q_pe.clamp(0, 60);
                q_pe = (q_pe as f64 * (1.0 - 0.5 * (a[0].a[0].frac_rep as f64 + a[1].a[0].frac_rep as f64)) + 0.499) as i32;

                let mut q_se = [0_i32; 2];
                if o > score_un {
                    for i in 0..2_usize {
                        let zi = usize::try_from(z[i]).expect("z");
                        let secondary = a[i].a[zi].secondary;
                        let secondary_score = if secondary >= 0 {
                            Some(a[i].a[usize::try_from(secondary).expect("secondary")].score)
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
                    q_se[0] = if q_se[0] > q_pe { q_se[0] } else if q_pe < q_se[0] + 40 { q_pe } else { q_se[0] + 40 };
                    q_se[1] = if q_se[1] > q_pe { q_se[1] } else if q_pe < q_se[1] + 40 { q_pe } else { q_se[1] + 40 };
                    extra_flag |= 2;
                    for i in 0..2_usize {
                        let c = &a[i].a[usize::try_from(z[i]).expect("z")];
                        q_se[i] = q_se[i].min(raw_mapq(c.score - c.csub, opt.a));
                    }
                } else {
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

                for i in 0..2_usize {
                    let zi = usize::try_from(z[i]).expect("z");
                    let k = a[i].a[zi].secondary_all;
                    if k >= 0 && k < n_pri[i] {
                        for j in 0..a[i].n {
                            if a[i].a[j].secondary_all == k || j == usize::try_from(k).expect("k") {
                                a[i].a[j].secondary_all = z[i];
                            }
                        }
                        a[i].a[zi].secondary_all = -1;
                    }
                }

                let xa = if (opt.flag & MEM_F_ALL) == 0 {
                    [
                        mem_gen_alt(opt, bns, pac, &a[0], s[0].l_seq, s[0].seq.as_deref().unwrap_or("")),
                        mem_gen_alt(opt, bns, pac, &a[1], s[1].l_seq, s[1].seq.as_deref().unwrap_or("")),
                    ]
                } else {
                    [vec![None; a[0].n], vec![None; a[1].n]]
                };

                for i in 0..2_usize {
                    let zi = usize::try_from(z[i]).expect("z");
                    h[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq.as_deref().unwrap_or(""), Some(&a[i].a[zi]));
                    h[i].mapq = u32::try_from(q_se[i]).expect("mapq");
                    h[i].flag |= (0x40 << i) | extra_flag;
                    h[i].XA = xa[i][zi].clone();
                    aa[i].push(h[i].clone());
                    if usize::try_from(n_pri[i]).expect("n_pri") < a[i].n {
                        let p = &a[i].a[usize::try_from(n_pri[i]).expect("n_pri")];
                        if p.score >= opt.T && p.secondary < 0 && p.is_alt != 0 {
                            g[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq.as_deref().unwrap_or(""), Some(p));
                            g[i].flag |= 0x800 | (0x40 << i) | extra_flag;
                            g[i].XA = xa[i][usize::try_from(n_pri[i]).expect("n_pri")].clone();
                            aa[i].push(g[i].clone());
                        }
                    }
                }

                for i in 0..aa[0].len() {
                    mem_aln2sam(opt, bns, &mut str_, &s[0], i32::try_from(aa[0].len()).expect("len"), &aa[0], i32::try_from(i).expect("i"), Some(&h[1]));
                }
                s[0].sam = Some(str_.as_str().to_owned());
                str_.l = 0;
                if !str_.s.is_empty() {
                    str_.s[0] = 0;
                }
                for i in 0..aa[1].len() {
                    mem_aln2sam(opt, bns, &mut str_, &s[1], i32::try_from(aa[1].len()).expect("len"), &aa[1], i32::try_from(i).expect("i"), Some(&h[0]));
                }
                s[1].sam = Some(str_.as_str().to_owned());
                assert_eq!(s[0].name, s[1].name, "paired reads have different names");
                return n;
            }
        }
    }

    for i in 0..2_usize {
        let mut which = -1_i32;
        if a[i].n > 0 {
            if a[i].a[0].score >= opt.T {
                which = 0;
            } else if usize::try_from(n_pri[i]).expect("n_pri") < a[i].n
                && a[i].a[usize::try_from(n_pri[i]).expect("n_pri")].score >= opt.T
            {
                which = n_pri[i];
            }
        }
        h[i] = if which >= 0 {
            mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq.as_deref().unwrap_or(""), Some(&a[i].a[usize::try_from(which).expect("which")]))
        } else {
            mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq.as_deref().unwrap_or(""), None)
        };
    }
    if (opt.flag & MEM_F_NOPAIRING) == 0 && h[0].rid == h[1].rid && h[0].rid >= 0 {
        let mut dist = 0_i64;
        let d = usize::try_from(mem_infer_dir(bns.l_pac, a[0].a[0].rb, a[1].a[0].rb, &mut dist)).expect("dir");
        if pes[d].failed == 0 && dist >= i64::from(pes[d].low) && dist <= i64::from(pes[d].high) {
            extra_flag |= 2;
        }
    }
    mem_reg2sam(opt, bns, pac, &mut s[0], &mut a[0], 0x41 | extra_flag, Some(&h[1]));
    mem_reg2sam(opt, bns, pac, &mut s[1], &mut a[1], 0x81 | extra_flag, Some(&h[0]));
    assert_eq!(s[0].name, s[1].name, "paired reads have different names");
    n
}

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
    for i in 0..ma.n {
        let mut dist = 0_i64;
        let r = usize::try_from(mem_infer_dir(l_pac, a.rb, ma.a[i].rb, &mut dist)).expect("dir");
        if dist >= i64::from(pes[r].low) && dist <= i64::from(pes[r].high) {
            skip[r] = 1;
        }
    }
    if skip.iter().copied().sum::<i32>() == 4 {
        while mmc.seqPairArrayAux[tid].len() < usize::try_from(gcnt + 4).expect("gcnt") {
            mmc.seqPairArrayAux[tid].push(SeqPair { id: -1, ..Default::default() });
        }
        for r in 0..4_usize {
            mmc.seqPairArrayAux[tid][usize::try_from(gcnt).expect("gcnt") + r].id = -1;
        }
        return pcnt;
    }

    let l_ms_usize = usize::try_from(l_ms).expect("l_ms");
    for r in 0..4_usize {
        if skip[r] != 0 {
            while mmc.seqPairArrayAux[tid].len() < usize::try_from(gcnt + r as i32 + 1).expect("gcnt") {
                mmc.seqPairArrayAux[tid].push(SeqPair { id: -1, ..Default::default() });
            }
            mmc.seqPairArrayAux[tid][usize::try_from(gcnt).expect("gcnt") + r].id = -1;
            continue;
        }
        let is_rev = ((r >> 1) != (r & 1)) as i32;
        let is_larger = ((r >> 1) == 0) as i32;
        let seq = if is_rev != 0 {
            let mut rev = vec![0_u8; l_ms_usize];
            for i in 0..l_ms_usize {
                rev[l_ms_usize - 1 - i] = if ms[i] < 4 { 3 - ms[i] } else { 4 };
            }
            rev
        } else {
            ms[..l_ms_usize].to_vec()
        };

        let mut rb;
        let mut re;
        if is_rev == 0 {
            rb = if is_larger != 0 { a.rb + i64::from(pes[r].low) } else { a.rb - i64::from(pes[r].high) };
            re = if is_larger != 0 { a.rb + i64::from(pes[r].high) } else { a.rb - i64::from(pes[r].low) };
            re += i64::from(l_ms);
        } else {
            rb = if is_larger != 0 { a.rb + i64::from(pes[r].low) } else { a.rb - i64::from(pes[r].high) };
            rb -= i64::from(l_ms);
            re = if is_larger != 0 { a.rb + i64::from(pes[r].high) } else { a.rb - i64::from(pes[r].low) };
        }
        if rb < 0 {
            rb = 0;
        }
        if re > (l_pac << 1) {
            re = l_pac << 1;
        }
        let mut rid = -1_i32;
        let ref_seq = if rb < re {
            let mid = (rb + re) >> 1;
            Some(bns_fetch_seq(bns, pac, &mut rb, mid, &mut re, &mut rid))
        } else {
            None
        };
        if a.rid == rid && re - rb >= i64::from(opt.min_seed_len) {
            let xtra = KSW_XSUBO | KSW_XSTART | if l_ms * opt.a < 250 { KSW_XBYTE } else { 0 } | (opt.min_seed_len * opt.a);
            let (ref_offset, qer_offset) = if pcnt != 0 {
                let prev = mmc.seqPairArrayLeft128[tid][usize::try_from(pcnt - 1).expect("pcnt")];
                (
                    prev.idr + prev.len1,
                    prev.idq + prev.len2,
                )
            } else {
                (0_i32, 0_i32)
            };
            let len1 = i32::try_from(re - rb).expect("len1");
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
            while mmc.seqPairArrayAux[tid].len() < usize::try_from(gcnt + r as i32 + 1).expect("gcnt") {
                mmc.seqPairArrayAux[tid].push(SeqPair { id: -1, ..Default::default() });
            }
            mmc.seqPairArrayAux[tid][usize::try_from(gcnt).expect("gcnt") + r].id = pcnt;
            let ref_start = usize::try_from(ref_offset).expect("ref_offset");
            let qer_start = usize::try_from(qer_offset).expect("qer_offset");
            let ref_end = ref_start + usize::try_from(len1).expect("len1");
            let qer_end = qer_start + usize::try_from(len2).expect("len2");
            if mmc.seqBufLeftRef[tid].len() < ref_end {
                mmc.seqBufLeftRef[tid].resize(ref_end, 0);
            }
            if mmc.seqBufLeftQer[tid].len() < qer_end {
                mmc.seqBufLeftQer[tid].resize(qer_end, 0);
            }
            mmc.seqBufLeftRef[tid][ref_start..ref_end]
                .copy_from_slice(ref_seq.as_ref().expect("reference seq"));
            mmc.seqBufLeftQer[tid][qer_start..qer_end].copy_from_slice(&seq);
            let sp_index = usize::try_from(pcnt).expect("pcnt");
            if mmc.seqPairArrayLeft128[tid].len() <= sp_index {
                mmc.seqPairArrayLeft128[tid].resize(sp_index + 1, SeqPair::default());
            }
            mmc.seqPairArrayLeft128[tid][sp_index] = sp;
            pcnt += 1;
        } else {
            while mmc.seqPairArrayAux[tid].len() < usize::try_from(gcnt + r as i32 + 1).expect("gcnt") {
                mmc.seqPairArrayAux[tid].push(SeqPair { id: -1, ..Default::default() });
            }
            mmc.seqPairArrayAux[tid][usize::try_from(gcnt).expect("gcnt") + r].id = -1;
        }
    }
    pcnt
}

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
    gar: &[i32],
    _mmc: &mut mem_cache,
) -> i32 {
    let l_pac = bns.l_pac;
    let mut skip = [0_i32; 4];
    for r in 0..4 {
        skip[r] = if pes[r].failed != 0 { 1 } else { 0 };
    }
    for i in 0..ma.n {
        let mut dist = 0_i64;
        let r = usize::try_from(mem_infer_dir(l_pac, a.rb, ma.a[i].rb, &mut dist)).expect("dir");
        if dist >= i64::from(pes[r].low) && dist <= i64::from(pes[r].high) {
            skip[r] = 1;
        }
    }
    if skip.iter().copied().sum::<i32>() == 4 {
        return 0;
    }

    let mut n = 0_i32;
    let l_ms_usize = usize::try_from(l_ms).expect("l_ms");
    for r in 0..4_usize {
        if skip[r] != 0 {
            continue;
        }
        let is_rev = ((r >> 1) != (r & 1)) as i32;
        let is_larger = ((r >> 1) == 0) as i32;
        let seq = if is_rev != 0 {
            let mut rev = vec![0_u8; l_ms_usize];
            for i in 0..l_ms_usize {
                rev[l_ms_usize - 1 - i] = if ms[i] < 4 { 3 - ms[i] } else { 4 };
            }
            rev
        } else {
            ms[..l_ms_usize].to_vec()
        };
        let mut rb;
        let mut re;
        if is_rev == 0 {
            rb = if is_larger != 0 { a.rb + i64::from(pes[r].low) } else { a.rb - i64::from(pes[r].high) };
            re = if is_larger != 0 { a.rb + i64::from(pes[r].high) } else { a.rb - i64::from(pes[r].low) };
            re += i64::from(l_ms);
        } else {
            rb = if is_larger != 0 { a.rb + i64::from(pes[r].low) } else { a.rb - i64::from(pes[r].high) };
            rb -= i64::from(l_ms);
            re = if is_larger != 0 { a.rb + i64::from(pes[r].high) } else { a.rb - i64::from(pes[r].low) };
        }
        if rb < 0 {
            rb = 0;
        }
        if re > (l_pac << 1) {
            re = l_pac << 1;
        }
        let mut rid = -1_i32;
        let ref_seq = if rb < re {
            let mid = (rb + re) >> 1;
            Some(bns_fetch_seq(bns, pac, &mut rb, mid, &mut re, &mut rid))
        } else {
            None
        };
        if a.rid == rid && re - rb >= i64::from(opt.min_seed_len) {
            let xtra = KSW_XSUBO | KSW_XSTART | if l_ms * opt.a < 250 { KSW_XBYTE } else { 0 } | (opt.min_seed_len * opt.a);
            let index = gar.get(usize::try_from(gcnt).expect("gcnt") + r).copied().unwrap_or(-1);
            let aln = if index < 0 {
                let mut query = seq.clone();
                let mut target = ref_seq.expect("reference seq");
                ksw_align2(
                    l_ms,
                    &mut query,
                    i32::try_from(re - rb).expect("ref span"),
                    &mut target,
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
                myaln[usize::try_from(index).expect("index")]
            };
            if aln.score >= opt.min_seed_len && aln.qb >= 0 {
                let b = mem_alnreg_t {
                    rid: a.rid,
                    is_alt: a.is_alt,
                    qb: if is_rev != 0 { l_ms - (aln.qe + 1) } else { aln.qb },
                    qe: if is_rev != 0 { l_ms - aln.qb } else { aln.qe + 1 },
                    rb: if is_rev != 0 { (l_pac << 1) - (rb + i64::from(aln.te) + 1) } else { rb + i64::from(aln.tb) },
                    re: if is_rev != 0 { (l_pac << 1) - (rb + i64::from(aln.tb)) } else { rb + i64::from(aln.te) + 1 },
                    score: aln.score,
                    csub: aln.score2,
                    secondary: -1,
                    seedcov: i32::try_from(((if is_rev != 0 { (l_pac << 1) - (rb + i64::from(aln.tb)) } else { rb + i64::from(aln.te) + 1 })
                        - (if is_rev != 0 { (l_pac << 1) - (rb + i64::from(aln.te) + 1) } else { rb + i64::from(aln.tb) }))
                        .min(i64::from(if is_rev != 0 { l_ms - aln.qb } else { aln.qe + 1 } - if is_rev != 0 { l_ms - (aln.qe + 1) } else { aln.qb }))
                        >> 1)
                    .expect("seedcov"),
                    ..Default::default()
                };
                ma.a.push(b);
                ma.n = ma.a.len();
                ma.m = ma.a.len();

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
            ma.n = usize::try_from(mem_sort_dedup_patch(opt, None, None, None, i32::try_from(ma.n).expect("n"), &mut ma.a)).expect("n");
            ma.m = ma.a.len();
        }
    }
    n
}

#[cfg(test)]
mod tests {
    use super::{
        cal_sub, mem_infer_dir, mem_matesw, mem_pair, mem_pestat, mem_sam_pe,
        mem_sam_pe_batch, mem_sam_pe_batch_post, mem_sam_pe_batch_pre, revseq,
    };
    use crate::generated::bwamem_cpp::mem_opt_init;
    use crate::generated::bntseq_h::{bntann1_t, bntseq_t};
    use crate::generated::bwa_h::bseq1_t;
    use crate::generated::bwamem_h::{mem_alnreg_t, mem_alnreg_v, mem_cache, mem_pestat_t};
    use crate::generated::ksw_h::kswr_t;
    use crate::generated::bandedswa_h::SeqPair;

    fn pack_seq(seq: &[u8]) -> Vec<u8> {
        let mut pac = vec![0_u8; (seq.len() + 3) / 4];
        for (i, &base) in seq.iter().enumerate() {
            let shift = (((!i64::try_from(i).expect("i")) & 3) << 1) as u8;
            pac[i >> 2] |= base << shift;
        }
        pac
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
                mem_alnreg_t { qb: 0, qe: 10, score: 20, ..Default::default() },
                mem_alnreg_t { qb: 2, qe: 9, score: 15, ..Default::default() },
            ],
        };
        assert_eq!(cal_sub(&opt, &r), 15);

        let r = mem_alnreg_v {
            n: 1,
            m: 1,
            a: vec![mem_alnreg_t { qb: 0, qe: 10, score: 20, ..Default::default() }],
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
            anns: vec![bntann1_t { offset: 0, len: 1000, name: "chr1".into(), ..Default::default() }],
            ..Default::default()
        };
        let pac = pack_seq(&vec![0_u8; 1000]);
        let pes = [
            mem_pestat_t { low: 10, high: 30, failed: 0, avg: 20.0, std: 1.0 },
            mem_pestat_t { failed: 1, ..Default::default() },
            mem_pestat_t { failed: 1, ..Default::default() },
            mem_pestat_t { failed: 1, ..Default::default() },
        ];
        let a = mem_alnreg_t { rid: 0, rb: 100, re: 110, ..Default::default() };
        let mut ma = mem_alnreg_v {
            n: 1,
            m: 1,
            a: vec![mem_alnreg_t { rid: 0, rb: 120, re: 130, score: 10, ..Default::default() }],
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
            mem_pestat_t { low: 2, high: 6, failed: 0, avg: 4.0, std: 1.0 },
            mem_pestat_t { failed: 1, ..Default::default() },
            mem_pestat_t { failed: 1, ..Default::default() },
            mem_pestat_t { failed: 1, ..Default::default() },
        ];
        let a = mem_alnreg_t { rid: 0, rb: 0, re: 4, is_alt: 0, ..Default::default() };
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
            anns: vec![bntann1_t { offset: 0, len: 1000, name: "chr1".into(), ..Default::default() }],
            ..Default::default()
        };
        let pes = [
            mem_pestat_t { failed: 1, ..Default::default() },
            mem_pestat_t { low: 15, high: 30, failed: 0, avg: 20.0, std: 3.0 },
            mem_pestat_t { failed: 1, ..Default::default() },
            mem_pestat_t { failed: 1, ..Default::default() },
        ];
        let s = [bseq1_t::default(), bseq1_t::default()];
        let mut regs = [
            mem_alnreg_v {
                n: 2,
                m: 2,
                a: vec![
                    mem_alnreg_t { rid: 0, rb: 1000 + 100, score: 40, ..Default::default() },
                    mem_alnreg_t { rid: 0, rb: 1000 + 150, score: 20, ..Default::default() },
                ],
            },
            mem_alnreg_v {
                n: 2,
                m: 2,
                a: vec![
                    mem_alnreg_t { rid: 0, rb: 880, score: 39, ..Default::default() },
                    mem_alnreg_t { rid: 0, rb: 840, score: 15, ..Default::default() },
                ],
            },
        ];
        let n_pri = [2, 2];
        let mut sub = -1;
        let mut n_sub = -1;
        let mut z = [-1, -1];
        let ret = mem_pair(&opt, &bns, &[], &pes, &s, &mut regs, 7, &mut sub, &mut n_sub, &mut z, &n_pri);
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
            anns: vec![bntann1_t { offset: 0, len: 1000, name: "chr1".into(), ..Default::default() }],
            ..Default::default()
        };
        let pes = [mem_pestat_t { failed: 1, ..Default::default() }; 4];
        let s = [bseq1_t::default(), bseq1_t::default()];
        let mut regs = [
            mem_alnreg_v { n: 1, m: 1, a: vec![mem_alnreg_t { rid: 0, rb: 100, score: 20, ..Default::default() }] },
            mem_alnreg_v { n: 1, m: 1, a: vec![mem_alnreg_t { rid: 0, rb: 200, score: 20, ..Default::default() }] },
        ];
        let mut sub = -1;
        let mut n_sub = -1;
        let mut z = [-1, -1];
        let ret = mem_pair(&opt, &bns, &[], &pes, &s, &mut regs, 1, &mut sub, &mut n_sub, &mut z, &[1, 1]);
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
            anns: vec![bntann1_t { offset: 0, len: 8, name: "chr1".into(), ..Default::default() }],
            ..Default::default()
        };
        let pac = pack_seq(&ref_seq);
        let pes = [mem_pestat_t::default(); 4];
        let mut s = [
            bseq1_t { l_seq: 4, name: Some("pair1".into()), seq: Some("ACGT".into()), qual: Some("IIII".into()), ..Default::default() },
            bseq1_t { l_seq: 4, name: Some("pair1".into()), seq: Some("ACGT".into()), qual: Some("JJJJ".into()), ..Default::default() },
        ];
        let mut a = [
            mem_alnreg_v { n: 1, m: 1, a: vec![mem_alnreg_t { rid: 0, rb: 0, re: 4, qb: 0, qe: 4, score: 4, truesc: 4, w: 100, ..Default::default() }] },
            mem_alnreg_v { n: 1, m: 1, a: vec![mem_alnreg_t { rid: 0, rb: 4, re: 8, qb: 0, qe: 4, score: 4, truesc: 4, w: 100, ..Default::default() }] },
        ];
        let rescued = mem_sam_pe(&opt, &bns, &pac, &pes, 1, &mut s, &mut a);
        assert_eq!(rescued, 0);
        assert!(s[0].sam.as_deref().expect("sam").starts_with("pair1\t65\t"), "{}", s[0].sam.as_deref().unwrap());
        assert!(s[1].sam.as_deref().expect("sam").starts_with("pair1\t129\t"), "{}", s[1].sam.as_deref().unwrap());
    }

    #[test]
    fn mem_sam_pe_writes_preferred_proper_pair_records() {
        let mut opt = (*mem_opt_init()).clone();
        opt.T = 1;
        let ref_seq = vec![0_u8, 1, 2, 3, 0, 1, 2, 3];
        let bns = bntseq_t {
            l_pac: 8,
            n_seqs: 1,
            anns: vec![bntann1_t { offset: 0, len: 8, name: "chr1".into(), ..Default::default() }],
            ..Default::default()
        };
        let pac = pack_seq(&ref_seq);
        let pes = [
            mem_pestat_t { failed: 1, ..Default::default() },
            mem_pestat_t { low: 1, high: 10, failed: 0, avg: 3.0, std: 1.0 },
            mem_pestat_t { failed: 1, ..Default::default() },
            mem_pestat_t { failed: 1, ..Default::default() },
        ];
        let mut s = [
            bseq1_t { l_seq: 4, name: Some("pair2".into()), seq: Some("ACGT".into()), qual: Some("IIII".into()), ..Default::default() },
            bseq1_t { l_seq: 4, name: Some("pair2".into()), seq: Some("ACGT".into()), qual: Some("JJJJ".into()), ..Default::default() },
        ];
        let mut a = [
            mem_alnreg_v { n: 1, m: 1, a: vec![mem_alnreg_t { rid: 0, rb: 8 + 4, re: 8 + 8, qb: 0, qe: 4, score: 40, truesc: 4, w: 100, ..Default::default() }] },
            mem_alnreg_v { n: 1, m: 1, a: vec![mem_alnreg_t { rid: 0, rb: 0, re: 4, qb: 0, qe: 4, score: 39, truesc: 4, w: 100, ..Default::default() }] },
        ];
        let rescued = mem_sam_pe(&opt, &bns, &pac, &pes, 2, &mut s, &mut a);
        assert_eq!(rescued, 0);
        let sam0 = s[0].sam.as_deref().expect("sam0");
        let sam1 = s[1].sam.as_deref().expect("sam1");
        assert!(sam0.starts_with("pair2\t99\t") || sam0.starts_with("pair2\t83\t"), "{sam0}");
        assert!(sam1.starts_with("pair2\t147\t") || sam1.starts_with("pair2\t163\t"), "{sam1}");
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
            mem_pestat_t { low: 3, high: 8, failed: 0, avg: 5.0, std: 1.0 },
            mem_pestat_t { failed: 1, ..Default::default() },
            mem_pestat_t { failed: 1, ..Default::default() },
            mem_pestat_t { failed: 1, ..Default::default() },
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
        mem_sam_pe_batch(&opt, &mut mmc, pcnt, 0, &mut aln, max_ref_len, max_qer_len, 0);
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
            mem_pestat_t { low: 1, high: 8, failed: 0, avg: 4.0, std: 1.0 },
            mem_pestat_t { failed: 1, ..Default::default() },
            mem_pestat_t { failed: 1, ..Default::default() },
            mem_pestat_t { failed: 1, ..Default::default() },
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
        mem_sam_pe_batch(&opt, &mut mmc, pcnt, 0, &mut aln, max_ref_len, max_qer_len, 0);
        assert!(aln.iter().take(usize::try_from(pcnt).expect("pcnt")).any(|x| x.score > 0));
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
