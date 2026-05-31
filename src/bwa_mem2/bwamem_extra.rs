#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Port of `bwa-mem2/src/bwamem_extra.cpp`.

use crate::bwa_mem2::bntseq::bntseq_t;
use crate::bwa_mem2::bwamem::{mem_alnreg_t, mem_alnreg_v, mem_opt_t, mem_reg2aln};
use crate::bwa_mem2::bwt::{bwt_t, bwtintv_v};
use crate::bwa_mem2::kstring::{kputc, kputl, kputs, kputw, kstring_t};

// --- bwamem_extra.cpp ---

// SMEM iterator interface.

#[doc = "Original struct: __smem_i (bwa-mem2/src/bwamem_extra.cpp)"]
#[derive(Debug, Default, Clone)]
pub struct __smem_i {
    pub bwt: bwt_t,
    pub query: Vec<u8>,
    pub start: i32,
    pub len: i32,
    pub min_intv: i32,
    pub max_len: i32,
    pub max_intv: u64,
    /// matches; to be returned by `smem_next()`
    pub matches: Box<bwtintv_v>,
    /// sub-matches inside the longest match; temporary
    pub sub: Box<bwtintv_v>,
    /// temporary arrays
    pub tmpvec: [Box<bwtintv_v>; 2],
}

/// Initialise an SMEM iterator bound to a BWT, with default `min_intv=1`,
/// `max_len=INT_MAX` and `max_intv=0`. The query is set separately via
/// [`smem_set_query`].
#[doc = "Original function: smem_itr_init:52"]
pub fn smem_itr_init(bwt: &bwt_t) -> __smem_i {
    __smem_i {
        bwt: *bwt,
        query: Vec::new(),
        start: 0,
        len: 0,
        min_intv: 1,
        max_len: i32::MAX,
        max_intv: 0,
        matches: Box::new(bwtintv_v::default()),
        sub: Box::new(bwtintv_v::default()),
        tmpvec: [
            Box::new(bwtintv_v::default()),
            Box::new(bwtintv_v::default()),
        ],
    }
}

/// Destroy an SMEM iterator. In upstream C++ this frees the per-iterator
/// scratch `bwtintv_v` buffers; the Rust port relies on `Drop` and only needs
/// to consume the value.
#[doc = "Original function: smem_itr_destroy:67"]
pub fn smem_itr_destroy(_itr: __smem_i) {}

/// Bind a 2-bit-encoded query of length `len` to the iterator and rewind the
/// scan cursor to position 0.
#[doc = "Original function: smem_set_query:76"]
pub fn smem_set_query(itr: &mut __smem_i, len: i32, query: &[u8]) {
    itr.query = query[..usize::try_from(len).expect("len")].to_vec();
    itr.start = 0;
    itr.len = len;
}

/// Override the iterator's SMEM filtering thresholds (minimum SA interval
/// size, maximum match length, maximum SA interval).
#[doc = "Original function: smem_config:83"]
pub fn smem_config(itr: &mut __smem_i, min_intv: i32, max_len: i32, max_intv: u64) {
    itr.min_intv = min_intv;
    itr.max_len = max_len;
    itr.max_intv = max_intv;
}

#[doc = "Original function: smem_next:90"]
pub(crate) fn smem_next(itr: &mut __smem_i) -> Option<&bwtintv_v> {
    itr.tmpvec[0].n = 0;
    itr.tmpvec[1].n = 0;
    itr.matches.n = 0;
    itr.sub.n = 0;
    if itr.start >= itr.len || itr.start < 0 {
        return None;
    }
    // skip ambiguous bases
    while itr.start < itr.len && itr.query[usize::try_from(itr.start).expect("start")] > 3 {
        itr.start += 1;
    }
    if itr.start == itr.len {
        return None;
    }
    // search for SMEM via bwt_smem1a (upstream call disabled inside #if 0).
    panic!("smem_next is inside #if 0 in upstream bwa-mem2/src/bwamem_extra.cpp")
}

// Extra functions.

/// The difference from `mem_align1_core()` is that this routine: 1) calls
/// `mem_mark_primary_se()`; 2) does not modify the input sequence.
#[doc = "Original function: mem_align1:107"]
pub(crate) fn mem_align1(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
    _arg2: crate::support::Opaque,
    _arg3: crate::support::Opaque,
    _arg4: crate::support::Opaque,
    _arg5: crate::support::Opaque,
) -> crate::support::Opaque {
    panic!("mem_align1 is inside #if 0 in upstream bwa-mem2/src/bwamem_extra.cpp")
}

/// Return the index of the primary alignment that hit `i` should be reported
/// under in the `XA` tag, or `-1` if it shouldn't be reported. Requires the
/// alignment's score to exceed the primary's by at least `XA_drop_ratio`.
#[doc = "Original function: get_pri_idx:122"]
pub fn get_pri_idx(XA_drop_ratio: f64, a: &[mem_alnreg_t], i: i32) -> i32 {
    let i = i as usize;
    let k = a[i].secondary_all;
    if k >= 0 && (a[i].score as f64) >= (a[k as usize].score as f64) * XA_drop_ratio {
        k
    } else {
        -1
    }
}

thread_local! {
    static MEM_GEN_ALT_SCRATCH: std::cell::RefCell<(Vec<i32>, Vec<bool>, kstring_t, Vec<i32>)> =
        std::cell::RefCell::new((Vec::new(), Vec::new(), kstring_t::default(), Vec::new()));
}

/// Build the per-primary `XA:Z:` SAM-tag strings for each alignment in `a`.
///
/// Returns a vector with one entry per element of `a`: `Some(xa)` for primary
/// hits whose secondaries qualify, `None` otherwise. Each `XA` entry is a
/// `;`-terminated list of `chr,pos,cigar,NM` tuples summarising the
/// secondaries that hit the same primary.
///
/// ONLY works after `mem_mark_primary_se()`.
///
/// (Upstream comment: "Okay, returning strings is bad, but this has happened
/// a lot elsewhere. If I have time, I need serious code cleanup.")
#[doc = "Original function: mem_gen_alt:130"]
pub fn mem_gen_alt(
    opt: &mem_opt_t,
    bns: &bntseq_t,
    pac: &[u8],
    a: &mem_alnreg_v,
    l_query: i32,
    query: &str,
) -> Vec<Option<String>> {
    let (cnt_owned, has_alt_owned, str_owned, pri_idx_owned) =
        MEM_GEN_ALT_SCRATCH.with(|c| std::mem::take(&mut *c.borrow_mut()));
    let mut cnt = cnt_owned;
    let mut has_alt = has_alt_owned;
    let mut pri_idx = pri_idx_owned;
    cnt.clear();
    cnt.resize(a.n, 0_i32);
    has_alt.clear();
    has_alt.resize(a.n, false);
    pri_idx.clear();
    pri_idx.reserve(a.n);
    let mut tot = 0_i32;
    let drop_ratio = opt.XA_drop_ratio as f64;
    for i in 0..a.n {
        let r = get_pri_idx(drop_ratio, &a.a, i as i32);
        pri_idx.push(r);
        if r >= 0 {
            let r = r as usize;
            cnt[r] += 1;
            tot += 1;
            if a.a[i].is_alt != 0 {
                has_alt[r] = true;
            }
        }
    }
    if tot == 0 {
        MEM_GEN_ALT_SCRATCH.with(|c| {
            let mut scratch = c.borrow_mut();
            scratch.0 = cnt;
            scratch.1 = has_alt;
            scratch.2 = str_owned;
            scratch.3 = pri_idx;
        });
        return vec![None; a.n];
    }

    // Build XA strings directly into per-primary kstring_t entries — avoids the prior
    // str_ scratch + kputsn copy round-trip. The unused str_ stays in scratch for tot>0
    // re-entry parity.
    let mut aln = vec![kstring_t::default(); a.n];
    for i in 0..a.n {
        let r = pri_idx[i];
        if r < 0 {
            continue;
        }
        let r = r as usize;
        if cnt[r] > opt.max_XA_hits_alt || (!has_alt[r] && cnt[r] > opt.max_XA_hits) {
            continue;
        }
        let t = mem_reg2aln(opt, bns, pac, l_query, query, Some(&a.a[i]));
        let dst = &mut aln[r];
        kputs(&bns.anns[t.rid as usize].name, dst);
        kputc(i32::from(b','), dst);
        kputc(i32::from(if t.is_rev != 0 { b'-' } else { b'+' }), dst);
        kputl(t.pos + 1, dst);
        kputc(i32::from(b','), dst);
        const MIDSHN_LUT: [u8; 6] = *b"MIDSHN";
        for &cigar in t.cigar.iter().take(t.n_cigar as usize) {
            kputw((cigar >> 4) as i32, dst);
            let op_idx = (cigar & 0xf) as usize;
            kputc(
                i32::from(unsafe { *MIDSHN_LUT.get_unchecked(op_idx.min(5)) }),
                dst,
            );
        }
        kputc(i32::from(b','), dst);
        kputw(t.NM as i32, dst);
        kputc(i32::from(b';'), dst);
    }

    // SAFETY: every byte written above came from kput* helpers (ASCII) or from chromosome
    // names from bns.anns (which are stored as Rust String, hence valid UTF-8). Skip the
    // O(n) validation pass.
    let result: Vec<Option<String>> = aln
        .into_iter()
        .map(|s| {
            if s.l == 0 {
                None
            } else {
                let mut bytes = s.s;
                bytes.truncate(s.l);
                Some(unsafe { String::from_utf8_unchecked(bytes) })
            }
        })
        .collect();
    MEM_GEN_ALT_SCRATCH.with(|c| {
        let mut scratch = c.borrow_mut();
        scratch.0 = cnt;
        scratch.1 = has_alt;
        scratch.2 = str_owned;
        scratch.3 = pri_idx;
    });
    result
}

#[cfg(test)]
mod tests {
    use super::{
        get_pri_idx, mem_align1, mem_gen_alt, smem_config, smem_itr_destroy, smem_itr_init,
        smem_next, smem_set_query,
    };
    use crate::bwa_mem2::bntseq::{bntann1_t, bntseq_t};
    use crate::bwa_mem2::bwamem::mem_opt_init;
    use crate::bwa_mem2::bwamem::{mem_alnreg_t, mem_alnreg_v};
    use crate::bwa_mem2::bwt::bwt_t;

    fn pack_seq(seq: &[u8]) -> Vec<u8> {
        let mut pac = vec![0_u8; (seq.len() + 3) / 4];
        for (i, &base) in seq.iter().enumerate() {
            let shift = (((!i64::try_from(i).expect("i")) & 3) << 1) as u8;
            pac[i >> 2] |= base << shift;
        }
        pac
    }

    #[test]
    fn get_pri_idx_requires_secondary_all_and_drop_ratio_pass() {
        let regs = vec![
            mem_alnreg_t {
                score: 10,
                ..Default::default()
            },
            mem_alnreg_t {
                score: 8,
                secondary_all: 0,
                ..Default::default()
            },
            mem_alnreg_t {
                score: 4,
                secondary_all: 0,
                ..Default::default()
            },
        ];
        assert_eq!(get_pri_idx(0.7, &regs, 1), 0);
        assert_eq!(get_pri_idx(0.7, &regs, 2), -1);
    }

    #[test]
    fn mem_gen_alt_builds_xa_strings_for_primary_hits() {
        let mut opt = (*mem_opt_init()).clone();
        opt.XA_drop_ratio = 0.5;
        opt.max_XA_hits = 5;
        opt.max_XA_hits_alt = 5;
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
        let pac = pack_seq(&[0, 1, 2, 3, 0, 1, 2, 3]);
        let regs = mem_alnreg_v {
            n: 2,
            m: 2,
            a: vec![
                mem_alnreg_t {
                    rb: 0,
                    re: 4,
                    qb: 0,
                    qe: 4,
                    rid: 0,
                    score: 10,
                    truesc: 4,
                    w: 100,
                    ..Default::default()
                },
                mem_alnreg_t {
                    rb: 4,
                    re: 8,
                    qb: 0,
                    qe: 4,
                    rid: 0,
                    score: 8,
                    truesc: 4,
                    secondary_all: 0,
                    w: 100,
                    ..Default::default()
                },
            ],
        };

        let xa = mem_gen_alt(&opt, &bns, &pac, &regs, 4, "ACGT");
        assert_eq!(xa.len(), 2);
        let primary = xa[0].as_deref().expect("primary xa");
        assert!(primary.contains("chr1,+5,4M,0;"), "{primary}");
        assert!(xa[1].is_none());
    }

    #[test]
    fn smem_iterator_init_and_configure_tracks_state() {
        let bwt = bwt_t::default();
        let mut itr = smem_itr_init(&bwt);
        assert_eq!(itr.min_intv, 1);
        assert_eq!(itr.max_len, i32::MAX);
        assert_eq!(itr.max_intv, 0);
        smem_set_query(&mut itr, 5, &[0, 1, 4, 2, 3]);
        assert_eq!(itr.query, vec![0, 1, 4, 2, 3]);
        assert_eq!(itr.start, 0);
        assert_eq!(itr.len, 5);
        smem_config(&mut itr, 7, 99, 1234);
        assert_eq!(itr.min_intv, 7);
        assert_eq!(itr.max_len, 99);
        assert_eq!(itr.max_intv, 1234);
        smem_itr_destroy(itr);
    }

    #[test]
    fn smem_next_returns_none_once_query_is_exhausted() {
        let bwt = bwt_t::default();
        let mut itr = smem_itr_init(&bwt);
        smem_set_query(&mut itr, 3, &[4, 4, 4]);
        assert!(smem_next(&mut itr).is_none());
        assert_eq!(itr.start, 3);
    }

    #[test]
    #[should_panic(expected = "#if 0")]
    fn smem_next_panics_for_upstream_disabled_active_search_path() {
        let bwt = bwt_t::default();
        let mut itr = smem_itr_init(&bwt);
        smem_set_query(&mut itr, 1, &[0]);
        let _ = smem_next(&mut itr);
    }

    #[test]
    #[should_panic(expected = "#if 0")]
    fn mem_align1_panics_because_upstream_excludes_it() {
        let _ = mem_align1(0, 0, 0, 0, 0, 0);
    }
}
