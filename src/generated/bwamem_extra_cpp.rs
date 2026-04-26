#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Generated scaffold for `bwa-mem2/src/bwamem_extra.cpp`.

use crate::generated::bntseq_h::bntseq_t;
use crate::generated::bwamem_cpp::mem_reg2aln;
use crate::generated::bwamem_h::{mem_alnreg_t, mem_alnreg_v, mem_opt_t};
use crate::generated::bwt_h::{bwt_t, bwtintv_v};
use crate::generated::kstring_h::{kputc, kputl, kputs, kputsn, kputw, kstring_t};

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
    pub matches: Box<bwtintv_v>,
    pub sub: Box<bwtintv_v>,
    pub tmpvec: [Box<bwtintv_v>; 2],
}

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

#[doc = "Original function: smem_itr_destroy:67"]
pub fn smem_itr_destroy(_itr: __smem_i) {}

#[doc = "Original function: smem_set_query:76"]
pub fn smem_set_query(itr: &mut __smem_i, len: i32, query: &[u8]) {
    itr.query = query[..usize::try_from(len).expect("len")].to_vec();
    itr.start = 0;
    itr.len = len;
}

#[doc = "Original function: smem_config:83"]
pub fn smem_config(itr: &mut __smem_i, min_intv: i32, max_len: i32, max_intv: u64) {
    itr.min_intv = min_intv;
    itr.max_len = max_len;
    itr.max_intv = max_intv;
}

#[doc = "Original function: smem_next:90"]
pub fn smem_next(itr: &mut __smem_i) -> Option<&bwtintv_v> {
    itr.tmpvec[0].n = 0;
    itr.tmpvec[1].n = 0;
    itr.matches.n = 0;
    itr.sub.n = 0;
    if itr.start >= itr.len || itr.start < 0 {
        return None;
    }
    while itr.start < itr.len && itr.query[usize::try_from(itr.start).expect("start")] > 3 {
        itr.start += 1;
    }
    if itr.start == itr.len {
        return None;
    }
    todo!("smem_next requires bwt_smem1a, which is not translated yet")
}

#[doc = "Original function: mem_align1:107"]
pub fn mem_align1(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
    _arg2: crate::support::Opaque,
    _arg3: crate::support::Opaque,
    _arg4: crate::support::Opaque,
    _arg5: crate::support::Opaque,
) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("mem_align1")
}

#[doc = "Original function: get_pri_idx:122"]
pub fn get_pri_idx(XA_drop_ratio: f64, a: &[mem_alnreg_t], i: i32) -> i32 {
    let i = usize::try_from(i).expect("i");
    let k = a[i].secondary_all;
    if k >= 0
        && (a[i].score as f64) >= (a[usize::try_from(k).expect("k")].score as f64) * XA_drop_ratio
    {
        k
    } else {
        -1
    }
}

#[doc = "Original function: mem_gen_alt:130"]
pub fn mem_gen_alt(
    opt: &mem_opt_t,
    bns: &bntseq_t,
    pac: &[u8],
    a: &mem_alnreg_v,
    l_query: i32,
    query: &str,
) -> Vec<Option<String>> {
    let mut cnt = vec![0_i32; a.n];
    let mut has_alt = vec![false; a.n];
    let mut tot = 0_i32;
    for i in 0..a.n {
        let r = get_pri_idx(opt.XA_drop_ratio as f64, &a.a, i32::try_from(i).expect("i"));
        if r >= 0 {
            let r = usize::try_from(r).expect("r");
            cnt[r] += 1;
            tot += 1;
            if a.a[i].is_alt != 0 {
                has_alt[r] = true;
            }
        }
    }
    if tot == 0 {
        return vec![None; a.n];
    }

    let mut aln = vec![kstring_t::default(); a.n];
    let mut str_ = kstring_t::default();
    for i in 0..a.n {
        let r = get_pri_idx(opt.XA_drop_ratio as f64, &a.a, i32::try_from(i).expect("i"));
        if r < 0 {
            continue;
        }
        let r = usize::try_from(r).expect("r");
        if cnt[r] > opt.max_XA_hits_alt || (!has_alt[r] && cnt[r] > opt.max_XA_hits) {
            continue;
        }
        let t = mem_reg2aln(opt, bns, pac, l_query, query, Some(&a.a[i]));
        str_.l = 0;
        kputs(
            &bns.anns[usize::try_from(t.rid).expect("rid")].name,
            &mut str_,
        );
        kputc(i32::from(b','), &mut str_);
        kputc(
            i32::from(if t.is_rev != 0 { b'-' } else { b'+' }),
            &mut str_,
        );
        kputl(t.pos + 1, &mut str_);
        kputc(i32::from(b','), &mut str_);
        for &cigar in t
            .cigar
            .iter()
            .take(usize::try_from(t.n_cigar).expect("n_cigar"))
        {
            kputw(i32::try_from(cigar >> 4).expect("len"), &mut str_);
            kputc(
                i32::from(b"MIDSHN"[usize::try_from(cigar & 0xf).expect("op")]),
                &mut str_,
            );
        }
        kputc(i32::from(b','), &mut str_);
        kputw(i32::try_from(t.NM).expect("NM"), &mut str_);
        kputc(i32::from(b';'), &mut str_);
        kputsn(
            str_.as_bytes(),
            i32::try_from(str_.l).expect("str len"),
            &mut aln[r],
        );
    }

    aln.into_iter()
        .map(|s| {
            if s.l == 0 {
                None
            } else {
                let mut bytes = s.s;
                bytes.truncate(s.l);
                Some(String::from_utf8(bytes).expect("XA tag contains invalid UTF-8"))
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::{
        get_pri_idx, mem_gen_alt, smem_config, smem_itr_destroy, smem_itr_init, smem_next,
        smem_set_query,
    };
    use crate::generated::bntseq_h::{bntann1_t, bntseq_t};
    use crate::generated::bwamem_cpp::mem_opt_init;
    use crate::generated::bwamem_h::{mem_alnreg_t, mem_alnreg_v};
    use crate::generated::bwt_h::bwt_t;

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
}
