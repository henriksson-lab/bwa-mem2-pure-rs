#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Generated scaffold for `bwa-mem2/src/bwa.cpp`.

use std::sync::atomic::{AtomicI32, Ordering};
use std::sync::{LazyLock, Mutex};

use crate::generated::bntseq_cpp::bns_get_seq_into;
use crate::generated::bntseq_h::bntseq_t;
use crate::generated::bwa_h::bseq1_t;
use crate::generated::kseq_h::{kseq_read, kseq_t};
use crate::generated::kstring_h::{kputc, kputw, kstring_t};
use crate::generated::ksw_cpp::ksw_global2;
use crate::generated::utils_cpp::{err_fputc, err_fputs, ErrFile};

// Thread-local rseq scratch for bwa_gen_cigar2 (allocated once per ~100K alignments otherwise).
thread_local! {
    static BWA_GEN_CIGAR_RSEQ: std::cell::RefCell<Vec<u8>> =
        const { std::cell::RefCell::new(Vec::new()) };
}

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct bwa_cigar_t {
    pub cigar: Vec<u32>,
    pub md: Vec<u8>,
}

pub static BWA_VERBOSE: AtomicI32 = AtomicI32::new(3);
pub static BWA_RG_ID: LazyLock<Mutex<String>> = LazyLock::new(|| Mutex::new(String::new()));
pub static BWA_PG: LazyLock<Mutex<Option<String>>> = LazyLock::new(|| Mutex::new(None));

#[doc = "Original function: trim_readno:62"]
pub fn trim_readno(s: &mut kstring_t) {
    if s.l > 2 && s.s[s.l - 2] == b'/' && s.s[s.l - 1].is_ascii_digit() {
        s.l -= 2;
        if s.s.len() <= s.l {
            s.s.resize(s.l + 1, 0);
            s.m = s.m.max(s.s.len());
        }
        s.s[s.l] = 0;
    }
}

#[doc = "Original function: kseq2bseq1:68"]
pub fn kseq2bseq1(ks: &kseq_t, s: &mut bseq1_t) {
    s.name = Some(ks.name.as_str().to_string().into_boxed_str());
    s.comment = if ks.comment.l != 0 {
        Some(ks.comment.as_str().to_string().into_boxed_str())
    } else {
        None
    };
    s.seq = Some(ks.seq.as_str().to_string().into_boxed_str());
    s.qual = if ks.qual.l != 0 {
        Some(ks.qual.as_str().to_string().into_boxed_str())
    } else {
        None
    };
    s.l_seq = i32::try_from(ks.seq.l).expect("sequence length");
}

#[doc = "Original function: bseq_read:78"]
pub fn bseq_read(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
    _arg2: crate::support::Opaque,
    _arg3: crate::support::Opaque,
    _arg4: crate::support::Opaque,
    _arg5: crate::support::Opaque,
    _arg6: crate::support::Opaque,
) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("bseq_read")
}

#[doc = "Original function: bseq_read_orig:170"]
pub fn bseq_read_orig(
    chunk_size: i64,
    n_: &mut i32,
    ks1: &mut kseq_t,
    ks2: Option<&mut kseq_t>,
    s: &mut i64,
) -> Vec<bseq1_t> {
    let mut size = 0_i64;
    let mut seqs = Vec::new();
    let mut ks2 = ks2;
    while kseq_read(ks1) >= 0 {
        if let Some(ks2_) = ks2.as_deref_mut() {
            if kseq_read(ks2_) < 0 {
                eprintln!("[W::bseq_read_orig] the 2nd file has fewer sequences.");
                break;
            }
        }

        trim_readno(&mut ks1.name);
        let mut seq = bseq1_t::default();
        kseq2bseq1(ks1, &mut seq);
        seq.id = i32::try_from(seqs.len()).expect("id");
        size += i64::from(seq.l_seq);
        seqs.push(seq);

        if let Some(ks2_) = ks2.as_deref_mut() {
            trim_readno(&mut ks2_.name);
            let mut seq = bseq1_t::default();
            kseq2bseq1(ks2_, &mut seq);
            seq.id = i32::try_from(seqs.len()).expect("id");
            size += i64::from(seq.l_seq);
            seqs.push(seq);
        }

        if size >= chunk_size && (seqs.len() & 1) == 0 {
            break;
        }
    }
    if size == 0 {
        if let Some(ks2_) = ks2.as_deref_mut() {
            if kseq_read(ks2_) >= 0 {
                eprintln!("[W::bseq_read_orig] the 1st file has fewer sequences.");
            }
        }
    }
    *n_ = i32::try_from(seqs.len()).expect("n_seqs");
    *s = size;
    seqs
}

#[doc = "Original function: bseq_read_one_fasta_file:218"]
pub fn bseq_read_one_fasta_file(
    chunk_size: i64,
    n_: &mut i32,
    text: &str,
    s: &mut i64,
) -> Vec<bseq1_t> {
    let mut ks = kseq_t::from_text(text);
    bseq_read_orig(chunk_size, n_, &mut ks, None, s)
}

#[doc = "Original function: bseq_classify:226"]
pub fn bseq_classify(n: i32, seqs: &[bseq1_t], m: &mut [i32; 2], sep: &mut [Vec<bseq1_t>; 2]) {
    let limit = usize::try_from(n.max(0)).expect("n");
    let seqs = &seqs[..limit.min(seqs.len())];
    let mut single = Vec::new();
    let mut paired = Vec::new();

    if seqs.is_empty() {
        m[0] = 0;
        m[1] = 0;
        sep[0].clear();
        sep[1].clear();
        return;
    }

    let mut has_last = true;
    for i in 1..seqs.len() {
        if has_last {
            if seqs[i].name == seqs[i - 1].name {
                paired.push(seqs[i - 1].clone());
                paired.push(seqs[i].clone());
                has_last = false;
            } else {
                single.push(seqs[i - 1].clone());
            }
        } else {
            has_last = true;
        }
    }
    if has_last {
        single.push(seqs[seqs.len() - 1].clone());
    }

    m[0] = i32::try_from(single.len()).expect("single len");
    m[1] = i32::try_from(paired.len()).expect("paired len");
    sep[0] = single;
    sep[1] = paired;
}

#[doc = "Original function: bwa_fill_scmat:248"]
#[inline]
pub fn bwa_fill_scmat(a: i32, b: i32, mat: &mut [i8; 25]) {
    let mut k = 0_usize;
    for i in 0..4 {
        for j in 0..4 {
            mat[k] = if i == j { a as i8 } else { -(b as i8) };
            k += 1;
        }
        mat[k] = -1;
        k += 1;
    }
    while k < 25 {
        mat[k] = -1;
        k += 1;
    }
}

#[doc = "Original function: bwa_gen_cigar2:260"]
pub fn bwa_gen_cigar2(
    mat: &[i8; 25],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    w_: i32,
    l_pac: i64,
    pac: &[u8],
    l_query: i32,
    query: &mut [u8],
    rb: i64,
    re: i64,
    score: &mut i32,
    n_cigar: Option<&mut i32>,
    NM: Option<&mut i32>,
) -> Option<bwa_cigar_t> {
    let mut n_cigar_slot = n_cigar;
    let mut nm_slot = NM;
    if let Some(n_cigar) = n_cigar_slot.as_deref_mut() {
        *n_cigar = 0;
    }
    if let Some(nm) = nm_slot.as_deref_mut() {
        *nm = -1;
    }
    if l_query <= 0 || rb >= re || (rb < l_pac && re > l_pac) {
        return None;
    }

    let mut rlen = 0_i64;
    // Take the thread-local rseq buffer (capacity reused across ~100K calls).
    let mut rseq = BWA_GEN_CIGAR_RSEQ.with(|cell| std::mem::take(&mut *cell.borrow_mut()));
    if !bns_get_seq_into(l_pac, pac, rb, re, &mut rlen, &mut rseq) || re - rb != rlen {
        BWA_GEN_CIGAR_RSEQ.with(|cell| *cell.borrow_mut() = rseq);
        return None;
    }

    let l_query_usize = usize::try_from(l_query).expect("l_query");
    if rb >= l_pac {
        query[..l_query_usize].reverse();
        rseq.reverse();
    }

    let mut out_cigar: Option<Vec<u32>> = None;
    if i64::from(l_query) == re - rb && w_ == 0 {
        if n_cigar_slot.is_some() {
            out_cigar = Some(vec![(u32::try_from(l_query).expect("l_query") << 4)]);
            if let Some(n_cigar) = n_cigar_slot.as_deref_mut() {
                *n_cigar = 1;
            }
        }
        *score = 0;
        for i in 0..l_query_usize {
            *score += i32::from(mat[usize::from(rseq[i]) * 5 + usize::from(query[i])]);
        }
    } else {
        let half_query = (l_query + 1) >> 1;
        let max_ins =
            (((half_query * i32::from(mat[0]) - o_ins) as f64) / e_ins as f64 + 1.0) as i32;
        let max_del =
            (((half_query * i32::from(mat[0]) - o_del) as f64) / e_del as f64 + 1.0) as i32;
        let mut max_gap = max_ins.max(max_del);
        max_gap = max_gap.max(1);
        let mut w =
            (max_gap + i32::try_from((rlen - i64::from(l_query)).abs()).expect("band delta") + 1)
                >> 1;
        w = w.min(w_);
        let min_w = i32::try_from((rlen - i64::from(l_query)).abs()).expect("min_w") + 3;
        w = w.max(min_w);

        let mut n_cigar_tmp = 0_i32;
        let mut cigar_tmp = Vec::new();
        *score = ksw_global2(
            l_query,
            &query[..l_query_usize],
            i32::try_from(rlen).expect("rlen"),
            &rseq,
            5,
            &mat[..],
            o_del,
            e_del,
            o_ins,
            e_ins,
            w,
            if n_cigar_slot.is_some() {
                Some(&mut n_cigar_tmp)
            } else {
                None
            },
            if n_cigar_slot.is_some() {
                Some(&mut cigar_tmp)
            } else {
                None
            },
        );
        if n_cigar_slot.is_some() {
            if let Some(n_cigar) = n_cigar_slot.as_deref_mut() {
                *n_cigar = n_cigar_tmp;
            }
            out_cigar = Some(cigar_tmp);
        }
    }

    let mut md = Vec::new();
    if nm_slot.is_some() && n_cigar_slot.is_some() {
        let mut str = kstring_t::default();
        let int2base = if rb < l_pac { b"ACGTN" } else { b"TGCAN" };
        let cigar = out_cigar.as_ref().expect("cigar available with n_cigar");
        let mut x = 0_usize;
        let mut y = 0_usize;
        let mut u = 0_i32;
        let mut n_mm = 0_i32;
        let mut n_gap = 0_i32;
        for (k, cigar_op) in cigar.iter().enumerate() {
            let op = cigar_op & 0xf;
            let len = usize::try_from(cigar_op >> 4).expect("cigar len");
            if op == 0 {
                for i in 0..len {
                    if query[x + i] != rseq[y + i] {
                        kputw(u, &mut str);
                        kputc(i32::from(int2base[usize::from(rseq[y + i])]), &mut str);
                        n_mm += 1;
                        u = 0;
                    } else {
                        u += 1;
                    }
                }
                x += len;
                y += len;
            } else if op == 2 {
                if k > 0 && k + 1 < cigar.len() {
                    kputw(u, &mut str);
                    kputc(i32::from(b'^'), &mut str);
                    for i in 0..len {
                        kputc(i32::from(int2base[usize::from(rseq[y + i])]), &mut str);
                    }
                    u = 0;
                    n_gap += i32::try_from(len).expect("gap len");
                }
                y += len;
            } else if op == 1 {
                x += len;
                n_gap += i32::try_from(len).expect("gap len");
            }
        }
        kputw(u, &mut str);
        if let Some(nm) = nm_slot.as_deref_mut() {
            *nm = n_mm + n_gap;
        }
        // Move the kstring buffer instead of cloning bytes.
        str.s.truncate(str.l);
        md = str.s;
    }

    if rb >= l_pac {
        query[..l_query_usize].reverse();
    }

    // Restore the rseq buffer to the thread-local for the next call's reuse.
    BWA_GEN_CIGAR_RSEQ.with(|cell| *cell.borrow_mut() = rseq);

    out_cigar.map(|cigar| bwa_cigar_t { cigar, md })
}

#[doc = "Original function: bwa_gen_cigar:349"]
pub fn bwa_gen_cigar(
    mat: &[i8; 25],
    q: i32,
    r: i32,
    w_: i32,
    l_pac: i64,
    pac: &[u8],
    l_query: i32,
    query: &mut [u8],
    rb: i64,
    re: i64,
    score: &mut i32,
    n_cigar: Option<&mut i32>,
    NM: Option<&mut i32>,
) -> Option<bwa_cigar_t> {
    bwa_gen_cigar2(
        mat, q, r, q, r, w_, l_pac, pac, l_query, query, rb, re, score, n_cigar, NM,
    )
}

#[doc = "Original function: bwa_idx_infer_prefix:358"]
pub fn bwa_idx_infer_prefix(_arg0: crate::support::Opaque) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("bwa_idx_infer_prefix")
}

#[doc = "Original function: bwa_idx_load_bwt:384"]
pub fn bwa_idx_load_bwt(_arg0: crate::support::Opaque) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("bwa_idx_load_bwt")
}

#[doc = "Original function: bwa_idx_load_from_disk:402"]
pub fn bwa_idx_load_from_disk(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("bwa_idx_load_from_disk")
}

#[doc = "Original function: bwa_idx_load:433"]
pub fn bwa_idx_load(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("bwa_idx_load")
}

#[doc = "Original function: bwa_idx_destroy:438"]
pub fn bwa_idx_destroy(_arg0: crate::support::Opaque) {
    crate::support::stub::<()>("bwa_idx_destroy")
}

#[doc = "Original function: bwa_mem2idx:452"]
pub fn bwa_mem2idx(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
    _arg2: crate::support::Opaque,
) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("bwa_mem2idx")
}

#[doc = "Original function: bwa_idx2mem:477"]
pub fn bwa_idx2mem(_arg0: crate::support::Opaque) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("bwa_idx2mem")
}

#[doc = "Original function: bwa_print_sam_hdr:523"]
pub fn bwa_print_sam_hdr(bns: &bntseq_t, hdr_line: Option<&str>, fp: &mut ErrFile) {
    let mut n_sq = 0_i32;
    if let Some(hdr_line) = hdr_line {
        let mut p = 0_usize;
        while let Some(found) = hdr_line[p..].find("@SQ\t") {
            let idx = p + found;
            if idx == 0 || hdr_line.as_bytes()[idx - 1] == b'\n' {
                n_sq += 1;
            }
            p = idx + 4;
        }
    }
    if n_sq == 0 {
        for ann in &bns.anns {
            let line = format!("@SQ\tSN:{}\tLN:{}", ann.name, ann.len);
            err_fputs(&line, fp);
            if ann.is_alt != 0 {
                err_fputs("\tAH:*\n", fp);
            } else {
                err_fputc(i32::from(b'\n'), fp);
            }
        }
    } else if n_sq != bns.n_seqs && BWA_VERBOSE.load(Ordering::Relaxed) >= 2 {
        eprintln!(
            "[W::bwa_print_sam_hdr] {} @SQ lines provided with -H; {} sequences in the index. Continue anyway.",
            n_sq, bns.n_seqs
        );
    }

    if let Some(hdr_line) = hdr_line {
        err_fputs(hdr_line, fp);
        err_fputs("\n", fp);
    }
    if let Some(pg) = BWA_PG.lock().expect("bwa_pg lock").as_ref() {
        err_fputs(pg, fp);
    }
}

#[doc = "Original function: bwa_escape:567"]
pub fn bwa_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut chars = s.chars();
    while let Some(ch) = chars.next() {
        if ch == '\\' {
            match chars.next() {
                Some('t') => out.push('\t'),
                Some('n') => out.push('\n'),
                Some('r') => out.push('\r'),
                Some('\\') => out.push('\\'),
                Some(other) => {
                    out.push('\\');
                    out.push(other);
                }
                None => out.push('\\'),
            }
        } else {
            out.push(ch);
        }
    }
    out
}

#[doc = "Original function: bwa_set_rg:583"]
pub fn bwa_set_rg(s: &str) -> Option<String> {
    *BWA_RG_ID.lock().expect("bwa_rg_id lock") = String::new();
    if !s.starts_with("@RG") {
        if BWA_VERBOSE.load(Ordering::Relaxed) >= 1 {
            eprintln!("[E::bwa_set_rg] the read group line is not started with @RG");
        }
        return None;
    }
    let rg_line = bwa_escape(s);
    let Some(id_idx) = rg_line.find("\tID:") else {
        if BWA_VERBOSE.load(Ordering::Relaxed) >= 1 {
            eprintln!("[E::bwa_set_rg] no ID at the read group line");
        }
        return None;
    };
    let id_start = id_idx + 4;
    let id_end = rg_line[id_start..]
        .find(['\t', '\n'])
        .map(|x| id_start + x)
        .unwrap_or(rg_line.len());
    let id = &rg_line[id_start..id_end];
    if id.len() + 1 > 256 {
        if BWA_VERBOSE.load(Ordering::Relaxed) >= 1 {
            eprintln!("[E::bwa_set_rg] @RG:ID is longer than 255 characters");
        }
        return None;
    }
    *BWA_RG_ID.lock().expect("bwa_rg_id lock") = id.to_string();
    Some(rg_line)
}

#[doc = "Original function: bwa_insert_header:612"]
pub fn bwa_insert_header(s: Option<&str>, hdr: Option<String>) -> Option<String> {
    let s = s?;
    if !s.starts_with('@') {
        return hdr;
    }
    // C++ bwa.cpp:623 escapes only the newly-appended portion (`bwa_escape(hdr + len)`),
    // preserving any earlier escape pass on the existing header.
    let escaped_new = bwa_escape(s);
    Some(if let Some(mut hdr) = hdr {
        hdr.push('\n');
        hdr.push_str(&escaped_new);
        hdr
    } else {
        escaped_new
    })
}

#[cfg(test)]
mod tests {
    use super::{
        bseq_classify, bseq_read_one_fasta_file, bseq_read_orig, bwa_escape, bwa_fill_scmat,
        bwa_gen_cigar, bwa_gen_cigar2, bwa_insert_header, bwa_print_sam_hdr, bwa_set_rg,
        kseq2bseq1, trim_readno, BWA_PG, BWA_RG_ID, BWA_VERBOSE,
    };
    use crate::generated::bntseq_h::{bntann1_t, bntseq_t};
    use crate::generated::bwa_h::bseq1_t;
    use crate::generated::kseq_h::{kseq_read, kseq_t};
    use crate::generated::kstring_h::kstring_t;
    use crate::generated::utils_cpp::{err_fclose, err_xopen_core};
    use std::fs;
    use std::sync::atomic::Ordering;

    fn pack_seq(seq: &[u8]) -> Vec<u8> {
        let mut pac = vec![0_u8; (seq.len() + 3) / 4];
        for (i, &base) in seq.iter().enumerate() {
            let shift = (((!i64::try_from(i).expect("i")) & 3) << 1) as u8;
            pac[i >> 2] |= base << shift;
        }
        pac
    }

    fn temp_path(label: &str) -> std::path::PathBuf {
        let mut p = std::env::temp_dir();
        let stamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos();
        p.push(format!("bwa_mem2_rs_bwa_{label}_{stamp}.sam"));
        p
    }

    #[test]
    fn trim_readno_removes_numeric_slash_suffix_only() {
        let mut s = kstring_t {
            l: 6,
            m: 8,
            s: b"read/1\0x".to_vec(),
        };
        trim_readno(&mut s);
        assert_eq!(s.as_str(), "read");

        let mut unchanged = kstring_t {
            l: 6,
            m: 7,
            s: b"read/x\0".to_vec(),
        };
        trim_readno(&mut unchanged);
        assert_eq!(unchanged.as_str(), "read/x");
    }

    #[test]
    fn bseq_classify_splits_singletons_and_pairs_by_adjacent_name() {
        let seqs = vec![
            bseq1_t {
                name: Some("solo0".into()),
                id: 0,
                ..Default::default()
            },
            bseq1_t {
                name: Some("pair1".into()),
                id: 1,
                ..Default::default()
            },
            bseq1_t {
                name: Some("pair1".into()),
                id: 2,
                ..Default::default()
            },
            bseq1_t {
                name: Some("solo3".into()),
                id: 3,
                ..Default::default()
            },
            bseq1_t {
                name: Some("pair4".into()),
                id: 4,
                ..Default::default()
            },
            bseq1_t {
                name: Some("pair4".into()),
                id: 5,
                ..Default::default()
            },
            bseq1_t {
                name: Some("solo6".into()),
                id: 6,
                ..Default::default()
            },
        ];
        let mut counts = [0, 0];
        let mut sep = [Vec::new(), Vec::new()];
        bseq_classify(
            i32::try_from(seqs.len()).expect("len"),
            &seqs,
            &mut counts,
            &mut sep,
        );

        assert_eq!(counts, [3, 4]);
        assert_eq!(
            sep[0]
                .iter()
                .map(|s| s.name.as_deref().unwrap())
                .collect::<Vec<_>>(),
            vec!["solo0", "solo3", "solo6"]
        );
        assert_eq!(
            sep[1].iter().map(|s| s.id).collect::<Vec<_>>(),
            vec![1, 2, 4, 5]
        );
    }

    #[test]
    fn kseq2bseq1_copies_current_record_fields() {
        let mut ks = kseq_t::from_text("@read0 desc\nACGT\n+\nIIII\n");
        assert_eq!(kseq_read(&mut ks), 4);
        let mut out = bseq1_t::default();
        kseq2bseq1(&ks, &mut out);
        assert_eq!(out.name.as_deref(), Some("read0"));
        assert_eq!(out.comment.as_deref(), Some("desc"));
        assert_eq!(out.seq.as_deref(), Some("ACGT"));
        assert_eq!(out.qual.as_deref(), Some("IIII"));
        assert_eq!(out.l_seq, 4);
    }

    #[test]
    fn bseq_read_orig_preserves_pair_chunking_and_trimmed_names() {
        let mut ks1 = kseq_t::from_text("@pairA/1\nACGT\n+\nIIII\n@pairB/1\nTT\n+\nJJ\n");
        let mut ks2 = kseq_t::from_text("@pairA/2\nTGCA\n+\nHHHH\n@pairB/2\nAA\n+\nKK\n");
        let mut n = 0;
        let mut size = 0;
        let seqs = bseq_read_orig(5, &mut n, &mut ks1, Some(&mut ks2), &mut size);
        assert_eq!(n, 2);
        assert_eq!(size, 8);
        assert_eq!(seqs.len(), 2);
        assert_eq!(seqs[0].name.as_deref(), Some("pairA"));
        assert_eq!(seqs[1].name.as_deref(), Some("pairA"));
        assert_eq!(seqs[0].id, 0);
        assert_eq!(seqs[1].id, 1);
    }

    #[test]
    fn bseq_read_one_fasta_file_reads_all_records_from_text() {
        let mut n = 0;
        let mut size = 0;
        let seqs = bseq_read_one_fasta_file(100, &mut n, ">r0\nACGT\n>r1 note\nTT\n", &mut size);
        assert_eq!(n, 2);
        assert_eq!(size, 6);
        assert_eq!(seqs[0].name.as_deref(), Some("r0"));
        assert_eq!(seqs[0].comment, None);
        assert_eq!(seqs[1].name.as_deref(), Some("r1"));
        assert_eq!(seqs[1].comment.as_deref(), Some("note"));
        assert_eq!(seqs[1].seq.as_deref(), Some("TT"));
    }

    #[test]
    fn bwa_gen_cigar2_handles_ungapped_forward_alignment_and_md() {
        let mut mat = [0_i8; 25];
        bwa_fill_scmat(1, 4, &mut mat);
        let pac = pack_seq(&[0, 1, 2, 3]);
        let mut query = vec![0_u8, 1, 2, 3];
        let mut score = 0;
        let mut n_cigar = 0;
        let mut nm = -1;
        let result = bwa_gen_cigar2(
            &mat,
            6,
            1,
            6,
            1,
            0,
            4,
            &pac,
            4,
            &mut query,
            0,
            4,
            &mut score,
            Some(&mut n_cigar),
            Some(&mut nm),
        )
        .expect("cigar");
        assert_eq!(score, 4);
        assert_eq!(n_cigar, 1);
        assert_eq!(nm, 0);
        assert_eq!(result.cigar, vec![4 << 4]);
        assert_eq!(result.md, b"4");
        assert_eq!(query, vec![0, 1, 2, 3]);
    }

    #[test]
    fn bwa_gen_cigar2_handles_mismatch_md_and_nm() {
        let mut mat = [0_i8; 25];
        bwa_fill_scmat(1, 4, &mut mat);
        let pac = pack_seq(&[0, 1, 2, 3]);
        let mut query = vec![0_u8, 1, 3, 3];
        let mut score = 0;
        let mut n_cigar = 0;
        let mut nm = -1;
        let result = bwa_gen_cigar2(
            &mat,
            6,
            1,
            6,
            1,
            0,
            4,
            &pac,
            4,
            &mut query,
            0,
            4,
            &mut score,
            Some(&mut n_cigar),
            Some(&mut nm),
        )
        .expect("cigar");
        assert_eq!(score, -1);
        assert_eq!(n_cigar, 1);
        assert_eq!(nm, 1);
        assert_eq!(result.md, b"2G1");
    }

    #[test]
    fn bwa_gen_cigar_wrapper_matches_symmetric_gap_penalties() {
        let mut mat = [0_i8; 25];
        bwa_fill_scmat(2, 3, &mut mat);
        let pac = pack_seq(&[0, 1, 2]);
        let mut query = vec![0_u8, 1, 2];
        let mut score = 0;
        let mut n_cigar = 0;
        let result = bwa_gen_cigar(
            &mat,
            5,
            1,
            0,
            3,
            &pac,
            3,
            &mut query,
            0,
            3,
            &mut score,
            Some(&mut n_cigar),
            None,
        )
        .expect("cigar");
        assert_eq!(score, 6);
        assert_eq!(n_cigar, 1);
        assert_eq!(result.cigar, vec![3 << 4]);
    }

    #[test]
    fn bwa_escape_unescapes_control_sequences() {
        assert_eq!(
            bwa_escape(r"@RG\tID:foo\nSM:bar\\baz"),
            "@RG\tID:foo\nSM:bar\\baz"
        );
    }

    #[test]
    fn bwa_set_rg_extracts_id_and_returns_escaped_line() {
        BWA_VERBOSE.store(3, Ordering::Relaxed);
        let rg = bwa_set_rg(r"@RG\tID:foo\tSM:bar").expect("rg");
        assert_eq!(rg, "@RG\tID:foo\tSM:bar");
        assert_eq!(&*BWA_RG_ID.lock().expect("lock"), "foo");
    }

    #[test]
    fn bwa_insert_header_appends_and_escapes() {
        let hdr = bwa_insert_header(Some(r"@HD\tVN:1.6"), None).expect("hdr");
        let hdr = bwa_insert_header(Some(r"@RG\tID:foo"), Some(hdr)).expect("hdr2");
        assert_eq!(hdr, "@HD\tVN:1.6\n@RG\tID:foo");
    }

    #[test]
    fn bwa_print_sam_hdr_emits_sq_lines_header_and_pg() {
        let path = temp_path("hdr");
        let mut fp = err_xopen_core("test", path.to_str().expect("utf8"), "w");
        *BWA_PG.lock().expect("lock") = Some("@PG\tID:bwa-mem2\n".to_string());
        BWA_VERBOSE.store(3, Ordering::Relaxed);
        let bns = bntseq_t {
            n_seqs: 2,
            anns: vec![
                bntann1_t {
                    name: "chr1".into(),
                    len: 10,
                    is_alt: 0,
                    ..Default::default()
                },
                bntann1_t {
                    name: "alt1".into(),
                    len: 7,
                    is_alt: 1,
                    ..Default::default()
                },
            ],
            ..Default::default()
        };
        bwa_print_sam_hdr(&bns, Some("@HD\tVN:1.6"), &mut fp);
        assert_eq!(err_fclose(fp), 0);
        let text = fs::read_to_string(&path).expect("read");
        assert!(text.contains("@SQ\tSN:chr1\tLN:10\n"));
        assert!(text.contains("@SQ\tSN:alt1\tLN:7\tAH:*\n"));
        assert!(text.contains("@HD\tVN:1.6\n"));
        assert!(text.contains("@PG\tID:bwa-mem2\n"));
        let _ = fs::remove_file(path);
        *BWA_PG.lock().expect("lock") = None;
    }
}
