#![allow(dead_code, non_snake_case, non_camel_case_types, non_upper_case_globals)]

//! Generated scaffold for `bwa-mem2/src/bntseq.cpp`.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

use flate2::read::MultiGzDecoder;

use crate::generated::bntseq_h::{bntamb1_t, bntann1_t};
use crate::generated::bntseq_h::{bns_depos, bntseq_t};

fn parse_i64(text: &str, ctx: &str) -> i64 {
    text.parse().unwrap_or_else(|_| panic!("parse error reading {ctx}"))
}

fn parse_i32(text: &str, ctx: &str) -> i32 {
    text.parse().unwrap_or_else(|_| panic!("parse error reading {ctx}"))
}

fn parse_u32(text: &str, ctx: &str) -> u32 {
    text.parse().unwrap_or_else(|_| panic!("parse error reading {ctx}"))
}

fn next_line(lines: &mut impl Iterator<Item = std::io::Result<String>>, ctx: &str) -> String {
    lines.next()
        .unwrap_or_else(|| panic!("Error reading {ctx} : Unexpected end of file"))
        .unwrap_or_else(|e| panic!("Error reading {ctx} : {e}"))
}

fn split_header<'a>(line: &'a str, ctx: &str) -> (&'a str, &'a str, &'a str) {
    let first_end = line.find(char::is_whitespace).unwrap_or_else(|| panic!("parse error reading {ctx}"));
    let first = &line[..first_end];
    let rest = line[first_end..].trim_start();
    let second_end = rest.find(char::is_whitespace).unwrap_or(rest.len());
    let second = &rest[..second_end];
    if second.is_empty() {
        panic!("parse error reading {ctx}");
    }
    let tail = &rest[second_end..];
    (first, second, tail)
}

#[inline(always)]
fn get_pac(pac: &[u8], l: i64) -> u8 {
    let idx = (l >> 2) as usize;
    let shift = ((!l) & 3) << 1;
    // Caller (bns_get_seq) bounds-checks beg/end against l_pac before iterating, so idx is in
    // range. Skip the bounds check to keep this inlined; bns_get_seq is in the chain-extension
    // hot path (showed 1.2% in profile).
    unsafe { (*pac.get_unchecked(idx) >> shift) & 3 }
}

fn set_pac(pac: &mut [u8], l: i64, c: u8) {
    let idx = usize::try_from(l >> 2).expect("pac index overflow");
    let shift = ((!l) & 3) << 1;
    pac[idx] |= c << shift;
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SequenceRecord {
    name: String,
    comment: Option<String>,
    seq: Vec<u8>,
}

fn nt4(c: u8) -> u8 {
    match c {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' | b'U' | b'u' => 3,
        b'-' => 5,
        _ => 4,
    }
}

fn parse_fasta<R: BufRead>(reader: R) -> Vec<SequenceRecord> {
    let mut records = Vec::new();
    let mut current: Option<SequenceRecord> = None;
    for line in reader.lines() {
        let line = line.expect("read fasta");
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(record) = current.take() {
                records.push(record);
            }
            let mut parts = rest.splitn(2, char::is_whitespace);
            let name = parts.next().unwrap_or("").to_string();
            let comment = parts
                .next()
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .map(ToOwned::to_owned);
            current = Some(SequenceRecord { name, comment, seq: Vec::new() });
        } else if !line.is_empty() {
            current
                .as_mut()
                .expect("FASTA sequence before header")
                .seq
                .extend_from_slice(line.trim().as_bytes());
        }
    }
    if let Some(record) = current {
        records.push(record);
    }
    records
}

#[doc = "Original function: bns_dump:73"]
pub fn bns_dump(bns: &bntseq_t, prefix: &str) {
    let ann_name = format!("{prefix}.ann");
    let amb_name = format!("{prefix}.amb");

    {
        let fp = File::create(&ann_name).unwrap_or_else(|e| panic!("fail to open file '{ann_name}' : {e}"));
        let mut w = BufWriter::new(fp);
        writeln!(w, "{} {} {}", bns.l_pac, bns.n_seqs, bns.seed).expect("write ann header");
        for ann in &bns.anns {
            write!(w, "{} {}", ann.gi, ann.name).expect("write ann name");
            if !ann.anno.is_empty() {
                writeln!(w, " {}", ann.anno).expect("write ann comment");
            } else {
                writeln!(w).expect("write ann newline");
            }
            writeln!(w, "{} {} {}", ann.offset, ann.len, ann.n_ambs).expect("write ann body");
        }
        w.flush().expect("flush ann");
    }

    {
        let fp = File::create(&amb_name).unwrap_or_else(|e| panic!("fail to open file '{amb_name}' : {e}"));
        let mut w = BufWriter::new(fp);
        writeln!(w, "{} {} {}", bns.l_pac, bns.n_seqs, bns.n_holes).expect("write amb header");
        for amb in &bns.ambs {
            writeln!(w, "{} {} {}", amb.offset, amb.len, amb.amb as char).expect("write amb body");
        }
        w.flush().expect("flush amb");
    }
}

#[doc = "Original function: bns_restore_core:106"]
pub fn bns_restore_core(ann_filename: &str, amb_filename: &str, pac_filename: &str) -> bntseq_t {
    let ann_file = File::open(ann_filename)
        .unwrap_or_else(|e| panic!("fail to open file '{ann_filename}' : {e}"));
    let mut ann_lines = BufReader::new(ann_file).lines();
    let header = next_line(&mut ann_lines, ann_filename);
    let mut fields = header.split_whitespace();
    let l_pac = parse_i64(fields.next().unwrap_or(""), ann_filename);
    let n_seqs = parse_i32(fields.next().unwrap_or(""), ann_filename);
    let seed = parse_u32(fields.next().unwrap_or(""), ann_filename);

    let mut bns = bntseq_t {
        l_pac,
        n_seqs,
        seed,
        anns: Vec::with_capacity(n_seqs as usize),
        ..Default::default()
    };

    for _ in 0..n_seqs {
        let seq_header = next_line(&mut ann_lines, ann_filename);
        let (gi, name, tail) = split_header(&seq_header, ann_filename);
        let anno = if tail.len() > 1 && tail != " (null)" {
            tail[1..].to_string()
        } else {
            String::new()
        };
        let metrics = next_line(&mut ann_lines, ann_filename);
        let mut metrics_fields = metrics.split_whitespace();
        let offset = parse_i64(metrics_fields.next().unwrap_or(""), ann_filename);
        let len = parse_i32(metrics_fields.next().unwrap_or(""), ann_filename);
        let n_ambs = parse_i32(metrics_fields.next().unwrap_or(""), ann_filename);
        bns.anns.push(bntann1_t {
            offset,
            len,
            n_ambs,
            gi: parse_u32(gi, ann_filename),
            is_alt: 0,
            name: name.to_string(),
            anno,
        });
    }

    let amb_file = File::open(amb_filename)
        .unwrap_or_else(|e| panic!("fail to open file '{amb_filename}' : {e}"));
    let mut amb_lines = BufReader::new(amb_file).lines();
    let amb_header = next_line(&mut amb_lines, amb_filename);
    let mut amb_fields = amb_header.split_whitespace();
    let amb_l_pac = parse_i64(amb_fields.next().unwrap_or(""), amb_filename);
    let amb_n_seqs = parse_i32(amb_fields.next().unwrap_or(""), amb_filename);
    let n_holes = parse_i32(amb_fields.next().unwrap_or(""), amb_filename);
    assert!(
        amb_l_pac == bns.l_pac && amb_n_seqs == bns.n_seqs,
        "inconsistent .ann and .amb files."
    );
    bns.n_holes = n_holes;
    bns.ambs = Vec::with_capacity(n_holes as usize);
    for _ in 0..n_holes {
        let line = next_line(&mut amb_lines, amb_filename);
        let mut fields = line.split_whitespace();
        let offset = parse_i64(fields.next().unwrap_or(""), amb_filename);
        let len = parse_i32(fields.next().unwrap_or(""), amb_filename);
        let amb = fields
            .next()
            .unwrap_or_else(|| panic!("parse error reading {amb_filename}"))
            .as_bytes()[0];
        bns.ambs.push(bntamb1_t { offset, len, amb });
    }

    bns.fp_pac = Some(
        File::open(pac_filename)
            .unwrap_or_else(|e| panic!("fail to open file '{pac_filename}' : {e}"))
    );
    bns
}

#[doc = "Original function: bns_restore:188"]
pub fn bns_restore(prefix: &str) -> bntseq_t {
    let ann_filename = format!("{prefix}.ann");
    let amb_filename = format!("{prefix}.amb");
    let pac_filename = format!("{prefix}.pac");
    let alt_filename = format!("{prefix}.alt");
    let mut bns = bns_restore_core(&ann_filename, &amb_filename, &pac_filename);
    if Path::new(&alt_filename).is_file() {
        let file = File::open(&alt_filename)
            .unwrap_or_else(|e| panic!("fail to open file '{alt_filename}' : {e}"));
        let mut names = HashMap::new();
        for (i, ann) in bns.anns.iter().enumerate() {
            names.insert(ann.name.clone(), i);
        }
        for line in BufReader::new(file).lines() {
            let line = line.unwrap_or_else(|e| panic!("Error reading {alt_filename} : {e}"));
            let token = line
                .split(['\t', '\n', '\r', ' '])
                .find(|part| !part.is_empty())
                .unwrap_or("");
            if token.starts_with('@') || token.is_empty() {
                continue;
            }
            if let Some(&idx) = names.get(token) {
                bns.anns[idx].is_alt = 1;
            }
        }
    }
    bns
}

#[doc = "Original function: bns_destroy:230"]
pub fn bns_destroy(bns: &mut bntseq_t) {
    bns.fp_pac = None;
    bns.ambs.clear();
    bns.anns.clear();
    bns.n_seqs = 0;
    bns.n_holes = 0;
}

#[doc = "Original function: add1:249"]
fn add1(
    seq: &SequenceRecord,
    bns: &mut bntseq_t,
    pac: &mut Vec<u8>,
    m_pac: &mut i64,
) {
    let offset = bns.anns.last().map(|p| p.offset + i64::from(p.len)).unwrap_or(0);
    let ann_index = bns.anns.len();
    bns.anns.push(bntann1_t {
        name: seq.name.clone(),
        anno: seq.comment.clone().unwrap_or_else(|| "(null)".to_string()),
        gi: 0,
        len: i32::try_from(seq.seq.len()).expect("sequence too long"),
        offset,
        n_ambs: 0,
        is_alt: 0,
    });

    let mut lasts = 0_u8;
    let mut rng_state = bns.seed as u64;
    for (i, &base) in seq.seq.iter().enumerate() {
        let mut c = nt4(base);
        if c >= 4 {
            if lasts == base {
                if let Some(last) = bns.ambs.last_mut() {
                    last.len += 1;
                }
            } else {
                bns.ambs.push(bntamb1_t {
                    len: 1,
                    offset: offset + i as i64,
                    amb: base,
                });
                bns.anns[ann_index].n_ambs += 1;
                bns.n_holes += 1;
            }
        }
        lasts = base;

        if c >= 4 {
            rng_state = rng_state.wrapping_mul(0x5deece66d).wrapping_add(0xb);
            c = ((rng_state >> 16) & 3) as u8;
        }
        if bns.l_pac == *m_pac {
            *m_pac <<= 1;
            let new_len = usize::try_from(*m_pac / 4).expect("pac resize overflow");
            pac.resize(new_len, 0);
        }
        set_pac(pac, bns.l_pac, c);
        bns.l_pac += 1;
    }

    bns.n_seqs += 1;
}

#[doc = "Original function: bns_fasta2bntseq:298"]
pub fn bns_fasta2bntseq<R: BufRead>(reader: R, prefix: &str, for_only: i32) -> i64 {
    let records = parse_fasta(reader);
    let mut bns = bntseq_t {
        seed: 11,
        ..Default::default()
    };
    let mut m_pac = 0x10000_i64;
    let mut pac = vec![0_u8; usize::try_from(m_pac / 4).expect("initial pac size")];

    for record in &records {
        add1(record, &mut bns, &mut pac, &mut m_pac);
    }

    if for_only == 0 {
        m_pac = ((bns.l_pac * 2 + 3) / 4) * 4;
        pac.resize(usize::try_from(m_pac / 4).expect("reverse pac resize"), 0);
        let original_l_pac = bns.l_pac;
        let mut l = original_l_pac - 1;
        while l >= 0 {
            let base = get_pac(&pac, l);
            set_pac(&mut pac, bns.l_pac, 3 - base);
            bns.l_pac += 1;
            if l == 0 {
                break;
            }
            l -= 1;
        }
    }

    let ret = bns.l_pac;
    let pac_name = format!("{prefix}.pac");
    {
        let fp = File::create(&pac_name).unwrap_or_else(|e| panic!("fail to open file '{pac_name}' : {e}"));
        let mut w = BufWriter::new(fp);
        let bytes_to_write = usize::try_from(bns.l_pac >> 2).expect("pac write size")
            + if (bns.l_pac & 3) == 0 { 0 } else { 1 };
        w.write_all(&pac[..bytes_to_write]).expect("write pac");
        if bns.l_pac % 4 == 0 {
            w.write_all(&[0]).expect("write pac pad");
        }
        w.write_all(&[u8::try_from(bns.l_pac % 4).expect("tail byte")]).expect("write pac tail");
        w.flush().expect("flush pac");
    }

    bns_dump(&bns, prefix);
    ret
}

#[doc = "Original function: bwa_fa2pac:359"]
pub fn bwa_fa2pac(argv: &[String]) -> i32 {
    let mut for_only = 0;
    let mut positional = Vec::new();
    for arg in &argv[1..] {
        if arg == "-f" {
            for_only = 1;
        } else {
            positional.push(arg.clone());
        }
    }
    if positional.is_empty() {
        eprintln!("Usage: bwa fa2pac [-f] <in.fasta> [<out.prefix>]");
        return 1;
    }
    let input = &positional[0];
    let prefix = positional.get(1).unwrap_or(input);
    let file = File::open(input).unwrap_or_else(|e| panic!("fail to open file '{input}' : {e}"));
    let reader: Box<dyn Read> = if input.ends_with(".gz") {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let reader = BufReader::new(reader);
    let _ = bns_fasta2bntseq(reader, prefix, for_only);
    0
}

#[doc = "Original function: bns_pos2rid:378"]
pub fn bns_pos2rid(bns: &bntseq_t, pos_f: i64) -> i32 {
    if pos_f >= bns.l_pac {
        return -1;
    }
    let mut left = 0_i32;
    let mut mid = 0_i32;
    let mut right = bns.n_seqs;
    while left < right {
        mid = (left + right) >> 1;
        let mid_usize = usize::try_from(mid).expect("mid");
        if pos_f >= bns.anns[mid_usize].offset {
            if mid == bns.n_seqs - 1 {
                break;
            }
            let next = usize::try_from(mid + 1).expect("mid+1");
            if pos_f < bns.anns[next].offset {
                break;
            }
            left = mid + 1;
        } else {
            right = mid;
        }
    }
    mid
}

#[doc = "Original function: bns_intv2rid:394"]
pub fn bns_intv2rid(bns: &bntseq_t, rb: i64, re: i64) -> i32 {
    let mut is_rev = 0;
    if rb < bns.l_pac && re > bns.l_pac {
        return -2;
    }
    assert!(rb <= re);
    let rid_b = bns_pos2rid(bns, bns_depos(bns, rb, &mut is_rev));
    let rid_e = if rb < re {
        bns_pos2rid(bns, bns_depos(bns, re - 1, &mut is_rev))
    } else {
        rid_b
    };
    if rid_b == rid_e { rid_b } else { -1 }
}

#[doc = "Original function: bns_cnt_ambi:404"]
pub fn bns_cnt_ambi(bns: &bntseq_t, pos_f: i64, len: i32, ref_id: Option<&mut i32>) -> i32 {
    if let Some(ref_id) = ref_id {
        *ref_id = bns_pos2rid(bns, pos_f);
    }
    let mut left = 0_i32;
    let mut right = bns.n_holes;
    let mut nn = 0_i32;
    while left < right {
        let mid = (left + right) >> 1;
        let amb = &bns.ambs[usize::try_from(mid).expect("mid")];
        if pos_f >= amb.offset + i64::from(amb.len) {
            left = mid + 1;
        } else if pos_f + i64::from(len) <= amb.offset {
            right = mid;
        } else {
            if pos_f >= amb.offset {
                let overlap_end = (amb.offset + i64::from(amb.len)).min(pos_f + i64::from(len));
                nn += i32::try_from(overlap_end - pos_f).expect("overlap");
            } else {
                let overlap_end = (amb.offset + i64::from(amb.len)).min(pos_f + i64::from(len));
                nn += i32::try_from(overlap_end - amb.offset).expect("overlap");
            }
            break;
        }
    }
    nn
}

#[doc = "Original function: bns_get_seq:427"]
pub fn bns_get_seq(l_pac: i64, pac: &[u8], mut beg: i64, mut end: i64, len: &mut i64) -> Option<Vec<u8>> {
    if end < beg {
        std::mem::swap(&mut beg, &mut end);
    }
    if end > l_pac << 1 {
        end = l_pac << 1;
    }
    if beg < 0 {
        beg = 0;
    }
    if beg >= l_pac || end <= l_pac {
        *len = end - beg;
        let mut seq = Vec::with_capacity(usize::try_from(end - beg + 64).expect("seq capacity"));
        if beg >= l_pac {
            let beg_f = (l_pac << 1) - 1 - end;
            let end_f = (l_pac << 1) - 1 - beg;
            let mut k = end_f;
            while k > beg_f {
                seq.push(3 - get_pac(pac, k));
                k -= 1;
            }
        } else {
            let mut k = beg;
            while k < end {
                seq.push(get_pac(pac, k));
                k += 1;
            }
        }
        Some(seq)
    } else {
        *len = 0;
        None
    }
}

#[doc = "Original function: bns_fetch_seq:453"]
pub fn bns_fetch_seq(
    bns: &bntseq_t,
    pac: &[u8],
    beg: &mut i64,
    mid: i64,
    end: &mut i64,
    rid: &mut i32,
) -> Vec<u8> {
    if *end < *beg {
        std::mem::swap(end, beg);
    }
    assert!(*beg <= mid && mid < *end);

    let mut is_rev = 0;
    *rid = bns_pos2rid(bns, bns_depos(bns, mid, &mut is_rev));
    let ann = &bns.anns[usize::try_from(*rid).expect("rid")];
    let mut far_beg = ann.offset;
    let mut far_end = far_beg + i64::from(ann.len);
    if is_rev != 0 {
        let tmp = far_beg;
        far_beg = (bns.l_pac << 1) - far_end;
        far_end = (bns.l_pac << 1) - tmp;
    }
    *beg = (*beg).max(far_beg);
    *end = (*end).min(far_end);

    let mut len = 0;
    let seq = bns_get_seq(bns.l_pac, pac, *beg, *end, &mut len)
        .unwrap_or_else(|| panic!("[E::bns_fetch_seq] begin={}, mid={}, end={}", *beg, mid, *end));
    assert_eq!(*end - *beg, len);
    seq
}

#[cfg(test)]
mod tests {
    use std::fs::{self, File};
    use std::io::{BufReader, Write};
    use std::path::{Path, PathBuf};
    use std::time::{SystemTime, UNIX_EPOCH};

    use crate::generated::bntseq_h::{bntamb1_t, bntann1_t, bntseq_t};
    use super::{
        add1, bns_cnt_ambi, bns_destroy, bns_dump, bns_fasta2bntseq, bns_fetch_seq, bns_get_seq,
        bns_intv2rid, bns_pos2rid, bns_restore, bns_restore_core, bwa_fa2pac, parse_fasta,
        SequenceRecord,
    };

    fn temp_prefix(name: &str) -> PathBuf {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_{name}_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        dir.join("ref")
    }

    fn sample_bns() -> bntseq_t {
        bntseq_t {
            l_pac: 8,
            n_seqs: 2,
            anns: vec![
                bntann1_t { offset: 0, len: 4, name: "chr1".into(), ..Default::default() },
                bntann1_t { offset: 4, len: 4, name: "chr2".into(), ..Default::default() },
            ],
            n_holes: 1,
            ambs: vec![bntamb1_t { offset: 2, len: 3, amb: b'N' }],
            ..Default::default()
        }
    }

    fn sample_pac() -> Vec<u8> {
        vec![(0 << 6) | (1 << 4) | (2 << 2) | 3, (1 << 6) | (0 << 4) | (3 << 2) | 2]
    }

    #[test]
    fn fasta_parser_preserves_name_comment_and_sequence() {
        let records = parse_fasta(BufReader::new(&b">chr1 comment one\nAC\nTG\n>chr2\nNN\n"[..]));
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, "chr1");
        assert_eq!(records[0].comment.as_deref(), Some("comment one"));
        assert_eq!(records[0].seq, b"ACTG");
        assert_eq!(records[1].name, "chr2");
        assert_eq!(records[1].comment, None);
    }

    #[test]
    fn bns_destroy_clears_owned_state() {
        let mut bns = bntseq_t {
            n_seqs: 2,
            n_holes: 1,
            anns: vec![bntann1_t::default(), bntann1_t::default()],
            ambs: vec![Default::default()],
            ..Default::default()
        };
        bns_destroy(&mut bns);
        assert_eq!(bns.n_seqs, 0);
        assert_eq!(bns.n_holes, 0);
        assert!(bns.anns.is_empty());
        assert!(bns.ambs.is_empty());
        assert!(bns.fp_pac.is_none());
    }

    #[test]
    fn bns_restore_core_reads_ann_amb_and_pac() {
        let prefix = temp_prefix("restore_core");
        fs::write(
            prefix.with_extension("ann"),
            "12 2 7\n1 chr1 comment one\n0 5 0\n2 chr2 (null)\n5 7 1\n",
        )
        .expect("write ann");
        fs::write(prefix.with_extension("amb"), "12 2 1\n9 3 N\n").expect("write amb");
        let mut pac = File::create(prefix.with_extension("pac")).expect("write pac");
        pac.write_all(&[1, 2, 3, 4]).expect("pac bytes");

        let bns = bns_restore_core(
            prefix.with_extension("ann").to_str().expect("utf8"),
            prefix.with_extension("amb").to_str().expect("utf8"),
            prefix.with_extension("pac").to_str().expect("utf8"),
        );

        assert_eq!(bns.l_pac, 12);
        assert_eq!(bns.n_seqs, 2);
        assert_eq!(bns.seed, 7);
        assert_eq!(bns.anns[0].name, "chr1");
        assert_eq!(bns.anns[0].anno, "comment one");
        assert_eq!(bns.anns[1].anno, "");
        assert_eq!(bns.ambs.len(), 1);
        assert!(bns.fp_pac.is_some());
        fs::remove_dir_all(prefix.parent().expect("temp dir")).expect("cleanup");
    }

    #[test]
    fn bns_restore_marks_alt_contigs() {
        let prefix = temp_prefix("restore");
        fs::write(
            prefix.with_extension("ann"),
            "12 2 7\n1 chr1 comment one\n0 5 0\n2 chr2 comment two\n5 7 1\n",
        )
        .expect("write ann");
        fs::write(prefix.with_extension("amb"), "12 2 0\n").expect("write amb");
        let mut pac = File::create(prefix.with_extension("pac")).expect("write pac");
        pac.write_all(&[1, 2, 3, 4]).expect("pac bytes");
        fs::write(prefix.with_extension("alt"), "@header\nchr2\tfoo\n").expect("write alt");

        let bns = bns_restore(prefix.to_str().expect("utf8"));
        assert_eq!(bns.anns[0].is_alt, 0);
        assert_eq!(bns.anns[1].is_alt, 1);
        fs::remove_dir_all(prefix.parent().expect("temp dir")).expect("cleanup");
    }

    #[test]
    fn add1_builds_ann_amb_and_pac() {
        let mut bns = bntseq_t { seed: 11, ..Default::default() };
        let mut pac = vec![0_u8; 8];
        let mut m_pac = 32_i64;
        let rec = SequenceRecord {
            name: "chr1".into(),
            comment: Some("comment".into()),
            seq: b"ACNN".to_vec(),
        };
        add1(&rec, &mut bns, &mut pac, &mut m_pac);
        assert_eq!(bns.n_seqs, 1);
        assert_eq!(bns.anns[0].name, "chr1");
        assert_eq!(bns.anns[0].anno, "comment");
        assert_eq!(bns.anns[0].n_ambs, 1);
        assert_eq!(bns.ambs.len(), 1);
        assert_eq!(bns.ambs[0].offset, 2);
        assert_eq!(bns.l_pac, 4);
    }

    #[test]
    fn bns_dump_round_trips_restore_core() {
        let prefix = temp_prefix("dump");
        let pac_path = prefix.with_extension("pac");
        let mut pac_file = File::create(&pac_path).expect("create pac");
        pac_file.write_all(&[1, 2, 3, 4]).expect("pac bytes");

        let bns = bntseq_t {
            l_pac: 12,
            n_seqs: 2,
            seed: 7,
            anns: vec![
                bntann1_t { gi: 1, name: "chr1".into(), anno: "comment one".into(), offset: 0, len: 5, n_ambs: 0, ..Default::default() },
                bntann1_t { gi: 2, name: "chr2".into(), anno: String::new(), offset: 5, len: 7, n_ambs: 1, ..Default::default() },
            ],
            n_holes: 1,
            ambs: vec![bntamb1_t { offset: 9, len: 3, amb: b'N' }],
            fp_pac: Some(File::open(&pac_path).expect("open pac")),
        };
        bns_dump(&bns, prefix.to_str().expect("utf8"));
        let restored = bns_restore_core(
            prefix.with_extension("ann").to_str().expect("utf8"),
            prefix.with_extension("amb").to_str().expect("utf8"),
            prefix.with_extension("pac").to_str().expect("utf8"),
        );
        assert_eq!(restored.l_pac, 12);
        assert_eq!(restored.anns[0].anno, "comment one");
        assert_eq!(restored.anns[1].anno, "");
        assert_eq!(restored.ambs[0].amb, b'N');
        fs::remove_dir_all(prefix.parent().expect("temp dir")).expect("cleanup");
    }

    #[test]
    fn bns_fasta2bntseq_writes_index_files() {
        let prefix = temp_prefix("fa2pac");
        let fasta = b">chr1 comment one\nACGT\n>chr2\nNN\n";
        let l_pac = bns_fasta2bntseq(BufReader::new(&fasta[..]), prefix.to_str().expect("utf8"), 1);
        assert_eq!(l_pac, 6);
        assert!(prefix.with_extension("pac").is_file());
        let restored = bns_restore_core(
            prefix.with_extension("ann").to_str().expect("utf8"),
            prefix.with_extension("amb").to_str().expect("utf8"),
            prefix.with_extension("pac").to_str().expect("utf8"),
        );
        assert_eq!(restored.n_seqs, 2);
        assert_eq!(restored.anns[0].name, "chr1");
        assert_eq!(restored.anns[1].n_ambs, 1);
        fs::remove_dir_all(prefix.parent().expect("temp dir")).expect("cleanup");
    }

    #[test]
    fn bwa_fa2pac_accepts_plain_fasta_args() {
        let prefix = temp_prefix("argv");
        let fasta_path = prefix.with_extension("fa");
        fs::write(&fasta_path, b">chr1\nACGT\n").expect("write fasta");
        let prefix_out = prefix.with_extension("out");
        let argv = vec![
            "bwa".to_string(),
            fasta_path.to_str().expect("utf8").to_string(),
            prefix_out.to_str().expect("utf8").to_string(),
        ];
        assert_eq!(bwa_fa2pac(&argv), 0);
        assert!(Path::new(&format!("{}.pac", prefix_out.display())).is_file());
        fs::remove_dir_all(prefix.parent().expect("temp dir")).expect("cleanup");
    }

    #[test]
    fn bns_pos2rid_binary_searches_reference_blocks() {
        let bns = sample_bns();
        assert_eq!(bns_pos2rid(&bns, 0), 0);
        assert_eq!(bns_pos2rid(&bns, 3), 0);
        assert_eq!(bns_pos2rid(&bns, 4), 1);
        assert_eq!(bns_pos2rid(&bns, 8), -1);
    }

    #[test]
    fn bns_intv2rid_handles_boundary_cases() {
        let bns = sample_bns();
        assert_eq!(bns_intv2rid(&bns, 1, 3), 0);
        assert_eq!(bns_intv2rid(&bns, 3, 5), -1);
        assert_eq!(bns_intv2rid(&bns, 7, 9), -2);
    }

    #[test]
    fn bns_cnt_ambi_counts_overlap_and_sets_ref_id() {
        let bns = sample_bns();
        let mut rid = -1;
        assert_eq!(bns_cnt_ambi(&bns, 1, 4, Some(&mut rid)), 3);
        assert_eq!(rid, 0);
        assert_eq!(bns_cnt_ambi(&bns, 5, 2, None), 0);
    }

    #[test]
    fn bns_get_seq_handles_forward_reverse_and_boundary() {
        let pac = sample_pac();
        let mut len = -1;
        assert_eq!(bns_get_seq(8, &pac, 1, 4, &mut len), Some(vec![1, 2, 3]));
        assert_eq!(len, 3);
        assert_eq!(bns_get_seq(8, &pac, 9, 12, &mut len), Some(vec![0, 3, 2]));
        assert_eq!(bns_get_seq(8, &pac, 7, 9, &mut len), None);
        assert_eq!(len, 0);
    }

    #[test]
    fn bns_fetch_seq_clamps_to_reference_interval() {
        let bns = sample_bns();
        let pac = sample_pac();
        let mut beg = 2;
        let mut end = 7;
        let mut rid = -1;
        let seq = bns_fetch_seq(&bns, &pac, &mut beg, 5, &mut end, &mut rid);
        assert_eq!(rid, 1);
        assert_eq!((beg, end), (4, 7));
        assert_eq!(seq, vec![1, 0, 3]);
    }
}
