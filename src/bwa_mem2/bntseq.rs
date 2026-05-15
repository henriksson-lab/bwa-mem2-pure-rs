#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Port of `bwa-mem2/src/bntseq.h` + `bwa-mem2/src/bntseq.cpp`.

use std::fs::File;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;
use flate2::read::MultiGzDecoder;

// --- bntseq.h ---

#[doc = "Original struct: bntann1_t (bwa-mem2/src/bntseq.h)"]
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct bntann1_t {
    pub offset: i64,
    pub len: i32,
    pub n_ambs: i32,
    pub gi: u32,
    pub is_alt: i32,
    pub name: String,
    pub anno: String,
}

#[doc = "Original struct: bntamb1_t (bwa-mem2/src/bntseq.h)"]
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct bntamb1_t {
    pub offset: i64,
    pub len: i32,
    pub amb: u8,
}

#[doc = "Original struct: bntseq_t (bwa-mem2/src/bntseq.h)"]
#[derive(Debug, Default)]
pub struct bntseq_t {
    pub l_pac: i64,
    pub n_seqs: i32,
    pub seed: u32,
    pub anns: Vec<bntann1_t>,
    pub n_holes: i32,
    pub ambs: Vec<bntamb1_t>,
    pub fp_pac: Option<File>,
}

/// Map a packed position to its forward-strand coordinate.
///
/// If `pos >= bns.l_pac` the position is on the reverse strand; this
/// returns the reflected forward-strand position and sets `is_rev = 1`.
/// Otherwise the position is returned unchanged with `is_rev = 0`.
#[doc = "Original function: bns_depos:87"]
#[inline]
pub fn bns_depos(bns: &bntseq_t, pos: i64, is_rev: &mut i32) -> i64 {
    if pos >= bns.l_pac {
        *is_rev = 1;
        (bns.l_pac << 1) - 1 - pos
    } else {
        *is_rev = 0;
        pos
    }
}

#[cfg(test)]
mod tests_h {
    use super::{bns_depos, bntseq_t};

    #[test]
    fn bns_depos_tracks_forward_and_reverse_coordinates() {
        let bns = bntseq_t {
            l_pac: 10,
            ..Default::default()
        };
        let mut is_rev = -1;
        assert_eq!(bns_depos(&bns, 3, &mut is_rev), 3);
        assert_eq!(is_rev, 0);
        assert_eq!(bns_depos(&bns, 12, &mut is_rev), 7);
        assert_eq!(is_rev, 1);
    }
}

// --- bntseq.cpp ---

fn parse_i64(text: &str, ctx: &str) -> i64 {
    text.parse()
        .unwrap_or_else(|_| panic!("parse error reading {ctx}"))
}

fn parse_i32(text: &str, ctx: &str) -> i32 {
    text.parse()
        .unwrap_or_else(|_| panic!("parse error reading {ctx}"))
}

fn parse_u32(text: &str, ctx: &str) -> u32 {
    text.parse()
        .unwrap_or_else(|_| panic!("parse error reading {ctx}"))
}

fn next_line(lines: &mut impl Iterator<Item = std::io::Result<String>>, ctx: &str) -> String {
    lines
        .next()
        .unwrap_or_else(|| panic!("Error reading {ctx} : Unexpected end of file"))
        .unwrap_or_else(|e| panic!("Error reading {ctx} : {e}"))
}

fn split_header<'a>(line: &'a str, ctx: &str) -> (&'a str, &'a str, &'a str) {
    let first_end = line
        .find(char::is_whitespace)
        .unwrap_or_else(|| panic!("parse error reading {ctx}"));
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

#[inline]
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
            current = Some(SequenceRecord {
                name,
                comment,
                seq: Vec::new(),
            });
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

/// Write the `.ann` and `.amb` companion files for a `bntseq_t`.
///
/// `.ann` lists each contig's gi, name, comment, offset, length, and
/// ambiguous-base count; `.amb` lists the ambiguous-base runs. Used
/// by `bns_fasta2bntseq` at the end of `bwa index`.
#[doc = "Original function: bns_dump:73"]
pub fn bns_dump(bns: &bntseq_t, prefix: &str) {
    let ann_name = format!("{prefix}.ann");
    let amb_name = format!("{prefix}.amb");

    // dump .ann
    {
        let fp = File::create(&ann_name)
            .unwrap_or_else(|e| panic!("fail to open file '{ann_name}' : {e}"));
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

    // dump .amb
    {
        let fp = File::create(&amb_name)
            .unwrap_or_else(|e| panic!("fail to open file '{amb_name}' : {e}"));
        let mut w = BufWriter::new(fp);
        writeln!(w, "{} {} {}", bns.l_pac, bns.n_seqs, bns.n_holes).expect("write amb header");
        for amb in &bns.ambs {
            writeln!(w, "{} {} {}", amb.offset, amb.len, amb.amb as char).expect("write amb body");
        }
        w.flush().expect("flush amb");
    }
}

/// Reconstruct a `bntseq_t` from explicit `.ann`/`.amb`/`.pac` paths.
///
/// Parses contig metadata and ambiguous-base runs, opens the `.pac`
/// for later random access, and asserts that the `.ann` and `.amb`
/// headers agree on `l_pac` and `n_seqs`.
#[doc = "Original function: bns_restore_core:106"]
pub fn bns_restore_core(ann_filename: &str, amb_filename: &str, pac_filename: &str) -> bntseq_t {
    // read .ann
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
        // read gi and sequence name; remainder of header line is the FASTA comment ("anno"),
        // skipping the leading space. " (null)" placeholder collapses to an empty anno.
        let seq_header = next_line(&mut ann_lines, ann_filename);
        let (gi, name, tail) = split_header(&seq_header, ann_filename);
        let anno = if tail.len() > 1 && tail != " (null)" {
            tail[1..].to_string()
        } else {
            String::new()
        };
        // read the rest: offset, len, n_ambs
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

    // read .amb
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

    // open .pac
    bns.fp_pac = Some(
        File::open(pac_filename)
            .unwrap_or_else(|e| panic!("fail to open file '{pac_filename}' : {e}")),
    );
    bns
}

/// Restore a `bntseq_t` from `prefix.{ann,amb,pac}` and optional `prefix.alt`.
///
/// Wraps `bns_restore_core`; when `prefix.alt` exists, each contig
/// whose name appears in the file is marked `is_alt = 1`.
#[doc = "Original function: bns_restore:188"]
pub fn bns_restore(prefix: &str) -> bntseq_t {
    let ann_filename = format!("{prefix}.ann");
    let amb_filename = format!("{prefix}.amb");
    let pac_filename = format!("{prefix}.pac");
    let alt_filename = format!("{prefix}.alt");
    let mut bns = bns_restore_core(&ann_filename, &amb_filename, &pac_filename);
    // read .alt file if present: mark contigs whose names appear there as ALT
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

/// Release the `.pac` file handle and clear annotation/ambiguity tables.
///
/// Mirrors C's `bns_destroy(bntseq_t *)` minus the `free(bns)` — the
/// caller still owns the struct in Rust.
#[doc = "Original function: bns_destroy:230"]
pub fn bns_destroy(bns: &mut bntseq_t) {
    bns.fp_pac = None;
    bns.ambs.clear();
    bns.anns.clear();
    bns.n_seqs = 0;
    bns.n_holes = 0;
}

/// Append one FASTA record into the growing `bntseq_t` and `.pac` buffer.
///
/// Translates ASCII bases via `nt4`, opens/extends ambiguous-base
/// (`N`) runs as `bntamb1_t`, packs 2-bit bases into `pac` (doubling
/// `m_pac` when full), and replaces `N`s with a random base drawn
/// from a shared LCG state (`rng_state`) that matches glibc
/// `srand48(bns.seed) ; lrand48()` semantics in upstream.
#[doc = "Original function: add1:249"]
fn add1(
    seq: &SequenceRecord,
    bns: &mut bntseq_t,
    pac: &mut Vec<u8>,
    m_pac: &mut i64,
    rng_state: &mut u64,
) {
    let offset = bns
        .anns
        .last()
        .map(|p| p.offset + i64::from(p.len))
        .unwrap_or(0);
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
    // C++ bntseq.cpp:315 calls srand48(bns.seed) ONCE in bns_fasta2bntseq; lrand48() state is then
    // shared across all add1 calls (one per record). Caller initializes `rng_state` once and
    // threads it in here so the state persists across records.
    for (i, &base) in seq.seq.iter().enumerate() {
        let mut c = nt4(base);
        if c >= 4 {
            // N (ambiguous base)
            if lasts == base {
                // contiguous N — extend the current hole rather than opening a new one
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

        // fill buffer: pack 2-bit base into the pac, replacing N with a random base.
        if c >= 4 {
            // Step LCG mod 2^48, then take bits 17..18 (== `(state >> 17) & 3`) — what
            // `lrand48() & 3` returns. The mask after multiply trims to 48 bits.
            *rng_state =
                rng_state.wrapping_mul(0x5deece66d).wrapping_add(0xb) & 0x0000_FFFF_FFFF_FFFF;
            c = ((*rng_state >> 17) & 3) as u8;
        }
        if bns.l_pac == *m_pac {
            // double the pac size
            *m_pac <<= 1;
            let new_len = usize::try_from(*m_pac / 4).expect("pac resize overflow");
            pac.resize(new_len, 0);
        }
        set_pac(pac, bns.l_pac, c);
        bns.l_pac += 1;
    }

    bns.n_seqs += 1;
}

/// Convert a FASTA stream to packed `.pac` + `.ann` + `.amb` indices.
///
/// Calls `add1` per record, optionally appends the reverse complement
/// of the concatenated forward strand (controlled by `for_only`), and
/// writes the final `.pac` with the upstream tail byte convention.
/// Returns the total packed length `l_pac` after the optional rc append.
///
/// # Arguments
/// * `reader` - FASTA input (gzip is handled by the caller)
/// * `prefix` - output filename prefix for `.pac`/`.ann`/`.amb`
/// * `for_only` - if nonzero, skip the reverse-complement half
#[doc = "Original function: bns_fasta2bntseq:298"]
pub fn bns_fasta2bntseq<R: BufRead>(reader: R, prefix: &str, for_only: i32) -> i64 {
    let records = parse_fasta(reader);
    let mut bns = bntseq_t {
        // fixed seed for random generator (matches C++ bntseq.cpp:314)
        seed: 11,
        ..Default::default()
    };
    let mut m_pac = 0x10000_i64;
    let mut pac = vec![0_u8; usize::try_from(m_pac / 4).expect("initial pac size")];

    // Mirror glibc srand48(bns.seed): X = ((seed & 0xFFFFFFFF) << 16) | 0x330e (48-bit state).
    // C++ bntseq.cpp:315 calls this once before the per-record add1 loop.
    let mut rng_state: u64 = (((bns.seed as u32) as u64) << 16) | 0x330e_u64;
    for record in &records {
        add1(record, &mut bns, &mut pac, &mut m_pac, &mut rng_state);
    }

    if for_only == 0 {
        // add the reverse complemented sequence after the forward strand
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
    // finalize .pac file
    {
        let fp = File::create(&pac_name)
            .unwrap_or_else(|e| panic!("fail to open file '{pac_name}' : {e}"));
        let mut w = BufWriter::new(fp);
        let bytes_to_write = usize::try_from(bns.l_pac >> 2).expect("pac write size")
            + if (bns.l_pac & 3) == 0 { 0 } else { 1 };
        w.write_all(&pac[..bytes_to_write]).expect("write pac");
        // the following codes make the pac file size always (l_pac/4+1+1)
        if bns.l_pac % 4 == 0 {
            w.write_all(&[0]).expect("write pac pad");
        }
        w.write_all(&[u8::try_from(bns.l_pac % 4).expect("tail byte")])
            .expect("write pac tail");
        // close .pac file
        w.flush().expect("flush pac");
    }

    bns_dump(&bns, prefix);
    ret
}

/// CLI driver for the `fa2pac` subcommand.
///
/// Parses `-f` (forward-only) plus 1-2 positional arguments
/// (`<in.fasta> [<out.prefix>]`), opens the input (transparently
/// handling `.gz`), and invokes `bns_fasta2bntseq`.
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

/// Look up the contig id (`rid`) bracketing a forward-strand position.
///
/// Returns `-1` if `pos_f >= bns.l_pac`. Otherwise binary-searches
/// `bns.anns` for the contig whose `[offset, offset + len)` brackets
/// `pos_f`.
#[doc = "Original function: bns_pos2rid:378"]
#[inline(always)]
pub fn bns_pos2rid(bns: &bntseq_t, pos_f: i64) -> i32 {
    if pos_f >= bns.l_pac {
        return -1;
    }
    let mut left = 0_i32;
    let mut mid = 0_i32;
    let mut right = bns.n_seqs;
    // binary search for the rid whose contig brackets pos_f
    while left < right {
        mid = (left + right) >> 1;
        let mid_usize = mid as usize;
        if pos_f >= bns.anns[mid_usize].offset {
            if mid == bns.n_seqs - 1 {
                break;
            }
            if pos_f < bns.anns[mid_usize + 1].offset {
                break; // bracketed
            }
            left = mid + 1;
        } else {
            right = mid;
        }
    }
    mid
}

/// Resolve a `[rb, re)` reference interval to a single `rid`.
///
/// Returns `-2` if the interval straddles the forward/reverse boundary,
/// the shared `rid` if both endpoints map to the same contig, or `-1`
/// if they fall on different contigs.
#[doc = "Original function: bns_intv2rid:394"]
#[inline(always)]
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
    if rid_b == rid_e {
        rid_b
    } else {
        -1
    }
}

/// Count ambiguous bases overlapping a `[pos_f, pos_f + len)` window.
///
/// Binary-searches `bns.ambs` for the first overlapping run and
/// returns the size of the intersection. Optionally writes the
/// owning contig's `rid` into `ref_id` via `bns_pos2rid`.
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
        let amb = &bns.ambs[mid as usize];
        if pos_f >= amb.offset + i64::from(amb.len) {
            left = mid + 1;
        } else if pos_f + i64::from(len) <= amb.offset {
            right = mid;
        } else {
            // overlap with an ambiguous-base run
            if pos_f >= amb.offset {
                let overlap_end = (amb.offset + i64::from(amb.len)).min(pos_f + i64::from(len));
                nn += (overlap_end - pos_f) as i32;
            } else {
                let overlap_end = (amb.offset + i64::from(amb.len)).min(pos_f + i64::from(len));
                nn += (overlap_end - amb.offset) as i32;
            }
            break;
        }
    }
    nn
}

/// Decode the 2-bit reference slice `[beg, end)` into a fresh byte vec.
///
/// Forward intervals walk `pac` directly; intervals on the reverse
/// half reflect back into the forward strand and emit
/// complemented bases. Intervals straddling the forward/reverse
/// boundary return `None`. Writes the decoded length into `len`.
#[doc = "Original function: bns_get_seq:427"]
pub fn bns_get_seq(l_pac: i64, pac: &[u8], beg: i64, end: i64, len: &mut i64) -> Option<Vec<u8>> {
    let mut seq = Vec::new();
    if bns_get_seq_into(l_pac, pac, beg, end, len, &mut seq) {
        Some(seq)
    } else {
        None
    }
}

/// Variant of `bns_get_seq` that writes into a caller-provided buffer (cleared first).
/// Returns true if the requested range was valid; false otherwise.
pub fn bns_get_seq_into(
    l_pac: i64,
    pac: &[u8],
    mut beg: i64,
    mut end: i64,
    len: &mut i64,
    out: &mut Vec<u8>,
) -> bool {
    // if end is smaller, swap
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
        out.clear();
        out.reserve((end - beg + 64) as usize);
        if beg >= l_pac {
            // reverse strand: mirror the requested interval back into the forward half of the
            // pac and emit base complements as we walk backward.
            let beg_f = (l_pac << 1) - 1 - end;
            let end_f = (l_pac << 1) - 1 - beg;
            // Walk backward over (beg_f, end_f] complementing each base.
            decode_pac_reverse_complement_into(pac, beg_f + 1, end_f + 1, out);
        } else {
            // forward strand: walk [beg, end).
            decode_pac_forward_into(pac, beg, end, out);
        }
        true
    } else {
        // bridging the forward-reverse boundary: return nothing
        *len = 0;
        false
    }
}

/// Decode 4 bases per byte for a forward range [beg, end) and append to `out`.
/// pac is the bit-packed reference (4 bases per byte, bit-position 6,4,2,0 for offsets 0,1,2,3).
#[inline]
fn decode_pac_forward_into(pac: &[u8], beg: i64, end: i64, out: &mut Vec<u8>) {
    let mut k = beg;
    // Prefix: bring k to a 4-byte boundary
    while k < end && (k & 3) != 0 {
        out.push(get_pac(pac, k));
        k += 1;
    }
    // Middle: read each byte once, emit 4 bases via extend_from_slice (single bounds check).
    while k + 4 <= end {
        let byte = unsafe { *pac.get_unchecked((k >> 2) as usize) };
        out.extend_from_slice(&[
            (byte >> 6) & 3,
            (byte >> 4) & 3,
            (byte >> 2) & 3,
            byte & 3,
        ]);
        k += 4;
    }
    // Suffix
    while k < end {
        out.push(get_pac(pac, k));
        k += 1;
    }
}

/// Decode 4 bases per byte for a reverse-complement walk over (beg, end] reversed and append.
/// We emit `3 - base` for each k from end-1 down to beg.
#[inline]
fn decode_pac_reverse_complement_into(pac: &[u8], beg: i64, end: i64, out: &mut Vec<u8>) {
    // Walk k from end-1 down to beg, outputting 3 - get_pac(pac, k).
    let mut k = end - 1;
    // Suffix (high-k side): align k+1 down to a 4-byte boundary
    while k >= beg && (k & 3) != 3 {
        out.push(3 - get_pac(pac, k));
        k -= 1;
    }
    // Middle: read each byte once, emit 4 complemented bases in reverse via extend_from_slice.
    while k - 3 >= beg {
        let byte = unsafe { *pac.get_unchecked((k >> 2) as usize) };
        // For this byte, k corresponds to offset 3 (lowest 2 bits), then offsets 2, 1, 0.
        out.extend_from_slice(&[
            3 - (byte & 3),
            3 - ((byte >> 2) & 3),
            3 - ((byte >> 4) & 3),
            3 - ((byte >> 6) & 3),
        ]);
        k -= 4;
    }
    // Prefix (low-k side)
    while k >= beg {
        out.push(3 - get_pac(pac, k));
        k -= 1;
    }
}

/// Fetch the reference slice around `mid`, clipped to the owning contig.
///
/// Determines the contig containing `mid` via `bns_pos2rid`, clamps
/// the caller-supplied `[beg, end)` to that contig's range
/// (flipping for reverse-strand `mid`), and decodes the resulting
/// window. `beg`, `end`, and `rid` are updated in place.
#[doc = "Original function: bns_fetch_seq:453"]
pub fn bns_fetch_seq(
    bns: &bntseq_t,
    pac: &[u8],
    beg: &mut i64,
    mid: i64,
    end: &mut i64,
    rid: &mut i32,
) -> Vec<u8> {
    let mut out = Vec::new();
    bns_fetch_seq_into(bns, pac, beg, mid, end, rid, &mut out);
    out
}

/// Variant of `bns_fetch_seq` that writes into a caller-provided buffer.
pub fn bns_fetch_seq_into(
    bns: &bntseq_t,
    pac: &[u8],
    beg: &mut i64,
    mid: i64,
    end: &mut i64,
    rid: &mut i32,
    out: &mut Vec<u8>,
) {
    // if end is smaller, swap
    if *end < *beg {
        std::mem::swap(end, beg);
    }
    assert!(*beg <= mid && mid < *end);

    let mut is_rev = 0;
    *rid = bns_pos2rid(bns, bns_depos(bns, mid, &mut is_rev));
    let ann = &bns.anns[*rid as usize];
    let mut far_beg = ann.offset;
    let mut far_end = far_beg + i64::from(ann.len);
    if is_rev != 0 {
        // flip to the reverse strand
        let tmp = far_beg;
        far_beg = (bns.l_pac << 1) - far_end;
        far_end = (bns.l_pac << 1) - tmp;
    }
    *beg = (*beg).max(far_beg);
    *end = (*end).min(far_end);

    let mut len = 0;
    let ok = bns_get_seq_into(bns.l_pac, pac, *beg, *end, &mut len, out);
    if !ok {
        panic!(
            "[E::bns_fetch_seq] begin={}, mid={}, end={}",
            *beg, mid, *end
        );
    }
    assert_eq!(*end - *beg, len);
}

#[cfg(test)]
mod tests {
    use std::fs::{self, File};
    use std::io::{BufReader, Write};
    use std::path::{Path, PathBuf};
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::{
        add1, bns_cnt_ambi, bns_destroy, bns_dump, bns_fasta2bntseq, bns_fetch_seq, bns_get_seq,
        bns_intv2rid, bns_pos2rid, bns_restore, bns_restore_core, bwa_fa2pac, parse_fasta,
        SequenceRecord,
    };
    use crate::bwa_mem2::bntseq::{bntamb1_t, bntann1_t, bntseq_t};

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
                bntann1_t {
                    offset: 0,
                    len: 4,
                    name: "chr1".into(),
                    ..Default::default()
                },
                bntann1_t {
                    offset: 4,
                    len: 4,
                    name: "chr2".into(),
                    ..Default::default()
                },
            ],
            n_holes: 1,
            ambs: vec![bntamb1_t {
                offset: 2,
                len: 3,
                amb: b'N',
            }],
            ..Default::default()
        }
    }

    fn sample_pac() -> Vec<u8> {
        vec![
            (0 << 6) | (1 << 4) | (2 << 2) | 3,
            (1 << 6) | (0 << 4) | (3 << 2) | 2,
        ]
    }

    #[test]
    fn fasta_parser_preserves_name_comment_and_sequence() {
        let records = parse_fasta(BufReader::new(
            &b">chr1 comment one\nAC\nTG\n>chr2\nNN\n"[..],
        ));
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
        let mut bns = bntseq_t {
            seed: 11,
            ..Default::default()
        };
        let mut pac = vec![0_u8; 8];
        let mut m_pac = 32_i64;
        let rec = SequenceRecord {
            name: "chr1".into(),
            comment: Some("comment".into()),
            seq: b"ACNN".to_vec(),
        };
        let mut rng_state: u64 = (((bns.seed as u32) as u64) << 16) | 0x330e_u64;
        add1(&rec, &mut bns, &mut pac, &mut m_pac, &mut rng_state);
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
                bntann1_t {
                    gi: 1,
                    name: "chr1".into(),
                    anno: "comment one".into(),
                    offset: 0,
                    len: 5,
                    n_ambs: 0,
                    ..Default::default()
                },
                bntann1_t {
                    gi: 2,
                    name: "chr2".into(),
                    anno: String::new(),
                    offset: 5,
                    len: 7,
                    n_ambs: 1,
                    ..Default::default()
                },
            ],
            n_holes: 1,
            ambs: vec![bntamb1_t {
                offset: 9,
                len: 3,
                amb: b'N',
            }],
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
        let l_pac = bns_fasta2bntseq(
            BufReader::new(&fasta[..]),
            prefix.to_str().expect("utf8"),
            1,
        );
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
