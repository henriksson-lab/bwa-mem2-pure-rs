#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Port of `bwa-mem2/src/bwa.h` + `bwa-mem2/src/bwa.cpp`.

use std::sync::atomic::{AtomicBool, AtomicI32, Ordering};
use std::sync::{LazyLock, Mutex, RwLock};
use crate::bwa_mem2::bntseq::{bns_get_seq_into, bntseq_t};
use crate::bwa_mem2::kseq::{kseq_read, kseq_t};
use crate::bwa_mem2::kstring::{kputc, kputw, kstring_t};
use crate::bwa_mem2::ksw::ksw_global2;
use crate::bwa_mem2::utils::{ErrFile, err_fputc, err_fputs};

// --- bwa.h ---

#[doc = "Original struct: bseq1_t (bwa-mem2/src/bwa.h)"]
#[derive(Debug, Default, Clone)]
pub struct bseq1_t {
    pub l_seq: i32,
    pub id: i32,
    pub name: Option<String>,
    pub comment: Option<String>,
    pub seq: Option<String>,
    pub qual: Option<String>,
    pub sam: Option<String>,
    pub seq_nt4: Vec<u8>,
}

#[doc = "Original struct: bwaidx_t (bwa-mem2/src/bwa.h)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct bwaidx_t {
    pub _opaque: (),
}

// --- bwa.cpp ---

// Thread-local rseq scratch for bwa_gen_cigar2 (allocated once per ~100K alignments otherwise).
thread_local! {
    static BWA_GEN_CIGAR_RSEQ: std::cell::RefCell<Vec<u8>> =
        const { std::cell::RefCell::new(Vec::new()) };
    // Pool of reusable cigar/md buffers. bwa_gen_cigar2 takes from these pools and returns the
    // filled buffers as bwa_cigar_t fields (no per-call clone). After SAM emission, mem_reg2sam
    // returns each mem_aln_t's cigar/md back to the pool. Per 700K-read run with ~1M alignments
    // that's ~2M alloc/free cycles eliminated.
    pub(crate) static CIGAR_OPS_POOL: std::cell::RefCell<Vec<Vec<u32>>> =
        const { std::cell::RefCell::new(Vec::new()) };
    pub(crate) static CIGAR_MD_POOL: std::cell::RefCell<Vec<Vec<u8>>> =
        const { std::cell::RefCell::new(Vec::new()) };
}

#[inline]
pub(crate) fn take_cigar_buf() -> Vec<u32> {
    CIGAR_OPS_POOL.with(|p| p.borrow_mut().pop().unwrap_or_default())
}

#[inline]
pub(crate) fn return_cigar_buf(mut buf: Vec<u32>) {
    buf.clear();
    CIGAR_OPS_POOL.with(|p| p.borrow_mut().push(buf));
}

#[inline]
pub(crate) fn take_md_buf() -> Vec<u8> {
    CIGAR_MD_POOL.with(|p| p.borrow_mut().pop().unwrap_or_default())
}

#[inline]
pub(crate) fn return_md_buf(mut buf: Vec<u8>) {
    buf.clear();
    CIGAR_MD_POOL.with(|p| p.borrow_mut().push(buf));
}

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct bwa_cigar_t {
    pub cigar: Vec<u32>,
    pub md: Vec<u8>,
}

pub static BWA_VERBOSE: AtomicI32 = AtomicI32::new(3);
// RwLock instead of Mutex: read-mostly access pattern (set once at startup via bwa_set_rg,
// then read once per SAM record × 700K records × N threads). Mutex serialized all readers via
// lock cmpxchg (was ~13% of mem_aln2sam in -t 4 profile). RwLock allows concurrent shared reads.
pub static BWA_RG_ID: LazyLock<RwLock<String>> = LazyLock::new(|| RwLock::new(String::new()));
// Hot-path fast skip: when -R is not passed (the common case), every per-record SAM emission
// would otherwise acquire the RwLock just to discover rg_id is empty. This atomic flag lets us
// branch-out early with a single relaxed atomic load.
pub static BWA_RG_ID_NONEMPTY: AtomicBool = AtomicBool::new(false);
pub static BWA_PG: LazyLock<Mutex<Option<String>>> = LazyLock::new(|| Mutex::new(None));

/// Strip a trailing `/1` or `/2` read-number suffix from a FASTQ name in-place.
///
/// Mutates the kstring_t by shortening its logical length (and writing a NUL terminator) when
/// the last two bytes are `/` followed by an ASCII digit.
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

/// Move the current kseq record (name / comment / seq / qual) into a `bseq1_t`.
///
/// The kstring buffers in `ks` are taken (not copied) and reinterpreted as `String`s; bytes came
/// from FASTQ (printable ASCII per spec), so UTF-8 validation is skipped. `ks`'s fields are left
/// empty for the next `kseq_read` to refill.
///
/// TODO (upstream): it would be better to allocate one chunk of memory, but probably it does not
/// matter in practice.
#[doc = "Original function: kseq2bseq1:68"]
pub fn kseq2bseq1(ks: &mut kseq_t, s: &mut bseq1_t) {
    // Swap kstring buffers into bseq1_t instead of copying. The Vec<u8> in each kstring_t holds
    // the bytes plus a trailing 0 sentinel — truncate to logical length .l before reinterpreting
    // as String. Bytes came from FASTQ (printable ASCII by spec), so from_utf8_unchecked is safe.
    // ks's kstring buffers are left empty (.l = 0 stays consistent with empty .s); next
    // kseq_read pushes bytes into them via read_until.
    fn take_kstring_as_string(ks: &mut kstring_t) -> Option<String> {
        if ks.l == 0 {
            return None;
        }
        let mut bytes = std::mem::take(&mut ks.s);
        bytes.truncate(ks.l);
        ks.l = 0;
        ks.m = 0;
        Some(unsafe { String::from_utf8_unchecked(bytes) })
    }
    s.l_seq = ks.seq.l as i32;
    s.name = take_kstring_as_string(&mut ks.name);
    s.comment = take_kstring_as_string(&mut ks.comment);
    s.seq = take_kstring_as_string(&mut ks.seq);
    s.qual = take_kstring_as_string(&mut ks.qual);
}

/// Customized FASTA/Q chunk reader for MPI processing (upstream).
///
/// Reads chunks of records into `bseq1_t`. Currently a stub in this port.
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

/// Read a chunk of records from one or two FASTQ streams into `bseq1_t`s.
///
/// Stops once accumulated sequence length meets `chunk_size` (rounded up to an even count when
/// paired). Emits warnings if one mate file is shorter than the other. Writes the per-record
/// count to `n_` and the total sequence-byte count to `s`.
#[doc = "Original function: bseq_read_orig:170"]
pub fn bseq_read_orig(
    chunk_size: i64,
    n_: &mut i32,
    ks1: &mut kseq_t,
    ks2: Option<&mut kseq_t>,
    s: &mut i64,
) -> Vec<bseq1_t> {
    let mut size = 0_i64;
    // Pre-reserve seqs to a typical batch's record count: chunk_size bytes / ~200 bytes/record
    // ~= 100K records. Avoids the doubling-realloc churn during batch fill.
    let initial_cap = ((chunk_size / 200).max(1024) as usize).min(200_000);
    let mut seqs = Vec::with_capacity(initial_cap);
    let mut ks2 = ks2;
    while kseq_read(ks1) >= 0 {
        if let Some(ks2_) = ks2.as_deref_mut() {
            // the 2nd file has fewer reads
            if kseq_read(ks2_) < 0 {
                eprintln!("[W::bseq_read_orig] the 2nd file has fewer sequences.");
                break;
            }
        }

        trim_readno(&mut ks1.name);
        let mut seq = bseq1_t::default();
        kseq2bseq1(ks1, &mut seq);
        seq.id = seqs.len() as i32;
        size += i64::from(seq.l_seq);
        seqs.push(seq);

        if let Some(ks2_) = ks2.as_deref_mut() {
            trim_readno(&mut ks2_.name);
            let mut seq = bseq1_t::default();
            kseq2bseq1(ks2_, &mut seq);
            seq.id = seqs.len() as i32;
            size += i64::from(seq.l_seq);
            seqs.push(seq);
        }

        if size >= chunk_size && (seqs.len() & 1) == 0 {
            break;
        }
    }
    // test if the 2nd file is finished
    if size == 0 {
        if let Some(ks2_) = ks2.as_deref_mut() {
            if kseq_read(ks2_) >= 0 {
                eprintln!("[W::bseq_read_orig] the 1st file has fewer sequences.");
            }
        }
    }
    *n_ = seqs.len() as i32;
    *s = size;
    seqs
}

/// Convenience wrapper that reads a chunk from a single FASTA/Q text source.
///
/// Builds a transient `kseq_t` from `text` and delegates to `bseq_read_orig`.
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

/// Split a sequence batch into singletons and paired-end pairs by adjacent matching names.
///
/// Walks `seqs` and groups consecutive reads with identical names into `sep[1]` (paired);
/// everything else goes into `sep[0]` (singletons). Counts are written to `m[0]`/`m[1]`.
#[doc = "Original function: bseq_classify:226"]
pub fn bseq_classify(n: i32, seqs: &[bseq1_t], m: &mut [i32; 2], sep: &mut [Vec<bseq1_t>; 2]) {
    let limit = n.max(0) as usize;
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

    m[0] = single.len() as i32;
    m[1] = paired.len() as i32;
    sep[0] = single;
    sep[1] = paired;
}

/// Fill the 5x5 ACGTN scoring matrix.
///
/// Sets diagonal of the ACGT submatrix to `+a` (match) and off-diagonals to `-b` (mismatch).
/// Row and column 4 (the ambiguous base N) are filled with -1.
#[doc = "Original function: bwa_fill_scmat:248"]
#[inline]
pub fn bwa_fill_scmat(a: i32, b: i32, mat: &mut [i8; 25]) {
    let mut k = 0_usize;
    for i in 0..4 {
        for j in 0..4 {
            mat[k] = if i == j { a as i8 } else { -(b as i8) };
            k += 1;
        }
        // ambiguous base
        mat[k] = -1;
        k += 1;
    }
    // DEFAULT AMBIG row (N vs anything)
    while k < 25 {
        mat[k] = -1;
        k += 1;
    }
}

/// Generate the CIGAR when the alignment end points are known.
///
/// Runs banded global Needleman-Wunsch (`ksw_global2`) between the query and the reference
/// window `[rb, re)`, optionally producing the cigar ops and the NM/MD tag values. Rejects
/// alignments bridging the forward/reverse strand boundary. When `rb >= l_pac` (reverse strand),
/// query and reference are reversed for the DP so indels land at the leftmost position.
///
/// # Arguments
/// * `mat` - 5x5 ACGTN scoring matrix (see `bwa_fill_scmat`)
/// * `o_del`, `e_del`, `o_ins`, `e_ins` - gap open / extend penalties for del / ins
/// * `w_` - user-supplied band-width cap
/// * `l_pac` - packed-reference length (single-strand)
/// * `pac` - packed reference (2-bit)
/// * `l_query`, `query` - query length and 0-3 encoded query bytes
/// * `rb`, `re` - reference window (concatenated fwd|rev coordinates)
/// * `score` - written with the alignment score
/// * `n_cigar` - if `Some`, populated with the number of CIGAR ops emitted
/// * `NM` - if `Some`, populated with NM (mismatches + gap bases)
///
/// # Returns
/// `Some(bwa_cigar_t)` with the CIGAR ops and the MD string when `n_cigar` was requested, or
/// `None` for rejected windows.
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
    // reject if negative length or bridging the forward and reverse strand
    if l_query <= 0 || rb >= re || (rb < l_pac && re > l_pac) {
        return None;
    }

    let mut rlen = 0_i64;
    // Take the thread-local rseq buffer (capacity reused across ~100K calls).
    let mut rseq = BWA_GEN_CIGAR_RSEQ.with(|cell| std::mem::take(&mut *cell.borrow_mut()));
    if !bns_get_seq_into(l_pac, pac, rb, re, &mut rlen, &mut rseq) || re - rb != rlen {
        // possible if out of range
        BWA_GEN_CIGAR_RSEQ.with(|cell| *cell.borrow_mut() = rseq);
        return None;
    }

    let l_query_usize = l_query as usize;
    if rb >= l_pac {
        // then reverse both query and rseq; this is to ensure indels to be placed at the leftmost
        // position
        query[..l_query_usize].reverse();
        rseq.reverse();
    }

    let mut out_cigar: Option<Vec<u32>> = None;
    // no gap; no need to do DP. UPDATE (upstream): we come to this block now... FIXME: due to an
    // issue in mem_reg2aln(), we never come to this block. This does not affect accuracy, but it
    // hurts performance.
    if i64::from(l_query) == re - rb && w_ == 0 {
        if n_cigar_slot.is_some() {
            out_cigar = Some(vec![(l_query as u32) << 4]);
            if let Some(n_cigar) = n_cigar_slot.as_deref_mut() {
                *n_cigar = 1;
            }
        }
        // Hoist score into a local so the compiler doesn't reload from memory each iter.
        // Use unchecked indexing on mat (size 25, max index = 4*5+4 = 24) and rseq/query
        // (bounded by l_query_usize); per-byte load + LUT lookup is the bottleneck.
        let mut acc = 0_i32;
        unsafe {
            for i in 0..l_query_usize {
                let ref_b = *rseq.get_unchecked(i);
                let qer_b = *query.get_unchecked(i);
                let idx = usize::from(ref_b) * 5 + usize::from(qer_b);
                acc += i32::from(*mat.get_unchecked(idx));
            }
        }
        *score = acc;
    } else {
        // set the band-width: budget the max gap derivable from the query's half-score at the
        // configured gap-open/extend penalties, then take the tighter of the geometric band
        // (around |rlen-l_query|) and the user-supplied w_, with a floor of |rlen-l_query|+3.
        let half_query = (l_query + 1) >> 1;
        let max_ins =
            (((half_query * i32::from(mat[0]) - o_ins) as f64) / e_ins as f64 + 1.0) as i32;
        let max_del =
            (((half_query * i32::from(mat[0]) - o_del) as f64) / e_del as f64 + 1.0) as i32;
        let mut max_gap = max_ins.max(max_del);
        max_gap = max_gap.max(1);
        let band_delta = (rlen - i64::from(l_query)).abs() as i32;
        let mut w = (max_gap + band_delta + 1) >> 1;
        w = w.min(w_);
        let min_w = band_delta + 3;
        w = w.max(min_w);
        // NW alignment

        let mut n_cigar_tmp = 0_i32;
        // Take a cigar buffer from the pool. The buffer (with prior capacity) becomes the
        // returned bwa_cigar_t.cigar — no clone. mem_reg2sam returns it to the pool after SAM
        // emission.
        let mut cigar_tmp = take_cigar_buf();
        *score = ksw_global2(
            l_query,
            &query[..l_query_usize],
            rlen as i32,
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
        } else {
            return_cigar_buf(cigar_tmp);
        }
    }

    let mut md = take_md_buf();
    // compute NM and MD (append MD to CIGAR buffer in the C++; here we build md into a separate
    // kstring backed by a pooled Vec<u8>).
    if nm_slot.is_some() && n_cigar_slot.is_some() {
        // Take an md buffer from the pool, wrap in a temp kstring_t for kput* helpers, then
        // unwrap and truncate to logical length for the returned bwa_cigar_t.md. The kstring
        // invariant requires s.len() == m; the pool buf has cap but len=0, so resize up so
        // ensure_allocated sees a valid initial state.
        let mut buf = std::mem::take(&mut md);
        let cap = buf.capacity();
        if buf.len() < cap {
            buf.resize(cap, 0);
        }
        let mut str = kstring_t {
            l: 0,
            m: buf.len(),
            s: buf,
        };
        // Const-size LUTs let the compiler elide bounds checks on the per-base index.
        // rseq is 2-bit-encoded (0..=4); using min(4) to keep get_unchecked sound.
        const ACGTN_LUT: [u8; 5] = *b"ACGTN";
        const TGCAN_LUT: [u8; 5] = *b"TGCAN";
        let int2base: &[u8; 5] = if rb < l_pac { &ACGTN_LUT } else { &TGCAN_LUT };
        let cigar = out_cigar.as_ref().expect("cigar available with n_cigar");
        let mut x = 0_usize;
        let mut y = 0_usize;
        let mut u = 0_i32;
        let mut n_mm = 0_i32;
        let mut n_gap = 0_i32;
        for (k, cigar_op) in cigar.iter().enumerate() {
            let op = cigar_op & 0xf;
            let len = (cigar_op >> 4) as usize;
            // op codes: 0=match, 1=insertion, 2=deletion
            if op == 0 {
                for i in 0..len {
                    if query[x + i] != rseq[y + i] {
                        kputw(u, &mut str);
                        let base_idx = usize::from(rseq[y + i].min(4));
                        kputc(
                            i32::from(unsafe { *int2base.get_unchecked(base_idx) }),
                            &mut str,
                        );
                        n_mm += 1;
                        u = 0;
                    } else {
                        u += 1;
                    }
                }
                x += len;
                y += len;
            } else if op == 2 {
                // deletion: don't emit the leading "^..." block if D is the first or the last CIGAR
                if k > 0 && k + 1 < cigar.len() {
                    kputw(u, &mut str);
                    kputc(i32::from(b'^'), &mut str);
                    for i in 0..len {
                        let base_idx = usize::from(rseq[y + i].min(4));
                        kputc(
                            i32::from(unsafe { *int2base.get_unchecked(base_idx) }),
                            &mut str,
                        );
                    }
                    u = 0;
                    n_gap += len as i32;
                }
                y += len;
            } else if op == 1 {
                // insertion
                x += len;
                n_gap += len as i32;
            }
        }
        kputw(u, &mut str);
        if let Some(nm) = nm_slot.as_deref_mut() {
            *nm = n_mm + n_gap;
        }
        md = str.s;
        md.truncate(str.l);
    } else {
        return_md_buf(md);
        md = Vec::new();
    }

    if rb >= l_pac {
        // reverse back query (rseq is local scratch, no need to restore order)
        query[..l_query_usize].reverse();
    }

    // Restore the rseq buffer to the thread-local for the next call's reuse.
    BWA_GEN_CIGAR_RSEQ.with(|cell| *cell.borrow_mut() = rseq);

    out_cigar.map(|cigar| bwa_cigar_t { cigar, md })
}

/// Convenience wrapper around `bwa_gen_cigar2` using a single gap-open/extend pair `q`/`r` for
/// both deletions and insertions.
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

/// Emit the SAM header to `fp`.
///
/// If `hdr_line` is `None` (or contains no `@SQ`), one `@SQ\tSN:<name>\tLN:<len>` line is
/// emitted per reference sequence (with `\tAH:*` for ALT contigs). When `hdr_line` is provided
/// it is appended verbatim, followed by the `@PG` line stashed in `BWA_PG` (if any).
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

/// Replace C-style backslash escapes (`\t`, `\n`, `\r`, `\\`) with the literal byte.
///
/// Used to unescape user-supplied `@RG` lines so command-line input like `\t` produces a real
/// tab. Unrecognized escapes are passed through unchanged.
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

/// Parse the user-supplied `@RG` line, set `BWA_RG_ID`, and return the escaped line.
///
/// Validates that `s` begins with `@RG`, extracts the `ID:` token (at most 255 characters), and
/// stores it in `BWA_RG_ID` for use by `mem_aln2sam`. Returns `None` on any validation failure;
/// the escaped header line on success.
#[doc = "Original function: bwa_set_rg:583"]
pub fn bwa_set_rg(s: &str) -> Option<String> {
    *BWA_RG_ID.write().expect("bwa_rg_id lock") = String::new();
    BWA_RG_ID_NONEMPTY.store(false, Ordering::Relaxed);
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
    *BWA_RG_ID.write().expect("bwa_rg_id lock") = id.to_string();
    BWA_RG_ID_NONEMPTY.store(true, Ordering::Relaxed);
    Some(rg_line)
}

/// Append a user-supplied SAM header line to the accumulated header.
///
/// Skips lines that do not start with `@`. Only the newly-appended portion is run through
/// `bwa_escape`, preserving any earlier escape pass on the existing header (C++ bwa.cpp:623).
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
    use crate::bwa_mem2::bntseq::{bntann1_t, bntseq_t};
    use crate::bwa_mem2::bwa::bseq1_t;
    use crate::bwa_mem2::kseq::{kseq_read, kseq_t};
    use crate::bwa_mem2::kstring::kstring_t;
    use crate::bwa_mem2::utils::{err_fclose, err_xopen_core};
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
        kseq2bseq1(&mut ks, &mut out);
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
        assert_eq!(&*BWA_RG_ID.read().expect("lock"), "foo");
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
