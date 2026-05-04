#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Generated scaffold for `bwa-mem2/src/bwamem.cpp`.

use std::cell::Cell;
use std::collections::HashMap;
use std::sync::{Arc, LazyLock, Mutex};

use rayon::{scope, ThreadPoolBuilder};

use crate::generated::bandedswa_cpp::BandedPairWiseSW;
use crate::generated::bandedswa_h::SeqPair;
use crate::generated::bntseq_cpp::bns_intv2rid;
use crate::generated::bntseq_cpp::bns_pos2rid;
use crate::generated::bntseq_h::bns_depos;
use crate::generated::bntseq_h::bntseq_t;
use crate::generated::bwa_cpp::{bwa_fill_scmat, bwa_gen_cigar2, BWA_VERBOSE};
use crate::generated::bwa_h::bseq1_t;
use crate::generated::bwamem_extra_cpp::mem_gen_alt;
use crate::generated::bwamem_h::{
    mem_aln_t, mem_alnreg_t, mem_alnreg_v, mem_cache, mem_chain_t, mem_chain_v, mem_opt_t,
    mem_pestat_t, mem_seed_t, smem_aux_t, worker_t,
};
use crate::generated::bwamem_pair_cpp::{
    mem_pestat, mem_sam_pe_batch, mem_sam_pe_batch_post, mem_sam_pe_batch_pre,
};
use crate::generated::fmi_search_cpp::FMI_search;
use crate::generated::fmi_search_h::SMEM;
use crate::generated::kstring_cpp::ksprintf;
use crate::generated::kstring_h::{kputc, kputl, kputs, kputsn, kputw, ks_resize, kstring_t};
use crate::generated::ksw_cpp::ksw_align2;
use crate::generated::ksw_h::{kswr_t, KSW_XBYTE, KSW_XSTART};
use crate::generated::macro_h::{BATCH_SIZE, SEEDS_PER_CHAIN};
use crate::generated::utils_h::hash_64;

const PATCH_MAX_R_BW: f64 = 0.05;
const PATCH_MIN_SC_RATIO: f64 = 0.90;
const MEM_SHORT_EXT: i32 = 50;
const MEM_SHORT_LEN: i32 = 200;
const MEM_HSP_COEF: f64 = 1.1;
const MEM_MINSC_COEF: f64 = 5.5;
const MEM_SEEDSW_COEF: f64 = 0.05;
const MAX_SEQ_LEN8: usize = 128;
const MAX_SEQ_LEN16: usize = 32768;
const MAX_BAND_TRY: i32 = 2;
const H0_: i32 = -99;
const AVG_SEEDS_PER_READ: usize = 64;
const MEM_MAPQ_COEF: f64 = 30.0;
const MEM_MAPQ_MAX: i32 = 60;
const MEM_F_PE: i32 = 0x2;
const MEM_F_PRIMARY5: i32 = 0x800;
const MEM_F_ALL: i32 = 0x8;
const MEM_F_NO_MULTI: i32 = 0x10;
const MEM_F_KEEP_SUPP_MAPQ: i32 = 0x1000;
const MEM_F_REF_HDR: i32 = 0x100;
const MEM_F_SOFTCLIP: i32 = 0x200;

static COMPUTE_POOLS: LazyLock<Mutex<HashMap<usize, Arc<rayon::ThreadPool>>>> =
    LazyLock::new(|| Mutex::new(HashMap::new()));

// 256-byte lookup table for ASCII-to-2-bit encoding. Replaces the per-char match in seq_to_nt4
// — a single load + index in the inner loop instead of a chained-compare branch ladder. The
// match is also auto-vectorizable but at the cost of a wider scalar op count per byte; the LUT
// is auto-vectorized by LLVM with VPGATHERDD on AVX2/AVX-512 targets.
const NT4_LUT: [u8; 256] = {
    let mut t = [4_u8; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t[b'U' as usize] = 3;
    t[b'u' as usize] = 3;
    t
};

#[inline]
fn seq_to_nt4(seq: &str) -> Vec<u8> {
    seq.as_bytes()
        .iter()
        .map(|&c| NT4_LUT[c as usize])
        .collect()
}

#[inline]
fn clear_vec_with_cap<T>(v: &mut Vec<T>, keep_cap: usize) {
    v.clear();
    if v.capacity() > keep_cap {
        v.shrink_to(keep_cap);
    }
}

const KEEP_REG_CAP: usize = 16;
const KEEP_CHAIN_CAP: usize = 16;

// Per-chain seeds-buffer pool. Each new mem_chain_t starts with a Vec::with_capacity(1) for its
// seeds; without pooling that's ~25K alloc/free cycles per chunk × 35 chunks per 700K-read run.
// The pool drains seeds Vecs out of mem_chain_t before drop, preserving their capacity for reuse.
thread_local! {
    static CHAIN_SEEDS_POOL: std::cell::RefCell<Vec<Vec<mem_seed_t>>> =
        const { std::cell::RefCell::new(Vec::new()) };
}

#[inline]
fn take_seeds_buf() -> Vec<mem_seed_t> {
    // Start fresh buffers at capacity 8: typical chains hold ~5-10 seeds. With cap=1 the
    // doubling growth (1→2→4→8→16) costs 4 extra reallocs per cold chain. Pool reuse means
    // most calls hit the pop() path; this only affects the warm-up batch.
    CHAIN_SEEDS_POOL.with(|c| {
        c.borrow_mut().pop().unwrap_or_else(|| Vec::with_capacity(8))
    })
}

#[inline]
fn return_seeds_buf(mut buf: Vec<mem_seed_t>) {
    buf.clear();
    CHAIN_SEEDS_POOL.with(|c| c.borrow_mut().push(buf));
}

#[inline]
fn clear_chain_a_pooled(v: &mut Vec<mem_chain_t>, keep_cap: usize) {
    for c in v.drain(..) {
        return_seeds_buf(c.seeds);
    }
    if v.capacity() > keep_cap {
        v.shrink_to(keep_cap);
    }
}

fn debug_rss(label: &str) {
    if std::env::var_os("BWA_DEBUG_RSS").is_none() {
        return;
    }
    if let Ok(status) = std::fs::read_to_string("/proc/self/status") {
        let rss = status
            .lines()
            .find(|line| line.starts_with("VmRSS:"))
            .unwrap_or("VmRSS:\t?");
        let hwm = status
            .lines()
            .find(|line| line.starts_with("VmHWM:"))
            .unwrap_or("VmHWM:\t?");
        eprintln!("[rss::{label}] {rss}; {hwm}");
    }
}

#[inline]
fn trim_allocator_rss() {
    #[cfg(all(target_os = "linux", target_env = "gnu"))]
    unsafe {
        libc::malloc_trim(0);
    }
}

#[inline]
pub fn nt4_cow(seq: &crate::generated::bwa_h::bseq1_t) -> std::borrow::Cow<'_, [u8]> {
    if !seq.seq_nt4.is_empty() || seq.l_seq == 0 {
        std::borrow::Cow::Borrowed(&seq.seq_nt4)
    } else {
        std::borrow::Cow::Owned(seq_to_nt4(seq.seq.as_deref().unwrap_or("")))
    }
}

fn pac_to_reference_layout(l_pac: i64, pac: &[u8]) -> Vec<u8> {
    let l_pac_usize = usize::try_from(l_pac).expect("l_pac");
    let mut forward = vec![0_u8; l_pac_usize];
    for (i, base) in forward.iter_mut().enumerate() {
        let shift = (((!(i as i64)) & 3) << 1) as u8;
        *base = (pac[i >> 2] >> shift) & 3;
    }
    let mut ref_string = forward.clone();
    ref_string.extend(forward.iter().rev().map(|&b| if b < 4 { 3 - b } else { b }));
    ref_string
}

fn debug_trace_read(name: Option<&str>) -> bool {
    static TARGET: std::sync::OnceLock<Option<String>> = std::sync::OnceLock::new();
    let Some(target) = TARGET
        .get_or_init(|| std::env::var("BWA_MEM2_RS_TRACE_READ").ok())
        .as_deref()
    else {
        return false;
    };
    !target.is_empty() && name == Some(target)
}

fn debug_timings() -> bool {
    static FLAG: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *FLAG.get_or_init(|| std::env::var_os("BWA_MEM2_RS_TIMINGS").is_some())
}

#[derive(Clone, Copy)]
struct WorkerPtr(*mut worker_t);

unsafe impl Send for WorkerPtr {}
unsafe impl Sync for WorkerPtr {}

impl WorkerPtr {
    unsafe fn run(self, func: fn(&mut worker_t, i32, i32, i32), start: i32, len: i32, tid: i32) {
        func(&mut *self.0, start, len, tid);
    }
}

fn compute_pool(nthreads: usize) -> Arc<rayon::ThreadPool> {
    let mut pools = COMPUTE_POOLS.lock().expect("compute pool cache");
    pools
        .entry(nthreads)
        .or_insert_with(|| {
            Arc::new(
                ThreadPoolBuilder::new()
                    .num_threads(nthreads)
                    .build()
                    .expect("rayon thread pool"),
            )
        })
        .clone()
}

pub fn with_current_rayon_pool<R>(f: impl FnOnce() -> R) -> R {
    struct Reset(bool);
    impl Drop for Reset {
        fn drop(&mut self) {
            USE_CURRENT_RAYON_POOL.with(|flag| flag.set(self.0));
        }
    }

    let previous = USE_CURRENT_RAYON_POOL.with(|flag| {
        let previous = flag.get();
        flag.set(true);
        previous
    });
    let _reset = Reset(previous);
    f()
}

fn run_worker_chunks_parallel(
    worker: &mut worker_t,
    n: i32,
    func: fn(&mut worker_t, i32, i32, i32),
) {
    let n_usize = n.max(0) as usize;
    if n_usize == 0 {
        return;
    }
    let nthreads = worker.nthreads.max(1) as usize;
    let chunk_size = BATCH_SIZE;
    let nchunks = n_usize.div_ceil(chunk_size);
    if nthreads <= 1 || nchunks <= 1 {
        for chunk_idx in 0..nchunks {
            let start = chunk_idx * chunk_size;
            let end = (start + chunk_size).min(n_usize);
            func(worker, start as i32, (end - start) as i32, 0);
        }
        return;
    }

    let worker_ptr = WorkerPtr(worker as *mut worker_t);
    let run_chunks = || {
        scope(|s| {
            for tid in 0..nthreads {
                s.spawn(move |_| {
                    let mut chunk_idx = tid;
                    while chunk_idx < nchunks {
                        let start = chunk_idx * chunk_size;
                        let end = (start + chunk_size).min(n_usize);
                        unsafe {
                            // Chunks are disjoint in read-space, and `tid` selects distinct
                            // per-thread scratch lanes in `mem_cache`, matching the original kt_for boundary.
                            worker_ptr.run(
                                func,
                                start as i32,
                                (end - start) as i32,
                                tid as i32,
                            );
                        }
                        chunk_idx += nthreads;
                    }
                });
            }
        });
    };
    if USE_CURRENT_RAYON_POOL.with(Cell::get) {
        run_chunks();
    } else {
        let pool = compute_pool(nthreads);
        pool.install(run_chunks);
    }
}

fn ensure_mem_cache_thread_slots(mmc: &mut mem_cache, nthreads: usize) {
    mmc.seqBufLeftRef.resize_with(nthreads, Vec::new);
    mmc.seqBufLeftQer.resize_with(nthreads, Vec::new);
    mmc.seqBufRightRef.resize_with(nthreads, Vec::new);
    mmc.seqBufRightQer.resize_with(nthreads, Vec::new);
    mmc.wsize_buf_ref.resize(nthreads, 0);
    mmc.wsize_buf_qer.resize(nthreads, 0);

    mmc.seqPairArrayAux.resize_with(nthreads, Vec::new);
    mmc.seqPairArrayLeft128.resize_with(nthreads, Vec::new);
    mmc.seqPairArrayRight128.resize_with(nthreads, Vec::new);
    mmc.wsize.resize(nthreads, 0);

    mmc.wsize_mem.resize(nthreads, 0);
    mmc.wsize_mem_s.resize(nthreads, 0);
    mmc.wsize_mem_r.resize(nthreads, 0);
    mmc.matchArray.resize_with(nthreads, Vec::new);
    mmc.min_intv_ar.resize_with(nthreads, Vec::new);
    mmc.query_pos_ar.resize_with(nthreads, Vec::new);
    mmc.enc_qdb.resize_with(nthreads, Vec::new);
    mmc.rid.resize_with(nthreads, Vec::new);
    mmc.lim.resize_with(nthreads, || vec![0; BATCH_SIZE + 32]);
}

fn run_worker_bwt_aln_chunks_serial(worker: &mut worker_t, n: i32) {
    let n_usize = n.max(0) as usize;
    if n_usize == 0 {
        return;
    }
    let chunk_size = BATCH_SIZE;
    let nchunks = n_usize.div_ceil(chunk_size);
    for chunk_idx in 0..nchunks {
        let start = chunk_idx * chunk_size;
        let end = (start + chunk_size).min(n_usize);
        let start_i32 = start as i32;
        let len_i32 = (end - start) as i32;
        worker_bwt(worker, start_i32, len_i32, 0);
        worker_aln(worker, start_i32, len_i32, 0);
        trim_allocator_rss();
    }
}

// Fused parallel runner: each thread does bwt+aln on its chunks back-to-back, eliminating
// the rayon scope-join barrier between phases. Equivalent to two run_worker_chunks_parallel
// calls but avoids the inter-phase barrier overhead.
fn run_worker_bwt_aln_chunks_parallel(worker: &mut worker_t, n: i32) {
    let n_usize = n.max(0) as usize;
    if n_usize == 0 {
        return;
    }
    let nthreads = worker.nthreads.max(1) as usize;
    let chunk_size = BATCH_SIZE;
    let nchunks = n_usize.div_ceil(chunk_size);
    if nthreads <= 1 || nchunks <= 1 {
        for chunk_idx in 0..nchunks {
            let start = chunk_idx * chunk_size;
            let end = (start + chunk_size).min(n_usize);
            let start_i32 = start as i32;
            let len_i32 = (end - start) as i32;
            worker_bwt(worker, start_i32, len_i32, 0);
            worker_aln(worker, start_i32, len_i32, 0);
        }
        return;
    }

    let worker_ptr = WorkerPtr(worker as *mut worker_t);
    let run_chunks = || {
        scope(|s| {
            for tid in 0..nthreads {
                s.spawn(move |_| {
                    let mut chunk_idx = tid;
                    while chunk_idx < nchunks {
                        let start = chunk_idx * chunk_size;
                        let end = (start + chunk_size).min(n_usize);
                        unsafe {
                            // Same chunk on same thread for bwt + aln — preserves scratch reuse.
                            worker_ptr.run(
                                worker_bwt,
                                start as i32,
                                (end - start) as i32,
                                tid as i32,
                            );
                            worker_ptr.run(
                                worker_aln,
                                start as i32,
                                (end - start) as i32,
                                tid as i32,
                            );
                        }
                        chunk_idx += nthreads;
                    }
                });
            }
        });
    };
    if USE_CURRENT_RAYON_POOL.with(Cell::get) {
        run_chunks();
    } else {
        let pool = compute_pool(nthreads);
        pool.install(run_chunks);
    }
}

#[inline(always)]
fn chn_beg(ch: &mem_chain_t) -> i32 {
    ch.seeds[0].qbeg
}

#[inline(always)]
fn chn_end(ch: &mem_chain_t) -> i32 {
    let last = &ch.seeds[(ch.n - 1) as usize];
    last.qbeg + last.len
}

#[inline]
fn recompute_seedcov_for_aln(chain: &mem_chain_t, a: &mut mem_alnreg_t) {
    if a.rb == H0_ as i64 || a.qb == H0_ || a.qe == H0_ || a.re == H0_ as i64 {
        return;
    }
    a.seedcov = 0;
    for t in chain.seeds.iter().take(chain.n as usize) {
        if t.qbeg >= a.qb
            && t.qbeg + t.len <= a.qe
            && t.rbeg >= a.rb
            && t.rbeg + i64::from(t.len) <= a.re
        {
            a.seedcov += t.len;
        }
    }
}

fn stable_alnreg_hash(id: i64, p: &mem_alnreg_t) -> u64 {
    let mut x = hash_64(u64::from_ne_bytes(id.to_ne_bytes()));
    x = hash_64(x ^ u64::from_ne_bytes(p.rb.to_ne_bytes()));
    x = hash_64(x ^ u64::from_ne_bytes(p.re.to_ne_bytes()));
    x = hash_64(x ^ u64::from_ne_bytes(i64::from(p.qb).to_ne_bytes()));
    x = hash_64(x ^ u64::from_ne_bytes(i64::from(p.qe).to_ne_bytes()));
    x = hash_64(x ^ u64::from_ne_bytes(i64::from(p.rid).to_ne_bytes()));
    x = hash_64(x ^ u64::from_ne_bytes(i64::from(p.score).to_ne_bytes()));
    x = hash_64(x ^ u64::from_ne_bytes(i64::from(p.seedcov).to_ne_bytes()));
    hash_64(x ^ u64::from(p.is_alt))
}

#[doc = "Original function: cal_max_gap:66"]
#[inline]
pub fn cal_max_gap(opt: &mem_opt_t, qlen: i32) -> i32 {
    let l_del = (((qlen * opt.a - opt.o_del) as f64) / opt.e_del as f64 + 1.0) as i32;
    let l_ins = (((qlen * opt.a - opt.o_ins) as f64) / opt.e_ins as f64 + 1.0) as i32;
    let mut l = if l_del > l_ins { l_del } else { l_ins };
    l = if l > 1 { l } else { 1 };
    if l < (opt.w << 1) {
        l
    } else {
        opt.w << 1
    }
}

#[doc = "Original function: smem_aux_init:80"]
pub fn smem_aux_init() -> Vec<smem_aux_t> {
    let mut a = Vec::with_capacity(BATCH_SIZE);
    for _ in 0..BATCH_SIZE {
        a.push(smem_aux_t::default());
    }
    a
}

#[doc = "Original function: smem_aux_destroy:93"]
pub fn smem_aux_destroy(a: &mut Vec<smem_aux_t>) {
    for item in a.iter_mut() {
        item.tmpv[0].a.clear();
        item.tmpv[1].a.clear();
        item.mem.a.clear();
        item.mem1.a.clear();
    }
    a.clear();
}

#[doc = "Original function: mem_opt_init:107"]
pub fn mem_opt_init() -> Box<mem_opt_t> {
    let mut o = Box::new(mem_opt_t::default());
    o.flag = 0;
    o.a = 1;
    o.b = 4;
    o.o_del = 6;
    o.o_ins = 6;
    o.e_del = 1;
    o.e_ins = 1;
    o.w = 100;
    o.T = 30;
    o.zdrop = 100;
    o.pen_unpaired = 17;
    o.pen_clip5 = 5;
    o.pen_clip3 = 5;
    o.max_mem_intv = 20;
    o.min_seed_len = 19;
    o.split_width = 10;
    o.max_occ = 500;
    o.max_chain_gap = 10000;
    o.max_ins = 10000;
    o.mask_level = 0.50;
    o.drop_ratio = 0.50;
    o.XA_drop_ratio = 0.80;
    o.split_factor = 1.5;
    o.chunk_size = 10_000_000;
    o.n_threads = 1;
    o.max_XA_hits = 5;
    o.max_XA_hits_alt = 200;
    o.max_matesw = 50;
    o.mask_level_redun = 0.95;
    o.min_chain_weight = 0;
    o.max_chain_extend = 1 << 30;
    o.mapQ_coef_len = 50.0;
    o.mapQ_coef_fac = o.mapQ_coef_len.ln() as i32;
    bwa_fill_scmat(o.a, o.b, &mut o.mat);
    o
}

#[doc = "Original function: sort_alnreg_re:162"]
pub fn sort_alnreg_re(_n: i32, a: &mut [mem_alnreg_t]) {
    ks_introsort_mem_ars2(a);
}

#[doc = "Original function: sort_alnreg_score:166"]
pub fn sort_alnreg_score(_n: i32, a: &mut [mem_alnreg_t]) {
    ks_introsort_mem_ars(a);
}

#[inline(always)]
fn mem_ars2_lt(lhs: &mem_alnreg_t, rhs: &mem_alnreg_t) -> bool {
    lhs.re < rhs.re
}

#[inline(always)]
fn mem_ars_lt(lhs: &mem_alnreg_t, rhs: &mem_alnreg_t) -> bool {
    lhs.score > rhs.score
        || (lhs.score == rhs.score && (lhs.rb < rhs.rb || (lhs.rb == rhs.rb && lhs.qb < rhs.qb)))
}

fn mem_ars2_insertion_sort(a: &mut [mem_alnreg_t]) {
    for i in 1..a.len() {
        let mut j = i;
        while j > 0 && mem_ars2_lt(&a[j], &a[j - 1]) {
            a.swap(j, j - 1);
            j -= 1;
        }
    }
}

fn mem_ars_insertion_sort(a: &mut [mem_alnreg_t]) {
    for i in 1..a.len() {
        let mut j = i;
        while j > 0 && mem_ars_lt(&a[j], &a[j - 1]) {
            a.swap(j, j - 1);
            j -= 1;
        }
    }
}

fn mem_ars2_combsort(a: &mut [mem_alnreg_t]) {
    const SHRINK_FACTOR: f64 = 1.2473309501039787;
    let mut gap = a.len();
    let mut do_swap = true;
    while do_swap || gap > 2 {
        if gap > 2 {
            gap = (gap as f64 / SHRINK_FACTOR) as usize;
            if gap == 9 || gap == 10 {
                gap = 11;
            }
        }
        do_swap = false;
        for i in 0..a.len().saturating_sub(gap) {
            let j = i + gap;
            if mem_ars2_lt(&a[j], &a[i]) {
                a.swap(i, j);
                do_swap = true;
            }
        }
    }
    if gap != 1 {
        mem_ars2_insertion_sort(a);
    }
}

fn mem_ars_combsort(a: &mut [mem_alnreg_t]) {
    const SHRINK_FACTOR: f64 = 1.2473309501039787;
    let mut gap = a.len();
    let mut do_swap = true;
    while do_swap || gap > 2 {
        if gap > 2 {
            gap = (gap as f64 / SHRINK_FACTOR) as usize;
            if gap == 9 || gap == 10 {
                gap = 11;
            }
        }
        do_swap = false;
        for i in 0..a.len().saturating_sub(gap) {
            let j = i + gap;
            if mem_ars_lt(&a[j], &a[i]) {
                a.swap(i, j);
                do_swap = true;
            }
        }
    }
    if gap != 1 {
        mem_ars_insertion_sort(a);
    }
}

fn ks_introsort_mem_ars2(a: &mut [mem_alnreg_t]) {
    if a.len() < 2 {
        return;
    }
    if a.len() == 2 {
        if mem_ars2_lt(&a[1], &a[0]) {
            a.swap(0, 1);
        }
        return;
    }

    let mut d = 2usize;
    while (1usize << d) < a.len() {
        d += 1;
    }
    d <<= 1;

    let mut stack = [(0_usize, 0_usize, 0_usize); 128];
    let mut stack_len = 0_usize;
    let mut s = 0usize;
    let mut t = a.len() - 1;

    loop {
        if s < t {
            d -= 1;
            if d == 0 {
                mem_ars2_combsort(&mut a[s..=t]);
                t = s;
                continue;
            }

            let mut i = s;
            let mut j = t;
            let mut k = i + ((j - i) >> 1) + 1;
            if mem_ars2_lt(&a[k], &a[i]) {
                if mem_ars2_lt(&a[k], &a[j]) {
                    k = j;
                }
            } else {
                k = if mem_ars2_lt(&a[j], &a[i]) { i } else { j };
            }
            let rp = a[k].clone();
            if k != t {
                a.swap(k, t);
            }

            loop {
                loop {
                    i += 1;
                    if !mem_ars2_lt(&a[i], &rp) {
                        break;
                    }
                }
                loop {
                    j -= 1;
                    if i > j || !mem_ars2_lt(&rp, &a[j]) {
                        break;
                    }
                }
                if j <= i {
                    break;
                }
                a.swap(i, j);
            }
            a.swap(i, t);

            if i - s > t - i {
                if i - s > 16 {
                    debug_assert!(stack_len < stack.len());
                    stack[stack_len] = (s, i - 1, d);
                    stack_len += 1;
                }
                s = if t - i > 16 { i + 1 } else { t };
            } else {
                if t - i > 16 {
                    debug_assert!(stack_len < stack.len());
                    stack[stack_len] = (i + 1, t, d);
                    stack_len += 1;
                }
                t = if i - s > 16 { i - 1 } else { s };
            }
        } else if stack_len != 0 {
            stack_len -= 1;
            let (left, right, depth) = stack[stack_len];
            s = left;
            t = right;
            d = depth;
        } else {
            mem_ars2_insertion_sort(a);
            return;
        }
    }
}

fn ks_introsort_mem_ars(a: &mut [mem_alnreg_t]) {
    if a.len() < 2 {
        return;
    }
    if a.len() == 2 {
        if mem_ars_lt(&a[1], &a[0]) {
            a.swap(0, 1);
        }
        return;
    }

    let mut d = 2usize;
    while (1usize << d) < a.len() {
        d += 1;
    }
    d <<= 1;

    let mut stack = [(0_usize, 0_usize, 0_usize); 128];
    let mut stack_len = 0_usize;
    let mut s = 0usize;
    let mut t = a.len() - 1;

    loop {
        if s < t {
            d -= 1;
            if d == 0 {
                mem_ars_combsort(&mut a[s..=t]);
                t = s;
                continue;
            }

            let mut i = s;
            let mut j = t;
            let mut k = i + ((j - i) >> 1) + 1;
            if mem_ars_lt(&a[k], &a[i]) {
                if mem_ars_lt(&a[k], &a[j]) {
                    k = j;
                }
            } else {
                k = if mem_ars_lt(&a[j], &a[i]) { i } else { j };
            }
            let rp = a[k].clone();
            if k != t {
                a.swap(k, t);
            }

            loop {
                loop {
                    i += 1;
                    if !mem_ars_lt(&a[i], &rp) {
                        break;
                    }
                }
                loop {
                    j -= 1;
                    if i > j || !mem_ars_lt(&rp, &a[j]) {
                        break;
                    }
                }
                if j <= i {
                    break;
                }
                a.swap(i, j);
            }
            a.swap(i, t);

            if i - s > t - i {
                if i - s > 16 {
                    debug_assert!(stack_len < stack.len());
                    stack[stack_len] = (s, i - 1, d);
                    stack_len += 1;
                }
                s = if t - i > 16 { i + 1 } else { t };
            } else {
                if t - i > 16 {
                    debug_assert!(stack_len < stack.len());
                    stack[stack_len] = (i + 1, t, d);
                    stack_len += 1;
                }
                t = if i - s > 16 { i - 1 } else { s };
            }
        } else if stack_len != 0 {
            stack_len -= 1;
            let (left, right, depth) = stack[stack_len];
            s = left;
            t = right;
            d = depth;
        } else {
            mem_ars_insertion_sort(a);
            return;
        }
    }
}

#[doc = "Original function: mem_patch_reg:175"]
pub fn mem_patch_reg(
    opt: &mem_opt_t,
    bns: Option<&bntseq_t>,
    pac: Option<&[u8]>,
    query: Option<&mut [u8]>,
    a: &mem_alnreg_t,
    b: &mem_alnreg_t,
    w_out: &mut i32,
) -> i32 {
    let (Some(bns), Some(pac), Some(query)) = (bns, pac, query) else {
        return 0;
    };
    assert!(a.rid == b.rid && a.rb <= b.rb);
    if a.rb < bns.l_pac && b.rb >= bns.l_pac {
        return 0;
    }
    if a.qb >= b.qb || a.qe >= b.qe || a.re >= b.re {
        return 0;
    }

    let mut w = ((a.re - b.rb) - i64::from(a.qe - b.qb)) as i32;
    w = w.abs();
    let mut r = (a.re - b.rb) as f64 / (b.re - a.rb) as f64
        - f64::from(a.qe - b.qb) / f64::from(b.qe - a.qb);
    r = r.abs();

    if a.re < b.rb || a.qe < b.qb {
        if w > (opt.w << 1) || r >= PATCH_MAX_R_BW {
            return 0;
        }
    } else if w > (opt.w << 2) || r >= PATCH_MAX_R_BW * 2.0 {
        return 0;
    }

    w += a.w + b.w;
    w = w.min(opt.w << 2);

    let mut score = 0;
    let query_beg = a.qb as usize;
    let query_end = b.qe as usize;
    let _ = bwa_gen_cigar2(
        &opt.mat,
        opt.o_del,
        opt.e_del,
        opt.o_ins,
        opt.e_ins,
        w,
        bns.l_pac,
        pac,
        b.qe - a.qb,
        &mut query[query_beg..query_end],
        a.rb,
        b.re,
        &mut score,
        None,
        None,
    );

    let q_s = (((b.qe - a.qb) as f64 / f64::from((b.qe - b.qb) + (a.qe - a.qb)))
        * f64::from(b.score + a.score)
        + 0.499) as i32;
    let r_s = ((((b.re - a.rb) as f64) / ((b.re - b.rb) + (a.re - a.rb)) as f64)
        * f64::from(b.score + a.score)
        + 0.499) as i32;
    if score as f64 / f64::from(q_s.max(r_s)) < PATCH_MIN_SC_RATIO {
        return 0;
    }
    *w_out = w;
    score
}

#[doc = "Original function: mem_dedup_patch:239"]
pub fn mem_dedup_patch(
    opt: &mem_opt_t,
    bns: Option<&bntseq_t>,
    pac: Option<&[u8]>,
    query: Option<&mut [u8]>,
    mut n: i32,
    a: &mut Vec<mem_alnreg_t>,
) -> i32 {
    if n <= 1 {
        return n;
    }
    let n_us = n as usize;
    for reg in a.iter_mut().take(n_us) {
        reg.n_comp = 1;
    }

    let mut query = query;
    for i in 1..n_us {
        // Read p fields fresh: C++ uses a pointer `p = &a[i]` so post-merge updates to
        // `a[i].rb`/`a[i].qb` are visible inside the inner j-loop and widen its reach left.
        if a[i].rid != a[i - 1].rid || a[i].rb >= a[i - 1].re + i64::from(opt.max_chain_gap) {
            continue;
        }
        let mut j = i;
        while j > 0 {
            j -= 1;
            if a[i].rid != a[j].rid || a[i].rb >= a[j].re + i64::from(opt.max_chain_gap) {
                break;
            }
            if a[j].qe == a[j].qb {
                continue;
            }
            let p_qb = a[i].qb;
            let p_qe = a[i].qe;
            let p_rb = a[i].rb;
            let p_re = a[i].re;
            let q = a[j];
            let or_ = q.re - p_rb;
            let oq = if q.qb < p_qb {
                q.qe - p_qb
            } else {
                p_qe - q.qb
            };
            let mr = (q.re - q.rb).min(p_re - p_rb);
            let mq = i64::from((q.qe - q.qb).min(p_qe - p_qb));
            if (or_ as f32) > opt.mask_level_redun * mr as f32
                && (oq as f32) > opt.mask_level_redun * mq as f32
            {
                if a[i].score < a[j].score {
                    a[i].qe = a[i].qb;
                    break;
                } else {
                    a[j].qe = a[j].qb;
                }
            } else if a[j].rb < a[i].rb {
                let mut w = 0;
                let score =
                    mem_patch_reg(opt, bns, pac, query.as_deref_mut(), &a[j], &a[i], &mut w);
                if score > 0 {
                    let q_n_comp = a[j].n_comp;
                    let q_seedcov = a[j].seedcov;
                    let q_sub = a[j].sub;
                    let q_csub = a[j].csub;
                    let q_qb = a[j].qb;
                    let q_rb = a[j].rb;
                    a[i].n_comp += q_n_comp + 1;
                    a[i].seedcov = a[i].seedcov.max(q_seedcov);
                    a[i].sub = a[i].sub.max(q_sub);
                    a[i].csub = a[i].csub.max(q_csub);
                    a[i].qb = q_qb;
                    a[i].rb = q_rb;
                    a[i].truesc = score;
                    a[i].score = score;
                    a[i].w = w;
                    a[j].qb = a[j].qe;
                }
            }
        }
    }

    let mut m = 0_usize;
    for i in 0..n_us {
        if a[i].qe > a[i].qb {
            if m != i {
                a[m] = a[i];
            }
            m += 1;
        }
    }
    a.truncate(m);
    n = m as i32;
    n
}

#[doc = "Original function: mem_sort_dedup_patch:292"]
pub fn mem_sort_dedup_patch(
    opt: &mem_opt_t,
    bns: Option<&bntseq_t>,
    pac: Option<&[u8]>,
    query: Option<&mut [u8]>,
    n: i32,
    a: &mut Vec<mem_alnreg_t>,
) -> i32 {
    if n <= 1 {
        return n;
    }
    sort_alnreg_re(n, a.as_mut_slice());
    let n = mem_dedup_patch(opt, bns, pac, query, n, a);
    sort_alnreg_score(n, a.as_mut_slice());
    let n_us = n as usize;
    for i in 1..n_us {
        if a[i].score == a[i - 1].score && a[i].rb == a[i - 1].rb && a[i].qb == a[i - 1].qb {
            a[i].qe = a[i].qb;
        }
    }
    let mut m = if n > 0 { 1_usize } else { 0 };
    for i in 1..n_us {
        if a[i].qe > a[i].qb {
            if m != i {
                a[m] = a[i];
            }
            m += 1;
        }
    }
    a.truncate(m);
    m as i32
}

#[doc = "Original function: test_and_merge:357"]
#[inline]
pub fn test_and_merge(
    opt: &mem_opt_t,
    l_pac: i64,
    c: &mut mem_chain_t,
    p: &mem_seed_t,
    seed_rid: i32,
    _tid: i32,
) -> i32 {
    let last = &c.seeds[(c.n - 1) as usize];
    let qend = last.qbeg + last.len;
    let rend = last.rbeg + i64::from(last.len);

    if seed_rid != c.rid {
        return 0;
    }
    if p.qbeg >= c.seeds[0].qbeg
        && p.qbeg + p.len <= qend
        && p.rbeg >= c.seeds[0].rbeg
        && p.rbeg + i64::from(p.len) <= rend
    {
        return 1;
    }

    if (last.rbeg < l_pac || c.seeds[0].rbeg < l_pac) && p.rbeg >= l_pac {
        return 0;
    }

    let x = p.qbeg - last.qbeg;
    let y = p.rbeg - last.rbeg;
    if y >= 0
        && i64::from(x) - y <= i64::from(opt.w)
        && y - i64::from(x) <= i64::from(opt.w)
        && x - last.len < opt.max_chain_gap
        && y - i64::from(last.len) < i64::from(opt.max_chain_gap)
    {
        if c.n == c.m {
            c.m <<= 1;
            let want = c.m as usize;
            if c.seeds.capacity() < want {
                c.seeds.reserve(want - c.seeds.capacity());
            }
        }
        let used = c.n as usize;
        if used == c.seeds.len() {
            c.seeds.push(*p);
        } else {
            c.seeds[used] = *p;
        }
        c.n += 1;
        return 1;
    }
    0
}

#[doc = "Original function: mem_seed_sw:401"]
#[inline]
pub fn mem_seed_sw(
    opt: &mem_opt_t,
    bns: &bntseq_t,
    pac: &[u8],
    l_query: i32,
    query: &[u8],
    s: &mem_seed_t,
) -> i32 {
    let l_pac = bns.l_pac;
    if s.len >= MEM_SHORT_LEN {
        return -1;
    }
    let mut qb = s.qbeg;
    let mut qe = s.qbeg + s.len;
    let mut rb = s.rbeg;
    let mut re = s.rbeg + i64::from(s.len);
    let mid = (rb + re) >> 1;
    qb = (qb - MEM_SHORT_EXT).max(0);
    qe = (qe + MEM_SHORT_EXT).min(l_query);
    rb = (rb - i64::from(MEM_SHORT_EXT)).max(0);
    re = (re + i64::from(MEM_SHORT_EXT)).min(l_pac << 1);
    if rb < l_pac && l_pac < re {
        if mid < l_pac {
            re = l_pac;
        } else {
            rb = l_pac;
        }
    }
    if qe - qb >= MEM_SHORT_LEN || re - rb >= i64::from(MEM_SHORT_LEN) {
        return -1;
    }

    let mut rid = -1;
    let mut fetch_rb = rb;
    let mut fetch_re = re;
    MEM_SEED_SW_SCRATCH.with(|cell| {
        let mut buf = cell.borrow_mut();
        let (rseq, qseq) = &mut *buf;
        crate::generated::bntseq_cpp::bns_fetch_seq_into(
            bns,
            pac,
            &mut fetch_rb,
            mid,
            &mut fetch_re,
            &mut rid,
            rseq,
        );
        qseq.clear();
        // qb/qe are i32 ≥ 0 (clamped via .max(0)/.min(l_query) above); fits in usize.
        qseq.extend_from_slice(&query[qb as usize..qe as usize]);
        let rseq_len = rseq.len() as i32;
        let x = ksw_align2(
            qe - qb,
            qseq,
            rseq_len,
            rseq,
            5,
            &opt.mat,
            opt.o_del,
            opt.e_del,
            opt.o_ins,
            opt.e_ins,
            KSW_XSTART,
            None,
        );
        x.score
    })
}

// Thread-local scratch for mem_seed_sw rseq + qseq (called ~1.25M times per 50K-read run).
thread_local! {
    static MEM_SEED_SW_SCRATCH: std::cell::RefCell<(Vec<u8>, Vec<u8>)> =
        const { std::cell::RefCell::new((Vec::new(), Vec::new())) };
}

// Thread-local scratch for the per-chunk `aln` Vec in worker_sam (PE path). Sized to pcnt+256
// (kswr_t = 28 bytes). Avoids zeroing 30-60KB on every chunk; the Vec grows and is reused across
// chunks.
thread_local! {
    static WORKER_SAM_ALN_SCRATCH: std::cell::RefCell<Vec<crate::generated::ksw_h::kswr_t>> =
        const { std::cell::RefCell::new(Vec::new()) };
}

// Thread-local scratch for mem_mark_primary_se. `z` is the working set of primary candidates; `rank`
// is the post-sort permutation. Both are sized by per-read n_alnreg (~3-10 typically).
thread_local! {
    static MEM_MARK_PRIMARY_SE_SCRATCH: std::cell::RefCell<(Vec<usize>, Vec<usize>)> =
        const { std::cell::RefCell::new((Vec::new(), Vec::new())) };
}

// Thread-local scratch for mem_reg2sam: `aa` is the per-read accumulator of mem_aln_t's that get
// emitted as SAM. Pooling the outer Vec saves the per-read alloc; the inner mem_aln_t's still
// allocate via mem_reg2aln.
thread_local! {
    static MEM_REG2SAM_AA_SCRATCH: std::cell::RefCell<Vec<crate::generated::bwamem_h::mem_aln_t>> =
        const { std::cell::RefCell::new(Vec::new()) };
}

#[doc = "Original function: mem_chain_weight:429"]
#[inline]
pub fn mem_chain_weight(c: &mem_chain_t) -> i32 {
    let n_usize = c.n as usize; // c.n is i32 ≥ 0, fits in usize.
    let mut end = 0_i64;
    let mut w = 0_i64;
    for s in c.seeds.iter().take(n_usize) {
        if i64::from(s.qbeg) >= end {
            w += i64::from(s.len);
        } else if i64::from(s.qbeg + s.len) > end {
            w += i64::from(s.qbeg + s.len) - end;
        }
        end = end.max(i64::from(s.qbeg + s.len));
    }
    let tmp = w;
    end = 0;
    w = 0;
    for s in c.seeds.iter().take(n_usize) {
        if s.rbeg >= end {
            w += i64::from(s.len);
        } else if s.rbeg + i64::from(s.len) > end {
            w += s.rbeg + i64::from(s.len) - end;
        }
        end = end.max(s.rbeg + i64::from(s.len));
    }
    let w = w.min(tmp);
    if w < (1_i64 << 30) {
        w as i32
    } else {
        (1 << 30) - 1
    }
}

fn mem_flt_lt(lhs: &mem_chain_t, rhs: &mem_chain_t) -> bool {
    lhs.w > rhs.w
}

fn mem_flt_insertion_sort(a: &mut [mem_chain_t]) {
    for i in 1..a.len() {
        let mut j = i;
        while j > 0 && mem_flt_lt(&a[j], &a[j - 1]) {
            a.swap(j, j - 1);
            j -= 1;
        }
    }
}

fn mem_flt_combsort(a: &mut [mem_chain_t]) {
    const SHRINK_FACTOR: f64 = 1.2473309501039787;
    let mut gap = a.len();
    let mut do_swap = true;
    while do_swap || gap > 2 {
        if gap > 2 {
            gap = (gap as f64 / SHRINK_FACTOR) as usize;
            if gap == 9 || gap == 10 {
                gap = 11;
            }
        }
        do_swap = false;
        for i in 0..a.len().saturating_sub(gap) {
            let j = i + gap;
            if mem_flt_lt(&a[j], &a[i]) {
                a.swap(i, j);
                do_swap = true;
            }
        }
    }
    if gap != 1 {
        mem_flt_insertion_sort(a);
    }
}

fn ks_introsort_mem_flt(a: &mut [mem_chain_t]) {
    if a.len() < 2 {
        return;
    }
    if a.len() == 2 {
        if mem_flt_lt(&a[1], &a[0]) {
            a.swap(0, 1);
        }
        return;
    }

    let mut d = 2usize;
    while (1usize << d) < a.len() {
        d += 1;
    }
    d <<= 1;

    let mut stack = [(0_usize, 0_usize, 0_usize); 128];
    let mut stack_len = 0_usize;
    let mut s = 0usize;
    let mut t = a.len() - 1;

    loop {
        if s < t {
            d -= 1;
            if d == 0 {
                mem_flt_combsort(&mut a[s..=t]);
                t = s;
                continue;
            }

            let mut i = s;
            let mut j = t;
            let mut k = i + ((j - i) >> 1) + 1;
            if mem_flt_lt(&a[k], &a[i]) {
                if mem_flt_lt(&a[k], &a[j]) {
                    k = j;
                }
            } else {
                k = if mem_flt_lt(&a[j], &a[i]) { i } else { j };
            }
            let rp = a[k].clone();
            if k != t {
                a.swap(k, t);
            }

            loop {
                loop {
                    i += 1;
                    if !mem_flt_lt(&a[i], &rp) {
                        break;
                    }
                }
                loop {
                    j -= 1;
                    if i > j || !mem_flt_lt(&rp, &a[j]) {
                        break;
                    }
                }
                if j <= i {
                    break;
                }
                a.swap(i, j);
            }
            a.swap(i, t);

            if i - s > t - i {
                if i - s > 16 {
                    debug_assert!(stack_len < stack.len());
                    stack[stack_len] = (s, i - 1, d);
                    stack_len += 1;
                }
                s = if t - i > 16 { i + 1 } else { t };
            } else {
                if t - i > 16 {
                    debug_assert!(stack_len < stack.len());
                    stack[stack_len] = (i + 1, t, d);
                    stack_len += 1;
                }
                t = if i - s > 16 { i - 1 } else { s };
            }
        } else if stack_len != 0 {
            stack_len -= 1;
            let (left, right, depth) = stack[stack_len];
            s = left;
            t = right;
            d = depth;
        } else {
            mem_flt_insertion_sort(a);
            return;
        }
    }
}

#[doc = "Original function: mem_print_chain:450"]
pub fn mem_print_chain(bns: &bntseq_t, chn: &[mem_chain_t]) -> String {
    let mut out = String::new();
    for (i, p) in chn.iter().enumerate() {
        out.push_str(&format!(
            "* Found CHAIN({i}): n={}; weight={}",
            p.n,
            mem_chain_weight(p)
        ));
        for seed in p.seeds.iter().take(usize::try_from(p.n).expect("p.n")) {
            let mut is_rev = 0;
            let mut pos = bns_depos(bns, seed.rbeg, &mut is_rev);
            if is_rev != 0 {
                pos -= i64::from(seed.len - 1);
            }
            let ann = &bns.anns[usize::try_from(p.rid).expect("rid")];
            let strand = if is_rev != 0 { '-' } else { '+' };
            out.push_str(&format!(
                "\t{};{};{},{}({}:{strand}{})",
                seed.score,
                seed.len,
                seed.qbeg,
                seed.rbeg,
                ann.name,
                pos - ann.offset + 1
            ));
        }
        out.push('\n');
    }
    out
}

// Thread-local kept-seeds buffer for mem_flt_chained_seeds (called ~250K times per 50K-read run).
thread_local! {
    static MEM_FLT_CHAINED_KEPT: std::cell::RefCell<Vec<crate::generated::bwamem_h::mem_seed_t>> =
        const { std::cell::RefCell::new(Vec::new()) };

    // Reused across mem_chain_flt calls (~50K per 50K-read run). filtered/ranges/result/chains/
    // group are all small (chains-per-read typically <20) but per-call allocs add up.
    static MEM_CHAIN_FLT_SCRATCH: std::cell::RefCell<(
        Vec<mem_chain_t>,
        Vec<(usize, usize)>,
        Vec<mem_chain_t>,
        Vec<usize>,
        Vec<mem_chain_t>,
    )> = std::cell::RefCell::new((Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new()));
}

#[doc = "Original function: mem_flt_chained_seeds:472"]
pub fn mem_flt_chained_seeds(
    opt: &mem_opt_t,
    bns: &bntseq_t,
    pac: &[u8],
    seq_: &[bseq1_t],
    n_chn: i32,
    a: &mut [mem_chain_t],
) {
    for c in a.iter_mut().take(n_chn as usize) {
        let seqid = c.seqid as usize;
        let query_cow = nt4_cow(&seq_[seqid]);
        let query: &[u8] = &query_cow;
        let l_query = seq_[seqid].l_seq;

        let min_l = if opt.min_chain_weight != 0 {
            MEM_HSP_COEF * f64::from(opt.min_chain_weight)
        } else {
            MEM_MINSC_COEF * f64::from(l_query).ln()
        };
        let min_hsp_score = (f64::from(opt.a) * min_l + 0.499) as i32;
        if min_l > MEM_SEEDSW_COEF * f64::from(l_query) {
            continue;
        }

        let c_n_usize = c.n as usize; // c.n is i32, fits in usize for typical reads.
        MEM_FLT_CHAINED_KEPT.with(|cell| {
            let mut kept = cell.borrow_mut();
            kept.clear();
            kept.reserve(c_n_usize);
            for mut s in c.seeds.iter().take(c_n_usize).copied() {
                s.score = mem_seed_sw(opt, bns, pac, l_query, query, &s);
                if s.score < 0 || s.score >= min_hsp_score {
                    s.score = if s.score < 0 { s.len * opt.a } else { s.score };
                    kept.push(s);
                }
            }
            c.n = kept.len() as i32;
            for (dst, src) in c.seeds.iter_mut().zip(kept.iter().copied()) {
                *dst = src;
            }
        });
    }
}

#[doc = "Original function: mem_chain_flt:506"]
pub fn mem_chain_flt(opt: &mem_opt_t, n_chn_: i32, a_: &mut Vec<mem_chain_t>, _tid: i32) -> i32 {
    if n_chn_ == 0 {
        return 0;
    }

    MEM_CHAIN_FLT_SCRATCH.with(|cell| {
        let mut buf = cell.borrow_mut();
        let (filtered, ranges, result, chains, group) = &mut *buf;
        filtered.clear();
        ranges.clear();
        result.clear();

        let take_n = (n_chn_ as usize).min(a_.len());
        for mut c in a_.drain(..take_n) {
            c.first = -1;
            c.kept = 0;
            c.w = mem_chain_weight(&c) as u32;
            if (c.w as i32) >= opt.min_chain_weight {
                filtered.push(c);
            } else {
                // Chain is being dropped — pool its seeds Vec for reuse.
                return_seeds_buf(c.seeds);
            }
        }
        a_.clear();
        if filtered.is_empty() {
            return 0;
        }

        let mut start = 0_usize;
        let mut pseqid = filtered[0].seqid;
        for i in 1..filtered.len() {
            if filtered[i].seqid != pseqid {
                ranges.push((start, i));
                start = i;
            }
            pseqid = filtered[i].seqid;
        }
        ranges.push((start, filtered.len()));

        for &(start, end) in ranges.iter() {
            // Move chains out of filtered into group (avoids cloning seeds Vec).
            // Filtered slots become default mem_chain_t (empty seeds Vec, no alloc).
            group.clear();
            group.extend(
                filtered[start..end]
                    .iter_mut()
                    .map(|c| std::mem::take(c)),
            );
            ks_introsort_mem_flt(group);

            chains.clear();
            chains.push(0_usize);
            group[0].kept = 3;
            for i in 1..group.len() {
                let mut large_ovlp = false;
                let mut admitted = true;
                for &j in chains.iter() {
                    let b_max = chn_beg(&group[j]).max(chn_beg(&group[i]));
                    let e_min = chn_end(&group[j]).min(chn_end(&group[i]));
                    if e_min > b_max && (group[j].is_alt == 0 || group[i].is_alt != 0) {
                        let li = chn_end(&group[i]) - chn_beg(&group[i]);
                        let lj = chn_end(&group[j]) - chn_beg(&group[j]);
                        let min_l = li.min(lj);
                        if ((e_min - b_max) as f32) >= (min_l as f32) * opt.mask_level
                            && min_l < opt.max_chain_gap
                        {
                            large_ovlp = true;
                            if group[j].first < 0 {
                                group[j].first = i as i32;
                            }
                            if (group[i].w as f32) < (group[j].w as f32) * opt.drop_ratio
                                && (group[j].w - group[i].w) as i32 >= (opt.min_seed_len << 1)
                            {
                                admitted = false;
                                break;
                            }
                        }
                    }
                }
                if admitted {
                    chains.push(i);
                    group[i].kept = if large_ovlp { 2 } else { 3 };
                }
            }

            for &idx in chains.iter() {
                let first = group[idx].first;
                if first >= 0 {
                    group[first as usize].kept = 1;
                }
            }

            // Match C++ bwamem.cpp:597-603: when max_chain_extend triggers a break, the
            // index of the breaking chain is the start of the zero-out loop (NOT i+1).
            let mut limited = 0_i32;
            let mut i = 0_usize;
            let mut hit_limit = false;
            while i < group.len() {
                if group[i].kept != 0 && group[i].kept != 3 {
                    limited += 1;
                    if limited >= opt.max_chain_extend {
                        hit_limit = true;
                        break;
                    }
                }
                i += 1;
            }
            let tail_start = if hit_limit { i } else { group.len() };
            for c in group.iter_mut().skip(tail_start) {
                if c.kept < 3 {
                    c.kept = 0;
                }
            }

            // Drain group; kept chains move into result, dropped chains return their seeds
            // Vec to the pool for reuse instead of being freed.
            for c in group.drain(..) {
                if c.kept != 0 {
                    result.push(c);
                } else {
                    return_seeds_buf(c.seeds);
                }
            }
        }

        a_.append(result);
        a_.len() as i32
    })
}

#[doc = "Original function: kv_push:564"]
pub fn kv_push__L564(_arg0: crate::support::Opaque, _arg1: crate::support::Opaque) {
    crate::support::stub::<()>("kv_push")
}

#[doc = "Original function: mem_collect_smem:626"]
pub fn mem_collect_smem(
    fmi: &FMI_search,
    opt: &mem_opt_t,
    seq_: &[bseq1_t],
    nseq: i32,
    match_array: &mut Vec<SMEM>,
    min_intv_ar: &mut Vec<i32>,
    query_pos_ar: &mut Vec<i16>,
    enc_qdb: &mut Vec<u8>,
    rid: &mut Vec<i32>,
    mmc: &mut mem_cache,
    tot_smem: &mut i64,
    tid: i32,
) {
    let nseq_usize = nseq as usize;
    let split_len = (opt.min_seed_len as f32 * opt.split_factor + 0.499) as i32;
    let mut num_smem1 = 0_i64;
    let mut num_smem2 = 0_i64;
    let mut num_smem3 = 0_i64;
    let mut max_readlength = -1_i32;
    let mut query_cum_len_ar =
        MEM_COLLECT_QUERY_CUM_LEN.with(|c| std::mem::take(&mut *c.borrow_mut()));
    query_cum_len_ar.clear();
    query_cum_len_ar.resize(nseq_usize, 0);

    min_intv_ar.resize(nseq_usize.max(min_intv_ar.len()), 0);
    rid.resize(nseq_usize.max(rid.len()), 0);

    enc_qdb.clear();
    for (l, seq) in seq_.iter().take(nseq_usize).enumerate() {
        min_intv_ar[l] = 1;
        let nt4 = nt4_cow(seq);
        enc_qdb.extend_from_slice(&nt4[..seq.l_seq as usize]);
        rid[l] = l as i32;
    }

    if nseq_usize > 0 {
        max_readlength = seq_[0].l_seq;
        query_cum_len_ar[0] = 0;
    }
    for i in 1..nseq_usize {
        query_cum_len_ar[i] = query_cum_len_ar[i - 1] + seq_[i - 1].l_seq;
        if max_readlength < seq_[i].l_seq {
            max_readlength = seq_[i].l_seq;
        }
    }

    match_array.clear();
    fmi.getSMEMsAllPosOneThread(
        enc_qdb,
        &mut min_intv_ar[..nseq_usize],
        &mut rid[..nseq_usize],
        nseq,
        nseq,
        seq_,
        &query_cum_len_ar,
        max_readlength,
        opt.min_seed_len,
        match_array,
        &mut num_smem1,
    );

    let mut pos = 0_usize;
    let mut mem_lim = 0_i64;
    for i in 0..(num_smem1 as usize) {
        let p = match_array[i];
        let start = p.m as i32;
        let end = p.n as i32 + 1;
        if end - start < split_len || p.s > i64::from(opt.split_width) {
            continue;
        }
        let r = p.rid as usize;
        if rid.len() <= pos {
            rid.push(0);
        }
        if query_pos_ar.len() <= pos {
            query_pos_ar.push(0);
        }
        if min_intv_ar.len() <= pos {
            min_intv_ar.push(0);
        }
        rid[pos] = p.rid as i32;
        query_pos_ar[pos] = ((end + start) >> 1) as i16;
        assert!(i32::from(query_pos_ar[pos]) < seq_[r].l_seq);
        min_intv_ar[pos] = (p.s + 1) as i32;
        pos += 1;
        mem_lim += i64::from(end - start);
    }

    let tid_usize = tid as usize;
    if mmc.wsize_mem.len() <= tid_usize {
        mmc.wsize_mem.resize(tid_usize + 1, 0);
    }
    if mmc.wsize_mem[tid_usize] < mem_lim + num_smem1 {
        let want = (mem_lim + num_smem1).max(0) as usize;
        let extra = want.saturating_sub(match_array.capacity());
        match_array.reserve(extra);
        mmc.wsize_mem[tid_usize] = match_array.capacity() as i64;
    }

    let original_len = match_array.len();
    let pos_i32 = pos as i32;
    fmi.getSMEMsOnePosOneThread(
        enc_qdb,
        &mut query_pos_ar[..pos],
        &mut min_intv_ar[..pos],
        &mut rid[..pos],
        pos_i32,
        pos_i32,
        seq_,
        &query_cum_len_ar,
        max_readlength,
        opt.min_seed_len,
        match_array,
        &mut num_smem2,
    );
    num_smem2 = match_array.len().saturating_sub(original_len) as i64;

    if opt.max_mem_intv > 0 {
        let want =
            num_smem2 + (enc_qdb.len() as f64 / f64::from(opt.min_seed_len + 1)) as i64 + num_smem1;
        if mmc.wsize_mem[tid_usize] < want {
            let want = want.max(0) as usize;
            let extra = want.saturating_sub(match_array.capacity());
            match_array.reserve(extra);
            mmc.wsize_mem[tid_usize] = match_array.capacity() as i64;
        }
        for value in min_intv_ar.iter_mut().take(nseq_usize) {
            *value = opt.max_mem_intv as i32;
        }
        num_smem3 = fmi.bwtSeedStrategyAllPosOneThread(
            enc_qdb,
            &min_intv_ar[..nseq_usize],
            nseq,
            seq_,
            &query_cum_len_ar,
            opt.min_seed_len + 1,
            match_array,
        );
    }

    *tot_smem = num_smem1 + num_smem2 + num_smem3;
    fmi.sortSMEMs(match_array, &[*tot_smem], nseq, seq_[0].l_seq, 1);

    let tot_smem_usize = *tot_smem as usize;
    let mut smem_ptr = 0_usize;
    while smem_ptr < tot_smem_usize {
        let rid_value = match_array[smem_ptr].rid;
        let mut pos = smem_ptr;
        while pos + 1 < tot_smem_usize && match_array[pos + 1].rid == rid_value {
            pos += 1;
        }
        match_array[smem_ptr..=pos].sort_by_key(|s| (s.m, s.n));
        smem_ptr = pos + 1;
    }
    MEM_COLLECT_QUERY_CUM_LEN.with(|c| *c.borrow_mut() = query_cum_len_ar);
}

#[doc = "Original function: mem_chain_seeds:806"]
pub fn mem_chain_seeds(
    fmi: &FMI_search,
    opt: &mem_opt_t,
    bns: &bntseq_t,
    seq_: &[bseq1_t],
    nseq: i32,
    tid: i32,
    chain_ar: &mut [mem_chain_v],
    _seedBuf: &mut [mem_seed_t],
    _seedBufSize: i64,
    matchArray: &[SMEM],
    num_smem: i64,
) {
    // C++ uses kb_init(chn, KB_DEFAULT_SIZE + 8). With mem_chain_t's layout this
    // gives B-tree degree 5; duplicate-key lookup depends on the node split shape.
    const CHAIN_TREE_T: usize = 5;

    struct ChainTreeNode {
        is_internal: bool,
        key_len: usize,
        keys: [mem_chain_t; 2 * CHAIN_TREE_T - 1],
        child_len: usize,
        children: [usize; 2 * CHAIN_TREE_T],
    }

    impl ChainTreeNode {
        #[inline]
        fn new(is_internal: bool) -> Self {
            Self {
                is_internal,
                key_len: 0,
                keys: std::array::from_fn(|_| mem_chain_t::default()),
                child_len: 0,
                children: [0; 2 * CHAIN_TREE_T],
            }
        }

        #[inline]
        fn keys(&self) -> &[mem_chain_t] {
            &self.keys[..self.key_len]
        }

        #[inline]
        fn key_mut(&mut self, idx: usize) -> &mut mem_chain_t {
            &mut self.keys[idx]
        }

        #[inline]
        fn insert_key(&mut self, idx: usize, key: mem_chain_t) {
            debug_assert!(self.key_len < self.keys.len());
            debug_assert!(idx <= self.key_len);
            // Shift keys[idx..key_len] right by 1, leaving the default value (from keys[key_len])
            // at position idx, then overwrite. rotate_right is O(n) like the prior mem::take loop
            // but avoids the per-element take + assign overhead.
            self.keys[idx..=self.key_len].rotate_right(1);
            self.keys[idx] = key;
            self.key_len += 1;
        }

        #[inline]
        fn insert_child(&mut self, idx: usize, child: usize) {
            debug_assert!(self.child_len < self.children.len());
            debug_assert!(idx <= self.child_len);
            // children is [usize; ..] (Copy), copy_within is a single memmove.
            self.children.copy_within(idx..self.child_len, idx + 1);
            self.children[idx] = child;
            self.child_len += 1;
        }
    }

    struct ChainTree {
        nodes: Vec<ChainTreeNode>,
        root: usize,
        len: usize,
    }

    impl ChainTree {
        #[inline]
        fn new_node(is_internal: bool) -> ChainTreeNode {
            ChainTreeNode::new(is_internal)
        }

        fn new() -> Self {
            Self::with_capacity(16)
        }

        fn with_capacity(capacity: usize) -> Self {
            let mut nodes = Vec::with_capacity(capacity.max(1));
            nodes.push(Self::new_node(false));
            Self {
                nodes,
                root: 0,
                len: 0,
            }
        }

        fn reset(&mut self, capacity: usize) {
            self.nodes.clear();
            let want = capacity.max(1);
            if self.nodes.capacity() < want {
                self.nodes.reserve(want - self.nodes.capacity());
            }
            self.nodes.push(Self::new_node(false));
            self.root = 0;
            self.len = 0;
        }

        #[inline]
        fn len(&self) -> usize {
            self.len
        }

        #[inline]
        fn getp_aux(keys: &[mem_chain_t], pos: i64) -> (isize, i32) {
            if keys.is_empty() {
                return (-1, 0);
            }
            let mut begin = 0_usize;
            let mut end = keys.len();
            while begin < end {
                let mid = (begin + end) >> 1;
                if keys[mid].pos < pos {
                    begin = mid + 1;
                } else {
                    end = mid;
                }
            }
            if begin == keys.len() {
                return ((keys.len() - 1) as isize, 1);
            }
            let r = if pos < keys[begin].pos {
                -1
            } else if keys[begin].pos < pos {
                1
            } else {
                0
            };
            if r < 0 {
                (begin as isize - 1, r)
            } else {
                (begin as isize, r)
            }
        }

        #[inline]
        fn interval_index(&self, pos: i64) -> Option<(usize, usize)> {
            let mut x = self.root;
            let mut lower = None;
            loop {
                let node = &self.nodes[x];
                let (i, r) = Self::getp_aux(node.keys(), pos);
                if i >= 0 && r == 0 {
                    return Some((x, i as usize));
                }
                if i >= 0 {
                    lower = Some((x, i as usize));
                }
                if !node.is_internal {
                    return lower;
                }
                x = node.children[(i + 1) as usize];
            }
        }

        #[inline]
        fn interval_mut(&mut self, pos: i64) -> Option<&mut mem_chain_t> {
            let (node_idx, key_idx) = self.interval_index(pos)?;
            self.nodes
                .get_mut(node_idx)
                .map(|node| node.key_mut(key_idx))
        }

        fn split_child(&mut self, parent_idx: usize, child_slot: usize) {
            let child_idx = self.nodes[parent_idx].children[child_slot];
            let z_idx = self.nodes.len();
            let (median, z_node) = {
                let child = &mut self.nodes[child_idx];
                debug_assert_eq!(child.key_len, 2 * CHAIN_TREE_T - 1);
                let mut z_node = ChainTreeNode::new(child.is_internal);
                for idx in 0..(CHAIN_TREE_T - 1) {
                    z_node.keys[idx] = std::mem::take(&mut child.keys[CHAIN_TREE_T + idx]);
                }
                z_node.key_len = CHAIN_TREE_T - 1;
                let median = std::mem::take(&mut child.keys[CHAIN_TREE_T - 1]);
                child.key_len = CHAIN_TREE_T - 1;
                if child.is_internal {
                    for idx in 0..CHAIN_TREE_T {
                        z_node.children[idx] = child.children[CHAIN_TREE_T + idx];
                    }
                    z_node.child_len = CHAIN_TREE_T;
                    child.child_len = CHAIN_TREE_T;
                }
                (median, z_node)
            };
            self.nodes.push(z_node);
            let parent = &mut self.nodes[parent_idx];
            parent.insert_child(child_slot + 1, z_idx);
            parent.insert_key(child_slot, median);
        }

        #[inline]
        fn insert(&mut self, key: mem_chain_t) {
            if self.nodes[self.root].key_len == 2 * CHAIN_TREE_T - 1 {
                let old_root = self.root;
                let new_root = self.nodes.len();
                let mut new_node = Self::new_node(true);
                new_node.children[0] = old_root;
                new_node.child_len = 1;
                self.nodes.push(new_node);
                self.root = new_root;
                self.split_child(new_root, 0);
            }
            self.insert_nonfull(self.root, key);
            self.len += 1;
        }

        #[inline]
        fn insert_nonfull(&mut self, node_idx: usize, key: mem_chain_t) {
            if !self.nodes[node_idx].is_internal {
                let (i, _) = Self::getp_aux(self.nodes[node_idx].keys(), key.pos);
                self.nodes[node_idx].insert_key((i + 1) as usize, key);
                return;
            }

            let (i_raw, _) = Self::getp_aux(self.nodes[node_idx].keys(), key.pos);
            let mut child_slot = (i_raw + 1) as usize;
            let child_idx = self.nodes[node_idx].children[child_slot];
            if self.nodes[child_idx].key_len == 2 * CHAIN_TREE_T - 1 {
                self.split_child(node_idx, child_slot);
                if key.pos > self.nodes[node_idx].keys[child_slot].pos {
                    child_slot += 1;
                }
            }
            let next_idx = self.nodes[node_idx].children[child_slot];
            self.insert_nonfull(next_idx, key);
        }

        fn traverse_into(&mut self, mut out: Vec<mem_chain_t>) -> Vec<mem_chain_t> {
            fn visit(tree: &mut ChainTree, node_idx: usize, out: &mut Vec<mem_chain_t>) {
                let is_internal = tree.nodes[node_idx].is_internal;
                let key_count = tree.nodes[node_idx].key_len;
                let children = tree.nodes[node_idx].children;
                for i in 0..key_count {
                    if is_internal {
                        visit(tree, children[i], out);
                    }
                    out.push(std::mem::take(&mut tree.nodes[node_idx].keys[i]));
                }
                if is_internal {
                    visit(tree, children[key_count], out);
                }
            }

            out.clear();
            let extra = self.len().saturating_sub(out.capacity());
            out.reserve(extra);
            let root = self.root;
            visit(self, root, &mut out);
            out
        }
    }

    let mut pos = 0_i64;
    let mut smem_ptr = 0_usize;
    let l_pac = bns.l_pac;
    let mut smem_buf_size = 6000_i64;
    let mut sa_coord = MEM_CHAIN_SEEDS_SA_COORD.with(|c| std::mem::take(&mut *c.borrow_mut()));
    let initial_sa_cap = (opt.max_occ as usize) * (smem_buf_size as usize);
    if sa_coord.capacity() < initial_sa_cap {
        sa_coord.reserve(initial_sa_cap - sa_coord.capacity());
    }
    let mut chains = ChainTree::new();

    for chain in chain_ar
        .iter_mut()
        .take(nseq as usize)
    {
        clear_chain_a_pooled(&mut chain.a, KEEP_CHAIN_CAP);
        chain.n = 0;
        chain.m = 0;
        chain.cc = 0;
    }

    for l in 0..nseq as usize {
        if pos >= num_smem - 1 {
            break;
        }
        if smem_ptr >= matchArray.len() || (matchArray[smem_ptr].rid as usize) > l {
            continue;
        }
        if seq_[l].l_seq < opt.min_seed_len {
            continue;
        }
        assert_eq!(matchArray[smem_ptr].rid as usize, l);

        let chain = &mut chain_ar[l];
        let chain_out = std::mem::take(&mut chain.a);
        let mut b = 0_i32;
        let mut e = 0_i32;
        let mut l_rep = 0_i32;
        pos = smem_ptr as i64 - 1;
        loop {
            pos += 1;
            let pos_us = pos as usize; // pos >= 0, fits in usize.
            let p = &matchArray[pos_us];
            let sb = p.m as i32; // u32 → i32 ok for typical positions.
            let se = p.n as i32 + 1;
            if p.s > i64::from(opt.max_occ) {
                if sb > e {
                    l_rep += e - b;
                    b = sb;
                    e = se;
                } else {
                    e = e.max(se);
                }
            }
            if pos >= num_smem - 1
                || matchArray[pos_us].rid != matchArray[pos_us + 1].rid
            {
                break;
            }
        }
        l_rep += e - b;

        let smem_count = pos - smem_ptr as i64 + 1;
        if smem_count >= smem_buf_size {
            smem_buf_size *= 2;
        }
        sa_coord.clear();
        let mut cnt_ = 0_i64;
        let mut id = 0_i64;
        fmi.get_sa_entries_prefetch(
            &matchArray[smem_ptr..=pos as usize],
            &mut sa_coord,
            &mut cnt_,
            smem_count,
            opt.max_occ,
            tid,
            &mut id,
        );
        let node_cap = ((cnt_.max(1) as usize) / CHAIN_TREE_T).saturating_add(2);
        chains.reset(node_cap);

        let mut mypos = 0_usize;
        for p in &matchArray[smem_ptr..=pos as usize] {
            // p.m and p.n are u32, fit in i32 for typical positions.
            let p_m_i32 = p.m as i32;
            let p_n_i32 = p.n as i32;
            let slen = p_n_i32 + 1 - p_m_i32;
            let step = if p.s > i64::from(opt.max_occ) {
                p.s / i64::from(opt.max_occ)
            } else {
                1
            };

            let mut count = 0_i32;
            let mut k = 0_i64;
            while k < p.s && count < opt.max_occ {
                let rbeg = sa_coord[mypos];
                mypos += 1;
                let s = mem_seed_t {
                    rbeg,
                    qbeg: p_m_i32,
                    score: slen,
                    len: slen,
                    ..Default::default()
                };
                let mut tmp = mem_chain_t {
                    pos: rbeg,
                    ..Default::default()
                };
                let rid = bns_intv2rid(bns, s.rbeg, s.rbeg + i64::from(s.len));
                if rid >= 0 {
                    let mut to_add = true;
                    if chains.len() != 0 {
                        if let Some(lower) = chains.interval_mut(tmp.pos) {
                            if test_and_merge(opt, l_pac, lower, &s, rid, tid) != 0 {
                                to_add = false;
                            }
                        }
                    }
                    if to_add {
                        tmp.n = 1;
                        tmp.m = SEEDS_PER_CHAIN as i32;
                        tmp.seeds = take_seeds_buf();
                        tmp.seeds.push(s);
                        tmp.rid = rid;
                        tmp.seqid = l as i32;
                        tmp.is_alt = u32::from(bns.anns[rid as usize].is_alt != 0);
                        chains.insert(tmp);
                    }
                }
                k += step;
                count += 1;
            }
        }

        smem_ptr = (pos + 1) as usize;
        let mut chains = chains.traverse_into(chain_out);
        for c in &mut chains {
            c.frac_rep = l_rep as f32 / seq_[l].l_seq as f32;
        }
        chain.n = chains.len();
        chain.m = chains.len();
        chain.a = chains;
    }
    MEM_CHAIN_SEEDS_SA_COORD.with(|c| *c.borrow_mut() = sa_coord);
}

#[doc = "Original function: mem_kernel1_core:976"]
pub fn mem_kernel1_core(
    fmi: &FMI_search,
    opt: &mem_opt_t,
    seq_: &[bseq1_t],
    nseq: i32,
    chain_ar: &mut [mem_chain_v],
    seedBuf: &mut [mem_seed_t],
    seedBufSize: i64,
    mmc: &mut mem_cache,
    tid: i32,
) -> i32 {
    let tid_usize = tid as usize;
    if mmc.matchArray.len() <= tid_usize {
        mmc.matchArray.resize_with(tid_usize + 1, Vec::new);
    }
    if mmc.min_intv_ar.len() <= tid_usize {
        mmc.min_intv_ar.resize_with(tid_usize + 1, Vec::new);
    }
    if mmc.rid.len() <= tid_usize {
        mmc.rid.resize_with(tid_usize + 1, Vec::new);
    }
    if mmc.query_pos_ar.len() <= tid_usize {
        mmc.query_pos_ar.resize_with(tid_usize + 1, Vec::new);
    }
    if mmc.enc_qdb.len() <= tid_usize {
        mmc.enc_qdb.resize_with(tid_usize + 1, Vec::new);
    }
    if mmc.wsize_mem.len() <= tid_usize {
        mmc.wsize_mem.resize(tid_usize + 1, 0);
    }
    if mmc.wsize_mem_s.len() <= tid_usize {
        mmc.wsize_mem_s.resize(tid_usize + 1, 0);
    }
    if mmc.wsize_mem_r.len() <= tid_usize {
        mmc.wsize_mem_r.resize(tid_usize + 1, 0);
    }

    let mut tot_len = 0_i64;
    for seq in seq_.iter().take(nseq as usize) {
        tot_len += i64::from(seq.l_seq);
    }

    if tot_len >= mmc.wsize_mem[tid_usize] {
        mmc.wsize_mem[tid_usize] = tot_len;
        mmc.wsize_mem_s[tid_usize] = tot_len;
        mmc.wsize_mem_r[tid_usize] = tot_len;
        let cap = tot_len as usize;
        let extra = cap.saturating_sub(mmc.matchArray[tid_usize].capacity());
        mmc.matchArray[tid_usize].reserve(extra);
        let min_intv_extra = cap.saturating_sub(mmc.min_intv_ar[tid_usize].capacity());
        let query_pos_extra = cap.saturating_sub(mmc.query_pos_ar[tid_usize].capacity());
        let enc_qdb_extra = cap.saturating_sub(mmc.enc_qdb[tid_usize].capacity());
        let rid_extra = cap.saturating_sub(mmc.rid[tid_usize].capacity());
        mmc.min_intv_ar[tid_usize].reserve(min_intv_extra);
        mmc.query_pos_ar[tid_usize].reserve(query_pos_extra);
        mmc.enc_qdb[tid_usize].reserve(enc_qdb_extra);
        mmc.rid[tid_usize].reserve(rid_extra);
    }

    let mut num_smem = 0_i64;
    let mut matchArray = std::mem::take(&mut mmc.matchArray[tid_usize]);
    let mut min_intv_ar = std::mem::take(&mut mmc.min_intv_ar[tid_usize]);
    let mut query_pos_ar = std::mem::take(&mut mmc.query_pos_ar[tid_usize]);
    let mut enc_qdb = std::mem::take(&mut mmc.enc_qdb[tid_usize]);
    let mut rid = std::mem::take(&mut mmc.rid[tid_usize]);
    min_intv_ar.clear();
    query_pos_ar.clear();
    enc_qdb.clear();
    rid.clear();

    mem_collect_smem(
        fmi,
        opt,
        seq_,
        nseq,
        &mut matchArray,
        &mut min_intv_ar,
        &mut query_pos_ar,
        &mut enc_qdb,
        &mut rid,
        mmc,
        &mut num_smem,
        tid,
    );
    debug_rss("kernel1_after_collect_smem");

    mem_chain_seeds(
        fmi,
        opt,
        fmi.base.idx.bns.as_ref().expect("loaded bns"),
        seq_,
        nseq,
        tid,
        chain_ar,
        seedBuf,
        seedBufSize,
        &matchArray,
        num_smem,
    );
    debug_rss("kernel1_after_chain_seeds");

    for chn in chain_ar
        .iter_mut()
        .take(nseq as usize)
    {
        // chn.n is usize, fits in i32 for chain counts.
        chn.n = mem_chain_flt(opt, chn.n as i32, &mut chn.a, tid) as usize;
        chn.m = chn.a.len();
    }
    debug_rss("kernel1_after_chain_flt");

    let bns = fmi.base.idx.bns.as_ref().expect("loaded bns");
    let pac = &fmi.base.idx.pac;
    for chn in chain_ar
        .iter_mut()
        .take(nseq as usize)
    {
        mem_flt_chained_seeds(opt, bns, pac, seq_, chn.n as i32, &mut chn.a);
    }
    debug_rss("kernel1_after_flt_chained");

    mmc.matchArray[tid_usize] = matchArray;
    mmc.min_intv_ar[tid_usize] = min_intv_ar;
    mmc.query_pos_ar[tid_usize] = query_pos_ar;
    mmc.enc_qdb[tid_usize] = enc_qdb;
    mmc.rid[tid_usize] = rid;
    1
}

#[doc = "Original function: mem_kernel2_core:1093"]
pub fn mem_kernel2_core(
    fmi: &FMI_search,
    opt: &mem_opt_t,
    seq_: &[bseq1_t],
    regs: &mut [mem_alnreg_v],
    nseq: i32,
    chain_ar: &mut [mem_chain_v],
    mmc: &mut mem_cache,
    ref_string: &[u8],
    tid: i32,
) -> i32 {
    for reg in regs.iter_mut().take(nseq as usize) {
        clear_vec_with_cap(&mut reg.a, KEEP_REG_CAP);
        reg.n = 0;
        reg.m = 0;
    }

    let bns = fmi.base.idx.bns.as_ref().expect("loaded bns");
    let pac = &fmi.base.idx.pac;
    mem_chain2aln_across_reads_V2(
        opt, bns, pac, seq_, nseq, chain_ar, regs, mmc, ref_string, tid,
    );
    for (l, reg) in regs
        .iter()
        .take(nseq as usize)
        .enumerate()
    {
        if debug_trace_read(seq_[l].name.as_deref()) {
            eprintln!(
                "[trace::kernel2:after_chain2aln] read={} regs={:?}",
                seq_[l].name.as_deref().unwrap_or(""),
                reg.a
                    .iter()
                    .take(reg.n.min(16))
                    .map(|r| (
                        r.score,
                        r.sub,
                        r.csub,
                        r.qb,
                        r.qe,
                        r.rb,
                        r.re,
                        r.secondary,
                        r.secondary_all,
                        r.c
                    ))
                    .collect::<Vec<_>>(),
            );
        }
    }

    for chain in chain_ar
        .iter_mut()
        .take(nseq as usize)
    {
        clear_chain_a_pooled(&mut chain.a, KEEP_CHAIN_CAP);
        chain.n = 0;
        chain.m = 0;
    }

    for reg in regs.iter_mut().take(nseq as usize) {
        let mut m = 0_usize;
        for i in 0..reg.n {
            if reg.a[i].qe > reg.a[i].qb {
                if m != i {
                    reg.a[m] = reg.a[i];
                }
                m += 1;
            }
        }
        reg.n = m;
        reg.a.truncate(m);
    }

    // Hoisted to amortize the per-read encoded_query allocation across nseq iterations.
    // mem_sort_dedup_patch needs &mut [u8] (bwa_gen_cigar2 reverses in place + restores),
    // so we copy seq_nt4 into this buffer per read but keep the capacity.
    let mut encoded_query: Vec<u8> = Vec::new();
    for (l, reg) in regs
        .iter_mut()
        .take(nseq as usize)
        .enumerate()
    {
        let nt4 = nt4_cow(&seq_[l]);
        encoded_query.clear();
        encoded_query.extend_from_slice(&nt4);
        reg.n = mem_sort_dedup_patch(
            opt,
            Some(bns),
            Some(pac),
            Some(&mut encoded_query),
            reg.n as i32,
            &mut reg.a,
        ) as usize;
        reg.m = reg.a.len();
        if debug_trace_read(seq_[l].name.as_deref()) {
            eprintln!(
                "[trace::kernel2:after_dedup] read={} regs={:?}",
                seq_[l].name.as_deref().unwrap_or(""),
                reg.a
                    .iter()
                    .take(reg.n.min(16))
                    .map(|r| (
                        r.score,
                        r.sub,
                        r.csub,
                        r.qb,
                        r.qe,
                        r.rb,
                        r.re,
                        r.secondary,
                        r.secondary_all,
                        r.c
                    ))
                    .collect::<Vec<_>>(),
            );
        }
    }

    for reg in regs.iter_mut().take(nseq as usize) {
        for p in reg.a.iter_mut().take(reg.n) {
            if p.rid >= 0 && bns.anns[p.rid as usize].is_alt != 0 {
                p.is_alt = 1;
            }
        }
    }

    1
}

#[doc = "Original function: worker_aln:1175"]
pub fn worker_aln(data: &mut worker_t, seq_id: i32, batch_size: i32, tid: i32) {
    let opt = data.opt.as_deref().expect("worker opt");
    let fmi = data.fmi.as_ref().expect("worker fmi");
    let seq_start = seq_id as usize;
    let batch = batch_size as usize;
    mem_kernel2_core(
        fmi,
        opt,
        &data.seqs[seq_start..seq_start + batch],
        &mut data.regs[seq_start..seq_start + batch],
        batch_size,
        &mut data.chain_ar[seq_start..seq_start + batch],
        &mut data.mmc,
        &data.ref_string,
        tid,
    );
}

#[doc = "Original function: worker_bwt:1193"]
pub fn worker_bwt(data: &mut worker_t, seq_id: i32, batch_size: i32, tid: i32) {
    let opt = data.opt.as_deref().expect("worker opt");
    let fmi = data.fmi.as_ref().expect("worker fmi");
    let seq_start = seq_id as usize;
    let batch = batch_size as usize;
    for seq in data.seqs[seq_start..seq_start + batch].iter_mut() {
        if seq.seq_nt4.is_empty() {
            if let Some(text) = seq.seq.as_deref() {
                seq.seq_nt4 = seq_to_nt4(text);
            }
        }
    }
    mem_kernel1_core(
        fmi,
        opt,
        &data.seqs[seq_start..seq_start + batch],
        batch_size,
        &mut data.chain_ar[seq_start..seq_start + batch],
        &mut [],
        data.seedBufSize,
        &mut data.mmc,
        tid,
    );
}

#[doc = "Original function: sort_classify:1216"]
pub fn sort_classify(mmc: &mut mem_cache, pcnt: i64, tid: i32) -> i64 {
    let tid_usize = tid as usize;
    let count = pcnt as usize;
    let seqPairArray = &mut mmc.seqPairArrayLeft128[tid_usize];
    let seqPairArrayAux = &mut mmc.seqPairArrayRight128[tid_usize];

    let mut pos8 = 0_usize;
    let mut pos16 = 0_usize;
    for i in 0..count {
        let s = seqPairArray[i];
        let xtra = s.h0;
        let size = if (xtra & KSW_XBYTE) != 0 { 1 } else { 2 };
        if size == 1 {
            seqPairArray[pos8] = seqPairArray[i];
            pos8 += 1;
        } else {
            if seqPairArrayAux.len() <= pos16 {
                seqPairArrayAux.push(s);
            } else {
                seqPairArrayAux[pos16] = s;
            }
            pos16 += 1;
        }
    }
    assert_eq!(pos8 + pos16, count);
    for i in pos8..count {
        seqPairArray[i] = seqPairArrayAux[i - pos8];
    }
    pos8 as i64
}

#[doc = "Original function: worker_sam:1245"]
pub fn worker_sam(data: &mut worker_t, seqid: i32, batch_size: i32, tid: i32) {
    if data.opt.as_ref().expect("worker opt missing").flag & MEM_F_PE != 0 {
        let start = seqid as usize;
        let end = start + batch_size as usize;
        assert_eq!(
            (end - start) & 1,
            0,
            "paired-end worker_sam requires an even batch size"
        );
        let tid_usize = tid as usize;
        let opt = data.opt.as_ref().expect("worker opt missing");
        let fmi = data.fmi.as_ref().expect("worker fmi missing");
        let bns = fmi.base.idx.bns.as_ref().expect("fmi bns missing");
        let pac = &fmi.base.idx.pac;

        let mut pcnt = 0_i32;
        let mut maxRefLen = 0_i32;
        let mut maxQerLen = 0_i32;
        let mut gcnt = 0_i32;
        let mut pos = start >> 1;
        let timing = debug_timings();
        let t0 = std::time::Instant::now();
        // Pre-reserve seqPairArrayAux for the upper bound: 4 entries per pair × pairs in batch.
        // The mid-batch grow loops in mem_sam_pe_batch_pre/mem_matesw_batch_pre push one entry at
        // a time; reserving up front avoids the doubling-realloc churn during warmup.
        // Also pre-reserve seqPairArrayLeft128 for at most 4 entries per pair (one per direction
        // when matesw rescue is admitted). Same reasoning — the per-pair resize() in
        // mem_matesw_batch_pre at line 1787 grows by 1 each call.
        {
            let target_aux = ((end - start) >> 1) * 4;
            let aux = &mut data.mmc.seqPairArrayAux[tid_usize];
            if aux.capacity() < target_aux {
                aux.reserve(target_aux - aux.capacity());
            }
            let left = &mut data.mmc.seqPairArrayLeft128[tid_usize];
            if left.capacity() < target_aux {
                left.reserve(target_aux - left.capacity());
            }
            // Estimate seqBufLeftRef/Qer capacity: target_aux entries × MAX_SEQ_LEN_REF (~256)
            // bytes per ref entry, MAX_SEQ_LEN_QER (~128) per query entry. Pre-reserving here
            // prevents the per-pair resize() in mem_matesw_batch_pre (lines 1777/1780) from
            // doubling-realloc during warmup batches.
            let target_ref = target_aux * 256;
            let target_qer = target_aux * 128;
            let ref_buf = &mut data.mmc.seqBufLeftRef[tid_usize];
            if ref_buf.capacity() < target_ref {
                ref_buf.reserve(target_ref - ref_buf.capacity());
            }
            let qer_buf = &mut data.mmc.seqBufLeftQer[tid_usize];
            if qer_buf.capacity() < target_qer {
                qer_buf.reserve(target_qer - qer_buf.capacity());
            }
        }
        for i in (start..end).step_by(2) {
            let pair_id = ((data.n_processed >> 1) + pos as i64) as u64;
            pos += 1;
            let seq_pair: &mut [bseq1_t; 2] = (&mut data.seqs[i..i + 2])
                .try_into()
                .expect("paired seq slice");
            let reg_pair: &mut [mem_alnreg_v; 2] = (&mut data.regs[i..i + 2])
                .try_into()
                .expect("paired reg slice");
            mem_sam_pe_batch_pre(
                opt,
                bns,
                pac,
                &data.pes,
                pair_id,
                seq_pair,
                reg_pair,
                &mut data.mmc,
                &mut pcnt,
                &mut gcnt,
                &mut maxRefLen,
                &mut maxQerLen,
                tid_usize,
            );
        }
        let t1 = std::time::Instant::now();

        let pcnt8 = sort_classify(&mut data.mmc, i64::from(pcnt), tid) as i32;
        let aln_len = (pcnt + 256) as usize;
        WORKER_SAM_ALN_SCRATCH.with(|cell| {
            let mut aln = cell.borrow_mut();
            aln.clear();
            aln.resize(aln_len, kswr_t::default());
            mem_sam_pe_batch(
                opt,
                &mut data.mmc,
                pcnt,
                pcnt8,
                &mut aln,
                maxRefLen,
                maxQerLen,
                tid_usize,
            );
        });
        let t2 = std::time::Instant::now();

        gcnt = 0;
        pos = start >> 1;
        WORKER_SAM_ALN_SCRATCH.with(|cell| {
            let aln = cell.borrow();
            for i in (start..end).step_by(2) {
                let pair_id = ((data.n_processed >> 1) + pos as i64) as u64;
                pos += 1;
                let seq_pair: &mut [bseq1_t; 2] = (&mut data.seqs[i..i + 2])
                    .try_into()
                    .expect("paired seq slice");
                let reg_pair: &mut [mem_alnreg_v; 2] = (&mut data.regs[i..i + 2])
                    .try_into()
                    .expect("paired reg slice");
                mem_sam_pe_batch_post(
                    opt,
                    bns,
                    pac,
                    &data.pes,
                    pair_id,
                    seq_pair,
                    reg_pair,
                    &aln,
                    &mut data.mmc,
                    &mut gcnt,
                    tid_usize,
                );
                for regs in reg_pair.iter_mut() {
                    clear_vec_with_cap(&mut regs.a, KEEP_REG_CAP);
                    regs.n = 0;
                    regs.m = 0;
                }
                seq_pair[0].name = None;
                seq_pair[0].comment = None;
                seq_pair[0].seq = None;
                seq_pair[0].qual = None;
                seq_pair[1].name = None;
                seq_pair[1].comment = None;
                seq_pair[1].seq = None;
                seq_pair[1].qual = None;
                seq_pair[0].seq_nt4 = Vec::new();
                seq_pair[1].seq_nt4 = Vec::new();
            }
        });
        if timing {
            let t3 = std::time::Instant::now();
            eprintln!(
                "[timing::worker_sam] seqid={} n={} tid={} pcnt={} pre={:.3}s batch_sw={:.3}s post={:.3}s total={:.3}s",
                seqid,
                batch_size,
                tid,
                pcnt,
                (t1 - t0).as_secs_f64(),
                (t2 - t1).as_secs_f64(),
                (t3 - t2).as_secs_f64(),
                (t3 - t0).as_secs_f64(),
            );
        }
    } else {
        let opt = data.opt.as_ref().expect("worker opt missing");
        let fmi = data.fmi.as_ref().expect("worker fmi missing");
        let start = seqid as usize;
        let end = start + batch_size as usize;
        for i in start..end {
            mem_mark_primary_se(
                opt,
                data.regs[i].n as i32,
                &mut data.regs[i].a,
                data.n_processed + i as i64,
            );
            if (opt.flag & MEM_F_PRIMARY5) != 0 {
                mem_reorder_primary5(opt.T, &mut data.regs[i]);
            }
            mem_reg2sam(
                opt,
                fmi.base.idx.bns.as_ref().expect("fmi bns missing"),
                &fmi.base.idx.pac,
                &mut data.seqs[i],
                &mut data.regs[i],
                0,
                None,
            );
            clear_vec_with_cap(&mut data.regs[i].a, KEEP_REG_CAP);
            data.regs[i].n = 0;
            data.regs[i].m = 0;
            data.seqs[i].name = None;
            data.seqs[i].comment = None;
            data.seqs[i].seq = None;
            data.seqs[i].qual = None;
            data.seqs[i].seq_nt4 = Vec::new();
        }
    }
}

#[doc = "Original function: mem_process_seqs:1338"]
pub fn mem_process_seqs(
    opt: &mut mem_opt_t,
    n_processed: i64,
    n: i32,
    seqs: &mut Vec<bseq1_t>,
    pes0: Option<&[mem_pestat_t; 4]>,
    w: &mut worker_t,
) {
    w.opt = Some(Box::new(opt.clone()));
    w.n_processed = n_processed;
    w.nthreads = i16::try_from(opt.n_threads.max(1)).expect("nthreads");
    ensure_mem_cache_thread_slots(
        &mut w.mmc,
        w.nthreads.max(1) as usize,
    );
    w.seqs = std::mem::take(seqs);
    let n_usize = n.max(0) as usize;
    if w.regs.len() < n_usize {
        w.regs.resize(n_usize, mem_alnreg_v::default());
    }
    if w.chain_ar.len() < n_usize {
        w.chain_ar.resize(n_usize, mem_chain_v::default());
    }
    w.seedBuf.clear();
    w.seedBufSize = 0;
    w.nreads = n;
    if w.ref_string.is_empty() {
        let fmi = w.fmi.as_ref().expect("worker fmi missing");
        let bns = fmi.base.idx.bns.as_ref().expect("fmi bns missing");
        w.ref_string = pac_to_reference_layout(bns.l_pac, &fmi.base.idx.pac);
    }

    let n_ = n;
    if opt.flag & MEM_F_PE != 0 {
        let timing = debug_timings();
        let t0 = std::time::Instant::now();
        let t1 = if w.nthreads <= 1 {
            run_worker_bwt_aln_chunks_serial(w, n_);
            trim_allocator_rss();
            debug_rss("after_worker_bwt_aln_serial");
            std::time::Instant::now()
        } else {
            run_worker_bwt_aln_chunks_parallel(w, n_);
            trim_allocator_rss();
            debug_rss("after_worker_bwt_aln");
            std::time::Instant::now()
        };
        let t2 = std::time::Instant::now();
        if let Some(pes0) = pes0 {
            w.pes = *pes0;
        } else {
            let fmi = w.fmi.as_ref().expect("worker fmi missing");
            let bns = fmi.base.idx.bns.as_ref().expect("fmi bns missing");
            mem_pestat(opt, bns.l_pac, n, &w.regs, &mut w.pes);
            if BWA_VERBOSE.load(std::sync::atomic::Ordering::Relaxed) >= 4 {
                eprintln!(
                    "[dbg::mem_process_seqs] inferred_pes n={} t={} => \
FF(low={},high={},failed={},avg={:.2},std={:.2}) \
FR(low={},high={},failed={},avg={:.2},std={:.2}) \
RF(low={},high={},failed={},avg={:.2},std={:.2}) \
RR(low={},high={},failed={},avg={:.2},std={:.2})",
                    n,
                    opt.n_threads,
                    w.pes[0].low,
                    w.pes[0].high,
                    w.pes[0].failed,
                    w.pes[0].avg,
                    w.pes[0].std,
                    w.pes[1].low,
                    w.pes[1].high,
                    w.pes[1].failed,
                    w.pes[1].avg,
                    w.pes[1].std,
                    w.pes[2].low,
                    w.pes[2].high,
                    w.pes[2].failed,
                    w.pes[2].avg,
                    w.pes[2].std,
                    w.pes[3].low,
                    w.pes[3].high,
                    w.pes[3].failed,
                    w.pes[3].avg,
                    w.pes[3].std,
                );
            }
        }
        let t3 = std::time::Instant::now();
        run_worker_chunks_parallel(w, n_, worker_sam);
        trim_allocator_rss();
        debug_rss("after_worker_sam");
        if timing {
            let t4 = std::time::Instant::now();
            eprintln!(
                "[timing::mem_process_seqs] n={} t={} bwt={:.3}s aln={:.3}s pestat={:.3}s sam={:.3}s total={:.3}s",
                n_,
                opt.n_threads,
                (t1 - t0).as_secs_f64(),
                (t2 - t1).as_secs_f64(),
                (t3 - t2).as_secs_f64(),
                (t4 - t3).as_secs_f64(),
                (t4 - t0).as_secs_f64(),
            );
        }
    } else {
        let timing = debug_timings();
        let t0 = std::time::Instant::now();
        let (t1, t2) = if w.nthreads <= 1 {
            run_worker_bwt_aln_chunks_serial(w, n_);
            trim_allocator_rss();
            debug_rss("after_worker_bwt_aln_serial");
            let t = std::time::Instant::now();
            (t, t)
        } else {
            run_worker_bwt_aln_chunks_parallel(w, n_);
            trim_allocator_rss();
            debug_rss("after_worker_bwt_aln");
            let t = std::time::Instant::now();
            (t, t)
        };
        run_worker_chunks_parallel(w, n_, worker_sam);
        trim_allocator_rss();
        debug_rss("after_worker_sam");
        if timing {
            let t3 = std::time::Instant::now();
            eprintln!(
                "[timing::mem_process_seqs] n={} t={} bwt+aln={:.3}s sam={:.3}s total={:.3}s",
                n_,
                opt.n_threads,
                (t1 - t0).as_secs_f64(),
                (t3 - t2).as_secs_f64(),
                (t3 - t0).as_secs_f64(),
            );
        }
    }
    *seqs = std::mem::take(&mut w.seqs);
}

#[doc = "Original function: mem_mark_primary_se_core:1392"]
pub fn mem_mark_primary_se_core(
    opt: &mem_opt_t,
    n: i32,
    a: &mut [mem_alnreg_t],
    z: &mut Vec<usize>,
) {
    let mut tmp = opt.a + opt.b;
    tmp = tmp.max(opt.o_del + opt.e_del);
    tmp = tmp.max(opt.o_ins + opt.e_ins);
    z.clear();
    z.push(0);
    let n_us = n as usize;
    for i in 1..n_us {
        let mut k = 0_usize;
        while k < z.len() {
            let j = z[k];
            let b_max = a[j].qb.max(a[i].qb);
            let e_min = a[j].qe.min(a[i].qe);
            if e_min > b_max {
                let min_l = (a[i].qe - a[i].qb).min(a[j].qe - a[j].qb);
                if ((e_min - b_max) as f32) >= (min_l as f32) * opt.mask_level {
                    if a[j].sub == 0 {
                        a[j].sub = a[i].score;
                    }
                    if a[j].score - a[i].score <= tmp && (a[j].is_alt != 0 || a[i].is_alt == 0) {
                        a[j].sub_n += 1;
                    }
                    break;
                }
            }
            k += 1;
        }
        if k == z.len() {
            z.push(i);
        } else {
            a[i].secondary = z[k] as i32;
        }
    }
}

#[doc = "Original function: kv_push:1399"]
pub fn kv_push__L1399(_arg0: crate::support::Opaque, _arg1: crate::support::Opaque) {
    crate::support::stub::<()>("kv_push")
}

#[doc = "Original function: mem_mark_primary_se:1420"]
pub fn mem_mark_primary_se(opt: &mem_opt_t, n: i32, a: &mut [mem_alnreg_t], id: i64) -> i32 {
    if n == 0 {
        return 0;
    }
    let mut n_pri = 0_i32;
    let n_us = n as usize;
    for (i, p) in a.iter_mut().take(n_us).enumerate() {
        p.sub = 0;
        p.alt_sc = 0;
        p.secondary = -1;
        p.secondary_all = -1;
        p.hash = hash_64(u64::from_ne_bytes((id + i as i64).to_ne_bytes()));
        if p.is_alt == 0 {
            n_pri += 1;
        }
    }
    a[..n_us].sort_by(|lhs, rhs| {
        rhs.score
            .cmp(&lhs.score)
            .then_with(|| lhs.is_alt.cmp(&rhs.is_alt))
            .then_with(|| lhs.hash.cmp(&rhs.hash))
    });
    MEM_MARK_PRIMARY_SE_SCRATCH.with(|cell| {
        let mut scratch = cell.borrow_mut();
        let (z, rank) = &mut *scratch;
        mem_mark_primary_se_core(opt, n, a, z);
        for i in 0..n_us {
            // C++ bwamem.cpp:1437-1438: assigns alt_sc unconditionally when the secondary points
            // to an ALT hit (regardless of whether sec.score is zero). Match that structure.
            let secondary = a[i].secondary;
            let assign_alt_sc =
                a[i].is_alt == 0 && secondary >= 0 && a[secondary as usize].is_alt != 0;
            let alt_score = if assign_alt_sc {
                a[secondary as usize].score
            } else {
                0
            };
            let p = &mut a[i];
            p.secondary_all = i as i32;
            if assign_alt_sc {
                p.alt_sc = alt_score;
            }
        }
        if n_pri >= 0 && n_pri < n {
            if n_pri > 0 {
                a[..n_us].sort_by(|lhs, rhs| {
                    lhs.is_alt
                        .cmp(&rhs.is_alt)
                        .then_with(|| rhs.score.cmp(&lhs.score))
                        .then_with(|| lhs.hash.cmp(&rhs.hash))
                });
            }
            rank.clear();
            rank.resize(n_us, 0_usize);
            for (i, p) in a.iter().take(n_us).enumerate() {
                rank[p.secondary_all as usize] = i;
            }
            for i in 0..n_us {
                if a[i].secondary >= 0 {
                    a[i].secondary_all = rank[a[i].secondary as usize] as i32;
                    if a[i].is_alt != 0 {
                        a[i].secondary = i32::MAX;
                    }
                } else {
                    a[i].secondary_all = -1;
                }
            }
            if n_pri > 0 {
                let n_pri_us = n_pri as usize;
                for p in a.iter_mut().take(n_pri_us) {
                    p.sub = 0;
                    p.secondary = -1;
                }
                z.clear();
                mem_mark_primary_se_core(opt, n_pri, &mut a[..n_pri_us], z);
            }
        } else {
            for p in a.iter_mut().take(n_us) {
                p.secondary_all = p.secondary;
            }
        }
    });
    n_pri
}

#[doc = "Original function: mem_approx_mapq_se:1470"]
#[inline]
pub fn mem_approx_mapq_se(opt: &mem_opt_t, a: &mem_alnreg_t) -> i32 {
    let mut sub = if a.sub != 0 {
        a.sub
    } else {
        opt.min_seed_len * opt.a
    };
    sub = sub.max(a.csub);
    if sub >= a.score {
        return 0;
    }
    let l = (a.qe - a.qb).max((a.re - a.rb) as i32);
    let identity = 1.0 - f64::from(l * opt.a - a.score) / f64::from(opt.a + opt.b) / f64::from(l);
    let mut mapq = if a.score == 0 {
        0
    } else if opt.mapQ_coef_len > 0.0 {
        let mut tmp = if (l as f32) < opt.mapQ_coef_len {
            1.0
        } else {
            f64::from(opt.mapQ_coef_fac) / f64::from(l).ln()
        };
        tmp *= identity * identity;
        (6.02 * f64::from(a.score - sub) / f64::from(opt.a) * tmp * tmp + 0.499) as i32
    } else {
        let mut mapq = (MEM_MAPQ_COEF
            * (1.0 - f64::from(sub) / f64::from(a.score))
            * f64::from(a.seedcov).ln()
            + 0.499) as i32;
        if identity < 0.95 {
            mapq = (f64::from(mapq) * identity * identity + 0.499) as i32;
        }
        mapq
    };
    if a.sub_n > 0 {
        mapq -= (4.343 * f64::from(a.sub_n + 1).ln() + 0.499) as i32;
    }
    mapq = mapq.clamp(0, MEM_MAPQ_MAX);
    (f64::from(mapq) * (1.0 - f64::from(a.frac_rep)) + 0.499) as i32
}

#[doc = "Original function: mem_reorder_primary5:1496"]
#[inline]
pub fn mem_reorder_primary5(T: i32, a: &mut mem_alnreg_v) {
    let mut n_pri = 0_i32;
    let mut left_st = i32::MAX;
    let mut left_k = None;
    for (k, p) in a.a.iter().take(a.n).enumerate() {
        if p.secondary < 0 && p.is_alt == 0 && p.score >= T {
            n_pri += 1;
            if p.qb < left_st {
                left_st = p.qb;
                left_k = Some(k);
            }
        }
    }
    if n_pri <= 1 {
        return;
    }
    let left_k = left_k.expect("left_k");
    if left_k == 0 {
        return;
    }
    a.a.swap(0, left_k);
    let left_k_i32 = left_k as i32;
    for p in a.a.iter_mut().take(a.n).skip(1) {
        if p.secondary == 0 {
            p.secondary = left_k_i32;
        } else if p.secondary == left_k_i32 {
            p.secondary = 0;
        }
        if p.secondary_all == 0 {
            p.secondary_all = left_k_i32;
        } else if p.secondary_all == left_k_i32 {
            p.secondary_all = 0;
        }
    }
}

#[doc = "Original function: mem_reg2sam:1521"]
pub fn mem_reg2sam(
    opt: &mem_opt_t,
    bns: &bntseq_t,
    pac: &[u8],
    s: &mut bseq1_t,
    a: &mut mem_alnreg_v,
    extra_flag: i32,
    m: Option<&mem_aln_t>,
) {
    let mut xa = if (opt.flag & MEM_F_ALL) == 0 {
        mem_gen_alt(opt, bns, pac, a, s.l_seq, s.seq.as_deref().unwrap_or(""))
    } else {
        vec![None; a.n]
    };
    // Take aa from thread-local; restored cleared at function end.
    let mut aa = MEM_REG2SAM_AA_SCRATCH.with(|c| std::mem::take(&mut *c.borrow_mut()));
    // Pre-size for a typical SAM record (~500B) to skip the early grow-allocs that show up in
    // mem_aln2sam's per-record kputw/kputs calls. The kstring buffer is moved into s.sam at the
    // end, so we pay one malloc per record either way — this just makes it the right size.
    let mut str_ = kstring_t::default();
    ks_resize(&mut str_, 512);
    let mut l = 0_i32;

    for k in 0..a.n {
        let p = &a.a[k];
        if p.score < opt.T {
            continue;
        }
        if p.secondary >= 0 && (p.is_alt != 0 || (opt.flag & MEM_F_ALL) == 0) {
            continue;
        }
        if p.secondary >= 0
            && (p.secondary as usize) < a.a.len()
            && (p.score as f32)
                < (a.a[p.secondary as usize].score as f32) * opt.drop_ratio
        {
            continue;
        }

        let mut q = mem_reg2aln(
            opt,
            bns,
            pac,
            s.l_seq,
            s.seq.as_deref().unwrap_or(""),
            Some(p),
        );
        assert!(q.rid >= 0, "mem_reg2sam expected mapped alignment");
        // Take xa[k] by value (xa not used after this loop) — avoids the String clone.
        q.XA = std::mem::take(&mut xa[k]);
        q.flag |= extra_flag;
        if p.secondary >= 0 {
            q.sub = -1;
        }
        if l != 0 && p.secondary < 0 {
            q.flag |= if (opt.flag & MEM_F_NO_MULTI) != 0 {
                0x10000
            } else {
                0x800
            };
        }
        if (opt.flag & MEM_F_KEEP_SUPP_MAPQ) == 0
            && l != 0
            && p.is_alt == 0
            && !aa.is_empty()
            && q.mapq > aa[0].mapq
        {
            q.mapq = aa[0].mapq;
        }
        aa.push(q);
        l += 1;
    }

    if aa.is_empty() {
        let mut t = mem_reg2aln(opt, bns, pac, s.l_seq, s.seq.as_deref().unwrap_or(""), None);
        t.flag |= extra_flag;
        mem_aln2sam(opt, bns, &mut str_, s, 1, &[t], 0, m);
    } else {
        let n = aa.len() as i32;
        for k in 0..aa.len() {
            mem_aln2sam(opt, bns, &mut str_, s, n, &aa, k as i32, m);
        }
    }

    // Move the kstring buffer into the String to avoid the as_str().to_owned() clone.
    // SAFETY: every byte written to `str_` comes from kput* helpers (ASCII) or from FASTQ
    // sequence/quality, which the SAM/FASTQ formats define as printable ASCII. Skipping the
    // UTF-8 validation pass shaves ~hundreds of MB of byte scans on a 700K-read run.
    let mut bytes = str_.s;
    bytes.truncate(str_.l);
    s.sam = Some(unsafe { String::from_utf8_unchecked(bytes) });
    // Drain mem_aln_ts and return their cigar/md buffers to the global pools (preserving
    // capacity for next call's mem_reg2aln → bwa_gen_cigar2). XA String is just dropped.
    for aln in aa.drain(..) {
        crate::generated::bwa_cpp::return_cigar_buf(aln.cigar);
        crate::generated::bwa_cpp::return_md_buf(aln.md);
    }
    MEM_REG2SAM_AA_SCRATCH.with(|c| *c.borrow_mut() = aa);
}

#[doc = "Original function: add_cigar:1579"]
#[inline]
pub fn add_cigar(
    opt: &mem_opt_t,
    p: &mem_aln_t,
    str_: &mut crate::generated::kstring_h::kstring_t,
    which: i32,
) {
    add_cigar_fields(opt, p.n_cigar, &p.cigar, p.is_alt, str_, which);
}

#[inline]
pub fn add_cigar_fields(
    opt: &mem_opt_t,
    n_cigar: i32,
    cigar: &[u32],
    is_alt: u32,
    str_: &mut crate::generated::kstring_h::kstring_t,
    which: i32,
) {
    if n_cigar > 0 {
        const MIDSH_LUT: [u8; 5] = *b"MIDSH";
        for &c_op in cigar.iter().take(n_cigar as usize) {
            let mut c = (c_op & 0xf) as i32;
            if (opt.flag & 0x200) == 0 && is_alt == 0 && (c == 3 || c == 4) {
                c = if which != 0 { 4 } else { 3 };
            }
            kputw((c_op >> 4) as i32, str_);
            // c is 0..=4 (CIGAR ops M/I/D/S/H). Const-size LUT for bounds-check elision.
            kputc(
                i32::from(unsafe { *MIDSH_LUT.get_unchecked((c as usize).min(4)) }),
                str_,
            );
        }
    } else {
        kputc(i32::from(b'*'), str_);
    }
}

#[doc = "Original function: mem_aln2sam:1592"]
pub fn mem_aln2sam(
    opt: &mem_opt_t,
    bns: &bntseq_t,
    str_: &mut kstring_t,
    s: &bseq1_t,
    n: i32,
    list: &[mem_aln_t],
    which: i32,
    m_: Option<&mem_aln_t>,
) {
    fn cigar_end_clip(cigar: u32) -> bool {
        let op = cigar & 0xf;
        op == 3 || op == 4
    }

    let which_usize = which as usize;
    // Avoid cloning p (cigar/md Vecs + Option<String> XA per call). Track mutated Copy fields in
    // shadow vars; read non-mutated heap fields (cigar, md, XA) directly from list[which_usize].
    let p_ref = &list[which_usize];
    let mut p_flag = p_ref.flag;
    let mut p_rid = p_ref.rid;
    let mut p_pos = p_ref.pos;
    let mut p_is_rev = p_ref.is_rev;
    let mut p_n_cigar = p_ref.n_cigar;
    // Avoid cloning the mate's cigar/md/XA. Track mutated mate fields in shadow vars; read
    // non-mutated fields (cigar, md) directly from m_ when needed.
    let mut mate_rid = m_.map_or(-1, |m| m.rid);
    let mut mate_pos = m_.map_or(0, |m| m.pos);
    let mut mate_is_rev = m_.map_or(0, |m| m.is_rev);
    let mut mate_n_cigar = m_.map_or(0, |m| m.n_cigar);

    p_flag |= if m_.is_some() { 0x1 } else { 0 };
    p_flag |= if p_rid < 0 { 0x4 } else { 0 };
    p_flag |= if m_.is_some() && mate_rid < 0 { 0x8 } else { 0 };
    if p_rid < 0 && m_.is_some() && mate_rid >= 0 {
        p_rid = mate_rid;
        p_pos = mate_pos;
        p_is_rev = mate_is_rev;
        p_n_cigar = 0;
    }
    if m_.is_some() && mate_rid < 0 && p_rid >= 0 {
        mate_rid = p_rid;
        mate_pos = p_pos;
        mate_is_rev = p_is_rev;
        mate_n_cigar = 0;
    }
    p_flag |= if p_is_rev != 0 { 0x10 } else { 0 };
    p_flag |= if m_.is_some() && mate_is_rev != 0 {
        0x20
    } else {
        0
    };

    let qname = s.name.as_deref().expect("read name missing");
    let l_seq = s.l_seq as usize;
    let seq_text = s.seq.as_deref().unwrap_or("");
    let qual_text = s.qual.as_deref();

    ks_resize(
        str_,
        str_.l + l_seq + qname.len() + qual_text.map_or(0, str::len) + 20,
    );
    kputsn(qname.as_bytes(), qname.len() as i32, str_);
    kputc(i32::from(b'\t'), str_);
    kputw(
        (p_flag & 0xffff) | if (p_flag & 0x10000) != 0 { 0x100 } else { 0 },
        str_,
    );
    kputc(i32::from(b'\t'), str_);

    if p_rid >= 0 {
        let ann = &bns.anns[p_rid as usize];
        kputs(&ann.name, str_);
        kputc(i32::from(b'\t'), str_);
        kputl(p_pos + 1, str_);
        kputc(i32::from(b'\t'), str_);
        kputw(p_ref.mapq as i32, str_);
        kputc(i32::from(b'\t'), str_);
        add_cigar_fields(opt, p_n_cigar, &p_ref.cigar, p_ref.is_alt, str_, which);
    } else {
        kputsn(b"*\t0\t0\t*", 7, str_);
    }
    kputc(i32::from(b'\t'), str_);

    if m_.is_some() {
        if mate_rid >= 0 {
            if p_rid == mate_rid {
                kputc(i32::from(b'='), str_);
            } else {
                kputs(&bns.anns[mate_rid as usize].name, str_);
            }
            kputc(i32::from(b'\t'), str_);
            kputl(mate_pos + 1, str_);
            kputc(i32::from(b'\t'), str_);
            if p_rid == mate_rid {
                let p0 = p_pos
                    + if p_is_rev != 0 && p_n_cigar > 0 {
                        i64::from(get_rlen(p_n_cigar, &p_ref.cigar) - 1)
                    } else {
                        0
                    };
                // Note: when mate_n_cigar was zeroed by the mutation path, mate.cigar isn't read
                // (the n_cigar==0 branch below catches it). Otherwise the original mate.cigar is
                // intact from m_.
                let p1 = mate_pos
                    + if mate_is_rev != 0 && mate_n_cigar > 0 {
                        i64::from(get_rlen(mate_n_cigar, &m_.expect("mate present").cigar) - 1)
                    } else {
                        0
                    };
                if mate_n_cigar == 0 || p_n_cigar == 0 {
                    kputc(i32::from(b'0'), str_);
                } else {
                    let tlen = -(p0 - p1
                        + if p0 > p1 {
                            1
                        } else if p0 < p1 {
                            -1
                        } else {
                            0
                        });
                    kputl(tlen, str_);
                }
            } else {
                kputc(i32::from(b'0'), str_);
            }
        } else {
            kputsn(b"*\t0\t0", 5, str_);
        }
    } else {
        kputsn(b"*\t0\t0", 5, str_);
    }
    kputc(i32::from(b'\t'), str_);

    if (p_flag & 0x100) != 0 {
        kputsn(b"*\t*", 3, str_);
    } else {
        let mut qb = 0_usize;
        let mut qe = l_seq;
        if p_n_cigar > 0 && which != 0 && (opt.flag & MEM_F_SOFTCLIP) == 0 && p_ref.is_alt == 0 {
            if cigar_end_clip(p_ref.cigar[0]) {
                qb += (p_ref.cigar[0] >> 4) as usize;
            }
            if cigar_end_clip(*p_ref.cigar.last().expect("non-empty cigar")) {
                qe -= (p_ref.cigar.last().expect("cigar") >> 4) as usize;
            }
        }
        // Reuse s.seq_nt4 if populated (worker_bwt fills it); otherwise encode on the fly.
        let seq_nt4_owned: Option<Vec<u8>>;
        let seq_nt4: &[u8] = if !s.seq_nt4.is_empty() {
            seq_nt4_owned = None;
            &s.seq_nt4
        } else {
            seq_nt4_owned = Some(seq_to_nt4(seq_text));
            seq_nt4_owned.as_deref().expect("seq_nt4 encoded")
        };
        let _ = &seq_nt4_owned; // suppress unused warning when populated path is taken
        ks_resize(str_, str_.l + (qe - qb) + 1);
        // Fixed-size LUTs so the compiler elides bounds checks on the inner per-byte index.
        // nt4-encoded input is 0..=4 by construction (worker_bwt's seq_to_nt4); SAFETY contract
        // is the same as the upstream C++ which indexes `"ACGTN"` directly.
        const FWD_LUT: [u8; 5] = *b"ACGTN";
        const REV_LUT: [u8; 5] = *b"TGCAN";
        if p_is_rev == 0 {
            let dst = str_.l;
            let src = &seq_nt4[qb..qe];
            let n = src.len();
            let out = &mut str_.s[dst..dst + n];
            for (o, &base) in out.iter_mut().zip(src.iter()) {
                *o = unsafe { *FWD_LUT.get_unchecked(usize::from(base.min(4))) };
            }
            str_.l += n;
        } else {
            let mut rq_b = qb;
            let mut rq_e = qe;
            if p_n_cigar > 0 && which != 0 && (opt.flag & MEM_F_SOFTCLIP) == 0 && p_ref.is_alt == 0
            {
                rq_e = l_seq;
                rq_b = 0;
                if cigar_end_clip(p_ref.cigar[0]) {
                    rq_e -= (p_ref.cigar[0] >> 4) as usize;
                }
                if cigar_end_clip(*p_ref.cigar.last().expect("non-empty cigar")) {
                    rq_b += (p_ref.cigar.last().expect("cigar") >> 4) as usize;
                }
            }
            let dst = str_.l;
            let src = &seq_nt4[rq_b..rq_e];
            let n = src.len();
            let out = &mut str_.s[dst..dst + n];
            for (o, &base) in out.iter_mut().zip(src.iter().rev()) {
                *o = unsafe { *REV_LUT.get_unchecked(usize::from(base.min(4))) };
            }
            str_.l += n;
        }
        str_.s[str_.l] = 0;
        kputc(i32::from(b'\t'), str_);
        if let Some(qual) = qual_text {
            ks_resize(str_, str_.l + (qe - qb) + 1);
            if p_is_rev == 0 {
                let dst = str_.l;
                let src = &qual.as_bytes()[qb..qe];
                let n = src.len();
                str_.s[dst..dst + n].copy_from_slice(src);
                str_.l += n;
            } else {
                let mut rq_b = qb;
                let mut rq_e = qe;
                if p_n_cigar > 0
                    && which != 0
                    && (opt.flag & MEM_F_SOFTCLIP) == 0
                    && p_ref.is_alt == 0
                {
                    rq_e = l_seq;
                    rq_b = 0;
                    if cigar_end_clip(p_ref.cigar[0]) {
                        rq_e -= (p_ref.cigar[0] >> 4) as usize;
                    }
                    if cigar_end_clip(*p_ref.cigar.last().expect("non-empty cigar")) {
                        rq_b += (p_ref.cigar.last().expect("cigar") >> 4) as usize;
                    }
                }
                let dst = str_.l;
                let src = &qual.as_bytes()[rq_b..rq_e];
                let n = src.len();
                let out = &mut str_.s[dst..dst + n];
                for (o, &ch) in out.iter_mut().zip(src.iter().rev()) {
                    *o = ch;
                }
                str_.l += n;
            }
            str_.s[str_.l] = 0;
        } else {
            kputc(i32::from(b'*'), str_);
        }
    }

    if p_n_cigar > 0 {
        kputsn(b"\tNM:i:", 6, str_);
        kputw(p_ref.NM as i32, str_);
        kputsn(b"\tMD:Z:", 6, str_);
        // p_ref.md is bwa_gen_cigar2's output: ASCII digits, uppercase letters, '^' separators.
        // Always valid UTF-8; skip the per-record str_::from_utf8 validation.
        kputsn(&p_ref.md, p_ref.md.len() as i32, str_);
    }
    if let Some(mate) = m_ {
        // Use shadow mate_n_cigar (may be 0 from the mutation path even if original > 0).
        if mate_n_cigar > 0 {
            kputsn(b"\tMC:Z:", 6, str_);
            add_cigar(opt, mate, str_, which);
        }
    }
    if p_ref.score >= 0 {
        kputsn(b"\tAS:i:", 6, str_);
        kputw(p_ref.score, str_);
    }
    if p_ref.sub >= 0 {
        kputsn(b"\tXS:i:", 6, str_);
        kputw(p_ref.sub, str_);
    }
    if crate::generated::bwa_cpp::BWA_RG_ID_NONEMPTY
        .load(std::sync::atomic::Ordering::Relaxed)
    {
        let rg_id = crate::generated::bwa_cpp::BWA_RG_ID
            .read()
            .expect("bwa_rg_id lock");
        if !rg_id.is_empty() {
            kputsn(b"\tRG:Z:", 6, str_);
            kputs(rg_id.as_str(), str_);
        }
    }
    if (p_flag & 0x100) == 0 {
        let n_usize = n as usize;
        let has_other_primary = list
            .iter()
            .enumerate()
            .take(n_usize)
            .any(|(i, r)| i != which_usize && (r.flag & 0x100) == 0);
        if has_other_primary {
            kputsn(b"\tSA:Z:", 6, str_);
            for (i, r) in list.iter().enumerate().take(n_usize) {
                if i == which_usize || (r.flag & 0x100) != 0 {
                    continue;
                }
                kputs(&bns.anns[r.rid as usize].name, str_);
                kputc(i32::from(b','), str_);
                kputl(r.pos + 1, str_);
                kputc(i32::from(b','), str_);
                kputc(i32::from(if r.is_rev != 0 { b'-' } else { b'+' }), str_);
                kputc(i32::from(b','), str_);
                const SA_MIDSH_LUT: [u8; 5] = *b"MIDSH";
                for &cigar in r.cigar.iter().take(r.n_cigar as usize) {
                    kputw((cigar >> 4) as i32, str_);
                    let op_idx = (cigar & 0xf) as usize;
                    kputc(
                        i32::from(unsafe { *SA_MIDSH_LUT.get_unchecked(op_idx.min(4)) }),
                        str_,
                    );
                }
                kputc(i32::from(b','), str_);
                kputw(r.mapq as i32, str_);
                kputc(i32::from(b','), str_);
                kputw(r.NM as i32, str_);
                kputc(i32::from(b';'), str_);
            }
        }
        if p_ref.alt_sc > 0 {
            let _ = ksprintf(
                str_,
                format_args!("\tpa:f:{:.3}", (p_ref.score as f64) / (p_ref.alt_sc as f64)),
            );
        }
    }
    if let Some(xa) = p_ref.XA.as_deref() {
        kputsn(b"\tXA:Z:", 6, str_);
        kputs(xa, str_);
    }
    if let Some(comment) = s.comment.as_deref() {
        kputc(i32::from(b'\t'), str_);
        kputs(comment, str_);
    }
    if (opt.flag & MEM_F_REF_HDR) != 0 && p_rid >= 0 {
        let anno = &bns.anns[p_rid as usize].anno;
        if !anno.is_empty() {
            kputsn(b"\tXR:Z:", 6, str_);
            let tmp = str_.l;
            kputs(anno, str_);
            for ch in &mut str_.s[tmp..str_.l] {
                if *ch == b'\t' {
                    *ch = b' ';
                }
            }
        }
    }
    kputc(i32::from(b'\n'), str_);
}

#[doc = "Original function: mem_reg2aln:1732"]
pub fn mem_reg2aln(
    opt: &mem_opt_t,
    bns: &bntseq_t,
    pac: &[u8],
    l_query: i32,
    query_: &str,
    ar: Option<&mem_alnreg_t>,
) -> mem_aln_t {
    let mut a = mem_aln_t::default();
    let Some(ar) = ar else {
        a.rid = -1;
        a.pos = -1;
        a.flag |= 0x4;
        return a;
    };
    if ar.rb < 0 || ar.re < 0 {
        a.rid = -1;
        a.pos = -1;
        a.flag |= 0x4;
        return a;
    }

    let qb = ar.qb;
    let qe = ar.qe;
    let rb = ar.rb;
    let re = ar.re;
    // The full nt4-encoded query is only needed for the qb..qe slice below.
    // Skip the full seq_to_nt4 allocation by encoding the subslice directly.
    debug_assert_eq!(query_.len(), l_query as usize);

    a.mapq = if ar.secondary < 0 {
        mem_approx_mapq_se(opt, ar) as u32
    } else {
        0
    };
    if ar.secondary >= 0 {
        a.flag |= 0x100;
    }

    let ref_span = (re - rb) as i32;
    let tmp = infer_bw(qe - qb, ref_span, ar.truesc, opt.a, opt.o_del, opt.e_del);
    let mut w2 = infer_bw(qe - qb, ref_span, ar.truesc, opt.a, opt.o_ins, opt.e_ins);
    if w2 < tmp {
        w2 = tmp;
    }
    if w2 > opt.w {
        w2 = w2.min(ar.w);
    }

    let mut score = 0_i32;
    let mut nm = 0_i32;
    let mut last_sc = -(1 << 30);
    let cigar_result;
    let mut tries = 0_i32;
    // Hoist query_slice — bwa_gen_cigar2 does an in-place reverse + restore, so the buffer
    // is unchanged across retries. Saves the per-iter to_vec().
    let qb_us = qb as usize;
    let qe_us = qe as usize;
    let mut query_slice: Vec<u8> = query_.as_bytes()[qb_us..qe_us]
        .iter()
        .map(|&c| NT4_LUT[c as usize])
        .collect();
    loop {
        w2 = w2.min(opt.w << 2);
        let mut n_cigar = 0_i32;
        let next_cigar_result = bwa_gen_cigar2(
            &opt.mat,
            opt.o_del,
            opt.e_del,
            opt.o_ins,
            opt.e_ins,
            w2,
            bns.l_pac,
            pac,
            qe - qb,
            &mut query_slice,
            rb,
            re,
            &mut score,
            Some(&mut n_cigar),
            Some(&mut nm),
        );
        if score == last_sc || w2 == (opt.w << 2) || tries >= 2 || score >= ar.truesc - opt.a {
            a.n_cigar = n_cigar;
            cigar_result = next_cigar_result;
            break;
        }
        last_sc = score;
        w2 <<= 1;
        tries += 1;
    }

    let cigar_result = cigar_result.expect("mem_reg2aln requires CIGAR");
    a.cigar = cigar_result.cigar;
    a.md = cigar_result.md;
    a.n_cigar = a.cigar.len() as i32;
    a.NM = nm as u32;

    let mut is_rev = 0_i32;
    let mut pos = bns_depos(bns, if rb < bns.l_pac { rb } else { re - 1 }, &mut is_rev);
    a.is_rev = is_rev as u32;

    if a.n_cigar > 0 {
        if (a.cigar[0] & 0xf) == 2 {
            pos += i64::from(a.cigar[0] >> 4);
            a.cigar.remove(0);
            a.n_cigar -= 1;
        } else if (a.cigar[(a.n_cigar - 1) as usize] & 0xf) == 2 {
            a.cigar.pop();
            a.n_cigar -= 1;
        }
    }

    if qb != 0 || qe != l_query {
        let clip5 = if is_rev != 0 { l_query - qe } else { qb };
        let clip3 = if is_rev != 0 { qb } else { l_query - qe };
        if clip5 != 0 {
            a.cigar.insert(0, ((clip5 as u32) << 4) | 3);
            a.n_cigar += 1;
        }
        if clip3 != 0 {
            a.cigar.push(((clip3 as u32) << 4) | 3);
            a.n_cigar += 1;
        }
    }

    a.rid = bns_pos2rid(bns, pos);
    assert_eq!(a.rid, ar.rid);
    a.pos = pos - bns.anns[a.rid as usize].offset;
    a.score = ar.score;
    a.sub = ar.sub.max(ar.csub);
    a.is_alt = ar.is_alt;
    a.alt_sc = ar.alt_sc;
    a
}

#[doc = "Original function: infer_bw:1811"]
#[inline]
pub fn infer_bw(l1: i32, l2: i32, score: i32, a: i32, q: i32, r: i32) -> i32 {
    if l1 == l2 && l1 * a - score < (q + r - a) << 1 {
        return 0;
    }
    let mut w = (((l1.min(l2) * a - score - q) as f64) / r as f64 + 2.0) as i32;
    if w < (l1 - l2).abs() {
        w = (l1 - l2).abs();
    }
    w
}

#[doc = "Original function: get_rlen:1820"]
#[inline]
pub fn get_rlen(n_cigar: i32, cigar: &[u32]) -> i32 {
    let mut l = 0_i32;
    for &op in cigar.iter().take(n_cigar.max(0) as usize) {
        let code = op & 0xf;
        if code == 0 || code == 2 {
            l += (op >> 4) as i32;
        }
    }
    l
}

#[doc = "Original function: _mm_realloc:1834"]
pub fn mm_realloc<T: Copy + Default>(ptr: &[T], csize: i64, nsize: i64, _dsize: i16) -> Vec<T> {
    if nsize <= csize {
        return ptr.to_vec();
    }
    let mut nptr = vec![T::default(); usize::try_from(nsize).expect("nsize")];
    let copy_len = usize::try_from(csize).expect("csize").min(ptr.len());
    nptr[..copy_len].copy_from_slice(&ptr[..copy_len]);
    nptr
}

#[doc = "Original function: bns_get_seq_v2:1851"]
#[inline]
pub fn bns_get_seq_v2<'a>(
    l_pac: i64,
    _pac: &[u8],
    mut beg: i64,
    mut end: i64,
    len: &mut i64,
    ref_string: &'a [u8],
    _seqb: Option<&'a mut [u8]>,
) -> Option<&'a [u8]> {
    if end < beg {
        std::mem::swap(&mut beg, &mut end);
    }
    if end > (l_pac << 1) {
        end = l_pac << 1;
    }
    if beg < 0 {
        beg = 0;
    }
    if !(beg >= l_pac || end <= l_pac) {
        *len = 0;
        return None;
    }

    *len = end - beg;
    let beg_usize = beg as usize;
    let end_usize = end as usize;
    let seq = ref_string.get(beg_usize..end_usize);
    assert!(seq.is_some(), "reference window outside ref_string");
    seq
}

#[doc = "Original function: bns_fetch_seq_v2:1890"]
#[inline]
pub fn bns_fetch_seq_v2<'a>(
    bns: &bntseq_t,
    pac: &[u8],
    beg: &mut i64,
    mid: i64,
    end: &mut i64,
    rid: &mut i32,
    ref_string: &'a [u8],
    seqb: Option<&'a mut [u8]>,
) -> &'a [u8] {
    if *end < *beg {
        std::mem::swap(beg, end);
    }
    assert!(*beg <= mid && mid < *end);

    let mut is_rev = 0;
    *rid = crate::generated::bntseq_cpp::bns_pos2rid(bns, bns_depos(bns, mid, &mut is_rev));
    let ann = &bns.anns[*rid as usize];
    let mut far_beg = ann.offset;
    let mut far_end = far_beg + i64::from(ann.len);
    if is_rev != 0 {
        let tmp = far_beg;
        far_beg = (bns.l_pac << 1) - far_end;
        far_end = (bns.l_pac << 1) - tmp;
    }
    *beg = (*beg).max(far_beg);
    *end = (*end).min(far_end);

    let mut len = 0_i64;
    let seq = bns_get_seq_v2(bns.l_pac, pac, *beg, *end, &mut len, ref_string, seqb);
    assert!(seq.is_some() && *end - *beg == len);
    seq.expect("sequence window")
}

#[doc = "Original function: sortPairsLenExt:1926"]
pub fn sortPairsLenExt(
    pairArray: &mut [SeqPair],
    count: i32,
    tempArray: &mut [SeqPair],
    hist: &mut [i32],
    numPairs128: &mut i32,
    numPairs16: &mut i32,
    numPairs1: &mut i32,
    score_a: i32,
) {
    *numPairs128 = 0;
    *numPairs16 = 0;
    *numPairs1 = 0;

    let count = count as usize;
    assert!(pairArray.len() >= count);
    assert!(tempArray.len() >= count);
    assert!(hist.len() > MAX_SEQ_LEN8 + MAX_SEQ_LEN16);

    hist[..=MAX_SEQ_LEN8 + MAX_SEQ_LEN16].fill(0);
    let (hist8, tail) = hist.split_at_mut(MAX_SEQ_LEN8);
    let (hist16, hist3) = tail.split_at_mut(MAX_SEQ_LEN16);
    #[cfg(debug_assertions)]
    let mut arr = vec![0_i32; count];

    let max_len8 = MAX_SEQ_LEN8 as i32;
    let max_len16 = MAX_SEQ_LEN16 as i32;
    for sp in pairArray.iter().take(count).copied() {
        let minval = sp.h0 + sp.len1.min(sp.len2) * score_a;
        if sp.len1 < max_len8 && sp.len2 < max_len8 && minval < max_len8 {
            hist8[minval as usize] += 1;
        } else if sp.len1 < max_len16 && sp.len2 < max_len16 && minval < max_len16 {
            hist16[minval as usize] += 1;
        } else {
            hist3[0] += 1;
        }
    }

    let mut cumul_sum = 0_i32;
    for item in hist8.iter_mut() {
        let cur = *item;
        *item = cumul_sum;
        cumul_sum += cur;
    }
    for item in hist16.iter_mut() {
        let cur = *item;
        *item = cumul_sum;
        cumul_sum += cur;
    }
    hist3[0] = cumul_sum;

    for (_i, sp) in pairArray.iter().take(count).copied().enumerate() {
        let minval = sp.h0 + sp.len1.min(sp.len2) * score_a;
        if sp.len1 < max_len8 && sp.len2 < max_len8 && minval < max_len8 {
            let mv_us = minval as usize;
            let pos = hist8[mv_us] as usize;
            tempArray[pos] = sp;
            hist8[mv_us] += 1;
            *numPairs128 += 1;
            #[cfg(debug_assertions)]
            {
                assert_eq!(arr[pos], 0, "sortPairsLenExt repeated 128 slot");
                arr[pos] = 1;
            }
        } else if sp.len1 < max_len16 && sp.len2 < max_len16 && minval < max_len16 {
            let mv_us = minval as usize;
            let pos = hist16[mv_us] as usize;
            tempArray[pos] = sp;
            hist16[mv_us] += 1;
            *numPairs16 += 1;
            #[cfg(debug_assertions)]
            {
                assert_eq!(arr[pos], 0, "sortPairsLenExt repeated 16 slot");
                arr[pos] = (_i + 1) as i32;
            }
        } else {
            let pos = hist3[0] as usize;
            tempArray[pos] = sp;
            hist3[0] += 1;
            #[cfg(debug_assertions)]
            {
                arr[pos] = (_i + 1) as i32;
            }
            *numPairs1 += 1;
        }
    }

    pairArray[..count].copy_from_slice(&tempArray[..count]);
}

#[doc = "Original function: sortPairsLen:2025"]
#[inline]
pub fn sortPairsLen(
    pairArray: &mut [SeqPair],
    count: i32,
    tempArray: &mut [SeqPair],
    hist: &mut [i32],
) {
    let count = count as usize;
    assert!(pairArray.len() >= count);
    assert!(tempArray.len() >= count);
    assert!(hist.len() > MAX_SEQ_LEN16);

    hist[..=MAX_SEQ_LEN16].fill(0);

    // First pass: validate sp.len1 fits and bump bucket counts. After this pass, sum of all
    // hist[..MAX_SEQ_LEN16+1] entries equals count, so the cumsum below produces positions in
    // [0, count). That bounds the unchecked tempArray[pos] writes in the third pass.
    let pair_ptr = pairArray.as_ptr();
    let temp_ptr = tempArray.as_mut_ptr();
    let hist_ptr = hist.as_mut_ptr();
    for sp in pairArray.iter().take(count).copied() {
        let idx = sp.len1 as usize;
        assert!(idx <= MAX_SEQ_LEN16);
        unsafe { *hist_ptr.add(idx) += 1 };
    }
    let mut cumul_sum = 0_i32;
    for item in hist.iter_mut().take(MAX_SEQ_LEN16 + 1) {
        let cur = *item;
        *item = cumul_sum;
        cumul_sum += cur;
    }
    for k in 0..count {
        let sp = unsafe { *pair_ptr.add(k) };
        let idx = sp.len1 as usize;
        let pos = unsafe { *hist_ptr.add(idx) } as usize;
        unsafe {
            *temp_ptr.add(pos) = sp;
            *hist_ptr.add(idx) += 1;
        }
    }
    pairArray[..count].copy_from_slice(&tempArray[..count]);
}

#[derive(Default)]
struct Chain2alnScratch {
    left_pairs: Vec<SeqPair>,
    right_pairs: Vec<SeqPair>,
    left_ref_buf: Vec<u8>,
    right_ref_buf: Vec<u8>,
    left_qer_buf: Vec<u8>,
    right_qer_buf: Vec<u8>,
    hist: Vec<i32>,
    temp_pairs: Vec<SeqPair>,
    lim: Vec<i32>,
    sorted_seed_ranges: Vec<Vec<(usize, usize)>>,
    sorted_seed_indices: Vec<u64>,
    work_pairs: Vec<SeqPair>,
}

thread_local! {
    static CHAIN2ALN_SCRATCH: std::cell::RefCell<Chain2alnScratch> =
        std::cell::RefCell::new(Chain2alnScratch::default());
    static MEM_COLLECT_QUERY_CUM_LEN: std::cell::RefCell<Vec<i32>> =
        const { std::cell::RefCell::new(Vec::new()) };
    static MEM_CHAIN_SEEDS_SA_COORD: std::cell::RefCell<Vec<i64>> =
        const { std::cell::RefCell::new(Vec::new()) };
    static USE_CURRENT_RAYON_POOL: Cell<bool> = const { Cell::new(false) };
}

#[doc = "Original function: mem_chain2aln_across_reads_V2:2069"]
pub fn mem_chain2aln_across_reads_V2(
    opt: &mem_opt_t,
    bns: &bntseq_t,
    pac: &[u8],
    seq_: &[bseq1_t],
    nseq: i32,
    chain_ar: &mut [mem_chain_v],
    av_v: &mut [mem_alnreg_v],
    _mmc: &mut mem_cache,
    ref_string: &[u8],
    _tid: i32,
) {
    let l_pac = bns.l_pac;
    // Take chain2aln scratch out of thread-local; restored at function end. Reuses capacity
    // across the ~11 chain2aln calls per program — these buffers grow to ~600KB each per call.
    let scratch = CHAIN2ALN_SCRATCH.with(|c| std::mem::take(&mut *c.borrow_mut()));
    let Chain2alnScratch {
        mut left_pairs,
        mut right_pairs,
        mut left_ref_buf,
        mut right_ref_buf,
        mut left_qer_buf,
        mut right_qer_buf,
        mut hist,
        mut temp_pairs,
        mut lim,
        mut sorted_seed_ranges,
        mut sorted_seed_indices,
        mut work_pairs,
    } = scratch;
    left_pairs.clear();
    right_pairs.clear();
    left_ref_buf.clear();
    right_ref_buf.clear();
    left_qer_buf.clear();
    right_qer_buf.clear();
    let nseq_usize = nseq as usize;
    sorted_seed_ranges.resize_with(nseq_usize, Vec::new);
    for per_read in sorted_seed_ranges.iter_mut().take(nseq_usize) {
        per_read.clear();
    }
    sorted_seed_indices.clear();

    for l in 0..nseq_usize {
        // Cow::Borrowed in production (seq_nt4 pre-populated by worker_bwt).
        let query_cow = nt4_cow(&seq_[l]);
        let query: &[u8] = &query_cow;
        let l_query = seq_[l].l_seq;
        let trace_this = debug_trace_read(seq_[l].name.as_deref());
        let chn = &mut chain_ar[l];
        let av = &mut av_v[l];
        av.n = 0;
        av.m = chn.a.iter().map(|c| c.n as usize).sum();
        // Reuse the previous Vec's capacity instead of dropping + reallocating per read.
        av.a.clear();
        av.a.resize(av.m, mem_alnreg_t::default());
        sorted_seed_ranges[l].clear();

        for (j, c) in chn.a.iter_mut().enumerate().take(chn.n) {
            if c.n == 0 {
                sorted_seed_ranges[l].push((sorted_seed_indices.len(), 0));
                continue;
            }

            let mut rmax0 = l_pac << 1;
            let mut rmax1 = 0_i64;
            let seed_order_start = sorted_seed_indices.len();
            let seed_count = c.n as usize; // c.n is i32 ≥ 0.
            sorted_seed_indices.extend((0..seed_count).map(|idx| {
                let score = c.seeds[idx].score as u32;
                (u64::from(score) << 32) | (idx as u64)
            }));
            sorted_seed_indices[seed_order_start..seed_order_start + seed_count].sort_unstable();
            sorted_seed_ranges[l].push((seed_order_start, seed_count));
            let seed_order = &sorted_seed_indices[seed_order_start..seed_order_start + seed_count];

            for t in c.seeds.iter().take(seed_count) {
                let b = t.rbeg - i64::from(t.qbeg + cal_max_gap(opt, t.qbeg));
                let e = t.rbeg
                    + i64::from(t.len)
                    + i64::from(
                        (l_query - t.qbeg - t.len) + cal_max_gap(opt, l_query - t.qbeg - t.len),
                    );
                rmax0 = rmax0.min(b);
                rmax1 = rmax1.max(e);
            }
            rmax0 = rmax0.max(0);
            rmax1 = rmax1.min(l_pac << 1);
            if rmax0 < l_pac && l_pac < rmax1 {
                if c.seeds[0].rbeg < l_pac {
                    rmax1 = l_pac;
                } else {
                    rmax0 = l_pac;
                }
            }

            let mut fetch_beg = rmax0;
            let mut fetch_end = rmax1;
            let mut rid = 0;
            let rseq = bns_fetch_seq_v2(
                bns,
                pac,
                &mut fetch_beg,
                c.seeds[0].rbeg,
                &mut fetch_end,
                &mut rid,
                ref_string,
                None,
            );
            assert_eq!(c.rid, rid);

            for &seed_key in seed_order.iter().rev() {
                let seed_idx = seed_key as u32;
                let s = &mut c.seeds[seed_idx as usize];
                let aln_idx = av.n;
                av.n += 1;
                let a = &mut av.a[aln_idx];
                s.aln = aln_idx as i32;
                if trace_this {
                    eprintln!(
                        "[trace::chain2aln:init] read={} chain={} aln={} seed_idx={} qbeg={} len={} rbeg={} rmax0={} rmax1={}",
                        seq_[l].name.as_deref().unwrap_or(""),
                        j,
                        aln_idx,
                        seed_idx,
                        s.qbeg,
                        s.len,
                        s.rbeg,
                        fetch_beg,
                        fetch_end
                    );
                }
                *a = mem_alnreg_t {
                    w: opt.w,
                    score: -1,
                    truesc: -1,
                    rid: c.rid,
                    frac_rep: c.frac_rep,
                    seedlen0: s.len,
                    c: j,
                    rb: i64::from(H0_),
                    re: i64::from(H0_),
                    qb: H0_,
                    qe: H0_,
                    ..Default::default()
                };

                if s.qbeg > 0 {
                    let sp = SeqPair {
                        h0: s.len * opt.a,
                        seqid: l as i32,
                        regid: aln_idx as i32,
                        idq: left_qer_buf.len() as i32,
                        idr: left_ref_buf.len() as i32,
                        len2: s.qbeg,
                        len1: (s.rbeg - fetch_beg) as i32,
                        ..Default::default()
                    };
                    let qbeg_us = s.qbeg as usize;
                    left_qer_buf.extend(query[..qbeg_us].iter().rev().copied());
                    let left_len = sp.len1 as usize;
                    left_ref_buf.extend(rseq[..left_len].iter().rev().copied());
                    left_pairs.push(sp);
                    a.qb = s.qbeg;
                    a.rb = s.rbeg;
                } else {
                    a.score = s.len * opt.a;
                    a.truesc = a.score;
                    a.qb = 0;
                    a.rb = s.rbeg;
                }

                if s.qbeg + s.len != l_query {
                    let qe = s.qbeg + s.len;
                    let re = (s.rbeg + i64::from(s.len) - fetch_beg) as usize;
                    let sp = SeqPair {
                        h0: H0_,
                        seqid: l as i32,
                        regid: aln_idx as i32,
                        len2: l_query - qe,
                        len1: (fetch_end - fetch_beg) as i32 - re as i32,
                        idq: right_qer_buf.len() as i32,
                        idr: right_ref_buf.len() as i32,
                        ..Default::default()
                    };
                    let qe_us = qe as usize;
                    let len2_us = sp.len2 as usize;
                    right_qer_buf.extend_from_slice(&query[qe_us..qe_us + len2_us]);
                    let len1_us = sp.len1 as usize;
                    right_ref_buf.extend_from_slice(&rseq[re..re + len1_us]);
                    right_pairs.push(sp);
                    a.qe = qe;
                    a.re = fetch_beg + re as i64;
                    if trace_this {
                        eprintln!(
                            "[trace::chain2aln:right_seed] read={} chain={} aln={} len1={} len2={} qe={} re={}",
                            seq_[l].name.as_deref().unwrap_or(""),
                            j,
                            aln_idx,
                            sp.len1,
                            sp.len2,
                            a.qe,
                            a.re
                        );
                    }
                } else {
                    a.qe = l_query;
                    a.re = s.rbeg + i64::from(s.len);
                    if a.rb != i64::from(H0_) && a.qb != H0_ {
                        recompute_seedcov_for_aln(c, a);
                    }
                }
            }
        }
    }

    hist.clear();
    hist.resize(MAX_SEQ_LEN8 + MAX_SEQ_LEN16 + 1, 0);
    let max_pairs = left_pairs.len().max(right_pairs.len()).max(1);
    temp_pairs.clear();
    temp_pairs.resize(max_pairs, SeqPair::default());
    let bswLeft = BandedPairWiseSW::ctor(
        opt.o_del,
        opt.e_del,
        opt.o_ins,
        opt.e_ins,
        opt.zdrop,
        opt.pen_clip5,
        &opt.mat,
        opt.a,
        opt.b,
        1,
    );
    let bswRight = BandedPairWiseSW::ctor(
        opt.o_del,
        opt.e_del,
        opt.o_ins,
        opt.e_ins,
        opt.zdrop,
        opt.pen_clip3,
        &opt.mat,
        opt.a,
        opt.b,
        1,
    );

    if !left_pairs.is_empty() {
        let left_len = left_pairs.len();
        let mut num128 = 0;
        let mut num16 = 0;
        let mut num1 = 0;
        sortPairsLenExt(
            &mut left_pairs,
            left_len as i32,
            &mut temp_pairs[..left_len],
            &mut hist,
            &mut num128,
            &mut num16,
            &mut num1,
            opt.a,
        );
        // Bucket scalar (num1): pairs at left_pairs[num128 + num16 ..]. Clone only that slice
        // (was: full left_pairs.clone() — wasted memory by carrying the unused larger-pair prefix).
        let scalar_off = (num128 + num16) as usize;
        let scalar_len = num1 as usize;
        work_pairs.clear();
        work_pairs.extend_from_slice(&left_pairs[scalar_off..scalar_off + scalar_len]);
        let mut nump = work_pairs.len();
        for i in 0..MAX_BAND_TRY {
            let w = opt.w << i;
            bswLeft.scalarBandedSWAWrapper(
                &mut work_pairs[..nump],
                &left_ref_buf,
                &left_qer_buf,
                nump as i32,
                1,
                w,
            );
            let mut keep = 0_usize;
            for idx in 0..nump {
                let sp = work_pairs[idx];
                let a = &mut av_v[sp.seqid as usize].a
                    [sp.regid as usize];
                let prev = a.score;
                a.score = sp.score;
                if a.score == prev || sp.max_off < (w >> 1) + (w >> 2) || i + 1 == MAX_BAND_TRY {
                    if sp.gscore <= 0 || sp.gscore <= a.score - opt.pen_clip5 {
                        a.qb -= sp.qle;
                        a.rb -= i64::from(sp.tle);
                        a.truesc = a.score;
                    } else {
                        a.qb = 0;
                        a.rb -= i64::from(sp.gtle);
                        a.truesc = sp.gscore;
                    }
                    a.w = a.w.max(w);
                    let chain = &chain_ar[sp.seqid as usize].a[a.c];
                    recompute_seedcov_for_aln(chain, a);
                    if debug_trace_read(
                        seq_[sp.seqid as usize]
                            .name
                            .as_deref(),
                    ) {
                        eprintln!(
                            "[trace::chain2aln:left_done] read={} aln={} score={} qle={} tle={} gtle={} gscore={} qb={} rb={} qe={} re={}",
                            seq_[sp.seqid as usize].name.as_deref().unwrap_or(""),
                            sp.regid,
                            a.score,
                            sp.qle,
                            sp.tle,
                            sp.gtle,
                            sp.gscore,
                            a.qb,
                            a.rb,
                            a.qe,
                            a.re
                        );
                    }
                } else {
                    work_pairs[keep] = sp;
                    keep += 1;
                }
            }
            nump = keep;
            if nump == 0 {
                break;
            }
        }

        // Bucket i16 (num16): pairs at left_pairs[num128 .. num128 + num16].
        let i16_off = num128 as usize;
        let i16_len = num16 as usize;
        work_pairs.clear();
        work_pairs.extend_from_slice(&left_pairs[i16_off..i16_off + i16_len]);
        let mut nump = work_pairs.len();
        for i in 0..MAX_BAND_TRY {
            let w = opt.w << i;
            sortPairsLen(
                &mut work_pairs[..nump],
                nump as i32,
                &mut temp_pairs[..nump],
                &mut hist,
            );
            bswLeft.getScores16(
                &mut work_pairs[..nump],
                &left_ref_buf,
                &left_qer_buf,
                nump as i32,
                1,
                w,
            );
            let mut keep = 0_usize;
            for idx in 0..nump {
                let sp = work_pairs[idx];
                let a = &mut av_v[sp.seqid as usize].a
                    [sp.regid as usize];
                let prev = a.score;
                a.score = sp.score;
                if a.score == prev || sp.max_off < (w >> 1) + (w >> 2) || i + 1 == MAX_BAND_TRY {
                    if sp.gscore <= 0 || sp.gscore <= a.score - opt.pen_clip5 {
                        a.qb -= sp.qle;
                        a.rb -= i64::from(sp.tle);
                        a.truesc = a.score;
                    } else {
                        a.qb = 0;
                        a.rb -= i64::from(sp.gtle);
                        a.truesc = sp.gscore;
                    }
                    a.w = a.w.max(w);
                    let chain = &chain_ar[sp.seqid as usize].a[a.c];
                    recompute_seedcov_for_aln(chain, a);
                } else {
                    work_pairs[keep] = sp;
                    keep += 1;
                }
            }
            nump = keep;
            if nump == 0 {
                break;
            }
        }

        // Bucket u8 (num128): pairs at left_pairs[0 .. num128].
        let u8_len = num128 as usize;
        work_pairs.clear();
        work_pairs.extend_from_slice(&left_pairs[..u8_len]);
        let mut nump = work_pairs.len();
        for i in 0..MAX_BAND_TRY {
            let w = opt.w << i;
            sortPairsLen(
                &mut work_pairs[..nump],
                nump as i32,
                &mut temp_pairs[..nump],
                &mut hist,
            );
            bswLeft.getScores8(
                &mut work_pairs[..nump],
                &left_ref_buf,
                &left_qer_buf,
                nump as i32,
                1,
                w,
            );
            let mut keep = 0_usize;
            for idx in 0..nump {
                let sp = work_pairs[idx];
                let a = &mut av_v[sp.seqid as usize].a
                    [sp.regid as usize];
                let prev = a.score;
                a.score = sp.score;
                if a.score == prev || sp.max_off < (w >> 1) + (w >> 2) || i + 1 == MAX_BAND_TRY {
                    if sp.gscore <= 0 || sp.gscore <= a.score - opt.pen_clip5 {
                        a.qb -= sp.qle;
                        a.rb -= i64::from(sp.tle);
                        a.truesc = a.score;
                    } else {
                        a.qb = 0;
                        a.rb -= i64::from(sp.gtle);
                        a.truesc = sp.gscore;
                    }
                    a.w = a.w.max(w);
                    let chain = &chain_ar[sp.seqid as usize].a[a.c];
                    recompute_seedcov_for_aln(chain, a);
                } else {
                    work_pairs[keep] = sp;
                    keep += 1;
                }
            }
            nump = keep;
            if nump == 0 {
                break;
            }
        }
    }

    for sp in &mut right_pairs {
        let a = &av_v[sp.seqid as usize].a
            [sp.regid as usize];
        sp.h0 = a.score;
    }
    if !right_pairs.is_empty() {
        let right_len = right_pairs.len();
        let mut num128 = 0;
        let mut num16 = 0;
        let mut num1 = 0;
        sortPairsLenExt(
            &mut right_pairs,
            right_len as i32,
            &mut temp_pairs[..right_len],
            &mut hist,
            &mut num128,
            &mut num16,
            &mut num1,
            opt.a,
        );

        let mut process_right_bucket = |offset: usize, initial: usize, dispatch: u8| {
            // dispatch: 0 = scalar (num1 bucket), 1 = i16 SIMD (getScores16), 2 = u8 SIMD (getScores8)
            work_pairs.clear();
            work_pairs.extend_from_slice(&right_pairs[offset..offset + initial]);
            let mut nump = work_pairs.len();
            for i in 0..MAX_BAND_TRY {
                let w = opt.w << i;
                if dispatch != 0 {
                    sortPairsLen(
                        &mut work_pairs[..nump],
                        nump as i32,
                        &mut temp_pairs[..nump],
                        &mut hist,
                    );
                }
                match dispatch {
                    0 => bswRight.scalarBandedSWAWrapper(
                        &mut work_pairs[..nump],
                        &right_ref_buf,
                        &right_qer_buf,
                        nump as i32,
                        1,
                        w,
                    ),
                    1 => bswRight.getScores16(
                        &mut work_pairs[..nump],
                        &right_ref_buf,
                        &right_qer_buf,
                        nump as i32,
                        1,
                        w,
                    ),
                    2 => bswRight.getScores8(
                        &mut work_pairs[..nump],
                        &right_ref_buf,
                        &right_qer_buf,
                        nump as i32,
                        1,
                        w,
                    ),
                    _ => unreachable!(),
                }
                let mut keep = 0_usize;
                for idx in 0..nump {
                    let sp = work_pairs[idx];
                    let a = &mut av_v[sp.seqid as usize].a
                        [sp.regid as usize];
                    let prev = a.score;
                    a.score = sp.score;
                    if a.score == prev || sp.max_off < (w >> 1) + (w >> 2) || i + 1 == MAX_BAND_TRY
                    {
                        if sp.gscore <= 0 || sp.gscore <= a.score - opt.pen_clip3 {
                            a.qe += sp.qle;
                            a.re += i64::from(sp.tle);
                            a.truesc += a.score - sp.h0;
                        } else {
                            a.qe = seq_[sp.seqid as usize].l_seq;
                            a.re += i64::from(sp.gtle);
                            a.truesc += sp.gscore - sp.h0;
                        }
                        a.w = a.w.max(w);
                        let chain = &chain_ar[sp.seqid as usize].a[a.c];
                        recompute_seedcov_for_aln(chain, a);
                        if debug_trace_read(
                            seq_[sp.seqid as usize]
                                .name
                                .as_deref(),
                        ) {
                            eprintln!(
                            "[trace::chain2aln:right_done] read={} aln={} score={} h0={} qle={} tle={} gtle={} gscore={} qb={} rb={} qe={} re={}",
                            seq_[sp.seqid as usize].name.as_deref().unwrap_or(""),
                            sp.regid,
                            a.score,
                            sp.h0,
                            sp.qle,
                            sp.tle,
                            sp.gtle,
                            sp.gscore,
                            a.qb,
                            a.rb,
                            a.qe,
                            a.re
                        );
                        }
                    } else {
                        work_pairs[keep] = sp;
                        keep += 1;
                    }
                }
                nump = keep;
                if nump == 0 {
                    break;
                }
            }
        };

        process_right_bucket((num128 + num16) as usize, num1 as usize, 0);
        process_right_bucket(num128 as usize, num16 as usize, 1);
        process_right_bucket(0, num128 as usize, 2);
    }

    lim.clear();
    lim.resize(nseq_usize, 0);
    for l in 0..nseq_usize {
        let l_query = seq_[l].l_seq;
        let chn = &chain_ar[l];
        let av = &mut av_v[l];
        for (j, c) in chn.a.iter().enumerate().take(chn.n) {
            let (seed_order_start, seed_count) = sorted_seed_ranges[l][j];
            let srt2 = &mut sorted_seed_indices[seed_order_start..seed_order_start + seed_count];
            for k in (0..srt2.len()).rev() {
                let seed_key = srt2[k];
                if seed_key == u64::MAX {
                    continue;
                }
                let seed_idx = seed_key as u32;
                let s = &c.seeds[seed_idx as usize];
                let mut v = 0_i32;
                for p in av.a.iter().take(av.n) {
                    if v >= lim[l] {
                        break;
                    }
                    if p.qb == -1 && p.qe == -1 {
                        continue;
                    }
                    if s.rbeg < p.rb
                        || s.rbeg + i64::from(s.len) > p.re
                        || s.qbeg < p.qb
                        || s.qbeg + s.len > p.qe
                    {
                        v += 1;
                        continue;
                    }
                    if f64::from(s.len - p.seedlen0) > 0.1 * f64::from(l_query) {
                        v += 1;
                        continue;
                    }
                    let mut qd = s.qbeg - p.qb;
                    let mut rd = (s.rbeg - p.rb) as i32;
                    let mut max_gap = cal_max_gap(opt, qd.min(rd));
                    let mut w = max_gap.min(p.w);
                    if qd - rd < w && rd - qd < w {
                        break;
                    }
                    qd = p.qe - (s.qbeg + s.len);
                    rd = (p.re - (s.rbeg + i64::from(s.len))) as i32;
                    max_gap = cal_max_gap(opt, qd.min(rd));
                    w = max_gap.min(p.w);
                    if qd - rd < w && rd - qd < w {
                        break;
                    }
                    v += 1;
                }
                if v < lim[l] {
                    let mut overlap = false;
                    for v_idx in k + 1..srt2.len() {
                        let next_seed_key = srt2[v_idx];
                        if next_seed_key == u64::MAX {
                            continue;
                        }
                        let next_seed_idx = next_seed_key as u32;
                        let t = &c.seeds[next_seed_idx as usize];
                        if f64::from(t.len) < f64::from(s.len) * 0.95 {
                            continue;
                        }
                        if s.qbeg <= t.qbeg
                            && s.qbeg + s.len - t.qbeg >= s.len >> 2
                            && t.qbeg - s.qbeg != (t.rbeg - s.rbeg) as i32
                        {
                            overlap = true;
                            break;
                        }
                        if t.qbeg <= s.qbeg
                            && t.qbeg + t.len - s.qbeg >= s.len >> 2
                            && s.qbeg - t.qbeg != (s.rbeg - t.rbeg) as i32
                        {
                            overlap = true;
                            break;
                        }
                    }
                    if !overlap {
                        let aln_idx = s.aln as usize;
                        if debug_trace_read(seq_[l].name.as_deref()) {
                            eprintln!(
                                "[trace::chain2aln:prune] read={} chain={} seed_idx={} seed_aln={} lim={} seed=({},{},{}) reg_before={:?}",
                                seq_[l].name.as_deref().unwrap_or(""),
                                j,
                                seed_idx,
                                s.aln,
                                lim[l],
                                s.qbeg,
                                s.len,
                                s.rbeg,
                                av.a.get(usize::try_from(s.aln).expect("s.aln")).map(|r| (r.score, r.qb, r.qe, r.rb, r.re, r.seedlen0)),
                            );
                        }
                        let ar = &mut av.a[aln_idx];
                        ar.qb = -1;
                        ar.qe = -1;
                        srt2[k] = u64::MAX;
                        continue;
                    }
                }
                if debug_trace_read(seq_[l].name.as_deref()) {
                    eprintln!(
                        "[trace::chain2aln:keep] read={} chain={} seed_idx={} seed_aln={} lim={} seed=({},{},{})",
                        seq_[l].name.as_deref().unwrap_or(""),
                        j,
                        seed_idx,
                        s.aln,
                        lim[l],
                        s.qbeg,
                        s.len,
                        s.rbeg,
                    );
                }
                lim[l] += 1;
            }
        }
    }

    // Restore the scratch for the next chain2aln call's reuse.
    CHAIN2ALN_SCRATCH.with(|c| {
        *c.borrow_mut() = Chain2alnScratch {
            left_pairs,
            right_pairs,
            left_ref_buf,
            right_ref_buf,
            left_qer_buf,
            right_qer_buf,
            hist,
            temp_pairs,
            lim,
            sorted_seed_ranges,
            sorted_seed_indices,
            work_pairs,
        };
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::generated::bntseq_cpp::bns_fasta2bntseq;
    use crate::generated::bntseq_h::{bntann1_t, bntseq_t};
    use crate::generated::bwa_cpp::bseq_read_orig;
    use crate::generated::bwa_h::bseq1_t;
    use crate::generated::fastmap_cpp::memoryAlloc;
    use crate::generated::fastmap_h::ktp_aux_t;
    use crate::generated::fmi_search_cpp::FMI_search;
    use crate::generated::kseq_h::kseq_t;
    use crate::generated::macro_h::BATCH_SIZE;
    use std::fs;
    use std::io::Cursor;
    use std::path::Path;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn revcomp_nt4(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|&b| match b {
                0..=3 => 3 - b,
                _ => b,
            })
            .collect()
    }

    fn pack_seq(seq: &[u8]) -> Vec<u8> {
        let mut pac = vec![0_u8; (seq.len() + 3) / 4];
        for (i, &base) in seq.iter().enumerate() {
            let shift = (((!i64::try_from(i).expect("i")) & 3) << 1) as u8;
            pac[i >> 2] |= base << shift;
        }
        pac
    }

    fn single_seed_chain(seed: mem_seed_t, rid: i32) -> mem_chain_t {
        mem_chain_t {
            n: 1,
            m: 1,
            rid,
            seeds: vec![seed],
            ..Default::default()
        }
    }

    fn infer_real_fixture_pes(n_threads: i32) -> [mem_pestat_t; 4] {
        let prefix = ".tmp/real_bench/ecoli_rel606";
        let r1 = "external/tutorial_bwa_small/SRR2584863_1.trim.sub.fastq";
        let r2 = "external/tutorial_bwa_small/SRR2584863_2.trim.sub.fastq";
        assert!(
            Path::new(&format!("{prefix}.bwt.2bit.64")).exists(),
            "missing real-fixture index at {prefix}.bwt.2bit.64"
        );
        assert!(Path::new(r1).exists(), "missing real-fixture reads at {r1}");
        assert!(Path::new(r2).exists(), "missing real-fixture reads at {r2}");

        let r1_text = fs::read_to_string(r1).expect("read r1");
        let r2_text = fs::read_to_string(r2).expect("read r2");
        let mut ks1 = kseq_t::from_text(&r1_text);
        let mut ks2 = kseq_t::from_text(&r2_text);
        let mut n = 0_i32;
        let mut total_bases = 0_i64;
        let mut seqs = bseq_read_orig(1 << 20, &mut n, &mut ks1, Some(&mut ks2), &mut total_bases);
        assert!(n > 0);

        let mut fmi = FMI_search::ctor(prefix);
        fmi.load_index();

        let mut opt = (*mem_opt_init()).clone();
        opt.flag |= MEM_F_PE;
        opt.n_threads = n_threads;

        let aux = ktp_aux_t {
            opt: Some(Box::new(opt.clone())),
            ..Default::default()
        };
        let mut worker = worker_t {
            fmi: Some(fmi),
            ..Default::default()
        };
        memoryAlloc(&aux, &mut worker, n, n_threads);
        worker.opt = Some(Box::new(opt.clone()));
        worker.n_processed = 0;
        worker.nreads = n;
        worker.seqs = std::mem::take(&mut seqs);
        let l_pac = {
            let fmi_ref = worker.fmi.as_ref().expect("worker fmi");
            let bns = fmi_ref.base.idx.bns.as_ref().expect("bns");
            worker.ref_string = pac_to_reference_layout(bns.l_pac, &fmi_ref.base.idx.pac);
            bns.l_pac
        };

        run_worker_chunks_parallel(&mut worker, n, worker_bwt);
        run_worker_chunks_parallel(&mut worker, n, worker_aln);
        mem_pestat(&opt, l_pac, n, &worker.regs, &mut worker.pes);
        worker.pes
    }

    fn load_real_fixture_pair(target: &str) -> (Vec<bseq1_t>, FMI_search, mem_opt_t) {
        let prefix = ".tmp/real_bench/ecoli_rel606";
        let r1 = "external/tutorial_bwa_small/SRR2584863_1.trim.sub.fastq";
        let r2 = "external/tutorial_bwa_small/SRR2584863_2.trim.sub.fastq";
        assert!(
            Path::new(&format!("{prefix}.bwt.2bit.64")).exists(),
            "missing real-fixture index at {prefix}.bwt.2bit.64"
        );
        assert!(Path::new(r1).exists(), "missing real-fixture reads at {r1}");
        assert!(Path::new(r2).exists(), "missing real-fixture reads at {r2}");

        let r1_text = fs::read_to_string(r1).expect("read r1");
        let r2_text = fs::read_to_string(r2).expect("read r2");
        let mut ks1 = kseq_t::from_text(&r1_text);
        let mut ks2 = kseq_t::from_text(&r2_text);
        let mut total_bases = 0_i64;
        let found_pair = loop {
            let mut n = 0_i32;
            let seqs = bseq_read_orig(1 << 20, &mut n, &mut ks1, Some(&mut ks2), &mut total_bases);
            assert!(n > 0, "target {target} not found in fixture");
            if let Some(idx) = seqs.iter().position(|s| s.name.as_deref() == Some(target)) {
                let pair_start = idx & !1;
                break vec![seqs[pair_start].clone(), seqs[pair_start + 1].clone()];
            }
        };

        let mut fmi = FMI_search::ctor(prefix);
        fmi.load_index();

        let mut opt = (*mem_opt_init()).clone();
        opt.flag |= MEM_F_PE;
        opt.n_threads = 2;

        (found_pair, fmi, opt)
    }

    #[test]
    fn cal_max_gap_matches_original_formula_and_caps_to_band() {
        let mut opt = mem_opt_t::default();
        opt.a = 1;
        opt.o_del = 6;
        opt.e_del = 1;
        opt.o_ins = 6;
        opt.e_ins = 1;
        opt.w = 100;
        assert_eq!(cal_max_gap(&opt, 10), 5);
        assert_eq!(cal_max_gap(&opt, 1000), 200);
    }

    #[test]
    #[ignore]
    fn real_fixture_pestat_matches_across_threads() {
        let pes1 = infer_real_fixture_pes(1);
        let pes2 = infer_real_fixture_pes(2);
        eprintln!("pes_t1={pes1:?}");
        eprintln!("pes_t2={pes2:?}");
        for i in 0..4 {
            assert_eq!(pes1[i].low, pes2[i].low, "orientation {i} low");
            assert_eq!(pes1[i].high, pes2[i].high, "orientation {i} high");
            assert_eq!(pes1[i].failed, pes2[i].failed, "orientation {i} failed");
            assert!(
                (pes1[i].avg - pes2[i].avg).abs() < 1e-9,
                "orientation {i} avg"
            );
            assert!(
                (pes1[i].std - pes2[i].std).abs() < 1e-9,
                "orientation {i} std"
            );
        }
    }

    #[test]
    #[ignore]
    fn real_fixture_106000_keeps_secondary_after_kernel2() {
        let (seqs, fmi, opt) = load_real_fixture_pair("SRR2584863.106000");
        let aux = ktp_aux_t {
            opt: Some(Box::new(opt.clone())),
            ..Default::default()
        };
        let mut worker = worker_t {
            fmi: Some(fmi),
            ..Default::default()
        };
        memoryAlloc(&aux, &mut worker, 2, 2);
        worker.opt = Some(Box::new(opt.clone()));
        worker.nreads = 2;
        worker.seqs = seqs;

        {
            let fmi_ref = worker.fmi.as_ref().expect("worker fmi");
            let bns = fmi_ref.base.idx.bns.as_ref().expect("bns");
            worker.ref_string = pac_to_reference_layout(bns.l_pac, &fmi_ref.base.idx.pac);
        }

        worker_bwt(&mut worker, 0, 2, 0);
        worker_aln(&mut worker, 0, 2, 0);

        let regs = &worker.regs[0].a[..worker.regs[0].n];
        assert!(
            regs.iter()
                .any(|r| r.score == 150 && r.qb == 0 && r.qe == 150),
            "missing primary in {:?}",
            regs.iter()
                .map(|r| (r.score, r.qb, r.qe, r.rb, r.re))
                .collect::<Vec<_>>()
        );
        assert!(
            regs.iter()
                .any(|r| r.score == 21 && r.qb == 120 && r.qe == 141),
            "missing expected secondary in {:?}",
            regs.iter()
                .map(|r| (r.score, r.qb, r.qe, r.rb, r.re))
                .collect::<Vec<_>>()
        );
    }

    #[test]
    #[ignore]
    fn real_fixture_207026_kernel2_scores_match_expected_shape() {
        let (seqs, fmi, opt) = load_real_fixture_pair("SRR2584863.207026");
        let aux = ktp_aux_t {
            opt: Some(Box::new(opt.clone())),
            ..Default::default()
        };
        let mut worker = worker_t {
            fmi: Some(fmi),
            ..Default::default()
        };
        memoryAlloc(&aux, &mut worker, 2, 2);
        worker.opt = Some(Box::new(opt.clone()));
        worker.nreads = 2;
        worker.seqs = seqs;

        {
            let fmi_ref = worker.fmi.as_ref().expect("worker fmi");
            let bns = fmi_ref.base.idx.bns.as_ref().expect("bns");
            worker.ref_string = pac_to_reference_layout(bns.l_pac, &fmi_ref.base.idx.pac);
        }

        worker_bwt(&mut worker, 0, 2, 0);
        worker_aln(&mut worker, 0, 2, 0);

        let regs = &worker.regs[0].a[..worker.regs[0].n];
        let scores: Vec<i32> = regs.iter().map(|r| r.score).collect();
        assert!(scores.contains(&150), "missing primary: {scores:?}");
        assert!(scores.contains(&108), "missing 108 alternative: {scores:?}");
        assert!(
            !scores.contains(&126),
            "unexpected 126 alternative survived kernel2: {:?}",
            regs.iter()
                .map(|r| (r.score, r.qb, r.qe, r.rb, r.re, r.seedcov))
                .collect::<Vec<_>>()
        );
    }

    #[test]
    #[ignore]
    fn real_fixture_207026_patch_reg_for_126_vs_150() {
        let (seqs, fmi, opt) = load_real_fixture_pair("SRR2584863.207026");
        let aux = ktp_aux_t {
            opt: Some(Box::new(opt.clone())),
            ..Default::default()
        };
        let mut worker = worker_t {
            fmi: Some(fmi),
            ..Default::default()
        };
        memoryAlloc(&aux, &mut worker, 2, 2);
        worker.opt = Some(Box::new(opt.clone()));
        worker.nreads = 2;
        worker.seqs = seqs;

        worker_bwt(&mut worker, 0, 2, 0);
        let (bns, pac) = {
            let fmi_ref = worker.fmi.as_ref().expect("worker fmi");
            (
                fmi_ref.base.idx.bns.as_ref().expect("bns"),
                fmi_ref.base.idx.pac.clone(),
            )
        };
        let regs = vec![
            mem_alnreg_t {
                score: 126,
                qb: 0,
                qe: 150,
                rb: 3_319_983,
                re: 3_320_134,
                ..Default::default()
            },
            mem_alnreg_t {
                score: 150,
                qb: 0,
                qe: 150,
                rb: 3_320_075,
                re: 3_320_225,
                ..Default::default()
            },
        ];
        let query_text = worker.seqs[0].seq.as_deref().expect("seq");
        let mut query = seq_to_nt4(query_text);
        let mut w = 0_i32;
        let score = mem_patch_reg(
            &opt,
            Some(bns),
            Some(&pac),
            Some(&mut query),
            &regs[0],
            &regs[1],
            &mut w,
        );
        eprintln!("patch_score={score} w={w}");
        assert!(
            score <= 0,
            "unexpectedly patchable in Rust with score={score} w={w}"
        );
    }

    #[test]
    #[ignore]
    fn real_fixture_207026_chain_shape_after_worker_bwt() {
        let (seqs, fmi, opt) = load_real_fixture_pair("SRR2584863.207026");
        let aux = ktp_aux_t {
            opt: Some(Box::new(opt.clone())),
            ..Default::default()
        };
        let mut worker = worker_t {
            fmi: Some(fmi),
            ..Default::default()
        };
        memoryAlloc(&aux, &mut worker, 2, 2);
        worker.opt = Some(Box::new(opt.clone()));
        worker.nreads = 2;
        worker.seqs = seqs;
        worker_bwt(&mut worker, 0, 2, 0);

        let chains = &worker.chain_ar[0].a[..worker.chain_ar[0].n];
        eprintln!(
            "chains={:?}",
            chains
                .iter()
                .map(|c| {
                    (
                        c.pos,
                        c.w,
                        c.kept,
                        c.n,
                        c.seeds
                            .iter()
                            .take(usize::try_from(c.n).expect("c.n"))
                            .map(|s| (s.qbeg, s.len, s.rbeg, s.score))
                            .collect::<Vec<_>>(),
                    )
                })
                .collect::<Vec<_>>()
        );
    }

    #[test]
    #[ignore]
    fn real_fixture_207026_chain_counts_by_kernel1_stage() {
        let (seqs, fmi, opt) = load_real_fixture_pair("SRR2584863.207026");
        let mut match_array = Vec::new();
        let mut min_intv_ar = Vec::new();
        let mut query_pos_ar = Vec::new();
        let mut enc_qdb = Vec::new();
        let mut rid = Vec::new();
        let mut mmc = mem_cache::default();
        let mut num_smem = 0_i64;
        let mut seed_buf = vec![mem_seed_t::default(); 4096];
        let seed_buf_len = seed_buf.len() as i64;
        let mut chain_ar = vec![mem_chain_v::default(); 2];

        mem_collect_smem(
            &fmi,
            &opt,
            &seqs,
            2,
            &mut match_array,
            &mut min_intv_ar,
            &mut query_pos_ar,
            &mut enc_qdb,
            &mut rid,
            &mut mmc,
            &mut num_smem,
            0,
        );

        mem_chain_seeds(
            &fmi,
            &opt,
            fmi.base.idx.bns.as_ref().expect("bns"),
            &seqs,
            2,
            0,
            &mut chain_ar,
            &mut seed_buf,
            seed_buf_len,
            &match_array,
            num_smem,
        );
        eprintln!(
            "after_chain_seeds={:?}",
            chain_ar[0]
                .a
                .iter()
                .map(|c| (c.pos, c.w, c.n))
                .collect::<Vec<_>>()
        );

        chain_ar[0].n = usize::try_from(mem_chain_flt(
            &opt,
            i32::try_from(chain_ar[0].n).expect("n"),
            &mut chain_ar[0].a,
            0,
        ))
        .expect("filtered n");
        eprintln!(
            "after_chain_flt={:?}",
            chain_ar[0]
                .a
                .iter()
                .take(chain_ar[0].n)
                .map(|c| (c.pos, c.w, c.kept, c.n))
                .collect::<Vec<_>>()
        );

        mem_flt_chained_seeds(
            &opt,
            fmi.base.idx.bns.as_ref().expect("bns"),
            &fmi.base.idx.pac,
            &seqs,
            i32::try_from(chain_ar[0].n).expect("n"),
            &mut chain_ar[0].a,
        );
        eprintln!(
            "after_flt_chained={:?}",
            chain_ar[0]
                .a
                .iter()
                .take(chain_ar[0].n)
                .map(|c| {
                    (
                        c.pos,
                        c.w,
                        c.n,
                        c.seeds
                            .iter()
                            .take(usize::try_from(c.n).expect("c.n"))
                            .map(|s| (s.qbeg, s.len, s.rbeg, s.score))
                            .collect::<Vec<_>>(),
                    )
                })
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn mm_realloc_grows_and_preserves_prefix() {
        let grown = mm_realloc(&[1_u8, 2, 3], 3, 6, 1);
        assert_eq!(&grown[..3], &[1, 2, 3]);
        assert_eq!(grown.len(), 6);
    }

    #[test]
    fn mm_realloc_does_not_shrink() {
        let original = [4_u8, 5, 6];
        let same = mm_realloc(&original, 3, 2, 1);
        assert_eq!(same, original);
    }

    #[test]
    fn smem_aux_init_allocates_batch_size_empty_vectors() {
        let aux = smem_aux_init();
        assert_eq!(aux.len(), BATCH_SIZE);
        assert!(aux
            .iter()
            .all(|x| x.mem.a.is_empty() && x.mem1.a.is_empty()));
        assert!(aux
            .iter()
            .all(|x| x.tmpv[0].a.is_empty() && x.tmpv[1].a.is_empty()));
    }

    #[test]
    fn smem_aux_destroy_clears_owned_memory() {
        let mut aux = smem_aux_init();
        aux[0].mem.a.push(Default::default());
        aux[0].tmpv[0].a.push(Default::default());
        smem_aux_destroy(&mut aux);
        assert!(aux.is_empty());
    }

    #[test]
    fn mem_opt_init_sets_expected_defaults() {
        let opt = mem_opt_init();
        assert_eq!(opt.a, 1);
        assert_eq!(opt.b, 4);
        assert_eq!(opt.w, 100);
        assert_eq!(opt.min_seed_len, 19);
        assert_eq!(opt.max_occ, 500);
        assert_eq!(opt.n_threads, 1);
        assert_eq!(opt.max_XA_hits, 5);
        assert_eq!(opt.max_XA_hits_alt, 200);
        assert_eq!(opt.mat[0], 1);
        assert_eq!(opt.mat[1], -4);
        assert_eq!(opt.mat[4], -1);
        assert_eq!(opt.mat[24], -1);
    }

    #[test]
    fn sort_alnreg_re_orders_by_reference_end() {
        let mut regs = vec![
            mem_alnreg_t {
                re: 30,
                ..Default::default()
            },
            mem_alnreg_t {
                re: 10,
                ..Default::default()
            },
            mem_alnreg_t {
                re: 20,
                ..Default::default()
            },
        ];
        sort_alnreg_re(regs.len() as i32, &mut regs);
        assert_eq!(
            regs.iter().map(|x| x.re).collect::<Vec<_>>(),
            vec![10, 20, 30]
        );
    }

    #[test]
    fn sort_alnreg_score_orders_by_score_then_rb_then_qb() {
        let mut regs = vec![
            mem_alnreg_t {
                score: 20,
                rb: 50,
                qb: 10,
                ..Default::default()
            },
            mem_alnreg_t {
                score: 30,
                rb: 40,
                qb: 20,
                ..Default::default()
            },
            mem_alnreg_t {
                score: 30,
                rb: 40,
                qb: 10,
                ..Default::default()
            },
            mem_alnreg_t {
                score: 30,
                rb: 30,
                qb: 30,
                ..Default::default()
            },
        ];
        sort_alnreg_score(regs.len() as i32, &mut regs);
        assert_eq!(
            regs.iter()
                .map(|x| (x.score, x.rb, x.qb))
                .collect::<Vec<_>>(),
            vec![(30, 30, 30), (30, 40, 10), (30, 40, 20), (20, 50, 10)]
        );
    }

    #[test]
    fn sort_alnreg_helpers_handle_larger_partitions() {
        let mut regs: Vec<mem_alnreg_t> = (0..80)
            .rev()
            .map(|i| mem_alnreg_t {
                rid: 0,
                rb: i64::from((i * 37) % 80),
                re: i64::from((i * 53) % 80),
                qb: (i * 11) % 17,
                score: (i * 7) % 23,
                ..Default::default()
            })
            .collect();

        sort_alnreg_re(regs.len() as i32, &mut regs);
        for i in 1..regs.len() {
            assert!(regs[i - 1].re <= regs[i].re);
        }

        sort_alnreg_score(regs.len() as i32, &mut regs);
        for i in 1..regs.len() {
            assert!(
                regs[i - 1].score > regs[i].score
                    || (regs[i - 1].score == regs[i].score
                        && (regs[i - 1].rb < regs[i].rb
                            || (regs[i - 1].rb == regs[i].rb && regs[i - 1].qb <= regs[i].qb)))
            );
        }
    }

    #[test]
    fn mem_patch_reg_rejects_cross_strand_hits() {
        let opt = mem_opt_init();
        let bns = bntseq_t {
            l_pac: 8,
            anns: vec![bntann1_t {
                name: "chr1".into(),
                len: 8,
                ..Default::default()
            }],
            ..Default::default()
        };
        let pac = pack_seq(&[0, 1, 2, 3, 0, 1, 2, 3]);
        let mut query = vec![0_u8, 1, 2, 3];
        let a = mem_alnreg_t {
            rid: 0,
            rb: 2,
            re: 4,
            qb: 0,
            qe: 2,
            score: 2,
            w: 0,
            ..Default::default()
        };
        let b = mem_alnreg_t {
            rid: 0,
            rb: 9,
            re: 11,
            qb: 2,
            qe: 4,
            score: 2,
            w: 0,
            ..Default::default()
        };
        let mut w = -1;
        assert_eq!(
            mem_patch_reg(
                &opt,
                Some(&bns),
                Some(&pac),
                Some(&mut query),
                &a,
                &b,
                &mut w
            ),
            0
        );
        assert_eq!(w, -1);
    }

    #[test]
    fn mem_patch_reg_scores_colinear_merge() {
        let mut opt = mem_opt_init();
        opt.w = 20;
        let bns = bntseq_t {
            l_pac: 8,
            anns: vec![bntann1_t {
                name: "chr1".into(),
                len: 8,
                ..Default::default()
            }],
            ..Default::default()
        };
        let pac = pack_seq(&[0, 1, 2, 3, 0, 1, 2, 3]);
        let mut query = vec![0_u8, 1, 2, 3];
        let a = mem_alnreg_t {
            rid: 0,
            rb: 0,
            re: 2,
            qb: 0,
            qe: 2,
            score: 2,
            w: 0,
            ..Default::default()
        };
        let b = mem_alnreg_t {
            rid: 0,
            rb: 2,
            re: 4,
            qb: 2,
            qe: 4,
            score: 2,
            w: 0,
            ..Default::default()
        };
        let mut w = -1;
        let score = mem_patch_reg(
            &opt,
            Some(&bns),
            Some(&pac),
            Some(&mut query),
            &a,
            &b,
            &mut w,
        );
        assert_eq!(score, 4);
        assert_eq!(w, 0);
        assert_eq!(query, vec![0, 1, 2, 3]);
    }

    #[test]
    fn mem_dedup_patch_masks_redundant_lower_scoring_hit() {
        let mut opt = mem_opt_init();
        opt.mask_level_redun = 0.90;
        let mut regs = vec![
            mem_alnreg_t {
                rid: 0,
                rb: 0,
                re: 10,
                qb: 0,
                qe: 10,
                score: 20,
                ..Default::default()
            },
            mem_alnreg_t {
                rid: 0,
                rb: 1,
                re: 9,
                qb: 1,
                qe: 9,
                score: 10,
                ..Default::default()
            },
        ];
        let n = mem_dedup_patch(&opt, None, None, None, regs.len() as i32, &mut regs);
        assert_eq!(n, 1);
        assert_eq!(regs.len(), 1);
        assert_eq!(regs[0].score, 20);
        assert_eq!(regs[0].qb, 0);
        assert_eq!(regs[0].qe, 10);
    }

    #[test]
    fn mem_dedup_patch_merges_patchable_hits() {
        let mut opt = mem_opt_init();
        opt.w = 20;
        let bns = bntseq_t {
            l_pac: 8,
            anns: vec![bntann1_t {
                name: "chr1".into(),
                len: 8,
                ..Default::default()
            }],
            ..Default::default()
        };
        let pac = pack_seq(&[0, 1, 2, 3, 0, 1, 2, 3]);
        let mut query = vec![0_u8, 1, 2, 3];
        let mut regs = vec![
            mem_alnreg_t {
                rid: 0,
                rb: 0,
                re: 2,
                qb: 0,
                qe: 2,
                score: 2,
                seedcov: 1,
                sub: 3,
                csub: 4,
                ..Default::default()
            },
            mem_alnreg_t {
                rid: 0,
                rb: 2,
                re: 4,
                qb: 2,
                qe: 4,
                score: 2,
                seedcov: 5,
                sub: 1,
                csub: 2,
                ..Default::default()
            },
        ];
        let n = mem_dedup_patch(
            &opt,
            Some(&bns),
            Some(&pac),
            Some(&mut query),
            regs.len() as i32,
            &mut regs,
        );
        assert_eq!(n, 1);
        assert_eq!(regs.len(), 1);
        assert_eq!(regs[0].rb, 0);
        assert_eq!(regs[0].qb, 0);
        assert_eq!(regs[0].score, 4);
        assert_eq!(regs[0].truesc, 4);
        assert_eq!(regs[0].seedcov, 5);
        assert_eq!(regs[0].sub, 3);
        assert_eq!(regs[0].csub, 4);
        assert_eq!(regs[0].n_comp, 3);
    }

    #[test]
    fn mem_sort_dedup_patch_sorts_and_drops_identical_hits() {
        let opt = mem_opt_init();
        let mut regs = vec![
            mem_alnreg_t {
                rid: 0,
                rb: 20,
                re: 30,
                qb: 5,
                qe: 15,
                score: 7,
                ..Default::default()
            },
            mem_alnreg_t {
                rid: 0,
                rb: 10,
                re: 20,
                qb: 0,
                qe: 10,
                score: 7,
                ..Default::default()
            },
            mem_alnreg_t {
                rid: 0,
                rb: 10,
                re: 18,
                qb: 0,
                qe: 8,
                score: 7,
                ..Default::default()
            },
        ];
        let n = mem_sort_dedup_patch(&opt, None, None, None, regs.len() as i32, &mut regs);
        assert_eq!(n, 2);
        assert_eq!(regs.len(), 2);
        assert_eq!(regs[0].rb, 10);
        assert_eq!(regs[0].qb, 0);
        assert_eq!(regs[1].rb, 20);
        assert_eq!(regs[1].qb, 5);
    }

    #[test]
    fn mem_sort_dedup_patch_keeps_short_secondary_when_longer_overlapping_hit_is_redundant() {
        let opt = mem_opt_init();
        let mut regs = vec![
            mem_alnreg_t {
                rid: 0,
                score: 150,
                qb: 0,
                qe: 150,
                rb: 3_949_610,
                re: 3_949_760,
                ..Default::default()
            },
            mem_alnreg_t {
                rid: 0,
                score: 21,
                qb: 120,
                qe: 141,
                rb: 3_949_748,
                re: 3_949_769,
                ..Default::default()
            },
            mem_alnreg_t {
                rid: 0,
                score: 129,
                qb: 0,
                qe: 147,
                rb: 3_949_610,
                re: 3_949_769,
                ..Default::default()
            },
            mem_alnreg_t {
                rid: 0,
                score: 150,
                qb: -1,
                qe: -1,
                rb: 3_949_610,
                re: 3_949_760,
                ..Default::default()
            },
        ];

        let n = mem_sort_dedup_patch(&opt, None, None, None, regs.len() as i32, &mut regs);
        assert_eq!(n, 2);
        assert_eq!(regs[0].score, 150);
        assert_eq!(regs[1].score, 21);
        assert_eq!(regs[1].qb, 120);
        assert_eq!(regs[1].qe, 141);
    }

    #[test]
    fn test_and_merge_accepts_contained_seed_without_growing_chain() {
        let opt = mem_opt_init();
        let first = mem_seed_t {
            rbeg: 100,
            qbeg: 10,
            len: 20,
            ..Default::default()
        };
        let mut chain = single_seed_chain(first, 3);
        let contained = mem_seed_t {
            rbeg: 105,
            qbeg: 12,
            len: 5,
            ..Default::default()
        };
        assert_eq!(test_and_merge(&opt, 1000, &mut chain, &contained, 3, 0), 1);
        assert_eq!(chain.n, 1);
        assert_eq!(chain.seeds.len(), 1);
    }

    #[test]
    fn test_and_merge_rejects_seed_on_different_strand() {
        let opt = mem_opt_init();
        let first = mem_seed_t {
            rbeg: 100,
            qbeg: 10,
            len: 20,
            ..Default::default()
        };
        let mut chain = single_seed_chain(first, 1);
        let seed = mem_seed_t {
            rbeg: 1200,
            qbeg: 35,
            len: 10,
            ..Default::default()
        };
        assert_eq!(test_and_merge(&opt, 1000, &mut chain, &seed, 1, 0), 0);
        assert_eq!(chain.n, 1);
    }

    #[test]
    fn test_and_merge_grows_chain_and_resizes_seed_storage() {
        let mut opt = mem_opt_init();
        opt.w = 10;
        opt.max_chain_gap = 100;
        let first = mem_seed_t {
            rbeg: 100,
            qbeg: 10,
            len: 20,
            ..Default::default()
        };
        let mut chain = single_seed_chain(first, 2);
        let next = mem_seed_t {
            rbeg: 130,
            qbeg: 40,
            len: 10,
            ..Default::default()
        };
        assert_eq!(test_and_merge(&opt, 1000, &mut chain, &next, 2, 0), 1);
        assert_eq!(chain.n, 2);
        assert_eq!(chain.m, 2);
        assert_eq!(chain.seeds[1], next);
    }

    #[test]
    fn mem_seed_sw_returns_negative_one_for_long_seed_or_window() {
        let opt = mem_opt_init();
        let bns = bntseq_t {
            l_pac: 300,
            anns: vec![bntann1_t {
                name: "chr1".into(),
                len: 300,
                ..Default::default()
            }],
            ..Default::default()
        };
        let pac = pack_seq(&vec![0_u8; 300]);
        let query = vec![0_u8; 250];
        let long_seed = mem_seed_t {
            rbeg: 10,
            qbeg: 10,
            len: 200,
            ..Default::default()
        };
        assert_eq!(mem_seed_sw(&opt, &bns, &pac, 250, &query, &long_seed), -1);
    }

    #[test]
    fn mem_seed_sw_scores_short_exact_match_window() {
        let opt = mem_opt_init();
        let ref_seq = [0_u8, 1, 2, 3, 0, 1, 2, 3];
        let bns = bntseq_t {
            l_pac: 8,
            anns: vec![bntann1_t {
                name: "chr1".into(),
                len: 8,
                ..Default::default()
            }],
            ..Default::default()
        };
        let pac = pack_seq(&ref_seq);
        let query = vec![3_u8, 3, 0, 1, 2, 3, 0];
        let seed = mem_seed_t {
            rbeg: 2,
            qbeg: 2,
            len: 2,
            ..Default::default()
        };
        assert_eq!(
            mem_seed_sw(
                &opt,
                &bns,
                &pac,
                i32::try_from(query.len()).expect("qlen"),
                &query,
                &seed
            ),
            5
        );
    }

    #[test]
    fn mem_chain_weight_uses_minimum_non_overlapping_span_on_query_and_ref() {
        let chain = mem_chain_t {
            n: 3,
            seeds: vec![
                mem_seed_t {
                    qbeg: 0,
                    rbeg: 10,
                    len: 5,
                    ..Default::default()
                },
                mem_seed_t {
                    qbeg: 3,
                    rbeg: 20,
                    len: 5,
                    ..Default::default()
                },
                mem_seed_t {
                    qbeg: 10,
                    rbeg: 23,
                    len: 4,
                    ..Default::default()
                },
            ],
            ..Default::default()
        };
        assert_eq!(mem_chain_weight(&chain), 12);
    }

    #[test]
    fn mem_print_chain_formats_seed_coordinates() {
        let bns = bntseq_t {
            l_pac: 100,
            anns: vec![bntann1_t {
                name: "chr1".into(),
                offset: 10,
                len: 90,
                ..Default::default()
            }],
            ..Default::default()
        };
        let chain = mem_chain_t {
            rid: 0,
            n: 1,
            seeds: vec![mem_seed_t {
                score: 7,
                len: 5,
                qbeg: 3,
                rbeg: 20,
                ..Default::default()
            }],
            ..Default::default()
        };
        assert_eq!(
            mem_print_chain(&bns, &[chain]),
            "* Found CHAIN(0): n=1; weight=5\t7;5;3,20(chr1:+11)\n"
        );
    }

    #[test]
    fn mem_flt_chained_seeds_keeps_high_scoring_and_fallback_seeds() {
        let mut opt = mem_opt_init();
        opt.min_chain_weight = 1;
        let bns = bntseq_t {
            l_pac: 32,
            anns: vec![bntann1_t {
                name: "chr1".into(),
                len: 32,
                ..Default::default()
            }],
            ..Default::default()
        };
        let pac = pack_seq(&[
            0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0,
            1, 2, 3,
        ]);
        let seqs = vec![bseq1_t {
            l_seq: 32,
            seq: Some("TTACGTACGTACGTACGTACGTTTTTTTTTTT".into()),
            ..Default::default()
        }];
        let mut chains = vec![mem_chain_t {
            seqid: 0,
            n: 2,
            seeds: vec![
                mem_seed_t {
                    rbeg: 2,
                    qbeg: 2,
                    len: 20,
                    ..Default::default()
                },
                mem_seed_t {
                    rbeg: 0,
                    qbeg: 0,
                    len: 201,
                    score: 0,
                    ..Default::default()
                },
            ],
            ..Default::default()
        }];
        mem_flt_chained_seeds(&opt, &bns, &pac, &seqs, 1, &mut chains);
        assert_eq!(chains[0].n, 2);
        assert!(chains[0].seeds[0].score > 0);
        assert_eq!(chains[0].seeds[1].score, 201);
    }

    #[test]
    fn mem_chain_flt_drops_low_weight_chains() {
        let mut opt = mem_opt_init();
        opt.min_chain_weight = 10;
        let mut chains = vec![
            mem_chain_t {
                seqid: 0,
                n: 1,
                seeds: vec![mem_seed_t {
                    qbeg: 0,
                    rbeg: 0,
                    len: 5,
                    ..Default::default()
                }],
                ..Default::default()
            },
            mem_chain_t {
                seqid: 0,
                n: 1,
                seeds: vec![mem_seed_t {
                    qbeg: 10,
                    rbeg: 10,
                    len: 12,
                    ..Default::default()
                }],
                ..Default::default()
            },
        ];
        let n = mem_chain_flt(&opt, 2, &mut chains, 0);
        assert_eq!(n, 1);
        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0].w, 12);
    }

    #[test]
    fn mem_chain_flt_marks_overlapped_chain_as_shadowed() {
        let mut opt = mem_opt_init();
        opt.min_chain_weight = 1;
        opt.mask_level = 0.5;
        opt.drop_ratio = 0.5;
        opt.min_seed_len = 2;
        opt.max_chain_extend = 10;
        let mut chains = vec![
            mem_chain_t {
                seqid: 0,
                n: 2,
                seeds: vec![
                    mem_seed_t {
                        qbeg: 0,
                        rbeg: 0,
                        len: 10,
                        ..Default::default()
                    },
                    mem_seed_t {
                        qbeg: 10,
                        rbeg: 10,
                        len: 10,
                        ..Default::default()
                    },
                ],
                ..Default::default()
            },
            mem_chain_t {
                seqid: 0,
                n: 1,
                seeds: vec![mem_seed_t {
                    qbeg: 5,
                    rbeg: 5,
                    len: 10,
                    ..Default::default()
                }],
                ..Default::default()
            },
        ];
        let n = mem_chain_flt(&opt, 2, &mut chains, 0);
        assert_eq!(n, 2);
        assert_eq!(chains[0].kept, 3);
        assert_eq!(chains[1].kept, 1);
        assert_eq!(chains[0].first, 1);
    }

    #[test]
    fn mem_chain_flt_matches_cpp_equal_weight_order_for_read68_shape() {
        let mut chains = vec![
            mem_chain_t {
                pos: 2635834,
                w: 20,
                ..Default::default()
            },
            mem_chain_t {
                pos: 3529141,
                w: 20,
                ..Default::default()
            },
            mem_chain_t {
                pos: 3671498,
                w: 20,
                ..Default::default()
            },
            mem_chain_t {
                pos: 4425228,
                w: 20,
                ..Default::default()
            },
            mem_chain_t {
                pos: 5736348,
                w: 131,
                ..Default::default()
            },
        ];
        ks_introsort_mem_flt(&mut chains);
        let positions: Vec<i64> = chains.iter().map(|c| c.pos).collect();
        assert_eq!(positions, vec![5736348, 4425228, 2635834, 3529141, 3671498]);
    }

    #[test]
    fn mem_collect_smem_collects_and_sorts_real_index_matches() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_collect_smem_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("ref");
        let fasta = b">chr1\nACGTACGTACGT\n";
        let l_pac = bns_fasta2bntseq(
            Cursor::new(fasta.as_slice()),
            prefix.to_str().expect("utf8"),
            1,
        );
        assert_eq!(l_pac, 12);

        let mut fmi = FMI_search::ctor(prefix.to_str().expect("utf8"));
        assert_eq!(fmi.build_index(), 0);
        fmi.load_index();

        let seqs = vec![
            bseq1_t {
                l_seq: 6,
                seq: Some("ACGTAC".into()),
                ..Default::default()
            },
            bseq1_t {
                l_seq: 4,
                seq: Some("GTAC".into()),
                ..Default::default()
            },
        ];
        let mut opt = mem_opt_init();
        opt.min_seed_len = 2;
        opt.split_factor = 1.5;
        opt.split_width = 10;
        opt.max_mem_intv = 20;
        let mut match_array = Vec::new();
        let mut min_intv_ar = Vec::new();
        let mut query_pos_ar = Vec::new();
        let mut enc_qdb = Vec::new();
        let mut rid = Vec::new();
        let mut mmc = mem_cache::default();
        let mut tot_smem = 0_i64;

        mem_collect_smem(
            &fmi,
            &opt,
            &seqs,
            i32::try_from(seqs.len()).expect("nseq"),
            &mut match_array,
            &mut min_intv_ar,
            &mut query_pos_ar,
            &mut enc_qdb,
            &mut rid,
            &mut mmc,
            &mut tot_smem,
            0,
        );

        assert!(tot_smem > 0);
        assert_eq!(
            usize::try_from(tot_smem).expect("tot_smem"),
            match_array.len()
        );
        for window in match_array.windows(2) {
            let a = window[0];
            let b = window[1];
            if a.rid == b.rid {
                assert!((a.m, a.n) <= (b.m, b.n));
            }
        }

        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn mem_chain_seeds_builds_non_empty_chains_from_real_smems() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_chain_seeds_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("ref");
        let fasta = b">chr1\nACGTACGTACGT\n";
        bns_fasta2bntseq(
            Cursor::new(fasta.as_slice()),
            prefix.to_str().expect("utf8"),
            1,
        );
        let bns = crate::generated::bntseq_cpp::bns_restore(prefix.to_str().expect("utf8"));
        let pac = fs::read(prefix.with_extension("pac")).expect("read pac");

        let mut fmi = FMI_search::ctor(prefix.to_str().expect("utf8"));
        assert_eq!(fmi.build_index(), 0);
        fmi.load_index();

        let seqs = vec![
            bseq1_t {
                l_seq: 6,
                seq: Some("ACGTAC".into()),
                ..Default::default()
            },
            bseq1_t {
                l_seq: 4,
                seq: Some("GTAC".into()),
                ..Default::default()
            },
        ];
        let mut opt = mem_opt_init();
        opt.min_seed_len = 2;
        opt.split_factor = 1.5;
        opt.split_width = 10;
        opt.max_mem_intv = 20;
        let mut match_array = Vec::new();
        let mut min_intv_ar = Vec::new();
        let mut query_pos_ar = Vec::new();
        let mut enc_qdb = Vec::new();
        let mut rid = Vec::new();
        let mut mmc = mem_cache::default();
        let mut tot_smem = 0_i64;
        mem_collect_smem(
            &fmi,
            &opt,
            &seqs,
            i32::try_from(seqs.len()).expect("nseq"),
            &mut match_array,
            &mut min_intv_ar,
            &mut query_pos_ar,
            &mut enc_qdb,
            &mut rid,
            &mut mmc,
            &mut tot_smem,
            0,
        );

        let mut chain_ar = vec![mem_chain_v::default(); seqs.len()];
        let mut seed_buf = vec![mem_seed_t::default(); 128];
        mem_chain_seeds(
            &fmi,
            &opt,
            &bns,
            &seqs,
            i32::try_from(seqs.len()).expect("nseq"),
            0,
            &mut chain_ar,
            &mut seed_buf,
            128,
            &match_array,
            tot_smem,
        );

        assert!(chain_ar.iter().any(|c| c.n > 0));
        for chains in &chain_ar {
            for c in &chains.a {
                assert!(c.n > 0);
                assert_eq!(
                    usize::try_from(c.n).expect("c.n"),
                    c.seeds
                        .iter()
                        .take(usize::try_from(c.n).expect("n"))
                        .count()
                );
                assert!(c.rid >= 0);
            }
        }

        drop(pac);
        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn mem_kernel1_core_runs_full_kernel1_pipeline_on_real_index_data() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_kernel1_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("ref");
        let fasta = b">chr1\nACGTACGTACGT\n";
        bns_fasta2bntseq(
            Cursor::new(fasta.as_slice()),
            prefix.to_str().expect("utf8"),
            1,
        );

        let mut fmi = FMI_search::ctor(prefix.to_str().expect("utf8"));
        assert_eq!(fmi.build_index(), 0);
        fmi.load_index();

        let seqs = vec![
            bseq1_t {
                l_seq: 6,
                seq: Some("ACGTAC".into()),
                ..Default::default()
            },
            bseq1_t {
                l_seq: 4,
                seq: Some("GTAC".into()),
                ..Default::default()
            },
        ];
        let mut opt = mem_opt_init();
        opt.min_seed_len = 2;
        opt.split_factor = 1.5;
        opt.split_width = 10;
        opt.max_mem_intv = 20;
        opt.min_chain_weight = 1;

        let mut chain_ar = vec![mem_chain_v::default(); seqs.len()];
        let mut seed_buf = vec![mem_seed_t::default(); 256];
        let seed_buf_len = i64::try_from(seed_buf.len()).expect("seed buf len");
        let mut mmc = mem_cache::default();

        let ret = mem_kernel1_core(
            &fmi,
            &opt,
            &seqs,
            i32::try_from(seqs.len()).expect("nseq"),
            &mut chain_ar,
            &mut seed_buf,
            seed_buf_len,
            &mut mmc,
            0,
        );

        assert_eq!(ret, 1);
        assert!(!mmc.matchArray.is_empty());
        assert!(!mmc.enc_qdb.is_empty());
        assert!(chain_ar.iter().any(|c| c.n > 0));
        for chains in &chain_ar {
            assert_eq!(chains.n, chains.a.len());
            for c in &chains.a {
                assert!(c.n >= 0);
                assert_eq!(c.rid >= 0, c.n > 0);
            }
        }

        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn bns_get_seq_v2_returns_forward_and_reverse_windows_from_ref_string_layout() {
        let forward = vec![0_u8, 1, 2, 3];
        let mut ref_string = forward.clone();
        ref_string.extend_from_slice(&revcomp_nt4(&forward));
        let mut len = -1_i64;

        let seq = bns_get_seq_v2(4, &[], 1, 3, &mut len, &ref_string, None).expect("forward seq");
        assert_eq!(seq, &[1, 2]);
        assert_eq!(len, 2);

        let seq = bns_get_seq_v2(4, &[], 4, 7, &mut len, &ref_string, None).expect("reverse seq");
        assert_eq!(seq, &[0, 1, 2]);
        assert_eq!(len, 3);

        let seq = bns_get_seq_v2(4, &[], 3, 5, &mut len, &ref_string, None);
        assert!(seq.is_none());
        assert_eq!(len, 0);
    }

    #[test]
    fn bns_fetch_seq_v2_clips_to_reference_interval_and_returns_expected_rid() {
        let bns = bntseq_t {
            l_pac: 8,
            n_seqs: 2,
            anns: vec![
                bntann1_t {
                    offset: 0,
                    len: 4,
                    name: "chr1".into(),
                    anno: String::new(),
                    ..Default::default()
                },
                bntann1_t {
                    offset: 4,
                    len: 4,
                    name: "chr2".into(),
                    anno: String::new(),
                    ..Default::default()
                },
            ],
            ..Default::default()
        };
        let forward = vec![0_u8, 1, 2, 3, 3, 2, 1, 0];
        let mut ref_string = forward.clone();
        ref_string.extend_from_slice(&revcomp_nt4(&forward));

        let mut beg = 2_i64;
        let mut end = 7_i64;
        let mut rid = -1;
        let seq = bns_fetch_seq_v2(
            &bns,
            &[],
            &mut beg,
            5,
            &mut end,
            &mut rid,
            &ref_string,
            None,
        );
        assert_eq!(rid, 1);
        assert_eq!((beg, end), (4, 7));
        assert_eq!(seq, &[3, 2, 1]);

        let mut beg = 9_i64;
        let mut end = 15_i64;
        let mut rid = -1;
        let seq = bns_fetch_seq_v2(
            &bns,
            &[],
            &mut beg,
            10,
            &mut end,
            &mut rid,
            &ref_string,
            None,
        );
        assert_eq!(rid, 1);
        assert_eq!((beg, end), (9, 12));
        assert_eq!(seq, &ref_string[9..12]);
    }

    #[test]
    fn sort_pairs_len_ext_partitions_into_8_16_and_scalar_buckets() {
        let mut pairs = vec![
            SeqPair {
                id: 10,
                len1: 5,
                len2: 7,
                h0: 3,
                ..Default::default()
            },
            SeqPair {
                id: 20,
                len1: 200,
                len2: 150,
                h0: 10,
                ..Default::default()
            },
            SeqPair {
                id: 30,
                len1: 40_000,
                len2: 2,
                h0: 1,
                ..Default::default()
            },
            SeqPair {
                id: 11,
                len1: 5,
                len2: 6,
                h0: 2,
                ..Default::default()
            },
        ];
        let mut temp = vec![SeqPair::default(); pairs.len()];
        let mut hist = vec![0_i32; MAX_SEQ_LEN8 + MAX_SEQ_LEN16 + 1];
        let mut n128 = 0;
        let mut n16 = 0;
        let mut n1 = 0;
        let count = pairs.len() as i32;

        sortPairsLenExt(
            &mut pairs, count, &mut temp, &mut hist, &mut n128, &mut n16, &mut n1, 1,
        );

        assert_eq!((n128, n16, n1), (2, 1, 1));
        assert_eq!(
            pairs.iter().map(|p| p.id).collect::<Vec<_>>(),
            vec![11, 10, 20, 30]
        );
    }

    #[test]
    fn sort_pairs_len_orders_by_len1_stably() {
        let mut pairs = vec![
            SeqPair {
                id: 1,
                len1: 9,
                ..Default::default()
            },
            SeqPair {
                id: 2,
                len1: 3,
                ..Default::default()
            },
            SeqPair {
                id: 3,
                len1: 9,
                ..Default::default()
            },
            SeqPair {
                id: 4,
                len1: 1,
                ..Default::default()
            },
        ];
        let mut temp = vec![SeqPair::default(); pairs.len()];
        let mut hist = vec![0_i32; MAX_SEQ_LEN16 + 1];
        let count = pairs.len() as i32;

        sortPairsLen(&mut pairs, count, &mut temp, &mut hist);

        assert_eq!(
            pairs.iter().map(|p| (p.len1, p.id)).collect::<Vec<_>>(),
            vec![(1, 4), (3, 2), (9, 1), (9, 3)]
        );
    }

    #[test]
    fn mem_chain2aln_across_reads_v2_extends_internal_seed_to_full_match() {
        let forward = vec![0_u8, 1, 2, 3, 0, 1, 2, 3];
        let mut ref_string = forward.clone();
        ref_string.extend_from_slice(&revcomp_nt4(&forward));
        let bns = bntseq_t {
            l_pac: 8,
            n_seqs: 1,
            anns: vec![bntann1_t {
                offset: 0,
                len: 8,
                name: "chr1".into(),
                anno: String::new(),
                ..Default::default()
            }],
            ..Default::default()
        };
        let seqs = vec![bseq1_t {
            l_seq: 5,
            seq: Some("CGTAC".into()),
            ..Default::default()
        }];
        let chain = mem_chain_t {
            seqid: 0,
            rid: 0,
            n: 1,
            m: 1,
            seeds: vec![mem_seed_t {
                qbeg: 1,
                rbeg: 2,
                len: 3,
                score: 3,
                ..Default::default()
            }],
            ..Default::default()
        };
        let mut chain_ar = vec![mem_chain_v {
            n: 1,
            m: 1,
            a: vec![chain],
            ..Default::default()
        }];
        let mut av_v = vec![mem_alnreg_v::default()];
        let mut mmc = mem_cache::default();
        let opt = mem_opt_init();

        mem_chain2aln_across_reads_V2(
            &opt,
            &bns,
            &[],
            &seqs,
            1,
            &mut chain_ar,
            &mut av_v,
            &mut mmc,
            &ref_string,
            0,
        );

        assert_eq!(av_v[0].n, 1);
        let aln = &av_v[0].a[0];
        assert_eq!((aln.qb, aln.qe), (0, 5));
        assert_eq!((aln.rb, aln.re), (1, 6));
        assert!(aln.score > 0);
        assert!(aln.truesc > 0);
        assert_eq!(aln.seedcov, 3);
    }

    #[test]
    fn mem_chain2aln_across_reads_v2_keeps_distinct_chain_alignment_when_later_chain_has_multiple_seeds(
    ) {
        let query = "TGATGGAAGTGATCCGTGATATGGGCGTAGAAAAAACCGTTAGTCTGGCTGGAAGGTGTTTCATCGGTACTGTTCGGTGTAGCCGCATTCTCGACCTTTTATGTGATTGTGGTGTGGATGCCCAAATATGCGATGGCTTTTGCTGGTATG";
        let query_nt4 = seq_to_nt4(query);

        let mut forward = vec![0_u8; 320];
        forward[100..217].copy_from_slice(&query_nt4[33..150]);
        forward[200..241].copy_from_slice(&query_nt4[0..41]);

        let bns = bntseq_t {
            l_pac: i64::try_from(forward.len()).expect("len"),
            n_seqs: 1,
            anns: vec![bntann1_t {
                offset: 0,
                len: i32::try_from(forward.len()).expect("len"),
                name: "chr1".into(),
                anno: String::new(),
                ..Default::default()
            }],
            ..Default::default()
        };
        let pac = pack_seq(&forward);
        let mut ref_string = forward.clone();
        ref_string.extend_from_slice(&revcomp_nt4(&forward));

        let seqs = vec![bseq1_t {
            l_seq: i32::try_from(query.len()).expect("query len"),
            seq: Some(query.to_string().into()),
            ..Default::default()
        }];
        let chain0 = mem_chain_t {
            seqid: 0,
            rid: 0,
            n: 1,
            m: 1,
            seeds: vec![mem_seed_t {
                qbeg: 33,
                rbeg: 100,
                len: 117,
                score: 117,
                ..Default::default()
            }],
            ..Default::default()
        };
        let chain1 = mem_chain_t {
            seqid: 0,
            rid: 0,
            n: 2,
            m: 2,
            seeds: vec![
                mem_seed_t {
                    qbeg: 0,
                    rbeg: 200,
                    len: 20,
                    score: 20,
                    ..Default::default()
                },
                mem_seed_t {
                    qbeg: 0,
                    rbeg: 200,
                    len: 41,
                    score: 41,
                    ..Default::default()
                },
            ],
            ..Default::default()
        };
        let mut chain_ar = vec![mem_chain_v {
            n: 2,
            m: 2,
            a: vec![chain0, chain1],
            ..Default::default()
        }];
        let mut av_v = vec![mem_alnreg_v::default()];
        let mut mmc = mem_cache::default();
        let opt = mem_opt_init();

        mem_chain2aln_across_reads_V2(
            &opt,
            &bns,
            &pac,
            &seqs,
            1,
            &mut chain_ar,
            &mut av_v,
            &mut mmc,
            &ref_string,
            0,
        );

        let chain_ids: Vec<usize> = av_v[0].a.iter().take(av_v[0].n).map(|aln| aln.c).collect();
        assert!(
            chain_ids.contains(&0),
            "missing chain 0 alignment: {:?}",
            av_v[0].a
        );
        assert!(
            chain_ids.contains(&1),
            "missing chain 1 alignment: {:?}",
            av_v[0].a
        );
    }

    #[test]
    fn mem_kernel2_core_runs_alignment_pipeline_and_marks_alt_hits() {
        let forward = vec![0_u8, 1, 2, 3, 0, 1, 2, 3];
        let mut ref_string = forward.clone();
        ref_string.extend_from_slice(&revcomp_nt4(&forward));

        let mut fmi = FMI_search::ctor("kernel2-test");
        fmi.base.idx.bns = Some(bntseq_t {
            l_pac: 8,
            n_seqs: 1,
            anns: vec![bntann1_t {
                offset: 0,
                len: 8,
                is_alt: 1,
                name: "chr1".into(),
                anno: String::new(),
                ..Default::default()
            }],
            ..Default::default()
        });
        fmi.base.idx.pac = vec![];

        let seqs = vec![bseq1_t {
            l_seq: 5,
            seq: Some("CGTAC".into()),
            ..Default::default()
        }];
        let chain = mem_chain_t {
            seqid: 0,
            rid: 0,
            n: 1,
            m: 1,
            seeds: vec![mem_seed_t {
                qbeg: 1,
                rbeg: 2,
                len: 3,
                score: 3,
                ..Default::default()
            }],
            ..Default::default()
        };
        let mut chain_ar = vec![mem_chain_v {
            n: 1,
            m: 1,
            a: vec![chain],
            ..Default::default()
        }];
        let mut regs = vec![mem_alnreg_v::default()];
        let mut mmc = mem_cache::default();
        let opt = mem_opt_init();

        let ret = mem_kernel2_core(
            &fmi,
            &opt,
            &seqs,
            &mut regs,
            1,
            &mut chain_ar,
            &mut mmc,
            &ref_string,
            0,
        );

        assert_eq!(ret, 1);
        assert_eq!(chain_ar[0].n, 0);
        assert_eq!(regs[0].n, 1);
        let aln = &regs[0].a[0];
        assert_eq!((aln.qb, aln.qe), (0, 5));
        assert_eq!((aln.rb, aln.re), (1, 6));
        assert_eq!(aln.is_alt, 1);
    }

    #[test]
    fn sort_classify_partitions_byte_and_word_sw_requests() {
        let mut mmc = mem_cache::default();
        mmc.seqPairArrayLeft128.push(vec![
            SeqPair {
                id: 1,
                h0: KSW_XBYTE,
                ..Default::default()
            },
            SeqPair {
                id: 2,
                h0: 0,
                ..Default::default()
            },
            SeqPair {
                id: 3,
                h0: KSW_XBYTE | 7,
                ..Default::default()
            },
        ]);
        mmc.seqPairArrayRight128.push(vec![SeqPair::default(); 3]);

        let pos8 = sort_classify(&mut mmc, 3, 0);

        assert_eq!(pos8, 2);
        assert_eq!(
            mmc.seqPairArrayLeft128[0]
                .iter()
                .map(|p| p.id)
                .collect::<Vec<_>>(),
            vec![1, 3, 2]
        );
    }

    #[test]
    fn worker_bwt_populates_chains_for_batch_window() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_worker_bwt_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("ref");
        let fasta = b">chr1\nACGTACGTACGT\n";
        bns_fasta2bntseq(
            Cursor::new(fasta.as_slice()),
            prefix.to_str().expect("utf8"),
            1,
        );
        let mut fmi = FMI_search::ctor(prefix.to_str().expect("utf8"));
        assert_eq!(fmi.build_index(), 0);
        fmi.load_index();

        let seqs = vec![
            bseq1_t {
                l_seq: 6,
                seq: Some("ACGTAC".into()),
                ..Default::default()
            },
            bseq1_t {
                l_seq: 4,
                seq: Some("GTAC".into()),
                ..Default::default()
            },
        ];
        let mut opt = mem_opt_init();
        opt.min_seed_len = 2;
        opt.split_factor = 1.5;
        opt.split_width = 10;
        opt.max_mem_intv = 20;
        let mut worker = worker_t {
            opt: Some(opt),
            seqs,
            chain_ar: vec![mem_chain_v::default(); 2],
            seedBuf: vec![mem_seed_t::default(); 2 * AVG_SEEDS_PER_READ],
            seedBufSize: i64::try_from(2 * AVG_SEEDS_PER_READ).expect("seedBufSize"),
            nreads: 2,
            fmi: Some(fmi),
            ..Default::default()
        };

        worker_bwt(&mut worker, 0, 2, 0);

        assert!(worker.chain_ar.iter().any(|c| c.n > 0));
        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn worker_aln_populates_regs_for_batch_window() {
        let forward = vec![0_u8, 1, 2, 3, 0, 1, 2, 3];
        let mut ref_string = forward.clone();
        ref_string.extend_from_slice(&revcomp_nt4(&forward));

        let mut fmi = FMI_search::ctor("worker-aln-test");
        fmi.base.idx.bns = Some(bntseq_t {
            l_pac: 8,
            n_seqs: 1,
            anns: vec![bntann1_t {
                offset: 0,
                len: 8,
                is_alt: 0,
                name: "chr1".into(),
                anno: String::new(),
                ..Default::default()
            }],
            ..Default::default()
        });
        fmi.base.idx.pac = vec![];

        let seqs = vec![bseq1_t {
            l_seq: 5,
            seq: Some("CGTAC".into()),
            ..Default::default()
        }];
        let chain = mem_chain_t {
            seqid: 0,
            rid: 0,
            n: 1,
            m: 1,
            seeds: vec![mem_seed_t {
                qbeg: 1,
                rbeg: 2,
                len: 3,
                score: 3,
                ..Default::default()
            }],
            ..Default::default()
        };
        let opt = mem_opt_init();
        let mut worker = worker_t {
            opt: Some(opt),
            seqs,
            regs: vec![mem_alnreg_v::default()],
            chain_ar: vec![mem_chain_v {
                n: 1,
                m: 1,
                a: vec![chain],
                ..Default::default()
            }],
            ref_string,
            fmi: Some(fmi),
            ..Default::default()
        };

        worker_aln(&mut worker, 0, 1, 0);

        assert_eq!(worker.regs[0].n, 1);
        assert_eq!((worker.regs[0].a[0].qb, worker.regs[0].a[0].qe), (0, 5));
    }

    #[test]
    fn mem_mark_primary_se_core_marks_overlapping_lower_hit_secondary() {
        let mut opt = mem_opt_init();
        opt.mask_level = 0.5;
        let mut regs = vec![
            mem_alnreg_t {
                qb: 0,
                qe: 10,
                score: 20,
                ..Default::default()
            },
            mem_alnreg_t {
                qb: 2,
                qe: 9,
                score: 18,
                ..Default::default()
            },
            mem_alnreg_t {
                qb: 20,
                qe: 30,
                score: 15,
                ..Default::default()
            },
        ];
        let mut z = Vec::new();
        mem_mark_primary_se_core(&opt, regs.len() as i32, &mut regs, &mut z);
        assert_eq!(regs[1].secondary, 0);
        assert_eq!(regs[0].sub, 18);
        assert_eq!(z, vec![0, 2]);
    }

    #[test]
    fn mem_mark_primary_se_sorts_and_counts_primary_hits() {
        let opt = mem_opt_init();
        let mut regs = vec![
            mem_alnreg_t {
                score: 15,
                qb: 0,
                qe: 10,
                is_alt: 0,
                ..Default::default()
            },
            mem_alnreg_t {
                score: 12,
                qb: 1,
                qe: 9,
                is_alt: 1,
                ..Default::default()
            },
            mem_alnreg_t {
                score: 14,
                qb: 20,
                qe: 30,
                is_alt: 0,
                ..Default::default()
            },
        ];
        let n_pri = mem_mark_primary_se(&opt, regs.len() as i32, &mut regs, 42);
        assert_eq!(n_pri, 2);
        assert_eq!(regs[0].is_alt, 0);
        assert_eq!(regs[1].is_alt, 0);
        assert_eq!(regs[2].is_alt, 1);
        assert!(regs.iter().all(|r| r.hash != 0));
    }

    #[test]
    fn mem_mark_primary_se_uses_cpp_id_plus_index_hash_tiebreak() {
        let opt = mem_opt_init();
        let regs = vec![
            mem_alnreg_t {
                score: 150,
                qb: 0,
                qe: 150,
                rb: 2254222,
                re: 2254372,
                rid: 0,
                is_alt: 0,
                seedcov: 42,
                ..Default::default()
            },
            mem_alnreg_t {
                score: 29,
                qb: 0,
                qe: 29,
                rb: 2254585,
                re: 2254614,
                rid: 0,
                is_alt: 0,
                seedcov: 12,
                ..Default::default()
            },
            mem_alnreg_t {
                score: 29,
                qb: 0,
                qe: 29,
                rb: 2255037,
                re: 2255066,
                rid: 0,
                is_alt: 0,
                seedcov: 12,
                ..Default::default()
            },
        ];
        let mut regs_a = regs.clone();
        let mut regs_b = vec![regs[0], regs[2], regs[1]];
        let n_pri_a = mem_mark_primary_se(&opt, regs_a.len() as i32, &mut regs_a, 7);
        let n_pri_b = mem_mark_primary_se(&opt, regs_b.len() as i32, &mut regs_b, 7);
        assert_eq!(n_pri_a, n_pri_b);
        assert_eq!(regs_a[0].rb, 2254222);
        assert_eq!(regs_b[0].rb, 2254222);
        let tail_a: Vec<i64> = regs_a[1..].iter().map(|r| r.rb).collect();
        let tail_b: Vec<i64> = regs_b[1..].iter().map(|r| r.rb).collect();
        assert_ne!(tail_a, tail_b);
    }

    #[test]
    fn mem_approx_mapq_se_clamps_and_scales_by_repeat_fraction() {
        let mut opt = mem_opt_init();
        opt.mapQ_coef_len = 50.0;
        opt.mapQ_coef_fac = (50.0_f32.ln()) as i32;
        let reg = mem_alnreg_t {
            score: 100,
            sub: 80,
            csub: 0,
            qb: 0,
            qe: 100,
            rb: 0,
            re: 100,
            seedcov: 50,
            sub_n: 2,
            frac_rep: 0.25,
            ..Default::default()
        };
        let mapq = mem_approx_mapq_se(&opt, &reg);
        assert!((0..=60).contains(&mapq));
        assert!(mapq > 0);
    }

    #[test]
    fn mem_reorder_primary5_moves_leftmost_primary_to_front() {
        let mut regs = mem_alnreg_v {
            n: 3,
            m: 3,
            a: vec![
                mem_alnreg_t {
                    qb: 30,
                    score: 50,
                    secondary: -1,
                    is_alt: 0,
                    ..Default::default()
                },
                mem_alnreg_t {
                    qb: 10,
                    score: 45,
                    secondary: -1,
                    is_alt: 0,
                    ..Default::default()
                },
                mem_alnreg_t {
                    qb: 20,
                    score: 40,
                    secondary: 0,
                    secondary_all: 0,
                    is_alt: 0,
                    ..Default::default()
                },
            ],
        };
        mem_reorder_primary5(30, &mut regs);
        assert_eq!(regs.a[0].qb, 10);
        assert_eq!(regs.a[2].secondary, 1);
        assert_eq!(regs.a[2].secondary_all, 1);
    }

    #[test]
    fn get_rlen_counts_match_and_delete_ops_only() {
        let cigar = vec![
            (5_u32 << 4) | 0,
            (3_u32 << 4) | 1,
            (2_u32 << 4) | 2,
            (4_u32 << 4) | 4,
        ];
        assert_eq!(get_rlen(cigar.len() as i32, &cigar), 7);
    }

    #[test]
    fn add_cigar_formats_ops_and_rewrites_softclip_for_supplementary_primary_assembly() {
        let opt = mem_opt_init();
        let aln = mem_aln_t {
            n_cigar: 3,
            cigar: vec![(2_u32 << 4) | 4, (5_u32 << 4) | 0, (1_u32 << 4) | 3],
            ..Default::default()
        };
        let mut ks = crate::generated::kstring_h::kstring_t::default();
        add_cigar(&opt, &aln, &mut ks, 1);
        assert_eq!(ks.as_str(), "2H5M1H");

        let mut soft_opt = mem_opt_init();
        soft_opt.flag |= 0x200;
        let mut ks = crate::generated::kstring_h::kstring_t::default();
        add_cigar(&soft_opt, &aln, &mut ks, 1);
        assert_eq!(ks.as_str(), "2H5M1S");
    }

    #[test]
    fn infer_bw_returns_zero_when_equal_lengths_need_no_gap_band() {
        assert_eq!(infer_bw(100, 100, 98, 1, 6, 1), 0);
    }

    #[test]
    fn infer_bw_is_at_least_length_difference() {
        assert!(infer_bw(100, 90, 70, 1, 6, 1) >= 10);
    }

    #[test]
    fn mem_reg2aln_returns_unmapped_record_for_missing_region() {
        let opt = mem_opt_init();
        let bns = bntseq_t::default();
        let aln = mem_reg2aln(&opt, &bns, &[], 4, "ACGT", None);
        assert_eq!(aln.rid, -1);
        assert_eq!(aln.pos, -1);
        assert_ne!(aln.flag & 0x4, 0);
    }

    #[test]
    fn mem_reg2aln_generates_clipped_cigar_md_and_reference_position() {
        let opt = mem_opt_init();
        let forward = vec![0_u8, 1, 2, 3];
        let pac = pack_seq(&forward);
        let bns = bntseq_t {
            l_pac: 4,
            n_seqs: 1,
            anns: vec![crate::generated::bntseq_h::bntann1_t {
                offset: 0,
                len: 4,
                name: "chr1".into(),
                ..Default::default()
            }],
            ..Default::default()
        };
        let reg = mem_alnreg_t {
            rb: 0,
            re: 4,
            qb: 1,
            qe: 5,
            rid: 0,
            score: 4,
            truesc: 4,
            sub: 2,
            csub: 1,
            secondary: -1,
            w: 100,
            ..Default::default()
        };

        let aln = mem_reg2aln(&opt, &bns, &pac, 6, "TACGTN", Some(&reg));
        assert_eq!(aln.rid, 0);
        assert_eq!(aln.pos, 0);
        assert_eq!(aln.is_rev, 0);
        assert_eq!(aln.n_cigar, 3);
        assert_eq!(
            aln.cigar,
            vec![(1_u32 << 4) | 3, (4_u32 << 4), (1_u32 << 4) | 3]
        );
        assert_eq!(std::str::from_utf8(&aln.md).expect("md utf8"), "4");
        assert_eq!(aln.NM, 0);
        assert_eq!(aln.score, 4);
        assert_eq!(aln.sub, 2);
    }

    #[test]
    fn mem_aln2sam_formats_primary_paired_record_with_tags() {
        let opt = mem_opt_init();
        let bns = bntseq_t {
            l_pac: 8,
            n_seqs: 1,
            anns: vec![crate::generated::bntseq_h::bntann1_t {
                offset: 0,
                len: 8,
                name: "chr1".into(),
                anno: "anno\tfield".into(),
                ..Default::default()
            }],
            ..Default::default()
        };
        let s = bseq1_t {
            l_seq: 4,
            name: Some("read1".into()),
            seq: Some("ACGT".into()),
            qual: Some("!!!!".into()),
            comment: Some("CO:Z:note".into()),
            ..Default::default()
        };
        let list = vec![
            mem_aln_t {
                rid: 0,
                pos: 2,
                mapq: 60,
                n_cigar: 1,
                cigar: vec![4_u32 << 4],
                md: b"4".to_vec(),
                NM: 0,
                score: 4,
                sub: 1,
                alt_sc: 8,
                ..Default::default()
            },
            mem_aln_t {
                rid: 0,
                pos: 5,
                mapq: 20,
                n_cigar: 1,
                cigar: vec![4_u32 << 4],
                md: b"4".to_vec(),
                NM: 0,
                score: 4,
                flag: 0x800,
                ..Default::default()
            },
        ];
        let mate = mem_aln_t {
            rid: 0,
            pos: 6,
            is_rev: 1,
            n_cigar: 1,
            cigar: vec![4_u32 << 4],
            ..Default::default()
        };
        let mut out = kstring_t::default();
        let mut opt = (*opt).clone();
        opt.flag |= MEM_F_REF_HDR;
        mem_aln2sam(
            &opt,
            &bns,
            &mut out,
            &s,
            list.len() as i32,
            &list,
            0,
            Some(&mate),
        );
        let sam = out.as_str();
        assert!(
            sam.starts_with("read1\t33\tchr1\t3\t60\t4M\t=\t7\t8\tACGT\t!!!!"),
            "{sam}"
        );
        assert!(
            sam.contains("\tNM:i:0\tMD:Z:4\tMC:Z:4M\tAS:i:4\tXS:i:1"),
            "{sam}"
        );
        assert!(sam.contains("\tSA:Z:chr1,6,+,4M,20,0;"), "{sam}");
        assert!(sam.contains("\tpa:f:0.500"), "{sam}");
        assert!(sam.contains("\tCO:Z:note"), "{sam}");
        assert!(sam.contains("\tXR:Z:anno field"), "{sam}");
    }

    #[test]
    fn mem_aln2sam_emits_reverse_supplementary_sequence_and_hard_clipped_cigar() {
        let opt = mem_opt_init();
        let bns = bntseq_t {
            l_pac: 8,
            n_seqs: 1,
            anns: vec![crate::generated::bntseq_h::bntann1_t {
                offset: 0,
                len: 8,
                name: "chr1".into(),
                ..Default::default()
            }],
            ..Default::default()
        };
        let s = bseq1_t {
            l_seq: 6,
            name: Some("read2".into()),
            seq: Some("AACCGT".into()),
            qual: Some("123456".into()),
            ..Default::default()
        };
        let list = vec![
            mem_aln_t {
                rid: 0,
                pos: 0,
                mapq: 60,
                n_cigar: 1,
                cigar: vec![6_u32 << 4],
                md: b"6".to_vec(),
                NM: 0,
                score: 6,
                ..Default::default()
            },
            mem_aln_t {
                rid: 0,
                pos: 1,
                is_rev: 1,
                mapq: 30,
                n_cigar: 3,
                cigar: vec![(1_u32 << 4) | 3, (4_u32 << 4), (1_u32 << 4) | 3],
                md: b"4".to_vec(),
                NM: 0,
                score: 4,
                ..Default::default()
            },
        ];
        let mut out = kstring_t::default();
        mem_aln2sam(&opt, &bns, &mut out, &s, list.len() as i32, &list, 1, None);
        let sam = out.as_str();
        assert!(
            sam.starts_with("read2\t16\tchr1\t2\t30\t1H4M1H\t*\t0\t0\tCGGT\t5432"),
            "{sam}"
        );
        assert!(sam.contains("\tNM:i:0\tMD:Z:4\tAS:i:4\tXS:i:0"), "{sam}");
    }

    #[test]
    fn mem_reg2sam_writes_unmapped_record_when_no_region_passes_threshold() {
        let opt = mem_opt_init();
        let bns = bntseq_t {
            anns: vec![crate::generated::bntseq_h::bntann1_t {
                offset: 0,
                len: 4,
                name: "chr1".into(),
                ..Default::default()
            }],
            ..Default::default()
        };
        let mut s = bseq1_t {
            l_seq: 4,
            name: Some("read3".into()),
            seq: Some("ACGT".into()),
            qual: Some("!!!!".into()),
            ..Default::default()
        };
        let mut regs = mem_alnreg_v {
            n: 1,
            m: 1,
            // secondary_all=-1 matches the value mem_mark_primary_se sets before mem_reg2sam runs;
            // the default 0 would point get_pri_idx at this entry (a self-secondary), which never happens in real flow.
            a: vec![mem_alnreg_t {
                score: 5,
                secondary: -1,
                secondary_all: -1,
                ..Default::default()
            }],
        };

        mem_reg2sam(&opt, &bns, &[], &mut s, &mut regs, 0, None);
        let sam = s.sam.as_deref().expect("sam");
        assert!(
            sam.starts_with("read3\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t!!!!"),
            "{sam}"
        );
    }

    #[test]
    fn mem_reg2sam_writes_multiple_alignments_and_marks_supplementary() {
        let mut opt = (*mem_opt_init()).clone();
        opt.flag |= MEM_F_ALL;
        opt.T = 1;
        let forward = vec![0_u8, 1, 2, 3, 0, 1, 2, 3];
        let pac = pack_seq(&forward);
        let bns = bntseq_t {
            l_pac: 8,
            n_seqs: 1,
            anns: vec![crate::generated::bntseq_h::bntann1_t {
                offset: 0,
                len: 8,
                name: "chr1".into(),
                ..Default::default()
            }],
            ..Default::default()
        };
        let mut s = bseq1_t {
            l_seq: 4,
            name: Some("read4".into()),
            seq: Some("ACGT".into()),
            qual: Some("####".into()),
            ..Default::default()
        };
        let mut regs = mem_alnreg_v {
            n: 2,
            m: 2,
            a: vec![
                mem_alnreg_t {
                    rb: 0,
                    re: 4,
                    qb: 0,
                    qe: 4,
                    rid: 0,
                    score: 4,
                    truesc: 4,
                    sub: 1,
                    csub: 1,
                    secondary: -1,
                    w: 100,
                    ..Default::default()
                },
                mem_alnreg_t {
                    rb: 4,
                    re: 8,
                    qb: 0,
                    qe: 4,
                    rid: 0,
                    score: 4,
                    truesc: 4,
                    sub: 1,
                    csub: 1,
                    secondary: -1,
                    w: 100,
                    ..Default::default()
                },
            ],
        };

        mem_reg2sam(&opt, &bns, &pac, &mut s, &mut regs, 0, None);
        let sam = s.sam.as_deref().expect("sam");
        let lines: Vec<&str> = sam.trim_end().split('\n').collect();
        assert_eq!(lines.len(), 2, "{sam}");
        assert!(lines[0].starts_with("read4\t0\tchr1\t1\t"), "{sam}");
        assert!(lines[1].starts_with("read4\t2048\tchr1\t5\t"), "{sam}");
    }

    #[test]
    fn worker_sam_marks_primary_and_writes_sam_for_batch_window() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_worker_sam_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("ref");
        let fasta = b">chr1\nACGTACGTACGT\n";
        bns_fasta2bntseq(
            Cursor::new(fasta.as_slice()),
            prefix.to_str().expect("utf8"),
            1,
        );

        let mut fmi = FMI_search::ctor(prefix.to_str().expect("utf8"));
        assert_eq!(fmi.build_index(), 0);
        fmi.load_index();

        let seqs = vec![bseq1_t {
            l_seq: 6,
            name: Some("read5".into()),
            seq: Some("ACGTAC".into()),
            qual: Some("IIIIII".into()),
            ..Default::default()
        }];
        let mut opt = (*mem_opt_init()).clone();
        opt.min_seed_len = 2;
        opt.split_factor = 1.5;
        opt.split_width = 10;
        opt.max_mem_intv = 20;
        opt.min_chain_weight = 1;
        opt.T = 1;

        let mut worker = worker_t {
            opt: Some(Box::new(opt.clone())),
            seqs,
            chain_ar: vec![mem_chain_v::default(); 1],
            regs: vec![mem_alnreg_v::default(); 1],
            seedBuf: vec![mem_seed_t::default(); 256],
            seedBufSize: 256,
            mmc: mem_cache::default(),
            ref_string: pac_to_reference_layout(12, &fmi.base.idx.pac),
            fmi: Some(fmi),
            ..Default::default()
        };

        worker_bwt(&mut worker, 0, 1, 0);
        worker_aln(&mut worker, 0, 1, 0);
        assert!(worker.regs[0].n > 0);
        worker_sam(&mut worker, 0, 1, 0);
        let sam = worker.seqs[0].sam.as_deref().expect("sam");
        assert!(sam.starts_with("read5\t"), "{sam}");
        assert!(sam.contains("\tchr1\t"), "{sam}");

        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn mem_process_seqs_runs_single_end_pipeline_end_to_end() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_process_seqs_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("ref");
        let fasta = b">chr1\nACGTACGTACGT\n";
        bns_fasta2bntseq(
            Cursor::new(fasta.as_slice()),
            prefix.to_str().expect("utf8"),
            1,
        );

        let mut fmi = FMI_search::ctor(prefix.to_str().expect("utf8"));
        assert_eq!(fmi.build_index(), 0);
        fmi.load_index();

        let mut seqs = vec![bseq1_t {
            l_seq: 6,
            name: Some("read6".into()),
            seq: Some("ACGTAC".into()),
            qual: Some("IIIIII".into()),
            ..Default::default()
        }];
        let mut opt = (*mem_opt_init()).clone();
        opt.min_seed_len = 2;
        opt.split_factor = 1.5;
        opt.split_width = 10;
        opt.max_mem_intv = 20;
        opt.min_chain_weight = 1;
        opt.T = 1;

        let mut worker = worker_t {
            chain_ar: vec![mem_chain_v::default(); 1],
            regs: vec![mem_alnreg_v::default(); 1],
            seedBuf: vec![mem_seed_t::default(); 256],
            seedBufSize: 256,
            mmc: mem_cache::default(),
            fmi: Some(fmi),
            ..Default::default()
        };

        mem_process_seqs(&mut opt, 0, 1, &mut seqs, None, &mut worker);
        let sam = seqs[0].sam.as_deref().expect("sam");
        assert!(sam.starts_with("read6\t"), "{sam}");
        assert!(sam.contains("\tchr1\t"), "{sam}");

        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn mem_process_seqs_paired_end_copies_pestat_before_stub_boundary() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_process_seqs_pe_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("ref");
        let fasta = b">chr1\nACGTACGTACGT\n";
        bns_fasta2bntseq(
            Cursor::new(fasta.as_slice()),
            prefix.to_str().expect("utf8"),
            1,
        );

        let mut fmi = FMI_search::ctor(prefix.to_str().expect("utf8"));
        assert_eq!(fmi.build_index(), 0);
        fmi.load_index();

        let mut seqs = vec![
            bseq1_t {
                l_seq: 6,
                name: Some("read7".into()),
                seq: Some("ACGTAC".into()),
                qual: Some("IIIIII".into()),
                ..Default::default()
            },
            bseq1_t {
                l_seq: 6,
                name: Some("read7".into()),
                seq: Some("GTACGT".into()),
                qual: Some("IIIIII".into()),
                ..Default::default()
            },
        ];
        let mut opt = (*mem_opt_init()).clone();
        opt.flag |= MEM_F_PE;
        opt.min_seed_len = 2;
        opt.split_factor = 1.5;
        opt.split_width = 10;
        opt.max_mem_intv = 20;
        opt.min_chain_weight = 1;
        opt.T = 1;

        let pes0 = [
            mem_pestat_t {
                low: 10,
                high: 20,
                failed: 0,
                avg: 15.0,
                std: 1.5,
            },
            mem_pestat_t {
                low: 0,
                high: 0,
                failed: 1,
                avg: 0.0,
                std: 0.0,
            },
            mem_pestat_t {
                low: 0,
                high: 0,
                failed: 1,
                avg: 0.0,
                std: 0.0,
            },
            mem_pestat_t {
                low: 0,
                high: 0,
                failed: 1,
                avg: 0.0,
                std: 0.0,
            },
        ];

        let aux = ktp_aux_t {
            opt: Some(Box::new(opt.clone())),
            ..Default::default()
        };
        let mut worker = worker_t {
            fmi: Some(fmi),
            ..Default::default()
        };
        memoryAlloc(&aux, &mut worker, 2, opt.n_threads);

        mem_process_seqs(&mut opt, 0, 2, &mut seqs, Some(&pes0), &mut worker);
        assert_eq!(worker.pes[0].low, 10);
        assert_eq!(worker.pes[0].high, 20);
        assert_eq!(worker.pes[0].failed, 0);
        assert!(seqs[0]
            .sam
            .as_deref()
            .is_some_and(|sam| sam.contains("read7")));
        assert!(seqs[1]
            .sam
            .as_deref()
            .is_some_and(|sam| sam.contains("read7")));

        fs::remove_dir_all(&dir).expect("cleanup");
    }
}
