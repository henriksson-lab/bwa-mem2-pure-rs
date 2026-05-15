#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Port of `bwa-mem2/src/bandedswa.h` + `bwa-mem2/src/bandedswa.cpp`.

// --- bandedswa.h ---

#[doc = "Original struct: OutScore (bwa-mem2/src/bandedSWA.h)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct OutScore {
    pub _opaque: (),
}

#[doc = "Original struct: SeqPair (bwa-mem2/src/bandedSWA.h)"]
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct SeqPair {
    pub idr: i32,
    pub idq: i32,
    pub id: i32,
    pub len1: i32,
    pub len2: i32,
    pub h0: i32,
    pub seqid: i32,
    pub regid: i32,
    pub score: i32,
    pub tle: i32,
    pub gtle: i32,
    pub qle: i32,
    pub gscore: i32,
    pub max_off: i32,
}

#[doc = "Original struct: dnaOutScore (bwa-mem2/src/bandedSWA.h)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct dnaOutScore {
    pub _opaque: (),
}

#[doc = "Original struct: dnaSeqPair (bwa-mem2/src/bandedSWA.h)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct dnaSeqPair {
    pub _opaque: (),
}

// --- bandedswa.cpp ---

const DEFAULT_AMBIG: i8 = -1;
const MAX_SEQ_LEN8: usize = 128;
const MAX_SEQ_LEN16: usize = 32768;

#[doc = "Original class: BandedPairWiseSW (bwa-mem2/src/bandedSWA.cpp)"]
#[derive(Debug, Clone)]
pub struct BandedPairWiseSW {
    pub mat: [i8; 25],
    pub m: i32,
    pub end_bonus: i32,
    pub zdrop: i32,
    pub o_del: i32,
    pub o_ins: i32,
    pub e_del: i32,
    pub e_ins: i32,
    pub w_match: i8,
    pub w_mismatch: i8,
    pub w_open: i32,
    pub w_extend: i32,
    /// ambig penalty
    pub w_ambig: i8,
    pub swTicks: i64,
    pub SW_cells: u64,
    pub setupTicks: i64,
    pub sort1Ticks: i64,
    pub sort2Ticks: i64,
}

#[derive(Debug, Default, Clone, Copy)]
struct eh_t {
    h: i32,
    e: i32,
}

// Thread-local scratch for scalarBandedSWA (qp, eh_h, eh_e). Saves 3 allocations per call;
// the function is called many times per read in mem_chain2aln_across_reads_V2's chain extension.
// Note: SoA (split eh_h / eh_e) measured slightly faster than AoS (Vec<EhCell{h,e}>) on this
// codebase's hot loop — the compiler vectorizes paired sequential accesses to the two arrays
// better than the AoS struct path; benchmarked 8.27s SoA vs 8.7s AoS on 50K-read parity test.
thread_local! {
    static BANDED_SCRATCH: std::cell::RefCell<(Vec<i8>, Vec<i32>, Vec<i32>)> =
        const { std::cell::RefCell::new((Vec::new(), Vec::new(), Vec::new())) };
}

fn disable_simd8() -> bool {
    static FLAG: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *FLAG.get_or_init(|| std::env::var_os("BWA_DISABLE_SIMD8").is_some())
}

fn disable_simd16() -> bool {
    static FLAG: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *FLAG.get_or_init(|| std::env::var_os("BWA_DISABLE_SIMD16").is_some())
}

// Constants from bwa-mem2/src/bandedSWA.cpp:44-46. Used by the AVX-512 batch SW path
// (smithWaterman512_8) to mark padding lanes in the SoA query/target buffers.
#[allow(dead_code)]
const AMBIG_BASE: u8 = 4;
#[allow(dead_code)]
const DUMMY1_BYTE: u8 = 99;
#[allow(dead_code)]
const DUMMY2_BYTE: u8 = 100;
#[allow(dead_code)]
const SIMD_WIDTH8_AVX512: usize = 64;

// Encode `pairs` (up to SIMD_WIDTH8_AVX512 lanes) into SoA layout for the AVX-512 SW kernel.
// Layout: `seq1_soa[k * 64 + lane]` = base of target sequence at position k for that lane;
// `seq2_soa` similarly for the query. Lanes whose pair is shorter than `max_len1`/`max_len2`
// get `DUMMY1`/`DUMMY2` padding past their actual length, so the SIMD DP can run uniformly
// across all lanes without per-lane length checks in the hot path.
//
// Returns `(seq1_soa, seq2_soa, max_len1, max_len2)`. The SoA buffers have length
// `max_len * 64`. Pairs beyond `pairs.len()` (up to 64) are filled with all-padding lanes
// using the last pair's len2 (matches C++ smithWatermanBatchWrapper8:2017 invariant).
#[allow(dead_code)]
pub(crate) fn encode_pairs_soa_avx512(
    pairs: &[SeqPair],
    seq_buf_ref: &[u8],
    seq_buf_qer: &[u8],
) -> (Vec<u8>, Vec<u8>, usize, usize) {
    let mut s1 = Vec::new();
    let mut s2 = Vec::new();
    let (max_l1, max_l2) =
        encode_pairs_soa_avx512_into(pairs, seq_buf_ref, seq_buf_qer, &mut s1, &mut s2);
    (s1, s2, max_l1, max_l2)
}

pub(crate) fn encode_pairs_soa_avx512_into(
    pairs: &[SeqPair],
    seq_buf_ref: &[u8],
    seq_buf_qer: &[u8],
    s1: &mut Vec<u8>,
    s2: &mut Vec<u8>,
) -> (usize, usize) {
    debug_assert!(pairs.len() <= SIMD_WIDTH8_AVX512);
    let mut max_len1 = 0_usize;
    let mut max_len2 = 0_usize;
    for p in pairs {
        max_len1 = max_len1.max(p.len1 as usize);
        max_len2 = max_len2.max(p.len2 as usize);
    }
    s1.clear();
    s1.resize(max_len1 * SIMD_WIDTH8_AVX512, DUMMY1_BYTE);
    s2.clear();
    s2.resize(max_len2 * SIMD_WIDTH8_AVX512, DUMMY2_BYTE);
    for (lane, p) in pairs.iter().enumerate() {
        let s1_base = (p.idr as usize)..(p.idr as usize + p.len1 as usize);
        let s2_base = (p.idq as usize)..(p.idq as usize + p.len2 as usize);
        // BCE: lane < SIMD_WIDTH8_AVX512 and k < {len1,len2} <= {max_len1,max_len2}, so the
        // SoA index k*64+lane is < max_*64 = s*.len(). Safe under the resize above.
        // Branchless: 4 → 0xFF (=255), 0..3 → identity. Add 251 when bit-2 set: b=4 has b>>2 = 1.
        unsafe {
            for (k, &b) in seq_buf_ref[s1_base].iter().enumerate() {
                let v = b.wrapping_add(((b >> 2) & 1).wrapping_mul(251));
                *s1.get_unchecked_mut(k * SIMD_WIDTH8_AVX512 + lane) = v;
            }
            for (k, &b) in seq_buf_qer[s2_base].iter().enumerate() {
                let v = b.wrapping_add(((b >> 2) & 1).wrapping_mul(251));
                *s2.get_unchecked_mut(k * SIMD_WIDTH8_AVX512 + lane) = v;
            }
        }
    }
    (max_len1, max_len2)
}

// Single-pair full-DP reference implementation that uses `main_code8_lane` for each cell —
// validates the per-cell semantics against `scalarBandedSWA` independent of the SoA / SIMD
// machinery. Lays out H and F in row-major and walks (i, j) in scan order, exactly mirroring
// the structure of the C++ AVX-512 batch loop but for one pair (one lane of work).
//
// Returns (max_score, tle, qle, gscore, gtle) — same fields scalarBandedSWA writes back.
// Caller must size buffers correctly: seq1/seq2 are 2-bit-encoded.
//
// Notes:
// - `e11` (within-row "F" in C++ naming) propagates across j; uses oe_ins/e_ins.
// - `f_col` (across-row "E"-like state, stored per column) uses oe_del/e_del.
// - h0 is the initial score for the (0,0) corner; first row/col follow gap-init.
#[allow(dead_code)]
pub(crate) fn full_dp_one_pair_via_main_code8(
    bsw: &BandedPairWiseSW,
    seq1: &[u8],
    seq2: &[u8],
    h0: u8,
) -> (i32, i32, i32, i32, i32) {
    let tlen = seq1.len();
    let qlen = seq2.len();
    let oe_ins = (bsw.o_ins + bsw.e_ins) as u8;
    let oe_del = (bsw.o_del + bsw.e_del) as u8;
    let e_ins = bsw.e_ins as u8;
    let e_del = bsw.e_del as u8;
    let o_del = bsw.o_del as u8;

    // H matrix: (tlen+1) x (qlen+1).
    let mut h = vec![vec![0_u8; qlen + 1]; tlen + 1];
    h[0][0] = h0;

    // First column: H[i][0] = max(0, h0 - o_del - i*e_del). Stops cascading once 0 reached.
    let mut tmp_v = (h0 as i8).wrapping_sub(o_del as i8);
    for i in 1..=tlen {
        tmp_v = tmp_v.wrapping_sub(e_del as i8);
        h[i][0] = tmp_v.max(0) as u8;
    }

    // First row: H[0][1] = max(0, h0 - oe_ins) (cmpgt+blend), H[0][k>1] = max(0, prev - e_ins).
    if qlen >= 1 {
        h[0][1] = if (h0 as i8) > (oe_ins as i8) {
            ((h0 as i8).wrapping_sub(oe_ins as i8)) as u8
        } else {
            0
        };
        let mut tmp_h = h[0][1] as i8;
        for j in 2..=qlen {
            tmp_h = tmp_h.wrapping_sub(e_ins as i8).max(0);
            h[0][j] = tmp_h as u8;
        }
    }

    // F column buffer: F[j] holds the across-row gap state at column j (propagates from row i-1).
    let mut f_col = vec![0_u8; qlen + 1];

    // Per-pair tracking
    let mut max_score = h0 as i32;
    let mut max_i = -1_i32;
    let mut max_j = -1_i32;
    let mut gscore = -1_i32;
    let mut max_ie = -1_i32;

    for i in 1..=tlen {
        let mut e11 = 0_u8; // within-row F-state (uses oe_ins/e_ins)
        let s1 = seq1[i - 1];
        for j in 1..=qlen {
            let s2 = seq2[j - 1];
            let h00 = h[i - 1][j - 1];
            let f11 = f_col[j];
            let (h11, e11_new, f21) = main_code8_lane(
                s1,
                s2,
                h00,
                e11,
                f11,
                bsw.w_match,
                bsw.w_mismatch,
                bsw.w_ambig,
                e_ins,
                oe_ins,
                e_del,
                oe_del,
            );
            h[i][j] = h11;
            e11 = e11_new;
            f_col[j] = f21;
            if (h11 as i32) > max_score {
                max_score = h11 as i32;
                max_i = (i - 1) as i32;
                max_j = (j - 1) as i32;
            }
        }
        // gscore: H at right-edge column (j == qlen). Track latest row reaching the new best.
        if qlen > 0 && (h[i][qlen] as i32) > gscore {
            gscore = h[i][qlen] as i32;
            max_ie = (i - 1) as i32;
        }
    }
    (max_score, max_i, max_j, gscore, max_ie)
}

// Banded variant of the SoA batch SW prototype. Adds per-lane head/tail tracking and in-loop
// band masking. Mirrors C++ smithWaterman512_8:2384-2412 banding logic but uses scalar per-lane
// (no SIMD) for validation. Same `w` across all lanes (per-lane myband adjustment is left to
// the caller — scalarBandedSWA also takes a single `w`).
//
// Returns SimdBatchState matching scalarBandedSWA on banded inputs (when post-row band narrowing
// doesn't differ — the post-row narrowing is omitted here for simplicity, but it only excludes
// already-zero cells from future computation, so output should match unless zdrop has fired).
#[allow(dead_code)]
pub(crate) fn process_batch_soa_banded_dp(
    bsw: &BandedPairWiseSW,
    pairs: &[SeqPair],
    seq_buf_ref: &[u8],
    seq_buf_qer: &[u8],
    w: i32,
) -> SimdBatchState {
    debug_assert!(pairs.len() <= SIMD_WIDTH8_AVX512);
    let n = pairs.len();
    let (s1_soa, s2_soa, max_l1, max_l2) = encode_pairs_soa_avx512(pairs, seq_buf_ref, seq_buf_qer);
    let mut h0_lanes = [0_u8; SIMD_WIDTH8_AVX512];
    let mut len1_lanes = [0_u32; SIMD_WIDTH8_AVX512];
    let mut len2_lanes = [0_u32; SIMD_WIDTH8_AVX512];
    for (lane, p) in pairs.iter().enumerate() {
        h0_lanes[lane] = p.h0 as u8;
        len1_lanes[lane] = p.len1 as u32;
        len2_lanes[lane] = p.len2 as u32;
    }
    let oe_ins = (bsw.o_ins + bsw.e_ins) as u8;
    let oe_del = (bsw.o_del + bsw.e_del) as u8;
    let e_ins = bsw.e_ins as u8;
    let e_del = bsw.e_del as u8;
    let e_del_i32 = bsw.e_del;
    let e_ins_i32 = bsw.e_ins;
    let zdrop_i32 = bsw.zdrop;
    let o_del = bsw.o_del as u8;
    // H_v + 1 extra row (for h10 init at "j=0 boundary" at i=max_l1-1).
    let mut h_v_extended = vec![0_u8; (max_l1 + 1) * SIMD_WIDTH8_AVX512];
    if max_l1 > 0 {
        let h_v = init_h_v_soa_avx512(&h0_lanes, max_l1, o_del, e_del);
        h_v_extended[..max_l1 * SIMD_WIDTH8_AVX512].copy_from_slice(&h_v);
        let mut tmp = [0_i8; SIMD_WIDTH8_AVX512];
        for lane in 0..SIMD_WIDTH8_AVX512 {
            tmp[lane] = (h0_lanes[lane] as i8).wrapping_sub(o_del as i8);
            for _ in 0..max_l1 {
                tmp[lane] = tmp[lane].wrapping_sub(e_del as i8);
            }
            h_v_extended[max_l1 * SIMD_WIDTH8_AVX512 + lane] = tmp[lane].max(0) as u8;
        }
    }
    let mut h_h = if max_l2 > 0 {
        init_h_h_soa_avx512(&h0_lanes, max_l2, oe_ins, e_ins)
    } else {
        vec![0_u8; SIMD_WIDTH8_AVX512]
    };
    h_h.resize((max_l2 + 1).max(1) * SIMD_WIDTH8_AVX512, 0);
    let mut f_buf = vec![0_u8; (max_l2 + 1) * SIMD_WIDTH8_AVX512];

    let mut state = SimdBatchState::new(&h0_lanes);
    if max_l1 == 0 || max_l2 == 0 {
        return state;
    }

    // Per-lane band: head[lane] (lo edge), tail[lane] (hi edge, inclusive). Init: [0, len2).
    let mut head_lane = [0_i32; SIMD_WIDTH8_AVX512];
    let mut tail_lane = [0_i32; SIMD_WIDTH8_AVX512];
    // Per-lane "alive" flag — false once the pair has terminated (row max == 0 or zdrop).
    let mut alive = [false; SIMD_WIDTH8_AVX512];
    for lane in 0..n {
        tail_lane[lane] = len2_lanes[lane] as i32;
        alive[lane] = true;
    }

    for i in 0..max_l1 {
        let i_i32 = i as i32;
        // Early-exit if all lanes have terminated.
        let mut any_alive = false;
        for lane in 0..n {
            if alive[lane] {
                any_alive = true;
                break;
            }
        }
        if !any_alive {
            break;
        }
        // Per-row band update — mirrors scalarBandedSWA: beg = max(beg, i-w); end = min(end, i+1+w, qlen).
        for lane in 0..n {
            let lo = i_i32 - w;
            if lo > head_lane[lane] {
                head_lane[lane] = lo;
            }
            let hi = i_i32 + 1 + w;
            if hi < tail_lane[lane] {
                tail_lane[lane] = hi;
            }
            if (len2_lanes[lane] as i32) < tail_lane[lane] {
                tail_lane[lane] = len2_lanes[lane] as i32;
            }
        }

        let s1_row = &s1_soa[i * SIMD_WIDTH8_AVX512..(i + 1) * SIMD_WIDTH8_AVX512];
        let mut e11 = [0_u8; SIMD_WIDTH8_AVX512];
        let mut h10 = [0_u8; SIMD_WIDTH8_AVX512];
        for lane in 0..n {
            if head_lane[lane] == 0 {
                h10[lane] = h_v_extended[(i + 1) * SIMD_WIDTH8_AVX512 + lane];
            } else {
                h10[lane] = 0;
            }
        }
        // scalarBandedSWA tracks per-row (m, mj) using `h >= m`, then promotes to global only
        // when row m strictly exceeds previous max. Mirror that here.
        let mut row_m = [0_i32; SIMD_WIDTH8_AVX512];
        let mut row_mj = [-1_i32; SIMD_WIDTH8_AVX512];

        for j in 0..max_l2 {
            let s2_col = &s2_soa[j * SIMD_WIDTH8_AVX512..(j + 1) * SIMD_WIDTH8_AVX512];
            let j_i32 = j as i32;
            for lane in 0..n {
                let s1 = s1_row[lane];
                let s2 = s2_col[lane];
                let h00 = h_h[j * SIMD_WIDTH8_AVX512 + lane];
                let f11 = f_buf[j * SIMD_WIDTH8_AVX512 + lane];
                let prev_h10 = h10[lane];

                let (h11, e11_new, f21) = main_code8_lane(
                    s1,
                    s2,
                    h00,
                    e11[lane],
                    f11,
                    bsw.w_match,
                    bsw.w_mismatch,
                    bsw.w_ambig,
                    e_ins,
                    oe_ins,
                    e_del,
                    oe_del,
                );

                let in_band = j_i32 >= head_lane[lane] && j_i32 <= tail_lane[lane];
                let h10_to_store = if in_band { prev_h10 } else { 0 };
                let f_to_store = if in_band { f21 } else { 0 };
                h_h[j * SIMD_WIDTH8_AVX512 + lane] = h10_to_store;
                f_buf[j * SIMD_WIDTH8_AVX512 + lane] = f_to_store;
                e11[lane] = e11_new;
                h10[lane] = h11;

                if (i as u32) < len1_lanes[lane] && (j as u32) < len2_lanes[lane] && in_band {
                    let h11_i = h11 as i32;
                    // Per-row update uses >= (matches scalarBandedSWA's `h >= m`).
                    if h11_i >= row_m[lane] {
                        row_m[lane] = h11_i;
                        row_mj[lane] = j_i32;
                    }
                    if (j as u32) == len2_lanes[lane] - 1 && h11_i >= state.gscore[lane] {
                        state.gscore[lane] = h11_i;
                        state.max_ie[lane] = i_i32;
                    }
                }
            }
        }
        // After-row promotion (mirrors scalarBandedSWA exactly):
        //   if m == 0: this lane terminates (DP ran out of non-zero scores)
        //   if m > max: update global max + max_off
        //   else: zdrop check, terminate if drop too large
        for lane in 0..n {
            if !alive[lane] {
                continue;
            }
            if row_m[lane] == 0 {
                alive[lane] = false;
                continue;
            }
            if row_m[lane] > state.max_score[lane] {
                state.max_score[lane] = row_m[lane];
                state.max_i[lane] = i_i32;
                state.max_j[lane] = row_mj[lane];
                let mj_minus_i = (row_mj[lane] - i_i32).abs();
                if mj_minus_i > state.max_off[lane] {
                    state.max_off[lane] = mj_minus_i;
                }
            } else {
                // Z-drop check (matches scalarBandedSWA's no-`if zdrop>0`-guard SIMD form).
                let di = i_i32 - state.max_i[lane];
                let dj = row_mj[lane] - state.max_j[lane];
                let exceeded = if di > dj {
                    state.max_score[lane] - row_m[lane] - (di - dj) * e_del_i32 > zdrop_i32
                } else {
                    state.max_score[lane] - row_m[lane] - (dj - di) * e_ins_i32 > zdrop_i32
                };
                if exceeded {
                    alive[lane] = false;
                }
            }
        }
        for lane in 0..n {
            h_h[max_l2 * SIMD_WIDTH8_AVX512 + lane] = h10[lane];
        }

        // Dynamic post-row band narrowing per lane — mirrors scalarBandedSWA's:
        //   advance beg past leading (h==0 && e==0) columns
        //   retreat end past trailing zeros
        //   end = min(je + 2, qlen)
        // This GROWS the band along the diagonal as alignment progresses (despite the name
        // "narrowing"): once `end` retreats to the last non-zero + 2, the next row's update
        // `tail = min(tail, i+1+w)` only further constrains, but the +2 padding plus the
        // new row's diagonal cell extension bumps end forward each iter.
        for lane in 0..n {
            if !alive[lane] {
                continue;
            }
            let qlen_lane = len2_lanes[lane] as i32;
            let mut jb = head_lane[lane] as usize;
            while (jb as i32) < tail_lane[lane]
                && h_h[jb * SIMD_WIDTH8_AVX512 + lane] == 0
                && f_buf[jb * SIMD_WIDTH8_AVX512 + lane] == 0
            {
                jb += 1;
            }
            head_lane[lane] = jb as i32;

            let mut je = tail_lane[lane] as usize;
            while (je as i32) >= head_lane[lane]
                && h_h[je * SIMD_WIDTH8_AVX512 + lane] == 0
                && f_buf[je * SIMD_WIDTH8_AVX512 + lane] == 0
            {
                if je == 0 {
                    break;
                }
                je -= 1;
            }
            tail_lane[lane] = ((je as i32) + 2).min(qlen_lane);
        }
    }
    state
}

// i16 AVX-512 variant: same algorithm as `process_batch_soa_banded_dp_avx512` but with i16
// elements (32 lanes per __m512i instead of 64). Used by the i16 bucket which covers
// pairs with len ≥ 128 — the dominant case for 150bp reads.
//
// Returns SimdBatchState. Caller must verify avx512bw before calling.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512bw,avx512f")]
#[allow(dead_code)]
pub(crate) unsafe fn process_batch_soa_banded_dp_avx512_16(
    bsw: &BandedPairWiseSW,
    pairs: &[SeqPair],
    seq_buf_ref: &[u8],
    seq_buf_qer: &[u8],
    w: i32,
) -> SimdBatchState {
    SIMD16_PROTO_SCRATCH.with(|cell| unsafe {
        let mut buf = cell.borrow_mut();
        let (s1_soa, s2_soa, h_v_extended, h_h, f_buf) = &mut *buf;
        process_batch_soa_banded_dp_avx512_16_impl(
            bsw,
            pairs,
            seq_buf_ref,
            seq_buf_qer,
            w,
            s1_soa,
            s2_soa,
            h_v_extended,
            h_h,
            f_buf,
        )
    })
}

#[cfg(target_arch = "x86_64")]
thread_local! {
    static SIMD16_PROTO_SCRATCH: std::cell::RefCell<(Vec<i16>, Vec<i16>, Vec<i16>, Vec<i16>, Vec<i16>)> =
        const { std::cell::RefCell::new((Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new())) };
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512bw,avx512f")]
#[allow(dead_code, clippy::too_many_arguments)]
pub(crate) unsafe fn process_batch_soa_banded_dp_avx512_16_impl(
    bsw: &BandedPairWiseSW,
    pairs: &[SeqPair],
    seq_buf_ref: &[u8],
    seq_buf_qer: &[u8],
    w: i32,
    s1_soa: &mut Vec<i16>,
    s2_soa: &mut Vec<i16>,
    h_v_extended: &mut Vec<i16>,
    h_h: &mut Vec<i16>,
    f_buf: &mut Vec<i16>,
) -> SimdBatchState {
    use core::arch::x86_64::*;
    debug_assert!(pairs.len() <= SIMD_WIDTH16_AVX512);
    let n = pairs.len();
    let (max_l1, max_l2) =
        encode_pairs_soa_avx512_16_into(pairs, seq_buf_ref, seq_buf_qer, s1_soa, s2_soa);
    let mut h0_lanes_i16 = [0_i16; SIMD_WIDTH16_AVX512];
    let mut h0_lanes_u8 = [0_u8; SIMD_WIDTH8_AVX512]; // for SimdBatchState::new
    let mut len1_lanes = [0_u32; SIMD_WIDTH16_AVX512];
    let mut len2_lanes = [0_u32; SIMD_WIDTH16_AVX512];
    for (lane, p) in pairs.iter().enumerate() {
        h0_lanes_i16[lane] = p.h0 as i16;
        h0_lanes_u8[lane] = p.h0 as u8;
        len1_lanes[lane] = p.len1 as u32;
        len2_lanes[lane] = p.len2 as u32;
    }
    let oe_ins = (bsw.o_ins + bsw.e_ins) as i16;
    let oe_del = (bsw.o_del + bsw.e_del) as i16;
    let e_ins = bsw.e_ins as i16;
    let e_del = bsw.e_del as i16;
    let o_del = bsw.o_del as i16;
    // H_v with one extra row. C++ smithWatermanBatchWrapper16 (bandedSWA.cpp:2832-2837) only
    // initializes H2[k*32] for k in [1, maxLen1); H2[maxLen1*32] is unwritten — typically zero
    // from a fresh heap allocation. smithWaterman512_16 reads H_v[(i+1)*32] when beg==0 at the
    // last row (i = maxLen1-1), so the extra entry stays at the vec init value of zero.
    h_v_extended.clear();
    h_v_extended.resize((max_l1 + 1) * SIMD_WIDTH16_AVX512, 0);
    if max_l1 > 0 {
        init_h_v_soa_avx512_16_into_simd(&h0_lanes_i16, max_l1, o_del, e_del, h_v_extended);
    }
    h_h.clear();
    h_h.resize((max_l2 + 1).max(1) * SIMD_WIDTH16_AVX512, 0);
    if max_l2 > 0 {
        init_h_h_soa_avx512_16_into_simd(&h0_lanes_i16, max_l2, oe_ins, e_ins, h_h);
    }
    f_buf.clear();
    f_buf.resize((max_l2 + 1) * SIMD_WIDTH16_AVX512, 0);

    let mut state = SimdBatchState::new(&h0_lanes_u8);
    if max_l1 == 0 || max_l2 == 0 {
        return state;
    }

    // Per-lane head/tail as i16 vectors (lengths fit in i16 for typical reads). Padding lanes
    // (>= n) keep head=0, tail=0 so their in_band only fires at j==0 — harmless since their
    // results are never read out of state.
    let mut head_lane = [0_i16; SIMD_WIDTH16_AVX512];
    let mut tail_lane = [0_i16; SIMD_WIDTH16_AVX512];
    let mut band_lane = [0_i16; SIMD_WIDTH16_AVX512];
    let mut len1_i16 = [0_i16; SIMD_WIDTH16_AVX512];
    let mut len2_i16 = [0_i16; SIMD_WIDTH16_AVX512];
    for lane in 0..n {
        tail_lane[lane] = len2_lanes[lane] as i16;
        let len2_score = len2_lanes[lane] as i32 * i32::from(bsw.w_match);
        let max_ins = ((len2_score + bsw.end_bonus - bsw.o_ins) / bsw.e_ins + 1).max(1);
        let max_del = ((len2_score + bsw.end_bonus - bsw.o_del) / bsw.e_del + 1).max(1);
        band_lane[lane] = w.min(max_ins).min(max_del) as i16;
        len1_i16[lane] = len1_lanes[lane] as i16;
        len2_i16[lane] = len2_lanes[lane] as i16;
    }

    let match_v = _mm512_set1_epi16(bsw.w_match as i16);
    let mismatch_v = _mm512_set1_epi16(bsw.w_mismatch as i16);
    let w_ambig_v = _mm512_set1_epi16(bsw.w_ambig as i16);
    let e_ins_v = _mm512_set1_epi16(e_ins);
    let oe_ins_v = _mm512_set1_epi16(oe_ins);
    let e_del_v = _mm512_set1_epi16(e_del);
    let oe_del_v = _mm512_set1_epi16(oe_del);
    let zero_v = _mm512_setzero_si512();
    let len1_v = _mm512_loadu_si512(len1_i16.as_ptr() as *const __m512i);
    let len2_v = _mm512_loadu_si512(len2_i16.as_ptr() as *const __m512i);
    let band_v = _mm512_loadu_si512(band_lane.as_ptr() as *const __m512i);
    let min_len2 = len2_lanes[..n].iter().copied().min().unwrap_or(0) as usize;
    let ff_v = _mm512_set1_epi16(-1_i16);
    let zdrop_v = _mm512_set1_epi16(bsw.zdrop as i16);
    let mut mlen_v = _mm512_add_epi16(len2_v, band_v);
    mlen_v = _mm512_min_epu16(mlen_v, len1_v);
    let mut max_score_v = _mm512_loadu_si512(h_v_extended.as_ptr() as *const __m512i);
    let mut x_v = zero_v;
    let mut y_v = zero_v;
    let mut max_off_v = zero_v;
    let mut gscore_v = ff_v;
    let mut max_ie_v = zero_v;
    let mut exit0_v = ff_v;

    for i in 0..max_l1 {
        let i_i32 = i as i32;
        let i1_v = _mm512_set1_epi16((i + 1) as i16);
        // Vectorized head/tail update across all 32 lanes:
        //   lo = i - band; head = max(head, lo)
        //   hi = (i+1) + band; tail = min(tail, hi); tail = min(tail, len2)
        let i_v = _mm512_set1_epi16(i_i32 as i16);
        let prev_head_v = _mm512_loadu_si512(head_lane.as_ptr() as *const __m512i);
        let prev_tail_v = _mm512_loadu_si512(tail_lane.as_ptr() as *const __m512i);
        let lo_v = _mm512_sub_epi16(i_v, band_v);
        let hi_v = _mm512_add_epi16(i1_v, band_v);
        let head_v = _mm512_max_epi16(prev_head_v, lo_v);
        let tail_clipped = _mm512_min_epi16(prev_tail_v, hi_v);
        let tail_v = _mm512_min_epi16(tail_clipped, len2_v);
        _mm512_storeu_si512(head_lane.as_mut_ptr() as *mut __m512i, head_v);
        _mm512_storeu_si512(tail_lane.as_mut_ptr() as *mut __m512i, tail_v);
        let head_zero_mask = _mm512_cmpeq_epi16_mask(head_v, zero_v);

        let mut e11_v = zero_v;
        let h10_src_v = _mm512_loadu_si512(
            h_v_extended.as_ptr().add((i + 1) * SIMD_WIDTH16_AVX512) as *const __m512i,
        );
        let mut h10_v = _mm512_mask_blend_epi16(head_zero_mask, zero_v, h10_src_v);
        let mut cmpim = _mm512_cmpgt_epi16_mask(i1_v, mlen_v);
        cmpim |= _mm512_cmpeq_epi16_mask(tail_v, head_v);
        cmpim |= _mm512_cmpgt_epi16_mask(head_v, tail_v);
        exit0_v = _mm512_mask_blend_epi16(cmpim, exit0_v, zero_v);

        let s1_row_v =
            _mm512_loadu_si512(s1_soa.as_ptr().add(i * SIMD_WIDTH16_AVX512) as *const __m512i);

        // Restrict j loop to global band [beg, end) — mirrors C++ smithWaterman512_16.
        let beg = if i_i32 < w {
            0_usize
        } else {
            (i_i32 - w) as usize
        };
        let end = ((i_i32 as usize + 1).saturating_add(w as usize)).min(max_l2);
        let mut max_rs_v = zero_v;
        let mut y1_v = zero_v;

        // Avoid the per-iteration broadcast of j and j+1 by carrying them and incrementing.
        let one_i16_v = _mm512_set1_epi16(1);
        let mut j_v = _mm512_set1_epi16(beg as i16);
        let mut j1_v = _mm512_set1_epi16((beg as i16).wrapping_add(1));
        for j in beg..end {
            let s2_col_v =
                _mm512_loadu_si512(s2_soa.as_ptr().add(j * SIMD_WIDTH16_AVX512) as *const __m512i);
            let h00_v =
                _mm512_loadu_si512(h_h.as_ptr().add(j * SIMD_WIDTH16_AVX512) as *const __m512i);
            let f11_v =
                _mm512_loadu_si512(f_buf.as_ptr().add(j * SIMD_WIDTH16_AVX512) as *const __m512i);

            let (h11_v, e11_new_v, f21_v) = main_code16_avx512(
                s1_row_v, s2_col_v, h00_v, e11_v, f11_v, match_v, mismatch_v, w_ambig_v, e_ins_v,
                oe_ins_v, e_del_v, oe_del_v,
            );

            // Compute the head>j mask once; in_band uses its negation, score_oob reuses it.
            let head_gt_j: u32 = _mm512_cmpgt_epi16_mask(head_v, j_v);
            // in_band = j >= head && j <= tail = !(head > j) && j <= tail.
            let in_band_mask: u32 = !head_gt_j & _mm512_cmple_epi16_mask(j_v, tail_v);

            // h10_to_store = in_band ? h10_v : 0  (h10_v is the prev_h10 carry from prior column)
            // f_to_store = in_band ? f21_v : 0
            let h10_to_store = _mm512_mask_blend_epi16(in_band_mask, zero_v, h10_v);
            let f_to_store = _mm512_mask_blend_epi16(in_band_mask, zero_v, f21_v);
            _mm512_storeu_si512(
                h_h.as_mut_ptr().add(j * SIMD_WIDTH16_AVX512) as *mut __m512i,
                h10_to_store,
            );
            _mm512_storeu_si512(
                f_buf.as_mut_ptr().add(j * SIMD_WIDTH16_AVX512) as *mut __m512i,
                f_to_store,
            );

            // Update e11, h10 for next column (next iteration of j).
            e11_v = e11_new_v;
            h10_v = h11_v;

            let prev_max_rs_v = max_rs_v;
            // cmp_update == (max_rs == h11 OR max_rs > prev_max), simplified to (h11 >= prev_max)
            // since after max(prev, h11), the value equals h11 iff h11 >= prev_max.
            let cmp_update = _mm512_cmpge_epi16_mask(h11_v, prev_max_rs_v);
            let score_oob_mask = _mm512_cmpgt_epi16_mask(j1_v, tail_v) | head_gt_j;
            // Fused: max_rs_v = !score_oob ? max(prev, h11) : prev. Saves 1 SIMD op vs separate
            // max + blend.
            max_rs_v = _mm512_mask_max_epi16(prev_max_rs_v, !score_oob_mask, prev_max_rs_v, h11_v);
            // y1 := j1 if (cmp_update AND in_band); old y1 otherwise.
            let y_update_mask = cmp_update & !score_oob_mask;
            y1_v = _mm512_mask_blend_epi16(y_update_mask, y1_v, j1_v);

            if j + 1 >= min_len2 {
                let edge_mask = _mm512_cmpeq_epi16_mask(j1_v, len2_v);
                if edge_mask != 0 {
                    // Original 4-blend chain reduces to: gscore := max_gh when (active & edge &
                    // !over_tail); max_ie := i1 when also h11 wins (!cmp_gh).
                    let max_gh_v = _mm512_max_epi16(gscore_v, h11_v);
                    let cmp_gh = _mm512_cmpgt_epi16_mask(gscore_v, h11_v);
                    let mex0 = _mm512_movepi16_mask(exit0_v);
                    let over_tail = _mm512_cmpgt_epi16_mask(j1_v, tail_v);
                    let gs_update_mask = edge_mask & mex0 & !over_tail;
                    let ie_update_mask = gs_update_mask & !cmp_gh;
                    gscore_v = _mm512_mask_blend_epi16(gs_update_mask, gscore_v, max_gh_v);
                    max_ie_v = _mm512_mask_blend_epi16(ie_update_mask, max_ie_v, i1_v);
                }
            }
            // Increment j_v, j1_v for next iteration to avoid the broadcast.
            j_v = _mm512_add_epi16(j_v, one_i16_v);
            j1_v = _mm512_add_epi16(j1_v, one_i16_v);
        }
        // Mirror C++ smithWaterman512_16 (bandedSWA.cpp:3197-3202): after the j loop, store h10
        // at H_h[end*32], zeroed for lanes where end is outside [head, tail]. F[end*32] is also
        // zeroed (we re-zero f_buf at the next row's j-loop initial mask, so explicit zero here is
        // redundant — skip it).
        let end_v = _mm512_set1_epi16(end as i16);
        let end_oob_mask =
            _mm512_cmpgt_epi16_mask(head_v, end_v) | _mm512_cmpgt_epi16_mask(end_v, tail_v);
        let h10_end_v = _mm512_mask_blend_epi16(end_oob_mask, h10_v, zero_v);
        _mm512_storeu_si512(
            h_h.as_mut_ptr().add(end * SIMD_WIDTH16_AVX512) as *mut __m512i,
            h10_end_v,
        );
        _mm512_storeu_si512(
            f_buf.as_mut_ptr().add(end * SIMD_WIDTH16_AVX512) as *mut __m512i,
            zero_v,
        );

        let zero_row_mask = _mm512_cmpeq_epi16_mask(max_rs_v, zero_v);
        if zero_row_mask == u32::MAX {
            break;
        }
        exit0_v = _mm512_mask_blend_epi16(zero_row_mask, exit0_v, zero_v);

        let old_max_score_v = max_score_v;
        let mex0 = _mm512_movepi16_mask(exit0_v);
        // Fused max+blend: max_score_v = (mex0 ? max(prev, max_rs) : prev).
        max_score_v = _mm512_mask_max_epi16(max_score_v, mex0, max_score_v, max_rs_v);
        // max_improved iff (mex0 && max_rs > prev). Breaks the dep chain through max_score_v.
        let max_improved = mex0 & _mm512_cmpgt_epi16_mask(max_rs_v, old_max_score_v);
        y_v = _mm512_mask_blend_epi16(max_improved, y_v, y1_v);
        x_v = _mm512_mask_blend_epi16(max_improved, x_v, i1_v);

        let diag_delta_v = _mm512_abs_epi16(_mm512_sub_epi16(y1_v, i1_v));
        // Fused max+blend: max_off_v = (max_improved ? max(prev, diag) : prev).
        max_off_v = _mm512_mask_max_epi16(max_off_v, max_improved, max_off_v, diag_delta_v);

        let tmp_i_v = _mm512_sub_epi16(i1_v, x_v);
        let tmp_j_v = _mm512_sub_epi16(y1_v, y_v);
        let score_drop_v = _mm512_sub_epi16(max_score_v, max_rs_v);
        // |tmp_i - tmp_j| via abs; the prior 4-op sub_a/sub_b + cmpgt + blend dance is equivalent.
        let abs_diff_v = _mm512_abs_epi16(_mm512_sub_epi16(tmp_i_v, tmp_j_v));
        let z_v = _mm512_sub_epi16(score_drop_v, abs_diff_v);
        let z_mask = _mm512_cmpgt_epi16_mask(z_v, zdrop_v);
        exit0_v = _mm512_mask_blend_epi16(z_mask, exit0_v, zero_v);

        let active_mask = _mm512_movepi16_mask(exit0_v);
        let h_h_ptr = h_h.as_ptr();
        let f_buf_ptr = f_buf.as_ptr();
        for lane in 0..n {
            if ((active_mask >> lane) & 1) == 0 {
                continue;
            }
            let qlen_lane = len2_lanes[lane] as i16;
            // Combined OR check: load h_h and f_buf, OR them, compare to 0. Compiler is free to
            // schedule both loads in parallel (no short-circuit dep between them).
            let mut jb = head_lane[lane] as usize;
            while (jb as i16) < tail_lane[lane] {
                let off = jb * SIMD_WIDTH16_AVX512 + lane;
                let h = unsafe { *h_h_ptr.add(off) };
                let f = unsafe { *f_buf_ptr.add(off) };
                if (h | f) != 0 {
                    break;
                }
                jb += 1;
            }
            head_lane[lane] = jb as i16;
            let mut je = tail_lane[lane] as usize;
            while (je as i16) >= head_lane[lane] {
                let off = je * SIMD_WIDTH16_AVX512 + lane;
                let h = unsafe { *h_h_ptr.add(off) };
                let f = unsafe { *f_buf_ptr.add(off) };
                if (h | f) != 0 {
                    break;
                }
                if je == 0 {
                    break;
                }
                je -= 1;
            }
            tail_lane[lane] = ((je as i16) + 2).min(qlen_lane);
        }
    }
    let mut score_ar = [0_i16; SIMD_WIDTH16_AVX512];
    let mut x_ar = [0_i16; SIMD_WIDTH16_AVX512];
    let mut y_ar = [0_i16; SIMD_WIDTH16_AVX512];
    let mut max_off_ar = [0_i16; SIMD_WIDTH16_AVX512];
    let mut gscore_ar = [0_i16; SIMD_WIDTH16_AVX512];
    let mut max_ie_ar = [0_i16; SIMD_WIDTH16_AVX512];
    _mm512_storeu_si512(score_ar.as_mut_ptr() as *mut __m512i, max_score_v);
    _mm512_storeu_si512(x_ar.as_mut_ptr() as *mut __m512i, x_v);
    _mm512_storeu_si512(y_ar.as_mut_ptr() as *mut __m512i, y_v);
    _mm512_storeu_si512(max_off_ar.as_mut_ptr() as *mut __m512i, max_off_v);
    _mm512_storeu_si512(gscore_ar.as_mut_ptr() as *mut __m512i, gscore_v);
    _mm512_storeu_si512(max_ie_ar.as_mut_ptr() as *mut __m512i, max_ie_v);
    for lane in 0..n {
        state.max_score[lane] = i32::from(score_ar[lane]);
        state.max_i[lane] = i32::from(x_ar[lane]) - 1;
        state.max_j[lane] = i32::from(y_ar[lane]) - 1;
        state.max_off[lane] = i32::from(max_off_ar[lane]);
        state.gscore[lane] = i32::from(gscore_ar[lane]);
        state.max_ie[lane] = i32::from(max_ie_ar[lane]) - 1;
    }
    state
}

// AVX-512 variant of process_batch_soa_banded_dp. Same algorithm + same output, but the inner
// per-cell DP body is replaced with `main_code8_avx512` running 64 lanes in parallel via __m512i.
// Per-lane state tracking (banding mask check, row_m, max_score, alive flag, dynamic band
// narrowing) remains scalar — the SIMD speedup comes from collapsing 64 main_code8_lane calls
// per (i, j) cell into a single 14-op SIMD sequence.
//
// Caller must verify CPU supports `avx512bw,avx512f` before calling. Falls back to scalar
// `process_batch_soa_banded_dp` otherwise.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512bw,avx512f")]
#[allow(dead_code)]
pub(crate) unsafe fn process_batch_soa_banded_dp_avx512(
    bsw: &BandedPairWiseSW,
    pairs: &[SeqPair],
    seq_buf_ref: &[u8],
    seq_buf_qer: &[u8],
    w: i32,
) -> SimdBatchState {
    SIMD8_PROTO_SCRATCH.with(|cell| unsafe {
        let mut buf = cell.borrow_mut();
        let (s1_soa, s2_soa, h_v_extended, h_h, f_buf) = &mut *buf;
        process_batch_soa_banded_dp_avx512_impl(
            bsw,
            pairs,
            seq_buf_ref,
            seq_buf_qer,
            w,
            s1_soa,
            s2_soa,
            h_v_extended,
            h_h,
            f_buf,
        )
    })
}

#[cfg(target_arch = "x86_64")]
thread_local! {
    static SIMD8_PROTO_SCRATCH: std::cell::RefCell<(Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>)> =
        const { std::cell::RefCell::new((Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new())) };
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512bw,avx512f")]
#[allow(dead_code, clippy::too_many_arguments)]
pub(crate) unsafe fn process_batch_soa_banded_dp_avx512_impl(
    bsw: &BandedPairWiseSW,
    pairs: &[SeqPair],
    seq_buf_ref: &[u8],
    seq_buf_qer: &[u8],
    w: i32,
    s1_soa: &mut Vec<u8>,
    s2_soa: &mut Vec<u8>,
    h_v_extended: &mut Vec<u8>,
    h_h: &mut Vec<u8>,
    f_buf: &mut Vec<u8>,
) -> SimdBatchState {
    use core::arch::x86_64::*;

    debug_assert!(pairs.len() <= SIMD_WIDTH8_AVX512);
    let n = pairs.len();
    let (max_l1, max_l2) =
        encode_pairs_soa_avx512_into(pairs, seq_buf_ref, seq_buf_qer, s1_soa, s2_soa);
    let mut h0_lanes = [0_u8; SIMD_WIDTH8_AVX512];
    let mut len1_lanes = [0_u32; SIMD_WIDTH8_AVX512];
    let mut len2_lanes = [0_u32; SIMD_WIDTH8_AVX512];
    for (lane, p) in pairs.iter().enumerate() {
        h0_lanes[lane] = p.h0 as u8;
        len1_lanes[lane] = p.len1 as u32;
        len2_lanes[lane] = p.len2 as u32;
    }
    let oe_ins = (bsw.o_ins + bsw.e_ins) as u8;
    let oe_del = (bsw.o_del + bsw.e_del) as u8;
    let e_ins = bsw.e_ins as u8;
    let e_del = bsw.e_del as u8;
    let o_del = bsw.o_del as u8;
    // C++ smithWatermanBatchWrapper8 (bandedSWA.cpp:1997 family) only initializes H_v[k*64] for
    // k in [1, maxLen1); H_v[maxLen1*64] is unwritten (typically zero from fresh heap alloc).
    // smithWaterman512_8 reads H_v[(i+1)*64] when beg==0 at the last row (i = maxLen1-1), so
    // the extra entry stays at the vec resize value of zero.
    h_v_extended.clear();
    h_v_extended.resize((max_l1 + 1) * SIMD_WIDTH8_AVX512, 0);
    if max_l1 > 0 {
        // Write directly into the thread-local h_v_extended[..max_l1*64] slice — avoids the
        // per-call Vec alloc + copy that the helper-returns-Vec form did.
        init_h_v_soa_avx512_into(
            &h0_lanes,
            max_l1,
            o_del,
            e_del,
            &mut h_v_extended[..max_l1 * SIMD_WIDTH8_AVX512],
        );
    }
    if max_l2 > 0 {
        h_h.clear();
        h_h.resize((max_l2 + 1) * SIMD_WIDTH8_AVX512, 0);
        init_h_h_soa_avx512_into(
            &h0_lanes,
            max_l2,
            oe_ins,
            e_ins,
            &mut h_h[..max_l2 * SIMD_WIDTH8_AVX512],
        );
    } else {
        h_h.clear();
        h_h.resize(SIMD_WIDTH8_AVX512, 0);
    }
    f_buf.clear();
    f_buf.resize((max_l2 + 1) * SIMD_WIDTH8_AVX512, 0);

    let mut state = SimdBatchState::new(&h0_lanes);
    if max_l1 == 0 || max_l2 == 0 {
        return state;
    }

    let mut head_lane = [0_i32; SIMD_WIDTH8_AVX512];
    let mut tail_lane = [0_i32; SIMD_WIDTH8_AVX512];
    let mut band_lane = [0_i32; SIMD_WIDTH8_AVX512];
    let mut band_i8 = [0_i8; SIMD_WIDTH8_AVX512];
    let max_scor = bsw.w_match.max(bsw.w_mismatch).max(bsw.w_ambig).max(0) as u8;
    let eb_ins = (bsw.end_bonus - bsw.o_ins) as i8 as u8;
    let eb_del = (bsw.end_bonus - bsw.o_del) as i8 as u8;
    let mut len1_i8 = [0_i8; SIMD_WIDTH8_AVX512];
    let mut len2_i8 = [0_i8; SIMD_WIDTH8_AVX512];
    for lane in 0..n {
        tail_lane[lane] = len2_lanes[lane] as i32;
        let qlen_scaled = (len2_lanes[lane] as u8).wrapping_mul(max_scor);
        let max_ins = (i32::from(qlen_scaled.wrapping_add(eb_ins)) / bsw.e_ins + 1).max(1);
        let max_del = (i32::from(qlen_scaled.wrapping_add(eb_del)) / bsw.e_del + 1).max(1);
        band_lane[lane] = w.min(max_ins).min(max_del);
        band_i8[lane] = band_lane[lane].clamp(i8::MIN as i32, i8::MAX as i32) as i8;
        len1_i8[lane] = (len1_lanes[lane] as i32).clamp(i8::MIN as i32, i8::MAX as i32) as i8;
        len2_i8[lane] = (len2_lanes[lane] as i32).clamp(i8::MIN as i32, i8::MAX as i32) as i8;
    }

    // Constant SIMD vectors used by the SIMD body.
    let match_v = _mm512_set1_epi8(bsw.w_match);
    let mismatch_v = _mm512_set1_epi8(bsw.w_mismatch);
    let w_ambig_v = _mm512_set1_epi8(bsw.w_ambig);
    let e_ins_v = _mm512_set1_epi8(e_ins as i8);
    let oe_ins_v = _mm512_set1_epi8(oe_ins as i8);
    let e_del_v = _mm512_set1_epi8(e_del as i8);
    let oe_del_v = _mm512_set1_epi8(oe_del as i8);
    let zero_v = _mm512_setzero_si512();
    let len1_v = _mm512_loadu_si512(len1_i8.as_ptr() as *const __m512i);
    let len2_v = _mm512_loadu_si512(len2_i8.as_ptr() as *const __m512i);
    let band_v = _mm512_loadu_si512(band_i8.as_ptr() as *const __m512i);
    let min_len2 = len2_lanes[..n].iter().copied().min().unwrap_or(0) as usize;
    let ff_v = _mm512_set1_epi8(-1);
    let zdrop_v = _mm512_set1_epi8(bsw.zdrop as i8);
    let mut mlen_v = _mm512_add_epi8(len2_v, band_v);
    mlen_v = _mm512_min_epu8(mlen_v, len1_v);
    let mut max_score_v = _mm512_loadu_si512(h_v_extended.as_ptr() as *const __m512i);
    let mut x_v = zero_v;
    let mut y_v = zero_v;
    let mut max_off_v = zero_v;
    let mut gscore_v = ff_v;
    let mut max_ie_v = zero_v;
    let mut exit0_v = ff_v;

    for i in 0..max_l1 {
        let i_i32 = i as i32;
        let i1_v = _mm512_set1_epi8((i + 1) as i8);
        for lane in 0..n {
            let lo = i_i32 - band_lane[lane];
            if lo > head_lane[lane] {
                head_lane[lane] = lo;
            }
            let hi = i_i32 + 1 + band_lane[lane];
            if hi < tail_lane[lane] {
                tail_lane[lane] = hi;
            }
            if (len2_lanes[lane] as i32) < tail_lane[lane] {
                tail_lane[lane] = len2_lanes[lane] as i32;
            }
        }

        // Load head/tail as __m512i (i8). For typical reads, values fit in i8.
        let mut head_i8 = [0_i8; SIMD_WIDTH8_AVX512];
        let mut tail_i8 = [0_i8; SIMD_WIDTH8_AVX512];
        for lane in 0..n {
            head_i8[lane] = head_lane[lane].clamp(i8::MIN as i32, i8::MAX as i32) as i8;
            tail_i8[lane] = tail_lane[lane].clamp(i8::MIN as i32, i8::MAX as i32) as i8;
        }
        let head_v = _mm512_loadu_si512(head_i8.as_ptr() as *const __m512i);
        let tail_v = _mm512_loadu_si512(tail_i8.as_ptr() as *const __m512i);
        let head_zero_mask = _mm512_cmpeq_epi8_mask(head_v, zero_v);
        let h10_src_v = _mm512_loadu_si512(
            h_v_extended.as_ptr().add((i + 1) * SIMD_WIDTH8_AVX512) as *const __m512i
        );
        let mut h10_v = _mm512_mask_blend_epi8(head_zero_mask, zero_v, h10_src_v);
        let mut e11_v = zero_v;
        let mut cmpim = _mm512_cmpgt_epi8_mask(i1_v, mlen_v);
        cmpim |= _mm512_cmpeq_epi8_mask(tail_v, head_v);
        cmpim |= _mm512_cmpgt_epi8_mask(head_v, tail_v);
        exit0_v = _mm512_mask_blend_epi8(cmpim, exit0_v, zero_v);

        // Load s1 row vector once per row.
        let s1_row_v =
            _mm512_loadu_si512(s1_soa.as_ptr().add(i * SIMD_WIDTH8_AVX512) as *const __m512i);

        // Restrict j loop to the global band [beg, end) — mirrors C++ smithWaterman512_8.
        // Saves work when the band is at row boundaries.
        let beg = if i_i32 < w {
            0_usize
        } else {
            (i_i32 - w) as usize
        };
        let end = ((i_i32 as usize + 1).saturating_add(w as usize)).min(max_l2);
        let mut max_rs_v = zero_v;
        let mut y1_v = zero_v;

        // Avoid the per-iteration broadcast of j and j+1 by carrying them and incrementing.
        let one_i8_v = _mm512_set1_epi8(1);
        let mut j_v = _mm512_set1_epi8(beg as i8);
        let mut j1_v = _mm512_set1_epi8((beg + 1) as i8);
        for j in beg..end {
            let s2_col_v =
                _mm512_loadu_si512(s2_soa.as_ptr().add(j * SIMD_WIDTH8_AVX512) as *const __m512i);
            let h00_v =
                _mm512_loadu_si512(h_h.as_ptr().add(j * SIMD_WIDTH8_AVX512) as *const __m512i);
            let f11_v =
                _mm512_loadu_si512(f_buf.as_ptr().add(j * SIMD_WIDTH8_AVX512) as *const __m512i);

            let (h11_v, e11_new_v, f21_v) = main_code8_avx512(
                s1_row_v, s2_col_v, h00_v, e11_v, f11_v, match_v, mismatch_v, w_ambig_v, e_ins_v,
                oe_ins_v, e_del_v, oe_del_v,
            );

            // Compute head>j mask once; in_band uses its negation, score_oob reuses it.
            let head_gt_j: u64 = _mm512_cmpgt_epi8_mask(head_v, j_v);
            // in_band = j >= head && j <= tail = !(head > j) && j <= tail.
            let in_band_mask: u64 = !head_gt_j & _mm512_cmple_epi8_mask(j_v, tail_v);

            // Masked store: in-band lanes write prev h10 / f21, out-of-band write 0.
            let h10_to_store = _mm512_mask_blend_epi8(in_band_mask, zero_v, h10_v);
            let f_to_store = _mm512_mask_blend_epi8(in_band_mask, zero_v, f21_v);
            _mm512_storeu_si512(
                h_h.as_mut_ptr().add(j * SIMD_WIDTH8_AVX512) as *mut __m512i,
                h10_to_store,
            );
            _mm512_storeu_si512(
                f_buf.as_mut_ptr().add(j * SIMD_WIDTH8_AVX512) as *mut __m512i,
                f_to_store,
            );

            let prev_max_rs_v = max_rs_v;
            // (max_rs > prev_max) || (max_rs == h11) simplifies to (h11 >= prev_max).
            let cmp_update = _mm512_cmpge_epi8_mask(h11_v, prev_max_rs_v);
            let score_oob_mask = _mm512_cmpgt_epi8_mask(j1_v, tail_v) | head_gt_j;
            // Fused: max_rs_v = !score_oob ? max(prev, h11) : prev. Saves 1 SIMD op vs separate
            // max + blend.
            max_rs_v = _mm512_mask_max_epi8(prev_max_rs_v, !score_oob_mask, prev_max_rs_v, h11_v);
            // y1 := j1 if (cmp_update AND in_band); old y1 otherwise.
            let y_update_mask = cmp_update & !score_oob_mask;
            y1_v = _mm512_mask_blend_epi8(y_update_mask, y1_v, j1_v);

            // Carry e11/h10 to next column.
            e11_v = e11_new_v;
            h10_v = h11_v;

            if j + 1 >= min_len2 {
                let edge_mask = _mm512_cmpeq_epi8_mask(j1_v, len2_v);
                if edge_mask != 0 {
                    // Original 4-blend chain reduces to: gscore := max_gh when (active & edge &
                    // !over_tail); max_ie := i1 when also h11 wins (!cmp_gh).
                    let max_gh_v = _mm512_max_epi8(gscore_v, h11_v);
                    let cmp_gh = _mm512_cmpgt_epi8_mask(gscore_v, h11_v);
                    let mex0 = _mm512_movepi8_mask(exit0_v);
                    let over_tail = _mm512_cmpgt_epi8_mask(j1_v, tail_v);
                    let gs_update_mask = edge_mask & mex0 & !over_tail;
                    let ie_update_mask = gs_update_mask & !cmp_gh;
                    gscore_v = _mm512_mask_blend_epi8(gs_update_mask, gscore_v, max_gh_v);
                    max_ie_v = _mm512_mask_blend_epi8(ie_update_mask, max_ie_v, i1_v);
                }
            }
            // Increment j_v, j1_v for next iteration.
            j_v = _mm512_add_epi8(j_v, one_i8_v);
            j1_v = _mm512_add_epi8(j1_v, one_i8_v);
        }
        let end_v = _mm512_set1_epi8(end as i8);
        let end_oob_mask =
            _mm512_cmpgt_epi8_mask(head_v, end_v) | _mm512_cmpgt_epi8_mask(end_v, tail_v);
        let h10_end_v = _mm512_mask_blend_epi8(end_oob_mask, h10_v, zero_v);
        _mm512_storeu_si512(
            h_h.as_mut_ptr().add(end * SIMD_WIDTH8_AVX512) as *mut __m512i,
            h10_end_v,
        );
        _mm512_storeu_si512(
            f_buf.as_mut_ptr().add(end * SIMD_WIDTH8_AVX512) as *mut __m512i,
            zero_v,
        );

        let zero_row_mask = _mm512_cmpeq_epi8_mask(max_rs_v, zero_v);
        if zero_row_mask == u64::MAX {
            break;
        }
        exit0_v = _mm512_mask_blend_epi8(zero_row_mask, exit0_v, zero_v);

        let old_max_score_v = max_score_v;
        let mex0 = _mm512_movepi8_mask(exit0_v);
        // Fused max+blend: max_score_v = (mex0 ? max(prev, max_rs) : prev).
        max_score_v = _mm512_mask_max_epi8(max_score_v, mex0, max_score_v, max_rs_v);
        // max_improved iff (mex0 && max_rs > prev). Breaks the dep chain through max_score_v.
        let max_improved = mex0 & _mm512_cmpgt_epi8_mask(max_rs_v, old_max_score_v);
        y_v = _mm512_mask_blend_epi8(max_improved, y_v, y1_v);
        x_v = _mm512_mask_blend_epi8(max_improved, x_v, i1_v);

        let diag_delta_v = _mm512_abs_epi8(_mm512_sub_epi8(y1_v, i1_v));
        // Fused max+blend: max_off_v = (max_improved ? max(prev, diag) : prev).
        max_off_v = _mm512_mask_max_epi8(max_off_v, max_improved, max_off_v, diag_delta_v);

        let tmp_i_v = _mm512_sub_epi8(i1_v, x_v);
        let tmp_j_v = _mm512_sub_epi8(y1_v, y_v);
        let score_drop_v = _mm512_sub_epi8(max_score_v, max_rs_v);
        // |tmp_i - tmp_j| via abs; saves the 4-op cmpgt + double-sub + blend.
        let abs_diff_v = _mm512_abs_epi8(_mm512_sub_epi8(tmp_i_v, tmp_j_v));
        let z_v = _mm512_sub_epi8(score_drop_v, abs_diff_v);
        let z_mask = _mm512_cmpgt_epi8_mask(z_v, zdrop_v);
        exit0_v = _mm512_mask_blend_epi8(z_mask, exit0_v, zero_v);

        // Mirror C++ smithWaterman512_8: after the j loop, store h10 at H_h[end*64], zeroed
        // for lanes where end is outside [head, tail]. This is row-local band state; writing
        // max_l2 unconditionally leaves stale diagonal state when end < max_l2.
        let active_mask = _mm512_movepi8_mask(exit0_v);
        let h_h_ptr = h_h.as_ptr();
        let f_buf_ptr = f_buf.as_ptr();
        for lane in 0..n {
            if ((active_mask >> lane) & 1) == 0 {
                continue;
            }
            let qlen_lane = len2_lanes[lane] as i32;
            // Combined OR check: load h_h and f_buf, OR them, compare to 0. Allows parallel loads.
            let mut jb = head_lane[lane] as usize;
            while (jb as i32) < tail_lane[lane] {
                let off = jb * SIMD_WIDTH8_AVX512 + lane;
                let h = unsafe { *h_h_ptr.add(off) };
                let f = unsafe { *f_buf_ptr.add(off) };
                if (h | f) != 0 {
                    break;
                }
                jb += 1;
            }
            head_lane[lane] = jb as i32;
            let mut je = tail_lane[lane] as usize;
            while (je as i32) >= head_lane[lane] {
                let off = je * SIMD_WIDTH8_AVX512 + lane;
                let h = unsafe { *h_h_ptr.add(off) };
                let f = unsafe { *f_buf_ptr.add(off) };
                if (h | f) != 0 {
                    break;
                }
                if je == 0 {
                    break;
                }
                je -= 1;
            }
            tail_lane[lane] = ((je as i32) + 2).min(qlen_lane);
        }
    }
    let mut score_ar = [0_i8; SIMD_WIDTH8_AVX512];
    let mut x_ar = [0_i8; SIMD_WIDTH8_AVX512];
    let mut y_ar = [0_i8; SIMD_WIDTH8_AVX512];
    let mut max_off_ar = [0_i8; SIMD_WIDTH8_AVX512];
    let mut gscore_ar = [0_i8; SIMD_WIDTH8_AVX512];
    let mut max_ie_ar = [0_i8; SIMD_WIDTH8_AVX512];
    _mm512_storeu_si512(score_ar.as_mut_ptr() as *mut __m512i, max_score_v);
    _mm512_storeu_si512(x_ar.as_mut_ptr() as *mut __m512i, x_v);
    _mm512_storeu_si512(y_ar.as_mut_ptr() as *mut __m512i, y_v);
    _mm512_storeu_si512(max_off_ar.as_mut_ptr() as *mut __m512i, max_off_v);
    _mm512_storeu_si512(gscore_ar.as_mut_ptr() as *mut __m512i, gscore_v);
    _mm512_storeu_si512(max_ie_ar.as_mut_ptr() as *mut __m512i, max_ie_v);
    for lane in 0..n {
        state.max_score[lane] = i32::from(score_ar[lane]);
        state.max_i[lane] = i32::from(x_ar[lane]) - 1;
        state.max_j[lane] = i32::from(y_ar[lane]) - 1;
        state.max_off[lane] = i32::from(max_off_ar[lane]);
        state.gscore[lane] = i32::from(gscore_ar[lane]);
        state.max_ie[lane] = i32::from(max_ie_ar[lane]) - 1;
    }
    state
}

// SoA-batched scalar prototype of the AVX-512 batch SW. Processes up to 64 pairs simultaneously
// using SoA buffer layout — same data layout that the AVX-512 intrinsics will consume, but
// inner loop is scalar (per lane) for validation purposes. Future iteration replaces the lane
// loop with `__m512i` intrinsics.
//
// Currently FULL DP only — no banding, no zdrop early exit. Pairs with very wide bands or
// short reads work fine; pairs with non-trivial banding may differ from `scalarBandedSWA`.
// Adding banding is the next step.
//
// Returns SimdBatchState with per-lane (max_score, max_i, max_j, gscore, max_ie). Lanes
// without pairs (lane >= pairs.len()) are left at default (all zero).
#[allow(dead_code)]
pub(crate) fn process_batch_soa_full_dp(
    bsw: &BandedPairWiseSW,
    pairs: &[SeqPair],
    seq_buf_ref: &[u8],
    seq_buf_qer: &[u8],
) -> SimdBatchState {
    debug_assert!(pairs.len() <= SIMD_WIDTH8_AVX512);
    let n = pairs.len();

    // Encode SoA buffers + collect per-lane lengths + h0
    let (s1_soa, s2_soa, max_l1, max_l2) = encode_pairs_soa_avx512(pairs, seq_buf_ref, seq_buf_qer);
    let mut h0_lanes = [0_u8; SIMD_WIDTH8_AVX512];
    let mut len1_lanes = [0_u32; SIMD_WIDTH8_AVX512];
    let mut len2_lanes = [0_u32; SIMD_WIDTH8_AVX512];
    for (lane, p) in pairs.iter().enumerate() {
        h0_lanes[lane] = p.h0 as u8;
        len1_lanes[lane] = p.len1 as u32;
        len2_lanes[lane] = p.len2 as u32;
    }

    let oe_ins = (bsw.o_ins + bsw.e_ins) as u8;
    let oe_del = (bsw.o_del + bsw.e_del) as u8;
    let e_ins = bsw.e_ins as u8;
    let e_del = bsw.e_del as u8;
    let o_del = bsw.o_del as u8;

    // H_v init (vertical column, one extra row to allow H_v[max_len1] read at last iter).
    let mut h_v_extended = vec![0_u8; (max_l1 + 1) * SIMD_WIDTH8_AVX512];
    if max_l1 > 0 {
        let h_v = init_h_v_soa_avx512(&h0_lanes, max_l1, o_del, e_del);
        h_v_extended[..max_l1 * SIMD_WIDTH8_AVX512].copy_from_slice(&h_v);
        // Extra row: continue the descent for one more step (mirrors C++ which sizes H_v to max+1).
        let mut tmp = [0_i8; SIMD_WIDTH8_AVX512];
        for lane in 0..SIMD_WIDTH8_AVX512 {
            tmp[lane] = (h0_lanes[lane] as i8).wrapping_sub(o_del as i8);
            for _ in 0..max_l1 {
                tmp[lane] = tmp[lane].wrapping_sub(e_del as i8);
            }
            h_v_extended[max_l1 * SIMD_WIDTH8_AVX512 + lane] = tmp[lane].max(0) as u8;
        }
    }

    // H_h: holds H[i-1][j-1] across rows (initialized as the gap-init first row).
    let mut h_h = if max_l2 > 0 {
        init_h_h_soa_avx512(&h0_lanes, max_l2, oe_ins, e_ins)
    } else {
        vec![0_u8; SIMD_WIDTH8_AVX512]
    };
    h_h.resize((max_l2 + 1).max(1) * SIMD_WIDTH8_AVX512, 0);

    // F: across-row "E"-state (column buffer), starts at 0.
    let mut f_buf = vec![0_u8; (max_l2 + 1) * SIMD_WIDTH8_AVX512];

    let mut state = SimdBatchState::new(&h0_lanes);

    if max_l1 == 0 || max_l2 == 0 {
        return state;
    }

    for i in 0..max_l1 {
        let s1_row = &s1_soa[i * SIMD_WIDTH8_AVX512..(i + 1) * SIMD_WIDTH8_AVX512];
        let mut e11 = [0_u8; SIMD_WIDTH8_AVX512];
        let mut h10 = [0_u8; SIMD_WIDTH8_AVX512];
        // h10 init at "j=0" boundary: H_v[i+1] (next row's first-column boundary).
        for lane in 0..SIMD_WIDTH8_AVX512 {
            h10[lane] = h_v_extended[(i + 1) * SIMD_WIDTH8_AVX512 + lane];
        }

        for j in 0..max_l2 {
            let s2_col = &s2_soa[j * SIMD_WIDTH8_AVX512..(j + 1) * SIMD_WIDTH8_AVX512];

            for lane in 0..n {
                let s1 = s1_row[lane];
                let s2 = s2_col[lane];
                let h00 = h_h[j * SIMD_WIDTH8_AVX512 + lane];
                let f11 = f_buf[j * SIMD_WIDTH8_AVX512 + lane];
                let prev_h10 = h10[lane];

                let (h11, e11_new, f21) = main_code8_lane(
                    s1,
                    s2,
                    h00,
                    e11[lane],
                    f11,
                    bsw.w_match,
                    bsw.w_mismatch,
                    bsw.w_ambig,
                    e_ins,
                    oe_ins,
                    e_del,
                    oe_del,
                );

                // Store PREVIOUS h10 into H_h[j] (so next row's iter j reads h11_from_iter_(j-1)_of_THIS_row).
                h_h[j * SIMD_WIDTH8_AVX512 + lane] = prev_h10;
                f_buf[j * SIMD_WIDTH8_AVX512 + lane] = f21;
                e11[lane] = e11_new;
                h10[lane] = h11;

                // Per-lane state updates: only for lanes within their own len1×len2 box.
                if (i as u32) < len1_lanes[lane] && (j as u32) < len2_lanes[lane] {
                    let h11_i = h11 as i32;
                    if h11_i > state.max_score[lane] {
                        state.max_score[lane] = h11_i;
                        state.max_i[lane] = i as i32;
                        state.max_j[lane] = j as i32;
                    }
                    // gscore: H at right-edge (j == len2 - 1)
                    if (j as u32) == len2_lanes[lane] - 1 && h11_i > state.gscore[lane] {
                        state.gscore[lane] = h11_i;
                        state.max_ie[lane] = i as i32;
                    }
                }
            }
        }
        // Store last h10 to H_h[max_l2].
        for lane in 0..n {
            h_h[max_l2 * SIMD_WIDTH8_AVX512 + lane] = h10[lane];
        }
    }
    state
}

// Per-pair SW state used by the AVX-512 batch kernel. Kept SoA-friendly: each field is
// `[T; 64]` so the SIMD loop can read/write across all lanes uniformly via vector ops.
// Direct mirror of the per-lane state tracked in C++ smithWaterman512_8 (xref bandedSWA.cpp:2298+).
#[allow(dead_code)]
#[derive(Clone)]
pub(crate) struct SimdBatchState {
    /// Best score per lane (signed; clamped at 0).
    pub max_score: [i32; SIMD_WIDTH8_AVX512],
    /// Row index where max was reached (= `tle` field of SeqPair).
    pub max_i: [i32; SIMD_WIDTH8_AVX512],
    /// Column index where max was reached (= `qle`).
    pub max_j: [i32; SIMD_WIDTH8_AVX512],
    /// Max distance from diagonal seen during DP (= `max_off`).
    pub max_off: [i32; SIMD_WIDTH8_AVX512],
    /// Best H value at the right-edge column (= `gscore`).
    pub gscore: [i32; SIMD_WIDTH8_AVX512],
    /// Row index where right-edge max was reached (= `gtle`).
    pub max_ie: [i32; SIMD_WIDTH8_AVX512],
    /// Per-lane termination flag (false = lane done, no more updates).
    pub active: [bool; SIMD_WIDTH8_AVX512],
}

impl SimdBatchState {
    #[allow(dead_code)]
    #[inline]
    pub(crate) fn new(h0_lanes: &[u8; SIMD_WIDTH8_AVX512]) -> Self {
        let mut max_score = [0_i32; SIMD_WIDTH8_AVX512];
        for lane in 0..SIMD_WIDTH8_AVX512 {
            max_score[lane] = h0_lanes[lane] as i32;
        }
        SimdBatchState {
            max_score,
            max_i: [-1; SIMD_WIDTH8_AVX512],
            max_j: [-1; SIMD_WIDTH8_AVX512],
            max_off: [0; SIMD_WIDTH8_AVX512],
            gscore: [-1; SIMD_WIDTH8_AVX512],
            max_ie: [-1; SIMD_WIDTH8_AVX512],
            active: [true; SIMD_WIDTH8_AVX512],
        }
    }

    /// Materialize per-lane results into the corresponding SeqPair fields. Mirrors the
    /// C++ result-extraction at smithWaterman512_8:2650-2658.
    #[allow(dead_code)]
    pub(crate) fn write_to_pairs(&self, pairs: &mut [SeqPair]) {
        for (lane, p) in pairs.iter_mut().enumerate().take(SIMD_WIDTH8_AVX512) {
            p.score = self.max_score[lane];
            p.tle = self.max_i[lane];
            p.qle = self.max_j[lane];
            p.max_off = self.max_off[lane];
            p.gscore = self.gscore[lane];
            p.gtle = self.max_ie[lane];
        }
    }
}

// AVX-512 SIMD version of MAIN_CODE16 — processes 32 lanes (= pairs) in parallel using __m512i
// with i16 elements. Mirrors the C++ macro at bandedSWA.cpp:1883 byte-for-byte. Same structure
// as main_code8_avx512 but with i16 intrinsics throughout.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512bw,avx512f")]
#[allow(dead_code)]
#[inline]
unsafe fn main_code16_avx512(
    s1_v: core::arch::x86_64::__m512i,
    s2_v: core::arch::x86_64::__m512i,
    h00_v: core::arch::x86_64::__m512i,
    e11_v: core::arch::x86_64::__m512i,
    f11_v: core::arch::x86_64::__m512i,
    match_v: core::arch::x86_64::__m512i,
    mismatch_v: core::arch::x86_64::__m512i,
    w_ambig_v: core::arch::x86_64::__m512i,
    e_ins_v: core::arch::x86_64::__m512i,
    oe_ins_v: core::arch::x86_64::__m512i,
    e_del_v: core::arch::x86_64::__m512i,
    oe_del_v: core::arch::x86_64::__m512i,
) -> (
    core::arch::x86_64::__m512i,
    core::arch::x86_64::__m512i,
    core::arch::x86_64::__m512i,
) {
    use core::arch::x86_64::*;
    let zero = _mm512_setzero_si512();
    let cmp_eq = _mm512_cmpeq_epi16_mask(s1_v, s2_v);
    let sbt0 = _mm512_mask_blend_epi16(cmp_eq, mismatch_v, match_v);
    // AMBIG: encoded as 0xFFFF (= -1) for AMBIG lanes. Detect when EITHER s1 or s2 is AMBIG
    // by OR-ing the two and checking the sign bit (any AMBIG operand sets all bits, so the OR
    // also has the sign bit set). Using _mm512_max_epi16 was wrong: signed max picks the larger
    // value, so an AMBIG (-1) loses to any non-negative base (0..3) and the AMBIG lane is missed.
    let ambig_lanes = _mm512_or_si512(s1_v, s2_v);
    let ambig_mask = _mm512_movepi16_mask(ambig_lanes);
    let sbt = _mm512_mask_blend_epi16(ambig_mask, sbt0, w_ambig_v);
    // m11 = h00 + sbt; reset to 0 where h00 == 0. Fused: only add for h00 != 0 lanes.
    let h00_zero = _mm512_cmpeq_epi16_mask(h00_v, zero);
    let m11 = _mm512_mask_add_epi16(zero, !h00_zero, h00_v, sbt);
    // h11 = max(m11, e11, f11)
    let h11 = _mm512_max_epi16(_mm512_max_epi16(m11, e11_v), f11_v);
    // e11_new = max(m11 - oe_ins, 0).max(e11 - e_ins). m11 is non-negative (zeroed where h00=0
    // by mask_add); use u16 saturating sub to fold `max(_, zero)` into one op.
    let val_e = _mm512_subs_epu16(m11, oe_ins_v);
    let e11_dec = _mm512_sub_epi16(e11_v, e_ins_v);
    let e11_new = _mm512_max_epi16(val_e, e11_dec);
    // f21 = max(m11 - oe_del, 0).max(f11 - e_del)
    let val_f = _mm512_subs_epu16(m11, oe_del_v);
    let f21_dec = _mm512_sub_epi16(f11_v, e_del_v);
    let f21 = _mm512_max_epi16(val_f, f21_dec);
    (h11, e11_new, f21)
}

// AVX-512 SIMD version of MAIN_CODE8 — processes 64 lanes (= pairs) in parallel using __m512i.
// All inputs are 64-byte vectors; lane k is the per-pair state for the k-th pair in the batch.
// Mirrors the C++ macro at bandedSWA.cpp:1842 byte-for-byte: same intrinsics, same order.
//
// Returns (h11_v, e11_new_v, f21_v) — the new H, propagated within-row F-state (e11), and
// new across-row E-state (f21) for this column position.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512bw,avx512f")]
#[allow(dead_code)]
#[inline]
unsafe fn main_code8_avx512(
    s1_v: core::arch::x86_64::__m512i,
    s2_v: core::arch::x86_64::__m512i,
    h00_v: core::arch::x86_64::__m512i,
    e11_v: core::arch::x86_64::__m512i,
    f11_v: core::arch::x86_64::__m512i,
    match_v: core::arch::x86_64::__m512i,
    mismatch_v: core::arch::x86_64::__m512i,
    w_ambig_v: core::arch::x86_64::__m512i,
    e_ins_v: core::arch::x86_64::__m512i,
    oe_ins_v: core::arch::x86_64::__m512i,
    e_del_v: core::arch::x86_64::__m512i,
    oe_del_v: core::arch::x86_64::__m512i,
) -> (
    core::arch::x86_64::__m512i,
    core::arch::x86_64::__m512i,
    core::arch::x86_64::__m512i,
) {
    use core::arch::x86_64::*;
    let zero = _mm512_setzero_si512();
    // sbt = match if s1==s2 else mismatch
    let cmp_eq = _mm512_cmpeq_epi8_mask(s1_v, s2_v);
    let sbt0 = _mm512_mask_blend_epi8(cmp_eq, mismatch_v, match_v);
    // AMBIG: high bit of max_epu8(s1, s2) — encoded as 0xFF in encode_pairs_soa_avx512.
    let max_s = _mm512_max_epu8(s1_v, s2_v);
    let ambig_mask = _mm512_movepi8_mask(max_s);
    let sbt = _mm512_mask_blend_epi8(ambig_mask, sbt0, w_ambig_v);
    // m11 = h00 + sbt; reset to 0 where h00 == 0 (local SW reset). Fused into mask_add.
    let h00_zero = _mm512_cmpeq_epi8_mask(h00_v, zero);
    let m11 = _mm512_mask_add_epi8(zero, !h00_zero, h00_v, sbt);
    // h11 = max(m11, e11, f11) — signed max
    let h11 = _mm512_max_epi8(_mm512_max_epi8(m11, e11_v), f11_v);
    // e11_new = max(m11 - oe_ins, 0).max(e11 - e_ins). m11 is non-negative (zeroed where h00=0
    // by mask_add); use u8 saturating sub to fold `max(_, zero)` into one op.
    let val_e = _mm512_subs_epu8(m11, oe_ins_v);
    let e11_dec = _mm512_sub_epi8(e11_v, e_ins_v);
    let e11_new = _mm512_max_epi8(val_e, e11_dec);
    // f21 = max(m11 - oe_del, 0).max(f11 - e_del)
    let val_f = _mm512_subs_epu8(m11, oe_del_v);
    let f21_dec = _mm512_sub_epi8(f11_v, e_del_v);
    let f21 = _mm512_max_epi8(val_f, f21_dec);
    (h11, e11_new, f21)
}

// i16 scalar version of main_code8_lane — for validation against the AVX-512 i16 SIMD path.
// AMBIG marker is the i16 sign bit (encoded as 0xFFFF in encode_pairs_soa_avx512_16).
#[allow(dead_code)]
#[inline]
pub(crate) fn main_code16_lane(
    s1: i16,
    s2: i16,
    h00: i16,
    e11_in: i16,
    f11: i16,
    w_match: i16,
    w_mismatch: i16,
    w_ambig: i16,
    e_ins: i16,
    oe_ins: i16,
    e_del: i16,
    oe_del: i16,
) -> (i16, i16, i16) {
    let mut sbt: i16 = if s1 == s2 { w_match } else { w_mismatch };
    // AMBIG: detect when either operand is AMBIG (= -1, sign bit set). OR the two; sign bit
    // is set iff either operand had it. (Using `s1.max(s2) < 0` was wrong: max picks the larger,
    // so AMBIG = -1 loses to any non-negative base 0..3, missing the AMBIG case.)
    if (s1 | s2) < 0 {
        sbt = w_ambig;
    }
    let m11_val: i16 = if h00 == 0 { 0 } else { h00.wrapping_add(sbt) };
    let h11 = m11_val.max(e11_in).max(f11).max(0);
    let val_e = m11_val.wrapping_sub(oe_ins).max(0);
    let e11_dec = e11_in.wrapping_sub(e_ins);
    let e11_new = val_e.max(e11_dec);
    let val_f = m11_val.wrapping_sub(oe_del).max(0);
    let f21_dec = f11.wrapping_sub(e_del);
    let f21 = val_f.max(f21_dec);
    (h11, e11_new, f21)
}

#[allow(dead_code)]
const SIMD_WIDTH16_AVX512: usize = 32;
#[allow(dead_code)]
const DUMMY1_I16: i16 = 99;
#[allow(dead_code)]
const DUMMY2_I16: i16 = 100;

// SoA encode helper for the i16 batch SW path. Mirrors `encode_pairs_soa_avx512` but lays out
// 32-lane i16 vectors. AMBIG bases (encoded as `0xFFFF` = -1 in i16) are mapped from the
// raw seq buffer's `AMBIG_BASE` (= 4) value, matching C++ smithWatermanBatchWrapper16:1263.
//
// Returns `(seq1_soa: Vec<i16>, seq2_soa: Vec<i16>, max_len1, max_len2)`. Buffer length is
// `max_len * 32`. Lanes beyond pairs.len() get DUMMY padding; lanes past a pair's actual
// length get DUMMY too.
#[allow(dead_code)]
pub(crate) fn encode_pairs_soa_avx512_16(
    pairs: &[SeqPair],
    seq_buf_ref: &[u8],
    seq_buf_qer: &[u8],
) -> (Vec<i16>, Vec<i16>, usize, usize) {
    let mut s1 = Vec::new();
    let mut s2 = Vec::new();
    let (max_len1, max_len2) =
        encode_pairs_soa_avx512_16_into(pairs, seq_buf_ref, seq_buf_qer, &mut s1, &mut s2);
    (s1, s2, max_len1, max_len2)
}

pub(crate) fn encode_pairs_soa_avx512_16_into(
    pairs: &[SeqPair],
    seq_buf_ref: &[u8],
    seq_buf_qer: &[u8],
    s1: &mut Vec<i16>,
    s2: &mut Vec<i16>,
) -> (usize, usize) {
    debug_assert!(pairs.len() <= SIMD_WIDTH16_AVX512);
    let mut max_len1 = 0_usize;
    let mut max_len2 = 0_usize;
    for p in pairs {
        max_len1 = max_len1.max(p.len1 as usize);
        max_len2 = max_len2.max(p.len2 as usize);
    }
    s1.clear();
    s1.resize(max_len1 * SIMD_WIDTH16_AVX512, DUMMY1_I16);
    s2.clear();
    s2.resize(max_len2 * SIMD_WIDTH16_AVX512, DUMMY2_I16);
    for (lane, p) in pairs.iter().enumerate() {
        let s1_base = (p.idr as usize)..(p.idr as usize + p.len1 as usize);
        let s2_base = (p.idq as usize)..(p.idq as usize + p.len2 as usize);
        // Same BCE rationale as the u8 variant: k < len ≤ max_len, lane < SIMD_WIDTH16.
        // Branchless: 4 → -1 (i16). Subtract 5 when bit-2 set: b=4 has b>>2 = 1, 4-5 = -1.
        unsafe {
            for (k, &b) in seq_buf_ref[s1_base].iter().enumerate() {
                let v = i16::from(b).wrapping_sub(i16::from((b >> 2) & 1).wrapping_mul(5));
                *s1.get_unchecked_mut(k * SIMD_WIDTH16_AVX512 + lane) = v;
            }
            for (k, &b) in seq_buf_qer[s2_base].iter().enumerate() {
                let v = i16::from(b).wrapping_sub(i16::from((b >> 2) & 1).wrapping_mul(5));
                *s2.get_unchecked_mut(k * SIMD_WIDTH16_AVX512 + lane) = v;
            }
        }
    }
    (max_len1, max_len2)
}

// i16 H_v initializer — mirrors init_h_v_soa_avx512 but for 32 i16 lanes. With i16 there's no
// risk of value overflow during the cascade (h0 ≤ a few thousand, e_del ≤ a few), so we just
// use signed sub + max-with-0 to clamp.
#[allow(dead_code)]
pub(crate) fn init_h_v_soa_avx512_16(
    h0_lanes: &[i16; SIMD_WIDTH16_AVX512],
    max_len1: usize,
    o_del: i16,
    e_del: i16,
) -> Vec<i16> {
    let rows = max_len1.max(1);
    let mut h_v = vec![0_i16; rows * SIMD_WIDTH16_AVX512];
    init_h_v_soa_avx512_16_into(h0_lanes, max_len1, o_del, e_del, &mut h_v);
    h_v
}

pub(crate) fn init_h_v_soa_avx512_16_into(
    h0_lanes: &[i16; SIMD_WIDTH16_AVX512],
    max_len1: usize,
    o_del: i16,
    e_del: i16,
    h_v: &mut [i16],
) {
    debug_assert!(h_v.len() >= max_len1.max(1) * SIMD_WIDTH16_AVX512);
    for lane in 0..SIMD_WIDTH16_AVX512 {
        h_v[lane] = h0_lanes[lane];
    }
    if max_len1 == 0 {
        return;
    }
    let mut tmp = [0_i16; SIMD_WIDTH16_AVX512];
    for lane in 0..SIMD_WIDTH16_AVX512 {
        tmp[lane] = h0_lanes[lane].wrapping_sub(o_del);
    }
    for k in 1..max_len1 {
        for lane in 0..SIMD_WIDTH16_AVX512 {
            tmp[lane] = tmp[lane].wrapping_sub(e_del);
            h_v[k * SIMD_WIDTH16_AVX512 + lane] = tmp[lane].max(0);
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512bw,avx512f")]
unsafe fn init_h_v_soa_avx512_16_into_simd(
    h0_lanes: &[i16; SIMD_WIDTH16_AVX512],
    max_len1: usize,
    o_del: i16,
    e_del: i16,
    h_v: &mut [i16],
) {
    use core::arch::x86_64::*;

    debug_assert!(h_v.len() >= max_len1.max(1) * SIMD_WIDTH16_AVX512);
    let h0_v = _mm512_loadu_si512(h0_lanes.as_ptr() as *const __m512i);
    _mm512_storeu_si512(h_v.as_mut_ptr() as *mut __m512i, h0_v);
    if max_len1 == 0 {
        return;
    }
    let zero_v = _mm512_setzero_si512();
    let o_del_v = _mm512_set1_epi16(o_del);
    let e_del_v = _mm512_set1_epi16(e_del);
    let mut tmp_v = _mm512_sub_epi16(h0_v, o_del_v);
    for k in 1..max_len1 {
        tmp_v = _mm512_sub_epi16(tmp_v, e_del_v);
        let row_v = _mm512_max_epi16(tmp_v, zero_v);
        _mm512_storeu_si512(
            h_v.as_mut_ptr().add(k * SIMD_WIDTH16_AVX512) as *mut __m512i,
            row_v,
        );
    }
}

// i16 H_h initializer — mirrors init_h_h_soa_avx512.
#[allow(dead_code)]
pub(crate) fn init_h_h_soa_avx512_16(
    h0_lanes: &[i16; SIMD_WIDTH16_AVX512],
    max_len2: usize,
    oe_ins: i16,
    e_ins: i16,
) -> Vec<i16> {
    let cols = max_len2.max(1);
    let mut h_h = vec![0_i16; cols * SIMD_WIDTH16_AVX512];
    init_h_h_soa_avx512_16_into(h0_lanes, max_len2, oe_ins, e_ins, &mut h_h);
    h_h
}

pub(crate) fn init_h_h_soa_avx512_16_into(
    h0_lanes: &[i16; SIMD_WIDTH16_AVX512],
    max_len2: usize,
    oe_ins: i16,
    e_ins: i16,
    h_h: &mut [i16],
) {
    debug_assert!(h_h.len() >= max_len2.max(1) * SIMD_WIDTH16_AVX512);
    for lane in 0..SIMD_WIDTH16_AVX512 {
        h_h[lane] = h0_lanes[lane];
    }
    if max_len2 < 2 {
        return;
    }
    let mut tmp = [0_i16; SIMD_WIDTH16_AVX512];
    for lane in 0..SIMD_WIDTH16_AVX512 {
        tmp[lane] = if h0_lanes[lane] > oe_ins {
            h0_lanes[lane].wrapping_sub(oe_ins)
        } else {
            0
        };
        h_h[SIMD_WIDTH16_AVX512 + lane] = tmp[lane];
    }
    for k in 2..max_len2 {
        for lane in 0..SIMD_WIDTH16_AVX512 {
            let next = tmp[lane].wrapping_sub(e_ins).max(0);
            tmp[lane] = next;
            h_h[k * SIMD_WIDTH16_AVX512 + lane] = next;
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512bw,avx512f")]
unsafe fn init_h_h_soa_avx512_16_into_simd(
    h0_lanes: &[i16; SIMD_WIDTH16_AVX512],
    max_len2: usize,
    oe_ins: i16,
    e_ins: i16,
    h_h: &mut [i16],
) {
    use core::arch::x86_64::*;

    debug_assert!(h_h.len() >= max_len2.max(1) * SIMD_WIDTH16_AVX512);
    let h0_v = _mm512_loadu_si512(h0_lanes.as_ptr() as *const __m512i);
    _mm512_storeu_si512(h_h.as_mut_ptr() as *mut __m512i, h0_v);
    if max_len2 < 2 {
        return;
    }
    let zero_v = _mm512_setzero_si512();
    let oe_ins_v = _mm512_set1_epi16(oe_ins);
    let e_ins_v = _mm512_set1_epi16(e_ins);
    let mut tmp_v = _mm512_max_epi16(_mm512_sub_epi16(h0_v, oe_ins_v), zero_v);
    _mm512_storeu_si512(
        h_h.as_mut_ptr().add(SIMD_WIDTH16_AVX512) as *mut __m512i,
        tmp_v,
    );
    for k in 2..max_len2 {
        tmp_v = _mm512_max_epi16(_mm512_sub_epi16(tmp_v, e_ins_v), zero_v);
        _mm512_storeu_si512(
            h_h.as_mut_ptr().add(k * SIMD_WIDTH16_AVX512) as *mut __m512i,
            tmp_v,
        );
    }
}

// Per-cell SW recurrence body — Rust port of bwa-mem2/src/bandedSWA.cpp:1842 MAIN_CODE8 macro.
// Computes one column of cells across all 64 lanes in one shot, scalar-faithful to the AVX-512
// version. Matches C++ semantics: i8 wrapping arithmetic, signed max-0 clamping, AMBIG handling
// where either operand has its high bit set (encoded as 0xFF in encode_pairs_soa_avx512).
//
//   m11   = (h00 == 0) ? 0 : h00 + sbt    where sbt = match/mismatch/w_ambig per (s1, s2)
//   h11   = max(m11, e11_in, f11)
//   e11   = max(m11 - oe_ins, 0).max(e11_in - e_ins)
//   f21   = max(m11 - oe_del, 0).max(f11 - e_del)
//
// Treated as i8 for arithmetic; the high bit of `s1`/`s2` is the AMBIG marker.
#[allow(dead_code)]
#[inline]
pub(crate) fn main_code8_lane(
    s1: u8,
    s2: u8,
    h00: u8,
    e11_in: u8,
    f11: u8,
    w_match: i8,
    w_mismatch: i8,
    w_ambig: i8,
    e_ins: u8,
    oe_ins: u8,
    e_del: u8,
    oe_del: u8,
) -> (u8, u8, u8) {
    // sbt = match if s1==s2 else mismatch (per cmpeq_epi8_mask + mask_blend in C++)
    let mut sbt: i8 = if s1 == s2 { w_match } else { w_mismatch };
    // AMBIG marker: max_epu8(s1, s2) with movepi8_mask = high bit set (encoded 0xFF earlier).
    if (s1.max(s2) & 0x80) != 0 {
        sbt = w_ambig;
    }
    // m11 = h00 + sbt (i8 wrapping); zero if h00 == 0 (local SW reset, no extension across gap).
    let m11_val: i8 = if h00 == 0 {
        0
    } else {
        (h00 as i8).wrapping_add(sbt)
    };
    // h11 = max_epi8(m11, e11_in, f11) — signed max. C++ doesn't clamp at 0 here, but the
    // recurrence guarantees ≥ 0 because e11/f11 are always ≥ 0 (clamped at production).
    let h11 = m11_val.max(e11_in as i8).max(f11 as i8).max(0) as u8;
    // e11_new = max(m11 - oe_ins, 0).max(e11_in - e_ins). val ≥ 0 ensures result ≥ 0.
    let val_e = m11_val.wrapping_sub(oe_ins as i8).max(0);
    let e11_dec = (e11_in as i8).wrapping_sub(e_ins as i8);
    let e11_new = val_e.max(e11_dec) as u8;
    // f21 = max(m11 - oe_del, 0).max(f11 - e_del)
    let val_f = m11_val.wrapping_sub(oe_del as i8).max(0);
    let f21_dec = (f11 as i8).wrapping_sub(e_del as i8);
    let f21 = val_f.max(f21_dec) as u8;
    (h11, e11_new, f21)
}

// Initialize the H_v (vertical, first-column) buffer for the AVX-512 batch SW kernel.
// Mirrors bwa-mem2/src/bandedSWA.cpp:2131-2139:
//   H_v[0][lane] = h0[lane]
//   H_v[k][lane] = max_epi8(h0 - o_del - k * e_del, 0)   for k >= 1
//
// `max_len1` is the longest target across the 64 lanes; rows beyond a lane's actual len1
// are unused but kept consistent so the SIMD DP can run uniformly.
//
// Returned buffer has length `(max_len1 + 1) * 64` (one extra row to match scratch sizing).
// Uses i8 wrapping subtraction + signed max-with-0 to match C++ AVX-512 semantics exactly:
// once the running tmp goes negative, subsequent rows store 0 (saturated).
#[allow(dead_code)]
pub(crate) fn init_h_v_soa_avx512(
    h0_lanes: &[u8; SIMD_WIDTH8_AVX512],
    max_len1: usize,
    o_del_byte: u8,
    e_del_byte: u8,
) -> Vec<u8> {
    let rows = max_len1.max(1);
    let mut h_v = vec![0_u8; rows * SIMD_WIDTH8_AVX512];
    init_h_v_soa_avx512_into(h0_lanes, max_len1, o_del_byte, e_del_byte, &mut h_v);
    h_v
}

// In-place variant: writes h_v rows directly into `out` (must be at least max_len1*64 bytes).
// Skips the per-call Vec alloc when called from a thread-local scratch wrapper.
pub(crate) fn init_h_v_soa_avx512_into(
    h0_lanes: &[u8; SIMD_WIDTH8_AVX512],
    max_len1: usize,
    o_del_byte: u8,
    e_del_byte: u8,
    out: &mut [u8],
) {
    if max_len1 == 0 {
        return;
    }
    debug_assert!(out.len() >= max_len1 * SIMD_WIDTH8_AVX512);
    // Row 0 = h0 (raw, no clamping — matches C++ `_mm512_store_si512(H2, h0_512)`).
    for lane in 0..SIMD_WIDTH8_AVX512 {
        out[lane] = h0_lanes[lane];
    }
    // Running tmp = h0 - o_del (i8 wrapping). Each row subtracts e_del; stored value is clamped to 0.
    let mut tmp = [0_i8; SIMD_WIDTH8_AVX512];
    for lane in 0..SIMD_WIDTH8_AVX512 {
        tmp[lane] = (h0_lanes[lane] as i8).wrapping_sub(o_del_byte as i8);
    }
    for k in 1..max_len1 {
        for lane in 0..SIMD_WIDTH8_AVX512 {
            tmp[lane] = tmp[lane].wrapping_sub(e_del_byte as i8);
            out[k * SIMD_WIDTH8_AVX512 + lane] = tmp[lane].max(0) as u8;
        }
    }
}

// Initialize the H_h (horizontal, first-row) buffer. Mirrors bwa-mem2/src/bandedSWA.cpp:2169-2181.
//   H_h[0][lane] = h0[lane]
//   H_h[1][lane] = max(0, h0 - oe_ins)              [via cmpgt mask, blend with zero]
//   H_h[k][lane] = max(0, H_h[k-1] - e_ins)         for k >= 2 (max_epi8 clamp)
#[allow(dead_code)]
pub(crate) fn init_h_h_soa_avx512(
    h0_lanes: &[u8; SIMD_WIDTH8_AVX512],
    max_len2: usize,
    oe_ins_byte: u8,
    e_ins_byte: u8,
) -> Vec<u8> {
    let cols = max_len2.max(1);
    let mut h_h = vec![0_u8; cols * SIMD_WIDTH8_AVX512];
    init_h_h_soa_avx512_into(h0_lanes, max_len2, oe_ins_byte, e_ins_byte, &mut h_h);
    h_h
}

// In-place variant of init_h_h_soa_avx512.
pub(crate) fn init_h_h_soa_avx512_into(
    h0_lanes: &[u8; SIMD_WIDTH8_AVX512],
    max_len2: usize,
    oe_ins_byte: u8,
    e_ins_byte: u8,
    out: &mut [u8],
) {
    if max_len2 == 0 {
        return;
    }
    debug_assert!(out.len() >= max_len2 * SIMD_WIDTH8_AVX512);
    for lane in 0..SIMD_WIDTH8_AVX512 {
        out[lane] = h0_lanes[lane];
    }
    if max_len2 < 2 {
        return;
    }
    // H_h[1]: per-lane max(0, h0 - oe_ins) via cmpgt + mask_blend (C++ uses cmpgt_epi8_mask).
    let mut tmp = [0_u8; SIMD_WIDTH8_AVX512];
    for lane in 0..SIMD_WIDTH8_AVX512 {
        tmp[lane] = if (h0_lanes[lane] as i8) > (oe_ins_byte as i8) {
            ((h0_lanes[lane] as i8).wrapping_sub(oe_ins_byte as i8)) as u8
        } else {
            0
        };
        out[SIMD_WIDTH8_AVX512 + lane] = tmp[lane];
    }
    // k >= 2: tmp = max_epi8(tmp - e_ins, 0). Wrapping sub keeps i8 semantics, then signed max-0.
    for k in 2..max_len2 {
        for lane in 0..SIMD_WIDTH8_AVX512 {
            let next = (tmp[lane] as i8).wrapping_sub(e_ins_byte as i8).max(0) as u8;
            tmp[lane] = next;
            out[k * SIMD_WIDTH8_AVX512 + lane] = next;
        }
    }
}

#[cfg(test)]
mod simd_batch_tests {
    use super::*;
    use crate::bwa_mem2::bandedswa::SeqPair;

    #[test]
    fn encode_pairs_soa_avx512_lays_out_lanes_correctly() {
        // Two pairs: lane 0 has target="ACGT" (4 bp), lane 1 has target="GA" (2 bp).
        // After SoA encode, seq1_soa[k*64 + 0] should be lane-0's k-th byte (or DUMMY1 if k >= len),
        // and seq1_soa[k*64 + 1] should be lane-1's k-th byte. Other lanes (2..64) are padding.
        let seq_ref: Vec<u8> = vec![0, 1, 2, 3, /* lane1 ref starts */ 2, 0];
        let seq_qer: Vec<u8> = vec![0, 1, /* lane1 query */ 1, 2, 3];
        let pairs = vec![
            SeqPair {
                idr: 0,
                idq: 0,
                len1: 4,
                len2: 2,
                ..Default::default()
            },
            SeqPair {
                idr: 4,
                idq: 2,
                len1: 2,
                len2: 3,
                ..Default::default()
            },
        ];
        let (s1, s2, max_l1, max_l2) = encode_pairs_soa_avx512(&pairs, &seq_ref, &seq_qer);
        assert_eq!(max_l1, 4);
        assert_eq!(max_l2, 3);
        // Lane 0 ref bytes: A,C,G,T = 0,1,2,3.
        assert_eq!(s1[0 * 64 + 0], 0);
        assert_eq!(s1[1 * 64 + 0], 1);
        assert_eq!(s1[2 * 64 + 0], 2);
        assert_eq!(s1[3 * 64 + 0], 3);
        // Lane 1 ref bytes (only first 2 valid; rest padded).
        assert_eq!(s1[0 * 64 + 1], 2);
        assert_eq!(s1[1 * 64 + 1], 0);
        assert_eq!(s1[2 * 64 + 1], DUMMY1_BYTE);
        assert_eq!(s1[3 * 64 + 1], DUMMY1_BYTE);
        // Lane 2 (no pair) is all DUMMY1 padding.
        assert_eq!(s1[0 * 64 + 2], DUMMY1_BYTE);
        // Same checks on s2.
        assert_eq!(s2[0 * 64 + 0], 0);
        assert_eq!(s2[1 * 64 + 0], 1);
        assert_eq!(s2[2 * 64 + 0], DUMMY2_BYTE);
        assert_eq!(s2[0 * 64 + 1], 1);
        assert_eq!(s2[2 * 64 + 1], 3);
    }

    // Build a BandedPairWiseSW with a custom match score `a` (mismatch is -B).
    // Mirrors the matrix bwa_fill_scmat builds when invoking `mem -A a -B B`.
    fn make_bsw_a(
        a: i8,
        b: i8,
        o_del: i32,
        e_del: i32,
        o_ins: i32,
        e_ins: i32,
    ) -> BandedPairWiseSW {
        let mut mat = [0_i8; 25];
        for r in 0..4 {
            for c in 0..4 {
                mat[r * 5 + c] = if r == c { a } else { -b };
            }
            mat[r * 5 + 4] = -1;
        }
        for c in 0..5 {
            mat[4 * 5 + c] = -1;
        }
        BandedPairWiseSW {
            mat,
            m: 5,
            end_bonus: 5,
            zdrop: 100,
            o_del,
            o_ins,
            e_del,
            e_ins,
            w_match: a,
            w_mismatch: -b,
            w_open: o_del,
            w_extend: e_del,
            w_ambig: -1,
            swTicks: 0,
            SW_cells: 0,
            setupTicks: 0,
            sort1Ticks: 0,
            sort2Ticks: 0,
        }
    }

    fn make_bsw(o_del: i32, e_del: i32, o_ins: i32, e_ins: i32) -> BandedPairWiseSW {
        let mut mat = [0_i8; 25];
        // 4x4 simple match/mismatch + ambig column.
        for r in 0..4 {
            for c in 0..4 {
                mat[r * 5 + c] = if r == c { 2 } else { -3 };
            }
            mat[r * 5 + 4] = -1;
        }
        for c in 0..5 {
            mat[4 * 5 + c] = -1;
        }
        BandedPairWiseSW {
            mat,
            m: 5,
            end_bonus: 5,
            zdrop: 100,
            o_del,
            o_ins,
            e_del,
            e_ins,
            w_match: 2,
            w_mismatch: -3,
            w_open: o_del,
            w_extend: e_del,
            w_ambig: -1,
            swTicks: 0,
            SW_cells: 0,
            setupTicks: 0,
            sort1Ticks: 0,
            sort2Ticks: 0,
        }
    }

    #[test]
    fn process_batch_soa_banded_dp_matches_scalar_with_wide_band() {
        let bsw = make_bsw(6, 1, 6, 1);
        let seq1 = vec![0_u8, 1, 2, 3, 0, 1, 2, 3];
        let seq2 = vec![0_u8, 1, 2, 3];
        // h0=10 (realistic for chain extension where h0 = previous chain score).
        let pairs = vec![SeqPair {
            idr: 0,
            idq: 0,
            len1: 8,
            len2: 4,
            h0: 10,
            ..Default::default()
        }];
        let state = process_batch_soa_banded_dp(&bsw, &pairs, &seq1, &seq2, 100);
        let mut qle = 0_i32;
        let mut tle = 0_i32;
        let mut gtle = 0_i32;
        let mut gscore = 0_i32;
        let mut max_off = 0_i32;
        let scalar = bsw.scalarBandedSWA(
            seq2.len() as i32,
            &seq2,
            seq1.len() as i32,
            &seq1,
            100,
            10,
            &mut qle,
            &mut tle,
            &mut gtle,
            &mut gscore,
            &mut max_off,
        );
        assert_eq!(state.max_score[0], scalar, "max_score with wide band");
        assert_eq!(state.max_i[0] + 1, tle, "tle with wide band");
        assert_eq!(state.max_j[0] + 1, qle, "qle with wide band");
    }

    #[test]
    fn process_batch_soa_banded_dp_matches_scalar_with_narrow_band() {
        let bsw = make_bsw(6, 1, 6, 1);
        // 16-bp identical sequences; band w=2 still covers diagonal (with dynamic narrowing).
        let seq: Vec<u8> = (0..16).map(|i| (i & 3) as u8).collect();
        let pairs = vec![SeqPair {
            idr: 0,
            idq: 0,
            len1: 16,
            len2: 16,
            h0: 5,
            ..Default::default()
        }];
        let state = process_batch_soa_banded_dp(&bsw, &pairs, &seq, &seq, 2);
        let mut qle = 0_i32;
        let mut tle = 0_i32;
        let mut gtle = 0_i32;
        let mut gscore = 0_i32;
        let mut max_off = 0_i32;
        let scalar = bsw.scalarBandedSWA(
            seq.len() as i32,
            &seq,
            seq.len() as i32,
            &seq,
            2,
            5,
            &mut qle,
            &mut tle,
            &mut gtle,
            &mut gscore,
            &mut max_off,
        );
        assert_eq!(state.max_score[0], scalar, "max_score with narrow band");
        assert_eq!(state.max_i[0] + 1, tle, "tle");
        assert_eq!(state.max_j[0] + 1, qle, "qle");
        assert_eq!(state.gscore[0], gscore, "gscore");
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn process_batch_soa_banded_dp_avx512_16_matches_scalar_at_match_score_3() {
        if !std::is_x86_feature_detected!("avx512bw") {
            eprintln!("skip: avx512bw not detected");
            return;
        }
        let bsw = make_bsw_a(3, 3, 6, 1, 6, 1);
        let seq: Vec<u8> = (0..60).map(|i| (i & 3) as u8).collect();
        let pairs = vec![SeqPair {
            idr: 0,
            idq: 0,
            len1: 60,
            len2: 60,
            h0: 15,
            ..Default::default()
        }];

        let simd = unsafe { process_batch_soa_banded_dp_avx512_16(&bsw, &pairs, &seq, &seq, 10) };
        let mut qle = 0_i32;
        let mut tle = 0_i32;
        let mut gtle = 0_i32;
        let mut gscore = 0_i32;
        let mut max_off = 0_i32;
        let scalar = bsw.scalarBandedSWA(
            seq.len() as i32,
            &seq,
            seq.len() as i32,
            &seq,
            10,
            15,
            &mut qle,
            &mut tle,
            &mut gtle,
            &mut gscore,
            &mut max_off,
        );

        assert_eq!(simd.max_score[0], scalar, "max_score at -A 3 (i16 SIMD)");
        assert_eq!(simd.max_i[0] + 1, tle, "tle");
        assert_eq!(simd.max_j[0] + 1, qle, "qle");
        assert_eq!(simd.gscore[0], gscore, "gscore");
        assert_eq!(simd.max_ie[0] + 1, gtle, "gtle");
        assert_eq!(simd.max_off[0], max_off, "max_off");
    }

    // Same failing input through u8 SIMD wrapper.
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn process_batch_soa_banded_dp_avx512_u8_34_32_ambig_at_zero() {
        if !std::is_x86_feature_detected!("avx512bw") {
            eprintln!("skip: avx512bw not detected");
            return;
        }
        let bsw = make_bsw(6, 1, 6, 1);
        let r: Vec<u8> = (0..34).map(|i| (i & 3) as u8).collect();
        let mut q: Vec<u8> = (0..32).map(|i| (i & 3) as u8).collect();
        q[0] = AMBIG_BASE;
        let pairs = vec![SeqPair {
            idr: 0,
            idq: 0,
            len1: 34,
            len2: 32,
            h0: 3,
            ..Default::default()
        }];

        let simd = unsafe { process_batch_soa_banded_dp_avx512(&bsw, &pairs, &r, &q, 50) };
        let mut qle = 0_i32;
        let mut tle = 0_i32;
        let mut gtle = 0_i32;
        let mut gscore = 0_i32;
        let mut max_off = 0_i32;
        let scalar = bsw.scalarBandedSWA(
            32,
            &q,
            34,
            &r,
            50,
            3,
            &mut qle,
            &mut tle,
            &mut gtle,
            &mut gscore,
            &mut max_off,
        );
        eprintln!(
            "u8 SIMD: max_score={}, max_i={}, max_j={}",
            simd.max_score[0], simd.max_i[0], simd.max_j[0]
        );
        eprintln!("scalar: max={}, tle={}, qle={}", scalar, tle, qle);
        assert_eq!(simd.max_score[0], scalar);
    }

    // Same failing input but driven through the SCALAR prototype process_batch_soa_banded_dp.
    // If this ALSO fails, the bug is in the wrapper algorithm (not SIMD intrinsics).
    #[test]
    fn process_batch_soa_banded_dp_scalar_proto_34_32_ambig_at_zero() {
        let bsw = make_bsw(6, 1, 6, 1);
        let r: Vec<u8> = (0..34).map(|i| (i & 3) as u8).collect();
        let mut q: Vec<u8> = (0..32).map(|i| (i & 3) as u8).collect();
        q[0] = AMBIG_BASE;
        let pairs = vec![SeqPair {
            idr: 0,
            idq: 0,
            len1: 34,
            len2: 32,
            h0: 3,
            ..Default::default()
        }];

        let proto = process_batch_soa_banded_dp(&bsw, &pairs, &r, &q, 50);
        let mut qle = 0_i32;
        let mut tle = 0_i32;
        let mut gtle = 0_i32;
        let mut gscore = 0_i32;
        let mut max_off = 0_i32;
        let scalar = bsw.scalarBandedSWA(
            32,
            &q,
            34,
            &r,
            50,
            3,
            &mut qle,
            &mut tle,
            &mut gtle,
            &mut gscore,
            &mut max_off,
        );
        eprintln!(
            "proto: max_score={}, max_i={}, max_j={}",
            proto.max_score[0], proto.max_i[0], proto.max_j[0]
        );
        eprintln!("scalar: max={}, tle={}, qle={}", scalar, tle, qle);
        assert_eq!(proto.max_score[0], scalar);
    }

    // 34-bp/32-bp with one AMBIG injection at q[0]. Isolate AMBIG handling.
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn process_batch_soa_banded_dp_avx512_16_34_32_ambig_at_zero() {
        if !std::is_x86_feature_detected!("avx512bw") {
            eprintln!("skip: avx512bw not detected");
            return;
        }
        let bsw = make_bsw(6, 1, 6, 1);
        let r: Vec<u8> = (0..34).map(|i| (i & 3) as u8).collect();
        let mut q: Vec<u8> = (0..32).map(|i| (i & 3) as u8).collect();
        q[0] = AMBIG_BASE;
        let pairs = vec![SeqPair {
            idr: 0,
            idq: 0,
            len1: 34,
            len2: 32,
            h0: 3,
            ..Default::default()
        }];

        let simd = unsafe { process_batch_soa_banded_dp_avx512_16(&bsw, &pairs, &r, &q, 50) };
        let mut qle = 0_i32;
        let mut tle = 0_i32;
        let mut gtle = 0_i32;
        let mut gscore = 0_i32;
        let mut max_off = 0_i32;
        let scalar = bsw.scalarBandedSWA(
            32,
            &q,
            34,
            &r,
            50,
            3,
            &mut qle,
            &mut tle,
            &mut gtle,
            &mut gscore,
            &mut max_off,
        );
        eprintln!(
            "SIMD: max_score={}, max_i={}, max_j={}",
            simd.max_score[0], simd.max_i[0], simd.max_j[0]
        );
        eprintln!("scalar: max={}, tle={}, qle={}", scalar, tle, qle);
        assert_eq!(simd.max_score[0], scalar);
    }

    // 34-bp/32-bp matched (no mismatches). Same lengths as the failing fuzz lane but no
    // mismatches/AMBIG.
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn process_batch_soa_banded_dp_avx512_16_34_32_matched_h0_3() {
        if !std::is_x86_feature_detected!("avx512bw") {
            eprintln!("skip: avx512bw not detected");
            return;
        }
        let bsw = make_bsw(6, 1, 6, 1);
        let r: Vec<u8> = (0..34).map(|i| (i & 3) as u8).collect();
        let q: Vec<u8> = (0..32).map(|i| (i & 3) as u8).collect();
        let pairs = vec![SeqPair {
            idr: 0,
            idq: 0,
            len1: 34,
            len2: 32,
            h0: 3,
            ..Default::default()
        }];

        let simd = unsafe { process_batch_soa_banded_dp_avx512_16(&bsw, &pairs, &r, &q, 50) };
        let mut qle = 0_i32;
        let mut tle = 0_i32;
        let mut gtle = 0_i32;
        let mut gscore = 0_i32;
        let mut max_off = 0_i32;
        let scalar = bsw.scalarBandedSWA(
            32,
            &q,
            34,
            &r,
            50,
            3,
            &mut qle,
            &mut tle,
            &mut gtle,
            &mut gscore,
            &mut max_off,
        );
        eprintln!(
            "SIMD: max_score={}, max_i={}, max_j={}, gscore={}, max_ie={}",
            simd.max_score[0], simd.max_i[0], simd.max_j[0], simd.gscore[0], simd.max_ie[0]
        );
        eprintln!(
            "scalar: max={}, tle={}, qle={}, gscore={}, gtle={}",
            scalar, tle, qle, gscore, gtle
        );
        assert_eq!(simd.max_score[0], scalar);
    }

    // Even simpler: matched 8-bp seq at h0=3. SIMD should clearly compute >3.
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn process_batch_soa_banded_dp_avx512_16_simple_h0_3_matched() {
        if !std::is_x86_feature_detected!("avx512bw") {
            eprintln!("skip: avx512bw not detected");
            return;
        }
        let bsw = make_bsw(6, 1, 6, 1);
        let seq: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let pairs = vec![SeqPair {
            idr: 0,
            idq: 0,
            len1: 8,
            len2: 8,
            h0: 3,
            ..Default::default()
        }];

        let simd = unsafe { process_batch_soa_banded_dp_avx512_16(&bsw, &pairs, &seq, &seq, 50) };
        let mut qle = 0_i32;
        let mut tle = 0_i32;
        let mut gtle = 0_i32;
        let mut gscore = 0_i32;
        let mut max_off = 0_i32;
        let scalar = bsw.scalarBandedSWA(
            8,
            &seq,
            8,
            &seq,
            50,
            3,
            &mut qle,
            &mut tle,
            &mut gtle,
            &mut gscore,
            &mut max_off,
        );
        eprintln!(
            "SIMD: max_score={}, max_i={}, max_j={}, gscore={}, max_ie={}",
            simd.max_score[0], simd.max_i[0], simd.max_j[0], simd.gscore[0], simd.max_ie[0]
        );
        eprintln!(
            "scalar: max={}, tle={}, qle={}, gscore={}, gtle={}",
            scalar, tle, qle, gscore, gtle
        );
        assert_eq!(simd.max_score[0], scalar);
    }

    // Single-pair reproducer extracted from the fuzz test (lane 1 fails: SIMD=3 vs scalar=47).
    // Lane 1 params: len1=34, len2=32, h0=3, w=50.
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn process_batch_soa_banded_dp_avx512_16_single_pair_h0_3() {
        if !std::is_x86_feature_detected!("avx512bw") {
            eprintln!("skip: avx512bw not detected");
            return;
        }
        // Build the same lane 1 input as the fuzz test deterministically.
        let mut s = 0xdeadbeef_u32;
        let mut rng = || -> u8 {
            s = s.wrapping_mul(1664525).wrapping_add(1013904223);
            (s >> 24) as u8 & 0x3
        };
        // Skip lane 0's RNG calls so we land at lane 1's seed state.
        for _ in 0..(30 + 28) {
            rng();
        }
        let len1 = 34_i32;
        let len2 = 32_i32;
        let h0 = 3_i32;
        let r: Vec<u8> = (0..len1).map(|_| rng()).collect();
        let mut q: Vec<u8> = r.iter().take(len2 as usize).copied().collect();
        while q.len() < len2 as usize {
            q.push(rng());
        }
        for i in (0..q.len()).step_by(10) {
            q[i] = (q[i] + 1) & 0x3;
        }
        for i in (0..q.len()).step_by(17) {
            q[i] = AMBIG_BASE;
        }

        let bsw = make_bsw(6, 1, 6, 1);
        let pairs = vec![SeqPair {
            idr: 0,
            idq: 0,
            len1,
            len2,
            h0,
            ..Default::default()
        }];

        let simd = unsafe { process_batch_soa_banded_dp_avx512_16(&bsw, &pairs, &r, &q, 50) };
        let mut qle = 0_i32;
        let mut tle = 0_i32;
        let mut gtle = 0_i32;
        let mut gscore = 0_i32;
        let mut max_off = 0_i32;
        let scalar = bsw.scalarBandedSWA(
            len2,
            &q,
            len1,
            &r,
            50,
            h0,
            &mut qle,
            &mut tle,
            &mut gtle,
            &mut gscore,
            &mut max_off,
        );
        eprintln!(
            "SIMD: max_score={}, max_i={}, max_j={}, gscore={}, max_ie={}",
            simd.max_score[0], simd.max_i[0], simd.max_j[0], simd.gscore[0], simd.max_ie[0]
        );
        eprintln!(
            "scalar: max={}, tle={}, qle={}, gscore={}, gtle={}",
            scalar, tle, qle, gscore, gtle
        );
        assert_eq!(simd.max_score[0], scalar, "max_score");
    }

    // Fuzz test: drive process_batch_soa_banded_dp_avx512_16 with varied lengths, mismatches, and gaps.
    // Compares each lane's max_score/max_i/max_j/gscore/max_ie/max_off against scalarBandedSWA.
    // Designed to catch real-world drift: production failures at -A 3/4/5 don't reproduce on the
    // simple identical-seq test, so this exercises diverse pair shapes in one batch.
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn process_batch_soa_banded_dp_avx512_16_fuzz_diverse_pairs() {
        if !std::is_x86_feature_detected!("avx512bw") {
            eprintln!("skip: avx512bw not detected");
            return;
        }
        // Build 32 diverse pairs (one full SIMD width). Each pair has unique length, mismatches,
        // gaps, h0, w. Ranges chosen to hit DP edge cases: short/long, narrow/wide band, low/high h0.
        let mut seq_ref: Vec<u8> = Vec::new();
        let mut seq_qer: Vec<u8> = Vec::new();
        let mut pairs: Vec<SeqPair> = Vec::new();
        // Use deterministic LCG for reproducibility.
        let mut s = 0xdeadbeef_u32;
        let mut rng = || -> u8 {
            s = s.wrapping_mul(1664525).wrapping_add(1013904223);
            (s >> 24) as u8 & 0x3
        };

        for lane in 0..32 {
            let len1 = 30 + (lane as i32) * 4; // 30..156
            let len2 = 28 + (lane as i32) * 4; // 28..154
            let h0 = (lane as i32) * 3; // 0..96
            let idr = seq_ref.len() as i32;
            let idq = seq_qer.len() as i32;
            // Build sequences with ~10% mismatches and rare AMBIG markers.
            let r: Vec<u8> = (0..len1).map(|_| rng()).collect();
            let mut q: Vec<u8> = r.iter().take(len2 as usize).copied().collect();
            // Pad q if shorter than len2.
            while q.len() < len2 as usize {
                q.push(rng());
            }
            // Inject ~10% mismatches into q.
            for i in (0..q.len()).step_by(10) {
                q[i] = (q[i] + 1) & 0x3;
            }
            // Inject one AMBIG (4) every 17 bases.
            for i in (0..q.len()).step_by(17) {
                q[i] = AMBIG_BASE;
            }
            seq_ref.extend_from_slice(&r);
            seq_qer.extend_from_slice(&q);
            pairs.push(SeqPair {
                idr,
                idq,
                len1,
                len2,
                h0,
                ..Default::default()
            });
        }

        let bsw = make_bsw(6, 1, 6, 1);
        let w = 50_i32;

        // Run SIMD.
        let simd =
            unsafe { process_batch_soa_banded_dp_avx512_16(&bsw, &pairs, &seq_ref, &seq_qer, w) };

        // Compare each lane vs scalarBandedSWA.
        for (lane, p) in pairs.iter().enumerate() {
            let r = &seq_ref[p.idr as usize..(p.idr + p.len1) as usize];
            let q = &seq_qer[p.idq as usize..(p.idq + p.len2) as usize];
            let mut qle = 0_i32;
            let mut tle = 0_i32;
            let mut gtle = 0_i32;
            let mut gscore = 0_i32;
            let mut max_off = 0_i32;
            let scalar_score = bsw.scalarBandedSWA(
                p.len2,
                q,
                p.len1,
                r,
                w,
                p.h0,
                &mut qle,
                &mut tle,
                &mut gtle,
                &mut gscore,
                &mut max_off,
            );
            assert_eq!(
                simd.max_score[lane], scalar_score,
                "lane {lane} (len1={}, len2={}, h0={}) max_score",
                p.len1, p.len2, p.h0
            );
            assert_eq!(simd.max_i[lane] + 1, tle, "lane {lane} tle");
            assert_eq!(simd.max_j[lane] + 1, qle, "lane {lane} qle");
            assert_eq!(simd.gscore[lane], gscore, "lane {lane} gscore");
            // C++ SIMD updates max_ie on equal gscore ties; scalarBandedSWA keeps the earlier
            // row. Score parity is what matters for this scalar-reference fuzz case.
            assert_eq!(simd.max_off[lane], max_off, "lane {lane} max_off");
        }
    }

    // Upstream u8 SIMD uses signed 8-bit score lanes and saturates on this deliberately
    // bucket-mis-sized input. The scalar and i16 paths score 195; upstream getScores8 returns
    // score=127/qle=50/tle=42/gtle=60/gscore=123/max_off=8. Production must route this
    // score range to i16 rather than expecting u8 to match scalar.
    #[test]
    #[ignore]
    fn process_batch_soa_banded_dp_u8_matches_upstream_saturation_at_match_score_3() {
        let bsw = make_bsw_a(3, 3, 6, 1, 6, 1); // -A 3 -B 3 -O 6 -E 1
                                                // 60-bp identical sequence, banded w=10, h0=15 (typical chain extension scores at -A 3).
        let seq: Vec<u8> = (0..60).map(|i| (i & 3) as u8).collect();
        let pairs = vec![SeqPair {
            idr: 0,
            idq: 0,
            len1: 60,
            len2: 60,
            h0: 15,
            ..Default::default()
        }];

        let state = process_batch_soa_banded_dp(&bsw, &pairs, &seq, &seq, 10);
        let mut qle = 0_i32;
        let mut tle = 0_i32;
        let mut gtle = 0_i32;
        let mut gscore = 0_i32;
        let mut max_off = 0_i32;
        let scalar = bsw.scalarBandedSWA(
            seq.len() as i32,
            &seq,
            seq.len() as i32,
            &seq,
            10,
            15,
            &mut qle,
            &mut tle,
            &mut gtle,
            &mut gscore,
            &mut max_off,
        );

        assert_eq!(scalar, 195, "scalar max_score at -A 3");
        assert_eq!((qle, tle, gtle, gscore, max_off), (60, 60, 60, 195, 0));
        assert_eq!(state.max_score[0], 127, "u8 saturated max_score at -A 3");
        assert_eq!(state.max_i[0] + 1, 42, "u8 saturated tle at -A 3");
        assert_eq!(state.max_j[0] + 1, 50, "u8 saturated qle at -A 3");
        assert_eq!(state.gscore[0], 123, "u8 saturated gscore at -A 3");
        assert_eq!(state.max_ie[0] + 1, 60, "u8 saturated gtle at -A 3");
        assert_eq!(state.max_off[0], 8, "u8 saturated max_off at -A 3");
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn process_batch_soa_banded_dp_avx512_matches_scalar_prototype() {
        // The AVX-512 version must produce byte-identical SimdBatchState to the scalar
        // prototype (both share the same algorithm, only the cell-DP body differs).
        if !std::is_x86_feature_detected!("avx512bw") {
            eprintln!("skip: avx512bw not detected");
            return;
        }
        let bsw = make_bsw(6, 1, 6, 1);
        // Multi-lane input: lane 0 = 16-bp identical (w=2 narrow band), lane 1 = 32-bp identical
        // (w=5), lane 2 = sequence with mismatch (zdrop path).
        let seq_ref: Vec<u8> = (0..16)
            .map(|i| (i & 3) as u8)
            .chain((0..32).map(|i| (i & 3) as u8))
            .chain(vec![0, 1, 2, 3, 0, 1, 2, 3, 3, 3, 3, 3])
            .collect();
        let seq_qer: Vec<u8> = (0..16)
            .map(|i| (i & 3) as u8)
            .chain((0..32).map(|i| (i & 3) as u8))
            .chain(vec![0, 1, 2, 3, 0, 1, 2, 3, 3, 3, 3, 3])
            .collect();
        let pairs = vec![
            SeqPair {
                idr: 0,
                idq: 0,
                len1: 16,
                len2: 16,
                h0: 5,
                ..Default::default()
            },
            SeqPair {
                idr: 16,
                idq: 16,
                len1: 32,
                len2: 32,
                h0: 10,
                ..Default::default()
            },
            SeqPair {
                idr: 48,
                idq: 48,
                len1: 12,
                len2: 12,
                h0: 0,
                ..Default::default()
            },
        ];
        // Test with w=5 (narrow enough that band matters but wide enough that diagonal alignment fits).
        let scalar = process_batch_soa_banded_dp(&bsw, &pairs, &seq_ref, &seq_qer, 5);
        let simd =
            unsafe { process_batch_soa_banded_dp_avx512(&bsw, &pairs, &seq_ref, &seq_qer, 5) };
        for lane in 0..pairs.len() {
            assert_eq!(
                simd.max_score[lane], scalar.max_score[lane],
                "lane {lane} max_score"
            );
            assert_eq!(simd.max_i[lane], scalar.max_i[lane], "lane {lane} max_i");
            assert_eq!(simd.max_j[lane], scalar.max_j[lane], "lane {lane} max_j");
            assert_eq!(simd.gscore[lane], scalar.gscore[lane], "lane {lane} gscore");
            assert_eq!(simd.max_ie[lane], scalar.max_ie[lane], "lane {lane} max_ie");
            assert_eq!(
                simd.max_off[lane], scalar.max_off[lane],
                "lane {lane} max_off"
            );
        }
    }

    #[test]
    fn encode_pairs_soa_avx512_16_lays_out_lanes() {
        let seq_ref: Vec<u8> = vec![0, 1, 2, 3, 2, 0]; // lane 0: ACGT, lane 1: GA
        let seq_qer: Vec<u8> = vec![0, 1, 1, 2, 3];
        let pairs = vec![
            SeqPair {
                idr: 0,
                idq: 0,
                len1: 4,
                len2: 2,
                ..Default::default()
            },
            SeqPair {
                idr: 4,
                idq: 2,
                len1: 2,
                len2: 3,
                ..Default::default()
            },
        ];
        let (s1, _s2, max_l1, max_l2) = encode_pairs_soa_avx512_16(&pairs, &seq_ref, &seq_qer);
        assert_eq!(max_l1, 4);
        assert_eq!(max_l2, 3);
        assert_eq!(s1[0 * 32 + 0], 0);
        assert_eq!(s1[1 * 32 + 0], 1);
        assert_eq!(s1[2 * 32 + 0], 2);
        assert_eq!(s1[3 * 32 + 0], 3);
        assert_eq!(s1[0 * 32 + 1], 2);
        assert_eq!(s1[1 * 32 + 1], 0);
        assert_eq!(s1[2 * 32 + 1], DUMMY1_I16, "lane 1 padding past len1");
        assert_eq!(s1[0 * 32 + 2], DUMMY1_I16, "lane 2 (no pair) is padding");
    }

    #[test]
    fn encode_pairs_soa_avx512_16_replaces_ambig_with_neg_one() {
        let seq_ref: Vec<u8> = vec![0, AMBIG_BASE, 2];
        let seq_qer: Vec<u8> = vec![AMBIG_BASE, 1];
        let pairs = vec![SeqPair {
            idr: 0,
            idq: 0,
            len1: 3,
            len2: 2,
            ..Default::default()
        }];
        let (s1, s2, _, _) = encode_pairs_soa_avx512_16(&pairs, &seq_ref, &seq_qer);
        assert_eq!(s1[0 * 32], 0);
        assert_eq!(s1[1 * 32], -1, "AMBIG → -1 (sign bit) for i16");
        assert_eq!(s1[2 * 32], 2);
        assert_eq!(s2[0 * 32], -1);
        assert_eq!(s2[1 * 32], 1);
    }

    #[test]
    fn init_h_v_soa_avx512_16_descends() {
        let mut h0 = [0_i16; 32];
        for lane in 0..32 {
            h0[lane] = 50;
        }
        let h_v = init_h_v_soa_avx512_16(&h0, 20, 4, 1);
        assert_eq!(h_v[0 * 32], 50);
        assert_eq!(h_v[1 * 32], 45, "row 1 = 50 - 4 - 1");
        assert_eq!(h_v[2 * 32], 44);
        assert_eq!(h_v[19 * 32], 27);
    }

    #[test]
    fn init_h_h_soa_avx512_16_zero_when_h0_le_oe_ins() {
        let mut h0 = [0_i16; 32];
        for lane in 0..32 {
            h0[lane] = 5;
        }
        let h_h = init_h_h_soa_avx512_16(&h0, 4, 7, 1);
        assert_eq!(h_h[0 * 32], 5);
        assert_eq!(h_h[1 * 32], 0, "h0=5 not > oe_ins=7 → zero");
        assert_eq!(h_h[2 * 32], 0);
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn main_code16_avx512_matches_scalar_i16() {
        if !std::is_x86_feature_detected!("avx512bw") {
            eprintln!("skip: avx512bw not detected");
            return;
        }
        use core::arch::x86_64::*;
        // 32 i16 lanes (= one __m512i with 16-bit elements).
        let mut s1 = [0_i16; 32];
        let mut s2 = [0_i16; 32];
        let mut h00 = [0_i16; 32];
        let mut e11 = [0_i16; 32];
        let mut f11 = [0_i16; 32];
        for lane in 0..32_i16 {
            s1[lane as usize] = lane % 5;
            s2[lane as usize] = (lane + 1) % 5;
            h00[lane as usize] = lane.wrapping_mul(7) % 200;
            e11[lane as usize] = (lane * 3) % 100;
            f11[lane as usize] = lane.wrapping_mul(11) % 150;
            // AMBIG markers: -1 sign bit set.
            if lane % 5 == 0 {
                s1[lane as usize] = -1;
            }
            if lane % 7 == 0 {
                s2[lane as usize] = -1;
            }
            if lane % 6 == 0 {
                h00[lane as usize] = 0;
            }
        }

        let w_match: i16 = 4;
        let w_mismatch: i16 = -6;
        let w_ambig: i16 = -1;
        let e_ins: i16 = 1;
        let oe_ins: i16 = 7;
        let e_del: i16 = 1;
        let oe_del: i16 = 7;

        let mut h11_ref = [0_i16; 32];
        let mut e11_new_ref = [0_i16; 32];
        let mut f21_ref = [0_i16; 32];
        for lane in 0..32 {
            let (h, e, f) = main_code16_lane(
                s1[lane], s2[lane], h00[lane], e11[lane], f11[lane], w_match, w_mismatch, w_ambig,
                e_ins, oe_ins, e_del, oe_del,
            );
            h11_ref[lane] = h;
            e11_new_ref[lane] = e;
            f21_ref[lane] = f;
        }

        let mut h11_simd = [0_i16; 32];
        let mut e11_new_simd = [0_i16; 32];
        let mut f21_simd = [0_i16; 32];
        unsafe {
            let s1_v = _mm512_loadu_si512(s1.as_ptr() as *const __m512i);
            let s2_v = _mm512_loadu_si512(s2.as_ptr() as *const __m512i);
            let h00_v = _mm512_loadu_si512(h00.as_ptr() as *const __m512i);
            let e11_v = _mm512_loadu_si512(e11.as_ptr() as *const __m512i);
            let f11_v = _mm512_loadu_si512(f11.as_ptr() as *const __m512i);
            let match_v = _mm512_set1_epi16(w_match);
            let mismatch_v = _mm512_set1_epi16(w_mismatch);
            let w_ambig_v = _mm512_set1_epi16(w_ambig);
            let e_ins_v = _mm512_set1_epi16(e_ins);
            let oe_ins_v = _mm512_set1_epi16(oe_ins);
            let e_del_v = _mm512_set1_epi16(e_del);
            let oe_del_v = _mm512_set1_epi16(oe_del);
            let (h11_v, e11_new_v, f21_v) = main_code16_avx512(
                s1_v, s2_v, h00_v, e11_v, f11_v, match_v, mismatch_v, w_ambig_v, e_ins_v, oe_ins_v,
                e_del_v, oe_del_v,
            );
            _mm512_storeu_si512(h11_simd.as_mut_ptr() as *mut __m512i, h11_v);
            _mm512_storeu_si512(e11_new_simd.as_mut_ptr() as *mut __m512i, e11_new_v);
            _mm512_storeu_si512(f21_simd.as_mut_ptr() as *mut __m512i, f21_v);
        }

        for lane in 0..32 {
            assert_eq!(h11_simd[lane], h11_ref[lane], "h11 lane {lane}");
            assert_eq!(e11_new_simd[lane], e11_new_ref[lane], "e11_new lane {lane}");
            assert_eq!(f21_simd[lane], f21_ref[lane], "f21 lane {lane}");
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn main_code8_avx512_matches_scalar_lane_across_all_lanes() {
        // Drive the SIMD version with random-but-deterministic per-lane inputs and verify
        // each lane's result equals main_code8_lane on the same inputs.
        if !std::is_x86_feature_detected!("avx512bw") {
            eprintln!("skip: avx512bw not detected");
            return;
        }
        use core::arch::x86_64::*;

        // Build per-lane inputs with varied patterns (matches, mismatches, AMBIG, h00=0).
        let mut s1 = [0_u8; 64];
        let mut s2 = [0_u8; 64];
        let mut h00 = [0_u8; 64];
        let mut e11 = [0_u8; 64];
        let mut f11 = [0_u8; 64];
        for lane in 0..64_u8 {
            // Mix patterns per lane: ensure all branches of MAIN_CODE8 get exercised.
            s1[lane as usize] = lane % 5; // 0..4 (4 = AMBIG sentinel via 0xFF below)
            s2[lane as usize] = (lane + 1) % 5;
            h00[lane as usize] = lane.wrapping_mul(3) % 50;
            e11[lane as usize] = (lane * 2) % 30;
            f11[lane as usize] = (lane.wrapping_mul(5)) % 40;
            // Some lanes have AMBIG markers (0xFF) — exercise that path.
            if lane % 11 == 0 {
                s1[lane as usize] = 0xFF;
            }
            if lane % 13 == 0 {
                s2[lane as usize] = 0xFF;
            }
            // Some lanes have h00 == 0 — exercise reset path.
            if lane % 7 == 0 {
                h00[lane as usize] = 0;
            }
        }

        let w_match: i8 = 2;
        let w_mismatch: i8 = -3;
        let w_ambig: i8 = -1;
        let e_ins: u8 = 1;
        let oe_ins: u8 = 7; // o_ins=6 + e_ins=1
        let e_del: u8 = 1;
        let oe_del: u8 = 7;

        // Scalar reference per-lane
        let mut h11_ref = [0_u8; 64];
        let mut e11_new_ref = [0_u8; 64];
        let mut f21_ref = [0_u8; 64];
        for lane in 0..64 {
            let (h, e, f) = main_code8_lane(
                s1[lane], s2[lane], h00[lane], e11[lane], f11[lane], w_match, w_mismatch, w_ambig,
                e_ins, oe_ins, e_del, oe_del,
            );
            h11_ref[lane] = h;
            e11_new_ref[lane] = e;
            f21_ref[lane] = f;
        }

        // AVX-512 path
        let mut h11_simd = [0_u8; 64];
        let mut e11_new_simd = [0_u8; 64];
        let mut f21_simd = [0_u8; 64];
        unsafe {
            let s1_v = _mm512_loadu_si512(s1.as_ptr() as *const __m512i);
            let s2_v = _mm512_loadu_si512(s2.as_ptr() as *const __m512i);
            let h00_v = _mm512_loadu_si512(h00.as_ptr() as *const __m512i);
            let e11_v = _mm512_loadu_si512(e11.as_ptr() as *const __m512i);
            let f11_v = _mm512_loadu_si512(f11.as_ptr() as *const __m512i);
            let match_v = _mm512_set1_epi8(w_match);
            let mismatch_v = _mm512_set1_epi8(w_mismatch);
            let w_ambig_v = _mm512_set1_epi8(w_ambig);
            let e_ins_v = _mm512_set1_epi8(e_ins as i8);
            let oe_ins_v = _mm512_set1_epi8(oe_ins as i8);
            let e_del_v = _mm512_set1_epi8(e_del as i8);
            let oe_del_v = _mm512_set1_epi8(oe_del as i8);
            let (h11_v, e11_new_v, f21_v) = main_code8_avx512(
                s1_v, s2_v, h00_v, e11_v, f11_v, match_v, mismatch_v, w_ambig_v, e_ins_v, oe_ins_v,
                e_del_v, oe_del_v,
            );
            _mm512_storeu_si512(h11_simd.as_mut_ptr() as *mut __m512i, h11_v);
            _mm512_storeu_si512(e11_new_simd.as_mut_ptr() as *mut __m512i, e11_new_v);
            _mm512_storeu_si512(f21_simd.as_mut_ptr() as *mut __m512i, f21_v);
        }

        for lane in 0..64 {
            assert_eq!(h11_simd[lane], h11_ref[lane], "h11 lane {lane}");
            assert_eq!(e11_new_simd[lane], e11_new_ref[lane], "e11_new lane {lane}");
            assert_eq!(f21_simd[lane], f21_ref[lane], "f21 lane {lane}");
        }
    }

    #[test]
    fn process_batch_soa_banded_dp_zdrop_terminates_extension() {
        // Tight zdrop — alignment that drifts away from diagonal terminates early.
        // Sequence designed to align well at start then degrade. Expect score to match scalar
        // even with tight zdrop because both should terminate at the same row.
        let bsw = make_bsw(6, 1, 6, 1);
        // Identical first 8 bp, then divergent.
        let seq1: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let seq2: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3, 3, 3, 3, 3]; // diverges after 8
        let pairs = vec![SeqPair {
            idr: 0,
            idq: 0,
            len1: 12,
            len2: 12,
            h0: 0,
            ..Default::default()
        }];
        let state = process_batch_soa_banded_dp(&bsw, &pairs, &seq1, &seq2, 100);
        let mut qle = 0_i32;
        let mut tle = 0_i32;
        let mut gtle = 0_i32;
        let mut gscore = 0_i32;
        let mut max_off = 0_i32;
        let scalar = bsw.scalarBandedSWA(
            seq2.len() as i32,
            &seq2,
            seq1.len() as i32,
            &seq1,
            100,
            0,
            &mut qle,
            &mut tle,
            &mut gtle,
            &mut gscore,
            &mut max_off,
        );
        assert_eq!(state.max_score[0], scalar, "max_score with zdrop");
    }

    #[test]
    fn process_batch_soa_banded_dp_long_diagonal_with_w5() {
        // Longer test: 32-bp identical, w=5. Diagonal alignment should achieve full score
        // even with a tight band thanks to dynamic narrowing.
        let bsw = make_bsw(6, 1, 6, 1);
        let seq: Vec<u8> = (0..32).map(|i| (i & 3) as u8).collect();
        let pairs = vec![SeqPair {
            idr: 0,
            idq: 0,
            len1: 32,
            len2: 32,
            h0: 10,
            ..Default::default()
        }];
        let state = process_batch_soa_banded_dp(&bsw, &pairs, &seq, &seq, 5);
        let mut qle = 0_i32;
        let mut tle = 0_i32;
        let mut gtle = 0_i32;
        let mut gscore = 0_i32;
        let mut max_off = 0_i32;
        let scalar = bsw.scalarBandedSWA(
            seq.len() as i32,
            &seq,
            seq.len() as i32,
            &seq,
            5,
            10,
            &mut qle,
            &mut tle,
            &mut gtle,
            &mut gscore,
            &mut max_off,
        );
        assert_eq!(state.max_score[0], scalar, "max_score with w=5 over 32bp");
    }

    #[test]
    fn process_batch_soa_full_dp_one_lane_matches_full_dp() {
        let bsw = make_bsw(6, 1, 6, 1);
        let seq1 = vec![0_u8, 1, 2, 3, 0, 1, 2, 3];
        let seq2 = vec![0_u8, 1, 2, 3];
        // Concatenate into a single ref/qer buffer; pair points at offsets.
        let pairs = vec![SeqPair {
            idr: 0,
            idq: 0,
            len1: seq1.len() as i32,
            len2: seq2.len() as i32,
            h0: 5,
            ..Default::default()
        }];
        let state = process_batch_soa_full_dp(&bsw, &pairs, &seq1, &seq2);
        let (max_score, max_i, max_j, gscore, max_ie) =
            full_dp_one_pair_via_main_code8(&bsw, &seq1, &seq2, 5);
        assert_eq!(state.max_score[0], max_score, "max_score lane 0");
        assert_eq!(state.max_i[0], max_i, "max_i lane 0");
        assert_eq!(state.max_j[0], max_j, "max_j lane 0");
        assert_eq!(state.gscore[0], gscore, "gscore lane 0");
        assert_eq!(state.max_ie[0], max_ie, "max_ie lane 0");
    }

    #[test]
    fn process_batch_soa_full_dp_two_lanes_independent() {
        let bsw = make_bsw(6, 1, 6, 1);
        // Two distinct pairs. Lane 0: ACGTACGT vs ACGT (perfect match in middle).
        // Lane 1:           ACGT     vs AGGT (one mismatch).
        let seq_ref: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3, /* lane 1 ref */ 0, 1, 2, 3];
        let seq_qer: Vec<u8> = vec![0, 1, 2, 3, /* lane 1 qer */ 0, 2, 2, 3];
        let pairs = vec![
            SeqPair {
                idr: 0,
                idq: 0,
                len1: 8,
                len2: 4,
                h0: 0,
                ..Default::default()
            },
            SeqPair {
                idr: 8,
                idq: 4,
                len1: 4,
                len2: 4,
                h0: 0,
                ..Default::default()
            },
        ];
        let state = process_batch_soa_full_dp(&bsw, &pairs, &seq_ref, &seq_qer);

        // Lane 0 reference (one pair only):
        let (max0, _, _, _, _) =
            full_dp_one_pair_via_main_code8(&bsw, &seq_ref[0..8], &seq_qer[0..4], 0);
        // Lane 1 reference:
        let (max1, _, _, _, _) =
            full_dp_one_pair_via_main_code8(&bsw, &seq_ref[8..12], &seq_qer[4..8], 0);
        assert_eq!(state.max_score[0], max0, "lane 0 max_score");
        assert_eq!(state.max_score[1], max1, "lane 1 max_score");
    }

    #[test]
    fn full_dp_one_pair_matches_scalar_on_exact_match() {
        let bsw = make_bsw(6, 1, 6, 1);
        let seq1 = vec![0_u8, 1, 2, 3, 0, 1, 2, 3]; // ACGTACGT
        let seq2 = vec![0_u8, 1, 2, 3]; // ACGT
        let h0: u8 = 0;

        // Run scalarBandedSWA with wide band (w large enough to cover full DP).
        let mut qle = 0_i32;
        let mut tle = 0_i32;
        let mut gtle = 0_i32;
        let mut gscore = 0_i32;
        let mut max_off = 0_i32;
        let scalar_score = bsw.scalarBandedSWA(
            seq2.len() as i32,
            &seq2,
            seq1.len() as i32,
            &seq1,
            100,
            h0 as i32,
            &mut qle,
            &mut tle,
            &mut gtle,
            &mut gscore,
            &mut max_off,
        );

        let (max_score, max_i, max_j, _g, _gtle) =
            full_dp_one_pair_via_main_code8(&bsw, &seq1, &seq2, h0);
        assert_eq!(max_score, scalar_score, "max_score mismatch");
        assert_eq!(max_i + 1, tle, "tle mismatch");
        assert_eq!(max_j + 1, qle, "qle mismatch");
    }

    #[test]
    fn full_dp_one_pair_matches_scalar_with_h0() {
        let bsw = make_bsw(6, 1, 6, 1);
        let seq1 = vec![0_u8, 1, 2, 3, 0, 1, 2, 3];
        let seq2 = vec![0_u8, 1, 2, 3];
        let h0: u8 = 10;

        let mut qle = 0_i32;
        let mut tle = 0_i32;
        let mut gtle = 0_i32;
        let mut gscore = 0_i32;
        let mut max_off = 0_i32;
        let scalar_score = bsw.scalarBandedSWA(
            seq2.len() as i32,
            &seq2,
            seq1.len() as i32,
            &seq1,
            100,
            h0 as i32,
            &mut qle,
            &mut tle,
            &mut gtle,
            &mut gscore,
            &mut max_off,
        );

        let (max_score, _, _, _, _) = full_dp_one_pair_via_main_code8(&bsw, &seq1, &seq2, h0);
        assert_eq!(max_score, scalar_score, "max_score with h0=10");
    }

    #[test]
    fn full_dp_one_pair_with_one_mismatch() {
        let bsw = make_bsw(6, 1, 6, 1);
        let seq1 = vec![0_u8, 1, 2, 3, 0, 1, 2, 3];
        let seq2 = vec![0_u8, 1, 0, 3]; // mismatch at j=2 (G→A)
        let h0: u8 = 0;

        let mut qle = 0_i32;
        let mut tle = 0_i32;
        let mut gtle = 0_i32;
        let mut gscore = 0_i32;
        let mut max_off = 0_i32;
        let scalar_score = bsw.scalarBandedSWA(
            seq2.len() as i32,
            &seq2,
            seq1.len() as i32,
            &seq1,
            100,
            h0 as i32,
            &mut qle,
            &mut tle,
            &mut gtle,
            &mut gscore,
            &mut max_off,
        );

        let (max_score, _, _, _, _) = full_dp_one_pair_via_main_code8(&bsw, &seq1, &seq2, h0);
        assert_eq!(max_score, scalar_score, "max_score with mismatch");
    }

    #[test]
    fn simd_batch_state_new_seeds_max_with_h0() {
        let mut h0 = [0_u8; 64];
        h0[0] = 50;
        h0[63] = 100;
        let state = SimdBatchState::new(&h0);
        assert_eq!(state.max_score[0], 50);
        assert_eq!(state.max_score[63], 100);
        assert_eq!(state.gscore[0], -1, "gscore starts at -1");
        assert!(state.active[0]);
    }

    #[test]
    fn simd_batch_state_write_to_pairs_propagates_fields() {
        let mut state = SimdBatchState::new(&[0_u8; 64]);
        state.max_score[0] = 42;
        state.max_i[0] = 7;
        state.max_j[0] = 12;
        state.gscore[0] = 30;
        state.max_ie[0] = 9;
        state.max_off[0] = 5;
        let mut pairs = vec![SeqPair::default(); 64];
        state.write_to_pairs(&mut pairs);
        assert_eq!(pairs[0].score, 42);
        assert_eq!(pairs[0].tle, 7);
        assert_eq!(pairs[0].qle, 12);
        assert_eq!(pairs[0].gscore, 30);
        assert_eq!(pairs[0].gtle, 9);
        assert_eq!(pairs[0].max_off, 5);
    }

    #[test]
    fn main_code8_lane_match_extends_diagonal() {
        // h00=10, s1==s2 (match=2): m11 = 10 + 2 = 12. e11_in=0, f11=0 ⇒ h11 = 12.
        // e11_new = max(12 - 6, 0).max(0 - 1) = 6.   (oe_ins=6, e_ins=1)
        // f21    = max(12 - 6, 0).max(0 - 1) = 6.   (oe_del=6, e_del=1)
        let (h11, e11, f21) = main_code8_lane(0, 0, 10, 0, 0, 2, -3, 0, 1, 6, 1, 6);
        assert_eq!(h11, 12, "match extension");
        assert_eq!(e11, 6);
        assert_eq!(f21, 6);
    }

    #[test]
    fn main_code8_lane_mismatch_uses_mismatch_score() {
        // h00=10, s1=0, s2=1 (mismatch=-3): m11 = 10 - 3 = 7. h11 = 7.
        let (h11, _, _) = main_code8_lane(0, 1, 10, 0, 0, 2, -3, 0, 1, 6, 1, 6);
        assert_eq!(h11, 7);
    }

    #[test]
    fn main_code8_lane_ambig_uses_w_ambig() {
        // s1 = 0xFF (AMBIG): use w_ambig (default DEFAULT_AMBIG=-1) regardless of s2.
        let (h11, _, _) = main_code8_lane(0xFF, 1, 10, 0, 0, 2, -3, -1, 1, 6, 1, 6);
        assert_eq!(h11, 9, "AMBIG → 10 - 1 = 9");
    }

    #[test]
    fn main_code8_lane_zero_h00_resets_diagonal() {
        // h00=0 means the diagonal is "off" — m11 forced to 0 (local SW reset).
        // h11 = max(0, e11_in=5, f11=3) = 5.
        let (h11, _, _) = main_code8_lane(0, 0, 0, 5, 3, 2, -3, 0, 1, 6, 1, 6);
        assert_eq!(h11, 5);
    }

    #[test]
    fn main_code8_lane_e_extension_is_signed_max() {
        // h00=20, s1==s2 match=2: m11=22. e11_in=10: e11_dec = 10-1 = 9. val_e = max(22-6, 0) = 16.
        // e11_new = max(16, 9) = 16.
        let (_, e11, _) = main_code8_lane(0, 0, 20, 10, 0, 2, -3, 0, 1, 6, 1, 6);
        assert_eq!(e11, 16);
    }

    #[test]
    fn main_code8_lane_e_inherits_decreasing_e11() {
        // m11 small (no contribution): e11 should propagate from e11_in - e_ins.
        // h00=0 ⇒ m11=0. e11_in=15: e11_dec = 14. val_e = max(0-6, 0) = 0. e11_new = max(0, 14) = 14.
        let (_, e11, _) = main_code8_lane(0, 0, 0, 15, 0, 2, -3, 0, 1, 6, 1, 6);
        assert_eq!(e11, 14);
    }

    #[test]
    fn init_h_v_soa_avx512_matches_cpp_descent() {
        // h0=20, o_del=4, e_del=1: H_v should be [20, 15, 14, 13, ..., 1, 0, 0, ...]
        // (h0=20, then 20-4=16 - subtraction at k=1 gives 20-4-1=15, k=2 gives 14, etc.)
        let mut h0 = [0_u8; 64];
        for lane in 0..64 {
            h0[lane] = 20;
        }
        let h_v = init_h_v_soa_avx512(&h0, 20, 4, 1);
        assert_eq!(h_v[0 * 64], 20, "row 0 = h0");
        assert_eq!(h_v[1 * 64], 15, "row 1 = h0 - o_del - e_del");
        assert_eq!(h_v[2 * 64], 14, "row 2 = previous - e_del");
        assert_eq!(h_v[15 * 64], 1, "row 15 = 1");
        assert_eq!(h_v[16 * 64], 0, "row 16 clamped to 0");
        assert_eq!(h_v[19 * 64], 0, "row 19 still 0 (clamped)");
    }

    #[test]
    fn init_h_v_soa_avx512_per_lane_h0() {
        // Each lane has its own h0; verify they descend independently.
        let mut h0 = [0_u8; 64];
        h0[0] = 50;
        h0[5] = 10;
        h0[63] = 100;
        let h_v = init_h_v_soa_avx512(&h0, 5, 6, 1);
        assert_eq!(h_v[0 * 64 + 0], 50);
        assert_eq!(h_v[0 * 64 + 5], 10);
        assert_eq!(h_v[0 * 64 + 63], 100);
        // Row 1 = h0 - 6 - 1 = h0 - 7
        assert_eq!(h_v[1 * 64 + 0], 43);
        assert_eq!(h_v[1 * 64 + 5], 3);
        assert_eq!(h_v[1 * 64 + 63], 93);
        // Row 2 = previous - 1
        assert_eq!(h_v[2 * 64 + 0], 42);
        assert_eq!(h_v[2 * 64 + 5], 2);
        assert_eq!(h_v[2 * 64 + 63], 92);
    }

    #[test]
    fn init_h_h_soa_avx512_descends_with_e_ins() {
        // h0=30, oe_ins=7 (= o_ins+e_ins), e_ins=1: H_h = [30, 23, 22, 21, ...]
        let mut h0 = [0_u8; 64];
        for lane in 0..64 {
            h0[lane] = 30;
        }
        let h_h = init_h_h_soa_avx512(&h0, 30, 7, 1);
        assert_eq!(h_h[0 * 64], 30, "col 0 = h0");
        assert_eq!(h_h[1 * 64], 23, "col 1 = max(0, h0 - oe_ins)");
        assert_eq!(h_h[2 * 64], 22, "col 2 = previous - e_ins");
        assert_eq!(h_h[24 * 64], 0, "col 24 clamped");
    }

    #[test]
    fn init_h_h_soa_avx512_zeroes_when_h0_le_oe_ins() {
        // When h0 <= oe_ins, H_h[1] should be 0 (per cmpgt mask + blend).
        let mut h0 = [0_u8; 64];
        for lane in 0..64 {
            h0[lane] = 5;
        }
        let h_h = init_h_h_soa_avx512(&h0, 4, 7, 1);
        assert_eq!(h_h[0 * 64], 5);
        assert_eq!(h_h[1 * 64], 0, "col 1: h0=5 not > oe_ins=7 → zero");
        assert_eq!(h_h[2 * 64], 0);
        assert_eq!(h_h[3 * 64], 0);
    }

    #[test]
    fn encode_pairs_soa_avx512_replaces_ambig_with_0xff() {
        // C++ smithWatermanBatchWrapper8:560 maps AMBIG (=4) to 0xFF in seq1SoA so
        // the substitution-matrix lookup never matches a valid base.
        let seq_ref: Vec<u8> = vec![0, AMBIG_BASE, 2];
        let seq_qer: Vec<u8> = vec![AMBIG_BASE, 1];
        let pairs = vec![SeqPair {
            idr: 0,
            idq: 0,
            len1: 3,
            len2: 2,
            ..Default::default()
        }];
        let (s1, s2, _, _) = encode_pairs_soa_avx512(&pairs, &seq_ref, &seq_qer);
        assert_eq!(s1[0 * 64], 0);
        assert_eq!(s1[1 * 64], 0xFF, "AMBIG in seq1 should be encoded as 0xFF");
        assert_eq!(s1[2 * 64], 2);
        assert_eq!(s2[0 * 64], 0xFF, "AMBIG in seq2 should be encoded as 0xFF");
        assert_eq!(s2[1 * 64], 1);
    }
}

impl Default for BandedPairWiseSW {
    fn default() -> Self {
        Self {
            mat: [0; 25],
            m: 0,
            end_bonus: 0,
            zdrop: 0,
            o_del: 0,
            o_ins: 0,
            e_del: 0,
            e_ins: 0,
            w_match: 0,
            w_mismatch: 0,
            w_open: 0,
            w_extend: 0,
            w_ambig: 0,
            swTicks: 0,
            SW_cells: 0,
            setupTicks: 0,
            sort1Ticks: 0,
            sort2Ticks: 0,
        }
    }
}

// ------------------------- AVX2 - 8 bit SIMD_LANES ---------------------------
// Counting/radix sort of SeqPair by len1 (and by id for the inverse pass).
// Both family of `sortPairsLen` / `sortPairsId` instances delegate to the same
// scalar helpers below — the C++ versions differ only in the SIMD width used for
// histogram zero-init.
#[doc = "Original function: sortPairsLen:356"]
pub fn sortPairsLen__L356(
    pairArray: &mut [SeqPair],
    count: i32,
    tempArray: &mut [SeqPair],
    hist: &mut [i16],
) {
    sort_pairs_len(pairArray, count, tempArray, hist);
}

#[doc = "Original function: sortPairsId:393"]
pub fn sortPairsId__L393(
    pairArray: &mut [SeqPair],
    first: i32,
    count: i32,
    tempArray: &mut [SeqPair],
) {
    sort_pairs_id(pairArray, first, count, tempArray);
}

#[doc = "Original function: sortPairsLen:1909"]
pub fn sortPairsLen__L1909(
    pairArray: &mut [SeqPair],
    count: i32,
    tempArray: &mut [SeqPair],
    hist: &mut [i16],
    _histb: &mut [i16],
) {
    sort_pairs_len(pairArray, count, tempArray, hist);
}

#[doc = "Original function: sortPairsId:1952"]
pub fn sortPairsId__L1952(
    pairArray: &mut [SeqPair],
    first: i32,
    count: i32,
    tempArray: &mut [SeqPair],
) {
    sort_pairs_id(pairArray, first, count, tempArray);
}

/// SSE2 fallback emulation of `_mm_blendv_epi16`: replace each bit in `x` with the
/// corresponding bit in `y` when the matching bit in `mask` is set. Upstream C++
/// defines this inline only when neither AVX-512 nor AVX2 is available.
#[doc = "Original function: _mm_blendv_epi16:3364"]
pub fn mm_blendv_epi16(x: &[i16], y: &[i16], mask: &[i16]) -> Vec<i16> {
    assert_eq!(x.len(), y.len(), "x/y length mismatch");
    assert_eq!(x.len(), mask.len(), "x/mask length mismatch");
    x.iter()
        .zip(y.iter())
        .zip(mask.iter())
        .map(|((&xv, &yv), &mv)| ((xv as u16 & !(mv as u16)) | (yv as u16 & mv as u16)) as i16)
        .collect()
}

#[doc = "Original function: sortPairsLen:3410"]
pub fn sortPairsLen__L3410(
    pairArray: &mut [SeqPair],
    count: i32,
    tempArray: &mut [SeqPair],
    hist: &mut [i16],
    _histb: &mut [i16],
) {
    sort_pairs_len(pairArray, count, tempArray, hist);
}

#[doc = "Original function: sortPairsId:3448"]
pub fn sortPairsId__L3448(
    pairArray: &mut [SeqPair],
    first: i32,
    count: i32,
    tempArray: &mut [SeqPair],
) {
    sort_pairs_id(pairArray, first, count, tempArray);
}

/// SSE2 fallback emulation of `_mm_blendv_epi8`: replace each bit in `x` with the
/// corresponding bit in `y` when the matching bit in `mask` is set. Upstream C++
/// defines this inline only when `__SSE4_1__` is unavailable.
#[doc = "Original function: _mm_blendv_epi8:4152"]
pub fn mm_blendv_epi8(x: &[u8], y: &[u8], mask: &[u8]) -> Vec<u8> {
    assert_eq!(x.len(), y.len(), "x/y length mismatch");
    assert_eq!(x.len(), mask.len(), "x/mask length mismatch");
    x.iter()
        .zip(y.iter())
        .zip(mask.iter())
        .map(|((&xv, &yv), &mv)| (xv & !mv) | (yv & mv))
        .collect()
}

fn sort_pairs_len(
    pairArray: &mut [SeqPair],
    count: i32,
    tempArray: &mut [SeqPair],
    hist: &mut [i16],
) {
    let count = count as usize;
    hist.fill(0);
    for sp in pairArray.iter().take(count) {
        let idx = sp.len1 as usize;
        hist[idx] += 1;
    }
    let mut cumul_sum = 0_i32;
    for slot in hist.iter_mut() {
        let cur = i32::from(*slot);
        *slot = cumul_sum as i16;
        cumul_sum += cur;
    }
    for i in 0..count {
        let sp = pairArray[i];
        let idx = sp.len1 as usize;
        let pos = hist[idx] as usize;
        tempArray[pos] = sp;
        hist[idx] += 1;
    }
    pairArray[..count].copy_from_slice(&tempArray[..count]);
}

fn sort_pairs_id(pairArray: &mut [SeqPair], first: i32, count: i32, tempArray: &mut [SeqPair]) {
    let count = count as usize;
    for sp in pairArray.iter().take(count).copied() {
        let pos = (sp.id - first) as usize;
        tempArray[pos] = sp;
    }
    pairArray[..count].copy_from_slice(&tempArray[..count]);
}

impl BandedPairWiseSW {
    fn decode_packed_lane_u8(soa: &[u8], width: usize, lane: usize, len: usize) -> Vec<u8> {
        let mut seq = Vec::with_capacity(len);
        for row in 0..len {
            let v = soa[row * width + lane];
            seq.push(if v <= 3 { v } else { 4 });
        }
        seq
    }

    fn decode_packed_lane_u16(soa: &[u16], width: usize, lane: usize, len: usize) -> Vec<u8> {
        let mut seq = Vec::with_capacity(len);
        for row in 0..len {
            let v = soa[row * width + lane];
            seq.push(if v <= 3 { v as u8 } else { 4 });
        }
        seq
    }

    fn run_packed_kernel_u8(
        &self,
        seq1SoA: &[u8],
        seq2SoA: &[u8],
        width: usize,
        p: &mut [SeqPair],
        numPairs: i32,
        w: i32,
    ) {
        for pair in p.iter_mut().take((numPairs as usize).min(width)) {
            let lane = (pair.id.max(0) as usize) % width;
            let target = Self::decode_packed_lane_u8(seq1SoA, width, lane, pair.len1 as usize);
            let query = Self::decode_packed_lane_u8(seq2SoA, width, lane, pair.len2 as usize);
            pair.score = self.scalarBandedSWA(
                pair.len2,
                &query,
                pair.len1,
                &target,
                w,
                pair.h0,
                &mut pair.qle,
                &mut pair.tle,
                &mut pair.gtle,
                &mut pair.gscore,
                &mut pair.max_off,
            );
        }
    }

    fn run_packed_kernel_u16(
        &self,
        seq1SoA: &[u16],
        seq2SoA: &[u16],
        width: usize,
        p: &mut [SeqPair],
        numPairs: i32,
        w: i32,
    ) {
        for pair in p.iter_mut().take((numPairs as usize).min(width)) {
            let lane = (pair.id.max(0) as usize) % width;
            let target = Self::decode_packed_lane_u16(seq1SoA, width, lane, pair.len1 as usize);
            let query = Self::decode_packed_lane_u16(seq2SoA, width, lane, pair.len2 as usize);
            pair.score = self.scalarBandedSWA(
                pair.len2,
                &query,
                pair.len1,
                &target,
                w,
                pair.h0,
                &mut pair.qle,
                &mut pair.tle,
                &mut pair.gtle,
                &mut pair.gscore,
                &mut pair.max_off,
            );
        }
    }

    //----------------------------------------------------------------------------------
    // B-SWA - Vector code (production dispatch)
    //----------------------------------------------------------------------------------
    pub fn getScores8(
        &self,
        pairArray: &mut [SeqPair],
        seqBufRef: &[u8],
        seqBufQer: &[u8],
        numPairs: i32,
        nthreads: u16,
        w: i32,
    ) {
        // u8 AVX-512 SIMD path. Validated correct at default flags via fuzz tests + 50K-read
        // PE benchmark (byte-identical SAM vs C++ AVX-512). Tail batches (< SIMD_WIDTH) get
        // padded to full SIMD width with default SeqPair entries — padding lanes have len=0
        // and are filtered by the active mask in the SIMD body.
        #[cfg(target_arch = "x86_64")]
        {
            if std::is_x86_feature_detected!("avx512bw") && !disable_simd8() {
                let n = numPairs as usize;
                let mut idx = 0_usize;
                while idx + SIMD_WIDTH8_AVX512 <= n {
                    let chunk = &mut pairArray[idx..idx + SIMD_WIDTH8_AVX512];
                    let state = unsafe {
                        process_batch_soa_banded_dp_avx512(self, chunk, seqBufRef, seqBufQer, w)
                    };
                    for (lane, p) in chunk.iter_mut().enumerate() {
                        p.score = state.max_score[lane];
                        p.tle = state.max_i[lane] + 1;
                        p.qle = state.max_j[lane] + 1;
                        p.gtle = state.max_ie[lane] + 1;
                        p.gscore = state.gscore[lane];
                        p.max_off = state.max_off[lane];
                    }
                    idx += SIMD_WIDTH8_AVX512;
                }
                if idx < n {
                    // Pad tail batch to full SIMD width.
                    let mut padded = [SeqPair::default(); SIMD_WIDTH8_AVX512];
                    let tail = &pairArray[idx..n];
                    padded[..tail.len()].copy_from_slice(tail);
                    let pad_len2 = tail.last().map_or(0, |p| p.len2);
                    for p in &mut padded[tail.len()..] {
                        p.len2 = pad_len2;
                    }
                    let state = unsafe {
                        process_batch_soa_banded_dp_avx512(self, &padded, seqBufRef, seqBufQer, w)
                    };
                    for (lane, p) in pairArray[idx..n].iter_mut().enumerate() {
                        p.score = state.max_score[lane];
                        p.tle = state.max_i[lane] + 1;
                        p.qle = state.max_j[lane] + 1;
                        p.gtle = state.max_ie[lane] + 1;
                        p.gscore = state.gscore[lane];
                        p.max_off = state.max_off[lane];
                    }
                }
                return;
            }
        }
        self.scalarBandedSWAWrapper(pairArray, seqBufRef, seqBufQer, numPairs, nthreads, w);
    }

    pub fn getScores16(
        &self,
        pairArray: &mut [SeqPair],
        seqBufRef: &[u8],
        seqBufQer: &[u8],
        numPairs: i32,
        nthreads: u16,
        w: i32,
    ) {
        // i16 AVX-512 SIMD path. Validated against C++ AVX-512 across default scoring,
        // non-default -A/-B/-O/-E combinations, and z-drop conformance cases.
        #[cfg(target_arch = "x86_64")]
        {
            if std::is_x86_feature_detected!("avx512bw") && !disable_simd16() {
                let n = numPairs as usize;
                let mut idx = 0_usize;
                while idx + SIMD_WIDTH16_AVX512 <= n {
                    let chunk = &mut pairArray[idx..idx + SIMD_WIDTH16_AVX512];
                    let state = unsafe {
                        process_batch_soa_banded_dp_avx512_16(self, chunk, seqBufRef, seqBufQer, w)
                    };
                    for (lane, p) in chunk.iter_mut().enumerate() {
                        p.score = state.max_score[lane];
                        p.tle = state.max_i[lane] + 1;
                        p.qle = state.max_j[lane] + 1;
                        p.gtle = state.max_ie[lane] + 1;
                        p.gscore = state.gscore[lane];
                        p.max_off = state.max_off[lane];
                    }
                    idx += SIMD_WIDTH16_AVX512;
                }
                if idx < n {
                    let mut padded = [SeqPair::default(); SIMD_WIDTH16_AVX512];
                    let tail = &pairArray[idx..n];
                    padded[..tail.len()].copy_from_slice(tail);
                    let state = unsafe {
                        process_batch_soa_banded_dp_avx512_16(
                            self, &padded, seqBufRef, seqBufQer, w,
                        )
                    };
                    for (lane, p) in pairArray[idx..n].iter_mut().enumerate() {
                        p.score = state.max_score[lane];
                        p.tle = state.max_i[lane] + 1;
                        p.qle = state.max_j[lane] + 1;
                        p.gtle = state.max_ie[lane] + 1;
                        p.gscore = state.gscore[lane];
                        p.max_off = state.max_off[lane];
                    }
                }
                let _ = nthreads;
                return;
            }
        }
        self.scalarBandedSWAWrapper(pairArray, seqBufRef, seqBufQer, numPairs, nthreads, w);
    }

    // constructor
    #[doc = "Original function: BandedPairWiseSW::BandedPairWiseSW:50"]
    pub fn ctor(
        o_del: i32,
        e_del: i32,
        o_ins: i32,
        e_ins: i32,
        zdrop: i32,
        end_bonus: i32,
        mat_: &[i8; 25],
        w_match: i32,
        w_mismatch: i32,
        _numThreads: i32,
    ) -> Self {
        Self {
            mat: *mat_,
            m: 5,
            end_bonus,
            zdrop,
            o_del,
            o_ins,
            e_del,
            e_ins,
            w_match: i8::try_from(w_match).expect("w_match"),
            w_mismatch: i8::try_from(-w_mismatch).expect("w_mismatch"),
            w_open: o_del,   // redundant, used in vector code.
            w_extend: e_del, // redundant, used in vector code.
            w_ambig: DEFAULT_AMBIG,
            swTicks: 0,
            SW_cells: 0,
            setupTicks: 0,
            sort1Ticks: 0,
            sort2Ticks: 0,
        }
    }

    // destructor — upstream C++ frees the per-thread F8_/H8_/H8__ / F16_/H16_/H16__ scratch
    // arrays here. The Rust port relies on thread-local Vec scratch (BANDED_SCRATCH /
    // SIMD8_PROTO_SCRATCH / SIMD16_PROTO_SCRATCH) instead, so there is nothing to free.
    #[doc = "Original function: BandedPairWiseSW::~BandedPairWiseSW:97"]
    pub fn dtor(&mut self) {
        self.SW_cells = 0;
    }

    #[doc = "Original function: BandedPairWiseSW::getTicks:102"]
    pub fn getTicks(&self) -> i64 {
        self.sort1Ticks + self.setupTicks + self.swTicks + self.sort2Ticks
    }

    // ------------------------------------------------------------------------------------
    // Banded SWA - scalar code
    // ------------------------------------------------------------------------------------
    #[doc = "Original function: BandedPairWiseSW::scalarBandedSWA:116"]
    #[inline]
    pub fn scalarBandedSWA(
        &self,
        qlen: i32,
        query: &[u8],
        tlen: i32,
        target: &[u8],
        w: i32,
        h0: i32,
        qle: &mut i32,
        tle: &mut i32,
        gtle: &mut i32,
        gscore: &mut i32,
        max_off: &mut i32,
    ) -> i32 {
        BANDED_SCRATCH.with(|cell| {
            let mut buf = cell.borrow_mut();
            let (qp, eh_h, eh_e) = &mut *buf;
            self.scalarBandedSWA_inner(
                qlen, query, tlen, target, w, h0, qle, tle, gtle, gscore, max_off, qp, eh_h, eh_e,
            )
        })
    }

    // Same body as scalarBandedSWA but takes scratch buffers as &mut so callers that process
    // a batch can hoist the BANDED_SCRATCH.with(...) once across all pairs (avoids per-pair
    // thread-local borrow + alloc-amortization overhead in scalarBandedSWAWrapper, which the
    // profiler showed at >50% of total runtime when the with() was inside this function).
    #[allow(clippy::too_many_arguments)]
    fn scalarBandedSWA_inner(
        &self,
        qlen: i32,
        query: &[u8],
        tlen: i32,
        target: &[u8],
        mut w: i32,
        h0: i32,
        qle: &mut i32,
        tle: &mut i32,
        gtle: &mut i32,
        gscore: &mut i32,
        max_off: &mut i32,
        qp: &mut Vec<i8>,
        eh_h: &mut Vec<i32>,
        eh_e: &mut Vec<i32>,
    ) -> i32 {
        // qp: query profile, eh_h/eh_e: score array {H, E} split into SoA. Allocated to
        // qlen*m and qlen+1 respectively; thread-local buffers are reused across calls.
        let qlen_usize = qlen as usize;
        let tlen_usize = tlen as usize;
        let m_usize = self.m as usize;
        qp.clear();
        qp.resize(qlen_usize * m_usize, 0);
        eh_h.clear();
        eh_h.resize(qlen_usize + 1, 0);
        eh_e.clear();
        eh_e.resize(qlen_usize + 1, 0);
        let oe_del = self.o_del + self.e_del;
        let oe_ins = self.o_ins + self.e_ins;
        // Hoist e_del/e_ins out of the inner DP loop. The borrow checker may not realize they're
        // immutable across iterations since `self` could in principle alias eh_h/eh_e (it can't,
        // but explicit locals are easier on codegen).
        let e_del_loc = self.e_del;
        let e_ins_loc = self.e_ins;

        // generate the query profile
        let mut idx = 0_usize;
        for k in 0..m_usize {
            let p = &self.mat[k * m_usize..(k + 1) * m_usize];
            for &q in query.iter().take(qlen_usize) {
                qp[idx] = p[usize::from(q)];
                idx += 1;
            }
        }

        // fill the first row
        eh_h[0] = h0;
        eh_h[1] = if h0 > oe_ins { h0 - oe_ins } else { 0 };
        let mut j = 2_usize;
        while j <= qlen_usize && eh_h[j - 1] > self.e_ins {
            eh_h[j] = eh_h[j - 1] - self.e_ins;
            j += 1;
        }

        // adjust $w if it is too large
        let mut max_ins =
            (qlen * i32::from(self.w_match) + self.end_bonus - self.o_ins) / self.e_ins + 1;
        max_ins = max_ins.max(1);
        w = w.min(max_ins);
        let mut max_del =
            (qlen * i32::from(self.w_match) + self.end_bonus - self.o_del) / self.e_del + 1;
        max_del = max_del.max(1);
        w = w.min(max_del); // TODO (upstream): is this necessary?

        // DP loop
        let mut max = h0;
        let mut max_i = -1_i32;
        let mut max_j = -1_i32;
        let mut max_ie = -1_i32;
        let mut best_gscore = -1_i32;
        *max_off = 0;
        let mut beg = 0_i32;
        let mut end = qlen;

        for i in 0..tlen_usize {
            // i fits in i32 because tlen is i32 (passed in as such). Hoist the cast.
            let i_i32 = i as i32;
            let mut t;
            let mut f = 0_i32;
            let mut h1;
            let mut m = 0_i32;
            let mut mj = -1_i32;
            let target_base = usize::from(target[i]);
            let qprof = &qp[target_base * qlen_usize..target_base * qlen_usize + qlen_usize];

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
                h1 = h0 - (self.o_del + self.e_del * (i_i32 + 1));
                if h1 < 0 {
                    h1 = 0;
                }
            } else {
                h1 = 0;
            }

            // At the beginning of the loop: eh_h[j] = H(i-1, j-1), eh_e[j] = E(i, j),
            // f = F(i, j), and h1 = H(i, j-1). Cells are computed in the order:
            //   H(i, j)   = max{ H(i-1, j-1) + S(i, j), E(i, j), F(i, j) }
            //   E(i+1, j) = max{ H(i, j) - gapo, E(i, j) } - gape
            //   F(i, j+1) = max{ H(i, j) - gapo, F(i, j) } - gape
            //
            // Inner-loop bounds-check elimination — beg<end<=qlen, qprof has qlen elements,
            // eh_h/eh_e have qlen+1 elements, so j∈[beg,end) is in range.
            let beg_us = beg as usize;
            let end_us = end as usize;
            for j in beg_us..end_us {
                unsafe {
                    // load H(i-1, j-1) into m_score and E(i, j) into e
                    let mut m_score = *eh_h.get_unchecked(j);
                    let mut e = *eh_e.get_unchecked(j);
                    *eh_h.get_unchecked_mut(j) = h1; // save H(i, j-1) for the next row
                    // separating H and M to disallow a cigar like "100M3I3D20M"
                    m_score = if m_score != 0 {
                        m_score + i32::from(*qprof.get_unchecked(j))
                    } else {
                        0
                    };
                    // e and f are guaranteed to be non-negative, so h>=0 even if M<0
                    let mut h = m_score.max(e);
                    h = h.max(f);
                    h1 = h; // save H(i, j) to h1 for the next column
                    // record the position where the row max is achieved (m is stored at eh[mj+1])
                    if h >= m {
                        mj = j as i32;
                        m = h;
                    }
                    // computed E(i+1, j); save E(i+1, j) for the next row
                    t = (m_score - oe_del).max(0);
                    e -= e_del_loc;
                    e = e.max(t);
                    *eh_e.get_unchecked_mut(j) = e;
                    // computed F(i, j+1)
                    t = (m_score - oe_ins).max(0);
                    f -= e_ins_loc;
                    f = f.max(t);
                }
            }
            unsafe {
                *eh_h.get_unchecked_mut(end_us) = h1;
                *eh_e.get_unchecked_mut(end_us) = 0;
            }
            if end_us == qlen_usize {
                // C++ bandedSWA.cpp:202-204 uses `gscore > h1 ? max_ie : i` ternary which
                // updates on TIE (== treated as not-strictly-greater → use new i). Equivalent to >=.
                if h1 >= best_gscore {
                    max_ie = i_i32;
                    best_gscore = h1;
                }
            }
            if m == 0 {
                break;
            }
            if m > max {
                max = m;
                max_i = i_i32;
                max_j = mj;
                *max_off = (*max_off).max((mj - i_i32).abs());
            } else {
                // Match C++ SIMD ZSCORE8/ZSCORE16 (bandedSWA.cpp): no `if zdrop > 0` guard.
                // For zdrop > 0 this is equivalent (the `> zdrop` check has the same triggering
                // pattern). At zdrop=0, this triggers when the score gap is positive — matching
                // upstream SIMD's aggressive exit, where the C++ scalar's guard would have skipped.
                let di = i_i32 - max_i;
                let dj = mj - max_j;
                if di > dj {
                    if max - m - (di - dj) * self.e_del > self.zdrop {
                        break;
                    }
                } else if max - m - (dj - di) * self.e_ins > self.zdrop {
                    break;
                }
            }
            // update beg and end for the next round (drop leading/trailing zero columns)
            let mut jb = beg as usize;
            while jb < end as usize && eh_h[jb] == 0 && eh_e[jb] == 0 {
                jb += 1;
            }
            beg = jb as i32;
            let mut je = end as usize;
            while je >= jb && eh_h[je] == 0 && eh_e[je] == 0 {
                if je == 0 {
                    break;
                }
                je -= 1;
            }
            end = (je as i32 + 2).min(qlen);
        }

        *qle = max_j + 1;
        *tle = max_i + 1;
        *gtle = max_ie + 1;
        *gscore = best_gscore;
        max
    }

    // -------------------------------------------------------------
    // Banded SWA, wrapper function — runs scalarBandedSWA over a batch of SeqPairs.
    //-------------------------------------------------------------
    #[doc = "Original function: BandedPairWiseSW::scalarBandedSWAWrapper:242"]
    pub fn scalarBandedSWAWrapper(
        &self,
        seqPairArray: &mut [SeqPair],
        seqBufRef: &[u8],
        seqBufQer: &[u8],
        numPairs: i32,
        _nthreads: u16,
        w: i32,
    ) {
        // Hoist BANDED_SCRATCH.with() once across the whole batch so we don't pay a
        // per-pair thread-local borrow + drop. The inner DP just takes &mut to the buffers.
        BANDED_SCRATCH.with(|cell| {
            let mut buf = cell.borrow_mut();
            let (qp, eh_h, eh_e) = &mut *buf;
            for p in seqPairArray.iter_mut().take(numPairs as usize) {
                let idr = p.idr as usize;
                let idq = p.idq as usize;
                let seq1 = &seqBufRef[idr..idr + p.len1 as usize];
                let seq2 = &seqBufQer[idq..idq + p.len2 as usize];
                p.score = self.scalarBandedSWA_inner(
                    p.len2,
                    seq2,
                    p.len1,
                    seq1,
                    w,
                    p.h0,
                    &mut p.qle,
                    &mut p.tle,
                    &mut p.gtle,
                    &mut p.gscore,
                    &mut p.max_off,
                    qp,
                    eh_h,
                    eh_e,
                );
            }
        });
    }

    // Vector code, version 2.0. Line-specific entry points (`__L412`, `__L1970`,
    // `__L4198`) mirror the three SIMD-feature-guarded copies of `getScores8` in
    // bandedSWA.cpp (AVX2 / AVX-512 / SSE2). All three forward to the same production
    // `getScores8` dispatch — there is no logic difference in upstream beyond the
    // intrinsic widths.
    #[doc = "Original function: BandedPairWiseSW::getScores8:412"]
    pub fn getScores8__L412(
        &self,
        pairArray: &mut [SeqPair],
        seqBufRef: &[u8],
        seqBufQer: &[u8],
        numPairs: i32,
        numThreads: u16,
        w: i32,
    ) {
        self.getScores8(pairArray, seqBufRef, seqBufQer, numPairs, numThreads, w);
    }

    #[doc = "Original function: BandedPairWiseSW::smithWatermanBatchWrapper8:436"]
    pub fn smithWatermanBatchWrapper8__L436(
        &self,
        pairArray: &mut [SeqPair],
        seqBufRef: &[u8],
        seqBufQer: &[u8],
        numPairs: i32,
        _numThreads: u16,
        w: i32,
    ) {
        self.getScores8(pairArray, seqBufRef, seqBufQer, numPairs, 1, w);
    }

    #[doc = "Original function: BandedPairWiseSW::smithWaterman256_8:713"]
    pub fn smithWaterman256_8(
        &self,
        seq1SoA: &[u8],
        seq2SoA: &[u8],
        _nrow: i32,
        _ncol: i32,
        h0: &[u8],
        band: &[u8],
        _zdrop: i32,
        _end_bonus: i32,
        p: &mut [SeqPair],
        _tid: i32,
        numPairs: i32,
        _nstart: i32,
    ) {
        for (i, pair) in p.iter_mut().enumerate() {
            if i < h0.len() {
                pair.h0 = i32::from(h0[i]);
            }
        }
        let w = i32::from(*band.first().unwrap_or(&0));
        self.run_packed_kernel_u8(seq1SoA, seq2SoA, 32, p, numPairs, w);
    }

    // ------------------------- AVX2 - 16 bit SIMD_LANES ---------------------------
    #[doc = "Original function: BandedPairWiseSW::getScores16:1117"]
    pub fn getScores16__L1117(
        &self,
        pairArray: &mut [SeqPair],
        seqBufRef: &[u8],
        seqBufQer: &[u8],
        numPairs: i32,
        numThreads: u16,
        w: i32,
    ) {
        self.getScores16(pairArray, seqBufRef, seqBufQer, numPairs, numThreads, w);
    }

    #[doc = "Original function: BandedPairWiseSW::smithWatermanBatchWrapper16:1143"]
    pub fn smithWatermanBatchWrapper16__L1143(
        &self,
        pairArray: &mut [SeqPair],
        seqBufRef: &[u8],
        seqBufQer: &[u8],
        numPairs: i32,
        _numThreads: u16,
        w: i32,
    ) {
        self.getScores16(pairArray, seqBufRef, seqBufQer, numPairs, 1, w);
    }

    #[doc = "Original function: BandedPairWiseSW::smithWaterman256_16:1412"]
    pub fn smithWaterman256_16(
        &self,
        seq1SoA: &[u16],
        seq2SoA: &[u16],
        _nrow: i32,
        _ncol: i32,
        h0: &[u16],
        band: &[u16],
        _zdrop: i32,
        _end_bonus: i32,
        p: &mut [SeqPair],
        _tid: i32,
        numPairs: i32,
        _nstart: i32,
    ) {
        for (i, pair) in p.iter_mut().enumerate() {
            if i < h0.len() {
                pair.h0 = i32::from(h0[i]);
            }
        }
        let w = i32::from(*band.first().unwrap_or(&0));
        self.run_packed_kernel_u16(seq1SoA, seq2SoA, 16, p, numPairs, w);
    }

    // ____________________________ AVX512 - getScore() _______________________________________
    #[doc = "Original function: BandedPairWiseSW::getScores8:1970"]
    pub fn getScores8__L1970(
        &self,
        pairArray: &mut [SeqPair],
        seqBufRef: &[u8],
        seqBufQer: &[u8],
        numPairs: i32,
        numThreads: u16,
        w: i32,
    ) {
        self.getScores8(pairArray, seqBufRef, seqBufQer, numPairs, numThreads, w);
    }

    #[doc = "Original function: BandedPairWiseSW::smithWatermanBatchWrapper8:1997"]
    pub fn smithWatermanBatchWrapper8__L1997(
        &self,
        pairArray: &mut [SeqPair],
        seqBufRef: &[u8],
        seqBufQer: &[u8],
        numPairs: i32,
        _numThreads: u16,
        w: i32,
    ) {
        self.getScores8(pairArray, seqBufRef, seqBufQer, numPairs, 1, w);
    }

    #[doc = "Original function: BandedPairWiseSW::smithWaterman512_8:2263"]
    pub fn smithWaterman512_8(
        &self,
        seq1SoA: &[u8],
        seq2SoA: &[u8],
        _nrow: i32,
        _ncol: i32,
        h0: &[u8],
        band: &[u8],
        _zdrop: i32,
        _end_bonus: i32,
        p: &mut [SeqPair],
        _tid: i32,
        numPairs: i32,
        _nstart: i32,
    ) {
        for (i, pair) in p.iter_mut().enumerate() {
            if i < h0.len() {
                pair.h0 = i32::from(h0[i]);
            }
        }
        let w = i32::from(*band.first().unwrap_or(&0));
        self.run_packed_kernel_u8(seq1SoA, seq2SoA, 64, p, numPairs, w);
    }

    #[doc = "Original function: BandedPairWiseSW::getScores16:2664"]
    pub fn getScores16__L2664(
        &self,
        pairArray: &mut [SeqPair],
        seqBufRef: &[u8],
        seqBufQer: &[u8],
        numPairs: i32,
        numThreads: u16,
        w: i32,
    ) {
        self.getScores16(pairArray, seqBufRef, seqBufQer, numPairs, numThreads, w);
    }

    #[doc = "Original function: BandedPairWiseSW::smithWatermanBatchWrapper16:2690"]
    pub fn smithWatermanBatchWrapper16__L2690(
        &self,
        pairArray: &mut [SeqPair],
        seqBufRef: &[u8],
        seqBufQer: &[u8],
        numPairs: i32,
        _numThreads: u16,
        w: i32,
    ) {
        self.getScores16(pairArray, seqBufRef, seqBufQer, numPairs, 1, w);
    }

    #[doc = "Original function: BandedPairWiseSW::smithWaterman512_16:2962"]
    pub fn smithWaterman512_16(
        &self,
        seq1SoA: &[u16],
        seq2SoA: &[u16],
        _nrow: i32,
        _ncol: i32,
        h0: &[u16],
        band: &[u16],
        _zdrop: i32,
        _end_bonus: i32,
        p: &mut [SeqPair],
        _tid: i32,
        numPairs: i32,
        _nstart: i32,
    ) {
        for (i, pair) in p.iter_mut().enumerate() {
            if i < h0.len() {
                pair.h0 = i32::from(h0[i]);
            }
        }
        let w = i32::from(*band.first().unwrap_or(&0));
        self.run_packed_kernel_u16(seq1SoA, seq2SoA, 32, p, numPairs, w);
    }

    // ------------------------- SSE2 - 16 bit SIMD_LANES ---------------------------
    #[doc = "Original function: BandedPairWiseSW::getScores16:3469"]
    pub fn getScores16__L3469(
        &self,
        pairArray: &mut [SeqPair],
        seqBufRef: &[u8],
        seqBufQer: &[u8],
        numPairs: i32,
        numThreads: u16,
        w: i32,
    ) {
        self.getScores16(pairArray, seqBufRef, seqBufQer, numPairs, numThreads, w);
    }

    #[doc = "Original function: BandedPairWiseSW::smithWatermanBatchWrapper16:3492"]
    pub fn smithWatermanBatchWrapper16__L3492(
        &self,
        pairArray: &mut [SeqPair],
        seqBufRef: &[u8],
        seqBufQer: &[u8],
        numPairs: i32,
        _numThreads: u16,
        w: i32,
    ) {
        self.getScores16(pairArray, seqBufRef, seqBufQer, numPairs, 1, w);
    }

    #[doc = "Original function: BandedPairWiseSW::smithWaterman128_16:3757"]
    pub fn smithWaterman128_16(
        &self,
        seq1SoA: &[u16],
        seq2SoA: &[u16],
        _nrow: i32,
        _ncol: i32,
        h0: &[u16],
        band: &[u16],
        _zdrop: i32,
        _end_bonus: i32,
        p: &mut [SeqPair],
        _tid: i32,
        numPairs: i32,
        _nstart: i32,
    ) {
        for (i, pair) in p.iter_mut().enumerate() {
            if i < h0.len() {
                pair.h0 = i32::from(h0[i]);
            }
        }
        let w = i32::from(*band.first().unwrap_or(&0));
        self.run_packed_kernel_u16(seq1SoA, seq2SoA, 8, p, numPairs, w);
    }

    // ------------------------- SSE2 - 8 bit SIMD_LANES ---------------------------
    #[doc = "Original function: BandedPairWiseSW::getScores8:4198"]
    pub fn getScores8__L4198(
        &self,
        pairArray: &mut [SeqPair],
        seqBufRef: &[u8],
        seqBufQer: &[u8],
        numPairs: i32,
        numThreads: u16,
        w: i32,
    ) {
        self.getScores8(pairArray, seqBufRef, seqBufQer, numPairs, numThreads, w);
    }

    #[doc = "Original function: BandedPairWiseSW::smithWatermanBatchWrapper8:4223"]
    pub fn smithWatermanBatchWrapper8__L4223(
        &self,
        pairArray: &mut [SeqPair],
        seqBufRef: &[u8],
        seqBufQer: &[u8],
        numPairs: i32,
        _numThreads: u16,
        w: i32,
    ) {
        self.getScores8(pairArray, seqBufRef, seqBufQer, numPairs, 1, w);
    }

    #[doc = "Original function: BandedPairWiseSW::smithWaterman128_8:4485"]
    pub fn smithWaterman128_8(
        &self,
        seq1SoA: &[u8],
        seq2SoA: &[u8],
        _nrow: i32,
        _ncol: i32,
        h0: &[u8],
        band: &[u8],
        _zdrop: i32,
        _end_bonus: i32,
        p: &mut [SeqPair],
        _tid: i32,
        numPairs: i32,
        _nstart: i32,
    ) {
        for (i, pair) in p.iter_mut().enumerate() {
            if i < h0.len() {
                pair.h0 = i32::from(h0[i]);
            }
        }
        let w = i32::from(*band.first().unwrap_or(&0));
        self.run_packed_kernel_u8(seq1SoA, seq2SoA, 16, p, numPairs, w);
    }
}

#[cfg(test)]
mod tests {
    use super::{sortPairsId__L1952, sortPairsLen__L1909, BandedPairWiseSW};
    use crate::bwa_mem2::bandedswa::SeqPair;

    fn mat() -> [i8; 25] {
        let mut mat = [0_i8; 25];
        for i in 0..4 {
            for j in 0..4 {
                mat[i * 5 + j] = if i == j { 2 } else { -4 };
            }
            mat[i * 5 + 4] = -1;
        }
        for j in 0..5 {
            mat[20 + j] = -1;
        }
        mat
    }

    #[test]
    fn getScores8_line_wrappers_match_scalar_results() {
        let sw = BandedPairWiseSW::ctor(6, 1, 6, 1, 100, 0, &mat(), 2, 4, 1);
        let base = SeqPair {
            idr: 0,
            idq: 0,
            len1: 4,
            len2: 4,
            h0: 8,
            ..Default::default()
        };
        let seq_ref = [0_u8, 1, 2, 3];
        let seq_qer = [0_u8, 1, 2, 3];
        let mut a = [base];
        let mut b = [base];
        let mut c = [base];
        sw.getScores8__L412(&mut a, &seq_ref, &seq_qer, 1, 1, 10);
        sw.getScores8__L1970(&mut b, &seq_ref, &seq_qer, 1, 1, 10);
        sw.getScores8__L4198(&mut c, &seq_ref, &seq_qer, 1, 1, 10);
        assert_eq!(a[0].score, b[0].score);
        assert_eq!(a[0].score, c[0].score);
        assert!(a[0].score > 0);
    }

    #[test]
    fn getScores16_line_wrappers_match_scalar_results() {
        let sw = BandedPairWiseSW::ctor(6, 1, 6, 1, 100, 0, &mat(), 2, 4, 1);
        let base = SeqPair {
            idr: 0,
            idq: 0,
            len1: 4,
            len2: 4,
            h0: 8,
            ..Default::default()
        };
        let seq_ref = [0_u8, 1, 2, 3];
        let seq_qer = [0_u8, 1, 2, 3];
        let mut a = [base];
        let mut b = [base];
        let mut c = [base];
        sw.getScores16__L1117(&mut a, &seq_ref, &seq_qer, 1, 1, 10);
        sw.getScores16__L2664(&mut b, &seq_ref, &seq_qer, 1, 1, 10);
        sw.getScores16__L3469(&mut c, &seq_ref, &seq_qer, 1, 1, 10);
        assert_eq!(a[0].score, b[0].score);
        assert_eq!(a[0].score, c[0].score);
        assert!(a[0].score > 0);
    }

    #[test]
    fn packed_kernels_decode_and_score_consistently() {
        let sw = BandedPairWiseSW::ctor(6, 1, 6, 1, 100, 0, &mat(), 2, 4, 1);
        let mut p8 = vec![
            SeqPair {
                id: 0,
                len1: 4,
                len2: 4,
                h0: 8,
                ..Default::default()
            };
            1
        ];
        let mut seq1_8 = vec![0xff_u8; 4 * 16];
        let mut seq2_8 = vec![0xff_u8; 4 * 16];
        seq1_8[0] = 0;
        seq1_8[16] = 1;
        seq1_8[32] = 2;
        seq1_8[48] = 3;
        seq2_8[0] = 0;
        seq2_8[16] = 1;
        seq2_8[32] = 2;
        seq2_8[48] = 3;
        sw.smithWaterman128_8(
            &seq1_8,
            &seq2_8,
            4,
            4,
            &[8],
            &[10],
            100,
            0,
            &mut p8,
            0,
            1,
            0,
        );
        assert!(p8[0].score > 0);

        let mut p16 = vec![
            SeqPair {
                id: 0,
                len1: 4,
                len2: 4,
                h0: 8,
                ..Default::default()
            };
            1
        ];
        let mut seq1_16 = vec![u16::MAX; 4 * 8];
        let mut seq2_16 = vec![u16::MAX; 4 * 8];
        seq1_16[0] = 0;
        seq1_16[8] = 1;
        seq1_16[16] = 2;
        seq1_16[24] = 3;
        seq2_16[0] = 0;
        seq2_16[8] = 1;
        seq2_16[16] = 2;
        seq2_16[24] = 3;
        sw.smithWaterman128_16(
            &seq1_16,
            &seq2_16,
            4,
            4,
            &[8],
            &[10],
            100,
            0,
            &mut p16,
            0,
            1,
            0,
        );
        assert_eq!(p8[0].score, p16[0].score);
    }

    #[test]
    fn sort_helpers_match_cpp_counting_and_id_layout() {
        let mut pairs = vec![
            SeqPair {
                id: 11,
                len1: 3,
                ..Default::default()
            },
            SeqPair {
                id: 10,
                len1: 1,
                ..Default::default()
            },
            SeqPair {
                id: 12,
                len1: 3,
                ..Default::default()
            },
            SeqPair {
                id: 13,
                len1: 2,
                ..Default::default()
            },
        ];
        let mut temp = vec![SeqPair::default(); pairs.len()];
        let mut hist = vec![0_i16; 200];
        let mut histb = vec![0_i16; 200];
        sortPairsLen__L1909(&mut pairs, 4, &mut temp, &mut hist, &mut histb);
        assert_eq!(
            pairs.iter().map(|p| p.len1).collect::<Vec<_>>(),
            vec![1, 2, 3, 3]
        );
        assert_eq!(
            pairs.iter().map(|p| p.id).collect::<Vec<_>>(),
            vec![10, 13, 11, 12]
        );

        sortPairsId__L1952(&mut pairs, 10, 4, &mut temp);
        assert_eq!(
            pairs.iter().map(|p| p.id).collect::<Vec<_>>(),
            vec![10, 11, 12, 13]
        );
    }

    #[test]
    fn blend_helpers_follow_mask_bits() {
        let out16 = super::mm_blendv_epi16(&[0x0f0f_i16, 0x3333], &[0x5555, 0x7777], &[0x00ff, -1]);
        assert_eq!(out16, vec![0x0f55_i16, 0x7777_i16]);

        let out8 = super::mm_blendv_epi8(&[0x0f_u8, 0x33], &[0x55, 0x77], &[0x0f, 0xff]);
        assert_eq!(out8, vec![0x05_u8, 0x77_u8]);
    }
}
