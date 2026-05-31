#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Port of `bwa-mem2/src/fmi_search.h` + `bwa-mem2/src/fmi_search.cpp`.

use crate::bwa_mem2::bwa::bseq1_t;
use crate::bwa_mem2::read_index_ele::{indexEle, BWA_IDX_ALL};
use crate::bwa_mem2::sais::sais_suffixes_i64_upstream_port_mapped;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};

// --- fmi_search.h ---

#[doc = "Original struct: SMEM (bwa-mem2/src/FMI_search.h)"]
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct SMEM {
    pub rid: u32,
    pub m: u32,
    pub n: u32,
    pub k: i64,
    pub l: i64,
    pub s: i64,
}

#[doc = "Original struct: checkpoint_occ_scalar (bwa-mem2/src/FMI_search.h)"]
#[derive(Debug, Default, Clone, Copy)]
pub struct checkpoint_occ_scalar {
    pub _opaque: (),
}

#[doc = "Original struct: smem_struct (bwa-mem2/src/FMI_search.h)"]
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct smem_struct {
    pub rid: u32,
    pub m: u32,
    pub n: u32,
    pub k: i64,
    pub l: i64,
    pub s: i64,
}

#[doc = "Original function: _mm_countbits_64:61"]
pub(crate) fn mm_countbits_64(_arg0: crate::support::Opaque) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("_mm_countbits_64")
}

// --- fmi_search.cpp ---

thread_local! {
    // Reused across all calls to getSMEMsOnePosOneThread on a given Rayon worker.
    // C++ allocates SMEM prev[max_readlength] on the stack per call (kept tiny because the
    // compiler bounds maxLen). Rust has no equivalent stack-vec; per-call vec! adds significant
    // alloc overhead in the seed-finding inner loop. Reuse via thread-local Vec.
    static SMEM_PREV_SCRATCH: std::cell::RefCell<Vec<SMEM>> =
        const { std::cell::RefCell::new(Vec::new()) };

    // Reused across getSMEMsAllPosOneThread calls. The original C++ allocates this temporary with
    // _mm_malloc/_mm_free per call; in Rust that showed up as allocator work in the SMEM hotspot.
    static SMEM_QUERY_POS_SCRATCH: std::cell::RefCell<Vec<i16>> =
        const { std::cell::RefCell::new(Vec::new()) };

    // Reused across get_sa_entries_prefetch calls (per-chain in mem_chain_seeds).
    // Holds (pos_ar, map_ar) for the prefetch-bookkeeping pass.
    static SA_PREFETCH_SCRATCH: std::cell::RefCell<(Vec<i64>, Vec<usize>)> =
        const { std::cell::RefCell::new((Vec::new(), Vec::new())) };
}

const DUMMY_CHAR: u8 = 6;
const CP_BLOCK_SIZE: usize = 64;
const CP_FILENAME_SUFFIX: &str = ".bwt.2bit.64";
const CP_SHIFT: usize = 6;
const CP_MASK: usize = 63;
const SA_COMPX: usize = 3;
const SA_COMPX_MASK: usize = 0x7;

#[repr(C)]
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct CP_OCC {
    pub cp_count: [i64; 4],
    pub one_hot_bwt_str: [u64; 4],
}

#[derive(Debug, Default)]
pub struct FMI_search {
    pub base: indexEle,
    pub file_name: String,
    pub reference_seq_len: i64,
    pub sentinel_index: i64,
    pub count: [i64; 5],
    pub sa_ls_word: Vec<u32>,
    pub sa_ms_byte: Vec<i8>,
    pub cp_occ: Vec<CP_OCC>,
    pub one_hot_mask_array: Vec<u64>,
}

fn read_i64<R: Read>(reader: &mut R) -> i64 {
    let mut buf = [0_u8; 8];
    reader.read_exact(&mut buf).expect("read i64");
    i64::from_le_bytes(buf)
}

fn read_u64<R: Read>(reader: &mut R) -> u64 {
    let mut buf = [0_u8; 8];
    reader.read_exact(&mut buf).expect("read u64");
    u64::from_le_bytes(buf)
}

fn read_u32<R: Read>(reader: &mut R) -> u32 {
    let mut buf = [0_u8; 4];
    reader.read_exact(&mut buf).expect("read u32");
    u32::from_le_bytes(buf)
}

fn write_i8_slice<W: Write>(writer: &mut W, values: &[i8]) {
    let bytes = unsafe { std::slice::from_raw_parts(values.as_ptr().cast::<u8>(), values.len()) };
    writer.write_all(bytes).expect("write i8 slice");
}

#[cfg(target_endian = "little")]
fn write_u32_slice_le<W: Write>(writer: &mut W, values: &[u32]) {
    let bytes = unsafe {
        std::slice::from_raw_parts(values.as_ptr().cast::<u8>(), std::mem::size_of_val(values))
    };
    writer.write_all(bytes).expect("write u32 slice");
}

#[cfg(not(target_endian = "little"))]
fn write_u32_slice_le<W: Write>(writer: &mut W, values: &[u32]) {
    for value in values {
        writer
            .write_all(&value.to_le_bytes())
            .expect("write u32 slice");
    }
}

#[cfg(target_endian = "little")]
fn write_cp_occ_slice_le<W: Write>(writer: &mut W, values: &[CP_OCC]) {
    let bytes = unsafe {
        std::slice::from_raw_parts(values.as_ptr().cast::<u8>(), std::mem::size_of_val(values))
    };
    writer.write_all(bytes).expect("write cp_occ slice");
}

#[cfg(not(target_endian = "little"))]
fn write_cp_occ_slice_le<W: Write>(writer: &mut W, values: &[CP_OCC]) {
    for cpo in values {
        for value in cpo.cp_count {
            writer
                .write_all(&value.to_le_bytes())
                .expect("write cp_count");
        }
        for value in cpo.one_hot_bwt_str {
            writer
                .write_all(&value.to_le_bytes())
                .expect("write one_hot");
        }
    }
}

#[inline(always)]
fn countbits_64(x: u64) -> i32 {
    x.count_ones() as i32
}

#[inline]
fn decode_bwt_base(cp_occ: &[CP_OCC], sp: i64) -> u8 {
    let occ_id = (sp >> CP_SHIFT) as usize;
    let y = CP_BLOCK_SIZE - ((sp & CP_MASK as i64) as usize) - 1;
    debug_assert!(occ_id < cp_occ.len());
    let one_hot_bwt_str = unsafe { &cp_occ.get_unchecked(occ_id).one_hot_bwt_str };
    // Branchless: pack the y-th bit of each lane into a 4-bit nibble, then trailing_zeros gives
    // the base index (0..3). Sentinel value 4 when no lane bit is set. Avoids the chain of 4
    // mispredict-prone branches in the if/else-if version.
    let bit0 = ((one_hot_bwt_str[0] >> y) & 1) as u32;
    let bit1 = ((one_hot_bwt_str[1] >> y) & 1) as u32;
    let bit2 = ((one_hot_bwt_str[2] >> y) & 1) as u32;
    let bit3 = ((one_hot_bwt_str[3] >> y) & 1) as u32;
    let packed = bit0 | (bit1 << 1) | (bit2 << 2) | (bit3 << 3);
    if packed == 0 {
        4
    } else {
        packed.trailing_zeros() as u8
    }
}

#[inline(always)]
fn occ(cp_occ: &[CP_OCC], one_hot_mask_array: &[u64], sp: i64, b: u8) -> i64 {
    // sp is a BWT position bounded by reference_seq_len; cp_occ is sized to fit all
    // (sp >> CP_SHIFT) values, and y_sp is masked to CP_MASK so it indexes one_hot_mask_array
    // (which is exactly 64 elements). b ∈ 0..4 indexes cp_count (size 4) and one_hot_bwt_str.
    let occ_id_sp = (sp >> CP_SHIFT) as usize;
    let y_sp = (sp & CP_MASK as i64) as usize;
    let b_us = usize::from(b);
    unsafe {
        let entry = cp_occ.get_unchecked(occ_id_sp);
        let mut occ_sp = *entry.cp_count.get_unchecked(b_us);
        let one_hot_bwt_str_c_sp = *entry.one_hot_bwt_str.get_unchecked(b_us);
        let match_mask_sp = one_hot_bwt_str_c_sp & *one_hot_mask_array.get_unchecked(y_sp);
        occ_sp += i64::from(countbits_64(match_mask_sp));
        occ_sp
    }
}

/// `qsort`-style three-way comparator on `SMEM` records.
///
/// Orders primarily by ascending `rid`, then ascending `m` (left edge), then
/// descending `n` (right edge). Returns `-1`, `0`, or `1` like the C++
/// `compare_smem`. Used by `sortSMEMs` to stabilise per-thread SMEM slices.
#[doc = "Original function: compare_smem:987"]
#[inline]
pub fn compare_smem(a: &SMEM, b: &SMEM) -> i32 {
    match compare_smem_ord(a, b) {
        std::cmp::Ordering::Less => -1,
        std::cmp::Ordering::Greater => 1,
        std::cmp::Ordering::Equal => 0,
    }
}

#[inline]
pub fn compare_smem_ord(a: &SMEM, b: &SMEM) -> std::cmp::Ordering {
    a.rid
        .cmp(&b.rid)
        .then_with(|| a.m.cmp(&b.m))
        .then_with(|| b.n.cmp(&a.n))
}

impl FMI_search {
    /// Construct an `FMI_search` whose on-disk artefacts share the prefix `fname`.
    ///
    /// All FM-index buffers start empty; populate them via [`load_index`] for
    /// reads or [`build_index`] for index construction.
    #[doc = "Original function: FMI_search::FMI_search:44"]
    pub fn ctor(fname: &str) -> Self {
        Self {
            base: indexEle::ctor(),
            file_name: fname.to_string(),
            reference_seq_len: 0,
            sentinel_index: 0,
            count: [0; 5],
            sa_ls_word: Vec::new(),
            sa_ms_byte: Vec::new(),
            cp_occ: Vec::new(),
            one_hot_mask_array: Vec::new(),
        }
    }

    /// Release every owned FM-index buffer and the underlying `indexEle`.
    #[doc = "Original function: FMI_search::~FMI_search:58"]
    pub fn dtor(&mut self) {
        self.sa_ms_byte.clear();
        self.sa_ls_word.clear();
        self.cp_occ.clear();
        self.one_hot_mask_array.clear();
        self.base.dtor();
    }

    /// Return the unpacked length of the reference encoded in `fn_pac`.
    ///
    /// The `.pac` file stores 2-bit-packed bases; the final byte records the
    /// trailing-base count (0..3) so that the unpacked length is
    /// `(file_size - 1) * 4 + tail_byte`.
    ///
    /// # Arguments
    /// * `fn_pac` - path to the `.pac` file.
    #[doc = "Original function: FMI_search::pac_seq_len:70"]
    pub fn pac_seq_len(&self, fn_pac: &str) -> i64 {
        let mut fp =
            File::open(fn_pac).unwrap_or_else(|e| panic!("fail to open file '{fn_pac}' : {e}"));
        fp.seek(SeekFrom::End(-1)).expect("seek pac tail");
        let pac_len = fp.stream_position().expect("ftell pac") as i64;
        let mut c = [0_u8; 1];
        fp.read_exact(&mut c).expect("read pac tail");
        (pac_len - 1) * 4 + i64::from(c[0])
    }

    /// Decode the 2-bit-packed reference at `fn_pac` and append its forward
    /// strand followed by the reverse-complement strand to `reference_seq`.
    ///
    /// Used by index construction to materialise the concatenated
    /// forward+reverse reference that the suffix-array builder consumes.
    ///
    /// # Arguments
    /// * `fn_pac` - path to the `.pac` file.
    /// * `reference_seq` - destination string; the function appends ASCII
    ///   bases (`A`/`C`/`G`/`T`) to it without clearing existing contents.
    #[doc = "Original function: FMI_search::pac2nt:83"]
    pub fn pac2nt(&self, fn_pac: &str, reference_seq: &mut String) {
        let seq_len = self.pac_seq_len(fn_pac);
        assert!(seq_len > 0);
        let mut fp =
            File::open(fn_pac).unwrap_or_else(|e| panic!("fail to open file '{fn_pac}' : {e}"));
        let pac_size = (seq_len >> 2) + if (seq_len & 3) == 0 { 0 } else { 1 };
        let mut buf = vec![0_u8; usize::try_from(pac_size).expect("pac size")];
        fp.read_exact(&mut buf).expect("read pac");

        for i in 0..seq_len {
            let nt = (buf[usize::try_from(i >> 2).expect("idx")] >> ((3 - (i & 3)) << 1)) & 3;
            reference_seq.push(match nt {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                3 => 'T',
                _ => panic!("ERROR! Value of nt is not in 0,1,2,3!"),
            });
        }

        for i in (0..usize::try_from(seq_len).expect("seq len")).rev() {
            let c = reference_seq.as_bytes()[i] as char;
            reference_seq.push(match c {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                'T' => 'A',
                _ => unreachable!(),
            });
        }
    }

    /// Write the checkpointed FM-index (`<ref>.bwt.2bit.64`) from an
    /// already-computed `i64` suffix array.
    ///
    /// Materialises the BWT, builds the per-`CP_BLOCK_SIZE` checkpointed
    /// occurrence table, downsamples the suffix array by `SA_COMPX`, and
    /// updates `self` with the same in-memory copies for subsequent queries.
    ///
    /// # Arguments
    /// * `ref_file_name` - index prefix; the output is `<ref_file_name>.bwt.2bit.64`.
    /// * `binary_seq` - reference encoded as 0..3 base codes.
    /// * `ref_seq_len` - length of `binary_seq` excluding the appended `$` sentinel.
    /// * `sa_bwt` - suffix array of length `ref_seq_len + 1`.
    /// * `count` - cumulative base counts (size 5) used to derive `self.count`.
    ///
    /// # Returns
    /// Always `0` on success; matches the C++ `int` return convention.
    #[doc = "Original function: FMI_search::build_fm_index:144"]
    pub fn build_fm_index(
        &mut self,
        ref_file_name: &str,
        binary_seq: &[u8],
        ref_seq_len: i64,
        sa_bwt: &[i64],
        count: &[i64; 5],
    ) -> i32 {
        let outname = format!("{ref_file_name}{CP_FILENAME_SUFFIX}");
        let mut outstream = BufWriter::with_capacity(
            1 << 20,
            File::create(&outname)
                .unwrap_or_else(|e| panic!("fail to open file '{outname}' : {e}")),
        );

        let ref_seq_len_with_sentinel = ref_seq_len + 1;
        outstream
            .write_all(&ref_seq_len_with_sentinel.to_le_bytes())
            .expect("write ref len");
        for value in count {
            outstream
                .write_all(&value.to_le_bytes())
                .expect("write count");
        }

        let ref_seq_len_aligned =
            (((usize::try_from(ref_seq_len_with_sentinel).expect("len") + CP_BLOCK_SIZE - 1)
                / CP_BLOCK_SIZE)
                * CP_BLOCK_SIZE) as i64;
        let mut bwt = vec![DUMMY_CHAR; usize::try_from(ref_seq_len_aligned).expect("aligned len")];
        let mut sentinel_index = -1_i64;
        for i in 0..ref_seq_len_with_sentinel {
            let sai = sa_bwt[usize::try_from(i).expect("sa idx")];
            if sai == 0 {
                bwt[usize::try_from(i).expect("bwt idx")] = 4;
                sentinel_index = i;
            } else {
                let c = binary_seq[usize::try_from(sai - 1).expect("binary idx")];
                assert!(c < 4, "ERROR! i = {i}, c = {c}");
                bwt[usize::try_from(i).expect("bwt idx")] = c;
            }
        }

        // create checkpointed occ
        let cp_occ_size =
            (usize::try_from(ref_seq_len_with_sentinel).expect("cp len") >> CP_SHIFT) + 1;
        let mut cp_occ = vec![CP_OCC::default(); cp_occ_size];
        let mut cp_count = [0_i64; 16];
        for i in 0..usize::try_from(ref_seq_len_with_sentinel).expect("loop len") {
            if (i & CP_MASK) == 0 {
                let mut cpo = CP_OCC::default();
                cpo.cp_count.copy_from_slice(&cp_count[..4]);
                for j in 0..CP_BLOCK_SIZE {
                    for slot in &mut cpo.one_hot_bwt_str {
                        *slot <<= 1;
                    }
                    let c = bwt[i + j];
                    if c < 4 {
                        cpo.one_hot_bwt_str[usize::from(c)] += 1;
                    }
                }
                cp_occ[i >> CP_SHIFT] = cpo;
            }
            cp_count[usize::from(bwt[i])] += 1;
        }
        write_cp_occ_slice_le(&mut outstream, &cp_occ);
        self.cp_occ = cp_occ.clone();

        let sa_len = (usize::try_from(ref_seq_len_with_sentinel).expect("sa len") >> SA_COMPX) + 1;
        let mut sa_ls_word = vec![0_u32; sa_len];
        let mut sa_ms_byte = vec![0_i8; sa_len];
        let mut pos = 0_usize;
        for (i, value) in sa_bwt.iter().enumerate() {
            if (i & SA_COMPX_MASK) == 0 {
                sa_ls_word[pos] = (*value & 0xffff_ffff) as u32;
                sa_ms_byte[pos] = ((*value >> 32) & 0xff) as i8;
                pos += 1;
            }
        }
        write_i8_slice(&mut outstream, &sa_ms_byte);
        write_u32_slice_le(&mut outstream, &sa_ls_word);
        outstream
            .write_all(&sentinel_index.to_le_bytes())
            .expect("write sentinel");

        self.reference_seq_len = ref_seq_len_with_sentinel;
        self.sentinel_index = sentinel_index;
        self.count = *count;
        self.sa_ms_byte = sa_ms_byte;
        self.sa_ls_word = sa_ls_word;
        0
    }

    fn build_fm_index_from_suffixes<T>(
        &mut self,
        ref_file_name: &str,
        binary_seq: &[u8],
        ref_seq_len: i64,
        suffixes: Vec<T>,
        count: &[i64; 5],
    ) -> i32
    where
        T: Copy + Into<i64>,
    {
        let outname = format!("{ref_file_name}{CP_FILENAME_SUFFIX}");
        let mut outstream = BufWriter::with_capacity(
            1 << 20,
            File::create(&outname)
                .unwrap_or_else(|e| panic!("fail to open file '{outname}' : {e}")),
        );

        let ref_seq_len_with_sentinel = ref_seq_len + 1;
        outstream
            .write_all(&ref_seq_len_with_sentinel.to_le_bytes())
            .expect("write ref len");
        for value in count {
            outstream
                .write_all(&value.to_le_bytes())
                .expect("write count");
        }

        let ref_seq_len_aligned =
            (((usize::try_from(ref_seq_len_with_sentinel).expect("len") + CP_BLOCK_SIZE - 1)
                / CP_BLOCK_SIZE)
                * CP_BLOCK_SIZE) as i64;
        let ref_seq_len_with_sentinel_usize =
            usize::try_from(ref_seq_len_with_sentinel).expect("sa len");
        let mut bwt = vec![DUMMY_CHAR; usize::try_from(ref_seq_len_aligned).expect("aligned len")];
        let sa_len = (ref_seq_len_with_sentinel_usize >> SA_COMPX) + 1;
        let mut sa_ls_word = vec![0_u32; sa_len];
        let mut sa_ms_byte = vec![0_i8; sa_len];
        let mut sentinel_index = -1_i64;
        let mut pos = 0_usize;
        for i in 0..ref_seq_len_with_sentinel_usize {
            let sai = if i == 0 {
                ref_seq_len
            } else {
                suffixes[i - 1].into()
            };
            if (i & SA_COMPX_MASK) == 0 {
                sa_ls_word[pos] = (sai & 0xffff_ffff) as u32;
                sa_ms_byte[pos] = ((sai >> 32) & 0xff) as i8;
                pos += 1;
            };
            if sai == 0 {
                bwt[i] = 4;
                sentinel_index = i as i64;
            } else {
                let c = binary_seq[usize::try_from(sai - 1).expect("binary idx")];
                assert!(c < 4, "ERROR! i = {i}, c = {c}");
                bwt[i] = c;
            }
        }
        drop(suffixes);

        // create checkpointed occ
        let cp_occ_size =
            (usize::try_from(ref_seq_len_with_sentinel).expect("cp len") >> CP_SHIFT) + 1;
        let mut cp_occ = vec![CP_OCC::default(); cp_occ_size];
        let mut cp_count = [0_i64; 16];
        for i in 0..usize::try_from(ref_seq_len_with_sentinel).expect("loop len") {
            if (i & CP_MASK) == 0 {
                let mut cpo = CP_OCC::default();
                cpo.cp_count.copy_from_slice(&cp_count[..4]);
                for j in 0..CP_BLOCK_SIZE {
                    for slot in &mut cpo.one_hot_bwt_str {
                        *slot <<= 1;
                    }
                    let c = bwt[i + j];
                    if c < 4 {
                        cpo.one_hot_bwt_str[usize::from(c)] += 1;
                    }
                }
                cp_occ[i >> CP_SHIFT] = cpo;
            }
            cp_count[usize::from(bwt[i])] += 1;
        }
        write_cp_occ_slice_le(&mut outstream, &cp_occ);
        self.cp_occ = cp_occ;

        write_i8_slice(&mut outstream, &sa_ms_byte);
        write_u32_slice_le(&mut outstream, &sa_ls_word);
        outstream
            .write_all(&sentinel_index.to_le_bytes())
            .expect("write sentinel");

        self.reference_seq_len = ref_seq_len_with_sentinel;
        self.sentinel_index = sentinel_index;
        self.count = *count;
        self.sa_ms_byte = sa_ms_byte;
        self.sa_ls_word = sa_ls_word;
        0
    }

    /// Build the full FM-index for the `.pac` reference at `self.file_name`.
    ///
    /// Decodes the packed reference via [`pac2nt`], 2-bit-encodes it to a
    /// `<prefix>.0123` file, computes the cumulative base counts, runs SA-IS
    /// to obtain the suffix array (falling back to a naive comparator-sort
    /// when the reference exceeds `i32::MAX` bases), and delegates the
    /// checkpointed-FM emit to [`build_fm_index`]-internal helpers.
    ///
    /// # Returns
    /// Always `0` on success.
    #[doc = "Original function: FMI_search::build_index:306"]
    pub fn build_index(&mut self) -> i32 {
        let prefix = self.file_name.clone();
        let pac_file_name = format!("{prefix}.pac");
        let mut reference_seq = String::new();
        self.pac2nt(&pac_file_name, &mut reference_seq);
        let pac_len = i64::try_from(reference_seq.len()).expect("reference seq len");

        let mut binary_ref_seq = Vec::with_capacity(usize::try_from(pac_len).expect("cap"));
        let mut count = [0_i64; 16];
        for c in reference_seq.bytes() {
            let v = match c {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 4,
            };
            binary_ref_seq.push(v);
            if v < 4 {
                count[usize::from(v)] += 1;
            }
        }
        count[4] = count[0] + count[1] + count[2] + count[3];
        count[3] = count[0] + count[1] + count[2];
        count[2] = count[0] + count[1];
        count[1] = count[0];
        count[0] = 0;

        let binary_ref_name = format!("{prefix}.0123");
        std::fs::write(&binary_ref_name, &binary_ref_seq).expect("write binary ref");

        let counts5 = [count[0], count[1], count[2], count[3], count[4]];
        let suffixes = sais_suffixes_i64_upstream_port_mapped(reference_seq.as_bytes());
        drop(reference_seq);
        self.build_fm_index_from_suffixes(&prefix, &binary_ref_seq, pac_len, suffixes, &counts5)
    }

    /// Load the checkpointed FM-index from `<self.file_name>.bwt.2bit.64`
    /// into in-memory buffers, then load the companion `.bns`/`.pac` files
    /// via the embedded `indexEle` loader.
    ///
    /// Initialises `one_hot_mask_array`, the per-`CP_BLOCK_SIZE` `cp_occ`
    /// table, the downsampled suffix array (`sa_ms_byte` / `sa_ls_word`),
    /// and the sentinel index. Panics on any I/O failure.
    #[doc = "Original function: FMI_search::load_index:384"]
    pub fn load_index(&mut self) {
        self.one_hot_mask_array = vec![0_u64; 64];
        self.one_hot_mask_array[0] = 0;
        let base = 0x8000_0000_0000_0000_u64;
        self.one_hot_mask_array[1] = base;
        for i in 2..64 {
            self.one_hot_mask_array[i] = (self.one_hot_mask_array[i - 1] >> 1) | base;
        }

        let ref_file_name = self.file_name.clone();
        let cp_file_name = format!("{ref_file_name}{CP_FILENAME_SUFFIX}");
        // Read the BWT and FM index of the reference sequence.
        // 64-byte read_exact calls in the cp_occ loop become one syscall per call without buffering
        // — wrap in BufReader to batch them. C++ uses fread (buffered) so this matches semantics.
        let mut cpstream = BufReader::with_capacity(
            1 << 20,
            File::open(&cp_file_name)
                .unwrap_or_else(|e| panic!("ERROR! Unable to open the file: {cp_file_name}: {e}")),
        );

        eprintln!("* Index file found. Loading index from {cp_file_name}");
        self.reference_seq_len = read_i64(&mut cpstream);
        assert!(self.reference_seq_len > 0);
        assert!(self.reference_seq_len <= 0x7fff_ffff_ff);
        eprintln!(
            "* Reference seq len for bi-index = {}",
            self.reference_seq_len
        );

        // create checkpointed occ
        let cp_occ_size =
            usize::try_from((self.reference_seq_len >> CP_SHIFT) + 1).expect("cp occ size");
        for slot in &mut self.count {
            *slot = read_i64(&mut cpstream);
        }

        self.cp_occ = Vec::with_capacity(cp_occ_size);
        for _ in 0..cp_occ_size {
            let mut cp = CP_OCC::default();
            for value in &mut cp.cp_count {
                *value = read_i64(&mut cpstream);
            }
            for value in &mut cp.one_hot_bwt_str {
                *value = read_u64(&mut cpstream);
            }
            self.cp_occ.push(cp);
        }
        // update read count structure
        for value in &mut self.count {
            *value += 1;
        }

        let reference_seq_len_ =
            usize::try_from((self.reference_seq_len >> SA_COMPX) + 1).expect("compressed sa len");
        let mut sa_ms_raw = vec![0_u8; reference_seq_len_];
        cpstream.read_exact(&mut sa_ms_raw).expect("read sa ms");
        self.sa_ms_byte = sa_ms_raw.into_iter().map(|v| v as i8).collect();

        self.sa_ls_word = Vec::with_capacity(reference_seq_len_);
        for _ in 0..reference_seq_len_ {
            self.sa_ls_word.push(read_u32(&mut cpstream));
        }

        self.sentinel_index = read_i64(&mut cpstream);
        eprintln!("* sentinel-index: {}", self.sentinel_index);

        eprintln!("* Count:");
        for (i, count) in self.count.iter().enumerate() {
            eprintln!("{i},\t{count}");
        }
        eprintln!();
        eprintln!(
            "* Reading other elements of the index from files {}",
            ref_file_name
        );
        self.base.bwa_idx_load_ele(&ref_file_name, BWA_IDX_ALL);
        eprintln!("* Done reading Index!!");
    }

    /// Advance one SMEM-search step per active read on the calling thread.
    ///
    /// For each read referenced by `rid_array[i]`, starts at
    /// `query_pos_array[i]`, runs the forward extension followed by the
    /// backward extension, appends qualifying SMEMs to `match_array`, and
    /// writes the new position to `query_pos_array[i]` so that
    /// [`getSMEMsAllPosOneThread`] can iterate until all reads are exhausted.
    ///
    /// # Arguments
    /// * `enc_qdb` - flat 0..3 encoded read buffer indexed via `query_cum_len_ar`.
    /// * `query_pos_array` - current search head per active read; updated in-place.
    /// * `min_intv_array` - minimum BWT interval size per active read.
    /// * `rid_array` - read index per active slot.
    /// * `numReads` - number of currently active slots (length cap on the arrays).
    /// * `seq_` - read records, used for sequence lengths.
    /// * `query_cum_len_ar` - cumulative offsets of each read into `enc_qdb`.
    /// * `max_readlength` - upper bound on read length; sizes the prev-array scratch.
    /// * `minSeedLen` - minimum seed length to emit.
    /// * `match_array` - destination SMEMs; new entries are pushed.
    /// * `num_total_smem` - running total of emitted SMEMs; incremented in-place.
    #[doc = "Original function: FMI_search::getSMEMsOnePosOneThread:496"]
    pub fn getSMEMsOnePosOneThread(
        &self,
        enc_qdb: &[u8],
        query_pos_array: &mut [i16],
        min_intv_array: &mut [i32],
        rid_array: &mut [i32],
        numReads: i32,
        _batch_size: i32,
        seq_: &[bseq1_t],
        query_cum_len_ar: &[i32],
        max_readlength: i32,
        minSeedLen: i32,
        match_array: &mut Vec<SMEM>,
        num_total_smem: &mut i64,
    ) {
        let mut total = *num_total_smem;
        SMEM_PREV_SCRATCH.with(|cell| {
            let mut prev_array = cell.borrow_mut();
            let need = max_readlength.max(0) as usize;
            if prev_array.len() < need {
                prev_array.resize(need, SMEM::default());
            }
            self.getSMEMsOnePosOneThread_inner(
                enc_qdb,
                query_pos_array,
                min_intv_array,
                rid_array,
                numReads,
                _batch_size,
                seq_,
                query_cum_len_ar,
                max_readlength,
                minSeedLen,
                match_array,
                num_total_smem,
                &mut total,
                &mut prev_array,
            );
        });
        *num_total_smem = total;
    }

    #[allow(clippy::too_many_arguments)]
    fn getSMEMsOnePosOneThread_inner(
        &self,
        enc_qdb: &[u8],
        query_pos_array: &mut [i16],
        min_intv_array: &mut [i32],
        rid_array: &mut [i32],
        numReads: i32,
        _batch_size: i32,
        seq_: &[bseq1_t],
        query_cum_len_ar: &[i32],
        max_readlength: i32,
        minSeedLen: i32,
        match_array: &mut Vec<SMEM>,
        num_total_smem: &mut i64,
        total: &mut i64,
        prev_array: &mut [SMEM],
    ) {
        let _ = max_readlength;
        let _ = num_total_smem;
        let limit = usize::min(
            numReads as usize,
            usize::min(
                query_pos_array.len(),
                usize::min(min_intv_array.len(), rid_array.len()),
            ),
        );

        let min_seed_len_u32 = minSeedLen as u32;
        // Perform SMEM for original reads
        for i in 0..limit {
            let x = query_pos_array[i] as usize;
            let rid = rid_array[i] as usize;
            let mut next_x = (x + 1) as i32;
            let readlength = seq_[rid].l_seq as usize;
            let offset = query_cum_len_ar[rid] as usize;
            let min_intv = i64::from(min_intv_array[i]);
            let mut a = enc_qdb[offset + x];

            if a < 4 {
                let mut smem = SMEM {
                    rid: rid as u32,
                    m: x as u32,
                    n: x as u32,
                    k: self.count[usize::from(a)],
                    l: self.count[usize::from(3 - a)],
                    s: self.count[usize::from(a + 1)] - self.count[usize::from(a)],
                };
                let mut num_prev = 0_usize;

                for j in (x + 1)..readlength {
                    a = enc_qdb[offset + j];
                    // j < readlength which fits in i32/u32 for typical sequencing reads.
                    next_x = (j + 1) as i32;
                    if a < 4 {
                        // Forward extension is backward extension with the BWT of reverse complement
                        let mut smem_ = smem;
                        smem_.k = smem.l;
                        smem_.l = smem.k;
                        let new_smem_ = self.backwardExt(smem_, 3 - a);
                        let mut new_smem = new_smem_;
                        new_smem.k = new_smem_.l;
                        new_smem.l = new_smem_.k;
                        new_smem.n = j as u32;

                        let s_neq_mask = if new_smem.s != smem.s {
                            1_usize
                        } else {
                            0_usize
                        };
                        prev_array[num_prev] = smem;
                        num_prev += s_neq_mask;
                        if new_smem.s < min_intv {
                            next_x = j as i32;
                            break;
                        }
                        smem = new_smem;
                        // Mirror upstream's ENABLE_PREFETCH: warm cp_occ for next-iter k/l reads.
                        // Next iteration will read cp_occ[smem.l>>CP_SHIFT] and
                        // cp_occ[(smem.l+smem.s)>>CP_SHIFT] (after the k/l swap inside backwardExt).
                        // _mm_prefetch on invalid addresses is a silent no-op (no fault) so the
                        // bounds checks are unnecessary; smem.l is already within valid BWT range.
                        #[cfg(target_arch = "x86_64")]
                        unsafe {
                            use core::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
                            let k_idx = (smem.l >> CP_SHIFT) as usize;
                            let l_idx = ((smem.l + smem.s) >> CP_SHIFT) as usize;
                            _mm_prefetch(self.cp_occ.as_ptr().add(k_idx) as *const i8, _MM_HINT_T0);
                            _mm_prefetch(self.cp_occ.as_ptr().add(l_idx) as *const i8, _MM_HINT_T0);
                        }
                    } else {
                        break;
                    }
                }

                if smem.s >= min_intv {
                    prev_array[num_prev] = smem;
                    num_prev += 1;
                }

                prev_array[..num_prev].reverse();

                // Backward search
                for j in (0..x).rev() {
                    let mut num_curr = 0_usize;
                    let mut curr_s = -1_i64;
                    a = enc_qdb[offset + j];
                    if a > 3 {
                        break;
                    }

                    let mut p = 0_usize;
                    let j_u32 = j as u32; // j < x < readlength fits in u32.
                                          // Each backwardExt call reads cp_occ[smem.k>>CP_SHIFT] and
                                          // cp_occ[(smem.k+smem.s)>>CP_SHIFT]. The next iter's smem is prev_array[p+1]
                                          // (or [0] when first while breaks into second), known up-front. Warming those
                                          // lines hides the DRAM latency that dominates the SMEM hot loop.
                    #[cfg(target_arch = "x86_64")]
                    {
                        use core::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
                        if num_prev > 0 {
                            let smem0 = prev_array[0];
                            let k_idx = (smem0.k >> CP_SHIFT) as usize;
                            let l_idx = ((smem0.k + smem0.s) >> CP_SHIFT) as usize;
                            unsafe {
                                _mm_prefetch(
                                    self.cp_occ.as_ptr().add(k_idx) as *const i8,
                                    _MM_HINT_T0,
                                );
                                _mm_prefetch(
                                    self.cp_occ.as_ptr().add(l_idx) as *const i8,
                                    _MM_HINT_T0,
                                );
                            }
                        }
                    }

                    while p < num_prev {
                        let smem = prev_array[p];
                        // Prefetch next iter's cp_occ (whichever path we take, p+1 is consumed
                        // either by this loop or the next).
                        #[cfg(target_arch = "x86_64")]
                        if p + 1 < num_prev {
                            use core::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
                            let next_smem = prev_array[p + 1];
                            let k_idx = (next_smem.k >> CP_SHIFT) as usize;
                            let l_idx = ((next_smem.k + next_smem.s) >> CP_SHIFT) as usize;
                            unsafe {
                                _mm_prefetch(
                                    self.cp_occ.as_ptr().add(k_idx) as *const i8,
                                    _MM_HINT_T0,
                                );
                                _mm_prefetch(
                                    self.cp_occ.as_ptr().add(l_idx) as *const i8,
                                    _MM_HINT_T0,
                                );
                            }
                        }
                        let mut new_smem = self.backwardExt(smem, a);
                        new_smem.m = j_u32;

                        if (new_smem.s < min_intv) && (smem.n - smem.m + 1 >= min_seed_len_u32) {
                            match_array.push(smem);
                            *total += 1;
                            break;
                        }
                        if (new_smem.s >= min_intv) && (new_smem.s != curr_s) {
                            curr_s = new_smem.s;
                            prev_array[num_curr] = new_smem;
                            num_curr += 1;
                            p += 1;
                            break;
                        }
                        p += 1;
                    }

                    while p < num_prev {
                        let smem = prev_array[p];
                        #[cfg(target_arch = "x86_64")]
                        if p + 1 < num_prev {
                            use core::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
                            let next_smem = prev_array[p + 1];
                            let k_idx = (next_smem.k >> CP_SHIFT) as usize;
                            let l_idx = ((next_smem.k + next_smem.s) >> CP_SHIFT) as usize;
                            unsafe {
                                _mm_prefetch(
                                    self.cp_occ.as_ptr().add(k_idx) as *const i8,
                                    _MM_HINT_T0,
                                );
                                _mm_prefetch(
                                    self.cp_occ.as_ptr().add(l_idx) as *const i8,
                                    _MM_HINT_T0,
                                );
                            }
                        }
                        let mut new_smem = self.backwardExt(smem, a);
                        new_smem.m = j_u32;
                        if (new_smem.s >= min_intv) && (new_smem.s != curr_s) {
                            curr_s = new_smem.s;
                            prev_array[num_curr] = new_smem;
                            num_curr += 1;
                        }
                        p += 1;
                    }

                    num_prev = num_curr;
                    if num_curr == 0 {
                        break;
                    }
                }

                if num_prev != 0 {
                    let smem = prev_array[0];
                    if smem.n - smem.m + 1 >= min_seed_len_u32 {
                        match_array.push(smem);
                        *total += 1;
                    }
                }
            }

            query_pos_array[i] = next_x as i16;
        }
    }

    /// Drive [`getSMEMsOnePosOneThread`] until every read on this thread is
    /// fully consumed.
    ///
    /// Maintains a compacting active-list: after each one-position pass the
    /// arrays are repacked to keep only reads whose search head has not yet
    /// reached `l_seq`. Resets `*num_total_smem` to zero before counting.
    ///
    /// # Arguments
    /// * `enc_qdb` - flat 0..3 encoded read buffer.
    /// * `min_intv_array` - per-read minimum BWT interval.
    /// * `rid_array` - per-read identifiers (compacted in-place).
    /// * `numReads` - initial number of reads.
    /// * `batch_size` - forwarded to the per-step driver.
    /// * `seq_` - read records (for sequence length lookup).
    /// * `query_cum_len_ar` - cumulative read offsets into `enc_qdb`.
    /// * `max_readlength` - upper bound on read length.
    /// * `minSeedLen` - minimum seed length.
    /// * `match_array` - SMEM output.
    /// * `num_total_smem` - total emitted SMEM count.
    #[doc = "Original function: FMI_search::getSMEMsAllPosOneThread:672"]
    pub fn getSMEMsAllPosOneThread(
        &self,
        enc_qdb: &[u8],
        min_intv_array: &mut [i32],
        rid_array: &mut [i32],
        numReads: i32,
        batch_size: i32,
        seq_: &[bseq1_t],
        query_cum_len_ar: &[i32],
        max_readlength: i32,
        minSeedLen: i32,
        match_array: &mut Vec<SMEM>,
        num_total_smem: &mut i64,
    ) {
        let num_reads = numReads as usize;
        let mut num_active = numReads;
        *num_total_smem = 0;

        SMEM_QUERY_POS_SCRATCH.with(|cell| {
            let mut query_pos_array = cell.borrow_mut();
            query_pos_array.clear();
            query_pos_array.resize(num_reads, 0);

            loop {
                let mut tail = 0_usize;
                for head in 0..(num_active as usize) {
                    let readlength = seq_[rid_array[head] as usize].l_seq;
                    if i32::from(query_pos_array[head]) < readlength {
                        rid_array[tail] = rid_array[head];
                        query_pos_array[tail] = query_pos_array[head];
                        min_intv_array[tail] = min_intv_array[head];
                        tail += 1;
                    }
                }

                self.getSMEMsOnePosOneThread(
                    enc_qdb,
                    &mut query_pos_array[..tail],
                    &mut min_intv_array[..tail],
                    &mut rid_array[..tail],
                    tail as i32,
                    batch_size,
                    seq_,
                    query_cum_len_ar,
                    max_readlength,
                    minSeedLen,
                    match_array,
                    num_total_smem,
                );
                num_active = tail as i32;
                if num_active <= 0 {
                    break;
                }
            }
        });
    }

    /// Forward-only seed enumeration used by the alternate seeding strategy.
    ///
    /// For each read, scans left-to-right and accepts the first SMEM whose
    /// interval drops below `max_intv_array[i]` while its length is still at
    /// least `minSeedLen`; the seed is appended to `match_array` and the scan
    /// advances past it.
    ///
    /// # Arguments
    /// * `enc_qdb` - flat 0..3 encoded read buffer.
    /// * `max_intv_array` - per-read maximum BWT interval at which to accept a seed.
    /// * `numReads` - number of reads.
    /// * `seq_` - read records.
    /// * `query_cum_len_ar` - cumulative offsets of each read into `enc_qdb`.
    /// * `minSeedLen` - minimum seed length to emit.
    /// * `match_array` - seed output.
    ///
    /// # Returns
    /// Number of seeds appended.
    #[doc = "Original function: FMI_search::bwtSeedStrategyAllPosOneThread:726"]
    pub fn bwtSeedStrategyAllPosOneThread(
        &self,
        enc_qdb: &[u8],
        max_intv_array: &[i32],
        numReads: i32,
        seq_: &[bseq1_t],
        query_cum_len_ar: &[i32],
        minSeedLen: i32,
        match_array: &mut Vec<SMEM>,
    ) -> i64 {
        let mut num_total_seed = 0_i64;
        let limit = usize::min(
            numReads as usize,
            usize::min(
                seq_.len(),
                usize::min(max_intv_array.len(), query_cum_len_ar.len()),
            ),
        );

        for i in 0..limit {
            let readlength = seq_[i].l_seq;
            let mut x = 0_i16;
            while x < readlength as i16 {
                let mut next_x = x + 1;
                // Forward search
                let mut smem = SMEM {
                    rid: i as u32,
                    m: x as u32,
                    n: x as u32,
                    ..Default::default()
                };

                let offset = query_cum_len_ar[i] as usize;
                let a = enc_qdb[offset + x as usize];
                if a < 4 {
                    smem.k = self.count[usize::from(a)];
                    smem.l = self.count[usize::from(3 - a)];
                    smem.s = self.count[usize::from(a + 1)] - self.count[usize::from(a)];

                    for j in (i32::from(x) + 1)..readlength {
                        next_x = (j + 1) as i16;
                        let a = enc_qdb[offset + j as usize];
                        if a < 4 {
                            // Forward extension is backward extension with the BWT of reverse complement
                            let mut smem_ = smem;
                            smem_.k = smem.l;
                            smem_.l = smem.k;
                            let new_smem_ = self.backwardExt(smem_, 3 - a);
                            let mut new_smem = new_smem_;
                            new_smem.k = new_smem_.l;
                            new_smem.l = new_smem_.k;
                            new_smem.n = j as u32;
                            smem = new_smem;

                            if (smem.s < i64::from(max_intv_array[i]))
                                && ((smem.n as i32 - smem.m as i32 + 1) >= minSeedLen)
                            {
                                if smem.s > 0 {
                                    match_array.push(smem);
                                    num_total_seed += 1;
                                }
                                break;
                            }
                        } else {
                            break;
                        }
                    }
                }
                x = next_x;
            }
        }
        num_total_seed
    }

    /// Sequential per-thread SMEM enumerator for fixed-width reads.
    ///
    /// Partitions the `numReads` reads across `nthreads` quota-balanced slices
    /// (executed serially in this port) and runs the bidirectional SMEM walk
    /// over each read. Each thread writes its SMEMs into `match_array` at a
    /// reserved `tid * readlength` offset and reports the count via
    /// `num_total_smem[tid]`.
    ///
    /// # Arguments
    /// * `enc_qdb` - flat 0..3 encoded read buffer (all reads share `readlength`).
    /// * `numReads` - number of reads to process.
    /// * `readlength` - per-read length (uniform across the batch).
    /// * `minSeedLen` - minimum seed length to keep.
    /// * `nthreads` - thread-slot count used to partition the work.
    /// * `match_array` - destination for emitted SMEMs.
    /// * `num_total_smem` - per-thread SMEM counts (resized to `nthreads`).
    #[doc = "Original function: FMI_search::getSMEMs:815"]
    pub fn getSMEMs(
        &self,
        enc_qdb: &[u8],
        numReads: i32,
        _batch_size: i32,
        readlength: i32,
        minSeedLen: i32,
        nthreads: i32,
        match_array: &mut Vec<SMEM>,
        num_total_smem: &mut Vec<i64>,
    ) {
        let nthreads_usize = usize::try_from(nthreads).expect("nthreads");
        let readlength_usize = usize::try_from(readlength).expect("readlength");
        let mut prev_array = vec![SMEM::default(); nthreads_usize * readlength_usize];
        let mut curr_array = vec![SMEM::default(); nthreads_usize * readlength_usize];

        if num_total_smem.len() < nthreads_usize {
            num_total_smem.resize(nthreads_usize, 0);
        }

        for tid in 0..nthreads_usize {
            num_total_smem[tid] = 0;
            let prev_base = tid * readlength_usize;
            let curr_base = tid * readlength_usize;

            let per_thread_quota = (numReads + (nthreads - 1)) / nthreads;
            let first = i32::try_from(tid).expect("tid") * per_thread_quota;
            let mut last = (i32::try_from(tid).expect("tid") + 1) * per_thread_quota;
            if last > numReads {
                last = numReads;
            }
            let match_base = usize::try_from(first * readlength).expect("match base");

            // Perform SMEM for original reads
            for i in first..last {
                let mut x = readlength - 1;
                let mut num_prev = 0_usize;
                let mut num_smem = 0_usize;

                while x >= 0 {
                    // Forward search
                    let mut smem = SMEM {
                        rid: u32::try_from(i).expect("rid"),
                        m: u32::try_from(x).expect("m"),
                        n: u32::try_from(x).expect("n"),
                        ..Default::default()
                    };
                    let a = enc_qdb[usize::try_from(i * readlength + x).expect("enc idx")];

                    if a > 3 {
                        x -= 1;
                        continue;
                    }
                    smem.k = self.count[usize::from(a)];
                    smem.l = self.count[usize::from(3 - a)];
                    smem.s = self.count[usize::from(a + 1)] - self.count[usize::from(a)];

                    for j in (x + 1)..readlength {
                        let a = enc_qdb[usize::try_from(i * readlength + j).expect("enc idx")];
                        if a < 4 {
                            // Forward extension is backward extension with the BWT of reverse complement
                            let mut smem_ = smem;
                            smem_.k = smem.l;
                            smem_.l = smem.k;
                            let new_smem_ = self.backwardExt(smem_, 3 - a);
                            let mut new_smem = new_smem_;
                            new_smem.k = new_smem_.l;
                            new_smem.l = new_smem_.k;
                            new_smem.n = u32::try_from(j).expect("n");

                            if new_smem.s != smem.s {
                                prev_array[prev_base + num_prev] = smem;
                                num_prev += 1;
                            }
                            smem = new_smem;
                            if new_smem.s == 0 {
                                break;
                            }
                        } else {
                            prev_array[prev_base + num_prev] = smem;
                            num_prev += 1;
                            break;
                        }
                    }
                    if smem.s != 0 {
                        prev_array[prev_base + num_prev] = smem;
                        num_prev += 1;
                    }

                    prev_array[prev_base..prev_base + num_prev].reverse();
                    let mut use_prev = true;
                    let mut next_x = x - 1;
                    let mut cur_j = readlength;

                    // Backward search
                    for j in (0..x).rev() {
                        let mut num_curr = 0_usize;
                        let mut curr_s = -1_i64;
                        let a = enc_qdb[usize::try_from(i * readlength + j).expect("enc idx")];
                        if a > 3 {
                            next_x = j - 1;
                            break;
                        }

                        for p in 0..num_prev {
                            let smem = if use_prev {
                                prev_array[prev_base + p]
                            } else {
                                curr_array[curr_base + p]
                            };
                            let mut new_smem = self.backwardExt(smem, a);
                            new_smem.m = u32::try_from(j).expect("m");

                            if new_smem.s == 0 {
                                if (num_curr == 0) && (j < cur_j) {
                                    cur_j = j;
                                    if (i32::try_from(smem.n).expect("n")
                                        - i32::try_from(smem.m).expect("m")
                                        + 1)
                                        >= minSeedLen
                                    {
                                        let target = match_base
                                            + usize::try_from(num_total_smem[tid])
                                                .expect("num_total")
                                            + num_smem;
                                        if match_array.len() <= target {
                                            match_array.resize(target + 1, SMEM::default());
                                        }
                                        match_array[target] = smem;
                                        num_smem += 1;
                                    }
                                }
                            }
                            if (new_smem.s != 0) && (new_smem.s != curr_s) {
                                curr_s = new_smem.s;
                                if use_prev {
                                    curr_array[curr_base + num_curr] = new_smem;
                                } else {
                                    prev_array[prev_base + num_curr] = new_smem;
                                }
                                num_curr += 1;
                            }
                        }

                        use_prev = !use_prev;
                        num_prev = num_curr;
                        if num_curr == 0 {
                            next_x = j;
                            break;
                        } else {
                            next_x = j - 1;
                        }
                    }

                    if num_prev != 0 {
                        let smem = if use_prev {
                            prev_array[prev_base]
                        } else {
                            curr_array[curr_base]
                        };
                        if (i32::try_from(smem.n).expect("n") - i32::try_from(smem.m).expect("m")
                            + 1)
                            >= minSeedLen
                        {
                            let target = match_base
                                + usize::try_from(num_total_smem[tid]).expect("num_total")
                                + num_smem;
                            if match_array.len() <= target {
                                match_array.resize(target + 1, SMEM::default());
                            }
                            match_array[target] = smem;
                            num_smem += 1;
                        }
                        num_prev = 0;
                    }
                    x = next_x;
                }

                num_total_smem[tid] += i64::try_from(num_smem).expect("num_smem");
            }
        }
    }

    /// Sort each thread's SMEM slice in-place using [`compare_smem_ord`].
    ///
    /// # Arguments
    /// * `match_array` - flat per-thread layout produced by [`getSMEMs`].
    /// * `num_total_smem` - per-thread emitted counts (slice length).
    /// * `num_reads` / `readlength` - same partitioning as in [`getSMEMs`],
    ///   used to recover each thread's base offset.
    /// * `nthreads` - number of thread slices.
    #[doc = "Original function: FMI_search::sortSMEMs:1008"]
    pub fn sortSMEMs(
        &self,
        match_array: &mut [SMEM],
        num_total_smem: &[i64],
        num_reads: i32,
        readlength: i32,
        nthreads: i32,
    ) {
        let per_thread_quota = (num_reads + (nthreads - 1)) / nthreads;
        for tid in 0..nthreads as usize {
            let first = tid as i32 * per_thread_quota;
            let base = (first * readlength) as usize;
            let len = num_total_smem[tid] as usize;
            let end = base + len;
            match_array[base..end].sort_by(compare_smem_ord);
        }
    }

    /// One step of bidirectional FM-index extension to the left by base `a`.
    ///
    /// Computes the new `(k, l, s)` triplet for all four candidate bases from
    /// the checkpointed occurrence table, then writes the chosen base's
    /// result back into `smem`. The forward strand of the index handles
    /// rightward extension by passing the `3 - a` complement and swapping
    /// `k`/`l` around the call (mirroring the C++ implementation).
    ///
    /// # Arguments
    /// * `smem` - current SMEM interval; consumed and returned with updated `k`/`l`/`s`.
    /// * `a` - base to extend with (0..3); other values fall through to the `3` branch.
    #[doc = "Original function: FMI_search::backwardExt:1025"]
    #[inline(always)]
    pub fn backwardExt(&self, mut smem: SMEM, a: u8) -> SMEM {
        // Hoist the per-position lookup that's invariant across the 4 base loop. `entry_sp` /
        // `entry_ep` and `mask_sp` / `mask_ep` only depend on `sp` / `ep`, not on `b`.
        let sp = smem.k;
        let ep = smem.k + smem.s;
        let occ_id_sp = (sp >> CP_SHIFT) as usize;
        let y_sp = (sp & CP_MASK as i64) as usize;
        let occ_id_ep = (ep >> CP_SHIFT) as usize;
        let y_ep = (ep & CP_MASK as i64) as usize;
        let (k0, s0, k1, s1, k2, s2, k3, s3) = unsafe {
            let entry_sp = self.cp_occ.get_unchecked(occ_id_sp);
            let entry_ep = self.cp_occ.get_unchecked(occ_id_ep);
            let mask_sp = *self.one_hot_mask_array.get_unchecked(y_sp);
            let mask_ep = *self.one_hot_mask_array.get_unchecked(y_ep);

            let one_hot_sp0 = *entry_sp.one_hot_bwt_str.get_unchecked(0);
            let one_hot_ep0 = *entry_ep.one_hot_bwt_str.get_unchecked(0);
            let occ_sp0 = *entry_sp.cp_count.get_unchecked(0)
                + i64::from(countbits_64(one_hot_sp0 & mask_sp));
            let occ_ep0 = *entry_ep.cp_count.get_unchecked(0)
                + i64::from(countbits_64(one_hot_ep0 & mask_ep));

            let one_hot_sp1 = *entry_sp.one_hot_bwt_str.get_unchecked(1);
            let one_hot_ep1 = *entry_ep.one_hot_bwt_str.get_unchecked(1);
            let occ_sp1 = *entry_sp.cp_count.get_unchecked(1)
                + i64::from(countbits_64(one_hot_sp1 & mask_sp));
            let occ_ep1 = *entry_ep.cp_count.get_unchecked(1)
                + i64::from(countbits_64(one_hot_ep1 & mask_ep));

            let one_hot_sp2 = *entry_sp.one_hot_bwt_str.get_unchecked(2);
            let one_hot_ep2 = *entry_ep.one_hot_bwt_str.get_unchecked(2);
            let occ_sp2 = *entry_sp.cp_count.get_unchecked(2)
                + i64::from(countbits_64(one_hot_sp2 & mask_sp));
            let occ_ep2 = *entry_ep.cp_count.get_unchecked(2)
                + i64::from(countbits_64(one_hot_ep2 & mask_ep));

            let one_hot_sp3 = *entry_sp.one_hot_bwt_str.get_unchecked(3);
            let one_hot_ep3 = *entry_ep.one_hot_bwt_str.get_unchecked(3);
            let occ_sp3 = *entry_sp.cp_count.get_unchecked(3)
                + i64::from(countbits_64(one_hot_sp3 & mask_sp));
            let occ_ep3 = *entry_ep.cp_count.get_unchecked(3)
                + i64::from(countbits_64(one_hot_ep3 & mask_ep));

            (
                *self.count.get_unchecked(0) + occ_sp0,
                occ_ep0 - occ_sp0,
                *self.count.get_unchecked(1) + occ_sp1,
                occ_ep1 - occ_sp1,
                *self.count.get_unchecked(2) + occ_sp2,
                occ_ep2 - occ_sp2,
                *self.count.get_unchecked(3) + occ_sp3,
                occ_ep3 - occ_sp3,
            )
        };

        let mut sentinel_offset = 0_i64;
        if smem.k <= self.sentinel_index && (smem.k + smem.s) > self.sentinel_index {
            sentinel_offset = 1;
        }
        let l3 = smem.l + sentinel_offset;
        let l2 = l3 + s3;
        let l1 = l2 + s2;
        let l0 = l1 + s1;

        match a {
            0 => {
                smem.k = k0;
                smem.l = l0;
                smem.s = s0;
            }
            1 => {
                smem.k = k1;
                smem.l = l1;
                smem.s = s1;
            }
            2 => {
                smem.k = k2;
                smem.l = l2;
                smem.s = s2;
            }
            _ => {
                smem.k = k3;
                smem.l = l3;
                smem.s = s3;
            }
        }
        smem
    }

    /// Read the suffix-array value at BWT position `pos` from the uncompressed
    /// `(ms_byte, ls_word)` arrays (i.e. SA-compression factor 1).
    ///
    /// Reassembles the 40-bit SA entry as `(ms_byte << 32) | ls_word`.
    #[doc = "Original function: FMI_search::get_sa_entry:1054"]
    #[inline]
    pub fn get_sa_entry(&self, pos: i64) -> i64 {
        // BWT positions are non-negative and bounded by sa_ms_byte.len() == sa_ls_word.len()
        // by index construction. Caller-side correctness is enforced; skip the runtime checks.
        debug_assert!(pos >= 0 && (pos as usize) < self.sa_ms_byte.len());
        debug_assert!((pos as usize) < self.sa_ls_word.len());
        let pos_us = pos as usize;
        let ms = unsafe { *self.sa_ms_byte.get_unchecked(pos_us) };
        let ls = unsafe { *self.sa_ls_word.get_unchecked(pos_us) };
        (i64::from(ms) << 32) + i64::from(ls)
    }

    /// Bulk-resolve SA entries for an array of BWT positions.
    ///
    /// Writes `coord_array[i] = get_sa_entry(pos_array[i])` for the first
    /// `count` positions; grows `coord_array` to fit when necessary. Matches
    /// the uncompressed SA overload of `FMI_search::get_sa_entries`.
    ///
    /// # Arguments
    /// * `pos_array` - input BWT positions.
    /// * `coord_array` - output SA entries.
    /// * `count` - number of positions to process.
    #[doc = "Original function: FMI_search::get_sa_entries:1062"]
    pub fn get_sa_entries__L1062(
        &self,
        pos_array: &[i64],
        coord_array: &mut Vec<i64>,
        count: u32,
        _nthreads: i32,
    ) {
        let limit = usize::min(pos_array.len(), count as usize);
        if coord_array.len() < limit {
            coord_array.resize(limit, 0);
        }
        for (i, &pos) in pos_array[..limit].iter().enumerate() {
            coord_array[i] = self.get_sa_entry(pos);
        }
    }

    /// Expand each SMEM's BWT interval into up to `max_occ` SA-derived genomic
    /// coordinates using uncompressed [`get_sa_entry`] lookups.
    ///
    /// When the interval is larger than `max_occ`, samples it at a stride of
    /// `smem.s / max_occ`. Records the emitted count per SMEM in
    /// `coord_count_array`.
    ///
    /// # Arguments
    /// * `smem_array` - SMEMs whose intervals to expand.
    /// * `coord_array` - destination coordinates; grown as needed.
    /// * `coord_count_array` - per-SMEM emitted counts.
    /// * `count` - number of SMEMs to process.
    /// * `max_occ` - cap on coordinates emitted per SMEM.
    #[doc = "Original function: FMI_search::get_sa_entries:1077"]
    pub fn get_sa_entries__L1077(
        &self,
        smem_array: &[SMEM],
        coord_array: &mut Vec<i64>,
        coord_count_array: &mut Vec<i32>,
        count: u32,
        max_occ: i32,
    ) {
        let mut total_coord_count = 0_usize;
        let limit = usize::min(smem_array.len(), count as usize);
        if coord_count_array.len() < limit {
            coord_count_array.resize(limit, 0);
        }
        for (i, smem) in smem_array[..limit].iter().enumerate() {
            let mut c = 0_i32;
            let hi = smem.k + smem.s;
            let step = if smem.s > i64::from(max_occ) {
                smem.s / i64::from(max_occ)
            } else {
                1
            };
            let mut j = smem.k;
            while j < hi && c < max_occ {
                let sa_entry = self.get_sa_entry(j);
                let target = total_coord_count + c as usize;
                if coord_array.len() <= target {
                    coord_array.push(sa_entry);
                } else {
                    coord_array[target] = sa_entry;
                }
                j += step;
                c += 1;
            }
            coord_count_array[i] = c;
            total_coord_count += c as usize;
        }
    }

    // sa_compression
    /// Resolve a single SA entry against the SA-compressed index
    /// (compression factor `SA_COMPX`).
    ///
    /// Positions aligned to the `SA_COMPX_MASK` lattice are read directly
    /// from `sa_ms_byte`/`sa_ls_word`; otherwise the BWT is walked one
    /// position at a time until a compressed slot is hit, accumulating an
    /// offset that is added to the stored entry. The sentinel base aborts
    /// the walk with offset only.
    ///
    /// # Arguments
    /// * `pos` - BWT position to resolve.
    /// * `_tid` - thread id (unused in the Rust port; preserved for parity).
    #[doc = "Original function: FMI_search::get_sa_entry_compressed:1103"]
    #[inline]
    pub fn get_sa_entry_compressed(&self, pos: i64, _tid: i32) -> i64 {
        if (pos & SA_COMPX_MASK as i64) == 0 {
            let idx = (pos >> SA_COMPX) as usize;
            debug_assert!(idx < self.sa_ms_byte.len() && idx < self.sa_ls_word.len());
            unsafe {
                let ms = *self.sa_ms_byte.get_unchecked(idx);
                let ls = *self.sa_ls_word.get_unchecked(idx);
                return (i64::from(ms) << 32) + i64::from(ls);
            }
        }

        let mut offset = 0_i64;
        let mut sp = pos;
        loop {
            let b = decode_bwt_base(&self.cp_occ, sp);
            if b == 4 {
                return offset;
            }
            let occ_sp = occ(&self.cp_occ, &self.one_hot_mask_array, sp, b);
            sp = self.count[usize::from(b)] + occ_sp;
            offset += 1;
            if (sp & SA_COMPX_MASK as i64) == 0 {
                break;
            }
        }

        let idx = (sp >> SA_COMPX) as usize;
        debug_assert!(idx < self.sa_ms_byte.len() && idx < self.sa_ls_word.len());
        unsafe {
            let ms = *self.sa_ms_byte.get_unchecked(idx);
            let ls = *self.sa_ls_word.get_unchecked(idx);
            (i64::from(ms) << 32) + i64::from(ls) + offset
        }
    }

    /// Expand each SMEM's BWT interval into up to `max_occ` genomic
    /// coordinates using the SA-compressed lookup
    /// ([`get_sa_entry_compressed`]).
    ///
    /// Same striding behaviour as the uncompressed variant; aggregates the
    /// total coordinate count into `*coord_count_array`.
    ///
    /// # Arguments
    /// * `smem_array` - input SMEMs.
    /// * `coord_array` - destination coordinates.
    /// * `coord_count_array` - running coordinate total; incremented by the count.
    /// * `count` - number of SMEMs to process.
    /// * `max_occ` - cap on coordinates emitted per SMEM.
    /// * `tid` - thread id passed through to the compressed walker.
    #[doc = "Original function: FMI_search::get_sa_entries:1177"]
    pub fn get_sa_entries__L1177(
        &self,
        smem_array: &[SMEM],
        coord_array: &mut Vec<i64>,
        coord_count_array: &mut i32,
        count: u32,
        max_occ: i32,
        tid: i32,
    ) {
        let mut total_coord_count = 0_usize;
        let limit = usize::min(smem_array.len(), count as usize);
        if coord_array.is_empty() {
            coord_array.reserve(limit);
        }
        for smem in &smem_array[..limit] {
            let mut c = 0_i32;
            let hi = smem.k + smem.s;
            let step = if smem.s > i64::from(max_occ) {
                smem.s / i64::from(max_occ)
            } else {
                1
            };
            let mut j = smem.k;
            while j < hi && c < max_occ {
                let sa_entry = self.get_sa_entry_compressed(j, tid);
                let target = total_coord_count + c as usize;
                if coord_array.len() <= target {
                    coord_array.push(sa_entry);
                } else {
                    coord_array[target] = sa_entry;
                }
                j += step;
                c += 1;
            }
            *coord_count_array += c;
            total_coord_count += c as usize;
        }
    }

    // SA_COMPRESSION w/ PREFETCH
    /// Take one step of the SA-compressed BWT walk and report whether the
    /// final SA entry is now known.
    ///
    /// If `pos` is aligned to a compressed SA slot, fills `*sa_entry` and
    /// returns `1`. Otherwise decodes the BWT base, advances `pos`, and
    /// returns `1` only when the next position is itself a compressed slot
    /// (whose value is recorded with the accumulated `*offset` applied); it
    /// returns `0` to signal that the caller must continue walking.
    ///
    /// # Arguments
    /// * `pos` - current BWT position.
    /// * `sa_entry` - resolved SA entry or next walk position, written in-place.
    /// * `offset` - running offset accumulated during the walk.
    ///
    /// # Returns
    /// `1` when `*sa_entry` holds a final SA value, `0` when the walk must continue.
    #[doc = "Original function: FMI_search::call_one_step:1202"]
    #[inline]
    pub fn call_one_step(&self, pos: i64, sa_entry: &mut i64, offset: &mut i64) -> i64 {
        if (pos & SA_COMPX_MASK as i64) == 0 {
            let idx = (pos >> SA_COMPX) as usize;
            debug_assert!(idx < self.sa_ms_byte.len() && idx < self.sa_ls_word.len());
            unsafe {
                let ms = *self.sa_ms_byte.get_unchecked(idx);
                let ls = *self.sa_ls_word.get_unchecked(idx);
                *sa_entry = (i64::from(ms) << 32) + i64::from(ls);
            }
            return 1;
        }

        let sp = pos;
        let b = decode_bwt_base(&self.cp_occ, sp);
        if b == 4 {
            *sa_entry = 0;
            return 1;
        }

        let occ_sp = occ(&self.cp_occ, &self.one_hot_mask_array, sp, b);
        let next_sp = self.count[usize::from(b)] + occ_sp;
        *offset += 1;
        if (next_sp & SA_COMPX_MASK as i64) == 0 {
            let idx = (next_sp >> SA_COMPX) as usize;
            debug_assert!(idx < self.sa_ms_byte.len() && idx < self.sa_ls_word.len());
            unsafe {
                let ms = *self.sa_ms_byte.get_unchecked(idx);
                let ls = *self.sa_ls_word.get_unchecked(idx);
                *sa_entry = (i64::from(ms) << 32) + i64::from(ls) + *offset;
            }
            1
        } else {
            *sa_entry = next_sp;
            0
        }
    }

    /// Prefetch-aware variant of [`get_sa_entries__L1177`] that drives
    /// [`call_one_step`] manually so the next walk's `cp_occ` line can be
    /// warmed under the current step's latency.
    ///
    /// Buffers the per-SMEM positions and their destination indices in
    /// thread-local scratch, then walks them in order, prefetching the next
    /// position's `cp_occ` entry on `x86_64`. Updates `*coord_count_array`
    /// with the total emitted count and accumulates the number of walks
    /// driven into `*id_`.
    ///
    /// # Arguments
    /// * `smem_array` - SMEMs to resolve.
    /// * `coord_array` - destination genomic coordinates.
    /// * `coord_count_array` - running coordinate total.
    /// * `count` - number of SMEMs to process.
    /// * `max_occ` - cap on coordinates emitted per SMEM.
    /// * `_tid` - thread id, unused in the Rust port.
    /// * `id_` - counter incremented by the number of walks started.
    #[doc = "Original function: FMI_search::get_sa_entries_prefetch:1257"]
    pub fn get_sa_entries_prefetch(
        &self,
        smem_array: &[SMEM],
        coord_array: &mut Vec<i64>,
        coord_count_array: &mut i64,
        count: i64,
        max_occ: i32,
        _tid: i32,
        id_: &mut i64,
    ) {
        SA_PREFETCH_SCRATCH.with(|cell| {
            let mut buf = cell.borrow_mut();
            let (pos_ar, map_ar) = &mut *buf;
            pos_ar.clear();
            map_ar.clear();
            let mut total_coord_count = 0_usize;
            let limit = usize::min(smem_array.len(), count as usize);

            for smem in &smem_array[..limit] {
                let mut c = 0_i32;
                let hi = smem.k + smem.s;
                let step = if smem.s > i64::from(max_occ) {
                    smem.s / i64::from(max_occ)
                } else {
                    1
                };
                let mut j = smem.k;
                while j < hi && c < max_occ {
                    pos_ar.push(j);
                    map_ar.push(total_coord_count + c as usize);
                    j += step;
                    c += 1;
                }
                *coord_count_array += i64::from(c);
                total_coord_count += c as usize;
            }

            *id_ += pos_ar.len() as i64;
            if coord_array.len() < total_coord_count {
                coord_array.resize(total_coord_count, 0);
            }

            for (idx, &pos) in pos_ar.iter().enumerate() {
                // Prefetch cp_occ for the NEXT pos_ar entry while we walk the current one.
                // Each call_one_step does cp_occ[pos>>CP_SHIFT] if the position isn't at a
                // compressed SA slot — i.e., 7 of 8 cases. The walks are independent, so
                // pre-loading the next pos's cp_occ entry hides DRAM latency.
                #[cfg(target_arch = "x86_64")]
                if let Some(&next_pos) = pos_ar.get(idx + 1) {
                    use core::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
                    let next_idx = (next_pos >> CP_SHIFT) as usize;
                    unsafe {
                        _mm_prefetch(self.cp_occ.as_ptr().add(next_idx) as *const i8, _MM_HINT_T0);
                    }
                }
                let mut working = pos;
                let mut sp = 0_i64;
                let mut offset = 0_i64;
                loop {
                    let quit = self.call_one_step(working, &mut sp, &mut offset);
                    if quit != 0 {
                        coord_array[map_ar[idx]] = sp;
                        break;
                    }
                    working = sp;
                }
            }
        });
    }
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::io::Cursor;
    use std::io::Write;
    use std::time::{SystemTime, UNIX_EPOCH};

    use crate::bwa_mem2::bntseq::bns_fasta2bntseq;
    use crate::bwa_mem2::bwa::bseq1_t;

    use super::FMI_search;
    use crate::bwa_mem2::fmi_search::SMEM;

    fn reference_and_sa(prefix: &std::path::Path) -> (String, Vec<i64>) {
        let fmi = FMI_search::ctor(prefix.to_str().expect("utf8"));
        let mut reference_seq = String::new();
        fmi.pac2nt(
            prefix.with_extension("pac").to_str().expect("utf8"),
            &mut reference_seq,
        );
        let mut suffixes: Vec<i64> =
            (0..i64::try_from(reference_seq.len()).expect("len")).collect();
        suffixes.sort_by(|&a, &b| {
            let sa = &reference_seq.as_bytes()[usize::try_from(a).expect("a")..];
            let sb = &reference_seq.as_bytes()[usize::try_from(b).expect("b")..];
            sa.cmp(sb)
        });
        let mut sa_bwt = Vec::with_capacity(suffixes.len() + 1);
        sa_bwt.push(i64::try_from(reference_seq.len()).expect("len"));
        sa_bwt.extend(suffixes);
        (reference_seq, sa_bwt)
    }

    fn expected_prefetch_sa_entry(fmi: &FMI_search, pos: i64) -> i64 {
        let mut offset = 0_i64;
        let mut working = pos;
        loop {
            if (working & 0x7) == 0 {
                let idx = usize::try_from(working >> 3).expect("compressed idx");
                let mut sa_entry = i64::from(fmi.sa_ms_byte[idx]);
                sa_entry <<= 32;
                sa_entry += i64::from(fmi.sa_ls_word[idx]);
                return sa_entry;
            }

            let b = super::decode_bwt_base(&fmi.cp_occ, working);
            if b == 4 {
                return 0;
            }
            let occ_sp = super::occ(&fmi.cp_occ, &fmi.one_hot_mask_array, working, b);
            working = fmi.count[usize::from(b)] + occ_sp;
            offset += 1;
            if (working & 0x7) == 0 {
                let idx = usize::try_from(working >> 3).expect("compressed idx");
                let mut sa_entry = i64::from(fmi.sa_ms_byte[idx]);
                sa_entry <<= 32;
                sa_entry += i64::from(fmi.sa_ls_word[idx]);
                sa_entry += offset;
                return sa_entry;
            }
        }
    }

    fn bwt_from_sa(reference_seq: &str, sa_bwt: &[i64]) -> Vec<u8> {
        let mut bwt = Vec::with_capacity(sa_bwt.len());
        for &sai in sa_bwt {
            if sai == 0 {
                bwt.push(4);
            } else {
                let c = reference_seq.as_bytes()[usize::try_from(sai - 1).expect("sai")];
                bwt.push(match c {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    _ => panic!("unexpected base"),
                });
            }
        }
        bwt
    }

    fn explicit_backward_ext(fmi: &FMI_search, bwt: &[u8], smem: SMEM, a: u8) -> SMEM {
        let mut k = [0_i64; 4];
        let mut l = [0_i64; 4];
        let mut s = [0_i64; 4];
        let sp = usize::try_from(smem.k).expect("sp");
        let ep = usize::try_from(smem.k + smem.s).expect("ep");
        for b in 0..4_u8 {
            let occ_sp = bwt[..sp].iter().filter(|&&x| x == b).count() as i64;
            let occ_ep = bwt[..ep].iter().filter(|&&x| x == b).count() as i64;
            k[usize::from(b)] = fmi.count[usize::from(b)] + occ_sp;
            s[usize::from(b)] = occ_ep - occ_sp;
        }

        let sentinel_offset =
            if smem.k <= fmi.sentinel_index && (smem.k + smem.s) > fmi.sentinel_index {
                1
            } else {
                0
            };
        l[3] = smem.l + sentinel_offset;
        l[2] = l[3] + s[3];
        l[1] = l[2] + s[2];
        l[0] = l[1] + s[1];

        SMEM {
            rid: smem.rid,
            m: smem.m,
            n: smem.n,
            k: k[usize::from(a)],
            l: l[usize::from(a)],
            s: s[usize::from(a)],
        }
    }

    fn explicit_seed_strategy(
        fmi: &FMI_search,
        bwt: &[u8],
        enc_qdb: &[u8],
        max_intv_array: &[i32],
        seqs: &[bseq1_t],
        query_cum_len_ar: &[i32],
        min_seed_len: i32,
    ) -> Vec<SMEM> {
        let mut out = Vec::new();
        for i in 0..seqs.len() {
            let readlength = seqs[i].l_seq;
            let mut x = 0_i32;
            while x < readlength {
                let mut next_x = x + 1;
                let mut smem = SMEM {
                    rid: u32::try_from(i).expect("rid"),
                    m: u32::try_from(x).expect("m"),
                    n: u32::try_from(x).expect("n"),
                    ..Default::default()
                };
                let offset = usize::try_from(query_cum_len_ar[i]).expect("offset");
                let a = enc_qdb[offset + usize::try_from(x).expect("x")];
                if a < 4 {
                    smem.k = fmi.count[usize::from(a)];
                    smem.l = fmi.count[usize::from(3 - a)];
                    smem.s = fmi.count[usize::from(a + 1)] - fmi.count[usize::from(a)];

                    for j in (x + 1)..readlength {
                        next_x = j + 1;
                        let a = enc_qdb[offset + usize::try_from(j).expect("j")];
                        if a < 4 {
                            let mut smem_ = smem;
                            smem_.k = smem.l;
                            smem_.l = smem.k;
                            let new_smem_ = explicit_backward_ext(fmi, bwt, smem_, 3 - a);
                            let mut new_smem = new_smem_;
                            new_smem.k = new_smem_.l;
                            new_smem.l = new_smem_.k;
                            new_smem.n = u32::try_from(j).expect("n");
                            smem = new_smem;
                            if (smem.s < i64::from(max_intv_array[i]))
                                && ((i32::try_from(smem.n).expect("n")
                                    - i32::try_from(smem.m).expect("m")
                                    + 1)
                                    >= min_seed_len)
                            {
                                if smem.s > 0 {
                                    out.push(smem);
                                }
                                break;
                            }
                        } else {
                            break;
                        }
                    }
                }
                x = next_x;
            }
        }
        out
    }

    #[test]
    fn compare_smem_orders_by_rid_then_m_then_reverse_n() {
        let a = SMEM {
            rid: 1,
            m: 5,
            n: 10,
            ..Default::default()
        };
        let b = SMEM {
            rid: 2,
            m: 1,
            n: 1,
            ..Default::default()
        };
        let c = SMEM {
            rid: 1,
            m: 6,
            n: 1,
            ..Default::default()
        };
        let d = SMEM {
            rid: 1,
            m: 5,
            n: 8,
            ..Default::default()
        };
        let e = SMEM {
            rid: 1,
            m: 5,
            n: 10,
            ..Default::default()
        };

        assert_eq!(super::compare_smem(&a, &b), -1);
        assert_eq!(super::compare_smem(&b, &a), 1);
        assert_eq!(super::compare_smem(&a, &c), -1);
        assert_eq!(super::compare_smem(&c, &a), 1);
        assert_eq!(super::compare_smem(&a, &d), -1);
        assert_eq!(super::compare_smem(&d, &a), 1);
        assert_eq!(super::compare_smem(&a, &e), 0);
    }

    #[test]
    fn sort_smems_sorts_each_thread_slice_only() {
        let fmi = FMI_search::default();
        let mut entries = vec![
            SMEM {
                rid: 1,
                m: 8,
                n: 2,
                ..Default::default()
            },
            SMEM {
                rid: 0,
                m: 9,
                n: 1,
                ..Default::default()
            },
            SMEM {
                rid: 0,
                m: 9,
                n: 3,
                ..Default::default()
            },
            SMEM {
                rid: 3,
                m: 4,
                n: 1,
                ..Default::default()
            },
            SMEM {
                rid: 2,
                m: 7,
                n: 5,
                ..Default::default()
            },
            SMEM {
                rid: 2,
                m: 7,
                n: 8,
                ..Default::default()
            },
            SMEM {
                rid: 99,
                m: 99,
                n: 99,
                ..Default::default()
            },
            SMEM {
                rid: 88,
                m: 88,
                n: 88,
                ..Default::default()
            },
        ];
        let counts = vec![3_i64, 3_i64];
        fmi.sortSMEMs(&mut entries, &counts, 4, 2, 2);

        assert_eq!(
            &entries[0..3]
                .iter()
                .map(|s| (s.rid, s.m, s.n))
                .collect::<Vec<_>>(),
            &vec![(0, 9, 3), (0, 9, 1), (1, 8, 2)]
        );
        assert_eq!(
            &entries[4..7]
                .iter()
                .map(|s| (s.rid, s.m, s.n))
                .collect::<Vec<_>>(),
            &vec![(2, 7, 8), (2, 7, 5), (99, 99, 99)]
        );
        assert_eq!((entries[3].rid, entries[7].rid), (3, 88));
    }

    #[test]
    fn build_index_creates_expected_marker_file() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_fmi_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("ref");
        let mut pac = std::fs::File::create(dir.join("ref.pac")).expect("create pac");
        pac.write_all(&[(0 << 6) | (1 << 4) | (2 << 2) | 3, 0, 0])
            .expect("write pac");
        let mut fmi = FMI_search::ctor(prefix.to_str().expect("utf8"));
        assert_eq!(
            fmi.pac_seq_len(dir.join("ref.pac").to_str().expect("utf8")),
            4
        );
        let mut seq = String::new();
        fmi.pac2nt(dir.join("ref.pac").to_str().expect("utf8"), &mut seq);
        assert_eq!(seq, "ACGTACGT");
        assert_eq!(fmi.build_index(), 0);
        assert!(dir.join("ref.bwt.2bit.64").is_file());
        assert!(dir.join("ref.0123").is_file());
        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn load_index_round_trips_generated_index_and_reference_metadata() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_fmi_load_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("ref");
        let fasta = b">chr1\nACGTAC\n>chr2\nTTAA\n";
        let l_pac = bns_fasta2bntseq(
            Cursor::new(fasta.as_slice()),
            prefix.to_str().expect("utf8"),
            1,
        );
        assert_eq!(l_pac, 10);

        let mut built = FMI_search::ctor(prefix.to_str().expect("utf8"));
        assert_eq!(built.build_index(), 0);

        let mut loaded = FMI_search::ctor(prefix.to_str().expect("utf8"));
        loaded.load_index();

        assert_eq!(loaded.reference_seq_len, 21);
        assert_eq!(loaded.sentinel_index, built.sentinel_index);
        assert_eq!(loaded.count, built.count.map(|x| x + 1));
        assert_eq!(loaded.cp_occ, built.cp_occ);
        assert_eq!(loaded.sa_ls_word, built.sa_ls_word);
        assert_eq!(loaded.sa_ms_byte, built.sa_ms_byte);
        assert_eq!(loaded.one_hot_mask_array.len(), 64);
        assert_eq!(loaded.base.idx.pac.len(), 3);
        let bns = loaded.base.idx.bns.as_ref().expect("loaded bns");
        assert_eq!(bns.n_seqs, 2);
        assert_eq!(bns.anns[0].name, "chr1");
        assert_eq!(bns.anns[1].name, "chr2");

        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn get_sa_entry_compressed_matches_expected_suffix_array_positions() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_fmi_sa_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("ref");
        let fasta = b">chr1\nACGTAC\n>chr2\nTTAA\n";
        let _ = bns_fasta2bntseq(
            Cursor::new(fasta.as_slice()),
            prefix.to_str().expect("utf8"),
            1,
        );

        let mut built = FMI_search::ctor(prefix.to_str().expect("utf8"));
        assert_eq!(built.build_index(), 0);

        let (_reference_seq, sa_bwt) = reference_and_sa(&prefix);
        let mut loaded = FMI_search::ctor(prefix.to_str().expect("utf8"));
        loaded.load_index();

        for (pos, expected) in sa_bwt.iter().enumerate() {
            assert_eq!(
                loaded.get_sa_entry_compressed(i64::try_from(pos).expect("pos"), 0),
                *expected,
                "compressed sa mismatch at pos {pos}"
            );
        }

        for pos in (0..sa_bwt.len()).step_by(8) {
            assert_eq!(
                loaded.get_sa_entry(i64::try_from(pos >> 3).expect("compressed idx")),
                sa_bwt[pos],
                "direct compressed slot mismatch at pos {pos}"
            );
        }

        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn build_index_sorts_ambiguous_reference_by_ascii_sequence() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_fmi_ambiguous_sa_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("ref");
        let fasta = b">chr1\nTNGCATNA\n";
        let _ = bns_fasta2bntseq(
            Cursor::new(fasta.as_slice()),
            prefix.to_str().expect("utf8"),
            1,
        );

        let mut built = FMI_search::ctor(prefix.to_str().expect("utf8"));
        assert_eq!(built.build_index(), 0);

        let (_reference_seq, sa_bwt) = reference_and_sa(&prefix);
        let mut loaded = FMI_search::ctor(prefix.to_str().expect("utf8"));
        loaded.load_index();

        for (pos, expected) in sa_bwt.iter().enumerate() {
            assert_eq!(
                loaded.get_sa_entry_compressed(i64::try_from(pos).expect("pos"), 0),
                *expected,
                "ambiguous reference SA mismatch at BWT pos {pos}"
            );
        }

        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn call_one_step_walks_until_a_compressed_sa_slot() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_fmi_step_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("ref");
        let fasta = b">chr1\nACGTAC\n>chr2\nTTAA\n";
        let _ = bns_fasta2bntseq(
            Cursor::new(fasta.as_slice()),
            prefix.to_str().expect("utf8"),
            1,
        );

        let mut fmi = FMI_search::ctor(prefix.to_str().expect("utf8"));
        assert_eq!(fmi.build_index(), 0);
        fmi.load_index();

        for pos in 0..usize::try_from(fmi.reference_seq_len).expect("ref len") {
            let mut working = i64::try_from(pos).expect("pos");
            let mut sa_entry = 0_i64;
            let mut offset = 0_i64;
            loop {
                let quit = fmi.call_one_step(working, &mut sa_entry, &mut offset);
                if quit != 0 {
                    assert_eq!(
                        sa_entry,
                        expected_prefetch_sa_entry(&fmi, i64::try_from(pos).expect("pos")),
                        "call_one_step mismatch at pos {pos}"
                    );
                    break;
                }
                working = sa_entry;
            }
        }

        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn prefetch_and_direct_compressed_entry_collection_match() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_fmi_collect_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("ref");
        let fasta = b">chr1\nACGTAC\n>chr2\nTTAA\n";
        let _ = bns_fasta2bntseq(
            Cursor::new(fasta.as_slice()),
            prefix.to_str().expect("utf8"),
            1,
        );

        let mut fmi = FMI_search::ctor(prefix.to_str().expect("utf8"));
        assert_eq!(fmi.build_index(), 0);
        fmi.load_index();

        let smems = vec![
            SMEM {
                k: 0,
                s: 5,
                ..Default::default()
            },
            SMEM {
                k: 6,
                s: 7,
                ..Default::default()
            },
        ];

        let mut direct = Vec::new();
        let mut direct_count = 0_i32;
        fmi.get_sa_entries__L1177(
            &smems,
            &mut direct,
            &mut direct_count,
            smems.len() as u32,
            3,
            0,
        );

        let mut prefetched = Vec::new();
        let mut prefetched_count = 0_i64;
        let mut id = 0_i64;
        fmi.get_sa_entries_prefetch(
            &smems,
            &mut prefetched,
            &mut prefetched_count,
            smems.len() as i64,
            3,
            0,
            &mut id,
        );

        let expected_prefetched: Vec<i64> = vec![
            expected_prefetch_sa_entry(&fmi, 0),
            expected_prefetch_sa_entry(&fmi, 1),
            expected_prefetch_sa_entry(&fmi, 2),
            expected_prefetch_sa_entry(&fmi, 6),
            expected_prefetch_sa_entry(&fmi, 8),
            expected_prefetch_sa_entry(&fmi, 10),
        ];

        assert_eq!(direct, vec![20, 12, 8, 13, 17, 5]);
        assert_eq!(prefetched, expected_prefetched);
        assert_eq!(prefetched_count, i64::from(direct_count));
        assert_eq!(id, prefetched_count);

        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn get_sa_entries_direct_positions_match_scalar_lookup() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_fmi_entries_pos_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("ref");
        let fasta = b">chr1\nACGTAC\n>chr2\nTTAA\n";
        let _ = bns_fasta2bntseq(
            Cursor::new(fasta.as_slice()),
            prefix.to_str().expect("utf8"),
            1,
        );

        let mut fmi = FMI_search::ctor(prefix.to_str().expect("utf8"));
        assert_eq!(fmi.build_index(), 0);
        fmi.load_index();

        let positions = vec![0_i64, 1, 2];
        let mut coords = Vec::new();
        fmi.get_sa_entries__L1062(&positions, &mut coords, positions.len() as u32, 1);
        let expected: Vec<i64> = positions.iter().map(|&p| fmi.get_sa_entry(p)).collect();
        assert_eq!(coords, expected);

        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn backward_ext_matches_bwt_interval_math_for_real_index() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_fmi_backward_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("ref");
        let fasta = b">chr1\nACGTAC\n>chr2\nTTAA\n";
        let _ = bns_fasta2bntseq(
            Cursor::new(fasta.as_slice()),
            prefix.to_str().expect("utf8"),
            1,
        );

        let mut fmi = FMI_search::ctor(prefix.to_str().expect("utf8"));
        assert_eq!(fmi.build_index(), 0);
        fmi.load_index();

        let (reference_seq, sa_bwt) = reference_and_sa(&prefix);
        let bwt = bwt_from_sa(&reference_seq, &sa_bwt);

        for seed_base in 0_u8..4 {
            let smem = SMEM {
                k: fmi.count[usize::from(seed_base)],
                l: fmi.count[usize::from(3 - seed_base)],
                s: fmi.count[usize::from(seed_base) + 1] - fmi.count[usize::from(seed_base)],
                ..Default::default()
            };

            for ext_base in 0_u8..4 {
                let actual = fmi.backwardExt(smem, ext_base);

                let sp = usize::try_from(smem.k).expect("sp");
                let ep = usize::try_from(smem.k + smem.s).expect("ep");
                let mut exp_k = [0_i64; 4];
                let mut exp_s = [0_i64; 4];
                for b in 0..4_u8 {
                    let occ_sp = bwt[..sp].iter().filter(|&&x| x == b).count() as i64;
                    let occ_ep = bwt[..ep].iter().filter(|&&x| x == b).count() as i64;
                    exp_k[usize::from(b)] = fmi.count[usize::from(b)] + occ_sp;
                    exp_s[usize::from(b)] = occ_ep - occ_sp;
                }

                let mut exp_l = [0_i64; 4];
                let sentinel_offset =
                    if smem.k <= fmi.sentinel_index && (smem.k + smem.s) > fmi.sentinel_index {
                        1
                    } else {
                        0
                    };
                exp_l[3] = smem.l + sentinel_offset;
                exp_l[2] = exp_l[3] + exp_s[3];
                exp_l[1] = exp_l[2] + exp_s[2];
                exp_l[0] = exp_l[1] + exp_s[1];

                assert_eq!(actual.k, exp_k[usize::from(ext_base)]);
                assert_eq!(actual.l, exp_l[usize::from(ext_base)]);
                assert_eq!(actual.s, exp_s[usize::from(ext_base)]);
            }
        }

        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn bwt_seed_strategy_matches_explicit_interval_walk() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_fmi_seed_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("ref");
        let fasta = b">chr1\nACGTAC\n>chr2\nTTAA\n";
        let _ = bns_fasta2bntseq(
            Cursor::new(fasta.as_slice()),
            prefix.to_str().expect("utf8"),
            1,
        );

        let mut fmi = FMI_search::ctor(prefix.to_str().expect("utf8"));
        assert_eq!(fmi.build_index(), 0);
        fmi.load_index();

        let (reference_seq, sa_bwt) = reference_and_sa(&prefix);
        let bwt = bwt_from_sa(&reference_seq, &sa_bwt);
        let enc_qdb = vec![0_u8, 1, 2, 3, 0, 1];
        let seqs = vec![bseq1_t {
            l_seq: 6,
            ..Default::default()
        }];
        let max_intv = vec![2_i32];
        let query_offsets = vec![0_i32];

        let mut actual = Vec::new();
        let num = fmi.bwtSeedStrategyAllPosOneThread(
            &enc_qdb,
            &max_intv,
            1,
            &seqs,
            &query_offsets,
            2,
            &mut actual,
        );

        let expected =
            explicit_seed_strategy(&fmi, &bwt, &enc_qdb, &max_intv, &seqs, &query_offsets, 2);
        assert_eq!(num, i64::try_from(expected.len()).expect("expected len"));
        assert_eq!(actual, expected);

        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn get_smems_one_pos_advances_queries_and_emits_matches() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_fmi_onepos_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("ref");
        let fasta = b">chr1\nACGTAC\n>chr2\nTTAA\n";
        let _ = bns_fasta2bntseq(
            Cursor::new(fasta.as_slice()),
            prefix.to_str().expect("utf8"),
            1,
        );

        let mut fmi = FMI_search::ctor(prefix.to_str().expect("utf8"));
        assert_eq!(fmi.build_index(), 0);
        fmi.load_index();

        let enc_qdb = vec![0_u8, 1, 2, 3, 0, 1];
        let seqs = vec![bseq1_t {
            l_seq: 6,
            ..Default::default()
        }];
        let query_offsets = vec![0_i32];
        let mut query_pos = vec![0_i16];
        let mut min_intv = vec![2_i32];
        let mut rid = vec![0_i32];
        let mut matches = Vec::new();
        let mut total = 0_i64;

        fmi.getSMEMsOnePosOneThread(
            &enc_qdb,
            &mut query_pos,
            &mut min_intv,
            &mut rid,
            1,
            1,
            &seqs,
            &query_offsets,
            6,
            2,
            &mut matches,
            &mut total,
        );

        assert_eq!(total, i64::try_from(matches.len()).expect("matches len"));
        assert!(!matches.is_empty());
        assert!(query_pos[0] > 0);
        for smem in &matches {
            assert_eq!(smem.rid, 0);
            assert!(smem.n >= smem.m);
            assert!(smem.s > 0);
        }

        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn get_smems_all_pos_matches_manual_one_pos_driver() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_fmi_allpos_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("ref");
        let fasta = b">chr1\nACGTAC\n>chr2\nTTAA\n";
        let _ = bns_fasta2bntseq(
            Cursor::new(fasta.as_slice()),
            prefix.to_str().expect("utf8"),
            1,
        );

        let mut fmi = FMI_search::ctor(prefix.to_str().expect("utf8"));
        assert_eq!(fmi.build_index(), 0);
        fmi.load_index();

        let enc_qdb = vec![0_u8, 1, 2, 3, 0, 1, 3, 3, 0, 0];
        let seqs = vec![
            bseq1_t {
                l_seq: 6,
                ..Default::default()
            },
            bseq1_t {
                l_seq: 4,
                ..Default::default()
            },
        ];
        let query_offsets = vec![0_i32, 6_i32];

        let mut all_min_intv = vec![2_i32, 2_i32];
        let mut all_rid = vec![0_i32, 1_i32];
        let mut all_matches = Vec::new();
        let mut all_total = 0_i64;
        fmi.getSMEMsAllPosOneThread(
            &enc_qdb,
            &mut all_min_intv,
            &mut all_rid,
            2,
            2,
            &seqs,
            &query_offsets,
            6,
            2,
            &mut all_matches,
            &mut all_total,
        );

        let mut manual_query_pos = vec![0_i16, 0_i16];
        let mut manual_min_intv = vec![2_i32, 2_i32];
        let mut manual_rid = vec![0_i32, 1_i32];
        let mut manual_total = 0_i64;
        let mut manual_matches = Vec::new();
        let mut num_active = 2_i32;
        loop {
            let mut tail = 0_usize;
            for head in 0..usize::try_from(num_active).expect("num_active") {
                let readlength = seqs[usize::try_from(manual_rid[head]).expect("rid")].l_seq;
                if i32::from(manual_query_pos[head]) < readlength {
                    manual_rid[tail] = manual_rid[head];
                    manual_query_pos[tail] = manual_query_pos[head];
                    manual_min_intv[tail] = manual_min_intv[head];
                    tail += 1;
                }
            }
            fmi.getSMEMsOnePosOneThread(
                &enc_qdb,
                &mut manual_query_pos[..tail],
                &mut manual_min_intv[..tail],
                &mut manual_rid[..tail],
                i32::try_from(tail).expect("tail"),
                2,
                &seqs,
                &query_offsets,
                6,
                2,
                &mut manual_matches,
                &mut manual_total,
            );
            num_active = i32::try_from(tail).expect("tail");
            if num_active <= 0 {
                break;
            }
        }

        assert_eq!(all_total, manual_total);
        assert_eq!(all_matches, manual_matches);

        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn get_smems_records_per_thread_counts_and_layout() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_fmi_getsmems_layout_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("ref");
        let fasta = b">chr1\nACGTAC\n>chr2\nTTAA\n";
        let _ = bns_fasta2bntseq(
            Cursor::new(fasta.as_slice()),
            prefix.to_str().expect("utf8"),
            1,
        );

        let mut fmi = FMI_search::ctor(prefix.to_str().expect("utf8"));
        assert_eq!(fmi.build_index(), 0);
        fmi.load_index();

        let readlength = 6_i32;
        let enc_qdb = vec![0_u8, 1, 2, 3, 0, 1, 3, 3, 0, 0, 0, 0];
        let mut matches = Vec::new();
        let mut totals = Vec::new();
        fmi.getSMEMs(&enc_qdb, 2, 2, readlength, 2, 2, &mut matches, &mut totals);

        assert_eq!(totals.len(), 2);
        let sum: i64 = totals.iter().sum();
        let mut populated = 0_i64;
        for tid in 0..2_usize {
            let per_thread_quota = (2 + (2 - 1)) / 2;
            let first = i32::try_from(tid).expect("tid") * per_thread_quota;
            let base = usize::try_from(first * readlength).expect("base");
            for idx in 0..usize::try_from(totals[tid]).expect("tid count") {
                assert!(matches[base + idx].s >= 0);
                populated += 1;
            }
        }
        assert_eq!(sum, populated);

        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn get_smems_is_stable_across_thread_counts_after_sort() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_fmi_getsmems_threads_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("ref");
        let fasta = b">chr1\nACGTAC\n>chr2\nTTAA\n";
        let _ = bns_fasta2bntseq(
            Cursor::new(fasta.as_slice()),
            prefix.to_str().expect("utf8"),
            1,
        );

        let mut fmi = FMI_search::ctor(prefix.to_str().expect("utf8"));
        assert_eq!(fmi.build_index(), 0);
        fmi.load_index();

        let readlength = 6_i32;
        let enc_qdb = vec![0_u8, 1, 2, 3, 0, 1, 3, 3, 0, 0, 0, 0];

        let mut matches1 = Vec::new();
        let mut totals1 = Vec::new();
        fmi.getSMEMs(
            &enc_qdb,
            2,
            2,
            readlength,
            2,
            1,
            &mut matches1,
            &mut totals1,
        );
        fmi.sortSMEMs(&mut matches1, &totals1, 2, readlength, 1);
        let total1: usize = totals1
            .iter()
            .map(|&x| usize::try_from(x).expect("x"))
            .sum();
        let out1 = matches1[..total1].to_vec();

        let mut matches2 = Vec::new();
        let mut totals2 = Vec::new();
        fmi.getSMEMs(
            &enc_qdb,
            2,
            2,
            readlength,
            2,
            2,
            &mut matches2,
            &mut totals2,
        );
        fmi.sortSMEMs(&mut matches2, &totals2, 2, readlength, 2);
        let mut out2 = Vec::new();
        for tid in 0..2_usize {
            let per_thread_quota = (2 + (2 - 1)) / 2;
            let first = i32::try_from(tid).expect("tid") * per_thread_quota;
            let base = usize::try_from(first * readlength).expect("base");
            let len = usize::try_from(totals2[tid]).expect("len");
            out2.extend_from_slice(&matches2[base..base + len]);
        }
        out2.sort_by(|a, b| super::compare_smem(a, b).cmp(&0));

        let mut out1_sorted = out1;
        out1_sorted.sort_by(|a, b| super::compare_smem(a, b).cmp(&0));
        assert_eq!(out1_sorted, out2);

        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn get_sa_entries_plain_smem_ranges_match_scalar_lookup() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_fmi_entries_smem_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("ref");
        let fasta = b">chr1\nACGTAC\n>chr2\nTTAA\n";
        let _ = bns_fasta2bntseq(
            Cursor::new(fasta.as_slice()),
            prefix.to_str().expect("utf8"),
            1,
        );

        let mut fmi = FMI_search::ctor(prefix.to_str().expect("utf8"));
        assert_eq!(fmi.build_index(), 0);
        fmi.load_index();

        let smems = vec![
            SMEM {
                k: 0,
                s: 2,
                ..Default::default()
            },
            SMEM {
                k: 2,
                s: 1,
                ..Default::default()
            },
        ];
        let mut coords = Vec::new();
        let mut coord_counts = Vec::new();
        fmi.get_sa_entries__L1077(
            &smems,
            &mut coords,
            &mut coord_counts,
            smems.len() as u32,
            3,
        );

        assert_eq!(coord_counts, vec![2, 1]);
        let expected = vec![
            fmi.get_sa_entry(0),
            fmi.get_sa_entry(1),
            fmi.get_sa_entry(2),
        ];
        assert_eq!(coords, expected);

        fs::remove_dir_all(&dir).expect("cleanup");
    }
}
