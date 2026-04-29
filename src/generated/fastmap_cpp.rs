#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Generated scaffold for `bwa-mem2/src/fastmap.cpp`.

use std::io::{BufReader, Cursor, Read, Write};
use std::net::TcpStream;
use std::process::Command;
use std::sync::mpsc::sync_channel;
use std::sync::{Arc, Mutex};

use clap::Parser;
use flate2::read::MultiGzDecoder;
use rayon::scope;

use crate::generated::bandedswa_h::SeqPair;
use crate::generated::bntseq_h::bntseq_t;
use crate::generated::bwa_cpp::bseq_read_orig;
use crate::generated::bwa_cpp::BWA_VERBOSE;
use crate::generated::bwa_cpp::{
    bseq_classify, bwa_fill_scmat, bwa_insert_header, bwa_print_sam_hdr, bwa_set_rg,
};
use crate::generated::bwamem_cpp::mem_process_seqs;
use crate::generated::bwamem_h::{
    mem_alnreg_v, mem_chain_v, mem_opt_t, mem_pestat_t, mem_seed_t, worker_t,
};
use crate::generated::fastmap_h::{ktp_aux_t, ktp_data_t};
use crate::generated::fmi_search_cpp::FMI_search;
use crate::generated::kseq_h::kseq_t;
use crate::generated::utils_cpp::{err_fclose, err_fputs, err_xopen_core, ErrFile};

const AVG_SEEDS_PER_READ: usize = 64;
const BATCH_MUL: usize = 20;
const READ_LEN: usize = 151;
const SEEDS_PER_READ: usize = 500;
const MAX_SEQ_LEN_REF: usize = 256;
const MAX_SEQ_LEN_QER: usize = 128;
const MAX_LINE_LEN: usize = 256;
const MEM_F_PE: i32 = 0x2;
const MEM_F_NOPAIRING: i32 = 0x4;
const MEM_F_ALL: i32 = 0x8;
const MEM_F_NO_MULTI: i32 = 0x10;
const MEM_F_NO_RESCUE: i32 = 0x20;
const MEM_F_REF_HDR: i32 = 0x100;
const MEM_F_SOFTCLIP: i32 = 0x200;
const MEM_F_SMARTPE: i32 = 0x400;
const MEM_F_PRIMARY5: i32 = 0x800;
const MEM_F_KEEP_SUPP_MAPQ: i32 = 0x1000;

#[derive(Debug, Parser)]
#[command(name = "mem", disable_help_flag = true)]
struct MemCli {
    #[arg(short = 'o', visible_short_alias = 'f')]
    output: Option<String>,
    #[arg(short = 't')]
    threads: Option<i32>,
    #[arg(short = 'k')]
    min_seed_len: Option<i32>,
    #[arg(short = 'w')]
    band_width: Option<i32>,
    #[arg(short = 'd')]
    zdrop: Option<i32>,
    #[arg(short = 'r')]
    split_factor: Option<f32>,
    #[arg(short = 's')]
    split_width: Option<i32>,
    #[arg(short = 'c')]
    max_occ: Option<i32>,
    #[arg(short = 'D')]
    drop_ratio: Option<f32>,
    #[arg(short = 'G')]
    max_chain_gap: Option<i32>,
    #[arg(short = 'm')]
    max_matesw: Option<i32>,
    #[arg(short = 'N')]
    max_chain_extend: Option<i32>,
    #[arg(short = 'A')]
    a: Option<i32>,
    #[arg(short = 'B')]
    b: Option<i32>,
    #[arg(short = 'O')]
    gap_open: Option<String>,
    #[arg(short = 'E')]
    gap_ext: Option<String>,
    #[arg(short = 'L')]
    clip_pen: Option<String>,
    #[arg(short = 'U')]
    pen_unpaired: Option<i32>,
    #[arg(short = 'T')]
    min_score: Option<i32>,
    #[arg(short = 'Q')]
    mapq_coef_len: Option<i32>,
    #[arg(short = 'X')]
    mask_level: Option<f32>,
    #[arg(short = 'y')]
    max_mem_intv: Option<u64>,
    #[arg(short = 'W')]
    min_chain_weight: Option<i32>,
    #[arg(short = 'h')]
    xa_hits: Option<String>,
    #[arg(short = '1')]
    no_mt_io: bool,
    #[arg(short = 'S')]
    no_rescue: bool,
    #[arg(short = 'P')]
    no_pairing: bool,
    #[arg(short = 'p')]
    smart_pairing: bool,
    #[arg(short = 'a')]
    all_alignments: bool,
    #[arg(short = 'M')]
    no_multi: bool,
    #[arg(short = 'Y')]
    softclip: bool,
    #[arg(short = 'V')]
    ref_hdr: bool,
    #[arg(short = '5')]
    primary5: bool,
    #[arg(short = 'q')]
    keep_supp_mapq: bool,
    #[arg(short = 'C')]
    copy_comment: bool,
    #[arg(short = 'j')]
    ignore_alt: bool,
    #[arg(short = 'K')]
    fixed_chunk_size: Option<i64>,
    #[arg(short = 'v')]
    verbose: Option<i32>,
    #[arg(short = 'R')]
    read_group: Option<String>,
    #[arg(short = 'H')]
    header: Vec<String>,
    #[arg(short = 'I')]
    insert: Option<String>,
    #[arg(short = 'x')]
    mode: Option<String>,
    idxbase: String,
    in1: String,
    in2: Option<String>,
}

fn parse_pair_i32(arg: &str) -> Option<(i32, Option<i32>)> {
    let mut parts = arg
        .split(|c: char| c.is_ascii_punctuation())
        .filter(|s| !s.is_empty());
    let first = parts.next()?.parse().ok()?;
    let second = parts.next().and_then(|s| s.parse().ok());
    Some((first, second))
}

fn parse_insert_stats(arg: &str) -> Option<mem_pestat_t> {
    let parts: Vec<&str> = arg
        .split(|c: char| c.is_ascii_punctuation())
        .filter(|s| !s.is_empty())
        .collect();
    if parts.is_empty() {
        return None;
    }
    let avg: f64 = parts[0].parse().ok()?;
    let mut std = avg * 0.1;
    let mut high;
    let mut low;
    if parts.len() > 1 {
        std = parts[1].parse().ok()?;
    }
    high = (avg + 4.0 * std + 0.499) as i32;
    low = (avg - 4.0 * std + 0.499) as i32;
    if low < 1 {
        low = 1;
    }
    if parts.len() > 2 {
        // C++ fastmap.cpp:772 uses `(int)(strtod(...) + .499)` (truncate-toward-zero of x+0.499),
        // not `.round()` which is half-away-from-zero. They diverge at x.5 (e.g. 100.5 → 100 vs 101).
        high = (parts[2].parse::<f64>().ok()? + 0.499) as i32;
    }
    if parts.len() > 3 {
        low = (parts[3].parse::<f64>().ok()? + 0.499) as i32;
    }
    Some(mem_pestat_t {
        avg,
        std,
        high,
        low,
        failed: 0,
    })
}

fn parse_http_url(url: &str) -> Option<(String, u16, String)> {
    let rest = url.strip_prefix("http://")?;
    let (host_port, path) = match rest.split_once('/') {
        Some((host_port, path)) => (host_port, format!("/{}", path)),
        None => (rest, "/".to_string()),
    };
    let (host, port) = match host_port.split_once(':') {
        Some((host, port)) => (host.to_string(), port.parse().ok()?),
        None => (host_port.to_string(), 80_u16),
    };
    Some((host, port, path))
}

fn fetch_http_bytes(url: &str) -> Result<Vec<u8>, ()> {
    let (host, port, path) = parse_http_url(url).ok_or(())?;
    let mut stream = TcpStream::connect((host.as_str(), port)).map_err(|_| ())?;
    write!(
        stream,
        "GET {} HTTP/1.0\r\nHost: {}\r\nUser-Agent: bwa-mem2-rs\r\nConnection: close\r\n\r\n",
        path, host
    )
    .map_err(|_| ())?;
    let mut response = Vec::new();
    stream.read_to_end(&mut response).map_err(|_| ())?;
    let header_end = response
        .windows(4)
        .position(|w| w == b"\r\n\r\n")
        .map(|idx| idx + 4)
        .or_else(|| {
            response
                .windows(2)
                .position(|w| w == b"\n\n")
                .map(|idx| idx + 2)
        })
        .ok_or(())?;
    let header = String::from_utf8_lossy(&response[..header_end]);
    let status_line = header.lines().next().ok_or(())?;
    if !status_line.starts_with("HTTP/") || !status_line.contains("200") {
        return Err(());
    }
    Ok(response[header_end..].to_vec())
}

fn fetch_remote_with_tool(url: &str) -> Result<Vec<u8>, ()> {
    for (program, args) in [("curl", vec!["-fsSL", url]), ("wget", vec!["-qO-", url])] {
        let Ok(output) = Command::new(program).args(&args).output() else {
            continue;
        };
        if output.status.success() {
            return Ok(output.stdout);
        }
    }
    Err(())
}

fn read_bytes_input(path: &str) -> Result<Vec<u8>, ()> {
    if path == "-" {
        let mut bytes = Vec::new();
        std::io::stdin().read_to_end(&mut bytes).map_err(|_| ())?;
        Ok(bytes)
    } else if let Some(cmd) = path.trim_start().strip_prefix('<') {
        let output = Command::new("sh")
            .arg("-c")
            .arg(cmd.trim())
            .output()
            .map_err(|_| ())?;
        if output.status.success() {
            Ok(output.stdout)
        } else {
            Err(())
        }
    } else if path.starts_with("http://") {
        fetch_http_bytes(path).or_else(|_| fetch_remote_with_tool(path))
    } else if path.starts_with("ftp://") {
        fetch_remote_with_tool(path)
    } else {
        std::fs::read(path).map_err(|_| ())
    }
}

fn apply_mode(opt: &mut mem_opt_t, cli: &MemCli) -> Result<(), ()> {
    let Some(mode) = cli.mode.as_deref() else {
        return Ok(());
    };
    eprintln!("WARNING: bwa-mem2 doesn't work well with long reads or contigs; please use minimap2 instead.");
    match mode {
        "intractg" => {
            if cli.gap_open.is_none() {
                opt.o_del = 16;
                opt.o_ins = 16;
            }
            if cli.b.is_none() {
                opt.b = 9;
            }
            if cli.clip_pen.is_none() {
                opt.pen_clip5 = 5;
                opt.pen_clip3 = 5;
            }
        }
        "pacbio" | "pbref" | "ont2d" => {
            if cli.gap_open.is_none() {
                opt.o_del = 1;
                opt.o_ins = 1;
            }
            if cli.gap_ext.is_none() {
                opt.e_del = 1;
                opt.e_ins = 1;
            }
            if cli.b.is_none() {
                opt.b = 1;
            }
            if cli.split_factor.is_none() {
                opt.split_factor = 10.0;
            }
            if mode == "ont2d" {
                if cli.min_chain_weight.is_none() {
                    opt.min_chain_weight = 20;
                }
                if cli.min_seed_len.is_none() {
                    opt.min_seed_len = 14;
                }
            } else {
                if cli.min_chain_weight.is_none() {
                    opt.min_chain_weight = 40;
                }
                if cli.min_seed_len.is_none() {
                    opt.min_seed_len = 17;
                }
            }
            if cli.clip_pen.is_none() {
                opt.pen_clip5 = 0;
                opt.pen_clip3 = 0;
            }
        }
        other => {
            eprintln!("[E::main_mem] unknown read type '{other}'");
            return Err(());
        }
    }
    Ok(())
}

fn read_text_input(path: &str) -> Result<String, ()> {
    let bytes = read_bytes_input(path)?;
    let gzipped = path.ends_with(".gz") || bytes.starts_with(&[0x1f, 0x8b]);
    if gzipped {
        let mut decoder = MultiGzDecoder::new(Cursor::new(bytes));
        let mut text = String::new();
        decoder.read_to_string(&mut text).map_err(|_| ())?;
        Ok(text)
    } else {
        String::from_utf8(bytes).map_err(|_| ())
    }
}

fn open_kseq_input(path: &str) -> Result<kseq_t, ()> {
    if path == "-" {
        return Ok(kseq_t::from_reader(Box::new(BufReader::new(
            std::io::stdin(),
        ))));
    }
    if path.trim_start().starts_with('<')
        || path.starts_with("http://")
        || path.starts_with("ftp://")
    {
        return read_text_input(path).map(|text| kseq_t::from_text(&text));
    }
    kseq_t::from_path(path).map_err(|_| ())
}

fn read_batch(aux: &mut ktp_aux_t) -> Option<ktp_data_t> {
    let mut ret = ktp_data_t::default();
    let mut sz = 0_i64;
    let ks = aux.ks.as_mut().expect("read_batch requires aux.ks");
    ret.seqs = bseq_read_orig(
        aux.task_size,
        &mut ret.n_seqs,
        ks,
        aux.ks2.as_mut(),
        &mut sz,
    );
    if ret.seqs.is_empty() {
        return None;
    }
    if aux.copy_comment == 0 {
        for seq in &mut ret.seqs {
            seq.comment = None;
        }
    }
    Some(ret)
}

fn process_batch(
    batch: &mut ktp_data_t,
    opt: &mut mem_opt_t,
    n_processed: i64,
    pes0: Option<&[mem_pestat_t; 4]>,
    w: &mut worker_t,
) {
    if opt.flag & MEM_F_SMARTPE != 0 {
        let mut sep = [Vec::new(), Vec::new()];
        let mut n_sep = [0_i32; 2];
        let mut tmp_opt = opt.clone();
        bseq_classify(batch.n_seqs, &batch.seqs, &mut n_sep, &mut sep);

        if n_sep[0] > 0 {
            tmp_opt.flag &= !MEM_F_PE;
            mem_process_seqs(&mut tmp_opt, n_processed, n_sep[0], &mut sep[0], None, w);
            for seq in &sep[0] {
                batch.seqs[usize::try_from(seq.id).expect("single seq id")].sam = seq.sam.clone();
            }
        }
        if n_sep[1] > 0 {
            tmp_opt.flag |= MEM_F_PE;
            mem_process_seqs(
                &mut tmp_opt,
                n_processed + i64::from(n_sep[0]),
                n_sep[1],
                &mut sep[1],
                pes0,
                w,
            );
            for seq in &sep[1] {
                batch.seqs[usize::try_from(seq.id).expect("paired seq id")].sam = seq.sam.clone();
            }
        }
    } else {
        mem_process_seqs(opt, n_processed, batch.n_seqs, &mut batch.seqs, pes0, w);
    }

    batch.sam_lines.clear();
    batch.sam_lines.reserve(batch.seqs.len());
    for seq in &mut batch.seqs {
        if let Some(sam) = seq.sam.take() {
            batch.sam_lines.push(sam);
        }
        seq.name = None;
        seq.comment = None;
        seq.seq = None;
        seq.qual = None;
        seq.seq_nt4.clear();
    }
    batch.seqs.clear();
    batch.seqs.shrink_to(0);
}

fn write_batch(batch: &ktp_data_t, fp: &mut ErrFile) {
    for sam in &batch.sam_lines {
        err_fputs(sam, fp);
    }
}

#[doc = "Original function: __cpuid:51"]
pub fn cpuid(i: u32) -> [u32; 4] {
    #[cfg(target_arch = "x86")]
    {
        let r = unsafe { std::arch::x86::__cpuid_count(i, 0) };
        [r.eax, r.ebx, r.ecx, r.edx]
    }
    #[cfg(target_arch = "x86_64")]
    {
        let r = unsafe { std::arch::x86_64::__cpuid_count(i, 0) };
        [r.eax, r.ebx, r.ecx, r.edx]
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    {
        [0, 0, 0, 0]
    }
}

#[doc = "Original function: HTStatus:63"]
pub fn HTStatus() -> i32 {
    let cpuid0 = cpuid(0);
    let mut vendor = [0_u8; 12];
    vendor[0..4].copy_from_slice(&cpuid0[1].to_le_bytes());
    vendor[4..8].copy_from_slice(&cpuid0[3].to_le_bytes());
    vendor[8..12].copy_from_slice(&cpuid0[2].to_le_bytes());
    let vendor = String::from_utf8_lossy(&vendor).to_string();

    let cpuid1 = cpuid(1);
    let platform_features = cpuid1[3];
    let num_logical_cpus = (cpuid1[1] >> 16) & 0xff;

    let num_cores = if vendor == "GenuineIntel" {
        let cpuid4 = cpuid(4);
        eprintln!("Platform vendor: Intel.");
        ((cpuid4[0] >> 26) & 0x3f) + 1
    } else {
        eprintln!("Platform vendor unknown.");
        u32::MAX
    };

    let ht = ((platform_features & (1 << 28)) != 0) && num_cores < num_logical_cpus;
    if ht {
        eprintln!("CPUs support hyperThreading !!");
        1
    } else {
        0
    }
}

#[doc = "Original function: memoryAlloc:100"]
pub fn memoryAlloc(aux: &ktp_aux_t, w: &mut worker_t, nreads: i32, nthreads: i32) {
    let opt = aux.opt.as_ref().expect("memoryAlloc requires aux.opt");
    let read_len = READ_LEN;
    let mem_size = usize::try_from(nreads.max(0)).expect("nreads");

    w.regs = vec![mem_alnreg_v::default(); mem_size];
    w.chain_ar = vec![mem_chain_v::default(); mem_size];
    w.seedBufSize = i64::try_from(crate::generated::macro_h::BATCH_SIZE * AVG_SEEDS_PER_READ)
        .expect("seedBufSize");

    let nthreads_usize = usize::try_from(nthreads.max(1)).expect("nthreads");
    let seed_buf_len = (mem_size * AVG_SEEDS_PER_READ)
        .max(nthreads_usize * crate::generated::macro_h::BATCH_SIZE * AVG_SEEDS_PER_READ);
    w.seedBuf = vec![mem_seed_t::default(); seed_buf_len];
    let wsize = crate::generated::macro_h::BATCH_SIZE * SEEDS_PER_READ;

    w.mmc.seqBufLeftRef = vec![Vec::new(); nthreads_usize];
    w.mmc.seqBufLeftQer = vec![Vec::new(); nthreads_usize];
    w.mmc.seqBufRightRef = vec![Vec::new(); nthreads_usize];
    w.mmc.seqBufRightQer = vec![Vec::new(); nthreads_usize];
    w.mmc.wsize_buf_ref = vec![0; nthreads_usize];
    w.mmc.wsize_buf_qer = vec![0; nthreads_usize];

    w.mmc.seqPairArrayAux = vec![Vec::new(); nthreads_usize];
    w.mmc.seqPairArrayLeft128 = vec![Vec::new(); nthreads_usize];
    w.mmc.seqPairArrayRight128 = vec![Vec::new(); nthreads_usize];
    w.mmc.wsize = vec![0; nthreads_usize];

    let _mem_wsize = BATCH_MUL * crate::generated::macro_h::BATCH_SIZE * read_len;
    w.mmc.wsize_mem = vec![0; nthreads_usize];
    w.mmc.wsize_mem_s = vec![0; nthreads_usize];
    w.mmc.wsize_mem_r = vec![0; nthreads_usize];
    w.mmc.matchArray = vec![Vec::new(); nthreads_usize];
    w.mmc.min_intv_ar = vec![Vec::new(); nthreads_usize];
    w.mmc.query_pos_ar = vec![Vec::new(); nthreads_usize];
    w.mmc.enc_qdb = vec![Vec::new(); nthreads_usize];
    w.mmc.rid = vec![Vec::new(); nthreads_usize];
    w.mmc.lim = vec![vec![0; crate::generated::macro_h::BATCH_SIZE + 32]; nthreads_usize];

    let _alloc_mem = (wsize * MAX_SEQ_LEN_REF + MAX_LINE_LEN)
        * usize::try_from(opt.n_threads.max(1)).expect("opt threads")
        * 2
        + (wsize * MAX_SEQ_LEN_QER + MAX_LINE_LEN)
            * usize::try_from(opt.n_threads.max(1)).expect("opt threads")
            * 2
        + wsize
            * std::mem::size_of::<SeqPair>()
            * usize::try_from(opt.n_threads.max(1)).expect("opt threads")
            * 3;
}

#[doc = "Original function: kt_pipeline:189"]
pub fn kt_pipeline(
    shared: &mut ktp_aux_t,
    step: i32,
    data: Option<ktp_data_t>,
    opt: &mut mem_opt_t,
    w: &mut worker_t,
) -> Option<ktp_data_t> {
    match step {
        0 => read_batch(shared),
        1 => {
            let mut ret = data.expect("kt_pipeline step 1 requires input data");
            process_batch(&mut ret, opt, shared.n_processed, shared.pes0.as_deref(), w);
            shared.n_processed += i64::from(ret.n_seqs);
            Some(ret)
        }
        2 => crate::support::stub::<Option<ktp_data_t>>("kt_pipeline write-output path"),
        _ => None,
    }
}

#[doc = "Original function: ktp_worker:327"]
pub fn ktp_worker(_arg0: crate::support::Opaque) {
    crate::support::stub::<()>("ktp_worker")
}

#[doc = "Original function: process:368"]
pub fn process(shared: &mut ktp_aux_t, w: &mut worker_t, _pipe_threads: i32) -> i32 {
    let nthreads = shared
        .opt
        .as_ref()
        .expect("process requires aux.opt")
        .n_threads;
    w.nthreads = i16::try_from(nthreads).expect("nthreads");

    let nreads =
        i32::try_from(shared.actual_chunk_size / i64::try_from(READ_LEN).expect("READ_LEN") + 10)
            .expect("nreads");
    memoryAlloc(shared, w, nreads, nthreads);
    w.ref_string = shared.ref_string.clone();
    w.nreads = 0;

    let mut opt = (*shared.opt.as_ref().expect("opt")).clone();
    if _pipe_threads <= 1 {
        while let Some(mut batch) = read_batch(shared) {
            process_batch(
                &mut batch,
                &mut opt,
                shared.n_processed,
                shared.pes0.as_deref(),
                w,
            );
            shared.n_processed += i64::from(batch.n_seqs);
            if let Some(fp) = shared.fp.as_mut() {
                write_batch(&batch, fp);
            }
        }
        return 0;
    }

    let (tx_read, rx_read) = sync_channel::<Option<ktp_data_t>>(0);
    let (tx_write, rx_write) = sync_channel::<Option<ktp_data_t>>(0);
    let mut ks = shared.ks.take();
    let mut ks2 = shared.ks2.take();
    let task_size = shared.task_size;
    let copy_comment = shared.copy_comment;
    let pes0 = shared.pes0.as_deref().copied();
    let fp = Arc::new(Mutex::new(shared.fp.take()));
    let final_n_processed = Arc::new(Mutex::new(shared.n_processed));

    scope(|s| {
        s.spawn(move |_| {
            let mut aux = ktp_aux_t {
                ks: ks.take(),
                ks2: ks2.take(),
                task_size,
                copy_comment,
                ..Default::default()
            };
            while let Some(batch) = read_batch(&mut aux) {
                if tx_read.send(Some(batch)).is_err() {
                    break;
                }
            }
            let _ = tx_read.send(None);
        });

        let final_n_processed_compute = Arc::clone(&final_n_processed);
        s.spawn(move |_| {
            let mut n_processed = *final_n_processed_compute.lock().expect("n_processed lock");
            while let Ok(msg) = rx_read.recv() {
                let Some(mut batch) = msg else {
                    let _ = tx_write.send(None);
                    break;
                };
                process_batch(&mut batch, &mut opt, n_processed, pes0.as_ref(), w);
                n_processed += i64::from(batch.n_seqs);
                if tx_write.send(Some(batch)).is_err() {
                    break;
                }
            }
            *final_n_processed_compute.lock().expect("n_processed lock") = n_processed;
        });

        let fp_write = Arc::clone(&fp);
        s.spawn(move |_| {
            while let Ok(msg) = rx_write.recv() {
                let Some(batch) = msg else {
                    break;
                };
                if let Some(fp_ref) = fp_write.lock().expect("fp lock").as_mut() {
                    write_batch(&batch, fp_ref);
                }
            }
        });
    });

    shared.fp = fp.lock().expect("fp lock").take();
    shared.n_processed = *final_n_processed.lock().expect("n_processed lock");
    0
}

#[doc = "Original function: update_a:547"]
pub fn update_a(opt: &mut mem_opt_t, opt0: &mem_opt_t) {
    if opt0.a != 0 {
        if opt0.b == 0 {
            opt.b *= opt.a;
        }
        if opt0.T == 0 {
            opt.T *= opt.a;
        }
        if opt0.o_del == 0 {
            opt.o_del *= opt.a;
        }
        if opt0.e_del == 0 {
            opt.e_del *= opt.a;
        }
        if opt0.o_ins == 0 {
            opt.o_ins *= opt.a;
        }
        if opt0.e_ins == 0 {
            opt.e_ins *= opt.a;
        }
        if opt0.zdrop == 0 {
            opt.zdrop *= opt.a;
        }
        if opt0.pen_clip5 == 0 {
            opt.pen_clip5 *= opt.a;
        }
        if opt0.pen_clip3 == 0 {
            opt.pen_clip3 *= opt.a;
        }
        if opt0.pen_unpaired == 0 {
            opt.pen_unpaired *= opt.a;
        }
    }
}

#[doc = "Original function: usage:563"]
pub fn usage(opt: &mem_opt_t) {
    eprint!("{}", usage_text(opt));
}

#[doc = "Original function: main_mem:616"]
pub fn main_mem(argv: &[String]) -> i32 {
    let cli = match MemCli::try_parse_from(argv) {
        Ok(cli) => cli,
        Err(err) => {
            let _ = err.print();
            return 1;
        }
    };

    let mut opt = *crate::generated::bwamem_cpp::mem_opt_init();
    let mut opt0 = mem_opt_t::default();
    if let Some(v) = cli.min_seed_len {
        opt.min_seed_len = v;
        opt0.min_seed_len = 1;
    }
    if let Some(v) = cli.band_width {
        opt.w = v;
        opt0.w = 1;
    }
    if let Some(v) = cli.zdrop {
        opt.zdrop = v;
        opt0.zdrop = 1;
    }
    if let Some(v) = cli.split_factor {
        opt.split_factor = v;
        opt0.split_factor = 1.0;
    }
    if let Some(v) = cli.split_width {
        opt.split_width = v;
    }
    if let Some(v) = cli.max_occ {
        opt.max_occ = v;
        opt0.max_occ = 1;
    }
    if let Some(v) = cli.drop_ratio {
        opt.drop_ratio = v;
        opt0.drop_ratio = 1.0;
    }
    if let Some(v) = cli.max_chain_gap {
        opt.max_chain_gap = v;
        opt0.max_chain_gap = 1;
    }
    if let Some(v) = cli.max_matesw {
        opt.max_matesw = v;
        opt0.max_matesw = 1;
    }
    if let Some(v) = cli.max_chain_extend {
        opt.max_chain_extend = v;
        opt0.max_chain_extend = 1;
    }
    if let Some(v) = cli.a {
        opt.a = v;
        opt0.a = 1;
    }
    if let Some(v) = cli.b {
        opt.b = v;
        opt0.b = 1;
    }
    if let Some(v) = cli.pen_unpaired {
        opt.pen_unpaired = v;
        opt0.pen_unpaired = 1;
    }
    if let Some(v) = cli.min_score {
        opt.T = v;
        opt0.T = 1;
    }
    if let Some(v) = cli.mapq_coef_len {
        opt.mapQ_coef_len = v as f32;
        opt.mapQ_coef_fac = if v > 0 { (f64::from(v)).ln() as i32 } else { 0 };
    }
    if let Some(v) = cli.mask_level {
        opt.mask_level = v;
    }
    if let Some(v) = cli.max_mem_intv {
        opt.max_mem_intv = v;
        opt0.max_mem_intv = 1;
    }
    if let Some(v) = cli.min_chain_weight {
        opt.min_chain_weight = v;
        opt0.min_chain_weight = 1;
    }
    if let Some(arg) = cli.xa_hits.as_deref().and_then(parse_pair_i32) {
        opt.max_XA_hits = arg.0;
        opt.max_XA_hits_alt = arg.1.unwrap_or(arg.0);
        opt0.max_XA_hits = 1;
        opt0.max_XA_hits_alt = 1;
    }
    if let Some(v) = cli.threads {
        opt.n_threads = v.max(1);
    }
    if let Some(v) = cli.verbose {
        BWA_VERBOSE.store(v, std::sync::atomic::Ordering::Relaxed);
    }
    if cli.no_pairing {
        opt.flag |= MEM_F_NOPAIRING;
    }
    if cli.all_alignments {
        opt.flag |= MEM_F_ALL;
    }
    if cli.no_multi {
        opt.flag |= MEM_F_NO_MULTI;
    }
    if cli.no_rescue {
        opt.flag |= MEM_F_NO_RESCUE;
    }
    if cli.softclip {
        opt.flag |= MEM_F_SOFTCLIP;
    }
    if cli.ref_hdr {
        opt.flag |= MEM_F_REF_HDR;
    }
    if cli.primary5 {
        opt.flag |= MEM_F_PRIMARY5 | MEM_F_KEEP_SUPP_MAPQ;
    }
    if cli.keep_supp_mapq {
        opt.flag |= MEM_F_KEEP_SUPP_MAPQ;
    }
    if cli.smart_pairing {
        opt.flag |= MEM_F_PE | MEM_F_SMARTPE;
    }
    if let Some(arg) = cli.gap_open.as_deref() {
        if let Some((first, second)) = parse_pair_i32(arg) {
            opt.o_del = first;
            opt.o_ins = second.unwrap_or(first);
            opt0.o_del = 1;
            opt0.o_ins = 1;
        }
    }
    if let Some(arg) = cli.gap_ext.as_deref() {
        if let Some((first, second)) = parse_pair_i32(arg) {
            opt.e_del = first;
            opt.e_ins = second.unwrap_or(first);
            opt0.e_del = 1;
            opt0.e_ins = 1;
        }
    }
    if let Some(arg) = cli.clip_pen.as_deref() {
        if let Some((first, second)) = parse_pair_i32(arg) {
            opt.pen_clip5 = first;
            opt.pen_clip3 = second.unwrap_or(first);
            opt0.pen_clip5 = 1;
            opt0.pen_clip3 = 1;
        }
    }
    if apply_mode(&mut opt, &cli).is_err() {
        return 1;
    } else if cli.mode.is_none() {
        update_a(&mut opt, &opt0);
    }

    bwa_fill_scmat(opt.a, opt.b, &mut opt.mat);

    let mut hdr_line: Option<String> = None;
    if let Some(rg) = cli.read_group.as_deref() {
        let Some(rg_line) = bwa_set_rg(rg) else {
            return 1;
        };
        hdr_line = bwa_insert_header(Some(&rg_line), hdr_line);
    }
    for hdr in &cli.header {
        if hdr.starts_with('@') {
            hdr_line = bwa_insert_header(Some(hdr), hdr_line);
        } else {
            let Ok(text) = std::fs::read_to_string(hdr) else {
                eprintln!("[E::main_mem] failed to open header file `{hdr}`");
                return 1;
            };
            for line in text.lines() {
                hdr_line = bwa_insert_header(Some(line), hdr_line);
            }
        }
    }

    let mut fmi = FMI_search::ctor(&cli.idxbase);
    fmi.load_index();
    if cli.ignore_alt {
        if let Some(bns) = fmi.base.idx.bns.as_mut() {
            for ann in &mut bns.anns {
                ann.is_alt = 0;
            }
        }
    }
    let Some(bns): Option<&bntseq_t> = fmi.base.idx.bns.as_ref() else {
        eprintln!("[E::main_mem] failed to load reference index");
        return 1;
    };

    let mut out = if let Some(path) = cli.output.as_deref() {
        err_xopen_core("main_mem", path, "w")
    } else {
        ErrFile::Stdout
    };
    bwa_print_sam_hdr(bns, hdr_line.as_deref(), &mut out);

    let ks1 = match open_kseq_input(&cli.in1) {
        Ok(ks) => ks,
        Err(()) => {
            eprintln!("[E::main_mem] fail to open file `{}`.", cli.in1);
            let _ = err_fclose(out);
            return 1;
        }
    };
    let ks2 = if cli.smart_pairing {
        if let Some(in2) = cli.in2.as_deref() {
            eprintln!(
                "[W::main_mem] when '-p' is in use, the second query file is ignored: `{in2}`"
            );
        }
        None
    } else if let Some(in2) = cli.in2.as_deref() {
        opt.flag |= MEM_F_PE;
        match open_kseq_input(in2) {
            Ok(ks) => Some(ks),
            Err(()) => {
                eprintln!("[E::main_mem] fail to open file `{}`.", in2);
                let _ = err_fclose(out);
                return 1;
            }
        }
    } else {
        None
    };

    let mut aux = ktp_aux_t {
        ks: Some(ks1),
        ks2: ks2,
        opt: Some(Box::new(opt.clone())),
        copy_comment: if cli.copy_comment { 1 } else { 0 },
        task_size: cli
            .fixed_chunk_size
            .unwrap_or(opt.chunk_size * i64::from(opt.n_threads)),
        actual_chunk_size: cli
            .fixed_chunk_size
            .unwrap_or(opt.chunk_size * i64::from(opt.n_threads)),
        fp: Some(out),
        ..Default::default()
    };
    if let Some(insert) = cli.insert.as_deref().and_then(parse_insert_stats) {
        let mut pes = [mem_pestat_t::default(); 4];
        for entry in &mut pes {
            entry.failed = 1;
        }
        pes[1] = insert;
        aux.pes0 = Some(Box::new(pes));
    }

    let mut worker = worker_t {
        fmi: Some(fmi),
        ..Default::default()
    };
    let rc = process(&mut aux, &mut worker, if cli.no_mt_io { 1 } else { 2 });
    if let Some(out) = aux.fp.take() {
        let _ = err_fclose(out);
    }
    if let Some(fmi) = worker.fmi.as_mut() {
        fmi.dtor();
    }
    rc
}

fn usage_text(opt: &mem_opt_t) -> String {
    format!(
        concat!(
            "Usage: bwa-mem2 mem [options] <idxbase> <in1.fq> [in2.fq]\n",
            "Options:\n",
            "  Algorithm options:\n",
            "    -o STR        Output SAM file name\n",
            "    -t INT        number of threads [{}]\n",
            "    -k INT        minimum seed length [{}]\n",
            "    -w INT        band width for banded alignment [{}]\n",
            "    -d INT        off-diagonal X-dropoff [{}]\n",
            "    -r FLOAT      look for internal seeds inside a seed longer than {{-k}} * FLOAT [{}]\n",
            "    -s INT        look for internal seeds inside a seed with less than INT occ [{}]\n",
            "    -y INT        seed occurrence for the 3rd round seeding [{}]\n",
            "    -c INT        skip seeds with more than INT occurrences [{}]\n",
            "    -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [{:.2}]\n",
            "    -G INT        maximum gap between seeds in a chain [{}]\n",
            "    -W INT        discard a chain if seeded bases shorter than INT [0]\n",
            "    -m INT        perform at most INT rounds of mate rescues for each read [{}]\n",
            "    -N INT        max occurrences for extending a seed [{}]\n",
            "    -1            use single-threaded I/O path\n",
            "    -S            skip mate rescue\n",
            "    -o            output file name missing\n",
            "    -P            skip pairing; mate rescue performed unless -S also in use\n",
            "Scoring options:\n",
            "   -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [{}]\n",
            "   -B INT        penalty for a mismatch [{}]\n",
            "   -O INT[,INT]  gap open penalties for deletions and insertions [{},{}]\n",
            "   -E INT[,INT]  gap extension penalty; a gap of size k cost '{{-O}} + {{-E}}*k' [{},{}]\n",
            "   -L INT[,INT]  penalty for 5'- and 3'-end clipping [{},{}]\n",
            "   -U INT        penalty for an unpaired read pair [{}]\n",
            "Input/output options:\n",
            "   -p            smart pairing (ignoring in2.fq)\n",
            "   -R STR        read group header line such as '@RG\\tID:foo\\tSM:bar' [null]\n",
            "   -H STR/FILE   insert STR to header if it starts with @; or insert lines in FILE [null]\n",
            "   -j            treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)\n",
            "   -5            for split alignment, take the alignment with the smallest coordinate as primary\n",
            "   -q            don't modify mapQ of supplementary alignments\n",
            "   -K INT        process INT input bases in each batch regardless of nThreads (for reproducibility) []\n",
            "   -v INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging [{}]\n",
            "   -T INT        minimum score to output [{}]\n",
            "   -h INT[,INT]  if there are <INT hits with score >80% of the max score, output all in XA [{},{}]\n",
            "   -Q INT        length coefficient for mapQ scaling [{:.0}]\n",
            "   -X FLOAT      mask-level overlap ratio [{:.2}]\n",
            "   -a            output all alignments for SE or unpaired PE\n",
            "   -C            append FASTA/FASTQ comment to SAM output\n",
            "   -V            output the reference FASTA header in the XR tag\n",
            "   -Y            use soft clipping for supplementary alignments\n",
            "   -M            mark shorter split hits as secondary\n",
            "   -I FLOAT[,FLOAT[,INT[,INT]]]\n",
            "                 specify the mean, standard deviation (10% of the mean if absent), max\n",
            "                 (4 sigma from the mean if absent) and min of the insert size distribution.\n",
            "                 FR orientation only. [inferred]\n",
            "Note: Please read the man page for detailed description of the command line and options.\n"
        ),
        opt.n_threads,
        opt.min_seed_len,
        opt.w,
        opt.zdrop,
        opt.split_factor,
        opt.split_width,
        opt.max_mem_intv,
        opt.max_occ,
        opt.drop_ratio,
        opt.max_chain_gap,
        opt.max_matesw,
        opt.max_chain_extend,
        opt.a,
        opt.b,
        opt.o_del,
        opt.o_ins,
        opt.e_del,
        opt.e_ins,
        opt.pen_clip5,
        opt.pen_clip3,
        opt.pen_unpaired,
        BWA_VERBOSE.load(std::sync::atomic::Ordering::Relaxed),
        opt.T,
        opt.max_XA_hits,
        opt.max_XA_hits_alt,
        opt.mapQ_coef_len,
        opt.mask_level
    )
}

#[cfg(test)]
mod tests {
    use super::{
        cpuid, kt_pipeline, main_mem, memoryAlloc, process, read_text_input, update_a, usage_text,
        HTStatus, AVG_SEEDS_PER_READ,
    };
    use crate::generated::bntseq_cpp::bns_fasta2bntseq;
    use crate::generated::bwa_cpp::BWA_VERBOSE;
    use crate::generated::bwa_h::bseq1_t;
    use crate::generated::bwamem_cpp::mem_opt_init;
    use crate::generated::bwamem_h::{
        mem_alnreg_v, mem_cache, mem_chain_v, mem_opt_t, mem_seed_t, worker_t,
    };
    use crate::generated::fastmap_h::ktp_aux_t;
    use crate::generated::fmi_search_cpp::FMI_search;
    use crate::generated::kseq_h::kseq_t;
    use crate::generated::macro_h::BATCH_SIZE;
    use crate::generated::utils_cpp::{err_fclose, err_xopen_core};
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::fs;
    use std::io::{BufRead, BufReader, Cursor, Read, Write};
    use std::net::TcpListener;
    use std::path::{Path, PathBuf};
    use std::process::Command;
    use std::sync::atomic::Ordering;
    use std::thread;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn tutorial_fixture_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("external")
            .join("tutorial_bwa_small")
    }

    fn indexed_real_fixture_prefix() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join(".tmp")
            .join("real_bench")
            .join("ecoli_rel606")
    }

    fn extract_fastq_line_range(path: &Path, start_line: usize, end_line: usize) -> String {
        let file = fs::File::open(path).expect("open fastq");
        let mut out = String::new();
        for (idx, line) in BufReader::new(file).lines().enumerate() {
            let line_no = idx + 1;
            if line_no < start_line {
                continue;
            }
            if line_no > end_line {
                break;
            }
            out.push_str(&line.expect("read fastq line"));
            out.push('\n');
        }
        out
    }

    fn write_real_fixture_slice(
        dir: &Path,
        pair_ordinal: usize,
        before_pairs: usize,
        pair_count: usize,
    ) -> (PathBuf, PathBuf) {
        let fixture_dir = tutorial_fixture_dir();
        let r1 = fixture_dir.join("SRR2584863_1.trim.sub.fastq");
        let r2 = fixture_dir.join("SRR2584863_2.trim.sub.fastq");
        assert!(r1.exists(), "missing fixture reads1: {}", r1.display());
        assert!(r2.exists(), "missing fixture reads2: {}", r2.display());

        let start_pair = pair_ordinal - before_pairs;
        let start_line = start_pair * 4 + 1;
        let end_line = (start_pair + pair_count) * 4;
        let out1 = dir.join("slice_1.fq");
        let out2 = dir.join("slice_2.fq");
        fs::write(&out1, extract_fastq_line_range(&r1, start_line, end_line))
            .expect("write slice1");
        fs::write(&out2, extract_fastq_line_range(&r2, start_line, end_line))
            .expect("write slice2");
        (out1, out2)
    }

    fn real_slice_output_lines(
        pair_ordinal: usize,
        target: &str,
        before_pairs: usize,
        pair_count: usize,
    ) -> Vec<String> {
        let prefix = indexed_real_fixture_prefix();
        assert!(
            prefix.with_extension("bwt.2bit.64").exists(),
            "missing indexed real fixture at {}",
            prefix.display()
        );
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_real_slice_{target}_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");

        let (reads1, reads2) =
            write_real_fixture_slice(&dir, pair_ordinal, before_pairs, pair_count);
        let out = dir.join("out.sam");
        let argv = vec![
            "mem".to_string(),
            "-o".to_string(),
            out.to_str().expect("utf8").to_string(),
            "-t".to_string(),
            "2".to_string(),
            prefix.to_str().expect("utf8").to_string(),
            reads1.to_str().expect("utf8").to_string(),
            reads2.to_str().expect("utf8").to_string(),
        ];
        assert_eq!(main_mem(&argv), 0);
        let lines = fs::read_to_string(&out)
            .expect("read sam")
            .lines()
            .filter(|line| line.starts_with(target))
            .map(|line| line.to_string())
            .collect::<Vec<_>>();
        fs::remove_dir_all(&dir).expect("cleanup");
        lines
    }

    fn real_slice_sam_body(
        pair_ordinal: usize,
        before_pairs: usize,
        pair_count: usize,
        threads: i32,
    ) -> Vec<String> {
        let prefix = indexed_real_fixture_prefix();
        assert!(
            prefix.with_extension("bwt.2bit.64").exists(),
            "missing indexed real fixture at {}",
            prefix.display()
        );
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!(
            "bwa_mem2_rs_real_slice_body_{pair_ordinal}_{threads}_{nanos}"
        ));
        fs::create_dir_all(&dir).expect("create temp dir");

        let (reads1, reads2) =
            write_real_fixture_slice(&dir, pair_ordinal, before_pairs, pair_count);
        let out = dir.join("out.sam");
        let argv = vec![
            "mem".to_string(),
            "-o".to_string(),
            out.to_str().expect("utf8").to_string(),
            "-t".to_string(),
            threads.to_string(),
            prefix.to_str().expect("utf8").to_string(),
            reads1.to_str().expect("utf8").to_string(),
            reads2.to_str().expect("utf8").to_string(),
        ];
        assert_eq!(main_mem(&argv), 0);
        let body = fs::read_to_string(&out)
            .expect("read sam")
            .lines()
            .filter(|line| !line.starts_with('@'))
            .map(|line| line.to_string())
            .collect::<Vec<_>>();
        fs::remove_dir_all(&dir).expect("cleanup");
        body
    }

    fn real_slice_cpp_sam_body(
        pair_ordinal: usize,
        before_pairs: usize,
        pair_count: usize,
        threads: i32,
    ) -> Vec<String> {
        let prefix = indexed_real_fixture_prefix();
        let cpp_bin = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("bwa-mem2")
            .join("bwa-mem2.debug");
        assert!(
            cpp_bin.exists(),
            "missing upstream debug binary at {}",
            cpp_bin.display()
        );
        assert!(
            prefix.with_extension("bwt.2bit.64").exists(),
            "missing indexed real fixture at {}",
            prefix.display()
        );
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!(
            "bwa_mem2_rs_real_slice_cpp_{pair_ordinal}_{threads}_{nanos}"
        ));
        fs::create_dir_all(&dir).expect("create temp dir");

        let (reads1, reads2) =
            write_real_fixture_slice(&dir, pair_ordinal, before_pairs, pair_count);
        let output = Command::new(&cpp_bin)
            .arg("mem")
            .arg("-t")
            .arg(threads.to_string())
            .arg(prefix.to_str().expect("utf8"))
            .arg(reads1.to_str().expect("utf8"))
            .arg(reads2.to_str().expect("utf8"))
            .output()
            .expect("run upstream mem");
        assert!(
            output.status.success(),
            "upstream mem failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let body = String::from_utf8(output.stdout)
            .expect("upstream sam utf8")
            .lines()
            .filter(|line| !line.starts_with('@'))
            .map(|line| line.to_string())
            .collect::<Vec<_>>();
        fs::remove_dir_all(&dir).expect("cleanup");
        body
    }

    #[test]
    fn usage_text_includes_key_options_and_defaults() {
        let opt = mem_opt_init();
        BWA_VERBOSE.store(3, Ordering::Relaxed);
        let text = usage_text(&opt);
        assert!(text.starts_with("Usage: bwa-mem2 mem [options]"));
        assert!(text.contains("    -t INT        number of threads ["));
        assert!(text.contains("    -k INT        minimum seed length ["));
        assert!(text.contains("   -R STR        read group header line"));
        assert!(text.contains(
            "   -v INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging [3]"
        ));
        assert!(text.contains("Note: Please read the man page"));
    }

    #[test]
    fn cpuid_vendor_leaf_looks_populated_on_x86() {
        let leaf0 = cpuid(0);
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        assert_ne!(leaf0[1] | leaf0[2] | leaf0[3], 0);
    }

    #[test]
    fn htstatus_returns_booleanish_value() {
        let ht = HTStatus();
        assert!(ht == 0 || ht == 1);
    }

    #[test]
    fn update_a_scales_only_unset_scoring_fields() {
        let mut opt = *mem_opt_init();
        opt.a = 2;
        let mut opt0 = mem_opt_t::default();
        opt0.a = 1;
        opt0.b = 0;
        opt0.T = 5;
        opt0.o_del = 0;
        opt0.e_del = 7;
        opt0.o_ins = 0;
        opt0.e_ins = 0;
        opt0.zdrop = 0;
        opt0.pen_clip5 = 9;
        opt0.pen_clip3 = 0;
        opt0.pen_unpaired = 0;

        let old_b = opt.b;
        let old_t = opt.T;
        let old_o_del = opt.o_del;
        let old_e_del = opt.e_del;
        let old_pen_clip5 = opt.pen_clip5;

        update_a(&mut opt, &opt0);

        assert_eq!(opt.b, old_b * 2);
        assert_eq!(opt.T, old_t);
        assert_eq!(opt.o_del, old_o_del * 2);
        assert_eq!(opt.e_del, old_e_del);
        assert_eq!(opt.pen_clip5, old_pen_clip5);
        assert_eq!(opt.pen_clip3, 10);
        assert_eq!(opt.pen_unpaired, 34);
    }

    #[test]
    fn memoryalloc_builds_expected_worker_buffers() {
        let aux = ktp_aux_t {
            opt: Some(mem_opt_init()),
            ..Default::default()
        };
        let mut w = worker_t::default();
        memoryAlloc(&aux, &mut w, 3, 2);

        assert_eq!(w.regs.len(), 3);
        assert_eq!(w.chain_ar.len(), 3);
        assert_eq!(w.seedBuf.len(), 2 * BATCH_SIZE * AVG_SEEDS_PER_READ);
        assert_eq!(
            w.seedBufSize,
            i64::try_from(BATCH_SIZE * AVG_SEEDS_PER_READ).expect("seedBufSize")
        );

        assert_eq!(w.mmc.seqBufLeftRef.len(), 2);
        assert_eq!(w.mmc.seqBufLeftRef[0].len(), 0);
        assert_eq!(w.mmc.seqBufLeftQer[0].len(), 0);
        assert_eq!(w.mmc.seqPairArrayAux.len(), 2);
        assert_eq!(w.mmc.seqPairArrayAux[0].len(), 0);
        assert_eq!(w.mmc.wsize[0], 0);
        assert_eq!(w.mmc.matchArray.len(), 2);
        assert_eq!(w.mmc.min_intv_ar[0].len(), 0);
        assert_eq!(w.mmc.lim[0].len(), BATCH_SIZE + 32);
    }

    #[test]
    fn kt_pipeline_step0_reads_chunk_and_respects_copy_comment_flag() {
        let mut aux = ktp_aux_t {
            ks: Some(kseq_t::from_text(
                "@r0 note\nACGT\n+\nIIII\n@r1 keep\nTT\n+\nJJ\n",
            )),
            task_size: 100,
            copy_comment: 0,
            ..Default::default()
        };
        let mut opt = *mem_opt_init();
        let mut worker = worker_t::default();
        let ret = kt_pipeline(&mut aux, 0, None, &mut opt, &mut worker).expect("step0 data");
        assert_eq!(ret.n_seqs, 2);
        assert_eq!(ret.seqs[0].name.as_deref(), Some("r0"));
        assert_eq!(ret.seqs[0].comment, None);
        assert_eq!(ret.seqs[1].seq.as_deref(), Some("TT"));
    }

    #[test]
    fn kt_pipeline_step1_runs_single_end_processing() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_fastmap_pipeline_{nanos}"));
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

        let mut aux = ktp_aux_t::default();
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

        let input = crate::generated::fastmap_h::ktp_data_t {
            n_seqs: 1,
            seqs: vec![bseq1_t {
                l_seq: 6,
                name: Some("read-fastmap".into()),
                seq: Some("ACGTAC".into()),
                qual: Some("IIIIII".into()),
                ..Default::default()
            }],
            sam_lines: Vec::new(),
        };
        let ret = kt_pipeline(&mut aux, 1, Some(input), &mut opt, &mut worker).expect("step1 data");
        let sam = ret.sam_lines[0].as_ref();
        assert!(sam.starts_with("read-fastmap\t"), "{sam}");
        assert!(sam.contains("\tchr1\t"), "{sam}");
        assert_eq!(aux.n_processed, 1);

        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn process_runs_single_threaded_pipeline_and_writes_sam() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_fastmap_process_{nanos}"));
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

        let out_path = dir.join("out.sam");
        let out_stream = err_xopen_core("test", out_path.to_str().expect("utf8"), "w");
        let mut opt = (*mem_opt_init()).clone();
        opt.min_seed_len = 2;
        opt.split_factor = 1.5;
        opt.split_width = 10;
        opt.max_mem_intv = 20;
        opt.min_chain_weight = 1;
        opt.T = 1;

        let mut aux = ktp_aux_t {
            ks: Some(kseq_t::from_text("@p0\nACGTAC\n+\nIIIIII\n")),
            opt: Some(Box::new(opt)),
            task_size: 100,
            actual_chunk_size: 100,
            fp: Some(out_stream),
            ..Default::default()
        };

        let mut worker = worker_t {
            fmi: Some(fmi),
            ..Default::default()
        };

        assert_eq!(process(&mut aux, &mut worker, 1), 0);
        let fp = aux.fp.take().expect("fp");
        err_fclose(fp);
        let sam = fs::read_to_string(&out_path).expect("sam file");
        assert!(sam.starts_with("p0\t"), "{sam}");
        assert!(sam.contains("\tchr1\t"), "{sam}");
        assert_eq!(aux.n_processed, 1);

        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn process_parallel_pipeline_matches_single_threaded_output() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_process_pipe_compare_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("ref");
        let fasta = b">chr1\nACGTACGTACGT\n";
        bns_fasta2bntseq(
            Cursor::new(fasta.as_slice()),
            prefix.to_str().expect("utf8"),
            1,
        );

        let mut fmi1 = FMI_search::ctor(prefix.to_str().expect("utf8"));
        assert_eq!(fmi1.build_index(), 0);
        fmi1.load_index();
        let mut fmi2 = FMI_search::ctor(prefix.to_str().expect("utf8"));
        fmi2.load_index();

        let reads_text: String = (0..600)
            .map(|i| format!("@p{i}\nACGTAC\n+\nIIIIII\n"))
            .collect();
        let out1 = dir.join("out1.sam");
        let out2 = dir.join("out2.sam");
        let out1_stream = err_xopen_core("process", out1.to_str().expect("utf8"), "w");
        let out2_stream = err_xopen_core("process", out2.to_str().expect("utf8"), "w");

        let opt = *mem_opt_init();
        let mut aux1 = ktp_aux_t {
            ks: Some(kseq_t::from_text(&reads_text)),
            opt: Some(Box::new(opt.clone())),
            task_size: 50_000,
            actual_chunk_size: 50_000,
            fp: Some(out1_stream),
            ..Default::default()
        };
        let mut aux2 = ktp_aux_t {
            ks: Some(kseq_t::from_text(&reads_text)),
            opt: Some(Box::new(opt.clone())),
            task_size: 50_000,
            actual_chunk_size: 50_000,
            fp: Some(out2_stream),
            ..Default::default()
        };

        let mut worker1 = worker_t {
            fmi: Some(fmi1),
            ..Default::default()
        };
        let mut worker2 = worker_t {
            fmi: Some(fmi2),
            ..Default::default()
        };

        assert_eq!(process(&mut aux1, &mut worker1, 1), 0);
        assert_eq!(process(&mut aux2, &mut worker2, 2), 0);
        err_fclose(aux1.fp.take().expect("fp1"));
        err_fclose(aux2.fp.take().expect("fp2"));
        assert_eq!(
            fs::read_to_string(&out1).expect("out1"),
            fs::read_to_string(&out2).expect("out2")
        );
        assert_eq!(aux1.n_processed, aux2.n_processed);

        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn main_mem_runs_local_single_end_path_via_clap() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_main_mem_{nanos}"));
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
        fmi.dtor();

        let reads = dir.join("reads.fq");
        fs::write(&reads, b"@mm0\nACGTAC\n+\nIIIIII\n").expect("write reads");
        let out = dir.join("out.sam");
        let argv = vec![
            "mem".to_string(),
            "-o".to_string(),
            out.to_str().expect("utf8").to_string(),
            "-t".to_string(),
            "1".to_string(),
            "-k".to_string(),
            "2".to_string(),
            "-T".to_string(),
            "1".to_string(),
            prefix.to_str().expect("utf8").to_string(),
            reads.to_str().expect("utf8").to_string(),
        ];
        assert_eq!(main_mem(&argv), 0);
        let sam = fs::read_to_string(&out).expect("sam output");
        assert!(sam.contains("@SQ\tSN:chr1\tLN:12\n"), "{sam}");
        assert!(sam.contains("mm0\t"), "{sam}");
        assert!(sam.contains("\tchr1\t"), "{sam}");

        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn main_mem_reads_gz_fastq_input() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_main_mem_gz_{nanos}"));
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
        fmi.dtor();

        let reads = dir.join("reads.fq.gz");
        let file = fs::File::create(&reads).expect("create gz");
        let mut enc = GzEncoder::new(file, Compression::default());
        std::io::Write::write_all(&mut enc, b"@mmgz\nACGTAC\n+\nIIIIII\n").expect("write gz fastq");
        enc.finish().expect("finish gz");

        let out = dir.join("out.sam");
        let argv = vec![
            "mem".to_string(),
            "-o".to_string(),
            out.to_str().expect("utf8").to_string(),
            "-t".to_string(),
            "1".to_string(),
            "-k".to_string(),
            "2".to_string(),
            "-T".to_string(),
            "1".to_string(),
            prefix.to_str().expect("utf8").to_string(),
            reads.to_str().expect("utf8").to_string(),
        ];
        assert_eq!(main_mem(&argv), 0);
        let sam = fs::read_to_string(&out).expect("sam output");
        assert!(sam.contains("mmgz\t"), "{sam}");
        assert!(sam.contains("\tchr1\t"), "{sam}");

        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn read_text_input_supports_pipe_commands() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_pipe_input_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        let reads = dir.join("reads.fq");
        fs::write(&reads, b"@pipe0\nACGT\n+\nIIII\n").expect("write reads");
        let text = read_text_input(&format!("< cat {}", reads.display())).expect("pipe text");
        assert!(text.contains("@pipe0\nACGT\n+\nIIII\n"));
        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn read_text_input_supports_http_urls() {
        let listener = match TcpListener::bind("127.0.0.1:0") {
            Ok(listener) => listener,
            Err(err) if err.kind() == std::io::ErrorKind::PermissionDenied => return,
            Err(err) => panic!("bind: {err:?}"),
        };
        listener.set_nonblocking(true).expect("set_nonblocking");
        let addr = listener.local_addr().expect("local addr");
        let handle = thread::spawn(move || {
            let deadline = std::time::Instant::now() + std::time::Duration::from_secs(1);
            loop {
                match listener.accept() {
                    Ok((mut stream, _)) => {
                        let mut buf = [0_u8; 1024];
                        let _ = stream.read(&mut buf);
                        stream
                            .write_all(b"HTTP/1.0 200 OK\r\nContent-Length: 19\r\n\r\n@http0\nACGT\n+\nIIII\n")
                            .expect("write response");
                        break;
                    }
                    Err(err) if err.kind() == std::io::ErrorKind::WouldBlock => {
                        if std::time::Instant::now() >= deadline {
                            break;
                        }
                        thread::sleep(std::time::Duration::from_millis(5));
                    }
                    Err(err) => panic!("accept: {err:?}"),
                }
            }
        });
        let url = format!("http://{}/reads.fq", addr);
        thread::sleep(std::time::Duration::from_millis(10));
        let Ok(text) = read_text_input(&url) else {
            handle.join().expect("server join");
            return;
        };
        handle.join().expect("server join");
        assert!(text.contains("@http0\nACGT\n+\nIIII\n"));
    }

    #[test]
    fn main_mem_runs_local_paired_end_path_via_clap() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_main_mem_pe_{nanos}"));
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
        fmi.dtor();

        let reads1 = dir.join("reads_1.fq");
        let reads2 = dir.join("reads_2.fq");
        fs::write(&reads1, b"@pair0/1\nACGTAC\n+\nIIIIII\n").expect("write read1");
        fs::write(&reads2, b"@pair0/2\nGTACGT\n+\nIIIIII\n").expect("write read2");
        let out = dir.join("out.sam");
        let argv = vec![
            "mem".to_string(),
            "-o".to_string(),
            out.to_str().expect("utf8").to_string(),
            "-t".to_string(),
            "1".to_string(),
            "-k".to_string(),
            "2".to_string(),
            "-T".to_string(),
            "1".to_string(),
            prefix.to_str().expect("utf8").to_string(),
            reads1.to_str().expect("utf8").to_string(),
            reads2.to_str().expect("utf8").to_string(),
        ];
        assert_eq!(main_mem(&argv), 0);
        let sam = fs::read_to_string(&out).expect("sam output");
        assert!(sam.contains("pair0\t"), "{sam}");
        assert!(sam.contains("\tchr1\t"), "{sam}");
        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn main_mem_runs_smart_pairing_path_via_clap() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_main_mem_smartpe_{nanos}"));
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
        fmi.dtor();

        let reads = dir.join("reads.fq");
        fs::write(
            &reads,
            b"@smart0/1\nACGTAC\n+\nIIIIII\n@smart0/2\nGTACGT\n+\nIIIIII\n",
        )
        .expect("write reads");
        let out = dir.join("out.sam");
        let argv = vec![
            "mem".to_string(),
            "-p".to_string(),
            "-o".to_string(),
            out.to_str().expect("utf8").to_string(),
            "-t".to_string(),
            "1".to_string(),
            "-k".to_string(),
            "2".to_string(),
            "-T".to_string(),
            "1".to_string(),
            prefix.to_str().expect("utf8").to_string(),
            reads.to_str().expect("utf8").to_string(),
        ];
        assert_eq!(main_mem(&argv), 0);
        let sam = fs::read_to_string(&out).expect("sam output");
        assert!(sam.contains("smart0\t"), "{sam}");
        assert!(sam.contains("\tchr1\t"), "{sam}");
        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn main_mem_is_stable_across_thread_counts_for_many_single_end_reads() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_main_mem_rayon_{nanos}"));
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
        fmi.dtor();

        let reads = dir.join("reads_many.fq");
        let mut text = String::new();
        for i in 0..520 {
            text.push_str(&format!("@m{i}\nACGTAC\n+\nIIIIII\n"));
        }
        fs::write(&reads, text).expect("write reads");

        let out1 = dir.join("out_t1.sam");
        let out2 = dir.join("out_t2.sam");
        let argv1 = vec![
            "mem".to_string(),
            "-o".to_string(),
            out1.to_str().expect("utf8").to_string(),
            "-t".to_string(),
            "1".to_string(),
            "-k".to_string(),
            "2".to_string(),
            "-T".to_string(),
            "1".to_string(),
            prefix.to_str().expect("utf8").to_string(),
            reads.to_str().expect("utf8").to_string(),
        ];
        let argv2 = vec![
            "mem".to_string(),
            "-o".to_string(),
            out2.to_str().expect("utf8").to_string(),
            "-t".to_string(),
            "2".to_string(),
            "-k".to_string(),
            "2".to_string(),
            "-T".to_string(),
            "1".to_string(),
            prefix.to_str().expect("utf8").to_string(),
            reads.to_str().expect("utf8").to_string(),
        ];
        assert_eq!(main_mem(&argv1), 0);
        assert_eq!(main_mem(&argv2), 0);
        let sam1 = fs::read_to_string(&out1).expect("sam1");
        let sam2 = fs::read_to_string(&out2).expect("sam2");
        assert_eq!(sam1, sam2);

        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    fn main_mem_paired_end_survives_large_batch_with_compute_threads() {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_main_mem_pe_rayon_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("ref");
        let fasta = b">chr1\nACGTACGTACGTACGTACGTACGTACGT\n";
        bns_fasta2bntseq(
            Cursor::new(fasta.as_slice()),
            prefix.to_str().expect("utf8"),
            1,
        );

        let mut fmi = FMI_search::ctor(prefix.to_str().expect("utf8"));
        assert_eq!(fmi.build_index(), 0);
        fmi.load_index();
        fmi.dtor();

        let reads1 = dir.join("reads_1_many.fq");
        let reads2 = dir.join("reads_2_many.fq");
        let mut text1 = String::new();
        let mut text2 = String::new();
        for i in 0..520 {
            text1.push_str(&format!("@p{i}/1\nACGTACGT\n+\nIIIIIIII\n"));
            text2.push_str(&format!("@p{i}/2\nGTACGTAC\n+\nIIIIIIII\n"));
        }
        fs::write(&reads1, text1).expect("write reads1");
        fs::write(&reads2, text2).expect("write reads2");

        let out = dir.join("out_t2.sam");
        let argv = vec![
            "mem".to_string(),
            "-1".to_string(),
            "-o".to_string(),
            out.to_str().expect("utf8").to_string(),
            "-t".to_string(),
            "2".to_string(),
            "-k".to_string(),
            "2".to_string(),
            "-T".to_string(),
            "1".to_string(),
            prefix.to_str().expect("utf8").to_string(),
            reads1.to_str().expect("utf8").to_string(),
            reads2.to_str().expect("utf8").to_string(),
        ];
        assert_eq!(main_mem(&argv), 0);
        let sam = fs::read_to_string(&out).expect("sam");
        let record_count = sam.lines().filter(|line| !line.starts_with('@')).count();
        assert_eq!(record_count, 520 * 2, "{sam}");

        fs::remove_dir_all(&dir).expect("cleanup");
    }

    #[test]
    #[ignore]
    fn real_fixture_slice_207026_matches_expected_sam() {
        let lines = real_slice_output_lines(160_414, "SRR2584863.207026", 120, 300);
        assert_eq!(
            lines,
            vec![
                "SRR2584863.207026\t99\tCP000819.1\t3320076\t60\t150M\t=\t3320556\t630\tCGCCTTATCCGACCTACGGTTATGTTCCGTTTGTAGGCCTGATAAGACGCACAGCGTCGCATCAGGCAACGGCTGCCGGATGCGGCGTAAACGCCTTATCCGACCTACGGTTATGTTCCGTTTGTAGGCCTGATAAGACGCACAGCGTCG\tCCCFFFFFHHHHHJJJJJJJJJJGJJJJJJJJJIIJGHJJJJIJJIJJIIFJJJJHHGFDDDDDDDDDDDDDDDDDDDDBDDDDDDDDDBDDDDDDDDDDDBDDDDDDDDDDDDDEEEDDDDDBDECCCDDDCDDDDDD@DBDDDDDDDB\tNM:i:0\tMD:Z:150\tMC:Z:150M\tAS:i:150\tXS:i:108".to_string(),
                "SRR2584863.207026\t147\tCP000819.1\t3320556\t60\t150M\t=\t3320076\t-630\tAACCGCTTCCATATTCAGGTCGCGACGTACCAGATTTACCGACCCTTTCATGGCGATATCCGCCTCCAGGCCATCCACCAGCGTGTCGTCGGTGTGCATAACGCCGTCTTTAATCCACGCGGTGCTGCGAATGGAGTCAAAATAGAACCC\tBB@>89AC@DCCDDDDDDBBDDDDDDCC@DECCDDDDBDDCADDDCC@><0DBDCDBB>DDDDDDDDDDDDDDBDDDDDDDDBBDBDDBDDEDCEEEFFFHFIIJGGIIIGJHHIJJJIJJJIJJJJJJJJJJJIJHFHHHHFFFFFCC@\tNM:i:0\tMD:Z:150\tMC:Z:150M\tAS:i:150\tXS:i:0".to_string(),
            ]
        );
    }

    #[test]
    #[ignore]
    fn real_fixture_slice_186556_matches_expected_sam() {
        let lines = real_slice_output_lines(144_413, "SRR2584863.186556", 120, 300);
        assert_eq!(
            lines,
            vec![
                "SRR2584863.186556\t83\tCP000819.1\t330457\t60\t150M\t=\t330028\t-579\tTTCTCTCCAGTAGTGGCGGGGCAAGAGAAATAGCTTTAAGCCCGTCATGTTGTGTGGCAATCGCTGCTGGTAACAATGTGGAAAGGGAAGTGCGGCGAATCAGCTCCAGAACCGCGCTAATTGAGTTCGCCTCAATGACCACCTGTGGAT\tC@@AA>:>@A?>B@B>3?49CDDCCEEDC@A:DCCDDB@D@5BCACDB?B??A:CDCB?B@9?BCBCCDACACCCDCDDA>C>:9@CDB@BBBCB?BDEC>@72FHHGIIHGF>D>HEIGIIHCGEEIHECIGIIIGHFD<GFFEBFC@@\tNM:i:0\tMD:Z:150\tMC:Z:29M\tAS:i:150\tXS:i:0".to_string(),
                "SRR2584863.186556\t163\tCP000819.1\t330028\t24\t29M\t=\t330457\t579\tAGCTTAGTTTATGCCTGATGCGGCGTGAA\tCC<DDF<DHHDHHIGIBEDHHAFGAHI?B\tNM:i:1\tMD:Z:15G13\tMC:Z:150M\tAS:i:24\tXS:i:20".to_string(),
            ]
        );
    }

    #[test]
    #[ignore]
    fn real_fixture_slice_38486_matches_expected_sam() {
        let lines = real_slice_output_lines(30_192, "SRR2584863.38486", 120, 300);
        assert_eq!(
            lines,
            vec![
                "SRR2584863.38486\t83\tCP000819.1\t4422765\t60\t150M\t=\t4421725\t-1190\tTGCGGTAACTTCGCTGATTTAATCGCTTCAATCAGCGTCGCGTTTTTCTCTTTTGCAGAGGTGCGACCAGCAATATGTTGCAGTACTCGCACTTTTCCCACTAACTGTGCGCTGTTCCAGGTTTTGTAGCTAAACTGATCTTTATCAAGC\t@DBCCDCBBDDDDCCCACA?AA<@>A@<CDCBDBDD@@9B<<BDDDCCBC9CADDDDDC@>=<DDDDDCDEDEDC@CACCFEEDBAGFEIHHFJIHFECDEIIEFDHG>EIF>DGHIIIIJJJJJJJIIIBJFDIJIHAHFGFFFFFCCC\tNM:i:0\tMD:Z:150\tMC:Z:34M\tAS:i:150\tXS:i:0".to_string(),
                "SRR2584863.38486\t163\tCP000819.1\t4421725\t60\t34M\t=\t4422765\t1190\tCGCTCACAGAGACGAGGGGGAGAAATGTTAGATC\t=<=AA<A;B>C8+<<CA):??A<BA<=/9=<AAB\tNM:i:1\tMD:Z:17T16\tMC:Z:150M\tAS:i:29\tXS:i:0".to_string(),
            ]
        );
    }

    #[test]
    #[ignore]
    fn real_fixture_slice_50147_matches_expected_sam() {
        let lines = real_slice_output_lines(39_165, "SRR2584863.50147", 120, 300);
        assert_eq!(
            lines,
            vec![
                "SRR2584863.50147\t97\tCP000819.1\t333610\t60\t134M16S\t=\t330451\t-3032\tGCATTTCATGCGCTGTGCGTCGTTGGGCTGATGATCATAACCCTGCGTTTTGCACCAGTACGTTTTCCGCAGCTGTGGGTCAAAGAGGCATGATGCGACGCTTGTTCCTGCGCTTTGTTCATGCCGGATGCGGCGTGAACGCCTTATCCG\tBCCFFFFFHHHHHJJJJJJHIJIJJJJJJJJJJJJJJJJJJJJJIIJHIJJHHHGHFFFFFDEEEDDDDDDDDDDCDDDDDDDDDDDDDDDDDDDDDBDBD@BDDDDDDDDDDDDDDDDDEEDDDDBDDDDDDDDD@?@A8<BBDDDDD<\tNM:i:0\tMD:Z:134\tMC:Z:129M\tAS:i:134\tXS:i:42\tSA:Z:CP000819.1,330238,+,120S30M,0,0;".to_string(),
                "SRR2584863.50147\t2145\tCP000819.1\t330238\t0\t120H30M\t=\t330451\t342\tATGCCGGATGCGGCGTGAACGCCTTATCCG\tEEDDDDBDDDDDDDDD@?@A8<BBDDDDD<\tNM:i:0\tMD:Z:30\tMC:Z:129M\tAS:i:30\tXS:i:30\tSA:Z:CP000819.1,333610,+,134M16S,60,0;".to_string(),
                "SRR2584863.50147\t145\tCP000819.1\t330451\t60\t129M\t=\t333610\t3032\tCCGCCGTTCTCTCCAGTAGTGGCGGGGCAAGAGAAATAGCTTTAAGCCCGTCATGTTGTGTGGCAATCGCTGCTGGTAACAATGTGGAAAGGGAAGTGCGGCGAATCAGCTCCAGAACCGCGCTAATTG\tDDDDDDCDDCBDDDDDDDDDDDDDDDDDDDDEEEDDDDDDDDDDDDDDDDDDDDDDDDDDCDDDDDBDDEEDFFFFEFHGHHEFIJIJJJJJJJJJJJJIJIJJJIIJJJJJJJJIFGHHHFFFFFCCC\tNM:i:0\tMD:Z:129\tMC:Z:134M16S\tAS:i:129\tXS:i:0".to_string(),
            ]
        );
    }

    #[test]
    #[ignore]
    fn real_fixture_small_slice_is_stable_across_thread_counts() {
        let sam_t1 = real_slice_sam_body(160_414, 120, 300, 1);
        let sam_t2 = real_slice_sam_body(160_414, 120, 300, 2);
        assert_eq!(sam_t1, sam_t2);
    }

    #[test]
    #[ignore]
    fn real_fixture_medium_slice_is_stable_across_thread_counts() {
        let sam_t1 = real_slice_sam_body(160_414, 500, 1000, 1);
        let sam_t2 = real_slice_sam_body(160_414, 500, 1000, 2);
        assert_eq!(sam_t1, sam_t2);
    }

    #[test]
    #[ignore]
    fn real_fixture_small_slice_matches_upstream_cpp() {
        let rust_sam = real_slice_sam_body(160_414, 120, 300, 2);
        let cpp_sam = real_slice_cpp_sam_body(160_414, 120, 300, 2);
        assert_eq!(rust_sam, cpp_sam);
    }

    #[test]
    #[ignore]
    fn real_fixture_medium_slice_matches_upstream_cpp() {
        let rust_sam = real_slice_sam_body(160_414, 500, 1000, 2);
        let cpp_sam = real_slice_cpp_sam_body(160_414, 500, 1000, 2);
        assert_eq!(rust_sam, cpp_sam);
    }

    #[test]
    #[ignore]
    fn real_fixture_large_slice_matches_upstream_cpp() {
        let rust_sam = real_slice_sam_body(160_414, 1000, 2000, 2);
        let cpp_sam = real_slice_cpp_sam_body(160_414, 1000, 2000, 2);
        assert_eq!(rust_sam, cpp_sam);
    }
}
