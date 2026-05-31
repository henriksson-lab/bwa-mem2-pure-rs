use bwa_mem2_pure_rs::bwa_mem2::bwa::bseq1_t;
use bwa_mem2_pure_rs::bwa_mem2::bwa::{bseq_read_orig, bwa_fill_scmat};
use bwa_mem2_pure_rs::bwa_mem2::bwamem::worker_t;
use bwa_mem2_pure_rs::bwa_mem2::bwamem::{mem_alnreg_v, mem_chain_v, mem_seed_t};
use bwa_mem2_pure_rs::bwa_mem2::bwamem::{
    mem_chain2aln_across_reads_v2, mem_mark_primary_se, mem_opt_init, sort_classify, worker_bwt,
};
use bwa_mem2_pure_rs::bwa_mem2::bwamem_pair::{
    mem_pair, mem_pestat, mem_sam_pe_batch, mem_sam_pe_batch_post, mem_sam_pe_batch_pre,
};
use bwa_mem2_pure_rs::bwa_mem2::fastmap::ktp_aux_t;
use bwa_mem2_pure_rs::bwa_mem2::fastmap::{memory_alloc, update_a};
use bwa_mem2_pure_rs::bwa_mem2::fmi_search::FMI_search;
use bwa_mem2_pure_rs::bwa_mem2::kseq::kseq_t;
use bwa_mem2_pure_rs::bwa_mem2::ksw::kswr_t;

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

fn main() {
    let mut args = std::env::args().skip(1);
    let prefix = args.next().expect("prefix");
    let r1 = args.next().expect("r1");
    let r2 = args.next().expect("r2");
    let target = args.next().expect("target read name");
    let n_threads: i32 = args
        .next()
        .as_deref()
        .unwrap_or("1")
        .parse()
        .expect("n_threads");
    let extra_args: Vec<String> = args.collect();

    let r1_text = std::fs::read_to_string(&r1).expect("read r1");
    let r2_text = std::fs::read_to_string(&r2).expect("read r2");
    let mut ks1 = kseq_t::from_text(&r1_text);
    let mut ks2 = kseq_t::from_text(&r2_text);

    let mut fmi = FMI_search::ctor(&prefix);
    fmi.load_index();

    let mut opt = (*mem_opt_init()).clone();
    let mut opt0 = opt.clone();
    opt.flag |= 0x2;
    opt.n_threads = n_threads.max(1);
    let mut i = 0;
    while i < extra_args.len() {
        match extra_args[i].as_str() {
            "-A" => {
                i += 1;
                opt.a = extra_args[i].parse().expect("-A value");
                opt0.a = 1;
            }
            "-d" => {
                i += 1;
                opt.zdrop = extra_args[i].parse().expect("-d value");
                opt0.zdrop = 1;
            }
            "-k" => {
                i += 1;
                opt.min_seed_len = extra_args[i].parse().expect("-k value");
                opt0.min_seed_len = 1;
            }
            other => panic!("unsupported dump flag {other}"),
        }
        i += 1;
    }
    update_a(&mut opt, &opt0);
    bwa_fill_scmat(opt.a, opt.b, &mut opt.mat);

    let task_size = opt.chunk_size * i64::from(opt.n_threads);
    let aux = ktp_aux_t {
        opt: Some(Box::new(opt.clone())),
        task_size,
        actual_chunk_size: task_size,
        ..Default::default()
    };

    let mut n = 0_i32;
    let mut total_bases = 0_i64;
    let mut batch_index = 0_i64;
    let mut n_processed = 0_i64;
    let seqs = loop {
        let seqs = bseq_read_orig(
            task_size,
            &mut n,
            &mut ks1,
            Some(&mut ks2),
            &mut total_bases,
        );
        assert!(n > 0, "target not found");
        if seqs
            .iter()
            .any(|s| s.name.as_deref() == Some(target.as_str()))
        {
            break seqs;
        }
        n_processed += i64::from(n);
        batch_index += 1;
    };

    let mut worker = worker_t {
        fmi: Some(fmi),
        ..Default::default()
    };
    memory_alloc(&aux, &mut worker, n, 1);
    worker.opt = Some(Box::new(opt.clone()));
    worker.seqs = seqs;
    worker.nreads = n;
    worker.n_processed = n_processed;
    worker.nthreads = i16::try_from(opt.n_threads.max(1)).expect("nthreads");
    let n_usize = usize::try_from(n).expect("n");
    worker.regs.resize(n_usize, mem_alnreg_v::default());
    worker.chain_ar.resize(n_usize, mem_chain_v::default());
    let want_seed_buf = usize::try_from(worker.nthreads.max(1))
        .expect("nthreads")
        .saturating_mul(bwa_mem2_pure_rs::bwa_mem2::r#macro::BATCH_SIZE)
        .saturating_mul(64);
    worker.seedBuf.resize(want_seed_buf, mem_seed_t::default());
    {
        let fmi_ref = worker.fmi.as_ref().expect("fmi");
        let bns = fmi_ref.base.idx.bns.as_ref().expect("bns");
        worker.ref_string = pac_to_reference_layout(bns.l_pac, &fmi_ref.base.idx.pac);
    }

    worker_bwt(&mut worker, 0, n, 0);
    let opt_ref = worker.opt.as_deref().expect("opt");
    let fmi_ref = worker.fmi.as_ref().expect("fmi");
    let bns_ref = fmi_ref.base.idx.bns.as_ref().expect("bns");
    mem_chain2aln_across_reads_v2(
        opt_ref,
        bns_ref,
        &fmi_ref.base.idx.pac,
        &worker.seqs,
        n,
        &mut worker.chain_ar,
        &mut worker.regs,
        &mut worker.mmc,
        &worker.ref_string,
        0,
    );

    mem_pestat(
        worker.opt.as_deref().expect("opt"),
        worker
            .fmi
            .as_ref()
            .expect("fmi")
            .base
            .idx
            .bns
            .as_ref()
            .expect("bns")
            .l_pac,
        n,
        &worker.regs,
        &mut worker.pes,
    );

    let target_idx = worker
        .seqs
        .iter()
        .position(|s| s.name.as_deref() == Some(target.as_str()))
        .expect("target not found");
    let pair_start = target_idx & !1;
    println!(
        "batch_index={batch_index} n_processed={n_processed} target_idx={target_idx} pair_start={pair_start} n={n} total_bases={total_bases}"
    );
    for side in 0..2usize {
        let seq_idx = pair_start + side;
        let seq = &worker.seqs[seq_idx];
        let regs = &worker.regs[seq_idx];
        println!(
            "side={side} name={:?} l_seq={} regs_n={}",
            seq.name, seq.l_seq, regs.n
        );
        for (i, reg) in regs.a.iter().take(regs.n.min(12)).enumerate() {
            println!(
                "  reg[{i}] c={} score={} truesc={} sub={} csub={} qb={} qe={} rb={} re={} rid={} sec={} sec_all={} seedcov={} frac_rep={:.6}",
                reg.c,
                reg.score,
                reg.truesc,
                reg.sub,
                reg.csub,
                reg.qb,
                reg.qe,
                reg.rb,
                reg.re,
                reg.rid,
                reg.secondary,
                reg.secondary_all,
                reg.seedcov,
                reg.frac_rep
            );
            if reg.c < worker.chain_ar[seq_idx].a.len() {
                let chain = &worker.chain_ar[seq_idx].a[reg.c];
                println!(
                    "    chain n={} rid={} w={} pos={} frac_rep={:.6}",
                    chain.n, chain.rid, chain.w, chain.pos, chain.frac_rep
                );
                for (k, seed) in chain
                    .seeds
                    .iter()
                    .take(usize::try_from(chain.n.min(8)).expect("chain.n"))
                    .enumerate()
                {
                    println!(
                        "      seed[{k}] qbeg={} len={} rbeg={} score={} aln={}",
                        seed.qbeg, seed.len, seed.rbeg, seed.score, seed.aln
                    );
                }
            }
        }
    }

    let mut pair_regs = [
        mem_alnreg_v {
            n: worker.regs[pair_start].n,
            m: worker.regs[pair_start].m,
            a: worker.regs[pair_start].a.clone(),
        },
        mem_alnreg_v {
            n: worker.regs[pair_start + 1].n,
            m: worker.regs[pair_start + 1].m,
            a: worker.regs[pair_start + 1].a.clone(),
        },
    ];
    let _pair_seqs = [
        worker.seqs[pair_start].clone(),
        worker.seqs[pair_start + 1].clone(),
    ];
    let mut n_pri = [0_i32; 2];
    n_pri[0] = mem_mark_primary_se(
        worker.opt.as_deref().expect("opt"),
        i32::try_from(pair_regs[0].n).expect("n"),
        &mut pair_regs[0].a,
        i64::try_from((pair_start as u64) << 1).expect("id"),
    );
    n_pri[1] = mem_mark_primary_se(
        worker.opt.as_deref().expect("opt"),
        i32::try_from(pair_regs[1].n).expect("n"),
        &mut pair_regs[1].a,
        i64::try_from(((pair_start as u64) << 1) | 1).expect("id"),
    );
    println!("after_primary n_pri={n_pri:?}");
    for side in 0..2usize {
        println!("post_primary side={side}");
        for (i, reg) in pair_regs[side]
            .a
            .iter()
            .take(pair_regs[side].n.min(12))
            .enumerate()
        {
            println!(
                "  reg[{i}] score={} sub={} csub={} qb={} qe={} rb={} re={} sec={} sec_all={} seedcov={}",
                reg.score, reg.sub, reg.csub, reg.qb, reg.qe, reg.rb, reg.re, reg.secondary, reg.secondary_all, reg.seedcov
            );
        }
    }
    let mut sub = 0_i32;
    let mut n_sub = 0_i32;
    let mut z = [0_i32; 2];
    let pair_score = mem_pair(
        worker.opt.as_deref().expect("opt"),
        worker
            .fmi
            .as_ref()
            .expect("fmi")
            .base
            .idx
            .bns
            .as_ref()
            .expect("bns"),
        &[],
        &worker.pes,
        &mut pair_regs,
        i32::try_from(pair_start >> 1).expect("pair id"),
        &mut sub,
        &mut n_sub,
        &mut z,
        &n_pri,
    );
    println!("mem_pair score={pair_score} sub={sub} n_sub={n_sub} z={z:?}");

    for (i, pe) in worker.pes.iter().enumerate() {
        println!(
            "pes[{i}] low={} high={} failed={} avg={:.6} std={:.6}",
            pe.low, pe.high, pe.failed, pe.avg, pe.std
        );
    }

    let mut sam_regs = worker.regs.clone();
    let mut pcnt = 0_i32;
    let mut gcnt = 0_i32;
    let mut max_ref_len = 0_i32;
    let mut max_qer_len = 0_i32;
    for i in (0..usize::try_from(n).expect("n")).step_by(2) {
        let pair_id = u64::try_from(i >> 1).expect("pair_id");
        let seq_pair: &[bseq1_t; 2] = (&worker.seqs[i..i + 2]).try_into().expect("seq_pair");
        let reg_pair: &[mem_alnreg_v; 2] = (&sam_regs[i..i + 2]).try_into().expect("reg_pair");
        mem_sam_pe_batch_pre(
            worker.opt.as_deref().expect("opt"),
            worker
                .fmi
                .as_ref()
                .expect("fmi")
                .base
                .idx
                .bns
                .as_ref()
                .expect("bns"),
            &worker.fmi.as_ref().expect("fmi").base.idx.pac,
            &worker.pes,
            pair_id,
            seq_pair,
            reg_pair,
            &mut worker.mmc,
            &mut pcnt,
            &mut gcnt,
            &mut max_ref_len,
            &mut max_qer_len,
            0,
        );
    }
    println!(
        "batch_pre pcnt={pcnt} gcnt={gcnt} max_ref_len={max_ref_len} max_qer_len={max_qer_len}"
    );
    let pcnt8 = i32::try_from(sort_classify(&mut worker.mmc, i64::from(pcnt), 0)).expect("pcnt8");
    let mut aln = vec![kswr_t::default(); usize::try_from(pcnt + 256).expect("aln len")];
    mem_sam_pe_batch(
        worker.opt.as_deref().expect("opt"),
        &mut worker.mmc,
        pcnt,
        pcnt8,
        &mut aln,
        max_ref_len,
        max_qer_len,
        0,
    );
    gcnt = 0;
    for i in (0..usize::try_from(n).expect("n")).step_by(2) {
        let pair_id = u64::try_from(i >> 1).expect("pair_id");
        let seq_pair: &mut [bseq1_t; 2] =
            (&mut worker.seqs[i..i + 2]).try_into().expect("seq_pair");
        let reg_pair: &mut [mem_alnreg_v; 2] =
            (&mut sam_regs[i..i + 2]).try_into().expect("reg_pair");
        if i == pair_start {
            println!("before_batch_post target pair gcnt={gcnt}");
            let gar_ids: Vec<i32> = worker.mmc.seqPairArrayAux[0]
                .iter()
                .map(|sp| sp.id)
                .skip(usize::try_from(gcnt).expect("gcnt"))
                .take(64)
                .collect();
            println!("  gar_ids={gar_ids:?}");
            for (off, idx) in gar_ids.iter().enumerate() {
                if *idx >= 0 {
                    let a = aln[usize::try_from(*idx).expect("idx")];
                    println!(
                        "  aln_slot off={} idx={} score={} score2={} qb={} qe={} tb={} te={}",
                        off, idx, a.score, a.score2, a.qb, a.qe, a.tb, a.te
                    );
                }
            }
            for side in 0..2usize {
                println!("  batch_pre side={side} regs_n={}", reg_pair[side].n);
                for (j, reg) in reg_pair[side]
                    .a
                    .iter()
                    .take(reg_pair[side].n.min(12))
                    .enumerate()
                {
                    println!(
                        "    reg[{j}] score={} sub={} csub={} qb={} qe={} rb={} re={} sec={} sec_all={} seedcov={}",
                        reg.score, reg.sub, reg.csub, reg.qb, reg.qe, reg.rb, reg.re, reg.secondary, reg.secondary_all, reg.seedcov
                    );
                }
            }
        }
        mem_sam_pe_batch_post(
            worker.opt.as_deref().expect("opt"),
            worker
                .fmi
                .as_ref()
                .expect("fmi")
                .base
                .idx
                .bns
                .as_ref()
                .expect("bns"),
            &worker.fmi.as_ref().expect("fmi").base.idx.pac,
            &worker.pes,
            pair_id,
            seq_pair,
            reg_pair,
            &aln,
            &mut worker.mmc,
            &mut gcnt,
            0,
        );
        if i == pair_start {
            println!("after_batch_post target pair gcnt={gcnt}");
            for side in 0..2usize {
                println!("  batch_post side={side} regs_n={}", reg_pair[side].n);
                for (j, reg) in reg_pair[side]
                    .a
                    .iter()
                    .take(reg_pair[side].n.min(16))
                    .enumerate()
                {
                    println!(
                        "    reg[{j}] score={} sub={} csub={} qb={} qe={} rb={} re={} sec={} sec_all={} seedcov={}",
                        reg.score, reg.sub, reg.csub, reg.qb, reg.qe, reg.rb, reg.re, reg.secondary, reg.secondary_all, reg.seedcov
                    );
                }
            }
            println!(
                "target SAM1={}",
                worker.seqs[i].sam.as_deref().unwrap_or("")
            );
            println!(
                "target SAM2={}",
                worker.seqs[i + 1].sam.as_deref().unwrap_or("")
            );
            break;
        }
    }
}
