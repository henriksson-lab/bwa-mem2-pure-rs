use bwa_mem2_rs::generated::bwa_cpp::bseq_read_orig;
use bwa_mem2_rs::generated::bwamem_cpp::{
    mem_chain_flt, mem_chain_seeds, mem_collect_smem, mem_flt_chained_seeds, mem_mark_primary_se,
    mem_opt_init, worker_aln, worker_bwt,
};
use bwa_mem2_rs::generated::bwamem_h::mem_chain_v;
use bwa_mem2_rs::generated::bwamem_h::worker_t;
use bwa_mem2_rs::generated::fastmap_cpp::memoryAlloc;
use bwa_mem2_rs::generated::fastmap_h::ktp_aux_t;
use bwa_mem2_rs::generated::fmi_search_cpp::FMI_search;
use bwa_mem2_rs::generated::kseq_h::kseq_t;

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
    let fq = args.next().expect("fq");

    let fq_text = std::fs::read_to_string(&fq).expect("read fq");
    let mut ks = kseq_t::from_text(&fq_text);
    let mut n = 0_i32;
    let mut total_bases = 0_i64;
    let seqs = bseq_read_orig(i64::MAX, &mut n, &mut ks, None, &mut total_bases);
    assert_eq!(n, 1);

    let mut fmi = FMI_search::ctor(&prefix);
    fmi.load_index();

    let opt = (*mem_opt_init()).clone();
    let aux = ktp_aux_t {
        opt: Some(Box::new(opt.clone())),
        ..Default::default()
    };
    let mut worker = worker_t {
        fmi: Some(fmi),
        ..Default::default()
    };
    memoryAlloc(&aux, &mut worker, 1, 1);
    worker.opt = Some(Box::new(opt));
    worker.n_processed = 0;
    worker.nreads = 1;
    worker.seqs = seqs;
    {
        let fmi_ref = worker.fmi.as_ref().expect("worker fmi");
        let bns = fmi_ref.base.idx.bns.as_ref().expect("bns");
        worker.ref_string = pac_to_reference_layout(bns.l_pac, &fmi_ref.base.idx.pac);
    }

    let mut tot_smem = 0_i64;
    let mut match_array = Vec::new();
    let mut min_intv_ar = Vec::new();
    let mut query_pos_ar = Vec::new();
    let mut enc_qdb = Vec::new();
    let mut rid = Vec::new();
    mem_collect_smem(
        worker.fmi.as_ref().expect("worker fmi"),
        worker.opt.as_deref().expect("opt"),
        &worker.seqs,
        1,
        &mut match_array,
        &mut min_intv_ar,
        &mut query_pos_ar,
        &mut enc_qdb,
        &mut rid,
        &mut worker.mmc,
        &mut tot_smem,
        0,
    );
    println!("smems n={tot_smem}");
    for (i, smem) in match_array
        .iter()
        .take(usize::try_from(tot_smem).unwrap_or(0))
        .enumerate()
    {
        println!(
            "smem[{i}] rid={} m={} n={} k={} s={}",
            smem.rid, smem.m, smem.n, smem.k, smem.s
        );
        let mut coords = Vec::new();
        let mut coord_count = 0_i64;
        let mut id = 0_i64;
        worker
            .fmi
            .as_ref()
            .expect("worker fmi")
            .get_sa_entries_prefetch(
                std::slice::from_ref(smem),
                &mut coords,
                &mut coord_count,
                1,
                worker.opt.as_deref().expect("opt").max_occ,
                0,
                &mut id,
            );
        println!("  sa_count={coord_count} sa={coords:?}");
    }

    let mut pre_chain_ar = vec![mem_chain_v::default(); 1];
    let seed_buf_size = i64::try_from(worker.seedBuf.len()).expect("seedBuf len");
    mem_chain_seeds(
        worker.fmi.as_ref().expect("worker fmi"),
        worker.opt.as_deref().expect("opt"),
        worker
            .fmi
            .as_ref()
            .expect("worker fmi")
            .base
            .idx
            .bns
            .as_ref()
            .expect("bns"),
        &worker.seqs,
        1,
        0,
        &mut pre_chain_ar,
        &mut worker.seedBuf,
        seed_buf_size,
        &match_array,
        tot_smem,
    );
    println!("pre_filter chains n={}", pre_chain_ar[0].n);
    for (i, ch) in pre_chain_ar[0].a.iter().take(pre_chain_ar[0].n).enumerate() {
        println!(
            "pre_chain[{i}] rid={} pos={} w={} kept={} first={} n={} cseed={} frac_rep={:.6}",
            ch.rid, ch.pos, ch.w, ch.kept, ch.first, ch.n, ch.cseed, ch.frac_rep
        );
        for (j, s) in ch
            .seeds
            .iter()
            .take(usize::try_from(ch.n).unwrap_or(0))
            .enumerate()
            .take(8)
        {
            println!(
                "  pre_seed[{j}] rbeg={} qbeg={} len={} score={} rid={}",
                s.rbeg, s.qbeg, s.len, s.score, ch.rid
            );
        }
    }

    let opt_ref = worker.opt.as_deref().expect("opt");
    let bns = worker
        .fmi
        .as_ref()
        .expect("worker fmi")
        .base
        .idx
        .bns
        .as_ref()
        .expect("bns");
    let pac = &worker.fmi.as_ref().expect("worker fmi").base.idx.pac;
    let mut post_chain_ar = pre_chain_ar.clone();
    for chn in post_chain_ar.iter_mut() {
        chn.n = usize::try_from(mem_chain_flt(opt_ref, chn.n as i32, &mut chn.a, 0))
            .expect("filtered chain count");
        chn.a.truncate(chn.n);
        mem_flt_chained_seeds(opt_ref, bns, pac, &worker.seqs, chn.n as i32, &mut chn.a);
    }
    println!("post_filter chains n={}", post_chain_ar[0].n);
    for (i, ch) in post_chain_ar[0]
        .a
        .iter()
        .take(post_chain_ar[0].n)
        .enumerate()
    {
        println!(
            "post_chain[{i}] rid={} pos={} w={} kept={} first={} n={} cseed={} frac_rep={:.6}",
            ch.rid, ch.pos, ch.w, ch.kept, ch.first, ch.n, ch.cseed, ch.frac_rep
        );
        for (j, s) in ch
            .seeds
            .iter()
            .take(usize::try_from(ch.n).unwrap_or(0))
            .enumerate()
            .take(8)
        {
            println!(
                "  post_seed[{j}] rbeg={} qbeg={} len={} score={} rid={}",
                s.rbeg, s.qbeg, s.len, s.score, ch.rid
            );
        }
    }

    worker_bwt(&mut worker, 0, 1, 0);
    println!("chains n={}", worker.chain_ar[0].n);
    for (i, ch) in worker.chain_ar[0]
        .a
        .iter()
        .take(worker.chain_ar[0].n)
        .enumerate()
    {
        println!(
            "chain[{i}] rid={} pos={} w={} kept={} first={} n={} cseed={} frac_rep={:.6}",
            ch.rid, ch.pos, ch.w, ch.kept, ch.first, ch.n, ch.cseed, ch.frac_rep
        );
        for (j, s) in ch
            .seeds
            .iter()
            .take(usize::try_from(ch.n).unwrap_or(0))
            .enumerate()
            .take(8)
        {
            println!(
                "  seed[{j}] rbeg={} qbeg={} len={} score={} rid={}",
                s.rbeg, s.qbeg, s.len, s.score, ch.rid
            );
        }
    }

    worker_aln(&mut worker, 0, 1, 0);
    println!("regs n={}", worker.regs[0].n);
    for (i, reg) in worker.regs[0].a.iter().take(worker.regs[0].n).enumerate() {
        println!(
            "reg[{i}] score={} qb={} qe={} rb={} re={} rid={} truesc={} sub={} csub={} seedcov={} secondary={} secondary_all={}",
            reg.score, reg.qb, reg.qe, reg.rb, reg.re, reg.rid, reg.truesc, reg.sub, reg.csub, reg.seedcov, reg.secondary, reg.secondary_all
        );
    }

    let n_pri = mem_mark_primary_se(
        worker.opt.as_deref().expect("opt"),
        i32::try_from(worker.regs[0].n).expect("n"),
        &mut worker.regs[0].a,
        0,
    );
    println!("after_primary n_pri={n_pri}");
    for (i, reg) in worker.regs[0].a.iter().take(worker.regs[0].n).enumerate() {
        println!(
            "post_reg[{i}] score={} qb={} qe={} rb={} re={} rid={} truesc={} sub={} csub={} seedcov={} secondary={} secondary_all={}",
            reg.score, reg.qb, reg.qe, reg.rb, reg.re, reg.rid, reg.truesc, reg.sub, reg.csub, reg.seedcov, reg.secondary, reg.secondary_all
        );
    }
}
