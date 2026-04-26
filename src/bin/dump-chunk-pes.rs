use bwa_mem2_rs::generated::bwa_cpp::bseq_read_orig;
use bwa_mem2_rs::generated::bwamem_cpp::mem_process_seqs;
use bwa_mem2_rs::generated::bwamem_cpp::mem_opt_init;
use bwa_mem2_rs::generated::bwamem_h::{mem_alnreg_v, mem_chain_v, mem_seed_t, worker_t};
use bwa_mem2_rs::generated::fastmap_cpp::memoryAlloc;
use bwa_mem2_rs::generated::fastmap_h::ktp_aux_t;
use bwa_mem2_rs::generated::fmi_search_cpp::FMI_search;
use bwa_mem2_rs::generated::kseq_h::kseq_t;

fn main() {
    let mut args = std::env::args().skip(1);
    let threads: i32 = args.next().expect("threads").parse().expect("threads");
    let prefix = args.next().expect("prefix");
    let r1 = args.next().expect("r1");
    let r2 = args.next().expect("r2");

    let r1_text = std::fs::read_to_string(&r1).expect("read r1");
    let r2_text = std::fs::read_to_string(&r2).expect("read r2");
    let mut ks1 = kseq_t::from_text(&r1_text);
    let mut ks2 = kseq_t::from_text(&r2_text);

    let mut fmi = FMI_search::ctor(&prefix);
    fmi.load_index();

    let mut opt = (*mem_opt_init()).clone();
    opt.flag |= 0x2;
    opt.n_threads = threads;

    let task_size = opt.chunk_size * i64::from(opt.n_threads);
    let aux = ktp_aux_t { opt: Some(Box::new(opt.clone())), task_size, actual_chunk_size: task_size, ..Default::default() };
    let nreads = i32::try_from(task_size / 151 + 10).expect("nreads");
    let mut worker = worker_t { fmi: Some(fmi), ..Default::default() };
    memoryAlloc(&aux, &mut worker, nreads, threads);

    let mut n_processed = 0_i64;
    let mut batch_id = 0_i32;
    loop {
        let mut n = 0_i32;
        let mut total_bases = 0_i64;
        let mut seqs = bseq_read_orig(task_size, &mut n, &mut ks1, Some(&mut ks2), &mut total_bases);
        if n == 0 {
            break;
        }
        if worker.nreads < n {
            worker.nreads = n;
            let nreads = usize::try_from(n).expect("nreads");
            worker.regs = vec![mem_alnreg_v::default(); nreads];
            worker.chain_ar = vec![mem_chain_v::default(); nreads];
            worker.seedBuf = vec![mem_seed_t::default(); nreads * 64];
        }
        mem_process_seqs(&mut opt, n_processed, n, &mut seqs, None, &mut worker);
        println!("batch={batch_id} n_seqs={n} total_bases={total_bases}");
        for (i, pe) in worker.pes.iter().enumerate() {
            println!(
                "  pes[{i}]: low={} high={} failed={} avg={:.10} std={:.10}",
                pe.low, pe.high, pe.failed, pe.avg, pe.std
            );
        }
        n_processed += i64::from(n);
        batch_id += 1;
    }
}
