use std::path::Path;
use std::sync::Arc;

use crate::bwa_mem2::bntseq::bntseq_t;
use crate::bwa_mem2::bwa::bseq1_t;
use crate::bwa_mem2::bwamem::{mem_opt_init, with_current_rayon_pool, MEM_F_PE};
use crate::bwa_mem2::bwamem::{mem_opt_t, worker_t};
use crate::bwa_mem2::fastmap::{ktp_aux_t, ktp_data_t, memory_alloc, process_batch};
use crate::bwa_mem2::fmi_search::FMI_search;
use crate::output::RunOutput;

pub type Result<T> = std::result::Result<T, String>;

pub struct MemReadPair<'a> {
    pub name: &'a str,
    pub r1: &'a [u8],
    pub q1: &'a [u8],
    pub r2: &'a [u8],
    pub q2: &'a [u8],
}

pub struct MemAligner {
    opt: mem_opt_t,
    worker: worker_t,
    n_processed: i64,
    rayon_pool: Option<Arc<rayon::ThreadPool>>,
}

pub struct MemAlignerBuilder<'a> {
    index_prefix: &'a Path,
    threads: usize,
    rayon_pool: Option<Arc<rayon::ThreadPool>>,
}

impl<'a> MemAlignerBuilder<'a> {
    pub fn threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }

    pub fn thread_pool(mut self, rayon_pool: Arc<rayon::ThreadPool>) -> Self {
        self.rayon_pool = Some(rayon_pool);
        self
    }

    pub fn build(self) -> Result<MemAligner> {
        MemAligner::new_inner(self.index_prefix, self.threads, self.rayon_pool)
    }
}

impl MemAligner {
    pub fn builder(index_prefix: &Path) -> MemAlignerBuilder<'_> {
        MemAlignerBuilder {
            index_prefix,
            threads: 1,
            rayon_pool: None,
        }
    }

    pub fn new(index_prefix: &Path, threads: usize) -> Result<Self> {
        Self::builder(index_prefix).threads(threads).build()
    }

    pub fn new_with_thread_pool(
        index_prefix: &Path,
        threads: usize,
        rayon_pool: Arc<rayon::ThreadPool>,
    ) -> Result<Self> {
        Self::builder(index_prefix)
            .threads(threads)
            .thread_pool(rayon_pool)
            .build()
    }

    fn new_inner(
        index_prefix: &Path,
        threads: usize,
        rayon_pool: Option<Arc<rayon::ThreadPool>>,
    ) -> Result<Self> {
        let prefix = index_prefix
            .to_str()
            .ok_or_else(|| format!("index path is not valid UTF-8: {}", index_prefix.display()))?;
        for suffix in [".bwt.2bit.64", ".ann", ".amb", ".pac"] {
            let path = format!("{prefix}{suffix}");
            if !Path::new(&path).is_file() {
                return Err(format!(
                    "failed to load bwa-mem2 index from {prefix}: missing {path}"
                ));
            }
        }
        let mut fmi = FMI_search::ctor(prefix);
        fmi.load_index();
        if fmi.base.idx.bns.is_none() {
            return Err(format!("failed to load bwa-mem2 index from {}", prefix));
        }

        let mut opt = *mem_opt_init();
        opt.n_threads = i32::try_from(threads.max(1))
            .map_err(|_| format!("thread count is too large: {}", threads))?;
        opt.flag |= MEM_F_PE;

        let mut worker = worker_t {
            fmi: Some(fmi),
            ..Default::default()
        };
        worker.nthreads = i16::try_from(opt.n_threads.max(1)).expect("nthreads");
        let actual_chunk_size = opt.chunk_size * i64::from(opt.n_threads.max(1));
        let nreads = i32::try_from(actual_chunk_size / 151 + 10)
            .map_err(|_| format!("estimated read capacity is too large: {actual_chunk_size}"))?;
        let aux = ktp_aux_t {
            opt: Some(Box::new(opt.clone())),
            task_size: actual_chunk_size,
            actual_chunk_size,
            ..Default::default()
        };
        memory_alloc(&aux, &mut worker, nreads, opt.n_threads);

        Ok(Self {
            opt,
            worker,
            n_processed: 0,
            rayon_pool,
        })
    }

    pub fn opt(&self) -> &mem_opt_t {
        &self.opt
    }

    pub fn opt_mut(&mut self) -> &mut mem_opt_t {
        &mut self.opt
    }

    pub fn sam_header(&self) -> Result<String> {
        let bns = self.bns()?;
        let mut out = String::new();
        for ann in &bns.anns {
            out.push_str("@SQ\tSN:");
            out.push_str(&ann.name);
            out.push_str("\tLN:");
            out.push_str(&ann.len.to_string());
            if ann.is_alt != 0 {
                out.push_str("\tAH:*");
            }
            out.push('\n');
        }
        out.push_str("@PG\tID:bwa-mem2-rs\tPN:bwa-mem2-rs\n");
        Ok(out)
    }

    pub fn align_pairs(&mut self, pairs: &[MemReadPair<'_>]) -> Result<Vec<String>> {
        let mut sam = Vec::new();
        self.align_pairs_into(pairs, |line| {
            sam.push(line.to_owned());
            Ok(())
        })?;
        Ok(sam)
    }

    /// User-approved deviation from CCC: this library callback API exposes generated SAM
    /// records without materializing a Vec, while leaving the translated core unchanged.
    pub fn align_pairs_into<F>(&mut self, pairs: &[MemReadPair<'_>], mut sink: F) -> Result<()>
    where
        F: FnMut(&str) -> Result<()>,
    {
        self.align_pairs_into_indexed(pairs, |_pair_index, line| sink(line))
    }

    /// User-approved deviation from CCC: the callback includes the input pair index so callers can
    /// keep metadata out of QNAMEs and recover one-read-to-many-output associations explicitly.
    pub fn align_pairs_into_indexed<F>(
        &mut self,
        pairs: &[MemReadPair<'_>],
        mut sink: F,
    ) -> Result<()>
    where
        F: FnMut(usize, &str) -> Result<()>,
    {
        if pairs.is_empty() {
            return Ok(());
        }

        let mut seqs = Vec::with_capacity(pairs.len() * 2);
        for pair in pairs {
            seqs.push(make_bseq(
                i32::try_from(seqs.len()).map_err(|_| "too many reads in batch".to_string())?,
                pair.name,
                pair.r1,
                pair.q1,
            )?);
            seqs.push(make_bseq(
                i32::try_from(seqs.len()).map_err(|_| "too many reads in batch".to_string())?,
                pair.name,
                pair.r2,
                pair.q2,
            )?);
        }

        let n = i32::try_from(seqs.len()).map_err(|_| "too many reads in batch".to_string())?;
        let mut batch = ktp_data_t {
            n_seqs: n,
            seqs,
            sam_lines: Vec::new(),
        };
        if let Some(pool) = self.rayon_pool.clone() {
            pool.install(|| {
                with_current_rayon_pool(|| {
                    process_batch(
                        &mut batch,
                        &mut self.opt,
                        self.n_processed,
                        None,
                        &mut self.worker,
                    );
                });
            });
        } else {
            process_batch(
                &mut batch,
                &mut self.opt,
                self.n_processed,
                None,
                &mut self.worker,
            );
        }
        self.n_processed += i64::from(n);

        for (read_index, line) in batch.sam_lines.drain(..).enumerate() {
            sink(read_index / 2, &line)?;
        }
        Ok(())
    }

    pub fn write_sam_for_pairs(
        &mut self,
        pairs: &[MemReadPair<'_>],
        output: &dyn RunOutput,
    ) -> Result<()> {
        for line in self.sam_header()?.lines() {
            output
                .stdout(format_args!("{line}"))
                .map_err(|e| format!("failed to write SAM header: {e}"))?;
        }
        self.align_pairs_into(pairs, |line| {
            output
                .stdout(format_args!("{}", line.trim_end_matches('\n')))
                .map_err(|e| format!("failed to write SAM record: {e}"))
        })?;
        Ok(())
    }

    fn bns(&self) -> Result<&bntseq_t> {
        self.worker
            .fmi
            .as_ref()
            .and_then(|fmi| fmi.base.idx.bns.as_ref())
            .ok_or_else(|| "bwa-mem2 index is not loaded".to_string())
    }
}

impl Drop for MemAligner {
    fn drop(&mut self) {
        if let Some(fmi) = self.worker.fmi.as_mut() {
            fmi.dtor();
        }
    }
}

/// User-approved deviation from CCC: construct translated `bseq1_t` records from in-memory
/// library inputs instead of the original FASTQ reader path.
fn make_bseq(id: i32, name: &str, seq: &[u8], qual: &[u8]) -> Result<bseq1_t> {
    if seq.len() != qual.len() {
        return Err(format!(
            "sequence/quality length mismatch for {}: {} != {}",
            name,
            seq.len(),
            qual.len()
        ));
    }
    let seq = std::str::from_utf8(seq)
        .map_err(|e| format!("read sequence for {} is not valid UTF-8: {}", name, e))?;
    let qual = std::str::from_utf8(qual)
        .map_err(|e| format!("read qualities for {} are not valid UTF-8: {}", name, e))?;
    Ok(bseq1_t {
        l_seq: i32::try_from(seq.len()).map_err(|_| format!("read is too long for {}", name))?,
        id,
        name: Some(name.to_string()),
        comment: None,
        seq: Some(seq.to_string()),
        qual: Some(qual.to_string()),
        sam: None,
    })
}

#[cfg(test)]
mod tests {
    use super::MemAligner;
    use super::MemAlignerBuilder;
    use super::MemReadPair;
    use crate::bwa_mem2::bntseq::bns_fasta2bntseq;
    use crate::bwa_mem2::fastmap::main_mem;
    use crate::bwa_mem2::fmi_search::FMI_search;
    use crate::output::{RunOutput, SharedWriterOutput};
    use rayon::ThreadPoolBuilder;
    use std::fs;
    use std::io::{Cursor, Write};
    use std::sync::Arc;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[test]
    fn builder_api_is_available() {
        let _builder: fn(&std::path::Path) -> MemAlignerBuilder<'_> = MemAligner::builder;
        let pool = Arc::new(
            ThreadPoolBuilder::new()
                .num_threads(1)
                .build()
                .expect("thread pool"),
        );
        let _builder = MemAligner::builder(std::path::Path::new("index"))
            .threads(1)
            .thread_pool(pool.clone());
        assert_eq!(pool.current_num_threads(), 1);
    }

    #[test]
    fn legacy_new_with_thread_pool_is_available() {
        let _ctor: fn(
            &std::path::Path,
            usize,
            Arc<rayon::ThreadPool>,
        ) -> super::Result<MemAligner> = MemAligner::new_with_thread_pool;
    }

    #[test]
    fn new_returns_err_for_missing_index_prefix() {
        let mut path = std::env::temp_dir();
        path.push(format!(
            "bwa_mem2_rs_missing_index_{}",
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .expect("clock")
                .as_nanos()
        ));
        let err = match MemAligner::new(&path, 1) {
            Ok(_) => panic!("missing index should return Err"),
            Err(err) => err,
        };
        assert!(err.contains("failed to load bwa-mem2 index"), "{err}");
    }

    #[test]
    fn output_capture_api_is_available() {
        let output = SharedWriterOutput::with_stream_labels(Vec::new());
        output.stdout(format_args!("@HD\tVN:1.6")).unwrap();
        output.stderr(format_args!("diagnostic")).unwrap();
        let text = String::from_utf8(output.into_inner().unwrap()).unwrap();
        assert_eq!(text, "[stdout] @HD\tVN:1.6\n[stderr] diagnostic\n");
    }

    #[test]
    fn builder_shared_pool_and_output_capture_write_sam() {
        let dir = unique_temp_dir("bwa_mem2_rs_mem_api");
        fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("ref");
        let prefix_str = prefix.to_str().expect("utf8 path");

        let reference = deterministic_reference(800);
        let fasta = format!(">chr1\n{reference}\n");
        let l_pac = bns_fasta2bntseq(Cursor::new(fasta.as_bytes()), prefix_str, 1);
        assert_eq!(l_pac, 800);
        let mut fmi = FMI_search::ctor(prefix_str);
        assert_eq!(fmi.build_index(), 0);
        fmi.dtor();

        let pool = Arc::new(
            ThreadPoolBuilder::new()
                .num_threads(1)
                .build()
                .expect("thread pool"),
        );
        let mut aligner = MemAligner::builder(&prefix)
            .threads(1)
            .thread_pool(pool)
            .build()
            .expect("load generated index");
        {
            let opt = aligner.opt_mut();
            opt.flag &= !super::MEM_F_PE;
            opt.min_seed_len = 2;
            opt.split_factor = 1.5;
            opt.split_width = 10;
            opt.max_mem_intv = 20;
            opt.min_chain_weight = 1;
            opt.T = 1;
        }
        let mut pairs = Vec::new();
        let qual = b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        let names: Vec<String> = (0..20).map(|i| format!("read-lib-{i}")).collect();
        for i in 0..20 {
            let start = 10 + i * 9;
            pairs.push(MemReadPair {
                name: names[i].as_str(),
                r1: reference[start..start + 60].as_bytes(),
                q1: qual,
                r2: reference[start + 140..start + 200].as_bytes(),
                q2: qual,
            });
        }
        let output = SharedWriterOutput::with_stream_labels(Vec::new());

        aligner
            .write_sam_for_pairs(&pairs, &output)
            .expect("write SAM");

        let text = String::from_utf8(output.into_inner().unwrap()).unwrap();
        assert!(text.contains("[stdout] @SQ\tSN:chr1\tLN:800"), "{text}");
        assert!(
            text.contains("[stdout] @PG\tID:bwa-mem2-rs\tPN:bwa-mem2-rs"),
            "{text}"
        );
        assert!(text.contains("[stdout] read-lib-0\t"), "{text}");
        assert!(text.contains("\tchr1\t"), "{text}");

        fs::remove_dir_all(dir).expect("remove temp dir");
    }

    #[test]
    fn mem_aligner_matches_stock_fastmap_process_for_paired_batch() {
        let dir = unique_temp_dir("bwa_mem2_rs_mem_api_parity");
        fs::create_dir_all(&dir).expect("create temp dir");
        let prefix = dir.join("ref");
        let prefix_str = prefix.to_str().expect("utf8 path");

        let reference = deterministic_reference(1000);
        let fasta = format!(">chr1\n{reference}\n");
        assert_eq!(
            bns_fasta2bntseq(Cursor::new(fasta.as_bytes()), prefix_str, 1),
            1000
        );
        let mut fmi = FMI_search::ctor(prefix_str);
        assert_eq!(fmi.build_index(), 0);
        fmi.dtor();

        let reads1 = dir.join("reads_1.fq");
        let reads2 = dir.join("reads_2.fq");
        let out = dir.join("stock.sam");
        let qual = b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        let names: Vec<String> = (0..12).map(|i| format!("read-parity-{i}")).collect();
        let mut pairs = Vec::new();
        {
            let mut fq1 = fs::File::create(&reads1).expect("create reads1");
            let mut fq2 = fs::File::create(&reads2).expect("create reads2");
            for i in 0..names.len() {
                let start = 20 + i * 13;
                let r1 = reference[start..start + 60].as_bytes();
                let r2 = reference[start + 180..start + 240].as_bytes();
                writeln!(
                    fq1,
                    "@{}\n{}\n+\n{}",
                    names[i],
                    std::str::from_utf8(r1).expect("r1 utf8"),
                    std::str::from_utf8(qual).expect("qual utf8")
                )
                .expect("write reads1");
                writeln!(
                    fq2,
                    "@{}\n{}\n+\n{}",
                    names[i],
                    std::str::from_utf8(r2).expect("r2 utf8"),
                    std::str::from_utf8(qual).expect("qual utf8")
                )
                .expect("write reads2");
                pairs.push(MemReadPair {
                    name: names[i].as_str(),
                    r1,
                    q1: qual,
                    r2,
                    q2: qual,
                });
            }
        }

        let argv = vec![
            "mem".to_string(),
            "-o".to_string(),
            out.to_str().expect("utf8").to_string(),
            "-t".to_string(),
            "1".to_string(),
            "-k".to_string(),
            "2".to_string(),
            "-W".to_string(),
            "1".to_string(),
            "-T".to_string(),
            "1".to_string(),
            prefix_str.to_string(),
            reads1.to_str().expect("utf8").to_string(),
            reads2.to_str().expect("utf8").to_string(),
        ];
        assert_eq!(main_mem(&argv), 0);
        let stock_body = fs::read_to_string(&out)
            .expect("read stock sam")
            .lines()
            .filter(|line| !line.starts_with('@'))
            .map(|line| line.to_string())
            .collect::<Vec<_>>();

        let mut aligner = MemAligner::new(&prefix, 1).expect("load generated index");
        {
            let opt = aligner.opt_mut();
            opt.min_seed_len = 2;
            opt.min_chain_weight = 1;
            opt.T = 1;
        }
        let embedded_body = aligner
            .align_pairs(&pairs)
            .expect("align pairs")
            .into_iter()
            .map(|line| line.trim_end_matches('\n').to_owned())
            .collect::<Vec<_>>();

        assert_eq!(embedded_body, stock_body);
        fs::remove_dir_all(dir).expect("remove temp dir");
    }

    fn unique_temp_dir(prefix: &str) -> std::path::PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system time")
            .as_nanos();
        std::env::temp_dir().join(format!("{prefix}_{}_{}", std::process::id(), nanos))
    }

    fn deterministic_reference(len: usize) -> String {
        const BASES: &[u8; 4] = b"ACGT";
        let mut out = String::with_capacity(len);
        let mut state = 0x1357_2468_u32;
        for _ in 0..len {
            state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            out.push(BASES[((state >> 29) & 3) as usize] as char);
        }
        out
    }
}
