use std::path::Path;
use std::sync::Arc;

use crate::generated::bntseq_h::bntseq_t;
use crate::generated::bwa_h::bseq1_t;
use crate::generated::bwamem_cpp::{mem_opt_init, mem_process_seqs, with_current_rayon_pool};
use crate::generated::bwamem_h::{mem_opt_t, worker_t};
use crate::generated::fmi_search_cpp::FMI_search;
use crate::output::RunOutput;

const MEM_F_PE: i32 = 0x2;

pub type Result<T> = std::result::Result<T, String>;

pub struct MemReadPair<'a> {
    pub name: String,
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
        let mut fmi = FMI_search::ctor(prefix);
        fmi.load_index();
        if fmi.base.idx.bns.is_none() {
            return Err(format!("failed to load bwa-mem2 index from {}", prefix));
        }

        let mut opt = *mem_opt_init();
        opt.n_threads = i32::try_from(threads.max(1))
            .map_err(|_| format!("thread count is too large: {}", threads))?;
        opt.flag |= MEM_F_PE;

        Ok(Self {
            opt,
            worker: worker_t {
                fmi: Some(fmi),
                ..Default::default()
            },
            n_processed: 0,
            rayon_pool,
        })
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
        if pairs.is_empty() {
            return Ok(Vec::new());
        }

        let mut seqs = Vec::with_capacity(pairs.len() * 2);
        for pair in pairs {
            seqs.push(make_bseq(
                i32::try_from(seqs.len()).map_err(|_| "too many reads in batch".to_string())?,
                &pair.name,
                pair.r1,
                pair.q1,
            )?);
            seqs.push(make_bseq(
                i32::try_from(seqs.len()).map_err(|_| "too many reads in batch".to_string())?,
                &pair.name,
                pair.r2,
                pair.q2,
            )?);
        }

        let n = i32::try_from(seqs.len()).map_err(|_| "too many reads in batch".to_string())?;
        if let Some(pool) = self.rayon_pool.clone() {
            pool.install(|| {
                with_current_rayon_pool(|| {
                    mem_process_seqs(
                        &mut self.opt,
                        self.n_processed,
                        n,
                        &mut seqs,
                        None,
                        &mut self.worker,
                    );
                });
            });
        } else {
            mem_process_seqs(
                &mut self.opt,
                self.n_processed,
                n,
                &mut seqs,
                None,
                &mut self.worker,
            );
        }
        self.n_processed += i64::from(n);

        let mut sam = Vec::with_capacity(self.worker.seqs.len());
        for seq in &mut self.worker.seqs {
            if let Some(line) = seq.sam.take() {
                sam.push(line.into_string());
            }
        }
        self.worker.seqs.clear();
        Ok(sam)
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
        for line in self.align_pairs(pairs)? {
            output
                .stdout(format_args!("{}", line.trim_end_matches('\n')))
                .map_err(|e| format!("failed to write SAM record: {e}"))?;
        }
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
        name: Some(name.to_string().into_boxed_str()),
        comment: None,
        seq: Some(seq.to_string().into_boxed_str()),
        qual: Some(qual.to_string().into_boxed_str()),
        sam: None,
        seq_nt4: Vec::new(),
    })
}

#[cfg(test)]
mod tests {
    use super::MemAligner;
    use super::MemAlignerBuilder;
    use crate::output::{RunOutput, SharedWriterOutput};
    use rayon::ThreadPoolBuilder;
    use std::sync::Arc;

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
    fn output_capture_api_is_available() {
        let output = SharedWriterOutput::with_stream_labels(Vec::new());
        output.stdout(format_args!("@HD\tVN:1.6")).unwrap();
        output.stderr(format_args!("diagnostic")).unwrap();
        let text = String::from_utf8(output.into_inner().unwrap()).unwrap();
        assert_eq!(text, "[stdout] @HD\tVN:1.6\n[stderr] diagnostic\n");
    }
}
