//! Reusable metagenomic gene predictor with parallel model evaluation.
//!
//! Caches the 50 metagenomic models and uses a rayon thread pool
//! to score qualifying models in parallel.

use std::alloc::{alloc_zeroed, handle_alloc_error, Layout};
use std::os::raw::c_int;
use std::sync::Arc;

use rayon::{prelude::*, ThreadPool};

use crate::types::{Gene, Node, Training, MAX_GENES, MAX_SEQ, NUM_META};
use super::convert::gene_to_predicted;
use super::encode::SequenceBuffer;
use super::types::{PredictedGene, ProdigalConfig, ProdigalError};

use super::predict::{sort_nodes, validate_config};

/// Recommended stack size for Rayon pools used by `MetaPredictor`.
///
/// The low-level translated prediction routines use deep call stacks, so the
/// default Rust thread stack can be too small for worker threads.
pub const META_PREDICTOR_STACK_SIZE: usize = 32 * 1024 * 1024; // 32 MB

use crate::node::{add_nodes, reset_node_scores, score_nodes, record_overlapping_starts};
use crate::dprog::{dprog, eliminate_bad_genes};
use crate::gene::{add_genes, tweak_final_starts, record_gene_data};

/// Reusable metagenomic gene predictor.
///
/// Pre-loads 50 metagenomic models and evaluates qualifying models
/// in parallel using a rayon thread pool with large stacks.
pub struct MetaPredictor {
    pool: Arc<ThreadPool>,
    models: Arc<Vec<Box<Training>>>,
    config: ProdigalConfig,
}

impl MetaPredictor {
    /// Create a new predictor with default config.
    pub fn new() -> Result<Self, ProdigalError> {
        Self::with_config(ProdigalConfig::default())
    }

    /// Create a new predictor with custom config.
    pub fn with_config(config: ProdigalConfig) -> Result<Self, ProdigalError> {
        validate_config(&config)?;

        let pool = rayon::ThreadPoolBuilder::new()
            .stack_size(META_PREDICTOR_STACK_SIZE)
            .build()
            .map_err(|e| ProdigalError::Io(std::io::Error::new(
                std::io::ErrorKind::Other, e.to_string(),
            )))?;

        Self::with_config_and_thread_pool(config, Arc::new(pool))
    }

    /// Create a new predictor with default config and a shared Rayon thread pool.
    ///
    /// The supplied pool is reused for parallel metagenomic model evaluation.
    /// For large inputs, build the pool with at least `META_PREDICTOR_STACK_SIZE`
    /// stack bytes per worker thread.
    pub fn with_thread_pool(pool: Arc<ThreadPool>) -> Result<Self, ProdigalError> {
        Self::with_config_and_thread_pool(ProdigalConfig::default(), pool)
    }

    /// Create a new predictor with custom config and a shared Rayon thread pool.
    ///
    /// The supplied pool is reused for parallel metagenomic model evaluation.
    /// For large inputs, build the pool with at least `META_PREDICTOR_STACK_SIZE`
    /// stack bytes per worker thread.
    pub fn with_config_and_thread_pool(
        config: ProdigalConfig,
        pool: Arc<ThreadPool>,
    ) -> Result<Self, ProdigalError> {
        validate_config(&config)?;
        let models = Arc::new(load_meta_models());

        Ok(MetaPredictor { pool, models, config })
    }

    /// Predict genes in the given sequence.
    pub fn predict(&self, seq: &[u8]) -> Result<Vec<PredictedGene>, ProdigalError> {
        if seq.is_empty() {
            return Err(ProdigalError::EmptySequence);
        }
        if seq.len() > MAX_SEQ {
            return Err(ProdigalError::SequenceTooLong {
                length: seq.len(),
                max: MAX_SEQ,
            });
        }

        let seq = seq.to_vec();
        let models = Arc::clone(&self.models);
        let config = self.config.clone();
        let pool = Arc::clone(&self.pool);

        std::thread::Builder::new()
            .stack_size(META_PREDICTOR_STACK_SIZE)
            .spawn(move || pool.install(|| predict_parallel(&seq, &models, &config)))
            .expect("failed to spawn worker thread")
            .join()
            .expect("worker thread panicked")
    }
}

fn load_meta_models() -> Vec<Box<Training>> {
    let mut models: Vec<Box<Training>> = Vec::with_capacity(NUM_META);
    for i in 0..NUM_META {
        unsafe {
            let layout = Layout::new::<Training>();
            let ptr = alloc_zeroed(layout) as *mut Training;
            if ptr.is_null() {
                handle_alloc_error(layout);
            }
            crate::training_data::load_metagenome(i, ptr);
            models.push(Box::from_raw(ptr));
        }
    }
    models
}

/// Nodes built for one translation table, shared across models in that group.
struct TransTableGroup {
    /// Indices into the models array for qualifying (GC-filtered) models.
    model_indices: Vec<usize>,
    /// Template node array (built once by add_nodes + sort).
    nodes: Vec<Node>,
    /// Number of valid nodes.
    nn: c_int,
}

fn predict_parallel(
    seq: &[u8],
    models: &[Box<Training>],
    config: &ProdigalConfig,
) -> Result<Vec<PredictedGene>, ProdigalError> {
    let closed = if config.closed_ends { 1 } else { 0 };

    let mut buf = SequenceBuffer::new();
    let (slen, gc) = unsafe { buf.encode(seq, config.mask_n_runs) };
    if slen == 0 {
        return Err(ProdigalError::EmptySequence);
    }
    buf.ensure_node_capacity(slen);

    // GC window for model selection
    let mut low = 0.88495 * gc - 0.0102337;
    if low > 0.65 { low = 0.65; }
    let mut high = 0.86596 * gc + 0.1131991;
    if high < 0.35 { high = 0.35; }

    // Phase 1: Build node arrays per translation table group,
    // filtering models by GC range.
    let mut groups: Vec<TransTableGroup> = Vec::new();
    let mut nn: c_int = 0;

    for i in 0..NUM_META {
        let need_rebuild = i == 0 || models[i].trans_table != models[i - 1].trans_table;

        if need_rebuild {
            // Build nodes for this translation table
            // We need a mutable Training pointer but add_nodes only reads from it
            let mut tinf_copy = (*models[i]).clone();
            unsafe {
                buf.clear_nodes(nn);
                nn = add_nodes(
                    buf.seq.as_mut_ptr(), buf.rseq.as_mut_ptr(), slen,
                    buf.nodes.as_mut_ptr(), closed,
                    buf.masks.as_mut_ptr(), buf.nmask,
                    &mut tinf_copy,
                );
            }
            sort_nodes(&mut buf.nodes[..nn as usize]);

            groups.push(TransTableGroup {
                model_indices: Vec::new(),
                nodes: buf.nodes[..nn as usize].to_vec(),
                nn,
            });
        }

        // GC filter
        if models[i].gc >= low && models[i].gc <= high {
            groups.last_mut().unwrap().model_indices.push(i);
        }
    }

    // Phase 2: Score all qualifying models in parallel.
    // Each task gets its own copy of the node array.
    // Safety: seq/rseq buffers are immutable during parallel scoring,
    // and the usize round-trip preserves the pointer value.
    let seq_addr = buf.seq.as_ptr() as usize;
    let rseq_addr = buf.rseq.as_ptr() as usize;

    struct ModelScore {
        phase: usize,
        score: f64,
    }

    let best = groups.par_iter().flat_map(|group| {
        group.model_indices.par_iter().map(|&model_idx| {
            let mut nodes = group.nodes.clone();
            let nn = group.nn;
            // Clone the Training struct for this thread (needed for mutable API)
            let mut tinf = (*models[model_idx]).clone();

            unsafe {
                reset_node_scores(nodes.as_mut_ptr(), nn);
                score_nodes(
                    seq_addr as *mut u8, rseq_addr as *mut u8, slen,
                    nodes.as_mut_ptr(), nn, &mut tinf, closed, 1,
                );
                record_overlapping_starts(nodes.as_mut_ptr(), nn, &mut tinf, 1);
                let ipath = dprog(nodes.as_mut_ptr(), nn, &mut tinf, 1);
                if ipath < 0 || ipath >= nn {
                    return ModelScore { phase: model_idx, score: f64::NEG_INFINITY };
                }
                ModelScore {
                    phase: model_idx,
                    score: nodes[ipath as usize].score,
                }
            }
        })
    })
    .reduce(
        || ModelScore { phase: 0, score: f64::NEG_INFINITY },
        |a, b| if a.score >= b.score { a } else { b },
    );

    if best.score == f64::NEG_INFINITY {
        return Ok(Vec::new());
    }

    // Phase 3: Re-run the best model to extract genes.
    // Need to rebuild nodes for the best model's translation table.
    let mut tinf = (*models[best.phase]).clone();

    // Find the group that contains the best model to get its nodes
    let best_group = groups.iter().find(|g| g.model_indices.contains(&best.phase)).unwrap();
    let mut nodes = best_group.nodes.clone();
    let nn = best_group.nn;
    let mut genes: Vec<Gene> = vec![unsafe { std::mem::zeroed() }; MAX_GENES];

    unsafe {
        reset_node_scores(nodes.as_mut_ptr(), nn);
        score_nodes(
            buf.seq.as_mut_ptr(), buf.rseq.as_mut_ptr(), slen,
            nodes.as_mut_ptr(), nn, &mut tinf, closed, 1,
        );
        record_overlapping_starts(nodes.as_mut_ptr(), nn, &mut tinf, 1);
        let ipath = dprog(nodes.as_mut_ptr(), nn, &mut tinf, 1);
        eliminate_bad_genes(nodes.as_mut_ptr(), ipath, &mut tinf);
        let ng = add_genes(genes.as_mut_ptr(), nodes.as_mut_ptr(), ipath);
        tweak_final_starts(genes.as_mut_ptr(), ng, nodes.as_mut_ptr(), nn, &mut tinf);
        record_gene_data(genes.as_mut_ptr(), ng, nodes.as_mut_ptr(), &mut tinf, 1);

        let mut result = Vec::with_capacity(ng as usize);
        for i in 0..ng {
            result.push(gene_to_predicted(&genes[i as usize], nodes.as_ptr(), &tinf));
        }
        Ok(result)
    }
}
