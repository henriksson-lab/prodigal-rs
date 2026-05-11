//! Reusable metagenomic gene predictor.
//!
//! Caches the 50 metagenomic models and uses a rayon thread pool
//! with a large stack for prediction.

use std::alloc::{alloc_zeroed, handle_alloc_error, Layout};
use std::os::raw::c_int;
use std::sync::Arc;

use rayon::prelude::*;
use rayon::ThreadPool;

use super::convert::gene_to_predicted;
use super::encode::SequenceBuffer;
use super::types::{PredictedGene, ProdigalConfig, ProdigalError};
use crate::types::{Gene, Node, Training, MAX_GENES, MAX_SEQ, NUM_META};

use super::predict::validate_config;

/// Recommended stack size for Rayon pools used by `MetaPredictor`.
///
/// The low-level translated prediction routines use deep call stacks, so the
/// default Rust thread stack can be too small for worker threads.
pub const META_PREDICTOR_STACK_SIZE: usize = 32 * 1024 * 1024; // 32 MB

use crate::dprog::{dprog, eliminate_bad_genes};
use crate::gene::{add_genes, record_gene_data, tweak_final_starts};
use crate::node::{add_nodes, record_overlapping_starts, reset_node_scores, score_nodes};

/// Reusable metagenomic gene predictor.
///
/// Pre-loads 50 metagenomic models and evaluates qualifying models with
/// Prodigal-compatible state preservation across adjacent model bins.
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
            .map_err(|e| {
                ProdigalError::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    e.to_string(),
                ))
            })?;

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

        Ok(MetaPredictor {
            pool,
            models,
            config,
        })
    }

    /// Predict genes in the given sequence.
    pub fn predict(&self, seq: &[u8]) -> Result<Vec<PredictedGene>, ProdigalError> {
        validate_sequence(seq)?;

        self.pool
            .install(|| predict_parallel(seq, &self.models, &self.config))
    }

    /// Predict genes in a batch of sequences, preserving input order.
    pub fn predict_batch<S>(&self, seqs: &[S]) -> Result<Vec<Vec<PredictedGene>>, ProdigalError>
    where
        S: AsRef<[u8]> + Sync,
    {
        for seq in seqs {
            validate_sequence(seq.as_ref())?;
        }

        self.pool.install(|| {
            seqs.par_iter()
                .map(|seq| predict_parallel(seq.as_ref(), &self.models, &self.config))
                .collect()
        })
    }
}

fn validate_sequence(seq: &[u8]) -> Result<(), ProdigalError> {
    if seq.is_empty() {
        return Err(ProdigalError::EmptySequence);
    }
    if seq.len() > MAX_SEQ {
        return Err(ProdigalError::SequenceTooLong {
            length: seq.len(),
            max: MAX_SEQ,
        });
    }
    Ok(())
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
    if low > 0.65 {
        low = 0.65;
    }
    let mut high = 0.86596 * gc + 0.1131991;
    if high < 0.35 {
        high = 0.35;
    }

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
                    buf.seq.as_mut_ptr(),
                    buf.rseq.as_mut_ptr(),
                    slen,
                    buf.nodes.as_mut_ptr(),
                    closed,
                    buf.masks.as_mut_ptr(),
                    buf.nmask,
                    &mut tinf_copy,
                );
            }
            buf.nodes[..nn as usize]
                .sort_unstable_by(|a, b| a.ndx.cmp(&b.ndx).then(b.strand.cmp(&a.strand)));

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

    let mut best_score = f64::NEG_INFINITY;
    let mut best_nodes = Vec::new();
    let mut best_genes: Vec<Gene> = Vec::new();
    let mut best_tinf: Option<Training> = None;

    unsafe {
        for group in groups.iter_mut() {
            for &model_idx in &group.model_indices {
                let mut tinf = (*models[model_idx]).clone();
                reset_node_scores(group.nodes.as_mut_ptr(), group.nn);
                score_nodes(
                    buf.seq.as_mut_ptr(),
                    buf.rseq.as_mut_ptr(),
                    slen,
                    group.nodes.as_mut_ptr(),
                    group.nn,
                    &mut tinf,
                    closed,
                    1,
                );
                record_overlapping_starts(group.nodes.as_mut_ptr(), group.nn, &mut tinf, 1);
                let ipath = dprog(group.nodes.as_mut_ptr(), group.nn, &mut tinf, 1);
                if ipath < 0 || ipath >= group.nn {
                    continue;
                }

                let score = group.nodes[ipath as usize].score;
                if score > best_score {
                    best_score = score;
                    eliminate_bad_genes(group.nodes.as_mut_ptr(), ipath, &mut tinf);

                    let mut genes: Vec<Gene> = vec![std::mem::zeroed(); MAX_GENES];
                    let ng = add_genes(genes.as_mut_ptr(), group.nodes.as_mut_ptr(), ipath);
                    tweak_final_starts(
                        genes.as_mut_ptr(),
                        ng,
                        group.nodes.as_mut_ptr(),
                        group.nn,
                        &mut tinf,
                    );
                    genes.truncate(ng as usize);

                    best_nodes = group.nodes.clone();
                    best_genes = genes;
                    best_tinf = Some(tinf);
                }
            }
        }
    }

    let Some(mut tinf) = best_tinf else {
        return Ok(Vec::new());
    };

    unsafe {
        record_gene_data(
            best_genes.as_mut_ptr(),
            best_genes.len() as c_int,
            best_nodes.as_mut_ptr(),
            &mut tinf,
            1,
        );

        let mut result = Vec::with_capacity(best_genes.len());
        for gene in &best_genes {
            result.push(gene_to_predicted(
                gene,
                best_nodes.as_ptr(),
                &tinf,
                slen as usize,
            ));
        }
        Ok(result)
    }
}
