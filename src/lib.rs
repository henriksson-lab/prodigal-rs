pub mod types;
pub mod reader;
pub mod bitmap;
pub mod training;
pub mod training_data;
pub mod sequence;
pub mod node;
pub mod dprog;
pub mod gene;
pub mod metagenomic;
pub mod pipeline;
pub mod api;

// Re-export the high-level API at crate root
pub use api::{
    predict, predict_meta, predict_meta_with, predict_with, train, train_with,
    MetaPredictor, PredictedGene, ProdigalConfig, ProdigalError, StartCodon, Strand, TrainingData,
};
