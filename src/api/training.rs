use crate::types::Training;
use super::types::ProdigalError;
use std::io::{Read, Write};
use std::path::Path;

/// Opaque wrapper around learned training parameters.
///
/// Can be saved/loaded to reuse training across runs.
pub struct TrainingData {
    pub(crate) inner: Box<Training>,
}

impl TrainingData {
    /// Load training data from a Prodigal binary training file.
    pub fn load(path: impl AsRef<Path>) -> Result<Self, ProdigalError> {
        let mut file = std::fs::File::open(path)?;
        let mut inner = Box::new(unsafe { std::mem::zeroed::<Training>() });
        let size = std::mem::size_of::<Training>();
        let buf = unsafe { std::slice::from_raw_parts_mut(&mut *inner as *mut Training as *mut u8, size) };
        file.read_exact(buf)?;
        Ok(TrainingData { inner })
    }

    /// Save training data to a file.
    pub fn save(&self, path: impl AsRef<Path>) -> Result<(), ProdigalError> {
        let mut file = std::fs::File::create(path)?;
        let size = std::mem::size_of::<Training>();
        let buf = unsafe { std::slice::from_raw_parts(&*self.inner as *const Training as *const u8, size) };
        file.write_all(buf)?;
        Ok(())
    }

    /// GC content this model was trained on.
    pub fn gc(&self) -> f64 {
        self.inner.gc
    }

    /// Translation table.
    pub fn translation_table(&self) -> i32 {
        self.inner.trans_table
    }

    /// Whether this model uses Shine-Dalgarno RBS scoring.
    pub fn uses_sd(&self) -> bool {
        self.inner.uses_sd != 0
    }
}

impl Clone for TrainingData {
    fn clone(&self) -> Self {
        TrainingData {
            inner: self.inner.clone(),
        }
    }
}
