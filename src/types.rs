use std::os::raw::{c_char, c_int};

// Constants from the C headers
pub const MAX_SEQ: usize = 32_000_000;
pub const MAX_LINE: usize = 10_000;
pub const MAX_MASKS: usize = 5_000;
pub const MAX_GENES: usize = 30_000;
pub const STT_NOD: usize = 100_000;
pub const MIN_GENE: c_int = 90;
pub const MIN_EDGE_GENE: c_int = 60;
pub const MAX_SAM_OVLP: c_int = 60;
pub const ST_WINDOW: c_int = 60;
pub const OPER_DIST: c_int = 60;
pub const EDGE_BONUS: f64 = 0.74;
pub const EDGE_UPS: f64 = -1.00;
pub const META_PEN: f64 = 7.5;
pub const NUM_META: usize = 50;
pub const WINDOW: usize = 120;
pub const MASK_SIZE: usize = 50;

// Node type constants
pub const ATG: c_int = 0;
pub const GTG: c_int = 1;
pub const TTG: c_int = 2;
pub const STOP: c_int = 3;

// dprog.h constants
pub const MAX_OPP_OVLP: c_int = 200;
pub const MAX_NODE_DIST: c_int = 500;

/// Masked region in a sequence (sequence.h)
#[repr(C)]
#[derive(Debug, Clone, Copy, Default)]
pub struct Mask {
    pub begin: c_int,
    pub end: c_int,
}

/// Training parameters learned from a genome (training.h)
#[repr(C)]
#[derive(Clone)]
pub struct Training {
    pub gc: f64,
    pub trans_table: c_int,
    pub st_wt: f64,
    pub bias: [f64; 3],
    pub type_wt: [f64; 3],
    pub uses_sd: c_int,
    pub rbs_wt: [f64; 28],
    pub ups_comp: [[f64; 4]; 32],
    pub mot_wt: [[[f64; 4096]; 4]; 4],
    pub no_mot: f64,
    pub gene_dc: [f64; 4096],
}

/// Upstream motif information (node.h)
#[repr(C)]
#[derive(Debug, Clone, Copy, Default)]
pub struct Motif {
    pub ndx: c_int,
    pub len: c_int,
    pub spacer: c_int,
    pub spacendx: c_int,
    pub score: f64,
}

/// A node in the dynamic programming graph (node.h)
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct Node {
    pub type_: c_int,
    pub edge: c_int,
    pub ndx: c_int,
    pub strand: c_int,
    pub stop_val: c_int,
    pub star_ptr: [c_int; 3],
    pub gc_bias: c_int,
    pub gc_score: [f64; 3],
    pub cscore: f64,
    pub gc_cont: f64,
    pub rbs: [c_int; 2],
    pub mot: Motif,
    pub uscore: f64,
    pub tscore: f64,
    pub rscore: f64,
    pub sscore: f64,
    pub traceb: c_int,
    pub tracef: c_int,
    pub ov_mark: c_int,
    pub score: f64,
    pub elim: c_int,
}

/// A predicted gene (gene.h)
#[repr(C)]
#[derive(Clone)]
pub struct Gene {
    pub begin: c_int,
    pub end: c_int,
    pub start_ndx: c_int,
    pub stop_ndx: c_int,
    pub gene_data: [c_char; 500],
    pub score_data: [c_char; 500],
}

/// A metagenomic model bin (metagenomic.h)
#[repr(C)]
pub struct MetagenomicBin {
    pub index: c_int,
    pub clusnum: c_int,
    pub desc: [c_char; 500],
    pub weight: f64,
    pub gc: f64,
    pub tinf: *mut Training,
}
