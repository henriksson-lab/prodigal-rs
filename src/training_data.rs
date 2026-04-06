use crate::types::Training;
use std::ptr;

const TRAINING_SIZE: usize = std::mem::size_of::<Training>();

macro_rules! metagenome_blob {
    ($n:literal) => {
        include_bytes!(concat!("../data/metagenome_", $n, ".bin"))
    };
}

static METAGENOME_BLOBS: [&[u8; TRAINING_SIZE]; 50] = [
    metagenome_blob!("0"),
    metagenome_blob!("1"),
    metagenome_blob!("2"),
    metagenome_blob!("3"),
    metagenome_blob!("4"),
    metagenome_blob!("5"),
    metagenome_blob!("6"),
    metagenome_blob!("7"),
    metagenome_blob!("8"),
    metagenome_blob!("9"),
    metagenome_blob!("10"),
    metagenome_blob!("11"),
    metagenome_blob!("12"),
    metagenome_blob!("13"),
    metagenome_blob!("14"),
    metagenome_blob!("15"),
    metagenome_blob!("16"),
    metagenome_blob!("17"),
    metagenome_blob!("18"),
    metagenome_blob!("19"),
    metagenome_blob!("20"),
    metagenome_blob!("21"),
    metagenome_blob!("22"),
    metagenome_blob!("23"),
    metagenome_blob!("24"),
    metagenome_blob!("25"),
    metagenome_blob!("26"),
    metagenome_blob!("27"),
    metagenome_blob!("28"),
    metagenome_blob!("29"),
    metagenome_blob!("30"),
    metagenome_blob!("31"),
    metagenome_blob!("32"),
    metagenome_blob!("33"),
    metagenome_blob!("34"),
    metagenome_blob!("35"),
    metagenome_blob!("36"),
    metagenome_blob!("37"),
    metagenome_blob!("38"),
    metagenome_blob!("39"),
    metagenome_blob!("40"),
    metagenome_blob!("41"),
    metagenome_blob!("42"),
    metagenome_blob!("43"),
    metagenome_blob!("44"),
    metagenome_blob!("45"),
    metagenome_blob!("46"),
    metagenome_blob!("47"),
    metagenome_blob!("48"),
    metagenome_blob!("49"),
];

/// Load metagenome training data for bin `n` (0..49) into the given pointer.
pub unsafe fn load_metagenome(n: usize, tptr: *mut Training) {
    ptr::copy_nonoverlapping(
        METAGENOME_BLOBS[n].as_ptr(),
        tptr as *mut u8,
        TRAINING_SIZE,
    );
}
