/*******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Copyright (C) 2007-2016 University of Tennessee / UT-Battelle

    Code Author:  Doug Hyatt

    Rust port of main.c — the pipeline orchestration module.
*******************************************************************************/

use std::ffi::CStr;
use std::os::raw::{c_char, c_int, c_void};
use std::os::unix::io::{FromRawFd, IntoRawFd};

use crate::types::{Gene, Mask, MetagenomicBin, Node, Training};
use crate::types::{MAX_GENES, MAX_LINE, MAX_MASKS, MAX_SEQ, NUM_META, STT_NOD};

const VERSION: &str = "2.6.3";
const DATE: &str = "February, 2016";
const VERSION_CSTR: &[u8] = b"2.6.3\0";
const MIN_SINGLE_GENOME: c_int = 20000;
const IDEAL_SINGLE_GENOME: c_int = 100000;

// ---------------------------------------------------------------------------
// Rust reader functions (reader.rs, replaces zlib gzopen/gzgets/gzseek/gzclose)
// ---------------------------------------------------------------------------
extern "C" {
    fn seq_reader_open(path: *const c_char) -> *mut c_void;
    fn seq_reader_close(file: *mut c_void) -> c_int;
    fn seq_reader_seek(
        file: *mut c_void,
        offset: i64,
        whence: c_int,
    ) -> i64;
}

// ---------------------------------------------------------------------------
// External C functions from the other Rust modules (all #[no_mangle] extern "C")
// ---------------------------------------------------------------------------
extern "C" {
    // sequence.rs
    fn read_seq_training(
        fp: *mut c_void,
        seq: *mut u8,
        useq: *mut u8,
        gc: *mut f64,
        do_mask: c_int,
        mlist: *mut Mask,
        nmask: *mut c_int,
    ) -> c_int;

    fn next_seq_multi(
        fp: *mut c_void,
        seq: *mut u8,
        useq: *mut u8,
        num_seq: *mut c_int,
        gc: *mut f64,
        do_mask: c_int,
        mlist: *mut Mask,
        nmask: *mut c_int,
        cur_header: *mut c_char,
        new_header: *mut c_char,
    ) -> c_int;

    fn rcom_seq(seq: *mut u8, rseq: *mut u8, useq: *mut u8, slen: c_int);

    fn calc_short_header(header: *mut c_char, short_header: *mut c_char, num_seq: c_int);

    fn calc_most_gc_frame(seq: *mut u8, slen: c_int) -> *mut c_int;

    // node.rs
    fn add_nodes(
        seq: *mut u8,
        rseq: *mut u8,
        slen: c_int,
        nodes: *mut Node,
        closed: c_int,
        mlist: *mut Mask,
        nmask: c_int,
        tinf: *mut Training,
    ) -> c_int;

    fn record_gc_bias(gc: *mut c_int, nodes: *mut Node, nn: c_int, tinf: *mut Training);

    fn calc_dicodon_gene(
        tinf: *mut Training,
        seq: *mut u8,
        rseq: *mut u8,
        slen: c_int,
        nodes: *mut Node,
        ipath: c_int,
    );

    fn raw_coding_score(
        seq: *mut u8,
        rseq: *mut u8,
        slen: c_int,
        nodes: *mut Node,
        nn: c_int,
        tinf: *mut Training,
    );

    fn rbs_score(
        seq: *mut u8,
        rseq: *mut u8,
        slen: c_int,
        nodes: *mut Node,
        nn: c_int,
        tinf: *mut Training,
    );

    fn train_starts_sd(
        seq: *mut u8,
        rseq: *mut u8,
        slen: c_int,
        nodes: *mut Node,
        nn: c_int,
        tinf: *mut Training,
    );

    fn train_starts_nonsd(
        seq: *mut u8,
        rseq: *mut u8,
        slen: c_int,
        nodes: *mut Node,
        nn: c_int,
        tinf: *mut Training,
    );

    fn determine_sd_usage(tinf: *mut Training);

    fn record_overlapping_starts(nodes: *mut Node, nn: c_int, tinf: *mut Training, stage: c_int);

    fn score_nodes(
        seq: *mut u8,
        rseq: *mut u8,
        slen: c_int,
        nodes: *mut Node,
        nn: c_int,
        tinf: *mut Training,
        closed: c_int,
        is_meta: c_int,
    );

    fn reset_node_scores(nod: *mut Node, nn: c_int);

    fn write_start_file(
        fh: c_int,
        nod: *mut Node,
        nn: c_int,
        tinf: *mut Training,
        sctr: c_int,
        slen: c_int,
        is_meta: c_int,
        mdesc: *mut c_char,
        version: *mut c_char,
        header: *mut c_char,
    );

    // dprog.rs
    fn dprog(nod: *mut Node, nn: c_int, tinf: *mut Training, final_: c_int) -> c_int;

    fn eliminate_bad_genes(nod: *mut Node, dbeg: c_int, tinf: *mut Training);

    // gene.rs
    fn add_genes(glist: *mut Gene, nod: *mut Node, dbeg: c_int) -> c_int;

    fn tweak_final_starts(
        genes: *mut Gene,
        ng: c_int,
        nod: *mut Node,
        nn: c_int,
        tinf: *mut Training,
    );

    fn record_gene_data(
        genes: *mut Gene,
        ng: c_int,
        nod: *mut Node,
        tinf: *mut Training,
        sctr: c_int,
    );

    fn print_genes(
        fp: c_int,
        genes: *mut Gene,
        ng: c_int,
        nod: *mut Node,
        slen: c_int,
        flag: c_int,
        sctr: c_int,
        is_meta: c_int,
        mdesc: *mut c_char,
        tinf: *mut Training,
        header: *mut c_char,
        short_hdr: *mut c_char,
        version: *mut c_char,
    );

    fn write_translations(
        fh: c_int,
        genes: *mut Gene,
        ng: c_int,
        nod: *mut Node,
        seq: *mut u8,
        rseq: *mut u8,
        useq: *mut u8,
        slen: c_int,
        tinf: *mut Training,
        sctr: c_int,
        short_hdr: *mut c_char,
    );

    fn write_nucleotide_seqs(
        fh: c_int,
        genes: *mut Gene,
        ng: c_int,
        nod: *mut Node,
        seq: *mut u8,
        rseq: *mut u8,
        useq: *mut u8,
        slen: c_int,
        tinf: *mut Training,
        sctr: c_int,
        short_hdr: *mut c_char,
    );

    // training.rs
    fn read_training_file(path: *const c_char, tinf: *mut Training) -> c_int;

    fn write_training_file(path: *const c_char, tinf: *mut Training) -> c_int;

    // metagenomic.rs
    fn initialize_metagenomic_bins(meta: *mut MetagenomicBin);
}

// ---------------------------------------------------------------------------
// Helper: sort a Node slice by (ndx, -strand) matching compare_nodes logic
// ---------------------------------------------------------------------------
unsafe fn sort_nodes(nodes: &mut [Node]) {
    nodes.sort_unstable_by(|a, b| {
        a.ndx.cmp(&b.ndx).then(b.strand.cmp(&a.strand))
    });
}

// ---------------------------------------------------------------------------
// Helper: print version and exit
// ---------------------------------------------------------------------------
fn version() -> ! {
    eprintln!("\nProdigal V{}: {}\n", VERSION, DATE);
    std::process::exit(0);
}

// ---------------------------------------------------------------------------
// Helper: print brief usage and exit
// ---------------------------------------------------------------------------
fn usage(msg: &str) -> ! {
    eprintln!("\n{}", msg);
    eprint!("\nUsage:  prodigal [-a trans_file] [-c] [-d nuc_file]");
    eprintln!(" [-f output_type]");
    eprint!("                 [-g tr_table] [-h] [-i input_file] [-m]");
    eprintln!(" [-n] [-o output_file]");
    eprint!("                 [-p mode] [-q] [-s start_file]");
    eprintln!(" [-t training_file] [-v]");
    eprintln!("\nDo 'prodigal -h' for more information.\n");
    std::process::exit(15);
}

// ---------------------------------------------------------------------------
// Helper: print full help and exit
// ---------------------------------------------------------------------------
fn help() -> ! {
    eprint!("\nUsage:  prodigal [-a trans_file] [-c] [-d nuc_file]");
    eprintln!(" [-f output_type]");
    eprint!("                 [-g tr_table] [-h] [-i input_file] [-m]");
    eprintln!(" [-n] [-o output_file]");
    eprint!("                 [-p mode] [-q] [-s start_file]");
    eprintln!(" [-t training_file] [-v]");
    eprint!("\n         -a:  Write protein translations to the selected ");
    eprintln!("file.");
    eprint!("         -c:  Closed ends.  Do not allow genes to run off ");
    eprintln!("edges.");
    eprint!("         -d:  Write nucleotide sequences of genes to the ");
    eprintln!("selected file.");
    eprint!("         -f:  Select output format (gbk, gff, or sco).  ");
    eprintln!("Default is gbk.");
    eprint!("         -g:  Specify a translation table to use (default");
    eprintln!(" 11).");
    eprintln!("         -h:  Print help menu and exit.");
    eprint!("         -i:  Specify FASTA/Genbank input file (default ");
    eprintln!("reads from stdin).");
    eprint!("         -m:  Treat runs of N as masked sequence; don't");
    eprintln!(" build genes across them.");
    eprint!("         -n:  Bypass Shine-Dalgarno trainer and force");
    eprintln!(" a full motif scan.");
    eprint!("         -o:  Specify output file (default writes to ");
    eprintln!("stdout).");
    eprint!("         -p:  Select procedure (single or meta).  Default");
    eprintln!(" is single.");
    eprintln!("         -q:  Run quietly (suppress normal stderr output).");
    eprint!("         -s:  Write all potential genes (with scores) to");
    eprintln!(" the selected file.");
    eprint!("         -t:  Write a training file (if none exists); ");
    eprintln!("otherwise, read and use");
    eprintln!("              the specified training file.");
    eprintln!("         -v:  Print version number and exit.\n");
    std::process::exit(0);
}

// ---------------------------------------------------------------------------
// Helper: copy stdin to a temp file (for piped input that needs rewinding)
// ---------------------------------------------------------------------------
fn copy_standard_input_to_file(path: &std::ffi::CStr, quiet: c_int) -> c_int {
    use std::io::{BufRead, Write};

    if quiet == 0 {
        eprint!("Piped input detected, copying stdin to a tmp file...");
    }

    let path_str = match path.to_str() {
        Ok(s) => s,
        Err(_) => return -1,
    };
    let file = match std::fs::File::create(path_str) {
        Ok(f) => f,
        Err(_) => return -1,
    };
    let mut writer = std::io::BufWriter::new(file);
    let stdin = std::io::stdin();
    let reader = stdin.lock();
    for line in reader.lines() {
        match line {
            Ok(l) => {
                if writeln!(writer, "{}", l).is_err() {
                    return -1;
                }
            }
            Err(_) => break,
        }
    }

    if quiet == 0 {
        eprintln!("done!");
        eprintln!("-------------------------------------");
    }
    0
}

// ---------------------------------------------------------------------------
// Main pipeline: replaces C main()
// ---------------------------------------------------------------------------
#[allow(unused_assignments)]
pub unsafe fn run_pipeline(args: &[String]) -> i32 {
    let argc = args.len() as c_int;

    // Convert args to C strings for compatibility (still needed for extern "C" calls)
    let c_args: Vec<std::ffi::CString> = args
        .iter()
        .map(|s| std::ffi::CString::new(s.as_str()).unwrap())
        .collect();
    let argv: Vec<*const c_char> = c_args.iter().map(|s| s.as_ptr()).collect();

    // Variable declarations
    let mut rv: c_int;
    let mut slen: c_int;
    let mut nn: c_int;
    let mut ng: c_int;
    let mut ipath: c_int;
    let gc_frame: *mut c_int;
    let mut do_training: c_int;
    let mut output: c_int;
    let mut max_phase: c_int;
    let mut closed: c_int;
    let mut do_mask: c_int;
    let mut nmask: c_int;
    let mut force_nonsd: c_int;
    let mut user_tt: c_int;
    let mut is_meta: c_int;
    let mut num_seq: c_int;
    let mut quiet: c_int;
    let mut piped: c_int;
    let mut max_slen: c_int;
    let mut max_score: f64;
    let mut gc: f64 = 0.0;
    let mut low: f64;
    let mut high: f64;

    // We still need raw C pointers for the extern "C" file path args
    let mut train_file: *const c_char = std::ptr::null();
    let mut start_file: *const c_char = std::ptr::null();
    let mut trans_file: *const c_char = std::ptr::null();
    let mut nuc_file: *const c_char = std::ptr::null();
    let mut input_file: *const c_char = std::ptr::null();
    let mut output_file: *const c_char = std::ptr::null();

    let input_copy_string = format!("tmp.prodigal.stdin.{}", std::process::id());
    let input_copy_cstr = std::ffi::CString::new(input_copy_string.as_str()).unwrap();

    let mut cur_header: [c_char; MAX_LINE] = [0; MAX_LINE];
    let mut new_header: [c_char; MAX_LINE] = [0; MAX_LINE];
    let mut short_header: [c_char; MAX_LINE] = [0; MAX_LINE];

    let mut output_ptr: c_int;
    let mut start_ptr: c_int;
    let mut trans_ptr: c_int;
    let mut nuc_ptr: c_int;
    let mut input_ptr: *mut c_void = std::ptr::null_mut();

    // Allocate memory using Vec
    let mut seq_vec: Vec<u8> = vec![0u8; MAX_SEQ / 4];
    let mut rseq_vec: Vec<u8> = vec![0u8; MAX_SEQ / 4];
    let mut useq_vec: Vec<u8> = vec![0u8; MAX_SEQ / 8];
    let mut nodes_vec: Vec<Node> = vec![unsafe { std::mem::zeroed() }; STT_NOD];
    let mut genes_vec: Vec<Gene> = vec![unsafe { std::mem::zeroed() }; MAX_GENES];

    let seq = seq_vec.as_mut_ptr();
    let rseq = rseq_vec.as_mut_ptr();
    let useq = useq_vec.as_mut_ptr();
    let mut nodes = nodes_vec.as_mut_ptr();
    let genes = genes_vec.as_mut_ptr();

    let mut tinf: Training = std::mem::zeroed();

    let mut meta: [MetagenomicBin; NUM_META] = std::mem::zeroed();
    let mut mlist: [Mask; MAX_MASKS] = [Mask { begin: 0, end: 0 }; MAX_MASKS];

    // Allocate Training structs for metagenomic bins using Box::into_raw
    // We collect the raw pointers so we can free them at the end.
    let mut meta_tinf_ptrs: Vec<*mut Training> = Vec::with_capacity(NUM_META);
    for i in 0..NUM_META {
        std::ptr::copy_nonoverlapping(
            b"None\0".as_ptr(),
            meta[i].desc.as_mut_ptr() as *mut u8,
            5,
        );
        let ptr = Box::into_raw(Box::new(std::mem::zeroed::<Training>()));
        meta[i].tinf = ptr;
        meta_tinf_ptrs.push(ptr);
    }

    // Initialize variables
    nn = 0;
    slen = 0;
    ipath = 0;
    ng = 0;
    nmask = 0;
    user_tt = 0;
    is_meta = 0;
    num_seq = 0;
    quiet = 0;
    max_phase = 0;
    max_score = -100.0;
    do_training = 0;
    piped = 0;
    max_slen = 0;
    output = 0;
    closed = 0;
    do_mask = 0;
    force_nonsd = 0;

    // Set up file descriptors (stdout = fd 1)
    let stdout_fd: c_int = 1;
    start_ptr = stdout_fd;
    trans_ptr = stdout_fd;
    nuc_ptr = stdout_fd;
    output_ptr = stdout_fd;

    // Set default training parameters
    tinf.st_wt = 4.35;
    tinf.trans_table = 11;

    // Parse command line arguments using Rust string matching
    let mut i: usize = 1;
    while i < argc as usize {
        let arg = args[i].as_str();

        // Check if this is a flag that requires a parameter but is the last arg
        if i == (argc - 1) as usize {
            match arg.to_lowercase().as_str() {
                "-t" | "-a" | "-g" | "-f" | "-s" | "-i" | "-o" | "-p" | "-d" => {
                    usage("-a/-f/-g/-i/-o/-p/-s options require parameters.");
                }
                _ => {}
            }
        }

        match arg.to_lowercase().as_str() {
            "-c" => {
                closed = 1;
            }
            "-q" => {
                quiet = 1;
            }
            "-m" => {
                do_mask = 1;
            }
            "-n" => {
                force_nonsd = 1;
            }
            "-h" => {
                help();
            }
            "-v" => {
                version();
            }
            "-a" => {
                trans_file = argv[i + 1];
                i += 1;
            }
            "-d" => {
                nuc_file = argv[i + 1];
                i += 1;
            }
            "-i" => {
                input_file = argv[i + 1];
                i += 1;
            }
            "-o" => {
                output_file = argv[i + 1];
                i += 1;
            }
            "-s" => {
                start_file = argv[i + 1];
                i += 1;
            }
            "-t" => {
                train_file = argv[i + 1];
                i += 1;
            }
            "-g" => {
                let val_str = args[i + 1].as_str();
                match val_str.parse::<c_int>() {
                    Ok(tt) => {
                        if tt < 1
                            || tt > 25
                            || tt == 7
                            || tt == 8
                            || (tt >= 17 && tt <= 20)
                        {
                            usage("Invalid translation table specified.");
                        }
                        tinf.trans_table = tt;
                        user_tt = tt;
                    }
                    Err(_) => {
                        usage("Invalid translation table specified.");
                    }
                }
                i += 1;
            }
            "-p" => {
                let val = args[i + 1].as_str();
                let first = val.as_bytes().first().copied().unwrap_or(0);
                if first == b'0' || first == b's' || first == b'S' {
                    is_meta = 0;
                } else if first == b'1' || first == b'm' || first == b'M' {
                    is_meta = 1;
                } else {
                    usage("Invalid meta/single genome type specified.");
                }
                i += 1;
            }
            "-f" => {
                let farg = args[i + 1].as_str();
                match farg.to_lowercase().as_str() {
                    "0" | "gbk" => output = 0,
                    "1" | "gca" => output = 1,
                    "2" | "sco" => output = 2,
                    "3" | "gff" => output = 3,
                    _ => {
                        // Also check if starts with 0/1/2/3
                        if farg.starts_with('0') {
                            output = 0;
                        } else if farg.starts_with('1') {
                            output = 1;
                        } else if farg.starts_with('2') {
                            output = 2;
                        } else if farg.starts_with('3') {
                            output = 3;
                        } else {
                            usage("Invalid output format specified.");
                        }
                    }
                }
                i += 1;
            }
            _ => {
                usage("Unknown option.");
            }
        }

        i += 1;
    }

    // Print header
    if quiet == 0 {
        eprintln!("-------------------------------------");
        eprintln!("PRODIGAL v{} [{}]         ", VERSION, DATE);
        eprintln!("Univ of Tenn / Oak Ridge National Lab");
        eprintln!("Doug Hyatt, Loren Hauser, et al.     ");
        eprintln!("-------------------------------------");
    }

    // Read in the training file (if specified)
    if !train_file.is_null() {
        if is_meta == 1 {
            eprint!("\nError: cannot specify metagenomic sequence with a");
            eprintln!(" training file.");
            return 2;
        }
        rv = read_training_file(train_file, &mut tinf);
        if rv == 1 {
            do_training = 1;
        } else {
            if force_nonsd == 1 {
                eprint!("\nError: cannot force non-SD finder with a training");
                eprintln!(" file already created!");
                return 3;
            }
            if quiet == 0 {
                let tf_str = CStr::from_ptr(train_file).to_string_lossy();
                eprint!("Reading in training data from file {}...", tf_str);
            }
            if user_tt > 0 && user_tt != tinf.trans_table {
                eprint!("\n\nWarning: user-specified translation table does");
                eprintln!("not match the one in the specified training file! \n");
            }
            if rv == -1 {
                eprintln!("\n\nError: training file did not read correctly!");
                return 4;
            }
            if quiet == 0 {
                eprintln!("done!");
                eprintln!("-------------------------------------");
            }
        }
    }

    // Determine where standard input is coming from
    if is_meta == 0 && train_file.is_null() && input_file.is_null() {
        // Use std::fs::metadata on /dev/stdin to check file type
        use std::os::unix::fs::FileTypeExt;
        match std::fs::metadata("/dev/stdin") {
            Err(_) => {
                eprintln!("\nError: can't fstat standard input.\n");
                return 5;
            }
            Ok(md) => {
                let ft = md.file_type();
                if ft.is_char_device() {
                    help();
                } else if ft.is_file() {
                    // do nothing
                } else if ft.is_fifo() {
                    piped = 1;
                    if copy_standard_input_to_file(&input_copy_cstr, quiet) == -1 {
                        eprintln!("\nError: can't copy stdin to file.\n");
                        return 5;
                    }
                    input_file = input_copy_cstr.as_ptr();
                }
            }
        }
    }

    // Check i/o files and prepare them for reading/writing
    if !input_file.is_null() {
        input_ptr = seq_reader_open(input_file);
        if input_ptr.is_null() {
            let f = CStr::from_ptr(input_file).to_string_lossy();
            eprintln!("\nError: can't open input file {}.\n", f);
            return 5;
        }
    }
    if input_ptr.is_null() {
        input_ptr = seq_reader_open(
            b"/dev/stdin\0".as_ptr() as *const c_char,
        );
        if input_ptr.is_null() {
            eprintln!("\nError: can't open stdin.\n");
            return 5;
        }
    }
    if !output_file.is_null() {
        let path = CStr::from_ptr(output_file).to_string_lossy();
        match std::fs::File::create(path.as_ref()) {
            Ok(f) => { output_ptr = f.into_raw_fd(); }
            Err(_) => { eprintln!("\nError: can't open output file {}.\n", path); return 6; }
        }
    }
    if !start_file.is_null() {
        let path = CStr::from_ptr(start_file).to_string_lossy();
        match std::fs::File::create(path.as_ref()) {
            Ok(f) => { start_ptr = f.into_raw_fd(); }
            Err(_) => { eprintln!("\nError: can't open start file {}.\n", path); return 7; }
        }
    }
    if !trans_file.is_null() {
        let path = CStr::from_ptr(trans_file).to_string_lossy();
        match std::fs::File::create(path.as_ref()) {
            Ok(f) => { trans_ptr = f.into_raw_fd(); }
            Err(_) => { eprintln!("\nError: can't open translation file {}.\n", path); return 8; }
        }
    }
    if !nuc_file.is_null() {
        let path = CStr::from_ptr(nuc_file).to_string_lossy();
        match std::fs::File::create(path.as_ref()) {
            Ok(f) => { nuc_ptr = f.into_raw_fd(); }
            Err(_) => { eprintln!("\nError: can't open gene nucleotide file {}.\n", path); return 16; }
        }
    }

    // =========================================================================
    // Single Genome Training
    // =========================================================================
    if is_meta == 0 && (do_training == 1 || (do_training == 0 && train_file.is_null())) {
        if quiet == 0 {
            eprintln!("Request:  Single Genome, Phase:  Training");
            eprint!("Reading in the sequence(s) to train...");
        }
        slen = read_seq_training(
            input_ptr,
            seq,
            useq,
            &mut tinf.gc,
            do_mask,
            mlist.as_mut_ptr(),
            &mut nmask,
        );
        if slen == 0 {
            eprint!("\n\nSequence read failed (file must be Fasta, ");
            eprintln!("Genbank, or EMBL format).\n");
            return 9;
        }
        if slen < MIN_SINGLE_GENOME {
            eprint!("\n\nError:  Sequence must be {}", MIN_SINGLE_GENOME);
            eprint!(" characters (only {} read).\n(Consider", slen);
            eprint!(" running with the -p meta option or finding");
            eprintln!(" more contigs from the same genome.)\n");
            return 10;
        }
        if slen < IDEAL_SINGLE_GENOME {
            eprint!("\n\nWarning:  ideally Prodigal should be given at");
            eprint!(" least {} bases for ", IDEAL_SINGLE_GENOME);
            eprint!("training.\nYou may get better results with the ");
            eprintln!("-p meta option.\n");
        }
        rcom_seq(seq, rseq, useq, slen);
        if quiet == 0 {
            eprintln!("{} bp seq created, {:.2} pct GC", slen, tinf.gc * 100.0);
        }

        // Find all potential starts and stops
        if quiet == 0 {
            eprint!("Locating all potential starts and stops...");
        }
        if slen > max_slen && slen > (STT_NOD as c_int) * 8 {
            let new_size = slen as usize / 8;
            nodes_vec.resize(new_size, std::mem::zeroed());
            nodes = nodes_vec.as_mut_ptr();
            max_slen = slen;
        }
        nn = add_nodes(
            seq,
            rseq,
            slen,
            nodes,
            closed,
            mlist.as_mut_ptr(),
            nmask,
            &mut tinf,
        );
        sort_nodes(&mut nodes_vec[..nn as usize]);
        if quiet == 0 {
            eprintln!("{} nodes", nn);
        }

        // Scan ORFs for GC bias
        if quiet == 0 {
            eprint!("Looking for GC bias in different frames...");
        }
        gc_frame = calc_most_gc_frame(seq, slen);
        if gc_frame.is_null() {
            eprintln!("Malloc failed on gc frame plot\n");
            return 11;
        }
        record_gc_bias(gc_frame, nodes, nn, &mut tinf);
        if quiet == 0 {
            eprintln!(
                "frame bias scores: {:.2} {:.2} {:.2}",
                tinf.bias[0], tinf.bias[1], tinf.bias[2]
            );
        }
        drop(Vec::from_raw_parts(gc_frame, slen as usize, slen as usize));

        // Initial DP with GC frame bias
        if quiet == 0 {
            eprint!("Building initial set of genes to train from...");
        }
        record_overlapping_starts(nodes, nn, &mut tinf, 0);
        ipath = dprog(nodes, nn, &mut tinf, 0);
        if quiet == 0 {
            eprintln!("done!");
        }

        // Gather dicodon statistics
        if quiet == 0 {
            eprint!("Creating coding model and scoring nodes...");
        }
        calc_dicodon_gene(&mut tinf, seq, rseq, slen, nodes, ipath);
        raw_coding_score(seq, rseq, slen, nodes, nn, &mut tinf);
        if quiet == 0 {
            eprintln!("done!");
        }

        // Determine SD usage and train starts
        if quiet == 0 {
            eprint!("Examining upstream regions and training starts...");
        }
        rbs_score(seq, rseq, slen, nodes, nn, &mut tinf);
        train_starts_sd(seq, rseq, slen, nodes, nn, &mut tinf);
        determine_sd_usage(&mut tinf);
        if force_nonsd == 1 {
            tinf.uses_sd = 0;
        }
        if tinf.uses_sd == 0 {
            train_starts_nonsd(seq, rseq, slen, nodes, nn, &mut tinf);
        }
        if quiet == 0 {
            eprintln!("done!");
        }

        // If training specified, write the training file and exit
        if do_training == 1 {
            if quiet == 0 {
                let tf_str = CStr::from_ptr(train_file).to_string_lossy();
                eprint!("Writing data to training file {}...", tf_str);
            }
            rv = write_training_file(train_file, &mut tinf);
            if rv != 0 {
                eprintln!("\nError: could not write training file!");
                return 12;
            } else {
                if quiet == 0 {
                    eprintln!("done!");
                }
                // Clean up and return
                // Vec memory freed automatically when dropped
                for ptr in &meta_tinf_ptrs {
                    drop(Box::from_raw(*ptr));
                }
                seq_reader_close(input_ptr);
                if output_ptr != stdout_fd {
                    drop(std::fs::File::from_raw_fd(output_ptr));
                }
                if start_ptr != stdout_fd {
                    drop(std::fs::File::from_raw_fd(start_ptr));
                }
                if trans_ptr != stdout_fd {
                    drop(std::fs::File::from_raw_fd(trans_ptr));
                }
                return 0;
            }
        }

        // Rewind input file
        if quiet == 0 {
            eprintln!("-------------------------------------");
        }
        if seq_reader_seek(input_ptr, 0, 0) == -1 {  // SEEK_SET = 0
            eprintln!("\nError: could not rewind input file.");
            return 13;
        }

        // Reset sequence/DP variables
        std::ptr::write_bytes(seq, 0, (slen as usize / 4 + 1).min(seq_vec.len()));
        std::ptr::write_bytes(rseq, 0, (slen as usize / 4 + 1).min(rseq_vec.len()));
        std::ptr::write_bytes(useq, 0, (slen as usize / 8 + 1).min(useq_vec.len()));
        std::ptr::write_bytes(nodes, 0, nn as usize);
        nn = 0;
        slen = 0;
        ipath = 0;
        nmask = 0;
    }
    // Initialize metagenomic bins
    else if is_meta == 1 {
        if quiet == 0 {
            eprintln!("Request:  Metagenomic, Phase:  Training");
            eprint!("Initializing training files...");
        }
        initialize_metagenomic_bins(meta.as_mut_ptr());
        if quiet == 0 {
            eprintln!("done!");
            eprintln!("-------------------------------------");
        }
    }

    // Print header for gene finding phase
    if quiet == 0 {
        if is_meta == 1 {
            eprintln!("Request:  Metagenomic, Phase:  Gene Finding");
        } else {
            eprintln!("Request:  Single Genome, Phase:  Gene Finding");
        }
    }

    // Read and process each sequence
    {
        let s = std::ffi::CString::new("Prodigal_Seq_1").unwrap();
        let bytes = s.as_bytes_with_nul();
        std::ptr::copy_nonoverlapping(bytes.as_ptr() as *const c_char, cur_header.as_mut_ptr(), bytes.len());
    }
    {
        let s = std::ffi::CString::new("Prodigal_Seq_2").unwrap();
        let bytes = s.as_bytes_with_nul();
        std::ptr::copy_nonoverlapping(bytes.as_ptr() as *const c_char, new_header.as_mut_ptr(), bytes.len());
    }

    loop {
        slen = next_seq_multi(
            input_ptr,
            seq,
            useq,
            &mut num_seq,
            &mut gc,
            do_mask,
            mlist.as_mut_ptr(),
            &mut nmask,
            cur_header.as_mut_ptr(),
            new_header.as_mut_ptr(),
        );
        if slen == -1 {
            break;
        }

        rcom_seq(seq, rseq, useq, slen);
        if slen == 0 {
            eprint!("\nSequence read failed (file must be Fasta, ");
            eprintln!("Genbank, or EMBL format).\n");
            return 14;
        }

        if quiet == 0 {
            eprint!("Finding genes in sequence #{} ({} bp)...", num_seq, slen);
        }

        // Reallocate if this is the biggest sequence we've seen
        if slen > max_slen && slen > (STT_NOD as c_int) * 8 {
            let new_size = slen as usize / 8;
            nodes_vec.resize(new_size, std::mem::zeroed());
            nodes = nodes_vec.as_mut_ptr();
            max_slen = slen;
        }

        // Calculate short header
        calc_short_header(cur_header.as_mut_ptr(), short_header.as_mut_ptr(), num_seq);

        if is_meta == 0 {
            // ---- Single Genome Version ----
            nn = add_nodes(
                seq,
                rseq,
                slen,
                nodes,
                closed,
                mlist.as_mut_ptr(),
                nmask,
                &mut tinf,
            );
            sort_nodes(&mut nodes_vec[..nn as usize]);

            score_nodes(seq, rseq, slen, nodes, nn, &mut tinf, closed, is_meta);
            if start_ptr != stdout_fd {
                write_start_file(
                    start_ptr,
                    nodes,
                    nn,
                    &mut tinf,
                    num_seq,
                    slen,
                    0,
                    std::ptr::null_mut(),
                    VERSION_CSTR.as_ptr() as *mut c_char,
                    cur_header.as_mut_ptr(),
                );
            }
            record_overlapping_starts(nodes, nn, &mut tinf, 1);
            ipath = dprog(nodes, nn, &mut tinf, 1);
            eliminate_bad_genes(nodes, ipath, &mut tinf);
            ng = add_genes(genes, nodes, ipath);
            tweak_final_starts(genes, ng, nodes, nn, &mut tinf);
            record_gene_data(genes, ng, nodes, &mut tinf, num_seq);
            if quiet == 0 {
                eprintln!("done!");
            }

            // Output the genes
            print_genes(
                output_ptr,
                genes,
                ng,
                nodes,
                slen,
                output,
                num_seq,
                0,
                std::ptr::null_mut(),
                &mut tinf,
                cur_header.as_mut_ptr(),
                short_header.as_mut_ptr(),
                VERSION_CSTR.as_ptr() as *mut c_char,
            );
            // fd writes go directly to kernel, no fflush needed
            if trans_ptr != stdout_fd {
                write_translations(
                    trans_ptr,
                    genes,
                    ng,
                    nodes,
                    seq,
                    rseq,
                    useq,
                    slen,
                    &mut tinf,
                    num_seq,
                    short_header.as_mut_ptr(),
                );
            }
            if nuc_ptr != stdout_fd {
                write_nucleotide_seqs(
                    nuc_ptr,
                    genes,
                    ng,
                    nodes,
                    seq,
                    rseq,
                    useq,
                    slen,
                    &mut tinf,
                    num_seq,
                    short_header.as_mut_ptr(),
                );
            }
        } else {
            // ---- Metagenomic Version ----
            low = 0.88495 * gc - 0.0102337;
            if low > 0.65 {
                low = 0.65;
            }
            high = 0.86596 * gc + 0.1131991;
            if high < 0.35 {
                high = 0.35;
            }

            max_score = -100.0;
            for mi in 0..NUM_META as c_int {
                if mi == 0
                    || (*meta[mi as usize].tinf).trans_table
                        != (*meta[(mi - 1) as usize].tinf).trans_table
                {
                    std::ptr::write_bytes(nodes, 0, nn as usize);
                    nn = add_nodes(
                        seq,
                        rseq,
                        slen,
                        nodes,
                        closed,
                        mlist.as_mut_ptr(),
                        nmask,
                        meta[mi as usize].tinf,
                    );
                    sort_nodes(&mut nodes_vec[..nn as usize]);
                }
                if (*meta[mi as usize].tinf).gc < low
                    || (*meta[mi as usize].tinf).gc > high
                {
                    continue;
                }
                reset_node_scores(nodes, nn);
                score_nodes(
                    seq,
                    rseq,
                    slen,
                    nodes,
                    nn,
                    meta[mi as usize].tinf,
                    closed,
                    is_meta,
                );
                record_overlapping_starts(nodes, nn, meta[mi as usize].tinf, 1);
                ipath = dprog(nodes, nn, meta[mi as usize].tinf, 1);
                if (*nodes.offset(ipath as isize)).score > max_score {
                    max_phase = mi;
                    max_score = (*nodes.offset(ipath as isize)).score;
                    eliminate_bad_genes(nodes, ipath, meta[mi as usize].tinf);
                    ng = add_genes(genes, nodes, ipath);
                    tweak_final_starts(genes, ng, nodes, nn, meta[mi as usize].tinf);
                    record_gene_data(genes, ng, nodes, meta[mi as usize].tinf, num_seq);
                }
            }

            // Recover the nodes for the best run
            std::ptr::write_bytes(nodes, 0, nn as usize);
            nn = add_nodes(
                seq,
                rseq,
                slen,
                nodes,
                closed,
                mlist.as_mut_ptr(),
                nmask,
                meta[max_phase as usize].tinf,
            );
            sort_nodes(&mut nodes_vec[..nn as usize]);
            score_nodes(
                seq,
                rseq,
                slen,
                nodes,
                nn,
                meta[max_phase as usize].tinf,
                closed,
                is_meta,
            );
            if start_ptr != stdout_fd {
                write_start_file(
                    start_ptr,
                    nodes,
                    nn,
                    meta[max_phase as usize].tinf,
                    num_seq,
                    slen,
                    1,
                    meta[max_phase as usize].desc.as_mut_ptr(),
                    VERSION_CSTR.as_ptr() as *mut c_char,
                    cur_header.as_mut_ptr(),
                );
            }

            if quiet == 0 {
                eprintln!("done!");
            }

            // Output the genes
            print_genes(
                output_ptr,
                genes,
                ng,
                nodes,
                slen,
                output,
                num_seq,
                1,
                meta[max_phase as usize].desc.as_mut_ptr(),
                meta[max_phase as usize].tinf,
                cur_header.as_mut_ptr(),
                short_header.as_mut_ptr(),
                VERSION_CSTR.as_ptr() as *mut c_char,
            );
            // fd writes go directly to kernel, no fflush needed
            if trans_ptr != stdout_fd {
                write_translations(
                    trans_ptr,
                    genes,
                    ng,
                    nodes,
                    seq,
                    rseq,
                    useq,
                    slen,
                    meta[max_phase as usize].tinf,
                    num_seq,
                    short_header.as_mut_ptr(),
                );
            }
            if nuc_ptr != stdout_fd {
                write_nucleotide_seqs(
                    nuc_ptr,
                    genes,
                    ng,
                    nodes,
                    seq,
                    rseq,
                    useq,
                    slen,
                    meta[max_phase as usize].tinf,
                    num_seq,
                    short_header.as_mut_ptr(),
                );
            }
        }

        // Reset sequence/DP variables
        std::ptr::write_bytes(seq, 0, (slen as usize / 4 + 1).min(seq_vec.len()));
        std::ptr::write_bytes(rseq, 0, (slen as usize / 4 + 1).min(rseq_vec.len()));
        std::ptr::write_bytes(useq, 0, (slen as usize / 8 + 1).min(useq_vec.len()));
        std::ptr::write_bytes(nodes, 0, nn as usize);
        nn = 0;
        slen = 0;
        ipath = 0;
        nmask = 0;
        std::ptr::copy_nonoverlapping(
            new_header.as_ptr(),
            cur_header.as_mut_ptr(),
            MAX_LINE,
        );
        {
            let s = format!("Prodigal_Seq_{}\n", num_seq + 1);
            let c = std::ffi::CString::new(s).unwrap();
            let bytes = c.as_bytes_with_nul();
            std::ptr::copy_nonoverlapping(
                bytes.as_ptr() as *const c_char,
                new_header.as_mut_ptr(),
                bytes.len().min(MAX_LINE),
            );
        }
    }

    if num_seq == 0 {
        eprintln!("\nError:  no input sequences to analyze.\n");
        return 18;
    }

    // Free metagenomic training data (allocated with Box)
    for ptr in &meta_tinf_ptrs {
        drop(Box::from_raw(*ptr));
    }

    // Close all filehandles
    seq_reader_close(input_ptr);
    if output_ptr != stdout_fd {
        drop(std::fs::File::from_raw_fd(output_ptr));
    }
    if start_ptr != stdout_fd {
        drop(std::fs::File::from_raw_fd(start_ptr));
    }
    if trans_ptr != stdout_fd {
        drop(std::fs::File::from_raw_fd(trans_ptr));
    }
    if nuc_ptr != stdout_fd {
        drop(std::fs::File::from_raw_fd(nuc_ptr));
    }

    // Remove tmp file
    if piped == 1 {
        let path_str = input_copy_cstr.to_str().unwrap_or("");
        if std::fs::remove_file(path_str).is_err() {
            eprintln!("Could not delete tmp file {}.", path_str);
            return 18;
        }
    }

    0
}
