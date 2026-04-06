/*******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Copyright (C) 2007-2016 University of Tennessee / UT-Battelle

    Code Author:  Doug Hyatt

    Rust port of main.c — the pipeline orchestration module.
*******************************************************************************/

use std::os::raw::{c_char, c_int, c_void};

use crate::types::{Gene, Mask, MetagenomicBin, Node, Training};
use crate::types::{MAX_GENES, MAX_LINE, MAX_MASKS, MAX_SEQ, NUM_META, STT_NOD};

const VERSION: &[u8] = b"2.6.3\0";
const DATE: &[u8] = b"February, 2016\0";
const MIN_SINGLE_GENOME: c_int = 20000;
const IDEAL_SINGLE_GENOME: c_int = 100000;

// ---------------------------------------------------------------------------
// External C functions from zlib
// ---------------------------------------------------------------------------
extern "C" {
    fn gzopen(path: *const c_char, mode: *const c_char) -> *mut c_void;
    fn gzclose(file: *mut c_void) -> c_int;
    fn gzseek(file: *mut c_void, offset: libc::c_long, whence: c_int) -> libc::c_long;
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

    fn compare_nodes(v1: *const c_void, v2: *const c_void) -> c_int;

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
        fh: *mut libc::FILE,
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
        fp: *mut libc::FILE,
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
        fh: *mut libc::FILE,
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
        fh: *mut libc::FILE,
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
// Helper: print version and exit
// ---------------------------------------------------------------------------
unsafe fn version() {
    libc::fprintf(
        libc::fdopen(2, b"w\0".as_ptr() as *const c_char),
        b"\nProdigal V%s: %s\n\n\0".as_ptr() as *const c_char,
        VERSION.as_ptr() as *const c_char,
        DATE.as_ptr() as *const c_char,
    );
    libc::exit(0);
}

// ---------------------------------------------------------------------------
// Helper: print brief usage and exit
// ---------------------------------------------------------------------------
unsafe fn usage(msg: *const c_char) {
    let stderr = libc::fdopen(2, b"w\0".as_ptr() as *const c_char);
    libc::fprintf(stderr, b"\n%s\n\0".as_ptr() as *const c_char, msg);
    libc::fprintf(
        stderr,
        b"\nUsage:  prodigal [-a trans_file] [-c] [-d nuc_file]\0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b" [-f output_type]\n\0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b"                 [-g tr_table] [-h] [-i input_file] [-m]\0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b" [-n] [-o output_file]\n\0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b"                 [-p mode] [-q] [-s start_file]\0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b" [-t training_file] [-v]\n\0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b"\nDo 'prodigal -h' for more information.\n\n\0".as_ptr() as *const c_char,
    );
    libc::exit(15);
}

// ---------------------------------------------------------------------------
// Helper: print full help and exit
// ---------------------------------------------------------------------------
unsafe fn help() {
    let stderr = libc::fdopen(2, b"w\0".as_ptr() as *const c_char);
    libc::fprintf(
        stderr,
        b"\nUsage:  prodigal [-a trans_file] [-c] [-d nuc_file]\0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b" [-f output_type]\n\0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b"                 [-g tr_table] [-h] [-i input_file] [-m]\0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b" [-n] [-o output_file]\n\0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b"                 [-p mode] [-q] [-s start_file]\0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b" [-t training_file] [-v]\n\0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b"\n         -a:  Write protein translations to the selected \0".as_ptr() as *const c_char,
    );
    libc::fprintf(stderr, b"file.\n\0".as_ptr() as *const c_char);
    libc::fprintf(
        stderr,
        b"         -c:  Closed ends.  Do not allow genes to run off \0".as_ptr() as *const c_char,
    );
    libc::fprintf(stderr, b"edges.\n\0".as_ptr() as *const c_char);
    libc::fprintf(
        stderr,
        b"         -d:  Write nucleotide sequences of genes to the \0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b"selected file.\n\0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b"         -f:  Select output format (gbk, gff, or sco).  \0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b"Default is gbk.\n\0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b"         -g:  Specify a translation table to use (default\0".as_ptr() as *const c_char,
    );
    libc::fprintf(stderr, b" 11).\n\0".as_ptr() as *const c_char);
    libc::fprintf(
        stderr,
        b"         -h:  Print help menu and exit.\n\0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b"         -i:  Specify FASTA/Genbank input file (default \0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b"reads from stdin).\n\0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b"         -m:  Treat runs of N as masked sequence; don't\0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b" build genes across them.\n\0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b"         -n:  Bypass Shine-Dalgarno trainer and force\0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b" a full motif scan.\n\0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b"         -o:  Specify output file (default writes to \0".as_ptr() as *const c_char,
    );
    libc::fprintf(stderr, b"stdout).\n\0".as_ptr() as *const c_char);
    libc::fprintf(
        stderr,
        b"         -p:  Select procedure (single or meta).  Default\0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b" is single.\n\0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b"         -q:  Run quietly (suppress normal stderr output).\n\0".as_ptr()
            as *const c_char,
    );
    libc::fprintf(
        stderr,
        b"         -s:  Write all potential genes (with scores) to\0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b" the selected file.\n\0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b"         -t:  Write a training file (if none exists); \0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b"otherwise, read and use\n\0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b"              the specified training file.\n\0".as_ptr() as *const c_char,
    );
    libc::fprintf(
        stderr,
        b"         -v:  Print version number and exit.\n\n\0".as_ptr() as *const c_char,
    );
    libc::exit(0);
}

// ---------------------------------------------------------------------------
// Helper: copy stdin to a temp file (for piped input that needs rewinding)
// ---------------------------------------------------------------------------
unsafe fn copy_standard_input_to_file(path: *const c_char, quiet: c_int) -> c_int {
    let mut line: [c_char; MAX_LINE + 1] = [0; MAX_LINE + 1];
    let stderr = libc::fdopen(2, b"w\0".as_ptr() as *const c_char);

    if quiet == 0 {
        libc::fprintf(
            stderr,
            b"Piped input detected, copying stdin to a tmp file...\0".as_ptr() as *const c_char,
        );
    }

    let wp = libc::fopen(path, b"w\0".as_ptr() as *const c_char);
    if wp.is_null() {
        return -1;
    }
    while !libc::fgets(
        line.as_mut_ptr(),
        MAX_LINE as c_int,
        libc::fdopen(0, b"r\0".as_ptr() as *const c_char),
    )
    .is_null()
    {
        libc::fprintf(wp, b"%s\0".as_ptr() as *const c_char, line.as_ptr());
    }
    libc::fclose(wp);

    if quiet == 0 {
        libc::fprintf(
            stderr,
            b"done!\n\0".as_ptr() as *const c_char,
        );
        libc::fprintf(
            stderr,
            b"-------------------------------------\n\0".as_ptr() as *const c_char,
        );
    }
    0
}

// ---------------------------------------------------------------------------
// Main pipeline: replaces C main()
// ---------------------------------------------------------------------------
pub unsafe fn run_pipeline(args: &[String]) -> i32 {
    let argc = args.len() as c_int;

    // Convert args to C strings for compatibility
    let c_args: Vec<std::ffi::CString> = args
        .iter()
        .map(|s| std::ffi::CString::new(s.as_str()).unwrap())
        .collect();
    let argv: Vec<*const c_char> = c_args.iter().map(|s| s.as_ptr()).collect();

    let stderr = libc::fdopen(2, b"w\0".as_ptr() as *const c_char);

    // Variable declarations
    let mut rv: c_int;
    let mut slen: c_int;
    let mut nn: c_int;
    let mut ng: c_int;
    let mut ipath: c_int;
    let mut gc_frame: *mut c_int;
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
    let mut fnum: c_int;
    let mut max_score: f64;
    let mut gc: f64 = 0.0;
    let mut low: f64;
    let mut high: f64;

    let mut train_file: *const c_char = std::ptr::null();
    let mut start_file: *const c_char = std::ptr::null();
    let mut trans_file: *const c_char = std::ptr::null();
    let mut nuc_file: *const c_char = std::ptr::null();
    let mut input_file: *const c_char = std::ptr::null();
    let mut output_file: *const c_char = std::ptr::null();

    let mut input_copy: [c_char; MAX_LINE] = [0; MAX_LINE];
    let mut cur_header: [c_char; MAX_LINE] = [0; MAX_LINE];
    let mut new_header: [c_char; MAX_LINE] = [0; MAX_LINE];
    let mut short_header: [c_char; MAX_LINE] = [0; MAX_LINE];

    let mut output_ptr: *mut libc::FILE;
    let mut start_ptr: *mut libc::FILE;
    let mut trans_ptr: *mut libc::FILE;
    let mut nuc_ptr: *mut libc::FILE;
    let mut input_ptr: *mut c_void = std::ptr::null_mut();

    let mut fbuf: libc::stat = std::mem::zeroed();

    // Allocate memory
    let seq = libc::malloc(MAX_SEQ / 4) as *mut u8;
    let rseq = libc::malloc(MAX_SEQ / 4) as *mut u8;
    let useq = libc::malloc(MAX_SEQ / 8) as *mut u8;
    let mut nodes = libc::malloc(STT_NOD * std::mem::size_of::<Node>()) as *mut Node;
    let genes = libc::malloc(MAX_GENES * std::mem::size_of::<Gene>()) as *mut Gene;

    if seq.is_null() || rseq.is_null() || nodes.is_null() || genes.is_null() {
        libc::fprintf(
            stderr,
            b"\nError: Malloc failed on sequence/orfs\n\n\0".as_ptr() as *const c_char,
        );
        return 1;
    }

    libc::memset(seq as *mut c_void, 0, MAX_SEQ / 4);
    libc::memset(rseq as *mut c_void, 0, MAX_SEQ / 4);
    libc::memset(useq as *mut c_void, 0, MAX_SEQ / 8);
    libc::memset(
        nodes as *mut c_void,
        0,
        STT_NOD * std::mem::size_of::<Node>(),
    );
    libc::memset(
        genes as *mut c_void,
        0,
        MAX_GENES * std::mem::size_of::<Gene>(),
    );

    let mut tinf: Training = std::mem::zeroed();

    let mut meta: [MetagenomicBin; NUM_META] = std::mem::zeroed();
    let mut mlist: [Mask; MAX_MASKS] = [Mask { begin: 0, end: 0 }; MAX_MASKS];

    for i in 0..NUM_META {
        libc::memset(
            &mut meta[i] as *mut MetagenomicBin as *mut c_void,
            0,
            std::mem::size_of::<MetagenomicBin>(),
        );
        libc::strcpy(
            meta[i].desc.as_mut_ptr(),
            b"None\0".as_ptr() as *const c_char,
        );
        meta[i].tinf = libc::malloc(std::mem::size_of::<Training>()) as *mut Training;
        if meta[i].tinf.is_null() {
            libc::fprintf(
                stderr,
                b"\nError: Malloc failed on training structure.\n\n\0".as_ptr() as *const c_char,
            );
            return 1;
        }
        libc::memset(
            meta[i].tinf as *mut c_void,
            0,
            std::mem::size_of::<Training>(),
        );
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

    let stdout = libc::fdopen(1, b"w\0".as_ptr() as *const c_char);
    start_ptr = stdout;
    trans_ptr = stdout;
    nuc_ptr = stdout;
    output_ptr = stdout;

    // Filename for input copy if needed
    let pid = libc::getpid();
    libc::sprintf(
        input_copy.as_mut_ptr(),
        b"tmp.prodigal.stdin.%d\0".as_ptr() as *const c_char,
        pid,
    );

    // Set default training parameters
    tinf.st_wt = 4.35;
    tinf.trans_table = 11;

    // Parse command line arguments
    let mut i: c_int = 1;
    while i < argc {
        let arg = argv[i as usize];

        // Check if this is a flag that requires a parameter but is the last arg
        if i == argc - 1
            && (libc::strcmp(arg, b"-t\0".as_ptr() as *const c_char) == 0
                || libc::strcmp(arg, b"-T\0".as_ptr() as *const c_char) == 0
                || libc::strcmp(arg, b"-a\0".as_ptr() as *const c_char) == 0
                || libc::strcmp(arg, b"-A\0".as_ptr() as *const c_char) == 0
                || libc::strcmp(arg, b"-g\0".as_ptr() as *const c_char) == 0
                || libc::strcmp(arg, b"-G\0".as_ptr() as *const c_char) == 0
                || libc::strcmp(arg, b"-f\0".as_ptr() as *const c_char) == 0
                || libc::strcmp(arg, b"-F\0".as_ptr() as *const c_char) == 0
                || libc::strcmp(arg, b"-s\0".as_ptr() as *const c_char) == 0
                || libc::strcmp(arg, b"-S\0".as_ptr() as *const c_char) == 0
                || libc::strcmp(arg, b"-i\0".as_ptr() as *const c_char) == 0
                || libc::strcmp(arg, b"-I\0".as_ptr() as *const c_char) == 0
                || libc::strcmp(arg, b"-o\0".as_ptr() as *const c_char) == 0
                || libc::strcmp(arg, b"-O\0".as_ptr() as *const c_char) == 0
                || libc::strcmp(arg, b"-p\0".as_ptr() as *const c_char) == 0
                || libc::strcmp(arg, b"-P\0".as_ptr() as *const c_char) == 0)
        {
            usage(
                b"-a/-f/-g/-i/-o/-p/-s options require parameters.\0".as_ptr() as *const c_char,
            );
        } else if libc::strcmp(arg, b"-c\0".as_ptr() as *const c_char) == 0
            || libc::strcmp(arg, b"-C\0".as_ptr() as *const c_char) == 0
        {
            closed = 1;
        } else if libc::strcmp(arg, b"-q\0".as_ptr() as *const c_char) == 0
            || libc::strcmp(arg, b"-Q\0".as_ptr() as *const c_char) == 0
        {
            quiet = 1;
        } else if libc::strcmp(arg, b"-m\0".as_ptr() as *const c_char) == 0
            || libc::strcmp(arg, b"-M\0".as_ptr() as *const c_char) == 0
        {
            do_mask = 1;
        } else if libc::strcmp(arg, b"-n\0".as_ptr() as *const c_char) == 0
            || libc::strcmp(arg, b"-N\0".as_ptr() as *const c_char) == 0
        {
            force_nonsd = 1;
        } else if libc::strcmp(arg, b"-h\0".as_ptr() as *const c_char) == 0
            || libc::strcmp(arg, b"-H\0".as_ptr() as *const c_char) == 0
        {
            help();
        } else if libc::strcmp(arg, b"-v\0".as_ptr() as *const c_char) == 0
            || libc::strcmp(arg, b"-V\0".as_ptr() as *const c_char) == 0
        {
            version();
        } else if libc::strcmp(arg, b"-a\0".as_ptr() as *const c_char) == 0
            || libc::strcmp(arg, b"-A\0".as_ptr() as *const c_char) == 0
        {
            trans_file = argv[(i + 1) as usize];
            i += 1;
        } else if libc::strcmp(arg, b"-d\0".as_ptr() as *const c_char) == 0
            || libc::strcmp(arg, b"-D\0".as_ptr() as *const c_char) == 0
        {
            nuc_file = argv[(i + 1) as usize];
            i += 1;
        } else if libc::strcmp(arg, b"-i\0".as_ptr() as *const c_char) == 0
            || libc::strcmp(arg, b"-I\0".as_ptr() as *const c_char) == 0
        {
            input_file = argv[(i + 1) as usize];
            i += 1;
        } else if libc::strcmp(arg, b"-o\0".as_ptr() as *const c_char) == 0
            || libc::strcmp(arg, b"-O\0".as_ptr() as *const c_char) == 0
        {
            output_file = argv[(i + 1) as usize];
            i += 1;
        } else if libc::strcmp(arg, b"-s\0".as_ptr() as *const c_char) == 0
            || libc::strcmp(arg, b"-S\0".as_ptr() as *const c_char) == 0
        {
            start_file = argv[(i + 1) as usize];
            i += 1;
        } else if libc::strcmp(arg, b"-t\0".as_ptr() as *const c_char) == 0
            || libc::strcmp(arg, b"-T\0".as_ptr() as *const c_char) == 0
        {
            train_file = argv[(i + 1) as usize];
            i += 1;
        } else if libc::strcmp(arg, b"-g\0".as_ptr() as *const c_char) == 0
            || libc::strcmp(arg, b"-G\0".as_ptr() as *const c_char) == 0
        {
            tinf.trans_table = libc::atoi(argv[(i + 1) as usize]);
            if tinf.trans_table < 1
                || tinf.trans_table > 25
                || tinf.trans_table == 7
                || tinf.trans_table == 8
                || (tinf.trans_table >= 17 && tinf.trans_table <= 20)
            {
                usage(
                    b"Invalid translation table specified.\0".as_ptr() as *const c_char,
                );
            }
            user_tt = tinf.trans_table;
            i += 1;
        } else if libc::strcmp(arg, b"-p\0".as_ptr() as *const c_char) == 0
            || libc::strcmp(arg, b"-P\0".as_ptr() as *const c_char) == 0
        {
            let p = *argv[(i + 1) as usize] as u8;
            if p == b'0' || p == b's' || p == b'S' {
                is_meta = 0;
            } else if p == b'1' || p == b'm' || p == b'M' {
                is_meta = 1;
            } else {
                usage(
                    b"Invalid meta/single genome type specified.\0".as_ptr() as *const c_char,
                );
            }
            i += 1;
        } else if libc::strcmp(arg, b"-f\0".as_ptr() as *const c_char) == 0
            || libc::strcmp(arg, b"-F\0".as_ptr() as *const c_char) == 0
        {
            let farg = argv[(i + 1) as usize];
            if libc::strncmp(farg, b"0\0".as_ptr() as *const c_char, 1) == 0
                || libc::strcmp(farg, b"gbk\0".as_ptr() as *const c_char) == 0
                || libc::strcmp(farg, b"GBK\0".as_ptr() as *const c_char) == 0
            {
                output = 0;
            } else if libc::strncmp(farg, b"1\0".as_ptr() as *const c_char, 1) == 0
                || libc::strcmp(farg, b"gca\0".as_ptr() as *const c_char) == 0
                || libc::strcmp(farg, b"GCA\0".as_ptr() as *const c_char) == 0
            {
                output = 1;
            } else if libc::strncmp(farg, b"2\0".as_ptr() as *const c_char, 1) == 0
                || libc::strcmp(farg, b"sco\0".as_ptr() as *const c_char) == 0
                || libc::strcmp(farg, b"SCO\0".as_ptr() as *const c_char) == 0
            {
                output = 2;
            } else if libc::strncmp(farg, b"3\0".as_ptr() as *const c_char, 1) == 0
                || libc::strcmp(farg, b"gff\0".as_ptr() as *const c_char) == 0
                || libc::strcmp(farg, b"GFF\0".as_ptr() as *const c_char) == 0
            {
                output = 3;
            } else {
                usage(
                    b"Invalid output format specified.\0".as_ptr() as *const c_char,
                );
            }
            i += 1;
        } else {
            usage(b"Unknown option.\0".as_ptr() as *const c_char);
        }

        i += 1;
    }

    // Print header
    if quiet == 0 {
        libc::fprintf(
            stderr,
            b"-------------------------------------\n\0".as_ptr() as *const c_char,
        );
        libc::fprintf(
            stderr,
            b"PRODIGAL v%s [%s]         \n\0".as_ptr() as *const c_char,
            VERSION.as_ptr() as *const c_char,
            DATE.as_ptr() as *const c_char,
        );
        libc::fprintf(
            stderr,
            b"Univ of Tenn / Oak Ridge National Lab\n\0".as_ptr() as *const c_char,
        );
        libc::fprintf(
            stderr,
            b"Doug Hyatt, Loren Hauser, et al.     \n\0".as_ptr() as *const c_char,
        );
        libc::fprintf(
            stderr,
            b"-------------------------------------\n\0".as_ptr() as *const c_char,
        );
    }

    // Read in the training file (if specified)
    if !train_file.is_null() {
        if is_meta == 1 {
            libc::fprintf(
                stderr,
                b"\nError: cannot specify metagenomic sequence with a\0".as_ptr()
                    as *const c_char,
            );
            libc::fprintf(
                stderr,
                b" training file.\n\0".as_ptr() as *const c_char,
            );
            return 2;
        }
        rv = read_training_file(train_file, &mut tinf);
        if rv == 1 {
            do_training = 1;
        } else {
            if force_nonsd == 1 {
                libc::fprintf(
                    stderr,
                    b"\nError: cannot force non-SD finder with a training\0".as_ptr()
                        as *const c_char,
                );
                libc::fprintf(
                    stderr,
                    b" file already created!\n\0".as_ptr() as *const c_char,
                );
                return 3;
            }
            if quiet == 0 {
                libc::fprintf(
                    stderr,
                    b"Reading in training data from file %s...\0".as_ptr() as *const c_char,
                    train_file,
                );
            }
            if user_tt > 0 && user_tt != tinf.trans_table {
                libc::fprintf(
                    stderr,
                    b"\n\nWarning: user-specified translation table does\0".as_ptr()
                        as *const c_char,
                );
                libc::fprintf(
                    stderr,
                    b"not match the one in the specified training file! \n\n\0".as_ptr()
                        as *const c_char,
                );
            }
            if rv == -1 {
                libc::fprintf(
                    stderr,
                    b"\n\nError: training file did not read correctly!\n\0".as_ptr()
                        as *const c_char,
                );
                return 4;
            }
            if quiet == 0 {
                libc::fprintf(
                    stderr,
                    b"done!\n\0".as_ptr() as *const c_char,
                );
                libc::fprintf(
                    stderr,
                    b"-------------------------------------\n\0".as_ptr() as *const c_char,
                );
            }
        }
    }

    // Determine where standard input is coming from
    if is_meta == 0 && train_file.is_null() && input_file.is_null() {
        fnum = libc::fileno(libc::fdopen(0, b"r\0".as_ptr() as *const c_char));
        if libc::fstat(fnum, &mut fbuf) == -1 {
            libc::fprintf(
                stderr,
                b"\nError: can't fstat standard input.\n\n\0".as_ptr() as *const c_char,
            );
            return 5;
        }
        if (fbuf.st_mode & libc::S_IFMT) == libc::S_IFCHR {
            help();
        } else if (fbuf.st_mode & libc::S_IFMT) == libc::S_IFREG {
            // do nothing
        } else if (fbuf.st_mode & libc::S_IFMT) == libc::S_IFIFO {
            piped = 1;
            if copy_standard_input_to_file(input_copy.as_ptr(), quiet) == -1 {
                libc::fprintf(
                    stderr,
                    b"\nError: can't copy stdin to file.\n\n\0".as_ptr() as *const c_char,
                );
                return 5;
            }
            input_file = input_copy.as_ptr();
        }
    }

    // Check i/o files and prepare them for reading/writing
    if !input_file.is_null() {
        input_ptr = gzopen(input_file, b"r\0".as_ptr() as *const c_char);
        if input_ptr.is_null() {
            libc::fprintf(
                stderr,
                b"\nError: can't open input file %s.\n\n\0".as_ptr() as *const c_char,
                input_file,
            );
            return 5;
        }
    }
    if input_ptr.is_null() {
        input_ptr = gzopen(
            b"/dev/stdin\0".as_ptr() as *const c_char,
            b"r\0".as_ptr() as *const c_char,
        );
        if input_ptr.is_null() {
            libc::fprintf(
                stderr,
                b"\nError: can't open input file %s.\n\n\0".as_ptr() as *const c_char,
                input_file,
            );
            return 5;
        }
    }
    if !output_file.is_null() {
        output_ptr = libc::fopen(output_file, b"w\0".as_ptr() as *const c_char);
        if output_ptr.is_null() {
            libc::fprintf(
                stderr,
                b"\nError: can't open output file %s.\n\n\0".as_ptr() as *const c_char,
                output_file,
            );
            return 6;
        }
    }
    if !start_file.is_null() {
        start_ptr = libc::fopen(start_file, b"w\0".as_ptr() as *const c_char);
        if start_ptr.is_null() {
            libc::fprintf(
                stderr,
                b"\nError: can't open start file %s.\n\n\0".as_ptr() as *const c_char,
                start_file,
            );
            return 7;
        }
    }
    if !trans_file.is_null() {
        trans_ptr = libc::fopen(trans_file, b"w\0".as_ptr() as *const c_char);
        if trans_ptr.is_null() {
            libc::fprintf(
                stderr,
                b"\nError: can't open translation file %s.\n\n\0".as_ptr() as *const c_char,
                trans_file,
            );
            return 8;
        }
    }
    if !nuc_file.is_null() {
        nuc_ptr = libc::fopen(nuc_file, b"w\0".as_ptr() as *const c_char);
        if nuc_ptr.is_null() {
            libc::fprintf(
                stderr,
                b"\nError: can't open gene nucleotide file %s.\n\n\0".as_ptr() as *const c_char,
                nuc_file,
            );
            return 16;
        }
    }

    // =========================================================================
    // Single Genome Training
    // =========================================================================
    if is_meta == 0 && (do_training == 1 || (do_training == 0 && train_file.is_null())) {
        if quiet == 0 {
            libc::fprintf(
                stderr,
                b"Request:  Single Genome, Phase:  Training\n\0".as_ptr() as *const c_char,
            );
            libc::fprintf(
                stderr,
                b"Reading in the sequence(s) to train...\0".as_ptr() as *const c_char,
            );
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
            libc::fprintf(
                stderr,
                b"\n\nSequence read failed (file must be Fasta, \0".as_ptr() as *const c_char,
            );
            libc::fprintf(
                stderr,
                b"Genbank, or EMBL format).\n\n\0".as_ptr() as *const c_char,
            );
            return 9;
        }
        if slen < MIN_SINGLE_GENOME {
            libc::fprintf(
                stderr,
                b"\n\nError:  Sequence must be %d\0".as_ptr() as *const c_char,
                MIN_SINGLE_GENOME,
            );
            libc::fprintf(
                stderr,
                b" characters (only %d read).\n(Consider\0".as_ptr() as *const c_char,
                slen,
            );
            libc::fprintf(
                stderr,
                b" running with the -p meta option or finding\0".as_ptr() as *const c_char,
            );
            libc::fprintf(
                stderr,
                b" more contigs from the same genome.)\n\n\0".as_ptr() as *const c_char,
            );
            return 10;
        }
        if slen < IDEAL_SINGLE_GENOME {
            libc::fprintf(
                stderr,
                b"\n\nWarning:  ideally Prodigal should be given at\0".as_ptr() as *const c_char,
            );
            libc::fprintf(
                stderr,
                b" least %d bases for \0".as_ptr() as *const c_char,
                IDEAL_SINGLE_GENOME,
            );
            libc::fprintf(
                stderr,
                b"training.\nYou may get better results with the \0".as_ptr() as *const c_char,
            );
            libc::fprintf(
                stderr,
                b"-p meta option.\n\n\0".as_ptr() as *const c_char,
            );
        }
        rcom_seq(seq, rseq, useq, slen);
        if quiet == 0 {
            libc::fprintf(
                stderr,
                b"%d bp seq created, %.2f pct GC\n\0".as_ptr() as *const c_char,
                slen,
                tinf.gc * 100.0,
            );
        }

        // Find all potential starts and stops
        if quiet == 0 {
            libc::fprintf(
                stderr,
                b"Locating all potential starts and stops...\0".as_ptr() as *const c_char,
            );
        }
        if slen > max_slen && slen > (STT_NOD as c_int) * 8 {
            nodes = libc::realloc(
                nodes as *mut c_void,
                (slen as usize / 8) * std::mem::size_of::<Node>(),
            ) as *mut Node;
            if nodes.is_null() {
                libc::fprintf(
                    stderr,
                    b"Realloc failed on nodes\n\n\0".as_ptr() as *const c_char,
                );
                return 11;
            }
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
        libc::qsort(
            nodes as *mut c_void,
            nn as libc::size_t,
            std::mem::size_of::<Node>(),
            Some(compare_nodes),
        );
        if quiet == 0 {
            libc::fprintf(
                stderr,
                b"%d nodes\n\0".as_ptr() as *const c_char,
                nn,
            );
        }

        // Scan ORFs for GC bias
        if quiet == 0 {
            libc::fprintf(
                stderr,
                b"Looking for GC bias in different frames...\0".as_ptr() as *const c_char,
            );
        }
        gc_frame = calc_most_gc_frame(seq, slen);
        if gc_frame.is_null() {
            libc::fprintf(
                stderr,
                b"Malloc failed on gc frame plot\n\n\0".as_ptr() as *const c_char,
            );
            return 11;
        }
        record_gc_bias(gc_frame, nodes, nn, &mut tinf);
        if quiet == 0 {
            libc::fprintf(
                stderr,
                b"frame bias scores: %.2f %.2f %.2f\n\0".as_ptr() as *const c_char,
                tinf.bias[0],
                tinf.bias[1],
                tinf.bias[2],
            );
        }
        libc::free(gc_frame as *mut c_void);

        // Initial DP with GC frame bias
        if quiet == 0 {
            libc::fprintf(
                stderr,
                b"Building initial set of genes to train from...\0".as_ptr() as *const c_char,
            );
        }
        record_overlapping_starts(nodes, nn, &mut tinf, 0);
        ipath = dprog(nodes, nn, &mut tinf, 0);
        if quiet == 0 {
            libc::fprintf(
                stderr,
                b"done!\n\0".as_ptr() as *const c_char,
            );
        }

        // Gather dicodon statistics
        if quiet == 0 {
            libc::fprintf(
                stderr,
                b"Creating coding model and scoring nodes...\0".as_ptr() as *const c_char,
            );
        }
        calc_dicodon_gene(&mut tinf, seq, rseq, slen, nodes, ipath);
        raw_coding_score(seq, rseq, slen, nodes, nn, &mut tinf);
        if quiet == 0 {
            libc::fprintf(
                stderr,
                b"done!\n\0".as_ptr() as *const c_char,
            );
        }

        // Determine SD usage and train starts
        if quiet == 0 {
            libc::fprintf(
                stderr,
                b"Examining upstream regions and training starts...\0".as_ptr() as *const c_char,
            );
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
            libc::fprintf(
                stderr,
                b"done!\n\0".as_ptr() as *const c_char,
            );
        }

        // If training specified, write the training file and exit
        if do_training == 1 {
            if quiet == 0 {
                libc::fprintf(
                    stderr,
                    b"Writing data to training file %s...\0".as_ptr() as *const c_char,
                    train_file,
                );
            }
            rv = write_training_file(train_file, &mut tinf);
            if rv != 0 {
                libc::fprintf(
                    stderr,
                    b"\nError: could not write training file!\n\0".as_ptr() as *const c_char,
                );
                return 12;
            } else {
                if quiet == 0 {
                    libc::fprintf(
                        stderr,
                        b"done!\n\0".as_ptr() as *const c_char,
                    );
                }
                // Clean up and return
                libc::free(seq as *mut c_void);
                libc::free(rseq as *mut c_void);
                libc::free(useq as *mut c_void);
                libc::free(nodes as *mut c_void);
                libc::free(genes as *mut c_void);
                for j in 0..NUM_META {
                    libc::free(meta[j].tinf as *mut c_void);
                }
                gzclose(input_ptr);
                if output_ptr != stdout {
                    libc::fclose(output_ptr);
                }
                if start_ptr != stdout {
                    libc::fclose(start_ptr);
                }
                if trans_ptr != stdout {
                    libc::fclose(trans_ptr);
                }
                return 0;
            }
        }

        // Rewind input file
        if quiet == 0 {
            libc::fprintf(
                stderr,
                b"-------------------------------------\n\0".as_ptr() as *const c_char,
            );
        }
        if gzseek(input_ptr, 0, libc::SEEK_SET) == -1 {
            libc::fprintf(
                stderr,
                b"\nError: could not rewind input file.\n\0".as_ptr() as *const c_char,
            );
            return 13;
        }

        // Reset sequence/DP variables
        libc::memset(seq as *mut c_void, 0, slen as usize / 4 + 1);
        libc::memset(rseq as *mut c_void, 0, slen as usize / 4 + 1);
        libc::memset(useq as *mut c_void, 0, slen as usize / 8 + 1);
        libc::memset(
            nodes as *mut c_void,
            0,
            nn as usize * std::mem::size_of::<Node>(),
        );
        nn = 0;
        slen = 0;
        ipath = 0;
        nmask = 0;
    }
    // Initialize metagenomic bins
    else if is_meta == 1 {
        if quiet == 0 {
            libc::fprintf(
                stderr,
                b"Request:  Metagenomic, Phase:  Training\n\0".as_ptr() as *const c_char,
            );
            libc::fprintf(
                stderr,
                b"Initializing training files...\0".as_ptr() as *const c_char,
            );
        }
        initialize_metagenomic_bins(meta.as_mut_ptr());
        if quiet == 0 {
            libc::fprintf(
                stderr,
                b"done!\n\0".as_ptr() as *const c_char,
            );
            libc::fprintf(
                stderr,
                b"-------------------------------------\n\0".as_ptr() as *const c_char,
            );
        }
    }

    // Print header for gene finding phase
    if quiet == 0 {
        if is_meta == 1 {
            libc::fprintf(
                stderr,
                b"Request:  Metagenomic, Phase:  Gene Finding\n\0".as_ptr() as *const c_char,
            );
        } else {
            libc::fprintf(
                stderr,
                b"Request:  Single Genome, Phase:  Gene Finding\n\0".as_ptr() as *const c_char,
            );
        }
    }

    // Read and process each sequence
    libc::sprintf(
        cur_header.as_mut_ptr(),
        b"Prodigal_Seq_1\0".as_ptr() as *const c_char,
    );
    libc::sprintf(
        new_header.as_mut_ptr(),
        b"Prodigal_Seq_2\0".as_ptr() as *const c_char,
    );

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
            libc::fprintf(
                stderr,
                b"\nSequence read failed (file must be Fasta, \0".as_ptr() as *const c_char,
            );
            libc::fprintf(
                stderr,
                b"Genbank, or EMBL format).\n\n\0".as_ptr() as *const c_char,
            );
            return 14;
        }

        if quiet == 0 {
            libc::fprintf(
                stderr,
                b"Finding genes in sequence #%d (%d bp)...\0".as_ptr() as *const c_char,
                num_seq,
                slen,
            );
        }

        // Reallocate if this is the biggest sequence we've seen
        if slen > max_slen && slen > (STT_NOD as c_int) * 8 {
            nodes = libc::realloc(
                nodes as *mut c_void,
                (slen as usize / 8) * std::mem::size_of::<Node>(),
            ) as *mut Node;
            if nodes.is_null() {
                libc::fprintf(
                    stderr,
                    b"Realloc failed on nodes\n\n\0".as_ptr() as *const c_char,
                );
                return 11;
            }
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
            libc::qsort(
                nodes as *mut c_void,
                nn as libc::size_t,
                std::mem::size_of::<Node>(),
                Some(compare_nodes),
            );

            score_nodes(seq, rseq, slen, nodes, nn, &mut tinf, closed, is_meta);
            if start_ptr != stdout {
                write_start_file(
                    start_ptr,
                    nodes,
                    nn,
                    &mut tinf,
                    num_seq,
                    slen,
                    0,
                    std::ptr::null_mut(),
                    VERSION.as_ptr() as *mut c_char,
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
                libc::fprintf(
                    stderr,
                    b"done!\n\0".as_ptr() as *const c_char,
                );
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
                VERSION.as_ptr() as *mut c_char,
            );
            libc::fflush(output_ptr);
            if trans_ptr != stdout {
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
            if nuc_ptr != stdout {
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
                    libc::memset(
                        nodes as *mut c_void,
                        0,
                        nn as usize * std::mem::size_of::<Node>(),
                    );
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
                    libc::qsort(
                        nodes as *mut c_void,
                        nn as libc::size_t,
                        std::mem::size_of::<Node>(),
                        Some(compare_nodes),
                    );
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
            libc::memset(
                nodes as *mut c_void,
                0,
                nn as usize * std::mem::size_of::<Node>(),
            );
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
            libc::qsort(
                nodes as *mut c_void,
                nn as libc::size_t,
                std::mem::size_of::<Node>(),
                Some(compare_nodes),
            );
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
            if start_ptr != stdout {
                write_start_file(
                    start_ptr,
                    nodes,
                    nn,
                    meta[max_phase as usize].tinf,
                    num_seq,
                    slen,
                    1,
                    meta[max_phase as usize].desc.as_mut_ptr(),
                    VERSION.as_ptr() as *mut c_char,
                    cur_header.as_mut_ptr(),
                );
            }

            if quiet == 0 {
                libc::fprintf(
                    stderr,
                    b"done!\n\0".as_ptr() as *const c_char,
                );
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
                VERSION.as_ptr() as *mut c_char,
            );
            libc::fflush(output_ptr);
            if trans_ptr != stdout {
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
            if nuc_ptr != stdout {
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
        libc::memset(seq as *mut c_void, 0, slen as usize / 4 + 1);
        libc::memset(rseq as *mut c_void, 0, slen as usize / 4 + 1);
        libc::memset(useq as *mut c_void, 0, slen as usize / 8 + 1);
        libc::memset(
            nodes as *mut c_void,
            0,
            nn as usize * std::mem::size_of::<Node>(),
        );
        nn = 0;
        slen = 0;
        ipath = 0;
        nmask = 0;
        libc::strcpy(cur_header.as_mut_ptr(), new_header.as_ptr());
        libc::sprintf(
            new_header.as_mut_ptr(),
            b"Prodigal_Seq_%d\n\0".as_ptr() as *const c_char,
            num_seq + 1,
        );
    }

    if num_seq == 0 {
        libc::fprintf(
            stderr,
            b"\nError:  no input sequences to analyze.\n\n\0".as_ptr() as *const c_char,
        );
        return 18;
    }

    // Free all memory
    libc::free(seq as *mut c_void);
    libc::free(rseq as *mut c_void);
    libc::free(useq as *mut c_void);
    libc::free(nodes as *mut c_void);
    libc::free(genes as *mut c_void);
    for j in 0..NUM_META {
        libc::free(meta[j].tinf as *mut c_void);
    }

    // Close all filehandles
    gzclose(input_ptr);
    if output_ptr != stdout {
        libc::fclose(output_ptr);
    }
    if start_ptr != stdout {
        libc::fclose(start_ptr);
    }
    if trans_ptr != stdout {
        libc::fclose(trans_ptr);
    }
    if nuc_ptr != stdout {
        libc::fclose(nuc_ptr);
    }

    // Remove tmp file
    if piped == 1 {
        if libc::remove(input_copy.as_ptr()) != 0 {
            libc::fprintf(
                stderr,
                b"Could not delete tmp file %s.\n\0".as_ptr() as *const c_char,
                input_copy.as_ptr(),
            );
            return 18;
        }
    }

    0
}
