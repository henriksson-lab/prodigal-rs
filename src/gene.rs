/*******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Copyright (C) 2007-2016 University of Tennessee / UT-Battelle

    Code Author:  Doug Hyatt

    Rust translation maintaining exact C ABI and output compatibility.
*******************************************************************************/

use std::os::raw::{c_char, c_double, c_int};

use crate::types::{Gene, Node, Training, MAX_GENES, MAX_LINE, MAX_SAM_OVLP, STOP};

extern "C" {
    fn intergenic_mod(n1: *mut Node, n2: *mut Node, tinf: *mut Training) -> f64;
    fn is_a(seq: *const u8, n: c_int) -> c_int;
    fn is_c(seq: *const u8, n: c_int) -> c_int;
    fn is_g(seq: *const u8, n: c_int) -> c_int;
    fn is_t(seq: *const u8, n: c_int) -> c_int;
    fn is_n(seq: *const u8, n: c_int) -> c_int;
    fn amino(seq: *const u8, n: c_int, tinf: *mut Training, is_init: c_int) -> c_char;
    fn mer_text(qt: *mut c_char, len: c_int, ndx: c_int);
    // fn start_text(st: *mut c_char, type_: c_int);

    static stderr: *mut libc::FILE;
}

/// Copies genes from the dynamic programming to a final array.
#[no_mangle]
pub unsafe extern "C" fn add_genes(
    glist: *mut Gene,
    nod: *mut Node,
    dbeg: c_int,
) -> c_int {
    if dbeg == -1 {
        return 0;
    }

    let mut path = dbeg;
    let mut ctr: c_int = 0;

    while (*nod.offset(path as isize)).traceb != -1 {
        path = (*nod.offset(path as isize)).traceb;
    }

    while path != -1 {
        let np = &*nod.offset(path as isize);
        if np.elim == 1 {
            path = np.tracef;
            continue;
        }
        if np.strand == 1 && np.type_ != STOP {
            (*glist.offset(ctr as isize)).begin = np.ndx + 1;
            (*glist.offset(ctr as isize)).start_ndx = path;
        }
        if np.strand == -1 && np.type_ == STOP {
            (*glist.offset(ctr as isize)).begin = np.ndx - 1;
            (*glist.offset(ctr as isize)).stop_ndx = path;
        }
        if np.strand == 1 && np.type_ == STOP {
            (*glist.offset(ctr as isize)).end = np.ndx + 3;
            (*glist.offset(ctr as isize)).stop_ndx = path;
            ctr += 1;
        }
        if np.strand == -1 && np.type_ != STOP {
            (*glist.offset(ctr as isize)).end = np.ndx + 1;
            (*glist.offset(ctr as isize)).start_ndx = path;
            ctr += 1;
        }
        path = (*nod.offset(path as isize)).tracef;
        if ctr == MAX_GENES as c_int {
            libc::fprintf(
                stderr,
                b"warning, max # of genes exceeded, truncating...\n\0".as_ptr() as *const c_char,
            );
            return ctr;
        }
    }
    ctr
}

/// Tweak final starts to improve gene predictions.
#[no_mangle]
pub unsafe extern "C" fn tweak_final_starts(
    genes: *mut Gene,
    ng: c_int,
    nod: *mut Node,
    nn: c_int,
    tinf: *mut Training,
) {
    let mut maxndx: [c_int; 2];
    let mut maxsc: [f64; 2];
    let mut maxigm: [f64; 2];

    for i in 0..ng {
        let ndx = (*genes.offset(i as isize)).start_ndx;
        let sc = (*nod.offset(ndx as isize)).sscore + (*nod.offset(ndx as isize)).cscore;
        let mut igm: f64 = 0.0;

        if i > 0
            && (*nod.offset(ndx as isize)).strand == 1
            && (*nod.offset((*genes.offset((i - 1) as isize)).start_ndx as isize)).strand == 1
        {
            igm = intergenic_mod(
                nod.offset((*genes.offset((i - 1) as isize)).stop_ndx as isize),
                nod.offset(ndx as isize),
                tinf,
            );
        }
        if i > 0
            && (*nod.offset(ndx as isize)).strand == 1
            && (*nod.offset((*genes.offset((i - 1) as isize)).start_ndx as isize)).strand == -1
        {
            igm = intergenic_mod(
                nod.offset((*genes.offset((i - 1) as isize)).start_ndx as isize),
                nod.offset(ndx as isize),
                tinf,
            );
        }
        if i < ng - 1
            && (*nod.offset(ndx as isize)).strand == -1
            && (*nod.offset((*genes.offset((i + 1) as isize)).start_ndx as isize)).strand == 1
        {
            igm = intergenic_mod(
                nod.offset(ndx as isize),
                nod.offset((*genes.offset((i + 1) as isize)).start_ndx as isize),
                tinf,
            );
        }
        if i < ng - 1
            && (*nod.offset(ndx as isize)).strand == -1
            && (*nod.offset((*genes.offset((i + 1) as isize)).start_ndx as isize)).strand == -1
        {
            igm = intergenic_mod(
                nod.offset(ndx as isize),
                nod.offset((*genes.offset((i + 1) as isize)).stop_ndx as isize),
                tinf,
            );
        }

        // Search upstream and downstream for the #2 and #3 scoring starts
        maxndx = [-1, -1];
        maxsc = [0.0, 0.0];
        maxigm = [0.0, 0.0];

        for j in (ndx - 100)..(ndx + 100) {
            if j < 0 || j >= nn || j == ndx {
                continue;
            }
            let nj = &*nod.offset(j as isize);
            if nj.type_ == STOP || nj.stop_val != (*nod.offset(ndx as isize)).stop_val {
                continue;
            }

            let mut tigm: f64 = 0.0;
            if i > 0
                && nj.strand == 1
                && (*nod.offset((*genes.offset((i - 1) as isize)).start_ndx as isize)).strand == 1
            {
                if (*nod.offset((*genes.offset((i - 1) as isize)).stop_ndx as isize)).ndx
                    - nj.ndx
                    > MAX_SAM_OVLP
                {
                    continue;
                }
                tigm = intergenic_mod(
                    nod.offset((*genes.offset((i - 1) as isize)).stop_ndx as isize),
                    nod.offset(j as isize),
                    tinf,
                );
            }
            if i > 0
                && nj.strand == 1
                && (*nod.offset((*genes.offset((i - 1) as isize)).start_ndx as isize)).strand
                    == -1
            {
                if (*nod.offset((*genes.offset((i - 1) as isize)).start_ndx as isize)).ndx
                    - nj.ndx
                    >= 0
                {
                    continue;
                }
                tigm = intergenic_mod(
                    nod.offset((*genes.offset((i - 1) as isize)).start_ndx as isize),
                    nod.offset(j as isize),
                    tinf,
                );
            }
            if i < ng - 1
                && nj.strand == -1
                && (*nod.offset((*genes.offset((i + 1) as isize)).start_ndx as isize)).strand == 1
            {
                if nj.ndx
                    - (*nod.offset((*genes.offset((i + 1) as isize)).start_ndx as isize)).ndx
                    >= 0
                {
                    continue;
                }
                tigm = intergenic_mod(
                    nod.offset(j as isize),
                    nod.offset((*genes.offset((i + 1) as isize)).start_ndx as isize),
                    tinf,
                );
            }
            if i < ng - 1
                && nj.strand == -1
                && (*nod.offset((*genes.offset((i + 1) as isize)).start_ndx as isize)).strand
                    == -1
            {
                if nj.ndx
                    - (*nod.offset((*genes.offset((i + 1) as isize)).stop_ndx as isize)).ndx
                    > MAX_SAM_OVLP
                {
                    continue;
                }
                tigm = intergenic_mod(
                    nod.offset(j as isize),
                    nod.offset((*genes.offset((i + 1) as isize)).stop_ndx as isize),
                    tinf,
                );
            }

            if maxndx[0] == -1 {
                maxndx[0] = j;
                maxsc[0] = nj.cscore + nj.sscore;
                maxigm[0] = tigm;
            } else if nj.cscore + nj.sscore + tigm > maxsc[0] {
                maxndx[1] = maxndx[0];
                maxsc[1] = maxsc[0];
                maxigm[1] = maxigm[0];
                maxndx[0] = j;
                maxsc[0] = nj.cscore + nj.sscore;
                maxigm[0] = tigm;
            } else if maxndx[1] == -1 || nj.cscore + nj.sscore + tigm > maxsc[1] {
                maxndx[1] = j;
                maxsc[1] = nj.cscore + nj.sscore;
                maxigm[1] = tigm;
            }
        }

        // Change the start if it's a TTG with better coding/RBS/upstream score
        // Also change the start if it's <=15bp but has better coding/RBS
        for j in 0..2 {
            let mndx = maxndx[j];
            if mndx == -1 {
                continue;
            }

            let nm = &*nod.offset(mndx as isize);
            let nn_ref = &*nod.offset(ndx as isize);

            // Start of less common type but with better coding, rbs, and
            // upstream.  Must be 18 or more bases away from original.
            if nm.tscore < nn_ref.tscore
                && maxsc[j] - nm.tscore >= sc - nn_ref.tscore + (*tinf).st_wt
                && nm.rscore > nn_ref.rscore
                && nm.uscore > nn_ref.uscore
                && nm.cscore > nn_ref.cscore
                && (nm.ndx - nn_ref.ndx).abs() > 15
            {
                maxsc[j] += nn_ref.tscore - nm.tscore;
            }
            // Close starts.  Ignore coding and see if start has better rbs
            // and type.
            else if (nm.ndx - nn_ref.ndx).abs() <= 15
                && nm.rscore + nm.tscore > nn_ref.rscore + nn_ref.tscore
                && nn_ref.edge == 0
                && nm.edge == 0
            {
                if nn_ref.cscore > nm.cscore {
                    maxsc[j] += nn_ref.cscore - nm.cscore;
                }
                if nn_ref.uscore > nm.uscore {
                    maxsc[j] += nn_ref.uscore - nm.uscore;
                }
                if igm > maxigm[j] {
                    maxsc[j] += igm - maxigm[j];
                }
            } else {
                maxsc[j] = -1000.0;
            }
        }

        // Change the gene coordinates to the new maximum.
        let mut mndx: c_int = -1;
        for j in 0..2i32 {
            if maxndx[j as usize] == -1 {
                continue;
            }
            if mndx == -1 && maxsc[j as usize] + maxigm[j as usize] > sc + igm {
                mndx = j;
            } else if mndx >= 0
                && maxsc[j as usize] + maxigm[j as usize]
                    > maxsc[mndx as usize] + maxigm[mndx as usize]
            {
                mndx = j;
            }
        }
        if mndx != -1
            && (*nod.offset(maxndx[mndx as usize] as isize)).strand == 1
        {
            (*genes.offset(i as isize)).start_ndx = maxndx[mndx as usize];
            (*genes.offset(i as isize)).begin =
                (*nod.offset(maxndx[mndx as usize] as isize)).ndx + 1;
        } else if mndx != -1
            && (*nod.offset(maxndx[mndx as usize] as isize)).strand == -1
        {
            (*genes.offset(i as isize)).start_ndx = maxndx[mndx as usize];
            (*genes.offset(i as isize)).end =
                (*nod.offset(maxndx[mndx as usize] as isize)).ndx + 1;
        }
    }
}

/// Record gene data strings into the gene_data and score_data fields.
#[no_mangle]
pub unsafe extern "C" fn record_gene_data(
    genes: *mut Gene,
    ng: c_int,
    nod: *mut Node,
    tinf: *mut Training,
    sctr: c_int,
) {
    // SD motif string tables
    let sd_string: [&[u8]; 28] = [
        b"None\0",
        b"GGA/GAG/AGG\0",
        b"3Base/5BMM\0",
        b"4Base/6BMM\0",
        b"AGxAG\0",
        b"AGxAG\0",
        b"GGA/GAG/AGG\0",
        b"GGxGG\0",
        b"GGxGG\0",
        b"AGxAG\0",
        b"AGGAG(G)/GGAGG\0",
        b"AGGA/GGAG/GAGG\0",
        b"AGGA/GGAG/GAGG\0",
        b"GGA/GAG/AGG\0",
        b"GGxGG\0",
        b"AGGA\0",
        b"GGAG/GAGG\0",
        b"AGxAGG/AGGxGG\0",
        b"AGxAGG/AGGxGG\0",
        b"AGxAGG/AGGxGG\0",
        b"AGGAG/GGAGG\0",
        b"AGGAG\0",
        b"AGGAG\0",
        b"GGAGG\0",
        b"GGAGG\0",
        b"AGGAGG\0",
        b"AGGAGG\0",
        b"AGGAGG\0",
    ];

    let sd_spacer: [&[u8]; 28] = [
        b"None\0",
        b"3-4bp\0",
        b"13-15bp\0",
        b"13-15bp\0",
        b"11-12bp\0",
        b"3-4bp\0",
        b"11-12bp\0",
        b"11-12bp\0",
        b"3-4bp\0",
        b"5-10bp\0",
        b"13-15bp\0",
        b"3-4bp\0",
        b"11-12bp\0",
        b"5-10bp\0",
        b"5-10bp\0",
        b"5-10bp\0",
        b"5-10bp\0",
        b"11-12bp\0",
        b"3-4bp\0",
        b"5-10bp\0",
        b"11-12bp\0",
        b"3-4bp\0",
        b"5-10bp\0",
        b"3-4bp\0",
        b"5-10bp\0",
        b"11-12bp\0",
        b"3-4bp\0",
        b"5-10bp\0",
    ];

    let type_string: [&[u8]; 4] = [b"ATG\0", b"GTG\0", b"TTG\0", b"Edge\0"];

    let mut buffer: [c_char; 500] = [0; 500];
    let mut qt: [c_char; 10] = [0; 10];

    for i in 0..ng {
        let gi = &mut *genes.offset(i as isize);
        let ndx = gi.start_ndx;
        let sndx = gi.stop_ndx;
        let n = &*nod.offset(ndx as isize);
        let sn = &*nod.offset(sndx as isize);

        // Record basic gene data
        let partial_left: c_int = if (n.edge == 1 && n.strand == 1)
            || (sn.edge == 1 && n.strand == -1)
        {
            1
        } else {
            0
        };
        let partial_right: c_int = if (sn.edge == 1 && n.strand == 1)
            || (n.edge == 1 && n.strand == -1)
        {
            1
        } else {
            0
        };
        let st_type: c_int = if n.edge == 1 { 3 } else { n.type_ };

        libc::sprintf(
            gi.gene_data.as_mut_ptr(),
            b"ID=%d_%d;partial=%d%d;start_type=%s;\0".as_ptr() as *const c_char,
            sctr,
            i + 1,
            partial_left,
            partial_right,
            type_string[st_type as usize].as_ptr() as *const c_char,
        );

        // Record rbs data
        let rbs1 = (*tinf).rbs_wt[n.rbs[0] as usize] * (*tinf).st_wt;
        let rbs2 = (*tinf).rbs_wt[n.rbs[1] as usize] * (*tinf).st_wt;

        if (*tinf).uses_sd == 1 {
            if rbs1 > rbs2 {
                libc::sprintf(
                    buffer.as_mut_ptr(),
                    b"rbs_motif=%s;rbs_spacer=%s\0".as_ptr() as *const c_char,
                    sd_string[n.rbs[0] as usize].as_ptr() as *const c_char,
                    sd_spacer[n.rbs[0] as usize].as_ptr() as *const c_char,
                );
                libc::strcat(
                    gi.gene_data.as_mut_ptr(),
                    buffer.as_ptr(),
                );
            } else {
                libc::sprintf(
                    buffer.as_mut_ptr(),
                    b"rbs_motif=%s;rbs_spacer=%s\0".as_ptr() as *const c_char,
                    sd_string[n.rbs[1] as usize].as_ptr() as *const c_char,
                    sd_spacer[n.rbs[1] as usize].as_ptr() as *const c_char,
                );
                libc::strcat(
                    gi.gene_data.as_mut_ptr(),
                    buffer.as_ptr(),
                );
            }
        } else {
            mer_text(qt.as_mut_ptr(), n.mot.len, n.mot.ndx);
            if (*tinf).no_mot > -0.5
                && rbs1 > rbs2
                && rbs1 > n.mot.score * (*tinf).st_wt
            {
                libc::sprintf(
                    buffer.as_mut_ptr(),
                    b"rbs_motif=%s;rbs_spacer=%s\0".as_ptr() as *const c_char,
                    sd_string[n.rbs[0] as usize].as_ptr() as *const c_char,
                    sd_spacer[n.rbs[0] as usize].as_ptr() as *const c_char,
                );
                libc::strcat(
                    gi.gene_data.as_mut_ptr(),
                    buffer.as_ptr(),
                );
            } else if (*tinf).no_mot > -0.5
                && rbs2 >= rbs1
                && rbs2 > n.mot.score * (*tinf).st_wt
            {
                libc::sprintf(
                    buffer.as_mut_ptr(),
                    b"rbs_motif=%s;rbs_spacer=%s\0".as_ptr() as *const c_char,
                    sd_string[n.rbs[1] as usize].as_ptr() as *const c_char,
                    sd_spacer[n.rbs[1] as usize].as_ptr() as *const c_char,
                );
                libc::strcat(
                    gi.gene_data.as_mut_ptr(),
                    buffer.as_ptr(),
                );
            } else if n.mot.len == 0 {
                libc::strcat(
                    gi.gene_data.as_mut_ptr(),
                    b"rbs_motif=None;rbs_spacer=None\0".as_ptr() as *const c_char,
                );
            } else {
                libc::sprintf(
                    buffer.as_mut_ptr(),
                    b"rbs_motif=%s;rbs_spacer=%dbp\0".as_ptr() as *const c_char,
                    qt.as_ptr() as *const c_char,
                    n.mot.spacer,
                );
                libc::strcat(
                    gi.gene_data.as_mut_ptr(),
                    buffer.as_ptr(),
                );
            }
        }
        libc::sprintf(
            buffer.as_mut_ptr(),
            b";gc_cont=%.3f\0".as_ptr() as *const c_char,
            n.gc_cont as c_double,
        );
        libc::strcat(
            gi.gene_data.as_mut_ptr(),
            buffer.as_ptr(),
        );

        // Record score data
        let confidence = calculate_confidence(n.cscore + n.sscore, (*tinf).st_wt);
        libc::sprintf(
            gi.score_data.as_mut_ptr(),
            b"conf=%.2f;score=%.2f;cscore=%.2f;sscore=%.2f;rscore=%.2f;uscore=%.2f;\0".as_ptr()
                as *const c_char,
            confidence as c_double,
            (n.cscore + n.sscore) as c_double,
            n.cscore as c_double,
            n.sscore as c_double,
            n.rscore as c_double,
            n.uscore as c_double,
        );

        libc::sprintf(
            buffer.as_mut_ptr(),
            b"tscore=%.2f;\0".as_ptr() as *const c_char,
            n.tscore as c_double,
        );
        libc::strcat(
            gi.score_data.as_mut_ptr(),
            buffer.as_ptr(),
        );
    }
}

/// Print the genes. 'flag' indicates which format to use.
#[no_mangle]
pub unsafe extern "C" fn print_genes(
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
) {
    let mut left: [c_char; 50] = [0; 50];
    let mut right: [c_char; 50] = [0; 50];
    let mut seq_data: [c_char; MAX_LINE * 2] = [0; MAX_LINE * 2];
    let mut run_data: [c_char; MAX_LINE] = [0; MAX_LINE];
    let mut buffer: [c_char; MAX_LINE] = [0; MAX_LINE];

    // Initialize sequence data
    libc::sprintf(
        seq_data.as_mut_ptr(),
        b"seqnum=%d;seqlen=%d;seqhdr=\"%s\"\0".as_ptr() as *const c_char,
        sctr,
        slen,
        header,
    );

    // Initialize run data string
    if is_meta == 0 {
        libc::sprintf(
            run_data.as_mut_ptr(),
            b"version=Prodigal.v%s;run_type=Single;\0".as_ptr() as *const c_char,
            version,
        );
        libc::strcat(
            run_data.as_mut_ptr(),
            b"model=\"Ab initio\";\0".as_ptr() as *const c_char,
        );
    } else {
        libc::sprintf(
            run_data.as_mut_ptr(),
            b"version=Prodigal.v%s;run_type=Metagenomic;\0".as_ptr() as *const c_char,
            version,
        );
        libc::sprintf(
            buffer.as_mut_ptr(),
            b"model=\"%s\";\0".as_ptr() as *const c_char,
            mdesc,
        );
        libc::strcat(run_data.as_mut_ptr(), buffer.as_ptr());
    }
    libc::sprintf(
        buffer.as_mut_ptr(),
        b"gc_cont=%.2f;transl_table=%d;uses_sd=%d\0".as_ptr() as *const c_char,
        ((*tinf).gc * 100.0) as c_double,
        (*tinf).trans_table,
        (*tinf).uses_sd,
    );
    libc::strcat(run_data.as_mut_ptr(), buffer.as_ptr());

    left[0] = 0;
    right[0] = 0;

    // Print the gff header once
    if flag == 3 && sctr == 1 {
        libc::fprintf(fp, b"##gff-version  3\n\0".as_ptr() as *const c_char);
    }

    // Print sequence/model information
    if flag == 0 {
        libc::fprintf(
            fp,
            b"DEFINITION  %s;%s\n\0".as_ptr() as *const c_char,
            seq_data.as_ptr(),
            run_data.as_ptr(),
        );
        libc::fprintf(
            fp,
            b"FEATURES             Location/Qualifiers\n\0".as_ptr() as *const c_char,
        );
    } else if flag != 1 {
        libc::fprintf(
            fp,
            b"# Sequence Data: %s\n\0".as_ptr() as *const c_char,
            seq_data.as_ptr(),
        );
        libc::fprintf(
            fp,
            b"# Model Data: %s\n\0".as_ptr() as *const c_char,
            run_data.as_ptr(),
        );
    }

    // Print the genes
    for i in 0..ng {
        let gi = &*genes.offset(i as isize);
        let ndx = gi.start_ndx;
        let sndx = gi.stop_ndx;
        let n = &*nod.offset(ndx as isize);
        let sn = &*nod.offset(sndx as isize);

        // Print the coordinates and data
        if n.strand == 1 {
            if n.edge == 1 {
                libc::sprintf(
                    left.as_mut_ptr(),
                    b"<%d\0".as_ptr() as *const c_char,
                    gi.begin,
                );
            } else {
                libc::sprintf(
                    left.as_mut_ptr(),
                    b"%d\0".as_ptr() as *const c_char,
                    gi.begin,
                );
            }
            if sn.edge == 1 {
                libc::sprintf(
                    right.as_mut_ptr(),
                    b">%d\0".as_ptr() as *const c_char,
                    gi.end,
                );
            } else {
                libc::sprintf(
                    right.as_mut_ptr(),
                    b"%d\0".as_ptr() as *const c_char,
                    gi.end,
                );
            }

            if flag == 0 {
                libc::fprintf(
                    fp,
                    b"     CDS             %s..%s\n\0".as_ptr() as *const c_char,
                    left.as_ptr(),
                    right.as_ptr(),
                );
                libc::fprintf(
                    fp,
                    b"                     \0".as_ptr() as *const c_char,
                );
                libc::fprintf(
                    fp,
                    b"/note=\"%s;%s\"\n\0".as_ptr() as *const c_char,
                    gi.gene_data.as_ptr(),
                    gi.score_data.as_ptr(),
                );
            }
            if flag == 1 {
                libc::fprintf(
                    fp,
                    b"gene_prodigal=%d|1|f|y|y|3|0|%d|%d|%d|%d|-1|-1|1.0\n\0".as_ptr()
                        as *const c_char,
                    i + 1,
                    gi.begin,
                    gi.end,
                    gi.begin,
                    gi.end,
                );
            }
            if flag == 2 {
                libc::fprintf(
                    fp,
                    b">%d_%d_%d_+\n\0".as_ptr() as *const c_char,
                    i + 1,
                    gi.begin,
                    gi.end,
                );
            }
            if flag == 3 {
                libc::fprintf(
                    fp,
                    b"%s\tProdigal_v%s\tCDS\t%d\t%d\t%.1f\t+\t0\t%s;%s\0".as_ptr()
                        as *const c_char,
                    short_hdr,
                    version,
                    gi.begin,
                    gi.end,
                    (n.cscore + n.sscore) as c_double,
                    gi.gene_data.as_ptr(),
                    gi.score_data.as_ptr(),
                );
                libc::fprintf(fp, b"\n\0".as_ptr() as *const c_char);
            }
        } else {
            if sn.edge == 1 {
                libc::sprintf(
                    left.as_mut_ptr(),
                    b"<%d\0".as_ptr() as *const c_char,
                    gi.begin,
                );
            } else {
                libc::sprintf(
                    left.as_mut_ptr(),
                    b"%d\0".as_ptr() as *const c_char,
                    gi.begin,
                );
            }
            if n.edge == 1 {
                libc::sprintf(
                    right.as_mut_ptr(),
                    b">%d\0".as_ptr() as *const c_char,
                    gi.end,
                );
            } else {
                libc::sprintf(
                    right.as_mut_ptr(),
                    b"%d\0".as_ptr() as *const c_char,
                    gi.end,
                );
            }

            if flag == 0 {
                libc::fprintf(
                    fp,
                    b"     CDS             complement(%s..%s)\n\0".as_ptr() as *const c_char,
                    left.as_ptr(),
                    right.as_ptr(),
                );
                libc::fprintf(
                    fp,
                    b"                     \0".as_ptr() as *const c_char,
                );
                libc::fprintf(
                    fp,
                    b"/note=\"%s;%s\"\n\0".as_ptr() as *const c_char,
                    gi.gene_data.as_ptr(),
                    gi.score_data.as_ptr(),
                );
            }
            if flag == 1 {
                libc::fprintf(
                    fp,
                    b"gene_prodigal=%d|1|r|y|y|3|0|%d|%d|%d|%d|-1|-1|1.0\n\0".as_ptr()
                        as *const c_char,
                    i + 1,
                    slen + 1 - gi.end,
                    slen + 1 - gi.begin,
                    slen + 1 - gi.end,
                    slen + 1 - gi.begin,
                );
            }
            if flag == 2 {
                libc::fprintf(
                    fp,
                    b">%d_%d_%d_-\n\0".as_ptr() as *const c_char,
                    i + 1,
                    gi.begin,
                    gi.end,
                );
            }
            if flag == 3 {
                libc::fprintf(
                    fp,
                    b"%s\tProdigal_v%s\tCDS\t%d\t%d\t%.1f\t-\t0\t%s;%s\0".as_ptr()
                        as *const c_char,
                    short_hdr,
                    version,
                    gi.begin,
                    gi.end,
                    (n.cscore + n.sscore) as c_double,
                    gi.gene_data.as_ptr(),
                    gi.score_data.as_ptr(),
                );
                libc::fprintf(fp, b"\n\0".as_ptr() as *const c_char);
            }
        }
    }

    // Footer
    if flag == 0 {
        libc::fprintf(fp, b"//\n\0".as_ptr() as *const c_char);
    }
}

/// Print the gene translations.
#[no_mangle]
pub unsafe extern "C" fn write_translations(
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
) {
    for i in 0..ng {
        let gi = &*genes.offset(i as isize);
        if (*nod.offset(gi.start_ndx as isize)).strand == 1 {
            libc::fprintf(
                fh,
                b">%s_%d # %d # %d # 1 # %s\n\0".as_ptr() as *const c_char,
                short_hdr,
                i + 1,
                gi.begin,
                gi.end,
                gi.gene_data.as_ptr(),
            );
            let mut j = gi.begin;
            while j < gi.end {
                if is_n(useq, j - 1) == 1 || is_n(useq, j) == 1 || is_n(useq, j + 1) == 1 {
                    libc::fprintf(fh, b"X\0".as_ptr() as *const c_char);
                } else {
                    let is_init = if j == gi.begin { 1 } else { 0 };
                    let edge_val = 1 - (*nod.offset(gi.start_ndx as isize)).edge;
                    libc::fprintf(
                        fh,
                        b"%c\0".as_ptr() as *const c_char,
                        amino(seq, j - 1, tinf, is_init & edge_val) as c_int,
                    );
                }
                if (j - gi.begin) % 180 == 177 {
                    libc::fprintf(fh, b"\n\0".as_ptr() as *const c_char);
                }
                j += 3;
            }
            if (j - gi.begin) % 180 != 0 {
                libc::fprintf(fh, b"\n\0".as_ptr() as *const c_char);
            }
        } else {
            libc::fprintf(
                fh,
                b">%s_%d # %d # %d # -1 # %s\n\0".as_ptr() as *const c_char,
                short_hdr,
                i + 1,
                gi.begin,
                gi.end,
                gi.gene_data.as_ptr(),
            );
            let mut j = slen + 1 - gi.end;
            while j < slen + 1 - gi.begin {
                if is_n(useq, slen - j) == 1
                    || is_n(useq, slen - 1 - j) == 1
                    || is_n(useq, slen - 2 - j) == 1
                {
                    libc::fprintf(fh, b"X\0".as_ptr() as *const c_char);
                } else {
                    let is_init = if j == slen + 1 - gi.end { 1 } else { 0 };
                    let edge_val = 1 - (*nod.offset(gi.start_ndx as isize)).edge;
                    libc::fprintf(
                        fh,
                        b"%c\0".as_ptr() as *const c_char,
                        amino(rseq, j - 1, tinf, is_init & edge_val) as c_int,
                    );
                }
                if (j - slen - 1 + gi.end) % 180 == 177 {
                    libc::fprintf(fh, b"\n\0".as_ptr() as *const c_char);
                }
                j += 3;
            }
            if (j - slen - 1 + gi.end) % 180 != 0 {
                libc::fprintf(fh, b"\n\0".as_ptr() as *const c_char);
            }
        }
    }
}

/// Print the gene nucleotide sequences.
#[no_mangle]
pub unsafe extern "C" fn write_nucleotide_seqs(
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
) {
    for i in 0..ng {
        let gi = &*genes.offset(i as isize);
        if (*nod.offset(gi.start_ndx as isize)).strand == 1 {
            libc::fprintf(
                fh,
                b">%s_%d # %d # %d # 1 # %s\n\0".as_ptr() as *const c_char,
                short_hdr,
                i + 1,
                gi.begin,
                gi.end,
                gi.gene_data.as_ptr(),
            );
            let mut j = gi.begin - 1;
            while j < gi.end {
                if is_a(seq, j) == 1 {
                    libc::fprintf(fh, b"A\0".as_ptr() as *const c_char);
                } else if is_t(seq, j) == 1 {
                    libc::fprintf(fh, b"T\0".as_ptr() as *const c_char);
                } else if is_g(seq, j) == 1 {
                    libc::fprintf(fh, b"G\0".as_ptr() as *const c_char);
                } else if is_c(seq, j) == 1 && is_n(useq, j) == 0 {
                    libc::fprintf(fh, b"C\0".as_ptr() as *const c_char);
                } else {
                    libc::fprintf(fh, b"N\0".as_ptr() as *const c_char);
                }
                if (j - gi.begin + 1) % 70 == 69 {
                    libc::fprintf(fh, b"\n\0".as_ptr() as *const c_char);
                }
                j += 1;
            }
            if (j - gi.begin + 1) % 70 != 0 {
                libc::fprintf(fh, b"\n\0".as_ptr() as *const c_char);
            }
        } else {
            libc::fprintf(
                fh,
                b">%s_%d # %d # %d # -1 # %s\n\0".as_ptr() as *const c_char,
                short_hdr,
                i + 1,
                gi.begin,
                gi.end,
                gi.gene_data.as_ptr(),
            );
            let mut j = slen - gi.end;
            while j < slen + 1 - gi.begin {
                if is_a(rseq, j) == 1 {
                    libc::fprintf(fh, b"A\0".as_ptr() as *const c_char);
                } else if is_t(rseq, j) == 1 {
                    libc::fprintf(fh, b"T\0".as_ptr() as *const c_char);
                } else if is_g(rseq, j) == 1 {
                    libc::fprintf(fh, b"G\0".as_ptr() as *const c_char);
                } else if is_c(rseq, j) == 1 && is_n(useq, slen - 1 - j) == 0 {
                    libc::fprintf(fh, b"C\0".as_ptr() as *const c_char);
                } else {
                    libc::fprintf(fh, b"N\0".as_ptr() as *const c_char);
                }
                if (j - slen + gi.end) % 70 == 69 {
                    libc::fprintf(fh, b"\n\0".as_ptr() as *const c_char);
                }
                j += 1;
            }
            if (j - slen + gi.end) % 70 != 0 {
                libc::fprintf(fh, b"\n\0".as_ptr() as *const c_char);
            }
        }
    }
}

/// Convert score to a percent confidence.
#[no_mangle]
pub unsafe extern "C" fn calculate_confidence(score: f64, start_weight: f64) -> f64 {
    let mut conf: f64;

    if score / start_weight < 41.0 {
        conf = (score / start_weight).exp();
        conf = (conf / (conf + 1.0)) * 100.0;
    } else {
        conf = 99.99;
    }
    if conf <= 50.00 {
        conf = 50.00;
    }
    conf
}
