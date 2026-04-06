/*******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Copyright (C) 2007-2016 University of Tennessee / UT-Battelle

    Code Author:  Doug Hyatt

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*******************************************************************************/

use std::os::raw::{c_char, c_double, c_int, c_uint, c_void};

use crate::types::{Mask, Training, MAX_LINE, MAX_MASKS, MAX_SEQ, MASK_SIZE, WINDOW};

extern "C" {
    fn seq_reader_gets(file: *mut c_void, buf: *mut c_char, len: c_int) -> *mut c_char;

    // bitmap functions from bitmap.rs
    fn test(bm: *const u8, ndx: c_int) -> u8;
    fn set(bm: *mut u8, ndx: c_int);
    fn toggle(bm: *mut u8, ndx: c_int);

    // C stderr stream
    static stderr: *mut libc::FILE;
}

/// Helper to get C stderr as *mut libc::FILE for fprintf
#[inline(always)]
unsafe fn c_stderr() -> *mut libc::FILE {
    stderr
}

/*******************************************************************************
  Read the sequence for training purposes.  If we encounter multiple
  sequences, we insert TTAATTAATTAA between each one to force stops in all
  six frames.  When we hit MAX_SEQ bp, we stop and return what we've got so
  far for training.  This routine reads in FASTA, and has a very 'loose'
  Genbank and Embl parser, but, to be safe, FASTA should generally be
  preferred.
*******************************************************************************/

#[no_mangle]
pub unsafe extern "C" fn read_seq_training(
    fp: *mut c_void,
    seq: *mut u8,
    useq: *mut u8,
    gc: *mut c_double,
    do_mask: c_int,
    mlist: *mut Mask,
    nm: *mut c_int,
) -> c_int {
    let mut line: [c_char; MAX_LINE + 1] = [0; MAX_LINE + 1];
    let mut hdr: c_int = 0;
    let mut fhdr: c_int = 0;
    let mut bctr: c_int = 0;
    let mut len: c_int = 0;
    let mut wrn: c_int = 0;
    let mut gc_cont: c_int = 0;
    let mut mask_beg: c_int = -1;
    let mut gapsize: c_uint = 0;

    line[MAX_LINE] = 0;
    while seq_reader_gets(fp, line.as_mut_ptr(), MAX_LINE as c_int) != std::ptr::null_mut() {
        if hdr == 0
            && *line.as_ptr().add(libc::strlen(line.as_ptr()) - 1) != b'\n' as c_char
            && wrn == 0
        {
            wrn = 1;
            libc::fprintf(
                c_stderr(),
                b"\n\nWarning: saw non-sequence line longer than \0".as_ptr() as *const c_char,
            );
            libc::fprintf(
                c_stderr(),
                b"%d chars, sequence might not be read \0".as_ptr() as *const c_char,
                MAX_LINE as c_int,
            );
            libc::fprintf(
                c_stderr(),
                b"correctly.\n\n\0".as_ptr() as *const c_char,
            );
        }
        if line[0] == b'>' as c_char
            || (line[0] == b'S' as c_char && line[1] == b'Q' as c_char)
            || (libc::strlen(line.as_ptr()) > 6
                && libc::strncmp(line.as_ptr(), b"ORIGIN\0".as_ptr() as *const c_char, 6) == 0)
        {
            hdr = 1;
            if fhdr > 0 {
                for i in 0u32..12 {
                    if i % 4 == 0 || i % 4 == 1 {
                        set(seq, bctr);
                        set(seq, bctr + 1);
                    }
                    bctr += 2;
                    len += 1;
                }
            }
            fhdr += 1;
        } else if hdr == 1 && (line[0] == b'/' as c_char && line[1] == b'/' as c_char) {
            hdr = 0;
        } else if hdr == 1 {
            if libc::strstr(line.as_ptr(), b"Expand\0".as_ptr() as *const c_char)
                != std::ptr::null_mut()
                && libc::strstr(line.as_ptr(), b"gap\0".as_ptr() as *const c_char)
                    != std::ptr::null_mut()
            {
                libc::sscanf(
                    libc::strstr(line.as_ptr(), b"gap\0".as_ptr() as *const c_char).add(4),
                    b"%u\0".as_ptr() as *const c_char,
                    &mut gapsize as *mut c_uint,
                );
                if gapsize < 1 || gapsize > MAX_LINE as c_uint {
                    libc::fprintf(
                        c_stderr(),
                        b"Error: gap size in gbk file can't exceed line\0".as_ptr()
                            as *const c_char,
                    );
                    libc::fprintf(
                        c_stderr(),
                        b" size.\n\0".as_ptr() as *const c_char,
                    );
                    libc::exit(51);
                }
                for i in 0..gapsize {
                    line[i as usize] = b'n' as c_char;
                }
                line[gapsize as usize] = 0;
            }
            let slen = libc::strlen(line.as_ptr());
            for i in 0..slen {
                if (line[i] as u8) < b'A' || (line[i] as u8) > b'z' {
                    continue;
                }
                if do_mask == 1
                    && mask_beg != -1
                    && line[i] != b'N' as c_char
                    && line[i] != b'n' as c_char
                {
                    if len - mask_beg >= MASK_SIZE as c_int {
                        if *nm == MAX_MASKS as c_int {
                            libc::fprintf(
                                c_stderr(),
                                b"Error: saw too many regions of 'N''s in the \0".as_ptr()
                                    as *const c_char,
                            );
                            libc::fprintf(
                                c_stderr(),
                                b"sequence.\n\0".as_ptr() as *const c_char,
                            );
                            libc::exit(52);
                        }
                        (*mlist.add(*nm as usize)).begin = mask_beg;
                        (*mlist.add(*nm as usize)).end = len - 1;
                        *nm += 1;
                    }
                    mask_beg = -1;
                }
                if do_mask == 1
                    && mask_beg == -1
                    && (line[i] == b'N' as c_char || line[i] == b'n' as c_char)
                {
                    mask_beg = len;
                }
                if line[i] == b'g' as c_char || line[i] == b'G' as c_char {
                    set(seq, bctr);
                    gc_cont += 1;
                } else if line[i] == b't' as c_char || line[i] == b'T' as c_char {
                    set(seq, bctr);
                    set(seq, bctr + 1);
                } else if line[i] == b'c' as c_char || line[i] == b'C' as c_char {
                    set(seq, bctr + 1);
                    gc_cont += 1;
                } else if line[i] != b'a' as c_char && line[i] != b'A' as c_char {
                    set(seq, bctr + 1);
                    set(useq, len);
                }
                bctr += 2;
                len += 1;
            }
        }
        if len + MAX_LINE as c_int >= MAX_SEQ as c_int {
            libc::fprintf(
                c_stderr(),
                b"\n\nWarning:  Sequence is long (max %d for training).\n\0".as_ptr()
                    as *const c_char,
                MAX_SEQ as c_int,
            );
            libc::fprintf(
                c_stderr(),
                b"Training on the first %d bases.\n\n\0".as_ptr() as *const c_char,
                MAX_SEQ as c_int,
            );
            break;
        }
    }
    if fhdr > 1 {
        for i in 0u32..12 {
            if i % 4 == 0 || i % 4 == 1 {
                set(seq, bctr);
                set(seq, bctr + 1);
            }
            bctr += 2;
            len += 1;
        }
    }
    *gc = gc_cont as c_double / len as c_double;
    len
}

/* This routine reads in the next sequence in a FASTA/GB/EMBL file */

#[no_mangle]
pub unsafe extern "C" fn next_seq_multi(
    fp: *mut c_void,
    seq: *mut u8,
    useq: *mut u8,
    sctr: *mut c_int,
    gc: *mut c_double,
    do_mask: c_int,
    mlist: *mut Mask,
    nm: *mut c_int,
    cur_hdr: *mut c_char,
    new_hdr: *mut c_char,
) -> c_int {
    let mut line: [c_char; MAX_LINE + 1] = [0; MAX_LINE + 1];
    let mut reading_seq: c_int = 0;
    let mut genbank_end: c_int = 0;
    let mut bctr: c_int = 0;
    let mut len: c_int = 0;
    let mut wrn: c_int = 0;
    let mut gc_cont: c_int = 0;
    let mut mask_beg: c_int = -1;
    let mut gapsize: c_uint = 0;

    libc::sprintf(
        new_hdr,
        b"Prodigal_Seq_%d\0".as_ptr() as *const c_char,
        *sctr + 2,
    );

    if *sctr > 0 {
        reading_seq = 1;
    }
    line[MAX_LINE] = 0;
    while seq_reader_gets(fp, line.as_mut_ptr(), MAX_LINE as c_int) != std::ptr::null_mut() {
        if reading_seq == 0
            && *line.as_ptr().add(libc::strlen(line.as_ptr()) - 1) != b'\n' as c_char
            && wrn == 0
        {
            wrn = 1;
            libc::fprintf(
                c_stderr(),
                b"\n\nWarning: saw non-sequence line longer than \0".as_ptr() as *const c_char,
            );
            libc::fprintf(
                c_stderr(),
                b"%d chars, sequence might not be read \0".as_ptr() as *const c_char,
                MAX_LINE as c_int,
            );
            libc::fprintf(
                c_stderr(),
                b"correctly.\n\n\0".as_ptr() as *const c_char,
            );
        }
        if libc::strlen(line.as_ptr()) > 10
            && libc::strncmp(line.as_ptr(), b"DEFINITION\0".as_ptr() as *const c_char, 10) == 0
        {
            if genbank_end == 0 {
                libc::strcpy(cur_hdr, line.as_ptr().add(12));
                *cur_hdr.add(libc::strlen(cur_hdr) - 1) = 0;
            } else {
                libc::strcpy(new_hdr, line.as_ptr().add(12));
                *new_hdr.add(libc::strlen(new_hdr) - 1) = 0;
            }
        }
        if line[0] == b'>' as c_char
            || (line[0] == b'S' as c_char && line[1] == b'Q' as c_char)
            || (libc::strlen(line.as_ptr()) > 6
                && libc::strncmp(line.as_ptr(), b"ORIGIN\0".as_ptr() as *const c_char, 6) == 0)
        {
            if reading_seq == 1 || genbank_end == 1 || *sctr > 0 {
                if line[0] == b'>' as c_char {
                    libc::strcpy(new_hdr, line.as_ptr().add(1));
                    *new_hdr.add(libc::strlen(new_hdr) - 1) = 0;
                }
                break;
            }
            if line[0] == b'>' as c_char {
                libc::strcpy(cur_hdr, line.as_ptr().add(1));
                *cur_hdr.add(libc::strlen(cur_hdr) - 1) = 0;
            }
            reading_seq = 1;
        } else if reading_seq == 1
            && (line[0] == b'/' as c_char && line[1] == b'/' as c_char)
        {
            reading_seq = 0;
            genbank_end = 1;
        } else if reading_seq == 1 {
            if libc::strstr(line.as_ptr(), b"Expand\0".as_ptr() as *const c_char)
                != std::ptr::null_mut()
                && libc::strstr(line.as_ptr(), b"gap\0".as_ptr() as *const c_char)
                    != std::ptr::null_mut()
            {
                libc::sscanf(
                    libc::strstr(line.as_ptr(), b"gap\0".as_ptr() as *const c_char).add(4),
                    b"%u\0".as_ptr() as *const c_char,
                    &mut gapsize as *mut c_uint,
                );
                if gapsize < 1 || gapsize > MAX_LINE as c_uint {
                    libc::fprintf(
                        c_stderr(),
                        b"Error: gap size in gbk file can't exceed line\0".as_ptr()
                            as *const c_char,
                    );
                    libc::fprintf(
                        c_stderr(),
                        b" size.\n\0".as_ptr() as *const c_char,
                    );
                    libc::exit(54);
                }
                for i in 0..gapsize {
                    line[i as usize] = b'n' as c_char;
                }
                line[gapsize as usize] = 0;
            }
            let slen = libc::strlen(line.as_ptr());
            for i in 0..slen {
                if (line[i] as u8) < b'A' || (line[i] as u8) > b'z' {
                    continue;
                }
                if do_mask == 1
                    && mask_beg != -1
                    && line[i] != b'N' as c_char
                    && line[i] != b'n' as c_char
                {
                    if len - mask_beg >= MASK_SIZE as c_int {
                        if *nm == MAX_MASKS as c_int {
                            libc::fprintf(
                                c_stderr(),
                                b"Error: saw too many regions of 'N''s in the \0".as_ptr()
                                    as *const c_char,
                            );
                            libc::fprintf(
                                c_stderr(),
                                b"sequence.\n\0".as_ptr() as *const c_char,
                            );
                            libc::exit(55);
                        }
                        (*mlist.add(*nm as usize)).begin = mask_beg;
                        (*mlist.add(*nm as usize)).end = len - 1;
                        *nm += 1;
                    }
                    mask_beg = -1;
                }
                if do_mask == 1
                    && mask_beg == -1
                    && (line[i] == b'N' as c_char || line[i] == b'n' as c_char)
                {
                    mask_beg = len;
                }
                if line[i] == b'g' as c_char || line[i] == b'G' as c_char {
                    set(seq, bctr);
                    gc_cont += 1;
                } else if line[i] == b't' as c_char || line[i] == b'T' as c_char {
                    set(seq, bctr);
                    set(seq, bctr + 1);
                } else if line[i] == b'c' as c_char || line[i] == b'C' as c_char {
                    set(seq, bctr + 1);
                    gc_cont += 1;
                } else if line[i] != b'a' as c_char && line[i] != b'A' as c_char {
                    set(seq, bctr + 1);
                    set(useq, len);
                }
                bctr += 2;
                len += 1;
            }
        }
        if len + MAX_LINE as c_int >= MAX_SEQ as c_int {
            libc::fprintf(
                c_stderr(),
                b"Sequence too long (max %d permitted).\n\0".as_ptr() as *const c_char,
                MAX_SEQ as c_int,
            );
            libc::exit(56);
        }
    }
    if len == 0 {
        return -1;
    }
    *gc = gc_cont as c_double / len as c_double;
    *sctr = *sctr + 1;
    len
}

/* Takes first word of header */
#[no_mangle]
pub unsafe extern "C" fn calc_short_header(
    header: *mut c_char,
    short_header: *mut c_char,
    sctr: c_int,
) {
    libc::strcpy(short_header, header);
    let hlen = libc::strlen(header);
    let mut i = 0usize;
    while i < hlen {
        if *header.add(i) == b' ' as c_char
            || *header.add(i) == b'\t' as c_char
            || *header.add(i) == b'\r' as c_char
            || *header.add(i) == b'\n' as c_char
        {
            libc::strncpy(short_header, header, i);
            *short_header.add(i) = 0;
            break;
        }
        i += 1;
    }
    if i == 0 {
        libc::sprintf(
            short_header,
            b"Prodigal_Seq_%d\0".as_ptr() as *const c_char,
            sctr,
        );
    }
}

/* Takes rseq and fills it up with the rev complement of seq */

#[no_mangle]
pub unsafe extern "C" fn rcom_seq(
    seq: *mut u8,
    rseq: *mut u8,
    useq: *mut u8,
    len: c_int,
) {
    let slen = len * 2;
    for i in 0..slen {
        if test(seq, i) == 0 {
            let offset = if i % 2 == 0 { -1 } else { 1 };
            set(rseq, slen - i - 1 + offset);
        }
    }
    for i in 0..len {
        if test(useq, i) == 1 {
            toggle(rseq, slen - 1 - i * 2);
            toggle(rseq, slen - 2 - i * 2);
        }
    }
}

/* Simple routines to say whether or not bases are */
/* a, c, t, g, starts, stops, etc. */

#[no_mangle]
pub unsafe extern "C" fn is_a(seq: *mut u8, n: c_int) -> c_int {
    let ndx = n * 2;
    if test(seq, ndx) == 1 || test(seq, ndx + 1) == 1 {
        return 0;
    }
    1
}

#[no_mangle]
pub unsafe extern "C" fn is_c(seq: *mut u8, n: c_int) -> c_int {
    let ndx = n * 2;
    if test(seq, ndx) == 1 || test(seq, ndx + 1) == 0 {
        return 0;
    }
    1
}

#[no_mangle]
pub unsafe extern "C" fn is_g(seq: *mut u8, n: c_int) -> c_int {
    let ndx = n * 2;
    if test(seq, ndx) == 0 || test(seq, ndx + 1) == 1 {
        return 0;
    }
    1
}

#[no_mangle]
pub unsafe extern "C" fn is_t(seq: *mut u8, n: c_int) -> c_int {
    let ndx = n * 2;
    if test(seq, ndx) == 0 || test(seq, ndx + 1) == 0 {
        return 0;
    }
    1
}

#[no_mangle]
pub unsafe extern "C" fn is_n(useq: *mut u8, n: c_int) -> c_int {
    if test(useq, n) == 0 {
        return 0;
    }
    1
}

#[no_mangle]
pub unsafe extern "C" fn is_stop(seq: *mut u8, n: c_int, tinf: *mut Training) -> c_int {
    let tt = (*tinf).trans_table;

    /* TAG */
    if is_t(seq, n) == 1 && is_a(seq, n + 1) == 1 && is_g(seq, n + 2) == 1 {
        if tt == 6 || tt == 15 || tt == 16 || tt == 22 {
            return 0;
        }
        return 1;
    }

    /* TGA */
    if is_t(seq, n) == 1 && is_g(seq, n + 1) == 1 && is_a(seq, n + 2) == 1 {
        if (tt >= 2 && tt <= 5)
            || tt == 9
            || tt == 10
            || tt == 13
            || tt == 14
            || tt == 21
            || tt == 25
        {
            return 0;
        }
        return 1;
    }

    /* TAA */
    if is_t(seq, n) == 1 && is_a(seq, n + 1) == 1 && is_a(seq, n + 2) == 1 {
        if tt == 6 || tt == 14 {
            return 0;
        }
        return 1;
    }

    /* Code 2 */
    if tt == 2 && is_a(seq, n) == 1 && is_g(seq, n + 1) == 1 && is_a(seq, n + 2) == 1 {
        return 1;
    }
    if tt == 2 && is_a(seq, n) == 1 && is_g(seq, n + 1) == 1 && is_g(seq, n + 2) == 1 {
        return 1;
    }

    /* Code 22 */
    if tt == 22 && is_t(seq, n) == 1 && is_c(seq, n + 1) == 1 && is_a(seq, n + 2) == 1 {
        return 1;
    }

    /* Code 23 */
    if tt == 23 && is_t(seq, n) == 1 && is_t(seq, n + 1) == 1 && is_a(seq, n + 2) == 1 {
        return 1;
    }

    0
}

#[no_mangle]
pub unsafe extern "C" fn is_start(seq: *mut u8, n: c_int, tinf: *mut Training) -> c_int {
    let tt = (*tinf).trans_table;

    /* ATG */
    if is_a(seq, n) == 1 && is_t(seq, n + 1) == 1 && is_g(seq, n + 2) == 1 {
        return 1;
    }

    /* Codes that only use ATG */
    if tt == 6 || tt == 10 || tt == 14 || tt == 15 || tt == 16 || tt == 22 {
        return 0;
    }

    /* GTG */
    if is_g(seq, n) == 1 && is_t(seq, n + 1) == 1 && is_g(seq, n + 2) == 1 {
        if tt == 1 || tt == 3 || tt == 12 || tt == 22 {
            return 0;
        }
        return 1;
    }

    /* TTG */
    if is_t(seq, n) == 1 && is_t(seq, n + 1) == 1 && is_g(seq, n + 2) == 1 {
        if tt < 4 || tt == 9 || (tt >= 21 && tt < 25) {
            return 0;
        }
        return 1;
    }

    /* We do not handle other initiation codons */
    0
}

#[no_mangle]
pub unsafe extern "C" fn is_atg(seq: *mut u8, n: c_int) -> c_int {
    if is_a(seq, n) == 0 || is_t(seq, n + 1) == 0 || is_g(seq, n + 2) == 0 {
        return 0;
    }
    1
}

#[no_mangle]
pub unsafe extern "C" fn is_gtg(seq: *mut u8, n: c_int) -> c_int {
    if is_g(seq, n) == 0 || is_t(seq, n + 1) == 0 || is_g(seq, n + 2) == 0 {
        return 0;
    }
    1
}

#[no_mangle]
pub unsafe extern "C" fn is_ttg(seq: *mut u8, n: c_int) -> c_int {
    if is_t(seq, n) == 0 || is_t(seq, n + 1) == 0 || is_g(seq, n + 2) == 0 {
        return 0;
    }
    1
}

#[no_mangle]
pub unsafe extern "C" fn is_gc(seq: *mut u8, n: c_int) -> c_int {
    let ndx = n * 2;
    if test(seq, ndx) != test(seq, ndx + 1) {
        return 1;
    }
    0
}

#[no_mangle]
pub unsafe extern "C" fn gc_content(seq: *mut u8, a: c_int, b: c_int) -> c_double {
    let mut sum: c_double = 0.0;
    let mut gc: c_double = 0.0;
    for i in a..=b {
        if is_g(seq, i) == 1 || is_c(seq, i) == 1 {
            gc += 1.0;
        }
        sum += 1.0;
    }
    gc / sum
}

/* Returns a single amino acid for this position */
#[no_mangle]
pub unsafe extern "C" fn amino(
    seq: *mut u8,
    n: c_int,
    tinf: *mut Training,
    is_init: c_int,
) -> c_char {
    let tt = (*tinf).trans_table;

    if is_stop(seq, n, tinf) == 1 {
        return b'*' as c_char;
    }
    if is_start(seq, n, tinf) == 1 && is_init == 1 {
        return b'M' as c_char;
    }
    if is_t(seq, n) == 1 && is_t(seq, n + 1) == 1 && is_t(seq, n + 2) == 1 {
        return b'F' as c_char;
    }
    if is_t(seq, n) == 1 && is_t(seq, n + 1) == 1 && is_c(seq, n + 2) == 1 {
        return b'F' as c_char;
    }
    if is_t(seq, n) == 1 && is_t(seq, n + 1) == 1 && is_a(seq, n + 2) == 1 {
        return b'L' as c_char;
    }
    if is_t(seq, n) == 1 && is_t(seq, n + 1) == 1 && is_g(seq, n + 2) == 1 {
        return b'L' as c_char;
    }
    if is_t(seq, n) == 1 && is_c(seq, n + 1) == 1 {
        return b'S' as c_char;
    }
    if is_t(seq, n) == 1 && is_a(seq, n + 1) == 1 && is_t(seq, n + 2) == 1 {
        return b'Y' as c_char;
    }
    if is_t(seq, n) == 1 && is_a(seq, n + 1) == 1 && is_c(seq, n + 2) == 1 {
        return b'Y' as c_char;
    }
    if is_t(seq, n) == 1 && is_a(seq, n + 1) == 1 && is_a(seq, n + 2) == 1 {
        if tt == 6 {
            return b'Q' as c_char;
        }
        if tt == 14 {
            return b'Y' as c_char;
        }
    }
    if is_t(seq, n) == 1 && is_a(seq, n + 1) == 1 && is_g(seq, n + 2) == 1 {
        if tt == 6 || tt == 15 {
            return b'Q' as c_char;
        }
        if tt == 22 {
            return b'L' as c_char;
        }
    }
    if is_t(seq, n) == 1 && is_g(seq, n + 1) == 1 && is_t(seq, n + 2) == 1 {
        return b'C' as c_char;
    }
    if is_t(seq, n) == 1 && is_g(seq, n + 1) == 1 && is_c(seq, n + 2) == 1 {
        return b'C' as c_char;
    }
    if is_t(seq, n) == 1 && is_g(seq, n + 1) == 1 && is_a(seq, n + 2) == 1 {
        if tt == 25 {
            return b'G' as c_char;
        } else {
            return b'W' as c_char;
        }
    }
    if is_t(seq, n) == 1 && is_g(seq, n + 1) == 1 && is_g(seq, n + 2) == 1 {
        return b'W' as c_char;
    }
    if is_c(seq, n) == 1 && is_t(seq, n + 1) == 1 && is_t(seq, n + 2) == 1 {
        if tt == 3 {
            return b'T' as c_char;
        }
        return b'L' as c_char;
    }
    if is_c(seq, n) == 1 && is_t(seq, n + 1) == 1 && is_c(seq, n + 2) == 1 {
        if tt == 3 {
            return b'T' as c_char;
        }
        return b'L' as c_char;
    }
    if is_c(seq, n) == 1 && is_t(seq, n + 1) == 1 && is_a(seq, n + 2) == 1 {
        if tt == 3 {
            return b'T' as c_char;
        }
        return b'L' as c_char;
    }
    if is_c(seq, n) == 1 && is_t(seq, n + 1) == 1 && is_g(seq, n + 2) == 1 {
        if tt == 3 {
            return b'T' as c_char;
        }
        if tt == 12 {
            return b'S' as c_char;
        }
        return b'L' as c_char;
    }
    if is_c(seq, n) == 1 && is_c(seq, n + 1) == 1 {
        return b'P' as c_char;
    }
    if is_c(seq, n) == 1 && is_a(seq, n + 1) == 1 && is_t(seq, n + 2) == 1 {
        return b'H' as c_char;
    }
    if is_c(seq, n) == 1 && is_a(seq, n + 1) == 1 && is_c(seq, n + 2) == 1 {
        return b'H' as c_char;
    }
    if is_c(seq, n) == 1 && is_a(seq, n + 1) == 1 && is_a(seq, n + 2) == 1 {
        return b'Q' as c_char;
    }
    if is_c(seq, n) == 1 && is_a(seq, n + 1) == 1 && is_g(seq, n + 2) == 1 {
        return b'Q' as c_char;
    }
    if is_c(seq, n) == 1 && is_g(seq, n + 1) == 1 {
        return b'R' as c_char;
    }
    if is_a(seq, n) == 1 && is_t(seq, n + 1) == 1 && is_t(seq, n + 2) == 1 {
        return b'I' as c_char;
    }
    if is_a(seq, n) == 1 && is_t(seq, n + 1) == 1 && is_c(seq, n + 2) == 1 {
        return b'I' as c_char;
    }
    if is_a(seq, n) == 1 && is_t(seq, n + 1) == 1 && is_a(seq, n + 2) == 1 {
        if tt == 2 || tt == 3 || tt == 5 || tt == 13 || tt == 21 {
            return b'M' as c_char;
        }
        return b'I' as c_char;
    }
    if is_a(seq, n) == 1 && is_t(seq, n + 1) == 1 && is_g(seq, n + 2) == 1 {
        return b'M' as c_char;
    }
    if is_a(seq, n) == 1 && is_c(seq, n + 1) == 1 {
        return b'T' as c_char;
    }
    if is_a(seq, n) == 1 && is_a(seq, n + 1) == 1 && is_t(seq, n + 2) == 1 {
        return b'N' as c_char;
    }
    if is_a(seq, n) == 1 && is_a(seq, n + 1) == 1 && is_c(seq, n + 2) == 1 {
        return b'N' as c_char;
    }
    if is_a(seq, n) == 1 && is_a(seq, n + 1) == 1 && is_a(seq, n + 2) == 1 {
        if tt == 9 || tt == 14 || tt == 21 {
            return b'N' as c_char;
        }
        return b'K' as c_char;
    }
    if is_a(seq, n) == 1 && is_a(seq, n + 1) == 1 && is_g(seq, n + 2) == 1 {
        return b'K' as c_char;
    }
    if is_a(seq, n) == 1 && is_g(seq, n + 1) == 1 && is_t(seq, n + 2) == 1 {
        return b'S' as c_char;
    }
    if is_a(seq, n) == 1 && is_g(seq, n + 1) == 1 && is_c(seq, n + 2) == 1 {
        return b'S' as c_char;
    }
    if is_a(seq, n) == 1
        && is_g(seq, n + 1) == 1
        && (is_a(seq, n + 2) == 1 || is_g(seq, n + 2) == 1)
    {
        if tt == 13 {
            return b'G' as c_char;
        }
        if tt == 5 || tt == 9 || tt == 14 || tt == 21 {
            return b'S' as c_char;
        }
        return b'R' as c_char;
    }
    if is_g(seq, n) == 1 && is_t(seq, n + 1) == 1 {
        return b'V' as c_char;
    }
    if is_g(seq, n) == 1 && is_c(seq, n + 1) == 1 {
        return b'A' as c_char;
    }
    if is_g(seq, n) == 1 && is_a(seq, n + 1) == 1 && is_t(seq, n + 2) == 1 {
        return b'D' as c_char;
    }
    if is_g(seq, n) == 1 && is_a(seq, n + 1) == 1 && is_c(seq, n + 2) == 1 {
        return b'D' as c_char;
    }
    if is_g(seq, n) == 1 && is_a(seq, n + 1) == 1 && is_a(seq, n + 2) == 1 {
        return b'E' as c_char;
    }
    if is_g(seq, n) == 1 && is_a(seq, n + 1) == 1 && is_g(seq, n + 2) == 1 {
        return b'E' as c_char;
    }
    if is_g(seq, n) == 1 && is_g(seq, n + 1) == 1 {
        return b'G' as c_char;
    }
    b'X' as c_char
}

/* Converts an amino acid letter to a numerical value */
#[no_mangle]
pub unsafe extern "C" fn amino_num(aa: c_char) -> c_int {
    let c = aa as u8;
    if c == b'a' || c == b'A' { return 0; }
    if c == b'c' || c == b'C' { return 1; }
    if c == b'd' || c == b'D' { return 2; }
    if c == b'e' || c == b'E' { return 3; }
    if c == b'f' || c == b'F' { return 4; }
    if c == b'g' || c == b'G' { return 5; }
    if c == b'h' || c == b'H' { return 6; }
    if c == b'i' || c == b'I' { return 7; }
    if c == b'k' || c == b'K' { return 8; }
    if c == b'l' || c == b'L' { return 9; }
    if c == b'm' || c == b'M' { return 10; }
    if c == b'n' || c == b'N' { return 11; }
    if c == b'p' || c == b'P' { return 12; }
    if c == b'q' || c == b'Q' { return 13; }
    if c == b'r' || c == b'R' { return 14; }
    if c == b's' || c == b'S' { return 15; }
    if c == b't' || c == b'T' { return 16; }
    if c == b'v' || c == b'V' { return 17; }
    if c == b'w' || c == b'W' { return 18; }
    if c == b'y' || c == b'Y' { return 19; }
    -1
}

/* Converts a numerical value to an amino acid letter */
#[no_mangle]
pub unsafe extern "C" fn amino_letter(num: c_int) -> c_char {
    if num == 0 { return b'A' as c_char; }
    if num == 1 { return b'C' as c_char; }
    if num == 2 { return b'D' as c_char; }
    if num == 3 { return b'E' as c_char; }
    if num == 4 { return b'F' as c_char; }
    if num == 5 { return b'G' as c_char; }
    if num == 6 { return b'H' as c_char; }
    if num == 7 { return b'I' as c_char; }
    if num == 8 { return b'K' as c_char; }
    if num == 9 { return b'L' as c_char; }
    if num == 10 { return b'M' as c_char; }
    if num == 11 { return b'N' as c_char; }
    if num == 12 { return b'P' as c_char; }
    if num == 13 { return b'Q' as c_char; }
    if num == 14 { return b'R' as c_char; }
    if num == 15 { return b'S' as c_char; }
    if num == 16 { return b'T' as c_char; }
    if num == 17 { return b'V' as c_char; }
    if num == 18 { return b'W' as c_char; }
    if num == 19 { return b'Y' as c_char; }
    b'X' as c_char
}

/* Returns the corresponding frame on the reverse strand */

#[no_mangle]
pub unsafe extern "C" fn rframe(fr: c_int, slen: c_int) -> c_int {
    let mut md = slen % 3 - 1;
    if md == 0 {
        md = 3;
    }
    md - fr
}

/* Simple 3-way max function */

#[no_mangle]
pub unsafe extern "C" fn max_fr(n1: c_int, n2: c_int, n3: c_int) -> c_int {
    if n1 > n2 {
        if n1 > n3 { 0 } else { 2 }
    } else {
        if n2 > n3 { 1 } else { 2 }
    }
}

/*******************************************************************************
  Creates a GC frame plot for a given sequence.  This is simply a string with
  the highest GC content frame for a window centered on position for every
  position in the sequence.
*******************************************************************************/

#[no_mangle]
pub unsafe extern "C" fn calc_most_gc_frame(seq: *mut u8, slen: c_int) -> *mut c_int {
    let gp: *mut c_int = libc::malloc((slen as usize) * std::mem::size_of::<c_double>()) as *mut c_int;
    let fwd: *mut c_int = libc::malloc((slen as usize) * std::mem::size_of::<c_int>()) as *mut c_int;
    let bwd: *mut c_int = libc::malloc((slen as usize) * std::mem::size_of::<c_int>()) as *mut c_int;
    let tot: *mut c_int = libc::malloc((slen as usize) * std::mem::size_of::<c_int>()) as *mut c_int;

    if fwd.is_null() || bwd.is_null() || gp.is_null() || tot.is_null() {
        return std::ptr::null_mut();
    }

    for i in 0..slen {
        *fwd.add(i as usize) = 0;
        *bwd.add(i as usize) = 0;
        *tot.add(i as usize) = 0;
        *gp.add(i as usize) = -1;
    }

    for i in 0..3 {
        let mut j = i;
        while j < slen {
            if j < 3 {
                *fwd.add(j as usize) = is_gc(seq, j);
            } else {
                *fwd.add(j as usize) = *fwd.add((j - 3) as usize) + is_gc(seq, j);
            }
            if j < 3 {
                *bwd.add((slen - j - 1) as usize) = is_gc(seq, slen - j - 1);
            } else {
                *bwd.add((slen - j - 1) as usize) =
                    *bwd.add((slen - j + 2) as usize) + is_gc(seq, slen - j - 1);
            }
            j += 1;
        }
    }

    for i in 0..slen {
        *tot.add(i as usize) = *fwd.add(i as usize) + *bwd.add(i as usize) - is_gc(seq, i);
        if i - WINDOW as c_int / 2 >= 0 {
            *tot.add(i as usize) -= *fwd.add((i - WINDOW as c_int / 2) as usize);
        }
        if i + WINDOW as c_int / 2 < slen {
            *tot.add(i as usize) -= *bwd.add((i + WINDOW as c_int / 2) as usize);
        }
    }
    libc::free(fwd as *mut c_void);
    libc::free(bwd as *mut c_void);

    let mut i = 0;
    while i < slen - 2 {
        let win = max_fr(
            *tot.add(i as usize),
            *tot.add((i + 1) as usize),
            *tot.add((i + 2) as usize),
        );
        for j in 0..3 {
            *gp.add((i + j) as usize) = win;
        }
        i += 3;
    }
    libc::free(tot as *mut c_void);
    gp
}

/* Converts a word of size len to a number */
#[no_mangle]
pub unsafe extern "C" fn mer_ndx(len: c_int, seq: *mut u8, pos: c_int) -> c_int {
    let mut ndx: c_int = 0;
    for i in 0..(2 * len) {
        ndx |= (test(seq, pos * 2 + i) as c_int) << i;
    }
    ndx
}

/* Gives a text string for a start */
#[no_mangle]
pub unsafe extern "C" fn start_text(st: *mut c_char, type_: c_int) {
    if type_ == 0 {
        *st.add(0) = b'A' as c_char;
    } else if type_ == 1 {
        *st.add(0) = b'G' as c_char;
    } else if type_ == 2 {
        *st.add(0) = b'T' as c_char;
    }
    *st.add(1) = b'T' as c_char;
    *st.add(2) = b'G' as c_char;
    *st.add(3) = 0;
}

/* Gives a text string for a mer of size 'len' (useful for outputting motifs) */
#[no_mangle]
pub unsafe extern "C" fn mer_text(qt: *mut c_char, len: c_int, ndx: c_int) {
    let letters: [c_char; 4] = [
        b'A' as c_char,
        b'G' as c_char,
        b'C' as c_char,
        b'T' as c_char,
    ];
    if len == 0 {
        libc::strcpy(qt, b"None\0".as_ptr() as *const c_char);
    } else {
        for i in 0..len {
            let mut val = (ndx & (1 << (2 * i))) + (ndx & (1 << (2 * i + 1)));
            val >>= i * 2;
            *qt.add(i as usize) = letters[val as usize];
        }
        *qt.add(len as usize) = 0;
    }
}

/* Builds a 'len'-mer background for whole sequence */
#[no_mangle]
pub unsafe extern "C" fn calc_mer_bg(
    len: c_int,
    seq: *mut u8,
    rseq: *mut u8,
    slen: c_int,
    bg: *mut c_double,
) {
    let mut glob: c_int = 0;
    let mut size: c_int = 1;

    for _i in 0..len {
        size *= 4;
    }
    let counts: *mut c_int = libc::malloc((size as usize) * std::mem::size_of::<c_int>()) as *mut c_int;
    for i in 0..size {
        *counts.add(i as usize) = 0;
    }
    for i in 0..(slen - len + 1) {
        *counts.add(mer_ndx(len, seq, i) as usize) += 1;
        *counts.add(mer_ndx(len, rseq, i) as usize) += 1;
        glob += 2;
    }
    for i in 0..size {
        *bg.add(i as usize) = (*counts.add(i as usize) as c_double * 1.0) / (glob as c_double * 1.0);
    }
    libc::free(counts as *mut c_void);
}

/*******************************************************************************
  Finds the highest-scoring region similar to AGGAGG in a given stretch of
  sequence upstream of a start.
*******************************************************************************/

#[no_mangle]
pub unsafe extern "C" fn shine_dalgarno_exact(
    seq: *mut u8,
    pos: c_int,
    start: c_int,
    rwt: *mut c_double,
) -> c_int {
    let mut cur_val: c_int = 0;
    let mut match_: [c_double; 6] = [0.0; 6];
    let mut cur_ctr: c_double;
    let mut dis_flag: c_double;

    let limit = imin(6, start - 4 - pos);
    for i in 0..6 {
        match_[i] = -10.0;
    }

    /* Compare the 6-base region to AGGAGG */
    for i in 0..limit {
        if pos + i >= 0 {
            if i % 3 == 0 && is_a(seq, pos + i) == 1 {
                match_[i as usize] = 2.0;
            } else if i % 3 != 0 && is_g(seq, pos + i) == 1 {
                match_[i as usize] = 3.0;
            }
        }
    }

    /* Find the maximally scoring motif */
    let mut max_val: c_int = 0;
    let mut i = limit;
    while i >= 3 {
        for j in 0..=(limit - i) {
            cur_ctr = -2.0;
            let mut mism: c_int = 0;
            for k in j..(j + i) {
                cur_ctr += match_[k as usize];
                if match_[k as usize] < 0.0 {
                    mism += 1;
                }
            }
            if mism > 0 {
                continue;
            }
            let rdis = start - (pos + j + i);
            if rdis < 5 && i < 5 {
                dis_flag = 2.0;
            } else if rdis < 5 && i >= 5 {
                dis_flag = 1.0;
            } else if rdis > 10 && rdis <= 12 && i < 5 {
                dis_flag = 1.0;
            } else if rdis > 10 && rdis <= 12 && i >= 5 {
                dis_flag = 2.0;
            } else if rdis >= 13 {
                dis_flag = 3.0;
            } else {
                dis_flag = 0.0;
            }
            if rdis > 15 || cur_ctr < 6.0 {
                continue;
            }

            /* Exact-Matching RBS Motifs */
            if cur_ctr < 6.0 {
                cur_val = 0;
            } else if cur_ctr == 6.0 && dis_flag == 2.0 {
                cur_val = 1;
            } else if cur_ctr == 6.0 && dis_flag == 3.0 {
                cur_val = 2;
            } else if cur_ctr == 8.0 && dis_flag == 3.0 {
                cur_val = 3;
            } else if cur_ctr == 9.0 && dis_flag == 3.0 {
                cur_val = 3;
            } else if cur_ctr == 6.0 && dis_flag == 1.0 {
                cur_val = 6;
            } else if cur_ctr == 11.0 && dis_flag == 3.0 {
                cur_val = 10;
            } else if cur_ctr == 12.0 && dis_flag == 3.0 {
                cur_val = 10;
            } else if cur_ctr == 14.0 && dis_flag == 3.0 {
                cur_val = 10;
            } else if cur_ctr == 8.0 && dis_flag == 2.0 {
                cur_val = 11;
            } else if cur_ctr == 9.0 && dis_flag == 2.0 {
                cur_val = 11;
            } else if cur_ctr == 8.0 && dis_flag == 1.0 {
                cur_val = 12;
            } else if cur_ctr == 9.0 && dis_flag == 1.0 {
                cur_val = 12;
            } else if cur_ctr == 6.0 && dis_flag == 0.0 {
                cur_val = 13;
            } else if cur_ctr == 8.0 && dis_flag == 0.0 {
                cur_val = 15;
            } else if cur_ctr == 9.0 && dis_flag == 0.0 {
                cur_val = 16;
            } else if cur_ctr == 11.0 && dis_flag == 2.0 {
                cur_val = 20;
            } else if cur_ctr == 11.0 && dis_flag == 1.0 {
                cur_val = 21;
            } else if cur_ctr == 11.0 && dis_flag == 0.0 {
                cur_val = 22;
            } else if cur_ctr == 12.0 && dis_flag == 2.0 {
                cur_val = 20;
            } else if cur_ctr == 12.0 && dis_flag == 1.0 {
                cur_val = 23;
            } else if cur_ctr == 12.0 && dis_flag == 0.0 {
                cur_val = 24;
            } else if cur_ctr == 14.0 && dis_flag == 2.0 {
                cur_val = 25;
            } else if cur_ctr == 14.0 && dis_flag == 1.0 {
                cur_val = 26;
            } else if cur_ctr == 14.0 && dis_flag == 0.0 {
                cur_val = 27;
            }

            if *rwt.add(cur_val as usize) < *rwt.add(max_val as usize) {
                continue;
            }
            if *rwt.add(cur_val as usize) == *rwt.add(max_val as usize) && cur_val < max_val {
                continue;
            }
            max_val = cur_val;
        }
        i -= 1;
    }

    max_val
}

/*******************************************************************************
  Finds the highest-scoring region similar to AGGAGG in a given stretch of
  sequence upstream of a start.  Only considers 5/6-mers with 1 mismatch.
*******************************************************************************/

#[no_mangle]
pub unsafe extern "C" fn shine_dalgarno_mm(
    seq: *mut u8,
    pos: c_int,
    start: c_int,
    rwt: *mut c_double,
) -> c_int {
    let mut cur_val: c_int = 0;
    let mut match_: [c_double; 6] = [0.0; 6];
    let mut cur_ctr: c_double;
    let mut dis_flag: c_double;

    let limit = imin(6, start - 4 - pos);
    for i in 0..6 {
        match_[i] = -10.0;
    }

    /* Compare the 6-base region to AGGAGG */
    for i in 0..limit {
        if pos + i >= 0 {
            if i % 3 == 0 {
                if is_a(seq, pos + i) == 1 {
                    match_[i as usize] = 2.0;
                } else {
                    match_[i as usize] = -3.0;
                }
            } else {
                if is_g(seq, pos + i) == 1 {
                    match_[i as usize] = 3.0;
                } else {
                    match_[i as usize] = -2.0;
                }
            }
        }
    }

    /* Find the maximally scoring motif */
    let mut max_val: c_int = 0;
    let mut i = limit;
    while i >= 5 {
        for j in 0..=(limit - i) {
            cur_ctr = -2.0;
            let mut mism: c_int = 0;
            for k in j..(j + i) {
                cur_ctr += match_[k as usize];
                if match_[k as usize] < 0.0 {
                    mism += 1;
                }
                if match_[k as usize] < 0.0 && (k <= j + 1 || k >= j + i - 2) {
                    cur_ctr -= 10.0;
                }
            }
            if mism != 1 {
                continue;
            }
            let rdis = start - (pos + j + i);
            if rdis < 5 {
                dis_flag = 1.0;
            } else if rdis > 10 && rdis <= 12 {
                dis_flag = 2.0;
            } else if rdis >= 13 {
                dis_flag = 3.0;
            } else {
                dis_flag = 0.0;
            }
            if rdis > 15 || cur_ctr < 6.0 {
                continue;
            }

            /* Single-Mismatch RBS Motifs */
            if cur_ctr < 6.0 {
                cur_val = 0;
            } else if cur_ctr == 6.0 && dis_flag == 3.0 {
                cur_val = 2;
            } else if cur_ctr == 7.0 && dis_flag == 3.0 {
                cur_val = 2;
            } else if cur_ctr == 9.0 && dis_flag == 3.0 {
                cur_val = 3;
            } else if cur_ctr == 6.0 && dis_flag == 2.0 {
                cur_val = 4;
            } else if cur_ctr == 6.0 && dis_flag == 1.0 {
                cur_val = 5;
            } else if cur_ctr == 6.0 && dis_flag == 0.0 {
                cur_val = 9;
            } else if cur_ctr == 7.0 && dis_flag == 2.0 {
                cur_val = 7;
            } else if cur_ctr == 7.0 && dis_flag == 1.0 {
                cur_val = 8;
            } else if cur_ctr == 7.0 && dis_flag == 0.0 {
                cur_val = 14;
            } else if cur_ctr == 9.0 && dis_flag == 2.0 {
                cur_val = 17;
            } else if cur_ctr == 9.0 && dis_flag == 1.0 {
                cur_val = 18;
            } else if cur_ctr == 9.0 && dis_flag == 0.0 {
                cur_val = 19;
            }

            if *rwt.add(cur_val as usize) < *rwt.add(max_val as usize) {
                continue;
            }
            if *rwt.add(cur_val as usize) == *rwt.add(max_val as usize) && cur_val < max_val {
                continue;
            }
            max_val = cur_val;
        }
        i -= 1;
    }

    max_val
}

/* Returns the minimum of two numbers */
#[no_mangle]
pub unsafe extern "C" fn imin(x: c_int, y: c_int) -> c_int {
    if x < y {
        x
    } else {
        y
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_imin() {
        unsafe {
            assert_eq!(imin(3, 5), 3);
            assert_eq!(imin(5, 3), 3);
            assert_eq!(imin(-1, 0), -1);
        }
    }

    #[test]
    fn test_amino_num_letter_roundtrip() {
        unsafe {
            let letters = b"ACDEFGHIKLMNPQRSTVWY";
            for (i, &c) in letters.iter().enumerate() {
                assert_eq!(amino_num(c as c_char), i as c_int);
                assert_eq!(amino_letter(i as c_int) as u8, c);
            }
            assert_eq!(amino_num(b'Z' as c_char), -1);
            assert_eq!(amino_letter(20) as u8, b'X');
        }
    }

    #[test]
    fn test_rframe() {
        unsafe {
            // slen=10: 10%3=1, md=1-1=0, md becomes 3 => 3-fr
            assert_eq!(rframe(0, 10), 3);
            assert_eq!(rframe(1, 10), 2);
        }
    }

    #[test]
    fn test_max_fr() {
        unsafe {
            assert_eq!(max_fr(3, 2, 1), 0);
            assert_eq!(max_fr(1, 3, 2), 1);
            assert_eq!(max_fr(1, 2, 3), 2);
            assert_eq!(max_fr(1, 1, 1), 2); // ties go to 2
        }
    }
}
