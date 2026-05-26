use std::fs::File;
use std::io::{self, BufRead, BufReader, Cursor, Read, Seek, SeekFrom};
use std::os::raw::{c_char, c_int, c_void};

use flate2::read::GzDecoder;

/// Opaque sequence reader that handles both plain and gzip-compressed files.
/// Passed through the pipeline as `*mut c_void`.
pub enum SeqReader {
    Plain(BufReader<File>),
    Gzip {
        path: String,
        inner: BufReader<GzDecoder<File>>,
    },
    Memory(BufReader<Cursor<Vec<u8>>>),
}

impl SeqReader {
    /// Open a file, auto-detecting gzip by magic bytes.
    pub fn open(path: &str) -> io::Result<Self> {
        let mut file = File::open(path)?;

        // Check gzip magic bytes (0x1f 0x8b)
        let mut magic = [0u8; 2];
        let n = file.read(&mut magic)?;
        file.seek(SeekFrom::Start(0))?;

        if n >= 2 && magic[0] == 0x1f && magic[1] == 0x8b {
            Ok(SeqReader::Gzip {
                path: path.to_string(),
                inner: BufReader::new(GzDecoder::new(file)),
            })
        } else {
            Ok(SeqReader::Plain(BufReader::new(file)))
        }
    }

    /// Slurp all of stdin into memory and wrap it in a seekable reader so the
    /// pipeline's rewind-and-rescan logic works on piped input.
    pub fn stdin() -> io::Result<Self> {
        let mut data = Vec::new();
        io::stdin().read_to_end(&mut data)?;
        Ok(SeqReader::Memory(BufReader::new(Cursor::new(data))))
    }

    fn inner_mut(&mut self) -> &mut dyn BufRead {
        match self {
            SeqReader::Plain(inner) => inner,
            SeqReader::Gzip { inner, .. } => inner,
            SeqReader::Memory(inner) => inner,
        }
    }

    /// Read a line (up to `max_len - 1` bytes) into the provided C buffer.
    /// Returns true if a line was read, false on EOF.
    /// Mimics gzgets: fills buf with up to max_len-1 chars + NUL, stops at newline.
    pub fn gets(&mut self, buf: *mut c_char, max_len: c_int) -> bool {
        let max = max_len as usize;
        if max == 0 {
            return false;
        }

        // Clear the buffer
        unsafe {
            *buf = 0;
        }

        let mut total = 0usize;

        // Copy directly from the buffered reader until newline or max-1.
        loop {
            let inner = self.inner_mut();
            let available = inner.fill_buf().unwrap_or(&[]);
            if available.is_empty() {
                break; // EOF
            }

            // Find newline or end of available data
            let limit = (max - 1 - total).min(available.len());
            let mut consumed = 0;
            for i in 0..limit {
                unsafe {
                    *(buf as *mut u8).add(total) = available[i];
                }
                consumed = i + 1;
                total += 1;
                if available[i] == b'\n' {
                    break;
                }
            }
            inner.consume(consumed);

            if consumed == 0 {
                break;
            }
            let last = unsafe { *(buf as *const u8).add(total - 1) };
            if total >= max - 1 || last == b'\n' {
                break;
            }
        }

        if total == 0 {
            return false;
        }

        // NUL terminator
        unsafe {
            *buf.add(total) = 0;
        }
        true
    }

    /// Seek to beginning of the stream.
    pub fn rewind(&mut self) -> bool {
        match self {
            SeqReader::Plain(inner) => inner.seek(SeekFrom::Start(0)).is_ok(),
            SeqReader::Gzip { path, inner } => match File::open(path) {
                Ok(file) => {
                    *inner = BufReader::new(GzDecoder::new(file));
                    true
                }
                Err(_) => false,
            },
            SeqReader::Memory(inner) => inner.seek(SeekFrom::Start(0)).is_ok(),
        }
    }
}

// ---------------------------------------------------------------------------
// C-compatible wrapper functions (replace gzopen/gzgets/gzseek/gzclose)
// ---------------------------------------------------------------------------

/// Open a file for reading. Returns an opaque pointer (like gzopen).
/// The caller must eventually call `seq_reader_close`.
pub unsafe fn seq_reader_open(path: *const c_char) -> *mut c_void {
    let c_str = std::ffi::CStr::from_ptr(path);
    let path_str = match c_str.to_str() {
        Ok(s) => s,
        Err(_) => return std::ptr::null_mut(),
    };
    match SeqReader::open(path_str) {
        Ok(reader) => Box::into_raw(Box::new(reader)) as *mut c_void,
        Err(_) => std::ptr::null_mut(),
    }
}

/// Read all stdin into a seekable in-memory reader.
pub fn seq_reader_open_stdin() -> *mut c_void {
    match SeqReader::stdin() {
        Ok(reader) => Box::into_raw(Box::new(reader)) as *mut c_void,
        Err(_) => std::ptr::null_mut(),
    }
}

/// Read a line from the reader (like gzgets). Returns null on EOF.
pub unsafe fn seq_reader_gets(
    handle: *mut c_void,
    buf: *mut c_char,
    max_len: c_int,
) -> *mut c_char {
    let reader = &mut *(handle as *mut SeqReader);
    if reader.gets(buf, max_len) {
        buf
    } else {
        std::ptr::null_mut()
    }
}

/// Seek to the beginning (like gzseek(fp, 0, SEEK_SET)). Returns 0 on success, -1 on failure.
pub unsafe fn seq_reader_seek(handle: *mut c_void, _offset: i64, _whence: c_int) -> i64 {
    let reader = &mut *(handle as *mut SeqReader);
    if reader.rewind() {
        0
    } else {
        -1
    }
}

/// Close the reader and free memory (like gzclose).
pub unsafe fn seq_reader_close(handle: *mut c_void) -> c_int {
    if !handle.is_null() {
        drop(Box::from_raw(handle as *mut SeqReader));
    }
    0
}
