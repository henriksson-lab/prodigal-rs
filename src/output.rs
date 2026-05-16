use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Write};
use std::os::raw::c_int;
use std::sync::{Mutex, OnceLock};

const STDOUT_HANDLE: c_int = 1;

struct OutputRegistry {
    next: c_int,
    writers: HashMap<c_int, Box<dyn Write + Send>>,
}

impl OutputRegistry {
    /// Create an empty registry whose first issued handle will be 3
    /// (reserving 0/1/2 for the standard streams).
    fn new() -> Self {
        Self {
            next: 3,
            writers: HashMap::new(),
        }
    }

    /// Store `writer` in the registry and return a freshly assigned handle.
    fn insert(&mut self, writer: Box<dyn Write + Send>) -> c_int {
        let handle = self.next;
        self.next += 1;
        self.writers.insert(handle, writer);
        handle
    }
}

/// Lazily initialize and return the process-wide output registry.
fn registry() -> &'static Mutex<OutputRegistry> {
    static REGISTRY: OnceLock<Mutex<OutputRegistry>> = OnceLock::new();
    REGISTRY.get_or_init(|| Mutex::new(OutputRegistry::new()))
}

/// Return the sentinel handle that routes writes to standard output.
pub fn stdout_handle() -> c_int {
    STDOUT_HANDLE
}

/// Create `path`, register a buffered writer for it, and return its handle.
pub fn create_file(path: &str) -> io::Result<c_int> {
    let file = File::create(path)?;
    let mut registry = registry().lock().expect("output registry poisoned");
    Ok(registry.insert(Box::new(io::BufWriter::new(file))))
}

/// Write `text` to the writer identified by `handle`, routing to stdout when
/// `handle` is the stdout sentinel and erroring on unknown handles.
pub fn write_to_handle(handle: c_int, text: &str) -> io::Result<()> {
    if handle == STDOUT_HANDLE {
        let mut stdout = io::stdout().lock();
        stdout.write_all(text.as_bytes())?;
        return Ok(());
    }

    let mut registry = registry().lock().expect("output registry poisoned");
    match registry.writers.get_mut(&handle) {
        Some(writer) => writer.write_all(text.as_bytes()),
        None => Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!("unknown output handle {handle}"),
        )),
    }
}

/// Flush and (for non-stdout handles) drop the writer associated with
/// `handle`; a no-op if the handle is unknown.
pub fn close_handle(handle: c_int) {
    if handle == STDOUT_HANDLE {
        let _ = io::stdout().lock().flush();
        return;
    }

    let mut registry = registry().lock().expect("output registry poisoned");
    if let Some(mut writer) = registry.writers.remove(&handle) {
        let _ = writer.flush();
    }
}
