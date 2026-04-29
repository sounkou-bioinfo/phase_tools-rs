#[cfg(target_arch = "x86_64")]
use libc::{c_char, c_int, c_void};
#[cfg(target_arch = "x86_64")]
use std::ffi::{CStr, CString};
#[cfg(target_arch = "x86_64")]
use std::ptr;

#[cfg(target_arch = "x86_64")]
#[allow(non_camel_case_types, non_snake_case, non_upper_case_globals)]
mod ffi {
    include!(concat!(env!("OUT_DIR"), "/fermi_lite_bindings.rs"));
}

#[derive(Debug, Clone)]
pub struct Unitig {
    pub seq: String,
    pub len: i32,
    pub supporting_reads: i32,
}

#[derive(Debug, Clone)]
pub struct AssemblyRead {
    pub seq: String,
    pub qual: Option<String>,
}

impl AssemblyRead {
    pub fn sequence<S: Into<String>>(seq: S) -> Self {
        Self {
            seq: seq.into(),
            qual: None,
        }
    }

    pub fn fastq<S: Into<String>, Q: Into<String>>(seq: S, qual: Q) -> Self {
        Self {
            seq: seq.into(),
            qual: Some(qual.into()),
        }
    }
}

#[derive(Debug, Clone)]
pub struct AssembleOptions {
    pub threads: i32,
    pub min_asm_overlap: i32,
    pub min_count: i32,
    pub max_count: i32,
    pub error_correction_k: i32,
}

impl Default for AssembleOptions {
    fn default() -> Self {
        Self {
            threads: 1,
            min_asm_overlap: 21,
            min_count: 1,
            max_count: 1000,
            // fermi-lite uses ec_k < 0 to skip error correction. This is more
            // predictable for tiny local adjudication windows and smoke tests.
            // Use --ec-k 0 in the CLI to enable fermi-lite's auto EC path,
            // which can consume FASTQ/BAM base qualities when supplied.
            error_correction_k: -1,
        }
    }
}

#[cfg(target_arch = "x86_64")]
fn malloc_c_string(value: &str, label: &str) -> Result<*mut c_char, String> {
    let c = CString::new(value).map_err(|_| format!("{label} contains NUL byte"))?;
    let bytes = c.as_bytes_with_nul();
    let ptr = unsafe { libc::malloc(bytes.len()) as *mut u8 };
    if ptr.is_null() {
        return Err(format!("malloc failed for fermi-lite {label}"));
    }
    unsafe {
        ptr.copy_from_nonoverlapping(bytes.as_ptr(), bytes.len());
    }
    Ok(ptr as *mut c_char)
}

#[cfg(target_arch = "x86_64")]
fn cleanup_partial(seqs_ptr: *mut ffi::bseq1_t, n_init: usize) {
    unsafe {
        for j in 0..n_init {
            let rec = seqs_ptr.add(j);
            if !(*rec).seq.is_null() {
                libc::free((*rec).seq as *mut c_void);
            }
            if !(*rec).qual.is_null() {
                libc::free((*rec).qual as *mut c_void);
            }
        }
        libc::free(seqs_ptr as *mut c_void);
    }
}

#[allow(dead_code)]
pub fn assemble_sequences<S: AsRef<str>>(
    sequences: &[S],
    options: &AssembleOptions,
) -> Result<Vec<Unitig>, String> {
    let reads = sequences
        .iter()
        .map(|seq| AssemblyRead::sequence(seq.as_ref()))
        .collect::<Vec<_>>();
    assemble_reads(&reads, options)
}

#[cfg(target_arch = "x86_64")]
pub fn assemble_reads(
    reads: &[AssemblyRead],
    options: &AssembleOptions,
) -> Result<Vec<Unitig>, String> {
    let clean = reads
        .iter()
        .map(|read| AssemblyRead {
            seq: read.seq.trim().to_ascii_uppercase(),
            qual: read.qual.as_ref().map(|q| q.trim().to_string()),
        })
        .filter(|read| !read.seq.is_empty())
        .collect::<Vec<_>>();
    if clean.is_empty() {
        return Ok(Vec::new());
    }
    if clean.len() > c_int::MAX as usize {
        return Err("too many sequences for fermi-lite".to_string());
    }

    let n = clean.len();
    let seqs_ptr =
        unsafe { libc::calloc(n, std::mem::size_of::<ffi::bseq1_t>()) as *mut ffi::bseq1_t };
    if seqs_ptr.is_null() {
        return Err("calloc failed for fermi-lite reads".to_string());
    }

    for (i, read) in clean.iter().enumerate() {
        let seq_ptr = match malloc_c_string(&read.seq, "sequence") {
            Ok(ptr) => ptr,
            Err(e) => {
                cleanup_partial(seqs_ptr, i);
                return Err(e);
            }
        };
        let qual_ptr = match &read.qual {
            Some(qual) => {
                if qual.len() != read.seq.len() {
                    cleanup_partial(seqs_ptr, i);
                    unsafe { libc::free(seq_ptr as *mut c_void) };
                    return Err(format!(
                        "quality length ({}) does not match sequence length ({})",
                        qual.len(),
                        read.seq.len()
                    ));
                }
                match malloc_c_string(qual, "quality string") {
                    Ok(ptr) => ptr,
                    Err(e) => {
                        cleanup_partial(seqs_ptr, i);
                        unsafe { libc::free(seq_ptr as *mut c_void) };
                        return Err(e);
                    }
                }
            }
            None => ptr::null_mut(),
        };
        unsafe {
            let rec = seqs_ptr.add(i);
            (*rec).l_seq = read.seq.len() as i32;
            (*rec).seq = seq_ptr;
            (*rec).qual = qual_ptr;
        }
    }

    let mut opt = unsafe {
        let mut opt: ffi::fml_opt_t = std::mem::zeroed();
        ffi::fml_opt_init(&mut opt);
        opt
    };
    opt.n_threads = options.threads.max(1);
    opt.min_asm_ovlp = options.min_asm_overlap.max(1);
    opt.min_cnt = options.min_count.max(1);
    opt.max_cnt = options.max_count.max(opt.min_cnt);
    opt.ec_k = options.error_correction_k;

    let mut n_utg = 0i32;
    let utg_ptr = unsafe {
        ffi::fm_verbose = 0;
        ffi::fml_assemble(&opt, n as c_int, seqs_ptr, &mut n_utg)
    };
    if utg_ptr.is_null() {
        return Ok(Vec::new());
    }
    if n_utg <= 0 {
        unsafe { ffi::fml_utg_destroy(n_utg, utg_ptr) };
        return Ok(Vec::new());
    }

    let mut out = Vec::new();
    for i in 0..n_utg as usize {
        let utg = unsafe { &*utg_ptr.add(i) };
        if utg.seq.is_null() {
            continue;
        }
        let seq = unsafe { CStr::from_ptr(utg.seq) }
            .to_string_lossy()
            .into_owned();
        out.push(Unitig {
            seq,
            len: utg.len,
            supporting_reads: utg.nsr,
        });
    }
    unsafe { ffi::fml_utg_destroy(n_utg, utg_ptr) };
    Ok(out)
}

#[cfg(not(target_arch = "x86_64"))]
pub fn assemble_reads(
    _reads: &[AssemblyRead],
    _options: &AssembleOptions,
) -> Result<Vec<Unitig>, String> {
    Err("fermi-lite assembly is currently enabled only on x86_64 targets because upstream ksw.c requires SSE2".to_string())
}
