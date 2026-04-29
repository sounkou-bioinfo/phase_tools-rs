// Fallback fermi-lite bindings used when build-time libclang is unavailable.
// build.rs prefers bindgen-generated bindings from vendor/fermi-lite/fml.h and
// copies this narrow ABI-compatible subset only as a portability fallback.

use std::os::raw::{c_char, c_float, c_int, c_void};

#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct bseq1_t {
    pub l_seq: i32,
    pub seq: *mut c_char,
    pub qual: *mut c_char,
}

#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct magopt_t {
    pub flag: c_int,
    pub min_ovlp: c_int,
    pub min_elen: c_int,
    pub min_ensr: c_int,
    pub min_insr: c_int,
    pub max_bdist: c_int,
    pub max_bdiff: c_int,
    pub max_bvtx: c_int,
    pub min_merge_len: c_int,
    pub trim_len: c_int,
    pub trim_depth: c_int,
    pub min_dratio1: c_float,
    pub max_bcov: c_float,
    pub max_bfrac: c_float,
}

#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct fml_opt_t {
    pub n_threads: c_int,
    pub ec_k: c_int,
    pub min_cnt: c_int,
    pub max_cnt: c_int,
    pub min_asm_ovlp: c_int,
    pub min_merge_len: c_int,
    pub mag_opt: magopt_t,
}

pub type fml_ovlp_t = c_void;

#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct fml_utg_t {
    pub len: i32,
    pub nsr: i32,
    pub seq: *mut c_char,
    pub cov: *mut c_char,
    pub n_ovlp: [i32; 2],
    pub ovlp: *mut fml_ovlp_t,
}

extern "C" {
    pub fn fml_opt_init(opt: *mut fml_opt_t);
    pub fn fml_assemble(
        opt: *const fml_opt_t,
        n_seqs: c_int,
        seqs: *mut bseq1_t,
        n_utg: *mut c_int,
    ) -> *mut fml_utg_t;
    pub fn fml_utg_destroy(n_utg: c_int, utg: *mut fml_utg_t);
    pub static mut fm_verbose: c_int;
}
