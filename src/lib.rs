#![forbid(unsafe_code)]
//! Difficult-locus target callers with explicit assay observability and
//! machine-checkable decision certificates.
//!
//! The public surface is intentionally finite. Phasing, local assembly, and
//! read-evidence algorithms may exist as internal kernels, but this crate does
//! not expose a grab bag of unrelated BAM/VCF utilities.

pub mod certificate;
pub mod digest;
pub mod hba;
pub mod model;
pub mod registry;
pub mod unum;
