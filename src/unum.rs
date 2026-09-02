use crate::digest::{sha256_file, sha256_named_files, to_hex};
use crate::model::{AssayProfile, NoCallReason, Target};
use crate::registry::require_observable;
use std::ffi::OsString;
use std::fmt;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnumBamMode {
    Alignment,
    NoAlignment,
}

impl fmt::Display for UnumBamMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Alignment => "alignment",
            Self::NoAlignment => "no-alignment",
        })
    }
}

impl std::str::FromStr for UnumBamMode {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value {
            "alignment" => Ok(Self::Alignment),
            "no-alignment" | "no_alignment" => Ok(Self::NoAlignment),
            _ => Err(format!(
                "unknown Unum BAM mode '{value}'; expected alignment or no-alignment"
            )),
        }
    }
}

#[derive(Debug, Clone)]
pub struct UnumRequest {
    pub target: Target,
    pub assay: AssayProfile,
    pub executable: PathBuf,
    pub bam: PathBuf,
    pub ref_seq_fasta: PathBuf,
    pub ref_coord_fasta: Option<PathBuf>,
    pub reference_genome: Option<PathBuf>,
    pub bam_mode: UnumBamMode,
    pub output_prefix: PathBuf,
    pub threads: u32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnumRun {
    pub backend_version: String,
    pub genotype_tsv: PathBuf,
    pub allele_vcf: Option<PathBuf>,
    pub call_count: u64,
    pub input_sha256: String,
    pub resource_sha256: String,
    pub output_sha256: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnumError(String);

impl UnumError {
    fn new(message: impl Into<String>) -> Self {
        Self(message.into())
    }
}

impl fmt::Display for UnumError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for UnumError {}

pub fn validate_request(request: &UnumRequest) -> Result<(), UnumError> {
    if !matches!(request.target, Target::Hla | Target::Kir) {
        return Err(UnumError::new(format!(
            "Unum backend is restricted to HLA/KIR, not {}",
            request.target
        )));
    }
    require_observable(request.target, request.assay)
        .map_err(|reason| UnumError::new(format!("target is not observable: {reason}")))?;
    if request.threads == 0 {
        return Err(UnumError::new("threads must be >= 1"));
    }
    if request.bam_mode == UnumBamMode::Alignment && request.ref_coord_fasta.is_none() {
        return Err(UnumError::new(
            "Unum alignment mode requires --ref-coord",
        ));
    }
    if request.output_prefix.as_os_str().is_empty() {
        return Err(UnumError::new("output prefix must not be empty"));
    }
    Ok(())
}

/// Run the maintained Unum/T1K HLA/KIR allele-typing backend.
///
/// No shell is involved. The wrapper owns input/resource hashing, normalized
/// call cardinality, and the certificate boundary; Unum owns read extraction,
/// k-mer filtering, banded alignment, abundance estimation, and allele calling.
pub fn run_unum(request: &UnumRequest) -> Result<UnumRun, UnumError> {
    validate_request(request)?;

    if let Some(parent) = request.output_prefix.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent).map_err(|error| {
                UnumError::new(format!(
                    "cannot create output directory '{}': {error}",
                    parent.display()
                ))
            })?;
        }
    }

    let backend_version = unum_version(&request.executable);

    let mut command = Command::new(&request.executable);
    command
        .arg("run")
        .arg("-b")
        .arg(&request.bam)
        .arg("-f")
        .arg(&request.ref_seq_fasta)
        .arg("--bam-mode")
        .arg(request.bam_mode.to_string())
        .arg("-o")
        .arg(&request.output_prefix)
        .arg("-t")
        .arg(request.threads.to_string());

    if let Some(path) = &request.ref_coord_fasta {
        command.arg("-c").arg(path);
    }
    if let Some(path) = &request.reference_genome {
        command.arg("-r").arg(path);
    }

    let output = command
        .output()
        .map_err(|error| UnumError::new(format!("failed to execute Unum: {error}")))?;
    if !output.status.success() {
        return Err(UnumError::new(format!(
            "Unum exited with {}:\n{}",
            output.status,
            String::from_utf8_lossy(&output.stderr).trim()
        )));
    }

    let genotype_tsv = suffixed_path(&request.output_prefix, "_genotype.tsv");
    if !genotype_tsv.is_file() {
        return Err(UnumError::new(format!(
            "Unum succeeded but did not create '{}'",
            genotype_tsv.display()
        )));
    }
    let allele_vcf_path = suffixed_path(&request.output_prefix, "_allele.vcf");
    let allele_vcf = allele_vcf_path.is_file().then_some(allele_vcf_path);

    let call_count = count_genotype_calls(&genotype_tsv)?;
    let input_sha256 = to_hex(
        &sha256_file(&request.bam)
            .map_err(|error| UnumError::new(format!("cannot hash BAM/CRAM: {error}")))?,
    );

    let mut resources = vec![("ref-seq", request.ref_seq_fasta.as_path())];
    if let Some(path) = &request.ref_coord_fasta {
        resources.push(("ref-coord", path.as_path()));
    }
    if let Some(path) = &request.reference_genome {
        resources.push(("reference-genome", path.as_path()));
    }
    let resource_sha256 = to_hex(
        &sha256_named_files(&resources)
            .map_err(|error| UnumError::new(format!("cannot hash Unum resources: {error}")))?,
    );

    let mut outputs = vec![("genotype-tsv", genotype_tsv.as_path())];
    if let Some(path) = &allele_vcf {
        outputs.push(("allele-vcf", path.as_path()));
    }
    let output_sha256 = to_hex(
        &sha256_named_files(&outputs)
            .map_err(|error| UnumError::new(format!("cannot hash Unum outputs: {error}")))?,
    );

    Ok(UnumRun {
        backend_version,
        genotype_tsv,
        allele_vcf,
        call_count,
        input_sha256,
        resource_sha256,
        output_sha256,
    })
}

fn unum_version(executable: &Path) -> String {
    Command::new(executable)
        .arg("--version")
        .output()
        .ok()
        .filter(|output| output.status.success())
        .and_then(|output| String::from_utf8(output.stdout).ok())
        .map(|value| sanitize_version(value.trim()))
        .filter(|value| !value.is_empty())
        .unwrap_or_else(|| "unknown".to_string())
}

fn sanitize_version(value: &str) -> String {
    value
        .chars()
        .map(|character| {
            if matches!(character, '\n' | '\r' | '\t' | '=') {
                ' '
            } else {
                character
            }
        })
        .collect::<String>()
        .trim()
        .to_string()
}

fn count_genotype_calls(path: &Path) -> Result<u64, UnumError> {
    let text = fs::read_to_string(path).map_err(|error| {
        UnumError::new(format!(
            "cannot read Unum genotype output '{}': {error}",
            path.display()
        ))
    })?;

    Ok(text
        .lines()
        .filter(|raw_line| {
            let line = raw_line.trim();
            if line.is_empty() || line.starts_with('#') {
                return false;
            }
            let first = line.split_whitespace().next().unwrap_or_default();
            !first.eq_ignore_ascii_case("gene")
        })
        .count() as u64)
}

fn suffixed_path(prefix: &Path, suffix: &str) -> PathBuf {
    let mut value = OsString::from(prefix.as_os_str());
    value.push(suffix);
    PathBuf::from(value)
}

pub fn no_call_reason_for_zero_calls(run: &UnumRun) -> Option<NoCallReason> {
    (run.call_count == 0).then_some(NoCallReason::InsufficientEvidence)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn request(target: Target) -> UnumRequest {
        UnumRequest {
            target,
            assay: AssayProfile::wgs(),
            executable: "unum".into(),
            bam: "sample.bam".into(),
            ref_seq_fasta: "alleles.fa".into(),
            ref_coord_fasta: Some("coords.fa".into()),
            reference_genome: None,
            bam_mode: UnumBamMode::Alignment,
            output_prefix: "out/sample".into(),
            threads: 2,
        }
    }

    #[test]
    fn only_hla_and_kir_are_accepted() {
        assert_eq!(validate_request(&request(Target::Hla)), Ok(()));
        assert_eq!(validate_request(&request(Target::Kir)), Ok(()));
        assert!(validate_request(&request(Target::Hba)).is_err());
    }

    #[test]
    fn alignment_mode_requires_coordinates() {
        let mut value = request(Target::Hla);
        value.ref_coord_fasta = None;
        assert!(validate_request(&value).is_err());
    }

    #[test]
    fn output_suffix_does_not_replace_extension() {
        assert_eq!(
            suffixed_path(Path::new("sample.name"), "_genotype.tsv"),
            PathBuf::from("sample.name_genotype.tsv")
        );
    }
}
