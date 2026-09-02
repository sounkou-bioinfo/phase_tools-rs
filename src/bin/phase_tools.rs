use phase_tools::certificate::{
    DecisionCertificate, SelectionWitness, PROOF_CONTRACT, REGISTRY_VERSION,
};
use phase_tools::digest::{sha256_bytes, sha256_file, to_hex};
use phase_tools::hba::{
    read_evidence, read_hypotheses, select_hba, HbaDecision, HbaHypothesis, HbaOutcome,
};
use phase_tools::model::{AssayKind, AssayProfile, CallStatus, NoCallReason, Target};
use phase_tools::registry::{require_observable, spec, TARGET_SPECS};
use phase_tools::unum::{run_unum, UnumBamMode, UnumRequest};
use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::path::{Path, PathBuf};
use std::str::FromStr;

fn main() {
    if let Err(error) = run() {
        eprintln!("phase-tools: {error}");
        std::process::exit(2);
    }
}

fn run() -> Result<(), String> {
    let arguments = env::args().skip(1).collect::<Vec<_>>();
    let Some((command, rest)) = arguments.split_first() else {
        print_help();
        return Ok(());
    };

    match command.as_str() {
        "targets" => command_targets(rest),
        "hba" => command_hba(rest),
        "unum" => command_unum(rest),
        "verify" => command_verify(rest),
        "proof-contract" => {
            command_proof_contract();
            Ok(())
        }
        "help" | "--help" | "-h" => {
            print_help();
            Ok(())
        }
        "version" | "--version" | "-V" => {
            println!("phase-tools {}", env!("CARGO_PKG_VERSION"));
            Ok(())
        }
        _ => Err(format!("unknown command '{command}'; run phase-tools help")),
    }
}

fn command_targets(arguments: &[String]) -> Result<(), String> {
    let options = Options::parse(arguments, &["--validated-enrichment"])?;
    options.ensure_allowed(&["--assay", "--validated-enrichment"])?;

    let profile = options
        .value("--assay")
        .map(|value| parse_assay_profile(value, options.flag("--validated-enrichment")))
        .transpose()?;

    println!("target\tfamily\tbackend\timplementation\twgs\twes\tevidence\tselected_profile");
    for entry in TARGET_SPECS {
        let evidence = entry
            .evidence
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(",");
        let selected = profile.map_or_else(
            || ".".to_string(),
            |assay| match require_observable(entry.target, assay) {
                Ok(()) => "observable".to_string(),
                Err(reason) => reason.to_string(),
            },
        );
        println!(
            "{}\t{}\t{}\t{}\tobservable\t{}\t{}\t{}",
            entry.target,
            entry.family,
            entry.backend,
            entry.implementation,
            entry.wes,
            evidence,
            selected
        );
    }

    Ok(())
}

fn command_hba(arguments: &[String]) -> Result<(), String> {
    let options = Options::parse(arguments, &["--validated-enrichment"])?;
    options.ensure_allowed(&[
        "--assay",
        "--validated-enrichment",
        "--evidence",
        "--hypotheses",
        "--min-margin",
        "--certificate",
    ])?;

    let assay = parse_assay_profile(
        options.required("--assay")?,
        options.flag("--validated-enrichment"),
    )?;
    let evidence_path = PathBuf::from(options.required("--evidence")?);
    let hypotheses_path = PathBuf::from(options.required("--hypotheses")?);
    let certificate_path = PathBuf::from(options.required("--certificate")?);
    let required_margin = options
        .value("--min-margin")
        .map(|value| parse_u64(value, "--min-margin"))
        .transpose()?
        .unwrap_or(10);

    let input_sha256 = hash_file(&evidence_path, "HBA evidence")?;
    let resource_sha256 = hash_file(&hypotheses_path, "HBA hypothesis resource")?;

    let (output, status, call_count, selection) = match require_observable(Target::Hba, assay) {
        Err(reason) => (
            format_hba_no_call(reason, None, None, &[]),
            CallStatus::NoCall(reason),
            0,
            None,
        ),
        Ok(()) => {
            let evidence = read_evidence(&evidence_path).map_err(|error| error.to_string())?;
            let hypotheses =
                read_hypotheses(&hypotheses_path).map_err(|error| error.to_string())?;
            let decision = select_hba(&evidence, &hypotheses, required_margin)
                .map_err(|error| error.to_string())?;
            certificate_parts_for_hba(&decision, &hypotheses, required_margin)?
        }
    };

    print!("{output}");
    let output_sha256 = to_hex(&sha256_bytes(output.as_bytes()));
    let certificate = DecisionCertificate {
        target: Target::Hba,
        assay,
        backend: spec(Target::Hba).backend.to_string(),
        backend_version: env!("CARGO_PKG_VERSION").to_string(),
        status,
        call_count,
        input_sha256,
        resource_sha256,
        output_sha256,
        selection,
    };
    validate_and_write_certificate(&certificate, &certificate_path)
}

fn certificate_parts_for_hba(
    decision: &HbaDecision,
    hypotheses: &[HbaHypothesis],
    required_margin: u64,
) -> Result<(String, CallStatus, u64, Option<SelectionWitness>), String> {
    match &decision.outcome {
        HbaOutcome::Called {
            winner_index,
            winner_penalty,
            runner_up_penalty,
            margin,
        } => {
            let hypothesis = hypotheses
                .get(*winner_index)
                .ok_or_else(|| "HBA winner index is outside the hypothesis resource".to_string())?;
            let output = format!(
                concat!(
                    "status\ttarget\thypothesis\tcall\thba1_cn\thba2_cn\t",
                    "penalty\trunner_up_penalty\tmargin\treason\n",
                    "called\tHBA\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t.\n"
                ),
                hypothesis.id,
                hypothesis.call,
                hypothesis.hba1_copy_number,
                hypothesis.hba2_copy_number,
                winner_penalty,
                runner_up_penalty,
                margin,
            );
            let selection = SelectionWitness {
                candidate_count: hypotheses.len() as u64,
                winner_index: *winner_index as u64,
                winner_penalty: *winner_penalty,
                runner_up_penalty: *runner_up_penalty,
                required_margin,
            };
            Ok((output, CallStatus::Called, 1, Some(selection)))
        }
        HbaOutcome::NoCall {
            reason,
            top_penalty,
            runner_up_penalty,
            missing_features,
        } => Ok((
            format_hba_no_call(*reason, *top_penalty, *runner_up_penalty, missing_features),
            CallStatus::NoCall(*reason),
            0,
            None,
        )),
    }
}

fn format_hba_no_call(
    reason: NoCallReason,
    top_penalty: Option<u64>,
    runner_up_penalty: Option<u64>,
    missing_features: &[String],
) -> String {
    let top = top_penalty.map_or_else(|| ".".to_string(), |value| value.to_string());
    let runner = runner_up_penalty.map_or_else(|| ".".to_string(), |value| value.to_string());
    let detail = if missing_features.is_empty() {
        reason.to_string()
    } else {
        format!("{reason}:{}", missing_features.join(","))
    };
    format!(
        concat!(
            "status\ttarget\thypothesis\tcall\thba1_cn\thba2_cn\t",
            "penalty\trunner_up_penalty\tmargin\treason\n",
            "no-call\tHBA\t.\t.\t.\t.\t{}\t{}\t.\t{}\n"
        ),
        top, runner, detail
    )
}

fn command_unum(arguments: &[String]) -> Result<(), String> {
    let options = Options::parse(arguments, &["--validated-enrichment"])?;
    options.ensure_allowed(&[
        "--target",
        "--assay",
        "--validated-enrichment",
        "--unum",
        "--bam",
        "--ref-seq",
        "--ref-coord",
        "--reference",
        "--bam-mode",
        "--output-prefix",
        "--threads",
        "--certificate",
    ])?;

    let target = options.required("--target")?.parse::<Target>()?;
    if !matches!(target, Target::Hla | Target::Kir) {
        return Err("--target for the unum command must be HLA or KIR".to_string());
    }
    let assay = parse_assay_profile(
        options.required("--assay")?,
        options.flag("--validated-enrichment"),
    )?;
    let request = UnumRequest {
        target,
        assay,
        executable: PathBuf::from(options.required("--unum")?),
        bam: PathBuf::from(options.required("--bam")?),
        ref_seq_fasta: PathBuf::from(options.required("--ref-seq")?),
        ref_coord_fasta: options.value("--ref-coord").map(PathBuf::from),
        reference_genome: options.value("--reference").map(PathBuf::from),
        bam_mode: options.required("--bam-mode")?.parse::<UnumBamMode>()?,
        output_prefix: PathBuf::from(options.required("--output-prefix")?),
        threads: options
            .value("--threads")
            .map(|value| parse_u32(value, "--threads"))
            .transpose()?
            .unwrap_or(1),
    };
    let certificate_path = PathBuf::from(options.required("--certificate")?);

    let backend_run = run_unum(&request).map_err(|error| error.to_string())?;
    let status = if backend_run.call_count == 0 {
        CallStatus::NoCall(NoCallReason::InsufficientEvidence)
    } else {
        CallStatus::Called
    };
    let certificate = DecisionCertificate {
        target,
        assay,
        backend: spec(target).backend.to_string(),
        backend_version: backend_run.backend_version.clone(),
        status,
        call_count: backend_run.call_count,
        input_sha256: backend_run.input_sha256,
        resource_sha256: backend_run.resource_sha256,
        output_sha256: backend_run.output_sha256,
        selection: None,
    };
    validate_and_write_certificate(&certificate, &certificate_path)?;

    println!("status\ttarget\tcalls\tgenotype_tsv\tallele_vcf\tcertificate");
    println!(
        "{}\t{}\t{}\t{}\t{}\t{}",
        match status {
            CallStatus::Called => "called",
            CallStatus::NoCall(_) => "no-call",
        },
        target,
        backend_run.call_count,
        backend_run.genotype_tsv.display(),
        backend_run
            .allele_vcf
            .as_deref()
            .map_or_else(|| ".".to_string(), |path| path.display().to_string()),
        certificate_path.display()
    );

    Ok(())
}

fn command_verify(arguments: &[String]) -> Result<(), String> {
    let options = Options::parse(arguments, &[])?;
    options.ensure_allowed(&["--certificate"])?;
    let path = PathBuf::from(options.required("--certificate")?);
    let certificate = DecisionCertificate::read(&path)?;
    match certificate.verify() {
        Ok(()) => {
            println!(
                "verified\ttarget={}\tstatus={}\tproof_contract={}",
                certificate.target,
                match certificate.status {
                    CallStatus::Called => "called",
                    CallStatus::NoCall(_) => "no-call",
                },
                PROOF_CONTRACT
            );
            Ok(())
        }
        Err(errors) => Err(format!(
            "certificate '{}' failed verification:\n- {}",
            path.display(),
            errors.join("\n- ")
        )),
    }
}

fn command_proof_contract() {
    println!("registry={REGISTRY_VERSION}");
    println!("proof_contract={PROOF_CONTRACT}");
    println!("lean_registry_theorem=PhaseTools.mem_allTargets");
    println!("lean_dragen_closure_theorem=PhaseTools.dragen45_closed");
    println!("lean_wes_theorem=PhaseTools.dragen45_wes_iff");
    println!("lean_selection_theorem=PhaseTools.verifySelection_sound");
}

fn validate_and_write_certificate(
    certificate: &DecisionCertificate,
    path: &Path,
) -> Result<(), String> {
    if let Err(errors) = certificate.verify() {
        return Err(format!(
            "internal certificate construction failed:\n- {}",
            errors.join("\n- ")
        ));
    }
    certificate.write(path)
}

fn parse_assay_profile(value: &str, validated_enrichment: bool) -> Result<AssayProfile, String> {
    let kind = AssayKind::from_str(value)?;
    if kind == AssayKind::Wgs && validated_enrichment {
        return Err("--validated-enrichment is meaningful only with --assay wes".to_string());
    }
    Ok(AssayProfile {
        kind,
        validated_target_enrichment: validated_enrichment,
    })
}

fn hash_file(path: &Path, description: &str) -> Result<String, String> {
    sha256_file(path)
        .map(|digest| to_hex(&digest))
        .map_err(|error| format!("cannot hash {description} '{}': {error}", path.display()))
}

fn parse_u64(value: &str, option: &str) -> Result<u64, String> {
    value
        .parse::<u64>()
        .map_err(|_| format!("{option} requires an unsigned integer, got '{value}'"))
}

fn parse_u32(value: &str, option: &str) -> Result<u32, String> {
    let parsed = value
        .parse::<u32>()
        .map_err(|_| format!("{option} requires a positive integer, got '{value}'"))?;
    if parsed == 0 {
        Err(format!("{option} must be >= 1"))
    } else {
        Ok(parsed)
    }
}

#[derive(Debug, Default)]
struct Options {
    values: BTreeMap<String, String>,
    flags: BTreeSet<String>,
}

impl Options {
    fn parse(arguments: &[String], flag_names: &[&str]) -> Result<Self, String> {
        let mut output = Self::default();
        let mut index = 0;

        while index < arguments.len() {
            let option = &arguments[index];
            if !option.starts_with("--") {
                return Err(format!("unexpected positional argument '{option}'"));
            }
            if flag_names.contains(&option.as_str()) {
                if !output.flags.insert(option.clone()) {
                    return Err(format!("duplicate option '{option}'"));
                }
                index += 1;
                continue;
            }

            let value = arguments
                .get(index + 1)
                .ok_or_else(|| format!("option '{option}' requires a value"))?;
            if value.starts_with("--") {
                return Err(format!("option '{option}' requires a value"));
            }
            if output
                .values
                .insert(option.clone(), value.clone())
                .is_some()
            {
                return Err(format!("duplicate option '{option}'"));
            }
            index += 2;
        }

        Ok(output)
    }

    fn ensure_allowed(&self, allowed: &[&str]) -> Result<(), String> {
        for key in self.values.keys().chain(self.flags.iter()) {
            if !allowed.contains(&key.as_str()) {
                return Err(format!("unknown option '{key}'"));
            }
        }
        Ok(())
    }

    fn required(&self, key: &str) -> Result<&str, String> {
        self.value(key)
            .ok_or_else(|| format!("missing required option '{key}'"))
    }

    fn value(&self, key: &str) -> Option<&str> {
        self.values.get(key).map(String::as_str)
    }

    fn flag(&self, key: &str) -> bool {
        self.flags.contains(key)
    }
}

fn print_help() {
    println!(
        r#"phase-tools: proof-carrying short-read callers for difficult loci

USAGE
  phase-tools targets [--assay wgs|wes] [--validated-enrichment]

  phase-tools hba
    --assay wgs|wes [--validated-enrichment]
    --evidence FILE
    --hypotheses FILE
    [--min-margin INTEGER]
    --certificate FILE

  phase-tools unum
    --target HLA|KIR
    --assay wgs|wes [--validated-enrichment]
    --unum PATH
    --bam BAM_OR_CRAM
    --ref-seq ALLELE_FASTA
    [--ref-coord COORD_FASTA]
    [--reference GENOME_FASTA]
    --bam-mode alignment|no-alignment
    --output-prefix PATH
    [--threads INTEGER]
    --certificate FILE

  phase-tools verify --certificate FILE
  phase-tools proof-contract

The HBA command consumes prepared integer evidence and a versioned hypothesis
catalogue. The Unum command delegates HLA/KIR allele inference to the maintained
Rust T1K port while this tool owns observability, hashing, normalized call
cardinality, and the decision certificate.
"#
    );
}
