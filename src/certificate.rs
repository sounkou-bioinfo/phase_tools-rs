use crate::digest::is_sha256_hex;
use crate::model::{
    AssayKind, AssayProfile, CallStatus, ImplementationStatus, NoCallReason, Target,
};
use crate::registry::{require_observable, spec};
use std::collections::BTreeMap;
use std::fmt;
use std::fs;
use std::path::Path;
use std::str::FromStr;

pub const CERTIFICATE_SCHEMA: &str = "phase-tools-certificate-v1";
pub const REGISTRY_VERSION: &str = "dragen-4.5-plus-hla-kir-v1";
pub const PROOF_CONTRACT: &str = "PhaseTools.Certificate.V1";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SelectionWitness {
    pub candidate_count: u64,
    pub winner_index: u64,
    pub winner_penalty: u64,
    pub runner_up_penalty: u64,
    pub required_margin: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DecisionCertificate {
    pub target: Target,
    pub assay: AssayProfile,
    pub backend: String,
    pub backend_version: String,
    pub status: CallStatus,
    pub call_count: u64,
    pub input_sha256: String,
    pub resource_sha256: String,
    pub output_sha256: String,
    pub selection: Option<SelectionWitness>,
}

impl DecisionCertificate {
    #[must_use]
    pub fn to_text(&self) -> String {
        let (status, reason) = match self.status {
            CallStatus::Called => ("called", ".".to_string()),
            CallStatus::NoCall(reason) => ("no-call", reason.to_string()),
        };
        let (
            candidate_count,
            winner_index,
            winner_penalty,
            runner_up_penalty,
            required_margin,
        ) = if let Some(selection) = self.selection {
            (
                selection.candidate_count.to_string(),
                selection.winner_index.to_string(),
                selection.winner_penalty.to_string(),
                selection.runner_up_penalty.to_string(),
                selection.required_margin.to_string(),
            )
        } else {
            (
                ".".to_string(),
                ".".to_string(),
                ".".to_string(),
                ".".to_string(),
                ".".to_string(),
            )
        };

        format!(
            concat!(
                "schema={}\n",
                "registry={}\n",
                "proof_contract={}\n",
                "target={}\n",
                "assay={}\n",
                "validated_target_enrichment={}\n",
                "backend={}\n",
                "backend_version={}\n",
                "status={}\n",
                "reason={}\n",
                "call_count={}\n",
                "input_sha256={}\n",
                "resource_sha256={}\n",
                "output_sha256={}\n",
                "candidate_count={}\n",
                "winner_index={}\n",
                "winner_penalty={}\n",
                "runner_up_penalty={}\n",
                "required_margin={}\n"
            ),
            CERTIFICATE_SCHEMA,
            REGISTRY_VERSION,
            PROOF_CONTRACT,
            self.target,
            self.assay.kind,
            self.assay.validated_target_enrichment,
            self.backend,
            self.backend_version,
            status,
            reason,
            self.call_count,
            self.input_sha256,
            self.resource_sha256,
            self.output_sha256,
            candidate_count,
            winner_index,
            winner_penalty,
            runner_up_penalty,
            required_margin,
        )
    }

    pub fn write(&self, path: &Path) -> Result<(), String> {
        fs::write(path, self.to_text())
            .map_err(|error| format!("cannot write certificate '{}': {error}", path.display()))
    }

    pub fn read(path: &Path) -> Result<Self, String> {
        let text = fs::read_to_string(path)
            .map_err(|error| format!("cannot read certificate '{}': {error}", path.display()))?;
        Self::parse(&text)
    }

    pub fn parse(text: &str) -> Result<Self, String> {
        let mut fields = BTreeMap::<String, String>::new();
        for (line_index, raw_line) in text.lines().enumerate() {
            let line = raw_line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let Some((key, value)) = line.split_once('=') else {
                return Err(format!(
                    "certificate line {} must be key=value",
                    line_index + 1
                ));
            };
            if key.is_empty()
                || value
                    .chars()
                    .any(|character| matches!(character, '\n' | '\r'))
            {
                return Err(format!("invalid certificate line {}", line_index + 1));
            }
            if fields.insert(key.to_string(), value.to_string()).is_some() {
                return Err(format!("duplicate certificate key '{key}'"));
            }
        }

        const ALLOWED: &[&str] = &[
            "schema",
            "registry",
            "proof_contract",
            "target",
            "assay",
            "validated_target_enrichment",
            "backend",
            "backend_version",
            "status",
            "reason",
            "call_count",
            "input_sha256",
            "resource_sha256",
            "output_sha256",
            "candidate_count",
            "winner_index",
            "winner_penalty",
            "runner_up_penalty",
            "required_margin",
        ];
        for key in fields.keys() {
            if !ALLOWED.contains(&key.as_str()) {
                return Err(format!("unknown certificate key '{key}'"));
            }
        }

        require_exact(&fields, "schema", CERTIFICATE_SCHEMA)?;
        require_exact(&fields, "registry", REGISTRY_VERSION)?;
        require_exact(&fields, "proof_contract", PROOF_CONTRACT)?;

        let target = required(&fields, "target")?.parse::<Target>()?;
        let kind = required(&fields, "assay")?.parse::<AssayKind>()?;
        let validated_target_enrichment =
            parse_bool(required(&fields, "validated_target_enrichment")?)?;
        let backend = required(&fields, "backend")?.to_string();
        let backend_version = required(&fields, "backend_version")?.to_string();
        let call_count = parse_u64(required(&fields, "call_count")?, "call_count")?;
        let input_sha256 = required(&fields, "input_sha256")?.to_string();
        let resource_sha256 = required(&fields, "resource_sha256")?.to_string();
        let output_sha256 = required(&fields, "output_sha256")?.to_string();

        let status = match required(&fields, "status")? {
            "called" => {
                require_exact(&fields, "reason", ".")?;
                CallStatus::Called
            }
            "no-call" => {
                let reason = required(&fields, "reason")?.parse::<NoCallReason>()?;
                CallStatus::NoCall(reason)
            }
            value => return Err(format!("invalid certificate status '{value}'")),
        };

        let selection_values = [
            required(&fields, "candidate_count")?,
            required(&fields, "winner_index")?,
            required(&fields, "winner_penalty")?,
            required(&fields, "runner_up_penalty")?,
            required(&fields, "required_margin")?,
        ];
        let all_missing = selection_values.iter().all(|value| *value == ".");
        let any_missing = selection_values.iter().any(|value| *value == ".");
        let selection = if all_missing {
            None
        } else if any_missing {
            return Err("selection witness fields must be all present or all '.'".to_string());
        } else {
            Some(SelectionWitness {
                candidate_count: parse_u64(selection_values[0], "candidate_count")?,
                winner_index: parse_u64(selection_values[1], "winner_index")?,
                winner_penalty: parse_u64(selection_values[2], "winner_penalty")?,
                runner_up_penalty: parse_u64(selection_values[3], "runner_up_penalty")?,
                required_margin: parse_u64(selection_values[4], "required_margin")?,
            })
        };

        Ok(Self {
            target,
            assay: AssayProfile {
                kind,
                validated_target_enrichment,
            },
            backend,
            backend_version,
            status,
            call_count,
            input_sha256,
            resource_sha256,
            output_sha256,
            selection,
        })
    }

    /// Verify the finite target/assay contract and the certificate witnesses.
    ///
    /// This checks decision consistency and content-addressing syntax. It does
    /// not establish biological truth; that requires independently validated
    /// resources, assay performance, and truth data.
    pub fn verify(&self) -> Result<(), Vec<String>> {
        let mut errors = Vec::new();
        let target_spec = spec(self.target);

        if self.assay.kind == AssayKind::Wgs && self.assay.validated_target_enrichment {
            errors.push(
                "validated_target_enrichment must be false for a WGS certificate".to_string(),
            );
        }
        if self.backend != target_spec.backend.to_string() {
            errors.push(format!(
                "backend '{}' does not match registry backend '{}'",
                self.backend, target_spec.backend
            ));
        }
        if self.backend_version.is_empty()
            || self
                .backend_version
                .chars()
                .any(|character| matches!(character, '\n' | '\r' | '\t' | '='))
        {
            errors.push("backend_version is empty or contains a forbidden character".to_string());
        }

        for (name, value) in [
            ("input_sha256", self.input_sha256.as_str()),
            ("resource_sha256", self.resource_sha256.as_str()),
            ("output_sha256", self.output_sha256.as_str()),
        ] {
            if !is_sha256_hex(value) {
                errors.push(format!("{name} is not canonical lowercase SHA-256"));
            }
        }

        match self.status {
            CallStatus::Called => {
                if self.call_count == 0 {
                    errors.push("called certificate must have call_count > 0".to_string());
                }
                if target_spec.implementation == ImplementationStatus::ContractOnly {
                    errors.push(
                        "contract-only target cannot issue a called certificate".to_string(),
                    );
                }
                if let Err(reason) = require_observable(self.target, self.assay) {
                    errors.push(format!(
                        "called certificate violates assay observability: {reason}"
                    ));
                }
            }
            CallStatus::NoCall(reason) => {
                if self.call_count != 0 {
                    errors.push("no-call certificate must have call_count = 0".to_string());
                }
                let observability = require_observable(self.target, self.assay);
                match reason {
                    NoCallReason::AssayNotObservable
                        if observability != Err(NoCallReason::AssayNotObservable) =>
                    {
                        errors.push(
                            "assay-not-observable reason contradicts the registry".to_string(),
                        );
                    }
                    NoCallReason::MissingValidatedEnrichment
                        if observability != Err(NoCallReason::MissingValidatedEnrichment) =>
                    {
                        errors.push(
                            "missing-validated-enrichment reason contradicts the registry"
                                .to_string(),
                        );
                    }
                    _ => {}
                }
            }
        }

        match (self.target, self.status, self.selection) {
            (Target::Hba, CallStatus::Called, Some(witness)) => {
                verify_selection(witness, &mut errors);
                if self.call_count != 1 {
                    errors.push("native HBA selection must emit exactly one call".to_string());
                }
            }
            (Target::Hba, CallStatus::Called, None) => {
                errors.push("called HBA certificate is missing its selection witness".to_string());
            }
            (_, _, Some(_)) => {
                errors.push(
                    "selection witness is currently defined only for called HBA results"
                        .to_string(),
                );
            }
            _ => {}
        }

        if errors.is_empty() {
            Ok(())
        } else {
            Err(errors)
        }
    }
}

fn verify_selection(witness: SelectionWitness, errors: &mut Vec<String>) {
    if witness.candidate_count == 0 {
        errors.push("selection candidate_count must be > 0".to_string());
    }
    if witness.winner_index >= witness.candidate_count {
        errors.push("selection winner_index is outside candidate_count".to_string());
    }
    match witness.winner_penalty.checked_add(witness.required_margin) {
        Some(bound) if bound <= witness.runner_up_penalty => {}
        Some(_) => errors.push(
            "selection margin witness does not separate winner from runner-up".to_string(),
        ),
        None => errors.push("selection margin witness overflows u64".to_string()),
    }
}

fn required<'a>(fields: &'a BTreeMap<String, String>, key: &str) -> Result<&'a str, String> {
    fields
        .get(key)
        .map(String::as_str)
        .ok_or_else(|| format!("missing certificate key '{key}'"))
}

fn require_exact(
    fields: &BTreeMap<String, String>,
    key: &str,
    expected: &str,
) -> Result<(), String> {
    let actual = required(fields, key)?;
    if actual == expected {
        Ok(())
    } else {
        Err(format!(
            "certificate key '{key}' is '{actual}', expected '{expected}'"
        ))
    }
}

fn parse_bool(value: &str) -> Result<bool, String> {
    match value {
        "true" => Ok(true),
        "false" => Ok(false),
        _ => Err(format!("invalid boolean '{value}'")),
    }
}

fn parse_u64(value: &str, name: &str) -> Result<u64, String> {
    value
        .parse::<u64>()
        .map_err(|_| format!("invalid unsigned integer for {name}: '{value}'"))
}

impl fmt::Display for DecisionCertificate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.to_text())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::digest::{sha256_bytes, to_hex};

    fn digest(label: &[u8]) -> String {
        to_hex(&sha256_bytes(label))
    }

    #[test]
    fn valid_hba_certificate_round_trips() {
        let certificate = DecisionCertificate {
            target: Target::Hba,
            assay: AssayProfile::wgs(),
            backend: "native-hba-v1".to_string(),
            backend_version: "0.2.0".to_string(),
            status: CallStatus::Called,
            call_count: 1,
            input_sha256: digest(b"input"),
            resource_sha256: digest(b"resource"),
            output_sha256: digest(b"output"),
            selection: Some(SelectionWitness {
                candidate_count: 3,
                winner_index: 1,
                winner_penalty: 10,
                runner_up_penalty: 30,
                required_margin: 20,
            }),
        };
        assert_eq!(certificate.verify(), Ok(()));

        let parsed = DecisionCertificate::parse(&certificate.to_text()).unwrap();
        assert_eq!(parsed, certificate);
        assert_eq!(parsed.verify(), Ok(()));
    }

    #[test]
    fn selection_margin_is_checked() {
        let mut certificate = DecisionCertificate {
            target: Target::Hba,
            assay: AssayProfile::wgs(),
            backend: "native-hba-v1".to_string(),
            backend_version: "0.2.0".to_string(),
            status: CallStatus::Called,
            call_count: 1,
            input_sha256: digest(b"input"),
            resource_sha256: digest(b"resource"),
            output_sha256: digest(b"output"),
            selection: Some(SelectionWitness {
                candidate_count: 2,
                winner_index: 0,
                winner_penalty: 20,
                runner_up_penalty: 30,
                required_margin: 11,
            }),
        };
        assert!(certificate.verify().is_err());

        certificate.selection.as_mut().unwrap().required_margin = 10;
        assert_eq!(certificate.verify(), Ok(()));
    }

    #[test]
    fn called_hba_requires_enriched_wes() {
        let certificate = DecisionCertificate {
            target: Target::Hba,
            assay: AssayProfile::wes(false),
            backend: "native-hba-v1".to_string(),
            backend_version: "0.2.0".to_string(),
            status: CallStatus::Called,
            call_count: 1,
            input_sha256: digest(b"input"),
            resource_sha256: digest(b"resource"),
            output_sha256: digest(b"output"),
            selection: Some(SelectionWitness {
                candidate_count: 2,
                winner_index: 0,
                winner_penalty: 0,
                runner_up_penalty: 10,
                required_margin: 10,
            }),
        };
        assert!(certificate.verify().is_err());
    }
}
