use crate::model::NoCallReason;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::fs;
use std::path::Path;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HbaEvidence {
    pub values: BTreeMap<String, u64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FeatureExpectation {
    pub feature: String,
    pub expected: u64,
    pub tolerance: u64,
    pub weight: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HbaHypothesis {
    pub id: String,
    pub call: String,
    pub hba1_copy_number: u8,
    pub hba2_copy_number: u8,
    pub prior_penalty: u64,
    pub expectations: Vec<FeatureExpectation>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HbaOutcome {
    Called {
        winner_index: usize,
        winner_penalty: u64,
        runner_up_penalty: u64,
        margin: u64,
    },
    NoCall {
        reason: NoCallReason,
        top_penalty: Option<u64>,
        runner_up_penalty: Option<u64>,
        missing_features: Vec<String>,
    },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HbaDecision {
    pub penalties: Vec<u64>,
    pub outcome: HbaOutcome,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HbaError(String);

impl HbaError {
    fn new(message: impl Into<String>) -> Self {
        Self(message.into())
    }
}

impl fmt::Display for HbaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for HbaError {}

pub fn read_evidence(path: &Path) -> Result<HbaEvidence, HbaError> {
    let text = fs::read_to_string(path)
        .map_err(|error| HbaError::new(format!("cannot read '{}': {error}", path.display())))?;
    parse_evidence(&text)
}

pub fn parse_evidence(text: &str) -> Result<HbaEvidence, HbaError> {
    let mut values = BTreeMap::new();

    for (line_index, raw_line) in text.lines().enumerate() {
        let line_number = line_index + 1;
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() != 2 {
            return Err(HbaError::new(format!(
                "evidence line {line_number} must contain feature and value"
            )));
        }
        if fields[0].eq_ignore_ascii_case("feature")
            && fields[1].eq_ignore_ascii_case("value")
        {
            continue;
        }

        validate_identifier(fields[0], "feature", line_number)?;
        let value = parse_u64(fields[1], "evidence value", line_number)?;
        if values.insert(fields[0].to_string(), value).is_some() {
            return Err(HbaError::new(format!(
                "duplicate evidence feature '{}' on line {line_number}",
                fields[0]
            )));
        }
    }

    if values.is_empty() {
        return Err(HbaError::new("evidence contains no feature rows"));
    }

    Ok(HbaEvidence { values })
}

pub fn read_hypotheses(path: &Path) -> Result<Vec<HbaHypothesis>, HbaError> {
    let text = fs::read_to_string(path)
        .map_err(|error| HbaError::new(format!("cannot read '{}': {error}", path.display())))?;
    parse_hypotheses(&text)
}

/// Parse one row per hypothesis/feature expectation.
///
/// Columns:
/// `hypothesis, call, hba1_cn, hba2_cn, prior_penalty, feature, expected,
/// tolerance, weight`.
pub fn parse_hypotheses(text: &str) -> Result<Vec<HbaHypothesis>, HbaError> {
    let mut hypotheses = Vec::<HbaHypothesis>::new();
    let mut indexes = BTreeMap::<String, usize>::new();
    let mut seen_features = BTreeSet::<(String, String)>::new();

    for (line_index, raw_line) in text.lines().enumerate() {
        let line_number = line_index + 1;
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() != 9 {
            return Err(HbaError::new(format!(
                "hypothesis line {line_number} must contain 9 tab-separated columns"
            )));
        }
        if fields[0].eq_ignore_ascii_case("hypothesis")
            && fields[5].eq_ignore_ascii_case("feature")
        {
            continue;
        }

        validate_identifier(fields[0], "hypothesis", line_number)?;
        if fields[1].is_empty() {
            return Err(HbaError::new(format!(
                "empty HBA call on line {line_number}"
            )));
        }
        validate_identifier(fields[5], "feature", line_number)?;

        let hba1_copy_number = parse_u8(fields[2], "hba1_cn", line_number)?;
        let hba2_copy_number = parse_u8(fields[3], "hba2_cn", line_number)?;
        let prior_penalty = parse_u64(fields[4], "prior_penalty", line_number)?;
        let expected = parse_u64(fields[6], "expected", line_number)?;
        let tolerance = parse_u64(fields[7], "tolerance", line_number)?;
        let weight = parse_u64(fields[8], "weight", line_number)?;
        if tolerance == 0 {
            return Err(HbaError::new(format!(
                "tolerance must be > 0 on line {line_number}"
            )));
        }
        if weight == 0 {
            return Err(HbaError::new(format!(
                "weight must be > 0 on line {line_number}"
            )));
        }

        let feature_key = (fields[0].to_string(), fields[5].to_string());
        if !seen_features.insert(feature_key) {
            return Err(HbaError::new(format!(
                "duplicate feature '{}' for hypothesis '{}' on line {line_number}",
                fields[5], fields[0]
            )));
        }

        let expectation = FeatureExpectation {
            feature: fields[5].to_string(),
            expected,
            tolerance,
            weight,
        };

        if let Some(&index) = indexes.get(fields[0]) {
            let hypothesis = &mut hypotheses[index];
            if hypothesis.call != fields[1]
                || hypothesis.hba1_copy_number != hba1_copy_number
                || hypothesis.hba2_copy_number != hba2_copy_number
                || hypothesis.prior_penalty != prior_penalty
            {
                return Err(HbaError::new(format!(
                    "inconsistent metadata for hypothesis '{}' on line {line_number}",
                    fields[0]
                )));
            }
            hypothesis.expectations.push(expectation);
        } else {
            let index = hypotheses.len();
            indexes.insert(fields[0].to_string(), index);
            hypotheses.push(HbaHypothesis {
                id: fields[0].to_string(),
                call: fields[1].to_string(),
                hba1_copy_number,
                hba2_copy_number,
                prior_penalty,
                expectations: vec![expectation],
            });
        }
    }

    if hypotheses.len() < 2 {
        return Err(HbaError::new(
            "HBA hypothesis resource must contain at least two hypotheses",
        ));
    }

    Ok(hypotheses)
}

/// Evaluate a resource-defined HBA hypothesis catalogue using integer penalties.
///
/// For each feature, the contribution is
/// `ceil(abs(observed - expected) / tolerance) * weight`. Lower is better.
/// Integer arithmetic makes the exact decision portable and certifiable.
pub fn select_hba(
    evidence: &HbaEvidence,
    hypotheses: &[HbaHypothesis],
    required_margin: u64,
) -> Result<HbaDecision, HbaError> {
    if hypotheses.len() < 2 {
        return Err(HbaError::new(
            "HBA selection requires at least two hypotheses",
        ));
    }

    let mut required_features = BTreeSet::new();
    for hypothesis in hypotheses {
        if hypothesis.expectations.is_empty() {
            return Err(HbaError::new(format!(
                "hypothesis '{}' has no feature expectations",
                hypothesis.id
            )));
        }
        for expectation in &hypothesis.expectations {
            required_features.insert(expectation.feature.clone());
        }
    }

    let missing_features = required_features
        .into_iter()
        .filter(|feature| !evidence.values.contains_key(feature))
        .collect::<Vec<_>>();
    if !missing_features.is_empty() {
        return Ok(HbaDecision {
            penalties: Vec::new(),
            outcome: HbaOutcome::NoCall {
                reason: NoCallReason::InsufficientEvidence,
                top_penalty: None,
                runner_up_penalty: None,
                missing_features,
            },
        });
    }

    let penalties = hypotheses
        .iter()
        .map(|hypothesis| score_hypothesis(evidence, hypothesis))
        .collect::<Result<Vec<_>, _>>()?;

    let winner_penalty = *penalties
        .iter()
        .min()
        .ok_or_else(|| HbaError::new("HBA selection received no hypotheses"))?;
    let winners = penalties
        .iter()
        .enumerate()
        .filter(|(_, penalty)| **penalty == winner_penalty)
        .map(|(index, _)| index)
        .collect::<Vec<_>>();

    if winners.len() != 1 {
        return Ok(HbaDecision {
            penalties,
            outcome: HbaOutcome::NoCall {
                reason: NoCallReason::AmbiguousTopScore,
                top_penalty: Some(winner_penalty),
                runner_up_penalty: Some(winner_penalty),
                missing_features: Vec::new(),
            },
        });
    }

    let winner_index = winners[0];
    let runner_up_penalty = penalties
        .iter()
        .enumerate()
        .filter(|(index, _)| *index != winner_index)
        .map(|(_, penalty)| *penalty)
        .min()
        .ok_or_else(|| HbaError::new("HBA selection has no runner-up hypothesis"))?;
    let margin = runner_up_penalty
        .checked_sub(winner_penalty)
        .ok_or_else(|| HbaError::new("internal HBA penalty ordering failure"))?;

    if margin < required_margin {
        return Ok(HbaDecision {
            penalties,
            outcome: HbaOutcome::NoCall {
                reason: NoCallReason::AmbiguousTopScore,
                top_penalty: Some(winner_penalty),
                runner_up_penalty: Some(runner_up_penalty),
                missing_features: Vec::new(),
            },
        });
    }

    Ok(HbaDecision {
        penalties,
        outcome: HbaOutcome::Called {
            winner_index,
            winner_penalty,
            runner_up_penalty,
            margin,
        },
    })
}

fn score_hypothesis(
    evidence: &HbaEvidence,
    hypothesis: &HbaHypothesis,
) -> Result<u64, HbaError> {
    let mut penalty = hypothesis.prior_penalty;

    for expectation in &hypothesis.expectations {
        let observed = evidence
            .values
            .get(&expectation.feature)
            .copied()
            .ok_or_else(|| {
                HbaError::new(format!(
                    "missing evidence feature '{}'",
                    expectation.feature
                ))
            })?;
        let delta = observed.abs_diff(expectation.expected);
        let remainder = if delta % expectation.tolerance == 0 {
            0
        } else {
            1
        };
        let units = delta / expectation.tolerance + remainder;
        let contribution = units.checked_mul(expectation.weight).ok_or_else(|| {
            HbaError::new(format!(
                "penalty overflow for hypothesis '{}' feature '{}'",
                hypothesis.id, expectation.feature
            ))
        })?;
        penalty = penalty.checked_add(contribution).ok_or_else(|| {
            HbaError::new(format!("penalty overflow for hypothesis '{}'", hypothesis.id))
        })?;
    }

    Ok(penalty)
}

fn validate_identifier(value: &str, name: &str, line_number: usize) -> Result<(), HbaError> {
    if value.is_empty()
        || value
            .bytes()
            .any(|byte| byte.is_ascii_whitespace() || byte == b'=')
    {
        Err(HbaError::new(format!(
            "invalid {name} '{value}' on line {line_number}"
        )))
    } else {
        Ok(())
    }
}

fn parse_u64(value: &str, name: &str, line_number: usize) -> Result<u64, HbaError> {
    value.parse::<u64>().map_err(|_| {
        HbaError::new(format!(
            "invalid unsigned integer for {name} on line {line_number}: '{value}'"
        ))
    })
}

fn parse_u8(value: &str, name: &str, line_number: usize) -> Result<u8, HbaError> {
    value.parse::<u8>().map_err(|_| {
        HbaError::new(format!(
            "invalid copy number for {name} on line {line_number}: '{value}'"
        ))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const HYPOTHESES: &str = "\
hypothesis\tcall\thba1_cn\thba2_cn\tprior_penalty\tfeature\texpected\ttolerance\tweight
normal\talpha-alpha/alpha-alpha\t2\t2\t0\thba_total_depth_milli\t4000\t100\t10
normal\talpha-alpha/alpha-alpha\t2\t2\t0\thba1_unique_depth_milli\t2000\t100\t10
single_deletion\t-alpha/alpha-alpha\t1\t2\t0\thba_total_depth_milli\t3000\t100\t10
single_deletion\t-alpha/alpha-alpha\t1\t2\t0\thba1_unique_depth_milli\t1000\t100\t10
";

    #[test]
    fn exact_normal_evidence_selects_normal() {
        let evidence = parse_evidence(
            "feature\tvalue\nhba_total_depth_milli\t4000\nhba1_unique_depth_milli\t2000\n",
        )
        .unwrap();
        let hypotheses = parse_hypotheses(HYPOTHESES).unwrap();
        let decision = select_hba(&evidence, &hypotheses, 10).unwrap();

        assert_eq!(
            decision.outcome,
            HbaOutcome::Called {
                winner_index: 0,
                winner_penalty: 0,
                runner_up_penalty: 200,
                margin: 200,
            }
        );
    }

    #[test]
    fn ties_are_no_calls() {
        let evidence = parse_evidence(
            "hba_total_depth_milli\t3500\nhba1_unique_depth_milli\t1500\n",
        )
        .unwrap();
        let hypotheses = parse_hypotheses(HYPOTHESES).unwrap();
        let decision = select_hba(&evidence, &hypotheses, 1).unwrap();

        assert!(matches!(
            decision.outcome,
            HbaOutcome::NoCall {
                reason: NoCallReason::AmbiguousTopScore,
                ..
            }
        ));
    }

    #[test]
    fn missing_features_are_explicit_no_calls() {
        let evidence = parse_evidence("hba_total_depth_milli\t4000\n").unwrap();
        let hypotheses = parse_hypotheses(HYPOTHESES).unwrap();
        let decision = select_hba(&evidence, &hypotheses, 10).unwrap();

        assert_eq!(
            decision.outcome,
            HbaOutcome::NoCall {
                reason: NoCallReason::InsufficientEvidence,
                top_penalty: None,
                runner_up_penalty: None,
                missing_features: vec!["hba1_unique_depth_milli".to_string()],
            }
        );
    }

    #[test]
    fn hypothesis_metadata_must_be_consistent() {
        let bad = HYPOTHESES.replace(
            "normal\talpha-alpha/alpha-alpha\t2\t2\t0\thba1",
            "normal\tbroken\t2\t2\t0\thba1",
        );
        assert!(parse_hypotheses(&bad).is_err());
    }
}
