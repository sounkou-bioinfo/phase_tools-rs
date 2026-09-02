use std::fmt;
use std::str::FromStr;

/// Coarse sequencing-assay class.
///
/// Observability is decided jointly from this class, the target contract, and
/// whether a validated target-enrichment profile was declared.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssayKind {
    Wgs,
    Wes,
}

impl fmt::Display for AssayKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Wgs => f.write_str("wgs"),
            Self::Wes => f.write_str("wes"),
        }
    }
}

impl FromStr for AssayKind {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.to_ascii_lowercase().as_str() {
            "wgs" | "genome" => Ok(Self::Wgs),
            "wes" | "exome" => Ok(Self::Wes),
            _ => Err(format!("unknown assay '{value}'; expected wgs or wes")),
        }
    }
}

/// Assay facts that are relevant to difficult-locus observability.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AssayProfile {
    pub kind: AssayKind,
    pub validated_target_enrichment: bool,
}

impl AssayProfile {
    #[must_use]
    pub const fn wgs() -> Self {
        Self {
            kind: AssayKind::Wgs,
            validated_target_enrichment: false,
        }
    }

    #[must_use]
    pub const fn wes(validated_target_enrichment: bool) -> Self {
        Self {
            kind: AssayKind::Wes,
            validated_target_enrichment,
        }
    }
}

/// Closed difficult-locus scope.
///
/// `ALL` is deliberately exhaustive: adding a target requires touching every
/// target-dependent match and the Lean registry proof.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Target {
    Hla,
    Kir,
    Cyp2b6,
    Cyp2d6,
    Cyp21a2,
    Gba,
    Hba,
    Lpa,
    Rh,
    Smn,
}

impl Target {
    pub const ALL: [Self; 10] = [
        Self::Hla,
        Self::Kir,
        Self::Cyp2b6,
        Self::Cyp2d6,
        Self::Cyp21a2,
        Self::Gba,
        Self::Hba,
        Self::Lpa,
        Self::Rh,
        Self::Smn,
    ];

    /// Illumina DRAGEN v4.5 Targeted Caller target names.
    ///
    /// HLA is a separate DRAGEN caller and KIR is supplied here through the
    /// Unum/T1K allele-typing lane.
    pub const DRAGEN_45_TARGETED: [Self; 8] = [
        Self::Cyp2b6,
        Self::Cyp2d6,
        Self::Cyp21a2,
        Self::Gba,
        Self::Hba,
        Self::Lpa,
        Self::Rh,
        Self::Smn,
    ];
}

impl fmt::Display for Target {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Hla => "HLA",
            Self::Kir => "KIR",
            Self::Cyp2b6 => "CYP2B6",
            Self::Cyp2d6 => "CYP2D6",
            Self::Cyp21a2 => "CYP21A2",
            Self::Gba => "GBA",
            Self::Hba => "HBA",
            Self::Lpa => "LPA",
            Self::Rh => "RH",
            Self::Smn => "SMN",
        })
    }
}

impl FromStr for Target {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let normalized = value
            .trim()
            .to_ascii_lowercase()
            .chars()
            .filter(|character| !matches!(character, '-' | '_'))
            .collect::<String>();
        match normalized.as_str() {
            "hla" => Ok(Self::Hla),
            "kir" => Ok(Self::Kir),
            "cyp2b6" => Ok(Self::Cyp2b6),
            "cyp2d6" => Ok(Self::Cyp2d6),
            "cyp21a2" => Ok(Self::Cyp21a2),
            "gba" => Ok(Self::Gba),
            "hba" | "alpha" | "alphaglobin" => Ok(Self::Hba),
            "lpa" => Ok(Self::Lpa),
            "rh" | "rhd" | "rhce" => Ok(Self::Rh),
            "smn" | "smn1" | "smn2" => Ok(Self::Smn),
            _ => Err(format!("unknown target '{value}'")),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CallerFamily {
    AlleleTyping,
    ParalogCopyNumber,
    RepeatCopyNumber,
    BloodGroup,
}

impl fmt::Display for CallerFamily {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::AlleleTyping => "allele-typing",
            Self::ParalogCopyNumber => "paralog-copy-number",
            Self::RepeatCopyNumber => "repeat-copy-number",
            Self::BloodGroup => "blood-group",
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Backend {
    Unum,
    NativeHba,
    NativePlanned,
}

impl fmt::Display for Backend {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Unum => "unum",
            Self::NativeHba => "native-hba-v1",
            Self::NativePlanned => "native-planned",
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ImplementationStatus {
    Runnable,
    SolverKernel,
    ContractOnly,
}

impl fmt::Display for ImplementationStatus {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Runnable => "runnable",
            Self::SolverKernel => "solver-kernel",
            Self::ContractOnly => "contract-only",
        })
    }
}

/// What a WES profile must provide before a target may be attempted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WesSupport {
    StandardCapture,
    ValidatedTargetEnrichment,
    Unsupported,
}

impl fmt::Display for WesSupport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::StandardCapture => "standard-capture",
            Self::ValidatedTargetEnrichment => "validated-enrichment",
            Self::Unsupported => "unsupported",
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvidenceKind {
    AlleleKmers,
    BandedAlignments,
    AlleleAbundance,
    UniqueDepth,
    TotalDepth,
    ParalogSites,
    JunctionReads,
    SmallVariants,
    RepeatSpanningReads,
    PhaseLinks,
}

impl fmt::Display for EvidenceKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::AlleleKmers => "allele-kmers",
            Self::BandedAlignments => "banded-alignments",
            Self::AlleleAbundance => "allele-abundance",
            Self::UniqueDepth => "unique-depth",
            Self::TotalDepth => "total-depth",
            Self::ParalogSites => "paralog-sites",
            Self::JunctionReads => "junction-reads",
            Self::SmallVariants => "small-variants",
            Self::RepeatSpanningReads => "repeat-spanning-reads",
            Self::PhaseLinks => "phase-links",
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NoCallReason {
    AssayNotObservable,
    MissingValidatedEnrichment,
    InsufficientEvidence,
    AmbiguousTopScore,
    BackendFailure,
    ResourceMismatch,
}

impl fmt::Display for NoCallReason {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::AssayNotObservable => "assay-not-observable",
            Self::MissingValidatedEnrichment => "missing-validated-enrichment",
            Self::InsufficientEvidence => "insufficient-evidence",
            Self::AmbiguousTopScore => "ambiguous-top-score",
            Self::BackendFailure => "backend-failure",
            Self::ResourceMismatch => "resource-mismatch",
        })
    }
}

impl FromStr for NoCallReason {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value {
            "assay-not-observable" => Ok(Self::AssayNotObservable),
            "missing-validated-enrichment" => Ok(Self::MissingValidatedEnrichment),
            "insufficient-evidence" => Ok(Self::InsufficientEvidence),
            "ambiguous-top-score" => Ok(Self::AmbiguousTopScore),
            "backend-failure" => Ok(Self::BackendFailure),
            "resource-mismatch" => Ok(Self::ResourceMismatch),
            _ => Err(format!("unknown no-call reason '{value}'")),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CallStatus {
    Called,
    NoCall(NoCallReason),
}
