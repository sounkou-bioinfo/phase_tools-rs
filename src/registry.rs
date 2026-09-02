use crate::model::{
    AssayKind, AssayProfile, Backend, CallerFamily, EvidenceKind, ImplementationStatus,
    NoCallReason, Target, WesSupport,
};

#[derive(Debug, Clone, Copy)]
pub struct TargetSpec {
    pub target: Target,
    pub family: CallerFamily,
    pub backend: Backend,
    pub implementation: ImplementationStatus,
    pub wes: WesSupport,
    pub evidence: &'static [EvidenceKind],
}

const HLA_KIR_EVIDENCE: &[EvidenceKind] = &[
    EvidenceKind::AlleleKmers,
    EvidenceKind::BandedAlignments,
    EvidenceKind::AlleleAbundance,
];

const PARALOG_EVIDENCE: &[EvidenceKind] = &[
    EvidenceKind::UniqueDepth,
    EvidenceKind::TotalDepth,
    EvidenceKind::ParalogSites,
    EvidenceKind::JunctionReads,
    EvidenceKind::SmallVariants,
    EvidenceKind::PhaseLinks,
];

const LPA_EVIDENCE: &[EvidenceKind] = &[
    EvidenceKind::UniqueDepth,
    EvidenceKind::RepeatSpanningReads,
    EvidenceKind::SmallVariants,
];

const RH_EVIDENCE: &[EvidenceKind] = &[
    EvidenceKind::UniqueDepth,
    EvidenceKind::TotalDepth,
    EvidenceKind::ParalogSites,
    EvidenceKind::JunctionReads,
    EvidenceKind::SmallVariants,
    EvidenceKind::PhaseLinks,
];

pub const TARGET_SPECS: [TargetSpec; 10] = [
    TargetSpec {
        target: Target::Hla,
        family: CallerFamily::AlleleTyping,
        backend: Backend::Unum,
        implementation: ImplementationStatus::Runnable,
        wes: WesSupport::StandardCapture,
        evidence: HLA_KIR_EVIDENCE,
    },
    TargetSpec {
        target: Target::Kir,
        family: CallerFamily::AlleleTyping,
        backend: Backend::Unum,
        implementation: ImplementationStatus::Runnable,
        wes: WesSupport::StandardCapture,
        evidence: HLA_KIR_EVIDENCE,
    },
    TargetSpec {
        target: Target::Cyp2b6,
        family: CallerFamily::ParalogCopyNumber,
        backend: Backend::NativePlanned,
        implementation: ImplementationStatus::ContractOnly,
        wes: WesSupport::Unsupported,
        evidence: PARALOG_EVIDENCE,
    },
    TargetSpec {
        target: Target::Cyp2d6,
        family: CallerFamily::ParalogCopyNumber,
        backend: Backend::NativePlanned,
        implementation: ImplementationStatus::ContractOnly,
        wes: WesSupport::Unsupported,
        evidence: PARALOG_EVIDENCE,
    },
    TargetSpec {
        target: Target::Cyp21a2,
        family: CallerFamily::ParalogCopyNumber,
        backend: Backend::NativePlanned,
        implementation: ImplementationStatus::ContractOnly,
        wes: WesSupport::Unsupported,
        evidence: PARALOG_EVIDENCE,
    },
    TargetSpec {
        target: Target::Gba,
        family: CallerFamily::ParalogCopyNumber,
        backend: Backend::NativePlanned,
        implementation: ImplementationStatus::ContractOnly,
        wes: WesSupport::Unsupported,
        evidence: PARALOG_EVIDENCE,
    },
    TargetSpec {
        target: Target::Hba,
        family: CallerFamily::ParalogCopyNumber,
        backend: Backend::NativeHba,
        implementation: ImplementationStatus::SolverKernel,
        wes: WesSupport::ValidatedTargetEnrichment,
        evidence: PARALOG_EVIDENCE,
    },
    TargetSpec {
        target: Target::Lpa,
        family: CallerFamily::RepeatCopyNumber,
        backend: Backend::NativePlanned,
        implementation: ImplementationStatus::ContractOnly,
        wes: WesSupport::Unsupported,
        evidence: LPA_EVIDENCE,
    },
    TargetSpec {
        target: Target::Rh,
        family: CallerFamily::BloodGroup,
        backend: Backend::NativePlanned,
        implementation: ImplementationStatus::ContractOnly,
        wes: WesSupport::Unsupported,
        evidence: RH_EVIDENCE,
    },
    TargetSpec {
        target: Target::Smn,
        family: CallerFamily::ParalogCopyNumber,
        backend: Backend::NativePlanned,
        implementation: ImplementationStatus::ContractOnly,
        wes: WesSupport::ValidatedTargetEnrichment,
        evidence: PARALOG_EVIDENCE,
    },
];

#[must_use]
pub const fn spec(target: Target) -> TargetSpec {
    match target {
        Target::Hla => TARGET_SPECS[0],
        Target::Kir => TARGET_SPECS[1],
        Target::Cyp2b6 => TARGET_SPECS[2],
        Target::Cyp2d6 => TARGET_SPECS[3],
        Target::Cyp21a2 => TARGET_SPECS[4],
        Target::Gba => TARGET_SPECS[5],
        Target::Hba => TARGET_SPECS[6],
        Target::Lpa => TARGET_SPECS[7],
        Target::Rh => TARGET_SPECS[8],
        Target::Smn => TARGET_SPECS[9],
    }
}

/// Return `Ok(())` only when the declared assay profile makes the target
/// observable enough to attempt. Evidence-level quality gates still apply.
pub fn require_observable(target: Target, assay: AssayProfile) -> Result<(), NoCallReason> {
    if assay.kind == AssayKind::Wgs {
        return Ok(());
    }

    match spec(target).wes {
        WesSupport::StandardCapture => Ok(()),
        WesSupport::ValidatedTargetEnrichment if assay.validated_target_enrichment => Ok(()),
        WesSupport::ValidatedTargetEnrichment => Err(NoCallReason::MissingValidatedEnrichment),
        WesSupport::Unsupported => Err(NoCallReason::AssayNotObservable),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn registry_is_exactly_the_closed_target_enum() {
        let registered = TARGET_SPECS
            .iter()
            .map(|entry| entry.target)
            .collect::<Vec<_>>();
        assert_eq!(registered.as_slice(), Target::ALL.as_slice());
    }

    #[test]
    fn dragen_targeted_set_is_closed_in_registry() {
        for target in Target::DRAGEN_45_TARGETED {
            assert_eq!(spec(target).target, target);
        }
    }

    #[test]
    fn all_targets_are_attemptable_on_wgs() {
        for target in Target::ALL {
            assert_eq!(require_observable(target, AssayProfile::wgs()), Ok(()));
        }
    }

    #[test]
    fn enriched_wes_expands_only_conditional_targets() {
        let plain = AssayProfile::wes(false);
        let enriched = AssayProfile::wes(true);

        assert_eq!(
            require_observable(Target::Hba, plain),
            Err(NoCallReason::MissingValidatedEnrichment)
        );
        assert_eq!(require_observable(Target::Hba, enriched), Ok(()));
        assert_eq!(
            require_observable(Target::Smn, plain),
            Err(NoCallReason::MissingValidatedEnrichment)
        );
        assert_eq!(require_observable(Target::Smn, enriched), Ok(()));

        assert_eq!(
            require_observable(Target::Cyp2d6, enriched),
            Err(NoCallReason::AssayNotObservable)
        );
    }
}
