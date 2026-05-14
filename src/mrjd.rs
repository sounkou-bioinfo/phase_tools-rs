//! Initial multi-region joint-detection kernels.
//!
//! This module is a small, auditable foundation for high-homology / multi-region
//! candidate evidence. It is not a DRAGEN-equivalent caller: it groups candidate
//! SNV evidence by region-group offset so downstream code can inspect homologous
//! support across multiple loci before a full posterior haplotype model exists.

use crate::io::fasta::Fai;
use rust_htslib::bam::{self, Read};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Region {
    pub group: String,
    pub chrom: String,
    pub start1: i64,
    pub end1: i64,
    pub copy: String,
}

impl Region {
    pub fn len(&self) -> i64 {
        self.end1 - self.start1 + 1
    }

    fn contains_pos0(&self, pos0: i64) -> bool {
        (self.start1 - 1) <= pos0 && pos0 < self.end1
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CandidateConfig {
    pub min_mapq: u8,
    pub min_baseq: u8,
    pub min_alt_count: u32,
    pub min_alt_fraction: f64,
}

impl Default for CandidateConfig {
    fn default() -> Self {
        Self {
            min_mapq: 0,
            min_baseq: 13,
            min_alt_count: 2,
            min_alt_fraction: 0.20,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mapq255Policy {
    DropWhenFiltering,
    KeepWhenFiltering,
}

impl Default for Mapq255Policy {
    fn default() -> Self {
        Self::DropWhenFiltering
    }
}

impl Mapq255Policy {
    fn keep_when_filtering(self) -> bool {
        matches!(self, Self::KeepWhenFiltering)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct RegionObservation {
    pub copy: String,
    pub chrom: String,
    pub pos1: i64,
    pub ref_base: u8,
    pub depth: u32,
    pub alt_count: u32,
    pub alt_fraction: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct JointCandidate {
    pub group: String,
    pub offset1: i64,
    pub alt_base: u8,
    pub alt_positive_depth: u32,
    pub alt_positive_alt_count: u32,
    pub regions_with_alt: usize,
    pub region_count: usize,
    pub observations: Vec<RegionObservation>,
}

fn is_regions_header(fields: &[&str]) -> bool {
    fields.len() >= 4
        && fields[0].eq_ignore_ascii_case("group")
        && fields[1].eq_ignore_ascii_case("chrom")
        && fields[2].eq_ignore_ascii_case("start")
        && fields[3].eq_ignore_ascii_case("end")
        && fields
            .get(4)
            .map_or(true, |value| value.eq_ignore_ascii_case("copy"))
}

pub fn read_regions_tsv(path: &str) -> Result<Vec<Region>, String> {
    let file = File::open(path).map_err(|e| format!("cannot open regions TSV '{path}': {e}"))?;
    let reader = BufReader::new(file);
    let mut regions = Vec::new();
    let mut first_record = true;
    for (line_no, line) in reader.lines().enumerate() {
        let line_no = line_no + 1;
        let line = line.map_err(|e| format!("failed reading regions TSV '{path}': {e}"))?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields = trimmed.split('\t').collect::<Vec<_>>();
        if first_record && is_regions_header(&fields) {
            first_record = false;
            continue;
        }
        first_record = false;
        if fields.len() < 4 || fields.len() > 5 {
            return Err(format!(
                "regions TSV line {line_no} must have 4 or 5 tab-separated fields: group chrom start end [copy]"
            ));
        }
        let group = fields[0].to_string();
        let chrom = fields[1].to_string();
        let start1 = fields[2]
            .parse::<i64>()
            .map_err(|_| format!("regions TSV line {line_no} has invalid start"))?;
        let end1 = fields[3]
            .parse::<i64>()
            .map_err(|_| format!("regions TSV line {line_no} has invalid end"))?;
        if group.is_empty() || chrom.is_empty() {
            return Err(format!(
                "regions TSV line {line_no} has empty group or chrom"
            ));
        }
        if start1 < 1 || end1 < start1 {
            return Err(format!(
                "regions TSV line {line_no} has invalid interval {chrom}:{start1}-{end1}"
            ));
        }
        let copy = fields
            .get(4)
            .filter(|value| !value.is_empty())
            .map(|value| (*value).to_string())
            .unwrap_or_else(|| format!("{}:{}-{}", chrom, start1, end1));
        regions.push(Region {
            group,
            chrom,
            start1,
            end1,
            copy,
        });
    }
    if regions.is_empty() {
        return Err(format!("regions TSV '{path}' did not contain any regions"));
    }
    Ok(regions)
}

pub fn validate_candidate_config(cfg: CandidateConfig) -> Result<(), String> {
    if cfg.min_alt_count == 0 {
        return Err("min_alt_count must be >= 1".to_string());
    }
    if !(0.0..=1.0).contains(&cfg.min_alt_fraction) {
        return Err("min_alt_fraction must be between 0 and 1".to_string());
    }
    Ok(())
}

fn base_index(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn base_from_index(index: usize) -> u8 {
    [b'A', b'C', b'G', b'T'][index]
}

fn region_count_for_offset(regions: &[Region], group: &str, offset1: i64) -> usize {
    regions
        .iter()
        .filter(|region| region.group == group && offset1 >= 1 && offset1 <= region.len())
        .count()
}

fn usable_record(record: &bam::Record, min_mapq: u8, keep_mapq_255: bool) -> bool {
    if record.is_unmapped()
        || record.is_quality_check_failed()
        || record.is_duplicate()
        || record.is_secondary()
        || record.is_supplementary()
    {
        return false;
    }
    if min_mapq == 0 {
        return true;
    }
    let mapq = record.mapq();
    if mapq == 255 {
        return keep_mapq_255;
    }
    mapq >= min_mapq
}

pub fn detect_snv_candidates(
    bam_path: &str,
    reference_path: &str,
    fai: &Fai,
    regions: &[Region],
    cfg: CandidateConfig,
    threads: usize,
) -> Result<Vec<JointCandidate>, String> {
    detect_snv_candidates_with_mapq255_policy(
        bam_path,
        reference_path,
        fai,
        regions,
        cfg,
        Mapq255Policy::default(),
        threads,
    )
}

pub fn detect_snv_candidates_with_mapq255_policy(
    bam_path: &str,
    reference_path: &str,
    fai: &Fai,
    regions: &[Region],
    cfg: CandidateConfig,
    mapq255_policy: Mapq255Policy,
    threads: usize,
) -> Result<Vec<JointCandidate>, String> {
    validate_candidate_config(cfg)?;
    if regions.is_empty() {
        return Ok(Vec::new());
    }
    let mut bam = bam::IndexedReader::from_path(bam_path)
        .map_err(|e| format!("cannot open indexed BAM/CRAM '{bam_path}': {e}"))?;
    bam.set_reference(reference_path)
        .map_err(|e| format!("failed to set CRAM reference '{reference_path}': {e}"))?;
    if threads > 1 {
        bam.set_threads(threads)
            .map_err(|e| format!("failed to enable BAM/CRAM threads for '{bam_path}': {e}"))?;
    }

    let mut grouped: BTreeMap<(String, i64, u8), Vec<RegionObservation>> = BTreeMap::new();
    for region in regions {
        let ref_seq = fai.fetch_seq(&region.chrom, region.start1, region.end1)?;
        if ref_seq.len() != region.len() as usize {
            return Err(format!(
                "reference interval {}:{}-{} returned length {}, expected {}",
                region.chrom,
                region.start1,
                region.end1,
                ref_seq.len(),
                region.len()
            ));
        }
        let start0 = region.start1 - 1;
        let end0 = region.end1;
        bam.fetch((region.chrom.as_bytes(), start0, end0))
            .map_err(|e| {
                format!(
                    "failed to fetch {}:{}-{} from '{bam_path}': {e}",
                    region.chrom, region.start1, region.end1
                )
            })?;
        for pileup in bam.pileup() {
            let pileup =
                pileup.map_err(|e| format!("failed to read pileup from '{bam_path}': {e}"))?;
            let pos0 = pileup.pos() as i64;
            if !region.contains_pos0(pos0) {
                continue;
            }
            let offset0 = pos0 - start0;
            let Some(&ref_base) = ref_seq.get(offset0 as usize) else {
                continue;
            };
            let Some(ref_index) = base_index(ref_base) else {
                continue;
            };
            let mut counts = [0u32; 4];
            for alignment in pileup.alignments() {
                let record = alignment.record();
                if !usable_record(&record, cfg.min_mapq, mapq255_policy.keep_when_filtering()) {
                    continue;
                }
                let Some(qpos) = alignment.qpos() else {
                    continue;
                };
                if record.qual().get(qpos).copied().unwrap_or(0) < cfg.min_baseq {
                    continue;
                }
                let base = record.seq()[qpos];
                if let Some(index) = base_index(base) {
                    counts[index] += 1;
                }
            }
            let depth = counts.iter().sum::<u32>();
            if depth == 0 {
                continue;
            }
            for (index, &alt_count) in counts.iter().enumerate() {
                if index == ref_index || alt_count < cfg.min_alt_count {
                    continue;
                }
                let alt_fraction = alt_count as f64 / depth as f64;
                if alt_fraction < cfg.min_alt_fraction {
                    continue;
                }
                let offset1 = offset0 + 1;
                grouped
                    .entry((region.group.clone(), offset1, base_from_index(index)))
                    .or_default()
                    .push(RegionObservation {
                        copy: region.copy.clone(),
                        chrom: region.chrom.clone(),
                        pos1: pos0 + 1,
                        ref_base: ref_base.to_ascii_uppercase(),
                        depth,
                        alt_count,
                        alt_fraction,
                    });
            }
        }
    }

    let mut out = Vec::new();
    for ((group, offset1, alt_base), observations) in grouped {
        let alt_positive_depth = observations.iter().map(|obs| obs.depth).sum();
        let alt_positive_alt_count = observations.iter().map(|obs| obs.alt_count).sum();
        let region_count = region_count_for_offset(regions, &group, offset1);
        out.push(JointCandidate {
            group,
            offset1,
            alt_base,
            alt_positive_depth,
            alt_positive_alt_count,
            regions_with_alt: observations.len(),
            region_count,
            observations,
        });
    }
    Ok(out)
}

fn format_observations(candidate: &JointCandidate) -> String {
    candidate
        .observations
        .iter()
        .map(|obs| {
            format!(
                "{}|{}:{}|{}|{}|{}|{:.6}",
                obs.copy,
                obs.chrom,
                obs.pos1,
                obs.ref_base as char,
                obs.depth,
                obs.alt_count,
                obs.alt_fraction
            )
        })
        .collect::<Vec<_>>()
        .join(";")
}

pub fn write_candidates_tsv<W: Write>(
    mut out: W,
    candidates: &[JointCandidate],
) -> Result<(), String> {
    writeln!(
        out,
        "group\toffset1\talt\talt_positive_depth\talt_positive_alt_count\tregions_with_alt\tregion_count\tper_region"
    )
    .map_err(|e| format!("failed to write output: {e}"))?;
    for candidate in candidates {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            candidate.group,
            candidate.offset1,
            candidate.alt_base as char,
            candidate.alt_positive_depth,
            candidate.alt_positive_alt_count,
            candidate.regions_with_alt,
            candidate.region_count,
            format_observations(candidate)
        )
        .map_err(|e| format!("failed to write output: {e}"))?;
    }
    Ok(())
}

fn vcf_escape(value: &str) -> String {
    let mut out = String::new();
    for byte in value.bytes() {
        match byte {
            b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' | b'.' | b'_' | b'-' | b':' => {
                out.push(byte as char);
            }
            _ => out.push_str(&format!("%{byte:02X}")),
        }
    }
    out
}

fn diagnostic_event_id(candidate: &JointCandidate) -> String {
    format!(
        "MRJD:{}:{}:{}",
        vcf_escape(&candidate.group),
        candidate.offset1,
        candidate.alt_base as char
    )
}

pub fn write_diagnostic_vcf<W: Write>(
    mut out: W,
    regions: &[Region],
    candidates: &[JointCandidate],
) -> Result<(), String> {
    writeln!(out, "##fileformat=VCFv4.3").map_err(|e| format!("failed to write VCF: {e}"))?;
    writeln!(out, "##source=phase_tools-rs_multi_region_joint_detect")
        .map_err(|e| format!("failed to write VCF: {e}"))?;
    let contigs = regions
        .iter()
        .map(|region| region.chrom.as_str())
        .collect::<BTreeSet<_>>();
    for contig in contigs {
        writeln!(out, "##contig=<ID={}>", vcf_escape(contig))
            .map_err(|e| format!("failed to write VCF: {e}"))?;
    }
    writeln!(out, "##INFO=<ID=EVENT,Number=1,Type=String,Description=\"Multi-region joint-detection event identifier shared by homologous observations\">").map_err(|e| format!("failed to write VCF: {e}"))?;
    writeln!(out, "##INFO=<ID=EVENTTYPE,Number=1,Type=String,Description=\"Event class; currently MRJD_SNV for diagnostic SNV evidence\">").map_err(|e| format!("failed to write VCF: {e}"))?;
    writeln!(out, "##INFO=<ID=MRJD_GROUP,Number=1,Type=String,Description=\"Region group from the multi-region manifest\">").map_err(|e| format!("failed to write VCF: {e}"))?;
    writeln!(out, "##INFO=<ID=MRJD_COPY,Number=1,Type=String,Description=\"Copy/region label from the multi-region manifest\">").map_err(|e| format!("failed to write VCF: {e}"))?;
    writeln!(out, "##INFO=<ID=MRJD_OFFSET,Number=1,Type=Integer,Description=\"1-based offset within the grouped region\">").map_err(|e| format!("failed to write VCF: {e}"))?;
    writeln!(out, "##INFO=<ID=MRJD_REGIONS_WITH_ALT,Number=1,Type=Integer,Description=\"Number of regions in this group with alt-positive evidence for the event\">").map_err(|e| format!("failed to write VCF: {e}"))?;
    writeln!(out, "##INFO=<ID=MRJD_REGION_COUNT,Number=1,Type=Integer,Description=\"Number of manifest regions in this group covering the event offset\">").map_err(|e| format!("failed to write VCF: {e}"))?;
    writeln!(out, "##INFO=<ID=MRJD_ALT_POSITIVE_DEPTH,Number=1,Type=Integer,Description=\"Summed depth across alt-positive region observations only\">").map_err(|e| format!("failed to write VCF: {e}"))?;
    writeln!(out, "##INFO=<ID=MRJD_ALT_POSITIVE_ALT_COUNT,Number=1,Type=Integer,Description=\"Summed alt count across alt-positive region observations only\">").map_err(|e| format!("failed to write VCF: {e}"))?;
    writeln!(out, "##INFO=<ID=MRJD_COPY_DEPTH,Number=1,Type=Integer,Description=\"Depth for this specific region observation\">").map_err(|e| format!("failed to write VCF: {e}"))?;
    writeln!(out, "##INFO=<ID=MRJD_COPY_ALT_COUNT,Number=1,Type=Integer,Description=\"Alt count for this specific region observation\">").map_err(|e| format!("failed to write VCF: {e}"))?;
    writeln!(out, "##INFO=<ID=MRJD_COPY_ALT_FRACTION,Number=1,Type=Float,Description=\"Alt fraction for this specific region observation\">").map_err(|e| format!("failed to write VCF: {e}"))?;
    writeln!(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
        .map_err(|e| format!("failed to write VCF: {e}"))?;

    let mut rows = Vec::new();
    for candidate in candidates {
        let event = diagnostic_event_id(candidate);
        for obs in &candidate.observations {
            rows.push((
                obs.chrom.clone(),
                obs.pos1,
                obs.ref_base as char,
                candidate.alt_base as char,
                event.clone(),
                candidate,
                obs,
            ));
        }
    }
    rows.sort_by(|a, b| {
        a.0.cmp(&b.0)
            .then_with(|| a.1.cmp(&b.1))
            .then_with(|| a.3.cmp(&b.3))
            .then_with(|| a.5.group.cmp(&b.5.group))
            .then_with(|| a.5.offset1.cmp(&b.5.offset1))
            .then_with(|| a.6.copy.cmp(&b.6.copy))
    });
    for (chrom, pos1, ref_base, alt_base, event, candidate, obs) in rows {
        let info = format!(
            "EVENT={};EVENTTYPE=MRJD_SNV;MRJD_GROUP={};MRJD_COPY={};MRJD_OFFSET={};MRJD_REGIONS_WITH_ALT={};MRJD_REGION_COUNT={};MRJD_ALT_POSITIVE_DEPTH={};MRJD_ALT_POSITIVE_ALT_COUNT={};MRJD_COPY_DEPTH={};MRJD_COPY_ALT_COUNT={};MRJD_COPY_ALT_FRACTION={:.6}",
            event,
            vcf_escape(&candidate.group),
            vcf_escape(&obs.copy),
            candidate.offset1,
            candidate.regions_with_alt,
            candidate.region_count,
            candidate.alt_positive_depth,
            candidate.alt_positive_alt_count,
            obs.depth,
            obs.alt_count,
            obs.alt_fraction
        );
        writeln!(
            out,
            "{}\t{}\t.\t{}\t{}\t.\tPASS\t{}",
            vcf_escape(&chrom),
            pos1,
            ref_base,
            alt_base,
            info
        )
        .map_err(|e| format!("failed to write VCF: {e}"))?;
    }
    out.flush()
        .map_err(|e| format!("failed to write VCF: {e}"))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_region_manifest() {
        let dir = std::env::temp_dir();
        let path = dir.join(format!("phase_tools_regions_{}.tsv", std::process::id()));
        std::fs::write(
            &path,
            "group\tchrom\tstart\tend\tcopy\nG1\tchr1\t10\t20\tcopy1\n",
        )
        .unwrap();
        let regions = read_regions_tsv(path.to_str().unwrap()).unwrap();
        std::fs::remove_file(&path).ok();
        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].group, "G1");
        assert_eq!(regions[0].len(), 11);
    }

    #[test]
    fn headerless_group_named_group_is_not_dropped() {
        let dir = std::env::temp_dir();
        let path = dir.join(format!(
            "phase_tools_regions_group_{}.tsv",
            std::process::id()
        ));
        std::fs::write(&path, "group\tchr1\t10\t20\tcopy1\n").unwrap();
        let regions = read_regions_tsv(path.to_str().unwrap()).unwrap();
        std::fs::remove_file(&path).ok();
        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].group, "group");
        assert_eq!(regions[0].chrom, "chr1");
    }

    #[test]
    fn vcf_escape_percent_encodes_reserved_characters() {
        assert_eq!(vcf_escape("G 1;%=x"), "G%201%3B%25%3Dx");
    }

    #[test]
    fn writes_candidate_tsv_from_library() {
        let candidates = vec![JointCandidate {
            group: "G1".to_string(),
            offset1: 6,
            alt_base: b'C',
            alt_positive_depth: 5,
            alt_positive_alt_count: 3,
            regions_with_alt: 1,
            region_count: 2,
            observations: vec![RegionObservation {
                copy: "copy1".to_string(),
                chrom: "chr1".to_string(),
                pos1: 15,
                ref_base: b'A',
                depth: 5,
                alt_count: 3,
                alt_fraction: 0.6,
            }],
        }];
        let mut out = Vec::new();
        write_candidates_tsv(&mut out, &candidates).unwrap();
        let text = String::from_utf8(out).unwrap();
        assert!(text.starts_with("group\toffset1\talt\talt_positive_depth"));
        assert!(text.contains("G1\t6\tC\t5\t3\t1\t2\tcopy1|chr1:15|A|5|3|0.600000"));
    }

    #[test]
    fn writes_diagnostic_vcf_from_library() {
        let regions = vec![Region {
            group: "G 1".to_string(),
            chrom: "chr1".to_string(),
            start1: 10,
            end1: 20,
            copy: "copy;1".to_string(),
        }];
        let candidates = vec![JointCandidate {
            group: "G 1".to_string(),
            offset1: 6,
            alt_base: b'C',
            alt_positive_depth: 5,
            alt_positive_alt_count: 3,
            regions_with_alt: 1,
            region_count: 1,
            observations: vec![RegionObservation {
                copy: "copy;1".to_string(),
                chrom: "chr1".to_string(),
                pos1: 15,
                ref_base: b'A',
                depth: 5,
                alt_count: 3,
                alt_fraction: 0.6,
            }],
        }];
        let mut out = Vec::new();
        write_diagnostic_vcf(&mut out, &regions, &candidates).unwrap();
        let text = String::from_utf8(out).unwrap();
        assert!(text.contains("##INFO=<ID=EVENT,"));
        assert!(text.contains("EVENT=MRJD:G%201:6:C"));
        assert!(text.contains("MRJD_COPY=copy%3B1"));
        assert!(text.contains("chr1\t15\t.\tA\tC\t.\tPASS"));
    }
}
