use phase_tools::io::fasta::Fai;
use phase_tools::mrjd::{
    detect_snv_candidates, read_regions_tsv, CandidateConfig, JointCandidate, Region,
};
use std::collections::BTreeSet;
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

fn usage() -> &'static str {
    "usage: multi_region_joint_detect --reference ref.fa --regions regions.tsv [options] reads.bam|reads.cram\n\n\
Initial multi-region SNV candidate evidence scanner. Regions sharing a group are\n\
compared by 1-based offset within each region, so homologous loci can be audited\n\
together before downstream joint genotyping. This is not a DRAGEN-equivalent\n\
caller and currently emits a TSV diagnostics table plus optional VCF diagnostics.\n\
Depth/count totals are over alt-positive region observations; ref-only or\n\
uncovered copies are counted in region_count but omitted from per_region. The\n\
VCF sidecar emits one record per alt-positive region observation. String INFO\n\
values are percent-escaped. Duplicate, QC-fail, secondary, and supplementary\n\
records are excluded. MAPQ 255 is excluded only when --min-mapq is greater than 0.\n\n\
regions.tsv columns:\n\
  group<TAB>chrom<TAB>start<TAB>end[<TAB>copy]\n\n\
options:\n\
  -r, --reference FILE       FASTA reference with .fai\n\
      --regions FILE         Region manifest TSV\n\
      --min-mapq N           Minimum read MAPQ (default: 0)\n\
      --min-baseq N          Minimum base quality (default: 13)\n\
      --min-alt-count N      Minimum per-region alt count (default: 2)\n\
      --min-alt-fraction F   Minimum per-region alt fraction (default: 0.20)\n\
  -@, --threads N            BAM/CRAM reader threads (default: 1)\n\
  -o, --output FILE          Output TSV path (default: stdout)\n\
      --vcf FILE             Write diagnostic VCF sidecar with EVENT/EVENTTYPE\n\
  -h, --help                 Show this help\n"
}

#[derive(Debug)]
struct Config {
    bam: String,
    reference: String,
    regions: String,
    candidate: CandidateConfig,
    threads: usize,
    output: Option<String>,
    vcf: Option<String>,
}

fn die(msg: &str) -> ! {
    eprintln!("error: {msg}");
    std::process::exit(1);
}

fn parse_u8(s: &str, name: &str) -> u8 {
    s.parse::<u8>()
        .unwrap_or_else(|_| die(&format!("{name} must be an integer between 0 and 255")))
}

fn parse_u32(s: &str, name: &str) -> u32 {
    s.parse::<u32>()
        .unwrap_or_else(|_| die(&format!("{name} must be an integer")))
}

fn parse_f64(s: &str, name: &str) -> f64 {
    s.parse::<f64>()
        .unwrap_or_else(|_| die(&format!("{name} must be a number")))
}

fn normalized_output_path(path: &str) -> Result<PathBuf, String> {
    let path = Path::new(path);
    let absolute = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir()
            .map_err(|e| format!("failed to get current directory: {e}"))?
            .join(path)
    };
    if absolute.exists() {
        return fs::canonicalize(&absolute)
            .map_err(|e| format!("failed to resolve '{}': {e}", absolute.display()));
    }
    let parent = absolute
        .parent()
        .ok_or_else(|| format!("invalid output path '{}'", absolute.display()))?;
    let file_name = absolute
        .file_name()
        .ok_or_else(|| format!("invalid output path '{}'", absolute.display()))?;
    let parent = fs::canonicalize(parent).map_err(|e| {
        format!(
            "failed to resolve output directory for '{}': {e}",
            absolute.display()
        )
    })?;
    Ok(parent.join(file_name))
}

fn validate_output_paths(cfg: &Config) -> Result<(), String> {
    let Some(vcf) = cfg.vcf.as_deref() else {
        return Ok(());
    };
    if let Some(output) = cfg.output.as_deref().filter(|path| *path != "-") {
        if normalized_output_path(output)? == normalized_output_path(vcf)? {
            return Err("--output and --vcf must be different paths".to_string());
        }
    }
    Ok(())
}

fn parse_args() -> Config {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.iter().any(|arg| arg == "-h" || arg == "--help") {
        print!("{}", usage());
        std::process::exit(0);
    }

    let mut reference = None;
    let mut regions = None;
    let mut candidate = CandidateConfig::default();
    let mut threads = 1usize;
    let mut output = None;
    let mut vcf = None;
    let mut positional = Vec::new();

    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "-r" | "--reference" => {
                i += 1;
                if i >= args.len() {
                    die("--reference requires an argument");
                }
                reference = Some(args[i].clone());
            }
            "--regions" => {
                i += 1;
                if i >= args.len() {
                    die("--regions requires an argument");
                }
                regions = Some(args[i].clone());
            }
            "--min-mapq" => {
                i += 1;
                if i >= args.len() {
                    die("--min-mapq requires an argument");
                }
                candidate.min_mapq = parse_u8(&args[i], "--min-mapq");
            }
            "--min-baseq" => {
                i += 1;
                if i >= args.len() {
                    die("--min-baseq requires an argument");
                }
                candidate.min_baseq = parse_u8(&args[i], "--min-baseq");
            }
            "--min-alt-count" => {
                i += 1;
                if i >= args.len() {
                    die("--min-alt-count requires an argument");
                }
                candidate.min_alt_count = parse_u32(&args[i], "--min-alt-count");
            }
            "--min-alt-fraction" => {
                i += 1;
                if i >= args.len() {
                    die("--min-alt-fraction requires an argument");
                }
                candidate.min_alt_fraction = parse_f64(&args[i], "--min-alt-fraction");
            }
            "-@" | "--threads" => {
                i += 1;
                if i >= args.len() {
                    die("--threads requires an argument");
                }
                threads = args[i]
                    .parse::<usize>()
                    .unwrap_or_else(|_| die("--threads must be an integer"));
                if threads == 0 {
                    die("--threads must be >= 1");
                }
            }
            "-o" | "--output" => {
                i += 1;
                if i >= args.len() {
                    die("--output requires an argument");
                }
                output = Some(args[i].clone());
            }
            "--vcf" => {
                i += 1;
                if i >= args.len() {
                    die("--vcf requires an argument");
                }
                if args[i] == "-" {
                    die("--vcf requires a file path; stdout is reserved for TSV output");
                }
                vcf = Some(args[i].clone());
            }
            x if x.starts_with('-') => die(&format!("unknown option: {x}")),
            _ => positional.push(args[i].clone()),
        }
        i += 1;
    }

    if positional.len() != 1 {
        die("expected exactly one BAM/CRAM input");
    }

    Config {
        bam: positional.remove(0),
        reference: reference.unwrap_or_else(|| die("--reference is required")),
        regions: regions.unwrap_or_else(|| die("--regions is required")),
        candidate,
        threads,
        output,
        vcf,
    }
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

fn write_candidates<W: Write>(mut out: W, candidates: &[JointCandidate]) -> Result<(), String> {
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

fn event_id(candidate: &JointCandidate) -> String {
    format!(
        "MRJD:{}:{}:{}",
        vcf_escape(&candidate.group),
        candidate.offset1,
        candidate.alt_base as char
    )
}

fn write_vcf_sidecar(
    path: &str,
    regions: &[Region],
    candidates: &[JointCandidate],
) -> Result<(), String> {
    let file = File::create(Path::new(path))
        .map_err(|e| format!("cannot create VCF sidecar '{path}': {e}"))?;
    let mut out = BufWriter::new(file);
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
        let event = event_id(candidate);
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

fn run() -> Result<(), String> {
    let cfg = parse_args();
    validate_output_paths(&cfg)?;
    let regions = read_regions_tsv(&cfg.regions)?;
    let fai = Fai::from_path(&cfg.reference)?;
    let candidates = detect_snv_candidates(
        &cfg.bam,
        &cfg.reference,
        &fai,
        &regions,
        cfg.candidate,
        cfg.threads,
    )?;
    if let Some(path) = cfg.vcf.as_deref() {
        write_vcf_sidecar(path, &regions, &candidates)?;
    }
    match cfg.output.as_deref() {
        None | Some("-") => {
            let stdout = io::stdout();
            write_candidates(BufWriter::new(stdout.lock()), &candidates)?;
        }
        Some(path) => {
            let file = File::create(Path::new(path))
                .map_err(|e| format!("cannot create output '{path}': {e}"))?;
            write_candidates(BufWriter::new(file), &candidates)?;
        }
    }
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        die(&e);
    }
}
