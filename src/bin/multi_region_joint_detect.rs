use phase_tools::io::fasta::Fai;
use phase_tools::io::vcf::infer_output_kind;
use phase_tools::mrjd::{
    detect_snv_candidates_with_mapq255_policy, read_regions_tsv, write_candidates_tsv,
    write_diagnostic_vcf, CandidateConfig, Mapq255Policy,
};
use std::fs::{self, File};
use std::io::{self, BufWriter};
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
records are excluded. MAPQ 255 is retained when no MAPQ threshold is requested\n\
and otherwise dropped as unknown unless --keep-mapq-255 is set.\n\n\
regions.tsv columns:\n\
  group<TAB>chrom<TAB>start<TAB>end[<TAB>copy]\n\n\
options:\n\
  -r, --reference FILE       FASTA reference with .fai\n\
      --regions FILE         Region manifest TSV\n\
      --min-mapq N           Minimum read MAPQ (default: 0)\n\
      --min-baseq N          Minimum base quality (default: 13)\n\
      --min-alt-count N      Minimum per-region alt count (default: 2)\n\
      --min-alt-fraction F   Minimum per-region alt fraction (default: 0.20)\n\
      --keep-mapq-255        Keep MAPQ 255 reads even when --min-mapq > 0\n\
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
    keep_mapq_255: bool,
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

fn validate_vcf_sidecar_path(path: &str) -> Result<(), String> {
    if infer_output_kind(Some(path)).is_indexable() {
        return Err(
            "--vcf currently writes plain diagnostic VCF; use a .vcf or non-indexable suffix, not .vcf.gz, .vcf.bgz, or .bcf".to_string(),
        );
    }
    Ok(())
}

fn validate_output_paths(cfg: &Config) -> Result<(), String> {
    let Some(vcf) = cfg.vcf.as_deref() else {
        return Ok(());
    };
    validate_vcf_sidecar_path(vcf)?;
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
    let mut keep_mapq_255 = false;
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
            "--keep-mapq-255" => {
                keep_mapq_255 = true;
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
        keep_mapq_255,
    }
}

fn run() -> Result<(), String> {
    let cfg = parse_args();
    validate_output_paths(&cfg)?;
    let regions = read_regions_tsv(&cfg.regions)?;
    let fai = Fai::from_path(&cfg.reference)?;
    let mapq255_policy = if cfg.keep_mapq_255 {
        Mapq255Policy::KeepWhenFiltering
    } else {
        Mapq255Policy::DropWhenFiltering
    };
    let candidates = detect_snv_candidates_with_mapq255_policy(
        &cfg.bam,
        &cfg.reference,
        &fai,
        &regions,
        cfg.candidate,
        mapq255_policy,
        cfg.threads,
    )?;
    if let Some(path) = cfg.vcf.as_deref() {
        let file = File::create(Path::new(path))
            .map_err(|e| format!("cannot create VCF sidecar '{path}': {e}"))?;
        write_diagnostic_vcf(BufWriter::new(file), &regions, &candidates)?;
    }
    match cfg.output.as_deref() {
        None | Some("-") => {
            let stdout = io::stdout();
            write_candidates_tsv(BufWriter::new(stdout.lock()), &candidates)?;
        }
        Some(path) => {
            let file = File::create(Path::new(path))
                .map_err(|e| format!("cannot create output '{path}': {e}"))?;
            write_candidates_tsv(BufWriter::new(file), &candidates)?;
        }
    }
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        die(&e);
    }
}
