use libc::c_void;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::{self, Read as BamRead};
use rust_htslib::htslib;
use std::ffi::CString;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

fn usage() -> &'static str {
    "usage: bam_ancestry --reference ref.fa --bam reads.bam --anchors ancestry.tsv [options]\n\n\
Experimental BAM/CRAM ancestry mixture probe. It counts REF/ALT bases at\n\
caller-supplied ancestry-informative SNV anchors, estimates observed ALT\n\
fractions, and fits a constrained least-squares mixture over reference\n\
population ALT frequencies. It applies no MAPQ/baseQ filter by default;\n\
optional thresholds are explicit. This is Summix-style in spirit, but it is not\n\
a full Summix replacement.\n\n\
Anchor TSV requires a header with columns: chrom, pos, ref, alt, then reference\n\
population ALT-frequency columns. Positions are 1-based. Use --populations to\n\
select/order population columns; otherwise all columns after alt are used. REF\n\
alleles are validated against the supplied FASTA.\n\n\
options:\n\
  -r, --reference FILE       Required FASTA reference (REF validation; CRAM decoding)\n\
      --bam FILE             Indexed BAM/CRAM read evidence\n\
      --anchors FILE         Anchor TSV with chrom,pos,ref,alt,popAF...\n\
      --populations LIST     Comma-separated population columns to use\n\
  -o, --output FILE          Output TSV (default: stdout)\n\
  -@, --threads N            htslib reader threads (default: 1)\n\
      --min-mapq N           Optional MAPQ cutoff (default: 0; no cutoff)\n\
      --min-baseq N          Optional baseQ cutoff (default: 0; no cutoff)\n\
      --min-observations N   Minimum REF+ALT observations for fitting (default: 1)\n\
      --include-duplicates   Include duplicate reads\n\
      --include-secondary    Include secondary alignments\n\
      --include-supplementary Include supplementary alignments\n\
  -h, --help                 Show this help\n"
}

#[derive(Debug)]
struct Config {
    reference: String,
    bam: String,
    anchors: String,
    populations: Option<Vec<String>>,
    output: Option<String>,
    threads: usize,
    min_mapq: u8,
    min_baseq: u8,
    min_observations: u64,
    include_duplicates: bool,
    include_secondary: bool,
    include_supplementary: bool,
}

#[derive(Debug, Clone)]
struct Anchor {
    chrom: String,
    pos: i64,
    ref_base: u8,
    alt_base: u8,
    pop_alt_af: Vec<f64>,
    source: String,
}

#[derive(Debug, Default, Clone)]
struct Counts {
    observations: u64,
    ref_count: u64,
    alt_count: u64,
    other_count: u64,
    ignored_count: u64,
}

#[derive(Debug)]
struct AnchorResult {
    anchor: Anchor,
    counts: Counts,
    observed_alt_fraction: Option<f64>,
    predicted_alt_fraction: Option<f64>,
    residual: Option<f64>,
    used_for_fit: bool,
}

struct Fai(*mut htslib::faidx_t);

impl Drop for Fai {
    fn drop(&mut self) {
        unsafe { htslib::fai_destroy(self.0) };
    }
}

fn die(msg: &str) -> ! {
    eprintln!("error: {msg}");
    std::process::exit(1);
}

fn parse_i64(s: &str, name: &str) -> i64 {
    s.parse::<i64>()
        .unwrap_or_else(|_| die(&format!("{name} must be an integer")))
}

fn parse_u8(s: &str, name: &str) -> u8 {
    let value = parse_i64(s, name);
    if !(0..=255).contains(&value) {
        die(&format!("{name} must be between 0 and 255"));
    }
    value as u8
}

fn parse_u64(s: &str, name: &str) -> u64 {
    let value = parse_i64(s, name);
    if value < 0 {
        die(&format!("{name} must be >= 0"));
    }
    value as u64
}

fn parse_population_list(s: &str) -> Vec<String> {
    let pops = s
        .split(',')
        .map(str::trim)
        .filter(|p| !p.is_empty())
        .map(ToOwned::to_owned)
        .collect::<Vec<_>>();
    if pops.is_empty() {
        die("--populations must contain at least one population name");
    }
    pops
}

fn parse_args() -> Config {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.iter().any(|a| a == "-h" || a == "--help") {
        print!("{}", usage());
        std::process::exit(0);
    }

    let mut reference = None;
    let mut bam = None;
    let mut anchors = None;
    let mut populations = None;
    let mut output = None;
    let mut threads = 1usize;
    let mut min_mapq = 0u8;
    let mut min_baseq = 0u8;
    let mut min_observations = 1u64;
    let mut include_duplicates = false;
    let mut include_secondary = false;
    let mut include_supplementary = false;
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
            "--bam" => {
                i += 1;
                if i >= args.len() {
                    die("--bam requires an argument");
                }
                bam = Some(args[i].clone());
            }
            "--anchors" => {
                i += 1;
                if i >= args.len() {
                    die("--anchors requires an argument");
                }
                anchors = Some(args[i].clone());
            }
            "--populations" => {
                i += 1;
                if i >= args.len() {
                    die("--populations requires an argument");
                }
                populations = Some(parse_population_list(&args[i]));
            }
            "-o" | "--output" => {
                i += 1;
                if i >= args.len() {
                    die("--output requires an argument");
                }
                output = Some(args[i].clone());
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
            "--min-mapq" => {
                i += 1;
                if i >= args.len() {
                    die("--min-mapq requires an argument");
                }
                min_mapq = parse_u8(&args[i], "--min-mapq");
            }
            "--min-baseq" => {
                i += 1;
                if i >= args.len() {
                    die("--min-baseq requires an argument");
                }
                min_baseq = parse_u8(&args[i], "--min-baseq");
            }
            "--min-observations" => {
                i += 1;
                if i >= args.len() {
                    die("--min-observations requires an argument");
                }
                min_observations = parse_u64(&args[i], "--min-observations");
            }
            "--include-duplicates" => include_duplicates = true,
            "--include-secondary" => include_secondary = true,
            "--include-supplementary" => include_supplementary = true,
            x if x.starts_with('-') => die(&format!("unknown option: {x}")),
            _ => positional.push(args[i].clone()),
        }
        i += 1;
    }
    if !positional.is_empty() {
        die("unexpected positional arguments; use --bam and --anchors");
    }

    Config {
        reference: reference.unwrap_or_else(|| die("--reference is required")),
        bam: bam.unwrap_or_else(|| die("--bam is required")),
        anchors: anchors.unwrap_or_else(|| die("--anchors is required")),
        populations,
        output,
        threads,
        min_mapq,
        min_baseq,
        min_observations,
        include_duplicates,
        include_secondary,
        include_supplementary,
    }
}

fn load_fai(path: &str) -> Result<Fai, String> {
    let c_path = CString::new(path.as_bytes()).map_err(|_| "reference path contains NUL")?;
    let fai = unsafe { htslib::fai_load(c_path.as_ptr()) };
    if fai.is_null() {
        Err(format!("cannot load FASTA index for '{path}'"))
    } else {
        Ok(Fai(fai))
    }
}

fn fetch_ref_base(fai: &Fai, chrom: &str, pos1: i64) -> Result<u8, String> {
    if pos1 < 1 {
        return Err(format!("invalid FASTA position {chrom}:{pos1}"));
    }
    let c_chrom = CString::new(chrom.as_bytes()).map_err(|_| "contig contains NUL")?;
    let mut len: htslib::hts_pos_t = 0;
    let ptr = unsafe {
        htslib::faidx_fetch_seq64(
            fai.0,
            c_chrom.as_ptr(),
            (pos1 - 1) as htslib::hts_pos_t,
            (pos1 - 1) as htslib::hts_pos_t,
            &mut len,
        )
    };
    if ptr.is_null() || len != 1 {
        unsafe {
            if !ptr.is_null() {
                libc::free(ptr as *mut c_void);
            }
        }
        return Err(format!("failed to fetch reference base {chrom}:{pos1}"));
    }
    let base = unsafe { *(ptr as *const u8) }.to_ascii_uppercase();
    unsafe { libc::free(ptr as *mut c_void) };
    Ok(base)
}

fn is_acgt(base: u8) -> bool {
    matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

fn parse_base(s: &str, name: &str, line_no: usize) -> Result<u8, String> {
    let bytes = s.as_bytes();
    if bytes.len() != 1 {
        return Err(format!(
            "{name} on anchor TSV line {line_no} must be one base"
        ));
    }
    let base = bytes[0].to_ascii_uppercase();
    if !is_acgt(base) {
        return Err(format!(
            "{name} on anchor TSV line {line_no} must be A/C/G/T"
        ));
    }
    Ok(base)
}

fn required_column(header: &[&str], names: &[&str]) -> Result<usize, String> {
    header
        .iter()
        .position(|col| names.iter().any(|name| col == name))
        .ok_or_else(|| {
            format!(
                "anchor TSV is missing required column '{}'; accepted aliases: {}",
                names[0],
                names.join(",")
            )
        })
}

fn read_anchors(
    path: &str,
    requested_pops: &Option<Vec<String>>,
) -> Result<(Vec<String>, Vec<Anchor>), String> {
    let file = File::open(path).map_err(|e| format!("cannot open anchor TSV '{path}': {e}"))?;
    let mut lines = BufReader::new(file).lines();
    let header_line = lines
        .next()
        .ok_or_else(|| format!("anchor TSV '{path}' is empty"))?
        .map_err(|e| format!("failed to read anchor TSV header: {e}"))?;
    let header = header_line.split('\t').collect::<Vec<_>>();
    let chrom_col = required_column(&header, &["chrom", "#chrom"])?;
    let pos_col = required_column(&header, &["pos", "position"])?;
    let ref_col = required_column(&header, &["ref", "ref_base"])?;
    let alt_col = required_column(&header, &["alt", "alt_base"])?;

    let populations = if let Some(pops) = requested_pops {
        pops.clone()
    } else {
        let first_pop_col = [chrom_col, pos_col, ref_col, alt_col]
            .into_iter()
            .max()
            .expect("non-empty columns")
            + 1;
        if header.len() <= first_pop_col {
            return Err(
                "anchor TSV must contain at least one population AF column after alt".to_string(),
            );
        }
        header[first_pop_col..]
            .iter()
            .map(|s| (*s).to_string())
            .collect::<Vec<_>>()
    };
    if populations.is_empty() {
        return Err("at least one population AF column is required".to_string());
    }
    let pop_cols = populations
        .iter()
        .map(|pop| {
            header
                .iter()
                .position(|col| col == pop)
                .ok_or_else(|| format!("population column '{pop}' not found in anchor TSV"))
        })
        .collect::<Result<Vec<_>, _>>()?;

    let mut anchors = Vec::new();
    for (idx, line) in lines.enumerate() {
        let line_no = idx + 2;
        let line = line.map_err(|e| format!("failed to read anchor TSV: {e}"))?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        let needed = [chrom_col, pos_col, ref_col, alt_col]
            .into_iter()
            .chain(pop_cols.iter().copied())
            .max()
            .expect("non-empty columns");
        if fields.len() <= needed {
            return Err(format!("anchor TSV line {line_no} has too few columns"));
        }
        let chrom = fields[chrom_col].to_string();
        let pos = fields[pos_col]
            .parse::<i64>()
            .map_err(|_| format!("invalid pos on anchor TSV line {line_no}"))?;
        if pos < 1 {
            return Err(format!("pos on anchor TSV line {line_no} must be >= 1"));
        }
        let ref_base = parse_base(fields[ref_col], "ref", line_no)?;
        let alt_base = parse_base(fields[alt_col], "alt", line_no)?;
        if ref_base == alt_base {
            return Err(format!(
                "ref and alt are identical on anchor TSV line {line_no}"
            ));
        }
        let mut pop_alt_af = Vec::with_capacity(pop_cols.len());
        for (pop, col) in populations.iter().zip(pop_cols.iter().copied()) {
            let value = fields[col]
                .parse::<f64>()
                .map_err(|_| format!("invalid {pop} ALT AF on anchor TSV line {line_no}"))?;
            if !value.is_finite() || !(0.0..=1.0).contains(&value) {
                return Err(format!(
                    "{pop} ALT AF on anchor TSV line {line_no} must be between 0 and 1"
                ));
            }
            pop_alt_af.push(value);
        }
        anchors.push(Anchor {
            chrom,
            pos,
            ref_base,
            alt_base,
            pop_alt_af,
            source: format!("TSV line {line_no}"),
        });
    }
    if anchors.is_empty() {
        Err("anchor TSV contained no anchors".to_string())
    } else {
        Ok((populations, anchors))
    }
}

fn validate_anchor_reference(anchors: &[Anchor], fai: &Fai) -> Result<(), String> {
    for anchor in anchors {
        let fasta_base = fetch_ref_base(fai, &anchor.chrom, anchor.pos)?;
        if !is_acgt(fasta_base) {
            return Err(format!(
                "reference FASTA base at {}:{} is not A/C/G/T",
                anchor.chrom, anchor.pos
            ));
        }
        if fasta_base != anchor.ref_base {
            return Err(format!(
                "anchor REF mismatch at {}:{} ({}): anchor REF={}, FASTA={}",
                anchor.chrom,
                anchor.pos,
                anchor.source,
                anchor.ref_base as char,
                fasta_base as char
            ));
        }
    }
    Ok(())
}

fn usable_record(record: &bam::Record, cfg: &Config) -> bool {
    !record.is_unmapped()
        && !record.is_quality_check_failed()
        && (cfg.min_mapq == 0 || (record.mapq() != 255 && record.mapq() >= cfg.min_mapq))
        && (cfg.include_duplicates || !record.is_duplicate())
        && (cfg.include_secondary || !record.is_secondary())
        && (cfg.include_supplementary || !record.is_supplementary())
}

fn base_at(record: &bam::Record, pos0: i64) -> Option<(u8, u8)> {
    let seq = record.seq();
    let qual = record.qual();
    let mut ref_pos = record.pos();
    let mut read_pos = 0usize;
    for cigar in record.cigar().iter() {
        match *cigar {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                for _ in 0..len {
                    if ref_pos == pos0 && read_pos < seq.len() {
                        return Some((
                            seq[read_pos].to_ascii_uppercase(),
                            qual.get(read_pos).copied().unwrap_or(255),
                        ));
                    }
                    ref_pos += 1;
                    read_pos += 1;
                }
            }
            Cigar::Ins(len) | Cigar::SoftClip(len) => read_pos += len as usize,
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                if pos0 >= ref_pos && pos0 < ref_pos + len as i64 {
                    return None;
                }
                ref_pos += len as i64;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    None
}

fn count_anchor(
    bam: &mut bam::IndexedReader,
    cfg: &Config,
    anchor: &Anchor,
) -> Result<Counts, String> {
    let start0 = anchor.pos - 1;
    bam.fetch((anchor.chrom.as_bytes(), start0, anchor.pos))
        .map_err(|e| {
            format!(
                "failed to fetch {}:{}-{}: {e}",
                anchor.chrom, anchor.pos, anchor.pos
            )
        })?;
    let mut counts = Counts::default();
    for result in bam.records() {
        let record = result.map_err(|e| format!("failed to read BAM/CRAM record: {e}"))?;
        if !usable_record(&record, cfg) {
            continue;
        }
        let Some((base, q)) = base_at(&record, start0) else {
            counts.ignored_count += 1;
            continue;
        };
        if q < cfg.min_baseq {
            counts.ignored_count += 1;
            continue;
        }
        let base = base.to_ascii_uppercase();
        if base == anchor.ref_base {
            counts.ref_count += 1;
        } else if base == anchor.alt_base {
            counts.alt_count += 1;
        } else if is_acgt(base) {
            counts.other_count += 1;
        } else {
            counts.ignored_count += 1;
            continue;
        }
        counts.observations += 1;
    }
    Ok(counts)
}

fn project_simplex(values: &[f64]) -> Vec<f64> {
    let n = values.len();
    if n == 1 {
        return vec![1.0];
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    let mut cssv = 0.0;
    let mut rho = 0usize;
    for (i, value) in sorted.iter().enumerate() {
        cssv += value;
        let theta = (cssv - 1.0) / (i + 1) as f64;
        if value - theta > 0.0 {
            rho = i + 1;
        }
    }
    let theta = (sorted.iter().take(rho).sum::<f64>() - 1.0) / rho as f64;
    values.iter().map(|v| (v - theta).max(0.0)).collect()
}

fn fit_mixture(matrix: &[Vec<f64>], y: &[f64], weights: &[f64], k: usize) -> Vec<f64> {
    if k == 1 {
        return vec![1.0];
    }
    let mut p = vec![1.0 / k as f64; k];
    let trace = matrix
        .iter()
        .zip(weights.iter().copied())
        .map(|(row, w)| w * row.iter().map(|v| v * v).sum::<f64>())
        .sum::<f64>();
    let step = if trace > 0.0 {
        1.0 / (2.0 * trace)
    } else {
        1.0
    };

    for _ in 0..10_000 {
        let mut grad = vec![0.0; k];
        for ((row, observed), weight) in matrix
            .iter()
            .zip(y.iter().copied())
            .zip(weights.iter().copied())
        {
            let predicted = row.iter().zip(p.iter()).map(|(a, b)| a * b).sum::<f64>();
            let residual = predicted - observed;
            for j in 0..k {
                grad[j] += 2.0 * weight * residual * row[j];
            }
        }
        let candidate = p
            .iter()
            .zip(grad.iter())
            .map(|(pj, gj)| pj - step * gj)
            .collect::<Vec<_>>();
        let next = project_simplex(&candidate);
        let max_delta = p
            .iter()
            .zip(next.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        p = next;
        if max_delta < 1e-12 {
            break;
        }
    }
    p
}

fn fraction(num: u64, den: u64) -> Option<f64> {
    if den == 0 {
        None
    } else {
        Some(num as f64 / den as f64)
    }
}

fn fmt_opt(value: Option<f64>) -> String {
    value
        .map(|v| format!("{v:.6}"))
        .unwrap_or_else(|| "NA".to_string())
}

fn open_output(path: &Option<String>) -> Result<Box<dyn Write>, String> {
    if let Some(path) = path {
        let file = File::create(path).map_err(|e| format!("cannot create output '{path}': {e}"))?;
        Ok(Box::new(BufWriter::new(file)))
    } else {
        Ok(Box::new(BufWriter::new(std::io::stdout())))
    }
}

fn run() -> Result<(), String> {
    let cfg = parse_args();
    let (populations, anchors) = read_anchors(&cfg.anchors, &cfg.populations)?;
    let fai = load_fai(&cfg.reference)?;
    validate_anchor_reference(&anchors, &fai)?;

    let mut bam = bam::IndexedReader::from_path(&cfg.bam)
        .map_err(|e| format!("cannot open indexed BAM/CRAM '{}': {e}", cfg.bam))?;
    bam.set_reference(&cfg.reference)
        .map_err(|e| format!("failed to set CRAM reference '{}': {e}", cfg.reference))?;
    if cfg.threads > 1 {
        bam.set_threads(cfg.threads)
            .map_err(|e| format!("failed to enable BAM/CRAM threads for '{}': {e}", cfg.bam))?;
    }

    let mut results = Vec::with_capacity(anchors.len());
    let mut matrix = Vec::new();
    let mut y = Vec::new();
    let mut weights = Vec::new();
    for anchor in anchors {
        let counts = count_anchor(&mut bam, &cfg, &anchor)?;
        let denom = counts.ref_count + counts.alt_count;
        let observed_alt_fraction = fraction(counts.alt_count, denom);
        let used_for_fit = denom >= cfg.min_observations && observed_alt_fraction.is_some();
        if used_for_fit {
            matrix.push(anchor.pop_alt_af.clone());
            y.push(observed_alt_fraction.expect("checked above"));
            weights.push(denom as f64);
        }
        results.push(AnchorResult {
            anchor,
            counts,
            observed_alt_fraction,
            predicted_alt_fraction: None,
            residual: None,
            used_for_fit,
        });
    }
    if matrix.is_empty() {
        return Err("no anchors had enough REF+ALT observations for ancestry fitting".to_string());
    }

    let proportions = fit_mixture(&matrix, &y, &weights, populations.len());
    for result in &mut results {
        if result.used_for_fit {
            let predicted = result
                .anchor
                .pop_alt_af
                .iter()
                .zip(proportions.iter())
                .map(|(a, p)| a * p)
                .sum::<f64>();
            result.predicted_alt_fraction = Some(predicted);
            result.residual = result
                .observed_alt_fraction
                .map(|observed| observed - predicted);
        }
    }
    let weighted_sse = results
        .iter()
        .filter(|r| r.used_for_fit)
        .map(|r| {
            let denom = (r.counts.ref_count + r.counts.alt_count) as f64;
            denom * r.residual.unwrap_or(0.0).powi(2)
        })
        .sum::<f64>();

    let mut out = open_output(&cfg.output)?;
    writeln!(out, "#anchors\t{}", results.len())
        .map_err(|e| format!("failed to write output: {e}"))?;
    writeln!(
        out,
        "#used_anchors\t{}",
        results.iter().filter(|r| r.used_for_fit).count()
    )
    .map_err(|e| format!("failed to write output: {e}"))?;
    writeln!(out, "#weighted_sse\t{weighted_sse:.6}")
        .map_err(|e| format!("failed to write output: {e}"))?;
    writeln!(out, "population\tproportion").map_err(|e| format!("failed to write output: {e}"))?;
    for (population, proportion) in populations.iter().zip(proportions.iter()) {
        writeln!(out, "{population}\t{proportion:.6}")
            .map_err(|e| format!("failed to write output: {e}"))?;
    }
    writeln!(
        out,
        "chrom\tpos\tref\talt\tused_for_fit\tobservations\tref_count\talt_count\tother_count\tignored_count\tobserved_alt_fraction\tpredicted_alt_fraction\tresidual\t{}",
        populations
            .iter()
            .map(|p| format!("{p}_alt_af"))
            .collect::<Vec<_>>()
            .join("\t")
    )
    .map_err(|e| format!("failed to write output: {e}"))?;
    for result in &results {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            result.anchor.chrom,
            result.anchor.pos,
            result.anchor.ref_base as char,
            result.anchor.alt_base as char,
            if result.used_for_fit { "1" } else { "0" },
            result.counts.observations,
            result.counts.ref_count,
            result.counts.alt_count,
            result.counts.other_count,
            result.counts.ignored_count,
            fmt_opt(result.observed_alt_fraction),
            fmt_opt(result.predicted_alt_fraction),
            fmt_opt(result.residual),
            result
                .anchor
                .pop_alt_af
                .iter()
                .map(|v| format!("{v:.6}"))
                .collect::<Vec<_>>()
                .join("\t")
        )
        .map_err(|e| format!("failed to write output: {e}"))?;
    }
    out.flush()
        .map_err(|e| format!("failed to flush output: {e}"))?;
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        die(&e);
    }
}
