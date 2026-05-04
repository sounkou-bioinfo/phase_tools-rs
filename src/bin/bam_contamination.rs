use libc::c_void;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::{self, Read as BamRead};
use rust_htslib::bcf::record::{GenotypeAllele, Numeric};
use rust_htslib::bcf::{self, Read as BcfRead};
use rust_htslib::htslib;
use std::ffi::CString;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

fn usage() -> &'static str {
    "usage: bam_contamination --reference ref.fa --bam reads.bam --anchors anchors.tsv|vcf|vcf.gz|vcf.bgz|bcf [options]\n\n\
Experimental anchor-site contamination probe for BAM/CRAM data. It counts read\n\
bases at caller-supplied anchor sites and reports raw reference infiltration at\n\
homozygous-alternate anchors, with an optional CHARR-like adjustment when a\n\
reference allele frequency is supplied. It applies no MAPQ/baseQ filter by\n\
default; optional thresholds are explicit.\n\n\
Anchor TSV requires a header with columns: chrom, pos, ref, alt, gt, optional\n\
ref_af. Positions are 1-based. Anchor VCF/VCF.GZ/VCF.BGZ/BCF input is also supported; GT is\n\
read from the selected sample, INFO/REF_AF is used when present, and INFO/AF is\n\
interpreted as ALT frequency when REF_AF is absent. GT currently supports\n\
biallelic 0/0, 0/1, 1/0, and 1/1 forms with / or |. REF alleles are validated\n\
against the supplied FASTA.\n\n\
options:\n\
  -r, --reference FILE       Required FASTA reference (REF validation; CRAM decoding)\n\
      --bam FILE             Indexed BAM/CRAM read evidence\n\
      --anchors FILE         Anchor TSV or VCF/VCF.GZ/VCF.BGZ/BCF SNV anchors\n\
      --sample NAME          Sample name for VCF/BCF anchors (default: first sample)\n\
  -o, --output FILE          Output TSV (default: stdout)\n\
  -@, --threads N            htslib reader threads (default: 1)\n\
      --min-mapq N           Optional MAPQ cutoff (default: 0; no cutoff)\n\
      --min-baseq N          Optional baseQ cutoff (default: 0; no cutoff)\n\
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
    sample: Option<String>,
    output: Option<String>,
    threads: usize,
    min_mapq: u8,
    min_baseq: u8,
    include_duplicates: bool,
    include_secondary: bool,
    include_supplementary: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GenotypeClass {
    HomRef,
    Het,
    HomAlt,
}

impl GenotypeClass {
    fn label(self) -> &'static str {
        match self {
            GenotypeClass::HomRef => "hom_ref",
            GenotypeClass::Het => "het",
            GenotypeClass::HomAlt => "hom_alt",
        }
    }
}

#[derive(Debug, Clone)]
struct Anchor {
    chrom: String,
    pos: i64,
    source: String,
    ref_base: u8,
    alt_base: u8,
    gt: String,
    class: GenotypeClass,
    ref_af: Option<f64>,
}

#[derive(Debug, Default)]
struct Counts {
    observations: u64,
    ref_count: u64,
    alt_count: u64,
    other_count: u64,
    ignored_count: u64,
}

#[derive(Debug, Default)]
struct Summary {
    anchors: u64,
    hom_alt_sites: u64,
    hom_alt_observations: u64,
    hom_alt_ref: u64,
    hom_alt_alt: u64,
    charr_sum: f64,
    charr_sites: u64,
}

struct Fai(*mut htslib::faidx_t);

impl Drop for Fai {
    fn drop(&mut self) {
        unsafe { htslib::fai_destroy(self.0) };
    }
}

impl Summary {
    fn add(&mut self, anchor: &Anchor, counts: &Counts) {
        self.anchors += 1;
        if anchor.class == GenotypeClass::HomAlt {
            self.hom_alt_sites += 1;
            let denom = counts.ref_count + counts.alt_count;
            self.hom_alt_observations += denom;
            self.hom_alt_ref += counts.ref_count;
            self.hom_alt_alt += counts.alt_count;
            if denom > 0 {
                let ref_balance = counts.ref_count as f64 / denom as f64;
                if let Some(ref_af) = anchor.ref_af {
                    if ref_af > 0.0 {
                        self.charr_sum += ref_balance / ref_af;
                        self.charr_sites += 1;
                    }
                }
            }
        }
    }

    fn mean_ref_balance(&self) -> String {
        let denom = self.hom_alt_ref + self.hom_alt_alt;
        if denom == 0 {
            "NA".to_string()
        } else {
            format!("{:.6}", self.hom_alt_ref as f64 / denom as f64)
        }
    }

    fn mean_charr_like(&self) -> String {
        if self.charr_sites == 0 {
            "NA".to_string()
        } else {
            format!("{:.6}", self.charr_sum / self.charr_sites as f64)
        }
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

fn parse_args() -> Config {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.iter().any(|a| a == "-h" || a == "--help") {
        print!("{}", usage());
        std::process::exit(0);
    }

    let mut reference = None;
    let mut bam = None;
    let mut anchors = None;
    let mut sample = None;
    let mut output = None;
    let mut threads = 1usize;
    let mut min_mapq = 0u8;
    let mut min_baseq = 0u8;
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
            "--sample" => {
                i += 1;
                if i >= args.len() {
                    die("--sample requires an argument");
                }
                sample = Some(args[i].clone());
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
        sample,
        output,
        threads,
        min_mapq,
        min_baseq,
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

fn parse_gt(gt: &str, line_no: usize) -> Result<GenotypeClass, String> {
    let sep = if gt.contains('|') { '|' } else { '/' };
    let parts = gt.split(sep).collect::<Vec<_>>();
    if parts.len() != 2 {
        return Err(format!(
            "gt on anchor TSV line {line_no} must be diploid biallelic"
        ));
    }
    let allele = |s: &str| match s {
        "0" => Ok(0u8),
        "1" => Ok(1u8),
        _ => Err(format!(
            "gt on anchor TSV line {line_no} currently supports only 0 and 1 alleles"
        )),
    };
    let a0 = allele(parts[0])?;
    let a1 = allele(parts[1])?;
    Ok(match (a0, a1) {
        (0, 0) => GenotypeClass::HomRef,
        (1, 1) => GenotypeClass::HomAlt,
        _ => GenotypeClass::Het,
    })
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

fn optional_column(header: &[&str], names: &[&str]) -> Option<usize> {
    header
        .iter()
        .position(|col| names.iter().any(|name| col == name))
}

fn read_tsv_anchors(path: &str) -> Result<Vec<Anchor>, String> {
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
    let gt_col = required_column(&header, &["gt", "genotype"])?;
    let ref_af_col = optional_column(&header, &["ref_af", "refAF", "ref_frequency"]);

    let mut anchors = Vec::new();
    for (idx, line) in lines.enumerate() {
        let line_no = idx + 2;
        let line = line.map_err(|e| format!("failed to read anchor TSV: {e}"))?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        let mut needed = [chrom_col, pos_col, ref_col, alt_col, gt_col]
            .into_iter()
            .max()
            .expect("non-empty columns");
        if let Some(col) = ref_af_col {
            needed = needed.max(col);
        }
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
        let gt = fields[gt_col].to_string();
        let class = parse_gt(&gt, line_no)?;
        let ref_af = if let Some(col) = ref_af_col {
            if fields[col].is_empty() || fields[col] == "." {
                None
            } else {
                let value = fields[col]
                    .parse::<f64>()
                    .map_err(|_| format!("invalid ref_af on anchor TSV line {line_no}"))?;
                if !(0.0..=1.0).contains(&value) {
                    return Err(format!(
                        "ref_af on anchor TSV line {line_no} must be between 0 and 1"
                    ));
                }
                Some(value)
            }
        } else {
            None
        };
        anchors.push(Anchor {
            chrom,
            pos,
            source: format!("TSV line {line_no}"),
            ref_base,
            alt_base,
            gt,
            class,
            ref_af,
        });
    }
    Ok(anchors)
}

fn genotype_index(allele: GenotypeAllele, context: &str) -> Result<u32, String> {
    allele
        .index()
        .ok_or_else(|| format!("missing GT allele at {context}"))
}

fn classify_vcf_gt(
    gt: &rust_htslib::bcf::record::Genotype,
    context: &str,
) -> Result<(GenotypeClass, String), String> {
    if gt.len() != 2 {
        return Err(format!("GT at {context} must be diploid"));
    }
    let a0 = genotype_index(gt[0], context)?;
    let a1 = genotype_index(gt[1], context)?;
    if a0 > 1 || a1 > 1 {
        return Err(format!(
            "GT at {context} currently supports only biallelic 0/1 anchors"
        ));
    }
    let sep = match gt[1] {
        GenotypeAllele::Phased(_) | GenotypeAllele::PhasedMissing => '|',
        GenotypeAllele::Unphased(_) | GenotypeAllele::UnphasedMissing => '/',
    };
    let class = match (a0, a1) {
        (0, 0) => GenotypeClass::HomRef,
        (1, 1) => GenotypeClass::HomAlt,
        _ => GenotypeClass::Het,
    };
    Ok((class, format!("{a0}{sep}{a1}")))
}

fn first_vcf_float(record: &bcf::Record, tag: &[u8], context: &str) -> Result<Option<f64>, String> {
    let value = match record.info(tag).float() {
        Ok(Some(values)) => values
            .first()
            .copied()
            .filter(|v| !v.is_missing() && v.is_finite())
            .map(|v| v as f64),
        Ok(None) => None,
        Err(_) => None,
    };
    if let Some(value) = value {
        if !value.is_finite() || !(0.0..=1.0).contains(&value) {
            let tag = String::from_utf8_lossy(tag);
            return Err(format!("INFO/{tag} at {context} must be between 0 and 1"));
        }
    }
    Ok(value)
}

fn vcf_ref_af(record: &bcf::Record, context: &str) -> Result<Option<f64>, String> {
    if let Some(ref_af) = first_vcf_float(record, b"REF_AF", context)? {
        return Ok(Some(ref_af));
    }
    if let Some(alt_af) = first_vcf_float(record, b"AF", context)? {
        return Ok(Some(1.0 - alt_af));
    }
    Ok(None)
}

fn read_vcf_anchors(cfg: &Config) -> Result<Vec<Anchor>, String> {
    let mut reader = bcf::Reader::from_path(&cfg.anchors)
        .map_err(|e| format!("cannot open anchor VCF/BCF '{}': {e}", cfg.anchors))?;
    if cfg.threads > 1 {
        reader.set_threads(cfg.threads).map_err(|e| {
            format!(
                "failed to enable VCF/BCF threads for '{}': {e}",
                cfg.anchors
            )
        })?;
    }
    let sample_idx = if let Some(sample) = &cfg.sample {
        reader
            .header()
            .sample_id(sample.as_bytes())
            .ok_or_else(|| format!("sample '{sample}' not found in anchor VCF/BCF"))?
    } else {
        if reader.header().sample_count() == 0 {
            return Err("anchor VCF/BCF must contain at least one sample with GT".to_string());
        }
        0
    };

    let mut anchors = Vec::new();
    for result in reader.records() {
        let record = result.map_err(|e| format!("failed to read anchor VCF/BCF: {e}"))?;
        let Some(rid) = record.rid() else {
            continue;
        };
        let chrom = String::from_utf8(
            record
                .header()
                .rid2name(rid)
                .map_err(|e| format!("failed to resolve VCF contig id {rid}: {e}"))?
                .to_vec(),
        )
        .map_err(|_| format!("VCF contig id {rid} is not valid UTF-8"))?;
        let pos = record.pos() + 1;
        let context = format!("VCF record {chrom}:{pos}");
        if pos < 1 {
            return Err(format!("{context} has invalid position"));
        }
        let alleles = record.alleles();
        if alleles.len() != 2 || alleles[0].len() != 1 || alleles[1].len() != 1 {
            continue;
        }
        let ref_base = alleles[0][0].to_ascii_uppercase();
        let alt_base = alleles[1][0].to_ascii_uppercase();
        if !is_acgt(ref_base) || !is_acgt(alt_base) || ref_base == alt_base {
            continue;
        }
        let genotypes = record
            .genotypes()
            .map_err(|e| format!("failed to read FORMAT/GT at {context}: {e}"))?;
        let gt = genotypes.get(sample_idx);
        let (class, gt) = classify_vcf_gt(&gt, &context)?;
        let ref_af = vcf_ref_af(&record, &context)?;
        anchors.push(Anchor {
            chrom,
            pos,
            source: context,
            ref_base,
            alt_base,
            gt,
            class,
            ref_af,
        });
    }
    if anchors.is_empty() {
        Err("anchor VCF/BCF contained no usable biallelic SNV anchors".to_string())
    } else {
        Ok(anchors)
    }
}

fn anchors_look_like_vcf(path: &str) -> bool {
    let lower = path.to_ascii_lowercase();
    lower.ends_with(".vcf")
        || lower.ends_with(".vcf.gz")
        || lower.ends_with(".vcf.bgz")
        || lower.ends_with(".bcf")
}

fn read_anchors(cfg: &Config) -> Result<Vec<Anchor>, String> {
    let anchors = if anchors_look_like_vcf(&cfg.anchors) {
        read_vcf_anchors(cfg)?
    } else {
        if cfg.sample.is_some() {
            return Err("--sample is only valid with VCF/BCF anchors".to_string());
        }
        read_tsv_anchors(&cfg.anchors)?
    };
    if anchors.is_empty() {
        Err("anchor input contained no anchors".to_string())
    } else {
        Ok(anchors)
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

fn fraction(num: u64, den: u64) -> String {
    if den == 0 {
        "NA".to_string()
    } else {
        format!("{:.6}", num as f64 / den as f64)
    }
}

fn charr_component(anchor: &Anchor, counts: &Counts) -> String {
    if anchor.class != GenotypeClass::HomAlt {
        return "NA".to_string();
    }
    let denom = counts.ref_count + counts.alt_count;
    let Some(ref_af) = anchor.ref_af else {
        return "NA".to_string();
    };
    if denom == 0 || ref_af <= 0.0 {
        "NA".to_string()
    } else {
        let ref_balance = counts.ref_count as f64 / denom as f64;
        format!("{:.6}", ref_balance / ref_af)
    }
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
    let anchors = read_anchors(&cfg)?;
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

    let mut out = open_output(&cfg.output)?;
    writeln!(out, "#anchors\t{}", anchors.len())
        .map_err(|e| format!("failed to write output: {e}"))?;
    writeln!(
        out,
        "chrom\tpos\tref\talt\tgt\tclass\tref_af\tobservations\tref_count\talt_count\tother_count\tignored_count\tref_fraction\talt_fraction\tother_fraction\tcharr_like_component"
    )
    .map_err(|e| format!("failed to write output header: {e}"))?;

    let mut summary = Summary::default();
    for anchor in &anchors {
        let counts = count_anchor(&mut bam, &cfg, anchor)?;
        summary.add(anchor, &counts);
        let total_with_other = counts.ref_count + counts.alt_count + counts.other_count;
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            anchor.chrom,
            anchor.pos,
            anchor.ref_base as char,
            anchor.alt_base as char,
            anchor.gt,
            anchor.class.label(),
            anchor
                .ref_af
                .map(|v| format!("{v:.6}"))
                .unwrap_or_else(|| "NA".to_string()),
            counts.observations,
            counts.ref_count,
            counts.alt_count,
            counts.other_count,
            counts.ignored_count,
            fraction(counts.ref_count, counts.ref_count + counts.alt_count),
            fraction(counts.alt_count, counts.ref_count + counts.alt_count),
            fraction(counts.other_count, total_with_other),
            charr_component(anchor, &counts)
        )
        .map_err(|e| format!("failed to write output: {e}"))?;
    }
    writeln!(out, "#hom_alt_sites\t{}", summary.hom_alt_sites)
        .map_err(|e| format!("failed to write output: {e}"))?;
    writeln!(
        out,
        "#hom_alt_observations\t{}",
        summary.hom_alt_observations
    )
    .map_err(|e| format!("failed to write output: {e}"))?;
    writeln!(
        out,
        "#hom_alt_mean_ref_balance\t{}",
        summary.mean_ref_balance()
    )
    .map_err(|e| format!("failed to write output: {e}"))?;
    writeln!(out, "#charr_like_sites\t{}", summary.charr_sites)
        .map_err(|e| format!("failed to write output: {e}"))?;
    writeln!(out, "#charr_like_mean\t{}", summary.mean_charr_like())
        .map_err(|e| format!("failed to write output: {e}"))?;
    out.flush()
        .map_err(|e| format!("failed to flush output: {e}"))?;
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        die(&e);
    }
}
