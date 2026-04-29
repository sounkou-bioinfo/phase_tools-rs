use libc::c_void;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::{self, Read};
use rust_htslib::htslib;
use std::collections::BTreeMap;
use std::ffi::CString;

fn usage() -> &'static str {
    "usage: bam_error_model --reference ref.fa [options] reads.bam|reads.cram\n\n\
Learn a simple empirical sequencing-error table from aligned reads by comparing\n\
BAM/CRAM bases to a FASTA reference. No MAPQ filter is applied by default; MAPQ\n\
is summarized as a covariate. Known variant sites are not masked yet, so real\n\
biological variants contribute to the mismatch rate unless callers restrict the\n\
regions accordingly.\n\n\
options:\n\
  -r, --reference FILE       FASTA reference with .fai\n\
      --region REG          Restrict to region CHR:START-END (1-based, repeatable)\n\
      --max-reads N         Stop after N usable reads\n\
      --min-mapq N          Optional MAPQ cutoff (default: 0; no cutoff)\n\
      --include-duplicates  Include duplicate reads\n\
      --include-secondary   Include secondary alignments\n\
      --include-supplementary Include supplementary alignments\n\
  -h, --help                Show this help\n"
}

#[derive(Debug)]
struct Config {
    bam: String,
    reference: String,
    regions: Vec<Region>,
    max_reads: Option<u64>,
    min_mapq: u8,
    include_duplicates: bool,
    include_secondary: bool,
    include_supplementary: bool,
}

#[derive(Debug, Clone)]
struct Region {
    chrom: String,
    start0: i64,
    end0: i64,
}

#[derive(Debug, Clone, Default)]
struct Stats {
    matches: u64,
    mismatches: u64,
    insertions: u64,
    deletions: u64,
}

impl Stats {
    fn observations(&self) -> u64 {
        self.matches + self.mismatches + self.insertions + self.deletions
    }

    fn errors(&self) -> u64 {
        self.mismatches + self.insertions + self.deletions
    }

    fn add_match(&mut self) {
        self.matches += 1;
    }

    fn add_mismatch(&mut self) {
        self.mismatches += 1;
    }

    fn add_insertion(&mut self) {
        self.insertions += 1;
    }

    fn add_deletion(&mut self) {
        self.deletions += 1;
    }
}

#[derive(Debug, Default)]
struct Model {
    reads: u64,
    overall: Stats,
    by_baseq: BTreeMap<String, Stats>,
    by_mapq: BTreeMap<String, Stats>,
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

fn parse_region(s: &str) -> Region {
    let (chrom, rest) = s
        .split_once(':')
        .unwrap_or_else(|| die("--region must be CHR:START-END"));
    let (start, end) = rest
        .replace(',', "")
        .split_once('-')
        .map(|(a, b)| (a.to_string(), b.to_string()))
        .unwrap_or_else(|| die("--region must be CHR:START-END"));
    let start1 = parse_i64(&start, "region start");
    let end1 = parse_i64(&end, "region end");
    if chrom.is_empty() || start1 < 1 || end1 < start1 {
        die("invalid --region coordinates");
    }
    Region {
        chrom: chrom.to_string(),
        start0: start1 - 1,
        end0: end1,
    }
}

fn parse_args() -> Config {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.iter().any(|a| a == "-h" || a == "--help") {
        print!("{}", usage());
        std::process::exit(0);
    }

    let mut reference = None;
    let mut regions = Vec::new();
    let mut max_reads = None;
    let mut min_mapq = 0u8;
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
            "--region" => {
                i += 1;
                if i >= args.len() {
                    die("--region requires an argument");
                }
                regions.push(parse_region(&args[i]));
            }
            "--max-reads" => {
                i += 1;
                if i >= args.len() {
                    die("--max-reads requires an argument");
                }
                let n = parse_i64(&args[i], "--max-reads");
                if n < 1 {
                    die("--max-reads must be >= 1");
                }
                max_reads = Some(n as u64);
            }
            "--min-mapq" => {
                i += 1;
                if i >= args.len() {
                    die("--min-mapq requires an argument");
                }
                let n = parse_i64(&args[i], "--min-mapq");
                if !(0..=255).contains(&n) {
                    die("--min-mapq must be between 0 and 255");
                }
                min_mapq = n as u8;
            }
            "--include-duplicates" => include_duplicates = true,
            "--include-secondary" => include_secondary = true,
            "--include-supplementary" => include_supplementary = true,
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
        regions,
        max_reads,
        min_mapq,
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

fn fetch_ref(fai: &Fai, chrom: &str, start0: i64, end0: i64) -> Result<Vec<u8>, String> {
    if end0 <= start0 {
        return Ok(Vec::new());
    }
    let c_chrom = CString::new(chrom.as_bytes()).map_err(|_| "contig contains NUL")?;
    let mut len: htslib::hts_pos_t = 0;
    let ptr = unsafe {
        htslib::faidx_fetch_seq64(
            fai.0,
            c_chrom.as_ptr(),
            start0 as htslib::hts_pos_t,
            (end0 - 1) as htslib::hts_pos_t,
            &mut len,
        )
    };
    let expected = end0 - start0;
    if ptr.is_null() || len != expected as htslib::hts_pos_t {
        unsafe {
            if !ptr.is_null() {
                libc::free(ptr as *mut c_void);
            }
        }
        return Err(format!(
            "failed to fetch reference {chrom}:{}-{} (got length {len})",
            start0 + 1,
            end0
        ));
    }
    let bytes = unsafe { std::slice::from_raw_parts(ptr as *const u8, len as usize) };
    let out = bytes
        .iter()
        .map(|b| b.to_ascii_uppercase())
        .collect::<Vec<_>>();
    unsafe { libc::free(ptr as *mut c_void) };
    Ok(out)
}

fn q_bin(q: u8) -> String {
    match q {
        255 => "unknown".to_string(),
        0..=9 => "00-09".to_string(),
        10..=19 => "10-19".to_string(),
        20..=29 => "20-29".to_string(),
        30..=39 => "30-39".to_string(),
        40..=49 => "40-49".to_string(),
        _ => "50+".to_string(),
    }
}

fn mapq_bin(q: u8) -> String {
    match q {
        255 => "unknown".to_string(),
        0 => "00".to_string(),
        1..=9 => "01-09".to_string(),
        10..=19 => "10-19".to_string(),
        20..=29 => "20-29".to_string(),
        30..=39 => "30-39".to_string(),
        40..=49 => "40-49".to_string(),
        50..=59 => "50-59".to_string(),
        _ => "60+".to_string(),
    }
}

fn is_acgt(b: u8) -> bool {
    matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

fn usable_record(record: &bam::Record, cfg: &Config) -> bool {
    !record.is_unmapped()
        && !record.is_quality_check_failed()
        && (cfg.min_mapq == 0 || (record.mapq() != 255 && record.mapq() >= cfg.min_mapq))
        && (cfg.include_duplicates || !record.is_duplicate())
        && (cfg.include_secondary || !record.is_secondary())
        && (cfg.include_supplementary || !record.is_supplementary())
}

fn add_match_or_mismatch(model: &mut Model, mapq_key: &str, baseq: u8, is_match: bool) {
    if is_match {
        model.overall.add_match();
        model
            .by_mapq
            .entry(mapq_key.to_string())
            .or_default()
            .add_match();
        model.by_baseq.entry(q_bin(baseq)).or_default().add_match();
    } else {
        model.overall.add_mismatch();
        model
            .by_mapq
            .entry(mapq_key.to_string())
            .or_default()
            .add_mismatch();
        model
            .by_baseq
            .entry(q_bin(baseq))
            .or_default()
            .add_mismatch();
    }
}

fn add_insertion(model: &mut Model, mapq_key: &str, baseq: u8) {
    model.overall.add_insertion();
    model
        .by_mapq
        .entry(mapq_key.to_string())
        .or_default()
        .add_insertion();
    model
        .by_baseq
        .entry(q_bin(baseq))
        .or_default()
        .add_insertion();
}

fn add_deletion(model: &mut Model, mapq_key: &str) {
    model.overall.add_deletion();
    model
        .by_mapq
        .entry(mapq_key.to_string())
        .or_default()
        .add_deletion();
}

fn in_clip(pos0: i64, clips: Option<&[(i64, i64)]>) -> bool {
    match clips {
        Some(clips) => clips
            .iter()
            .any(|&(start0, end0)| pos0 >= start0 && pos0 < end0),
        None => true,
    }
}

fn process_record(
    record: &bam::Record,
    chrom: &str,
    fai: &Fai,
    cfg: &Config,
    model: &mut Model,
    clips: Option<&[(i64, i64)]>,
) -> Result<bool, String> {
    if !usable_record(record, cfg) {
        return Ok(false);
    }
    let start0 = record.pos();
    let end0 = record.cigar().end_pos();
    if start0 < 0 || end0 <= start0 {
        return Ok(false);
    }
    let before_observations = model.overall.observations();
    let ref_seq = fetch_ref(fai, chrom, start0, end0)?;
    let seq = record.seq();
    let qual = record.qual();
    let mapq_key = mapq_bin(record.mapq());
    let mut ref_pos = start0;
    let mut read_pos = 0usize;
    let mut last_ref_pos: Option<i64> = None;

    for cigar in record.cigar().iter() {
        match *cigar {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                for _ in 0..len {
                    if read_pos < seq.len() {
                        let ref_idx = (ref_pos - start0) as usize;
                        if ref_idx < ref_seq.len() {
                            let rb = seq[read_pos].to_ascii_uppercase();
                            let fb = ref_seq[ref_idx];
                            if in_clip(ref_pos, clips) && is_acgt(rb) && is_acgt(fb) {
                                let q = qual.get(read_pos).copied().unwrap_or(255);
                                add_match_or_mismatch(model, &mapq_key, q, rb == fb);
                            }
                        }
                    }
                    last_ref_pos = Some(ref_pos);
                    ref_pos += 1;
                    read_pos += 1;
                }
            }
            Cigar::Ins(len) => {
                let anchor_pos = last_ref_pos.unwrap_or(ref_pos);
                for _ in 0..len {
                    if in_clip(anchor_pos, clips) && read_pos < seq.len() {
                        let rb = seq[read_pos].to_ascii_uppercase();
                        if is_acgt(rb) {
                            let q = qual.get(read_pos).copied().unwrap_or(255);
                            add_insertion(model, &mapq_key, q);
                        }
                    }
                    read_pos += 1;
                }
            }
            Cigar::Del(len) => {
                for _ in 0..len {
                    if in_clip(ref_pos, clips) {
                        add_deletion(model, &mapq_key);
                    }
                    last_ref_pos = Some(ref_pos);
                    ref_pos += 1;
                }
            }
            Cigar::RefSkip(len) => ref_pos += len as i64,
            Cigar::SoftClip(len) => read_pos += len as usize,
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    let contributed = model.overall.observations() > before_observations;
    if contributed {
        model.reads += 1;
    }
    Ok(contributed)
}

fn error_rate(stats: &Stats) -> String {
    let obs = stats.observations();
    if obs == 0 {
        "NA".to_string()
    } else {
        format!("{:.6}", stats.errors() as f64 / obs as f64)
    }
}

fn empirical_q(stats: &Stats) -> String {
    let obs = stats.observations();
    let err = stats.errors();
    if obs == 0 {
        "NA".to_string()
    } else if err == 0 {
        "inf".to_string()
    } else {
        let rate = err as f64 / obs as f64;
        format!("{:.3}", -10.0 * rate.log10())
    }
}

fn print_row(scope: &str, bin: &str, stats: &Stats) {
    println!(
        "{scope}\t{bin}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        stats.observations(),
        stats.matches,
        stats.mismatches,
        stats.insertions,
        stats.deletions,
        error_rate(stats),
        empirical_q(stats)
    );
}

fn print_model(model: &Model) {
    println!("#reads\t{}", model.reads);
    println!(
        "scope\tbin\tobservations\tmatches\tmismatches\tinsertions\tdeletions\terror_rate\tempirical_q"
    );
    print_row("overall", "all", &model.overall);
    for (bin, stats) in &model.by_baseq {
        print_row("baseq", bin, stats);
    }
    for (bin, stats) in &model.by_mapq {
        print_row("mapq", bin, stats);
    }
}

fn run_stream(cfg: &Config, fai: &Fai, model: &mut Model) -> Result<(), String> {
    let mut reader = bam::Reader::from_path(&cfg.bam)
        .map_err(|e| format!("cannot open BAM/CRAM '{}': {e}", cfg.bam))?;
    reader
        .set_reference(&cfg.reference)
        .map_err(|e| format!("failed to set CRAM reference '{}': {e}", cfg.reference))?;
    let header = bam::Header::from_template(reader.header());
    let header_view = bam::HeaderView::from_header(&header);
    for rec in reader.records() {
        let record = rec.map_err(|e| format!("failed to read BAM/CRAM record: {e}"))?;
        if let Some(max) = cfg.max_reads {
            if model.reads >= max {
                break;
            }
        }
        if record.tid() < 0 {
            continue;
        }
        let chrom = String::from_utf8_lossy(header_view.tid2name(record.tid() as u32)).into_owned();
        process_record(&record, &chrom, fai, cfg, model, None)?;
    }
    Ok(())
}

fn merge_regions(regions: &[Region]) -> Vec<Region> {
    let mut out = regions.to_vec();
    out.sort_by(|a, b| {
        a.chrom
            .cmp(&b.chrom)
            .then(a.start0.cmp(&b.start0))
            .then(a.end0.cmp(&b.end0))
    });
    let mut merged: Vec<Region> = Vec::new();
    for region in out {
        if let Some(last) = merged.last_mut() {
            if last.chrom == region.chrom && region.start0 <= last.end0 {
                last.end0 = last.end0.max(region.end0);
                continue;
            }
        }
        merged.push(region);
    }
    merged
}

fn grouped_regions(regions: &[Region]) -> BTreeMap<String, Vec<(i64, i64)>> {
    let mut grouped: BTreeMap<String, Vec<(i64, i64)>> = BTreeMap::new();
    for region in merge_regions(regions) {
        grouped
            .entry(region.chrom)
            .or_default()
            .push((region.start0, region.end0));
    }
    grouped
}

fn run_regions(cfg: &Config, fai: &Fai, model: &mut Model) -> Result<(), String> {
    let mut reader = bam::IndexedReader::from_path(&cfg.bam)
        .map_err(|e| format!("cannot open indexed BAM/CRAM '{}': {e}", cfg.bam))?;
    reader
        .set_reference(&cfg.reference)
        .map_err(|e| format!("failed to set CRAM reference '{}': {e}", cfg.reference))?;
    let header = bam::Header::from_template(reader.header());
    let header_view = bam::HeaderView::from_header(&header);

    for (chrom, clips) in grouped_regions(&cfg.regions) {
        let fetch_start = clips.iter().map(|&(start0, _)| start0).min().unwrap_or(0);
        let fetch_end = clips
            .iter()
            .map(|&(_, end0)| end0)
            .max()
            .unwrap_or(fetch_start);
        reader
            .fetch((chrom.as_bytes(), fetch_start, fetch_end))
            .map_err(|e| {
                format!(
                    "failed to fetch {chrom}:{}-{}: {e}",
                    fetch_start + 1,
                    fetch_end
                )
            })?;
        for rec in reader.records() {
            let record = rec.map_err(|e| format!("failed to read BAM/CRAM record: {e}"))?;
            if let Some(max) = cfg.max_reads {
                if model.reads >= max {
                    return Ok(());
                }
            }
            if record.tid() < 0 || !usable_record(&record, cfg) {
                continue;
            }
            let record_chrom =
                String::from_utf8_lossy(header_view.tid2name(record.tid() as u32)).into_owned();
            process_record(&record, &record_chrom, fai, cfg, model, Some(&clips))?;
        }
    }
    Ok(())
}

fn run() -> Result<(), String> {
    let cfg = parse_args();
    let fai = load_fai(&cfg.reference)?;
    let mut model = Model::default();
    if cfg.regions.is_empty() {
        run_stream(&cfg, &fai, &mut model)?;
    } else {
        run_regions(&cfg, &fai, &mut model)?;
    }
    print_model(&model);
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        die(&e);
    }
}
