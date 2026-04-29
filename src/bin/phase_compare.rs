use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::{self, Read};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};

const BCF_INT32_MISSING: i32 = i32::MIN;
const BCF_INT32_VECTOR_END: i32 = i32::MIN + 1;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct VariantKey {
    chrom: String,
    pos: i64, // 1-based
    ref_allele: String,
    alt_alleles: String,
}

#[derive(Debug, Clone)]
struct Site {
    key: VariantKey,
    alleles: [i32; 2],
    phased: bool,
    het: bool,
    ps: Option<i64>,
    selected_snv: bool,
    gt: String,
}

#[derive(Debug, Default, Clone)]
struct Stats {
    truth_records: u64,
    query_records: u64,
    common_records: u64,
    truth_only_records: u64,
    query_only_records: u64,
    common_het: u64,
    genotype_mismatches: u64,
    truth_phased_het: u64,
    query_phased_het: u64,
    both_phased_het_with_ps: u64,
    intersection_blocks: u64,
    intersection_variants: u64,
    assessed_pairs: u64,
    phase_match_pairs: u64,
    switch_errors: u64,
    blockwise_hamming: u64,
    blockwise_hamming_denominator: u64,
}

impl Stats {
    fn add(&mut self, other: &Stats) {
        self.truth_records += other.truth_records;
        self.query_records += other.query_records;
        self.common_records += other.common_records;
        self.truth_only_records += other.truth_only_records;
        self.query_only_records += other.query_only_records;
        self.common_het += other.common_het;
        self.genotype_mismatches += other.genotype_mismatches;
        self.truth_phased_het += other.truth_phased_het;
        self.query_phased_het += other.query_phased_het;
        self.both_phased_het_with_ps += other.both_phased_het_with_ps;
        self.intersection_blocks += other.intersection_blocks;
        self.intersection_variants += other.intersection_variants;
        self.assessed_pairs += other.assessed_pairs;
        self.phase_match_pairs += other.phase_match_pairs;
        self.switch_errors += other.switch_errors;
        self.blockwise_hamming += other.blockwise_hamming;
        self.blockwise_hamming_denominator += other.blockwise_hamming_denominator;
    }
}

#[derive(Debug)]
struct Config {
    truth: String,
    query: String,
    sample: Option<String>,
    truth_sample: Option<String>,
    query_sample: Option<String>,
    ignore_sample_name: bool,
    only_snvs: bool,
    threads: usize,
    switch_bed: Option<String>,
    pair_tsv: Option<String>,
    summary_tsv: Option<String>,
}

#[derive(Debug, Clone)]
struct BlockState {
    chrom: String,
    truth_ps: i64,
    query_ps: i64,
    prev_truth_site: Site,
    prev_query_site: Site,
    prev_orientation: u8,
    len: u64,
    switches: u64,
    orientation_counts: [u64; 2],
}

fn usage() -> &'static str {
    "usage: phase_compare [options] truth.vcf|bcf query.vcf|bcf\n\n\
Fast phase-concordance comparison for two VCF/BCF files. The tool compares\n\
exact shared variant records, diploid heterozygous GT phase, PS block\n\
intersections, pairwise phase matches, and switch errors. It does not perform\n\
generic haplotype variant-call matching like hap.py.\n\n\
options:\n\
  -s, --sample NAME          Sample name used in both files (default: first truth sample)\n\
      --truth-sample NAME    Sample name in truth file\n\
      --query-sample NAME    Sample name in query file\n\
      --ignore-sample-name   Use the first sample in the query if the truth sample name is absent\n\
      --only-snvs            Restrict to heterozygous selected SNV genotypes\n\
  -@, --threads N            htslib reader threads (default: 1)\n\
  -o, --report-prefix PREFIX Write PREFIX.summary.tsv, like hap.py's -o prefix\n\
      --summary-tsv FILE     Write summary TSV to FILE as well as stdout\n\
      --switch-bed FILE      Write switch-error intervals as BED\n\
      --switch-error-bed FILE Alias for --switch-bed, compatible with whatshap compare\n\
      --pair-tsv FILE        Write assessed adjacent-pair decisions\n\
      --tsv-pairwise FILE    Alias for --pair-tsv, compatible with whatshap compare\n\
  -r, --reference FILE       Accepted for hap.py-style compatibility; ignored\n\
      --engine NAME          Accepted for hap.py-style compatibility; ignored\n\
      --no-roc               Accepted for hap.py-style compatibility; ignored\n\
      --no-decompose         Accepted for hap.py-style compatibility; ignored\n\
      --names NAMES          Accepted for whatshap-compare compatibility; ignored\n\
      --tsv-multiway FILE    Accepted for whatshap-compare compatibility; ignored\n\
  -h, --help                 Show this help\n\n\
Output is a TSV summary with per-contig rows and a final TOTAL row.\n"
}

fn die(msg: &str) -> ! {
    eprintln!("error: {msg}");
    std::process::exit(1);
}

fn parse_args() -> Config {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.iter().any(|a| a == "-h" || a == "--help") {
        print!("{}", usage());
        std::process::exit(0);
    }

    let mut sample = None;
    let mut truth_sample = None;
    let mut query_sample = None;
    let mut ignore_sample_name = false;
    let mut only_snvs = false;
    let mut threads = 1usize;
    let mut switch_bed = None;
    let mut pair_tsv = None;
    let mut summary_tsv = None;
    let mut positional = Vec::new();

    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "-s" | "--sample" => {
                i += 1;
                if i >= args.len() {
                    die("--sample requires an argument");
                }
                sample = Some(args[i].clone());
            }
            "--truth-sample" => {
                i += 1;
                if i >= args.len() {
                    die("--truth-sample requires an argument");
                }
                truth_sample = Some(args[i].clone());
            }
            "--query-sample" => {
                i += 1;
                if i >= args.len() {
                    die("--query-sample requires an argument");
                }
                query_sample = Some(args[i].clone());
            }
            "--ignore-sample-name" => ignore_sample_name = true,
            "--only-snvs" => only_snvs = true,
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
            "-o" | "--report-prefix" => {
                i += 1;
                if i >= args.len() {
                    die("--report-prefix requires an argument");
                }
                summary_tsv = Some(format!("{}.summary.tsv", args[i]));
            }
            "--summary-tsv" => {
                i += 1;
                if i >= args.len() {
                    die("--summary-tsv requires an argument");
                }
                summary_tsv = Some(args[i].clone());
            }
            "--switch-bed" | "--switch-error-bed" => {
                i += 1;
                if i >= args.len() {
                    die("--switch-bed requires an argument");
                }
                switch_bed = Some(args[i].clone());
            }
            "--pair-tsv" | "--tsv-pairwise" => {
                i += 1;
                if i >= args.len() {
                    die("--pair-tsv requires an argument");
                }
                pair_tsv = Some(args[i].clone());
            }
            "-r" | "--reference" | "--engine" | "--names" | "--tsv-multiway" => {
                i += 1;
                if i >= args.len() {
                    die(&format!("{} requires an argument", args[i - 1]));
                }
            }
            "--no-roc" | "--no-decompose" => {}
            x if x.starts_with('-') => die(&format!("unknown option: {x}")),
            _ => positional.push(args[i].clone()),
        }
        i += 1;
    }

    if positional.len() != 2 {
        die("expected exactly two input VCF/BCF files");
    }
    if sample.is_some() && (truth_sample.is_some() || query_sample.is_some()) {
        die("use either --sample or --truth-sample/--query-sample, not both");
    }

    Config {
        truth: positional.remove(0),
        query: positional.remove(0),
        sample,
        truth_sample,
        query_sample,
        ignore_sample_name,
        only_snvs,
        threads,
        switch_bed,
        pair_tsv,
        summary_tsv,
    }
}

fn sample_names(reader: &bcf::Reader) -> Vec<String> {
    reader
        .header()
        .samples()
        .iter()
        .map(|s| String::from_utf8_lossy(s).into_owned())
        .collect()
}

fn find_sample_idx(
    reader: &bcf::Reader,
    wanted: Option<&str>,
    fallback_first: bool,
) -> Result<usize, String> {
    let samples = sample_names(reader);
    if samples.is_empty() {
        return Err("input has no samples".to_string());
    }
    if let Some(name) = wanted {
        if let Some(idx) = samples.iter().position(|s| s == name) {
            return Ok(idx);
        }
        if fallback_first {
            return Ok(0);
        }
        return Err(format!("sample not found: {name}"));
    }
    Ok(0)
}

fn allele_index(a: GenotypeAllele) -> Option<i32> {
    match a {
        GenotypeAllele::Unphased(i) | GenotypeAllele::Phased(i) => Some(i),
        GenotypeAllele::UnphasedMissing | GenotypeAllele::PhasedMissing => None,
    }
}

fn is_phased_after_first(a: GenotypeAllele) -> bool {
    matches!(a, GenotypeAllele::Phased(_) | GenotypeAllele::PhasedMissing)
}

fn get_sample_ps(record: &bcf::Record, sample_idx: usize) -> Option<i64> {
    let values = record.format(b"PS").integer().ok()?;
    let sample_values = values.get(sample_idx)?;
    let value = *sample_values.first()?;
    if value == BCF_INT32_MISSING || value == BCF_INT32_VECTOR_END {
        None
    } else {
        Some(value as i64)
    }
}

fn gt_string(alleles: [i32; 2], phased: bool) -> String {
    let sep = if phased { '|' } else { '/' };
    format!("{}{}{}", alleles[0], sep, alleles[1])
}

fn selected_snv(allele_strings: &[String], alleles: [i32; 2]) -> bool {
    alleles.iter().all(|&a| {
        a >= 0
            && allele_strings
                .get(a as usize)
                .map(|s| {
                    s.len() == 1
                        && s.bytes().all(|b| {
                            matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N')
                        })
                })
                .unwrap_or(false)
    })
}

fn site_from_record(
    record: &bcf::Record,
    header: &rust_htslib::bcf::header::HeaderView,
    sample_idx: usize,
) -> Result<Option<Site>, String> {
    let Some(rid) = record.rid() else {
        return Ok(None);
    };
    let chrom = String::from_utf8_lossy(
        header
            .rid2name(rid)
            .map_err(|e| format!("failed to resolve RID {rid}: {e}"))?,
    )
    .into_owned();
    let pos = record.pos() + 1;
    let allele_strings = record
        .alleles()
        .iter()
        .map(|a| String::from_utf8_lossy(a).into_owned())
        .collect::<Vec<_>>();
    if allele_strings.len() < 2 {
        return Ok(None);
    }
    let ref_allele = allele_strings[0].clone();
    let alt_alleles = allele_strings[1..].join(",");

    let genotypes = match record.genotypes() {
        Ok(g) => g,
        Err(_) => return Ok(None),
    };
    let gt = genotypes.get(sample_idx);
    if gt.len() != 2 {
        return Ok(None);
    }
    let Some(a0) = allele_index(gt[0]) else {
        return Ok(None);
    };
    let Some(a1) = allele_index(gt[1]) else {
        return Ok(None);
    };
    let phased = is_phased_after_first(gt[1]);
    let het = a0 != a1;
    let alleles = [a0, a1];
    let selected_snv = selected_snv(&allele_strings, alleles);
    let gt = gt_string(alleles, phased);
    let ps = get_sample_ps(record, sample_idx);

    Ok(Some(Site {
        key: VariantKey {
            chrom,
            pos,
            ref_allele,
            alt_alleles,
        },
        alleles,
        phased,
        het,
        ps,
        selected_snv,
        gt,
    }))
}

fn read_sites(
    path: &str,
    sample: Option<&str>,
    fallback_first: bool,
    threads: usize,
    only_snvs: bool,
) -> Result<Vec<Site>, String> {
    let mut reader =
        bcf::Reader::from_path(path).map_err(|e| format!("cannot open {path}: {e}"))?;
    if threads > 1 {
        reader
            .set_threads(threads)
            .map_err(|e| format!("failed to set threads for {path}: {e}"))?;
    }
    let sample_idx = find_sample_idx(&reader, sample, fallback_first)?;
    let header = reader.header().clone();
    let mut out = Vec::new();
    for result in reader.records() {
        let record = result.map_err(|e| format!("failed to read {path}: {e}"))?;
        if let Some(site) = site_from_record(&record, &header, sample_idx)? {
            if !only_snvs || site.selected_snv {
                out.push(site);
            }
        }
    }
    Ok(out)
}

fn same_unordered_gt(a: [i32; 2], b: [i32; 2]) -> bool {
    (a[0] == b[0] && a[1] == b[1]) || (a[0] == b[1] && a[1] == b[0])
}

fn orientation(truth: [i32; 2], query: [i32; 2]) -> Option<u8> {
    if truth[0] == query[0] && truth[1] == query[1] {
        Some(0)
    } else if truth[0] == query[1] && truth[1] == query[0] {
        Some(1)
    } else {
        None
    }
}

fn chrom_stats_mut<'a>(stats: &'a mut HashMap<String, Stats>, chrom: &str) -> &'a mut Stats {
    stats.entry(chrom.to_string()).or_default()
}

fn close_block(block: Option<BlockState>, stats: &mut HashMap<String, Stats>) {
    let Some(block) = block else {
        return;
    };
    if block.len < 2 {
        return;
    }
    let st = chrom_stats_mut(stats, &block.chrom);
    st.intersection_blocks += 1;
    st.intersection_variants += block.len;
    st.assessed_pairs += block.len - 1;
    st.switch_errors += block.switches;
    st.phase_match_pairs += (block.len - 1).saturating_sub(block.switches);
    st.blockwise_hamming += block.orientation_counts[0].min(block.orientation_counts[1]);
    st.blockwise_hamming_denominator += block.len;
}

fn write_pair_header(pair_writer: &mut Option<BufWriter<File>>) -> Result<(), String> {
    if let Some(w) = pair_writer.as_mut() {
        writeln!(
            w,
            "chrom\tprev_pos\tpos\ttruth_ps\tquery_ps\tprev_orientation\torientation\tstatus\tprev_gt_truth\tgt_truth\tprev_gt_query\tgt_query"
        )
        .map_err(|e| format!("failed to write pair TSV header: {e}"))?;
    }
    Ok(())
}

fn write_pair(
    pair_writer: &mut Option<BufWriter<File>>,
    prev_truth: &Site,
    current_truth: &Site,
    prev_query: &Site,
    current_query: &Site,
    truth_ps: i64,
    query_ps: i64,
    prev_orientation: u8,
    orientation: u8,
) -> Result<(), String> {
    if let Some(w) = pair_writer.as_mut() {
        let status = if prev_orientation == orientation {
            "match"
        } else {
            "switch"
        };
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            current_truth.key.chrom,
            prev_truth.key.pos,
            current_truth.key.pos,
            truth_ps,
            query_ps,
            prev_orientation,
            orientation,
            status,
            prev_truth.gt,
            current_truth.gt,
            prev_query.gt,
            current_query.gt
        )
        .map_err(|e| format!("failed to write pair TSV: {e}"))?;
    }
    Ok(())
}

fn write_switch(
    switch_writer: &mut Option<BufWriter<File>>,
    prev: &Site,
    current: &Site,
) -> Result<(), String> {
    if let Some(w) = switch_writer.as_mut() {
        let start0 = (prev.key.pos - 1).max(0);
        let end0 = current.key.pos.max(start0 + 1);
        writeln!(w, "{}\t{}\t{}", current.key.chrom, start0, end0)
            .map_err(|e| format!("failed to write switch BED: {e}"))?;
    }
    Ok(())
}

fn open_optional_writer(path: &Option<String>) -> Result<Option<BufWriter<File>>, String> {
    path.as_ref()
        .map(|p| {
            File::create(p)
                .map(BufWriter::new)
                .map_err(|e| format!("cannot create {p}: {e}"))
        })
        .transpose()
}

fn compare(cfg: &Config) -> Result<HashMap<String, Stats>, String> {
    let truth_probe = bcf::Reader::from_path(&cfg.truth)
        .map_err(|e| format!("cannot open {}: {e}", cfg.truth))?;
    let truth_samples = sample_names(&truth_probe);
    if truth_samples.is_empty() {
        return Err("truth input has no samples".to_string());
    }
    let truth_sample = cfg
        .truth_sample
        .as_deref()
        .or(cfg.sample.as_deref())
        .unwrap_or(&truth_samples[0])
        .to_string();
    drop(truth_probe);

    let query_sample = cfg
        .query_sample
        .as_deref()
        .or(cfg.sample.as_deref())
        .unwrap_or(&truth_sample)
        .to_string();

    let truth = read_sites(
        &cfg.truth,
        Some(&truth_sample),
        false,
        cfg.threads,
        cfg.only_snvs,
    )?;
    let query = read_sites(
        &cfg.query,
        Some(&query_sample),
        cfg.ignore_sample_name,
        cfg.threads,
        cfg.only_snvs,
    )?;

    let mut stats: HashMap<String, Stats> = HashMap::new();
    let mut truth_keys = HashSet::new();
    for site in &truth {
        truth_keys.insert(site.key.clone());
        chrom_stats_mut(&mut stats, &site.key.chrom).truth_records += 1;
    }

    let mut query_map: HashMap<VariantKey, Site> = HashMap::new();
    for site in query {
        chrom_stats_mut(&mut stats, &site.key.chrom).query_records += 1;
        query_map.insert(site.key.clone(), site);
    }

    for key in query_map.keys() {
        if !truth_keys.contains(key) {
            chrom_stats_mut(&mut stats, &key.chrom).query_only_records += 1;
        }
    }

    let mut switch_writer = open_optional_writer(&cfg.switch_bed)?;
    let mut pair_writer = open_optional_writer(&cfg.pair_tsv)?;
    write_pair_header(&mut pair_writer)?;

    let mut block: Option<BlockState> = None;
    for t in &truth {
        let Some(q) = query_map.get(&t.key) else {
            chrom_stats_mut(&mut stats, &t.key.chrom).truth_only_records += 1;
            // Non-common records are outside the intersection and do not break
            // phase blocks among later common records with the same PS values.
            continue;
        };

        {
            let st = chrom_stats_mut(&mut stats, &t.key.chrom);
            st.common_records += 1;
            if t.het && q.het {
                st.common_het += 1;
            }
            if t.het && t.phased {
                st.truth_phased_het += 1;
            }
            if q.het && q.phased {
                st.query_phased_het += 1;
            }
        }

        let assessable = if t.het && q.het {
            if !same_unordered_gt(t.alleles, q.alleles) {
                chrom_stats_mut(&mut stats, &t.key.chrom).genotype_mismatches += 1;
                None
            } else if t.phased && q.phased {
                match (t.ps, q.ps, orientation(t.alleles, q.alleles)) {
                    (Some(tp), Some(qp), Some(orient)) => {
                        chrom_stats_mut(&mut stats, &t.key.chrom).both_phased_het_with_ps += 1;
                        Some((tp, qp, orient))
                    }
                    _ => None,
                }
            } else {
                None
            }
        } else {
            None
        };

        let Some((truth_ps, query_ps, orient)) = assessable else {
            // Homozygous, unphased, missing-PS, and genotype-mismatched sites
            // are not assessed as switch transitions. They also do not by
            // themselves split a later intersection block; PS changes at the
            // next assessable site do that explicitly.
            continue;
        };

        let extend = block
            .as_ref()
            .map(|b| b.chrom == t.key.chrom && b.truth_ps == truth_ps && b.query_ps == query_ps)
            .unwrap_or(false);

        if !extend {
            close_block(block.take(), &mut stats);
            let mut orientation_counts = [0u64; 2];
            orientation_counts[orient as usize] += 1;
            block = Some(BlockState {
                chrom: t.key.chrom.clone(),
                truth_ps,
                query_ps,
                prev_truth_site: t.clone(),
                prev_query_site: q.clone(),
                prev_orientation: orient,
                len: 1,
                switches: 0,
                orientation_counts,
            });
            continue;
        }

        let b = block.as_mut().expect("extend implies a block");
        write_pair(
            &mut pair_writer,
            &b.prev_truth_site,
            t,
            &b.prev_query_site,
            q,
            truth_ps,
            query_ps,
            b.prev_orientation,
            orient,
        )?;
        if b.prev_orientation != orient {
            b.switches += 1;
            write_switch(&mut switch_writer, &b.prev_truth_site, t)?;
        }
        b.len += 1;
        b.orientation_counts[orient as usize] += 1;
        b.prev_truth_site = t.clone();
        b.prev_query_site = q.clone();
        b.prev_orientation = orient;
    }
    close_block(block.take(), &mut stats);

    Ok(stats)
}

fn rate(num: u64, den: u64) -> String {
    if den == 0 {
        "NA".to_string()
    } else {
        format!("{:.6}", num as f64 / den as f64)
    }
}

fn write_summary<W: Write>(stats: &HashMap<String, Stats>, out: &mut W) -> Result<(), String> {
    writeln!(
        out,
        "chrom\ttruth_records\tquery_records\tcommon_records\ttruth_only_records\tquery_only_records\tcommon_het\tgenotype_mismatches\ttruth_phased_het\tquery_phased_het\tboth_phased_het_with_ps\tintersection_blocks\tintersection_variants\tassessed_pairs\tphase_match_pairs\tswitch_errors\tswitch_rate\tblockwise_hamming\tblockwise_hamming_rate"
    )
    .map_err(|e| format!("failed to write summary header: {e}"))?;
    let mut chroms = stats.keys().cloned().collect::<Vec<_>>();
    chroms.sort_by(|a, b| {
        natural_chrom_key(a)
            .cmp(&natural_chrom_key(b))
            .then_with(|| a.cmp(b))
    });
    let mut total = Stats::default();
    for chrom in chroms {
        let st = &stats[&chrom];
        total.add(st);
        write_row(out, &chrom, st)?;
    }
    write_row(out, "TOTAL", &total)?;
    Ok(())
}

fn natural_chrom_key(chrom: &str) -> (u8, u32, String) {
    let c = chrom.strip_prefix("chr").unwrap_or(chrom);
    if let Ok(n) = c.parse::<u32>() {
        return (0, n, String::new());
    }
    match c {
        "X" => (1, 23, String::new()),
        "Y" => (1, 24, String::new()),
        "M" | "MT" => (1, 25, String::new()),
        _ => (2, u32::MAX, c.to_string()),
    }
}

fn write_row<W: Write>(out: &mut W, chrom: &str, st: &Stats) -> Result<(), String> {
    writeln!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        chrom,
        st.truth_records,
        st.query_records,
        st.common_records,
        st.truth_only_records,
        st.query_only_records,
        st.common_het,
        st.genotype_mismatches,
        st.truth_phased_het,
        st.query_phased_het,
        st.both_phased_het_with_ps,
        st.intersection_blocks,
        st.intersection_variants,
        st.assessed_pairs,
        st.phase_match_pairs,
        st.switch_errors,
        rate(st.switch_errors, st.assessed_pairs),
        st.blockwise_hamming,
        rate(st.blockwise_hamming, st.blockwise_hamming_denominator),
    )
    .map_err(|e| format!("failed to write summary row: {e}"))
}

fn main() {
    let cfg = parse_args();
    match compare(&cfg) {
        Ok(stats) => {
            let mut stdout = std::io::BufWriter::new(std::io::stdout());
            if let Err(e) = write_summary(&stats, &mut stdout) {
                die(&e);
            }
            if let Some(path) = cfg.summary_tsv.as_deref() {
                let mut out = File::create(path)
                    .map(BufWriter::new)
                    .unwrap_or_else(|e| die(&format!("cannot create {path}: {e}")));
                if let Err(e) = write_summary(&stats, &mut out) {
                    die(&e);
                }
            }
        }
        Err(e) => die(&e),
    }
}
