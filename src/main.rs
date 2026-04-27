use libc::c_void;
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::{self, Read};
use rust_htslib::htslib;
use std::cmp::Ordering;
use std::ffi::CString;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

const PS_MISSING: i64 = i64::MIN;
const BCF_INT32_MISSING: i32 = i32::MIN;
const BCF_INT32_VECTOR_END: i32 = i32::MIN + 1;

#[derive(Debug, Clone)]
struct Config {
    input_path: String,
    fasta_path: String,
    output_path: Option<String>,
    sample_name: Option<String>,
    max_gap: i64,
    min_variants: usize,
    no_ref_check: bool,
    no_header: bool,
    quiet: bool,
}

#[derive(Debug, Clone)]
struct HeaderInfo {
    sample: String,
    contigs: Vec<String>,
    rid_names: Vec<String>,
}

#[derive(Debug, Clone)]
struct Obs {
    rid: i32,
    hap: i32,
    ps: i64,
    pos: i64,
    end: i64,
    ref_allele: String,
    alt_allele: String,
    is_snv: bool,
}

#[derive(Debug, Clone)]
struct MnvCall {
    rid: i32,
    start: i64,
    end: i64,
    ref_seq: String,
    alt_seq: String,
    positions: String,
    nvars: usize,
    nsnps: usize,
    call_type: &'static str,
    hap_mask: i32,
    ps: i64,
}

#[derive(Debug, Clone, Default)]
struct Stats {
    records: u64,
    phased_records: u64,
    observations: u64,
    skipped_no_gt: u64,
    skipped_not_diploid: u64,
    skipped_missing_gt: u64,
    skipped_unphased: u64,
    skipped_ref: u64,
    skipped_unsupported_alt: u64,
    skipped_ref_allele: u64,
    emitted: u64,
}

struct Faidx(*mut htslib::faidx_t);
impl Drop for Faidx {
    fn drop(&mut self) {
        if !self.0.is_null() {
            unsafe { htslib::fai_destroy(self.0) };
        }
    }
}

fn die(msg: &str) -> ! {
    eprintln!("error: {msg}");
    std::process::exit(1);
}

fn print_usage<W: Write>(mut out: W) -> io::Result<()> {
    write!(
        out,
        concat!(
            "usage: phase_mnv -r ref.fa [options] input.vcf|input.bcf\n",
            "\n",
            "Build minimal merged MNV/complex records from phased variants in one sample.\n",
            "\n",
            "required:\n",
            "  -r, --reference FILE   Indexed or indexable FASTA reference\n",
            "\n",
            "options:\n",
            "  -s, --sample NAME      Sample to read (default: first sample)\n",
            "  -o, --output FILE      Output VCF path (default: stdout; plain text)\n",
            "  -g, --max-gap N        Allow up to N unchanged reference bases between\n",
            "                        phased variants when building one merged call (default: 0)\n",
            "      --min-vars N       Minimum source variants per emitted call (default: 2)\n",
            "      --min-snvs N       Alias for --min-vars\n",
            "      --no-ref-check     Do not fail when VCF REF differs from FASTA\n",
            "      --no-header        Suppress VCF header\n",
            "  -q, --quiet            Suppress summary on stderr\n",
            "  -h, --help             Show this help\n",
            "\n",
            "Notes:\n",
            "  * Only phased diploid GT (e.g. 0|1, 1|0, 1|1) is used. Unphased\n",
            "    genotypes and symbolic/breakend/non-DNA alleles are skipped.\n",
            "  * FORMAT/PS is honored when present; variants are only merged within the\n",
            "    same phase set. If PS is absent, the phase separator and proximity\n",
            "    define the merge block.\n",
            "  * With the default --max-gap 0, only adjacent phased variants are\n",
            "    merged. Pure SNV blocks are TYPE=MNV; blocks containing indels are\n",
            "    TYPE=COMPLEX.\n"
        )
    )
}

fn parse_i64(s: &str, name: &str) -> i64 {
    match s.parse::<i64>() {
        Ok(v) => v,
        Err(_) => {
            eprintln!("error: invalid {name}: {s}");
            std::process::exit(1);
        }
    }
}

fn parse_args() -> Config {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut fasta_path: Option<String> = None;
    let mut output_path: Option<String> = None;
    let mut sample_name: Option<String> = None;
    let mut max_gap = 0i64;
    let mut min_variants = 2usize;
    let mut no_ref_check = false;
    let mut no_header = false;
    let mut quiet = false;
    let mut positional: Vec<String> = Vec::new();

    let mut i = 0usize;
    while i < args.len() {
        let arg = &args[i];
        match arg.as_str() {
            "-h" | "--help" => {
                let _ = print_usage(io::stdout());
                std::process::exit(0);
            }
            "-r" | "--reference" => {
                i += 1;
                if i >= args.len() {
                    die("--reference requires an argument");
                }
                fasta_path = Some(args[i].clone());
            }
            "-s" | "--sample" => {
                i += 1;
                if i >= args.len() {
                    die("--sample requires an argument");
                }
                sample_name = Some(args[i].clone());
            }
            "-o" | "--output" => {
                i += 1;
                if i >= args.len() {
                    die("--output requires an argument");
                }
                output_path = Some(args[i].clone());
            }
            "-g" | "--max-gap" => {
                i += 1;
                if i >= args.len() {
                    die("--max-gap requires an argument");
                }
                max_gap = parse_i64(&args[i], "--max-gap");
            }
            "--min-vars" | "--min-snvs" => {
                i += 1;
                if i >= args.len() {
                    die("--min-vars requires an argument");
                }
                let value = parse_i64(&args[i], "--min-vars");
                if value < 0 {
                    die("--min-vars must be >= 2");
                }
                min_variants = value as usize;
            }
            "--no-ref-check" => no_ref_check = true,
            "--no-header" => no_header = true,
            "-q" | "--quiet" => quiet = true,
            _ if arg.starts_with('-') => {
                let _ = print_usage(io::stderr());
                die(&format!("unknown option: {arg}"));
            }
            _ => positional.push(arg.clone()),
        }
        i += 1;
    }

    if fasta_path.is_none() {
        let _ = print_usage(io::stderr());
        die("--reference is required");
    }
    if max_gap < 0 {
        die("--max-gap must be >= 0");
    }
    if min_variants < 2 {
        die("--min-vars must be >= 2");
    }
    if positional.len() != 1 {
        let _ = print_usage(io::stderr());
        die("exactly one input VCF/BCF is required");
    }

    Config {
        input_path: positional.remove(0),
        fasta_path: fasta_path.unwrap(),
        output_path,
        sample_name,
        max_gap,
        min_variants,
        no_ref_check,
        no_header,
        quiet,
    }
}

#[inline]
fn upbase(b: u8) -> u8 {
    b.to_ascii_uppercase()
}

fn is_dna_base(b: u8) -> bool {
    matches!(upbase(b), b'A' | b'C' | b'G' | b'T' | b'N')
}

fn is_symbolic_or_breakend(s: &[u8]) -> bool {
    s.is_empty() || s[0] == b'<' || s.contains(&b'[') || s.contains(&b']')
}

fn is_plain_dna_allele(s: &[u8]) -> bool {
    !is_symbolic_or_breakend(s) && s.iter().all(|&b| is_dna_base(b))
}

fn uppercase_ascii_string(s: &[u8]) -> String {
    s.iter().map(|&b| upbase(b) as char).collect()
}

fn is_snv_pair(ref_allele: &str, alt_allele: &str) -> bool {
    ref_allele.len() == 1 && alt_allele.len() == 1
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

fn get_sample_ps(record: &bcf::Record, sample_idx: usize) -> i64 {
    let values = match record.format(b"PS").integer() {
        Ok(v) => v,
        Err(_) => return PS_MISSING,
    };
    let Some(sample_values) = values.get(sample_idx) else {
        return PS_MISSING;
    };
    let Some(&value) = sample_values.first() else {
        return PS_MISSING;
    };
    if value == BCF_INT32_MISSING || value == BCF_INT32_VECTOR_END {
        PS_MISSING
    } else {
        value as i64
    }
}

fn collect_header_info(reader: &bcf::Reader, sample_idx: usize) -> Result<HeaderInfo, String> {
    let header = reader.header();
    let samples = header.samples();
    let sample = String::from_utf8_lossy(samples[sample_idx]).into_owned();

    // Keep the same compact contig output as the C implementation.  For normal
    // VCF headers this is the header/RID order used by htslib.
    let mut contigs = Vec::new();
    for rid in 0..header.contig_count() {
        let name = header
            .rid2name(rid)
            .map_err(|e| format!("failed to resolve contig {rid}: {e}"))?;
        contigs.push(String::from_utf8_lossy(name).into_owned());
    }
    let rid_names = contigs.clone();

    Ok(HeaderInfo {
        sample,
        contigs,
        rid_names,
    })
}

fn read_observations(cfg: &Config) -> Result<(HeaderInfo, Vec<Obs>, Stats), String> {
    let mut reader = bcf::Reader::from_path(&cfg.input_path)
        .map_err(|e| format!("cannot open input '{}': {}", cfg.input_path, e))?;
    let sample_count = reader.header().sample_count();
    if sample_count == 0 {
        return Err("input has no samples; phased GT is required".to_string());
    }
    let sample_idx = match cfg.sample_name.as_deref() {
        None => 0usize,
        Some(name) => reader.header().sample_id(name.as_bytes()).ok_or_else(|| {
            let available = reader
                .header()
                .samples()
                .iter()
                .map(|s| String::from_utf8_lossy(s).into_owned())
                .collect::<Vec<_>>()
                .join(" ");
            format!("sample '{name}' not found. Available samples: {available}")
        })?,
    };
    let header_info = collect_header_info(&reader, sample_idx)?;

    let mut obs = Vec::new();
    let mut st = Stats::default();
    let mut record = reader.empty_record();

    while let Some(result) = reader.read(&mut record) {
        result.map_err(|e| format!("failed to read VCF/BCF record: {e}"))?;
        st.records += 1;

        let genotypes = match record.genotypes() {
            Ok(g) => g,
            Err(_) => {
                st.skipped_no_gt += 1;
                continue;
            }
        };
        let gt = genotypes.get(sample_idx);
        if gt.len() != 2 {
            st.skipped_not_diploid += 1;
            continue;
        }
        if allele_index(gt[0]).is_none() || allele_index(gt[1]).is_none() {
            st.skipped_missing_gt += 1;
            continue;
        }
        if !is_phased_after_first(gt[1]) {
            st.skipped_unphased += 1;
            continue;
        }
        st.phased_records += 1;

        let alleles = record.alleles();
        if alleles.is_empty() {
            st.skipped_ref += 1;
            continue;
        }
        let ref_bytes = alleles[0];
        if !is_plain_dna_allele(ref_bytes) {
            st.skipped_ref += 1;
            continue;
        }

        let ps = get_sample_ps(&record, sample_idx);
        let pos1 = record.pos() + 1;
        let ref_allele = uppercase_ascii_string(ref_bytes);
        let ref_len = ref_allele.len() as i64;

        for hap in 0..2usize {
            let allele = allele_index(gt[hap]).unwrap();
            if allele == 0 {
                st.skipped_ref_allele += 1;
                continue;
            }
            if allele < 0 || allele as usize >= alleles.len() {
                st.skipped_unsupported_alt += 1;
                continue;
            }
            let alt_bytes = alleles[allele as usize];
            if !is_plain_dna_allele(alt_bytes) {
                st.skipped_unsupported_alt += 1;
                continue;
            }
            let alt_allele = uppercase_ascii_string(alt_bytes);
            if ref_allele.eq_ignore_ascii_case(&alt_allele) {
                st.skipped_unsupported_alt += 1;
                continue;
            }
            let is_snv = is_snv_pair(&ref_allele, &alt_allele);
            obs.push(Obs {
                rid: record.rid().unwrap_or(0) as i32,
                hap: hap as i32,
                ps,
                pos: pos1,
                end: pos1 + ref_len - 1,
                ref_allele: ref_allele.clone(),
                alt_allele,
                is_snv,
            });
            st.observations += 1;
        }
    }

    Ok((header_info, obs, st))
}

fn cmp_obs(a: &Obs, b: &Obs) -> Ordering {
    a.rid
        .cmp(&b.rid)
        .then_with(|| a.hap.cmp(&b.hap))
        .then_with(|| a.ps.cmp(&b.ps))
        .then_with(|| a.pos.cmp(&b.pos))
        .then_with(|| a.end.cmp(&b.end))
        .then_with(|| a.alt_allele.cmp(&b.alt_allele))
        .then_with(|| a.ref_allele.cmp(&b.ref_allele))
}

fn cmp_calls(a: &MnvCall, b: &MnvCall) -> Ordering {
    a.rid
        .cmp(&b.rid)
        .then_with(|| a.start.cmp(&b.start))
        .then_with(|| a.end.cmp(&b.end))
        .then_with(|| a.ref_seq.cmp(&b.ref_seq))
        .then_with(|| a.alt_seq.cmp(&b.alt_seq))
        .then_with(|| a.positions.cmp(&b.positions))
}

fn can_extend(a: &Obs, b: &Obs, max_gap: i64) -> bool {
    if a.rid != b.rid || a.hap != b.hap || a.ps != b.ps {
        return false;
    }
    if b.pos <= a.end {
        return false;
    }
    (b.pos - a.end - 1) <= max_gap
}

fn make_positions_string(obs: &[Obs]) -> String {
    let mut s = String::new();
    for (i, o) in obs.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        s.push_str(&o.pos.to_string());
    }
    s
}

fn fetch_seq(
    fai: *mut htslib::faidx_t,
    chrom: &str,
    start1: i64,
    end1: i64,
) -> Result<String, String> {
    let c_chrom = CString::new(chrom.as_bytes()).map_err(|_| "contig contains NUL".to_string())?;
    let mut len: htslib::hts_pos_t = 0;
    let ptr = unsafe {
        htslib::faidx_fetch_seq64(
            fai,
            c_chrom.as_ptr(),
            (start1 - 1) as htslib::hts_pos_t,
            (end1 - 1) as htslib::hts_pos_t,
            &mut len,
        )
    };
    let expected = end1 - start1 + 1;
    if ptr.is_null() || len != expected as htslib::hts_pos_t {
        unsafe {
            if !ptr.is_null() {
                libc::free(ptr as *mut c_void);
            }
        }
        return Err(format!(
            "failed to fetch reference {chrom}:{start1}-{end1} (got length {len})"
        ));
    }
    let bytes = unsafe { std::slice::from_raw_parts(ptr as *const u8, len as usize) };
    let s: String = bytes.iter().map(|&b| upbase(b) as char).collect();
    unsafe {
        libc::free(ptr as *mut c_void);
    }
    Ok(s)
}

fn fetch_left_base(fai: *mut htslib::faidx_t, chrom: &str, pos1: i64) -> Result<u8, String> {
    let seq = fetch_seq(fai, chrom, pos1, pos1)?;
    Ok(seq.as_bytes()[0])
}

// Normalize a biallelic VCF representation in-place using the left-aligned +
// parsimonious rules from:
//
//   Tan A, Abecasis GR, Kang HM. Unified representation of genetic variants.
//   Bioinformatics. 2015;31(13):2202-2204. doi:10.1093/bioinformatics/btv112
//
// Algorithm: repeatedly right-trim common suffixes with left-extension when an
// allele becomes empty, then left-trim common prefixes while all alleles remain
// non-empty. This makes our output left-aligned and parsimonious without an
// external `vt normalize`/`bcftools norm` pass.
fn normalize_biallelic(
    fai: *mut htslib::faidx_t,
    chrom: &str,
    pos1: &mut i64,
    ref_seq: &mut String,
    alt_seq: &mut String,
) -> Result<(), String> {
    let mut changed = true;
    while changed {
        changed = false;
        let rlen = ref_seq.len();
        let alen = alt_seq.len();

        if rlen == 0 || alen == 0 {
            if *pos1 <= 1 {
                return Err(format!(
                    "cannot left-extend variant at beginning of contig {chrom}"
                ));
            }
            let new_pos = *pos1 - 1;
            let base = fetch_left_base(fai, chrom, new_pos)? as char;
            ref_seq.insert(0, base);
            alt_seq.insert(0, base);
            *pos1 = new_pos;
            changed = true;
            continue;
        }

        let would_empty_at_contig_start = *pos1 == 1 && (rlen == 1 || alen == 1);
        let rlast = ref_seq.as_bytes()[rlen - 1].to_ascii_uppercase();
        let alast = alt_seq.as_bytes()[alen - 1].to_ascii_uppercase();
        if !would_empty_at_contig_start && rlast == alast {
            ref_seq.pop();
            alt_seq.pop();
            changed = true;
        }
    }

    while ref_seq.len() > 1
        && alt_seq.len() > 1
        && ref_seq.as_bytes()[0].eq_ignore_ascii_case(&alt_seq.as_bytes()[0])
    {
        ref_seq.remove(0);
        alt_seq.remove(0);
        *pos1 += 1;
    }

    Ok(())
}

fn add_call_from_block(
    cfg: &Config,
    header: &HeaderInfo,
    fai: *mut htslib::faidx_t,
    block: &[Obs],
    calls: &mut Vec<MnvCall>,
) -> Result<(), String> {
    let a = &block[0];
    let chrom = header
        .rid_names
        .get(a.rid as usize)
        .ok_or_else(|| format!("input record has invalid contig id {}", a.rid))?;
    let span_end = block.iter().map(|o| o.end).max().unwrap();
    let mut ref_seq = fetch_seq(fai, chrom, a.pos, span_end)
        .map_err(|e| format!("{e} from '{}'", cfg.fasta_path))?;
    let expected_len = span_end - a.pos + 1;
    let mut alt_seq = String::with_capacity(ref_seq.len() + 64);
    let mut cursor = a.pos;
    let mut nsnps = 0usize;
    let nvars = block.len();
    let mut all_snvs = true;

    for obs in block {
        if obs.pos < cursor {
            return Err(format!("overlapping phased records at {chrom}:{} on one haplotype; normalize/decompose input first", obs.pos));
        }
        let copy_start_off = (cursor - a.pos) as usize;
        let copy_len = (obs.pos - cursor) as usize;
        if copy_len > 0 {
            alt_seq.push_str(&ref_seq[copy_start_off..copy_start_off + copy_len]);
        }

        let off = (obs.pos - a.pos) as usize;
        let ref_len = obs.ref_allele.len();
        if off + ref_len > expected_len as usize {
            return Err("internal offset bug while building merged call".to_string());
        }
        if !cfg.no_ref_check {
            let fasta_piece = &ref_seq[off..off + ref_len];
            if !fasta_piece.eq_ignore_ascii_case(&obs.ref_allele) {
                return Err(format!(
                    "REF/FASTA mismatch at {chrom}:{} (VCF REF={} FASTA={}). Use --no-ref-check to ignore.",
                    obs.pos, obs.ref_allele, fasta_piece
                ));
            }
        }
        alt_seq.push_str(&obs.alt_allele);
        cursor = obs.end + 1;
        if obs.is_snv {
            nsnps += 1;
        } else {
            all_snvs = false;
        }
    }

    let tail_off = (cursor - a.pos) as usize;
    let tail_len = (span_end - cursor + 1).max(0) as usize;
    if tail_len > 0 {
        alt_seq.push_str(&ref_seq[tail_off..tail_off + tail_len]);
    }

    if ref_seq == alt_seq {
        return Ok(());
    }

    let mut norm_pos = a.pos;
    normalize_biallelic(fai, chrom, &mut norm_pos, &mut ref_seq, &mut alt_seq)?;

    calls.push(MnvCall {
        rid: a.rid,
        start: norm_pos,
        end: norm_pos + ref_seq.len() as i64 - 1,
        ref_seq,
        alt_seq,
        positions: make_positions_string(block),
        nvars,
        nsnps,
        call_type: if all_snvs { "MNV" } else { "COMPLEX" },
        hap_mask: 1 << a.hap,
        ps: a.ps,
    });
    Ok(())
}

fn build_calls(
    cfg: &Config,
    header: &HeaderInfo,
    fai: *mut htslib::faidx_t,
    obs: &mut [Obs],
) -> Result<Vec<MnvCall>, String> {
    let mut calls = Vec::new();
    if obs.is_empty() {
        return Ok(calls);
    }
    obs.sort_by(cmp_obs);
    let mut i = 0usize;
    while i < obs.len() {
        let mut j = i;
        while j + 1 < obs.len() && can_extend(&obs[j], &obs[j + 1], cfg.max_gap) {
            j += 1;
        }
        if j - i + 1 >= cfg.min_variants {
            add_call_from_block(cfg, header, fai, &obs[i..=j], &mut calls)?;
        }
        i = j + 1;
    }
    Ok(calls)
}

fn merge_duplicate_calls(calls: &mut Vec<MnvCall>) {
    if calls.is_empty() {
        return;
    }
    calls.sort_by(cmp_calls);
    let mut out: Vec<MnvCall> = Vec::with_capacity(calls.len());
    for call in calls.drain(..) {
        if let Some(last) = out.last_mut() {
            if last.rid == call.rid
                && last.start == call.start
                && last.end == call.end
                && last.ref_seq == call.ref_seq
                && last.alt_seq == call.alt_seq
                && last.positions == call.positions
            {
                last.hap_mask |= call.hap_mask;
                if last.ps != call.ps {
                    last.ps = PS_MISSING;
                }
                continue;
            }
        }
        out.push(call);
    }
    *calls = out;
}

fn gt_for_mask(mask: i32) -> &'static str {
    match mask & 3 {
        1 => "1|0",
        2 => "0|1",
        3 => "1|1",
        _ => "./.",
    }
}

fn haps_for_mask(mask: i32) -> &'static str {
    match mask & 3 {
        3 => "1,2",
        1 => "1",
        2 => "2",
        _ => ".",
    }
}

fn write_header<W: Write>(out: &mut W, cfg: &Config, header: &HeaderInfo) -> io::Result<()> {
    writeln!(out, "##fileformat=VCFv4.3")?;
    writeln!(out, "##source=phase_mnv")?;
    writeln!(
        out,
        "##phase_mnv_normalization=Tan2015_left_aligned_parsimonious"
    )?;
    writeln!(out, "##phase_mnv_normalization_citation=Tan_A_Abecasis_GR_Kang_HM_Bioinformatics_2015_31_13_2202_2204_doi_10.1093/bioinformatics/btv112")?;
    writeln!(out, "##phase_mnv_input={}", cfg.input_path)?;
    writeln!(out, "##reference={}", cfg.fasta_path)?;
    writeln!(out, "##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Merged call type: MNV for pure SNV blocks, COMPLEX when indels are included\">")?;
    writeln!(out, "##INFO=<ID=NVAR,Number=1,Type=Integer,Description=\"Number of phased source variants merged into this call\">")?;
    writeln!(out, "##INFO=<ID=NSNPS,Number=1,Type=Integer,Description=\"Number of source SNVs in this merged call\">")?;
    writeln!(out, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate of merged reference span\">")?;
    writeln!(out, "##INFO=<ID=SOURCE_POS,Number=.,Type=Integer,Description=\"Original source variant positions merged into this call\">")?;
    writeln!(out, "##INFO=<ID=HAPS,Number=.,Type=Integer,Description=\"One-based phased haplotypes carrying this merged call\">")?;
    writeln!(out, "##INFO=<ID=PS,Number=1,Type=Integer,Description=\"Phase set shared by merged variants, when present in input FORMAT/PS\">")?;
    writeln!(out, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotype for the constructed call in the selected sample\">")?;
    writeln!(out, "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set for the constructed call, or missing if absent/ambiguous\">")?;
    for contig in &header.contigs {
        writeln!(out, "##contig=<ID={contig}>")?;
    }
    writeln!(
        out,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
        header.sample
    )?;
    Ok(())
}

fn write_calls<W: Write>(
    out: &mut W,
    header: &HeaderInfo,
    calls: &[MnvCall],
    st: &mut Stats,
) -> io::Result<()> {
    for c in calls {
        let chrom = header
            .rid_names
            .get(c.rid as usize)
            .map(String::as_str)
            .unwrap_or(".");
        let gt = gt_for_mask(c.hap_mask);
        let haps = haps_for_mask(c.hap_mask);
        write!(
            out,
            "{}\t{}\t.\t{}\t{}\t.\tPASS\tTYPE={};NVAR={};NSNPS={};END={};SOURCE_POS={};HAPS={}",
            chrom,
            c.start,
            c.ref_seq,
            c.alt_seq,
            c.call_type,
            c.nvars,
            c.nsnps,
            c.end,
            c.positions,
            haps
        )?;
        if c.ps != PS_MISSING {
            write!(out, ";PS={}", c.ps)?;
        }
        let ps = if c.ps == PS_MISSING {
            ".".to_string()
        } else {
            c.ps.to_string()
        };
        writeln!(out, "\tGT:PS\t{}:{}", gt, ps)?;
        st.emitted += 1;
    }
    Ok(())
}

fn print_summary(cfg: &Config, st: &Stats, sample: &str) {
    if cfg.quiet {
        return;
    }
    eprintln!(
        "phase_mnv: sample={} records={} phased_records={} haplotype_variant_observations={} emitted_calls={}",
        sample, st.records, st.phased_records, st.observations, st.emitted
    );
    eprintln!(
        "phase_mnv: skipped no_gt={} non_diploid={} missing_gt={} unphased={} unsupported_ref={} unsupported_alt={} ref_hap_alleles={}",
        st.skipped_no_gt,
        st.skipped_not_diploid,
        st.skipped_missing_gt,
        st.skipped_unphased,
        st.skipped_ref,
        st.skipped_unsupported_alt,
        st.skipped_ref_allele
    );
}

fn run() -> Result<(), String> {
    let cfg = parse_args();
    let (header, mut obs, mut st) = read_observations(&cfg)?;

    let fasta = CString::new(cfg.fasta_path.as_bytes())
        .map_err(|_| "FASTA path contains NUL".to_string())?;
    let fai_ptr = unsafe { htslib::fai_load(fasta.as_ptr()) };
    if fai_ptr.is_null() {
        return Err(format!(
            "cannot load or create FASTA index for '{}'",
            cfg.fasta_path
        ));
    }
    let fai = Faidx(fai_ptr);

    let mut calls = build_calls(&cfg, &header, fai.0, &mut obs)?;
    merge_duplicate_calls(&mut calls);

    match cfg.output_path.as_deref() {
        None | Some("-") => {
            let stdout = io::stdout();
            let mut out = BufWriter::new(stdout.lock());
            if !cfg.no_header {
                write_header(&mut out, &cfg, &header).map_err(|e| e.to_string())?;
            }
            write_calls(&mut out, &header, &calls, &mut st).map_err(|e| e.to_string())?;
            out.flush().map_err(|e| e.to_string())?;
        }
        Some(path) => {
            let file = File::create(Path::new(path))
                .map_err(|e| format!("cannot open output '{}': {}", path, e))?;
            let mut out = BufWriter::new(file);
            if !cfg.no_header {
                write_header(&mut out, &cfg, &header).map_err(|e| e.to_string())?;
            }
            write_calls(&mut out, &header, &calls, &mut st).map_err(|e| e.to_string())?;
            out.flush().map_err(|e| e.to_string())?;
        }
    }

    print_summary(&cfg, &st, &header.sample);
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        die(&e);
    }
}
