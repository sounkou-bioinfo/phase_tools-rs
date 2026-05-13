use libc::c_void;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::{self, Read as BamRead};
use rust_htslib::bcf::header::{HeaderRecord, HeaderView, TagType};
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::{self, Read as BcfRead};
use rust_htslib::htslib;
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, HashSet};
use std::ffi::CString;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

const PS_MISSING: i64 = i64::MIN;
const BCF_INT32_MISSING: i32 = i32::MIN;
const BCF_INT32_VECTOR_END: i32 = i32::MIN + 1;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum UnsupportedAllelesPolicy {
    Skip,
    Fail,
}

impl UnsupportedAllelesPolicy {
    fn as_str(self) -> &'static str {
        match self {
            UnsupportedAllelesPolicy::Skip => "skip",
            UnsupportedAllelesPolicy::Fail => "fail",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OutputKind {
    PlainVcf,
    BgzfVcf,
    Bcf,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OutputIndexKind {
    Csi,
    Tbi,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum EmitMode {
    Mnv,
    Combined,
    AllSites,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PhaseAlgorithm {
    Mec,
    Greedy,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PhaseTag {
    Ps,
    Hp,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MnvAlgorithm {
    Proximity,
    NirvanaCodon,
}

impl MnvAlgorithm {
    fn as_str(self) -> &'static str {
        match self {
            MnvAlgorithm::Proximity => "proximity",
            MnvAlgorithm::NirvanaCodon => "nirvana-codon",
        }
    }
}

impl PhaseAlgorithm {
    fn as_str(self) -> &'static str {
        match self {
            PhaseAlgorithm::Mec => "mec",
            PhaseAlgorithm::Greedy => "greedy",
        }
    }
}

impl PhaseTag {
    fn as_str(self) -> &'static str {
        match self {
            PhaseTag::Ps => "PS",
            PhaseTag::Hp => "HP",
        }
    }
}

impl EmitMode {
    fn as_str(self) -> &'static str {
        match self {
            EmitMode::Mnv => "mnv",
            EmitMode::Combined => "combined",
            EmitMode::AllSites => "all-sites",
        }
    }
}

impl OutputIndexKind {
    fn as_str(self) -> &'static str {
        match self {
            OutputIndexKind::Csi => "csi",
            OutputIndexKind::Tbi => "tbi",
        }
    }

    fn min_shift(self) -> i32 {
        match self {
            OutputIndexKind::Csi => 14,
            OutputIndexKind::Tbi => 0,
        }
    }
}

impl OutputKind {
    fn as_str(self) -> &'static str {
        match self {
            OutputKind::PlainVcf => "vcf",
            OutputKind::BgzfVcf => "vcf.gz",
            OutputKind::Bcf => "bcf",
        }
    }
}

#[derive(Debug, Clone)]
struct Config {
    input_path: String,
    fasta_path: String,
    output_path: Option<String>,
    sample_name: Option<String>,
    phase_bam_path: Option<String>,
    phase_min_mapq: u8,
    phase_min_baseq: u8,
    phase_algorithm: PhaseAlgorithm,
    phase_tag: PhaseTag,
    phase_max_coverage: usize,
    phase_only_snvs: bool,
    phase_use_supplementary: bool,
    phase_ignore_read_groups: bool,
    phase_supplementary_distance: i64,
    phase_realign_overhang: i64,
    phase_read_list_path: Option<String>,
    threads: usize,
    emit_mode: EmitMode,
    output_index: Option<OutputIndexKind>,
    mnv_algorithm: MnvAlgorithm,
    codon_map_path: Option<String>,
    max_gap: i64,
    min_variants: usize,
    unsupported_alleles: UnsupportedAllelesPolicy,
    warn_on_n: bool,
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct SourceHaplotype {
    source_record: u64,
    rid: i32,
    pos: i64,
    hap: i32,
    allele: i32,
}

#[derive(Debug, Clone)]
struct Obs {
    rid: i32,
    hap: i32,
    ps: i64,
    pos: i64,
    end: i64,
    source_record: u64,
    source_allele: i32,
    ref_allele: String,
    alt_allele: String,
    is_snv: bool,
    codon_keys: Vec<String>,
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
    source_haplotypes: Vec<SourceHaplotype>,
}

#[derive(Debug, Clone, Default)]
struct Stats {
    records: u64,
    phased_records: u64,
    observations: u64,
    multiallelic_records: u64,
    observations_with_n: u64,
    bam_phase_candidates: u64,
    bam_phase_informative_reads: u64,
    bam_phase_components: u64,
    bam_phase_phased_variants: u64,
    bam_phase_unphased_variants: u64,
    bam_phase_conflicts: u64,
    bam_phase_selected_reads: u64,
    skipped_no_gt: u64,
    skipped_not_diploid: u64,
    skipped_missing_gt: u64,
    skipped_unphased: u64,
    skipped_ref: u64,
    skipped_unsupported_alt: u64,
    skipped_alt_out_of_range: u64,
    skipped_alt_symbolic_or_breakend: u64,
    skipped_alt_spanning_deletion: u64,
    skipped_alt_non_dna: u64,
    skipped_alt_same_as_ref: u64,
    skipped_ref_allele: u64,
    emitted: u64,
}

#[derive(Debug, Clone)]
struct PhaseCandidate {
    source_record: u64,
    rid: i32,
    chrom: String,
    pos: i64,
    end: i64,
    alleles: Vec<String>,
    gt_alleles: [i32; 2],
    input_order_alleles: [i32; 2],
    is_snv: bool,
}

#[derive(Debug, Clone, Copy)]
struct PhaseAssignment {
    ps: i64,
    rel: u8,
}

#[derive(Debug, Clone)]
struct ReadPhaseCall {
    calls: Vec<(usize, u8, u8)>,
    order: usize,
    name: String,
    mapq: u8,
}

#[derive(Debug, Clone)]
struct LocalReadCall {
    var: usize,
    allele: u8,
    qual: u8,
}

#[derive(Debug, Clone)]
struct LocalRead {
    calls: Vec<LocalReadCall>,
    first_var: usize,
    last_var: usize,
}

#[derive(Debug, Clone)]
struct CodonInterval {
    start: i64,
    end: i64,
    key: String,
}

#[derive(Debug, Clone, Default)]
struct CodonMap {
    by_chrom: HashMap<String, Vec<CodonInterval>>,
}

#[derive(Debug, Clone, Copy)]
struct MecTraceState {
    cost: u64,
    prev_mask: Option<u64>,
    orientation: u8,
    ambiguous: bool,
}

#[derive(Debug, Clone)]
struct MecTraceColumn {
    active_reads: Vec<usize>,
    states: HashMap<u64, MecTraceState>,
}

#[derive(Debug, Clone, Default)]
struct ReadEvent {
    base: Option<u8>,
    qual: Option<u8>,
    insertion_after: Vec<(u8, u8)>,
}

#[derive(Debug, Clone)]
struct Dsu {
    parent: Vec<usize>,
    rank: Vec<u8>,
    xor_to_parent: Vec<u8>,
}

impl Dsu {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            rank: vec![0; n],
            xor_to_parent: vec![0; n],
        }
    }

    fn find(&mut self, x: usize) -> (usize, u8) {
        let parent = self.parent[x];
        if parent == x {
            return (x, 0);
        }
        let (root, px) = self.find(parent);
        self.parent[x] = root;
        self.xor_to_parent[x] ^= px;
        (self.parent[x], self.xor_to_parent[x])
    }

    fn union(&mut self, a: usize, b: usize, parity: u8) -> bool {
        let (mut ra, xa) = self.find(a);
        let (mut rb, xb) = self.find(b);
        if ra == rb {
            return (xa ^ xb) == parity;
        }
        let x = xa ^ xb ^ parity;
        if self.rank[ra] < self.rank[rb] {
            std::mem::swap(&mut ra, &mut rb);
        }
        self.parent[rb] = ra;
        self.xor_to_parent[rb] = x;
        if self.rank[ra] == self.rank[rb] {
            self.rank[ra] += 1;
        }
        true
    }
}

struct Faidx(*mut htslib::faidx_t);
impl Drop for Faidx {
    fn drop(&mut self) {
        if !self.0.is_null() {
            unsafe { htslib::fai_destroy(self.0) };
        }
    }
}

struct BgzfWriter(*mut htslib::BGZF);

impl BgzfWriter {
    fn from_path(path: &str, threads: usize) -> Result<Self, String> {
        let c_path = CString::new(path.as_bytes())
            .map_err(|_| format!("output path contains NUL byte: {path}"))?;
        let mode = CString::new("w").expect("literal contains no NUL");
        let fp = unsafe { htslib::bgzf_open(c_path.as_ptr(), mode.as_ptr()) };
        if fp.is_null() {
            return Err(format!("cannot open BGZF output '{path}'"));
        }
        if threads > 1 && unsafe { htslib::bgzf_mt(fp, threads as i32, 256) } != 0 {
            unsafe { htslib::bgzf_close(fp) };
            return Err(format!("failed to enable BGZF output threads for '{path}'"));
        }
        Ok(Self(fp))
    }
}

impl Write for BgzfWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let n =
            unsafe { htslib::bgzf_write(self.0, buf.as_ptr() as *const c_void, buf.len() as _) };
        if n < 0 {
            Err(io::Error::new(io::ErrorKind::Other, "BGZF write failed"))
        } else {
            Ok(n as usize)
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        let ret = unsafe { htslib::bgzf_flush(self.0) };
        if ret != 0 {
            Err(io::Error::new(io::ErrorKind::Other, "BGZF flush failed"))
        } else {
            Ok(())
        }
    }
}

impl Drop for BgzfWriter {
    fn drop(&mut self) {
        if !self.0.is_null() {
            unsafe { htslib::bgzf_close(self.0) };
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
            "  -o, --output FILE      Output path (default: stdout). Format is inferred:\n",
            "                        .vcf = plain VCF, .vcf.gz/.vcf.bgz = BGZF VCF,\n",
            "                        .bcf = BCF; stdout defaults to plain VCF\n",
            "  -@, --threads N        Extra htslib/BGZF threads for decompression and\n",
            "                        compressed output (default: 1)\n",
            "      --write-index[=FMT]\n",
            "                        Build an index after writing -o output. FMT is csi\n",
            "                        (default) or tbi; requires .vcf.gz/.vcf.bgz/.bcf\n",
            "      --emit MODE        Output mode: mnv (default), combined, or all-sites.\n",
            "                        combined emits merged MNV/COMPLEX records plus\n",
            "                        selected-sample input variants not represented by a\n",
            "                        merge; all-sites preserves input records/header and\n",
            "                        updates GT/PS when used with --phase-from-bam\n",
            "  -g, --max-gap N        Allow up to N unchanged reference bases between\n",
            "                        phased variants when building one merged call (default: 0)\n",
            "      --mnv-algorithm MODE\n",
            "                        MNV construction: proximity (default) or\n",
            "                        nirvana-codon (SNV-only same-codon seed mode)\n",
            "      --codon-map FILE   BED-like codon map for --mnv-algorithm nirvana-codon:\n",
            "                        CHROM START0 END0 TRANSCRIPT CODON_ID [ignored...]\n",
            "      --min-vars N       Minimum source variants per emitted call (default: 2)\n",
            "      --min-snvs N       Alias for --min-vars\n",
            "      --unsupported-alleles MODE\n",
            "                        Selected unsupported allele policy: skip or fail\n",
            "                        (default: skip)\n",
            "      --phase-from-bam FILE\n",
            "                        Experimental Rust read-backed phasing from indexed BAM/CRAM\n",
            "                        before MNV construction; input GT phase/PS is ignored\n",
            "      --phase-algorithm MODE\n",
            "                        BAM phasing algorithm: mec or greedy (default: mec;\n",
            "                        whatshap is an alias for mec)\n",
            "      --tag TAG         Phasing tag for --emit all-sites and constructed output:\n",
            "                        PS (default) or HP\n",
            "      --only-snvs       Phase only biallelic SNV genotypes from BAM/CRAM\n",
            "      --output-read-list FILE\n",
            "                        Write selected BAM/CRAM reads used for MEC phasing\n",
            "      --ignore-read-groups\n",
            "                        Ignore BAM read-group sample names and use all reads\n",
            "      --use-supplementary\n",
            "                        Include supplementary alignments in BAM phasing\n",
            "      --supplementary-distance N\n",
            "                        Accepted for WhatsHap CLI compatibility; read grouping\n",
            "                        is currently by QNAME regardless of this distance\n",
            "      --phase-realign-overhang N\n",
            "                        REF/ALT realignment overhang for BAM allele detection\n",
            "                        (default: 10, matching WhatsHap reference mode)\n",
            "      --phase-max-coverage N\n",
            "                        Maximum selected read coverage per variant span for MEC\n",
            "                        phasing (default: 15; alias: --phase-internal-downsampling;\n",
            "                        inspired by WhatsHap --internal-downsampling)\n",
            "      --phase-min-mapq N  Minimum read MAPQ for --phase-from-bam (default: 20;\n",
            "                        aliases: --mapping-quality, --mapq)\n",
            "      --phase-min-baseq N Minimum base quality for --phase-from-bam (default: 13)\n",
            "      --warn-on-n        Warn when a selected REF/ALT allele contains N\n",
            "      --no-ref-check     Do not fail when VCF REF differs from FASTA\n",
            "      --no-header        Suppress VCF header\n",
            "  -q, --quiet            Suppress summary on stderr\n",
            "  -h, --help             Show this help\n",
            "\n",
            "Notes:\n",
            "  * Only phased diploid GT (e.g. 0|1, 1|0, 1|1, 1|2) is used.\n",
            "    Unphased, missing, and non-diploid genotypes are skipped.\n",
            "  * Multi-allelic input sites use the ALT allele selected by each\n",
            "    haplotype's GT allele index; unselected ALTs are ignored and output\n",
            "    remains biallelic. Example: GT 1|2 uses ALT1 on haplotype 1 and\n",
            "    ALT2 on haplotype 2.\n",
            "  * Symbolic, breakend, spanning-deletion '*', and non-DNA ALT alleles\n",
            "    are skipped by default and are currently not barriers; use\n",
            "    --unsupported-alleles fail to reject selected unsupported alleles.\n",
            "  * FORMAT/PS is honored when present; variants are only merged within the\n",
            "    same phase set. If PS is absent, the phase separator and proximity\n",
            "    define the merge block.\n",
            "  * --phase-from-bam is a Rust-only experimental phaser inspired by\n",
            "    WhatsHap's read-backed phasing model. The default mec algorithm solves\n",
            "    exact diploid single-sample MEC per connected component after WhatsHap-style\n",
            "    QNAME read-pair grouping, deterministic read selection, interval coverage\n",
            "    capping, bridge rescue, and blank-aware active-read DP; greedy keeps the\n",
            "    earlier pairwise parity heuristic. Phase comparisons should use switch or\n",
            "    pairwise phase metrics because a whole-block 0|1/1|0 orientation flip is\n",
            "    equivalent.\n",
            "  * With the default --max-gap 0 and --mnv-algorithm proximity, only\n",
            "    adjacent phased variants are merged. Pure SNV blocks are TYPE=MNV;\n",
            "    blocks containing indels are TYPE=COMPLEX. nirvana-codon mode only\n",
            "    recomposes phased SNVs sharing a codon key from --codon-map.\n",
            "  * Output format is inferred from -o/--output. BCF output always includes\n",
            "    a VCF/BCF header even if --no-header is set.\n",
            "  * --emit combined keeps the original VCF/BCF header for the selected\n",
            "    sample, appends phase_mnv metadata, and writes constructed MNV/COMPLEX\n",
            "    records plus residual selected-sample input records. Partial\n",
            "    multi-allelic residuals are allele-aware; multi-sample inputs are\n",
            "    subset to the selected sample, with common cohort INFO tags stripped\n",
            "    after subsetting. --phase-from-bam is not supported in combined mode.\n",
            "  * --emit all-sites keeps the original VCF/BCF header via htslib and\n",
            "    appends phase_mnv metadata instead of replacing it.\n",
            "  * --write-index builds a CSI sidecar by default after the output file is\n",
            "    closed. Use --write-index=tbi for a tabix/TBI index on BGZF VCF.\n",
            "    Indexing requires coordinate-sorted .vcf.gz/.vcf.bgz/.bcf output.\n",
            "  * Unless --quiet is set, summary stats go to stderr and include\n",
            "    input/reference/output (output=stdout for VCF stdout), settings,\n",
            "    skip counts, unsupported categories, and N counts.\n"
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

fn parse_unsupported_alleles_policy(s: &str) -> UnsupportedAllelesPolicy {
    match s {
        "skip" => UnsupportedAllelesPolicy::Skip,
        "fail" => UnsupportedAllelesPolicy::Fail,
        _ => die("--unsupported-alleles must be one of: skip, fail"),
    }
}

fn parse_emit_mode(s: &str) -> EmitMode {
    match s {
        "mnv" => EmitMode::Mnv,
        "combined" | "mnv-plus-unmerged" | "merged-plus-unmerged" => EmitMode::Combined,
        "all-sites" | "all" | "phased-vcf" => EmitMode::AllSites,
        _ => die("--emit must be one of: mnv, combined, all-sites"),
    }
}

fn parse_output_index_kind(s: &str) -> OutputIndexKind {
    match s {
        "csi" | "CSI" => OutputIndexKind::Csi,
        "tbi" | "TBI" | "tabix" => OutputIndexKind::Tbi,
        _ => die("--write-index must be csi or tbi"),
    }
}

fn parse_phase_algorithm(s: &str) -> PhaseAlgorithm {
    match s {
        "mec" | "whap" | "whatshap" => PhaseAlgorithm::Mec,
        "greedy" => PhaseAlgorithm::Greedy,
        "hapchat" | "heuristic" => die(
            "--phase-algorithm hapchat/heuristic is not implemented; use mec/whatshap or greedy",
        ),
        _ => die("--phase-algorithm must be one of: mec, whatshap, greedy"),
    }
}

fn parse_phase_tag(s: &str) -> PhaseTag {
    match s {
        "PS" | "ps" => PhaseTag::Ps,
        "HP" | "hp" => PhaseTag::Hp,
        _ => die("--tag must be one of: PS, HP"),
    }
}

fn parse_mnv_algorithm(s: &str) -> MnvAlgorithm {
    match s {
        "proximity" => MnvAlgorithm::Proximity,
        "nirvana-codon" | "codon" => MnvAlgorithm::NirvanaCodon,
        _ => die("--mnv-algorithm must be one of: proximity, nirvana-codon"),
    }
}

fn parse_args() -> Config {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut fasta_path: Option<String> = None;
    let mut output_path: Option<String> = None;
    let mut sample_name: Option<String> = None;
    let mut phase_bam_path: Option<String> = None;
    let mut phase_min_mapq = 20u8;
    let mut phase_min_baseq = 13u8;
    let mut phase_algorithm = PhaseAlgorithm::Mec;
    let mut phase_tag = PhaseTag::Ps;
    let mut phase_max_coverage = 15usize;
    let mut phase_only_snvs = false;
    let mut phase_use_supplementary = false;
    let mut phase_ignore_read_groups = false;
    let mut phase_supplementary_distance = 100_000i64;
    let mut phase_realign_overhang = 10i64;
    let mut phase_read_list_path: Option<String> = None;
    let mut threads = 1usize;
    let mut emit_mode = EmitMode::Mnv;
    let mut output_index: Option<OutputIndexKind> = None;
    let mut mnv_algorithm = MnvAlgorithm::Proximity;
    let mut codon_map_path: Option<String> = None;
    let mut max_gap = 0i64;
    let mut min_variants = 2usize;
    let mut unsupported_alleles = UnsupportedAllelesPolicy::Skip;
    let mut warn_on_n = false;
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
            "-@" | "--threads" => {
                i += 1;
                if i >= args.len() {
                    die("--threads requires an argument");
                }
                let value = parse_i64(&args[i], "--threads");
                if value < 1 {
                    die("--threads must be >= 1");
                }
                threads = value as usize;
            }
            "--emit" | "--emit-mode" => {
                i += 1;
                if i >= args.len() {
                    die("--emit requires an argument");
                }
                emit_mode = parse_emit_mode(&args[i]);
            }
            "--write-index" | "--index" => {
                output_index = Some(OutputIndexKind::Csi);
            }
            _ if arg.starts_with("--write-index=") => {
                let value = arg.trim_start_matches("--write-index=");
                output_index = Some(parse_output_index_kind(value));
            }
            _ if arg.starts_with("--index=") => {
                let value = arg.trim_start_matches("--index=");
                output_index = Some(parse_output_index_kind(value));
            }
            "--mnv-algorithm" => {
                i += 1;
                if i >= args.len() {
                    die("--mnv-algorithm requires an argument");
                }
                mnv_algorithm = parse_mnv_algorithm(&args[i]);
            }
            "--codon-map" => {
                i += 1;
                if i >= args.len() {
                    die("--codon-map requires an argument");
                }
                codon_map_path = Some(args[i].clone());
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
            "--unsupported-alleles" => {
                i += 1;
                if i >= args.len() {
                    die("--unsupported-alleles requires an argument");
                }
                unsupported_alleles = parse_unsupported_alleles_policy(&args[i]);
            }
            "--phase-from-bam" | "--bam" => {
                i += 1;
                if i >= args.len() {
                    die("--phase-from-bam requires an argument");
                }
                phase_bam_path = Some(args[i].clone());
            }
            "--phase-algorithm" | "--algorithm" => {
                i += 1;
                if i >= args.len() {
                    die("--phase-algorithm requires an argument");
                }
                phase_algorithm = parse_phase_algorithm(&args[i]);
            }
            "--tag" => {
                i += 1;
                if i >= args.len() {
                    die("--tag requires an argument");
                }
                phase_tag = parse_phase_tag(&args[i]);
            }
            "--only-snvs" => phase_only_snvs = true,
            "--output-read-list" => {
                i += 1;
                if i >= args.len() {
                    die("--output-read-list requires an argument");
                }
                phase_read_list_path = Some(args[i].clone());
            }
            "--ignore-read-groups" => phase_ignore_read_groups = true,
            "--use-supplementary" => phase_use_supplementary = true,
            "--supplementary-distance" => {
                i += 1;
                if i >= args.len() {
                    die("--supplementary-distance requires an argument");
                }
                let value = parse_i64(&args[i], "--supplementary-distance");
                if value < 0 {
                    die("--supplementary-distance must be >= 0");
                }
                phase_supplementary_distance = value;
            }
            "--phase-realign-overhang" | "--overhang" => {
                i += 1;
                if i >= args.len() {
                    die("--phase-realign-overhang requires an argument");
                }
                let value = parse_i64(&args[i], "--phase-realign-overhang");
                if !(0..=1000).contains(&value) {
                    die("--phase-realign-overhang must be between 0 and 1000");
                }
                phase_realign_overhang = value;
            }
            "--merge-reads" => die("--merge-reads is not implemented; phase_mnv_rs currently groups paired alignments by QNAME but does not implement WhatsHap ReadMerger"),
            "--error-rate" | "--maximum-error-rate" | "--threshold" | "--negative-threshold" => {
                i += 1;
                if i >= args.len() {
                    die(&format!("{} requires an argument", arg));
                }
                die("WhatsHap --merge-reads tuning options are not implemented")
            }
            "--ped" | "--genmap" | "--recombination-list" | "--changed-genotype-list" | "--default-gq" | "--gl-regularizer" | "--row-limit" | "--recombrate" => {
                i += 1;
                if i >= args.len() {
                    die(&format!("{} requires an argument", arg));
                }
                die("pedigree/PedMEC, heuristic row-limit, and genotype-distrust options are not implemented")
            }
            "--distrust-genotypes" | "--include-homozygous" | "--no-genetic-haplotyping" | "--use-ped-samples" => {
                die("pedigree/PedMEC and genotype-distrust options are not implemented")
            }
            "--chromosome" | "--exclude-chromosome" => {
                i += 1;
                if i >= args.len() {
                    die(&format!("{} requires an argument", arg));
                }
                die("WhatsHap-style chromosome include/exclude filters are not implemented; subset the input VCF/BCF upstream")
            }
            "--phase-max-coverage"
            | "--phase-internal-downsampling"
            | "--internal-downsampling" => {
                i += 1;
                if i >= args.len() {
                    die("--phase-max-coverage requires an argument");
                }
                let value = parse_i64(&args[i], "--phase-max-coverage");
                if !(1..=23).contains(&value) {
                    die("--phase-max-coverage must be between 1 and 23");
                }
                phase_max_coverage = value as usize;
            }
            "--phase-min-mapq" | "--mapping-quality" | "--mapq" => {
                i += 1;
                if i >= args.len() {
                    die("--phase-min-mapq requires an argument");
                }
                let value = parse_i64(&args[i], "--phase-min-mapq");
                if !(0..=255).contains(&value) {
                    die("--phase-min-mapq must be between 0 and 255");
                }
                phase_min_mapq = value as u8;
            }
            "--phase-min-baseq" => {
                i += 1;
                if i >= args.len() {
                    die("--phase-min-baseq requires an argument");
                }
                let value = parse_i64(&args[i], "--phase-min-baseq");
                if !(0..=255).contains(&value) {
                    die("--phase-min-baseq must be between 0 and 255");
                }
                phase_min_baseq = value as u8;
            }
            "--warn-on-n" | "--warm-on-n" => warn_on_n = true,
            "--no-ref-check" => no_ref_check = true,
            "--no-reference" => die("--no-reference is not implemented; phase_mnv_rs requires --reference"),
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
    if mnv_algorithm == MnvAlgorithm::NirvanaCodon && codon_map_path.is_none() {
        die("--mnv-algorithm nirvana-codon requires --codon-map");
    }

    Config {
        input_path: positional.remove(0),
        fasta_path: fasta_path.unwrap(),
        output_path,
        sample_name,
        phase_bam_path,
        phase_min_mapq,
        phase_min_baseq,
        phase_algorithm,
        phase_tag,
        phase_max_coverage,
        phase_only_snvs,
        phase_use_supplementary,
        phase_ignore_read_groups,
        phase_supplementary_distance,
        phase_realign_overhang,
        phase_read_list_path,
        threads,
        emit_mode,
        output_index,
        mnv_algorithm,
        codon_map_path,
        max_gap,
        min_variants,
        unsupported_alleles,
        warn_on_n,
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

fn contains_n_base(s: &[u8]) -> bool {
    s.iter().any(|&b| upbase(b) == b'N')
}

fn unsupported_alt_kind(s: &[u8]) -> &'static str {
    if s == b"*" {
        "spanning_deletion"
    } else if is_symbolic_or_breakend(s) {
        "symbolic_or_breakend"
    } else {
        "non_dna"
    }
}

fn phase_model_name(cfg: &Config) -> &'static str {
    match cfg.phase_algorithm {
        PhaseAlgorithm::Mec => "experimental_diploid_mec_dp_whatshap_style_selection",
        PhaseAlgorithm::Greedy => "experimental_read_linkage_greedy_parity",
    }
}

fn output_label(cfg: &Config) -> &str {
    match cfg.output_path.as_deref() {
        None | Some("-") => "stdout",
        Some(path) => path,
    }
}

fn infer_output_kind(cfg: &Config) -> OutputKind {
    match cfg.output_path.as_deref() {
        None | Some("-") => OutputKind::PlainVcf,
        Some(path) => {
            let lower = path.to_ascii_lowercase();
            if lower.ends_with(".bcf") {
                OutputKind::Bcf
            } else if lower.ends_with(".vcf.gz") || lower.ends_with(".vcf.bgz") {
                OutputKind::BgzfVcf
            } else {
                OutputKind::PlainVcf
            }
        }
    }
}

fn index_error_message(code: i32) -> &'static str {
    match code {
        -1 => "indexing failed",
        -2 => "opening output failed",
        -3 => "format not indexable",
        -4 => "failed to create or save the index",
        _ => "unknown indexing error",
    }
}

fn validate_index_request(cfg: &Config) -> Result<(), String> {
    let Some(index_kind) = cfg.output_index else {
        return Ok(());
    };
    let Some(path) = cfg.output_path.as_deref().filter(|p| *p != "-") else {
        return Err(
            "--write-index requires -o/--output with .vcf.gz, .vcf.bgz, or .bcf output".to_string(),
        );
    };
    match (infer_output_kind(cfg), index_kind) {
        (OutputKind::PlainVcf, _) => Err(format!(
            "--write-index cannot index plain VCF output '{path}'; use .vcf.gz, .vcf.bgz, or .bcf"
        )),
        (OutputKind::Bcf, OutputIndexKind::Tbi) => {
            Err("--write-index=tbi is only valid for BGZF VCF; BCF output requires CSI".to_string())
        }
        (OutputKind::BgzfVcf | OutputKind::Bcf, OutputIndexKind::Csi)
        | (OutputKind::BgzfVcf, OutputIndexKind::Tbi) => Ok(()),
    }
}

fn index_output_if_requested(cfg: &Config) -> Result<(), String> {
    let Some(index_kind) = cfg.output_index else {
        return Ok(());
    };
    let path = cfg
        .output_path
        .as_deref()
        .filter(|p| *p != "-")
        .ok_or_else(|| "--write-index requires file output".to_string())?;
    let c_path = CString::new(path.as_bytes())
        .map_err(|_| format!("output path contains NUL byte: {path}"))?;
    let threads = cfg.threads.min(i32::MAX as usize) as i32;
    let ret = unsafe {
        htslib::bcf_index_build3(
            c_path.as_ptr(),
            std::ptr::null(),
            index_kind.min_shift(),
            threads,
        )
    };
    if ret == 0 {
        Ok(())
    } else {
        Err(format!(
            "failed to build {} index for '{}': {} (htslib code {})",
            index_kind.as_str(),
            path,
            index_error_message(ret),
            ret
        ))
    }
}

impl CodonMap {
    fn keys_for(&self, chrom: &str, pos1: i64) -> Vec<String> {
        let Some(intervals) = self.by_chrom.get(chrom) else {
            return Vec::new();
        };
        intervals
            .iter()
            .filter(|iv| iv.start <= pos1 && pos1 <= iv.end)
            .map(|iv| iv.key.clone())
            .collect()
    }
}

fn load_codon_map(path: &str) -> Result<CodonMap, String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("cannot read codon map '{}': {}", path, e))?;
    let mut map = CodonMap::default();
    for (line_no, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields = line.split_whitespace().collect::<Vec<_>>();
        if fields.len() < 5 {
            return Err(format!(
                "invalid codon map '{}': line {} has fewer than 5 fields",
                path,
                line_no + 1
            ));
        }
        let start0 = fields[1].parse::<i64>().map_err(|_| {
            format!(
                "invalid codon map '{}': line {} START0 is not an integer",
                path,
                line_no + 1
            )
        })?;
        let end0 = fields[2].parse::<i64>().map_err(|_| {
            format!(
                "invalid codon map '{}': line {} END0 is not an integer",
                path,
                line_no + 1
            )
        })?;
        if start0 < 0 || end0 <= start0 {
            return Err(format!(
                "invalid codon map '{}': line {} requires 0 <= START0 < END0",
                path,
                line_no + 1
            ));
        }
        let key = format!("{}:{}", fields[3], fields[4]);
        map.by_chrom
            .entry(fields[0].to_string())
            .or_default()
            .push(CodonInterval {
                start: start0 + 1,
                end: end0,
                key,
            });
    }
    for intervals in map.by_chrom.values_mut() {
        intervals.sort_by_key(|iv| (iv.start, iv.end));
    }
    Ok(map)
}

fn load_optional_codon_map(cfg: &Config) -> Result<Option<CodonMap>, String> {
    cfg.codon_map_path
        .as_deref()
        .map(load_codon_map)
        .transpose()
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

    // Keep compact contig output. For normal VCF headers this is the header/RID
    // order used by htslib.
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

fn preflight_vcf_header(path: &str) -> Result<(), String> {
    let c_path = CString::new(path.as_bytes())
        .map_err(|_| format!("input path contains NUL byte: {path}"))?;
    let mode = CString::new("r").expect("literal contains no NUL");
    unsafe {
        let fp = htslib::hts_open(c_path.as_ptr(), mode.as_ptr());
        if fp.is_null() {
            return Err(format!("cannot open input '{path}'"));
        }
        let hdr = htslib::bcf_hdr_read(fp);
        if hdr.is_null() {
            htslib::hts_close(fp);
            return Err("failed to read VCF/BCF header".to_string());
        }
        htslib::bcf_hdr_destroy(hdr);
        if htslib::hts_close(fp) != 0 {
            return Err(format!("failed to close input '{path}' after header check"));
        }
    }
    Ok(())
}

fn is_read_usable_for_phasing(record: &bam::Record, cfg: &Config) -> bool {
    !record.is_unmapped()
        && !record.is_secondary()
        && (cfg.phase_use_supplementary || !record.is_supplementary())
        && !record.is_duplicate()
        && !record.is_quality_check_failed()
        && record.mapq() >= cfg.phase_min_mapq
}

fn bam_read_group_sample_map(header: &bam::HeaderView) -> HashMap<String, String> {
    let header = bam::Header::from_template(header);
    let header_map = header.to_hashmap();
    let mut rg_to_sample = HashMap::new();
    if let Some(read_groups) = header_map.get("RG") {
        for rg in read_groups {
            if let (Some(id), Some(sample)) = (rg.get("ID"), rg.get("SM")) {
                rg_to_sample.insert(id.clone(), sample.clone());
            }
        }
    }
    rg_to_sample
}

fn read_matches_requested_sample(
    record: &bam::Record,
    rg_to_sample: &HashMap<String, String>,
    sample: &str,
    ignore_read_groups: bool,
) -> bool {
    if ignore_read_groups || rg_to_sample.is_empty() {
        return true;
    }
    match record.aux(b"RG") {
        Ok(Aux::String(rg)) => rg_to_sample.get(rg).is_some_and(|sm| sm == sample),
        _ => false,
    }
}

fn merge_phase_call_lists(a: &mut Vec<(usize, u8, u8)>, b: &[(usize, u8, u8)]) {
    for &(idx, allele, qual) in b {
        match a
            .iter_mut()
            .find(|(existing_idx, _, _)| *existing_idx == idx)
        {
            Some((_, existing_allele, existing_qual)) => {
                if *existing_allele == allele {
                    *existing_qual = existing_qual.saturating_add(qual).min(60);
                } else if qual > *existing_qual {
                    *existing_allele = allele;
                    *existing_qual = qual;
                }
            }
            None => a.push((idx, allele, qual)),
        }
    }
    a.sort_by_key(|&(idx, _, _)| idx);
}

fn read_events(record: &bam::Record) -> HashMap<i64, ReadEvent> {
    let mut events: HashMap<i64, ReadEvent> = HashMap::new();
    let seq = record.seq();
    let qual = record.qual();
    let mut ref_pos = record.pos();
    let mut read_pos = 0usize;
    let mut last_ref_pos: Option<i64> = None;

    for cigar in record.cigar().iter() {
        match *cigar {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                for _ in 0..len {
                    if read_pos < seq.len() {
                        let base = upbase(seq[read_pos]);
                        let q = qual.get(read_pos).copied().unwrap_or(255);
                        let event = events.entry(ref_pos).or_default();
                        event.base = Some(base);
                        event.qual = Some(q);
                    }
                    last_ref_pos = Some(ref_pos);
                    ref_pos += 1;
                    read_pos += 1;
                }
            }
            Cigar::Ins(len) => {
                if let Some(anchor) = last_ref_pos {
                    let event = events.entry(anchor).or_default();
                    for _ in 0..len {
                        if read_pos < seq.len() {
                            let base = upbase(seq[read_pos]);
                            let q = qual.get(read_pos).copied().unwrap_or(255);
                            event.insertion_after.push((base, q));
                        }
                        read_pos += 1;
                    }
                } else {
                    read_pos += len as usize;
                }
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                for _ in 0..len {
                    events.entry(ref_pos).or_default();
                    last_ref_pos = Some(ref_pos);
                    ref_pos += 1;
                }
            }
            Cigar::SoftClip(len) => read_pos += len as usize,
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    events
}

fn phase_quality(q: u8) -> u8 {
    if q == 255 {
        60
    } else {
        q.clamp(1, 60)
    }
}

fn reference_end0(record: &bam::Record) -> i64 {
    let mut ref_pos = record.pos();
    for cigar in record.cigar().iter() {
        match *cigar {
            Cigar::Match(len)
            | Cigar::Equal(len)
            | Cigar::Diff(len)
            | Cigar::Del(len)
            | Cigar::RefSkip(len) => ref_pos += i64::from(len),
            Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    ref_pos
}

fn query_seq_for_ref_window(
    record: &bam::Record,
    start0: i64,
    end0: i64,
    min_baseq: u8,
) -> Option<(String, u8)> {
    if start0 >= end0 {
        return None;
    }
    let seq = record.seq();
    let qual = record.qual();
    let mut ref_pos = record.pos();
    let mut read_pos = 0usize;
    let mut last_ref_pos: Option<i64> = None;
    let mut query = Vec::new();
    let mut min_q = 255u8;
    let mut overlapped = false;

    for cigar in record.cigar().iter() {
        match *cigar {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                for _ in 0..len {
                    if read_pos < seq.len() && start0 <= ref_pos && ref_pos < end0 {
                        let q = qual.get(read_pos).copied().unwrap_or(255);
                        if q < min_baseq {
                            return None;
                        }
                        min_q = min_q.min(q);
                        query.push(upbase(seq[read_pos]));
                        overlapped = true;
                    }
                    last_ref_pos = Some(ref_pos);
                    ref_pos += 1;
                    read_pos += 1;
                }
            }
            Cigar::Ins(len) => {
                let keep = last_ref_pos.is_some_and(|anchor| start0 <= anchor && anchor < end0);
                for _ in 0..len {
                    if read_pos < seq.len() && keep {
                        let q = qual.get(read_pos).copied().unwrap_or(255);
                        if q < min_baseq {
                            return None;
                        }
                        min_q = min_q.min(q);
                        query.push(upbase(seq[read_pos]));
                    }
                    read_pos += 1;
                }
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                for _ in 0..len {
                    if start0 <= ref_pos && ref_pos < end0 {
                        overlapped = true;
                    }
                    last_ref_pos = Some(ref_pos);
                    ref_pos += 1;
                }
            }
            Cigar::SoftClip(len) => read_pos += len as usize,
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    if overlapped {
        Some((
            String::from_utf8_lossy(&query).into_owned(),
            phase_quality(min_q).min(30),
        ))
    } else {
        None
    }
}

fn edit_distance(a: &[u8], b: &[u8]) -> usize {
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    let mut curr = vec![0usize; b.len() + 1];
    for (i, &ac) in a.iter().enumerate() {
        curr[0] = i + 1;
        for (j, &bc) in b.iter().enumerate() {
            let cost = usize::from(ac.to_ascii_uppercase() != bc.to_ascii_uppercase());
            curr[j + 1] = (prev[j + 1] + 1).min(curr[j] + 1).min(prev[j] + cost);
        }
        std::mem::swap(&mut prev, &mut curr);
    }
    prev[b.len()]
}

fn realigned_read_allele_for_candidate(
    record: &bam::Record,
    fai: *mut htslib::faidx_t,
    candidate: &PhaseCandidate,
    min_baseq: u8,
    overhang: i64,
) -> Option<(String, u8)> {
    let ref_len = candidate.alleles.first()?.len() as i64;
    let var_start0 = candidate.pos - 1;
    let var_end0 = var_start0 + ref_len;
    let read_start0 = record.pos();
    let read_end0 = reference_end0(record);
    if read_start0 > var_start0 || read_end0 < var_end0 {
        return None;
    }
    let win_start0 = read_start0.max(var_start0.saturating_sub(overhang));
    let win_end0 = read_end0.min(var_end0.saturating_add(overhang));
    if win_start0 > var_start0 || win_end0 < var_end0 || win_start0 >= win_end0 {
        return None;
    }

    let ref_window = fetch_seq(fai, &candidate.chrom, win_start0 + 1, win_end0).ok()?;
    let left_len = (var_start0 - win_start0) as usize;
    let right_start = (var_end0 - win_start0) as usize;
    let left_pad = &ref_window[..left_len];
    let right_pad = &ref_window[right_start..];
    let (query, qual) = query_seq_for_ref_window(record, win_start0, win_end0, min_baseq)?;

    let allele_indices = candidate.input_order_alleles;
    let mut distances = Vec::new();
    for allele_index in allele_indices {
        if allele_index < 0 || allele_index as usize >= candidate.alleles.len() {
            return None;
        }
        let allele = format!(
            "{}{}{}",
            left_pad, candidate.alleles[allele_index as usize], right_pad
        );
        distances.push((
            allele_index,
            edit_distance(query.as_bytes(), allele.as_bytes()),
        ));
    }
    distances.sort_by_key(|&(_, distance)| distance);
    if distances.len() >= 2 && distances[0].1 < distances[1].1 {
        Some((candidate.alleles[distances[0].0 as usize].clone(), qual))
    } else {
        None
    }
}

fn read_allele_for_candidate(
    events: &HashMap<i64, ReadEvent>,
    candidate: &PhaseCandidate,
    min_baseq: u8,
) -> Option<(String, u8)> {
    let mut allele = Vec::new();
    let mut min_q = 255u8;
    let start0 = candidate.pos - 1;
    let ref_len = candidate.alleles.first()?.len() as i64;
    for pos0 in start0..start0 + ref_len {
        let event = events.get(&pos0)?;
        if let Some(base) = event.base {
            let q = event.qual.unwrap_or(0);
            if q < min_baseq {
                return None;
            }
            min_q = min_q.min(q);
            allele.push(base);
        }
        for &(base, q) in &event.insertion_after {
            if q < min_baseq {
                return None;
            }
            min_q = min_q.min(q);
            allele.push(base);
        }
    }
    if allele.is_empty() {
        None
    } else {
        Some((
            String::from_utf8_lossy(&allele).into_owned(),
            phase_quality(min_q),
        ))
    }
}

fn collect_read_phase_calls(
    record: &bam::Record,
    candidates: &[PhaseCandidate],
    candidates_by_pos0: &HashMap<i64, Vec<usize>>,
    fai: *mut htslib::faidx_t,
    min_baseq: u8,
    realign_overhang: i64,
) -> Vec<(usize, u8, u8)> {
    let events = read_events(record);
    if events.is_empty() {
        return Vec::new();
    }

    let mut calls = Vec::new();
    for &pos0 in events.keys() {
        let Some(indices) = candidates_by_pos0.get(&pos0) else {
            continue;
        };
        for &idx in indices {
            let candidate = &candidates[idx];
            let Some((read_allele, qual)) = realigned_read_allele_for_candidate(
                record,
                fai,
                candidate,
                min_baseq,
                realign_overhang,
            )
            .or_else(|| read_allele_for_candidate(&events, candidate, min_baseq)) else {
                continue;
            };
            let allele0 = &candidate.alleles[candidate.input_order_alleles[0] as usize];
            let allele1 = &candidate.alleles[candidate.input_order_alleles[1] as usize];
            if read_allele.eq_ignore_ascii_case(allele0) {
                calls.push((idx, 0, qual));
            } else if read_allele.eq_ignore_ascii_case(allele1) {
                calls.push((idx, 1, qual));
            }
        }
    }

    calls.sort_by_key(|&(idx, _, _)| idx);
    calls.dedup_by_key(|(idx, _, _)| *idx);
    calls
}

fn collect_bam_phase_read_calls(
    cfg: &Config,
    candidates: &[PhaseCandidate],
    sample: &str,
    st: &mut Stats,
) -> Result<Vec<ReadPhaseCall>, String> {
    let Some(bam_path) = cfg.phase_bam_path.as_deref() else {
        return Ok(Vec::new());
    };

    let mut bam = bam::IndexedReader::from_path(bam_path)
        .map_err(|e| format!("cannot open indexed BAM/CRAM '{bam_path}': {e}"))?;
    if cfg.threads > 1 {
        bam.set_threads(cfg.threads)
            .map_err(|e| format!("failed to enable BAM/CRAM threads for '{bam_path}': {e}"))?;
    }

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

    let rg_to_sample = bam_read_group_sample_map(bam.header());
    let mut by_chrom: HashMap<&str, Vec<usize>> = HashMap::new();
    for (idx, candidate) in candidates.iter().enumerate() {
        by_chrom
            .entry(candidate.chrom.as_str())
            .or_default()
            .push(idx);
    }

    let mut read_calls = Vec::new();
    let mut record = bam::Record::new();

    for (chrom, indices) in by_chrom {
        let beg = indices
            .iter()
            .map(|&idx| candidates[idx].pos - 1)
            .min()
            .unwrap_or(0)
            .max(0);
        let end = indices
            .iter()
            .map(|&idx| candidates[idx].end)
            .max()
            .unwrap_or(beg + 1);
        bam.fetch((chrom.as_bytes(), beg, end))
            .map_err(|e| format!("failed to fetch {chrom}:{beg}-{end} from '{bam_path}': {e}"))?;

        let mut candidates_by_pos0: HashMap<i64, Vec<usize>> = HashMap::new();
        for &idx in &indices {
            candidates_by_pos0
                .entry(candidates[idx].pos - 1)
                .or_default()
                .push(idx);
        }

        let mut grouped: HashMap<String, ReadPhaseCall> = HashMap::new();
        let mut group_order = Vec::new();
        while let Some(result) = bam.read(&mut record) {
            result.map_err(|e| format!("failed to read BAM/CRAM record from '{bam_path}': {e}"))?;
            if !is_read_usable_for_phasing(&record, cfg) {
                continue;
            }
            if !read_matches_requested_sample(
                &record,
                &rg_to_sample,
                sample,
                cfg.phase_ignore_read_groups,
            ) {
                continue;
            }
            let calls = collect_read_phase_calls(
                &record,
                candidates,
                &candidates_by_pos0,
                fai.0,
                cfg.phase_min_baseq,
                cfg.phase_realign_overhang,
            );
            if calls.is_empty() {
                continue;
            }
            let name = String::from_utf8_lossy(record.qname()).into_owned();
            match grouped.get_mut(&name) {
                Some(existing) => {
                    merge_phase_call_lists(&mut existing.calls, &calls);
                    existing.mapq = existing.mapq.max(record.mapq());
                }
                None => {
                    let order = group_order.len();
                    group_order.push(name.clone());
                    grouped.insert(
                        name.clone(),
                        ReadPhaseCall {
                            calls,
                            order,
                            name,
                            mapq: record.mapq(),
                        },
                    );
                }
            }
        }
        for name in group_order {
            if let Some(read) = grouped.remove(&name) {
                if read.calls.len() >= 2 {
                    st.bam_phase_informative_reads += 1;
                    read_calls.push(ReadPhaseCall {
                        order: read_calls.len(),
                        ..read
                    });
                }
            }
        }
    }

    Ok(read_calls)
}

fn assignments_from_dsu(
    candidates: &[PhaseCandidate],
    mut dsu: Dsu,
    st: &mut Stats,
) -> HashMap<usize, PhaseAssignment> {
    let mut components: HashMap<usize, Vec<usize>> = HashMap::new();
    for idx in 0..candidates.len() {
        let (root, _) = dsu.find(idx);
        components.entry(root).or_default().push(idx);
    }

    let mut assignments = HashMap::new();
    for members in components.values() {
        if members.len() < 2 {
            continue;
        }
        st.bam_phase_components += 1;
        let anchor = *members
            .iter()
            .min_by_key(|&&idx| (candidates[idx].pos, idx))
            .expect("non-empty component");
        let ps = candidates[anchor].pos;
        let (_, anchor_xor) = dsu.find(anchor);
        for &idx in members {
            let (_, x) = dsu.find(idx);
            assignments.insert(
                idx,
                PhaseAssignment {
                    ps,
                    rel: (x ^ anchor_xor) & 1,
                },
            );
        }
    }
    assignments
}

fn phase_candidates_greedy(
    candidates: &[PhaseCandidate],
    read_calls: &[ReadPhaseCall],
    st: &mut Stats,
) -> HashMap<usize, PhaseAssignment> {
    let mut pair_scores: HashMap<(usize, usize), i64> = HashMap::new();
    for read in read_calls {
        for i in 0..read.calls.len() {
            for j in i + 1..read.calls.len() {
                let (a, sa, qa) = read.calls[i];
                let (b, sb, qb) = read.calls[j];
                let key = if a < b { (a, b) } else { (b, a) };
                let parity = (sa ^ sb) & 1;
                let weight = i64::from(qa.min(qb));
                let score = pair_scores.entry(key).or_insert(0);
                if parity == 0 {
                    *score += weight;
                } else {
                    *score -= weight;
                }
            }
        }
    }

    let mut constraints: Vec<(i64, usize, usize, u8)> = pair_scores
        .into_iter()
        .filter_map(|((a, b), score)| {
            if score > 0 {
                Some((score, a, b, 0))
            } else if score < 0 {
                Some((-score, a, b, 1))
            } else {
                None
            }
        })
        .collect();
    constraints.sort_by(|x, y| {
        y.0.cmp(&x.0)
            .then_with(|| x.1.cmp(&y.1))
            .then_with(|| x.2.cmp(&y.2))
    });

    let mut dsu = Dsu::new(candidates.len());
    for (_weight, a, b, parity) in constraints {
        if !dsu.union(a, b, parity) {
            st.bam_phase_conflicts += 1;
        }
    }
    assignments_from_dsu(candidates, dsu, st)
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
struct ReadSelectionScore {
    new_score: i64,
    total_score: i64,
    min_qual: u8,
}

impl Ord for ReadSelectionScore {
    fn cmp(&self, other: &Self) -> Ordering {
        self.new_score
            .cmp(&other.new_score)
            .then_with(|| self.total_score.cmp(&other.total_score))
            .then_with(|| self.min_qual.cmp(&other.min_qual))
    }
}

impl PartialOrd for ReadSelectionScore {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
struct ReadSelectionHeapEntry {
    score: ReadSelectionScore,
    idx: usize,
    order: usize,
}

impl Ord for ReadSelectionHeapEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        self.score
            .cmp(&other.score)
            .then_with(|| other.order.cmp(&self.order))
            .then_with(|| other.idx.cmp(&self.idx))
    }
}

impl PartialOrd for ReadSelectionHeapEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn read_call_span(read: &ReadPhaseCall) -> Option<(usize, usize)> {
    let mut begin = usize::MAX;
    let mut end = 0usize;
    for &(var, _, _) in &read.calls {
        begin = begin.min(var);
        end = end.max(var + 1);
    }
    if begin == usize::MAX {
        None
    } else {
        Some((begin, end))
    }
}

fn read_selection_score(read: &ReadPhaseCall, covered: &[bool]) -> ReadSelectionScore {
    let Some((begin, end)) = read_call_span(read) else {
        return ReadSelectionScore {
            new_score: i64::MIN,
            total_score: i64::MIN,
            min_qual: 0,
        };
    };
    let span_len = end.saturating_sub(begin);
    let gaps = span_len.saturating_sub(read.calls.len());
    let new_observed = read
        .calls
        .iter()
        .filter(|&&(var, _, _)| var < covered.len() && !covered[var])
        .count();
    let min_qual = read.calls.iter().map(|&(_, _, q)| q).min().unwrap_or(0);
    ReadSelectionScore {
        new_score: new_observed as i64 - gaps as i64,
        total_score: read.calls.len() as i64 - gaps as i64,
        min_qual,
    }
}

fn read_selection_entry(
    read_calls: &[ReadPhaseCall],
    idx: usize,
    covered: &[bool],
) -> ReadSelectionHeapEntry {
    ReadSelectionHeapEntry {
        score: read_selection_score(&read_calls[idx], covered),
        idx,
        order: read_calls[idx].order,
    }
}

fn read_covers_new_observation(read: &ReadPhaseCall, covered: &[bool]) -> bool {
    read.calls
        .iter()
        .any(|&(var, _, _)| var < covered.len() && !covered[var])
}

fn read_span_violates_coverage(
    read: &ReadPhaseCall,
    coverage: &[usize],
    max_coverage: usize,
) -> bool {
    let Some((begin, end)) = read_call_span(read) else {
        return true;
    };
    if end > coverage.len() {
        return true;
    }
    coverage[begin..end].iter().any(|&cov| cov >= max_coverage)
}

fn add_read_span_coverage(read: &ReadPhaseCall, coverage: &mut [usize]) {
    if let Some((begin, end)) = read_call_span(read) {
        if end <= coverage.len() {
            for cov in &mut coverage[begin..end] {
                *cov += 1;
            }
        }
    }
}

fn mark_read_observations_covered(read: &ReadPhaseCall, covered: &mut [bool]) {
    for &(var, _, _) in &read.calls {
        if var < covered.len() {
            covered[var] = true;
        }
    }
}

fn union_read_variants(dsu: &mut Dsu, read: &ReadPhaseCall) {
    let Some(&(first, _, _)) = read.calls.first() else {
        return;
    };
    for &(var, _, _) in read.calls.iter().skip(1) {
        let _ = dsu.union(first, var, 0);
    }
}

fn read_component_count(dsu: &mut Dsu, read: &ReadPhaseCall) -> usize {
    let mut roots = Vec::new();
    for &(var, _, _) in &read.calls {
        let (root, _) = dsu.find(var);
        if !roots.contains(&root) {
            roots.push(root);
        }
    }
    roots.len()
}

fn select_read_slice(
    read_calls: &[ReadPhaseCall],
    variant_count: usize,
    max_coverage: usize,
    coverage: &mut [usize],
    undecided: &mut [bool],
) -> Vec<usize> {
    let mut covered = vec![false; variant_count];
    let mut heap = BinaryHeap::new();
    for idx in 0..read_calls.len() {
        if undecided[idx] {
            heap.push(read_selection_entry(read_calls, idx, &covered));
        }
    }

    let mut selected = Vec::new();
    while let Some(entry) = heap.pop() {
        if !undecided[entry.idx] {
            continue;
        }
        let current = read_selection_entry(read_calls, entry.idx, &covered);
        if current.score != entry.score {
            heap.push(current);
            continue;
        }
        if read_span_violates_coverage(&read_calls[entry.idx], coverage, max_coverage) {
            undecided[entry.idx] = false;
            continue;
        }
        if read_covers_new_observation(&read_calls[entry.idx], &covered) {
            undecided[entry.idx] = false;
            add_read_span_coverage(&read_calls[entry.idx], coverage);
            mark_read_observations_covered(&read_calls[entry.idx], &mut covered);
            selected.push(entry.idx);
        }
    }
    selected
}

fn select_bridging_reads(
    read_calls: &[ReadPhaseCall],
    variant_count: usize,
    max_coverage: usize,
    coverage: &mut [usize],
    undecided: &mut [bool],
    reads_in_slice: &[usize],
) -> Vec<usize> {
    if reads_in_slice.is_empty() {
        return Vec::new();
    }

    let mut dsu = Dsu::new(variant_count);
    for &idx in reads_in_slice {
        union_read_variants(&mut dsu, &read_calls[idx]);
    }

    let uncovered = vec![false; variant_count];
    let mut order: Vec<usize> = (0..read_calls.len())
        .filter(|&idx| undecided[idx])
        .collect();
    order.sort_by(|&a, &b| {
        read_selection_score(&read_calls[b], &uncovered)
            .cmp(&read_selection_score(&read_calls[a], &uncovered))
            .then_with(|| read_calls[a].order.cmp(&read_calls[b].order))
    });

    let mut bridged = Vec::new();
    for idx in order {
        if !undecided[idx] {
            continue;
        }
        if read_span_violates_coverage(&read_calls[idx], coverage, max_coverage) {
            undecided[idx] = false;
            continue;
        }
        if read_component_count(&mut dsu, &read_calls[idx]) < 2 {
            continue;
        }
        undecided[idx] = false;
        add_read_span_coverage(&read_calls[idx], coverage);
        union_read_variants(&mut dsu, &read_calls[idx]);
        bridged.push(idx);
    }
    bridged
}

fn select_read_calls_max_coverage(
    read_calls: &[ReadPhaseCall],
    variant_count: usize,
    max_coverage: usize,
) -> Vec<ReadPhaseCall> {
    let mut coverage = vec![0usize; variant_count];
    let mut undecided = vec![true; read_calls.len()];
    let mut selected_flag = vec![false; read_calls.len()];
    let mut selected_indices = Vec::new();

    while undecided.iter().any(|&x| x) {
        let undecided_before = undecided.iter().filter(|&&x| x).count();
        let reads_in_slice = select_read_slice(
            read_calls,
            variant_count,
            max_coverage,
            &mut coverage,
            &mut undecided,
        );
        for &idx in &reads_in_slice {
            if !selected_flag[idx] {
                selected_flag[idx] = true;
                selected_indices.push(idx);
            }
        }

        let bridging_reads = select_bridging_reads(
            read_calls,
            variant_count,
            max_coverage,
            &mut coverage,
            &mut undecided,
            &reads_in_slice,
        );
        for &idx in &bridging_reads {
            if !selected_flag[idx] {
                selected_flag[idx] = true;
                selected_indices.push(idx);
            }
        }

        let undecided_after = undecided.iter().filter(|&&x| x).count();
        if undecided_after == undecided_before {
            break;
        }
    }

    selected_indices.sort_by_key(|&idx| read_calls[idx].order);
    selected_indices
        .into_iter()
        .map(|idx| read_calls[idx].clone())
        .collect()
}

fn write_phase_read_list(path: &str, selected: &[ReadPhaseCall]) -> Result<(), String> {
    let file = File::create(path).map_err(|e| format!("cannot create read list '{path}': {e}"))?;
    let mut out = BufWriter::new(file);
    writeln!(
        out,
        "#read_name\torder\tmapq\tn_calls\tfirst_variant_index\tlast_variant_index\tcalls"
    )
    .map_err(|e| format!("failed to write read list '{path}': {e}"))?;
    for read in selected {
        let Some((begin, end)) = read_call_span(read) else {
            continue;
        };
        let calls = read
            .calls
            .iter()
            .map(|(idx, allele, qual)| format!("{idx}:{allele}:{qual}"))
            .collect::<Vec<_>>()
            .join(",");
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            read.name,
            read.order,
            read.mapq,
            read.calls.len(),
            begin,
            end.saturating_sub(1),
            calls
        )
        .map_err(|e| format!("failed to write read list '{path}': {e}"))?;
    }
    out.flush()
        .map_err(|e| format!("failed to flush read list '{path}': {e}"))?;
    Ok(())
}

fn project_mask(from_active: &[usize], to_pos: &HashMap<usize, usize>, from_mask: u64) -> u64 {
    let mut key = 0u64;
    for (from_bit, read_id) in from_active.iter().enumerate() {
        if let Some(&to_bit) = to_pos.get(read_id) {
            if (from_mask & (1u64 << from_bit)) != 0 {
                key |= 1u64 << to_bit;
            }
        }
    }
    key
}

fn mec_column_cost(
    mask: u64,
    observations: &[(usize, u8, u8)],
    active_pos: &HashMap<usize, usize>,
) -> (u64, u8, bool) {
    let mut costs = [0u64; 2];
    for orientation in 0..=1u8 {
        let mut cost = 0u64;
        for &(read_id, allele, qual) in observations {
            let bit = ((mask >> active_pos[&read_id]) & 1) as u8;
            let expected = orientation ^ bit;
            if allele != expected {
                cost += u64::from(qual);
            }
        }
        costs[orientation as usize] = cost;
    }
    if costs[0] <= costs[1] {
        (costs[0], 0, costs[0] == costs[1])
    } else {
        (costs[1], 1, false)
    }
}

fn solve_mec_component(members: &[usize], reads: &[ReadPhaseCall]) -> HashMap<usize, u8> {
    let mut sorted_members = members.to_vec();
    sorted_members.sort_unstable();
    let mut global_to_local = HashMap::new();
    for (local, &global) in sorted_members.iter().enumerate() {
        global_to_local.insert(global, local);
    }

    let mut local_reads = Vec::new();
    for read in reads {
        let mut calls = read
            .calls
            .iter()
            .filter_map(|&(global, allele, qual)| {
                global_to_local
                    .get(&global)
                    .map(|&var| LocalReadCall { var, allele, qual })
            })
            .collect::<Vec<_>>();
        if calls.is_empty() {
            continue;
        }
        calls.sort_by_key(|c| c.var);
        let first_var = calls.first().map(|c| c.var).unwrap_or(0);
        let last_var = calls.last().map(|c| c.var).unwrap_or(first_var);
        local_reads.push(LocalRead {
            calls,
            first_var,
            last_var,
        });
    }

    let mut by_var: Vec<Vec<(usize, u8, u8)>> = vec![Vec::new(); sorted_members.len()];
    for (read_id, read) in local_reads.iter().enumerate() {
        for call in &read.calls {
            by_var[call.var].push((read_id, call.allele, call.qual));
        }
    }

    let mut traces: Vec<MecTraceColumn> = Vec::new();
    for (var, observations) in by_var.iter().enumerate() {
        let active_reads = local_reads
            .iter()
            .enumerate()
            .filter_map(|(read_id, read)| {
                if read.first_var <= var && var <= read.last_var {
                    Some(read_id)
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        if active_reads.len() >= u64::BITS as usize {
            return HashMap::new();
        }
        let active_pos = active_reads
            .iter()
            .enumerate()
            .map(|(i, &read_id)| (read_id, i))
            .collect::<HashMap<_, _>>();
        let state_count = 1u64 << active_reads.len();
        let mut states = HashMap::new();

        if let Some(prev) = traces.last() {
            let mut prev_lookup: HashMap<u64, (u64, u64)> = HashMap::new();
            for (&prev_mask, prev_state) in &prev.states {
                let key = project_mask(&prev.active_reads, &active_pos, prev_mask);
                match prev_lookup.get_mut(&key) {
                    Some((best_cost, best_mask)) if prev_state.cost < *best_cost => {
                        *best_cost = prev_state.cost;
                        *best_mask = prev_mask;
                    }
                    None => {
                        prev_lookup.insert(key, (prev_state.cost, prev_mask));
                    }
                    _ => {}
                }
            }
            let shared_mask = active_reads
                .iter()
                .enumerate()
                .fold(0u64, |acc, (bit, read_id)| {
                    if prev.active_reads.contains(read_id) {
                        acc | (1u64 << bit)
                    } else {
                        acc
                    }
                });
            for mask in 0..state_count {
                let key = mask & shared_mask;
                let Some(&(prev_cost, prev_mask)) = prev_lookup.get(&key) else {
                    continue;
                };
                let (column_cost, orientation, ambiguous) =
                    mec_column_cost(mask, observations, &active_pos);
                states.insert(
                    mask,
                    MecTraceState {
                        cost: prev_cost + column_cost,
                        prev_mask: Some(prev_mask),
                        orientation,
                        ambiguous,
                    },
                );
            }
        } else {
            for mask in 0..state_count {
                let (column_cost, orientation, ambiguous) =
                    mec_column_cost(mask, observations, &active_pos);
                states.insert(
                    mask,
                    MecTraceState {
                        cost: column_cost,
                        prev_mask: None,
                        orientation,
                        ambiguous,
                    },
                );
            }
        }
        traces.push(MecTraceColumn {
            active_reads,
            states,
        });
    }

    let mut rel_by_local = vec![None; sorted_members.len()];
    let Some(last) = traces.last() else {
        return HashMap::new();
    };
    let Some((&mut_mask, _)) = last
        .states
        .iter()
        .min_by(|(ma, sa), (mb, sb)| sa.cost.cmp(&sb.cost).then_with(|| ma.cmp(mb)))
    else {
        return HashMap::new();
    };
    let mut mask = mut_mask;
    for col in (0..traces.len()).rev() {
        let state = traces[col].states[&mask];
        if !state.ambiguous {
            rel_by_local[col] = Some(state.orientation);
        }
        if let Some(prev_mask) = state.prev_mask {
            mask = prev_mask;
        }
    }

    sorted_members
        .into_iter()
        .enumerate()
        .filter_map(|(local, global)| rel_by_local[local].map(|rel| (global, rel)))
        .collect()
}

fn phase_candidates_mec(
    cfg: &Config,
    candidates: &[PhaseCandidate],
    read_calls: &[ReadPhaseCall],
    st: &mut Stats,
) -> Result<HashMap<usize, PhaseAssignment>, String> {
    let selected =
        select_read_calls_max_coverage(read_calls, candidates.len(), cfg.phase_max_coverage);
    st.bam_phase_selected_reads = selected.len() as u64;
    if let Some(path) = cfg.phase_read_list_path.as_deref() {
        write_phase_read_list(path, &selected)?;
    }

    let mut dsu = Dsu::new(candidates.len());
    for read in &selected {
        let first = read.calls[0].0;
        for &(idx, _, _) in read.calls.iter().skip(1) {
            let _ = dsu.union(first, idx, 0);
        }
    }

    let mut components: HashMap<usize, Vec<usize>> = HashMap::new();
    for read in &selected {
        for &(idx, _, _) in &read.calls {
            let (root, _) = dsu.find(idx);
            components.entry(root).or_default().push(idx);
        }
    }
    for members in components.values_mut() {
        members.sort_unstable();
        members.dedup();
    }

    let mut assignments = HashMap::new();
    for members in components.values() {
        if members.len() < 2 {
            continue;
        }
        st.bam_phase_components += 1;
        let rels = solve_mec_component(members, &selected);
        let anchor = *members
            .iter()
            .min_by_key(|&&idx| (candidates[idx].pos, idx))
            .expect("non-empty component");
        let ps = candidates[anchor].pos;
        let anchor_rel = rels.get(&anchor).copied().unwrap_or(0);
        for &idx in members {
            if let Some(rel) = rels.get(&idx).copied() {
                assignments.insert(
                    idx,
                    PhaseAssignment {
                        ps,
                        rel: (rel ^ anchor_rel) & 1,
                    },
                );
            }
        }
    }
    Ok(assignments)
}

fn phase_candidates_from_bam(
    cfg: &Config,
    candidates: &[PhaseCandidate],
    sample: &str,
    st: &mut Stats,
) -> Result<HashMap<usize, PhaseAssignment>, String> {
    if cfg.phase_bam_path.is_none() || candidates.is_empty() {
        return Ok(HashMap::new());
    }
    let read_calls = collect_bam_phase_read_calls(cfg, candidates, sample, st)?;
    let assignments = match cfg.phase_algorithm {
        PhaseAlgorithm::Mec => phase_candidates_mec(cfg, candidates, &read_calls, st)?,
        PhaseAlgorithm::Greedy => {
            if let Some(path) = cfg.phase_read_list_path.as_deref() {
                write_phase_read_list(path, &read_calls)?;
            }
            phase_candidates_greedy(candidates, &read_calls, st)
        }
    };
    st.bam_phase_phased_variants = assignments.len() as u64;
    st.bam_phase_unphased_variants = candidates.len().saturating_sub(assignments.len()) as u64;
    st.skipped_unphased += st.bam_phase_unphased_variants;
    Ok(assignments)
}

fn candidate_from_record(
    record: &bcf::Record,
    header_info: &HeaderInfo,
    sample_idx: usize,
    source_record: u64,
    st: &mut Stats,
    cfg: &Config,
) -> Result<Option<PhaseCandidate>, String> {
    if record.alleles().len() > 2 {
        st.multiallelic_records += 1;
    }
    let genotypes = match record.genotypes() {
        Ok(g) => g,
        Err(_) => {
            st.skipped_no_gt += 1;
            return Ok(None);
        }
    };
    let gt = genotypes.get(sample_idx);
    if gt.len() != 2 {
        st.skipped_not_diploid += 1;
        return Ok(None);
    }
    let Some(a0) = allele_index(gt[0]) else {
        st.skipped_missing_gt += 1;
        return Ok(None);
    };
    let Some(a1) = allele_index(gt[1]) else {
        st.skipped_missing_gt += 1;
        return Ok(None);
    };

    let alleles_raw = record.alleles();
    if alleles_raw.is_empty() {
        st.skipped_ref += 1;
        if cfg.unsupported_alleles == UnsupportedAllelesPolicy::Fail {
            return Err("unsupported empty REF allele".to_string());
        }
        return Ok(None);
    }

    let pos1 = record.pos() + 1;
    let rid = record.rid().unwrap_or(0) as usize;
    let chrom = header_info
        .rid_names
        .get(rid)
        .map(String::as_str)
        .unwrap_or(".");
    if !is_plain_dna_allele(alleles_raw[0]) {
        st.skipped_ref += 1;
        if cfg.unsupported_alleles == UnsupportedAllelesPolicy::Fail {
            let ref_allele = uppercase_ascii_string(alleles_raw[0]);
            return Err(format!(
                "unsupported REF allele at {chrom}:{pos1} REF={ref_allele}"
            ));
        }
        return Ok(None);
    }

    for allele in [a0, a1] {
        if allele < 0 || allele as usize >= alleles_raw.len() {
            st.skipped_unsupported_alt += 1;
            st.skipped_alt_out_of_range += 1;
            if cfg.unsupported_alleles == UnsupportedAllelesPolicy::Fail {
                return Err(format!(
                    "unsupported selected allele index at {chrom}:{pos1} allele={allele}"
                ));
            }
            return Ok(None);
        }
        if allele != 0 && !is_plain_dna_allele(alleles_raw[allele as usize]) {
            st.skipped_unsupported_alt += 1;
            let kind = unsupported_alt_kind(alleles_raw[allele as usize]);
            match kind {
                "symbolic_or_breakend" => st.skipped_alt_symbolic_or_breakend += 1,
                "spanning_deletion" => st.skipped_alt_spanning_deletion += 1,
                _ => st.skipped_alt_non_dna += 1,
            }
            if cfg.unsupported_alleles == UnsupportedAllelesPolicy::Fail {
                let alt_allele = uppercase_ascii_string(alleles_raw[allele as usize]);
                return Err(format!(
                    "unsupported selected allele at {chrom}:{pos1} ALT={alt_allele} kind={kind}"
                ));
            }
            return Ok(None);
        }
    }

    if a0 == a1 {
        st.skipped_unphased += 1;
        return Ok(None);
    }
    if !is_plain_dna_allele(alleles_raw[a0 as usize])
        || !is_plain_dna_allele(alleles_raw[a1 as usize])
    {
        return Ok(None);
    }

    let alleles = alleles_raw
        .iter()
        .map(|a| uppercase_ascii_string(a))
        .collect::<Vec<_>>();
    let ref_len = alleles[0].len() as i64;
    let is_snv = alleles[a0 as usize].len() == 1 && alleles[a1 as usize].len() == 1;
    if cfg.phase_only_snvs && !is_snv {
        st.skipped_unphased += 1;
        return Ok(None);
    }
    Ok(Some(PhaseCandidate {
        source_record,
        rid: record.rid().unwrap_or(0) as i32,
        chrom: chrom.to_string(),
        pos: pos1,
        end: pos1 + ref_len - 1,
        alleles,
        gt_alleles: [a0, a1],
        input_order_alleles: [a0, a1],
        is_snv,
    }))
}

fn add_observation_for_candidate(
    cfg: &Config,
    candidate: &PhaseCandidate,
    hap: usize,
    allele: i32,
    ps: i64,
    codon_map: Option<&CodonMap>,
    obs: &mut Vec<Obs>,
    st: &mut Stats,
) -> Result<(), String> {
    if allele == 0 {
        st.skipped_ref_allele += 1;
        return Ok(());
    }
    if allele < 0 || allele as usize >= candidate.alleles.len() {
        st.skipped_unsupported_alt += 1;
        st.skipped_alt_out_of_range += 1;
        if cfg.unsupported_alleles == UnsupportedAllelesPolicy::Fail {
            return Err(format!(
                "unsupported selected ALT allele index at {}:{} hap={} allele={allele}",
                candidate.chrom,
                candidate.pos,
                hap + 1
            ));
        }
        return Ok(());
    }

    let ref_allele = &candidate.alleles[0];
    let alt_allele = &candidate.alleles[allele as usize];
    if !is_plain_dna_allele(alt_allele.as_bytes()) {
        st.skipped_unsupported_alt += 1;
        let kind = unsupported_alt_kind(alt_allele.as_bytes());
        match kind {
            "symbolic_or_breakend" => st.skipped_alt_symbolic_or_breakend += 1,
            "spanning_deletion" => st.skipped_alt_spanning_deletion += 1,
            _ => st.skipped_alt_non_dna += 1,
        }
        if cfg.unsupported_alleles == UnsupportedAllelesPolicy::Fail {
            return Err(format!(
                "unsupported selected ALT allele at {}:{} hap={} ALT={} kind={kind}",
                candidate.chrom,
                candidate.pos,
                hap + 1,
                alt_allele
            ));
        }
        return Ok(());
    }
    if ref_allele.eq_ignore_ascii_case(alt_allele) {
        st.skipped_unsupported_alt += 1;
        st.skipped_alt_same_as_ref += 1;
        if cfg.unsupported_alleles == UnsupportedAllelesPolicy::Fail {
            return Err(format!(
                "unsupported selected ALT allele at {}:{} hap={} ALT equals REF ({alt_allele})",
                candidate.chrom,
                candidate.pos,
                hap + 1
            ));
        }
        return Ok(());
    }

    if contains_n_base(ref_allele.as_bytes()) || contains_n_base(alt_allele.as_bytes()) {
        st.observations_with_n += 1;
        if cfg.warn_on_n {
            eprintln!(
                "warning: N base in selected allele at {}:{} hap={} REF={} ALT={}",
                candidate.chrom,
                candidate.pos,
                hap + 1,
                ref_allele,
                alt_allele
            );
        }
    }

    let is_snv = candidate.is_snv && is_snv_pair(ref_allele, alt_allele);
    let codon_keys = if is_snv {
        codon_map
            .map(|m| m.keys_for(&candidate.chrom, candidate.pos))
            .unwrap_or_default()
    } else {
        Vec::new()
    };
    obs.push(Obs {
        rid: candidate.rid,
        hap: hap as i32,
        ps,
        pos: candidate.pos,
        end: candidate.end,
        source_record: candidate.source_record,
        source_allele: allele,
        ref_allele: ref_allele.clone(),
        alt_allele: alt_allele.clone(),
        is_snv,
        codon_keys,
    });
    st.observations += 1;
    Ok(())
}

fn read_observations_with_bam_phasing(
    cfg: &Config,
) -> Result<(HeaderInfo, Vec<Obs>, Stats), String> {
    preflight_vcf_header(&cfg.input_path)?;
    let mut reader = bcf::Reader::from_path(&cfg.input_path)
        .map_err(|e| format!("cannot open input '{}': {}", cfg.input_path, e))?;
    if cfg.threads > 1 {
        reader
            .set_threads(cfg.threads)
            .map_err(|e| format!("failed to enable VCF/BCF input threads: {e}"))?;
    }
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
    let codon_map = load_optional_codon_map(cfg)?;

    let mut candidates = Vec::new();
    let mut st = Stats::default();
    let mut record = reader.empty_record();

    while let Some(result) = reader.read(&mut record) {
        result.map_err(|e| format!("failed to read VCF/BCF record: {e}"))?;
        st.records += 1;
        let source_record = st.records - 1;
        if record.alleles().len() > 2 {
            st.multiallelic_records += 1;
        }

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
        let Some(a0) = allele_index(gt[0]) else {
            st.skipped_missing_gt += 1;
            continue;
        };
        let Some(a1) = allele_index(gt[1]) else {
            st.skipped_missing_gt += 1;
            continue;
        };

        let alleles_raw = record.alleles();
        if alleles_raw.is_empty() {
            st.skipped_ref += 1;
            if cfg.unsupported_alleles == UnsupportedAllelesPolicy::Fail {
                return Err("unsupported empty REF allele".to_string());
            }
            continue;
        }

        let pos1 = record.pos() + 1;
        let rid = record.rid().unwrap_or(0) as usize;
        let chrom = header_info
            .rid_names
            .get(rid)
            .map(String::as_str)
            .unwrap_or(".");
        if !is_plain_dna_allele(alleles_raw[0]) {
            st.skipped_ref += 1;
            if cfg.unsupported_alleles == UnsupportedAllelesPolicy::Fail {
                let ref_allele = uppercase_ascii_string(alleles_raw[0]);
                return Err(format!(
                    "unsupported REF allele at {chrom}:{pos1} REF={ref_allele}"
                ));
            }
            continue;
        }

        for allele in [a0, a1] {
            if allele < 0 || allele as usize >= alleles_raw.len() {
                st.skipped_unsupported_alt += 1;
                st.skipped_alt_out_of_range += 1;
                if cfg.unsupported_alleles == UnsupportedAllelesPolicy::Fail {
                    return Err(format!(
                        "unsupported selected allele index at {chrom}:{pos1} allele={allele}"
                    ));
                }
                continue;
            }
            if allele != 0 && !is_plain_dna_allele(alleles_raw[allele as usize]) {
                st.skipped_unsupported_alt += 1;
                let kind = unsupported_alt_kind(alleles_raw[allele as usize]);
                match kind {
                    "symbolic_or_breakend" => st.skipped_alt_symbolic_or_breakend += 1,
                    "spanning_deletion" => st.skipped_alt_spanning_deletion += 1,
                    _ => st.skipped_alt_non_dna += 1,
                }
                if cfg.unsupported_alleles == UnsupportedAllelesPolicy::Fail {
                    let alt_allele = uppercase_ascii_string(alleles_raw[allele as usize]);
                    return Err(format!(
                        "unsupported selected allele at {chrom}:{pos1} ALT={alt_allele} kind={kind}"
                    ));
                }
                continue;
            }
        }

        if a0 == a1 {
            st.skipped_unphased += 1;
            continue;
        }
        if a0 < 0 || a1 < 0 || a0 as usize >= alleles_raw.len() || a1 as usize >= alleles_raw.len()
        {
            continue;
        }
        if !is_plain_dna_allele(alleles_raw[a0 as usize])
            || !is_plain_dna_allele(alleles_raw[a1 as usize])
        {
            continue;
        }

        let alleles = alleles_raw
            .iter()
            .map(|a| uppercase_ascii_string(a))
            .collect::<Vec<_>>();
        let ref_len = alleles[0].len() as i64;
        let is_snv = alleles[a0 as usize].len() == 1 && alleles[a1 as usize].len() == 1;
        candidates.push(PhaseCandidate {
            source_record,
            rid: record.rid().unwrap_or(0) as i32,
            chrom: chrom.to_string(),
            pos: pos1,
            end: pos1 + ref_len - 1,
            alleles,
            gt_alleles: [a0, a1],
            input_order_alleles: [a0, a1],
            is_snv,
        });
    }

    st.bam_phase_candidates = candidates.len() as u64;
    let assignments = phase_candidates_from_bam(cfg, &candidates, &header_info.sample, &mut st)?;

    let mut obs = Vec::new();
    for (idx, candidate) in candidates.iter().enumerate() {
        let Some(assignment) = assignments.get(&idx).copied() else {
            continue;
        };
        st.phased_records += 1;
        let hap_alleles = if assignment.rel == 0 {
            candidate.gt_alleles
        } else {
            [candidate.gt_alleles[1], candidate.gt_alleles[0]]
        };
        add_observation_for_candidate(
            cfg,
            candidate,
            0,
            hap_alleles[0],
            assignment.ps,
            codon_map.as_ref(),
            &mut obs,
            &mut st,
        )?;
        add_observation_for_candidate(
            cfg,
            candidate,
            1,
            hap_alleles[1],
            assignment.ps,
            codon_map.as_ref(),
            &mut obs,
            &mut st,
        )?;
    }

    Ok((header_info, obs, st))
}

fn read_observations(cfg: &Config) -> Result<(HeaderInfo, Vec<Obs>, Stats), String> {
    if cfg.phase_bam_path.is_some() {
        return read_observations_with_bam_phasing(cfg);
    }
    preflight_vcf_header(&cfg.input_path)?;
    let mut reader = bcf::Reader::from_path(&cfg.input_path)
        .map_err(|e| format!("cannot open input '{}': {}", cfg.input_path, e))?;
    if cfg.threads > 1 {
        reader
            .set_threads(cfg.threads)
            .map_err(|e| format!("failed to enable VCF/BCF input threads: {e}"))?;
    }
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
    let codon_map = load_optional_codon_map(cfg)?;

    let mut obs = Vec::new();
    let mut st = Stats::default();
    let mut record = reader.empty_record();

    while let Some(result) = reader.read(&mut record) {
        result.map_err(|e| format!("failed to read VCF/BCF record: {e}"))?;
        st.records += 1;
        let source_record = st.records - 1;
        if record.alleles().len() > 2 {
            st.multiallelic_records += 1;
        }

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
            if cfg.unsupported_alleles == UnsupportedAllelesPolicy::Fail {
                return Err("unsupported empty REF allele".to_string());
            }
            continue;
        }
        let ref_bytes = alleles[0];
        let pos1 = record.pos() + 1;
        let rid = record.rid().unwrap_or(0) as usize;
        let chrom = header_info
            .rid_names
            .get(rid)
            .map(String::as_str)
            .unwrap_or(".");
        if !is_plain_dna_allele(ref_bytes) {
            st.skipped_ref += 1;
            if cfg.unsupported_alleles == UnsupportedAllelesPolicy::Fail {
                let ref_allele = uppercase_ascii_string(ref_bytes);
                return Err(format!(
                    "unsupported REF allele at {chrom}:{pos1} REF={ref_allele}"
                ));
            }
            continue;
        }

        let ps = get_sample_ps(&record, sample_idx);
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
                st.skipped_alt_out_of_range += 1;
                if cfg.unsupported_alleles == UnsupportedAllelesPolicy::Fail {
                    return Err(format!(
                        "unsupported selected ALT allele index at {chrom}:{pos1} hap={} allele={allele}",
                        hap + 1
                    ));
                }
                continue;
            }
            let alt_bytes = alleles[allele as usize];
            if !is_plain_dna_allele(alt_bytes) {
                st.skipped_unsupported_alt += 1;
                let kind = unsupported_alt_kind(alt_bytes);
                match kind {
                    "symbolic_or_breakend" => st.skipped_alt_symbolic_or_breakend += 1,
                    "spanning_deletion" => st.skipped_alt_spanning_deletion += 1,
                    _ => st.skipped_alt_non_dna += 1,
                }
                if cfg.unsupported_alleles == UnsupportedAllelesPolicy::Fail {
                    let alt_allele = uppercase_ascii_string(alt_bytes);
                    return Err(format!(
                        "unsupported selected ALT allele at {chrom}:{pos1} hap={} ALT={alt_allele} kind={kind}",
                        hap + 1
                    ));
                }
                continue;
            }
            let alt_allele = uppercase_ascii_string(alt_bytes);
            if ref_allele.eq_ignore_ascii_case(&alt_allele) {
                st.skipped_unsupported_alt += 1;
                st.skipped_alt_same_as_ref += 1;
                if cfg.unsupported_alleles == UnsupportedAllelesPolicy::Fail {
                    return Err(format!(
                        "unsupported selected ALT allele at {chrom}:{pos1} hap={} ALT equals REF ({alt_allele})",
                        hap + 1
                    ));
                }
                continue;
            }
            if contains_n_base(ref_bytes) || contains_n_base(alt_bytes) {
                st.observations_with_n += 1;
                if cfg.warn_on_n {
                    eprintln!(
                        "warning: N base in selected allele at {chrom}:{pos1} hap={} REF={} ALT={}",
                        hap + 1,
                        ref_allele,
                        alt_allele
                    );
                }
            }
            let is_snv = is_snv_pair(&ref_allele, &alt_allele);
            let codon_keys = if is_snv {
                codon_map
                    .as_ref()
                    .map(|m| m.keys_for(chrom, pos1))
                    .unwrap_or_default()
            } else {
                Vec::new()
            };
            obs.push(Obs {
                rid: record.rid().unwrap_or(0) as i32,
                hap: hap as i32,
                ps,
                pos: pos1,
                end: pos1 + ref_len - 1,
                source_record,
                source_allele: allele,
                ref_allele: ref_allele.clone(),
                alt_allele,
                is_snv,
                codon_keys,
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

fn source_haplotypes_from_obs(obs: &[Obs]) -> Vec<SourceHaplotype> {
    obs.iter()
        .map(|o| SourceHaplotype {
            source_record: o.source_record,
            rid: o.rid,
            pos: o.pos,
            hap: o.hap,
            allele: o.source_allele,
        })
        .collect()
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
        source_haplotypes: source_haplotypes_from_obs(block),
    });
    Ok(())
}

fn build_calls_proximity(
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

fn build_calls_nirvana_codon(
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

    let mut grouped: HashMap<(i32, i32, i64, String), Vec<Obs>> = HashMap::new();
    for o in obs.iter() {
        if !o.is_snv || o.codon_keys.is_empty() {
            continue;
        }
        for key in &o.codon_keys {
            grouped
                .entry((o.rid, o.hap, o.ps, key.clone()))
                .or_default()
                .push(o.clone());
        }
    }

    let mut blocks = grouped.into_values().collect::<Vec<_>>();
    blocks.sort_by(|a, b| {
        let aa = &a[0];
        let bb = &b[0];
        aa.rid
            .cmp(&bb.rid)
            .then_with(|| aa.hap.cmp(&bb.hap))
            .then_with(|| aa.ps.cmp(&bb.ps))
            .then_with(|| aa.pos.cmp(&bb.pos))
    });
    for mut block in blocks {
        block.sort_by(cmp_obs);
        block.dedup_by_key(|o| (o.rid, o.hap, o.ps, o.pos, o.end, o.alt_allele.clone()));
        if block.len() >= cfg.min_variants {
            add_call_from_block(cfg, header, fai, &block, &mut calls)?;
        }
    }
    Ok(calls)
}

fn build_calls(
    cfg: &Config,
    header: &HeaderInfo,
    fai: *mut htslib::faidx_t,
    obs: &mut [Obs],
) -> Result<Vec<MnvCall>, String> {
    match cfg.mnv_algorithm {
        MnvAlgorithm::Proximity => build_calls_proximity(cfg, header, fai, obs),
        MnvAlgorithm::NirvanaCodon => build_calls_nirvana_codon(cfg, header, fai, obs),
    }
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
                last.source_haplotypes.extend(call.source_haplotypes);
                last.source_haplotypes
                    .sort_by_key(|s| (s.source_record, s.rid, s.pos, s.hap, s.allele));
                last.source_haplotypes.dedup();
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

fn hp_for_mask(mask: i32, ps: i64) -> String {
    if ps == PS_MISSING {
        return ".".to_string();
    }
    let ps_label = ps;
    match mask & 3 {
        1 => format!("{ps_label}-2,{ps_label}-1"),
        2 => format!("{ps_label}-1,{ps_label}-2"),
        3 => format!("{ps_label}-2,{ps_label}-2"),
        _ => ".".to_string(),
    }
}

fn hp_for_assignment(assignment: PhaseAssignment) -> String {
    if assignment.rel == 0 {
        format!("{}-1,{}-2", assignment.ps, assignment.ps)
    } else {
        format!("{}-2,{}-1", assignment.ps, assignment.ps)
    }
}

fn push_output_header_records(h: &mut bcf::Header, cfg: &Config, header: &HeaderInfo) {
    h.push_record(b"##fileformat=VCFv4.3");
    h.push_record(b"##source=phase_mnv");
    h.push_record(b"##FILTER=<ID=PASS,Description=\"All filters passed\">");
    h.push_record(b"##phase_mnv_normalization=Tan2015_left_aligned_parsimonious");
    h.push_record(b"##phase_mnv_normalization_citation=Tan_A_Abecasis_GR_Kang_HM_Bioinformatics_2015_31_13_2202_2204_doi_10.1093/bioinformatics/btv112");
    h.push_record(format!("##phase_mnv_input={}", cfg.input_path).as_bytes());
    h.push_record(format!("##phase_mnv_mnv_algorithm={}", cfg.mnv_algorithm.as_str()).as_bytes());
    if let Some(path) = cfg.codon_map_path.as_deref() {
        h.push_record(format!("##phase_mnv_codon_map={path}").as_bytes());
    }
    if let Some(bam_path) = cfg.phase_bam_path.as_deref() {
        h.push_record(format!("##phase_mnv_phase_from_bam={bam_path}").as_bytes());
        h.push_record(format!("##phase_mnv_phase_model={}", phase_model_name(cfg)).as_bytes());
        h.push_record(format!("##phase_mnv_phase_tag={}", cfg.phase_tag.as_str()).as_bytes());
    }
    h.push_record(format!("##reference={}", cfg.fasta_path).as_bytes());
    h.push_record(b"##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Merged call type: MNV for pure SNV blocks, COMPLEX when indels are included\">");
    h.push_record(b"##INFO=<ID=NVAR,Number=1,Type=Integer,Description=\"Number of phased source variants merged into this call\">");
    h.push_record(b"##INFO=<ID=NSNPS,Number=1,Type=Integer,Description=\"Number of source SNVs in this merged call\">");
    h.push_record(b"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate of merged reference span\">");
    h.push_record(b"##INFO=<ID=SOURCE_POS,Number=.,Type=Integer,Description=\"Original source variant positions merged into this call\">");
    h.push_record(b"##INFO=<ID=HAPS,Number=.,Type=Integer,Description=\"One-based phased haplotypes carrying this merged call\">");
    h.push_record(b"##INFO=<ID=PS,Number=1,Type=Integer,Description=\"Phase set shared by merged variants, when present in input FORMAT/PS\">");
    h.push_record(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotype for the constructed call in the selected sample\">");
    match cfg.phase_tag {
        PhaseTag::Ps => h.push_record(b"##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set for the constructed call, or missing if absent/ambiguous\">"),
        PhaseTag::Hp => h.push_record(b"##FORMAT=<ID=HP,Number=.,Type=String,Description=\"Phasing haplotype identifier assigned by phase_mnv_rs\">"),
    };
    for contig in &header.contigs {
        h.push_record(format!("##contig=<ID={contig}>").as_bytes());
    }
    h.push_sample(header.sample.as_bytes());
}

fn make_bcf_header(cfg: &Config, header: &HeaderInfo) -> bcf::Header {
    let mut h = bcf::Header::new();
    push_output_header_records(&mut h, cfg, header);
    h
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
    writeln!(
        out,
        "##phase_mnv_mnv_algorithm={}",
        cfg.mnv_algorithm.as_str()
    )?;
    if let Some(path) = cfg.codon_map_path.as_deref() {
        writeln!(out, "##phase_mnv_codon_map={path}")?;
    }
    if let Some(bam_path) = cfg.phase_bam_path.as_deref() {
        writeln!(out, "##phase_mnv_phase_from_bam={bam_path}")?;
        writeln!(out, "##phase_mnv_phase_model={}", phase_model_name(cfg))?;
        writeln!(out, "##phase_mnv_phase_tag={}", cfg.phase_tag.as_str())?;
    }
    writeln!(out, "##reference={}", cfg.fasta_path)?;
    writeln!(out, "##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Merged call type: MNV for pure SNV blocks, COMPLEX when indels are included\">")?;
    writeln!(out, "##INFO=<ID=NVAR,Number=1,Type=Integer,Description=\"Number of phased source variants merged into this call\">")?;
    writeln!(out, "##INFO=<ID=NSNPS,Number=1,Type=Integer,Description=\"Number of source SNVs in this merged call\">")?;
    writeln!(out, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate of merged reference span\">")?;
    writeln!(out, "##INFO=<ID=SOURCE_POS,Number=.,Type=Integer,Description=\"Original source variant positions merged into this call\">")?;
    writeln!(out, "##INFO=<ID=HAPS,Number=.,Type=Integer,Description=\"One-based phased haplotypes carrying this merged call\">")?;
    writeln!(out, "##INFO=<ID=PS,Number=1,Type=Integer,Description=\"Phase set shared by merged variants, when present in input FORMAT/PS\">")?;
    writeln!(out, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotype for the constructed call in the selected sample\">")?;
    match cfg.phase_tag {
        PhaseTag::Ps => writeln!(out, "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set for the constructed call, or missing if absent/ambiguous\">")?,
        PhaseTag::Hp => writeln!(out, "##FORMAT=<ID=HP,Number=.,Type=String,Description=\"Phasing haplotype identifier assigned by phase_mnv_rs\">")?,
    }
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
    cfg: &Config,
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
        match cfg.phase_tag {
            PhaseTag::Ps => {
                let ps = if c.ps == PS_MISSING {
                    ".".to_string()
                } else {
                    c.ps.to_string()
                };
                writeln!(out, "\tGT:PS\t{}:{}", gt, ps)?;
            }
            PhaseTag::Hp => {
                let hp = hp_for_mask(c.hap_mask, c.ps);
                writeln!(out, "\tGT:HP\t{}:{}", gt, hp)?;
            }
        }
        st.emitted += 1;
    }
    Ok(())
}

fn parse_int_list(s: &str) -> Result<Vec<i32>, String> {
    if s.is_empty() || s == "." {
        return Ok(Vec::new());
    }
    s.split(',')
        .map(|x| {
            x.parse::<i32>()
                .map_err(|_| format!("internal error: invalid integer list value '{x}'"))
        })
        .collect()
}

fn missing_bcf_float() -> f32 {
    f32::from_bits(0x7F80_0001)
}

fn genotype_for_mask(mask: i32) -> [GenotypeAllele; 2] {
    match mask & 3 {
        1 => [GenotypeAllele::Unphased(1), GenotypeAllele::Phased(0)],
        2 => [GenotypeAllele::Unphased(0), GenotypeAllele::Phased(1)],
        3 => [GenotypeAllele::Unphased(1), GenotypeAllele::Phased(1)],
        _ => [
            GenotypeAllele::UnphasedMissing,
            GenotypeAllele::PhasedMissing,
        ],
    }
}

fn write_calls_bcf(
    writer: &mut bcf::Writer,
    cfg: &Config,
    calls: &[MnvCall],
    st: &mut Stats,
) -> Result<(), String> {
    for c in calls {
        let mut record = writer.empty_record();
        record.set_rid(Some(c.rid as u32));
        record.set_pos(c.start - 1);
        record
            .set_id(b".")
            .map_err(|e| format!("failed to set BCF ID: {e}"))?;
        record.set_qual(missing_bcf_float());
        record
            .push_filter(&b"PASS"[..])
            .map_err(|e| format!("failed to set BCF FILTER: {e}"))?;
        record
            .set_alleles(&[c.ref_seq.as_bytes(), c.alt_seq.as_bytes()])
            .map_err(|e| format!("failed to set BCF alleles: {e}"))?;
        record
            .push_info_string(b"TYPE", &[c.call_type.as_bytes()])
            .map_err(|e| format!("failed to set INFO/TYPE: {e}"))?;
        record
            .push_info_integer(b"NVAR", &[c.nvars as i32])
            .map_err(|e| format!("failed to set INFO/NVAR: {e}"))?;
        record
            .push_info_integer(b"NSNPS", &[c.nsnps as i32])
            .map_err(|e| format!("failed to set INFO/NSNPS: {e}"))?;
        record
            .push_info_integer(b"END", &[c.end as i32])
            .map_err(|e| format!("failed to set INFO/END: {e}"))?;
        let source_pos = parse_int_list(&c.positions)?;
        record
            .push_info_integer(b"SOURCE_POS", &source_pos)
            .map_err(|e| format!("failed to set INFO/SOURCE_POS: {e}"))?;
        let haps = parse_int_list(haps_for_mask(c.hap_mask))?;
        record
            .push_info_integer(b"HAPS", &haps)
            .map_err(|e| format!("failed to set INFO/HAPS: {e}"))?;
        if c.ps != PS_MISSING {
            record
                .push_info_integer(b"PS", &[c.ps as i32])
                .map_err(|e| format!("failed to set INFO/PS: {e}"))?;
        }
        record
            .push_genotypes(&genotype_for_mask(c.hap_mask))
            .map_err(|e| format!("failed to set FORMAT/GT: {e}"))?;
        match cfg.phase_tag {
            PhaseTag::Ps => {
                let ps_value = if c.ps == PS_MISSING {
                    BCF_INT32_MISSING
                } else {
                    c.ps as i32
                };
                record
                    .push_format_integer(b"PS", &[ps_value])
                    .map_err(|e| format!("failed to set FORMAT/PS: {e}"))?;
            }
            PhaseTag::Hp => {
                let hp = hp_for_mask(c.hap_mask, c.ps);
                record
                    .push_format_string(b"HP", &[hp.as_bytes()])
                    .map_err(|e| format!("failed to set FORMAT/HP: {e}"))?;
            }
        }
        writer
            .write(&record)
            .map_err(|e| format!("failed to write BCF record: {e}"))?;
        st.emitted += 1;
    }
    Ok(())
}

fn consumed_haplotypes_by_record(calls: &[MnvCall]) -> HashMap<u64, HashSet<(i32, i32)>> {
    let mut consumed: HashMap<u64, HashSet<(i32, i32)>> = HashMap::new();
    for call in calls {
        for source in &call.source_haplotypes {
            consumed
                .entry(source.source_record)
                .or_default()
                .insert((source.hap, source.allele));
        }
    }
    consumed
}

fn call_sorts_before_or_at_record(call: &MnvCall, record: &bcf::Record) -> bool {
    let record_rid = record.rid().map(|rid| rid as i32).unwrap_or(i32::MAX);
    let record_pos1 = record.pos() + 1;
    call.rid < record_rid || (call.rid == record_rid && call.start <= record_pos1)
}

fn genotype_has_alt(gt: &[GenotypeAllele]) -> bool {
    gt.iter()
        .filter_map(|&allele| allele_index(allele))
        .any(|allele| allele > 0)
}

fn clear_info_if_present(record: &mut bcf::Record, tag: &[u8]) -> Result<(), String> {
    let Ok((tag_type, _)) = record.header().info_type(tag) else {
        return Ok(());
    };
    match tag_type {
        TagType::Integer => record.clear_info_integer(tag),
        TagType::Float => record.clear_info_float(tag),
        TagType::String => record.clear_info_string(tag),
        TagType::Flag => record.clear_info_flag(tag),
    }
    .map_err(|e| {
        format!(
            "failed to clear INFO/{} after sample subsetting: {e}",
            String::from_utf8_lossy(tag)
        )
    })
}

fn clear_cohort_info_after_sample_subset(record: &mut bcf::Record) -> Result<(), String> {
    for tag in [
        b"AC".as_slice(),
        b"AN".as_slice(),
        b"AF".as_slice(),
        b"MLEAC".as_slice(),
        b"MLEAF".as_slice(),
        b"NS".as_slice(),
    ] {
        clear_info_if_present(record, tag)?;
    }
    Ok(())
}

fn write_one_input_record(
    writer: &mut bcf::Writer,
    mut record: bcf::Record,
    sample_idx: usize,
    clear_cohort_info: bool,
    st: &mut Stats,
) -> Result<(), String> {
    let selected_gt = record
        .genotypes()
        .ok()
        .map(|genotypes| genotypes.get(sample_idx).to_vec());
    if selected_gt
        .as_deref()
        .is_some_and(|gt| !genotype_has_alt(gt))
    {
        return Ok(());
    }

    writer.translate(&mut record);
    writer.subset(&mut record);
    if clear_cohort_info {
        clear_cohort_info_after_sample_subset(&mut record)?;
    }
    writer
        .write(&record)
        .map_err(|e| format!("failed to write unmerged input record: {e}"))?;
    st.emitted += 1;
    Ok(())
}

fn residual_genotype_allele(new_allele: Option<i32>, force_phased: bool) -> GenotypeAllele {
    match new_allele {
        Some(allele) if force_phased => GenotypeAllele::Phased(allele),
        Some(allele) => GenotypeAllele::Unphased(allele),
        None if force_phased => GenotypeAllele::PhasedMissing,
        None => GenotypeAllele::UnphasedMissing,
    }
}

fn write_residual_input_record(
    writer: &mut bcf::Writer,
    mut record: bcf::Record,
    sample_idx: usize,
    consumed_site: &HashSet<(i32, i32)>,
    st: &mut Stats,
) -> Result<(), String> {
    let allele_count = record.alleles().len();
    if allele_count <= 1 {
        return Ok(());
    }

    let original_gt = {
        let genotypes = record
            .genotypes()
            .map_err(|e| format!("failed to read GT for combined residual record: {e}"))?;
        genotypes.get(sample_idx).to_vec()
    };

    let mut consumed_alleles = HashSet::new();
    let mut residual_sample_alleles = HashSet::new();
    for (hap, gt_allele) in original_gt.iter().enumerate() {
        let Some(allele) = allele_index(*gt_allele) else {
            continue;
        };
        if allele <= 0 {
            continue;
        }
        if consumed_site.contains(&(hap as i32, allele)) {
            consumed_alleles.insert(allele);
        } else {
            residual_sample_alleles.insert(allele);
        }
    }
    for &(_, allele) in consumed_site {
        if allele > 0 {
            consumed_alleles.insert(allele);
        }
    }

    let mut remove = vec![false; allele_count];
    for allele in 1..allele_count {
        let allele_i = allele as i32;
        remove[allele] =
            consumed_alleles.contains(&allele_i) && !residual_sample_alleles.contains(&allele_i);
    }
    if remove.iter().skip(1).all(|&remove_allele| remove_allele) {
        return Ok(());
    }

    let mut old_to_new = vec![None; allele_count];
    old_to_new[0] = Some(0);
    let mut next_allele = 1i32;
    for allele in 1..allele_count {
        if !remove[allele] {
            old_to_new[allele] = Some(next_allele);
            next_allele += 1;
        }
    }

    let mut residual_gt = Vec::with_capacity(original_gt.len());
    for (hap, gt_allele) in original_gt.iter().enumerate() {
        let new_allele = match allele_index(*gt_allele) {
            Some(allele) if allele > 0 && consumed_site.contains(&(hap as i32, allele)) => Some(0),
            Some(allele) if allele >= 0 => old_to_new.get(allele as usize).copied().flatten(),
            Some(_) | None => None,
        };
        let force_phased = hap > 0 && is_phased_after_first(*gt_allele);
        residual_gt.push(residual_genotype_allele(new_allele, force_phased));
    }

    if !genotype_has_alt(&residual_gt) {
        return Ok(());
    }

    writer.translate(&mut record);
    writer.subset(&mut record);
    if remove.iter().any(|&remove_allele| remove_allele) {
        record
            .remove_alleles(&remove)
            .map_err(|e| format!("failed to remove consumed alleles from residual record: {e}"))?;
    }
    record
        .push_genotypes(&residual_gt)
        .map_err(|e| format!("failed to update residual combined GT: {e}"))?;
    clear_cohort_info_after_sample_subset(&mut record)?;
    writer
        .write(&record)
        .map_err(|e| format!("failed to write residual input record: {e}"))?;
    st.emitted += 1;
    Ok(())
}

fn validate_combined_request(cfg: &Config) -> Result<(), String> {
    if cfg.emit_mode != EmitMode::Combined {
        return Ok(());
    }
    if cfg.no_header {
        return Err(
            "--emit combined preserves the original VCF/BCF header; --no-header is not supported"
                .to_string(),
        );
    }
    if cfg.phase_bam_path.is_some() {
        return Err("--emit combined currently requires input-phased VCF/BCF; --phase-from-bam is supported with --emit mnv or --emit all-sites".to_string());
    }
    Ok(())
}

fn run_combined(
    cfg: &Config,
    header: &HeaderInfo,
    calls: &[MnvCall],
    st: &mut Stats,
) -> Result<(), String> {
    validate_combined_request(cfg)?;

    let consumed = consumed_haplotypes_by_record(calls);
    let mut reader = bcf::Reader::from_path(&cfg.input_path)
        .map_err(|e| format!("cannot open input '{}': {}", cfg.input_path, e))?;
    if cfg.threads > 1 {
        reader
            .set_threads(cfg.threads)
            .map_err(|e| format!("failed to enable VCF/BCF input threads: {e}"))?;
    }
    let clear_cohort_info = reader.header().sample_count() > 1;
    let sample_idx = reader
        .header()
        .sample_id(header.sample.as_bytes())
        .ok_or_else(|| {
            format!(
                "sample '{}' not found while writing combined output",
                header.sample
            )
        })?;
    let out_header = make_combined_header(cfg, reader.header(), &header.sample)?;
    let mut writer = open_htslib_writer(cfg, &out_header)?;
    let mut call_idx = 0usize;
    let mut record_ordinal = 0u64;
    let mut record = reader.empty_record();
    while let Some(result) = reader.read(&mut record) {
        result.map_err(|e| format!("failed to read VCF/BCF record: {e}"))?;
        let source_record = record_ordinal;
        record_ordinal += 1;
        while call_idx < calls.len() && call_sorts_before_or_at_record(&calls[call_idx], &record) {
            write_calls_bcf(&mut writer, cfg, &calls[call_idx..call_idx + 1], st)?;
            call_idx += 1;
        }
        if let Some(consumed_site) = consumed.get(&source_record) {
            write_residual_input_record(
                &mut writer,
                record.clone(),
                sample_idx,
                consumed_site,
                st,
            )?;
        } else {
            write_one_input_record(
                &mut writer,
                record.clone(),
                sample_idx,
                clear_cohort_info,
                st,
            )?;
        }
    }
    while call_idx < calls.len() {
        write_calls_bcf(&mut writer, cfg, &calls[call_idx..call_idx + 1], st)?;
        call_idx += 1;
    }
    drop(writer);
    index_output_if_requested(cfg)?;
    Ok(())
}

fn push_all_sites_header_records(h: &mut bcf::Header, cfg: &Config, input_header: &HeaderView) {
    h.push_record(b"##phase_mnv_emit_mode=all-sites");
    h.push_record(b"##phase_mnv_header_policy=preserve_input_header_and_append_phase_mnv_records");
    h.push_record(format!("##phase_mnv_input={}", cfg.input_path).as_bytes());
    h.push_record(format!("##phase_mnv_reference={}", cfg.fasta_path).as_bytes());
    if let Some(bam_path) = cfg.phase_bam_path.as_deref() {
        h.push_record(format!("##phase_mnv_phase_from_bam={bam_path}").as_bytes());
        h.push_record(format!("##phase_mnv_phase_model={}", phase_model_name(cfg)).as_bytes());
        h.push_record(format!("##phase_mnv_phase_tag={}", cfg.phase_tag.as_str()).as_bytes());
    }
    match cfg.phase_tag {
        PhaseTag::Ps => {
            if input_header.format_type(b"PS").is_err() {
                h.push_record(b"##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set assigned by phase_mnv_rs\">");
            }
        }
        PhaseTag::Hp => {
            if input_header.format_type(b"HP").is_err() {
                h.push_record(b"##FORMAT=<ID=HP,Number=.,Type=String,Description=\"Phasing haplotype identifier assigned by phase_mnv_rs\">");
            }
        }
    }
}

fn make_all_sites_header(cfg: &Config, input_header: &HeaderView) -> bcf::Header {
    let mut h = bcf::Header::from_template(input_header);
    push_all_sites_header_records(&mut h, cfg, input_header);
    h
}

fn header_has_info_record(input_header: &HeaderView, tag: &str) -> bool {
    input_header
        .header_records()
        .iter()
        .any(|record| match record {
            HeaderRecord::Info { values, .. } => values.get("ID").is_some_and(|id| id == tag),
            _ => false,
        })
}

fn header_has_format_record(input_header: &HeaderView, tag: &str) -> bool {
    input_header
        .header_records()
        .iter()
        .any(|record| match record {
            HeaderRecord::Format { values, .. } => values.get("ID").is_some_and(|id| id == tag),
            _ => false,
        })
}

fn header_has_filter_record(input_header: &HeaderView, tag: &str) -> bool {
    input_header
        .header_records()
        .iter()
        .any(|record| match record {
            HeaderRecord::Filter { values, .. } => values.get("ID").is_some_and(|id| id == tag),
            _ => false,
        })
}

fn push_info_record_if_absent(
    h: &mut bcf::Header,
    input_header: &HeaderView,
    tag: &str,
    record: &[u8],
) {
    if !header_has_info_record(input_header, tag) {
        h.push_record(record);
    }
}

fn push_format_record_if_absent(
    h: &mut bcf::Header,
    input_header: &HeaderView,
    tag: &str,
    record: &[u8],
) {
    if !header_has_format_record(input_header, tag) {
        h.push_record(record);
    }
}

fn push_filter_record_if_absent(
    h: &mut bcf::Header,
    input_header: &HeaderView,
    tag: &str,
    record: &[u8],
) {
    if !header_has_filter_record(input_header, tag) {
        h.push_record(record);
    }
}

fn push_combined_header_records(h: &mut bcf::Header, cfg: &Config, input_header: &HeaderView) {
    h.push_record(b"##phase_mnv_emit_mode=combined");
    h.push_record(b"##phase_mnv_header_policy=preserve_input_header_subset_selected_sample_and_replace_consumed_records");
    h.push_record(format!("##phase_mnv_input={}", cfg.input_path).as_bytes());
    h.push_record(format!("##phase_mnv_reference={}", cfg.fasta_path).as_bytes());
    h.push_record(format!("##phase_mnv_mnv_algorithm={}", cfg.mnv_algorithm.as_str()).as_bytes());
    if let Some(path) = cfg.codon_map_path.as_deref() {
        h.push_record(format!("##phase_mnv_codon_map={path}").as_bytes());
    }
    h.push_record(b"##phase_mnv_normalization=Tan2015_left_aligned_parsimonious");
    h.push_record(b"##phase_mnv_normalization_citation=Tan_A_Abecasis_GR_Kang_HM_Bioinformatics_2015_31_13_2202_2204_doi_10.1093/bioinformatics/btv112");
    push_filter_record_if_absent(
        h,
        input_header,
        "PASS",
        b"##FILTER=<ID=PASS,Description=\"All filters passed\">",
    );
    push_info_record_if_absent(h, input_header, "TYPE", b"##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Merged call type: MNV for pure SNV blocks, COMPLEX when indels are included\">");
    push_info_record_if_absent(h, input_header, "NVAR", b"##INFO=<ID=NVAR,Number=1,Type=Integer,Description=\"Number of phased source variants merged into this call\">");
    push_info_record_if_absent(h, input_header, "NSNPS", b"##INFO=<ID=NSNPS,Number=1,Type=Integer,Description=\"Number of source SNVs in this merged call\">");
    push_info_record_if_absent(h, input_header, "END", b"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate of merged reference span\">");
    push_info_record_if_absent(h, input_header, "SOURCE_POS", b"##INFO=<ID=SOURCE_POS,Number=.,Type=Integer,Description=\"Original source variant positions merged into this call\">");
    push_info_record_if_absent(h, input_header, "HAPS", b"##INFO=<ID=HAPS,Number=.,Type=Integer,Description=\"One-based phased haplotypes carrying this merged call\">");
    push_info_record_if_absent(h, input_header, "PS", b"##INFO=<ID=PS,Number=1,Type=Integer,Description=\"Phase set shared by merged variants, when present in input FORMAT/PS\">");
    push_format_record_if_absent(h, input_header, "GT", b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotype for the constructed call in the selected sample\">");
    match cfg.phase_tag {
        PhaseTag::Ps => push_format_record_if_absent(
            h,
            input_header,
            "PS",
            b"##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set for the constructed call, or missing if absent/ambiguous\">",
        ),
        PhaseTag::Hp => push_format_record_if_absent(
            h,
            input_header,
            "HP",
            b"##FORMAT=<ID=HP,Number=.,Type=String,Description=\"Phasing haplotype identifier assigned by phase_mnv_rs\">",
        ),
    };
}

fn make_combined_header(
    cfg: &Config,
    input_header: &HeaderView,
    sample: &str,
) -> Result<bcf::Header, String> {
    let mut h =
        bcf::Header::from_template_subset(input_header, &[sample.as_bytes()]).map_err(|e| {
            format!("failed to subset combined output header to sample '{sample}': {e}")
        })?;
    push_combined_header_records(&mut h, cfg, input_header);
    Ok(h)
}

fn open_htslib_writer(cfg: &Config, header: &bcf::Header) -> Result<bcf::Writer, String> {
    let mut writer = match (cfg.output_path.as_deref(), infer_output_kind(cfg)) {
        (None | Some("-"), OutputKind::PlainVcf) => {
            bcf::Writer::from_stdout(header, true, bcf::Format::Vcf)
                .map_err(|e| format!("cannot open VCF stdout: {e}"))?
        }
        (Some(path), OutputKind::PlainVcf) => {
            bcf::Writer::from_path(path, header, true, bcf::Format::Vcf)
                .map_err(|e| format!("cannot open VCF output '{}': {}", path, e))?
        }
        (Some(path), OutputKind::BgzfVcf) => {
            bcf::Writer::from_path(path, header, false, bcf::Format::Vcf)
                .map_err(|e| format!("cannot open BGZF VCF output '{}': {}", path, e))?
        }
        (Some(path), OutputKind::Bcf) => {
            bcf::Writer::from_path(path, header, false, bcf::Format::Bcf)
                .map_err(|e| format!("cannot open BCF output '{}': {}", path, e))?
        }
        (None, OutputKind::BgzfVcf | OutputKind::Bcf) => {
            return Err(
                "compressed VCF/BCF output to stdout requires explicit support; use -o FILE"
                    .to_string(),
            );
        }
    };
    if cfg.threads > 1 {
        writer
            .set_threads(cfg.threads)
            .map_err(|e| format!("failed to enable VCF/BCF output threads: {e}"))?;
    }
    Ok(writer)
}

fn update_record_phase(
    cfg: &Config,
    record: &mut bcf::Record,
    candidate: &PhaseCandidate,
    assignment: PhaseAssignment,
) -> Result<(), String> {
    let hap_alleles = if assignment.rel == 0 {
        candidate.gt_alleles
    } else {
        [candidate.gt_alleles[1], candidate.gt_alleles[0]]
    };
    let gt = [
        GenotypeAllele::Unphased(hap_alleles[0]),
        GenotypeAllele::Phased(hap_alleles[1]),
    ];
    record
        .push_genotypes(&gt)
        .map_err(|e| format!("failed to update FORMAT/GT: {e}"))?;
    match cfg.phase_tag {
        PhaseTag::Ps => record
            .push_format_integer(b"PS", &[assignment.ps as i32])
            .map_err(|e| format!("failed to update FORMAT/PS: {e}"))?,
        PhaseTag::Hp => {
            let hp = hp_for_assignment(assignment);
            record
                .push_format_string(b"HP", &[hp.as_bytes()])
                .map_err(|e| format!("failed to update FORMAT/HP: {e}"))?;
        }
    }
    Ok(())
}

fn run_all_sites(cfg: &Config) -> Result<(), String> {
    if cfg.no_header {
        return Err(
            "--emit all-sites preserves the original VCF/BCF header; --no-header is not supported"
                .to_string(),
        );
    }
    preflight_vcf_header(&cfg.input_path)?;

    let mut planning_reader = bcf::Reader::from_path(&cfg.input_path)
        .map_err(|e| format!("cannot open input '{}': {}", cfg.input_path, e))?;
    if cfg.threads > 1 {
        planning_reader
            .set_threads(cfg.threads)
            .map_err(|e| format!("failed to enable VCF/BCF input threads: {e}"))?;
    }
    let sample_count = planning_reader.header().sample_count();
    let sample_idx = match cfg.sample_name.as_deref() {
        None => 0usize,
        Some(name) => planning_reader
            .header()
            .sample_id(name.as_bytes())
            .ok_or_else(|| {
                let available = planning_reader
                    .header()
                    .samples()
                    .iter()
                    .map(|s| String::from_utf8_lossy(s).into_owned())
                    .collect::<Vec<_>>()
                    .join(" ");
                format!("sample '{name}' not found. Available samples: {available}")
            })?,
    };
    let sample_label = if sample_count == 0 {
        ".".to_string()
    } else {
        String::from_utf8_lossy(planning_reader.header().samples()[sample_idx]).into_owned()
    };
    let out_header = make_all_sites_header(cfg, planning_reader.header());

    let mut st = Stats::default();
    let mut assignments_by_index: HashMap<usize, PhaseAssignment> = HashMap::new();

    if cfg.phase_bam_path.is_some() {
        if sample_count == 0 {
            return Err(
                "--emit all-sites --phase-from-bam requires at least one sample".to_string(),
            );
        }
        if sample_count != 1 {
            return Err("--emit all-sites --phase-from-bam currently updates one-sample VCF/BCF inputs only; use --emit mnv for selected-sample MNV construction".to_string());
        }
        let header_info = collect_header_info(&planning_reader, sample_idx)?;
        let mut candidates = Vec::new();
        let mut record = planning_reader.empty_record();
        while let Some(result) = planning_reader.read(&mut record) {
            result.map_err(|e| format!("failed to read VCF/BCF record: {e}"))?;
            st.records += 1;
            let source_record = st.records - 1;
            if let Some(candidate) = candidate_from_record(
                &record,
                &header_info,
                sample_idx,
                source_record,
                &mut st,
                cfg,
            )? {
                candidates.push(candidate);
            }
        }
        st.bam_phase_candidates = candidates.len() as u64;
        assignments_by_index =
            phase_candidates_from_bam(cfg, &candidates, &header_info.sample, &mut st)?;
    }
    drop(planning_reader);

    let mut reader = bcf::Reader::from_path(&cfg.input_path)
        .map_err(|e| format!("cannot open input '{}': {}", cfg.input_path, e))?;
    if cfg.threads > 1 {
        reader
            .set_threads(cfg.threads)
            .map_err(|e| format!("failed to enable VCF/BCF input threads: {e}"))?;
    }
    let header_info = if cfg.phase_bam_path.is_some() {
        Some(collect_header_info(&reader, sample_idx)?)
    } else {
        None
    };
    let mut writer = open_htslib_writer(cfg, &out_header)?;
    let mut candidate_cursor = 0usize;
    let mut record_ordinal = 0u64;
    let mut record = reader.empty_record();
    while let Some(result) = reader.read(&mut record) {
        result.map_err(|e| format!("failed to read VCF/BCF record: {e}"))?;
        let source_record = record_ordinal;
        record_ordinal += 1;
        if cfg.phase_bam_path.is_none() {
            st.records += 1;
        }
        let maybe_assignment = if let Some(header_info) = header_info.as_ref() {
            let mut dummy = Stats::default();
            if let Some(candidate) = candidate_from_record(
                &record,
                header_info,
                sample_idx,
                source_record,
                &mut dummy,
                cfg,
            )? {
                let idx = candidate_cursor;
                candidate_cursor += 1;
                assignments_by_index
                    .get(&idx)
                    .copied()
                    .map(|assignment| (candidate, assignment))
            } else {
                None
            }
        } else {
            None
        };
        writer.translate(&mut record);
        if let Some((candidate, assignment)) = maybe_assignment {
            update_record_phase(cfg, &mut record, &candidate, assignment)?;
            st.phased_records += 1;
        }
        writer
            .write(&record)
            .map_err(|e| format!("failed to write VCF/BCF record: {e}"))?;
        st.emitted += 1;
    }

    drop(writer);
    index_output_if_requested(cfg)?;
    print_summary(cfg, &st, &sample_label);
    Ok(())
}

fn print_summary(cfg: &Config, st: &Stats, sample: &str) {
    if cfg.quiet {
        return;
    }
    eprintln!(
        "phase_mnv: input={} reference={} output={} sample={}",
        cfg.input_path,
        cfg.fasta_path,
        output_label(cfg),
        sample
    );
    eprintln!(
        "phase_mnv: settings max_gap={} min_vars={} unsupported_alleles={} warn_on_n={} no_ref_check={} no_header={} output_format={} threads={} emit={} mnv_algorithm={} write_index={}",
        cfg.max_gap,
        cfg.min_variants,
        cfg.unsupported_alleles.as_str(),
        cfg.warn_on_n,
        cfg.no_ref_check,
        cfg.no_header,
        infer_output_kind(cfg).as_str(),
        cfg.threads,
        cfg.emit_mode.as_str(),
        cfg.mnv_algorithm.as_str(),
        cfg.output_index
            .map(OutputIndexKind::as_str)
            .unwrap_or("none")
    );
    if let Some(bam_path) = cfg.phase_bam_path.as_deref() {
        eprintln!(
            "phase_mnv: bam_phase input={} algorithm={} tag={} max_coverage={} min_mapq={} min_baseq={} only_snvs={} use_supplementary={} ignore_read_groups={} supplementary_distance={} realign_overhang={} read_list={} candidates={} informative_reads={} selected_reads={} components={} phased_variants={} unphased_variants={} conflicts={}",
            bam_path,
            cfg.phase_algorithm.as_str(),
            cfg.phase_tag.as_str(),
            cfg.phase_max_coverage,
            cfg.phase_min_mapq,
            cfg.phase_min_baseq,
            cfg.phase_only_snvs,
            cfg.phase_use_supplementary,
            cfg.phase_ignore_read_groups,
            cfg.phase_supplementary_distance,
            cfg.phase_realign_overhang,
            cfg.phase_read_list_path.as_deref().unwrap_or("."),
            st.bam_phase_candidates,
            st.bam_phase_informative_reads,
            st.bam_phase_selected_reads,
            st.bam_phase_components,
            st.bam_phase_phased_variants,
            st.bam_phase_unphased_variants,
            st.bam_phase_conflicts
        );
    }
    eprintln!(
        "phase_mnv: records={} phased_records={} haplotype_variant_observations={} emitted_calls={}",
        st.records, st.phased_records, st.observations, st.emitted
    );
    eprintln!(
        "phase_mnv: skipped no_gt={} non_diploid={} missing_gt={} unphased={} ref_hap_alleles={}",
        st.skipped_no_gt,
        st.skipped_not_diploid,
        st.skipped_missing_gt,
        st.skipped_unphased,
        st.skipped_ref_allele
    );
    eprintln!(
        "phase_mnv: unsupported ref_non_dna={} alt_out_of_range={} alt_symbolic_or_breakend={} alt_spanning_deletion={} alt_non_dna={} alt_same_as_ref={} unsupported_alt_total={}",
        st.skipped_ref,
        st.skipped_alt_out_of_range,
        st.skipped_alt_symbolic_or_breakend,
        st.skipped_alt_spanning_deletion,
        st.skipped_alt_non_dna,
        st.skipped_alt_same_as_ref,
        st.skipped_unsupported_alt
    );
    eprintln!(
        "phase_mnv: multiallelic_records={} observations_with_n={}",
        st.multiallelic_records, st.observations_with_n
    );
}

fn run() -> Result<(), String> {
    let cfg = parse_args();
    validate_index_request(&cfg)?;
    validate_combined_request(&cfg)?;
    if cfg.emit_mode == EmitMode::AllSites {
        return run_all_sites(&cfg);
    }
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

    if cfg.emit_mode == EmitMode::Combined {
        run_combined(&cfg, &header, &calls, &mut st)?;
        print_summary(&cfg, &st, &header.sample);
        return Ok(());
    }

    match (cfg.output_path.as_deref(), infer_output_kind(&cfg)) {
        (None | Some("-"), OutputKind::PlainVcf) => {
            let stdout = io::stdout();
            let mut out = BufWriter::new(stdout.lock());
            if !cfg.no_header {
                write_header(&mut out, &cfg, &header).map_err(|e| e.to_string())?;
            }
            write_calls(&mut out, &cfg, &header, &calls, &mut st).map_err(|e| e.to_string())?;
            out.flush().map_err(|e| e.to_string())?;
        }
        (Some(path), OutputKind::PlainVcf) => {
            let file = File::create(Path::new(path))
                .map_err(|e| format!("cannot open output '{}': {}", path, e))?;
            let mut out = BufWriter::new(file);
            if !cfg.no_header {
                write_header(&mut out, &cfg, &header).map_err(|e| e.to_string())?;
            }
            write_calls(&mut out, &cfg, &header, &calls, &mut st).map_err(|e| e.to_string())?;
            out.flush().map_err(|e| e.to_string())?;
        }
        (Some(path), OutputKind::BgzfVcf) => {
            let mut out = BgzfWriter::from_path(path, cfg.threads)?;
            if !cfg.no_header {
                write_header(&mut out, &cfg, &header).map_err(|e| e.to_string())?;
            }
            write_calls(&mut out, &cfg, &header, &calls, &mut st).map_err(|e| e.to_string())?;
            out.flush().map_err(|e| e.to_string())?;
        }
        (Some(path), OutputKind::Bcf) => {
            let bcf_header = make_bcf_header(&cfg, &header);
            let mut writer = bcf::Writer::from_path(path, &bcf_header, false, bcf::Format::Bcf)
                .map_err(|e| format!("cannot open BCF output '{}': {}", path, e))?;
            if cfg.threads > 1 {
                writer
                    .set_threads(cfg.threads)
                    .map_err(|e| format!("failed to enable BCF output threads: {e}"))?;
            }
            write_calls_bcf(&mut writer, &cfg, &calls, &mut st)?;
        }
        (None, OutputKind::BgzfVcf | OutputKind::Bcf) => {
            return Err(
                "compressed VCF/BCF output to stdout requires explicit support; use -o FILE"
                    .to_string(),
            );
        }
    }

    index_output_if_requested(&cfg)?;
    print_summary(&cfg, &st, &header.sample);
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        die(&e);
    }
}
