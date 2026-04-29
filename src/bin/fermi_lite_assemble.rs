#[path = "../fermi_lite.rs"]
mod fermi_lite;

use fermi_lite::{assemble_sequences, AssembleOptions};
use std::io::{self, BufRead};

fn usage() -> &'static str {
    "usage: fermi_lite_assemble [options] [--seq SEQ ...]\n\n\
Small fermi-lite FFI smoke/utility binary. With --seq, assembles the supplied\n\
sequences. Without --seq, reads one plain sequence per non-empty stdin line,\n\
ignoring FASTA-style header lines. This is intended for local adjudication\n\
experiments, not as a full fermi-lite CLI replacement.\n\n\
options:\n\
      --seq SEQ              Add one input read/sequence\n\
  -@, --threads N            fermi-lite threads (default: 1)\n\
      --min-asm-ovlp N       minimum assembly overlap (default: 21)\n\
      --min-count N          minimum k-mer count threshold (default: 1)\n\
      --max-count N          maximum k-mer count threshold (default: 1000)\n\
      --ec-k N               error-correction k; negative disables EC (default: -1)\n\
  -h, --help                 Show this help\n"
}

fn die(msg: &str) -> ! {
    eprintln!("error: {msg}");
    std::process::exit(1);
}

fn read_stdin_sequences() -> Vec<String> {
    let mut seqs = Vec::new();
    let stdin = io::stdin();
    for line in stdin.lock().lines() {
        let line = line.unwrap_or_default();
        let s = line.trim();
        if s.is_empty() || s.starts_with('>') {
            continue;
        }
        seqs.push(s.to_string());
    }
    seqs
}

fn main() {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.iter().any(|a| a == "-h" || a == "--help") {
        print!("{}", usage());
        return;
    }

    let mut options = AssembleOptions::default();
    let mut seqs = Vec::<String>::new();
    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "--seq" => {
                i += 1;
                if i >= args.len() {
                    die("--seq requires an argument");
                }
                seqs.push(args[i].clone());
            }
            "-@" | "--threads" => {
                i += 1;
                if i >= args.len() {
                    die("--threads requires an argument");
                }
                options.threads = args[i]
                    .parse()
                    .unwrap_or_else(|_| die("--threads must be an integer"));
            }
            "--min-asm-ovlp" => {
                i += 1;
                if i >= args.len() {
                    die("--min-asm-ovlp requires an argument");
                }
                options.min_asm_overlap = args[i]
                    .parse()
                    .unwrap_or_else(|_| die("--min-asm-ovlp must be an integer"));
            }
            "--min-count" => {
                i += 1;
                if i >= args.len() {
                    die("--min-count requires an argument");
                }
                options.min_count = args[i]
                    .parse()
                    .unwrap_or_else(|_| die("--min-count must be an integer"));
            }
            "--max-count" => {
                i += 1;
                if i >= args.len() {
                    die("--max-count requires an argument");
                }
                options.max_count = args[i]
                    .parse()
                    .unwrap_or_else(|_| die("--max-count must be an integer"));
            }
            "--ec-k" => {
                i += 1;
                if i >= args.len() {
                    die("--ec-k requires an argument");
                }
                options.error_correction_k = args[i]
                    .parse()
                    .unwrap_or_else(|_| die("--ec-k must be an integer"));
            }
            x => die(&format!("unknown option: {x}")),
        }
        i += 1;
    }

    if seqs.is_empty() {
        seqs = read_stdin_sequences();
    }
    if seqs.is_empty() {
        die("no input sequences supplied");
    }

    let unitigs = assemble_sequences(&seqs, &options).unwrap_or_else(|e| die(&e));
    for (idx, unitig) in unitigs.iter().enumerate() {
        println!(
            ">utg{} len={} supporting_reads={}\n{}",
            idx + 1,
            unitig.len,
            unitig.supporting_reads,
            unitig.seq
        );
    }
}
