use libc::c_void;
use rust_htslib::bcf::{self, Format, Header, Read, Writer};
use rust_htslib::htslib;
use std::ffi::CString;
use std::ptr;

fn usage() -> &'static str {
    "usage: unphase_vcf [options] input.vcf|input.vcf.gz|input.bcf|-\n\n\
Write an unphased VCF stream from VCF/VCF.GZ/BCF input. GT separators are\n\
converted from phased to unphased, and FORMAT/PS plus FORMAT/PQ are removed\n\
by default. Other records, alleles, INFO fields, filters, and non-phase FORMAT\n\
values are preserved through htslib/rust-htslib.\n\n\
options:\n\
  -o, --output FILE       Output VCF path; .gz/.bgz writes BGZF (default: stdout)\n\
      --keep-phase-tags   Keep FORMAT/PS and FORMAT/PQ instead of removing them\n\
  -@, --threads N         htslib threads for compressed input/output (default: 1)\n\
  -h, --help              Show this help\n"
}

#[derive(Debug)]
struct Config {
    input: String,
    output: Option<String>,
    keep_phase_tags: bool,
    threads: usize,
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

    let mut output = None;
    let mut keep_phase_tags = false;
    let mut threads = 1usize;
    let mut positional = Vec::new();

    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "-o" | "--output" => {
                i += 1;
                if i >= args.len() {
                    die("--output requires an argument");
                }
                output = Some(args[i].clone());
            }
            "--keep-phase-tags" => keep_phase_tags = true,
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
            "-" => positional.push(args[i].clone()),
            x if x.starts_with('-') => die(&format!("unknown option: {x}")),
            _ => positional.push(args[i].clone()),
        }
        i += 1;
    }

    if positional.len() != 1 {
        die("expected exactly one input VCF/BCF file");
    }

    Config {
        input: positional.remove(0),
        output,
        keep_phase_tags,
        threads,
    }
}

fn make_header(input: &bcf::header::HeaderView, keep_phase_tags: bool) -> Result<Header, String> {
    let header = Header::from_template(input);
    if !keep_phase_tags {
        for tag in ["PS", "PQ"] {
            let c_tag = CString::new(tag).expect("literal has no NUL");
            unsafe {
                htslib::bcf_hdr_remove(header.inner, htslib::BCF_HL_FMT as i32, c_tag.as_ptr());
            }
        }
        let ret = unsafe { htslib::bcf_hdr_sync(header.inner) };
        if ret < 0 {
            return Err(
                "failed to sync output header after removing phase FORMAT tags".to_string(),
            );
        }
    }
    Ok(header)
}

fn unphase_gt(record: &mut bcf::Record) -> Result<(), String> {
    let c_gt = CString::new("GT").expect("literal has no NUL");
    let mut raw: *mut c_void = ptr::null_mut();
    let mut nraw = 0i32;
    let ret = unsafe {
        htslib::bcf_get_format_values(
            record.header().inner,
            record.inner,
            c_gt.as_ptr(),
            &mut raw,
            &mut nraw,
            htslib::BCF_HT_INT as i32,
        )
    };
    if ret == -1 || ret == -3 {
        return Ok(());
    }
    if ret < 0 {
        return Err(format!(
            "failed to read FORMAT/GT values (htslib code {ret})"
        ));
    }
    if raw.is_null() || ret == 0 {
        return Ok(());
    }

    let values = unsafe { std::slice::from_raw_parts_mut(raw as *mut i32, ret as usize) };
    for value in values.iter_mut() {
        if *value != htslib::bcf_int32_vector_end {
            *value &= !1;
        }
    }

    let update_ret = unsafe {
        htslib::bcf_update_format(
            record.header().inner,
            record.inner,
            c_gt.as_ptr(),
            raw as *const c_void,
            ret,
            htslib::BCF_HT_INT as i32,
        )
    };
    unsafe { libc::free(raw) };
    if update_ret < 0 {
        Err("failed to update FORMAT/GT".to_string())
    } else {
        Ok(())
    }
}

fn remove_format(record: &mut bcf::Record, tag: &str) -> Result<(), String> {
    let c_tag = CString::new(tag).expect("literal has no NUL");
    let ret = unsafe {
        htslib::bcf_update_format(
            record.header().inner,
            record.inner,
            c_tag.as_ptr(),
            ptr::null(),
            0,
            htslib::BCF_HT_INT as i32,
        )
    };
    if ret < 0 {
        Err(format!("failed to remove FORMAT/{tag}"))
    } else {
        Ok(())
    }
}

fn open_writer(cfg: &Config, header: &Header) -> Result<Writer, String> {
    match cfg.output.as_deref() {
        Some(path) if path != "-" => {
            let compressed = path.ends_with(".gz") || path.ends_with(".bgz");
            Writer::from_path(path, header, !compressed, Format::Vcf)
                .map_err(|e| format!("cannot open output VCF '{path}': {e}"))
        }
        _ => Writer::from_stdout(header, true, Format::Vcf)
            .map_err(|e| format!("cannot open stdout VCF writer: {e}")),
    }
}

fn run(cfg: &Config) -> Result<(), String> {
    let mut reader = if cfg.input == "-" {
        bcf::Reader::from_stdin().map_err(|e| format!("cannot open VCF/BCF from stdin: {e}"))?
    } else {
        bcf::Reader::from_path(&cfg.input)
            .map_err(|e| format!("cannot open input '{}': {e}", cfg.input))?
    };
    if cfg.threads > 1 {
        reader
            .set_threads(cfg.threads)
            .map_err(|e| format!("failed to enable input threads for '{}': {e}", cfg.input))?;
    }
    let header = make_header(reader.header(), cfg.keep_phase_tags)?;
    let mut writer = open_writer(cfg, &header)?;
    if cfg.threads > 1 {
        writer
            .set_threads(cfg.threads)
            .map_err(|e| format!("failed to enable output threads: {e}"))?;
    }

    for result in reader.records() {
        let mut record = result.map_err(|e| format!("failed to read '{}': {e}", cfg.input))?;
        unphase_gt(&mut record)?;
        if !cfg.keep_phase_tags {
            remove_format(&mut record, "PS")?;
            remove_format(&mut record, "PQ")?;
        }
        writer
            .write(&record)
            .map_err(|e| format!("failed to write VCF record: {e}"))?;
    }
    Ok(())
}

fn main() {
    let cfg = parse_args();
    if let Err(e) = run(&cfg) {
        die(&e);
    }
}
