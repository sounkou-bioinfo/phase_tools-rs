#!/usr/bin/env python3
"""Strip command/path-bearing metadata from VCF headers.

This is intended for local comparison artefacts where headers may otherwise
record private filesystem paths through producer command lines, especially
`##bcftools_*Command=...` and WhatsHap `##commandline=...` records.

The script supports plain VCF and VCF.GZ/VCF.BGZ input/output. It does not
support BCF; convert BCF to VCF first if header sanitisation is needed.
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path
from typing import TextIO

COMMAND_HEADER_RE = re.compile(r"^##[^=]*(?:Command|CommandLine)=")
BCFTOOLS_HEADER_RE = re.compile(r"^##bcftools_")
PATH_PREFIXES = (
    "##commandline=",
    "##phase_mnv_input=",
    "##phase_mnv_phase_from_bam=",
    "##reference=",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Strip command/path-bearing VCF header records."
    )
    parser.add_argument("input", help="Input VCF/VCF.GZ path, or '-' for stdin")
    parser.add_argument("output", help="Output VCF/VCF.GZ path, or '-' for stdout")
    parser.add_argument(
        "--keep-path-records",
        action="store_true",
        help="Only strip command-style records; keep phase_mnv/reference path records.",
    )
    return parser.parse_args()


def open_text(path: str, mode: str) -> TextIO:
    if path == "-":
        return sys.stdin if "r" in mode else sys.stdout
    if path.lower().endswith((".gz", ".bgz")):
        return gzip.open(path, mode + "t", encoding="utf-8")
    return open(path, mode, encoding="utf-8")


def drop_header(line: str, keep_path_records: bool) -> bool:
    if not line.startswith("##"):
        return False
    if BCFTOOLS_HEADER_RE.match(line):
        return True
    if COMMAND_HEADER_RE.match(line):
        return True
    if not keep_path_records and any(line.startswith(prefix) for prefix in PATH_PREFIXES):
        return True
    return False


def main() -> int:
    args = parse_args()
    if args.input != "-" and not Path(args.input).exists():
        sys.stderr.write(f"error: input does not exist: {args.input}\n")
        return 1

    with open_text(args.input, "r") as src, open_text(args.output, "w") as dst:
        for line in src:
            if drop_header(line, keep_path_records=args.keep_path_records):
                continue
            dst.write(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
