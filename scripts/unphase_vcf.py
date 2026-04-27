#!/usr/bin/env python3
"""Write an unphased VCF stream from VCF/VCF.GZ/BCF input.

The transformation is intentionally small and transparent:

* every GT value has phased separators ('|') replaced by unphased separators ('/');
* FORMAT/PS and FORMAT/PQ are removed by default, because they describe the
  phase state that is being discarded before read-backed re-phasing;
* all other records, INFO fields, alleles, filters, and sample-level FORMAT
  values are preserved.

BCF input requires `bcftools view -Ov`. Plain VCF and VCF.GZ input are read
with Python's standard library.
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable, Iterator, TextIO

PHASE_FORMAT_KEYS = {"PS", "PQ"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert phased GT fields in a VCF/BCF to unphased GT fields."
    )
    parser.add_argument(
        "vcf",
        help="Input VCF/VCF.GZ/BCF path, or '-' for VCF text on stdin.",
    )
    parser.add_argument(
        "--keep-phase-tags",
        action="store_true",
        help="Keep FORMAT/PS and FORMAT/PQ instead of removing them.",
    )
    return parser.parse_args()


class LineSource:
    def __init__(self, path: str):
        self.path = path
        self.proc: subprocess.Popen[str] | None = None
        self.handle: TextIO | None = None

    def __enter__(self) -> Iterable[str]:
        if self.path == "-":
            self.handle = sys.stdin
            return self.handle

        lower = self.path.lower()
        if lower.endswith(".vcf"):
            self.handle = open(self.path, "rt", encoding="utf-8")
            return self.handle
        if lower.endswith(".vcf.gz") or lower.endswith(".vcf.bgz"):
            self.handle = gzip.open(self.path, "rt", encoding="utf-8")
            return self.handle

        bcftools = shutil.which("bcftools")
        if bcftools is None:
            raise SystemExit(
                "error: BCF or non-standard compressed input requires bcftools in PATH"
            )
        self.proc = subprocess.Popen(
            [bcftools, "view", "-Ov", self.path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
        )
        assert self.proc.stdout is not None
        return self.proc.stdout

    def __exit__(self, exc_type, exc, tb) -> None:
        if self.handle is not None and self.handle is not sys.stdin:
            self.handle.close()
        if self.proc is not None:
            assert self.proc.stderr is not None
            stderr = self.proc.stderr.read()
            status = self.proc.wait()
            if status != 0 and exc_type is None:
                sys.stderr.write(stderr)
                raise SystemExit(f"error: bcftools view failed with status {status}")


def unphase_sample_value(format_keys: list[str], keep_indices: list[int], sample: str) -> str:
    if sample in {".", ""}:
        return sample
    values = sample.rstrip("\n").split(":")
    out: list[str] = []
    for idx in keep_indices:
        value = values[idx] if idx < len(values) else "."
        if format_keys[idx] == "GT":
            value = value.replace("|", "/")
        out.append(value)
    return ":".join(out) if out else "."


def transform_lines(lines: Iterable[str], drop_phase_tags: bool) -> Iterator[str]:
    for line in lines:
        if line.startswith("##FORMAT=<ID=") and drop_phase_tags:
            if any(line.startswith(f"##FORMAT=<ID={key},") for key in PHASE_FORMAT_KEYS):
                continue
            # Some producers emit minimal FORMAT header lines without a comma.
            if any(line.startswith(f"##FORMAT=<ID={key}>") for key in PHASE_FORMAT_KEYS):
                continue
        if line.startswith("#"):
            yield line
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 10:
            yield line
            continue

        format_keys = fields[8].split(":")
        if drop_phase_tags:
            keep_indices = [
                i for i, key in enumerate(format_keys) if key not in PHASE_FORMAT_KEYS
            ]
        else:
            keep_indices = list(range(len(format_keys)))

        if keep_indices != list(range(len(format_keys))):
            fields[8] = ":".join(format_keys[i] for i in keep_indices) or "."

        for col in range(9, len(fields)):
            fields[col] = unphase_sample_value(format_keys, keep_indices, fields[col])
        yield "\t".join(fields) + "\n"


def main() -> int:
    args = parse_args()
    if args.vcf != "-" and not Path(args.vcf).exists():
        sys.stderr.write(f"error: input does not exist: {args.vcf}\n")
        return 1
    with LineSource(args.vcf) as lines:
        for out_line in transform_lines(lines, drop_phase_tags=not args.keep_phase_tags):
            sys.stdout.write(out_line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
