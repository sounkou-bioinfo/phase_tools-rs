#!/usr/bin/env python3
"""Convert GTF CDS features into a BED-like codon map for phase_mnv_rs.

Output columns:

    CHROM  START0  END0  TRANSCRIPT  CODON_ID

Coordinates are 0-based half-open. Codons split by splice junctions are emitted
as multiple intervals with the same TRANSCRIPT/CODON_ID key.
"""

from __future__ import annotations

import argparse
import gzip
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable, TextIO

ATTR_RE = re.compile(r'([A-Za-z0-9_.:-]+)\s+"([^"]*)"')


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a phase_mnv_rs codon map from GTF CDS records.")
    parser.add_argument("gtf", help="Input GTF/GTF.GZ path, or '-' for stdin")
    parser.add_argument(
        "--transcript-attribute",
        default="transcript_id",
        help="GTF attribute used as transcript ID (default: transcript_id)",
    )
    parser.add_argument(
        "--feature",
        default="CDS",
        help="GTF feature type to use (default: CDS)",
    )
    parser.add_argument(
        "--positions-vcf",
        help="Optional VCF/VCF.GZ/BCF; only emit codon intervals that contain SNV positions from this file.",
    )
    return parser.parse_args()


def open_text(path: str) -> TextIO:
    if path == "-":
        return sys.stdin
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


class VariantLineSource:
    def __init__(self, path: str):
        self.path = path
        self.proc: subprocess.Popen[str] | None = None
        self.handle: TextIO | None = None

    def __enter__(self) -> Iterable[str]:
        lower = self.path.lower()
        if lower.endswith(".vcf"):
            self.handle = open(self.path, "rt", encoding="utf-8")
            return self.handle
        if lower.endswith((".vcf.gz", ".vcf.bgz")):
            self.handle = gzip.open(self.path, "rt", encoding="utf-8")
            return self.handle
        bcftools = shutil.which("bcftools")
        if bcftools is None:
            raise SystemExit("error: --positions-vcf BCF input requires bcftools in PATH")
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
        if self.handle is not None:
            self.handle.close()
        if self.proc is not None:
            assert self.proc.stderr is not None
            stderr = self.proc.stderr.read()
            status = self.proc.wait()
            if status != 0 and exc_type is None:
                sys.stderr.write(stderr)
                raise SystemExit(f"error: bcftools view failed with status {status}")


def load_snv_positions(path: str) -> dict[str, set[int]]:
    dna = set("ACGTNacgtn")
    positions: dict[str, set[int]] = defaultdict(set)
    if not Path(path).exists():
        raise SystemExit(f"error: --positions-vcf input does not exist: {path}")
    with VariantLineSource(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            chrom, pos_text, _id, ref, alts = fields[:5]
            if len(ref) != 1 or any(base not in dna for base in ref):
                continue
            if not any(len(alt) == 1 and all(base in dna for base in alt) for alt in alts.split(",")):
                continue
            try:
                pos0 = int(pos_text) - 1
            except ValueError:
                continue
            if pos0 >= 0:
                positions[chrom].add(pos0)
    return positions


def parse_attrs(text: str) -> dict[str, str]:
    attrs = {key: value for key, value in ATTR_RE.findall(text)}
    if attrs:
        return attrs
    # Fallback for less-common key=value GTF/GFF-like attributes.
    out: dict[str, str] = {}
    for item in text.rstrip(";").split(";"):
        item = item.strip()
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        out[key.strip()] = value.strip().strip('"')
    return out


def contiguous_runs(positions: list[int]) -> Iterable[tuple[int, int]]:
    if not positions:
        return
    positions.sort()
    start = prev = positions[0]
    for pos in positions[1:]:
        if pos == prev + 1:
            prev = pos
            continue
        yield start, prev + 1
        start = prev = pos
    yield start, prev + 1


def main() -> int:
    args = parse_args()
    if args.gtf != "-" and not Path(args.gtf).exists():
        sys.stderr.write(f"error: input does not exist: {args.gtf}\n")
        return 1

    positions = load_snv_positions(args.positions_vcf) if args.positions_vcf else None

    transcripts: dict[str, dict[str, object]] = defaultdict(lambda: {"segments": []})
    with open_text(args.gtf) as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != args.feature:
                continue
            attrs = parse_attrs(fields[8])
            transcript_id = attrs.get(args.transcript_attribute)
            if not transcript_id:
                continue
            try:
                start1 = int(fields[3])
                end1 = int(fields[4])
            except ValueError:
                sys.stderr.write(f"warning: skipping malformed coordinates at line {line_no}\n")
                continue
            if start1 < 1 or end1 < start1:
                sys.stderr.write(f"warning: skipping invalid interval at line {line_no}\n")
                continue
            chrom = fields[0]
            strand = fields[6]
            if strand not in {"+", "-"}:
                continue
            start0 = start1 - 1
            end0 = end1
            rec = transcripts[transcript_id]
            if "chrom" in rec and (rec["chrom"] != chrom or rec["strand"] != strand):
                sys.stderr.write(
                    f"warning: skipping transcript with inconsistent contig/strand: {transcript_id}\n"
                )
                rec["skip"] = True
                continue
            rec["chrom"] = chrom
            rec["strand"] = strand
            rec["segments"].append((start0, end0))

    for transcript_id in sorted(transcripts):
        rec = transcripts[transcript_id]
        if rec.get("skip"):
            continue
        segments = rec.get("segments", [])
        if not segments:
            continue
        chrom = str(rec["chrom"])
        strand = str(rec["strand"])
        ordered = sorted(segments, key=lambda x: x[0], reverse=(strand == "-"))
        codon_positions: list[int] = []
        codon_index = 1
        for start0, end0 in ordered:
            if strand == "+":
                iterator = range(start0, end0)
            else:
                iterator = range(end0 - 1, start0 - 1, -1)
            for pos in iterator:
                codon_positions.append(pos)
                if len(codon_positions) == 3:
                    for run_start, run_end in contiguous_runs(codon_positions):
                        if positions is not None:
                            chrom_positions = positions.get(chrom)
                            if not chrom_positions or not any(
                                pos in chrom_positions for pos in range(run_start, run_end)
                            ):
                                continue
                        print(f"{chrom}\t{run_start}\t{run_end}\t{transcript_id}\t{codon_index}")
                    codon_index += 1
                    codon_positions = []
        # Ignore incomplete terminal codons; they cannot seed a complete codon-level MNV.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
