#!/usr/bin/env python3
"""Prepare the two RPE-1 Multiome Perturb-seq differential-analysis results for IGVF upload.

The notebook writes CSVs with a pandas row index in an unnamed first column, and writes genomic
coordinates as floats because the intermediate frames carried NaNs. The portal wants every column
named and integer coordinates written as integers, and the Hs27 arm of this submission submitted
its equivalents as gzipped TSV. So this script does exactly four things, and nothing else:

  1. drops the unnamed index column
  2. rewrites coordinate columns as integers (asserted lossless: no blanks, every value integral)
  3. adds `peak_id` and `intended_target_name` to the differential-peak table
  4. writes tab-separated, gzip with mtime=0 so the md5 is reproducible

Step 3 needs a word, because both columns restore identifiers the notebook had and dropped.

`peak_id` is the snapatac2 feature name, `chr:start-end`. The notebook tests on the peak set from
`snap.tl.merge_peaks`, whose features are exactly those strings; it then splits each one into `chr`,
`start` and `end` and drops the original. Putting it back gives the table a single-column peak key,
which is how peaks are usually reported and what a reader needs to join rows back to the peak set.
It is reassembled from the three columns, so it carries no information they do not already hold.

`intended_target_name` is the perturbed gene. The notebook keys the table on `guide_id` and the
promoter coordinates but never writes the gene ID, because its readout columns are the peak's
coordinates rather than a gene. That leaves the element identified only positionally. The guide
library carries the mapping, so it is filled in here from `guide_id` — a lookup, not a
recomputation: the library's promoter coordinates are asserted to equal the ones the notebook
already wrote on every row, and the resulting column is asserted to agree with the expression
table's own `intended_target_name` for every guide present in both.

No rounding of statistics, no reordering, no filtering. Re-run it after any re-run of
`03_differential_analyses/multiome_data_organization.ipynb`.
"""
import argparse
import csv
import gzip
import hashlib
import os
from pathlib import Path

SRC = Path(os.environ.get("MULTIOME_RESULTS", Path(__file__).parent / ".." / "results")).resolve()
OUT = Path(__file__).parent / ".." / "results" / "upload"

INT_COLS = {"start", "end", "intended_target_start", "intended_target_end"}

GUIDES = Path(os.environ.get("MULTIOME_GUIDE_LIBRARY",
                             Path(__file__).parent / "multiome_paper_igvf_guides.csv"))

# source CSV -> (submission name, whether to add intended_target_name from the guide library)
FILES = {
    "multiome_paper_guide_effect_matrix.csv":
        ("rpe1_multiome_element_level_gex_results.tsv.gz", False),
    "multiome_paper_differential_peaks_by_guide.csv":
        ("rpe1_multiome_element_level_atac_results.tsv.gz", True),
}


def to_int(value, column, line):
    number = float(value)
    if not number.is_integer():
        raise ValueError(f"{column} line {line}: {value!r} is not an integer coordinate")
    return str(int(number))


def guide_targets(path):
    """guide_id -> (intended_target_name, chr, start, end) from the guide library."""
    targets = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            if not row["intended_target_name"]:
                continue          # the three non-targeting controls perturb no promoter
            targets[row["guide_id"]] = (
                row["intended_target_name"],
                row["intended_target_chr"],
                to_int(row["intended_target_start"], "intended_target_start", 0),
                to_int(row["intended_target_end"], "intended_target_end", 0),
            )
    return targets


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=OUT,
                    help="directory to write the submission copies into")
    out_dir = ap.parse_args().out
    out_dir.mkdir(parents=True, exist_ok=True)
    targets = guide_targets(GUIDES)
    for src_name, (out_name, add_target_name) in FILES.items():
        src, out = SRC / src_name, out_dir / out_name
        with open(src, newline="") as fh:
            reader = csv.reader(fh)
            header = next(reader)
            if header[0] != "":
                raise ValueError(f"{src_name}: expected an unnamed index column, got {header[0]!r}")
            header = header[1:]
            int_at = [i for i, c in enumerate(header) if c in INT_COLS]
            source_header = list(header)          # for error messages; header may gain a column
            if add_target_name:
                if "intended_target_name" in header:
                    raise ValueError(f"{src_name}: already has intended_target_name")
                # intended_target_name sits immediately before the coordinates it names, as in
                # the expression table; peak_id leads the peak columns it is assembled from
                insert_at = header.index("intended_target_chr")
                guide_at = header.index("guide_id")
                coord_at = [header.index("intended_target_" + c) for c in ("chr", "start", "end")]
                peak_at = [header.index(c) for c in ("chr", "start", "end")]
                peak_id_at = header.index("chr")
                header = header[:insert_at] + ["intended_target_name"] + header[insert_at:]
                header = header[:peak_id_at] + ["peak_id"] + header[peak_id_at:]
            rows = 0
            with gzip.GzipFile(out, "wb", mtime=0) as gz:
                gz.write(("\t".join(header) + "\n").encode())
                for line, row in enumerate(reader, start=2):
                    row = row[1:]
                    for i in int_at:
                        row[i] = to_int(row[i], source_header[i], line)
                    if add_target_name:
                        guide = row[guide_at]
                        if guide not in targets:
                            raise ValueError(f"{src_name} line {line}: {guide!r} not in {GUIDES.name}")
                        name, *coords = targets[guide]
                        written = [row[i] for i in coord_at]
                        if written != coords:
                            raise ValueError(
                                f"{src_name} line {line}: {guide!r} promoter is {written} in the "
                                f"results but {coords} in {GUIDES.name}")
                        row = row[:insert_at] + [name] + row[insert_at:]
                        chrom, start, end = (row[i] for i in peak_at)
                        row = row[:peak_id_at] + [f"{chrom}:{start}-{end}"] + row[peak_id_at:]
                    gz.write(("\t".join(row) + "\n").encode())
                    rows += 1
        md5 = hashlib.md5(out.read_bytes()).hexdigest()
        print(f"{out.name}\t{rows} rows\t{out.stat().st_size} bytes\tmd5 {md5}")


if __name__ == "__main__":
    main()
