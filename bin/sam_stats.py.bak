#!/usr/bin/env python3

import sys
from pathlib import Path
import pandas as pd

def parse_simple_stats(lines):
    """Parse SN lines from samtools stats"""
    sn = {}
    for l in lines:
        if l.startswith("SN"):
            parts = l.rstrip().split("\t")
            key = parts[1].rstrip(":")
            for p in parts[2:]:
                try:
                    sn[key] = float(p)
                    break
                except ValueError:
                    continue
    return sn


def main(stats_file, coverage_file):
    stats_file = Path(stats_file)
    coverage_file = Path(coverage_file)

    # -----------------------
    # Read samtools stats
    # -----------------------
    with open(stats_file) as f:
        stats_lines = [l for l in f if not l.startswith("#")]

    sn = parse_simple_stats(stats_lines)

    # -----------------------
    # Read samtools coverage
    # -----------------------
    cov = pd.read_csv(coverage_file, sep="\t")
    cov["length"] = cov["endpos"] - cov["startpos"] + 1

    genome_size = cov["length"].sum()
    covered_bases = cov["covbases"].sum()

    # Coverage-consistent mapped bases
    mapped_bases_cov = (cov["meandepth"] * cov["length"]).sum()

    # CIGAR-consistent mapped bases (preferred)
    mapped_bases_cigar = sn.get("bases mapped (cigar)", 0)

    # -----------------------
    # Metrics (RAD-correct)
    # -----------------------
    breadth_frac = covered_bases / genome_size if genome_size > 0 else 0
    breadth_pct = 100 * breadth_frac

    mapped_depth = mapped_bases_cov / covered_bases if covered_bases > 0 else 0
    global_depth = breadth_frac * mapped_depth

    dup_rate = 100 * sn["reads duplicated"] / sn["reads mapped"]
    
    map_fraction = 100 * sn["reads mapped"] / sn["raw total sequences"]
    
    mq0_rate = 100 * sn["reads MQ0"] / sn["reads mapped"]

    # -----------------------
    # Output (TSV)
    # -----------------------
    print("\t".join([
        "file",
        "genome_size",
        "mapped_bases_cigar",
        "mapped_reads_pct",
        "ref_covered_bases",
        "global_depth",
        "mapped_depth",
        "breadth_pct",
        "duplication_rate_pct",
        "mq0_rate_pct"
    ]))

    print("\t".join(map(str, [
        stats_file.name,
        f"{int(genome_size):,}",
        f"{int(mapped_bases_cigar):,}",
        f"{int(map_fraction):.2f}",
        f"{int(covered_bases):,}",
        f"{global_depth:.2f}",
        f"{mapped_depth:.2f}",
        f"{breadth_pct:.2f}",
        f"{dup_rate:.2f}",
        f"{mq0_rate:.2f}"
    ])))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(
            "USAGE:\n"
            "  radseq_samtools_qc.py <samtools.stats> <samtools.coverage.txt>\n"
        )

    main(sys.argv[1], sys.argv[2])
