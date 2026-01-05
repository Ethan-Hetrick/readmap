#!/usr/bin/env python3

import sys
from pathlib import Path
import pandas as pd


# -----------------------------
# samtools stats parser
# -----------------------------
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


# -----------------------------
# bcftools stats parser
# -----------------------------
def parse_bcftools_stats(lines):
    """
    Extract RAD-relevant metrics from bcftools stats
    """
    out = {
        "n_variants": None,
        "n_snps": None,
        "ts_tv": None,
        "mean_dp": None,
    }

    dp_sum = 0
    dp_sites = 0

    for l in lines:
        if l.startswith("SN"):
            parts = l.rstrip().split("\t")
            key = parts[2].rstrip(":")
            val = float(parts[3])

            if key == "number of records":
                out["n_variants"] = int(val)
            elif key == "number of SNPs":
                out["n_snps"] = int(val)

        elif l.startswith("TSTV"):
            parts = l.rstrip().split("\t")
            # ts/tv ratio is column 5 (0-based index 4)
            out["ts_tv"] = float(parts[4])

        elif l.startswith("DP"):
            parts = l.rstrip().split("\t")
            depth = int(parts[2])
            n_sites = int(parts[5])
            dp_sum += depth * n_sites
            dp_sites += n_sites

    if dp_sites > 0:
        out["mean_dp"] = dp_sum / dp_sites

    return out


# -----------------------------
# main
# -----------------------------
def main(stats_file, coverage_file, bcftools_stats_file=None):

    stats_file = Path(stats_file)
    coverage_file = Path(coverage_file)

    # -----------------------
    # samtools stats
    # -----------------------
    with open(stats_file) as f:
        stats_lines = [l for l in f if not l.startswith("#")]

    sn = parse_simple_stats(stats_lines)

    # -----------------------
    # samtools coverage
    # -----------------------
    cov = pd.read_csv(coverage_file, sep="\t")
    cov.columns = cov.columns.str.lstrip("#")
    cov["length"] = cov["endpos"] - cov["startpos"] + 1

    genome_size = cov["length"].sum()
    covered_bases = cov["covbases"].sum()

    mapped_bases_cov = (cov["meandepth"] * cov["length"]).sum()
    mapped_bases_cigar = sn.get("bases mapped (cigar)", 0)

    # -----------------------
    # RAD-correct metrics
    # -----------------------
    breadth_frac = covered_bases / genome_size if genome_size else 0
    breadth_pct = 100 * breadth_frac

    mapped_depth = mapped_bases_cov / covered_bases if covered_bases else 0
    global_depth = breadth_frac * mapped_depth

    dup_rate = (
        100 * sn["bases duplicated"] / mapped_bases_cigar
        if mapped_bases_cigar else 0
    )

    map_fraction = 100 * sn["reads mapped"] / sn["raw total sequences"]
    mq0_rate = 100 * sn["reads MQ0"] / sn["reads mapped"]

    # -----------------------
    # Optional bcftools stats
    # -----------------------
    bcf = {
        "mean_dp": "NA",
        "ts_tv": "NA",
        "n_variants": "NA",
        "n_snps": "NA",
    }

    if bcftools_stats_file:
        with open(bcftools_stats_file) as f:
            bcf_lines = [l for l in f if not l.startswith("#")]

        parsed = parse_bcftools_stats(bcf_lines)

        if parsed["mean_dp"] is not None:
            bcf["mean_dp"] = f"{parsed['mean_dp']:.2f}"
        if parsed["ts_tv"] is not None:
            bcf["ts_tv"] = f"{parsed['ts_tv']:.2f}"
        if parsed["n_variants"] is not None:
            bcf["n_variants"] = str(parsed["n_variants"])
        if parsed["n_snps"] is not None:
            bcf["n_snps"] = str(parsed["n_snps"])

    # -----------------------
    # Output
    # -----------------------
    header = [
        "file",
        "genome_size",
        "mapped_bases_cigar",
        "mapped_reads_pct",
        "ref_covered_bases",
        "global_depth",
        "mapped_depth",
        "breadth_pct",
        "duplication_rate_pct",
        "mq0_rate_pct",
        "mean_variant_dp",
        "ts_tv",
        "n_variants",
        "n_snps",
    ]

    values = [
        stats_file.name,
        f"{int(genome_size):,}",
        f"{int(mapped_bases_cigar):,}",
        f"{map_fraction:.2f}",
        f"{int(covered_bases):,}",
        f"{global_depth:.2f}",
        f"{mapped_depth:.2f}",
        f"{breadth_pct:.2f}",
        f"{dup_rate:.2f}",
        f"{mq0_rate:.2f}",
        bcf["mean_dp"],
        bcf["ts_tv"],
        bcf["n_variants"],
        bcf["n_snps"],
    ]

    print("\t".join(header))
    print("\t".join(values))


# -----------------------------
# CLI
# -----------------------------
if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        sys.exit(
            "USAGE:\n"
            "  sam_stats.py <samtools.stats> <samtools.coverage.txt> [bcftools.stats]\n"
        )

    main(*sys.argv[1:])