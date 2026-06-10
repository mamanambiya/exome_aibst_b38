#!/usr/bin/env python3
"""
Infer sample sex from chrX heterozygosity rates.

Males (hemizygous chrX) have very low het rates (<1-2%).
Females (diploid chrX) have normal het rates (>2-3%).

Input:  bcftools stats -s- output (stdin or file)
Output: TSV with columns: sample, sex, het_rate
        sex = 1 (male) or 2 (female)

Usage:
    bcftools stats -s- chrX.vcf.gz | python3 infer_sex_from_chrx.py > sex_map.txt
    python3 infer_sex_from_chrx.py --threshold 2.0 < bcftools_stats.txt > sex_map.txt
"""

import sys
import argparse


def main():
    parser = argparse.ArgumentParser(description="Infer sex from chrX het rates")
    parser.add_argument("input", nargs="?", default="-", help="bcftools stats output")
    parser.add_argument("--threshold", type=float, default=2.0,
                        help="Het rate %% cutoff: below=male, above=female (default: 2.0)")
    parser.add_argument("--output", "-o", default="-", help="Output file")
    args = parser.parse_args()

    infile = sys.stdin if args.input == "-" else open(args.input)
    outfile = sys.stdout if args.output == "-" else open(args.output, "w")

    samples = []
    for line in infile:
        if not line.startswith("PSC\t"):
            continue
        fields = line.strip().split("\t")
        sample = fields[2]
        hom_ref = int(fields[3])
        hom_alt = int(fields[4])
        het = int(fields[5])
        total = hom_ref + hom_alt + het
        het_rate = (het / total * 100) if total > 0 else 0
        sex = 1 if het_rate < args.threshold else 2
        samples.append((sample, sex, het_rate))

    n_male = sum(1 for _, s, _ in samples if s == 1)
    n_female = sum(1 for _, s, _ in samples if s == 2)
    print(f"# Inferred sex: {n_male} males, {n_female} females (threshold: {args.threshold}%)",
          file=sys.stderr)

    outfile.write("sample\tsex\thet_rate\n")
    for sample, sex, het_rate in samples:
        outfile.write(f"{sample}\t{sex}\t{het_rate:.2f}\n")

    if infile is not sys.stdin:
        infile.close()
    if outfile is not sys.stdout:
        outfile.close()


if __name__ == "__main__":
    main()
