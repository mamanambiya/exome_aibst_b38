#!/bin/bash
set -euo pipefail
INPUT=$1
OUTPUT=$2

# Strip PL/GQ/RNC, set missing GT to 0/0, then fix AD "." with awk
bcftools annotate -x FORMAT/PL,FORMAT/RNC,FORMAT/GQ "$INPUT" | \
    bcftools +setGT -- -t . -n 0 | \
    awk -F"\t" -v OFS="\t" '
    /^#/ {print; next}
    {
        split($9, fmt, ":")
        ad_idx = -1
        for (f in fmt) { if (fmt[f] == "AD") ad_idx = f }
        if (ad_idx > 0) {
            for (i = 10; i <= NF; i++) {
                split($i, vals, ":")
                if (vals[ad_idx] == ".") vals[ad_idx] = "0,0"
                $i = vals[1]
                for (v = 2; v <= length(vals); v++) $i = $i ":" vals[v]
            }
        }
        print
    }' | bgzip > "$OUTPUT"
tabix -p vcf "$OUTPUT"
