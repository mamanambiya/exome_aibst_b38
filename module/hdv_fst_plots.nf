/*
 * Additional HDV analysis processes:
 * - PGx-specific FST computation
 * - FST distribution plot
 * - HDV supplementary table generation
 */

process compute_pgx_fst {
    tag "pgx_fst_${dataset}"
    label "bigmem"
    publishDir "${params.outDir}/fst_pgx_specific/${dataset}", mode: 'copy'

    input:
        tuple val(dataset), file(vcf_file), file(pop_sample_files), file(pgx_bed)

    output:
        tuple val(dataset), file("pgx_pairwise_fst.tsv"), file("pgx_max_fst_per_variant.tsv")

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # Extract PGx-only variants
    bcftools view -R ${pgx_bed} ${vcf_file} -Oz -o pgx_only.vcf.gz
    bcftools index pgx_only.vcf.gz

    # Compute pairwise FST for all population pairs
    echo -e "POP1\\tPOP2\\tN1\\tN2\\tWeighted_FST\\tMean_FST" > pgx_pairwise_fst.tsv

    pop_files=(${pop_sample_files})
    for ((i=0; i<\${#pop_files[@]}; i++)); do
        pop1=\$(basename \${pop_files[i]} .samples.txt)
        n1=\$(wc -l < \${pop_files[i]})
        for ((j=i+1; j<\${#pop_files[@]}; j++)); do
            pop2=\$(basename \${pop_files[j]} .samples.txt)
            n2=\$(wc -l < \${pop_files[j]})
            if [ "\$n1" -gt 0 ] && [ "\$n2" -gt 0 ]; then
                output=\$(vcftools --gzvcf pgx_only.vcf.gz \\
                    --weir-fst-pop \${pop_files[i]} \\
                    --weir-fst-pop \${pop_files[j]} \\
                    --out fst_\${pop1}_\${pop2} 2>&1)
                weighted=\$(echo "\$output" | grep "weighted" | awk '{print \$NF}')
                mean=\$(echo "\$output" | grep "mean" | awk '{print \$NF}')
                echo -e "\${pop1}\\t\${pop2}\\t\${n1}\\t\${n2}\\t\${weighted}\\t\${mean}" >> pgx_pairwise_fst.tsv
            fi
        done
    done

    # Compute max FST per variant across all pairs
    python3 << 'PYEOF'
import os, glob

files = sorted(glob.glob("fst_*.weir.fst"))
max_fst = {}
max_pair = {}

for f in files:
    pair = os.path.basename(f).replace("fst_", "").replace(".weir.fst", "")
    with open(f) as fh:
        next(fh)
        for line in fh:
            parts = line.strip().split("\\t")
            if len(parts) < 3:
                continue
            chrom, pos, fst_val = parts[0], parts[1], parts[2]
            if fst_val in ("-nan", "nan"):
                continue
            try:
                fst = float(fst_val)
            except ValueError:
                continue
            key = (chrom, pos)
            if key not in max_fst or fst > max_fst[key]:
                max_fst[key] = fst
                max_pair[key] = pair

with open("pgx_max_fst_per_variant.tsv", "w") as out:
    out.write("CHROM\\tPOS\\tMAX_FST\\tMAX_PAIR\\n")
    for key in sorted(max_fst.keys(), key=lambda x: (x[0], int(x[1]))):
        out.write(f"{key[0]}\\t{key[1]}\\t{max_fst[key]:.6f}\\t{max_pair[key]}\\n")
PYEOF
    """
}

process plot_fst_distribution {
    tag "fst_dist_${dataset}"
    label "rplot"
    publishDir "${params.outDir}/REPORTS/fst/${dataset}", mode: 'copy'

    input:
        tuple val(dataset), file(max_fst_file)

    output:
        tuple val(dataset), file("fst_pgx_distribution.pdf"), file("fst_pgx_distribution.png")

    script:
    template "plot_fst_distribution.R"
}

process generate_hdv_supp_table {
    tag "hdv_supp_${dataset}"
    label "small"
    publishDir "${params.outDir}/REPORTS/hdv/${dataset}", mode: 'copy'

    input:
        tuple val(dataset), file(hdv_file)

    output:
        tuple val(dataset), file("hdv_supplementary_table.csv")

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # Generate supplementary table from HDV output
    # Add header and format for publication
    echo "Gene,rsID,Chromosome,Position,REF,ALT,Effect,Max_FST,Fold_Change,Fisher_P,ClinVar,PharmGKB,GWAS_Catalog" > hdv_supplementary_table.csv

    # Parse the combined HDV file (format depends on pipeline output)
    if [ -f ${hdv_file} ]; then
        tail -n+2 ${hdv_file} >> hdv_supplementary_table.csv
    fi
    """
}
