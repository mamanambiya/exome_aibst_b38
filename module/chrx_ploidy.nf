#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * chrX Ploidy Correction Module
 *
 * Handles sex chromosome-specific processing:
 * 1. Infer sample sex from chrX heterozygosity rates
 * 2. Fix male genotypes in non-PAR regions (diploid -> haploid)
 *
 * Uses bcftools +fixploidy with GRCh38 PAR coordinates.
 */


/*
 * Infer sex from chrX heterozygosity rates.
 * Runs once per dataset on the concatenated chrX VCF (after concat_chunks_vcf).
 *
 * Output: sex_map.txt (sample\tM/F)
 *         Compatible with bcftools +fixploidy -s format
 */
process infer_sex_chrx {
    tag "infer_sex_${dataset_name}"
    label "small"
    publishDir "${params.outDir}/chrx_sex_inference", mode: 'copy'

    input:
        tuple val(dataset_name), file(chrx_vcf), file(chrx_vcf_idx)

    output:
        tuple val(dataset_name), file(sex_map)

    script:
        sex_map = "${dataset_name}_sex_map.txt"
        threshold = params.chrx_het_threshold ?: 2.0
        """
        bcftools stats -s- ${chrx_vcf} | \
        awk -F'\\t' -v thresh=${threshold} '
            /^PSC/ {
                sample = \$3
                hom_ref = \$4
                hom_alt = \$5
                het = \$6
                total = hom_ref + hom_alt + het
                het_rate = (total > 0) ? het / total * 100 : 0
                sex = (het_rate < thresh) ? "M" : "F"
                print sample "\\t" sex
            }
        ' > ${sex_map}

        n_male=\$(grep -c 'M\$' ${sex_map} || true)
        n_female=\$(grep -c 'F\$' ${sex_map} || true)
        echo "Inferred sex: \${n_male} males, \${n_female} females (threshold: ${threshold}%)" >&2
        """
}

/*
 * Fix ploidy for chrX in males using bcftools +fixploidy.
 *
 * Uses GRCh38 ploidy definition:
 *   chrX non-PAR: males = haploid (1), females = diploid (2)
 *   chrX PAR1 (chrX:10001-2781479): all diploid (2)
 *   chrX PAR2 (chrX:155701383-156030895): all diploid (2)
 *
 * Matches fill_tags tuple format: (dataset_name, vcf, sample, chrm)
 */
process fix_chrx_ploidy {
    tag "fix_ploidy_${dataset_name}_${chrm}"
    label "small"

    input:
        tuple val(dataset_name), file(vcf_file), file(dataset_sample), val(chrm), file(sex_map)

    output:
        tuple val(dataset_name), file(vcf_out), file(dataset_sample), val(chrm)

    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_ploidy_fixed.vcf.gz"
        """
        # GRCh38 ploidy definition for chrX
        cat > ploidy.txt << 'PLOIDY'
chrX 1 10000 M 2
chrX 1 10000 F 2
chrX 10001 2781479 M 2
chrX 10001 2781479 F 2
chrX 2781480 155701382 M 1
chrX 2781480 155701382 F 2
chrX 155701383 156030895 M 2
chrX 155701383 156030895 F 2
PLOIDY

        bcftools +fixploidy ${vcf_file} -Oz -o ${vcf_out} -- \
            -p ploidy.txt \
            -s ${sex_map}

        bcftools index -t ${vcf_out}
        """
}
