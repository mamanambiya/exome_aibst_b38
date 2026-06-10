#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Pathogenicity Score Annotation Module
 *
 * Adds variant pathogenicity predictions to VCFs:
 * 1. CADD v1.7 (PHRED-scaled deleteriousness)
 * 2. AlphaMissense (protein structure-based pathogenicity)
 * 3. dbNSFP v4.5c (SIFT, MetaRNN, MutationTaster, GERP++, phastCons, phyloP)
 *
 * Slots into annotation chain after ClinVar, before clean_invalid_variants.
 * Tuple format: (dataset_name, vcf_file, dataset_sample, chrm, snpeff_database)
 */


/*
 * Annotate with CADD v1.7 PHRED scores.
 *
 * CADD files use bare chromosome names (1,2,...) while pipeline VCFs use chr prefix.
 * Handles this via: strip chr -> annotate -> restore chr.
 *
 * Adds INFO fields: CADD_RAW, CADD_PHRED
 */
process annotate_cadd {
    tag "cadd_${dataset_name}_${chrm}"
    label "medium"
    cpus 4

    input:
        tuple val(dataset_name), file(vcf_file), file(dataset_sample), val(chrm), val(snpeff_database)

    output:
        tuple val(dataset_name), file(vcf_out), file(dataset_sample), val(chrm), val(snpeff_database)

    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_cadd.vcf.gz"
        cadd_snvs = params.cadd_snvs
        cadd_indels = params.cadd_indels
        """
        # Header for CADD annotations
        cat > cadd.hdr << 'HDR'
##INFO=<ID=CADD_RAW,Number=A,Type=Float,Description="CADD v1.7 raw score">
##INFO=<ID=CADD_PHRED,Number=A,Type=Float,Description="CADD v1.7 PHRED-scaled score">
HDR

        # Chromosome name mapping: chr -> bare (for CADD lookup)
        for i in \$(seq 1 22) X Y M; do echo "chr\${i} \${i}"; done > chr_to_num.txt
        for i in \$(seq 1 22) X Y M; do echo "\${i} chr\${i}"; done > num_to_chr.txt

        tabix ${vcf_file}

        # Strip chr prefix -> annotate SNVs -> annotate indels -> restore chr prefix
        bcftools annotate --rename-chrs chr_to_num.txt ${vcf_file} -Ou | \\
        bcftools annotate -a ${cadd_snvs} -h cadd.hdr \\
            -c CHROM,POS,REF,ALT,CADD_RAW,CADD_PHRED -Ou | \\
        bcftools annotate -a ${cadd_indels} -h cadd.hdr \\
            -c CHROM,POS,REF,ALT,CADD_RAW,CADD_PHRED -Ou | \\
        bcftools annotate --rename-chrs num_to_chr.txt --threads ${task.cpus} -Oz -o ${vcf_out}

        tabix ${vcf_out}
        """
}


/*
 * Annotate with AlphaMissense pathogenicity predictions.
 *
 * AlphaMissense uses chr prefix — matches pipeline VCFs directly.
 *
 * Adds INFO fields: AM_PATHOGENICITY, AM_CLASS
 */
process annotate_alphamissense {
    tag "alphamissense_${dataset_name}_${chrm}"
    label "medium"
    cpus 4

    input:
        tuple val(dataset_name), file(vcf_file), file(dataset_sample), val(chrm), val(snpeff_database)

    output:
        tuple val(dataset_name), file(vcf_out), file(dataset_sample), val(chrm), val(snpeff_database)

    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_am.vcf.gz"
        am_db = params.alphamissense_annot ?: params.alphamissense
        """
        cat > am.hdr << 'HDR'
##INFO=<ID=AM_PATHOGENICITY,Number=A,Type=Float,Description="AlphaMissense pathogenicity score (0-1)">
##INFO=<ID=AM_CLASS,Number=A,Type=String,Description="AlphaMissense classification (likely_benign/ambiguous/likely_pathogenic)">
HDR

        tabix ${vcf_file}

        bcftools annotate -a ${am_db} -h am.hdr \\
            -c CHROM,POS,REF,ALT,-,-,-,-,AM_PATHOGENICITY,AM_CLASS \\
            --threads ${task.cpus} -Oz -o ${vcf_out} ${vcf_file}

        tabix ${vcf_out}
        """
}


/*
 * Annotate with dbNSFP v4.5c functional predictions.
 *
 * Uses SnpSift dbnsfp to add multiple prediction scores per variant.
 * Selected fields are most relevant for pharmacogenomics:
 *   SIFT, SIFT4G, MetaRNN, MutationTaster, FATHMM, PROVEAN, GERP++, phastCons, phyloP
 *
 * Adds INFO fields: dbNSFP_* prefix for each selected field
 */
process annotate_dbnsfp {
    tag "dbnsfp_${dataset_name}_${chrm}"
    label "medium"
    cpus 4

    input:
        tuple val(dataset_name), file(vcf_file), file(dataset_sample), val(chrm), val(snpeff_database)

    output:
        tuple val(dataset_name), file(vcf_out), file(dataset_sample), val(chrm), val(snpeff_database)

    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_dbnsfp.vcf.gz"
        dbnsfp_db = params.dbnsfp_db
        dbnsfp_fields = [
            'SIFT_pred',
            'SIFT4G_pred',
            'MetaRNN_pred',
            'MetaRNN_score',
            'GERP++_RS',
            'phastCons100way_vertebrate',
            'phyloP100way_vertebrate',
            'MutationTaster_pred',
            'FATHMM_pred',
            'PROVEAN_pred',
            'Interpro_domain',
            'Uniprot_acc'
        ].join(',')
        """
        set -o pipefail
        tabix ${vcf_file}

        SnpSift -Xmx6g dbnsfp \\
            -db ${dbnsfp_db} \\
            -f ${dbnsfp_fields} \\
            ${vcf_file} | \\
        bgzip -@ ${task.cpus} -c > ${vcf_out}

        tabix ${vcf_out}
        """
}
