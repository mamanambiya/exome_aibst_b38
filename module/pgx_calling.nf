#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.pgx_nonsv_genes = ['CYP3A5', 'NAT2', 'DPYD', 'SLCO1B1', 'CYP2C9', 'CYP2C19']
params.pgx_sv_genes = ['CYP2B6', 'CYP2D6', 'UGT1A1']
params.pypgx_container = '/cbio/users/mamana/singularity-containers/pypgx_0.25.0.sif'
params.bcftools_container = '/cbio/users/mamana/singularity-containers/depot.galaxyproject.org-singularity-bcftools-1.20--h8b25389_0.img'

process fix_vcf_for_pypgx {
    tag "fix_vcf"
    label "medium"
    container "${params.bcftools_container}"
    publishDir "${params.outDir}/pypgx_results", mode: 'copy'

    input:
        path(vcf)
        path(vcf_tbi)

    output:
        path("pypgx_input_fixed.vcf.gz"), emit: vcf
        path("pypgx_input_fixed.vcf.gz.tbi"), emit: tbi

    script:
        """
        fix_vcf_ad.sh ${vcf} pypgx_input_fixed.vcf.gz
        """
}

process pypgx_depth {
    tag "pypgx_depth"
    label "bigmem"
    container "${params.pypgx_container}"
    publishDir "${params.outDir}/pypgx_results", mode: 'copy'

    input:
        path(bam_files)
        path(bai_files)

    output:
        path("depth-of-coverage.zip"), emit: depth

    script:
        """
        export PYPGX_BUNDLE=/cbio/users/mamana/resources/reference-data/reference/pypgx-bundle
        pypgx prepare-depth-of-coverage \
            depth-of-coverage.zip \
            ${bam_files} \
            --assembly GRCh38
        """
}

process pypgx_control_stats {
    tag "pypgx_control"
    label "medium"
    container "${params.pypgx_container}"
    publishDir "${params.outDir}/pypgx_results", mode: 'copy'

    input:
        path(bam_files)
        path(bai_files)

    output:
        path("control-statistics-VDR.zip"), emit: control

    script:
        """
        export PYPGX_BUNDLE=/cbio/users/mamana/resources/reference-data/reference/pypgx-bundle
        pypgx compute-control-statistics \
            VDR \
            control-statistics-VDR.zip \
            ${bam_files} \
            --assembly GRCh38
        """
}

process pypgx_call_nonsv {
    tag "pypgx_${gene}"
    label "small"
    errorStrategy "ignore"
    container "${params.pypgx_container}"
    publishDir "${params.outDir}/pypgx_results/${gene}", mode: 'copy'

    input:
        each gene
        path(vcf)
        path(vcf_tbi)

    output:
        tuple val(gene), path("${gene}-results/"), emit: results

    script:
        """
        export PYPGX_BUNDLE=/cbio/users/mamana/resources/reference-data/reference/pypgx-bundle
        pypgx run-ngs-pipeline \
            ${gene} \
            ${gene}-results/ \
            --variants ${vcf} \
            --assembly GRCh38
        """
}

process pypgx_call_sv {
    tag "pypgx_${gene}"
    label "medium"
    errorStrategy "ignore"
    container "${params.pypgx_container}"
    publishDir "${params.outDir}/pypgx_results/${gene}", mode: 'copy'

    input:
        each gene
        path(vcf)
        path(vcf_tbi)
        path(depth)
        path(control)

    output:
        tuple val(gene), path("${gene}-results/"), emit: results

    script:
        """
        export PYPGX_BUNDLE=/cbio/users/mamana/resources/reference-data/reference/pypgx-bundle
        pypgx run-ngs-pipeline \
            ${gene} \
            ${gene}-results/ \
            --variants ${vcf} \
            --depth-of-coverage ${depth} \
            --control-statistics ${control} \
            --assembly GRCh38
        """
}

process pypgx_summary {
    tag "pypgx_summary"
    label "small"
    container "${params.pypgx_container}"
    publishDir "${params.outDir}/pypgx_results", mode: 'copy'

    input:
        path(result_dirs)

    output:
        path("pgx_star_allele_summary.tsv"), emit: summary

    script:
        """
        echo -e "Gene\\tSample\\tDiplotype\\tPhenotype\\tActivity_Score" > pgx_star_allele_summary.tsv

        for dir in ${result_dirs}; do
            gene=\$(basename \$dir | sed 's/-results//')
            if [ -f \$dir/results.tsv ]; then
                tail -n+2 \$dir/results.tsv | \
                awk -v g="\$gene" -F'\\t' '{print g "\\t" \$0}' >> pgx_star_allele_summary.tsv
            fi
        done
        """
}

workflow pgx_star_alleles {
    take:
        joint_vcf
        joint_vcf_tbi
        bam_ch
        fasta
        fasta_fai

    main:
        // Fix GLnexus VCF for PyPGx
        fix_vcf_for_pypgx(joint_vcf, joint_vcf_tbi)

        bams = bam_ch.map { id, bam, bai -> bam }.collect()
        bais = bam_ch.map { id, bam, bai -> bai }.collect()

        // Depth + control in parallel (no fasta needed)
        pypgx_depth(bams, bais)
        pypgx_control_stats(bams, bais)

        // Non-SV genes: 6 in parallel
        nonsv_genes = Channel.fromList(params.pgx_nonsv_genes)
        pypgx_call_nonsv(nonsv_genes, fix_vcf_for_pypgx.out.vcf, fix_vcf_for_pypgx.out.tbi)

        // SV genes: 3 in parallel, wait for depth + control
        sv_genes = Channel.fromList(params.pgx_sv_genes)
        pypgx_call_sv(
            sv_genes,
            fix_vcf_for_pypgx.out.vcf,
            fix_vcf_for_pypgx.out.tbi,
            pypgx_depth.out.depth,
            pypgx_control_stats.out.control
        )

        // Merge all
        all_results = pypgx_call_nonsv.out.results
            .mix(pypgx_call_sv.out.results)
            .map { gene, dir -> dir }
            .collect()

        pypgx_summary(all_results)

    emit:
        summary = pypgx_summary.out.summary
}
