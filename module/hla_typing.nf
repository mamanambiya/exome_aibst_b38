#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * HLA Typing Module
 *
 * Types HLA-A, HLA-B, HLA-C, HLA-DRB1 from WES BAM files.
 *
 * Uses OptiType (Class I: A/B/C, best WES accuracy ~98%) and
 * arcasHLA (Class I+II including DRB1, validated on RNA-seq but works on WES).
 *
 * At ~12x WES, two-field resolution (e.g., A*02:01) is reliable.
 * Three-field resolution should not be reported.
 */


/*
 * OptiType: HLA Class I typing (HLA-A, HLA-B, HLA-C)
 * Best-validated tool for WES data at two-field resolution.
 */
process optitype {
    tag "optitype_${sample_id}"
    label "medium"
    container 'docker://quay.io/biocontainers/optitype:1.3.5--hdfd78af_3'
    publishDir "${params.outDir}/hla_typing/optitype", mode: 'copy'

    input:
        tuple val(sample_id), file(bam), file(bai)

    output:
        tuple val(sample_id), path("${sample_id}_result.tsv"), emit: result

    script:
        """
        # Extract HLA reads from BAM using razers3 (bundled with OptiType)
        samtools fastq -@ 4 ${bam} \\
            -1 ${sample_id}_R1.fastq.gz \\
            -2 ${sample_id}_R2.fastq.gz \\
            -0 /dev/null -s /dev/null -n

        # Run OptiType
        OptiTypePipeline.py \\
            -i ${sample_id}_R1.fastq.gz ${sample_id}_R2.fastq.gz \\
            --dna \\
            -o . \\
            -p ${sample_id} \\
            --verbose

        # Rename output for consistency
        mv */result_*.tsv ${sample_id}_result.tsv 2>/dev/null || \\
        mv *_result.tsv ${sample_id}_result.tsv 2>/dev/null || true
        """
}


/*
 * arcasHLA: HLA Class I + II typing (HLA-A, B, C, DRB1)
 * Adds DRB1 which OptiType cannot type.
 */
process arcashla_extract {
    errorStrategy "ignore"
    tag "arcashla_extract_${sample_id}"
    label "medium"
    container 'docker://quay.io/biocontainers/arcas-hla:0.6.0--hdfd78af_2'
    containerOptions '-B /cbio/users/mamana/resources/arcashla_ref/dat:/usr/local/share/arcas-hla-0.6.0-2/dat'

    input:
        tuple val(sample_id), file(bam), file(bai)

    output:
        tuple val(sample_id), path("${sample_id}.extracted.1.fq.gz"), path("${sample_id}.extracted.2.fq.gz"), emit: reads

    script:
        """
        # Try arcasHLA extract first; if it fails (corrupted BAM records),
        # fall back to manual BAM-native extraction (avoids SAM text parsing)
        arcasHLA extract \\
            ${bam} \\
            -o . \\
            -t ${task.cpus} \\
            --unmapped \\
            -v 2>&1 | tee extract.log || true

        bam_stem=\$(basename ${bam} .bam)
        if [ "\${bam_stem}" != "${sample_id}" ]; then
            mv \${bam_stem}.extracted.1.fq.gz ${sample_id}.extracted.1.fq.gz 2>/dev/null || true
            mv \${bam_stem}.extracted.2.fq.gz ${sample_id}.extracted.2.fq.gz 2>/dev/null || true
        fi

        # Check if FASTQs are empty (indicates extract failure from corrupted reads)
        if [ ! -s ${sample_id}.extracted.1.fq.gz ] || [ \$(zcat ${sample_id}.extracted.1.fq.gz 2>/dev/null | head -1 | wc -c) -eq 0 ]; then
            echo "[fallback] arcasHLA extract produced empty FASTQs, using BAM-native extraction"
            # Extract chr6 + unmapped as BAM (bypasses SAM text parsing that chokes on bad records)
            samtools view -b -f 2 ${bam} chr6 > hla_chr6.bam 2>/dev/null || true
            samtools view -b -f 12 ${bam} > hla_unmapped.bam 2>/dev/null || true
            samtools merge -f hla_merged.bam hla_chr6.bam hla_unmapped.bam 2>/dev/null || true
            samtools sort -n -@ ${task.cpus} hla_merged.bam -o hla_sorted.bam 2>/dev/null || true
            bedtools bamtofastq -i hla_sorted.bam \\
                -fq ${sample_id}.extracted.1.fq \\
                -fq2 ${sample_id}.extracted.2.fq 2>/dev/null || true
            gzip -f ${sample_id}.extracted.1.fq ${sample_id}.extracted.2.fq
        fi
        """
}

process arcashla_genotype {
    errorStrategy "ignore"
    tag "arcashla_genotype_${sample_id}"
    label "small"
    container 'docker://quay.io/biocontainers/arcas-hla:0.6.0--hdfd78af_2'
    containerOptions '-B /cbio/users/mamana/resources/arcashla_ref/dat:/usr/local/share/arcas-hla-0.6.0-2/dat'
    publishDir "${params.outDir}/hla_typing/arcashla", mode: 'copy'

    input:
        tuple val(sample_id), file(fq1), file(fq2)

    output:
        tuple val(sample_id), path("${sample_id}.genotype.json"), emit: result

    script:
        """
        arcasHLA genotype \\
            ${fq1} ${fq2} \\
            -g A,B,C,DRB1 \\
            -o . \\
            -t ${task.cpus} \\
            -v
        """
}


/*
 * Merge OptiType + arcasHLA results into a consolidated table.
 * Uses OptiType for Class I (higher WES accuracy), arcasHLA for DRB1.
 */
process merge_hla_results {
    tag "merge_hla"
    label "small"
    publishDir "${params.outDir}/hla_typing", mode: 'copy'

    input:
        path(optitype_results)
        path(arcashla_results)

    output:
        path("hla_typing_summary.tsv"), emit: summary

    script:
        """
        python3 << 'PYEOF'
import json, glob, csv, os

# Parse OptiType results (Class I)
optitype = {}
for f in glob.glob("*_result.tsv"):
    sample = os.path.basename(f).replace("_result.tsv", "")
    with open(f) as fh:
        reader = csv.DictReader(fh, delimiter='\\t')
        for row in reader:
            optitype[sample] = {
                'HLA-A': [row.get('A1', ''), row.get('A2', '')],
                'HLA-B': [row.get('B1', ''), row.get('B2', '')],
                'HLA-C': [row.get('C1', ''), row.get('C2', '')],
            }

# Parse arcasHLA results (Class I + II)
arcas = {}
for f in glob.glob("*.genotype.json"):
    sample = os.path.basename(f).replace(".genotype.json", "")
    with open(f) as fh:
        data = json.load(fh)
        arcas[sample] = data

# Merge: OptiType for A/B/C (better WES accuracy), arcasHLA for DRB1
all_samples = sorted(set(list(optitype.keys()) + list(arcas.keys())))

with open("hla_typing_summary.tsv", "w") as out:
    out.write("Sample\\tHLA-A_1\\tHLA-A_2\\tHLA-B_1\\tHLA-B_2\\tHLA-C_1\\tHLA-C_2\\tHLA-DRB1_1\\tHLA-DRB1_2\\tSource_ClassI\\tSource_DRB1\\n")
    for sample in all_samples:
        row = [sample]
        # Class I from OptiType (preferred) or arcasHLA (fallback)
        if sample in optitype:
            for gene in ['HLA-A', 'HLA-B', 'HLA-C']:
                alleles = optitype[sample].get(gene, ['', ''])
                row.extend(alleles)
            source_i = "OptiType"
        elif sample in arcas:
            for gene in ['A', 'B', 'C']:
                alleles = arcas[sample].get(gene, ['', ''])
                row.extend(alleles if len(alleles) >= 2 else alleles + [''])
            source_i = "arcasHLA"
        else:
            row.extend([''] * 6)
            source_i = "NA"

        # DRB1 from arcasHLA only
        if sample in arcas and 'DRB1' in arcas[sample]:
            drb1 = arcas[sample]['DRB1']
            row.extend(drb1 if len(drb1) >= 2 else drb1 + [''])
            source_drb1 = "arcasHLA"
        else:
            row.extend(['', ''])
            source_drb1 = "NA"

        row.extend([source_i, source_drb1])
        out.write("\\t".join(str(x) for x in row) + "\\n")
PYEOF
        """
}


/*
 * Main workflow: HLA typing from BAMs
 */
workflow hla_typing {
    take:
        bam_ch  // channel of [sample_id, bam, bai]

    main:
        // OptiType for Class I (HLA-A, B, C)
        optitype(bam_ch)

        // arcasHLA for Class I + II (adds DRB1)
        arcashla_extract(bam_ch)
        arcashla_genotype(arcashla_extract.out.reads)

        // Merge results
        merge_hla_results(
            optitype.out.result.map { id, f -> f }.collect(),
            arcashla_genotype.out.result.map { id, f -> f }.collect()
        )

    emit:
        summary = merge_hla_results.out.summary
}
