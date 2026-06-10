#!/usr/bin/env nextflow
nextflow.enable.dsl=2

"""
Author: Mamana M., Claude AI Assistant
Affiliation: University of Cape Town
Aim: Advanced population structure and selection scan analyses (ADMIXTURE, iHS, XP-EHH, Tajima's D)
Date: 2026-01-26
"""

// ===== ADMIXTURE PROCESSES =====

process ped_to_bed {
    tag "ped_to_bed_${dataset}"
    label "medium"

    input:
        tuple val(dataset), path(ped_file), path(map_file), path(sample_file)
    output:
        tuple val(dataset), path("${base}.bed"), path("${base}.bim"), path("${base}.fam")
    script:
        base = ped_file.getSimpleName()
        """
        plink --file ${base} --make-bed --out ${base}
        """
}

process run_admixture {
    tag "admixture_${dataset}_K${k}"
    label "bigmem"
    publishDir "${params.outDir}/admixture/${dataset}", mode: 'copy'

    input:
        tuple val(dataset), path(bed_file), path(bim_file), path(fam_file), val(k)
    output:
        tuple val(dataset), val(k), path("${base}.${k}.Q"), path("${base}.${k}.P"), path("${base}.${k}.log")
    script:
        base = bed_file.getSimpleName()
        """
        admixture --cv -j${task.cpus} ${bed_file} ${k} | tee ${base}.${k}.log

        if [ ! -f "${base}.${k}.Q" ]; then
            echo "Error: ADMIXTURE failed to generate Q file for K=${k}"
            exit 1
        fi
        """
}

process extract_cv_error {
    tag "cv_error_${dataset}_K${k}"
    label "small"
    publishDir "${params.outDir}/admixture/${dataset}", mode: 'copy'

    input:
        tuple val(dataset), val(k), path(q_file), path(p_file), path(log_file)
    output:
        tuple val(dataset), val(k), path("${dataset}_K${k}_cv_error.txt")
    script:
        """
        # Extract CV error from log file
        grep "CV error" ${log_file} | awk '{print "${k}\\t" \$NF}' > ${dataset}_K${k}_cv_error.txt

        # Verify CV error was extracted
        if [ ! -s ${dataset}_K${k}_cv_error.txt ]; then
            echo "${k}\\tNA" > ${dataset}_K${k}_cv_error.txt
        fi
        """
}

process combine_cv_errors {
    tag "combine_cv_${dataset}"
    label "small"
    publishDir "${params.outDir}/admixture/${dataset}", mode: 'copy'

    input:
        tuple val(dataset), path(cv_files)
    output:
        tuple val(dataset), path("${dataset}_cv_errors_all_K.txt"), path("${dataset}_optimal_K.txt")
    script:
        """
        # Combine all CV errors
        echo -e "K\\tCV_Error" > ${dataset}_cv_errors_all_K.txt
        cat ${cv_files} | sort -n >> ${dataset}_cv_errors_all_K.txt

        # Find optimal K (lowest CV error)
        tail -n +2 ${dataset}_cv_errors_all_K.txt | \
            grep -v "NA" | \
            sort -k2 -n | \
            head -1 | \
            awk '{print "Optimal K: " \$1 "\\nCV Error: " \$2}' > ${dataset}_optimal_K.txt

        # If no valid CV errors, report error
        if [ ! -s ${dataset}_optimal_K.txt ]; then
            echo "Error: No valid CV errors found" > ${dataset}_optimal_K.txt
        fi
        """
}

process plot_admixture {
    tag "plot_admixture_${dataset}"
    label "bigmem"
    publishDir "${params.outDir}/admixture/${dataset}/plots", mode: 'copy'

    input:
        tuple val(dataset), path(bed_file), path(bim_file), path(fam_file), path(q_files), path(cv_summary)
    output:
        tuple val(dataset), path("${dataset}_barplot_K*.pdf"), path("${dataset}_cv_error_plot.pdf")
    script:
        base = bed_file.getSimpleName()
        """
        # Run R plotting script
        Rscript ${params.template_dir}/plot_admixture.R \
            ${dataset} \
            ${params.admixture_k_min} \
            ${params.admixture_k_max} \
            ${base}

        # Verify plots were created
        if [ ! -f ${dataset}_barplot_K${params.admixture_k_min}.pdf ]; then
            echo "Warning: ADMIXTURE plots may be incomplete"
        fi
        """
}

// ===== SELECTION SCAN PROCESSES =====

process phase_vcf {
    tag "phase_${pop}_chr${chr}"
    label "bigmem"
    publishDir "${params.outdir}/selection/phased/${pop}", mode: 'copy'

    input:
        tuple val(pop), val(chr), path(vcf_file)
    output:
        tuple val(pop), val(chr), path("${pop}_chr${chr}_phased.vcf.gz")
    script:
        """
        # Phase VCF using shapeit4
        genetic_map=${params.genetic_map_dir}/chr${chr}.b38.gmap.gz

        if [ ! -f "\${genetic_map}" ]; then
            echo "Error: Genetic map not found: \${genetic_map}"
            exit 1
        fi

        shapeit4 \
            --input ${vcf_file} \
            --map \${genetic_map} \
            --region ${chr} \
            --output ${pop}_chr${chr}_phased.vcf.gz \
            --thread ${task.cpus}

        # Index output
        tabix -f ${pop}_chr${chr}_phased.vcf.gz
        """
}

process vcf_to_selscan {
    tag "vcf_to_selscan_${pop}_chr${chr}"
    label "medium"
    publishDir "${params.outdir}/selection/selscan_input/${pop}", mode: 'copy'

    input:
        tuple val(pop), val(chr), path(phased_vcf)
    output:
        tuple val(pop), val(chr), path("${pop}_chr${chr}.hap.gz"), path("${pop}_chr${chr}.map.gz")
    script:
        """
        # Convert phased VCF to selscan format
        # Extract haplotypes (one column per haplotype)
        bcftools query -f '%POS[\\t%GT]\\n' ${phased_vcf} | \
            awk '{
                printf "%d", \$1;
                for(i=2; i<=NF; i++) {
                    gsub(/\\|/, "", \$i);
                    printf "\\t%s", \$i;
                }
                printf "\\n";
            }' | gzip > ${pop}_chr${chr}.hap.gz

        # Create map file (chr, rsid, genetic_distance, physical_position)
        bcftools query -f '%CHROM\\t%ID\\t0\\t%POS\\n' ${phased_vcf} | gzip > ${pop}_chr${chr}.map.gz
        """
}

process run_ihs {
    tag "ihs_${pop}_chr${chr}"
    label "bigmem"
    publishDir "${params.outdir}/selection/ihs/${pop}", mode: 'copy'

    input:
        tuple val(pop), val(chr), path(hap_file), path(map_file)
    output:
        tuple val(pop), val(chr), path("${pop}_chr${chr}.ihs.out")
    script:
        """
        # Run selscan iHS
        selscan --ihs \
            --hap <(zcat ${hap_file}) \
            --map <(zcat ${map_file}) \
            --out ${pop}_chr${chr} \
            --threads ${task.cpus}
        """
}

process normalize_ihs {
    tag "normalize_ihs_${pop}"
    label "medium"
    publishDir "${params.outdir}/selection/ihs/${pop}", mode: 'copy'

    input:
        tuple val(pop), path(ihs_files)
    output:
        tuple val(pop), path("${pop}_ihs_normalized.out")
    script:
        """
        # Combine all chromosome iHS files
        cat ${ihs_files} > ${pop}_ihs_all_chr.out

        # Normalize iHS scores
        norm --ihs --files ${pop}_ihs_all_chr.out --bp-win

        # Rename normalized output
        mv ${pop}_ihs_all_chr.out.norm ${pop}_ihs_normalized.out
        """
}

process identify_ihs_outliers {
    tag "ihs_outliers_${pop}"
    label "medium"
    publishDir "${params.outdir}/selection/ihs/${pop}", mode: 'copy'

    input:
        tuple val(pop), path(ihs_normalized)
    output:
        tuple val(pop), path("${pop}_ihs_outliers.txt"), path("${pop}_ihs_pgx_overlap.txt")
    script:
        threshold = params.ihs_threshold
        """
        # Identify outliers (|iHS| > threshold)
        awk -v thresh=${threshold} '\$6 > thresh || \$6 < -thresh' ${ihs_normalized} | \
            sort -k6 -rn > ${pop}_ihs_outliers.txt

        # Check for pharmacogene overlaps
        if [ -f "${params.pgx_bed}" ]; then
            bedtools intersect \
                -a <(awk '{print "chr"\$2"\\t"\$3"\\t"\$3+1"\\t"\$6}' ${pop}_ihs_outliers.txt) \
                -b ${params.pgx_bed} \
                -wa -wb > ${pop}_ihs_pgx_overlap.txt
        else
            touch ${pop}_ihs_pgx_overlap.txt
        fi

        # Report statistics
        echo "Total outliers (|iHS| > ${threshold}): \$(wc -l < ${pop}_ihs_outliers.txt)"
        if [ -f "${params.pgx_bed}" ]; then
            echo "Pharmacogene overlaps: \$(wc -l < ${pop}_ihs_pgx_overlap.txt)"
        fi
        """
}

process run_xpehh {
    tag "xpehh_${pop1}_${pop2}_chr${chr}"
    label "bigmem"
    publishDir "${params.outdir}/selection/xpehh/${pop1}_${pop2}", mode: 'copy'

    input:
        tuple val(pop1), val(pop2), val(chr), path(hap1), path(map1), path(hap2), path(map2)
    output:
        tuple val(pop1), val(pop2), val(chr), path("${pop1}_${pop2}_chr${chr}.xpehh.out")
    script:
        """
        # Run selscan XP-EHH
        selscan --xpehh \
            --hap <(zcat ${hap1}) \
            --map <(zcat ${map1}) \
            --hap-ref <(zcat ${hap2}) \
            --map-ref <(zcat ${map2}) \
            --out ${pop1}_${pop2}_chr${chr} \
            --threads ${task.cpus}
        """
}

process normalize_xpehh {
    tag "normalize_xpehh_${pop1}_${pop2}"
    label "medium"
    publishDir "${params.outdir}/selection/xpehh/${pop1}_${pop2}", mode: 'copy'

    input:
        tuple val(pop1), val(pop2), path(xpehh_files)
    output:
        tuple val(pop1), val(pop2), path("${pop1}_${pop2}_xpehh_normalized.out")
    script:
        """
        # Combine all chromosome XP-EHH files
        cat ${xpehh_files} > ${pop1}_${pop2}_xpehh_all_chr.out

        # Normalize XP-EHH scores
        norm --xpehh --files ${pop1}_${pop2}_xpehh_all_chr.out --bp-win

        # Rename normalized output
        mv ${pop1}_${pop2}_xpehh_all_chr.out.norm ${pop1}_${pop2}_xpehh_normalized.out
        """
}

process run_tajimas_d {
    tag "tajd_${pop}_chr${chr}"
    label "medium"
    publishDir "${params.outdir}/selection/tajimas_d/${pop}", mode: 'copy'

    input:
        tuple val(pop), val(chr), path(vcf_file)
    output:
        tuple val(pop), val(chr), path("${pop}_chr${chr}.Tajima.D")
    script:
        window_size = params.tajd_window_size
        """
        # Calculate Tajima's D using vcftools
        vcftools \
            --gzvcf ${vcf_file} \
            --TajimaD ${window_size} \
            --out ${pop}_chr${chr}

        # Verify output
        if [ ! -f ${pop}_chr${chr}.Tajima.D ]; then
            echo "Error: Tajima's D calculation failed"
            exit 1
        fi
        """
}

process combine_tajimas_d {
    tag "combine_tajd_${pop}"
    label "small"
    publishDir "${params.outdir}/selection/tajimas_d/${pop}", mode: 'copy'

    input:
        tuple val(pop), path(tajd_files)
    output:
        tuple val(pop), path("${pop}_tajimas_d_all_chr.txt")
    script:
        """
        # Combine all chromosome Tajima's D files
        head -1 ${tajd_files[0]} > ${pop}_tajimas_d_all_chr.txt
        for file in ${tajd_files}; do
            tail -n +2 \$file >> ${pop}_tajimas_d_all_chr.txt
        done

        # Report statistics
        echo "Total windows: \$(tail -n +2 ${pop}_tajimas_d_all_chr.txt | wc -l)"
        echo "Mean Tajima's D: \$(tail -n +2 ${pop}_tajimas_d_all_chr.txt | awk '{sum+=\$4; n++} END {print sum/n}')"
        """
}

process plot_selection_scans {
    tag "plot_selection_${pop}"
    label "bigmem"
    publishDir "${params.outdir}/selection/plots/${pop}", mode: 'copy'

    input:
        tuple val(pop), path(ihs_file), path(tajd_file)
    output:
        tuple val(pop), path("${pop}_ihs_manhattan.pdf"), path("${pop}_tajd_manhattan.pdf")
    script:
        """
        # Run R plotting script
        Rscript ${params.template_dir}/plot_selection_scans.R \
            ${pop} \
            ${ihs_file} \
            ${tajd_file} \
            ${params.pgx_bed} \
            .

        # Verify plots were created
        if [ ! -f ${pop}_ihs_manhattan.pdf ]; then
            echo "Warning: iHS Manhattan plot not created"
        fi
        if [ ! -f ${pop}_tajd_manhattan.pdf ]; then
            echo "Warning: Tajima's D Manhattan plot not created"
        fi
        """
}

// ===== WORKFLOWS =====

workflow admixture {
    take:
        ped_data  // Channel: [dataset, ped, map, sample] from vcf_to_plink1

    main:
        // Convert PED/MAP to binary PLINK BED/BIM/FAM
        ped_to_bed(ped_data)

        // Expand K values
        k_values = Channel.from(params.admixture_k_min..params.admixture_k_max)

        // Combine dataset with K values
        admixture_input = ped_to_bed.out
            .combine(k_values)

        // Run ADMIXTURE for each K
        run_admixture(admixture_input)

        // Extract CV errors
        extract_cv_error(run_admixture.out)

        // Group CV errors by dataset
        cv_errors_grouped = extract_cv_error.out
            .map { dataset, k, cv_file -> tuple(dataset, cv_file) }
            .groupTuple()

        // Combine CV errors
        combine_cv_errors(cv_errors_grouped)

        // Group Q files for plotting
        q_files_grouped = run_admixture.out
            .map { dataset, k, q_file, p_file, log_file -> tuple(dataset, q_file) }
            .groupTuple()

        // Combine with BED/BIM/FAM files for plotting
        plot_input = ped_to_bed.out
            .join(q_files_grouped)
            .join(combine_cv_errors.out.map { dataset, cv_all, opt_k -> tuple(dataset, cv_all) })

        // Generate plots
        plot_admixture(plot_input)

    emit:
        q_files = run_admixture.out
        optimal_k = combine_cv_errors.out
        plots = plot_admixture.out
}

workflow selection_scans {
    take:
        vcf_data  // Channel: [pop, chr, vcf]

    main:
        // Phase VCFs
        phase_vcf(vcf_data)

        // Convert to selscan format
        vcf_to_selscan(phase_vcf.out)

        // iHS Analysis
        run_ihs(vcf_to_selscan.out)

        // Group iHS results by population
        ihs_grouped = run_ihs.out
            .map { pop, chr, ihs_file -> tuple(pop, ihs_file) }
            .groupTuple()

        // Normalize iHS
        normalize_ihs(ihs_grouped)

        // Identify outliers
        identify_ihs_outliers(normalize_ihs.out)

        // XP-EHH Analysis (between population pairs)
        if (params.run_xpehh) {
            // Generate population pairs
            pop_pairs = vcf_to_selscan.out
                .map { pop, chr, hap, map -> tuple(pop, chr, hap, map) }
                .combine(vcf_to_selscan.out.map { pop, chr, hap, map -> tuple(pop, chr, hap, map) })
                .filter { pop1, chr1, hap1, map1, pop2, chr2, hap2, map2 ->
                    pop1 < pop2 && chr1 == chr2
                }
                .map { pop1, chr1, hap1, map1, pop2, chr2, hap2, map2 ->
                    tuple(pop1, pop2, chr1, hap1, map1, hap2, map2)
                }

            // Run XP-EHH
            run_xpehh(pop_pairs)

            // Group by population pair
            xpehh_grouped = run_xpehh.out
                .map { pop1, pop2, chr, xpehh_file -> tuple(pop1, pop2, xpehh_file) }
                .groupTuple()

            // Normalize XP-EHH
            normalize_xpehh(xpehh_grouped)
        }

        // Tajima's D Analysis
        run_tajimas_d(vcf_data)

        // Group Tajima's D by population
        tajd_grouped = run_tajimas_d.out
            .map { pop, chr, tajd_file -> tuple(pop, tajd_file) }
            .groupTuple()

        // Combine Tajima's D
        combine_tajimas_d(tajd_grouped)

        // Plotting
        plot_input = normalize_ihs.out
            .join(combine_tajimas_d.out)

        plot_selection_scans(plot_input)

    emit:
        phased_vcfs = phase_vcf.out
        ihs_results = normalize_ihs.out
        ihs_outliers = identify_ihs_outliers.out
        xpehh_results = params.run_xpehh ? normalize_xpehh.out : Channel.empty()
        tajd_results = combine_tajimas_d.out
        plots = plot_selection_scans.out
}
