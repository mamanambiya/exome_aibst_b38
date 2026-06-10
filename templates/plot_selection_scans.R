#!/usr/bin/env Rscript

# Selection Scan Plotting Script
# Generates Manhattan plots for iHS and Tajima's D
# Usage: Rscript plot_selection_scans.R <population> <ihs_file> <tajd_file> <pgx_bed> <output_dir>

library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    stop("Usage: Rscript plot_selection_scans.R <population> <ihs_file> <tajd_file> <pgx_bed> <output_dir>")
}

population <- args[1]
ihs_file <- args[2]
tajd_file <- args[3]
pgx_bed <- args[4]
output_dir <- args[5]

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat(paste("Processing selection scans for", population, "\n"))

# Read pharmacogene BED file
read_pgx_genes <- function(bed_file) {
    if (!file.exists(bed_file)) {
        warning(paste("Pharmacogene BED file not found:", bed_file))
        return(NULL)
    }
    pgx <- fread(bed_file, header = FALSE)
    colnames(pgx) <- c("chr", "start", "end", "gene", "score", "strand")
    pgx$chr <- as.numeric(gsub("chr", "", pgx$chr))
    return(pgx)
}

pgx_genes <- read_pgx_genes(pgx_bed)

# Function to calculate cumulative positions for Manhattan plot
calc_cumulative_pos <- function(data) {
    # Calculate chromosome lengths
    chr_lengths <- data %>%
        group_by(chr) %>%
        summarise(max_pos = max(pos), .groups = "drop") %>%
        arrange(chr) %>%
        mutate(
            cumulative = cumsum(as.numeric(max_pos)) - max_pos,
            chr_center = cumulative + max_pos / 2
        )

    # Add cumulative positions to data
    data <- data %>%
        left_join(chr_lengths %>% select(chr, cumulative), by = "chr") %>%
        mutate(pos_cum = pos + cumulative)

    return(list(data = data, chr_info = chr_lengths))
}

# Function to plot Manhattan plot
plot_manhattan <- function(data, chr_info, title, y_label, threshold = NULL,
                            pgx_genes = NULL, output_file) {
    # Color chromosomes alternately
    data$color <- ifelse(data$chr %% 2 == 0, "even", "odd")

    # Base plot
    p <- ggplot(data, aes(x = pos_cum, y = value, color = color)) +
        geom_point(alpha = 0.6, size = 0.8) +
        scale_color_manual(values = c("even" = "#4292C6", "odd" = "#084594")) +
        scale_x_continuous(
            breaks = chr_info$chr_center,
            labels = chr_info$chr
        ) +
        labs(
            title = title,
            x = "Chromosome",
            y = y_label
        ) +
        theme_minimal() +
        theme(
            legend.position = "none",
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.text.x = element_text(angle = 0, hjust = 0.5),
            plot.title = element_text(hjust = 0.5, face = "bold")
        )

    # Add threshold lines if provided
    if (!is.null(threshold)) {
        p <- p +
            geom_hline(yintercept = threshold, linetype = "dashed", color = "red", size = 0.5) +
            geom_hline(yintercept = -threshold, linetype = "dashed", color = "red", size = 0.5)
    }

    # Add pharmacogene positions if provided
    if (!is.null(pgx_genes)) {
        # Calculate cumulative positions for genes
        pgx_with_pos <- pgx_genes %>%
            left_join(chr_info %>% select(chr, cumulative), by = "chr") %>%
            mutate(
                gene_start_cum = start + cumulative,
                gene_end_cum = end + cumulative,
                gene_mid_cum = (gene_start_cum + gene_end_cum) / 2
            )

        # Add vertical lines for pharmacogenes
        p <- p +
            geom_segment(
                data = pgx_with_pos,
                aes(x = gene_mid_cum, xend = gene_mid_cum, y = min(data$value), yend = max(data$value)),
                color = "darkgreen", alpha = 0.3, size = 0.3, inherit.aes = FALSE
            )
    }

    ggsave(output_file, plot = p, width = 14, height = 6, dpi = 300)
    cat(paste("Saved:", output_file, "\n"))
}

# Process iHS results
if (file.exists(ihs_file)) {
    cat("\nProcessing iHS results...\n")

    # Read iHS data (normalized)
    ihs_data <- fread(ihs_file)
    colnames(ihs_data) <- c("snp_id", "chr", "pos", "freq", "ihs", "ihs_norm")

    # Calculate cumulative positions
    ihs_plot_data <- ihs_data %>%
        select(chr, pos, value = ihs_norm) %>%
        filter(!is.na(value) & !is.infinite(value))

    cumulative_result <- calc_cumulative_pos(ihs_plot_data)
    ihs_plot_data <- cumulative_result$data
    chr_info <- cumulative_result$chr_info

    # Identify outliers
    ihs_threshold <- 2.5
    outliers <- ihs_plot_data %>%
        filter(abs(value) > ihs_threshold) %>%
        arrange(desc(abs(value)))

    cat(paste("Found", nrow(outliers), "iHS outliers (|iHS| >", ihs_threshold, ")\n"))

    # Save outliers
    outlier_file <- file.path(output_dir, paste0(population, "_ihs_outliers.txt"))
    write.table(outliers %>% select(chr, pos, value),
                outlier_file, quote = FALSE, row.names = FALSE, sep = "\t")
    cat(paste("Saved outliers:", outlier_file, "\n"))

    # Check for pharmacogene overlaps
    if (!is.null(pgx_genes)) {
        pgx_overlaps <- outliers %>%
            rowwise() %>%
            filter(any(pgx_genes$chr == chr &
                       pgx_genes$start <= pos &
                       pgx_genes$end >= pos)) %>%
            ungroup()

        if (nrow(pgx_overlaps) > 0) {
            pgx_overlap_file <- file.path(output_dir, paste0(population, "_ihs_pgx_overlap.txt"))
            write.table(pgx_overlaps %>% select(chr, pos, value),
                        pgx_overlap_file, quote = FALSE, row.names = FALSE, sep = "\t")
            cat(paste("Found", nrow(pgx_overlaps), "iHS outliers overlapping pharmacogenes\n"))
            cat(paste("Saved overlaps:", pgx_overlap_file, "\n"))
        }
    }

    # Plot iHS Manhattan
    ihs_plot_file <- file.path(output_dir, paste0(population, "_ihs_manhattan.pdf"))
    plot_manhattan(
        ihs_plot_data,
        chr_info,
        title = paste0(population, " - iHS (Integrated Haplotype Score)"),
        y_label = "Normalized iHS",
        threshold = ihs_threshold,
        pgx_genes = pgx_genes,
        output_file = ihs_plot_file
    )
} else {
    cat("iHS file not found, skipping iHS plot\n")
}

# Process Tajima's D results
if (file.exists(tajd_file)) {
    cat("\nProcessing Tajima's D results...\n")

    # Read Tajima's D data
    tajd_data <- fread(tajd_file)
    colnames(tajd_data) <- c("chr", "bin_start", "n_snps", "tajd")

    # Calculate bin midpoint for plotting
    tajd_plot_data <- tajd_data %>%
        mutate(pos = bin_start + 25000) %>%  # Midpoint of 50kb window
        select(chr, pos, value = tajd) %>%
        filter(!is.na(value) & !is.infinite(value))

    # Calculate cumulative positions
    cumulative_result <- calc_cumulative_pos(tajd_plot_data)
    tajd_plot_data <- cumulative_result$data
    chr_info <- cumulative_result$chr_info

    # Identify extreme Tajima's D values
    tajd_threshold <- 2.0
    extreme_tajd <- tajd_plot_data %>%
        filter(abs(value) > tajd_threshold) %>%
        arrange(desc(abs(value)))

    cat(paste("Found", nrow(extreme_tajd), "extreme Tajima's D windows (|D| >", tajd_threshold, ")\n"))

    # Save extreme values
    extreme_file <- file.path(output_dir, paste0(population, "_tajd_extreme.txt"))
    write.table(extreme_tajd %>% select(chr, pos, value),
                extreme_file, quote = FALSE, row.names = FALSE, sep = "\t")
    cat(paste("Saved extreme values:", extreme_file, "\n"))

    # Plot Tajima's D Manhattan
    tajd_plot_file <- file.path(output_dir, paste0(population, "_tajd_manhattan.pdf"))
    plot_manhattan(
        tajd_plot_data,
        chr_info,
        title = paste0(population, " - Tajima's D"),
        y_label = "Tajima's D",
        threshold = tajd_threshold,
        pgx_genes = pgx_genes,
        output_file = tajd_plot_file
    )
} else {
    cat("Tajima's D file not found, skipping Tajima's D plot\n")
}

cat("\nSelection scan plotting complete!\n")
