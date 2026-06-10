#!/usr/bin/env Rscript

# ADMIXTURE Plotting Script
# Generates barplots for K=2 to K=max and CV error plot
# Usage: Rscript plot_admixture.R <dataset> <k_min> <k_max> <prefix>

library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    stop("Usage: Rscript plot_admixture.R <dataset> <k_min> <k_max> <prefix>")
}

dataset <- args[1]
k_min <- as.numeric(args[2])
k_max <- as.numeric(args[3])
prefix <- args[4]

# Read population file (.fam or .ped)
fam_file <- paste0(prefix, ".fam")
if (!file.exists(fam_file)) {
    fam_file <- paste0(prefix, ".ped")
}

fam <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)
colnames(fam)[1:2] <- c("pop", "ind")

# Function to read Q file for a given K
read_Q <- function(k, prefix) {
    q_file <- paste0(prefix, ".", k, ".Q")
    if (!file.exists(q_file)) {
        warning(paste("Q file not found:", q_file))
        return(NULL)
    }
    q <- read.table(q_file, header = FALSE)
    colnames(q) <- paste0("Anc", 1:k)
    q$ind <- fam$ind
    q$pop <- fam$pop
    q$K <- k
    return(q)
}

# Function to plot ADMIXTURE barplot for a given K
plot_admixture_barplot <- function(q_data, k, output_file) {
    # Reshape data for ggplot
    q_long <- q_data %>%
        select(-K) %>%
        gather(key = "Ancestry", value = "Proportion", starts_with("Anc"))

    # Order individuals by population
    q_long$ind <- factor(q_long$ind, levels = unique(q_long$ind))

    # Create color palette
    if (k <= 8) {
        colors <- brewer.pal(k, "Set2")
    } else {
        colors <- colorRampPalette(brewer.pal(8, "Set2"))(k)
    }

    # Create plot
    p <- ggplot(q_long, aes(x = ind, y = Proportion, fill = Ancestry)) +
        geom_bar(stat = "identity", width = 1) +
        facet_grid(~ pop, scales = "free_x", space = "free_x") +
        scale_fill_manual(values = colors) +
        labs(
            title = paste0(dataset, " - ADMIXTURE K=", k),
            x = "Individual",
            y = "Ancestry Proportion"
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.spacing = unit(0.1, "lines"),
            strip.text = element_text(angle = 90, hjust = 0),
            legend.position = "bottom"
        )

    ggsave(output_file, plot = p, width = 14, height = 6, dpi = 300)
    cat(paste("Saved:", output_file, "\n"))
}

# Function to read CV error
read_cv_error <- function(k, prefix) {
    log_file <- paste0(prefix, ".", k, ".log")
    if (!file.exists(log_file)) {
        warning(paste("Log file not found:", log_file))
        return(NA)
    }
    log_lines <- readLines(log_file)
    cv_line <- grep("CV error", log_lines, value = TRUE)
    if (length(cv_line) == 0) {
        return(NA)
    }
    cv_error <- as.numeric(sub(".*CV error \\(K=\\d+\\): ", "", cv_line[1]))
    return(cv_error)
}

# Main execution
cat(paste("Processing ADMIXTURE results for", dataset, "\n"))
cat(paste("K range:", k_min, "to", k_max, "\n"))

# Create output directory for plots
plot_dir <- paste0(dirname(prefix), "/plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Process each K value
cv_errors <- data.frame(K = k_min:k_max, CV_Error = NA)

for (k in k_min:k_max) {
    cat(paste("\nProcessing K =", k, "\n"))

    # Read Q file
    q_data <- read_Q(k, prefix)
    if (is.null(q_data)) {
        next
    }

    # Plot barplot
    output_file <- file.path(plot_dir, paste0(dataset, "_barplot_K", k, ".pdf"))
    plot_admixture_barplot(q_data, k, output_file)

    # Read CV error
    cv_error <- read_cv_error(k, prefix)
    cv_errors$CV_Error[cv_errors$K == k] <- cv_error
    cat(paste("CV error (K=", k, "):", cv_error, "\n"))
}

# Save CV errors
cv_file <- paste0(prefix, ".cv_errors_all_K.txt")
write.table(cv_errors, cv_file, quote = FALSE, row.names = FALSE, sep = "\t")
cat(paste("\nSaved CV errors:", cv_file, "\n"))

# Plot CV errors
cv_errors_clean <- cv_errors[!is.na(cv_errors$CV_Error), ]
if (nrow(cv_errors_clean) > 0) {
    optimal_k <- cv_errors_clean$K[which.min(cv_errors_clean$CV_Error)]
    min_cv <- min(cv_errors_clean$CV_Error)

    p_cv <- ggplot(cv_errors_clean, aes(x = K, y = CV_Error)) +
        geom_line(color = "blue", size = 1) +
        geom_point(color = "blue", size = 3) +
        geom_point(data = cv_errors_clean[cv_errors_clean$K == optimal_k, ],
                   aes(x = K, y = CV_Error), color = "red", size = 5) +
        labs(
            title = paste0(dataset, " - ADMIXTURE Cross-Validation Errors"),
            subtitle = paste0("Optimal K = ", optimal_k, " (CV error = ", round(min_cv, 4), ")"),
            x = "K (Number of Ancestral Populations)",
            y = "Cross-Validation Error"
        ) +
        scale_x_continuous(breaks = k_min:k_max) +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, color = "red")
        )

    cv_plot_file <- file.path(plot_dir, paste0(dataset, "_cv_error_plot.pdf"))
    ggsave(cv_plot_file, plot = p_cv, width = 10, height = 6, dpi = 300)
    cat(paste("\nSaved CV error plot:", cv_plot_file, "\n"))

    # Write optimal K to file
    optimal_file <- paste0(prefix, ".optimal_K.txt")
    writeLines(
        c(
            paste0("Dataset: ", dataset),
            paste0("Optimal K: ", optimal_k),
            paste0("CV Error: ", round(min_cv, 4)),
            "",
            "All K values and CV errors:",
            paste(cv_errors_clean$K, cv_errors_clean$CV_Error, sep = "\t")
        ),
        optimal_file
    )
    cat(paste("Saved optimal K:", optimal_file, "\n"))
    cat(paste("\n*** OPTIMAL K =", optimal_k, "***\n\n"))
} else {
    cat("\nWarning: No valid CV errors found\n")
}

cat("ADMIXTURE plotting complete!\n")
