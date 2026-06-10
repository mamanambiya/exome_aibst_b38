#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
library(scales)

fst_data <- read.delim("${max_fst_file}", header = TRUE)

fst_data\$category <- ifelse(fst_data\$MAX_FST > 0.25,
    "HDV candidate\n(FST > 0.25)", "Non-HDV")

n_hdv <- sum(fst_data\$MAX_FST > 0.25)
n_total <- nrow(fst_data)
pct_hdv <- round(n_hdv / n_total * 100, 1)

p <- ggplot(fst_data, aes(x = MAX_FST, fill = category)) +
    geom_histogram(bins = 100, alpha = 0.85, boundary = 0,
        color = "grey30", linewidth = 0.1) +
    geom_vline(xintercept = 0.25, linetype = "dashed",
        color = "red", linewidth = 0.8) +
    annotate("text", x = 0.27, y = Inf,
        label = expression(paste(F[ST], " = 0.25")),
        hjust = 0, vjust = 1.5, size = 3, color = "red",
        fontface = "italic") +
    annotate("text", x = 0.75, y = Inf,
        label = paste0("n = ", comma(n_hdv), "\n(", pct_hdv, "%)"),
        hjust = 0.5, vjust = 2, size = 3.5, color = "#d62728",
        fontface = "bold") +
    annotate("text", x = 0.05, y = Inf,
        label = paste0("n = ", comma(n_total - n_hdv), "\n(",
            round((n_total - n_hdv) / n_total * 100, 1), "%)"),
        hjust = 0.5, vjust = 2, size = 3.5, color = "#1f77b4",
        fontface = "bold") +
    scale_fill_manual(values = c(
        "Non-HDV" = "#1f77b4",
        "HDV candidate\n(FST > 0.25)" = "#d62728")) +
    scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(-0.05, 1.05)) +
    scale_y_continuous(labels = comma) +
    labs(
        x = expression(paste("Maximum pairwise ", F[ST],
            " across 12 AiBST populations")),
        y = "Number of pharmacogene variants",
        fill = NULL,
        title = expression(paste(F[ST],
            " distribution of pharmacogene variants")),
        subtitle = "Dashed line indicates HDV identification threshold"
    ) +
    theme_minimal(base_size = 11) +
    theme(
        legend.position = c(0.75, 0.6),
        legend.background = element_rect(fill = "white",
            color = "grey80", linewidth = 0.3),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey50"),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "grey30", linewidth = 0.3)
    )

ggsave("fst_pgx_distribution.pdf", p, width = 7, height = 4.5, dpi = 300)
ggsave("fst_pgx_distribution.png", p, width = 7, height = 4.5, dpi = 300)
