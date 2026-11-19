#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(DESeq2)
    library(dplyr)
    library(ggplot2)
    library(ggrepel)
    library(patchwork)
    library(pheatmap)
    library(readr)
    library(ragg)
    library(scales)
    library(tibble)
    library(tidyr)
})

options(readr.show_col_types = FALSE)

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(flag, default) {
    index <- match(flag, args)
    if (is.na(index)) {
        return(default)
    }
    if (index == length(args)) {
        stop("Missing value after ", flag)
    }
    args[[index + 1]]
}

results_dir <- arg_value("--results", "data/results/peak_analysis")
metadata_path <- arg_value(
    "--metadata", "code/clip_peak_analysis/deseq_samples.tsv"
)
figure_dir <- arg_value("--output", file.path(results_dir, "figures"))

alpha <- 0.05
effect_threshold <- 1
pseudocount <- 0.5
vst_minimum_total <- 10
heatmap_peak_count <- 75

required_files <- c(
    file.path(results_dir, "peak_count_matrix.tsv"),
    file.path(results_dir, "peaks.filtered_merged.tsv"),
    file.path(results_dir, "bam_library_sizes.tsv"),
    metadata_path
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
    stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

directories <- c(
    file.path(figure_dir, "QC"),
    file.path(figure_dir, "MSI2_328WT"),
    file.path(figure_dir, "MSI2_324WT"),
    file.path(figure_dir, "construct_comparison"),
    file.path(figure_dir, "peak_examples")
)
invisible(lapply(directories, dir.create, recursive = TRUE, showWarnings = FALSE))

theme_publication <- function(base_size = 10) {
    theme_classic(base_size = base_size) +
        theme(
            plot.title = element_text(face = "bold", size = rel(1.15)),
            plot.subtitle = element_text(color = "grey30"),
            plot.caption = element_text(color = "grey35", hjust = 0),
            axis.title = element_text(face = "bold"),
            strip.background = element_rect(fill = "grey94", color = NA),
            strip.text = element_text(face = "bold"),
            legend.position = "bottom",
            legend.title = element_text(face = "bold"),
            panel.grid.major.y = element_line(color = "grey90", linewidth = 0.25),
            panel.grid.minor = element_blank()
        )
}

assay_colors <- c(input = "#0072B2", IP = "#D55E00")
construct_shapes <- c(MSI2_328WT = 16, MSI2_324WT = 17)
pair_colors <- c(A = "#0072B2", B = "#D55E00")
significance_colors <- c(
    "Not tested (low count)" = "#D0D0D0",
    "Tested, not significant" = "#6B6B6B",
    "Significantly enriched" = "#D55E00",
    "Significantly depleted" = "#0072B2"
)
comparison_colors <- c(
    "Enriched in both" = "#7B3294",
    "Enriched only in 328WT" = "#E66101",
    "Enriched only in 324WT" = "#1B9E77",
    "Not significantly enriched" = "#B8B8B8"
)

save_plot <- function(plot, stem, width, height) {
    ggsave(
        paste0(stem, ".pdf"), plot = plot,
        device = grDevices::pdf, width = width, height = height,
        units = "in", bg = "white"
    )
    ggsave(
        paste0(stem, ".png"), plot = plot,
        device = ragg::agg_png, width = width, height = height,
        units = "in", dpi = 300, bg = "white"
    )
}

save_heatmap <- function(heatmap, stem, width, height) {
    grDevices::pdf(paste0(stem, ".pdf"), width = width, height = height)
    grid::grid.newpage()
    grid::grid.draw(heatmap$gtable)
    grDevices::dev.off()

    ragg::agg_png(
        paste0(stem, ".png"), width = width, height = height,
        units = "in", res = 300, background = "white"
    )
    grid::grid.newpage()
    grid::grid.draw(heatmap$gtable)
    grDevices::dev.off()
}

write_table <- function(x, path) {
    readr::write_tsv(as.data.frame(x), path, na = "NA")
}

metadata <- read_tsv(metadata_path) %>%
    mutate(
        assay = factor(assay, levels = c("input", "IP")),
        construct = factor(construct, levels = unique(construct)),
        comparison = factor(comparison, levels = c("MSI2_328WT", "MSI2_324WT")),
        pair = factor(pair, levels = c("A", "B")),
        short_sample = sub("NM2026_", "", sample)
    )

counts <- read_tsv(file.path(results_dir, "peak_count_matrix.tsv"))
peaks <- read_tsv(file.path(results_dir, "peaks.filtered_merged.tsv"))
library_sizes <- read_tsv(file.path(results_dir, "bam_library_sizes.tsv"))

stopifnot(
    identical(names(counts)[-1], metadata$sample),
    identical(counts$peak_id, peaks$peak_id),
    setequal(library_sizes$sample, metadata$sample)
)

count_matrix <- as.matrix(counts[, metadata$sample])
storage.mode(count_matrix) <- "integer"
rownames(count_matrix) <- counts$peak_id
if (any(count_matrix < 0L) || any(count_matrix != floor(count_matrix))) {
    stop("The peak count matrix must contain non-negative integers")
}

library_sizes <- library_sizes %>%
    slice(match(metadata$sample, sample))
full_library_factors <- library_sizes$primary_mapped_alignments /
    exp(mean(log(library_sizes$primary_mapped_alignments)))
names(full_library_factors) <- library_sizes$sample

# Sample-level VST uses informative peaks only. This affects visualization, not
# either existing DESeq2 fit or its unfiltered normalized-count output.
vst_keep <- rowSums(count_matrix) >= vst_minimum_total
vst_metadata <- as.data.frame(metadata)
rownames(vst_metadata) <- vst_metadata$sample
vst_metadata$construct <- factor(vst_metadata$comparison)
vst_metadata$assay <- factor(vst_metadata$assay)
dds_all <- DESeqDataSetFromMatrix(
    countData = count_matrix[vst_keep, , drop = FALSE],
    colData = vst_metadata,
    design = ~ construct + assay
)
sizeFactors(dds_all) <- full_library_factors[colnames(dds_all)]
vsd_all <- varianceStabilizingTransformation(dds_all, blind = TRUE)
vst_matrix <- assay(vsd_all)

write_table(
    data.frame(peak_id = rownames(vst_matrix), vst_matrix, check.names = FALSE),
    file.path(figure_dir, "QC", "VST_counts_for_QC.tsv")
)

parameter_table <- tibble(
    parameter = c(
        "adjusted_p_value_threshold", "effect_size_reference_threshold",
        "pair_ratio_pseudocount", "VST_minimum_total_count",
        "heatmap_target_peak_count", "omniCLIP_score_metric"
    ),
    value = c(alpha, effect_threshold, pseudocount, vst_minimum_total,
              heatmap_peak_count, "peak_score_max")
)
write_table(parameter_table, file.path(figure_dir, "analysis_parameters.tsv"))

# -----------------------------------------------------------------------------
# 1. Sample-level QC
# -----------------------------------------------------------------------------

pca <- prcomp(t(vst_matrix))
pca_variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_data <- as.data.frame(pca$x[, 1:2, drop = FALSE]) %>%
    rownames_to_column("sample") %>%
    left_join(metadata, by = "sample")
write_table(pca_data, file.path(figure_dir, "QC", "PCA_all_samples.tsv"))

pca_plot <- ggplot(
    pca_data,
    aes(PC1, PC2, color = assay, shape = comparison, label = short_sample)
) +
    geom_hline(yintercept = 0, color = "grey88", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey88", linewidth = 0.3) +
    geom_point(size = 3.2) +
    ggrepel::geom_text_repel(
        size = 3, color = "grey15", seed = 328324,
        box.padding = 0.45, point.padding = 0.3, max.overlaps = Inf
    ) +
    scale_color_manual(values = assay_colors, drop = FALSE) +
    scale_shape_manual(
        values = construct_shapes,
        labels = c(MSI2_328WT = "328WT", MSI2_324WT = "324WT"),
        drop = FALSE
    ) +
    labs(
        title = "All-sample VST peak-count PCA",
        subtitle = paste0(
            sum(vst_keep), " peaks with total count >= ", vst_minimum_total,
            "; full mapped-library size factors"
        ),
        x = sprintf("PC1 (%.1f%%)", pca_variance[1]),
        y = sprintf("PC2 (%.1f%%)", pca_variance[2]),
        color = "Assay", shape = "Construct"
    ) +
    theme_publication()
save_plot(
    pca_plot,
    file.path(figure_dir, "QC", "PCA_all_samples"),
    width = 7.4, height = 5.6
)

sample_cor <- cor(vst_matrix, method = "pearson")
write_table(
    data.frame(sample = rownames(sample_cor), sample_cor, check.names = FALSE),
    file.path(figure_dir, "QC", "sample_correlation_matrix.tsv")
)

annotation_col <- metadata %>%
    transmute(
        sample,
        Construct = factor(comparison, levels = c("MSI2_328WT", "MSI2_324WT")),
        Assay = assay,
        Pair = pair
    ) %>%
    column_to_rownames("sample")
annotation_colors <- list(
    Construct = c(MSI2_328WT = "#7B3294", MSI2_324WT = "#1B9E77"),
    Assay = assay_colors,
    Pair = pair_colors
)
cor_min <- floor(min(sample_cor) * 20) / 20
cor_breaks <- seq(cor_min, 1, length.out = 101)
cor_heatmap <- pheatmap(
    sample_cor,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    breaks = cor_breaks,
    annotation_col = annotation_col,
    annotation_row = annotation_col,
    annotation_colors = annotation_colors,
    border_color = NA,
    fontsize = 9,
    main = "Pearson correlation of VST peak counts",
    silent = TRUE
)
save_heatmap(
    cor_heatmap,
    file.path(figure_dir, "QC", "sample_correlation_heatmap"),
    width = 8.2, height = 7.4
)

size_factor_table <- bind_rows(lapply(
    c("MSI2_328WT", "MSI2_324WT"),
    function(comparison_name) {
        read_tsv(file.path(
            results_dir, "deseq2", comparison_name, "size_factors.tsv"
        ))
    }
))
sample_summary <- metadata %>%
    left_join(library_sizes, by = "sample") %>%
    left_join(size_factor_table %>% select(sample, size_factor, normalization),
              by = "sample") %>%
    mutate(
        assigned_to_fixed_peaks = as.numeric(colSums(count_matrix)[sample]),
        assigned_percent = 100 * assigned_to_fixed_peaks /
            primary_mapped_alignments,
        nonzero_peaks = as.numeric(colSums(count_matrix > 0)[sample])
    )
write_table(
    sample_summary,
    file.path(figure_dir, "QC", "sample_count_summary.tsv")
)

sample_summary_long <- sample_summary %>%
    select(sample, short_sample, assay, comparison,
           primary_mapped_alignments, assigned_to_fixed_peaks,
           assigned_percent, nonzero_peaks, size_factor) %>%
    pivot_longer(
        cols = c(primary_mapped_alignments, assigned_to_fixed_peaks,
                 assigned_percent, nonzero_peaks, size_factor),
        names_to = "metric", values_to = "value"
    ) %>%
    mutate(
        metric = factor(
            metric,
            levels = c(
                "primary_mapped_alignments", "assigned_to_fixed_peaks",
                "assigned_percent", "nonzero_peaks", "size_factor"
            ),
            labels = c(
                "Primary mapped alignments", "Counts assigned to fixed peaks",
                "Assigned to fixed peaks (%)", "Peaks with nonzero count",
                "DESeq2 size factor"
            )
        )
    )

count_summary_plot <- ggplot(
    sample_summary_long,
    aes(short_sample, value, fill = assay)
) +
    geom_col(width = 0.72) +
    facet_wrap(
        ~ metric,
        scales = "free_y",
        ncol = 3,
        axes = "all_x",
        axis.labels = "all_x"
    ) +
    scale_fill_manual(values = assay_colors, drop = FALSE) +
    scale_y_continuous(labels = label_number(big.mark = ","), expand = expansion(c(0, 0.08))) +
    labs(
        title = "Library and fixed-peak count summary",
        subtitle = "Sample labels omit the common NM2026_ prefix",
        x = "Sample", y = NULL, fill = "Assay"
    ) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(
    count_summary_plot,
    file.path(figure_dir, "QC", "library_count_summary"),
    width = 10.5, height = 6.8
)

# -----------------------------------------------------------------------------
# 2-4. Per-construct DESeq2, heatmaps, and pair reproducibility
# -----------------------------------------------------------------------------

build_construct <- function(comparison_name) {
    output_dir <- file.path(figure_dir, comparison_name)
    comparison_metadata <- metadata %>%
        filter(as.character(comparison) == comparison_name)
    sample_order <- comparison_metadata %>%
        arrange(pair, assay) %>%
        pull(sample)

    result <- read_tsv(file.path(
        results_dir, "deseq2", comparison_name, "IP_vs_input.all.tsv"
    ))
    normalized <- read_tsv(file.path(
        results_dir, "deseq2", comparison_name, "normalized_counts.tsv"
    ))
    stopifnot(all(sample_order %in% names(normalized)))

    pair_ratio_data <- normalized %>%
        select(peak_id) %>%
        left_join(peaks, by = "peak_id")
    for (pair_name in c("A", "B")) {
        ip_sample <- comparison_metadata %>%
            filter(pair == pair_name, assay == "IP") %>%
            pull(sample)
        input_sample <- comparison_metadata %>%
            filter(pair == pair_name, assay == "input") %>%
            pull(sample)
        if (length(ip_sample) != 1 || length(input_sample) != 1) {
            stop(comparison_name, " pair ", pair_name,
                 " is not one IP matched to one input")
        }
        pair_ratio_data[[paste0("pair_", pair_name, "_log2_ratio")]] <- log2(
            (normalized[[ip_sample]] + pseudocount) /
                (normalized[[input_sample]] + pseudocount)
        )
    }

    input_samples <- comparison_metadata %>%
        filter(assay == "input") %>% pull(sample)
    ip_samples <- comparison_metadata %>%
        filter(assay == "IP") %>% pull(sample)
    pair_ratio_data$mean_normalized_input <- rowMeans(
        as.matrix(normalized[, input_samples, drop = FALSE])
    )
    pair_ratio_data$mean_normalized_IP <- rowMeans(
        as.matrix(normalized[, ip_samples, drop = FALSE])
    )
    pair_ratio_data <- pair_ratio_data %>%
        left_join(
            result %>% select(
                peak_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
            ),
            by = "peak_id"
        ) %>%
        mutate(
            tested = peak_id %in% result$peak_id,
            significance = case_when(
                !tested ~ "Not tested (low count)",
                !is.na(padj) & padj < alpha & log2FoldChange > 0 ~
                    "Significantly enriched",
                !is.na(padj) & padj < alpha & log2FoldChange < 0 ~
                    "Significantly depleted",
                TRUE ~ "Tested, not significant"
            ),
            significance = factor(
                significance,
                levels = names(significance_colors)
            )
        )

    plot_data <- result %>%
        left_join(
            pair_ratio_data %>% select(
                peak_id, pair_A_log2_ratio, pair_B_log2_ratio,
                mean_normalized_input, mean_normalized_IP
            ),
            by = "peak_id"
        ) %>%
        mutate(
            significance = case_when(
                !is.na(padj) & padj < alpha & log2FoldChange > 0 ~
                    "Significantly enriched",
                !is.na(padj) & padj < alpha & log2FoldChange < 0 ~
                    "Significantly depleted",
                TRUE ~ "Tested, not significant"
            ),
            significance = factor(
                significance,
                levels = names(significance_colors)
            ),
            padj_for_plot = ifelse(is.na(padj), 1, pmax(padj, .Machine$double.xmin)),
            effect_at_least_one = abs(log2FoldChange) >= effect_threshold
        )
    write_table(plot_data, file.path(output_dir, "deseq_plot_data.tsv"))

    construct_label <- sub("MSI2_", "", comparison_name)
    ma_plot <- ggplot(plot_data, aes(baseMean, log2FoldChange, color = significance)) +
        geom_hline(yintercept = 0, color = "grey35", linewidth = 0.4) +
        geom_point(size = 1.15, alpha = 0.72) +
        scale_x_log10(labels = label_number()) +
        scale_color_manual(values = significance_colors, drop = TRUE) +
        labs(
            title = paste0(construct_label, " IP versus input MA plot"),
            subtitle = paste0(
                nrow(plot_data), " tested peaks; significance is padj < ", alpha
            ),
            x = "Mean normalized count (baseMean; log10 scale)",
            y = "IP/input log2 fold change", color = NULL
        ) +
        theme_publication()
    save_plot(
        ma_plot, file.path(output_dir, "MA_IP_vs_input"),
        width = 7.2, height = 5.6
    )

    volcano_labels <- plot_data %>%
        filter(!is.na(padj), padj < alpha, mean_normalized_input >= 1) %>%
        arrange(padj, desc(abs(log2FoldChange))) %>%
        slice_head(n = 8)
    volcano_plot <- ggplot(
        plot_data,
        aes(log2FoldChange, -log10(padj_for_plot), color = significance)
    ) +
        geom_hline(
            yintercept = -log10(alpha), color = "grey45",
            linetype = "dashed", linewidth = 0.4
        ) +
        geom_vline(
            xintercept = c(-effect_threshold, effect_threshold),
            color = "grey72", linetype = "dotted", linewidth = 0.35
        ) +
        geom_point(size = 1.15, alpha = 0.72) +
        ggrepel::geom_text_repel(
            data = volcano_labels,
            aes(label = peak_id),
            size = 2.6, color = "grey15", seed = 324328,
            min.segment.length = 0, max.overlaps = Inf,
            box.padding = 0.35, point.padding = 0.2
        ) +
        scale_color_manual(values = significance_colors, drop = TRUE) +
        labs(
            title = paste0(construct_label, " IP versus input volcano plot"),
            subtitle = paste0(
                "Dotted lines mark |log2FC| = ", effect_threshold,
                "; NA padj values are retained at y = 0"
            ),
            x = "IP/input log2 fold change",
            y = expression(-log[10](adjusted~italic(p))), color = NULL
        ) +
        theme_publication()
    save_plot(
        volcano_plot, file.path(output_dir, "volcano_IP_vs_input"),
        width = 7.2, height = 5.8
    )

    candidates <- pair_ratio_data %>%
        filter(
            tested, !is.na(padj), padj < alpha, log2FoldChange > 0,
            mean_normalized_input >= 1,
            pair_A_log2_ratio > 0, pair_B_log2_ratio > 0
        ) %>%
        mutate(
            rank_padj = rank(padj, ties.method = "min"),
            rank_effect = rank(-log2FoldChange, ties.method = "min"),
            selection_rank = rank_padj + rank_effect
        ) %>%
        arrange(selection_rank, padj, desc(log2FoldChange))
    if (nrow(candidates) < 25) {
        stop(comparison_name, " has too few reproducibly enriched heatmap candidates")
    }
    selected_heatmap_count <- min(heatmap_peak_count, nrow(candidates))
    selected_heatmap <- candidates %>%
        slice_head(n = selected_heatmap_count) %>%
        mutate(
            selection_rule = paste0(
                "padj < ", alpha,
                "; positive LFC and both pair ratios; mean input >= 1; ",
                "ranked by padj rank + LFC rank"
            )
        )
    write_table(selected_heatmap, file.path(output_dir, "heatmap_peaks.tsv"))

    heatmap_matrix <- normalized %>%
        filter(peak_id %in% selected_heatmap$peak_id) %>%
        slice(match(selected_heatmap$peak_id, peak_id)) %>%
        select(all_of(sample_order)) %>%
        as.matrix()
    rownames(heatmap_matrix) <- selected_heatmap$peak_id
    heatmap_log <- log2(heatmap_matrix + pseudocount)
    heatmap_scaled <- t(scale(t(heatmap_log)))
    heatmap_scaled[!is.finite(heatmap_scaled)] <- 0
    write_table(
        data.frame(peak_id = rownames(heatmap_scaled), heatmap_scaled,
                   check.names = FALSE),
        file.path(output_dir, "heatmap_values.row_z.tsv")
    )

    heatmap_annotation <- comparison_metadata %>%
        slice(match(sample_order, sample)) %>%
        transmute(sample, Assay = assay, Pair = pair) %>%
        column_to_rownames("sample")
    heatmap_col_labels <- comparison_metadata %>%
        slice(match(sample_order, sample)) %>%
        transmute(label = paste0(as.character(assay), " ", pair, "\n", short_sample)) %>%
        pull(label)
    enriched_heatmap <- pheatmap(
        heatmap_scaled,
        color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
        breaks = seq(-2.5, 2.5, length.out = 102),
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        show_rownames = FALSE,
        labels_col = heatmap_col_labels,
        annotation_col = heatmap_annotation,
        annotation_colors = annotation_colors[c("Assay", "Pair")],
        border_color = NA,
        fontsize = 9,
        main = paste0(
            construct_label, ": ", nrow(heatmap_scaled),
            " reproducibly enriched peaks (row z-score)"
        ),
        silent = TRUE
    )
    save_heatmap(
        enriched_heatmap,
        file.path(output_dir, "enriched_peak_heatmap"),
        width = 6.8, height = 8.2
    )

    tested_pair_data <- pair_ratio_data %>% filter(tested)
    pair_pearson <- cor(
        tested_pair_data$pair_A_log2_ratio,
        tested_pair_data$pair_B_log2_ratio,
        method = "pearson"
    )
    pair_spearman <- cor(
        tested_pair_data$pair_A_log2_ratio,
        tested_pair_data$pair_B_log2_ratio,
        method = "spearman"
    )
    pair_ratio_output <- pair_ratio_data %>%
        mutate(
            pair_ratio_pseudocount = pseudocount,
            correlation_population = ifelse(
                tested,
                "included in reported correlations",
                "shown but excluded from reported correlations"
            )
        )
    write_table(
        pair_ratio_output,
        file.path(output_dir, "pair_reproducibility.tsv")
    )

    pair_plot_data <- pair_ratio_data %>%
        arrange(significance)
    axis_limits <- range(c(
        pair_plot_data$pair_A_log2_ratio,
        pair_plot_data$pair_B_log2_ratio
    ), finite = TRUE)
    pair_plot <- ggplot(
        pair_plot_data,
        aes(pair_A_log2_ratio, pair_B_log2_ratio, color = significance)
    ) +
        geom_abline(slope = 1, intercept = 0, color = "grey55",
                    linetype = "dashed", linewidth = 0.4) +
        geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35) +
        geom_vline(xintercept = 0, color = "grey45", linewidth = 0.35) +
        geom_point(size = 0.7, alpha = 0.38) +
        coord_equal(xlim = axis_limits, ylim = axis_limits) +
        scale_color_manual(values = significance_colors, drop = TRUE) +
        annotate(
            "label", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.12,
            size = 3,
            label = sprintf(
                "Tested peaks: Pearson r = %.3f\nSpearman rho = %.3f; n = %s",
                pair_pearson, pair_spearman,
                scales::comma(nrow(tested_pair_data))
            ),
            linewidth = 0.2
        ) +
        labs(
            title = paste0(construct_label, " biological-pair reproducibility"),
            subtitle = paste0(
                "All ", scales::comma(nrow(pair_plot_data)),
                " peaks shown; pseudocount = ", pseudocount
            ),
            x = "Pair A log2(IP/input)",
            y = "Pair B log2(IP/input)", color = NULL,
            caption = "Descriptive pair agreement; this is not a replacement for the paired DESeq2 model."
        ) +
        theme_publication()
    save_plot(
        pair_plot, file.path(output_dir, "pair_reproducibility"),
        width = 7.2, height = 6.3
    )

    list(
        comparison = comparison_name,
        metadata = comparison_metadata,
        result = result,
        normalized = normalized,
        pair_data = pair_ratio_data,
        heatmap_peaks = selected_heatmap,
        pair_pearson = pair_pearson,
        pair_spearman = pair_spearman
    )
}

constructs <- list(
    MSI2_328WT = build_construct("MSI2_328WT"),
    MSI2_324WT = build_construct("MSI2_324WT")
)

# -----------------------------------------------------------------------------
# 5. Representative paired peak plots
# -----------------------------------------------------------------------------

select_peak_examples <- function(construct_data) {
    pool <- construct_data$pair_data %>%
        filter(tested) %>%
        mutate(
            minimum_pair_ratio = pmin(pair_A_log2_ratio, pair_B_log2_ratio),
            pair_disagreement = abs(pair_A_log2_ratio - pair_B_log2_ratio)
        )
    reproducible <- pool %>%
        filter(
            !is.na(padj), padj < alpha, log2FoldChange > 0,
            pair_A_log2_ratio > 0, pair_B_log2_ratio > 0,
            mean_normalized_input >= 1
        )

    selected <- tibble(peak_id = character())
    add_selection <- function(candidates, n, reason) {
        candidates <- candidates %>%
            filter(!peak_id %in% selected$peak_id) %>%
            slice_head(n = n) %>%
            mutate(selection_reason = reason)
        selected <<- bind_rows(selected, candidates)
    }

    add_selection(
        reproducible %>% arrange(desc(minimum_pair_ratio), padj),
        2, "Strong enrichment in both biological pairs"
    )
    add_selection(
        reproducible %>% arrange(padj, desc(log2FoldChange)),
        1, "High-confidence DESeq2 enrichment with nonzero input"
    )
    median_effect <- median(reproducible$log2FoldChange)
    add_selection(
        reproducible %>%
            mutate(distance_from_median_effect = abs(log2FoldChange - median_effect)) %>%
            arrange(distance_from_median_effect, padj),
        1, "Representative intermediate enriched peak"
    )
    add_selection(
        pool %>%
            filter(baseMean >= 5) %>%
            arrange(desc(pair_disagreement), desc(baseMean)),
        1, "Large pair-to-pair enrichment disagreement"
    )

    selected %>%
        slice_head(n = 5) %>%
        mutate(comparison = construct_data$comparison)
}

selected_examples <- bind_rows(lapply(constructs, select_peak_examples)) %>%
    select(
        comparison, peak_id, gene_ids, chromosome, start_0based,
        end_0based_exclusive, strand, peak_score_max, baseMean,
        log2FoldChange, padj, mean_normalized_input, mean_normalized_IP,
        pair_A_log2_ratio, pair_B_log2_ratio, selection_reason
    )
write_table(
    selected_examples,
    file.path(figure_dir, "peak_examples", "selected_peak_examples.tsv")
)

example_counts <- bind_rows(lapply(seq_len(nrow(selected_examples)), function(index) {
    selected <- selected_examples[index, ]
    construct_data <- constructs[[selected$comparison]]
    construct_data$normalized %>%
        filter(peak_id == selected$peak_id) %>%
        pivot_longer(-peak_id, names_to = "sample", values_to = "normalized_count") %>%
        inner_join(construct_data$metadata, by = "sample") %>%
        mutate(
            selection_reason = selected$selection_reason,
            panel_label = paste0(sub("MSI2_", "", selected$comparison),
                                 " | ", peak_id)
        )
})) %>%
    mutate(assay = factor(assay, levels = c("input", "IP")))
write_table(
    example_counts,
    file.path(figure_dir, "peak_examples", "paired_peak_counts.tsv")
)

paired_examples_plot <- ggplot(
    example_counts,
    aes(assay, normalized_count, group = pair, color = pair, shape = assay)
) +
    geom_line(linewidth = 0.65, alpha = 0.8) +
    geom_point(size = 2.4) +
    facet_wrap(~ panel_label, scales = "free_y", ncol = 5) +
    scale_color_manual(values = pair_colors, drop = FALSE) +
    scale_shape_manual(values = c(input = 1, IP = 16), drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.14))) +
    labs(
        title = "Representative peak counts preserve biological pairing",
        subtitle = "Five algorithmically selected examples per construct; free y-axis per peak",
        x = NULL, y = "Full-library normalized count",
        color = "Pair", shape = "Assay"
    ) +
    theme_publication(base_size = 9) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(
    paired_examples_plot,
    file.path(figure_dir, "peak_examples", "paired_peak_counts"),
    width = 12, height = 7.2
)

# -----------------------------------------------------------------------------
# 6-7. Direct construct comparison and descriptive enrichment categories
# -----------------------------------------------------------------------------

results_328 <- constructs$MSI2_328WT$result %>%
    select(
        peak_id, chromosome, start_0based, end_0based_exclusive, strand,
        gene_ids, peak_score_max, baseMean, log2FoldChange, padj
    ) %>%
    rename_with(~ paste0(.x, "_328"), -peak_id)
results_324 <- constructs$MSI2_324WT$result %>%
    select(peak_id, baseMean, log2FoldChange, padj) %>%
    rename_with(~ paste0(.x, "_324"), -peak_id)

construct_comparison <- inner_join(results_328, results_324, by = "peak_id") %>%
    mutate(
        enriched_328 = !is.na(padj_328) & padj_328 < alpha & log2FoldChange_328 > 0,
        enriched_324 = !is.na(padj_324) & padj_324 < alpha & log2FoldChange_324 > 0,
        enrichment_category = case_when(
            enriched_328 & enriched_324 ~ "Enriched in both",
            enriched_328 ~ "Enriched only in 328WT",
            enriched_324 ~ "Enriched only in 324WT",
            TRUE ~ "Not significantly enriched"
        ),
        enrichment_category = factor(
            enrichment_category,
            levels = names(comparison_colors)
        ),
        log2FC_difference_328_minus_324 =
            log2FoldChange_328 - log2FoldChange_324
    )
write_table(
    construct_comparison,
    file.path(
        figure_dir, "construct_comparison",
        "328_vs_324_peak_comparison.tsv"
    )
)

construct_pearson <- cor(
    construct_comparison$log2FoldChange_328,
    construct_comparison$log2FoldChange_324,
    method = "pearson"
)
construct_spearman <- cor(
    construct_comparison$log2FoldChange_328,
    construct_comparison$log2FoldChange_324,
    method = "spearman"
)
construct_labels <- construct_comparison %>%
    filter(baseMean_328 >= 5, baseMean_324 >= 5) %>%
    arrange(desc(abs(log2FC_difference_328_minus_324))) %>%
    slice_head(n = 8)
construct_limits <- range(c(
    construct_comparison$log2FoldChange_328,
    construct_comparison$log2FoldChange_324
), finite = TRUE)

construct_plot <- ggplot(
    construct_comparison,
    aes(log2FoldChange_328, log2FoldChange_324, color = enrichment_category)
) +
    geom_abline(slope = 1, intercept = 0, color = "grey55",
                linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35) +
    geom_vline(xintercept = 0, color = "grey45", linewidth = 0.35) +
    geom_point(size = 1.15, alpha = 0.68) +
    ggrepel::geom_text_repel(
        data = construct_labels,
        aes(label = peak_id),
        size = 2.6, color = "grey15", seed = 328324,
        min.segment.length = 0, max.overlaps = Inf,
        box.padding = 0.4, point.padding = 0.25
    ) +
    coord_equal(xlim = construct_limits, ylim = construct_limits) +
    scale_color_manual(values = comparison_colors, drop = FALSE) +
    annotate(
        "label", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.12,
        size = 3,
        label = sprintf(
            "Pearson r = %.3f\nSpearman rho = %.3f; n = %s",
            construct_pearson, construct_spearman,
            scales::comma(nrow(construct_comparison))
        ),
        linewidth = 0.2
    ) +
    labs(
        title = "IP enrichment effect sizes across MSI2 constructs",
        subtitle = "Only peaks passing the existing count filter in both separate DESeq2 analyses",
        x = "328WT IP/input log2 fold change",
        y = "324WT IP/input log2 fold change",
        color = NULL,
        caption = paste(
            "Categories are descriptive. Significance in one construct but not the other",
            "does not test a construct-by-assay difference."
        )
    ) +
    theme_publication()
save_plot(
    construct_plot,
    file.path(
        figure_dir, "construct_comparison",
        "328_vs_324_log2FC_scatter"
    ),
    width = 7.4, height = 6.5
)

category_summary <- construct_comparison %>%
    count(enrichment_category, .drop = FALSE, name = "peak_count") %>%
    mutate(percent_common_tested = 100 * peak_count / sum(peak_count))
write_table(
    category_summary,
    file.path(
        figure_dir, "construct_comparison", "enrichment_category_summary.tsv"
    )
)

category_plot <- ggplot(
    category_summary,
    aes(enrichment_category, peak_count, fill = enrichment_category)
) +
    geom_col(width = 0.68, show.legend = FALSE) +
    geom_text(
        aes(label = paste0(scales::comma(peak_count), " (",
                           sprintf("%.1f", percent_common_tested), "%)")),
        hjust = -0.08, size = 3.2
    ) +
    coord_flip(clip = "off") +
    scale_fill_manual(values = comparison_colors, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(
        title = "Descriptive enrichment categories across constructs",
        subtitle = paste0(
            scales::comma(nrow(construct_comparison)),
            " peaks tested in both separate DESeq2 analyses"
        ),
        x = NULL, y = "Number of peaks",
        caption = "Construct-only categories are not a formal test of construct differences."
    ) +
    theme_publication()
save_plot(
    category_plot,
    file.path(
        figure_dir, "construct_comparison", "enrichment_category_summary"
    ),
    width = 7.4, height = 4.8
)

# -----------------------------------------------------------------------------
# 8. Original omniCLIP score versus current enrichment
# -----------------------------------------------------------------------------

score_validation <- bind_rows(lapply(names(constructs), function(comparison_name) {
    constructs[[comparison_name]]$result %>%
        mutate(
            comparison = comparison_name,
            score_metric = "peak_score_max",
            score_quantile = factor(
                paste0("Q", dplyr::ntile(log10(peak_score_max), 5)),
                levels = paste0("Q", 1:5)
            ),
            significant_enriched = !is.na(padj) & padj < alpha &
                log2FoldChange > 0
        )
}))
write_table(
    score_validation,
    file.path(
        figure_dir, "construct_comparison", "omniCLIP_score_validation.tsv"
    )
)

score_correlations <- score_validation %>%
    group_by(comparison) %>%
    summarise(
        tested_peaks = n(),
        score_lfc_spearman = cor(
            peak_score_max, log2FoldChange, method = "spearman"
        ),
        .groups = "drop"
    )
write_table(
    score_correlations,
    file.path(
        figure_dir, "construct_comparison", "omniCLIP_score_correlations.tsv"
    )
)

for (comparison_name in names(constructs)) {
    score_data <- score_validation %>%
        filter(comparison == comparison_name) %>%
        mutate(
            significance = ifelse(
                significant_enriched,
                "Significantly enriched", "Other tested peak"
            )
        )
    rho <- score_correlations %>%
        filter(comparison == comparison_name) %>%
        pull(score_lfc_spearman)
    construct_label <- sub("MSI2_", "", comparison_name)

    score_scatter <- ggplot(
        score_data,
        aes(peak_score_max, log2FoldChange, color = significance)
    ) +
        geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35) +
        geom_point(size = 0.9, alpha = 0.5) +
        geom_smooth(
            method = "loess", formula = y ~ x, se = FALSE,
            linewidth = 0.8, color = "grey15"
        ) +
        scale_x_log10(labels = label_number()) +
        scale_color_manual(values = c(
            "Significantly enriched" = "#D55E00",
            "Other tested peak" = "#8A8A8A"
        )) +
        annotate(
            "label", x = Inf, y = -Inf, hjust = 1.08, vjust = -0.2,
            size = 3, label = sprintf("Spearman rho = %.3f", rho),
            linewidth = 0.2
        ) +
        labs(
            title = paste0(construct_label, ": original score versus new enrichment"),
            subtitle = "Merged-region maximum source-peak score; log10 x-axis",
            x = "Original omniCLIP peak score (peak_score_max)",
            y = "IP/input log2 fold change", color = NULL
        ) +
        theme_publication()

    score_boxplot <- ggplot(
        score_data,
        aes(score_quantile, log2FoldChange)
    ) +
        geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35) +
        geom_boxplot(
            width = 0.68, outlier.size = 0.5,
            fill = "grey85", color = "grey25"
        ) +
        labs(
            x = "Original omniCLIP score quantile (low to high)",
            y = "IP/input log2 fold change"
        ) +
        theme_publication() +
        theme(legend.position = "none")

    combined_score_plot <- score_scatter / score_boxplot +
        plot_layout(heights = c(2, 1))
    output_name <- ifelse(
        comparison_name == "MSI2_328WT",
        "omniCLIP_score_vs_328_enrichment",
        "omniCLIP_score_vs_324_enrichment"
    )
    save_plot(
        combined_score_plot,
        file.path(figure_dir, "construct_comparison", output_name),
        width = 7.3, height = 8.2
    )
}

# -----------------------------------------------------------------------------
# 9. Curated two-construct effect-size heatmap
# -----------------------------------------------------------------------------

shared_effects <- construct_comparison %>%
    filter(enriched_328, enriched_324) %>%
    mutate(selection_score = pmin(log2FoldChange_328, log2FoldChange_324)) %>%
    arrange(desc(selection_score)) %>%
    slice_head(n = 15) %>%
    mutate(selection_group = "Strongest common enrichment")

biased_328 <- construct_comparison %>%
    filter(!peak_id %in% shared_effects$peak_id, log2FoldChange_328 > 0) %>%
    arrange(desc(log2FC_difference_328_minus_324)) %>%
    slice_head(n = 15) %>%
    mutate(selection_group = "Largest 328WT-biased effect")

biased_324 <- construct_comparison %>%
    filter(
        !peak_id %in% c(shared_effects$peak_id, biased_328$peak_id),
        log2FoldChange_324 > 0
    ) %>%
    arrange(log2FC_difference_328_minus_324) %>%
    slice_head(n = 15) %>%
    mutate(selection_group = "Largest 324WT-biased effect")

effect_heatmap_peaks <- bind_rows(shared_effects, biased_328, biased_324) %>%
    mutate(selection_group = factor(
        selection_group,
        levels = c(
            "Strongest common enrichment",
            "Largest 328WT-biased effect",
            "Largest 324WT-biased effect"
        )
    ))
write_table(
    effect_heatmap_peaks,
    file.path(
        figure_dir, "construct_comparison", "combined_effect_size_heatmap_peaks.tsv"
    )
)

effect_matrix <- effect_heatmap_peaks %>%
    select(log2FoldChange_328, log2FoldChange_324) %>%
    as.matrix()
rownames(effect_matrix) <- effect_heatmap_peaks$peak_id
colnames(effect_matrix) <- c("328WT", "324WT")
effect_row_annotation <- data.frame(
    Selection = effect_heatmap_peaks$selection_group,
    row.names = effect_heatmap_peaks$peak_id
)
effect_limit <- max(abs(effect_matrix))
effect_breaks <- seq(-effect_limit, effect_limit, length.out = 102)
effect_heatmap <- pheatmap(
    effect_matrix,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
    breaks = effect_breaks,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = effect_row_annotation,
    annotation_names_row = FALSE,
    annotation_colors = list(Selection = c(
        "Strongest common enrichment" = "#7B3294",
        "Largest 328WT-biased effect" = "#E66101",
        "Largest 324WT-biased effect" = "#1B9E77"
    )),
    border_color = "grey88",
    fontsize = 8,
    fontsize_row = 6.5,
    angle_col = 0,
    main = "Curated IP/input effect sizes",
    silent = TRUE
)
save_heatmap(
    effect_heatmap,
    file.path(
        figure_dir, "construct_comparison", "combined_effect_size_heatmap"
    ),
    width = 7, height = 10.5
)

# Reproducibility record and concise terminal summary.
capture.output(sessionInfo(), file = file.path(figure_dir, "sessionInfo.txt"))

summary_table <- tibble(
    metric = c(
        "fixed_merged_peaks", "VST_QC_peaks", "328WT_tested_peaks",
        "328WT_significant", "328WT_significantly_enriched",
        "324WT_tested_peaks", "324WT_significant",
        "324WT_significantly_enriched", "common_tested_peaks",
        "significantly_enriched_in_both", "construct_LFC_Pearson",
        "construct_LFC_Spearman", "328WT_pair_ratio_Pearson",
        "328WT_pair_ratio_Spearman", "324WT_pair_ratio_Pearson",
        "324WT_pair_ratio_Spearman", "328WT_score_LFC_Spearman",
        "324WT_score_LFC_Spearman"
    ),
    value = c(
        nrow(count_matrix), sum(vst_keep), nrow(constructs$MSI2_328WT$result),
        sum(!is.na(constructs$MSI2_328WT$result$padj) &
                constructs$MSI2_328WT$result$padj < alpha),
        sum(!is.na(constructs$MSI2_328WT$result$padj) &
                constructs$MSI2_328WT$result$padj < alpha &
                constructs$MSI2_328WT$result$log2FoldChange > 0),
        nrow(constructs$MSI2_324WT$result),
        sum(!is.na(constructs$MSI2_324WT$result$padj) &
                constructs$MSI2_324WT$result$padj < alpha),
        sum(!is.na(constructs$MSI2_324WT$result$padj) &
                constructs$MSI2_324WT$result$padj < alpha &
                constructs$MSI2_324WT$result$log2FoldChange > 0),
        nrow(construct_comparison),
        sum(construct_comparison$enriched_328 & construct_comparison$enriched_324),
        construct_pearson, construct_spearman,
        constructs$MSI2_328WT$pair_pearson,
        constructs$MSI2_328WT$pair_spearman,
        constructs$MSI2_324WT$pair_pearson,
        constructs$MSI2_324WT$pair_spearman,
        score_correlations$score_lfc_spearman[
            score_correlations$comparison == "MSI2_328WT"
        ],
        score_correlations$score_lfc_spearman[
            score_correlations$comparison == "MSI2_324WT"
        ]
    )
)
write_table(summary_table, file.path(figure_dir, "visualization_summary.tsv"))

cat("Visualization workflow complete\n")
cat("Fixed merged peaks:", scales::comma(nrow(count_matrix)), "\n")
cat("328WT tested/significant/enriched:",
    nrow(constructs$MSI2_328WT$result), "/",
    summary_table$value[summary_table$metric == "328WT_significant"], "/",
    summary_table$value[summary_table$metric == "328WT_significantly_enriched"], "\n")
cat("324WT tested/significant/enriched:",
    nrow(constructs$MSI2_324WT$result), "/",
    summary_table$value[summary_table$metric == "324WT_significant"], "/",
    summary_table$value[summary_table$metric == "324WT_significantly_enriched"], "\n")
cat("Common tested / enriched in both:",
    nrow(construct_comparison), "/",
    sum(construct_comparison$enriched_328 & construct_comparison$enriched_324), "\n")
cat(sprintf("328WT vs 324WT LFC: Pearson %.3f; Spearman %.3f\n",
            construct_pearson, construct_spearman))
cat(sprintf("Pair A vs B: 328WT Pearson %.3f, Spearman %.3f; 324WT Pearson %.3f, Spearman %.3f\n",
            constructs$MSI2_328WT$pair_pearson,
            constructs$MSI2_328WT$pair_spearman,
            constructs$MSI2_324WT$pair_pearson,
            constructs$MSI2_324WT$pair_spearman))
cat(sprintf("Original score vs LFC: 328WT Spearman %.3f; 324WT Spearman %.3f\n",
            score_correlations$score_lfc_spearman[
                score_correlations$comparison == "MSI2_328WT"
            ],
            score_correlations$score_lfc_spearman[
                score_correlations$comparison == "MSI2_324WT"
            ]))
cat("Figures and plotting tables:", figure_dir, "\n")
