#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(DESeq2)
    library(dplyr)
    library(ggplot2)
    library(ggrepel)
    library(readr)
    library(ragg)
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
output_dir <- arg_value(
    "--output", file.path(results_dir, "figures", "all_nonzero_peaks")
)

pseudocount_cpm <- 0.5
kmeans_clusters <- 8L
kmeans_seed <- 328324L

required_files <- c(
    file.path(results_dir, "peak_count_matrix.tsv"),
    file.path(results_dir, "peaks.filtered_merged.tsv"),
    file.path(results_dir, "bam_library_sizes.tsv"),
    file.path(results_dir, "deseq2", "MSI2_328WT", "IP_vs_input.significant.tsv"),
    file.path(results_dir, "deseq2", "MSI2_324WT", "IP_vs_input.significant.tsv"),
    metadata_path
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
    stop("Missing required files: ", paste(missing_files, collapse = ", "))
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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
            panel.grid.major = element_line(color = "grey91", linewidth = 0.25),
            panel.grid.minor = element_blank()
        )
}

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

write_table <- function(x, path) {
    readr::write_tsv(as.data.frame(x), path, na = "NA")
}

metadata <- read_tsv(metadata_path) %>%
    mutate(
        comparison = factor(
            comparison, levels = c("MSI2_328WT", "MSI2_324WT")
        ),
        assay = factor(assay, levels = c("input", "IP")),
        pair = factor(pair, levels = c("A", "B")),
        library_label = factor(library_label, levels = c("SMI", "CLIP")),
        short_sample = sub("NM2026_", "", sample),
        plot_label = paste(short_sample, library_label, pair)
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

nonzero_keep <- rowSums(count_matrix) > 0L
if (sum(nonzero_keep) != 16046L) {
    stop(
        "Expected 16,046 peaks with any count, but found ",
        format(sum(nonzero_keep), big.mark = ",")
    )
}
nonzero_counts <- count_matrix[nonzero_keep, , drop = FALSE]
nonzero_peaks <- peaks[nonzero_keep, , drop = FALSE]

library_sizes <- library_sizes %>%
    slice(match(metadata$sample, sample))
if (!identical(library_sizes$sample, metadata$sample)) {
    stop("Library-size table could not be aligned to the sample metadata")
}
full_library_factors <- library_sizes$primary_mapped_alignments /
    exp(mean(log(library_sizes$primary_mapped_alignments)))
names(full_library_factors) <- library_sizes$sample

# The same full mapped-library totals are used for both CPM and VST size factors.
vst_metadata <- as.data.frame(metadata)
rownames(vst_metadata) <- vst_metadata$sample
vst_metadata$construct <- factor(vst_metadata$comparison)
dds <- DESeqDataSetFromMatrix(
    countData = nonzero_counts,
    colData = vst_metadata,
    design = ~ construct + assay
)
sizeFactors(dds) <- full_library_factors[colnames(dds)]
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
vst_matrix <- assay(vsd)

# -----------------------------------------------------------------------------
# PCA: every peak with at least one count in the eight-library matrix
# -----------------------------------------------------------------------------

pca <- prcomp(t(vst_matrix))
pca_variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_data <- as.data.frame(pca$x) %>%
    rownames_to_column("sample") %>%
    left_join(metadata, by = "sample")
write_table(
    pca_data,
    file.path(output_dir, "PCA_16046_nonzero_peaks.tsv")
)

assay_colors <- c(SMI = "#0072B2", CLIP = "#D55E00")
construct_shapes <- c(MSI2_328WT = 16, MSI2_324WT = 17)
pca_plot <- ggplot(
    pca_data,
    aes(
        PC1, PC2, color = library_label, shape = comparison,
        label = plot_label
    )
) +
    geom_hline(yintercept = 0, color = "grey88", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey88", linewidth = 0.3) +
    geom_point(size = 3.6) +
    ggrepel::geom_text_repel(
        size = 3, color = "grey15", seed = kmeans_seed,
        box.padding = 0.5, point.padding = 0.3, max.overlaps = Inf
    ) +
    scale_color_manual(values = assay_colors, drop = FALSE) +
    scale_shape_manual(
        values = construct_shapes,
        labels = c(MSI2_328WT = "328WT", MSI2_324WT = "324WT"),
        drop = FALSE
    ) +
    labs(
        title = "PCA of all peaks with any count",
        subtitle = paste0(
            "16,046 peaks; VST with full mapped-library size factors"
        ),
        x = sprintf("PC1 (%.1f%%)", pca_variance[1]),
        y = sprintf("PC2 (%.1f%%)", pca_variance[2]),
        color = "Library", shape = "Construct"
    ) +
    theme_publication()
save_plot(
    pca_plot,
    file.path(output_dir, "PCA_16046_nonzero_peaks"),
    width = 7.5, height = 5.6
)

# -----------------------------------------------------------------------------
# Heatmap: row-standardized VST counts for all 16,046 peaks
# -----------------------------------------------------------------------------

row_z <- t(scale(t(vst_matrix)))
row_z[!is.finite(row_z)] <- 0

set.seed(kmeans_seed)
kmeans_fit <- kmeans(
    row_z, centers = kmeans_clusters, nstart = 50, iter.max = 100
)
center_tree <- hclust(dist(kmeans_fit$centers), method = "complete")
cluster_order <- center_tree$order
cluster_remap <- setNames(seq_along(cluster_order), cluster_order)
ordered_cluster <- unname(cluster_remap[as.character(kmeans_fit$cluster)])

sample_for <- function(comparison_value, assay_value, pair_value) {
    matched <- metadata %>%
        filter(
            comparison == comparison_value,
            assay == assay_value,
            pair == pair_value
        ) %>%
        pull(sample)
    if (length(matched) != 1L) {
        stop(
            "Expected one sample for ", comparison_value, "/", assay_value,
            "/", pair_value, "; found ", length(matched)
        )
    }
    matched
}

s328_clip <- c(
    sample_for("MSI2_328WT", "IP", "A"),
    sample_for("MSI2_328WT", "IP", "B")
)
s328_smi <- c(
    sample_for("MSI2_328WT", "input", "A"),
    sample_for("MSI2_328WT", "input", "B")
)
s324_clip <- c(
    sample_for("MSI2_324WT", "IP", "A"),
    sample_for("MSI2_324WT", "IP", "B")
)
s324_smi <- c(
    sample_for("MSI2_324WT", "input", "A"),
    sample_for("MSI2_324WT", "input", "B")
)

heatmap_sort_score <-
    rowMeans(row_z[, s328_clip, drop = FALSE]) -
    rowMeans(row_z[, s328_smi, drop = FALSE]) +
    rowMeans(row_z[, s324_clip, drop = FALSE]) -
    rowMeans(row_z[, s324_smi, drop = FALSE])

heatmap_assignments <- tibble(
    peak_id = rownames(row_z),
    cluster = ordered_cluster,
    sort_score = heatmap_sort_score
) %>%
    arrange(cluster, desc(sort_score), peak_id) %>%
    group_by(cluster) %>%
    mutate(within_cluster_order = row_number()) %>%
    ungroup() %>%
    mutate(display_order = row_number())

cluster_sizes <- heatmap_assignments %>%
    count(cluster, name = "cluster_size") %>%
    mutate(cluster_label = paste0("K", cluster, " (n=", cluster_size, ")"))
heatmap_assignments <- heatmap_assignments %>%
    left_join(cluster_sizes, by = "cluster")
write_table(
    heatmap_assignments,
    file.path(output_dir, "heatmap_cluster_assignments.tsv")
)

heatmap_value_table <- data.frame(
    peak_id = rownames(row_z), row_z, check.names = FALSE
) %>%
    left_join(
        heatmap_assignments %>%
            select(
                peak_id, cluster, cluster_label, within_cluster_order,
                display_order
            ),
        by = "peak_id"
    ) %>%
    arrange(display_order) %>%
    select(
        peak_id, cluster, cluster_label, within_cluster_order, display_order,
        all_of(colnames(row_z))
    )
write_table(
    heatmap_value_table,
    file.path(output_dir, "heatmap_values.row_z.tsv")
)

heatmap_column_order <- c(
    s328_smi[1], s328_clip[1], s328_smi[2], s328_clip[2],
    s324_smi[1], s324_clip[1], s324_smi[2], s324_clip[2]
)
heatmap_column_labels <- c(
    "328 SMI A", "328 CLIP A", "328 SMI B", "328 CLIP B",
    "324 SMI A", "324 CLIP A", "324 SMI B", "324 CLIP B"
)

heatmap_long <- heatmap_value_table %>%
    select(
        peak_id, cluster_label, within_cluster_order,
        all_of(heatmap_column_order)
    ) %>%
    pivot_longer(
        cols = all_of(heatmap_column_order),
        names_to = "sample", values_to = "row_z"
    ) %>%
    mutate(
        sample = factor(
            sample,
            levels = heatmap_column_order,
            labels = heatmap_column_labels
        ),
        cluster_label = factor(
            cluster_label,
            levels = cluster_sizes$cluster_label
        )
    )

heatmap_plot <- ggplot(
    heatmap_long,
    aes(x = sample, y = within_cluster_order, fill = row_z)
) +
    geom_raster() +
    facet_grid(
        cluster_label ~ ., scales = "free_y", space = "free_y", switch = "y"
    ) +
    scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B",
        midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
        title = "Expression patterns across all peaks with any count",
        subtitle = paste0(
            "16,046 peaks; row-scaled VST counts grouped into ",
            kmeans_clusters, " k-means patterns"
        ),
        x = NULL, y = NULL, fill = "Row z-score\nof VST count"
    ) +
    theme_minimal(base_size = 9) +
    theme(
        plot.title = element_text(face = "bold", size = rel(1.15)),
        plot.subtitle = element_text(color = "grey30"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.spacing.y = grid::unit(0.8, "mm"),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, face = "bold", size = 7),
        legend.position = "right",
        legend.title = element_text(face = "bold")
    )
save_plot(
    heatmap_plot,
    file.path(output_dir, "heatmap_16046_nonzero_peaks"),
    width = 8.2, height = 9.2
)

# -----------------------------------------------------------------------------
# Descriptive CLIP/SMI CPM ratios and cross-construct scatter
# -----------------------------------------------------------------------------

library_millions <- library_sizes$primary_mapped_alignments / 1e6
names(library_millions) <- library_sizes$sample
cpm_matrix <- sweep(
    nonzero_counts, 2, library_millions[colnames(nonzero_counts)], "/"
)

ratio_metrics <- function(clip_samples, smi_samples, prefix) {
    clip_mean <- rowMeans(cpm_matrix[, clip_samples, drop = FALSE])
    smi_mean <- rowMeans(cpm_matrix[, smi_samples, drop = FALSE])
    primary_lfc <- log2(
        (clip_mean + pseudocount_cpm) / (smi_mean + pseudocount_cpm)
    )
    pair_a_lfc <- log2(
        (cpm_matrix[, clip_samples[1]] + pseudocount_cpm) /
            (cpm_matrix[, smi_samples[1]] + pseudocount_cpm)
    )
    pair_b_lfc <- log2(
        (cpm_matrix[, clip_samples[2]] + pseudocount_cpm) /
            (cpm_matrix[, smi_samples[2]] + pseudocount_cpm)
    )
    output <- data.frame(
        clip_mean, smi_mean, primary_lfc, pair_a_lfc, pair_b_lfc,
        mean_pair_lfc = (pair_a_lfc + pair_b_lfc) / 2,
        check.names = FALSE
    )
    names(output) <- paste0(
        prefix,
        c(
            "_mean_CLIP_CPM", "_mean_SMI_CPM", "_replicate_mean_LFC",
            "_pair_A_LFC", "_pair_B_LFC", "_mean_pair_LFC"
        )
    )
    output
}

ratio_328 <- ratio_metrics(s328_clip, s328_smi, "328WT")
ratio_324 <- ratio_metrics(s324_clip, s324_smi, "324WT")

sig_328 <- read_tsv(
    file.path(
        results_dir, "deseq2", "MSI2_328WT", "IP_vs_input.significant.tsv"
    )
) %>%
    filter(log2FoldChange > 0) %>%
    pull(peak_id)
sig_324 <- read_tsv(
    file.path(
        results_dir, "deseq2", "MSI2_324WT", "IP_vs_input.significant.tsv"
    )
) %>%
    filter(log2FoldChange > 0) %>%
    pull(peak_id)

ratio_table <- bind_cols(
    nonzero_peaks,
    as_tibble(ratio_328),
    as_tibble(ratio_324)
) %>%
    mutate(
        enriched_328WT = peak_id %in% sig_328,
        enriched_324WT = peak_id %in% sig_324,
        enrichment_category = case_when(
            enriched_328WT & enriched_324WT ~ "Enriched in both",
            enriched_328WT ~ "Enriched only in 328WT",
            enriched_324WT ~ "Enriched only in 324WT",
            TRUE ~ "Not significantly enriched"
        ),
        enrichment_category = factor(
            enrichment_category,
            levels = c(
                "Not significantly enriched", "Enriched only in 328WT",
                "Enriched only in 324WT", "Enriched in both"
            )
        )
    )
write_table(
    ratio_table,
    file.path(output_dir, "CLIP_SMI_pseudocount_LFC.tsv")
)

x_lfc <- ratio_table$`328WT_replicate_mean_LFC`
y_lfc <- ratio_table$`324WT_replicate_mean_LFC`
pearson_r <- cor(x_lfc, y_lfc, method = "pearson")
spearman_rho <- cor(x_lfc, y_lfc, method = "spearman")
scatter_limit <- ceiling(max(abs(c(x_lfc, y_lfc))))

category_colors <- c(
    "Not significantly enriched" = "#B8B8B8",
    "Enriched only in 328WT" = "#E66101",
    "Enriched only in 324WT" = "#1B9E77",
    "Enriched in both" = "#7B3294"
)
category_counts <- ratio_table %>%
    count(enrichment_category, .drop = FALSE) %>%
    mutate(legend_label = paste0(enrichment_category, " (n=", n, ")"))
category_labels <- setNames(
    category_counts$legend_label,
    as.character(category_counts$enrichment_category)
)

scatter_data <- ratio_table %>%
    arrange(enrichment_category != "Not significantly enriched")
scatter_plot <- ggplot(
    scatter_data,
    aes(
        x = `328WT_replicate_mean_LFC`,
        y = `324WT_replicate_mean_LFC`,
        color = enrichment_category
    )
) +
    geom_hline(yintercept = 0, color = "grey65", linewidth = 0.35) +
    geom_vline(xintercept = 0, color = "grey65", linewidth = 0.35) +
    geom_abline(
        slope = 1, intercept = 0, color = "grey40", linewidth = 0.5,
        linetype = "dashed"
    ) +
    geom_point(size = 0.8, alpha = 0.42) +
    scale_color_manual(
        values = category_colors, labels = category_labels, drop = FALSE
    ) +
    coord_equal(
        xlim = c(-scatter_limit, scatter_limit),
        ylim = c(-scatter_limit, scatter_limit), expand = FALSE
    ) +
    labs(
        title = "CLIP/SMI enrichment is compared across constructs",
        subtitle = sprintf(
            "16,046 peaks; Pearson r = %.2f, Spearman rho = %.2f",
            pearson_r, spearman_rho
        ),
        x = "328WT log2(CLIP / SMI)",
        y = "324WT log2(CLIP / SMI)",
        color = "Existing DESeq2 call",
        caption = paste0(
            "Each axis = log2[(mean CLIP CPM + 0.5)/(mean SMI CPM + 0.5)]. ",
            "CPM uses full mapped-library totals.\nColors show existing within-",
            "construct DESeq2 enrichment calls; this is not an interaction test."
        )
    ) +
    guides(
        color = guide_legend(
            nrow = 2, byrow = TRUE, title.position = "top"
        )
    ) +
    theme_publication() +
    theme(
        legend.text = element_text(size = 8),
        legend.box.just = "left"
    )
save_plot(
    scatter_plot,
    file.path(output_dir, "CLIP_SMI_pseudocount_LFC_scatter"),
    width = 7.4, height = 6.6
)

parameter_table <- tibble(
    parameter = c(
        "peak_inclusion_rule",
        "included_peak_count",
        "VST_size_factor_source",
        "heatmap_values",
        "heatmap_clustering",
        "heatmap_k",
        "kmeans_seed",
        "CPM_denominator",
        "LFC_pseudocount_CPM",
        "primary_LFC_formula",
        "replicate_averaging",
        "scatter_significance_context"
    ),
    value = c(
        "rowSums(raw counts across 8 libraries) > 0",
        sum(nonzero_keep),
        "primary mapped alignments, centered to geometric mean",
        "DESeq2 VST followed by per-peak z-score across 8 libraries",
        "k-means on row-z matrix; cluster order from center hclust",
        kmeans_clusters,
        kmeans_seed,
        "primary mapped alignments / 1e6",
        pseudocount_cpm,
        "log2((mean CLIP CPM + 0.5)/(mean SMI CPM + 0.5))",
        "two CLIP replicates averaged separately from two SMI replicates",
        "existing within-construct DESeq2 enriched calls; not interaction test"
    )
)
write_table(parameter_table, file.path(output_dir, "analysis_parameters.tsv"))

summary_table <- tibble(
    metric = c(
        "merged_peaks_total",
        "peaks_with_any_count",
        "peaks_zero_in_all_libraries",
        "PCA_PC1_percent_variance",
        "PCA_PC2_percent_variance",
        "LFC_Pearson_r_328WT_vs_324WT",
        "LFC_Spearman_rho_328WT_vs_324WT",
        "DESeq2_enriched_both",
        "DESeq2_enriched_only_328WT",
        "DESeq2_enriched_only_324WT",
        "DESeq2_not_significantly_enriched"
    ),
    value = c(
        nrow(count_matrix),
        sum(nonzero_keep),
        sum(!nonzero_keep),
        pca_variance[1],
        pca_variance[2],
        pearson_r,
        spearman_rho,
        sum(ratio_table$enrichment_category == "Enriched in both"),
        sum(ratio_table$enrichment_category == "Enriched only in 328WT"),
        sum(ratio_table$enrichment_category == "Enriched only in 324WT"),
        sum(ratio_table$enrichment_category == "Not significantly enriched")
    )
)
write_table(summary_table, file.path(output_dir, "summary.tsv"))

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

message(
    "Wrote all-nonzero-peak analysis for ", sum(nonzero_keep),
    " peaks to ", output_dir
)
