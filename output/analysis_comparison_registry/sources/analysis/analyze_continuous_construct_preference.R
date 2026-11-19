#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(readr)
    library(ragg)
    library(tibble)
    library(tidyr)
})

options(readr.show_col_types = FALSE)

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default) {
    index <- match(flag, args)
    if (is.na(index)) return(default)
    if (index == length(args)) stop("Missing value after ", flag)
    args[[index + 1]]
}

results_dir <- arg_value("--results", "data/results/peak_analysis")
lfc_path <- arg_value(
    "--lfc-table",
    file.path(
        results_dir, "figures", "all_nonzero_peaks",
        "CLIP_SMI_pseudocount_LFC.tsv"
    )
)
annotation_path <- arg_value(
    "--annotation-table",
    file.path(
        results_dir, "figures", "all_nonzero_peaks",
        "transcript_region_bias", "peak_transcript_region_annotations.tsv"
    )
)
count_path <- arg_value(
    "--count-matrix",
    file.path(results_dir, "peak_count_matrix.tsv")
)
output_dir <- arg_value(
    "--output",
    file.path(
        results_dir, "figures", "all_nonzero_peaks",
        "transcript_region_bias", "continuous_construct_preference"
    )
)

required_files <- c(lfc_path, annotation_path, count_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
    stop("Missing required files: ", paste(missing_files, collapse = ", "))
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

region_levels <- c("5' UTR", "CDS", "3' UTR", "Intron", "Other")
region_colors <- c(
    "5' UTR" = "#4C78A8",
    "CDS" = "#F58518",
    "3' UTR" = "#54A24B",
    "Intron" = "#B279A2",
    "Other" = "#79706E"
)

theme_publication <- function(base_size = 11) {
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
            panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
            panel.grid.minor = element_blank()
        )
}

save_plot <- function(plot, stem, width, height) {
    ggsave(
        paste0(stem, ".pdf"), plot = plot, device = grDevices::pdf,
        width = width, height = height, units = "in", bg = "white"
    )
    ggsave(
        paste0(stem, ".png"), plot = plot, device = ragg::agg_png,
        width = width, height = height, units = "in", dpi = 300,
        bg = "white"
    )
}

write_table <- function(x, filename) {
    write_tsv(as.data.frame(x), file.path(output_dir, filename), na = "NA")
}

lfc_peaks <- read_tsv(lfc_path)
annotations <- read_tsv(annotation_path)
count_matrix <- read_tsv(count_path)

required_lfc_columns <- c(
    "peak_id", "chromosome", "start_0based", "end_0based_exclusive",
    "strand", "gene_ids", "328WT_replicate_mean_LFC",
    "324WT_replicate_mean_LFC"
)
required_annotation_columns <- c(
    "peak_id", "transcript_region", "transcript_region_detail",
    "representative_source_peak_id", "representative_source_score"
)
if (!all(required_lfc_columns %in% names(lfc_peaks))) {
    stop(
        "LFC table lacks required columns: ",
        paste(setdiff(required_lfc_columns, names(lfc_peaks)), collapse = ", ")
    )
}
if (!all(required_annotation_columns %in% names(annotations))) {
    stop(
        "Annotation table lacks required columns: ",
        paste(
            setdiff(required_annotation_columns, names(annotations)),
            collapse = ", "
        )
    )
}

count_columns <- setdiff(names(count_matrix), "peak_id")
if (length(count_columns) != 8L) {
    stop("Expected eight sample columns in the raw peak-count matrix")
}
counted_peak_ids <- count_matrix %>%
    filter(rowSums(across(all_of(count_columns))) > 0) %>%
    pull(peak_id)
if (length(counted_peak_ids) != 16046L || anyDuplicated(counted_peak_ids)) {
    stop("Expected 16,046 unique peaks with at least one raw count")
}
if (nrow(lfc_peaks) != 16046L || anyDuplicated(lfc_peaks$peak_id)) {
    stop("Expected 16,046 unique rows in the existing LFC table")
}
if (!setequal(counted_peak_ids, lfc_peaks$peak_id)) {
    stop("The LFC table does not exactly match the construct-independent counted universe")
}
if (anyDuplicated(annotations$peak_id)) {
    stop("Transcript-region annotation contains duplicated peak IDs")
}

analysis_peaks <- lfc_peaks %>%
    inner_join(
        annotations %>% select(all_of(required_annotation_columns)),
        by = "peak_id"
    ) %>%
    mutate(
        transcript_region = factor(transcript_region, levels = region_levels),
        delta_delta =
            `324WT_replicate_mean_LFC` - `328WT_replicate_mean_LFC`,
        counted_universe_rule =
            "raw count > 0 in at least one of eight libraries"
    )

if (nrow(analysis_peaks) != 16046L || any(is.na(analysis_peaks$transcript_region))) {
    stop("All 16,046 counted peaks must have one simplified transcript-region annotation")
}
if (any(!is.finite(analysis_peaks$delta_delta))) {
    stop("All delta-delta values must be finite")
}

peak_output <- analysis_peaks %>%
    arrange(transcript_region, peak_id) %>%
    select(
        peak_id, chromosome, start_0based, end_0based_exclusive, strand,
        gene_ids, representative_source_peak_id, representative_source_score,
        transcript_region, transcript_region_detail,
        `324WT_replicate_mean_LFC`, `328WT_replicate_mean_LFC`, delta_delta,
        counted_universe_rule
    )
write_table(peak_output, "continuous_delta_delta_peaks.tsv")

region_summary <- analysis_peaks %>%
    group_by(transcript_region) %>%
    summarise(
        n_peaks = n(),
        median_delta_delta = median(delta_delta),
        q1_delta_delta = quantile(delta_delta, 0.25),
        q3_delta_delta = quantile(delta_delta, 0.75),
        mean_delta_delta = mean(delta_delta),
        sd_delta_delta = sd(delta_delta),
        percent_delta_positive = 100 * mean(delta_delta > 0),
        percent_delta_negative = 100 * mean(delta_delta < 0),
        .groups = "drop"
    ) %>%
    arrange(match(as.character(transcript_region), region_levels))
write_table(region_summary, "delta_delta_region_summary.tsv")

kruskal_result <- kruskal.test(delta_delta ~ transcript_region, data = analysis_peaks)
group_count <- nrow(region_summary)
epsilon_squared <- max(
    0,
    (unname(kruskal_result$statistic) - group_count + 1) /
        (nrow(analysis_peaks) - group_count)
)
overall_tests <- tibble(
    test = "Kruskal-Wallis rank-sum test",
    comparison = "delta_delta by five-category transcript-region annotation",
    peak_universe =
        "16,046 merged peaks with raw count > 0 in at least one of eight libraries",
    statistic = unname(kruskal_result$statistic),
    degrees_of_freedom = unname(kruskal_result$parameter),
    p_value = kruskal_result$p.value,
    epsilon_squared = epsilon_squared,
    interpretation_note =
        "Peak-level exploratory test; peaks from the same gene/transcript may not be independent"
)
write_table(overall_tests, "delta_delta_region_overall_tests.tsv")

cliffs_magnitude <- function(value) {
    absolute <- abs(value)
    case_when(
        absolute < 0.147 ~ "negligible",
        absolute < 0.330 ~ "small",
        absolute < 0.474 ~ "medium",
        TRUE ~ "large"
    )
}

region_pairs <- combn(region_levels, 2, simplify = FALSE)
pairwise_tests <- lapply(region_pairs, function(pair) {
    first <- analysis_peaks %>%
        filter(transcript_region == pair[[1]]) %>%
        pull(delta_delta)
    second <- analysis_peaks %>%
        filter(transcript_region == pair[[2]]) %>%
        pull(delta_delta)
    wilcoxon_result <- suppressWarnings(wilcox.test(
        first, second, alternative = "two.sided", exact = FALSE
    ))
    ks_result <- suppressWarnings(ks.test(
        first, second, alternative = "two.sided", exact = FALSE
    ))
    wilcoxon_u <- unname(wilcoxon_result$statistic)
    cliffs_delta <- 2 * wilcoxon_u / (length(first) * length(second)) - 1
    tibble(
        region_1 = pair[[1]],
        region_2 = pair[[2]],
        comparison = paste(pair[[1]], "vs", pair[[2]]),
        n_region_1 = length(first),
        n_region_2 = length(second),
        median_region_1 = median(first),
        median_region_2 = median(second),
        median_difference_region_1_minus_region_2 =
            median(first) - median(second),
        wilcoxon_rank_sum_W = wilcoxon_u,
        wilcoxon_p_value = wilcoxon_result$p.value,
        cliffs_delta_region_1_vs_region_2 = cliffs_delta,
        cliffs_delta_magnitude = cliffs_magnitude(cliffs_delta),
        ks_D = unname(ks_result$statistic),
        ks_p_value = ks_result$p.value
    )
}) %>%
    bind_rows() %>%
    mutate(
        wilcoxon_padj_BH = p.adjust(wilcoxon_p_value, method = "BH"),
        ks_padj_BH = p.adjust(ks_p_value, method = "BH"),
        nonindependence_note =
            "Peak-level exploratory comparison; within-gene/transcript dependence is not modeled"
    )
write_table(pairwise_tests, "delta_delta_region_pairwise_tests.tsv")

region_labels <- region_summary %>%
    mutate(
        display = sprintf(
            "%s  (n = %s; median = %.3f)",
            transcript_region, format(n_peaks, big.mark = ","),
            median_delta_delta
        )
    ) %>%
    select(transcript_region, display) %>%
    deframe()

distribution_plot <- ggplot(
    analysis_peaks,
    aes(x = transcript_region, y = delta_delta, fill = transcript_region)
) +
    geom_hline(yintercept = 0, color = "#9E2A2B", linewidth = 0.55) +
    geom_violin(
        scale = "width", trim = FALSE, color = "grey45", linewidth = 0.35,
        alpha = 0.75
    ) +
    geom_boxplot(
        width = 0.13, outlier.shape = NA, fill = "white", color = "grey20",
        linewidth = 0.4
    ) +
    stat_summary(
        fun = median, geom = "point", shape = 21, size = 2.2,
        fill = "black", color = "white", stroke = 0.4
    ) +
    scale_fill_manual(values = region_colors, guide = "none") +
    scale_x_discrete(labels = region_labels, limits = region_levels) +
    coord_flip() +
    labs(
        title = "Continuous 324WT-versus-328WT preference by transcript region",
        subtitle = paste0(
            "All 16,046 construct-independent counted peaks; delta-delta = ",
            "324WT log2(CLIP/SMI) - 328WT log2(CLIP/SMI)"
        ),
        x = NULL,
        y = "Delta-delta (positive favors 324WT; negative favors 328WT)",
        caption = paste0(
            "Red line marks delta-delta = 0. Boxplots show median and IQR; ",
            "labels report peak count and median. No delta-delta cutoff was applied."
        )
    ) +
    theme_publication() +
    theme(panel.grid.major.y = element_blank())
save_plot(
    distribution_plot,
    file.path(output_dir, "continuous_delta_delta_by_region"),
    width = 10.2, height = 6.2
)

joint_limits <- range(c(
    analysis_peaks$`324WT_replicate_mean_LFC`,
    analysis_peaks$`328WT_replicate_mean_LFC`
))
joint_padding <- diff(joint_limits) * 0.025
joint_limits <- joint_limits + c(-joint_padding, joint_padding)
facet_labels <- region_summary %>%
    mutate(
        display = paste0(
            transcript_region, "\n(n = ", format(n_peaks, big.mark = ","), ")"
        )
    ) %>%
    select(transcript_region, display) %>%
    deframe()

scatter_plot <- ggplot(
    analysis_peaks,
    aes(
        x = `324WT_replicate_mean_LFC`,
        y = `328WT_replicate_mean_LFC`
    )
) +
    geom_bin_2d(aes(fill = after_stat(ndensity)), bins = 38) +
    geom_abline(
        slope = 1, intercept = 0, linetype = "dashed", color = "grey15",
        linewidth = 0.55
    ) +
    facet_wrap(
        ~transcript_region, ncol = 3,
        labeller = as_labeller(facet_labels)
    ) +
    scale_fill_viridis_c(
        option = "C", limits = c(0, 1),
        name = "Relative bin density"
    ) +
    coord_equal(xlim = joint_limits, ylim = joint_limits, expand = FALSE) +
    labs(
        title = "Joint construct enrichment by transcript region",
        subtitle = paste0(
            "Common axes across facets; below y = x favors 324WT, above y = x ",
            "favors 328WT"
        ),
        x = "324WT log2(CLIP/SMI)",
        y = "328WT log2(CLIP/SMI)",
        caption = paste0(
            "Two-dimensional bins are colored by density normalized within each ",
            "facet. All 16,046 counted peaks are retained."
        )
    ) +
    theme_publication() +
    theme(
        legend.position = "right",
        panel.grid = element_blank()
    )
save_plot(
    scatter_plot,
    file.path(output_dir, "construct_enrichment_scatter_by_region"),
    width = 10.5, height = 8.2
)

focused_peaks <- analysis_peaks %>%
    filter(transcript_region %in% c("3' UTR", "Intron")) %>%
    droplevels()
focused_summary <- region_summary %>%
    filter(transcript_region %in% c("3' UTR", "Intron"))
focused_labels <- focused_summary %>%
    mutate(
        display = sprintf(
            "%s (n = %s; median = %.3f)",
            transcript_region, format(n_peaks, big.mark = ","),
            median_delta_delta
        )
    ) %>%
    select(transcript_region, display) %>%
    deframe()

density_plot <- ggplot(
    focused_peaks,
    aes(x = delta_delta, color = transcript_region, fill = transcript_region)
) +
    geom_vline(xintercept = 0, color = "#9E2A2B", linewidth = 0.55) +
    geom_density(alpha = 0.28, linewidth = 0.8, adjust = 1) +
    geom_vline(
        data = focused_summary,
        aes(xintercept = median_delta_delta, color = transcript_region),
        linetype = "dashed", linewidth = 0.75, show.legend = FALSE
    ) +
    scale_color_manual(values = region_colors, labels = focused_labels) +
    scale_fill_manual(values = region_colors, labels = focused_labels) +
    labs(
        title = "3' UTR and intron delta-delta densities",
        subtitle = "Identical 16,046-peak universe and enrichment calculation used in the five-region analysis",
        x = "Delta-delta (324WT LFC - 328WT LFC)",
        y = "Density",
        color = "Transcript region",
        fill = "Transcript region",
        caption = "Solid red line marks zero; colored dashed lines mark group medians."
    ) +
    theme_publication()
save_plot(
    density_plot,
    file.path(output_dir, "delta_delta_3UTR_vs_intron_density"),
    width = 9.2, height = 5.8
)

ecdf_plot <- ggplot(
    focused_peaks,
    aes(x = delta_delta, color = transcript_region)
) +
    geom_vline(xintercept = 0, color = "#9E2A2B", linewidth = 0.55) +
    stat_ecdf(linewidth = 0.9, geom = "step") +
    scale_color_manual(values = region_colors, labels = focused_labels) +
    labs(
        title = "3' UTR and intron delta-delta ECDFs",
        subtitle = "The ECDF compares the complete empirical distributions, not only their medians",
        x = "Delta-delta (324WT LFC - 328WT LFC)",
        y = "Cumulative fraction of peaks",
        color = "Transcript region",
        caption = "Solid red line marks delta-delta = 0."
    ) +
    theme_publication()
save_plot(
    ecdf_plot,
    file.path(output_dir, "delta_delta_3UTR_vs_intron_ecdf"),
    width = 9.2, height = 5.8
)

focused_test <- pairwise_tests %>%
    filter(
        (region_1 == "3' UTR" & region_2 == "Intron") |
            (region_1 == "Intron" & region_2 == "3' UTR")
    )
if (nrow(focused_test) != 1L) stop("Could not isolate the 3' UTR versus intron comparison")

readme_lines <- c(
    "# Continuous construct-preference analysis",
    "",
    "This directory contains the primary continuous analysis of whether relative",
    "324WT-versus-328WT CLIP enrichment differs by transcript-region annotation.",
    "It is separate from the older threshold-based delta-group analysis in the",
    "parent directory.",
    "",
    "## Peak universe and score",
    "",
    "The analysis includes all 16,046 merged peaks with a raw count greater than",
    "zero in at least one of the eight libraries. This filter is independent of",
    "delta-delta sign, preferred construct, and transcript-region annotation.",
    "No additional abundance or direction-dependent filter is applied.",
    "",
    "For each peak:",
    "",
    "`delta_delta = 324WT replicate-mean log2(CLIP/SMI) - 328WT replicate-mean log2(CLIP/SMI)`",
    "",
    "Each construct-specific enrichment is the existing CPM-normalized ratio",
    "`log2((mean CLIP CPM + 0.5)/(mean SMI CPM + 0.5))`. Positive delta-delta",
    "favors 324WT and negative delta-delta favors 328WT.",
    "",
    "## Statistical tests",
    "",
    paste0(
        "The overall five-region Kruskal-Wallis test has H = ",
        format(unname(kruskal_result$statistic), digits = 5), ", df = ",
        unname(kruskal_result$parameter), ", p = ",
        format(kruskal_result$p.value, scientific = TRUE, digits = 4),
        ", and epsilon-squared = ", format(epsilon_squared, digits = 4), "."
    ),
    "Pairwise Wilcoxon rank-sum tests compare locations and pairwise two-sample",
    "Kolmogorov-Smirnov tests compare complete distributions. Each family of ten",
    "pairwise p-values is adjusted separately using Benjamini-Hochberg.",
    "",
    paste0(
        "For 3' UTR versus intron, the median difference (3' UTR minus intron) is ",
        format(
            focused_test$median_difference_region_1_minus_region_2,
            digits = 4
        ), ", Cliff's delta is ",
        format(focused_test$cliffs_delta_region_1_vs_region_2, digits = 4),
        " (", focused_test$cliffs_delta_magnitude, "), Wilcoxon BH q = ",
        format(focused_test$wilcoxon_padj_BH, scientific = TRUE, digits = 4),
        ", and KS BH q = ",
        format(focused_test$ks_padj_BH, scientific = TRUE, digits = 4), "."
    ),
    "",
    "Small p-values should be interpreted together with median differences and",
    "rank-based effect sizes because the peak counts are large. These are",
    "peak-level exploratory tests; peaks from the same gene or transcript may not",
    "be statistically independent.",
    "",
    "## Files",
    "",
    "- `continuous_delta_delta_by_region.{pdf,png}`: violin/boxplot distributions",
    "  for all five transcript-region categories, with n and median in labels.",
    "- `construct_enrichment_scatter_by_region.{pdf,png}`: density-binned 324WT",
    "  versus 328WT enrichment facets with common axes and the y = x diagonal.",
    "- `delta_delta_3UTR_vs_intron_density.{pdf,png}`: focused density comparison.",
    "- `delta_delta_3UTR_vs_intron_ecdf.{pdf,png}`: focused ECDF comparison.",
    "- `continuous_delta_delta_peaks.tsv`: all peak-level values used in figures",
    "  and tests.",
    "- `delta_delta_region_summary.tsv`: n, median, quartiles, mean, SD, and sign",
    "  percentages by transcript region.",
    "- `delta_delta_region_overall_tests.tsv`: overall Kruskal-Wallis result and",
    "  epsilon-squared effect size.",
    "- `delta_delta_region_pairwise_tests.tsv`: pairwise Wilcoxon and KS tests, BH",
    "  q-values, median differences, and Cliff's delta.",
    "- `sessionInfo.txt`: R and package versions used for this run."
)
writeLines(readme_lines, file.path(output_dir, "README.md"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

message(
    "Continuous construct-preference analysis complete for ",
    format(nrow(analysis_peaks), big.mark = ","), " peaks; outputs written to ",
    output_dir
)
