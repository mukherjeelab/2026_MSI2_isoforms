#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(preprocessCore)
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
count_path <- arg_value(
    "--count-matrix", file.path(results_dir, "peak_count_matrix.tsv")
)
library_size_path <- arg_value(
    "--library-sizes", file.path(results_dir, "bam_library_sizes.tsv")
)
sample_path <- arg_value(
    "--sample-table", "code/clip_peak_analysis/deseq_samples.tsv"
)
annotation_path <- arg_value(
    "--annotation-table",
    file.path(
        results_dir, "figures", "all_nonzero_peaks",
        "transcript_region_bias", "peak_transcript_region_annotations.tsv"
    )
)
primary_delta_path <- arg_value(
    "--primary-delta-table",
    file.path(
        results_dir, "figures", "all_nonzero_peaks",
        "transcript_region_bias", "continuous_construct_preference",
        "continuous_delta_delta_peaks.tsv"
    )
)
output_dir <- arg_value(
    "--output",
    file.path(
        results_dir, "figures", "all_nonzero_peaks",
        "transcript_region_bias", "continuous_construct_preference",
        "quantile_normalized_log2rawcounts"
    )
)

log2_rawcount_pseudocount <- 0.5
region_levels <- c("5' UTR", "CDS", "3' UTR", "Intron", "Other")
region_colors <- c(
    "5' UTR" = "#4C78A8",
    "CDS" = "#F58518",
    "3' UTR" = "#54A24B",
    "Intron" = "#B279A2",
    "Other" = "#79706E"
)

required_files <- c(
    count_path, library_size_path, sample_path, annotation_path,
    primary_delta_path
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
    stop("Missing required files: ", paste(missing_files, collapse = ", "))
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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

count_matrix <- read_tsv(count_path)
library_sizes <- read_tsv(library_size_path)
samples <- read_tsv(sample_path)
annotations <- read_tsv(annotation_path)

sample_columns <- setdiff(names(count_matrix), "peak_id")
if (length(sample_columns) != 8L) {
    stop("Expected exactly eight sample columns in the count matrix")
}
if (!setequal(sample_columns, library_sizes$sample) ||
        !setequal(sample_columns, samples$sample)) {
    stop("Sample IDs differ among count, library-size, and metadata tables")
}
if (anyDuplicated(count_matrix$peak_id) || anyDuplicated(annotations$peak_id)) {
    stop("Peak IDs must be unique in count and annotation tables")
}

counted_peaks <- count_matrix %>%
    filter(rowSums(across(all_of(sample_columns))) > 0)
if (nrow(counted_peaks) != 16046L) {
    stop("Expected 16,046 peaks with at least one raw count")
}

sample_info <- samples %>%
    select(sample, assay, comparison, pair, construct, library_label) %>%
    left_join(library_sizes, by = "sample") %>%
    mutate(sample = factor(sample, levels = sample_columns)) %>%
    arrange(sample)
if (any(is.na(sample_info$primary_mapped_alignments)) ||
        nrow(sample_info) != 8L) {
    stop("All eight samples require full mapped-library sizes")
}

raw_matrix <- as.matrix(counted_peaks[, sample_columns])
storage.mode(raw_matrix) <- "double"
library_depths <- sample_info$primary_mapped_alignments
names(library_depths) <- as.character(sample_info$sample)
library_depths <- library_depths[sample_columns]

log2_rawcount_matrix <- log2(raw_matrix + log2_rawcount_pseudocount)
quantile_matrix <- preprocessCore::normalize.quantiles(log2_rawcount_matrix)
dimnames(quantile_matrix) <- dimnames(log2_rawcount_matrix)

# With sparse count data, many peaks are tied at zero within each library.
# preprocessCore maps ranks to a common reference distribution and assigns tied
# values their averaged reference rank. This preserves zero ties, so exact
# empirical marginals can differ when libraries have different zero fractions.
if (!all(dim(quantile_matrix) == dim(log2_rawcount_matrix)) ||
        any(!is.finite(quantile_matrix))) {
    stop("Quantile normalization returned invalid dimensions or values")
}

normalized_table <- as_tibble(quantile_matrix, .name_repair = "minimal") %>%
    setNames(paste0(sample_columns, "_QN_log2_rawcount")) %>%
    mutate(peak_id = counted_peaks$peak_id, .before = 1)
write_table(normalized_table, "quantile_normalized_log2_rawcount.tsv")

sample_summary <- bind_rows(lapply(sample_columns, function(sample_id) {
    metadata <- sample_info %>% filter(as.character(sample) == sample_id)
    tibble(
        sample = sample_id,
        comparison = metadata$comparison,
        assay = metadata$assay,
        pair = metadata$pair,
        primary_mapped_alignments = metadata$primary_mapped_alignments,
        pre_qn_q1 = quantile(log2_rawcount_matrix[, sample_id], 0.25),
        pre_qn_median = median(log2_rawcount_matrix[, sample_id]),
        pre_qn_mean = mean(log2_rawcount_matrix[, sample_id]),
        pre_qn_q3 = quantile(log2_rawcount_matrix[, sample_id], 0.75),
        post_qn_q1 = quantile(quantile_matrix[, sample_id], 0.25),
        post_qn_median = median(quantile_matrix[, sample_id]),
        post_qn_mean = mean(quantile_matrix[, sample_id]),
        post_qn_q3 = quantile(quantile_matrix[, sample_id], 0.75)
    )
}))
write_table(sample_summary, "normalization_sample_summary.tsv")

sample_sets <- list(
    `328WT_CLIP` = samples %>%
        filter(comparison == "MSI2_328WT", library_label == "CLIP") %>%
        pull(sample),
    `328WT_SMI` = samples %>%
        filter(comparison == "MSI2_328WT", library_label == "SMI") %>%
        pull(sample),
    `324WT_CLIP` = samples %>%
        filter(comparison == "MSI2_324WT", library_label == "CLIP") %>%
        pull(sample),
    `324WT_SMI` = samples %>%
        filter(comparison == "MSI2_324WT", library_label == "SMI") %>%
        pull(sample)
)
if (any(lengths(sample_sets) != 2L)) {
    stop("Expected two CLIP and two SMI libraries for each construct")
}

row_mean <- function(sample_ids) {
    rowMeans(quantile_matrix[, sample_ids, drop = FALSE])
}
enrichment_328 <- row_mean(sample_sets$`328WT_CLIP`) -
    row_mean(sample_sets$`328WT_SMI`)
enrichment_324 <- row_mean(sample_sets$`324WT_CLIP`) -
    row_mean(sample_sets$`324WT_SMI`)

analysis_peaks <- tibble(
    peak_id = counted_peaks$peak_id,
    `328WT_QN_log2_CLIP_minus_SMI` = enrichment_328,
    `324WT_QN_log2_CLIP_minus_SMI` = enrichment_324,
    delta_delta = enrichment_324 - enrichment_328
) %>%
    left_join(
        annotations %>%
            select(
                peak_id, representative_source_peak_id,
                representative_source_score, transcript_region,
                transcript_region_detail
            ),
        by = "peak_id"
    ) %>%
    mutate(
        transcript_region = factor(transcript_region, levels = region_levels),
        peak_universe_rule =
            "raw count > 0 in at least one of eight libraries",
        normalization = paste0(
            "raw counts; log2(count + ", log2_rawcount_pseudocount,
            "); quantile normalization across eight samples"
        )
    )
if (nrow(analysis_peaks) != 16046L ||
        any(is.na(analysis_peaks$transcript_region)) ||
        any(!is.finite(analysis_peaks$delta_delta))) {
    stop("The normalized continuous analysis must contain 16,046 annotated finite peaks")
}
write_table(
    analysis_peaks %>% arrange(transcript_region, peak_id),
    "continuous_delta_delta_peaks.tsv"
)

primary_delta <- read_tsv(primary_delta_path) %>%
    select(peak_id, primary_cpm_delta_delta = delta_delta)
normalization_comparison <- analysis_peaks %>%
    select(peak_id, quantile_normalized_delta_delta = delta_delta) %>%
    inner_join(primary_delta, by = "peak_id") %>%
    summarise(
        peak_count = n(),
        pearson_correlation = cor(
            primary_cpm_delta_delta, quantile_normalized_delta_delta
        ),
        spearman_correlation = cor(
            primary_cpm_delta_delta, quantile_normalized_delta_delta,
            method = "spearman"
        ),
        sign_agreement_fraction = mean(
            sign(primary_cpm_delta_delta) ==
                sign(quantile_normalized_delta_delta)
        ),
        primary_overall_median = median(primary_cpm_delta_delta),
        quantile_normalized_overall_median =
            median(quantile_normalized_delta_delta)
    )
if (normalization_comparison$peak_count != 16046L) {
    stop("Primary and quantile-normalized peak tables do not match")
}
write_table(
    normalization_comparison,
    "normalization_comparison_to_primary.tsv"
)

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
    comparison = "quantile-normalized delta_delta by transcript region",
    statistic = unname(kruskal_result$statistic),
    degrees_of_freedom = unname(kruskal_result$parameter),
    p_value = kruskal_result$p.value,
    epsilon_squared = epsilon_squared,
    peak_count = nrow(analysis_peaks),
    nonindependence_note =
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
            "%s  (n = %s; median = %.3f)", transcript_region,
            format(n_peaks, big.mark = ","), median_delta_delta
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
        title = "Quantile-normalized construct preference by transcript region",
        subtitle = paste0(
            "All 16,046 counted peaks; raw counts -> log2(count + ",
            log2_rawcount_pseudocount, ") -> quantile normalization across samples"
        ),
        x = NULL,
        y = "Delta-delta (positive favors 324WT; negative favors 328WT)",
        caption = paste0(
            "Delta-delta is the 324WT minus 328WT difference between mean ",
            "quantile-normalized log2 CLIP and SMI abundance. No direction-dependent filter."
        )
    ) +
    theme_publication() +
    theme(panel.grid.major.y = element_blank())
save_plot(
    distribution_plot,
    file.path(output_dir, "continuous_delta_delta_by_region"),
    width = 10.2, height = 6.2
)

focused_peaks <- analysis_peaks %>%
    filter(transcript_region %in% c("3' UTR", "Intron")) %>%
    droplevels()
focused_summary <- region_summary %>%
    filter(transcript_region %in% c("3' UTR", "Intron"))
focused_labels <- focused_summary %>%
    mutate(
        display = sprintf(
            "%s (n = %s; median = %.3f)", transcript_region,
            format(n_peaks, big.mark = ","), median_delta_delta
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
        title = "Quantile-normalized 3' UTR and intron delta-delta densities",
        subtitle = "All counted peaks; identical filtering and normalization to the five-region plot",
        x = "Delta-delta (324WT enrichment - 328WT enrichment)",
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
        title = "Quantile-normalized 3' UTR and intron delta-delta ECDFs",
        subtitle = "The ECDF compares the full normalized distributions",
        x = "Delta-delta (324WT enrichment - 328WT enrichment)",
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
    filter(region_1 == "3' UTR", region_2 == "Intron")
if (nrow(focused_test) != 1L) {
    stop("Could not isolate the 3' UTR versus intron result")
}

readme_lines <- c(
    "# Quantile-normalized log2-raw-count sensitivity analysis",
    "",
    "This directory contains a normalization sensitivity analysis of the same",
    "16,046 construct-independent counted peaks used by the primary continuous",
    "construct-preference workflow. It does not replace the original CPM-based",
    "analysis.",
    "",
    "## Normalization",
    "",
    "1. Subset the raw 43,794-peak matrix to peaks with a count above zero in at",
    "   least one of the eight libraries (16,046 peaks).",
    "2. Use raw peak counts directly, without CPM scaling.",
    paste0(
        "3. Transform each sample as `log2(count + ", log2_rawcount_pseudocount,
        ")`."
    ),
    "4. Quantile-normalize the eight transformed sample columns together using",
    "   `preprocessCore::normalize.quantiles`.",
    "5. For each construct, subtract the mean normalized SMI abundance from the",
    "   mean normalized CLIP abundance. Delta-delta is 324WT enrichment minus",
    "   328WT enrichment.",
    "",
    "This filter is independent of delta-delta sign, preferred construct, and",
    "transcript-region annotation. Quantile normalization maps sample ranks to a",
    "shared reference distribution. Tied zero-count peaks remain tied and receive",
    "their averaged reference rank, so sparse libraries do not have exactly",
    "identical empirical marginals after normalization. The procedure is still",
    "deliberately aggressive and should be interpreted as a sensitivity analysis",
    "rather than automatically as a superior normalization.",
    "",
    "## Results",
    "",
    paste0(
        "The overall Kruskal-Wallis test has H = ",
        format(unname(kruskal_result$statistic), digits = 5), ", df = ",
        unname(kruskal_result$parameter), ", p = ",
        format(kruskal_result$p.value, scientific = TRUE, digits = 4),
        ", and epsilon-squared = ", format(epsilon_squared, digits = 4), "."
    ),
    paste0(
        "For 3' UTR versus intron, the median difference is ",
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
    paste0(
        "Compared with the primary CPM-pseudocount delta-delta values, the",
        " quantile-normalized values have Pearson r = ",
        format(normalization_comparison$pearson_correlation, digits = 4),
        ", Spearman rho = ",
        format(normalization_comparison$spearman_correlation, digits = 4),
        ", and ",
        format(100 * normalization_comparison$sign_agreement_fraction, digits = 4),
        "% sign agreement."
    ),
    "",
    "Because sample sizes are large, p-values should be interpreted alongside",
    "median differences and Cliff's delta. These are peak-level exploratory tests",
    "and do not model dependence among peaks from the same gene or transcript.",
    "",
    "## Files",
    "",
    "- `quantile_normalized_log2_rawcount.tsv`: the normalized 16,046 by 8 matrix.",
    "- `normalization_sample_summary.tsv`: sample distributions before and after",
    "  quantile normalization.",
    "- `normalization_comparison_to_primary.tsv`: correlation, sign agreement,",
    "  and overall medians versus the primary CPM-pseudocount delta-delta values.",
    "- `continuous_delta_delta_peaks.tsv`: normalized construct enrichments,",
    "  delta-delta values, and transcript annotations for every eligible peak.",
    "- `continuous_delta_delta_by_region.{pdf,png}`: five-region continuous plot.",
    "- `delta_delta_3UTR_vs_intron_density.{pdf,png}` and",
    "  `delta_delta_3UTR_vs_intron_ecdf.{pdf,png}`: focused comparisons.",
    "- `delta_delta_region_summary.tsv`, `delta_delta_region_overall_tests.tsv`,",
    "  and `delta_delta_region_pairwise_tests.tsv`: summaries and exploratory tests.",
    "- `sessionInfo.txt`: R and package versions."
)
writeLines(readme_lines, file.path(output_dir, "README.md"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

message(
    "Quantile-normalized continuous construct-preference analysis complete for ",
    format(nrow(analysis_peaks), big.mark = ","), " peaks; outputs written to ",
    output_dir
)
