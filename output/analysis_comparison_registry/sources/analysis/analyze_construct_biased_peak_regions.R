#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(GenomicFeatures)
    library(GenomicRanges)
    library(ggplot2)
    library(IRanges)
    library(readr)
    library(ragg)
    library(tibble)
    library(tidyr)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
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
lfc_path <- arg_value(
    "--lfc-table",
    file.path(
        results_dir, "figures", "all_nonzero_peaks",
        "CLIP_SMI_pseudocount_LFC.tsv"
    )
)
output_dir <- arg_value(
    "--output",
    file.path(
        results_dir, "figures", "all_nonzero_peaks",
        "transcript_region_bias"
    )
)

strict_lfc_threshold <- 1
delta_lfc_threshold <- 1
lfc_pseudocount_cpm <- 0.5
chisq_monte_carlo_B <- 100000L
chisq_seed <- 324328L

required_files <- c(
    lfc_path,
    file.path(results_dir, "peak_members.tsv"),
    file.path(results_dir, "peaks.filtered_merged.tsv")
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
            panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
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

lfc_peaks <- read_tsv(lfc_path)
peak_members <- read_tsv(file.path(results_dir, "peak_members.tsv"))
merged_peaks <- read_tsv(file.path(results_dir, "peaks.filtered_merged.tsv"))

required_lfc_columns <- c(
    "peak_id", "328WT_replicate_mean_LFC", "324WT_replicate_mean_LFC",
    "gene_ids", "chromosome", "strand"
)
if (!all(required_lfc_columns %in% names(lfc_peaks))) {
    stop(
        "LFC table lacks required columns: ",
        paste(setdiff(required_lfc_columns, names(lfc_peaks)), collapse = ", ")
    )
}
if (nrow(lfc_peaks) != 16046L || anyDuplicated(lfc_peaks$peak_id)) {
    stop("Expected 16,046 unique nonzero peak rows in the LFC table")
}
if (nrow(peak_members) != 46108L || nrow(merged_peaks) != 43794L) {
    stop("Expected 46,108 source peaks and 43,794 merged peaks")
}

# For a merged region, use the summit from the highest-scoring contributing
# omniCLIP peak. This gives one biologically motivated annotation position per
# merged peak and avoids assigning a long peak by whichever feature covers the
# greatest number of bases.
representative_summits <- peak_members %>%
    arrange(peak_id, desc(source_score), source_line) %>%
    distinct(peak_id, .keep_all = TRUE) %>%
    transmute(
        peak_id,
        representative_source_peak_id = source_peak_id,
        representative_source_score = source_score,
        representative_gene_id = gene_id,
        summit_chromosome = chromosome,
        summit_start_0based = as.integer(summit_start),
        summit_end_0based_exclusive = as.integer(summit_end),
        summit_strand = strand
    )
if (nrow(representative_summits) != nrow(merged_peaks)) {
    stop("Could not identify one representative summit for every merged peak")
}

detail_region_levels <- c(
    "5' UTR", "CDS", "3' UTR", "Other exon/noncoding",
    "Ambiguous UTR/CDS", "Intron", "Intergenic"
)
region_levels <- c("5' UTR", "CDS", "3' UTR", "Intron", "Other")

annotation_cache_path <- file.path(
    output_dir, "peak_transcript_region_annotations.tsv"
)
source_annotation_cache_path <- file.path(
    output_dir, "source_peak_transcript_region_annotations.tsv"
)
required_annotation_columns <- c(
    "transcript_region_detail", "transcript_region", "overlaps_5UTR",
    "overlaps_CDS", "overlaps_3UTR", "overlaps_any_exon",
    "overlaps_transcript_span", "annotation_source"
)

cache_valid <- FALSE
if (file.exists(annotation_cache_path) &&
        file.exists(source_annotation_cache_path)) {
    cached_tables <- tryCatch(
        list(
            merged = read_tsv(annotation_cache_path),
            source = read_tsv(source_annotation_cache_path)
        ),
        error = function(e) NULL
    )
    if (!is.null(cached_tables)) {
        annotation_table <- cached_tables$merged
        source_annotation_table <- cached_tables$source
        cache_valid <-
            all(c("peak_id", required_annotation_columns) %in%
                names(annotation_table)) &&
            all(c("source_line", required_annotation_columns) %in%
                names(source_annotation_table)) &&
            nrow(annotation_table) == nrow(merged_peaks) &&
            nrow(source_annotation_table) == nrow(peak_members) &&
            !anyDuplicated(annotation_table$peak_id) &&
            !anyDuplicated(source_annotation_table$source_line) &&
            setequal(annotation_table$peak_id, merged_peaks$peak_id) &&
            setequal(source_annotation_table$source_line, peak_members$source_line)
    }
}

if (cache_valid) {
    annotation_table <- annotation_table %>%
        mutate(
            transcript_region_detail = factor(
                transcript_region_detail, levels = detail_region_levels
            ),
            transcript_region = factor(transcript_region, levels = region_levels)
        )
    source_annotation_table <- source_annotation_table %>%
        mutate(
            transcript_region_detail = factor(
                transcript_region_detail, levels = detail_region_levels
            ),
            transcript_region = factor(transcript_region, levels = region_levels)
        )
    message("Loaded validated transcript-region annotation caches")
} else {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    reduce_feature <- function(x) {
        reduce(unlist(x, use.names = FALSE), ignore.strand = FALSE)
    }
    five_utr_ranges <- reduce_feature(fiveUTRsByTranscript(txdb))
    three_utr_ranges <- reduce_feature(threeUTRsByTranscript(txdb))
    cds_ranges <- reduce(cds(txdb), ignore.strand = FALSE)
    exon_ranges <- reduce(exons(txdb), ignore.strand = FALSE)
    transcript_ranges <- reduce(transcripts(txdb), ignore.strand = FALSE)

    annotate_summits <- function(summit_table, id_column) {
        summit_ranges <- GRanges(
            seqnames = summit_table$summit_chromosome,
            ranges = IRanges(
                start = summit_table$summit_start_0based + 1L,
                end = summit_table$summit_end_0based_exclusive
            ),
            strand = summit_table$summit_strand
        )
        overlap_flag <- function(subject) {
            countOverlaps(summit_ranges, subject, ignore.strand = FALSE) > 0L
        }
        in_five_utr <- overlap_flag(five_utr_ranges)
        in_cds <- overlap_flag(cds_ranges)
        in_three_utr <- overlap_flag(three_utr_ranges)
        in_exon <- overlap_flag(exon_ranges)
        in_transcript <- overlap_flag(transcript_ranges)
        primary_exonic_hits <- rowSums(cbind(
            in_five_utr, in_cds, in_three_utr
        ))
        transcript_region_detail <- case_when(
            primary_exonic_hits > 1L ~ "Ambiguous UTR/CDS",
            in_five_utr ~ "5' UTR",
            in_cds ~ "CDS",
            in_three_utr ~ "3' UTR",
            in_exon ~ "Other exon/noncoding",
            in_transcript ~ "Intron",
            TRUE ~ "Intergenic"
        )
        simplified_region <- case_when(
            transcript_region_detail == "5' UTR" ~ "5' UTR",
            transcript_region_detail == "CDS" ~ "CDS",
            transcript_region_detail == "3' UTR" ~ "3' UTR",
            transcript_region_detail == "Intron" ~ "Intron",
            TRUE ~ "Other"
        )
        summit_table %>%
            mutate(
                transcript_region_detail = factor(
                    transcript_region_detail, levels = detail_region_levels
                ),
                transcript_region = factor(
                    simplified_region, levels = region_levels
                ),
                overlaps_5UTR = in_five_utr,
                overlaps_CDS = in_cds,
                overlaps_3UTR = in_three_utr,
                overlaps_any_exon = in_exon,
                overlaps_transcript_span = in_transcript,
                annotation_source =
                    "TxDb.Hsapiens.UCSC.hg38.knownGene 3.20.0",
                annotation_record_id = .data[[id_column]]
            )
    }

    annotation_table <- annotate_summits(
        representative_summits, "peak_id"
    ) %>%
        dplyr::select(-annotation_record_id)
    source_summits <- peak_members %>%
        transmute(
            source_line,
            source_peak_id,
            source_score,
            source_gene_id = gene_id,
            summit_chromosome = chromosome,
            summit_start_0based = as.integer(summit_start),
            summit_end_0based_exclusive = as.integer(summit_end),
            summit_strand = strand
        )
    source_annotation_table <- annotate_summits(
        source_summits, "source_line"
    ) %>%
        dplyr::select(-annotation_record_id)
    write_table(annotation_table, annotation_cache_path)
    write_table(source_annotation_table, source_annotation_cache_path)
    message("Built and cached transcript-region annotations")
}

analysis_peaks <- lfc_peaks %>%
    left_join(annotation_table, by = "peak_id") %>%
    mutate(
        delta_LFC_324_minus_328 =
            `324WT_replicate_mean_LFC` - `328WT_replicate_mean_LFC`,
        strict_quadrant_group = case_when(
            `324WT_replicate_mean_LFC` > strict_lfc_threshold &
                `328WT_replicate_mean_LFC` < -strict_lfc_threshold ~
                "324-biased",
            `328WT_replicate_mean_LFC` > strict_lfc_threshold &
                `324WT_replicate_mean_LFC` < -strict_lfc_threshold ~
                "328-biased",
            TRUE ~ "Other"
        ),
        delta_group = case_when(
            delta_LFC_324_minus_328 >= delta_lfc_threshold &
                `324WT_replicate_mean_LFC` > 0 ~
                "324-biased",
            delta_LFC_324_minus_328 <= -delta_lfc_threshold &
                `328WT_replicate_mean_LFC` > 0 ~
                "328-biased",
            abs(delta_LFC_324_minus_328) < delta_lfc_threshold &
                `324WT_replicate_mean_LFC` > 0 &
                `328WT_replicate_mean_LFC` > 0 ~
                "Similar-positive",
            TRUE ~ "Other"
        ),
        retained_for_continuous_delta =
            (delta_LFC_324_minus_328 >= 0 &
                `324WT_replicate_mean_LFC` > 0) |
            (delta_LFC_324_minus_328 <= 0 &
                `328WT_replicate_mean_LFC` > 0)
    )
if (any(is.na(analysis_peaks$transcript_region))) {
    stop("Transcript-region annotation is missing for one or more peaks")
}
write_table(
    analysis_peaks,
    file.path(output_dir, "construct_bias_peak_groups.tsv")
)

# -----------------------------------------------------------------------------
# Selection scatter
# -----------------------------------------------------------------------------

preference_colors <- c(
    "324-biased" = "#1B9E77",
    "328-biased" = "#E66101",
    "Similar-positive" = "#7570B3"
)
strict_counts <- analysis_peaks %>%
    count(strict_quadrant_group, name = "n")
strict_label <- setNames(
    paste0(
        strict_counts$strict_quadrant_group,
        " (n=", strict_counts$n, ")"
    ),
    strict_counts$strict_quadrant_group
)
scatter_limit <- ceiling(max(abs(c(
    analysis_peaks$`328WT_replicate_mean_LFC`,
    analysis_peaks$`324WT_replicate_mean_LFC`
))))

selection_scatter <- ggplot(
    analysis_peaks,
    aes(
        x = `328WT_replicate_mean_LFC`,
        y = `324WT_replicate_mean_LFC`
    )
) +
    geom_hline(yintercept = 0, color = "grey65", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey65", linewidth = 0.3) +
    geom_hline(
        yintercept = c(-strict_lfc_threshold, strict_lfc_threshold),
        color = "grey55", linewidth = 0.3, linetype = "dotted"
    ) +
    geom_vline(
        xintercept = c(-strict_lfc_threshold, strict_lfc_threshold),
        color = "grey55", linewidth = 0.3, linetype = "dotted"
    ) +
    geom_abline(
        slope = 1, intercept = 0, color = "grey40", linewidth = 0.45,
        linetype = "dashed"
    ) +
    geom_point(color = "#BDBDBD", alpha = 0.22, size = 0.65) +
    geom_point(
        data = analysis_peaks %>%
            filter(strict_quadrant_group != "Other") %>%
            mutate(
                strict_quadrant_group = factor(
                    strict_quadrant_group,
                    levels = c("324-biased", "328-biased")
                )
            ),
        aes(color = strict_quadrant_group),
        alpha = 0.75, size = 1.15
    ) +
    scale_color_manual(
        values = preference_colors,
        labels = strict_label[c("324-biased", "328-biased")],
        drop = FALSE
    ) +
    coord_equal(
        xlim = c(-scatter_limit, scatter_limit),
        ylim = c(-scatter_limit, scatter_limit),
        expand = FALSE
    ) +
    labs(
        title = "Strict opposite-enrichment peaks",
        subtitle = paste0(
            "Selected beyond +/-", strict_lfc_threshold,
            " on both axes; all 16,046 counted peaks shown"
        ),
        x = "328WT log2(CLIP / SMI)",
        y = "324WT log2(CLIP / SMI)",
        color = "Exploratory group",
        caption = paste0(
            "Upper-left peaks are 324-biased; lower-right peaks are 328-biased. ",
            "Cutoffs are descriptive, not statistical tests."
        )
    ) +
    theme_publication()
save_plot(
    selection_scatter,
    file.path(output_dir, "strict_quadrant_peak_selection"),
    width = 7.2, height = 6.5
)

# The requested delta grouping includes construct-biased peaks from every
# quadrant, provided that the preferred construct has positive enrichment.
# Similar-positive peaks lie within the delta band and are enriched in both.
delta_counts <- analysis_peaks %>%
    count(delta_group, name = "n")
delta_label <- setNames(
    paste0(delta_counts$delta_group, " (n=", delta_counts$n, ")"),
    delta_counts$delta_group
)
delta_selection_scatter <- ggplot(
    analysis_peaks,
    aes(
        x = `328WT_replicate_mean_LFC`,
        y = `324WT_replicate_mean_LFC`
    )
) +
    geom_hline(yintercept = 0, color = "grey65", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey65", linewidth = 0.3) +
    geom_abline(
        slope = 1, intercept = 0, color = "grey40", linewidth = 0.45,
        linetype = "dashed"
    ) +
    geom_abline(
        slope = 1, intercept = c(-delta_lfc_threshold, delta_lfc_threshold),
        color = "grey55", linewidth = 0.3, linetype = "dotted"
    ) +
    geom_point(color = "#BDBDBD", alpha = 0.2, size = 0.65) +
    geom_point(
        data = analysis_peaks %>%
            filter(delta_group != "Other") %>%
            mutate(
                delta_group = factor(
                    delta_group,
                    levels = c(
                        "324-biased", "328-biased", "Similar-positive"
                    )
                )
            ),
        aes(color = delta_group),
        alpha = 0.68, size = 0.95
    ) +
    scale_color_manual(
        values = preference_colors,
        labels = delta_label[
            c("324-biased", "328-biased", "Similar-positive")
        ],
        drop = FALSE
    ) +
    coord_equal(
        xlim = c(-scatter_limit, scatter_limit),
        ylim = c(-scatter_limit, scatter_limit),
        expand = FALSE
    ) +
    labs(
        title = "Delta-enrichment peak groups",
        subtitle = paste0(
            "Biased: |delta| >= ", delta_lfc_threshold,
            " and preferred LFC > 0\nSimilar-positive: |delta| < ",
            delta_lfc_threshold, " and both LFCs > 0"
        ),
        x = "328WT log2(CLIP / SMI)",
        y = "324WT log2(CLIP / SMI)",
        color = "Delta group",
        caption = paste0(
            "Parallel dotted lines mark delta = +/-", delta_lfc_threshold,
            ". Similar-positive peaks lie inside the band and above both axes."
        )
    ) +
    theme_publication()
save_plot(
    delta_selection_scatter,
    file.path(output_dir, "delta_group_peak_selection"),
    width = 7.2, height = 6.5
)

# -----------------------------------------------------------------------------
# Five-category transcript-region backgrounds and delta-group enrichment
# -----------------------------------------------------------------------------

delta_group_levels <- c("324-biased", "328-biased", "Similar-positive")
selected_delta_peaks <- analysis_peaks %>%
    filter(delta_group %in% delta_group_levels) %>%
    mutate(
        delta_group = factor(delta_group, levels = delta_group_levels),
        transcript_region = factor(transcript_region, levels = region_levels)
    )

make_background_counts <- function(data, background_id, background_label,
                                   peak_unit) {
    observed <- data %>%
        mutate(transcript_region = as.character(transcript_region)) %>%
        count(transcript_region, name = "background_region_n")
    crossing(transcript_region = region_levels) %>%
        left_join(observed, by = "transcript_region") %>%
        mutate(
            background_region_n = replace_na(background_region_n, 0L),
            background_id = background_id,
            background_label = background_label,
            peak_unit = peak_unit,
            background_total = nrow(data),
            background_fraction = background_region_n / background_total,
            background_percent = 100 * background_fraction,
            transcript_region = factor(transcript_region, levels = region_levels)
        ) %>%
        arrange(match(as.character(transcript_region), region_levels)) %>%
        dplyr::select(
            background_id, background_label, peak_unit, transcript_region,
            background_region_n, background_total, background_fraction,
            background_percent
        )
}

background_counts <- bind_rows(
    make_background_counts(
        source_annotation_table,
        "source_score_ge1",
        "46,108 source BED peaks (score >= 1)",
        "unmerged source peaks"
    ),
    make_background_counts(
        annotation_table,
        "merged_score_ge1",
        "43,794 merged peak regions",
        "merged peak regions"
    ),
    make_background_counts(
        analysis_peaks,
        "counted_nonzero",
        "16,046 merged peaks with any count",
        "merged peak regions"
    )
)
write_table(
    background_counts,
    file.path(output_dir, "transcript_region_background.tsv")
)

selected_totals <- selected_delta_peaks %>%
    count(delta_group, name = "selected_total", .drop = FALSE)
selected_counts <- selected_delta_peaks %>%
    count(delta_group, transcript_region, name = "selected_region_n", .drop = FALSE)

enrichment_table <- crossing(
    delta_group = delta_group_levels,
    transcript_region = region_levels
) %>%
    mutate(
        delta_group = factor(delta_group, levels = delta_group_levels),
        transcript_region = factor(transcript_region, levels = region_levels)
    ) %>%
    left_join(selected_totals, by = "delta_group") %>%
    left_join(selected_counts, by = c("delta_group", "transcript_region")) %>%
    mutate(selected_region_n = replace_na(selected_region_n, 0L)) %>%
    crossing(background_id = unique(background_counts$background_id)) %>%
    left_join(
        background_counts,
        by = c("background_id", "transcript_region")
    ) %>%
    mutate(
        selected_fraction = selected_region_n / selected_total,
        selected_percent = 100 * selected_fraction,
        fraction_fold_enrichment = selected_fraction / background_fraction,
        log2_fraction_enrichment = log2(fraction_fold_enrichment),
        percentage_point_difference = selected_percent - background_percent
    ) %>%
    arrange(
        match(as.character(background_id), c(
            "source_score_ge1", "merged_score_ge1", "counted_nonzero"
        )),
        match(as.character(delta_group), delta_group_levels),
        match(as.character(transcript_region), region_levels)
    )
write_table(
    enrichment_table,
    file.path(output_dir, "transcript_region_enrichment.tsv")
)

# -----------------------------------------------------------------------------
# Omnibus distribution tests and direct construct comparisons
# -----------------------------------------------------------------------------

chisq_summary <- function(contingency, comparison, groups) {
    asymptotic <- suppressWarnings(chisq.test(contingency, correct = FALSE))
    simulated <- suppressWarnings(chisq.test(
        contingency,
        simulate.p.value = TRUE,
        B = chisq_monte_carlo_B
    ))
    statistic <- unname(asymptotic$statistic)
    denominator <- sum(contingency) * (min(dim(contingency)) - 1)
    tibble(
        comparison = comparison,
        groups = paste(groups, collapse = " vs "),
        group_count = nrow(contingency),
        group_totals = paste(rowSums(contingency), collapse = ";"),
        region_category_count = ncol(contingency),
        pearson_chisq = statistic,
        degrees_of_freedom = unname(asymptotic$parameter),
        asymptotic_p_value = asymptotic$p.value,
        monte_carlo_p_value = simulated$p.value,
        monte_carlo_B = chisq_monte_carlo_B,
        min_expected_count = min(asymptotic$expected),
        expected_cells_below_5 = sum(asymptotic$expected < 5),
        cramers_v = sqrt(statistic / denominator)
    )
}

set.seed(chisq_seed)

# Primary matched-unit tests: each selected group versus all remaining merged
# peak regions. This retains the exact 43,794-row unit used by the LFC table.
merged_background_tests <- lapply(delta_group_levels, function(group_i) {
    test_data <- enrichment_table %>%
        filter(delta_group == group_i, background_id == "merged_score_ge1") %>%
        arrange(match(as.character(transcript_region), region_levels))
    row_tests <- lapply(seq_len(nrow(test_data)), function(i) {
        selected_region_n <- test_data$selected_region_n[[i]]
        selected_total <- test_data$selected_total[[i]]
        background_region_n <- test_data$background_region_n[[i]]
        background_total <- test_data$background_total[[i]]
        remaining_region_n <- background_region_n - selected_region_n
        remaining_not_region <-
            background_total - selected_total - remaining_region_n
        test <- fisher.test(matrix(
            c(
                selected_region_n,
                selected_total - selected_region_n,
                remaining_region_n,
                remaining_not_region
            ),
            nrow = 2, byrow = TRUE
        ))
        tibble(
            delta_group = group_i,
            transcript_region = test_data$transcript_region[[i]],
            selected_region_n = selected_region_n,
            selected_total = selected_total,
            selected_percent = test_data$selected_percent[[i]],
            background_region_n = background_region_n,
            background_total = background_total,
            background_percent = test_data$background_percent[[i]],
            percentage_point_difference =
                test_data$percentage_point_difference[[i]],
            fraction_fold_enrichment =
                test_data$fraction_fold_enrichment[[i]],
            log2_fraction_enrichment =
                test_data$log2_fraction_enrichment[[i]],
            fisher_odds_ratio = unname(test$estimate),
            fisher_p_value = test$p.value
        )
    }) %>%
        bind_rows() %>%
        mutate(fisher_padj_BH = p.adjust(fisher_p_value, method = "BH"))
    row_tests
}) %>%
    bind_rows()
write_table(
    merged_background_tests,
    file.path(output_dir, "delta_group_vs_merged_background_tests.tsv")
)

merged_overall_tests <- lapply(delta_group_levels, function(group_i) {
    test_data <- enrichment_table %>%
        filter(delta_group == group_i, background_id == "merged_score_ge1") %>%
        arrange(match(as.character(transcript_region), region_levels))
    contingency <- rbind(
        selected = test_data$selected_region_n,
        remaining_merged_peaks =
            test_data$background_region_n - test_data$selected_region_n
    )
    colnames(contingency) <- as.character(test_data$transcript_region)
    chisq_summary(
        contingency,
        "selected delta group vs remaining 43,794 merged peaks",
        c(group_i, "remaining merged peaks")
    ) %>%
        mutate(delta_group = group_i, .before = 1)
}) %>%
    bind_rows()
write_table(
    merged_overall_tests,
    file.path(output_dir, "overall_region_distribution_chisquare.tsv")
)

# Requested 46,108-source-peak sensitivity analysis. These tests treat the
# source-region proportions as fixed expected probabilities. They are labeled
# sensitivity analyses because selected rows are merged peaks, not source rows.
source_background_tests <- lapply(delta_group_levels, function(group_i) {
    test_data <- enrichment_table %>%
        filter(delta_group == group_i, background_id == "source_score_ge1") %>%
        arrange(match(as.character(transcript_region), region_levels))
    row_tests <- lapply(seq_len(nrow(test_data)), function(i) {
        test <- binom.test(
            test_data$selected_region_n[[i]],
            test_data$selected_total[[i]],
            p = test_data$background_fraction[[i]]
        )
        tibble(
            delta_group = group_i,
            transcript_region = test_data$transcript_region[[i]],
            selected_region_n = test_data$selected_region_n[[i]],
            selected_total = test_data$selected_total[[i]],
            selected_percent = test_data$selected_percent[[i]],
            background_region_n = test_data$background_region_n[[i]],
            background_total = test_data$background_total[[i]],
            background_percent = test_data$background_percent[[i]],
            percentage_point_difference =
                test_data$percentage_point_difference[[i]],
            fraction_fold_enrichment =
                test_data$fraction_fold_enrichment[[i]],
            log2_fraction_enrichment =
                test_data$log2_fraction_enrichment[[i]],
            binomial_p_value = test$p.value,
            binomial_conf_low_percent = 100 * test$conf.int[[1]],
            binomial_conf_high_percent = 100 * test$conf.int[[2]]
        )
    }) %>%
        bind_rows() %>%
        mutate(binomial_padj_BH = p.adjust(binomial_p_value, method = "BH"))
    row_tests
}) %>%
    bind_rows()
write_table(
    source_background_tests,
    file.path(output_dir, "delta_group_vs_source_background_tests.tsv")
)

source_overall_tests <- lapply(delta_group_levels, function(group_i) {
    test_data <- enrichment_table %>%
        filter(delta_group == group_i, background_id == "source_score_ge1") %>%
        arrange(match(as.character(transcript_region), region_levels))
    asymptotic <- chisq.test(
        x = test_data$selected_region_n,
        p = test_data$background_fraction,
        rescale.p = TRUE
    )
    simulated <- chisq.test(
        x = test_data$selected_region_n,
        p = test_data$background_fraction,
        rescale.p = TRUE,
        simulate.p.value = TRUE,
        B = chisq_monte_carlo_B
    )
    tibble(
        delta_group = group_i,
        comparison =
            "selected merged peaks vs fixed 46,108 source-peak proportions",
        selected_total = sum(test_data$selected_region_n),
        region_category_count = nrow(test_data),
        pearson_chisq = unname(asymptotic$statistic),
        degrees_of_freedom = unname(asymptotic$parameter),
        asymptotic_p_value = asymptotic$p.value,
        monte_carlo_p_value = simulated$p.value,
        monte_carlo_B = chisq_monte_carlo_B,
        min_expected_count = min(asymptotic$expected),
        expected_cells_below_5 = sum(asymptotic$expected < 5),
        cohen_w = sqrt(unname(asymptotic$statistic) /
            sum(test_data$selected_region_n))
    )
}) %>%
    bind_rows()
write_table(
    source_overall_tests,
    file.path(output_dir, "source_background_goodness_of_fit_chisquare.tsv")
)

# Direct comparison among the three mutually exclusive delta groups.
direct_counts <- crossing(
    delta_group = delta_group_levels,
    transcript_region = region_levels
) %>%
    left_join(
        selected_counts %>%
            mutate(
                delta_group = as.character(delta_group),
                transcript_region = as.character(transcript_region)
            ),
        by = c("delta_group", "transcript_region")
    ) %>%
    mutate(
        selected_region_n = replace_na(selected_region_n, 0L),
        delta_group = factor(delta_group, levels = delta_group_levels),
        transcript_region = factor(transcript_region, levels = region_levels)
    )
direct_contingency <- xtabs(
    selected_region_n ~ delta_group + transcript_region,
    data = direct_counts
)
three_group_overall <- chisq_summary(
    direct_contingency,
    "324-biased vs 328-biased vs Similar-positive",
    delta_group_levels
)
write_table(
    three_group_overall,
    file.path(output_dir, "three_group_overall_chisquare.tsv")
)

group_pairs <- combn(delta_group_levels, 2, simplify = FALSE)
pairwise_overall <- lapply(group_pairs, function(pair_i) {
    contingency <- direct_contingency[pair_i, , drop = FALSE]
    chisq_summary(
        contingency,
        paste(pair_i, collapse = " vs "),
        pair_i
    ) %>%
        mutate(
            group_1 = pair_i[[1]],
            group_2 = pair_i[[2]],
            .before = 1
        )
}) %>%
    bind_rows() %>%
    mutate(
        asymptotic_padj_BH = p.adjust(asymptotic_p_value, method = "BH"),
        monte_carlo_padj_BH = p.adjust(monte_carlo_p_value, method = "BH")
    )
write_table(
    pairwise_overall,
    file.path(output_dir, "pairwise_group_overall_chisquare.tsv")
)

pairwise_region_tests <- lapply(group_pairs, function(pair_i) {
    totals <- selected_totals %>%
        filter(delta_group %in% pair_i) %>%
        mutate(delta_group = as.character(delta_group))
    lapply(region_levels, function(region_i) {
        counts <- direct_counts %>%
            filter(
                delta_group %in% pair_i,
                transcript_region == region_i
            ) %>%
            mutate(delta_group = as.character(delta_group))
        n_1 <- counts$selected_region_n[counts$delta_group == pair_i[[1]]]
        n_2 <- counts$selected_region_n[counts$delta_group == pair_i[[2]]]
        total_1 <- totals$selected_total[totals$delta_group == pair_i[[1]]]
        total_2 <- totals$selected_total[totals$delta_group == pair_i[[2]]]
        test <- fisher.test(matrix(
            c(n_1, total_1 - n_1, n_2, total_2 - n_2),
            nrow = 2, byrow = TRUE
        ))
        fraction_1 <- n_1 / total_1
        fraction_2 <- n_2 / total_2
        tibble(
            group_1 = pair_i[[1]],
            group_2 = pair_i[[2]],
            transcript_region = region_i,
            group_1_region_n = n_1,
            group_1_total = total_1,
            group_1_percent = 100 * fraction_1,
            group_2_region_n = n_2,
            group_2_total = total_2,
            group_2_percent = 100 * fraction_2,
            percentage_point_difference_group_1_minus_group_2 =
                100 * (fraction_1 - fraction_2),
            fold_group_1_vs_group_2 = fraction_1 / fraction_2,
            log2_fold_group_1_vs_group_2 = log2(fraction_1 / fraction_2),
            fisher_odds_ratio = unname(test$estimate),
            fisher_p_value = test$p.value
        )
    }) %>%
        bind_rows() %>%
        mutate(fisher_padj_BH = p.adjust(fisher_p_value, method = "BH"))
}) %>%
    bind_rows()
write_table(
    pairwise_region_tests,
    file.path(output_dir, "pairwise_group_region_fisher.tsv")
)

# Keep legacy two-group filenames as exact subsets for downstream references.
write_table(
    pairwise_overall %>%
        filter(group_1 == "324-biased", group_2 == "328-biased"),
    file.path(output_dir, "direct_324_vs_328_overall_chisquare.tsv")
)
write_table(
    pairwise_region_tests %>%
        filter(group_1 == "324-biased", group_2 == "328-biased"),
    file.path(output_dir, "direct_324_vs_328_region_fisher.tsv")
)

enrichment_plot_data <- bind_rows(
    merged_background_tests %>%
        transmute(
            background_label = "43,794 merged peaks (matched units)",
            delta_group, transcript_region,
            log2_fraction_enrichment,
            adjusted_p_value = fisher_padj_BH
        ),
    source_background_tests %>%
        transmute(
            background_label = "46,108 source peaks (sensitivity)",
            delta_group, transcript_region,
            log2_fraction_enrichment,
            adjusted_p_value = binomial_padj_BH
        )
) %>%
    mutate(
        background_label = factor(
            background_label,
            levels = c(
                "43,794 merged peaks (matched units)",
                "46,108 source peaks (sensitivity)"
            )
        ),
        delta_group = factor(delta_group, levels = delta_group_levels),
        transcript_region = factor(transcript_region, levels = region_levels),
        significance_label = if_else(adjusted_p_value < 0.05, "*", ""),
        marker_y = log2_fraction_enrichment +
            if_else(log2_fraction_enrichment >= 0, 0.08, -0.08),
        marker_vjust = if_else(log2_fraction_enrichment >= 0, 0, 1)
    )
region_enrichment_plot <- ggplot(
    enrichment_plot_data,
    aes(
        x = transcript_region,
        y = log2_fraction_enrichment,
        fill = delta_group,
        group = delta_group
    )
) +
    geom_hline(yintercept = 0, color = "grey45", linewidth = 0.4) +
    geom_col(position = position_dodge(width = 0.82), width = 0.74) +
    geom_text(
        aes(y = marker_y, label = significance_label, vjust = marker_vjust),
        position = position_dodge(width = 0.82),
        color = "grey15", size = 4, show.legend = FALSE
    ) +
    facet_wrap(~background_label, ncol = 1) +
    scale_fill_manual(values = preference_colors, drop = FALSE) +
    labs(
        title = "Transcript-region enrichment of delta groups",
        subtitle = paste0(
            "Five categories; Similar-positive has |delta| < ",
            delta_lfc_threshold, " and positive enrichment in both constructs"
        ),
        x = NULL,
        y = "log2 region fold enrichment",
        fill = "Delta group",
        caption = paste0(
            "Stars mark BH-adjusted p < 0.05 (Fisher tests for the matched ",
            "merged background; exact binomial sensitivity tests for source ",
            "proportions).\nOther combines all remaining annotations."
        )
    ) +
    theme_publication() +
    theme(
        axis.text.x = element_text(angle = 28, hjust = 1),
        panel.grid.major.x = element_blank()
    )
save_plot(
    region_enrichment_plot,
    file.path(output_dir, "construct_bias_transcript_region_enrichment"),
    width = 8.4, height = 7.4
)

# -----------------------------------------------------------------------------
# Continuous delta-delta after the relaxed positive-preferred filter
# -----------------------------------------------------------------------------

delta_plot_peaks <- analysis_peaks %>%
    filter(retained_for_continuous_delta) %>%
    mutate(
        transcript_region = factor(transcript_region, levels = region_levels),
        delta_delta_324_minus_328 =
            `324WT_replicate_mean_LFC` - `328WT_replicate_mean_LFC`
    )

paired_test_for_region <- function(data) {
    if (nrow(data) < 2L) {
        return(tibble(
            wilcoxon_signed_rank_statistic = NA_real_,
            wilcoxon_p_value = NA_real_
        ))
    }
    test <- suppressWarnings(wilcox.test(
        data$`324WT_replicate_mean_LFC`,
        data$`328WT_replicate_mean_LFC`,
        paired = TRUE,
        exact = FALSE,
        alternative = "two.sided"
    ))
    tibble(
        wilcoxon_signed_rank_statistic = unname(test$statistic),
        wilcoxon_p_value = test$p.value
    )
}

regional_paired_tests <- lapply(region_levels, function(region_i) {
    region_data <- delta_plot_peaks %>%
        filter(transcript_region == region_i)
    paired_test_for_region(region_data) %>%
        mutate(transcript_region = region_i, .before = 1)
}) %>%
    bind_rows() %>%
    mutate(
        wilcoxon_padj_BH = p.adjust(wilcoxon_p_value, method = "BH")
    )

delta_summary <- delta_plot_peaks %>%
    group_by(transcript_region) %>%
    summarise(
        n = n(),
        median_324_LFC = median(`324WT_replicate_mean_LFC`),
        median_328_LFC = median(`328WT_replicate_mean_LFC`),
        mean_delta_delta = mean(delta_delta_324_minus_328),
        median_delta_delta = median(delta_delta_324_minus_328),
        q1_delta_delta = quantile(delta_delta_324_minus_328, 0.25),
        q3_delta_delta = quantile(delta_delta_324_minus_328, 0.75),
        percent_delta_positive = 100 * mean(delta_delta_324_minus_328 > 0),
        annotation_y = quantile(delta_delta_324_minus_328, 0.98) + 0.65,
        .groups = "drop"
    ) %>%
    left_join(regional_paired_tests, by = "transcript_region") %>%
    mutate(
        paired_test_label = case_when(
            is.na(wilcoxon_padj_BH) ~ "not tested",
            wilcoxon_padj_BH == 0 ~ "BH q < 1e-300",
            wilcoxon_padj_BH < 0.001 ~
                sprintf("BH q = %.1e", wilcoxon_padj_BH),
            TRUE ~ sprintf("BH q = %.3f", wilcoxon_padj_BH)
        )
    )
write_table(
    delta_summary,
    file.path(output_dir, "delta_LFC_by_transcript_region.tsv")
)
write_table(
    delta_plot_peaks %>%
        dplyr::select(
            peak_id, chromosome, start_0based, end_0based_exclusive,
            strand, gene_ids, transcript_region, transcript_region_detail,
            `324WT_replicate_mean_LFC`, `328WT_replicate_mean_LFC`,
            delta_delta_324_minus_328,
            retained_for_continuous_delta
        ),
    file.path(output_dir, "relaxed_filter_continuous_delta_delta_peaks.tsv")
)

overall_paired_test <- suppressWarnings(wilcox.test(
    delta_plot_peaks$`324WT_replicate_mean_LFC`,
    delta_plot_peaks$`328WT_replicate_mean_LFC`,
    paired = TRUE,
    exact = FALSE,
    alternative = "two.sided"
))
overall_paired_summary <- tibble(
    comparison = paste0(
        "324WT LFC vs 328WT LFC after retaining (delta >= 0 and 324WT ",
        "LFC > 0) or (delta <= 0 and 328WT LFC > 0)"
    ),
    peak_count = nrow(delta_plot_peaks),
    median_324_LFC = median(delta_plot_peaks$`324WT_replicate_mean_LFC`),
    median_328_LFC = median(delta_plot_peaks$`328WT_replicate_mean_LFC`),
    median_delta_delta = median(delta_plot_peaks$delta_delta_324_minus_328),
    wilcoxon_signed_rank_statistic = unname(overall_paired_test$statistic),
    wilcoxon_p_value = overall_paired_test$p.value
)
write_table(
    overall_paired_summary,
    file.path(output_dir, "relaxed_filter_overall_paired_test.tsv")
)

region_axis_labels <- setNames(
    paste0(delta_summary$transcript_region, "\n(n=", delta_summary$n, ")"),
    as.character(delta_summary$transcript_region)
)
density_regions <- delta_summary %>%
    filter(n >= 2L) %>%
    pull(transcript_region) %>%
    as.character()

delta_plot <- ggplot(
    delta_plot_peaks,
    aes(x = transcript_region, y = delta_delta_324_minus_328)
) +
    geom_hline(yintercept = 0, color = "grey45", linewidth = 0.4) +
    geom_violin(
        data = delta_plot_peaks %>%
            filter(transcript_region %in% density_regions),
        fill = "grey82", color = "grey45", linewidth = 0.35,
        scale = "width", trim = TRUE
    ) +
    geom_boxplot(
        data = delta_plot_peaks %>%
            filter(transcript_region %in% density_regions),
        width = 0.14, outlier.shape = NA, fill = "white",
        color = "grey20", linewidth = 0.4
    ) +
    geom_point(
        data = delta_plot_peaks %>%
            filter(!transcript_region %in% density_regions),
        color = "grey20", size = 2
    ) +
    geom_text(
        data = delta_summary,
        aes(
            x = transcript_region,
            y = annotation_y,
            label = paired_test_label
        ),
        inherit.aes = FALSE,
        size = 3.1,
        color = "grey20"
    ) +
    scale_x_discrete(limits = region_levels, labels = region_axis_labels) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.13))) +
    labs(
        title = "Continuous construct-preference delta by transcript region",
        subtitle = paste0(
            "324WT log2(CLIP/SMI) - 328WT log2(CLIP/SMI); ",
            "CPM pseudocount = ", lfc_pseudocount_cpm,
            "; paired Wilcoxon BH q shown"
        ),
        x = NULL,
        y = "Delta-delta: 324WT LFC - 328WT LFC",
        caption = paste0(
            "Retained ", nrow(delta_plot_peaks),
            " peaks satisfying (delta >= 0 and 324WT LFC > 0) or ",
            "(delta <= 0 and 328WT LFC > 0).\nPeaks depleted versus SMI in ",
            "both constructs were excluded. Positive delta-delta favors 324WT; ",
            "negative favors 328WT.\nOverall paired Wilcoxon p = ",
            format(overall_paired_test$p.value, digits = 2, scientific = TRUE),
            "."
        )
    ) +
    theme_publication() +
    theme(
        axis.text.x = element_text(angle = 28, hjust = 1),
        panel.grid.major.x = element_blank()
    )
save_plot(
    delta_plot,
    file.path(output_dir, "delta_LFC_by_transcript_region"),
    width = 8.4, height = 6.2
)

selection_summary <- bind_rows(
    analysis_peaks %>%
        count(strict_quadrant_group, name = "peak_count") %>%
        transmute(
            selection_rule = "strict_opposite_quadrants",
            group = strict_quadrant_group,
            peak_count
        ),
    analysis_peaks %>%
        count(delta_group, name = "peak_count") %>%
        transmute(
            selection_rule = "delta_biased_positive_preferred",
            group = delta_group,
            peak_count
        )
)
write_table(selection_summary, file.path(output_dir, "selection_summary.tsv"))

parameter_table <- tibble(
    parameter = c(
        "analysis_universe",
        "representative_position",
        "genome_build",
        "annotation_database",
        "strand_aware_annotation",
        "ambiguous_annotation_rule",
        "strict_324_rule",
        "strict_328_rule",
        "delta_324_rule",
        "delta_328_rule",
        "delta_similar_positive_rule",
        "delta_definition",
        "simplified_region_categories",
        "primary_enrichment_background",
        "primary_region_test",
        "source_sensitivity_background",
        "source_sensitivity_test",
        "direct_three_group_test",
        "direct_pairwise_test",
        "LFC_pseudocount_CPM",
        "delta_delta_plot_filter",
        "delta_delta_definition",
        "delta_delta_paired_test",
        "chisquare_monte_carlo_replicates",
        "chisquare_random_seed"
    ),
    value = c(
        "16,046 merged peaks with any count across eight libraries",
        "summit of highest-scoring source omniCLIP peak in each merged peak",
        "hg38",
        "TxDb.Hsapiens.UCSC.hg38.knownGene 3.20.0",
        "yes",
        "summit overlapping more than one of 5UTR/CDS/3UTR",
        "324 LFC > 1 and 328 LFC < -1",
        "328 LFC > 1 and 324 LFC < -1",
        "delta >= 1 and 324 LFC > 0",
        "delta <= -1 and 328 LFC > 0",
        "absolute delta < 1 and both 324 and 328 LFC > 0",
        "324WT replicate-mean LFC minus 328WT replicate-mean LFC",
        "5' UTR, CDS, 3' UTR, Intron, Other",
        "all 43,794 merged peak regions; same unit as selected peaks",
        "selected group vs remaining merged peaks: Pearson chi-square overall and two-sided Fisher per region with BH across five regions",
        "all 46,108 unmerged source BED peaks with score >= 1",
        "fixed source proportions: Pearson goodness-of-fit overall and exact binomial per region with BH across five regions; sensitivity analysis because peak units differ",
        paste0(
            "Pearson chi-square across 324-biased, 328-biased, and Similar-positive; Monte Carlo p with B=",
            chisq_monte_carlo_B
        ),
        "all three pairwise Pearson chi-square tests with BH; per-region two-sided Fisher tests with BH across five regions per pair",
        lfc_pseudocount_cpm,
        "retain peaks satisfying (delta >= 0 and 324WT LFC > 0) or (delta <= 0 and 328WT LFC > 0); no discrete group is assigned",
        "log2((mean 324WT CLIP CPM + 0.5)/(mean 324WT SMI CPM + 0.5)) minus log2((mean 328WT CLIP CPM + 0.5)/(mean 328WT SMI CPM + 0.5))",
        "paired two-sided Wilcoxon signed-rank test within each transcript region; BH across five regions; overall paired test also reported",
        chisq_monte_carlo_B,
        chisq_seed
    )
)
write_table(parameter_table, file.path(output_dir, "analysis_parameters.tsv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

message(
    "Annotated and compared transcript regions for ", nrow(analysis_peaks),
    " peaks; results written to ", output_dir
)
