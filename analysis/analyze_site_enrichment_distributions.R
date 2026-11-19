suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(patchwork)
    library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default) {
    i <- match(flag, args)
    if (is.na(i)) return(default)
    if (i == length(args)) stop("Missing value after ", flag)
    args[[i + 1L]]
}

input_dir <- arg_value("--input", "data/results")
peak_dir <- arg_value("--peak-input", "data/results/peak_analysis")
output_dir <- arg_value(
    "--output", "output/splicing_site_enrichment_distributions"
)
pseudocount_cpm <- as.numeric(arg_value("--pseudocount", "0.1"))
stopifnot(is.finite(pseudocount_cpm), pseudocount_cpm > 0)

results_dir <- file.path(output_dir, "results")
figures_dir <- file.path(output_dir, "figures")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_dir, "qc"), recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
    fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
}
write_tsv_gz <- function(x, path) {
    fwrite(
        as.data.table(x), path, sep = "\t", quote = FALSE,
        na = "NA", compress = "gzip"
    )
}
save_plot <- function(plot, stem, width = 11, height = 8) {
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
theme_publication <- function(base_size = 10) {
    theme_classic(base_size = base_size) +
        theme(
            plot.title = element_text(face = "bold", size = rel(1.1)),
            plot.subtitle = element_text(color = "grey30"),
            plot.caption = element_text(color = "grey35", hjust = 0),
            axis.title = element_text(face = "bold"),
            panel.grid.major.y = element_line(
                color = "grey90", linewidth = 0.25
            ),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = "grey94", color = NA),
            strip.text = element_text(face = "bold"),
            legend.position = "bottom",
            legend.title = element_text(face = "bold")
        )
}

construct_colors <- c(`324` = "#1B9E77", `328` = "#D95F02")
sign_colors <- c(Negative = "#4575B4", Zero = "#BDBDBD", Positive = "#D73027")

metadata <- fread(file.path(input_dir, "sample_metadata.tsv"))
metadata[, construct := as.character(construct)]
metadata[, library_millions := total_primary_mapped_reads / 1e6]
if (nrow(metadata) != 8L || anyDuplicated(metadata$sample_id)) {
    stop("Expected eight unique samples")
}
library_lookup <- setNames(metadata$library_millions, metadata$sample_id)

pair_map <- metadata[, {
    if (sum(assay == "CLIP") != 1L || sum(assay == "SMI") != 1L) {
        stop("Every pair must contain one CLIP and one SMI sample")
    }
    .(
        clip_sample = sample_id[assay == "CLIP"],
        smi_sample = sample_id[assay == "SMI"]
    )
}, by = .(
    construct, replicate = biological_replicate, pair_id
)]

sites_all <- fread(file.path(input_dir, "landmark_sites.tsv.gz"))
sites <- sites_all[
    retained_for_bin_counts == "yes" &
        site_type %in% c("splice_donor_5prime", "splice_acceptor_3prime")
]
sites[, primary_gene_id := sub(";.*$", "", gene_id)]
if (nrow(sites) != 201963L || anyDuplicated(sites$site_id)) {
    stop("Unexpected retained splice-site universe")
}

bin_counts <- fread(file.path(input_dir, "landmark_bin_counts.tsv.gz"))[
    site_type %in% c("splice_donor_5prime", "splice_acceptor_3prime")
]
if (any(bin_counts$raw_overlap_count <= 0L) ||
        anyDuplicated(bin_counts[, .(
            sample_id, site_id, relative_bin_start
        )])) {
    stop("Sparse bin counts must contain unique positive rows")
}
bin_counts[, normalized_count_cpm :=
    raw_overlap_count / library_lookup[sample_id]]

type_config <- list(
    splice_donor_5prime = list(
        label = "5' donor",
        regions = data.table(
            region = c("Exon side", "Intron side", "Full immediate"),
            region_start = c(-50L, 0L, -50L),
            region_end = c(0L, 50L, 50L)
        )
    ),
    splice_acceptor_3prime = list(
        label = "3' acceptor",
        regions = data.table(
            region = c("Intron side", "Exon side", "Full immediate"),
            region_start = c(-50L, 0L, -50L),
            region_end = c(0L, 50L, 50L)
        )
    )
)

pair_bin_enrichment <- function(counts) {
    out <- vector("list", nrow(pair_map))
    for (i in seq_len(nrow(pair_map))) {
        pair <- pair_map[i]
        clip <- counts[sample_id == pair$clip_sample, .(
            site_id, relative_bin_start, relative_bin_end,
            clip_cpm = normalized_count_cpm,
            clip_raw = raw_overlap_count
        )]
        smi <- counts[sample_id == pair$smi_sample, .(
            site_id, relative_bin_start, relative_bin_end,
            smi_cpm = normalized_count_cpm,
            smi_raw = raw_overlap_count
        )]
        joined <- merge(
            clip, smi,
            by = c("site_id", "relative_bin_start", "relative_bin_end"),
            all = TRUE
        )
        for (column in c("clip_cpm", "clip_raw", "smi_cpm", "smi_raw")) {
            set(joined, which(is.na(joined[[column]])), column, 0)
        }
        joined[, enrichment := log2(
            (clip_cpm + pseudocount_cpm) /
                (smi_cpm + pseudocount_cpm)
        )]
        joined[, `:=`(
            construct = pair$construct,
            replicate = pair$replicate,
            pair_id = pair$pair_id
        )]
        out[[i]] <- joined
    }
    rbindlist(out)
}

wilson_interval <- function(successes, n, level = 0.95) {
    if (n == 0) return(c(NA_real_, NA_real_))
    z <- qnorm(1 - (1 - level) / 2)
    p <- successes / n
    denominator <- 1 + z^2 / n
    center <- (p + z^2 / (2 * n)) / denominator
    half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denominator
    c(center - half, center + half)
}

site_values_all <- list()
replicate_values_all <- list()

for (type_name in names(type_config)) {
    config <- type_config[[type_name]]
    type_sites <- sites[site_type == type_name]
    type_counts <- bin_counts[site_type == type_name]
    pair_bins <- pair_bin_enrichment(type_counts)

    region_outputs <- vector("list", nrow(config$regions))
    replicate_outputs <- vector("list", nrow(config$regions))
    for (region_index in seq_len(nrow(config$regions))) {
        region_row <- config$regions[region_index]
        n_bins <- (region_row$region_end - region_row$region_start) / 10
        region_pair <- pair_bins[
            relative_bin_start >= region_row$region_start &
                relative_bin_end <= region_row$region_end,
            .(
                enrichment = sum(enrichment) / n_bins,
                clip_raw_pair = sum(clip_raw),
                smi_raw_pair = sum(smi_raw)
            ),
            by = .(site_id, construct, replicate, pair_id)
        ]

        pair_grid <- CJ(
            site_id = type_sites$site_id,
            construct = c("324", "328"),
            replicate = c("A", "B"),
            unique = TRUE
        )
        pair_grid <- merge(
            pair_grid,
            pair_map[, .(construct, replicate, pair_id)],
            by = c("construct", "replicate"), all.x = TRUE
        )
        region_pair <- merge(
            pair_grid, region_pair,
            by = c("site_id", "construct", "replicate", "pair_id"),
            all.x = TRUE
        )
        region_pair[is.na(enrichment), `:=`(
            enrichment = 0,
            clip_raw_pair = 0,
            smi_raw_pair = 0
        )]
        region_pair[, `:=`(
            site_type = type_name,
            landmark = config$label,
            region = region_row$region,
            region_start = region_row$region_start,
            region_end = region_row$region_end,
            measurable_both_pair = clip_raw_pair > 0 & smi_raw_pair > 0
        )]
        replicate_outputs[[region_index]] <- region_pair

        construct_site <- region_pair[, .(
            enrichment = mean(enrichment),
            clip_raw_construct = sum(clip_raw_pair),
            smi_raw_construct = sum(smi_raw_pair),
            measurable_both_same_pair_any = any(
                clip_raw_pair > 0 & smi_raw_pair > 0
            ),
            measurable_both_both_pairs = all(
                clip_raw_pair > 0 & smi_raw_pair > 0
            ),
            positive_in_both_replicates = all(enrichment > 0),
            positive_in_either_replicate = any(enrichment > 0)
        ), by = .(
            site_id, site_type, landmark, region,
            region_start, region_end, construct
        )]
        construct_site[, `:=`(
            measurable_clip = clip_raw_construct > 0,
            measurable_smi = smi_raw_construct > 0,
            measurable_both =
                clip_raw_construct > 0 & smi_raw_construct > 0,
            enrichment_sign = factor(
                fifelse(
                    enrichment > 0, "Positive",
                    fifelse(enrichment < 0, "Negative", "Zero")
                ),
                levels = c("Negative", "Zero", "Positive")
            )
        )]
        construct_site <- merge(
            construct_site,
            type_sites[, .(
                site_id, chromosome, genomic_position, strand,
                gene_id, gene_name, primary_gene_id
            )],
            by = "site_id", all.x = TRUE
        )
        region_outputs[[region_index]] <- construct_site
    }

    site_values_all[[type_name]] <- rbindlist(region_outputs)
    replicate_values_all[[type_name]] <- rbindlist(replicate_outputs)
}

site_values <- rbindlist(site_values_all)
replicate_values <- rbindlist(replicate_values_all)
write_tsv_gz(
    site_values,
    file.path(results_dir, "site_level_splice_enrichment.tsv.gz")
)
write_tsv_gz(
    site_values[enrichment > 0],
    file.path(results_dir, "positive_enriched_splice_sites.tsv.gz")
)

distribution_summary <- site_values[, {
    n_all <- .N
    n_positive <- sum(enrichment > 0)
    n_measurable <- sum(measurable_both)
    positive_ci_all <- wilson_interval(n_positive, n_all)
    positive_measurable <- sum(enrichment > 0 & measurable_both)
    positive_ci_measurable <- wilson_interval(
        positive_measurable, n_measurable
    )
    .(
        N_sites = n_all,
        N_measurable_CLIP = sum(measurable_clip),
        N_measurable_SMI = sum(measurable_smi),
        N_measurable_both = n_measurable,
        fraction_measurable_both = n_measurable / n_all,
        N_negative = sum(enrichment < 0),
        N_zero = sum(enrichment == 0),
        N_positive = n_positive,
        fraction_positive_all = n_positive / n_all,
        fraction_positive_all_ci_low = positive_ci_all[1L],
        fraction_positive_all_ci_high = positive_ci_all[2L],
        N_positive_measurable_both = positive_measurable,
        fraction_positive_measurable_both =
            positive_measurable / n_measurable,
        fraction_positive_measurable_ci_low =
            positive_ci_measurable[1L],
        fraction_positive_measurable_ci_high =
            positive_ci_measurable[2L],
        median_all_sites = median(enrichment),
        mean_all_sites = mean(enrichment),
        median_measurable_both = median(enrichment[measurable_both]),
        mean_measurable_both = mean(enrichment[measurable_both]),
        N_measurable_both_same_pair_any =
            sum(measurable_both_same_pair_any),
        median_measurable_both_same_pair_any = median(
            enrichment[measurable_both_same_pair_any]
        ),
        N_measurable_both_both_pairs = sum(measurable_both_both_pairs),
        median_measurable_both_both_pairs = median(
            enrichment[measurable_both_both_pairs]
        ),
        N_positive_in_both_replicates =
            sum(positive_in_both_replicates),
        fraction_positive_in_both_replicates =
            mean(positive_in_both_replicates),
        fraction_positive_in_either_replicate =
            mean(positive_in_either_replicate)
    )
}, by = .(
    site_type, landmark, region, region_start, region_end, construct
)]
write_tsv(
    distribution_summary,
    file.path(results_dir, "site_enrichment_distribution_summary.tsv")
)

replicate_summary <- replicate_values[, {
    n_positive <- sum(enrichment > 0)
    positive_ci <- wilson_interval(n_positive, .N)
    measurable_n <- sum(measurable_both_pair)
    .(
        N_sites = .N,
        N_measurable_both = measurable_n,
        fraction_positive_all = n_positive / .N,
        fraction_positive_ci_low = positive_ci[1L],
        fraction_positive_ci_high = positive_ci[2L],
        median_all_sites = median(enrichment),
        median_measurable_both =
            median(enrichment[measurable_both_pair]),
        fraction_positive_measurable_both =
            mean(enrichment[measurable_both_pair] > 0)
    )
}, by = .(
    site_type, landmark, region, region_start, region_end,
    construct, replicate, pair_id
)]
write_tsv(
    replicate_summary,
    file.path(results_dir, "replicate_site_enrichment_summary.tsv")
)

# Peak proximity uses all valid annotated splice sites, not the coverage-retained
# subset, to avoid conditioning the spatial null on read coverage.
splice_windows <- sites_all[
    site_type %in% c("splice_donor_5prime", "splice_acceptor_3prime"),
    .(
        chromosome, strand,
        start = pmax(0L, anchor_0based - 50L),
        end = anchor_0based + 50L,
        site_type
    )
]
splice_windows <- unique(splice_windows)
setkey(splice_windows, chromosome, strand, start, end)

peak_coordinates <- fread(file.path(peak_dir, "peaks.filtered_merged.tsv"))
peak_coordinates[, `:=`(
    peak_midpoint_0based = floor(
        (start_0based + end_0based_exclusive - 1) / 2
    ) |> as.integer(),
    primary_gene_id = sub(";.*$", "", gene_ids)
)]
peak_points <- peak_coordinates[, .(
    peak_id, chromosome, strand,
    start = peak_midpoint_0based,
    end = peak_midpoint_0based + 1L
)]
peak_overlap <- foverlaps(
    peak_points, splice_windows,
    by.x = c("chromosome", "strand", "start", "end"),
    by.y = c("chromosome", "strand", "start", "end"),
    type = "any", nomatch = 0L
)
peak_overlap <- unique(peak_overlap[, .(peak_id, site_type)])
peak_flags <- dcast(
    peak_overlap[, .(value = TRUE), by = .(peak_id, site_type)],
    peak_id ~ site_type, value.var = "value", fill = FALSE
)
setnames(
    peak_flags,
    c("splice_donor_5prime", "splice_acceptor_3prime"),
    c("within_50nt_donor", "within_50nt_acceptor"),
    skip_absent = TRUE
)
peak_annotation <- merge(
    peak_coordinates,
    peak_flags,
    by = "peak_id", all.x = TRUE
)
for (column in c("within_50nt_donor", "within_50nt_acceptor")) {
    if (!column %in% names(peak_annotation)) peak_annotation[, (column) := FALSE]
    peak_annotation[is.na(get(column)), (column) := FALSE]
}

tested_tables <- rbindlist(lapply(
    c(`324` = "MSI2_324WT", `328` = "MSI2_328WT"),
    function(directory) {
        table <- fread(file.path(
            peak_dir, "deseq2", directory, "IP_vs_input.all.tsv"
        ))
        table[, construct := sub("MSI2_", "", sub("WT$", "", directory))]
        table
    }
), use.names = TRUE, fill = TRUE)
tested_tables[, significant_enriched :=
    !is.na(padj) & padj < 0.05 & log2FoldChange > 0]
tested_peaks <- merge(
    tested_tables,
    peak_annotation[, .(
        peak_id, primary_gene_id, peak_midpoint_0based,
        within_50nt_donor, within_50nt_acceptor
    )],
    by = "peak_id", all.x = TRUE
)

peak_proximity_summary <- list()
peak_test_summary <- list()
for (construct_value in c("324", "328")) {
    construct_peaks <- tested_peaks[construct == construct_value]
    for (landmark_value in c("Donor", "Acceptor")) {
        flag_column <- if (landmark_value == "Donor") {
            "within_50nt_donor"
        } else {
            "within_50nt_acceptor"
        }
        construct_peaks[, proximal := get(flag_column)]
        for (group_value in c("Significant enriched", "Other tested")) {
            group_rows <- if (group_value == "Significant enriched") {
                construct_peaks[significant_enriched == TRUE]
            } else {
                construct_peaks[significant_enriched == FALSE]
            }
            n_proximal <- sum(group_rows$proximal)
            ci <- wilson_interval(n_proximal, nrow(group_rows))
            peak_proximity_summary[[length(peak_proximity_summary) + 1L]] <-
                data.table(
                    construct = construct_value,
                    landmark = landmark_value,
                    peak_group = group_value,
                    N_peaks = nrow(group_rows),
                    N_within_50nt = n_proximal,
                    fraction_within_50nt = n_proximal / nrow(group_rows),
                    fraction_ci_low = ci[1L],
                    fraction_ci_high = ci[2L]
                )
        }

        a <- sum(construct_peaks$significant_enriched & construct_peaks$proximal)
        b <- sum(construct_peaks$significant_enriched & !construct_peaks$proximal)
        c_value <- sum(!construct_peaks$significant_enriched & construct_peaks$proximal)
        d <- sum(!construct_peaks$significant_enriched & !construct_peaks$proximal)
        fisher <- fisher.test(matrix(c(a, b, c_value, d), nrow = 2L, byrow = TRUE))

        gene_counts <- construct_peaks[, .(
            a = sum(significant_enriched & proximal),
            b = sum(significant_enriched & !proximal),
            c = sum(!significant_enriched & proximal),
            d = sum(!significant_enriched & !proximal)
        ), by = primary_gene_id]
        informative <- gene_counts[
            (a + b) > 0 & (c + d) > 0 &
                (a + c) > 0 & (b + d) > 0
        ]
        if (nrow(informative) >= 2L) {
            strata <- array(NA_real_, dim = c(2L, 2L, nrow(informative)))
            for (stratum_index in seq_len(nrow(informative))) {
                row <- informative[stratum_index]
                strata[, , stratum_index] <- matrix(
                    c(row$a, row$b, row$c, row$d),
                    nrow = 2L, byrow = TRUE
                )
            }
            mh <- mantelhaen.test(strata, correct = FALSE)
            mh_or <- unname(mh$estimate)
            mh_p <- mh$p.value
            mh_low <- unname(mh$conf.int[1L])
            mh_high <- unname(mh$conf.int[2L])
        } else {
            mh_or <- mh_p <- mh_low <- mh_high <- NA_real_
        }
        peak_test_summary[[length(peak_test_summary) + 1L]] <- data.table(
            construct = construct_value,
            landmark = landmark_value,
            enriched_proximal = a,
            enriched_not_proximal = b,
            background_proximal = c_value,
            background_not_proximal = d,
            fisher_odds_ratio = unname(fisher$estimate),
            fisher_ci_low = unname(fisher$conf.int[1L]),
            fisher_ci_high = unname(fisher$conf.int[2L]),
            fisher_p = fisher$p.value,
            gene_stratified_informative_genes = nrow(informative),
            gene_stratified_common_odds_ratio = mh_or,
            gene_stratified_ci_low = mh_low,
            gene_stratified_ci_high = mh_high,
            gene_stratified_p = mh_p
        )
    }
}
peak_proximity_summary <- rbindlist(peak_proximity_summary)
peak_test_summary <- rbindlist(peak_test_summary)
write_tsv(
    tested_peaks,
    file.path(results_dir, "tested_peak_splice_proximity.tsv")
)
write_tsv(
    peak_proximity_summary,
    file.path(results_dir, "peak_splice_proximity_summary.tsv")
)
write_tsv(
    peak_test_summary,
    file.path(results_dir, "peak_splice_proximity_tests.tsv")
)

full_summary <- distribution_summary[region == "Full immediate"]
sign_plot_data <- full_summary[, .(
    construct, landmark,
    sign = factor(
        c("Negative", "Zero", "Positive"),
        levels = c("Negative", "Zero", "Positive")
    ),
    fraction = c(N_negative, N_zero, N_positive) / N_sites
), by = .I]

p_sign <- ggplot(
    sign_plot_data,
    aes(construct, fraction, fill = sign)
) +
    geom_col(width = 0.7) +
    facet_wrap(~ landmark) +
    scale_fill_manual(values = sign_colors) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
        title = "A. Sign of site-level enrichment",
        subtitle = "All retained sites; immediate -50 to +50 nt region",
        x = "Construct", y = "Fraction of sites", fill = "Enrichment"
    ) +
    theme_publication()

positive_plot_data <- rbindlist(list(
    full_summary[, .(
        construct, landmark, denominator = "All retained sites",
        fraction = fraction_positive_all,
        ci_low = fraction_positive_all_ci_low,
        ci_high = fraction_positive_all_ci_high
    )],
    full_summary[, .(
        construct, landmark, denominator = "Both assays measurable",
        fraction = fraction_positive_measurable_both,
        ci_low = fraction_positive_measurable_ci_low,
        ci_high = fraction_positive_measurable_ci_high
    )]
))
p_positive <- ggplot(
    positive_plot_data,
    aes(construct, fraction, color = construct)
) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.12) +
    geom_point(size = 2.5) +
    facet_grid(denominator ~ landmark, scales = "free_y") +
    scale_color_manual(values = construct_colors, guide = "none") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
        title = "B. Sites with positive enrichment",
        subtitle = "Wilson 95% confidence intervals",
        x = "Construct", y = "Fraction with log2(CLIP/SMI) > 0"
    ) +
    theme_publication(base_size = 9)

ecdf_data <- site_values[
    region == "Full immediate" & measurable_both
]
p_ecdf <- ggplot(
    ecdf_data,
    aes(enrichment, color = construct)
) +
    stat_ecdf(linewidth = 0.85) +
    geom_vline(xintercept = 0, color = "grey45", linetype = 2) +
    facet_wrap(~ landmark) +
    scale_color_manual(values = construct_colors) +
    coord_cartesian(xlim = c(-6, 6)) +
    labs(
        title = "C. Distribution among measurable sites",
        subtitle = "ECDF; x-axis limited to -6 to +6 for readability",
        x = "Site-level mean log2(CLIP/SMI)",
        y = "Cumulative fraction", color = "Construct"
    ) +
    theme_publication()

side_positive <- distribution_summary[
    region %in% c("Exon side", "Intron side")
]
p_side <- ggplot(
    side_positive,
    aes(region, fraction_positive_all, color = construct, group = construct)
) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.5) +
    geom_errorbar(
        aes(
            ymin = fraction_positive_all_ci_low,
            ymax = fraction_positive_all_ci_high
        ),
        width = 0.1
    ) +
    facet_wrap(~ landmark, scales = "free_x") +
    scale_color_manual(values = construct_colors) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
        title = "D. Positive enrichment by splice side",
        subtitle = "All retained sites; Wilson 95% confidence intervals",
        x = NULL, y = "Fraction positive", color = "Construct"
    ) +
    theme_publication()

distribution_figure <- (p_sign | p_positive) / (p_ecdf | p_side) +
    plot_annotation(
        title = "Site-level splice-landmark CLIP/SMI enrichment",
        subtitle = paste0(
            "Common primary-mapped CPM; p = ", pseudocount_cpm,
            " CPM; implicit zeros retained"
        ),
        caption = paste0(
            "Site values average 10-nt bin enrichment across the indicated ",
            "region and then average biological replicates."
        ),
        theme = theme(
            plot.title = element_text(face = "bold", size = 14),
            plot.subtitle = element_text(size = 10),
            plot.caption = element_text(size = 8, color = "grey35")
        )
    )
save_plot(
    distribution_figure,
    file.path(figures_dir, "splice_site_enrichment_distributions")
)

p_peak_fraction <- ggplot(
    peak_proximity_summary,
    aes(peak_group, fraction_within_50nt, color = peak_group)
) +
    geom_errorbar(
        aes(ymin = fraction_ci_low, ymax = fraction_ci_high),
        width = 0.12
    ) +
    geom_point(size = 2.7) +
    facet_grid(construct ~ landmark) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_color_manual(values = c(
        `Significant enriched` = "#7B3294",
        `Other tested` = "#737373"
    )) +
    labs(
        title = "A. Peak midpoints near splice boundaries",
        subtitle = "Same-strand midpoint within +/-50 nt; Wilson 95% CI",
        x = NULL, y = "Fraction of tested peaks", color = "Peak group"
    ) +
    theme_publication(base_size = 9) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))

overall_or <- peak_test_summary[, .(
    construct, landmark, test = "Overall Fisher",
    odds_ratio = fisher_odds_ratio,
    ci_low = fisher_ci_low, ci_high = fisher_ci_high,
    p_value = fisher_p
)]
stratified_or <- peak_test_summary[, .(
    construct, landmark, test = "Gene-stratified",
    odds_ratio = gene_stratified_common_odds_ratio,
    ci_low = gene_stratified_ci_low,
    ci_high = gene_stratified_ci_high,
    p_value = gene_stratified_p
)]
or_plot_data <- rbindlist(list(overall_or, stratified_or))
or_plot_data[
    !is.finite(odds_ratio) | odds_ratio <= 0 |
        !is.finite(ci_low) | ci_low <= 0 |
        !is.finite(ci_high) | ci_high <= 0,
    `:=`(odds_ratio = NA_real_, ci_low = NA_real_, ci_high = NA_real_)
]
p_peak_or <- ggplot(
    or_plot_data[!is.na(odds_ratio)],
    aes(test, odds_ratio, color = construct)
) +
    geom_hline(yintercept = 1, color = "grey45", linetype = 2) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.12) +
    geom_point(size = 2.7) +
    facet_wrap(~ landmark) +
    scale_color_manual(values = construct_colors) +
    scale_y_log10() +
    labs(
        title = "B. Enrichment of splice-proximal peaks",
        subtitle = "Odds ratios with 95% confidence intervals",
        x = NULL, y = "Odds ratio (log scale)", color = "Construct"
    ) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))

peak_figure <- (p_peak_fraction | p_peak_or) +
    plot_annotation(
        title = "Proximity of significant CLIP-enriched peaks to splice boundaries",
        subtitle = paste0(
            "Foreground: padj < 0.05 and log2FC > 0; background: other ",
            "DESeq2-tested fixed peaks"
        ),
        caption = paste0(
            "Boundary set uses all valid annotated donors/acceptors, not the ",
            "coverage-retained landmark subset. Gene-stratified estimates ",
            "control for primary gene ID where strata are informative."
        ),
        theme = theme(
            plot.title = element_text(face = "bold", size = 14),
            plot.subtitle = element_text(size = 10),
            plot.caption = element_text(size = 8, color = "grey35")
        )
    )
save_plot(
    peak_figure,
    file.path(figures_dir, "enriched_peak_splice_proximity"),
    width = 12, height = 6.5
)

parameters <- data.table(
    parameter = c(
        "site_universe", "site_level_enrichment",
        "construct_site_summary", "measurable_both_definition",
        "positive_definition", "normalization", "pseudocount_cpm",
        "peak_foreground", "peak_background", "peak_proximity",
        "peak_boundary_universe", "peak_stratified_sensitivity"
    ),
    value = c(
        "retained_for_bin_counts == yes",
        "mean of 10-nt bin log2((CLIP_CPM+p)/(SMI_CPM+p)) within region",
        "mean of biological replicate A and B site-level enrichments",
        "at least one raw regional CLIP count and one raw regional SMI count across construct replicates",
        "construct site-level enrichment > 0",
        "raw_count / primary_mapped_reads * 1e6",
        pseudocount_cpm,
        "DESeq2 padj < 0.05 and log2FoldChange > 0",
        "other peaks in the same construct-specific DESeq2 tested universe",
        "same-strand peak midpoint within 50 nt of a splice boundary",
        "all valid annotated donor or acceptor sites, independent of read coverage",
        "Mantel-Haenszel common odds ratio stratified by primary peak gene ID"
    )
)
write_tsv(parameters, file.path(results_dir, "analysis_parameters.tsv"))
writeLines(
    capture.output(sessionInfo()),
    file.path(results_dir, "qc", "sessionInfo.txt")
)

message("Site enrichment distribution analysis complete: ", normalizePath(output_dir))
