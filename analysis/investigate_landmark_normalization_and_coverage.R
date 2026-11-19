#!/usr/bin/env Rscript

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
old_output_dir <- arg_value("--old-output", "output/splicing_first_pass")
output_dir <- arg_value("--output", "output/splicing_common_cpm_investigation")
bootstrap_iterations <- as.integer(arg_value("--bootstrap", "1000"))
primary_pseudocount <- as.numeric(arg_value("--pseudocount", "0.1"))
sensitivity_pseudocount <- as.numeric(
    arg_value("--sensitivity-pseudocount", "0.5")
)
bootstrap_seed <- as.integer(arg_value("--seed", "328324"))

stopifnot(
    bootstrap_iterations >= 1000L,
    is.finite(primary_pseudocount), primary_pseudocount > 0,
    is.finite(sensitivity_pseudocount), sensitivity_pseudocount > 0
)

required_inputs <- file.path(input_dir, c(
    "sample_metadata.tsv", "landmark_sites.tsv.gz",
    "landmark_bin_counts.tsv.gz", "landmark_site_summary.tsv.gz",
    "landmark_metaprofile_raw.tsv.gz"
))
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
    stop("Missing required inputs: ", paste(missing_inputs, collapse = ", "))
}
if (!dir.exists(old_output_dir)) {
    stop("Previous first-pass output is required for normalization comparison")
}

results_dir <- file.path(output_dir, "results")
figures_dir <- file.path(output_dir, "figures")
result_subdirs <- c(
    "normalization", "metaprofiles", "replicate_regions",
    "side_of_splice", "coverage_dependence", "qc"
)
figure_subdirs <- c(
    "donors", "acceptors", "transcript_3prime_ends",
    "replicate_regions", "side_of_splice", "qc"
)
invisible(lapply(
    file.path(results_dir, result_subdirs), dir.create,
    recursive = TRUE, showWarnings = FALSE
))
invisible(lapply(
    file.path(figures_dir, figure_subdirs), dir.create,
    recursive = TRUE, showWarnings = FALSE
))

write_tsv <- function(x, path) {
    fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
}
write_tsv_gz <- function(x, path) {
    fwrite(
        as.data.table(x), path, sep = "\t", quote = FALSE,
        na = "NA", compress = "gzip"
    )
}
save_plot <- function(plot, stem, width = 7.2, height = 5.0) {
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
            plot.title = element_text(face = "bold", size = rel(1.12)),
            plot.subtitle = element_text(color = "grey30"),
            plot.caption = element_text(color = "grey35", hjust = 0),
            axis.title = element_text(face = "bold"),
            panel.grid.major.y = element_line(color = "grey90", linewidth = 0.25),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = "grey94", color = NA),
            strip.text = element_text(face = "bold"),
            legend.position = "bottom",
            legend.title = element_text(face = "bold")
        )
}

construct_colors <- c(`324` = "#1B9E77", `328` = "#D95F02")
normalization_colors <- c(
    `Previous construct-specific scale` = "#7A7A7A",
    `Common primary-mapped CPM` = "#6A3D9A"
)
weighting_colors <- c(
    `Site-weighted` = "#6A3D9A",
    `Read-weighted aggregate` = "#E69F00"
)
quartile_colors <- c(
    `Q1 lowest` = "#2166AC", `Q2` = "#67A9CF",
    `Q3` = "#EF8A62", `Q4 highest` = "#B2182B"
)

type_config <- list(
    splice_donor_5prime = list(
        slug = "donor", figure_dir = "donors", label = "5' splice donor",
        left = "EXON", right = "INTRON", expected_bins = 50L,
        regions = data.table(
            region = c("Exon side", "Intron side", "Full immediate"),
            start = c(-50L, 0L, -50L), end = c(0L, 50L, 50L)
        )
    ),
    splice_acceptor_3prime = list(
        slug = "acceptor", figure_dir = "acceptors",
        label = "3' splice acceptor", left = "INTRON", right = "EXON",
        expected_bins = 50L,
        regions = data.table(
            region = c("Intron side", "Exon side", "Full immediate"),
            start = c(-50L, 0L, -50L), end = c(0L, 50L, 50L)
        )
    ),
    annotated_transcript_3prime_end = list(
        slug = "transcript_3prime_end", figure_dir = "transcript_3prime_ends",
        label = "annotated transcript 3' end", left = "TRANSCRIPT",
        right = "DOWNSTREAM", expected_bins = 10L,
        regions = data.table(
            region = c("Upstream", "Downstream", "Full region"),
            start = c(-50L, 0L, -50L), end = c(0L, 50L, 50L)
        )
    )
)

metadata <- fread(file.path(input_dir, "sample_metadata.tsv"))
metadata[, construct := as.character(construct)]
if (nrow(metadata) != 8L || anyDuplicated(metadata$sample_id)) {
    stop("Expected eight unique samples")
}
if (!setequal(metadata$construct, c("324", "328")) ||
        !setequal(metadata$assay, c("CLIP", "SMI"))) {
    stop("Unexpected construct or assay labels")
}
metadata[, library_millions := total_primary_mapped_reads / 1e6]
metadata[, common_cpm_factor := library_millions]
normalization_qc <- metadata[, .(
    sample = sample_id, construct, assay,
    replicate = biological_replicate, pair_id,
    primary_mapped_reads = total_primary_mapped_reads,
    uniquely_mapped_reads = star_uniquely_mapped_reads,
    uniquely_mapped_percent = star_uniquely_mapped_percent,
    normalization_factor = common_cpm_factor,
    normalization_method = "CPM per primary mapped read",
    normalization_formula = "raw_count / primary_mapped_reads * 1e6"
)]
setorder(normalization_qc, construct, replicate, assay)
write_tsv(
    normalization_qc,
    file.path(results_dir, "normalization", "normalization_qc.tsv")
)

pair_map <- metadata[, {
    if (sum(assay == "CLIP") != 1L || sum(assay == "SMI") != 1L) {
        stop("Every pair must contain one CLIP and one SMI library")
    }
    .(
        clip_sample = sample_id[assay == "CLIP"],
        smi_sample = sample_id[assay == "SMI"]
    )
}, by = .(construct, biological_replicate, pair_id)]
if (nrow(pair_map) != 4L) stop("Expected four matched CLIP/SMI pairs")

sites_all <- fread(file.path(input_dir, "landmark_sites.tsv.gz"))
sites <- sites_all[retained_for_bin_counts == "yes"]
if (nrow(sites) != 237732L || anyDuplicated(sites$site_id)) {
    stop("Unexpected retained landmark universe")
}
sites[, primary_gene_id := sub(";.*$", "", gene_id)]
sites[, multiple_gene_ids := grepl(";", gene_id, fixed = TRUE)]
if (any(!nzchar(sites$primary_gene_id))) {
    stop("Every retained site needs a gene ID")
}

bin_counts <- fread(file.path(input_dir, "landmark_bin_counts.tsv.gz"))
if (any(bin_counts$raw_overlap_count <= 0L) ||
        anyDuplicated(bin_counts[, .(sample_id, site_id, relative_bin_start)])) {
    stop("Sparse bin table must contain unique positive rows")
}
library_lookup <- setNames(metadata$library_millions, metadata$sample_id)
bin_counts[, normalized_count_cpm := raw_overlap_count / library_lookup[sample_id]]

# Condition-independent coverage uses one window-overlap count per sample/site.
message("Reading condition-independent window coverage")
site_summary <- fread(
    file.path(input_dir, "landmark_site_summary.tsv.gz"),
    select = c("sample_id", "site_id", "site_type", "raw_overlap_count")
)
site_summary[, normalized_window_coverage_cpm :=
    raw_overlap_count / library_lookup[sample_id]]
coverage <- site_summary[, .(
    total_normalized_coverage_all_8_cpm = sum(normalized_window_coverage_cpm),
    total_raw_window_overlap_all_8 = sum(raw_overlap_count)
), by = .(site_id, site_type)]
rm(site_summary)
invisible(gc())
coverage <- coverage[sites[, .(
    site_id, primary_gene_id, multiple_gene_ids
)], on = "site_id", nomatch = 0]
if (nrow(coverage) != nrow(sites) ||
        any(coverage$total_normalized_coverage_all_8_cpm <= 0)) {
    stop("Coverage table does not match the retained-site universe")
}
quartile_labels <- c("Q1 lowest", "Q2", "Q3", "Q4 highest")
setorder(coverage, site_type, total_normalized_coverage_all_8_cpm, site_id)
coverage[, coverage_rank_fraction := seq_len(.N) / .N, by = site_type]
coverage[, coverage_quartile := factor(
    quartile_labels[pmin(4L, pmax(1L, ceiling(coverage_rank_fraction * 4)))],
    levels = quartile_labels
)]
coverage_thresholds <- coverage[, .(
    site_count = .N,
    gene_count = uniqueN(primary_gene_id),
    minimum_total_cpm = min(total_normalized_coverage_all_8_cpm),
    median_total_cpm = median(total_normalized_coverage_all_8_cpm),
    maximum_total_cpm = max(total_normalized_coverage_all_8_cpm),
    minimum_total_raw_window_overlap = min(total_raw_window_overlap_all_8),
    maximum_total_raw_window_overlap = max(total_raw_window_overlap_all_8)
), by = .(site_type, coverage_quartile)]
coverage_concentration <- coverage[, {
    coverage_values <- total_normalized_coverage_all_8_cpm
    top_n <- max(1L, ceiling(.N * 0.01))
    descending <- order(coverage_values, decreasing = TRUE)
    .(
        N_sites = .N,
        Q4_fraction_total_coverage = sum(
            coverage_values[coverage_quartile == "Q4 highest"]
        ) / sum(coverage_values),
        top_1pct_site_count = top_n,
        top_1pct_fraction_total_coverage = sum(
            coverage_values[descending[seq_len(top_n)]]
        ) / sum(coverage_values),
        maximum_to_median_coverage_ratio =
            max(coverage_values) / median(coverage_values)
    )
}, by = site_type]
write_tsv_gz(
    coverage,
    file.path(
        results_dir, "coverage_dependence",
        "condition_independent_site_coverage.tsv.gz"
    )
)
write_tsv(
    coverage_thresholds,
    file.path(results_dir, "coverage_dependence", "coverage_quartile_definitions.tsv")
)
write_tsv(
    coverage_concentration,
    file.path(results_dir, "coverage_dependence", "coverage_concentration_summary.tsv")
)

pair_enrichment <- function(counts, pseudocount) {
    out <- vector("list", nrow(pair_map))
    for (i in seq_len(nrow(pair_map))) {
        pm <- pair_map[i]
        clip <- counts[sample_id == pm$clip_sample, .(
            site_id, relative_bin_start, relative_bin_end,
            normalized_clip_cpm = normalized_count_cpm
        )]
        smi <- counts[sample_id == pm$smi_sample, .(
            site_id, relative_bin_start, relative_bin_end,
            normalized_smi_cpm = normalized_count_cpm
        )]
        joined <- merge(
            clip, smi,
            by = c("site_id", "relative_bin_start", "relative_bin_end"),
            all = TRUE
        )
        joined[is.na(normalized_clip_cpm), normalized_clip_cpm := 0]
        joined[is.na(normalized_smi_cpm), normalized_smi_cpm := 0]
        joined[, enrichment := log2(
            (normalized_clip_cpm + pseudocount) /
                (normalized_smi_cpm + pseudocount)
        )]
        joined[, `:=`(
            construct = pm$construct,
            replicate = pm$biological_replicate,
            pair_id = pm$pair_id
        )]
        out[[i]] <- joined
    }
    rbindlist(out)
}

profile_components <- function(pair_values, site_count, bins) {
    replicate_profile <- pair_values[, .(
        enrichment = sum(enrichment) / site_count
    ), by = .(construct, replicate, pair_id, relative_bin_start, relative_bin_end)]
    rep_grid <- CJ(
        construct = c("324", "328"), replicate = c("A", "B"),
        relative_bin_start = bins, unique = TRUE
    )
    rep_grid[, relative_bin_end := relative_bin_start + 10L]
    rep_grid <- merge(
        rep_grid,
        pair_map[, .(construct, replicate = biological_replicate, pair_id)],
        by = c("construct", "replicate"), all.x = TRUE
    )
    replicate_profile <- merge(
        rep_grid, replicate_profile,
        by = c(
            "construct", "replicate", "pair_id",
            "relative_bin_start", "relative_bin_end"
        ), all.x = TRUE
    )
    replicate_profile[is.na(enrichment), enrichment := 0]
    replicate_profile[, relative_position := relative_bin_start + 5]

    construct_sparse <- pair_values[, .(
        enrichment = sum(enrichment) / 2
    ), by = .(construct, site_id, relative_bin_start, relative_bin_end)]
    list(replicate = replicate_profile, construct_sparse = construct_sparse)
}

bootstrap_construct_profile <- function(
    construct_sparse, site_map, bins, seed_offset = 0L
) {
    gene_sizes <- site_map[, .(site_count = .N), by = primary_gene_id]
    genes <- gene_sizes$primary_gene_id
    mats <- list()
    for (construct_value in c("324", "328")) {
        gene_bin <- construct_sparse[construct == construct_value][
            site_map[, .(site_id, primary_gene_id)],
            on = "site_id", nomatch = 0
        ][, .(value = sum(enrichment)),
          by = .(primary_gene_id, relative_bin_start)]
        mat <- matrix(
            0, nrow = length(genes), ncol = length(bins),
            dimnames = list(genes, as.character(bins))
        )
        if (nrow(gene_bin)) {
            mat[cbind(
                match(gene_bin$primary_gene_id, genes),
                match(gene_bin$relative_bin_start, bins)
            )] <- gene_bin$value
        }
        mats[[construct_value]] <- mat
    }
    central_324 <- colSums(mats[["324"]]) / sum(gene_sizes$site_count)
    central_328 <- colSums(mats[["328"]]) / sum(gene_sizes$site_count)
    boot_324 <- matrix(NA_real_, bootstrap_iterations, length(bins))
    boot_328 <- matrix(NA_real_, bootstrap_iterations, length(bins))
    set.seed(bootstrap_seed + seed_offset)
    starts <- seq.int(1L, bootstrap_iterations, by = 100L)
    for (start in starts) {
        finish <- min(start + 99L, bootstrap_iterations)
        weights <- rmultinom(
            finish - start + 1L, size = length(genes),
            prob = rep(1 / length(genes), length(genes))
        )
        denominator <- as.numeric(crossprod(gene_sizes$site_count, weights))
        boot_324[start:finish, ] <- sweep(
            crossprod(weights, mats[["324"]]), 1L, denominator, "/"
        )
        boot_328[start:finish, ] <- sweep(
            crossprod(weights, mats[["328"]]), 1L, denominator, "/"
        )
    }
    boot_delta <- boot_324 - boot_328
    interval <- function(x) {
        t(apply(x, 2L, quantile, c(0.025, 0.975), names = FALSE))
    }
    ci324 <- interval(boot_324)
    ci328 <- interval(boot_328)
    cidd <- interval(boot_delta)
    data.table(
        relative_bin_start = bins,
        relative_bin_end = bins + 10L,
        relative_position = bins + 5,
        enrichment_324 = central_324,
        enrichment_324_ci_low = ci324[, 1L],
        enrichment_324_ci_high = ci324[, 2L],
        enrichment_328 = central_328,
        enrichment_328_ci_low = ci328[, 1L],
        enrichment_328_ci_high = ci328[, 2L],
        delta_delta = central_324 - central_328,
        delta_delta_ci_low = cidd[, 1L],
        delta_delta_ci_high = cidd[, 2L]
    )
}

delta_sparse_from_construct <- function(construct_sparse) {
    wide <- dcast(
        construct_sparse,
        site_id + relative_bin_start + relative_bin_end ~ construct,
        value.var = "enrichment", fill = 0
    )
    if (!all(c("324", "328") %in% names(wide))) {
        stop("Both constructs are required for delta-delta")
    }
    wide[, delta_delta := `324` - `328`]
    wide
}

site_region_values <- function(value_sparse, site_map, regions) {
    out <- vector("list", nrow(regions))
    for (i in seq_len(nrow(regions))) {
        region_row <- regions[i]
        n_bins <- (region_row$end - region_row$start) / 10
        sparse <- value_sparse[
            relative_bin_start >= region_row$start &
                relative_bin_end <= region_row$end,
            .(delta_delta = sum(delta_delta) / n_bins),
            by = site_id
        ]
        full <- merge(
            site_map[, .(
                site_id, gene_id, primary_gene_id, multiple_gene_ids
            )],
            sparse, by = "site_id", all.x = TRUE
        )
        full[is.na(delta_delta), delta_delta := 0]
        full[, `:=`(
            region = region_row$region,
            region_start = region_row$start,
            region_end = region_row$end
        )]
        out[[i]] <- full
    }
    rbindlist(out)
}

replicate_region_values <- function(pair_values, site_count, regions) {
    out <- vector("list", nrow(regions))
    for (i in seq_len(nrow(regions))) {
        region_row <- regions[i]
        n_bins <- (region_row$end - region_row$start) / 10
        values <- pair_values[
            relative_bin_start >= region_row$start &
                relative_bin_end <= region_row$end,
            .(regional_enrichment = sum(enrichment) / (site_count * n_bins)),
            by = .(construct, replicate, pair_id)
        ]
        values[, `:=`(
            region = region_row$region,
            region_start = region_row$start,
            region_end = region_row$end
        )]
        out[[i]] <- values
    }
    rbindlist(out)
}

bootstrap_cluster_scalar <- function(
    values, gene_col = "primary_gene_id", value_col = "delta_delta",
    seed_offset = 0L
) {
    gene_stats <- values[, .(
        value_sum = sum(get(value_col)), site_count = .N
    ), by = c(gene_col)]
    set.seed(bootstrap_seed + seed_offset)
    weights <- rmultinom(
        bootstrap_iterations, size = nrow(gene_stats),
        prob = rep(1 / nrow(gene_stats), nrow(gene_stats))
    )
    means <- as.numeric(crossprod(gene_stats$value_sum, weights)) /
        as.numeric(crossprod(gene_stats$site_count, weights))
    quantile(means, c(0.025, 0.975), names = FALSE)
}

paired_statistics <- function(x, y, unit_label, seed_offset = 0L) {
    difference <- x - y
    nonzero <- difference[is.finite(difference) & difference != 0]
    ranks <- rank(abs(nonzero), ties.method = "average")
    rank_biserial <- if (length(nonzero)) {
        sum(sign(nonzero) * ranks) / sum(ranks)
    } else {
        0
    }
    wilcox <- suppressWarnings(wilcox.test(
        x, y, paired = TRUE, exact = FALSE, correct = FALSE
    ))
    data.table(
        statistical_unit = unit_label,
        N = length(difference),
        mean_side_a = mean(x), median_side_a = median(x),
        mean_side_b = mean(y), median_side_b = median(y),
        mean_paired_difference = mean(difference),
        median_paired_difference = median(difference),
        paired_standardized_mean_difference = mean(difference) / sd(difference),
        paired_rank_biserial = rank_biserial,
        paired_wilcoxon_V = unname(wilcox$statistic),
        paired_wilcoxon_p = wilcox$p.value,
        nonzero_pair_count = length(nonzero)
    )
}

bootstrap_stratified_delta <- function(
    delta_sparse, site_map, bins, seed_offset = 0L
) {
    out <- vector("list", length(quartile_labels))
    for (i in seq_along(quartile_labels)) {
        stratum <- quartile_labels[[i]]
        map <- site_map[coverage_quartile == stratum]
        gene_sizes <- map[, .(site_count = .N), by = primary_gene_id]
        genes <- gene_sizes$primary_gene_id
        gene_bin <- delta_sparse[
            map[, .(site_id, primary_gene_id)], on = "site_id", nomatch = 0
        ][, .(value = sum(delta_delta)),
          by = .(primary_gene_id, relative_bin_start)]
        mat <- matrix(
            0, nrow = length(genes), ncol = length(bins),
            dimnames = list(genes, as.character(bins))
        )
        if (nrow(gene_bin)) {
            mat[cbind(
                match(gene_bin$primary_gene_id, genes),
                match(gene_bin$relative_bin_start, bins)
            )] <- gene_bin$value
        }
        central <- colSums(mat) / sum(gene_sizes$site_count)
        boot <- matrix(NA_real_, bootstrap_iterations, length(bins))
        set.seed(bootstrap_seed + seed_offset + i)
        starts <- seq.int(1L, bootstrap_iterations, by = 100L)
        for (start in starts) {
            finish <- min(start + 99L, bootstrap_iterations)
            weights <- rmultinom(
                finish - start + 1L, size = length(genes),
                prob = rep(1 / length(genes), length(genes))
            )
            denominator <- as.numeric(crossprod(gene_sizes$site_count, weights))
            boot[start:finish, ] <- sweep(
                crossprod(weights, mat), 1L, denominator, "/"
            )
        }
        ci <- t(apply(
            boot, 2L, quantile, c(0.025, 0.975), names = FALSE
        ))
        out[[i]] <- data.table(
            coverage_quartile = factor(stratum, levels = quartile_labels),
            relative_bin_start = bins,
            relative_bin_end = bins + 10L,
            relative_position = bins + 5,
            delta_delta = central,
            ci_low = ci[, 1L], ci_high = ci[, 2L],
            site_count = nrow(map), gene_count = length(genes)
        )
    }
    rbindlist(out)
}

aggregate_delta_profile <- function(counts, pseudocount) {
    aggregate <- counts[, .(
        sum_normalized_count_cpm = sum(normalized_count_cpm)
    ), by = .(sample_id, relative_bin_start, relative_bin_end)]
    aggregate <- merge(
        aggregate,
        metadata[, .(
            sample_id, construct, assay,
            replicate = biological_replicate, pair_id
        )],
        by = "sample_id"
    )
    pair <- dcast(
        aggregate,
        construct + replicate + pair_id + relative_bin_start +
            relative_bin_end ~ assay,
        value.var = "sum_normalized_count_cpm", fill = 0
    )
    pair[, enrichment := log2(
        (CLIP + pseudocount) / (SMI + pseudocount)
    )]
    construct <- pair[, .(enrichment = mean(enrichment)),
                      by = .(construct, relative_bin_start, relative_bin_end)]
    wide <- dcast(
        construct,
        relative_bin_start + relative_bin_end ~ construct,
        value.var = "enrichment", fill = 0
    )
    wide[, `:=`(
        relative_position = relative_bin_start + 5,
        delta_delta = `324` - `328`
    )]
    wide
}

aggregate_delta_profile_by_coverage <- function(
    counts, site_map, pseudocount
) {
    aggregate <- counts[
        site_map[, .(site_id, coverage_quartile)],
        on = "site_id", nomatch = 0
    ][, .(
        sum_normalized_count_cpm = sum(normalized_count_cpm)
    ), by = .(
        coverage_quartile, sample_id,
        relative_bin_start, relative_bin_end
    )]
    aggregate <- merge(
        aggregate,
        metadata[, .(
            sample_id, construct, assay,
            replicate = biological_replicate, pair_id
        )],
        by = "sample_id"
    )
    pair <- dcast(
        aggregate,
        coverage_quartile + construct + replicate + pair_id +
            relative_bin_start + relative_bin_end ~ assay,
        value.var = "sum_normalized_count_cpm", fill = 0
    )
    pair[, enrichment := log2(
        (CLIP + pseudocount) / (SMI + pseudocount)
    )]
    construct <- pair[, .(enrichment = mean(enrichment)), by = .(
        coverage_quartile, construct,
        relative_bin_start, relative_bin_end
    )]
    wide <- dcast(
        construct,
        coverage_quartile + relative_bin_start + relative_bin_end ~ construct,
        value.var = "enrichment", fill = 0
    )
    wide[, `:=`(
        relative_position = relative_bin_start + 5,
        delta_delta = `324` - `328`
    )]
    wide
}

add_side_labels <- function(plot, config, x_min, x_max) {
    plot + annotate(
        "text", x = c(x_min / 2, x_max / 2), y = Inf,
        label = c(config$left, config$right), vjust = 1.4,
        fontface = "bold", size = 3.1, color = "grey30"
    )
}

normalization_comparison <- list()
sensitivity_summary <- list()
replicate_region_all <- list()
side_summary_all <- list()
coverage_profile_all <- list()
coverage_region_all <- list()
coverage_read_profile_all <- list()
coverage_weighting_region_all <- list()
pattern_summary <- list()

for (type_index in seq_along(type_config)) {
    type_name <- names(type_config)[[type_index]]
    config <- type_config[[type_name]]
    message("Analyzing ", config$label)
    site_map <- sites[site_type == type_name, .(
        site_id, gene_id, primary_gene_id, multiple_gene_ids
    )]
    site_map <- coverage[site_type == type_name, .(
        site_id, coverage_quartile,
        total_normalized_coverage_all_8_cpm
    )][site_map, on = "site_id"]
    site_count <- nrow(site_map)
    counts <- bin_counts[site_type == type_name]
    bins <- sort(unique(counts$relative_bin_start))
    if (length(bins) != config$expected_bins) {
        stop("Unexpected bin count for ", type_name)
    }

    primary_pair <- pair_enrichment(counts, primary_pseudocount)
    primary_components <- profile_components(primary_pair, site_count, bins)
    common_profile <- bootstrap_construct_profile(
        primary_components$construct_sparse, site_map, bins,
        seed_offset = type_index * 100L
    )
    common_profile[, `:=`(
        site_type = type_name,
        retained_site_count = site_count,
        primary_gene_count = uniqueN(site_map$primary_gene_id),
        normalization = "common primary-mapped CPM",
        pseudocount_cpm = primary_pseudocount,
        bootstrap_iterations = bootstrap_iterations
    )]
    write_tsv(
        common_profile,
        file.path(
            results_dir, "metaprofiles",
            paste0(config$slug, "_common_cpm_site_weighted_metaprofile.tsv")
        )
    )
    write_tsv(
        primary_components$replicate,
        file.path(
            results_dir, "metaprofiles",
            paste0(config$slug, "_common_cpm_replicate_profiles.tsv")
        )
    )

    sensitivity_pair <- pair_enrichment(counts, sensitivity_pseudocount)
    sensitivity_components <- profile_components(
        sensitivity_pair, site_count, bins
    )
    sensitivity_profile <- bootstrap_construct_profile(
        sensitivity_components$construct_sparse, site_map, bins,
        seed_offset = type_index * 100L + 50L
    )
    write_tsv(
        sensitivity_profile,
        file.path(
            results_dir, "normalization",
            paste0(config$slug, "_common_cpm_p0.5_site_weighted_metaprofile.tsv")
        )
    )
    sensitivity_compare <- merge(
        common_profile[, .(
            relative_position,
            primary_delta_delta = delta_delta
        )],
        sensitivity_profile[, .(
            relative_position,
            sensitivity_delta_delta = delta_delta
        )],
        by = "relative_position"
    )
    sensitivity_summary[[type_name]] <- sensitivity_compare[, .(
        site_type = type_name,
        primary_pseudocount_cpm = primary_pseudocount,
        sensitivity_pseudocount_cpm = sensitivity_pseudocount,
        sign_change_bins = sum(
            sign(primary_delta_delta) != sign(sensitivity_delta_delta)
        ),
        pearson_r = cor(primary_delta_delta, sensitivity_delta_delta),
        spearman_rho = cor(
            primary_delta_delta, sensitivity_delta_delta,
            method = "spearman"
        ),
        primary_mean_delta_delta = mean(primary_delta_delta),
        sensitivity_mean_delta_delta = mean(sensitivity_delta_delta)
    )]

    old_path <- file.path(
        old_output_dir, "results", "metaprofiles",
        paste0(config$slug, "_site_weighted_metaprofile.tsv")
    )
    old_profile <- fread(old_path)
    comparison <- merge(
        old_profile[, .(
            relative_position,
            old_delta_delta = delta_delta,
            old_enrichment_324 = enrichment_324,
            old_enrichment_328 = enrichment_328
        )],
        common_profile[, .(
            relative_position,
            common_cpm_delta_delta = delta_delta,
            common_cpm_enrichment_324 = enrichment_324,
            common_cpm_enrichment_328 = enrichment_328
        )],
        by = "relative_position"
    )
    comparison[, `:=`(
        site_type = type_name,
        delta_sign_changed = sign(old_delta_delta) != sign(common_cpm_delta_delta),
        delta_absolute_change = abs(common_cpm_delta_delta - old_delta_delta)
    )]
    write_tsv(
        comparison,
        file.path(
            results_dir, "normalization",
            paste0(config$slug, "_old_vs_common_cpm_profiles.tsv")
        )
    )
    target_position <- if (type_name == "splice_donor_5prime") {
        -5
    } else if (type_name == "splice_acceptor_3prime") {
        15
    } else {
        NA_real_
    }
    normalization_comparison[[type_name]] <- comparison[, .(
        site_type = type_name,
        bin_count = .N,
        delta_delta_sign_change_bins = sum(delta_sign_changed),
        enrichment_324_sign_change_bins = sum(
            sign(old_enrichment_324) != sign(common_cpm_enrichment_324)
        ),
        enrichment_328_sign_change_bins = sum(
            sign(old_enrichment_328) != sign(common_cpm_enrichment_328)
        ),
        delta_profile_pearson_r = cor(old_delta_delta, common_cpm_delta_delta),
        delta_profile_spearman_rho = cor(
            old_delta_delta, common_cpm_delta_delta, method = "spearman"
        ),
        previous_mean_delta_delta = mean(old_delta_delta),
        common_cpm_mean_delta_delta = mean(common_cpm_delta_delta),
        mean_delta_delta_change = mean(common_cpm_delta_delta - old_delta_delta),
        previous_strongest_negative_position =
            relative_position[which.min(old_delta_delta)],
        common_cpm_strongest_negative_position =
            relative_position[which.min(common_cpm_delta_delta)],
        requested_extremum_position = target_position,
        requested_extremum_remains = if (is.na(target_position)) {
            NA
        } else {
            relative_position[which.min(common_cpm_delta_delta)] == target_position
        }
    )]

    delta_sparse <- delta_sparse_from_construct(
        primary_components$construct_sparse
    )
    region_sites <- site_region_values(delta_sparse, site_map, config$regions)
    replicate_regions <- replicate_region_values(
        primary_pair, site_count, config$regions
    )
    replicate_regions[, `:=`(
        site_type = type_name,
        landmark = config$label,
        sample_label = paste(construct, "replicate", replicate)
    )]
    replicate_region_all[[type_name]] <- replicate_regions
    write_tsv(
        replicate_regions,
        file.path(
            results_dir, "replicate_regions",
            paste0(config$slug, "_replicate_region_estimates.tsv")
        )
    )

    if (type_name %in% c(
        "splice_donor_5prime", "splice_acceptor_3prime"
    )) {
        side_a_name <- if (type_name == "splice_donor_5prime") {
            "Exon side"
        } else {
            "Exon side"
        }
        side_b_name <- if (type_name == "splice_donor_5prime") {
            "Intron side"
        } else {
            "Intron side"
        }
        side_a <- region_sites[region == side_a_name, .(
            site_id, gene_id, primary_gene_id, multiple_gene_ids,
            side_a_delta_delta = delta_delta
        )]
        side_b <- region_sites[region == side_b_name, .(
            site_id, side_b_delta_delta = delta_delta
        )]
        side_values <- merge(side_a, side_b, by = "site_id")
        side_values[, side_difference :=
            side_a_delta_delta - side_b_delta_delta]
        side_values[, `:=`(
            site_type = type_name,
            side_a = side_a_name,
            side_b = side_b_name,
            side_difference_definition = "exonic delta-delta - intronic delta-delta"
        )]
        write_tsv_gz(
            side_values,
            file.path(
                results_dir, "side_of_splice",
                paste0(config$slug, "_site_side_values.tsv.gz")
            )
        )
        gene_values <- side_values[, .(
            side_a_delta_delta = median(side_a_delta_delta),
            side_b_delta_delta = median(side_b_delta_delta),
            site_count = .N
        ), by = primary_gene_id]
        gene_values[, side_difference :=
            side_a_delta_delta - side_b_delta_delta]
        gene_values[, `:=`(
            site_type = type_name,
            side_a = side_a_name,
            side_b = side_b_name
        )]
        write_tsv(
            gene_values,
            file.path(
                results_dir, "side_of_splice",
                paste0(config$slug, "_gene_side_values.tsv")
            )
        )

        site_stats <- paired_statistics(
            side_values$side_a_delta_delta,
            side_values$side_b_delta_delta,
            "site"
        )
        site_ci <- bootstrap_cluster_scalar(
            side_values, value_col = "side_difference",
            seed_offset = type_index * 1000L + 1L
        )
        site_stats[, `:=`(
            N_sites = nrow(side_values),
            N_genes = uniqueN(side_values$primary_gene_id),
            mean_difference_bootstrap_ci_low = site_ci[1L],
            mean_difference_bootstrap_ci_high = site_ci[2L]
        )]
        gene_stats <- paired_statistics(
            gene_values$side_a_delta_delta,
            gene_values$side_b_delta_delta,
            "gene median"
        )
        set.seed(bootstrap_seed + type_index * 1000L + 2L)
        gene_boot <- replicate(
            bootstrap_iterations,
            mean(sample(gene_values$side_difference, replace = TRUE))
        )
        gene_stats[, `:=`(
            N_sites = nrow(side_values),
            N_genes = nrow(gene_values),
            mean_difference_bootstrap_ci_low = quantile(
                gene_boot, 0.025, names = FALSE
            ),
            mean_difference_bootstrap_ci_high = quantile(
                gene_boot, 0.975, names = FALSE
            )
        )]
        side_stats <- rbindlist(list(site_stats, gene_stats), fill = TRUE)
        side_stats[, `:=`(
            site_type = type_name,
            side_a = side_a_name,
            side_b = side_b_name,
            paired_difference = "exonic delta-delta - intronic delta-delta",
            negative_interpretation =
                "stronger 328-over-324 enrichment on the exonic side"
        )]
        side_summary_all[[type_name]] <- side_stats
        write_tsv(
            side_stats,
            file.path(
                results_dir, "side_of_splice",
                paste0(config$slug, "_side_summary.tsv")
            )
        )

        gene_long <- melt(
            gene_values,
            id.vars = c("primary_gene_id", "site_count"),
            measure.vars = c("side_a_delta_delta", "side_b_delta_delta"),
            variable.name = "side", value.name = "delta_delta"
        )
        gene_long[, side := factor(
            side,
            levels = c("side_a_delta_delta", "side_b_delta_delta"),
            labels = c(side_a_name, side_b_name)
        )]
        p_side_values <- ggplot(
            gene_long, aes(side, delta_delta, fill = side)
        ) +
            geom_violin(scale = "width", trim = TRUE, alpha = 0.65) +
            geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white") +
            geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
            scale_fill_manual(values = c("#4DAF4A", "#377EB8")) +
            labs(
                title = "Gene-level side values",
                x = NULL, y = "Gene median delta-delta"
            ) +
            theme_publication() +
            theme(legend.position = "none")
        p_side_difference <- ggplot(
            gene_values, aes(side_difference)
        ) +
            geom_histogram(bins = 60, fill = "#6A3D9A", color = "white") +
            geom_vline(xintercept = 0, color = "grey30", linetype = 2) +
            labs(
                title = "Gene-level paired side difference",
                x = "Exonic delta-delta - intronic delta-delta",
                y = "Genes"
            ) +
            theme_publication()
        p_side <- p_side_values + p_side_difference +
            plot_annotation(
                title = paste0(config$label, ": exonic versus intronic effect"),
                subtitle = paste0(
                    "Negative difference indicates stronger 328-over-324 enrichment on the exonic side; gene-level Wilcoxon p = ",
                    format(gene_stats$paired_wilcoxon_p, digits = 3, scientific = TRUE)
                )
            )
        save_plot(
            p_side,
            file.path(
                figures_dir, "side_of_splice",
                paste0(config$slug, "_exonic_vs_intronic")
            ),
            width = 9.5, height = 4.7
        )
    }

    stratified_profile <- bootstrap_stratified_delta(
        delta_sparse, site_map, bins,
        seed_offset = type_index * 10000L
    )
    stratified_profile[, site_type := type_name]
    coverage_profile_all[[type_name]] <- stratified_profile
    write_tsv(
        stratified_profile,
        file.path(
            results_dir, "coverage_dependence",
            paste0(config$slug, "_coverage_stratified_metaprofile.tsv")
        )
    )
    coverage_regions <- region_sites[
        site_map[, .(site_id, coverage_quartile)], on = "site_id"
    ][, {
        ci <- bootstrap_cluster_scalar(
            .SD, value_col = "delta_delta",
            seed_offset = type_index * 10000L + .GRP
        )
        .(
            site_count = .N,
            gene_count = uniqueN(primary_gene_id),
            mean_delta_delta = mean(delta_delta),
            median_delta_delta = median(delta_delta),
            mean_ci_low = ci[1L], mean_ci_high = ci[2L],
            proportion_negative = mean(delta_delta < 0),
            proportion_positive = mean(delta_delta > 0)
        )
    }, by = .(coverage_quartile, region, region_start, region_end)]
    coverage_regions[, site_type := type_name]
    coverage_region_all[[type_name]] <- coverage_regions
    write_tsv(
        coverage_regions,
        file.path(
            results_dir, "coverage_dependence",
            paste0(config$slug, "_coverage_stratified_regions.tsv")
        )
    )

    aggregate_profile <- aggregate_delta_profile(
        counts, primary_pseudocount
    )
    aggregate_profile[, site_type := type_name]
    write_tsv(
        aggregate_profile,
        file.path(
            results_dir, "metaprofiles",
            paste0(config$slug, "_common_cpm_read_weighted_metaprofile.tsv")
        )
    )

    coverage_read_profile <- aggregate_delta_profile_by_coverage(
        counts, site_map, primary_pseudocount
    )
    coverage_read_profile[, site_type := type_name]
    coverage_read_profile_all[[type_name]] <- coverage_read_profile
    write_tsv(
        coverage_read_profile,
        file.path(
            results_dir, "coverage_dependence",
            paste0(
                config$slug,
                "_coverage_stratified_read_weighted_metaprofile.tsv"
            )
        )
    )
    coverage_read_regions <- rbindlist(lapply(
        seq_len(nrow(config$regions)),
        function(i) {
            region_row <- config$regions[i]
            coverage_read_profile[
                relative_bin_start >= region_row$start &
                    relative_bin_end <= region_row$end,
                .(
                    mean_delta_delta = mean(delta_delta),
                    minimum_delta_delta = min(delta_delta),
                    maximum_delta_delta = max(delta_delta)
                ),
                by = coverage_quartile
            ][, `:=`(
                region = region_row$region,
                region_start = region_row$start,
                region_end = region_row$end
            )]
        }
    ))
    coverage_weighting_regions <- rbindlist(list(
        coverage_regions[, .(
            coverage_quartile, region, region_start, region_end,
            weighting = "Site-weighted",
            mean_delta_delta,
            mean_ci_low, mean_ci_high
        )],
        coverage_read_regions[, .(
            coverage_quartile, region, region_start, region_end,
            weighting = "Read-weighted aggregate",
            mean_delta_delta,
            mean_ci_low = NA_real_, mean_ci_high = NA_real_
        )]
    ))
    coverage_weighting_regions[, site_type := type_name]
    coverage_weighting_region_all[[type_name]] <- coverage_weighting_regions
    write_tsv(
        coverage_weighting_regions,
        file.path(
            results_dir, "coverage_dependence",
            paste0(config$slug, "_coverage_stratified_weighting_comparison.tsv")
        )
    )

    pattern_summary[[type_name]] <- data.table(
        site_type = type_name,
        site_weighted_mean_delta_delta = mean(common_profile$delta_delta),
        site_weighted_minimum_delta_delta = min(common_profile$delta_delta),
        site_weighted_minimum_position = common_profile[
            which.min(delta_delta), relative_position
        ],
        read_weighted_mean_delta_delta = mean(aggregate_profile$delta_delta),
        read_weighted_minimum_delta_delta = min(aggregate_profile$delta_delta),
        read_weighted_maximum_delta_delta = max(aggregate_profile$delta_delta)
    )

    x_min <- min(bins)
    x_max <- max(bins) + 10
    figure_base <- file.path(figures_dir, config$figure_dir)
    enrichment_long <- rbindlist(list(
        common_profile[, .(
            construct = "324", relative_position,
            enrichment = enrichment_324,
            ci_low = enrichment_324_ci_low,
            ci_high = enrichment_324_ci_high
        )],
        common_profile[, .(
            construct = "328", relative_position,
            enrichment = enrichment_328,
            ci_low = enrichment_328_ci_low,
            ci_high = enrichment_328_ci_high
        )]
    ))
    p_enrichment <- ggplot(
        enrichment_long,
        aes(relative_position, enrichment, color = construct, fill = construct)
    ) +
        geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
        geom_vline(xintercept = 0, color = "grey40", linetype = 2) +
        geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.16, color = NA) +
        geom_line(
            data = primary_components$replicate,
            aes(
                relative_position, enrichment,
                group = interaction(construct, replicate)
            ),
            alpha = 0.25, linewidth = 0.5
        ) +
        geom_line(linewidth = 1.05) +
        scale_color_manual(values = construct_colors) +
        scale_fill_manual(values = construct_colors) +
        scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
        labs(
            title = paste0(config$label, ": common-CPM site-weighted enrichment"),
            subtitle = paste0(
                "Mean site-level log2(CLIP/SMI); 95% gene-bootstrap CI; p = ",
                primary_pseudocount, " CPM"
            ),
            x = "Relative position (nt)", y = "Mean log2(CLIP / SMI)",
            color = "Construct", fill = "Construct",
            caption = paste0("Boundary: ", config$left, " | ", config$right)
        ) +
        theme_publication()
    p_enrichment <- add_side_labels(p_enrichment, config, x_min, x_max)
    save_plot(
        p_enrichment,
        file.path(figure_base, "common_cpm_site_weighted_enrichment")
    )

    p_delta <- ggplot(common_profile, aes(relative_position, delta_delta)) +
        geom_hline(yintercept = 0, color = "grey45", linewidth = 0.4) +
        geom_vline(xintercept = 0, color = "grey40", linetype = 2) +
        geom_ribbon(
            aes(ymin = delta_delta_ci_low, ymax = delta_delta_ci_high),
            fill = "#6A3D9A", alpha = 0.18
        ) +
        geom_line(color = "#6A3D9A", linewidth = 1.05) +
        scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
        labs(
            title = paste0(config$label, ": common-CPM site-weighted delta-delta"),
            subtitle = "Positive favors 324; negative favors 328; 95% gene-bootstrap CI",
            x = "Relative position (nt)",
            y = "324 enrichment - 328 enrichment",
            caption = paste0("Boundary: ", config$left, " | ", config$right)
        ) +
        theme_publication()
    p_delta <- add_side_labels(p_delta, config, x_min, x_max)
    save_plot(
        p_delta,
        file.path(figure_base, "common_cpm_site_weighted_delta_delta")
    )

    old_new_plot <- rbindlist(list(
        comparison[, .(
            relative_position, delta_delta = old_delta_delta,
            normalization = "Previous construct-specific scale"
        )],
        comparison[, .(
            relative_position, delta_delta = common_cpm_delta_delta,
            normalization = "Common primary-mapped CPM"
        )]
    ))
    p_old_new <- ggplot(
        old_new_plot,
        aes(relative_position, delta_delta, color = normalization)
    ) +
        geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
        geom_vline(xintercept = 0, color = "grey40", linetype = 2) +
        geom_line(linewidth = 1.0) +
        scale_color_manual(values = normalization_colors) +
        scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
        labs(
            title = paste0(config$label, ": normalization comparison"),
            subtitle = "Site-weighted delta-delta with p = 0.1 on each normalized-count scale",
            x = "Relative position (nt)",
            y = "324 enrichment - 328 enrichment", color = "Normalization",
            caption = paste0("Boundary: ", config$left, " | ", config$right)
        ) +
        theme_publication()
    p_old_new <- add_side_labels(p_old_new, config, x_min, x_max)
    save_plot(
        p_old_new,
        file.path(figure_base, "old_vs_common_cpm_delta_delta")
    )

    weighting_plot <- rbindlist(list(
        common_profile[, .(
            relative_position, delta_delta, weighting = "Site-weighted"
        )],
        aggregate_profile[, .(
            relative_position, delta_delta,
            weighting = "Read-weighted aggregate"
        )]
    ))
    p_weighting <- ggplot(
        weighting_plot,
        aes(relative_position, delta_delta, color = weighting)
    ) +
        geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
        geom_vline(xintercept = 0, color = "grey40", linetype = 2) +
        geom_line(linewidth = 1.0) +
        scale_color_manual(values = weighting_colors) +
        scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
        labs(
            title = paste0(config$label, ": site- versus read-weighted delta-delta"),
            subtitle = "Common primary-mapped CPM normalization",
            x = "Relative position (nt)",
            y = "324 enrichment - 328 enrichment", color = "Weighting",
            caption = paste0("Boundary: ", config$left, " | ", config$right)
        ) +
        theme_publication()
    p_weighting <- add_side_labels(p_weighting, config, x_min, x_max)
    save_plot(
        p_weighting,
        file.path(figure_base, "site_vs_read_weighted_delta_delta")
    )

    p_coverage_profile <- ggplot(
        stratified_profile,
        aes(
            relative_position, delta_delta,
            color = coverage_quartile, fill = coverage_quartile
        )
    ) +
        geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
        geom_vline(xintercept = 0, color = "grey40", linetype = 2) +
        geom_ribbon(
            aes(ymin = ci_low, ymax = ci_high),
            alpha = 0.10, color = NA
        ) +
        geom_line(linewidth = 0.95) +
        scale_color_manual(values = quartile_colors) +
        scale_fill_manual(values = quartile_colors) +
        scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
        labs(
            title = paste0(config$label, ": delta-delta by coverage quartile"),
            subtitle = "Quartiles use total normalized window coverage across all eight libraries",
            x = "Relative position (nt)",
            y = "324 enrichment - 328 enrichment",
            color = "Coverage", fill = "Coverage",
            caption = paste0("Boundary: ", config$left, " | ", config$right)
        ) +
        theme_publication()
    p_coverage_profile <- add_side_labels(
        p_coverage_profile, config, x_min, x_max
    )
    save_plot(
        p_coverage_profile,
        file.path(figure_base, "coverage_stratified_delta_delta_metaprofile")
    )

    coverage_regions[, coverage_quartile := factor(
        coverage_quartile, levels = quartile_labels
    )]
    p_coverage_region <- ggplot(
        coverage_weighting_regions,
        aes(
            coverage_quartile, mean_delta_delta,
            color = coverage_quartile, shape = weighting
        )
    ) +
        geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
        geom_errorbar(
            data = coverage_weighting_regions[weighting == "Site-weighted"],
            aes(ymin = mean_ci_low, ymax = mean_ci_high),
            width = 0.12, linewidth = 0.55
        ) +
        geom_point(size = 2.4) +
        facet_wrap(~ region, scales = "free_y") +
        scale_color_manual(values = quartile_colors, guide = "none") +
        scale_shape_manual(values = c(
            "Site-weighted" = 16,
            "Read-weighted aggregate" = 17
        )) +
        labs(
            title = paste0(config$label, ": regional delta-delta by coverage"),
            subtitle = paste0(
                "Site-weighted means have 95% gene-cluster bootstrap CI; ",
                "triangles show pooled-read estimates"
            ),
            x = "Condition-independent coverage quartile",
            y = "Mean regional delta-delta",
            shape = "Weighting"
        ) +
        theme_publication() +
        theme(axis.text.x = element_text(angle = 25, hjust = 1))
    save_plot(
        p_coverage_region,
        file.path(figure_base, "coverage_stratified_regional_delta_delta"),
        width = 8.4, height = 4.7
    )

    replicate_regions[, sample_label := factor(
        sample_label,
        levels = c(
            "324 replicate A", "324 replicate B",
            "328 replicate A", "328 replicate B"
        )
    )]
    p_replicate_region <- ggplot(
        replicate_regions,
        aes(sample_label, regional_enrichment, color = construct, group = construct)
    ) +
        geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
        geom_line(linewidth = 0.55, alpha = 0.7) +
        geom_point(size = 2.5) +
        facet_wrap(~ region, scales = "free_y") +
        scale_color_manual(values = construct_colors, guide = "none") +
        labs(
            title = paste0(config$label, ": replicate-level regional enrichment"),
            subtitle = "Each point is one independently matched CLIP/SMI biological replicate",
            x = NULL, y = "Mean site-weighted log2(CLIP / SMI)"
        ) +
        theme_publication() +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
    save_plot(
        p_replicate_region,
        file.path(
            figures_dir, "replicate_regions",
            paste0(config$slug, "_replicate_region_estimates")
        ),
        width = 8.5, height = 4.8
    )

    rm(
        counts, primary_pair, primary_components, sensitivity_pair,
        sensitivity_components, sensitivity_profile, delta_sparse,
        region_sites, aggregate_profile, coverage_read_profile,
        coverage_read_regions, coverage_weighting_regions
    )
    invisible(gc())
}

normalization_comparison <- rbindlist(normalization_comparison, fill = TRUE)
sensitivity_summary <- rbindlist(sensitivity_summary, fill = TRUE)
replicate_region_all <- rbindlist(replicate_region_all, fill = TRUE)
replicate_construct_contrasts <- dcast(
    replicate_region_all,
    site_type + landmark + region + region_start + region_end + replicate ~
        construct,
    value.var = "regional_enrichment"
)
replicate_construct_contrasts[, construct_difference_324_minus_328 :=
    `324` - `328`]
side_summary_all <- rbindlist(side_summary_all, fill = TRUE)
coverage_profile_all <- rbindlist(coverage_profile_all, fill = TRUE)
coverage_region_all <- rbindlist(coverage_region_all, fill = TRUE)
coverage_read_profile_all <- rbindlist(
    coverage_read_profile_all, fill = TRUE
)
coverage_weighting_region_all <- rbindlist(
    coverage_weighting_region_all, fill = TRUE
)
pattern_summary <- rbindlist(pattern_summary, fill = TRUE)

write_tsv(
    normalization_comparison,
    file.path(results_dir, "normalization", "old_vs_common_cpm_summary.tsv")
)
write_tsv(
    sensitivity_summary,
    file.path(results_dir, "normalization", "common_cpm_pseudocount_sensitivity.tsv")
)
write_tsv(
    replicate_region_all,
    file.path(results_dir, "replicate_regions", "all_replicate_region_estimates.tsv")
)
write_tsv(
    replicate_construct_contrasts,
    file.path(
        results_dir, "replicate_regions",
        "replicate_matched_construct_contrasts.tsv"
    )
)
write_tsv(
    side_summary_all,
    file.path(results_dir, "side_of_splice", "all_side_comparison_summary.tsv")
)
write_tsv(
    coverage_profile_all,
    file.path(results_dir, "coverage_dependence", "all_coverage_stratified_profiles.tsv")
)
write_tsv(
    coverage_region_all,
    file.path(results_dir, "coverage_dependence", "all_coverage_stratified_regions.tsv")
)
write_tsv(
    coverage_read_profile_all,
    file.path(
        results_dir, "coverage_dependence",
        "all_coverage_stratified_read_weighted_profiles.tsv"
    )
)
write_tsv(
    coverage_weighting_region_all,
    file.path(
        results_dir, "coverage_dependence",
        "all_coverage_stratified_weighting_comparison.tsv"
    )
)
write_tsv(
    pattern_summary,
    file.path(results_dir, "metaprofiles", "site_vs_read_weighted_summary.tsv")
)

parameters <- data.table(
    parameter = c(
        "normalization", "normalization_formula", "primary_pseudocount_cpm",
        "sensitivity_pseudocount_cpm", "implicit_zero_policy",
        "construct_enrichment_summary", "delta_delta_definition",
        "coverage_statistic", "coverage_strata", "bootstrap_unit",
        "bootstrap_iterations", "bootstrap_seed", "gene_assignment"
    ),
    value = c(
        "common CPM based on primary mapped library depth across all samples",
        "raw_count / primary_mapped_reads * 1e6",
        primary_pseudocount, sensitivity_pseudocount,
        "absent sparse sample/site/bin rows are exact zeros",
        "mean of two within-construct matched replicate enrichments",
        "enrichment_324 - enrichment_328",
        "sum of normalized whole-window overlap counts across all eight samples",
        "within-landmark empirical coverage quartiles",
        "primary gene ID cluster", bootstrap_iterations, bootstrap_seed,
        "first listed gene_id for multi-gene landmark annotations"
    )
)
write_tsv(
    parameters,
    file.path(results_dir, "normalization", "analysis_parameters.tsv")
)

readme <- c(
    "# MSI2 common-CPM normalization and discrepancy investigation",
    "",
    "This analysis preserves `output/splicing_first_pass/` and writes a separate common-scale reanalysis. No FASTQ or BAM file was accessed.",
    "",
    "## Common normalization",
    "",
    "Every sample uses `raw_count / primary_mapped_reads * 1e6`. CLIP/SMI enrichment is calculated within the true construct/replicate pair with a primary pseudocount of 0.1 CPM and a sensitivity pseudocount of 0.5 CPM.",
    "",
    "## Implicit zeros",
    "",
    "The sparse bin table stores positive rows only. Missing sample/site/bin combinations are exact zeros. Every site-weighted denominator includes the full condition-independent retained-site universe.",
    "",
    "## Replicate regions",
    "",
    "Regional enrichment is calculated independently for the four matched CLIP/SMI pairs. Region means average the specified bins within each site and then weight retained sites equally. Bins are not treated as biological replicates.",
    "",
    "## Side of splice",
    "",
    "Donor and acceptor site tables contain paired exonic and intronic regional delta-delta values. Side difference is exonic minus intronic delta-delta, so negative values indicate a stronger 328-over-324 effect on the exonic side. Gene values are medians across sites assigned to the primary gene ID. Paired Wilcoxon tests, paired standardized mean differences, rank-biserial effects, and gene-cluster bootstrap intervals are reported.",
    "",
    "## Coverage dependence",
    "",
    "Coverage is condition-independent: whole-window overlap counts are converted to CPM per sample and summed across all eight samples. Balanced quartiles are defined within landmark type without using construct direction or CLIP enrichment; site ID breaks exact coverage ties deterministically. Coverage-stratified metaprofiles and primary regional summaries use site-weighted delta-delta. A pooled-read comparator within each quartile diagnoses the aggregate-weighting reversal.",
    "",
    "## Limitations",
    "",
    "Only two biological replicates exist per construct. Gene-cluster bootstrap intervals quantify variation across genes and do not create additional experimental replication. This analysis measures CLIP enrichment near landmarks, not PSI or transcript-level differential splicing. Annotated transcript 3' ends are not validated cleavage or PAS coordinates.",
    "",
    "## Reproduction",
    "",
    "```bash",
    paste(
        "Rscript analysis/investigate_landmark_normalization_and_coverage.R",
        "--input data/results --old-output output/splicing_first_pass",
        "--output output/splicing_common_cpm_investigation",
        paste0("--bootstrap ", bootstrap_iterations),
        paste0("--pseudocount ", primary_pseudocount),
        paste0("--sensitivity-pseudocount ", sensitivity_pseudocount),
        paste0("--seed ", bootstrap_seed)
    ),
    "```"
)
writeLines(readme, file.path(output_dir, "README.md"))
writeLines(
    capture.output(sessionInfo()),
    file.path(results_dir, "qc", "sessionInfo.txt")
)
message("Common-CPM investigation complete: ", normalizePath(output_dir))
