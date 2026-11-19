#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ragg)
    library(scales)
})

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(flag, default) {
    i <- match(flag, args)
    if (is.na(i)) return(default)
    if (i == length(args)) stop("Missing value after ", flag)
    args[[i + 1L]]
}

input_dir <- arg_value("--input", "data/results")
output_dir <- arg_value("--output", "output/splicing_first_pass")
bootstrap_iterations <- as.integer(arg_value("--bootstrap", "1000"))
primary_pseudocount <- as.numeric(arg_value("--pseudocount", "0.1"))
sensitivity_pseudocount <- as.numeric(
    arg_value("--sensitivity-pseudocount", "0.5")
)
bootstrap_seed <- as.integer(arg_value("--seed", "328324"))

if (!is.finite(primary_pseudocount) || primary_pseudocount <= 0) {
    stop("Primary pseudocount must be positive")
}
if (!is.finite(sensitivity_pseudocount) || sensitivity_pseudocount <= 0) {
    stop("Sensitivity pseudocount must be positive")
}
if (bootstrap_iterations < 1000L) {
    stop("Use at least 1,000 gene-bootstrap iterations")
}

required_files <- file.path(
    input_dir,
    c(
        "sample_metadata.tsv", "landmark_sites.tsv.gz",
        "landmark_bin_counts.tsv.gz", "landmark_metaprofile_raw.tsv.gz",
        "landmark_extraction_QC.tsv", "README.md"
    )
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
    stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

results_dir <- file.path(output_dir, "results")
figures_dir <- file.path(output_dir, "figures")
result_subdirs <- c("normalized", "metaprofiles", "site_level", "gene_level", "qc")
figure_subdirs <- c("donors", "acceptors", "transcript_3prime_ends", "qc")
invisible(lapply(
    file.path(results_dir, result_subdirs),
    dir.create, recursive = TRUE, showWarnings = FALSE
))
invisible(lapply(
    file.path(figures_dir, figure_subdirs),
    dir.create, recursive = TRUE, showWarnings = FALSE
))

write_tsv <- function(x, path) {
    fwrite(as.data.table(x), path, sep = "\t", na = "NA", quote = FALSE)
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
group_colors <- c(
    `324 CLIP` = "#009E73", `324 SMI` = "#56B4E9",
    `328 CLIP` = "#D55E00", `328 SMI` = "#7B6BA8"
)

type_config <- list(
    splice_donor_5prime = list(
        slug = "donor", figure_dir = "donors", label = "5' splice donor",
        left = "EXON", right = "INTRON", expected_bins = 50L
    ),
    splice_acceptor_3prime = list(
        slug = "acceptor", figure_dir = "acceptors", label = "3' splice acceptor",
        left = "INTRON", right = "EXON", expected_bins = 50L
    ),
    annotated_transcript_3prime_end = list(
        slug = "transcript_3prime_end", figure_dir = "transcript_3prime_ends",
        label = "annotated transcript 3' end", left = "TRANSCRIPT",
        right = "DOWNSTREAM", expected_bins = 10L
    )
)

metadata <- fread(file.path(input_dir, "sample_metadata.tsv"))
required_metadata_columns <- c(
    "sample_id", "construct", "assay", "biological_replicate", "pair_id",
    "total_primary_mapped_reads", "star_uniquely_mapped_reads",
    "star_uniquely_mapped_percent", "existing_library_size_normalization_factor",
    "existing_normalization_scheme"
)
if (!all(required_metadata_columns %in% names(metadata))) {
    stop("sample_metadata.tsv is missing required columns")
}
metadata[, construct := as.character(construct)]
if (!setequal(metadata$construct, c("324", "328"))) {
    stop("Expected constructs 324 and 328")
}
if (!setequal(metadata$assay, c("CLIP", "SMI"))) {
    stop("Expected assays CLIP and SMI")
}
if (nrow(metadata) != 8L || anyDuplicated(metadata$sample_id)) {
    stop("Expected eight unique samples")
}
if (anyNA(metadata$existing_library_size_normalization_factor) ||
        any(metadata$existing_library_size_normalization_factor <= 0)) {
    stop("Existing positive library-size normalization factors are required")
}
if (!all(metadata$existing_normalization_scheme == "full_library")) {
    stop("Expected existing full_library normalization for every sample")
}

metadata[, normalization_reference_primary_mapped_reads :=
    exp(mean(log(total_primary_mapped_reads))), by = construct]
metadata[, recomputed_factor :=
    total_primary_mapped_reads / normalization_reference_primary_mapped_reads]
if (max(abs(
    metadata$recomputed_factor -
        metadata$existing_library_size_normalization_factor
)) > 1e-12) {
    stop(
        "Copied normalization factors do not match the construct-specific ",
        "primary-mapped geometric means used in the existing analyses"
    )
}
normalization_references <- unique(metadata[, .(
    construct, normalization_reference_primary_mapped_reads
)])
setorder(normalization_references, construct)

normalization_qc <- metadata[, .(
    sample = sample_id,
    construct,
    assay,
    replicate = biological_replicate,
    pair_id,
    primary_mapped_reads = total_primary_mapped_reads,
    uniquely_mapped_reads = star_uniquely_mapped_reads,
    uniquely_mapped_percent = star_uniquely_mapped_percent,
    normalization_factor = existing_library_size_normalization_factor,
    normalization_reference_primary_mapped_reads,
    normalization_method = paste0(
        "existing construct-specific full-library factor: primary mapped reads / ",
        "geometric mean within the four-library construct comparison"
    )
)]
setorder(normalization_qc, construct, replicate, assay)
write_tsv(
    normalization_qc,
    file.path(results_dir, "qc", "normalization_qc.tsv")
)

parameter_table <- data.table(
    parameter = c(
        "normalization_method", "normalization_formula",
        "normalization_reference_324", "normalization_reference_328",
        "primary_pseudocount_normalized_counts",
        "sensitivity_pseudocount_normalized_counts",
        "bootstrap_unit", "bootstrap_iterations", "bootstrap_seed",
        "construct_summary", "delta_delta_definition", "implicit_zero_policy"
    ),
    value = c(
        "existing project-specific construct-level full mapped-library factors",
        "normalized_count = raw_count / normalization_factor",
        sprintf(
            "%.12f primary mapped reads",
            normalization_references[construct == "324",
                                     normalization_reference_primary_mapped_reads]
        ),
        sprintf(
            "%.12f primary mapped reads",
            normalization_references[construct == "328",
                                     normalization_reference_primary_mapped_reads]
        ),
        primary_pseudocount, sensitivity_pseudocount,
        "primary gene ID; all sites assigned to the first listed gene ID",
        bootstrap_iterations, bootstrap_seed,
        "mean of two matched replicate-level CLIP/SMI log2 enrichments",
        "enrichment_324 - enrichment_328",
        "absent sample/site/bin sparse rows are exact zeros"
    )
)
write_tsv(
    parameter_table,
    file.path(results_dir, "normalized", "analysis_parameters.tsv")
)

pair_map <- metadata[, {
    if (sum(assay == "CLIP") != 1L || sum(assay == "SMI") != 1L) {
        stop("Every pair_id must have exactly one CLIP and one SMI sample")
    }
    .(
        clip_sample = sample_id[assay == "CLIP"],
        smi_sample = sample_id[assay == "SMI"]
    )
}, by = .(construct, biological_replicate, pair_id)]
if (nrow(pair_map) != 4L || any(pair_map[, .N, by = construct]$N != 2L)) {
    stop("Expected two matched CLIP/SMI pairs per construct")
}

sites_all <- fread(file.path(input_dir, "landmark_sites.tsv.gz"))
sites <- sites_all[retained_for_bin_counts == "yes"]
if (nrow(sites) == 0L || anyDuplicated(sites$site_id)) {
    stop("Retained landmark site table is empty or has duplicate site IDs")
}
sites[, primary_gene_id := sub(";.*$", "", gene_id)]
if (any(!nzchar(sites$primary_gene_id))) {
    stop("Every retained site must have a nonempty gene ID")
}
sites[, multiple_gene_ids := grepl(";", gene_id, fixed = TRUE)]

bin_counts <- fread(file.path(input_dir, "landmark_bin_counts.tsv.gz"))
if (any(bin_counts$raw_overlap_count <= 0L)) {
    stop("Sparse bin table must contain positive counts only")
}
if (!setequal(bin_counts$sample_id, metadata$sample_id)) {
    stop("Sparse bin table sample IDs do not match metadata")
}
if (anyDuplicated(bin_counts[, .(sample_id, site_id, relative_bin_start)])) {
    stop("Sparse bin table has duplicate sample/site/bin keys")
}
factor_lookup <- setNames(
    metadata$existing_library_size_normalization_factor,
    metadata$sample_id
)
bin_counts[, normalized_count :=
    raw_overlap_count / factor_lookup[sample_id]]

metaprofile_qc <- fread(file.path(input_dir, "landmark_metaprofile_raw.tsv.gz"))
reconciled <- bin_counts[, .(
    recomputed_sum_raw_overlap_count = sum(raw_overlap_count)
), by = .(sample_id, site_type, relative_bin_start, relative_bin_end)]
metaprofile_check <- merge(
    metaprofile_qc,
    reconciled,
    by = c("sample_id", "site_type", "relative_bin_start", "relative_bin_end"),
    all.x = TRUE
)
metaprofile_check[is.na(recomputed_sum_raw_overlap_count),
                  recomputed_sum_raw_overlap_count := 0]
metaprofile_check[, difference :=
    recomputed_sum_raw_overlap_count - sum_raw_overlap_count]
if (any(metaprofile_check$difference != 0)) {
    stop("Sparse bin sums do not reconcile to landmark_metaprofile_raw.tsv.gz")
}
write_tsv(
    metaprofile_check,
    file.path(results_dir, "qc", "aggregate_metaprofile_reconciliation.tsv")
)

pair_enrichment <- function(counts, pseudocount) {
    output <- vector("list", nrow(pair_map))
    for (i in seq_len(nrow(pair_map))) {
        pm <- pair_map[i]
        clip <- counts[sample_id == pm$clip_sample, .(
            site_id, relative_bin_start, relative_bin_end,
            normalized_clip = normalized_count
        )]
        smi <- counts[sample_id == pm$smi_sample, .(
            site_id, relative_bin_start, relative_bin_end,
            normalized_smi = normalized_count
        )]
        joined <- merge(
            clip, smi,
            by = c("site_id", "relative_bin_start", "relative_bin_end"),
            all = TRUE
        )
        joined[is.na(normalized_clip), normalized_clip := 0]
        joined[is.na(normalized_smi), normalized_smi := 0]
        joined[, enrichment := log2(
            (normalized_clip + pseudocount) /
                (normalized_smi + pseudocount)
        )]
        joined[, `:=`(
            construct = pm$construct,
            replicate = pm$biological_replicate,
            pair_id = pm$pair_id
        )]
        output[[i]] <- joined
    }
    rbindlist(output, use.names = TRUE)
}

profile_from_pair_values <- function(pair_values, site_count) {
    bins <- sort(unique(pair_values$relative_bin_start))
    replicate_profile <- pair_values[, .(
        enrichment = sum(enrichment) / site_count
    ), by = .(construct, replicate, pair_id, relative_bin_start, relative_bin_end)]
    full_replicate_grid <- CJ(
        construct = c("324", "328"),
        replicate = c("A", "B"),
        relative_bin_start = bins,
        unique = TRUE
    )
    full_replicate_grid[, relative_bin_end := relative_bin_start + 10L]
    full_replicate_grid <- merge(
        full_replicate_grid,
        pair_map[, .(construct, replicate = biological_replicate, pair_id)],
        by = c("construct", "replicate"), all.x = TRUE
    )
    replicate_profile <- merge(
        full_replicate_grid, replicate_profile,
        by = c(
            "construct", "replicate", "pair_id",
            "relative_bin_start", "relative_bin_end"
        ),
        all.x = TRUE
    )
    replicate_profile[is.na(enrichment), enrichment := 0]
    replicate_profile[, relative_position :=
        (relative_bin_start + relative_bin_end) / 2]

    construct_values <- pair_values[, .(
        enrichment = sum(enrichment) / 2
    ), by = .(construct, site_id, relative_bin_start, relative_bin_end)]
    construct_profile <- construct_values[, .(
        enrichment = sum(enrichment) / site_count
    ), by = .(construct, relative_bin_start, relative_bin_end)]
    full_construct_grid <- CJ(
        construct = c("324", "328"), relative_bin_start = bins,
        unique = TRUE
    )
    full_construct_grid[, relative_bin_end := relative_bin_start + 10L]
    construct_profile <- merge(
        full_construct_grid, construct_profile,
        by = c("construct", "relative_bin_start", "relative_bin_end"),
        all.x = TRUE
    )
    construct_profile[is.na(enrichment), enrichment := 0]
    construct_profile[, relative_position :=
        (relative_bin_start + relative_bin_end) / 2]

    list(
        replicate = replicate_profile,
        construct_sparse = construct_values,
        construct_profile = construct_profile
    )
}

bootstrap_profiles <- function(construct_sparse, site_map, bins) {
    gene_sizes <- site_map[, .(site_count = .N), by = primary_gene_id]
    genes <- gene_sizes$primary_gene_id
    n_gene <- length(genes)
    n_bin <- length(bins)
    mats <- list()
    for (construct_value in c("324", "328")) {
        gene_bin <- construct_sparse[construct == construct_value][
            site_map[, .(site_id, primary_gene_id)],
            on = "site_id", nomatch = 0
        ][, .(value = sum(enrichment)),
          by = .(primary_gene_id, relative_bin_start)]
        mat <- matrix(
            0, nrow = n_gene, ncol = n_bin,
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
    boot_324 <- matrix(NA_real_, bootstrap_iterations, n_bin)
    boot_328 <- matrix(NA_real_, bootstrap_iterations, n_bin)
    set.seed(bootstrap_seed + n_bin)
    batch_size <- 100L
    starts <- seq.int(1L, bootstrap_iterations, by = batch_size)
    uniform_prob <- rep(1 / n_gene, n_gene)
    for (start in starts) {
        finish <- min(start + batch_size - 1L, bootstrap_iterations)
        n_batch <- finish - start + 1L
        weights <- rmultinom(n_batch, size = n_gene, prob = uniform_prob)
        denominator <- as.numeric(crossprod(gene_sizes$site_count, weights))
        boot_324[start:finish, ] <- sweep(
            crossprod(weights, mats[["324"]]), 1L, denominator, "/"
        )
        boot_328[start:finish, ] <- sweep(
            crossprod(weights, mats[["328"]]), 1L, denominator, "/"
        )
    }
    boot_delta <- boot_324 - boot_328
    ci <- function(x) {
        t(apply(x, 2L, quantile, probs = c(0.025, 0.975), names = FALSE))
    }
    ci_324 <- ci(boot_324)
    ci_328 <- ci(boot_328)
    ci_delta <- ci(boot_delta)
    data.table(
        relative_bin_start = bins,
        enrichment_324 = central_324,
        enrichment_324_ci_low = ci_324[, 1L],
        enrichment_324_ci_high = ci_324[, 2L],
        enrichment_328 = central_328,
        enrichment_328_ci_low = ci_328[, 1L],
        enrichment_328_ci_high = ci_328[, 2L],
        delta_delta = central_324 - central_328,
        delta_delta_ci_low = ci_delta[, 1L],
        delta_delta_ci_high = ci_delta[, 2L],
        bootstrap_gene_count = n_gene,
        bootstrap_iterations = bootstrap_iterations
    )
}

aggregate_profiles <- function(counts, site_count, pseudocount) {
    bins <- sort(unique(counts$relative_bin_start))
    sample_grid <- CJ(
        sample_id = metadata$sample_id,
        relative_bin_start = bins,
        unique = TRUE
    )
    sample_grid[, relative_bin_end := relative_bin_start + 10L]
    aggregate <- counts[, .(
        sum_raw_overlap_count = sum(raw_overlap_count),
        sum_normalized_count = sum(normalized_count)
    ), by = .(sample_id, relative_bin_start, relative_bin_end)]
    aggregate <- merge(
        sample_grid, aggregate,
        by = c("sample_id", "relative_bin_start", "relative_bin_end"),
        all.x = TRUE
    )
    aggregate[is.na(sum_raw_overlap_count), sum_raw_overlap_count := 0]
    aggregate[is.na(sum_normalized_count), sum_normalized_count := 0]
    aggregate[, normalized_mean_per_retained_site :=
        sum_normalized_count / site_count]
    aggregate <- merge(
        aggregate,
        metadata[, .(
            sample_id, construct, assay,
            replicate = biological_replicate, pair_id
        )],
        by = "sample_id", all.x = TRUE
    )
    aggregate[, relative_position :=
        (relative_bin_start + relative_bin_end) / 2]

    pair_aggregate <- dcast(
        aggregate,
        construct + replicate + pair_id + relative_bin_start +
            relative_bin_end + relative_position ~ assay,
        value.var = "sum_normalized_count",
        fill = 0
    )
    pair_aggregate[, enrichment := log2(
        (CLIP + pseudocount) / (SMI + pseudocount)
    )]
    construct_aggregate <- pair_aggregate[, .(
        enrichment = mean(enrichment),
        replicate_min = min(enrichment),
        replicate_max = max(enrichment)
    ), by = .(construct, relative_bin_start, relative_bin_end, relative_position)]
    wide <- dcast(
        construct_aggregate,
        relative_bin_start + relative_bin_end + relative_position ~ construct,
        value.var = "enrichment"
    )
    setnames(wide, c("324", "328"), c("enrichment_324", "enrichment_328"))
    wide[, delta_delta := enrichment_324 - enrichment_328]
    list(
        coverage = aggregate,
        pair = pair_aggregate,
        construct = construct_aggregate,
        wide = wide
    )
}

add_side_labels <- function(plot, config, x_min, x_max) {
    plot +
        annotate(
            "text", x = c((x_min + 0) / 2, (0 + x_max) / 2), y = Inf,
            label = c(config$left, config$right), vjust = 1.4,
            fontface = "bold", size = 3.1, color = "grey30"
        )
}

site_gene_counts <- list()
replicate_consistency_summary <- list()
pseudocount_sensitivity <- list()
overall_pattern_summary <- list()

for (type_name in names(type_config)) {
    config <- type_config[[type_name]]
    message("Analyzing ", config$label)
    site_map <- sites[site_type == type_name, .(
        site_id, gene_id, primary_gene_id, multiple_gene_ids
    )]
    valid_site_count <- sites_all[site_type == type_name, .N]
    site_count <- nrow(site_map)
    counts <- bin_counts[site_type == type_name]
    bins <- sort(unique(counts$relative_bin_start))
    if (length(bins) != config$expected_bins) {
        stop("Unexpected bin count for ", site_type)
    }
    if (!all(counts$site_id %chin% site_map$site_id)) {
        stop("Bin counts contain a site outside the retained site universe")
    }

    primary_pair <- pair_enrichment(counts, primary_pseudocount)
    primary_profile <- profile_from_pair_values(primary_pair, site_count)
    boot_profile <- bootstrap_profiles(
        primary_profile$construct_sparse, site_map, bins
    )
    boot_profile[, `:=`(
        relative_bin_end = relative_bin_start + 10L,
        relative_position = relative_bin_start + 5,
        site_type = type_name,
        retained_site_count = site_count,
        primary_gene_count = uniqueN(site_map$primary_gene_id),
        pseudocount = primary_pseudocount
    )]
    setcolorder(boot_profile, c(
        "site_type", "relative_bin_start", "relative_bin_end",
        "relative_position", "retained_site_count", "primary_gene_count",
        "pseudocount", setdiff(names(boot_profile), c(
            "site_type", "relative_bin_start", "relative_bin_end",
            "relative_position", "retained_site_count", "primary_gene_count",
            "pseudocount"
        ))
    ))
    write_tsv(
        boot_profile,
        file.path(
            results_dir, "metaprofiles",
            paste0(config$slug, "_site_weighted_metaprofile.tsv")
        )
    )

    replicate_profile <- copy(primary_profile$replicate)
    replicate_profile[, `:=`(
        site_type = type_name,
        retained_site_count = site_count,
        pseudocount = primary_pseudocount
    )]
    write_tsv(
        replicate_profile,
        file.path(
            results_dir, "metaprofiles",
            paste0(config$slug, "_replicate_enrichment_profiles.tsv")
        )
    )

    sensitivity_pair <- pair_enrichment(counts, sensitivity_pseudocount)
    sensitivity_profile <- profile_from_pair_values(sensitivity_pair, site_count)
    sensitivity_wide <- dcast(
        sensitivity_profile$construct_profile,
        relative_bin_start + relative_bin_end + relative_position ~ construct,
        value.var = "enrichment", fill = 0
    )
    setnames(
        sensitivity_wide, c("324", "328"),
        c("sensitivity_enrichment_324", "sensitivity_enrichment_328")
    )
    sensitivity_wide[, sensitivity_delta_delta :=
        sensitivity_enrichment_324 - sensitivity_enrichment_328]
    sensitivity_table <- merge(
        boot_profile[, .(
            relative_bin_start, relative_bin_end, relative_position,
            primary_enrichment_324 = enrichment_324,
            primary_enrichment_328 = enrichment_328,
            primary_delta_delta = delta_delta
        )],
        sensitivity_wide,
        by = c("relative_bin_start", "relative_bin_end", "relative_position")
    )
    sensitivity_table[, `:=`(
        site_type = type_name,
        primary_pseudocount = primary_pseudocount,
        sensitivity_pseudocount = sensitivity_pseudocount,
        delta_absolute_difference = abs(
            primary_delta_delta - sensitivity_delta_delta
        ),
        delta_sign_agreement = sign(primary_delta_delta) ==
            sign(sensitivity_delta_delta)
    )]
    pseudocount_sensitivity[[type_name]] <- sensitivity_table

    aggregate_primary <- aggregate_profiles(
        counts, site_count, primary_pseudocount
    )
    aggregate_table <- merge(
        aggregate_primary$wide,
        dcast(
            aggregate_primary$construct,
            relative_bin_start + relative_bin_end + relative_position ~ construct,
            value.var = c("replicate_min", "replicate_max")
        ),
        by = c("relative_bin_start", "relative_bin_end", "relative_position")
    )
    aggregate_table[, `:=`(
        site_type = type_name,
        retained_site_count = site_count,
        pseudocount = primary_pseudocount,
        analysis_weighting = "read-weighted aggregate of summed normalized counts"
    )]
    write_tsv(
        aggregate_table,
        file.path(
            results_dir, "metaprofiles",
            paste0(config$slug, "_read_weighted_aggregate_metaprofile.tsv")
        )
    )
    coverage_table <- copy(aggregate_primary$coverage)
    coverage_table[, `:=`(
        site_type = type_name,
        retained_site_count = site_count
    )]
    write_tsv(
        coverage_table,
        file.path(
            results_dir, "normalized",
            paste0(config$slug, "_normalized_coverage_qc.tsv")
        )
    )

    site_gene_counts[[type_name]] <- data.table(
        site_type = type_name,
        valid_site_count = valid_site_count,
        retained_site_count = site_count,
        retained_fraction = site_count / valid_site_count,
        primary_gene_count = uniqueN(site_map$primary_gene_id),
        multi_gene_annotation_site_count = sum(site_map$multiple_gene_ids),
        bootstrap_site_count = site_count,
        bootstrap_gene_count = uniqueN(site_map$primary_gene_id),
        bootstrap_gene_assignment = "first listed gene_id"
    )

    replicate_wide <- dcast(
        replicate_profile,
        construct + relative_bin_start + relative_bin_end + relative_position ~
            replicate,
        value.var = "enrichment", fill = 0
    )
    replicate_stats <- replicate_wide[, .(
        pearson_r = cor(A, B, method = "pearson"),
        spearman_rho = cor(A, B, method = "spearman"),
        mean_absolute_difference = mean(abs(A - B)),
        root_mean_square_difference = sqrt(mean((A - B)^2))
    ), by = construct]
    replicate_stats[, `:=`(
        site_type = type_name,
        bin_count = length(bins),
        pseudocount = primary_pseudocount
    )]
    replicate_consistency_summary[[type_name]] <- replicate_stats
    replicate_wide <- merge(
        replicate_wide, replicate_stats[, .(construct, pearson_r, spearman_rho)],
        by = "construct"
    )
    write_tsv(
        replicate_wide,
        file.path(
            results_dir, "qc",
            paste0(config$slug, "_replicate_consistency.tsv")
        )
    )

    overall_pattern_summary[[type_name]] <- data.table(
        site_type = type_name,
        mean_delta_delta_across_bins = mean(boot_profile$delta_delta),
        median_delta_delta_across_bins = median(boot_profile$delta_delta),
        maximum_delta_delta = max(boot_profile$delta_delta),
        maximum_delta_position = boot_profile[
            which.max(delta_delta), relative_position
        ],
        minimum_delta_delta = min(boot_profile$delta_delta),
        minimum_delta_position = boot_profile[
            which.min(delta_delta), relative_position
        ],
        bins_ci_excludes_zero_positive = sum(
            boot_profile$delta_delta_ci_low > 0
        ),
        bins_ci_excludes_zero_negative = sum(
            boot_profile$delta_delta_ci_high < 0
        )
    )

    x_min <- min(boot_profile$relative_bin_start)
    x_max <- max(boot_profile$relative_bin_end)
    figure_base <- file.path(figures_dir, config$figure_dir)

    coverage_plot_data <- coverage_table[, .(
        normalized_mean = mean(normalized_mean_per_retained_site)
    ), by = .(construct, assay, relative_position)]
    coverage_plot_data[, group := paste(construct, assay)]
    coverage_table[, group := paste(construct, assay)]
    p_coverage <- ggplot() +
        geom_line(
            data = coverage_table,
            aes(
                relative_position, normalized_mean_per_retained_site,
                group = sample_id, color = group
            ),
            alpha = 0.28, linewidth = 0.55
        ) +
        geom_line(
            data = coverage_plot_data,
            aes(relative_position, normalized_mean, color = group),
            linewidth = 1.05
        ) +
        geom_vline(xintercept = 0, color = "grey40", linetype = 2) +
        scale_color_manual(values = group_colors) +
        scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
        labs(
            title = paste0(config$label, ": normalized coverage QC"),
            subtitle = "Bold lines are two-replicate means; faint lines are individual libraries",
            x = "Relative position (nt)",
            y = "Normalized mean count per retained site",
            color = "Group",
            caption = paste0("Boundary: ", config$left, " | ", config$right)
        ) +
        theme_publication()
    p_coverage <- add_side_labels(p_coverage, config, x_min, x_max)
    save_plot(
        p_coverage,
        file.path(figure_base, "raw_normalized_coverage_QC")
    )

    enrichment_long <- rbindlist(list(
        boot_profile[, .(
            construct = "324", relative_position,
            enrichment = enrichment_324,
            ci_low = enrichment_324_ci_low,
            ci_high = enrichment_324_ci_high
        )],
        boot_profile[, .(
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
            data = replicate_profile,
            aes(relative_position, enrichment, group = interaction(construct, replicate)),
            alpha = 0.25, linewidth = 0.5
        ) +
        geom_line(linewidth = 1.05) +
        scale_color_manual(values = construct_colors) +
        scale_fill_manual(values = construct_colors) +
        scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
        labs(
            title = paste0(config$label, ": site-weighted CLIP/SMI enrichment"),
            subtitle = paste0(
                "Mean site-level log2 enrichment; 95% gene-bootstrap CI; pseudocount ",
                primary_pseudocount
            ),
            x = "Relative position (nt)",
            y = "Mean log2(CLIP / SMI)",
            color = "Construct", fill = "Construct",
            caption = paste0("Boundary: ", config$left, " | ", config$right)
        ) +
        theme_publication()
    p_enrichment <- add_side_labels(p_enrichment, config, x_min, x_max)
    save_plot(
        p_enrichment,
        file.path(figure_base, "site_weighted_enrichment_metaprofile")
    )

    p_delta <- ggplot(
        boot_profile,
        aes(relative_position, delta_delta)
    ) +
        geom_hline(yintercept = 0, color = "grey45", linewidth = 0.4) +
        geom_vline(xintercept = 0, color = "grey40", linetype = 2) +
        geom_ribbon(
            aes(ymin = delta_delta_ci_low, ymax = delta_delta_ci_high),
            fill = "#6A3D9A", alpha = 0.18
        ) +
        geom_line(color = "#6A3D9A", linewidth = 1.05) +
        scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
        labs(
            title = paste0(config$label, ": site-weighted delta-delta"),
            subtitle = "Positive values indicate greater CLIP/SMI enrichment for 324; 95% gene-bootstrap CI",
            x = "Relative position (nt)",
            y = "324 enrichment - 328 enrichment",
            caption = paste0("Boundary: ", config$left, " | ", config$right)
        ) +
        theme_publication()
    p_delta <- add_side_labels(p_delta, config, x_min, x_max)
    save_plot(
        p_delta,
        file.path(figure_base, "site_weighted_delta_delta_metaprofile")
    )

    aggregate_long <- melt(
        aggregate_table,
        id.vars = c("relative_position"),
        measure.vars = c("enrichment_324", "enrichment_328"),
        variable.name = "construct", value.name = "enrichment"
    )
    aggregate_long[, construct := sub("enrichment_", "", construct)]
    p_aggregate_enrichment <- ggplot(
        aggregate_long,
        aes(relative_position, enrichment, color = construct)
    ) +
        geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
        geom_vline(xintercept = 0, color = "grey40", linetype = 2) +
        geom_line(linewidth = 1.05) +
        scale_color_manual(values = construct_colors) +
        scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
        labs(
            title = paste0(config$label, ": read-weighted aggregate enrichment"),
            subtitle = "Ratio of summed normalized coverage; secondary descriptive analysis",
            x = "Relative position (nt)", y = "Aggregate log2(CLIP / SMI)",
            color = "Construct",
            caption = paste0("Boundary: ", config$left, " | ", config$right)
        ) +
        theme_publication()
    p_aggregate_enrichment <- add_side_labels(
        p_aggregate_enrichment, config, x_min, x_max
    )
    save_plot(
        p_aggregate_enrichment,
        file.path(figure_base, "read_weighted_aggregate_enrichment_metaprofile")
    )

    p_aggregate_delta <- ggplot(
        aggregate_table,
        aes(relative_position, delta_delta)
    ) +
        geom_hline(yintercept = 0, color = "grey45", linewidth = 0.4) +
        geom_vline(xintercept = 0, color = "grey40", linetype = 2) +
        geom_line(color = "#6A3D9A", linewidth = 1.05) +
        scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
        labs(
            title = paste0(config$label, ": read-weighted aggregate delta-delta"),
            subtitle = "Positive values indicate greater aggregate CLIP/SMI enrichment for 324",
            x = "Relative position (nt)",
            y = "324 enrichment - 328 enrichment",
            caption = paste0("Boundary: ", config$left, " | ", config$right)
        ) +
        theme_publication()
    p_aggregate_delta <- add_side_labels(
        p_aggregate_delta, config, x_min, x_max
    )
    save_plot(
        p_aggregate_delta,
        file.path(figure_base, "read_weighted_aggregate_delta_delta_metaprofile")
    )

    stats_labels <- replicate_stats[, .(
        construct,
        label = sprintf("Pearson r = %.2f\nSpearman rho = %.2f", pearson_r, spearman_rho),
        x = -Inf, y = Inf
    )]
    p_replicate <- ggplot(
        replicate_wide,
        aes(A, B, color = construct)
    ) +
        geom_abline(slope = 1, intercept = 0, color = "grey55", linetype = 2) +
        geom_point(size = 1.8, alpha = 0.8) +
        geom_text(
            data = stats_labels,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = -0.05, vjust = 1.1, size = 3.1
        ) +
        facet_wrap(~ construct) +
        scale_color_manual(values = construct_colors, guide = "none") +
        labs(
            title = paste0(config$label, ": replicate profile consistency"),
            subtitle = "Each point is one 10-nt bin in the site-weighted enrichment profile",
            x = "Replicate A mean log2(CLIP / SMI)",
            y = "Replicate B mean log2(CLIP / SMI)"
        ) +
        coord_equal() +
        theme_publication()
    save_plot(
        p_replicate,
        file.path(figure_base, "replicate_consistency"),
        width = 7.2, height = 4.5
    )

    rm(
        counts, primary_pair, primary_profile, sensitivity_pair,
        sensitivity_profile, aggregate_primary
    )
    invisible(gc())
}

site_gene_counts <- rbindlist(site_gene_counts)
replicate_consistency_summary <- rbindlist(replicate_consistency_summary)
pseudocount_sensitivity <- rbindlist(pseudocount_sensitivity, fill = TRUE)
overall_pattern_summary <- rbindlist(overall_pattern_summary)

write_tsv(
    site_gene_counts,
    file.path(results_dir, "qc", "site_gene_counts.tsv")
)
write_tsv(
    replicate_consistency_summary,
    file.path(results_dir, "qc", "replicate_consistency_summary.tsv")
)
write_tsv(
    pseudocount_sensitivity,
    file.path(results_dir, "qc", "pseudocount_sensitivity.tsv")
)
write_tsv(
    overall_pattern_summary,
    file.path(results_dir, "metaprofiles", "first_pass_pattern_summary.tsv")
)

sensitivity_summary <- pseudocount_sensitivity[, .(
    delta_profile_pearson_r = cor(
        primary_delta_delta, sensitivity_delta_delta, method = "pearson"
    ),
    delta_profile_spearman_rho = cor(
        primary_delta_delta, sensitivity_delta_delta, method = "spearman"
    ),
    sign_agreement_fraction = mean(delta_sign_agreement),
    mean_absolute_delta_difference = mean(delta_absolute_difference),
    maximum_absolute_delta_difference = max(delta_absolute_difference)
), by = site_type]
write_tsv(
    sensitivity_summary,
    file.path(results_dir, "qc", "pseudocount_sensitivity_summary.tsv")
)

readme <- c(
    "# MSI2 landmark enrichment first-pass analysis",
    "",
    "This directory contains the requested local first-pass analysis of the completed BAM-derived landmark bundle. No FASTQ or BAM file was accessed.",
    "",
    "## Inputs",
    "",
    "- `data/results/sample_metadata.tsv`",
    "- `data/results/landmark_sites.tsv.gz`",
    "- `data/results/landmark_bin_counts.tsv.gz`",
    "- `data/results/landmark_metaprofile_raw.tsv.gz` (reconciliation only)",
    "- `data/results/landmark_extraction_QC.tsv`",
    "",
    "## Implicit zeros",
    "",
    "`landmark_bin_counts.tsv.gz` contains nonzero rows only. Every absent sample/site/bin combination was treated as an exact zero. Site-weighted means use the full retained-site denominator, including sites with zero coverage in a given sample/bin.",
    "",
    "## Normalization",
    "",
    paste0(
        "Counts were divided by the existing project-specific full-library factor copied into `sample_metadata.tsv`. Factors were originally computed within each four-library construct comparison as primary mapped reads divided by the comparison-specific geometric mean (324: ",
        sprintf("%.6f", normalization_references[construct == "324", normalization_reference_primary_mapped_reads]),
        "; 328: ",
        sprintf("%.6f", normalization_references[construct == "328", normalization_reference_primary_mapped_reads]),
        "). Unique-mapping percentages were not used for normalization."
    ),
    "",
    "## Matched enrichment and delta-delta",
    "",
    paste0("CLIP and SMI were matched within construct and biological replicate using `pair_id`. Replicate-level enrichment is `log2((normalized CLIP + ", primary_pseudocount, ") / (normalized SMI + ", primary_pseudocount, "))`. The pseudocount is applied after normalization. Construct profiles are the mean of the two replicate-level profiles. Delta-delta is `enrichment_324 - enrichment_328`; positive values favor 324 and negative values favor 328."),
    "",
    "## Site weighting and uncertainty",
    "",
    paste0("The primary profiles first calculate enrichment for each site and bin, then average over all retained sites. The 95% intervals use ", format(bootstrap_iterations, big.mark = ","), " cluster-bootstrap iterations over primary gene IDs. All sites are assigned to the first listed gene ID for bootstrap clustering; `site_gene_counts.tsv` reports multi-gene annotation counts."),
    "",
    "## Read-weighted aggregate analysis",
    "",
    "Secondary aggregate profiles sum normalized counts within each sample/bin, add the documented pseudocount to those summed values, then calculate matched CLIP/SMI enrichment. These profiles answer where reads are concentrated and are not the primary site-weighted result.",
    "",
    "## Pseudocount sensitivity",
    "",
    paste0("The primary pseudocount is ", primary_pseudocount, " normalized counts. The site-weighted analysis was repeated with ", sensitivity_pseudocount, "; results are in `results/qc/pseudocount_sensitivity.tsv` and its summary."),
    "",
    "## Orientation and terminology",
    "",
    "All positions are transcript-oriented. Negative values are upstream and positive values are downstream. Donor plots label EXON | INTRON, acceptor plots label INTRON | EXON, and transcript-end plots label TRANSCRIPT | DOWNSTREAM. Transcript ends are described only as annotated transcript 3' ends, not cleavage or polyadenylation sites.",
    "",
    "## Site filtering",
    "",
    "The primary universe is the condition-independent `retained_for_bin_counts == yes` set from extraction. No direction-dependent or CLIP-specific filter was added.",
    "",
    "## Limitations",
    "",
    "Only two biological replicates are available per construct/assay. Gene-bootstrap intervals describe uncertainty across annotated genes, not experimental-replicate uncertainty. Mean site-level log ratios are strongly shrunk toward zero by site/bin combinations with zero in both matched libraries, as intended by the full retained-site universe. The analysis addresses CLIP enrichment around landmarks; it does not estimate PSI or establish differential splicing. Annotated transcript 3' ends are not validated cleavage or PAS coordinates.",
    "",
    "## Reproduction",
    "",
    "Run:",
    "",
    "```bash",
    paste(
        "Rscript analysis/analyze_landmark_splicing_first_pass.R",
        "--input data/results --output output/splicing_first_pass",
        paste0("--bootstrap ", bootstrap_iterations),
        paste0("--pseudocount ", primary_pseudocount),
        paste0("--sensitivity-pseudocount ", sensitivity_pseudocount),
        paste0("--seed ", bootstrap_seed)
    ),
    "```"
)
writeLines(readme, file.path(output_dir, "README.md"))

writeLines(capture.output(sessionInfo()), file.path(results_dir, "qc", "sessionInfo.txt"))
message("First-pass analysis complete: ", normalizePath(output_dir))
