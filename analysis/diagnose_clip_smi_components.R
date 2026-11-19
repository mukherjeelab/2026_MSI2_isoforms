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
primary_dir <- arg_value(
    "--primary-output", "output/splicing_common_cpm_investigation"
)
output_dir <- arg_value(
    "--output", "output/splicing_clip_smi_diagnostic"
)
pseudocount_cpm <- as.numeric(arg_value("--pseudocount", "0.1"))
stopifnot(is.finite(pseudocount_cpm), pseudocount_cpm > 0)

required_inputs <- file.path(input_dir, c(
    "sample_metadata.tsv", "landmark_sites.tsv.gz",
    "landmark_bin_counts.tsv.gz"
))
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
    stop("Missing required inputs: ", paste(missing_inputs, collapse = ", "))
}
if (!dir.exists(primary_dir)) {
    stop("Existing common-CPM output is required and will be read unchanged")
}

results_dir <- file.path(output_dir, "results")
figures_dir <- file.path(output_dir, "figures")
invisible(lapply(
    c(
        results_dir,
        file.path(results_dir, "qc"),
        file.path(figures_dir, c(
            "donors", "acceptors", "transcript_3prime_ends", "qc"
        ))
    ),
    dir.create, recursive = TRUE, showWarnings = FALSE
))

write_tsv <- function(x, path) {
    fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
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

theme_publication <- function(base_size = 9.5) {
    theme_classic(base_size = base_size) +
        theme(
            plot.title = element_text(face = "bold", size = rel(1.08)),
            plot.subtitle = element_text(color = "grey30"),
            plot.caption = element_text(color = "grey35", hjust = 0),
            axis.title = element_text(face = "bold"),
            panel.grid.major.y = element_line(
                color = "grey90", linewidth = 0.25
            ),
            panel.grid.minor = element_blank(),
            legend.position = "bottom",
            legend.title = element_text(face = "bold")
        )
}

construct_colors <- c(`324` = "#1B9E77", `328` = "#D95F02")
replicate_types <- c(A = "solid", B = "22")

type_config <- list(
    splice_donor_5prime = list(
        slug = "donor", label = "5' splice donor",
        figure_dir = "donors", left = "EXON", right = "INTRON",
        expected_bins = 50L,
        regions = data.table(
            region = c("Exon side", "Intron side", "Full window"),
            start = c(-50L, 0L, -250L),
            end = c(0L, 50L, 250L)
        )
    ),
    splice_acceptor_3prime = list(
        slug = "acceptor", label = "3' splice acceptor",
        figure_dir = "acceptors", left = "INTRON", right = "EXON",
        expected_bins = 50L,
        regions = data.table(
            region = c("Intron side", "Exon side", "Full window"),
            start = c(-50L, 0L, -250L),
            end = c(0L, 50L, 250L)
        )
    ),
    annotated_transcript_3prime_end = list(
        slug = "transcript_3prime_end",
        label = "annotated transcript 3' end",
        figure_dir = "transcript_3prime_ends",
        left = "TRANSCRIPT", right = "DOWNSTREAM",
        expected_bins = 10L,
        regions = data.table(
            region = c("Upstream", "Downstream", "Full window"),
            start = c(-50L, 0L, -50L),
            end = c(0L, 50L, 50L)
        )
    )
)

metadata <- fread(file.path(input_dir, "sample_metadata.tsv"))
metadata[, construct := as.character(construct)]
metadata[, library_millions := total_primary_mapped_reads / 1e6]
if (nrow(metadata) != 8L || anyDuplicated(metadata$sample_id)) {
    stop("Expected eight unique samples")
}
if (!setequal(metadata$construct, c("324", "328")) ||
        !setequal(metadata$assay, c("CLIP", "SMI"))) {
    stop("Unexpected construct or assay labels")
}
if (any(metadata$library_millions <= 0)) {
    stop("Primary mapped library depths must be positive")
}

sites <- fread(file.path(input_dir, "landmark_sites.tsv.gz"))[
    retained_for_bin_counts == "yes"
]
if (nrow(sites) != 237732L || anyDuplicated(sites$site_id)) {
    stop("Unexpected retained-site universe")
}
sites[, primary_gene_id := sub(";.*$", "", gene_id)]

bin_counts <- fread(file.path(input_dir, "landmark_bin_counts.tsv.gz"))
if (any(bin_counts$raw_overlap_count <= 0L) ||
        anyDuplicated(bin_counts[, .(
            sample_id, site_id, relative_bin_start
        )])) {
    stop("Sparse bin table must contain unique positive rows")
}
library_lookup <- setNames(metadata$library_millions, metadata$sample_id)
bin_counts[, normalized_count_cpm :=
    raw_overlap_count / library_lookup[sample_id]]

make_sample_profiles <- function(counts, bins, site_count) {
    observed <- counts[, .(
        mean_site_weighted_cpm = sum(normalized_count_cpm) / site_count,
        positive_site_count = .N,
        nonzero_site_fraction = .N / site_count,
        log_component = log2(pseudocount_cpm) +
            sum(
                log2(normalized_count_cpm + pseudocount_cpm) -
                    log2(pseudocount_cpm)
            ) / site_count
    ), by = .(sample_id, relative_bin_start, relative_bin_end)]

    grid <- CJ(
        sample_id = metadata$sample_id,
        relative_bin_start = bins,
        unique = TRUE
    )
    grid[, relative_bin_end := relative_bin_start + 10L]
    profile <- merge(
        grid, observed,
        by = c("sample_id", "relative_bin_start", "relative_bin_end"),
        all.x = TRUE
    )
    profile[is.na(mean_site_weighted_cpm), `:=`(
        mean_site_weighted_cpm = 0,
        positive_site_count = 0,
        nonzero_site_fraction = 0,
        log_component = log2(pseudocount_cpm)
    )]
    profile <- merge(
        profile,
        metadata[, .(
            sample_id, construct, assay,
            replicate = biological_replicate, pair_id
        )],
        by = "sample_id", all.x = TRUE
    )
    profile[, relative_position := relative_bin_start + 5]
    profile[]
}

make_diagnostic_table <- function(
    sample_profile, primary_profile, type_name, site_count
) {
    construct_profile <- sample_profile[, .(
        mean_site_weighted_cpm = mean(mean_site_weighted_cpm),
        mean_nonzero_site_fraction = mean(nonzero_site_fraction),
        mean_log_component = mean(log_component)
    ), by = .(
        construct, assay, relative_bin_start,
        relative_bin_end, relative_position
    )]

    cpm_wide <- dcast(
        construct_profile,
        relative_bin_start + relative_bin_end + relative_position ~
            construct + assay,
        value.var = "mean_site_weighted_cpm"
    )
    setnames(
        cpm_wide,
        c("324_CLIP", "328_CLIP", "324_SMI", "328_SMI"),
        c(
            "mean_324_CLIP_CPM", "mean_328_CLIP_CPM",
            "mean_324_SMI_CPM", "mean_328_SMI_CPM"
        ),
        skip_absent = TRUE
    )

    occupancy_wide <- dcast(
        construct_profile,
        relative_bin_start + relative_bin_end + relative_position ~
            construct + assay,
        value.var = "mean_nonzero_site_fraction"
    )
    setnames(
        occupancy_wide,
        c("324_CLIP", "328_CLIP", "324_SMI", "328_SMI"),
        c(
            "mean_fraction_nonzero_324_CLIP",
            "mean_fraction_nonzero_328_CLIP",
            "mean_fraction_nonzero_324_SMI",
            "mean_fraction_nonzero_328_SMI"
        ),
        skip_absent = TRUE
    )

    log_wide <- dcast(
        construct_profile,
        relative_bin_start + relative_bin_end + relative_position ~
            construct + assay,
        value.var = "mean_log_component"
    )
    log_wide[, `:=`(
        clip_log_component_324_minus_328 = `324_CLIP` - `328_CLIP`,
        smi_log_component_324_minus_328 = `324_SMI` - `328_SMI`
    )]
    log_wide[, `:=`(
        clip_contribution_to_delta = clip_log_component_324_minus_328,
        smi_contribution_to_delta = -smi_log_component_324_minus_328,
        delta_delta_from_components =
            clip_log_component_324_minus_328 -
                smi_log_component_324_minus_328
    )]

    diagnostic <- merge(
        cpm_wide,
        occupancy_wide,
        by = c(
            "relative_bin_start", "relative_bin_end", "relative_position"
        )
    )
    diagnostic <- merge(
        diagnostic,
        log_wide[, .(
            relative_bin_start, relative_bin_end, relative_position,
            clip_log_component_324_minus_328,
            smi_log_component_324_minus_328,
            clip_contribution_to_delta,
            smi_contribution_to_delta,
            delta_delta_from_components
        )],
        by = c(
            "relative_bin_start", "relative_bin_end", "relative_position"
        )
    )
    diagnostic <- merge(
        diagnostic,
        primary_profile[, .(
            relative_bin_start, relative_bin_end, relative_position,
            enrichment_324, enrichment_328, delta_delta
        )],
        by = c(
            "relative_bin_start", "relative_bin_end", "relative_position"
        )
    )
    diagnostic[, `:=`(
        landmark_type = type_name,
        relative_bin = relative_position,
        `324_minus_328_CLIP_CPM` =
            mean_324_CLIP_CPM - mean_328_CLIP_CPM,
        `324_minus_328_SMI_CPM` =
            mean_324_SMI_CPM - mean_328_SMI_CPM,
        retained_site_count = site_count
    )]
    setcolorder(diagnostic, c(
        "landmark_type", "relative_bin", "relative_bin_start",
        "relative_bin_end", "retained_site_count",
        "mean_324_CLIP_CPM", "mean_328_CLIP_CPM",
        "324_minus_328_CLIP_CPM",
        "mean_324_SMI_CPM", "mean_328_SMI_CPM",
        "324_minus_328_SMI_CPM",
        "mean_fraction_nonzero_324_CLIP",
        "mean_fraction_nonzero_328_CLIP",
        "mean_fraction_nonzero_324_SMI",
        "mean_fraction_nonzero_328_SMI",
        "enrichment_324", "enrichment_328", "delta_delta",
        "clip_log_component_324_minus_328",
        "smi_log_component_324_minus_328",
        "clip_contribution_to_delta", "smi_contribution_to_delta",
        "delta_delta_from_components"
    ))
    list(construct = construct_profile, diagnostic = diagnostic)
}

add_side_labels <- function(plot, config, x_min, x_max) {
    plot + annotate(
        "text", x = c(x_min / 2, x_max / 2), y = Inf,
        label = c(config$left, config$right), vjust = 1.25,
        fontface = "bold", size = 3, color = "grey30"
    )
}

signal_panel <- function(
    construct_profile, sample_profile, assay_value,
    config, y_limit, panel_letter
) {
    construct_data <- construct_profile[assay == assay_value]
    replicate_data <- sample_profile[assay == assay_value]
    plot <- ggplot(
        construct_data,
        aes(relative_position, mean_site_weighted_cpm, color = construct)
    ) +
        geom_vline(xintercept = 0, color = "grey40", linetype = 2) +
        geom_line(
            data = replicate_data,
            aes(
                relative_position, mean_site_weighted_cpm,
                group = interaction(construct, replicate)
            ),
            linewidth = 0.45, alpha = 0.28
        ) +
        geom_line(linewidth = 1.05) +
        scale_color_manual(values = construct_colors) +
        coord_cartesian(ylim = c(0, y_limit)) +
        labs(
            title = paste0(panel_letter, ". ", assay_value, "-only coverage"),
            subtitle = paste0(
                "Observed ", assay_value,
                if (assay_value == "SMI") " RNA/background" else " reads"
            ),
            x = "Relative position (nt)",
            y = paste0("Mean ", assay_value, " CPM per retained site"),
            color = "Construct"
        ) +
        theme_publication()
    add_side_labels(
        plot, config,
        min(construct_data$relative_bin_start),
        max(construct_data$relative_bin_end)
    )
}

all_diagnostics <- list()
all_sample_profiles <- list()
all_regional <- list()
all_extrema <- list()

for (type_index in seq_along(type_config)) {
    type_name <- names(type_config)[[type_index]]
    config <- type_config[[type_name]]
    message("Diagnosing ", config$label)

    site_map <- sites[site_type == type_name, .(
        site_id, primary_gene_id
    )]
    site_count <- nrow(site_map)
    counts <- bin_counts[site_type == type_name]
    if (!all(counts$site_id %chin% site_map$site_id)) {
        stop("Bin count outside retained universe for ", type_name)
    }
    bins <- sort(unique(counts$relative_bin_start))
    if (length(bins) != config$expected_bins) {
        stop("Unexpected bin count for ", type_name)
    }

    primary_profile_path <- file.path(
        primary_dir, "results", "metaprofiles",
        paste0(config$slug, "_common_cpm_site_weighted_metaprofile.tsv")
    )
    primary_replicate_path <- file.path(
        primary_dir, "results", "metaprofiles",
        paste0(config$slug, "_common_cpm_replicate_profiles.tsv")
    )
    primary_profile <- fread(primary_profile_path)
    primary_replicate <- fread(primary_replicate_path)
    if (nrow(primary_profile) != length(bins) ||
            unique(primary_profile$pseudocount_cpm) != pseudocount_cpm) {
        stop("Primary profile mismatch for ", type_name)
    }

    sample_profile <- make_sample_profiles(counts, bins, site_count)
    sample_profile[, `:=`(
        landmark_type = type_name,
        retained_site_count = site_count
    )]
    components <- make_diagnostic_table(
        sample_profile, primary_profile, type_name, site_count
    )
    construct_profile <- components$construct
    diagnostic <- components$diagnostic

    identity_error <- max(abs(
        diagnostic$delta_delta - diagnostic$delta_delta_from_components
    ))
    if (!is.finite(identity_error) || identity_error > 1e-10) {
        stop("Component identity failed for ", type_name, ": ", identity_error)
    }

    write_tsv(
        diagnostic,
        file.path(
            results_dir, paste0(config$slug, "_component_diagnostic.tsv")
        )
    )

    replicate_cpm <- dcast(
        sample_profile,
        relative_bin_start + relative_bin_end + relative_position ~
            assay + construct + replicate,
        value.var = "mean_site_weighted_cpm"
    )
    replicate_enrichment <- dcast(
        primary_replicate,
        relative_bin_start + relative_bin_end + relative_position ~
            construct + replicate,
        value.var = "enrichment"
    )
    setnames(
        replicate_enrichment,
        c("324_A", "324_B", "328_A", "328_B"),
        paste0("enrichment_", c("324_A", "324_B", "328_A", "328_B"))
    )
    replicate_qc <- merge(
        replicate_cpm, replicate_enrichment,
        by = c(
            "relative_bin_start", "relative_bin_end", "relative_position"
        )
    )
    replicate_qc[, `:=`(
        landmark_type = type_name,
        delta_delta_A = enrichment_324_A - enrichment_328_A,
        delta_delta_B = enrichment_324_B - enrichment_328_B,
        retained_site_count = site_count
    )]
    write_tsv(
        replicate_qc,
        file.path(
            results_dir, paste0(config$slug, "_replicate_qc_profiles.tsv")
        )
    )

    regional <- rbindlist(lapply(
        seq_len(nrow(config$regions)),
        function(i) {
            region_row <- config$regions[i]
            diagnostic[
                relative_bin_start >= region_row$start &
                    relative_bin_end <= region_row$end,
                lapply(.SD, mean),
                .SDcols = c(
                    "mean_324_CLIP_CPM", "mean_328_CLIP_CPM",
                    "324_minus_328_CLIP_CPM",
                    "mean_324_SMI_CPM", "mean_328_SMI_CPM",
                    "324_minus_328_SMI_CPM",
                    "mean_fraction_nonzero_324_CLIP",
                    "mean_fraction_nonzero_328_CLIP",
                    "mean_fraction_nonzero_324_SMI",
                    "mean_fraction_nonzero_328_SMI",
                    "enrichment_324", "enrichment_328", "delta_delta",
                    "clip_contribution_to_delta",
                    "smi_contribution_to_delta"
                )
            ][, `:=`(
                landmark_type = type_name,
                region = region_row$region,
                region_start = region_row$start,
                region_end = region_row$end
            )]
        }
    ))
    setcolorder(regional, c(
        "landmark_type", "region", "region_start", "region_end",
        setdiff(names(regional), c(
            "landmark_type", "region", "region_start", "region_end"
        ))
    ))
    write_tsv(
        regional,
        file.path(
            results_dir,
            paste0(config$slug, "_regional_component_summary.tsv")
        )
    )

    extremum_position <- if (type_name == "splice_donor_5prime") {
        -5
    } else if (type_name == "splice_acceptor_3prime") {
        15
    } else {
        diagnostic[which.min(delta_delta), relative_bin]
    }
    extremum <- diagnostic[relative_bin == extremum_position]
    extremum[, `:=`(
        landmark_label = config$label,
        position_reason = if (type_name %in% c(
            "splice_donor_5prime", "splice_acceptor_3prime"
        )) "requested splice-site delta-delta minimum" else
            "observed delta-delta minimum"
    )]

    all_diagnostics[[type_name]] <- diagnostic
    all_sample_profiles[[type_name]] <- sample_profile
    all_regional[[type_name]] <- regional
    all_extrema[[type_name]] <- extremum

    signal_max <- max(
        components$construct$mean_site_weighted_cpm,
        sample_profile$mean_site_weighted_cpm
    ) * 1.04
    p_clip <- signal_panel(
        construct_profile, sample_profile, "CLIP",
        config, signal_max, "A"
    )
    p_smi <- signal_panel(
        construct_profile, sample_profile, "SMI",
        config, signal_max, "B"
    )

    enrichment_long <- rbindlist(list(
        primary_profile[, .(
            construct = "324", relative_position,
            value = enrichment_324,
            ci_low = enrichment_324_ci_low,
            ci_high = enrichment_324_ci_high
        )],
        primary_profile[, .(
            construct = "328", relative_position,
            value = enrichment_328,
            ci_low = enrichment_328_ci_low,
            ci_high = enrichment_328_ci_high
        )]
    ))
    p_enrichment <- ggplot(
        enrichment_long,
        aes(relative_position, value, color = construct, fill = construct)
    ) +
        geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
        geom_vline(xintercept = 0, color = "grey40", linetype = 2) +
        geom_ribbon(
            aes(ymin = ci_low, ymax = ci_high),
            alpha = 0.14, color = NA
        ) +
        geom_line(linewidth = 1.05) +
        scale_color_manual(values = construct_colors) +
        scale_fill_manual(values = construct_colors) +
        labs(
            title = "C. CLIP/SMI enrichment",
            subtitle = "Fixed primary profile; p = 0.1 CPM",
            x = "Relative position (nt)",
            y = "Mean log2(CLIP / SMI)",
            color = "Construct", fill = "Construct"
        ) +
        theme_publication()
    p_enrichment <- add_side_labels(
        p_enrichment, config,
        min(primary_profile$relative_bin_start),
        max(primary_profile$relative_bin_end)
    )

    p_delta <- ggplot(
        primary_profile,
        aes(relative_position, delta_delta)
    ) +
        geom_hline(yintercept = 0, color = "grey45", linewidth = 0.45) +
        geom_vline(xintercept = 0, color = "grey40", linetype = 2) +
        geom_ribbon(
            aes(ymin = delta_delta_ci_low, ymax = delta_delta_ci_high),
            fill = "#6A3D9A", alpha = 0.16
        ) +
        geom_line(color = "#6A3D9A", linewidth = 1.05) +
        labs(
            title = "D. Delta-delta enrichment",
            subtitle = "324 enrichment minus 328 enrichment",
            x = "Relative position (nt)", y = "324 - 328 enrichment"
        ) +
        theme_publication()
    p_delta <- add_side_labels(
        p_delta, config,
        min(primary_profile$relative_bin_start),
        max(primary_profile$relative_bin_end)
    )

    main_figure <- (p_clip | p_smi) / (p_enrichment | p_delta) +
        plot_annotation(
            title = paste0(config$label, ": CLIP and SMI components"),
            subtitle = paste0(
                "Site-weighted common primary-mapped CPM; ",
                format(site_count, big.mark = ","),
                " retained landmarks; faint lines are biological replicates"
            ),
            caption = paste0(
                "CLIP-only is observed read distribution; SMI-only is matched ",
                "RNA/background. CPM differences alone are not enrichment."
            ),
            theme = theme(
                plot.title = element_text(face = "bold", size = 14),
                plot.subtitle = element_text(size = 10),
                plot.caption = element_text(size = 8, color = "grey35")
            )
        )
    save_plot(
        main_figure,
        file.path(
            figures_dir, config$figure_dir,
            paste0(config$slug, "_clip_smi_enrichment_comparison")
        )
    )

    qc_signal_long <- sample_profile[, .(
        panel = assay,
        construct, replicate, relative_position,
        value = mean_site_weighted_cpm
    )]
    qc_enrichment_long <- primary_replicate[, .(
        panel = "CLIP/SMI enrichment",
        construct, replicate, relative_position,
        value = enrichment
    )]
    qc_delta_long <- rbindlist(list(
        replicate_qc[, .(
            panel = "Delta-delta", construct = "324 - 328",
            replicate = "A", relative_position,
            value = delta_delta_A
        )],
        replicate_qc[, .(
            panel = "Delta-delta", construct = "324 - 328",
            replicate = "B", relative_position,
            value = delta_delta_B
        )]
    ))
    qc_long <- rbindlist(list(
        qc_signal_long, qc_enrichment_long, qc_delta_long
    ), fill = TRUE)
    qc_long[, panel := factor(
        panel,
        levels = c("CLIP", "SMI", "CLIP/SMI enrichment", "Delta-delta")
    )]
    p_qc <- ggplot(
        qc_long,
        aes(
            relative_position, value,
            color = construct, linetype = replicate,
            group = interaction(construct, replicate)
        )
    ) +
        geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
        geom_vline(xintercept = 0, color = "grey40", linetype = 2) +
        geom_line(linewidth = 0.8) +
        facet_wrap(~ panel, ncol = 2, scales = "free_y") +
        scale_color_manual(
            values = c(construct_colors, `324 - 328` = "#6A3D9A")
        ) +
        scale_linetype_manual(values = replicate_types) +
        labs(
            title = paste0(config$label, ": replicate-level component QC"),
            subtitle = paste0(
                "A/B are biological replicates; delta panels compare ",
                "same-letter construct profiles"
            ),
            x = "Relative position (nt)", y = "Panel-specific value",
            color = "Series", linetype = "Replicate"
        ) +
        theme_publication() +
        theme(strip.text = element_text(face = "bold"))
    save_plot(
        p_qc,
        file.path(
            figures_dir, "qc",
            paste0(config$slug, "_replicate_component_qc")
        ),
        width = 10.5, height = 7.5
    )

    rm(
        counts, sample_profile, construct_profile, components,
        diagnostic, regional, extremum, primary_profile, primary_replicate,
        replicate_cpm, replicate_enrichment, replicate_qc,
        enrichment_long, qc_signal_long, qc_enrichment_long,
        qc_delta_long, qc_long
    )
    invisible(gc())
}

all_diagnostics <- rbindlist(all_diagnostics, fill = TRUE)
all_sample_profiles <- rbindlist(all_sample_profiles, fill = TRUE)
all_regional <- rbindlist(all_regional, fill = TRUE)
all_extrema <- rbindlist(all_extrema, fill = TRUE)

write_tsv(
    all_diagnostics,
    file.path(results_dir, "all_landmark_component_diagnostics.tsv")
)
write_tsv(
    all_sample_profiles,
    file.path(results_dir, "all_sample_clip_smi_profiles.tsv")
)
write_tsv(
    all_regional,
    file.path(results_dir, "all_regional_component_summaries.tsv")
)
write_tsv(
    all_extrema,
    file.path(results_dir, "landmark_extremum_component_diagnostics.tsv")
)

parameters <- data.table(
    parameter = c(
        "retained_site_filter", "normalization", "implicit_zero_policy",
        "site_weighting", "construct_summary", "pseudocount_cpm",
        "primary_enrichment_source", "component_identity",
        "clip_only_interpretation", "smi_only_interpretation"
    ),
    value = c(
        "sum raw whole-window overlap counts across all 8 samples >= 1",
        "raw_count / primary_mapped_reads * 1e6",
        "missing sparse sample/site/bin rows are exact zeros",
        "sum sample/bin CPM divided by all retained sites",
        "arithmetic mean of biological replicates A and B",
        pseudocount_cpm,
        file.path(
            primary_dir, "results/metaprofiles/",
            "*_common_cpm_site_weighted_metaprofile.tsv"
        ),
        "delta = CLIP log-component difference - SMI log-component difference",
        "observed CLIP read distribution, not standalone binding evidence",
        "matched RNA/background distribution, not binding"
    )
)
write_tsv(parameters, file.path(results_dir, "analysis_parameters.tsv"))

readme <- c(
    "# CLIP and SMI component diagnostic",
    "",
    "This additive analysis leaves the primary common-CPM CLIP/SMI analysis unchanged.",
    "",
    "CLIP-only panels show observed CLIP CPM per retained landmark. SMI-only panels show the matched RNA/background CPM distribution. Neither is called enrichment alone.",
    "",
    "The exact component decomposition uses site-level log terms with the same 0.1 CPM pseudocount as the primary analysis. The CLIP contribution plus the signed SMI contribution equals the unchanged delta-delta profile to numerical precision.",
    "",
    "Every sample/site/bin absent from the sparse count table is treated as zero, and every mean uses the full retained-site denominator. The retained filter is one or more raw whole-window overlaps across all eight samples combined.",
    "",
    "## Reproduction",
    "",
    "```bash",
    paste(
        "Rscript analysis/diagnose_clip_smi_components.R",
        "--input data/results",
        "--primary-output output/splicing_common_cpm_investigation",
        "--output output/splicing_clip_smi_diagnostic",
        paste0("--pseudocount ", pseudocount_cpm)
    ),
    "```"
)
writeLines(readme, file.path(output_dir, "README.md"))
writeLines(
    capture.output(sessionInfo()),
    file.path(results_dir, "qc", "sessionInfo.txt")
)

message("CLIP/SMI component diagnostic complete: ", normalizePath(output_dir))
