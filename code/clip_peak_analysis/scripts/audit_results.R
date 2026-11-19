#!/usr/bin/env Rscript

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
audit_dir <- arg_value("--output", file.path(results_dir, "audit"))
alpha <- 0.05
pseudocount <- 0.5

read_tsv <- function(path, ...) {
    if (!file.exists(path)) {
        stop("Missing required file: ", path)
    }
    read.delim(
        path,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        comment.char = "",
        ...
    )
}

assert <- function(condition, message) {
    if (!isTRUE(condition)) {
        stop(message)
    }
}

write_tsv <- function(x, path) {
    write.table(
        x, path, sep = "\t", quote = FALSE,
        row.names = FALSE, na = "NA"
    )
}

metadata <- read_tsv(metadata_path)
required_metadata <- c("sample", "assay", "comparison", "pair")
assert(
    all(required_metadata %in% names(metadata)),
    "Metadata is missing required columns"
)
assert(!anyDuplicated(metadata$sample), "Metadata sample names are not unique")

counts <- read_tsv(file.path(results_dir, "peak_count_matrix.tsv"))
assert(names(counts)[1] == "peak_id", "Count matrix must begin with peak_id")
assert(!anyDuplicated(counts$peak_id), "Count matrix peak IDs are not unique")
count_samples <- setdiff(names(counts), "peak_id")
assert(
    identical(count_samples, metadata$sample),
    "Count matrix columns are not in metadata sample order"
)
count_matrix <- as.matrix(counts[, count_samples, drop = FALSE])
storage.mode(count_matrix) <- "numeric"
assert(all(is.finite(count_matrix)), "Count matrix contains non-finite values")
assert(all(count_matrix >= 0), "Count matrix contains negative values")
assert(all(count_matrix == floor(count_matrix)), "Count matrix contains non-integers")
rownames(count_matrix) <- counts$peak_id

peaks <- read_tsv(file.path(results_dir, "peaks.filtered_merged.tsv"))
assert(!anyDuplicated(peaks$peak_id), "Merged peak IDs are not unique")
assert(
    identical(peaks$peak_id, counts$peak_id),
    "Merged peaks and count matrix rows differ"
)

library_sizes <- read_tsv(file.path(results_dir, "bam_library_sizes.tsv"))
assert(
    setequal(library_sizes$sample, metadata$sample),
    "Library-size samples differ from metadata"
)
library_sizes <- library_sizes[match(metadata$sample, library_sizes$sample), ]

feature_summary <- read_tsv(file.path(results_dir, "featurecounts.tsv.summary"))
feature_samples <- basename(names(feature_summary)[-1])
feature_samples <- sub("\\.mapped\\.sorted\\.bam$", "", feature_samples)
assert(
    setequal(feature_samples, metadata$sample),
    "featureCounts summary samples differ from metadata"
)
assigned_row <- which(feature_summary[[1]] == "Assigned")
assert(length(assigned_row) == 1, "featureCounts summary lacks one Assigned row")
assigned <- as.numeric(feature_summary[assigned_row, -1])
names(assigned) <- feature_samples
assigned <- assigned[metadata$sample]
assert(
    identical(unname(assigned), unname(colSums(count_matrix))),
    "Count-matrix column sums do not equal featureCounts Assigned totals"
)

size_factor_rows <- list()
for (comparison in unique(metadata$comparison)) {
    size_factor_path <- file.path(
        results_dir, "deseq2", comparison, "size_factors.tsv"
    )
    factors <- read_tsv(size_factor_path)
    factors$comparison <- comparison
    size_factor_rows[[comparison]] <- factors
}
size_factors <- do.call(rbind, size_factor_rows)
assert(
    setequal(size_factors$sample, metadata$sample),
    "DESeq2 size-factor samples differ from metadata"
)

sample_qc <- metadata
sample_qc$primary_mapped_alignments <- library_sizes$primary_mapped_alignments
sample_qc$assigned_to_fixed_peaks <- as.integer(assigned)
sample_qc$assigned_fraction <- assigned /
    library_sizes$primary_mapped_alignments
sample_qc$nonzero_peaks <- colSums(count_matrix > 0)
sample_qc$size_factor <- size_factors$size_factor[
    match(sample_qc$sample, size_factors$sample)
]
sample_qc$normalization <- size_factors$normalization[
    match(sample_qc$sample, size_factors$sample)
]

comparison_rows <- list()
results_by_comparison <- list()
for (comparison in unique(metadata$comparison)) {
    comparison_metadata <- metadata[metadata$comparison == comparison, ]
    samples <- comparison_metadata$sample
    result <- read_tsv(file.path(
        results_dir, "deseq2", comparison, "IP_vs_input.all.tsv"
    ))
    significant_file <- read_tsv(file.path(
        results_dir, "deseq2", comparison, "IP_vs_input.significant.tsv"
    ))
    normalized <- read_tsv(file.path(
        results_dir, "deseq2", comparison, "normalized_counts.tsv"
    ))

    expected_tested <- rownames(count_matrix)[rowSums(
        count_matrix[, samples, drop = FALSE]
    ) >= 10]
    assert(
        setequal(result$peak_id, expected_tested),
        paste(comparison, "tested peaks do not match the total-count filter")
    )
    expected_significant <- result$peak_id[
        !is.na(result$padj) & result$padj < alpha
    ]
    assert(
        identical(significant_file$peak_id, expected_significant),
        paste(comparison, "significant file does not match padj <", alpha)
    )
    assert(
        identical(normalized$peak_id, counts$peak_id),
        paste(comparison, "normalized counts do not contain all peaks in order")
    )
    assert(
        all(samples %in% names(normalized)),
        paste(comparison, "normalized counts are missing samples")
    )

    rownames(normalized) <- normalized$peak_id
    pairs <- unique(comparison_metadata$pair)
    assert(length(pairs) == 2, paste(comparison, "does not have exactly two pairs"))
    pair_ratios <- lapply(pairs, function(pair_name) {
        pair_metadata <- comparison_metadata[
            comparison_metadata$pair == pair_name, , drop = FALSE
        ]
        ip <- pair_metadata$sample[pair_metadata$assay == "IP"]
        input <- pair_metadata$sample[pair_metadata$assay == "input"]
        assert(
            length(ip) == 1 && length(input) == 1,
            paste(comparison, pair_name, "is not one IP matched to one input")
        )
        log2(
            (normalized[[ip]] + pseudocount) /
                (normalized[[input]] + pseudocount)
        )
    })
    names(pair_ratios) <- pairs
    tested_indices <- match(result$peak_id, normalized$peak_id)
    tested_cor_spearman <- cor(
        pair_ratios[[1]][tested_indices],
        pair_ratios[[2]][tested_indices],
        method = "spearman"
    )
    tested_cor_pearson <- cor(
        pair_ratios[[1]][tested_indices],
        pair_ratios[[2]][tested_indices],
        method = "pearson"
    )

    significant <- !is.na(result$padj) & result$padj < alpha
    enriched <- significant & result$log2FoldChange > 0
    depleted <- significant & result$log2FoldChange < 0
    enriched_indices <- match(result$peak_id[enriched], normalized$peak_id)
    consistent_enriched <- sum(
        pair_ratios[[1]][enriched_indices] > 0 &
            pair_ratios[[2]][enriched_indices] > 0
    )

    comparison_rows[[comparison]] <- data.frame(
        comparison = comparison,
        samples = nrow(comparison_metadata),
        tested_peaks = nrow(result),
        padj_nonmissing = sum(!is.na(result$padj)),
        significant_peaks = sum(significant),
        significantly_enriched = sum(enriched),
        significantly_depleted = sum(depleted),
        enriched_in_both_pairs = consistent_enriched,
        pair_ratio_pearson = tested_cor_pearson,
        pair_ratio_spearman = tested_cor_spearman,
        omniclip_score_lfc_spearman = cor(
            result$peak_score_max,
            result$log2FoldChange,
            method = "spearman"
        ),
        median_log2_fold_change = median(result$log2FoldChange),
        stringsAsFactors = FALSE
    )
    results_by_comparison[[comparison]] <- result
}

comparison_summary <- do.call(rbind, comparison_rows)
rownames(comparison_summary) <- NULL

comparison_names <- unique(metadata$comparison)
assert(
    length(comparison_names) == 2,
    "Construct comparison currently expects exactly two comparisons"
)
first <- results_by_comparison[[comparison_names[1]]][
    , c("peak_id", "log2FoldChange", "padj")
]
second <- results_by_comparison[[comparison_names[2]]][
    , c("peak_id", "log2FoldChange", "padj")
]
names(first)[-1] <- paste0(names(first)[-1], "_", comparison_names[1])
names(second)[-1] <- paste0(names(second)[-1], "_", comparison_names[2])
construct_comparison <- merge(first, second, by = "peak_id")

first_lfc <- construct_comparison[[paste0("log2FoldChange_", comparison_names[1])]]
second_lfc <- construct_comparison[[paste0("log2FoldChange_", comparison_names[2])]]
first_padj <- construct_comparison[[paste0("padj_", comparison_names[1])]]
second_padj <- construct_comparison[[paste0("padj_", comparison_names[2])]]
first_enriched <- !is.na(first_padj) & first_padj < alpha & first_lfc > 0
second_enriched <- !is.na(second_padj) & second_padj < alpha & second_lfc > 0
construct_comparison$enrichment_category <- ifelse(
    first_enriched & second_enriched,
    "enriched_in_both",
    ifelse(
        first_enriched,
        paste0("enriched_only_in_", comparison_names[1]),
        ifelse(
            second_enriched,
            paste0("enriched_only_in_", comparison_names[2]),
            "not_enriched_in_either"
        )
    )
)

dataset_summary <- data.frame(
    metric = c(
        "merged_peaks", "samples", "total_assigned_counts",
        "peaks_nonzero_in_any_sample", "peaks_zero_in_all_samples",
        "common_tested_peaks", "common_lfc_pearson",
        "common_lfc_spearman"
    ),
    value = c(
        nrow(count_matrix), ncol(count_matrix), sum(count_matrix),
        sum(rowSums(count_matrix) > 0), sum(rowSums(count_matrix) == 0),
        nrow(construct_comparison), cor(first_lfc, second_lfc),
        cor(first_lfc, second_lfc, method = "spearman")
    ),
    stringsAsFactors = FALSE
)

dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv(dataset_summary, file.path(audit_dir, "dataset_summary.tsv"))
write_tsv(sample_qc, file.path(audit_dir, "sample_qc.tsv"))
write_tsv(
    comparison_summary,
    file.path(audit_dir, "comparison_summary.tsv")
)
write_tsv(
    construct_comparison,
    file.path(audit_dir, "construct_comparison.common_tested.tsv")
)

message("All integrity checks passed.")
message("Audited ", nrow(count_matrix), " peaks across ", ncol(count_matrix), " samples.")
message("Wrote audit summaries to ", audit_dir)
