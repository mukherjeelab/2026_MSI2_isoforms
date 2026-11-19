#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name) {
    index <- match(name, args)
    if (is.na(index) || index == length(args)) {
        stop("Missing required argument: ", name)
    }
    args[[index + 1]]
}

counts_path <- get_arg("--counts")
peaks_path <- get_arg("--peaks")
metadata_path <- get_arg("--metadata")
library_sizes_path <- get_arg("--library-sizes")
comparison_column <- get_arg("--comparison-column")
comparison <- get_arg("--comparison")
normalized_path <- get_arg("--normalized")
all_results_path <- get_arg("--all-results")
significant_path <- get_arg("--significant")
size_factors_path <- get_arg("--size-factors")
qc_path <- get_arg("--qc")
session_path <- get_arg("--session")
design_text <- get_arg("--design")
reference_level <- get_arg("--reference")
test_level <- get_arg("--test")
minimum_total_count <- as.integer(get_arg("--minimum-total-count"))
normalization <- get_arg("--normalization")
size_factor_type <- get_arg("--size-factor-type")
alpha <- as.numeric(get_arg("--alpha"))

if (!normalization %in% c("full_library", "peak_rle")) {
    stop("normalization must be full_library or peak_rle")
}

all_metadata <- read.delim(
    metadata_path,
    comment.char = "#",
    stringsAsFactors = FALSE,
    check.names = FALSE
)
required_metadata <- c("sample", "assay")
if (!all(required_metadata %in% colnames(all_metadata))) {
    stop("Metadata must contain columns: ", paste(required_metadata, collapse = ", "))
}
if (!comparison_column %in% colnames(all_metadata)) {
    stop("Comparison column is missing from metadata: ", comparison_column)
}
if (anyDuplicated(all_metadata$sample)) {
    stop("Metadata sample names must be unique")
}
metadata <- all_metadata[
    all_metadata[[comparison_column]] == comparison, , drop = FALSE
]
if (!nrow(metadata)) {
    stop("No metadata rows found for comparison: ", comparison)
}
if (!all(c(reference_level, test_level) %in% metadata$assay)) {
    stop(
        "The assay column must contain both '", reference_level,
        "' and '", test_level, "'"
    )
}

library_sizes <- read.delim(
    library_sizes_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
)
required_library_columns <- c("sample", "primary_mapped_alignments")
if (!all(required_library_columns %in% colnames(library_sizes))) {
    stop(
        "Library-size table must contain columns: ",
        paste(required_library_columns, collapse = ", ")
    )
}
if (anyDuplicated(library_sizes$sample) ||
    !setequal(library_sizes$sample, all_metadata$sample)) {
    stop("Library-size samples must uniquely match all metadata samples")
}
library_sizes <- library_sizes[
    match(metadata$sample, library_sizes$sample), , drop = FALSE
]
if (any(!is.finite(library_sizes$primary_mapped_alignments)) ||
    any(library_sizes$primary_mapped_alignments <= 0)) {
    stop("All primary mapped-alignment counts must be positive numbers")
}

featurecounts <- read.delim(
    counts_path,
    comment.char = "#",
    stringsAsFactors = FALSE,
    check.names = FALSE
)
if (ncol(featurecounts) < 7) {
    stop("featureCounts output contains no sample columns")
}

count_columns <- seq.int(7, ncol(featurecounts))
observed_samples <- basename(colnames(featurecounts)[count_columns])
observed_samples <- sub("\\.mapped\\.sorted\\.bam$", "", observed_samples)
if (anyDuplicated(observed_samples)) {
    stop("Duplicate sample names derived from featureCounts BAM columns")
}
if (!all(metadata$sample %in% observed_samples)) {
    stop(
        "Not all comparison samples were found in featureCounts. BAMs: ",
        paste(observed_samples, collapse = ", "),
        "; comparison metadata: ", paste(metadata$sample, collapse = ", ")
    )
}

selected_count_columns <- count_columns[match(metadata$sample, observed_samples)]
count_data <- as.matrix(featurecounts[, selected_count_columns, drop = FALSE])
storage.mode(count_data) <- "integer"
rownames(count_data) <- featurecounts$Geneid
colnames(count_data) <- metadata$sample

peak_metadata <- read.delim(peaks_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!setequal(rownames(count_data), peak_metadata$peak_id)) {
    stop("Peak identifiers differ between featureCounts and merged peak metadata")
}
peak_metadata <- peak_metadata[match(rownames(count_data), peak_metadata$peak_id), ]

metadata$assay <- relevel(factor(metadata$assay), ref = reference_level)
for (column in setdiff(colnames(metadata), c("sample", "assay"))) {
    if (is.character(metadata[[column]])) {
        metadata[[column]] <- factor(metadata[[column]])
    }
}
rownames(metadata) <- metadata$sample

design_formula <- as.formula(design_text)
design_variables <- all.vars(design_formula)
missing_variables <- setdiff(design_variables, colnames(metadata))
if (length(missing_variables)) {
    stop("Design variables missing from metadata: ", paste(missing_variables, collapse = ", "))
}
model_matrix <- model.matrix(design_formula, data = metadata)
if (qr(model_matrix)$rank < ncol(model_matrix)) {
    stop("The DESeq2 design matrix is not full rank; check pair, batch, and assay columns")
}

keep <- rowSums(count_data) >= minimum_total_count
if (!any(keep)) {
    stop("No peaks pass minimum_total_count = ", minimum_total_count)
}

dds <- DESeqDataSetFromMatrix(
    countData = count_data[keep, , drop = FALSE],
    colData = metadata,
    design = design_formula
)
if (normalization == "full_library") {
    full_library_factors <- library_sizes$primary_mapped_alignments
    full_library_factors <- full_library_factors /
        exp(mean(log(full_library_factors)))
    names(full_library_factors) <- library_sizes$sample
    sizeFactors(dds) <- full_library_factors[colnames(dds)]
} else {
    dds <- estimateSizeFactors(dds, type = size_factor_type)
}
dds <- DESeq(dds)

size_factor_table <- data.frame(
    sample = colnames(dds),
    primary_mapped_alignments = library_sizes$primary_mapped_alignments[
        match(colnames(dds), library_sizes$sample)
    ],
    size_factor = unname(sizeFactors(dds)),
    normalization = normalization,
    check.names = FALSE
)
write.table(
    size_factor_table, size_factors_path,
    sep = "\t", quote = FALSE, row.names = FALSE
)

result <- results(
    dds,
    contrast = c("assay", test_level, reference_level),
    alpha = alpha
)
result_table <- data.frame(
    peak_id = rownames(result),
    as.data.frame(result),
    check.names = FALSE
)
result_table <- cbind(
    peak_metadata[match(result_table$peak_id, peak_metadata$peak_id), , drop = FALSE],
    result_table[, setdiff(colnames(result_table), "peak_id"), drop = FALSE]
)
result_table <- result_table[
    order(is.na(result_table$padj), result_table$padj, result_table$pvalue),
]

write.table(
    result_table, all_results_path,
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
significant <- result_table[!is.na(result_table$padj) & result_table$padj < alpha, ]
write.table(
    significant, significant_path,
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

normalized_all <- sweep(count_data, 2, sizeFactors(dds), "/")
write.table(
    data.frame(peak_id = rownames(normalized_all), normalized_all, check.names = FALSE),
    normalized_path,
    sep = "\t", quote = FALSE, row.names = FALSE
)

pdf(qc_path, width = 8, height = 7)
plotDispEsts(dds, main = "DESeq2 dispersion estimates")
plotMA(
    result,
    alpha = alpha,
    main = paste0(test_level, " vs ", reference_level)
)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
pca <- prcomp(t(assay(vsd)))
percent_variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
assay_colors <- as.integer(metadata$assay)
plot(
    pca$x[, 1], pca$x[, 2],
    col = assay_colors, pch = 19,
    xlab = sprintf("PC1 (%.1f%%)", percent_variance[1]),
    ylab = sprintf("PC2 (%.1f%%)", percent_variance[2]),
    main = "VST peak-count PCA"
)
text(pca$x[, 1], pca$x[, 2], labels = metadata$sample, pos = 3, cex = 0.7)
legend(
    "topright", legend = levels(metadata$assay),
    col = seq_along(levels(metadata$assay)), pch = 19
)
dev.off()

capture.output(sessionInfo(), file = session_path)

message("Samples: ", nrow(metadata))
message("Comparison: ", comparison)
message("Total merged peaks: ", nrow(count_data))
message("Peaks tested: ", sum(keep))
message("Normalization: ", normalization)
message("Significant peaks (padj < ", alpha, "): ", nrow(significant))
message("Positive log2FoldChange means enrichment in ", test_level)
