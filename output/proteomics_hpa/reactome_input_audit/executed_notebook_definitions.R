# Notebook fragment: clean-data
prot_df <- data_raw %>%
  transmute(
    gene = `Gene Symbol`,
    ensembl = `Ensembl Gene ID`,
    label = label,

    fc_324_gfp = `f_MSI2_324_WT /EGFP`,
    p_324_gfp  = `p_MSI2_324_WT /EGFP`,

    fc_328_gfp = `f_MSI2_328_WT /EGFP`,
    p_328_gfp  = `p_MSI2_328_WT /EGFP`,

    fc_40aa_gfp = `f_MSI2_328_40AA /EGFP`,
    p_40aa_gfp  = `p_MSI2_328_40AA /EGFP`,

    # original table has 40AA / 328
    fc_40aa_328 = `f_MSI2_328_40AA /MSI2_328_WT`,
    p_40aa_328  = `p_MSI2_328_40AA_MSI2_328_WT`,

    # derived direction: 328 / 40AA
    fc_328_40aa = 1 / `f_MSI2_328_40AA /MSI2_328_WT`,

    # original table has 328 / 324
    fc_328_324 = `f_MSI2_328_WT /MSI2_324_WT`,
    p_328_324  = `p_MSI2_328_WT /MSI2_324_WT`,

    # derived direction: 324 / 328
    fc_324_328 = 1 / `f_MSI2_328_WT /MSI2_324_WT`,

    EGFP_B = EGFP_Rep_B_n_n,
    EGFP_C = EGFP_Rep_C_n_n,

    WT324_A = FLAG_MSI2_324_WT_Rep_A_n_n,
    WT324_B = FLAG_MSI2_324_WT_Rep_B_n_n,
    WT324_C = FLAG_MSI2_324_WT_Rep_C_n_n,

    AA40_A = FLAG_MSI2_328_40AA_WT_Rep_A_n_n,
    AA40_B = FLAG_MSI2_328_40AA_WT_Rep_B_n_n,
    AA40_C = FLAG_MSI2_328_40AA_WT_Rep_C_n_n,

    WT328_A = FLAG_MSI2_328_WT_Rep_A_n_n,
    WT328_B = FLAG_MSI2_328_WT_Rep_B_n_n,
    WT328_C = FLAG_MSI2_328_WT_Rep_C_n_n
  ) %>%
  filter(!is.na(gene), gene != "")


# Notebook fragment: heatmap-functions
make_avg_kmeans_heatmap <- function(df,
                                    comparison = c("328_40AA", "324_328"),
                                    k = 3,
                                    seed = 123, #all eyes on me
                                    genes_to_label = NULL,
                                    title_text = NULL) {

  comparison <- match.arg(comparison)

  if (comparison == "328_40AA") {
    avg_df <- df %>%
      group_by(gene) %>%
      summarise(
        EGFP = mean(c(EGFP_B, EGFP_C), na.rm = TRUE),
        AA40 = mean(c(AA40_A, AA40_B, AA40_C), na.rm = TRUE),
        WT328 = mean(c(WT328_A, WT328_B, WT328_C), na.rm = TRUE),
        .groups = "drop"
      )

    condition_labels <- c("EGFP", "40AA", "328")

    if (is.null(title_text)) {
      title_text <- paste0("Full dataset averaged: EGFP vs 40AA vs 328, k = ", k)
    }
  }

  if (comparison == "324_328") {
    avg_df <- df %>%
      group_by(gene) %>%
      summarise(
        EGFP = mean(c(EGFP_B, EGFP_C), na.rm = TRUE),
        WT324 = mean(c(WT324_A, WT324_B, WT324_C), na.rm = TRUE),
        WT328 = mean(c(WT328_A, WT328_B, WT328_C), na.rm = TRUE),
        .groups = "drop"
      )

    condition_labels <- c("EGFP", "324", "328")

    if (is.null(title_text)) {
      title_text <- paste0("Full dataset averaged: EGFP vs 324 vs 328, k = ", k)
    }
  }

  heatmap_mat <- avg_df %>%
    column_to_rownames("gene") %>%
    as.matrix()

  heatmap_mat <- heatmap_mat[rowSums(is.na(heatmap_mat)) < ncol(heatmap_mat), ]

  heatmap_mat_scaled <- t(scale(t(heatmap_mat)))
  heatmap_mat_scaled <- heatmap_mat_scaled[complete.cases(heatmap_mat_scaled), ]

  set.seed(seed)
  km <- kmeans(heatmap_mat_scaled, centers = k)

  cluster_df <- tibble(
    gene = rownames(heatmap_mat_scaled),
    cluster = factor(km$cluster)
  )

  ordered_genes <- cluster_df %>%
    arrange(cluster, gene) %>%
    pull(gene)

  mat_ordered <- heatmap_mat_scaled[ordered_genes, ]
  clusters_ordered <- cluster_df$cluster[match(ordered_genes, cluster_df$gene)]

  gaps_row <- which(diff(as.numeric(clusters_ordered)) != 0)

  labels_row <- ifelse(
    rownames(mat_ordered) %in% genes_to_label,
    rownames(mat_ordered),
    ""
  )

  annotation_col <- data.frame(
    Condition = condition_labels
  )
  rownames(annotation_col) <- colnames(mat_ordered)

  annotation_row <- data.frame(
    Cluster = clusters_ordered
  )
  rownames(annotation_row) <- rownames(mat_ordered)

  heatmap_plot <- pheatmap(
    mat_ordered,
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    gaps_row = gaps_row,
    labels_row = labels_row,
    show_rownames = TRUE,
    fontsize_row = 9,
    scale = "none",
    main = title_text,
    silent = TRUE
  )

  list(
    avg_df = avg_df,
    matrix_scaled = heatmap_mat_scaled,
    cluster_df = cluster_df,
    plot = heatmap_plot
  )
}

# Notebook fragment: heatmaps-and-clusters
genes_to_label_all <- c("PABPC1", "IGF2BP3", "NCBP2", "MOV10", "EIF3B", "RPLP2", "RPL23A", "FAM98A", "EIF4G2", "YTHDF3",
  "EIF4A2", "DDX39A", "EXOSC2", "PAK2", "MAGEC1", "DPM1", "AIMP1", "DDX24", "WDR33")

genes_to_label_328 <- c("PABPC1", "IGF2BP3", "NCBP2", "MOV10", "EIF3B", "RPLP2", "RPL23A", "FAM98A", "EIF4G2", "YTHDF3")

heat_328_40aa <- make_avg_kmeans_heatmap(
  prot_df,
  comparison = "328_40AA",
  k = 3,
  genes_to_label = NULL
)

heat_328_40aa$plot

heat_324_328 <- make_avg_kmeans_heatmap(
  prot_df,
  comparison = "324_328",
  k = 3,
  genes_to_label = NULL
)

heat_324_328$plot

# Set these based on the heatmaps you like
cluster_of_interest_328 <- 3
cluster_of_interest_324 <- 2

genes_328_cluster <- heat_324_328$cluster_df %>%
  filter(cluster == cluster_of_interest_328) %>%
  pull(gene)

genes_324_cluster <- heat_324_328$cluster_df %>%
  filter(cluster == cluster_of_interest_324) %>%
  pull(gene)

genes_BOTH_cluster <- heat_324_328$cluster_df %>%
  filter(cluster == cluster_of_interest_324 | cluster == cluster_of_interest_328) %>%
  pull(gene)

# Notebook fragment: original-hit-definitions
filter_cluster_hits <- function(df, cluster_genes, fc_col, p_col,
                                fc_cutoff = 2, p_cutoff = 0.5) {
  df %>%
    filter(gene %in% cluster_genes) %>%
    filter(
      !is.na(.data[[fc_col]]),
      !is.na(.data[[p_col]]),
      .data[[fc_col]] >= fc_cutoff,
      .data[[p_col]] < p_cutoff
    ) %>%
    pull(gene) %>%
    unique()
}

genes_328_hits <- filter_cluster_hits(
  prot_df,
  cluster_genes = genes_328_cluster,
  fc_col = "fc_328_gfp",
  p_col = "p_328_gfp",
  fc_cutoff = 2,
  p_cutoff = 0.5
)

genes_324_hits <- filter_cluster_hits(
  prot_df,
  cluster_genes = genes_324_cluster,
  fc_col = "fc_324_gfp",
  p_col = "p_324_gfp",
  fc_cutoff = 2,
  p_cutoff = 0.5
)

filter_cluster_hits2 <- function(df,
                                cluster_genes,
                                fc_col,
                                p_col,
                                fc_cutoff = 2,
                                p_cutoff = 0.5,
                                direction = c("both", "up", "down"),
                                flip_direction = FALSE) {
  
  direction <- match.arg(direction)
  
  df %>%
    filter(gene %in% cluster_genes) %>%
    transmute(
      gene = gene,
      fc = .data[[fc_col]],
      pval = .data[[p_col]]
    ) %>%
    filter(
      !is.na(fc), !is.na(pval),
      is.finite(fc), is.finite(pval),
      fc > 0, pval > 0
    ) %>%
    mutate(
      log2FC = log2(fc),
      log2FC = if (flip_direction) -log2FC else log2FC
    ) %>%
    {
      out <- .
      
      if (direction == "both") {
        out <- out %>%
          filter(abs(log2FC) >= log2(fc_cutoff))
      }
      
      if (direction == "up") {
        out <- out %>%
          filter(log2FC >= log2(fc_cutoff))
      }
      
      if (direction == "down") {
        out <- out %>%
          filter(log2FC <= -log2(fc_cutoff))
      }
      
      out %>%
        filter(pval < p_cutoff)
    } %>%
    pull(gene) %>%
    unique()
}

genes_324v328_hits <- filter_cluster_hits2(
  prot_df,
  cluster_genes = genes_BOTH_cluster,
  fc_col = "fc_324_328",
  p_col = "p_328_324",
  fc_cutoff = 2,
  p_cutoff = 0.5,
  direction = "both"
)

split_cluster_by_fc <- function(df,
                                cluster_genes,
                                fc_col,
                                p_col,
                                fc_cutoff = 2,
                                p_cutoff = 0.5) {
  
  split_df <- df %>%
    filter(gene %in% cluster_genes) %>%
    transmute(
      gene = gene,
      fc = .data[[fc_col]],
      pval = .data[[p_col]]
    ) %>%
    filter(
      !is.na(gene), gene != "",
      !is.na(fc), !is.na(pval),
      is.finite(fc), is.finite(pval),
      fc > 0, pval > 0
    ) %>%
    mutate(
      log2FC = log2(fc),
      direction = case_when(
        log2FC >= log2(fc_cutoff) & pval < p_cutoff ~ "up_in_324",
        log2FC <= -log2(fc_cutoff) & pval < p_cutoff ~ "up_in_328",
        TRUE ~ "not_specific"
      )
    )
  
  list(
    table = split_df,
    up_in_324 = split_df %>% filter(direction == "up_in_324") %>% pull(gene) %>% unique(),
    up_in_328 = split_df %>% filter(direction == "up_in_328") %>% pull(gene) %>% unique(),
    not_specific = split_df %>% filter(direction == "not_specific") %>% pull(gene) %>% unique()
  )
}

split_324v328 <- split_cluster_by_fc(
  prot_df,
  cluster_genes = genes_BOTH_cluster,
  fc_col = "fc_324_328",
  p_col = "p_328_324",
  fc_cutoff = 2,
  p_cutoff = 0.5
)

genes_324_specific <- split_324v328$up_in_324
genes_328_specific <- split_324v328$up_in_328

length(genes_324_specific)
length(genes_328_specific)

length(genes_328_hits)
length(genes_324_hits)
length(genes_324v328_hits)

# Notebook fragment: original-entrez-mapping-function
map_symbols_to_entrez <- function(symbols) {
  bitr(
    unique(symbols),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  ) %>%
    distinct(SYMBOL, ENTREZID)
}


