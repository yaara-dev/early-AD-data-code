#!/usr/bin/env Rscript
# code/celltyping.R
# Seurat cell typing pipeline (publication-ready).
#
# Running from the repository root:
#   Rscript code/celltyping.R
#
# Outputs are written to the folder specified by `out_dir` (default: "output").

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
  library(tibble)
})

# -----------------------------
# Helpers
# -----------------------------
read_counts_csv <- function(path) {
  stopifnot(file.exists(path))
  m <- read.csv(path, row.names = 1, check.names = FALSE)
  m <- as.matrix(m)
  m <- as(m, "dgCMatrix")
  return(m)
}

# -----------------------------
# Main pipeline
# -----------------------------
run_celltyping_pipeline <- function(
  counts_csv,
  markers_csv,
  out_dir = "output",
  sample_id = "sample",
  qc_threshold = 40,
  resolution_values = seq(0.6, 1.4, by = 0.03),
  dim_values = 2:4,
  consensus_support_min = 0.20,
  run_pca_qc_plots = TRUE,
  run_jackstraw = TRUE
) {

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # --- Load counts
  counts_matrix <- read_counts_csv(counts_csv)

  # --- Create Seurat object
  seurat_object <- CreateSeuratObject(
    counts = counts_matrix,
    project = sample_id,
    assay = "RNA",
    min.cells = 1,
    min.features = 1
  )
  seurat_object <- AddMetaData(
    object = seurat_object,
    metadata = Matrix::colSums(counts_matrix),
    col.name = "counts_per_cell"
  )

  # --- QC filter
  num_cells_before <- ncol(counts_matrix)
  seurat_object <- subset(seurat_object, subset = counts_per_cell > qc_threshold)
  num_cells_after <- ncol(seurat_object)

  message("Number of cells before threshold: ", num_cells_before)
  message("Number of cells after threshold:  ", num_cells_after)

  all.genes <- rownames(seurat_object)

  # --- Standard preprocessing
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst")
  seurat_object <- ScaleData(seurat_object, features = all.genes)
  seurat_object <- RunPCA(seurat_object)

  # --- PCA QC plots
  if (isTRUE(run_pca_qc_plots)) {
    pdf(file.path(out_dir, paste0("PCA_QC__", sample_id, ".pdf")), width = 10, height = 8)
    print(VizDimLoadings(seurat_object, dims = 1:5, reduction = "pca"))
    print(DimHeatmap(seurat_object, dims = 1:10, cells = 80, balanced = TRUE))
    print(ElbowPlot(seurat_object))
    dev.off()
  }

  # --- JackStraw (can be slow)
  if (isTRUE(run_jackstraw)) {
    seurat_object <- JackStraw(object = seurat_object, dims = 15)
    seurat_object <- ScoreJackStraw(seurat_object, dims = 1:15)
    p_js <- JackStrawPlot(object = seurat_object, dims = 1:15) + theme(legend.position = "bottom")
    ggsave(
      filename = file.path(out_dir, paste0("JackStraw__", sample_id, ".pdf")),
      plot = p_js,
      width = 8,
      height = 5
    )
  }

  # --- Load marker catalog (expected columns: Marker, cell_type)
  gene_cell_type <- read.csv(markers_csv, stringsAsFactors = FALSE)
  stopifnot(all(c("Marker", "cell_type") %in% colnames(gene_cell_type)))

  gene_cell_type_counts <- gene_cell_type %>%
    group_by(cell_type) %>%
    summarise(num_marker_genes = n(), .groups = "drop")

  # --- Vote matrix across (resolution, dims)
  cell_type_counts <- data.frame(
    matrix(
      ncol = length(resolution_values) * length(dim_values),
      nrow = length(colnames(seurat_object))
    ),
    stringsAsFactors = FALSE
  )
  colnames(cell_type_counts) <- paste(
    rep(resolution_values, each = length(dim_values)),
    rep(dim_values, times = length(resolution_values)),
    sep = "_"
  )
  rownames(cell_type_counts) <- colnames(seurat_object)

  # --- Iterations
  for (resolution_n in resolution_values) {
    for (dim_n in dim_values) {

      message("Running clustering with resolution=", resolution_n, " dims=", dim_n)

      seurat_object <- FindNeighbors(seurat_object, dims = 1:dim_n)
      seurat_object <- FindClusters(seurat_object, resolution = resolution_n)
      seurat_object <- RunUMAP(seurat_object, dims = 1:dim_n, return.model = TRUE)

      cluster_markers <- FindAllMarkers(
        object = seurat_object,
        only.pos = TRUE,
        min.pct = 0.25,
        logfc.threshold = 0.25
      )

      filtered_markers <- cluster_markers %>%
        filter(gene %in% gene_cell_type$Marker, avg_log2FC > 0.2)

      matched_markers <- filtered_markers %>%
        inner_join(gene_cell_type, by = c("gene" = "Marker"))

      cluster_annotations <- matched_markers %>%
        group_by(cluster, cell_type) %>%
        summarise(total_log2FC = sum(avg_log2FC, na.rm = TRUE), .groups = "drop")

      cluster_marker_counts <- matched_markers %>%
        group_by(cluster, cell_type) %>%
        summarise(num_marker_genes_in_cluster = n(), .groups = "drop")

      cluster_annotations_with_counts <- cluster_annotations %>%
        left_join(cluster_marker_counts, by = c("cluster", "cell_type")) %>%
        left_join(gene_cell_type_counts, by = "cell_type") %>%
        mutate(adjusted_log2FC = (total_log2FC * num_marker_genes_in_cluster) / num_marker_genes)

      cluster_cell_types <- cluster_annotations_with_counts %>%
        group_by(cluster) %>%
        filter(adjusted_log2FC == max(adjusted_log2FC, na.rm = TRUE)) %>%
        summarise(dominant_cell_type = paste(unique(cell_type), collapse = "/"), .groups = "drop")

      cluster_cell_types_vec <- setNames(
        cluster_cell_types$dominant_cell_type,
        as.character(cluster_cell_types$cluster)
      )

      cell_type_annotations <- cluster_cell_types_vec[as.character(seurat_object$seurat_clusters)]
      names(cell_type_annotations) <- colnames(seurat_object)

      key <- paste(resolution_n, dim_n, sep = "_")
      cell_type_counts[, key] <- cell_type_annotations
    }
  }

  # --- Consensus robustness
  n_runs <- ncol(cell_type_counts)

  row_mode <- function(v) {
    v <- na.omit(as.character(v))
    if (length(v) == 0) return(NA_character_)
    tb <- sort(table(v), decreasing = TRUE)
    names(tb)[1]
  }
  row_support <- function(v, denom) {
    v <- na.omit(as.character(v))
    if (length(v) == 0) return(NA_real_)
    max(table(v)) / denom
  }

  final_label <- apply(cell_type_counts, 1, row_mode)
  support <- apply(cell_type_counts, 1, row_support, denom = n_runs)
  ignored_flag <- is.na(final_label) | (support < consensus_support_min)

  seurat_object$final_cell_type <- final_label
  seurat_object$consensus_support <- support
  seurat_object$ignored_by_consensus <- ignored_flag

  df_sup <- tibble(
    cell = colnames(seurat_object),
    support = support,
    ignored = ignored_flag,
    sample = sample_id
  )

  p_sup <- ggplot(df_sup, aes(support)) +
    geom_histogram(bins = 30) +
    geom_vline(xintercept = consensus_support_min, linetype = "dashed") +
    labs(
      x = "Fraction of runs supporting the final label",
      y = "Cells",
      title = paste0("Consensus support — ", sample_id)
    ) +
    theme_classic()

  ggsave(
    file.path(out_dir, paste0("support_hist__", sample_id, ".pdf")),
    p_sup,
    width = 5.5,
    height = 4
  )

  ignored_tbl <- df_sup %>%
    summarise(ignored = sum(ignored), total = n(), frac = ignored / total)

  write.csv(ignored_tbl, file.path(out_dir, "ignored_summary.csv"), row.names = FALSE)

  p_den <- ggplot(df_sup, aes(support)) +
    geom_density() +
    geom_vline(xintercept = consensus_support_min, linetype = "dashed") +
    labs(
      x = "Fraction of runs supporting the final label",
      y = "Density",
      title = paste0("Consensus support density — ", sample_id)
    ) +
    theme_classic()

  ggsave(
    file.path(out_dir, paste0("support_density__", sample_id, ".pdf")),
    p_den,
    width = 5.5,
    height = 4
  )

  # --- Final call table (kept)
  num_steps <- length(resolution_values)
  ignored_cells_count <- 0

  chosen_cell_types <- apply(cell_type_counts, 1, function(x) {
    freq_table <- table(x)
    valid_cell_types <- names(freq_table[freq_table > consensus_support_min * num_steps])

    if (length(valid_cell_types) > 0) {
      names(freq_table)[which.max(freq_table[valid_cell_types])]
    } else {
      ignored_cells_count <<- ignored_cells_count + 1
      NA_character_
    }
  })

  message("Number of ignored cells because iterations: ", ignored_cells_count)

  cell_type_mapping <- c(
    "oligodendrocytes" = 1,
    "excitatory neurons" = 2,
    "inhibitory neurons" = 3,
    "astrocytes" = 4,
    "endothelial" = 5,
    "activated microglia" = 6,
    "microglia" = 7,
    "vasculature" = 8,
    "NA" = NA
  )

  chosen_cell_types_numbers <- unname(cell_type_mapping[chosen_cell_types])

  cell_type_df <- data.frame(
    cell_index = colnames(seurat_object),
    cell_type_number = chosen_cell_types_numbers,
    cell_type = chosen_cell_types,
    consensus_support = support,
    ignored_by_consensus = ignored_flag,
    stringsAsFactors = FALSE
  )

  write.csv(cell_type_df, file.path(out_dir, paste0("cell_type__", sample_id, ".csv")), row.names = FALSE)

  p_umap <- DimPlot(seurat_object, reduction = "umap", group.by = "final_cell_type", label = TRUE) +
    ggtitle("UMAP colored by final cell type")
  ggsave(file.path(out_dir, paste0("umap__", sample_id, ".pdf")), p_umap, width = 7, height = 5)

  message("Outputs written to: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))

  return(list(seurat_object = seurat_object, cell_type_df = cell_type_df))
}

# -----------------------------
# Run input example when executed directly
# -----------------------------
if (sys.nframe() == 0) {
  run_celltyping_pipeline(
    counts_csv  = "input_example/input_counts.csv",
    markers_csv = "input_example/input_markers.csv",
    out_dir     = "output",
    sample_id   = "input_example",
    qc_threshold = 40,
    resolution_values = seq(0.6, 1.0, by = 0.1),
    dim_values = 2:3,
    consensus_support_min = 0.20,
    run_pca_qc_plots = TRUE,
    run_jackstraw = FALSE
  )
}
