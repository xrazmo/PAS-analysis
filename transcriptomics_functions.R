# transcriptomics_functions.R  — FIXED VERSION
# ─────────────────────────────────────────────────────────────────────────────
# Changes from original (all tagged with FIX-N matching review document):
#
#   FIX-2.1  run_DESeq2_sva: normalize counts before passing to sva()
#   FIX-2.2  run_DESeq2_sva: sv_df scoping fixed; no-SV branch now safe
#   FIX-2.3  assess_interactions: added "Shared (additive, no interaction)"
#             category for genes significant in BOTH phage and antibiotic
#   FIX-2.4  run_lfcShrink: combo contrast now built as a named numeric vector
#             (valid lfcShrink input) instead of list(character_vector)
#   FIX-3.1  run_gsea_base: KEGG downloaded once outside map() loop
#   FIX-3.2  plot_gsea_heatmap: identical() checks converted to stopifnot()
#   FIX-3.3  run_lfcShrink: interaction coefficients (containing "." or ":")
#             now routed to ashr; main effects remain on apeglm
#   FIX-3.4  assess_interactions: f_combo added to file-existence check
#   FIX-3.5  plot_venn_diagram: dead filter (pattern "ab4phage1" never matched)
#             removed
#   FIX-4.1  run_lfcShrink: library() calls inside function body removed
#   FIX-4.3  run_lfcShrink: mutate(coef=...) now guarded against NULL coef_name
#   FIX-4.4  run_gsea_base: ifelse() for dir.create() replaced with if()
#   FIX-4.5  plot_lfc_chicklets: !!filter_col / !!cols → all_of() / any_of()
#   FIX-4.6  plot_lfc_chicklets: height_df computed after dummy rows to keep
#             panel proportions correct; dummy rows constructed more robustly
#   FIX-4.7  plot_pathway_dotplot: typo "Entrichment" → "Enrichment"
#   MISC     set.seed(42) added before sva() and GSEA() for reproducibility
#   MISC     deprecated size= for line geoms → linewidth= throughout
# ─────────────────────────────────────────────────────────────────────────────

library(DESeq2)
library(sva)
library(tidyverse)
library(ggVennDiagram)
library(ggpubr)
library(ggchicklet)
library(ggh4x)
library(ggtext)
library(readr)
library(clusterProfiler)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(stringr)
library(RColorBrewer)


# ─────────────────────────────────────────────────────────────────────────────
# load_transcriptomics_data  — unchanged
# ─────────────────────────────────────────────────────────────────────────────
load_transcriptomics_data <- function(counts_path,
                                      metadata_path,
                                      annotation_path  = NULL,
                                      metadata_filter_col     = "comment",
                                      metadata_filter_pattern = "missing",
                                      sample_id_col = "sample_id",
                                      fname_col     = "fname") {

  message("Loading count data from: ", counts_path)
  counts <- read.table(counts_path, header = TRUE, row.names = "Geneid", sep = ",") %>%
    data.matrix()
  message("  -> Loaded ", nrow(counts), " genes across ", ncol(counts), " samples")

  message("Loading metadata from: ", metadata_path)
  if (grepl("\\.xlsx?$", metadata_path, ignore.case = TRUE)) {
    metadata <- readxl::read_xlsx(metadata_path)
  } else {
    metadata <- read.csv(metadata_path)
  }

  if (!is.null(metadata_filter_col) && metadata_filter_col %in% colnames(metadata)) {
    metadata <- metadata %>%
      filter(!grepl(metadata_filter_pattern, .data[[metadata_filter_col]])) %>%
      dplyr::select(-all_of(metadata_filter_col))
  }

  if (sample_id_col %in% colnames(metadata)) {
    metadata <- metadata %>% arrange(.data[[sample_id_col]])
  }

  if (fname_col %in% colnames(metadata)) {
    if (!all(colnames(counts) %in% metadata[[fname_col]])) {
      missing_samples <- setdiff(colnames(counts), metadata[[fname_col]])
      warning("Some count samples not found in metadata: ",
              paste(missing_samples, collapse = ", "))
    }
    metadata <- metadata[match(colnames(counts), metadata[[fname_col]]), ]
    stopifnot("Metadata and counts don't match!" =
                all(colnames(counts) == metadata[[fname_col]]))
  }

  rownames(metadata) <- metadata[[fname_col]]

  metadata <- metadata %>% mutate(
    time_point = factor(time_point, levels = c("1h", "4h", "24h")),
    phage      = factor(phage,   levels = c(0, 1)),
    ab_conc    = factor(ab_conc, levels = c(0, 4, 8))
  )

  annotation <- NULL
  if (!is.null(annotation_path)) {
    message("Loading annotation from: ", annotation_path)
    annotation <- read.csv(annotation_path)
    message("  -> Loaded annotations for ", nrow(annotation), " genes")
  }

  message("\n\u2713 Data loading complete!")
  return(list(counts = counts, metadata = metadata, annotation = annotation))
}


# ─────────────────────────────────────────────────────────────────────────────
# run_DESeq2_sva
# FIX-2.1  Normalize counts before sva() — raw integers invalid input for SVA
# FIX-2.2  sv_df scoping fixed: sv_terms computed inside both branches so
#           line 114-equivalent never references an undefined variable
# MISC     set.seed(42) before sva() for reproducibility
# ─────────────────────────────────────────────────────────────────────────────
run_DESeq2_sva <- function(counts, metadata,
                           formula_str = "~ ab_conc + phage",
                           time_point  = NULL,
                           n_svs       = NULL) {

  if (!is.null(time_point)) {
    message("\n\u2192 Subsetting for time point: ", time_point)
    metadata <- metadata %>%
      filter(as.character(.data[["time_point"]]) == !!time_point)
    counts <- counts[, metadata$fname]
  }

  # Design matrices for SVA
  mod  <- model.matrix(as.formula(formula_str), data = metadata)
  mod0 <- model.matrix(~1, data = metadata)

  # FIX-2.1: normalize + log-transform before SVA
  # Raw integer counts cause library-size variation to dominate SVs.
  message("Normalizing counts for SVA input...")
  dds_pre <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = metadata,
    design    = as.formula(formula_str)
  )
  dds_pre    <- estimateSizeFactors(dds_pre)
  norm_mat   <- log1p(counts(dds_pre, normalized = TRUE))

  message("Estimating surrogate variables...")
  set.seed(42)  # MISC: reproducibility
  svobj <- sva(norm_mat, mod, mod0)

  # FIX-2.2: define sv_terms in BOTH branches so it is always available below
  if (!is.null(svobj$sv) && is.matrix(svobj$sv) && ncol(svobj$sv) > 0) {
    sv_df  <- as.data.frame(svobj$sv)
    colnames(sv_df) <- paste0("SV", seq_len(ncol(sv_df)))
    metadata  <- cbind(metadata, sv_df)
    sv_terms  <- paste(colnames(sv_df), collapse = " + ")
    message("Detected ", ncol(svobj$sv), " surrogate variable(s). Added to metadata.")
  } else {
    sv_terms <- NULL
    message("No surrogate variables detected. Proceeding without SVs.")
  }

  # Build design formula — include SVs if detected
  rhs            <- if (!is.null(sv_terms)) paste(sv_terms, "+ ab_conc * phage")
                    else "ab_conc * phage"
  design_formula <- as.formula(paste("~", rhs))
  message("DESeq2 design: ", deparse(design_formula))

  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = metadata,
    design    = design_formula
  )
  dds <- DESeq(dds)

  return(list(dds = dds, svobj = svobj, metadata = metadata))
}


# ─────────────────────────────────────────────────────────────────────────────
# run_lfcShrink
# FIX-2.4  Combo contrast now built as a named numeric vector — valid for ashr
# FIX-3.3  Interaction coefficients (names containing "." or ":") routed to
#           ashr; main-effect coefficients remain on apeglm
# FIX-4.1  Removed library() calls inside function body
# FIX-4.3  mutate(coef = ...) guarded: only added when coef_name is not NULL
# ─────────────────────────────────────────────────────────────────────────────
run_lfcShrink <- function(dds,
                          coef_name  = NULL,
                          contrast_c = NULL,
                          LFC        = NULL,
                          padj       = NULL,
                          out_file   = NULL,
                          annotation = NULL) {

  if (!is.null(contrast_c)) {
    # FIX-2.4: build a proper named numeric contrast vector.
    # contrast_c is a character vector of coefficient names to sum (+1 each).
    # lfcShrink(type="ashr") accepts a numeric vector over all model coefs.
    all_coefs   <- resultsNames(dds)
    missing_c   <- setdiff(contrast_c, all_coefs)
    if (length(missing_c) > 0) {
      stop("Contrast coefficients not found in DESeq2 object: ",
           paste(missing_c, collapse = ", "),
           "\nAvailable: ", paste(all_coefs, collapse = ", "))
    }
    contr_vec              <- setNames(rep(0L, length(all_coefs)), all_coefs)
    contr_vec[contrast_c]  <- 1L
    res <- lfcShrink(dds, contrast = contr_vec, type = "ashr")

  } else {
    # Single named coefficient
    if (!(coef_name %in% resultsNames(dds))) {
      stop("Coefficient '", coef_name,
           "' not found in DESeq2 object. Available names: ",
           paste(resultsNames(dds), collapse = ", "))
    }

    # FIX-3.3: route interaction terms (coefficient names containing "." or ":")
    # to ashr; main effects use apeglm
    is_interaction <- grepl("[.:]", coef_name)
    shrink_type    <- if (is_interaction) "ashr" else "apeglm"
    message("  LFC shrinkage: coef='", coef_name,
            "' | type='", shrink_type, "'")
    res <- lfcShrink(dds, coef = coef_name, type = shrink_type)
  }

  res <- res %>%
    as.data.frame() %>%
    rownames_to_column("gene")

  # FIX-4.3: only add coef column when coef_name is not NULL
  if (!is.null(coef_name)) {
    res <- res %>% mutate(coef = coef_name)
  }

  if (!is.null(LFC) && !is.null(padj)) {
    res <- res %>% filter(abs(log2FoldChange) >= !!LFC & padj < !!padj)
  }

  if (!is.null(annotation)) {
    res <- res %>% merge(annotation, by.x = "gene", by.y = "ID", all.x = TRUE)
  }

  if (!is.null(out_file)) {
    write.csv(res, file = out_file, row.names = FALSE)
    message("Saved results to ", out_file)
  }

  return(res)
}


# ─────────────────────────────────────────────────────────────────────────────
# plot_venn_diagram
# FIX-3.5  Removed dead filter: grepl("ab4phage1", ...) never matched any
#           label — venn_list was always identical to deg_list
# ─────────────────────────────────────────────────────────────────────────────
plot_venn_diagram <- function(in_dir) {

  files <- c(
    file.path(in_dir, "all_genes_1h_phage0-ab4.csv"),
    file.path(in_dir, "all_genes_1h_phage1-ab0.csv"),
    file.path(in_dir, "all_genes_4h_phage0-ab4.csv"),
    file.path(in_dir, "all_genes_4h_phage1-ab0.csv"),
    file.path(in_dir, "all_genes_24h_phage0-ab4.csv"),
    file.path(in_dir, "all_genes_24h_phage1-ab0.csv")
  )

  missing_f <- files[!file.exists(files)]
  if (length(missing_f) > 0) {
    stop("Missing files:\n", paste(missing_f, collapse = "\n"))
  }

  deg_list <- map(files, ~ {
    read_csv(.x, show_col_types = FALSE) %>%
      filter(abs(log2FoldChange) >= 1 & padj < 0.05) %>%
      pull(gene) %>%
      unique()
  })

  labels <- c(
    "1h: ATM+",  "1h: Phage+",
    "4h: ATM+",  "4h: Phage+",
    "24h: ATM+", "24h: Phage+"
  )
  names(deg_list) <- labels

  # FIX-3.5: removed venn_list <- deg_list[!grepl("ab4phage1", ...)]
  # which matched nothing; use deg_list directly
  p <- ggVennDiagram(
    deg_list,
    label          = "count",
    category.names = names(deg_list),
    set_color      = c("#4CC9FE", "#133E87",
                       "#FFA09B", "#B82132",
                       "#00FF9C", "#347928")[seq_along(deg_list)],
    label_geom     = "text",
    label_alpha    = 0.8,
    edge_size      = 0.5,
    set_size       = 3,
    label_size     = 4
  ) +
    ggplot2::scale_fill_gradient2(low = "#ffffff", high = "#D84040") +
    theme(legend.position = "none", text = element_text(size = 5))

  # Tabular membership data
  labels2 <- c("AB_1h", "Phage_1h", "AB_4h", "Phage_4h", "AB_24h", "Phage_24h")
  names(deg_list) <- labels2
  all_genes <- unique(unlist(deg_list))

  membership <- data.frame(
    item   = all_genes,
    sapply(deg_list, function(x) all_genes %in% x)
  )
  membership$region <- apply(
    membership[-1], 1,
    function(x) paste(names(x)[x], collapse = " & ")
  )

  region_counts <- as.data.frame(table(membership$region))
  colnames(region_counts) <- c("groups", "count")

  return(list(fig = p, data = deg_list, count = region_counts))
}


# ─────────────────────────────────────────────────────────────────────────────
# plot_pcoa  — unchanged
# ─────────────────────────────────────────────────────────────────────────────
plot_pcoa <- function(counts, metadata) {

  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = metadata,
    design    = ~1
  )

  keep <- rowSums(counts(dds) >= 10) >= 3
  dds  <- dds[keep, ]
  dds  <- estimateSizeFactors(dds)

  vst_counts  <- vst(dds, blind = TRUE)
  dist_matrix <- dist(t(assay(vst_counts)), method = "euclidean")
  pcoa_result <- cmdscale(dist_matrix, eig = TRUE, k = 2)

  points <- pcoa_result$points %>%
    merge(metadata, by.x = "row.names", by.y = "fname") %>%
    mutate(gtiph = factor(
      paste0(time_point, ifelse(phage == 1, "/P+", "/P-")),
      levels = c("1h/P-", "1h/P+", "4h/P-", "4h/P+", "24h/P-", "24h/P+")
    ))

  eigenvalues          <- pcoa_result$eig
  positive_eigenvalues <- eigenvalues[eigenvalues > 0]
  variance_explained   <- round(100.0 * positive_eigenvalues /
                                  sum(positive_eigenvalues), 2)

  colors <- c(
    "1h/P-"  = "#4CC9FE", "1h/P+"  = "#133E87",
    "4h/P-"  = "#FFA09B", "4h/P+"  = "#B82132",
    "24h/P-" = "#00FF9C", "24h/P+" = "#347928"
  )

  p <- ggplot(data = points, aes(x = V1, y = V2)) +
    geom_point(aes(shape = ab_conc, fill = gtiph),
               size = 4, alpha = 0.8, stroke = 0.5) +
    stat_ellipse(
      aes(group = gtiph, color = gtiph, fill = gtiph),
      type      = "norm", level = 0.95, geom = "polygon",
      alpha     = 0.1, linewidth = 0.2, show.legend = FALSE
    ) +
    scale_shape_manual(name   = "Antibiotic concentration (mg/L)",
                       values = c(21, 25, 23)) +
    scale_fill_manual(name  = "Time point/Phage presence", values = colors) +
    scale_color_manual(name = "Time point/Phage presence", values = colors) +
    labs(
      x = paste0("PCoA1 [", variance_explained[1], "%]"),
      y = paste0("PCoA2 [", variance_explained[2], "%]")
    ) +
    theme_bw() +
    theme(
      text              = element_text(size = 10),
      axis.title        = element_text(size = 10),
      legend.text       = element_text(size = 9),
      legend.title      = element_text(size = 8, face = "bold"),
      legend.position   = "bottom",
      legend.box        = "vertical",
      legend.direction  = "vertical",
      legend.box.spacing = unit(0, "pt"),
      legend.margin     = margin(t = 3, r = 0, b = 0, l = 0, unit = "pt")
    ) +
    guides(
      fill  = guide_legend(ncol = 6, title.position = "top", title.hjust = 0.5,
                           override.aes = list(shape = 21, size = 3)),
      shape = guide_legend(ncol = 3, nrow = 1, title.position = "top",
                           title.hjust = 0.5,
                           override.aes = list(fill = "gray", size = 3))
    )

  return(list(fig = p, dds_obj = dds))
}


# ─────────────────────────────────────────────────────────────────────────────
# run_permanova
# Performs PERMANOVA (adonis2) on VST-transformed Euclidean distances to test
# significance of separation by time point, phage presence, and their
# interaction. Uses the same distance matrix as plot_pcoa() for consistency.
# Also runs betadisper to test homogeneity of within-group dispersions, which
# is a required assumption check for PERMANOVA interpretation.
# Pairwise PERMANOVA is performed for all time point combinations within each
# phage group to identify which transitions drive the overall separation.
# ─────────────────────────────────────────────────────────────────────────────
run_permanova <- function(counts,
                          metadata,
                          permutations = 999,
                          seed         = 42) {
  
  library(vegan)
  library(readr)
  # ── 1. Build VST distance matrix (identical to plot_pcoa) ──────────────────
  message("Building VST distance matrix...")
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = metadata,
    design    = ~1
  )
  
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds  <- dds[keep, ]
  dds  <- estimateSizeFactors(dds)
  
  vst_counts  <- vst(dds, blind = TRUE)
  vst_mat     <- t(assay(vst_counts))           # samples × genes
  dist_matrix <- dist(vst_mat, method = "euclidean")

  # Align metadata to distance matrix sample order
  meta <- as.data.frame(metadata)
  rownames(meta) <- meta$fname
  meta <- meta[rownames(as.matrix(dist_matrix)), ]
  meta$time_point <- factor(meta$time_point, levels = c("1h", "4h", "24h"))
  meta$phage      <- factor(meta$phage,      levels = c(0, 1),
                            labels = c("Phage-", "Phage+"))
  
  # ── 2. Global PERMANOVA: time_point * phage ─────────────────────────────────
  message("Running global PERMANOVA (time_point * phage)...")
  set.seed(seed)
  perm_global <- adonis2(
    dist_matrix ~ time_point * phage,
    data         = meta,
    permutations = permutations,
    by           = "terms"          # sequential SS — order matters
  )
  
  message("\n=== Global PERMANOVA results ===")
  print(perm_global)
  
  # ── 3. Homogeneity of dispersions (betadisper) ──────────────────────────────
  # Test per time_point×phage group — violation would qualify PERMANOVA results
  message("\nTesting dispersion homogeneity (betadisper)...")
  meta$group <- interaction(meta$time_point, meta$phage, sep = "_")
  
  bd      <- betadisper(dist_matrix, group = meta$group)
  bd_test <- permutest(bd, permutations = permutations)
  
  message("\n=== Betadisper permutation test ===")
  print(bd_test)
  
  # ── 4. Pairwise PERMANOVA: phage- vs phage+ within each time point ──────────
  message("\nRunning pairwise PERMANOVA: Phage- vs Phage+ within each time point...")
  
  tps <- levels(meta$time_point)
  
  pairwise_df <- purrr::map_dfr(tps, function(tp) {
    idx       <- which(meta$time_point == tp)
    sub_meta  <- meta[idx, ]
    sub_dist  <- as.dist(as.matrix(dist_matrix)[idx, idx])
    
    set.seed(seed)
    res <- adonis2(
      sub_dist ~ phage,
      data         = sub_meta,
      permutations = permutations,
      by           = "terms"
    )
    
    tibble::tibble(
      time_point = tp,
      comparison = "Phage- vs Phage+",
      R2         = round(res$R2[1], 4),
      F_value    = round(res$F[1],  3),
      p_value    = res$`Pr(>F)`[1]
    )
  }) %>%
    dplyr::mutate(p_adj = p.adjust(p_value, method = "BH"))
  
  message("\n=== Phage- vs Phage+ within each time point ===")
  print(as.data.frame(pairwise_df))
  
  # ── 5. Tidy global results table ────────────────────────────────────────────
  global_df <- as.data.frame(perm_global) %>%
    tibble::rownames_to_column("term") %>%
    dplyr::filter(!term %in% c("Residual", "Total")) %>%
    dplyr::rename(
      df      = Df,
      SS      = SumOfSqs,
      R2      = R2,
      F_value = F,
      p_value = `Pr(>F)`
    ) %>%
    dplyr::mutate(
      R2      = round(R2,      4),
      F_value = round(F_value, 3)
    )
  
  # ── 6. Return ────────────────────────────────────────────────────────────────
  invisible(list(
    global_permanova  = global_df,
    pairwise_permanova = pairwise_df,
    betadisper_test   = bd_test,
    betadisper_object = bd,
    dist_matrix       = dist_matrix
  ))
}

# ─────────────────────────────────────────────────────────────────────────────
# get_log_fold_vector  — unchanged
# Note: reads log2FoldChange from files saved by run_lfcShrink(), which are
# already shrunken LFCs — correct for GSEA ranking.
# ─────────────────────────────────────────────────────────────────────────────
get_log_fold_vector <- function(data) {
  data <- data %>%
    filter(!is.na(keggid)) %>%
    dplyr::select(keggid, gene, baseMean, log2FoldChange, lfcSE, pvalue, padj) %>%
    arrange(dplyr::desc(abs(log2FoldChange))) %>%
    dplyr::distinct(keggid, .keep_all = TRUE)

  lfc_vector          <- data$log2FoldChange
  names(lfc_vector)   <- data$keggid
  lfc_vector          <- sort(lfc_vector, decreasing = TRUE)
  return(lfc_vector)
}


# ─────────────────────────────────────────────────────────────────────────────
# run_gsea
# FIX-3.1  eco.pathways passed in as argument — download moved to run_gsea_base
# MISC     set.seed(42) added before GSEA() for reproducibility
# ─────────────────────────────────────────────────────────────────────────────
run_gsea <- function(data, eco.pathways) {

  lfc_vector <- get_log_fold_vector(data)

  set.seed(42)  # MISC: reproducibility of permutation testing
  gsea_results <- GSEA(
    geneList      = lfc_vector,
    minGSSize     = 10,
    maxGSSize     = 500,
    pvalueCutoff  = 1,
    eps           = 0,
    pAdjustMethod = "BH",
    TERM2GENE     = eco.pathways$KEGGPATHID2EXTID
  )

  gsea_result_df <- data.frame(gsea_results@result)
  gsea_result_df <- eco.pathways$KEGGPATHID2NAME %>%
    merge(gsea_result_df, by.y = "ID", by.x = "from", all.y = TRUE) %>%
    mutate(
      pathway = gsub(" - Escher.*", "", to),
      ID      = from,
      padj    = p.adjust
    ) %>%
    dplyr::select(ID, pathway, setSize, enrichmentScore, NES,
                  pvalue, padj, qvalue, rank, leading_edge, core_enrichment)

  return(gsea_result_df)
}


# ─────────────────────────────────────────────────────────────────────────────
# run_gsea_base
# FIX-3.1  download_KEGG called ONCE here and passed into run_gsea()
# FIX-4.4  ifelse() for dir.create() replaced with if()
# ─────────────────────────────────────────────────────────────────────────────
run_gsea_base <- function(in_dir, out_dir) {

  files <- c(
    file.path(in_dir, "all_genes_1h_phage0-ab4.csv"),
    file.path(in_dir, "all_genes_1h_phage1-ab0.csv"),
    file.path(in_dir, "combo_1h_abconc4.phage1.csv"),
    file.path(in_dir, "all_genes_4h_phage0-ab4.csv"),
    file.path(in_dir, "all_genes_4h_phage1-ab0.csv"),
    file.path(in_dir, "combo_4h_abconc4.phage1.csv"),
    file.path(in_dir, "all_genes_24h_phage0-ab4.csv"),
    file.path(in_dir, "all_genes_24h_phage1-ab0.csv"),
    file.path(in_dir, "combo_24h_abconc4.phage1.csv")
  )

  missing_f <- files[!file.exists(files)]
  if (length(missing_f) > 0) {
    stop("Missing files:\n", paste(missing_f, collapse = "\n"))
  }

  deg_list <- map(files, ~ read_csv(.x, show_col_types = FALSE))

  # FIX-3.1: download once, pass to each run_gsea() call
  message("Downloading KEGG pathways (eco)...")
  eco.pathways <- download_KEGG(species = "eco")

  gsea_list <- map(deg_list, ~ run_gsea(.x, eco.pathways))

  labels <- c(
    "ab_1h",    "phage_1h",    "combo_1h",
    "ab_4h",    "phage_4h",    "combo_4h",
    "ab_24h",   "phage_24h",   "combo_24h"
  )
  names(gsea_list) <- labels

  gsea_all <- gsea_list |>
    enframe(name = "group", value = "deg_data") |>
    unnest(deg_data) |>
    separate(
      group,
      into    = c("treatment", "time_point"),
      sep     = "_",
      remove  = TRUE,
      convert = TRUE
    ) |>
    mutate(treatment = case_when(
      treatment == "phage" ~ "Phage",
      treatment == "ab"    ~ "Antibiotic",
      treatment == "combo" ~ "Antibiotic + Phage",
      TRUE                 ~ treatment
    ))

  # FIX-4.4: standard if() instead of ifelse() for side effects
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  write.csv(gsea_all, file.path(out_dir, "gsea_all_3_cat.csv"))
}


# ─────────────────────────────────────────────────────────────────────────────
# plot_gsea_heatmap
# FIX-3.2  identical() checks now enforced with stopifnot()
# ─────────────────────────────────────────────────────────────────────────────
plot_gsea_heatmap <- function(gsea_all, pathway_category) {

  gsea_all <- gsea_all %>%
    filter(padj < 0.05) %>%
    merge(pathway_category, by = "pathway", all.x = TRUE) %>%
    # FIX-NA: pathways with no match in pathway_category get category = NA.
    # ComplexHeatmap col lists reject NA names → replace with "Other".
    mutate(category = ifelse(is.na(category) | category == "", "Other", category)) %>%
    arrange(category, NES) %>%
    mutate(porder = row_number())

  all_tp <- c("1h", "4h", "24h")
  all_tr <- c("Antibiotic", "Phage", "Antibiotic + Phage")
  row_id <- "pathway"

  pathways_ordered <- gsea_all %>%
    distinct(.data[[row_id]], porder, category) %>%
    arrange(porder) %>%
    pull(.data[[row_id]])

  df_full <- gsea_all %>%
    mutate(
      time_point = factor(time_point, levels = all_tp),
      treatment  = factor(treatment,  levels = all_tr),
      le_signal  = as.numeric(str_match(leading_edge, "signal=(\\d+)%")[, 2]),
      le_tags    = as.numeric(str_match(leading_edge, "tags=(\\d+)%")[, 2]),
      le_list    = as.numeric(str_match(leading_edge, "list=(\\d+)%")[, 2])
    ) %>%
    complete(
      !!rlang::sym(row_id) := pathways_ordered,
      treatment  = factor(all_tr, levels = all_tr),
      time_point = factor(all_tp, levels = all_tp)
    ) %>%
    left_join(
      gsea_all %>% distinct(.data[[row_id]], category),
      by = setNames(row_id, row_id)
    ) %>%
    # rename safely: drop .x (possibly NA from complete), keep .y from join
    dplyr::select(-any_of("category.x")) %>%
    rename(category = `category.y`) %>%
    # FIX-NA: complete() introduces new rows whose category is NA after the
    # left_join; replace to keep cat_color names clean.
    mutate(category = ifelse(is.na(category) | category == "", "Other", category))

  make_mat <- function(trt) {
    m <- df_full %>%
      filter(treatment == trt) %>%
      dplyr::select(!!rlang::sym(row_id), time_point, NES) %>%
      pivot_wider(names_from = time_point, values_from = NES) %>%
      arrange(match(.data[[row_id]], pathways_ordered))
    mat           <- as.matrix(m[, all_tp])
    rownames(mat) <- m[[row_id]]
    mat
  }

  make_mat_signal <- function(trt) {
    m <- df_full %>%
      filter(treatment == trt) %>%
      dplyr::select(!!rlang::sym(row_id), time_point, le_signal) %>%
      pivot_wider(names_from = time_point, values_from = le_signal) %>%
      arrange(match(.data[[row_id]], pathways_ordered))
    mat           <- as.matrix(m[, all_tp])
    rownames(mat) <- m[[row_id]]
    mat
  }

  mat_ab  <- make_mat("Antibiotic");       mat_ph  <- make_mat("Phage")
  mat_int <- make_mat("Antibiotic + Phage")
  mat     <- cbind(mat_ab, mat_ph, mat_int)

  mat_ab_sig  <- make_mat_signal("Antibiotic")
  mat_ph_sig  <- make_mat_signal("Phage")
  mat_int_sig <- make_mat_signal("Antibiotic + Phage")
  mat_signal  <- cbind(mat_ab_sig, mat_ph_sig, mat_int_sig)

  # FIX-3.2: enforce alignment — silent mismatch would cause wrong colors/sizes
  stopifnot(
    "NES and signal matrices have different row order" =
      identical(rownames(mat_signal), rownames(mat)),
    "NES and signal matrices have different column order" =
      identical(colnames(mat_signal), colnames(mat))
  )

  cap     <- quantile(mat_signal, 0.95, na.rm = TRUE)
  cap     <- ifelse(is.na(cap) || cap == 0, 1, cap)
  size01  <- pmin(mat_signal, cap) / cap

  row_category <- df_full %>%
    distinct(.data[[row_id]], category) %>%
    arrange(match(.data[[row_id]], pathways_ordered)) %>%
    pull(category)

  column_split    <- factor(
    rep(c("Antibiotic", "Phage", "Antibiotic + Phage"), each = length(all_tp)),
    levels = c("Antibiotic", "Phage", "Antibiotic + Phage")
  )
  colnames(mat)   <- rep(all_tp, times = 3)

  col_fun <- colorRamp2(
    c(min(mat, na.rm = TRUE), 0, max(mat, na.rm = TRUE)),
    c("#1746A2", "#F5F5F0", "#C40C0C")
  )

  cats      <- unique(na.omit(df_full$category))   # FIX-NA: drop any residual NAs
  n_pal     <- min(12, max(3, length(cats)))
  cat_color <- setNames(
    colorRampPalette(brewer.pal(n_pal, "Set3"))(length(cats)),
    cats
  )

  ha_row <- rowAnnotation(
    Category = row_category,
    col      = list(Category = cat_color),
    show_annotation_name = FALSE,
    annotation_legend_param = list(
      Category = list(
        direction      = "horizontal",
        title_position = "topcenter",
        ncol           = 3,
        title_gp       = grid::gpar(fontsize = 8, fontface = "bold"),
        labels_gp      = grid::gpar(fontsize = 8)
      )
    )
  )

  ha_top <- HeatmapAnnotation(
    Treatment = anno_block(
      labels    = levels(column_split),
      labels_gp = gpar(fontsize = 8, fontface = "bold"),
      gp        = gpar(fill = "#F9F8F6", col = NA)
    ),
    which                = "column",
    show_annotation_name = FALSE,
    height               = grid::unit(5, "mm")
  )

  ht <- Heatmap(
    mat,
    name             = "Normalized Enrichment Score",
    col              = col_fun,
    na_col           = "transparent",
    rect_gp          = gpar(col = NA, fill = NA),
    cluster_rows     = FALSE,
    cluster_columns  = FALSE,
    row_split        = row_category,
    column_split     = column_split,
    show_row_dend    = FALSE,
    show_column_dend = FALSE,
    row_names_side   = "left",
    row_names_gp     = gpar(fontsize = 9),
    column_names_gp  = gpar(fontsize = 9),
    column_title     = NULL,
    row_title        = NULL,
    row_gap          = unit(0, "mm"),
    column_gap       = unit(1, "mm"),
    top_annotation   = ha_top,
    border_gp        = gpar(col = "grey60", lwd = 0.5),
    heatmap_legend_param = list(
      direction      = "horizontal",
      title_position = "topcenter",
      title_gp       = grid::gpar(fontsize = 8, fontface = "bold")
    ),
    layer_fun = function(j, i, x, y, w, h, fill) {
      idx <- cbind(i, j)
      v   <- mat[idx]
      rr  <- size01[idx]
      ok  <- !is.na(v)
      if (!any(ok)) return()
      m     <- grid::unit.pmin(w, h)
      min_r <- unit(0.4, "mm")
      r     <- m * 0.5 * sqrt(rr)
      r     <- grid::unit.pmax(r, min_r)
      grid::grid.circle(
        x[ok], y[ok], r = r[ok],
        gp = grid::gpar(fill = col_fun(v[ok]), col = "#bdbdbd",
                        lwd = 0.4, alpha = 0.9)
      )
    }
  )

  ht_drawn <- draw(
    ht + ha_row,
    heatmap_legend_side  = "bottom",
    annotation_legend_side = "bottom",
    merge_legends        = TRUE,
    padding              = unit(c(0, 15, 0, 0), "mm")
  )

  ht_gg <- as.ggplot(function() draw(ht_drawn)) +
    theme(plot.margin = margin(0, 0, 0, 0))

  return(ht_gg)
}


# ─────────────────────────────────────────────────────────────────────────────
# plot_pathway_dotplot
# FIX-4.7  Typo "Entrichment" → "Enrichment"
# MISC     Deprecated size= → linewidth= for line geoms
# ─────────────────────────────────────────────────────────────────────────────
plot_pathway_dotplot <- function(gsea_all, category_sort = NULL) {

  gsea_all <- gsea_all %>%
    filter(padj < 0.05) %>%
    mutate(
      label      = paste0(pathway, " (", ID, ")"),
      time_point = factor(time_point, levels = c("1h", "4h", "24h")),
      treatment  = factor(treatment,
                          levels = c("Antibiotic", "Phage", "Antibiotic + Phage"))
    )

  if (is.null(category_sort)) {
    gsea_all <- gsea_all %>% arrange(treatment, NES) %>% mutate(porder = row_number())
  } else {
    gsea_all <- gsea_all %>%
      merge(category_sort, by = "pathway", all.x = TRUE) %>%
      arrange(category, NES) %>%
      mutate(porder = row_number())

    cat_blocks <- gsea_all %>%
      distinct(category, porder, pathway) %>%
      group_by(category) %>%
      summarise(
        ymin = min(porder) - 0.5,
        ymax = max(porder) + 0.5,
        ymid = (ymin + ymax) / 2,
        .groups = "drop"
      )
    print(cat_blocks)
  }

  pkegg <- ggplot(gsea_all, aes(y = reorder(label, porder), x = time_point)) +
    geom_point(aes(size = abs(NES), fill = NES),
               color = "#bdbdbd", shape = 21, alpha = 0.8) +
    facet_grid(category ~ treatment, scales = "free") +
    scale_fill_gradient2(
      low      = "#1746A2", mid = "#F5F5F0", high = "#C40C0C",
      midpoint = 0,
      name     = "Enrichment Direction"  # FIX-4.7: was "Entrichment Direction"
    ) +
    scale_size_continuous(range = c(2, 8), name = "Enrichment Strength") +
    scale_y_discrete(position = "right") +
    labs(x = "Time point", y = "KEGG Pathway") +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(
      axis.title            = element_text(size = 10),
      axis.text.x           = element_text(size = 10),
      axis.text.y           = element_text(size = 10),
      legend.text           = element_text(size = 9),
      legend.title          = element_text(size = 9),
      strip.background      = element_rect(fill = "white"),
      legend.box            = "horizontal",
      legend.direction      = "horizontal",
      legend.title.position = "top",
      legend.position       = "top",
      panel.grid.minor.x    = element_blank(),
      panel.grid.major.x    = element_line(linetype = "dashed",
                                            linewidth = 0.5)  # MISC: size → linewidth
    )

  return(pkegg)
}


# ─────────────────────────────────────────────────────────────────────────────
# assess_interactions
# FIX-2.3  Added "Shared (additive, no interaction)" for genes significant in
#           BOTH phage and antibiotic without a significant interaction term.
#           These previously fell silently into "No significant change".
# FIX-3.4  f_combo included in file-existence check
# ─────────────────────────────────────────────────────────────────────────────
assess_interactions <- function(
    in_dir,
    out_dir     = "./export/interaction/",
    time_points = c("1h", "4h", "24h"),
    annotation
) {
  stopifnot(dir.exists(in_dir))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  if (!is.data.frame(annotation) || !"ID" %in% colnames(annotation)) {
    stop("annotation must be a data frame with column 'ID'")
  }

  # Helper: robust gene-column extraction
  standardize_gene_col <- function(df) {
    gene_candidates <- c("gene", "Gene", "ID", "Row.names", "row.names", "X", "X1")
    found           <- gene_candidates[gene_candidates %in% colnames(df)]
    if (length(found) == 0) {
      stop("Could not find a gene column. Columns: ",
           paste(colnames(df), collapse = ", "))
    }
    df %>% dplyr::rename(gene = !!found[1])
  }

  prep_res <- function(df, prefix) {
    df      <- standardize_gene_col(df)
    req     <- c("gene", "log2FoldChange", "padj")
    missing_cols <- setdiff(req, colnames(df))
    if (length(missing_cols) > 0) {
      stop("Missing columns: ", paste(missing_cols, collapse = ", "),
           " (prefix=", prefix, ")")
    }
    df %>%
      dplyr::select(gene, log2FoldChange, padj) %>%
      dplyr::group_by(gene) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::rename(
        !!paste0("lfc_",  prefix) := log2FoldChange,
        !!paste0("padj_", prefix) := padj
      )
  }

  results_list <- list()

  for (t in time_points) {
    message("Processing time point: ", t)

    f_phage <- file.path(in_dir, paste0("all_genes_", t, "_phage1-ab0.csv"))
    f_ab    <- file.path(in_dir, paste0("all_genes_", t, "_phage0-ab4.csv"))
    f_combo <- file.path(in_dir, paste0("combo_",     t, "_abconc4.phage1.csv"))
    f_int   <- file.path(in_dir, paste0("all_genes_", t, "_abconc4.phage1.csv"))

    # FIX-3.4: f_combo now included in check
    if (!all(file.exists(c(f_phage, f_ab, f_int, f_combo)))) {
      warning("Skipping ", t, ": missing one or more required files")
      next
    }

    res_phage <- read.csv(f_phage, stringsAsFactors = FALSE, check.names = FALSE)
    res_ab    <- read.csv(f_ab,    stringsAsFactors = FALSE, check.names = FALSE)
    res_combo <- read.csv(f_combo, stringsAsFactors = FALSE, check.names = FALSE)
    res_int   <- read.csv(f_int,   stringsAsFactors = FALSE, check.names = FALSE)

    phage_short <- prep_res(res_phage, "phage")
    ab_short    <- prep_res(res_ab,    "ab")
    combo_short <- prep_res(res_combo, "combo")
    int_short   <- prep_res(res_int,   "interaction")

    res_all <- phage_short %>%
      full_join(ab_short,    by = "gene") %>%
      full_join(int_short,   by = "gene") %>%
      full_join(combo_short, by = "gene")

    # Diagnostics
    message("  N genes merged: ", nrow(res_all))
    message("  Duplicated genes: ", sum(duplicated(res_all$gene)))
    message("  NA padj rates (phage/ab/int): ",
            round(mean(is.na(res_all$padj_phage)),       3), " / ",
            round(mean(is.na(res_all$padj_ab)),          3), " / ",
            round(mean(is.na(res_all$padj_interaction)), 3))

    res_all <- res_all %>%
      mutate(
        sig_phage       = !is.na(padj_phage)       & padj_phage       < 0.05,
        sig_ab          = !is.na(padj_ab)           & padj_ab          < 0.05,
        sig_interaction = !is.na(padj_interaction)  & padj_interaction < 0.05,

        # FIX-2.3: genes significant in BOTH phage and antibiotic (no interaction)
        # were previously falling through to "No significant change" — fixed by
        # adding the "Shared" category as condition 2.
        response_class = case_when(
          sig_interaction                        ~ "Interaction-defined (Non-additive)",
          sig_ab & sig_phage & !sig_interaction  ~ "Shared (additive, no interaction)",
          sig_ab   & !sig_phage                  ~ "Antibiotic-driven (no interaction)",
          sig_phage & !sig_ab                    ~ "Phage-driven (no interaction)",
          TRUE                                   ~ "No significant change"
        )
      ) %>%
      left_join(annotation, by = c("gene" = "ID"))

    message("  Annotation hit rate: ",
            round(mean(res_all$gene %in% annotation$ID), 3))
    message("  Response class counts:")
    print(table(res_all$response_class))

    outfile <- file.path(out_dir, paste0("synergy_", t, ".csv"))
    write.csv(res_all, outfile, row.names = FALSE)
    results_list[[t]] <- res_all
  }

  df_results <- results_list |>
    tibble::enframe(name = "time_point", value = "response_data") |>
    tidyr::unnest(response_data)

  summary_responses <- df_results %>%
    group_by(time_point, response_class) %>%
    summarise(N = n(), .groups = "drop")

  write.csv(df_results,        file.path(out_dir, "combo_therapy_responses.csv"),         row.names = FALSE)
  write.csv(summary_responses, file.path(out_dir, "combo_therapy_responses_summary.csv"), row.names = FALSE)

  invisible(list(data = df_results, summary = summary_responses))
}


# ─────────────────────────────────────────────────────────────────────────────
# plot_lfc_chicklets
# FIX-4.5  !!filter_col / !!cols → all_of() / any_of()  (standard tidy eval)
# FIX-4.6  height_df computed BEFORE dummy rows are appended, then dummy rows
#           built more robustly without wiping structure via dummy[] <- NA
# MISC     Deprecated size= → linewidth= for geom_hline
# ─────────────────────────────────────────────────────────────────────────────
plot_lfc_chicklets <- function(
    data,
    response_class,
    gene_category  = NULL,
    filter_col     = c("time_point", "gene", "lfc_ab", "lfc_interaction", "gene_name"),
    label_font_pt  = 9,
    label_color    = "#435663") {

  df <- data %>%
    dplyr::filter(.data$response_class == !!response_class) %>%
    dplyr::rename(gene_name = `gene.y`)

  lvls       <- c("ab", "phage", "interaction")
  # FIX-4.5: intersect to avoid missing-column errors, then use all_of()
  filter_col <- intersect(filter_col, colnames(df))
  cols       <- intersect(c("lfc_ab", "lfc_phage", "lfc_interaction"), filter_col)

  df <- df %>%
    dplyr::select(dplyr::all_of(filter_col)) %>%   # FIX-4.5
    tidyr::pivot_longer(
      cols        = dplyr::any_of(cols),            # FIX-4.5
      names_to    = "group",
      values_to   = "LFC",
      names_prefix = "lfc_"
    ) %>%
    dplyr::filter(!is.na(gene_name)) %>%
    dplyr::mutate(group = factor(.data$group, levels = lvls))

  if (!is.null(gene_category)) {
    df <- df %>%
      merge(gene_category, by = "gene_name", all.x = TRUE) %>%
      arrange(functional_group) %>%
      mutate(
        gorder = dplyr::row_number(),
        label  = paste0(
          "***", xfun::html_escape(gene_name), "***",
          "<span style='font-size:", label_font_pt,
          "pt; color:", label_color, ";'> (",
          xfun::html_escape(functional_group), ")</span>"
        )
      )
  } else {
    df <- df %>%
      mutate(
        gorder = dplyr::row_number(),
        label  = paste0("***", xfun::html_escape(gene_name), "***")
      )
  }

  # FIX-4.6: compute height_df BEFORE adding dummy rows
  height_df <- df %>%
    dplyr::group_by(time_point) %>%
    dplyr::summarise(n_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
    dplyr::mutate(height = n_genes / min(n_genes))

  # FIX-4.6: build dummy rows properly — use scale_fill_manual(drop=FALSE)
  # to keep legend entries for missing levels; avoid wiping df structure.
  # One ghost row per missing level, placed in the first time_point facet.
  first_tp      <- as.character(height_df$time_point[1])
  present_lvls  <- unique(as.character(df$group))
  missing_lvls  <- setdiff(lvls, present_lvls)
  if (length(missing_lvls) > 0) {
    template <- df[1, ]
    dummy_rows <- purrr::map_dfr(missing_lvls, function(lv) {
      r             <- template
      r$LFC         <- NA_real_
      r$group       <- factor(lv, levels = lvls)
      r$time_point  <- first_tp
      r$gorder      <- NA_integer_
      r
    })
    df <- dplyr::bind_rows(df, dummy_rows)
  }

  p <- ggplot(data = df) +
    geom_chicklet(
      aes(x = reorder(label, gorder), y = LFC, color = group),
      position  = position_dodge2(reverse = TRUE, padding = 0.3),
      linewidth = 0.3, fill = "white", width = 0.6,      # MISC: size → linewidth
      radius    = grid::unit(2, "pt"), na.rm = TRUE
    ) +
    geom_chicklet(
      aes(x = reorder(label, gorder), y = LFC, color = group, fill = group),
      alpha     = 0.6,
      position  = position_dodge2(reverse = TRUE, padding = 0.3),
      width     = 0.6, linewidth = 0.3,                  # MISC: size → linewidth
      radius    = grid::unit(2, "pt"), na.rm = TRUE
    ) +
    geom_hline(yintercept = 0, color = "grey60",
               linewidth = 0.5, linetype = "dashed") +   # MISC: size → linewidth
    facet_wrap(~time_point, scales = "free_y", ncol = 1) +
    scale_fill_manual(
      values      = c(phage = "#347928", ab = "#FF7100", interaction = "#912BBC"),
      breaks      = c("ab", "phage", "interaction"),
      drop        = FALSE,
      labels      = c(ab = "Antibiotic Only", phage = "Phage Only",
                      interaction = "Interactions"),
      na.translate = FALSE,
      name        = "Response classes"
    ) +
    scale_color_manual(
      values = c(phage = "#347928", ab = "#FF7100", interaction = "#912BBC"),
      guide  = "none"
    ) +
    labs(
      x     = "Gene (functional group)",
      y     = "Log2(Fold Change)",
      title = response_class
    ) +
    theme_minimal() +
    theme(
      axis.title.y        = element_text(size = 10),
      axis.title.x        = element_text(size = 10),
      axis.text.y         = ggtext::element_markdown(size = 11),
      axis.text.y.left    = ggtext::element_markdown(size = 11),
      axis.text.y.right   = ggtext::element_markdown(size = 11),
      plot.title          = element_text(hjust = 0.5),
      legend.text         = element_text(size = 10),
      legend.title        = element_text(size = 10),
      legend.position     = "top",
      panel.grid.minor.x  = element_blank(),
      panel.grid.major.x  = element_line(linetype = "dashed",
                                          linewidth = 0.3)  # MISC
    ) +
    facetted_pos_scales(
      x = lapply(height_df$height, function(h) scale_x_discrete())
    ) +
    force_panelsizes(rows = height_df$height) +
    coord_flip()

  return(p)
}
