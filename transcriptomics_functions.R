# transcriptomics_analysis.R
library(DESeq2)
library(sva)
library(dplyr)
library(ggVennDiagram)
library(ggpubr)
# Function to load counts, metadata, and annotation
load_transcriptomics_data <- function(counts_path,
                                      metadata_path,
                                      annotation_path = NULL,
                                      metadata_filter_col = "comment",
                                      metadata_filter_pattern = "missing",
                                      sample_id_col = "sample_id",
                                      fname_col = "fname") {
  
  message("Loading count data from: ", counts_path)
  counts <- read.table(counts_path, header = TRUE, row.names = 'Geneid', sep = ',') %>%
    data.matrix()
  
  message("  -> Loaded ", nrow(counts), " genes across ", ncol(counts), " samples")
  
  # Load metadata
  message("Loading metadata from: ", metadata_path)
  if (grepl("\\.xlsx?$", metadata_path, ignore.case = TRUE)) {
    metadata <- readxl::read_xlsx(metadata_path)
  } else {
    metadata <- read.csv(metadata_path)
  }
  
  # Filter metadata
  if (!is.null(metadata_filter_col) && metadata_filter_col %in% colnames(metadata)) {
    metadata <- metadata %>%
      filter(!grepl(metadata_filter_pattern, .data[[metadata_filter_col]])) %>%
      dplyr::select(-all_of(metadata_filter_col))
  }
  
  # Arrange by sample ID
  if (sample_id_col %in% colnames(metadata)) {
    metadata <- metadata %>% arrange(.data[[sample_id_col]])
  }
  
  # Match metadata to counts
  if (fname_col %in% colnames(metadata)) {
    if (!all(colnames(counts) %in% metadata[[fname_col]])) {
      missing_samples <- setdiff(colnames(counts), metadata[[fname_col]])
      warning("Some count samples not found in metadata: ", paste(missing_samples, collapse = ", "))
    }
    metadata <- metadata[match(colnames(counts), metadata[[fname_col]]), ]
    stopifnot("Metadata and counts don't match!" = all(colnames(counts) == metadata[[fname_col]]))
  }
  
  rownames(metadata) <- metadata[[fname_col]]
  
  metadata <- metadata %>% mutate(
    time_point = factor(time_point, levels = c("1h", "4h", "24h")),
    phage      = factor(phage, levels = c(0,1)),
    ab_conc    = factor(ab_conc, levels = c(0,4,8))
  )
  
  # Load annotation (optional)
  annotation <- NULL
  if (!is.null(annotation_path)) {
    message("Loading annotation from: ", annotation_path)
    annotation <- read.csv(annotation_path)
    message("  -> Loaded annotations for ", nrow(annotation), " genes")
  }
  
  message("\n✓ Data loading complete!")
  return(list(counts = counts, metadata = metadata, annotation = annotation))
}

# Function to run DESeq2 with SVA per time point
run_DESeq2_sva <- function(counts, metadata, 
                           formula_str = "~ ab_conc + phage", 
                           time_point = NULL, n_svs = NULL) {
  
  if (!is.null(time_point)) {
    message("\n→ Subsetting for time point: ", time_point)
    metadata <- metadata %>% filter(as.character(.data[["time_point"]]) == !!time_point)
    counts <- counts[, metadata$fname]
  }
  
  # Design matrices for SVA
  mod <- model.matrix(as.formula(formula_str), data = metadata)
  mod0 <- model.matrix(~1, data = metadata)
  
  # Estimate surrogate variables
  message("Estimating surrogate variables...")
  
  svobj <- sva(counts, mod, mod0)
  # Check if any surrogate variables were detected
  if (!is.null(svobj$sv) && ncol(svobj$sv) > 0) {
    sv_df <- as.data.frame(svobj$sv)
    colnames(sv_df) <- paste0("SV", 1:ncol(sv_df))
    metadata <- cbind(metadata, sv_df)
    message("Detected ", ncol(svobj$sv), " surrogate variable(s). Added to metadata.")
  } else {
    message("No surrogate variables detected. Proceeding without SVs.")
  }
  
  # Build DESeq2 design including SVs
  sv_formula <- paste(colnames(sv_df), collapse = " + ")
  design_formula <- as.formula(paste("~", sv_formula, "+ ab_conc * phage"))
  message("DESeq2 design: ", deparse(design_formula))
  
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = design_formula)
  dds <- DESeq(dds)
  
  # Return
  return(list(dds = dds, svobj = svobj, metadata = metadata))
}

# Function to run lfcShrink and save results
run_lfcShrink <- function(dds, coef_name, LFC=NULL, padj=NULL, out_file = NULL,annotation=NULL) {
  library(DESeq2)
  library(dplyr)
  
  # Check if coefficient exists in DESeq2 object
  if(!(coef_name %in% resultsNames(dds))) {
    stop("Coefficient '", coef_name, "' not found in DESeq2 object. Available names: ", 
         paste(resultsNames(dds), collapse = ", "))
  }
  
  # Run lfcShrink
  res <- lfcShrink(dds, coef = coef_name, type = "apeglm") %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    mutate(coef = coef_name)  # store the contrast name
  
  if(!is.null(LFC) & !is.null(padj)){
    res <- res %>% filter(abs(log2FoldChange)	> !!LFC & padj < !!padj)
  }
  
  if(!is.null(annotation)){
    res <- res %>% merge(annotation, by.x="gene",by.y ="ID",all.x=T) 
  }
 
  
  # Optionally save to CSV
  if(!is.null(out_file)) {
    write.csv(res, file = out_file, row.names = FALSE)
    message("Saved results to ", out_file)
  }
  
  return(res)
}

plot_venn_diagram <- function(in_dir){
  
  files <- list.files(
    path = in_dir,
    pattern = "venn_degs*",
    full.names = TRUE
  )
  
  deg_list <- map(files, ~ {
    read_csv(.x, show_col_types = FALSE) %>%
      pull(gene) %>% unique()
  })
  
  labels <- c(
    "1h: ATM+",
    "1h: Phage+",
    
    "24h: ATM+",
    "24h: Phage+",
    
    "4h: ATM+",
    "4h: Phage+"
  )
  
  
  names(deg_list) <- labels
  venn_list <- deg_list[!grepl("ab4phage1", names(deg_list))]
  
 p<- ggVennDiagram(venn_list,
                label = "count",           # shows number in each region
                category.names = names(venn_list),
                set_color = c("#4CC9FE","#133E87",
                              "#00FF9C","#347928",
                              "#FFA09B", "#B82132")[1:length(venn_list)],
                label_geom = "text",
                label_alpha = 0.8,
                edge_size = 0.5) +
    ggplot2::scale_fill_gradient2(low = "#ffffff",high = "#D84040") +
    theme(legend.position = "none",text=element_text(size=10)) 
 
 # get tabular data
 labels2 <- c(
   "AB_1h",
   "Phage_1h",
   
   "AB_24h",
   "Phage_24h",
   
   "AB_4h",
   "Phage_4h"
 )
 names(deg_list) <- labels2
 all_genes <- unique(unlist(deg_list))
 
 membership <- data.frame(
   item = all_genes,
   sapply(deg_list, function(x) all_genes %in% x)
 )
 
 membership$region <- apply(
   membership[-1],
   1,
   function(x) paste(names(x)[x], collapse = " & ")
 )
 
 region_counts <- as.data.frame(table(membership$region))
 colnames(region_counts) <- c("groups", "count")
 
 return(list(fig=p, data=deg_list, count=region_counts))
}

plot_pcoa <- function(counts,metadata){
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = metadata,
    design    = ~ 1
  )
  
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  
  dds <- estimateSizeFactors(dds)
  
  #Option 3
  vst_counts <- vst(dds, blind = TRUE)  # Apply VST
  dist_matrix <- dist(t(assay(vst_counts)), method = "euclidean")
  
  pcoa_result <- cmdscale(dist_matrix, eig = TRUE, k = 2)  # k is the number of dimensions
  
  points <- pcoa_result$points %>% merge(metadata,by.x='row.names',by.y='fname') %>% 
    mutate(gtiph=factor(paste0(time_point,ifelse(phage==1,'/P+','/P-')),
                        levels=c('1h/P-','1h/P+','4h/P-','4h/P+','24h/P-','24h/P+')))
  
  eigenvalues <- pcoa_result$eig
  positive_eigenvalues <- eigenvalues[eigenvalues > 0]
  variance_explained <- round(100.0*positive_eigenvalues / sum(positive_eigenvalues),2)
  colors = c('1h/P-'="#4CC9FE",'1h/P+'="#133E87",
             '4h/P-'="#FFA09B",'4h/P+'="#B82132",
             '24h/P-'="#00FF9C",'24h/P+'="#347928")
  
  p <- ggplot(data=points,aes(x=V1, y=V2)) +
    geom_point(aes(shape=ab_conc,fill=gtiph),size=4,alpha=.8,stroke=0.5)+
    stat_ellipse( aes(group=gtiph, color=gtiph,fill=gtiph),
        
      type = "norm",          # "norm" = multivariate normal (most common)
      level = 0.95,           # 95% confidence ellipse
      geom = "polygon",       # filled ellipse
      alpha = 0.1,           # transparency of fill
      linewidth = 0.2,        # border thickness
      show.legend = FALSE     # hide duplicate legend entry
    ) +
    scale_shape_manual(name="Antibiotic concentration (mg/L)",values=c(21,25,23)) +
    scale_fill_manual(name="Time point/Phage presence", values=colors)+
    scale_color_manual(name="Time point/Phage presence", values=colors)+
    labs(y=paste0("PCoA2 [",variance_explained[2],"%]"),x=paste0("PCoA1 [",variance_explained[1],"%]"))+
    
    theme_bw()+
    theme(text=element_text(size=10),
          axis.title=element_text(size=10),
          legend.text=element_text(size=9),
          legend.title = element_text(size=9),
          legend.position="top", legend.box="vertical",legend.direction = "vertical",
          legend.box.spacing = unit(0, "pt"),
          legend.margin = margin(t = -10, r = 0, b = 0, l = 0, unit = "pt"))+
    guides(fill = guide_legend(ncol=6,title.position = "top",
                               title.hjust=0.5, override.aes = list(shape = 21,size=3)),
           shape = guide_legend(ncol=3,nrow=1,title.position = "top", 
                                title.hjust=0.5, override.aes = list(fill="gray",size=3)))
  
  return(list(fig=p, dds_obj=dds))
}

get_log_fold_vector<-function(data){
  data<- data  %>% filter(!is.na(keggid)) %>%
    dplyr::select(keggid,gene,baseMean,log2FoldChange,lfcSE,pvalue,padj) %>% 
    arrange(dplyr::desc(abs(log2FoldChange))) %>%   
    dplyr::distinct(keggid, .keep_all = TRUE)
  
  
  lfc_vector<- data$log2FoldChange 
  names(lfc_vector) <- data$keggid
  lfc_vector <- sort(lfc_vector, decreasing = TRUE)
  return(lfc_vector)
}

run_gsea<- function(data){
  
  lfc_vector <- get_log_fold_vector(data)
  
  eco.pathways <- download_KEGG(species="eco")
  
  gsea_results <- GSEA(
    geneList = lfc_vector, # Ordered ranked gene list
    minGSSize = 10, # Minimum gene set size
    maxGSSize = 500, # Maximum gene set set
    pvalueCutoff = 1, # p-value cutoff
    eps = 0, # Boundary for calculating the p value
    pAdjustMethod = "BH", # Benjamini-Hochberg correction
    TERM2GENE = eco.pathways$KEGGPATHID2EXTID
  )
  
  gsea_result_df <- data.frame(gsea_results@result) 
  
  gsea_result_df <- eco.pathways$KEGGPATHID2NAME %>%
    merge(gsea_result_df,by.y='ID',by.x='from',all.y=T) %>%
    mutate(pathway=gsub(" - Escher.*", "", to),ID=from,padj=p.adjust) %>%
    dplyr::select(ID,pathway,setSize,enrichmentScore,NES,pvalue,padj,qvalue,rank,leading_edge,core_enrichment)
  
  return(gsea_result_df)
}

run_gsea_venn <- function(in_dir){
  files <- list.files(
    path = in_dir,
    pattern = "all_genes*",
    full.names = TRUE
  )
  
  deg_list <- map(files, ~ {
    read_csv(.x, show_col_types = FALSE) 
  })
  
  gsea_list <- map(deg_list,~{run_gsea(.x)})
  labels <- c("ab_1h","phage_1h","ab_24h","phage_24h","ab_4h","phage_4h")
  names(gsea_list) <- labels
  
  gsea_all<- gsea_list |>
    enframe(name = "group", value = "deg_data") |>
    unnest(deg_data) |>
    separate(
      group,
      into = c("treatment", "time_point"),
      sep = "_",
      remove = TRUE,           # drop original group column
      convert = TRUE           # try to convert time_point to numeric if possible
    ) |>
    mutate(
      treatment = case_when(treatment == "phage" ~ "Phage", 
                            treatment == "ab"~ "Antibiotic", 
                            TRUE ~ treatment))
  
  write.csv(gsea_all,paste0(in_dir,"/gsea_all_venn.csv"))
}

plot_pathway_dotplot <- function(gsea_all){
  
  gsea_all <- gsea_all %>% filter(padj<0.05) %>% 
    mutate(label=paste0(pathway," (",ID,")"),
           time_point = factor(time_point, levels = c("1h", "4h", "24h"))) %>% 
    arrange(treatment,NES)%>% mutate(porder=row_number())
  
  pkegg <- ggplot(gsea_all, aes(y = reorder(label,porder), x = time_point)) +
    geom_point(aes(size  = abs(NES),fill=NES),color = "#bdbdbd",shape =21, alpha =0.8) +
    facet_grid(.~treatment)+
    scale_fill_gradient2(low = "#1746A2",mid = "#F5F5F0",high = "#C40C0C",
                         midpoint = 0,name = "Entrichment Direction") +
    scale_size_continuous(range = c(2, 8),name = "Enrichment Strength") +
    scale_y_discrete(position = "right")+
    labs(x="Time point",y="KEGG Pathway")+
    theme_bw()+
    theme(
          axis.title=element_text(size=10),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.text=element_text(size=9),
          legend.title = element_text(size=9),
          strip.background =element_rect(fill="white"),
          legend.box="horizontal",legend.direction = "horizontal",
          legend.title.position = "top",
          legend.position = 'top',
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_line(linetype='dashed',size=.5))
  return(pkegg)
}

###--------Replaced by KEGG GSEA
get_venn_region_deg_pathway <- function(
    venn_data,
    deg_1h,
    deg_4h,
    annotation,
    output_file = NULL
) {
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(KEGGREST)    # for keggGet
  library(clusterProfiler)  # for bitr_kegg
  
  message("→ Starting Venn region + KEGG pathway annotation...")
  
  # ────────────────────────────────────────────────
  # 1. Compute Venn regions (exclusive & intersection groups)
  # ────────────────────────────────────────────────
  message("  Computing Venn regions...")
  
  # Helper to compute exclusive sets
  exclusive <- function(target, others) {
    setdiff(target, Reduce(union, others))
  }
  
  # Define your regions generically (you can easily add/modify these)
  venn_regions <- list(
    early_only = exclusive(
      venn_data$Phage_1h,
      list(venn_data$Phage_4h, venn_data$AB_1h, venn_data$AB_4h, venn_data$AB_24h, venn_data$Phage_24h)
    ),
    sustained_1h_4h = exclusive(
      intersect(venn_data$Phage_1h, venn_data$Phage_4h),
      list(venn_data$AB_1h, venn_data$AB_4h, venn_data$AB_24h, venn_data$Phage_24h)
    ),
    late_emerging = exclusive(
      venn_data$Phage_4h,
      list(venn_data$Phage_1h, venn_data$AB_1h, venn_data$AB_4h, venn_data$AB_24h, venn_data$Phage_24h)
    )
    # Add more regions here if needed, e.g.:
    # phage_specific = exclusive(union(venn_data$Phage_1h, venn_data$Phage_4h), ...)
  )
  
  # Convert to tidy long format
  df_venn_genes <- venn_regions |>
    enframe(name = "venn_group", value = "gene") |>
    unnest(gene) |>
    left_join(annotation, by = c("gene" = "ID")) |>   # assuming ID is the gene column in annotation
    arrange(venn_group, gene)
  
  message("  → Found ", nrow(df_venn_genes), " gene-region assignments")
  
  # ────────────────────────────────────────────────
  # 2. Match DEGs (logFC, padj, etc.) to regions
  # ────────────────────────────────────────────────
  message("  Matching DE values to regions...")
  
  # Helper to filter DEGs for a region
  get_degs_for_group <- function(group_name, time_df) {
    genes_in_group <- df_venn_genes |>
      filter(venn_group == group_name) |>
      pull(gene)
    time_df |> filter(gene %in% genes_in_group)
  }
  
  lfc_matched <- list(
    early_only       = get_degs_for_group("early_only", deg_1h),
    sustained_1h     = get_degs_for_group("sustained_1h_4h", deg_1h),
    sustained_4h     = get_degs_for_group("sustained_1h_4h", deg_4h),
    late_emerging    = get_degs_for_group("late_emerging", deg_4h)
  )
  
  df_degs_venn <- lfc_matched |>
    enframe(name = "venn_group", value = "deg_data") |>
    unnest(deg_data)

  # ────────────────────────────────────────────────
  # 3. KEGG mapping & pathway info
  # ────────────────────────────────────────────────
  message("  Performing KEGG mapping and fetching pathway info...")
  
  unique_kegg_ids <- unique(df_degs_venn$keggid)
  
  kegg_map <- bitr_kegg(
    geneID   = unique_kegg_ids,
    fromType = "kegg",
    toType   = "Path",
    organism = "eco"
  )
  
  unique_paths <- unique(kegg_map$Path)
  
  message("  → Fetching info for ", length(unique_paths), " unique KEGG pathways...")
  
  # Fetch pathway info safely
  path_info_list <- map(
    unique_paths,
    ~ tryCatch(
      keggGet(.x)[[1]],
      error = function(e) {
        message("Warning: Failed to fetch ", .x, " – using NAs")
        list(NAME = NA, CLASS = NA, PATHWAY_MAP = NA)
      }
    )
  )
  
  path_info_df <- tibble(
    Path             = unique_paths,
    path_name        = map_chr(path_info_list, ~ paste(.x$NAME,        collapse = " | ") %||% NA_character_),
    path_class       = map_chr(path_info_list, ~ paste(.x$CLASS,       collapse = " | ") %||% NA_character_),
    path_description = map_chr(path_info_list, ~ paste(.x$PATHWAY_MAP, collapse = " | ") %||% NA_character_)
  )
  
  # Join pathway info to kegg_map
  kegg_map <- kegg_map |> left_join(path_info_df, by = "Path")
  
  # Final merge with DEG data
  df_final <- df_degs_venn |>
    left_join(kegg_map, by = c("keggid" = "kegg")) |>
    arrange(venn_group, gene)
  
  # ────────────────────────────────────────────────
  # 4. Save & return
  # ────────────────────────────────────────────────
  if (!is.null(output_file)) {
    message("  → Writing result to: ", output_file)
    write.csv(df_final, output_file, row.names = FALSE)
  } else {
    message("  → No output file specified; returning object only")
  }
  
  message("Done. Final table has ", nrow(df_final), " rows.")
  invisible(df_final)
}

get_pathways <- function(unique_kegg_ids){
  
  message("  Performing KEGG mapping and fetching pathway info...")
  
  kegg_map <- bitr_kegg(
    geneID   = unique_kegg_ids,
    fromType = "kegg",
    toType   = "Path",
    organism = "eco"
  )
  
  unique_paths <- unique(kegg_map$Path)
  
  message("  → Fetching info for ", length(unique_paths), " unique KEGG pathways...")
  
  # Fetch pathway info safely
  path_info_list <- map(
    unique_paths,
    ~ tryCatch(
      keggGet(.x)[[1]],
      error = function(e) {
        message("Warning: Failed to fetch ", .x, " – using NAs")
        list(NAME = NA, CLASS = NA, PATHWAY_MAP = NA)
      }
    )
  )
  
  path_info_df <- tibble(
    Path             = unique_paths,
    path_name        = map_chr(path_info_list, ~ paste(.x$NAME,        collapse = " | ") %||% NA_character_),
    path_class       = map_chr(path_info_list, ~ paste(.x$CLASS,       collapse = " | ") %||% NA_character_),
    path_description = map_chr(path_info_list, ~ paste(.x$PATHWAY_MAP, collapse = " | ") %||% NA_character_)
  )
  
  # Join pathway info to kegg_map
  kegg_map <- kegg_map |> left_join(path_info_df, by = "Path")
  
 return(kegg_map)
}
