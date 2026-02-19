# transcriptomics_analysis.R
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
run_lfcShrink <- function(dds, coef_name=NULL, contrast_c=NULL, LFC=NULL, padj=NULL, out_file = NULL,annotation=NULL) {
  library(DESeq2)
  library(dplyr)
  
  
  
  # Run lfcShrink
  if(!is.null(contrast_c)){
    res <- lfcShrink(dds, contrast = list(contrast_c), type = "ashr")
  }else{
    
    # Check if coefficient exists in DESeq2 object
    if(!(coef_name %in% resultsNames(dds))) {
      stop("Coefficient '", coef_name, "' not found in DESeq2 object. Available names: ", 
           paste(resultsNames(dds), collapse = ", "))
    }
    
    res <- lfcShrink(dds, coef = coef_name, type = "apeglm")
  }
  res <- res %>%
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
  
  files <- c(file.path(in_dir,"all_genes_1h_phage0-ab4.csv"),
             file.path(in_dir,"all_genes_1h_phage1-ab0.csv"),
            
             
             file.path(in_dir,"all_genes_4h_phage0-ab4.csv"),
             file.path(in_dir,"all_genes_4h_phage1-ab0.csv"),
             
             file.path(in_dir,"all_genes_24h_phage0-ab4.csv"),
             file.path(in_dir,"all_genes_24h_phage1-ab0.csv"))
  

  missing <- files[!file.exists(files)]
  if (length(missing) > 0) {
    stop("Missing files:\n", paste(missing, collapse = "\n"))
  }

  deg_list <- map(files, ~ {
    read_csv(.x, show_col_types = FALSE) %>%
      filter(abs(log2FoldChange)>1 & padj<0.05 ) %>%
      pull(gene) %>% unique()
  })
  
  labels <- c(
    "1h: ATM+",
    "1h: Phage+",
    
    "4h: ATM+",
    "4h: Phage+",
    
    "24h: ATM+",
    "24h: Phage+")
  
  
  names(deg_list) <- labels
  venn_list <- deg_list[!grepl("ab4phage1", names(deg_list))]
  
 p<- ggVennDiagram(venn_list,
                label = "count",           # shows number in each region
                category.names = names(venn_list),
                set_color = c("#4CC9FE","#133E87",
                              "#FFA09B", "#B82132",
                              "#00FF9C","#347928")[1:length(venn_list)],
                label_geom = "text",
                label_alpha = 0.8,
                edge_size = 0.5,
                set_size = 3,
                label_size = 4) +
    ggplot2::scale_fill_gradient2(low = "#ffffff",high = "#D84040") +
    theme(legend.position = "none",text=element_text(size=5)) 
 
 # get tabular data
 labels2 <- c(
   "AB_1h",
   "Phage_1h",
   
   "AB_4h",
   "Phage_4h",
   
   "AB_24h",
   "Phage_24h")
 
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
          legend.title = element_text(size=8,face="bold"),
          legend.position="bottom", legend.box="vertical",legend.direction = "vertical",
          legend.box.spacing = unit(0, "pt"),
          legend.margin = margin(t = 3, r = 0, b = 0, l = 0, unit = "pt"))+
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

run_gsea_base <- function(in_dir,out_dir){
  
  files <- c(file.path(in_dir,"all_genes_1h_phage0-ab4.csv"),
             file.path(in_dir,"all_genes_1h_phage1-ab0.csv"),
             file.path(in_dir,"combo_1h_abconc4.phage1.csv"),
             
             file.path(in_dir,"all_genes_4h_phage0-ab4.csv"),
             file.path(in_dir,"all_genes_4h_phage1-ab0.csv"),
             file.path(in_dir,"combo_4h_abconc4.phage1.csv"),
             
             file.path(in_dir,"all_genes_24h_phage0-ab4.csv"),
             file.path(in_dir,"all_genes_24h_phage1-ab0.csv"),
             file.path(in_dir,"combo_24h_abconc4.phage1.csv"))
  print(files)
  missing <- files[!file.exists(files)]
  if (length(missing) > 0) {
    stop("Missing files:\n", paste(missing, collapse = "\n"))
  }
  deg_list <- map(files, ~ {
    read_csv(.x, show_col_types = FALSE) 
  })
  
  gsea_list <- map(deg_list,~{run_gsea(.x)})
  labels <- c("ab_1h","phage_1h","combo_1h",
              "ab_4h","phage_4h","combo_4h",
              "ab_24h","phage_24h","combo_24h")
  
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
                            treatment == "combo"~ "Antibiotic + Phage", 
                            TRUE ~ treatment))
  
  ifelse(!dir.exists(file.path(out_dir)),
         dir.create(file.path(out_dir)),
         "Directory Exists")
  
  write.csv(gsea_all,paste0(out_dir,"/gsea_all_3_cat.csv"))
}

plot_gsea_heatmap <- function(gsea_all,pathway_category){
  
  gsea_all <- gsea_all %>% filter(padj<0.05) %>%
    merge(pathway_category,by="pathway",all.x=T) %>%
    arrange(category,NES) %>% 
    mutate(porder=row_number())
  
  
  all_tp <- c("1h","4h","24h")
  all_tr <- c("Antibiotic","Phage","Antibiotic + Phage")
  
  # choose your row id:
  row_id <- "pathway"   # or "label"
  
  
  pathways_ordered <- gsea_all %>%
    distinct(.data[[row_id]], porder, category) %>%
    arrange(porder) %>%
    pull(.data[[row_id]])
  
  df_full <- gsea_all %>%
    mutate(
      time_point = factor(time_point, levels = all_tp),
      treatment  = factor(treatment, levels = all_tr),
      le_signal = as.numeric(str_match(leading_edge, "signal=(\\d+)%")[,2]),
      le_tags   = as.numeric(str_match(leading_edge, "tags=(\\d+)%")[,2]),
      le_list   = as.numeric(str_match(leading_edge, "list=(\\d+)%")[,2])
    ) %>%
    # make sure every pathway has every treatment x time_point combination
    complete(
      !!rlang::sym(row_id) := pathways_ordered,
      treatment = factor(all_tr, levels = all_tr),
      time_point = factor(all_tp, levels = all_tp)
    ) %>%
    # bring category back (if complete created NAs)
    left_join(
      gsea_all %>% distinct(.data[[row_id]], category),
      by = setNames(row_id, row_id)
    )%>% rename(category=`category.y`)
  
  make_mat <- function(trt) {
    m <- df_full %>%
      filter(treatment == trt) %>%
      select(!!rlang::sym(row_id), time_point, NES) %>%
      pivot_wider(names_from = time_point, values_from = NES) %>%
      arrange(match(.data[[row_id]], pathways_ordered))
    
    mat <- as.matrix(m[, all_tp])
    rownames(mat) <- m[[row_id]]
    mat
  }
  make_mat_signal <- function(trt) {
    
    m <- df_full %>%
      filter(treatment == trt) %>%
      select(!!rlang::sym(row_id), time_point, le_signal) %>%
      tidyr::pivot_wider(
        names_from = time_point,
        values_from = le_signal
      ) %>%
      arrange(match(.data[[row_id]], pathways_ordered))
    
    mat <- as.matrix(m[, all_tp])
    rownames(mat) <- m[[row_id]]
    
    mat
  }
  
  mat_ab  <- make_mat("Antibiotic")
  mat_ph  <- make_mat("Phage")
  mat_int <- make_mat("Antibiotic + Phage")
  
  mat <- cbind(mat_ab, mat_ph, mat_int)         # 9 columns total
  
  mat_ab_sig  <- make_mat_signal("Antibiotic")
  mat_ph_sig  <- make_mat_signal("Phage")
  mat_int_sig <- make_mat_signal("Antibiotic + Phage")
  
  mat_signal <- cbind(mat_ab_sig, mat_ph_sig, mat_int_sig)
  
  
  identical(rownames(mat_signal), rownames(mat))
  identical(colnames(mat_signal), colnames(mat))
  
  # robust cap so one huge value doesn't dominate
  cap <- quantile(mat_signal, 0.95, na.rm = TRUE)
  cap <- ifelse(is.na(cap) || cap == 0, 1, cap)
  
  size01 <- pmin(mat_signal, cap) / cap
  
  
  row_category <- df_full %>%
    distinct(.data[[row_id]], category) %>%
    arrange(match(.data[[row_id]], pathways_ordered)) %>%
    pull(category)
  
  
  column_split <- factor(
    rep(c("Antibiotic","Phage","Antibiotic + Phage"), each = length(all_tp)),
    levels = c("Antibiotic","Phage","Antibiotic + Phage")
  )
  colnames(mat) <- rep(all_tp, times = 3)
  
  
  
  
  # color mapping for NES
  col_fun <- colorRamp2(
    c(min(mat, na.rm = TRUE), 0, max(mat, na.rm = TRUE)),
    c("#1746A2", "#F5F5F0", "#C40C0C")
  )
  cats <- unique(df_full$category)
  cat_color <- setNames(
    colorRampPalette(brewer.pal(12, "Set3"))(length(cats)),
    cats
  )
  
  # row annotation + split rows by category 
  ha_row <- rowAnnotation(
    Category = row_category,
    col = list(Category=cat_color),
    show_annotation_name = FALSE,
    annotation_legend_param = list(
      Category = list(
        direction = "horizontal",
        title_position = "topcenter",
        ncol=3,
        title_gp  = grid::gpar(fontsize = 8,fontface = "bold"),
        labels_gp = grid::gpar(fontsize = 8)
      )
    )
  )
  
  ha_top<- HeatmapAnnotation(
    Treatment = anno_block(
      labels = levels(column_split),
      labels_gp = gpar(fontsize = 8, fontface = "bold"),
      gp = gpar(fill = "#F9F8F6", col = NA)  # no colored blocks, text only
    ),
    which = "column",
    show_annotation_name = FALSE,
    height = grid::unit(5, "mm")
  )
  
  ht <- Heatmap(
    mat,
    name = "Normalize Enrichment Score",
    col = col_fun,
    na_col = "transparent",
    rect_gp = gpar(col = NA, fill = NA),      # no rectangles
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    row_split = row_category,                 # break by general category
    column_split = column_split,              # split into 3 treatment blocks
    
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 9),
    
    column_names_gp = gpar(fontsize = 9),
    column_title = NULL,
    row_title = NULL,
    
    row_gap = unit(0, "mm"),
    column_gap = unit(1, "mm"),
    
    top_annotation = ha_top,
    border_gp = gpar(col = "grey60", lwd = 0.5),
    heatmap_legend_param = list(
      direction = "horizontal",
      title_position = "topcenter",
      title_gp  = grid::gpar(fontsize = 8,fontface="bold")
    ),
    
    
    layer_fun = function(j, i, x, y, w, h, fill) {
      
      idx <- cbind(i, j)
      
      v  <- mat[idx]       # NES values for exactly those cells
      rr <- size01[idx]    # size values for exactly those cells
      
      ok <- !is.na(v)
      if (!any(ok)) return()
      
      m <- grid::unit.pmin(w, h)
      
      # radius: choose a max fraction of cell size
      min_r <- unit(0.4, "mm")
      r <- m * 0.5 * sqrt(rr)
      r <- grid::unit.pmax(r, min_r)
      
      grid::grid.circle(
        x[ok], y[ok],
        r = r[ok],
        gp = grid::gpar(
          fill = col_fun(v[ok]),     # signed NES -> color
          col  = "#bdbdbd",
          lwd  = 0.4,
          alpha = 0.9
        )
      )
    }
  )
  
  ht_drawn<-draw(
    ht + ha_row,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom",
    merge_legends = TRUE,          # <- important
    padding = unit(c(0, 15,0,0), "mm")
  )
  
  ht_gg <- as.ggplot(function() draw(ht_drawn))  +
    theme(plot.margin = margin(0, 0, 0, 0))
  
  return(ht_gg)
}
plot_pathway_dotplot <- function(gsea_all,category_sort=NULL){
 
  
  gsea_all <- gsea_all %>% filter(padj<0.05) %>% 
    mutate(label=paste0(pathway," (",ID,")"),
           time_point = factor(time_point, levels = c("1h", "4h", "24h")),
           treatment = factor(treatment,levels=c("Antibiotic","Phage","Antibiotic + Phage")))
  if(is.null(category_sort)){ 
   gsea_all <- gsea_all %>% arrange(treatment,NES)%>% mutate(porder=row_number())
  } else{
    
    gsea_all <- gsea_all %>% merge(category_sort,by="pathway",all.x=T) %>%
      arrange(category,NES) %>% 
      mutate(porder=row_number())

    cat_blocks <- gsea_all %>%
      distinct(category, porder, pathway) %>%
      group_by(category) %>%
      summarise(
        ymin = min(porder) - 0.5,
        ymax = max(porder) + 0.5,
        ymid = (ymin + ymax)/2,
        .groups = "drop"
      )
    print(cat_blocks)
  }
  
  pkegg <- ggplot(gsea_all, aes(y = reorder(label,porder), x = time_point)) +

    geom_point(aes(size  = abs(NES),fill=NES),color = "#bdbdbd",shape =21, alpha =0.8) +
    facet_grid(category~treatment,scales = "free")+
    scale_fill_gradient2(low = "#1746A2",mid = "#F5F5F0",high = "#C40C0C",
                         midpoint = 0,name = "Entrichment Direction") +
    scale_size_continuous(range = c(2, 8),name = "Enrichment Strength") +
    scale_y_discrete(position = "right")+
    labs(x="Time point",y="KEGG Pathway")+
    coord_cartesian(clip = "off") +    
    theme_minimal()+
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

assess_interactions <- function(
    in_dir,
    out_dir = "./export/interaction/",
    time_points = c("1h", "4h", "24h"),
    annotation
) {
  stopifnot(dir.exists(in_dir))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  if (!is.data.frame(annotation) || !"ID" %in% colnames(annotation)) {
    stop("annotation must be a data frame with column 'ID'")
  }
  
  # Helper: robust gene column extraction
  standardize_gene_col <- function(df) {
    # Common DESeq2 csv patterns: gene, Gene, rownames stored as X or Row.names
    gene_candidates <- c("gene", "Gene", "ID", "Row.names", "row.names", "X", "X1")
    found <- gene_candidates[gene_candidates %in% colnames(df)]
    if (length(found) == 0) {
      stop("Could not find a gene column in input. Columns: ", paste(colnames(df), collapse=", "))
    }
    gene_col <- found[1]
    df <- df %>% dplyr::rename(gene = !!gene_col)
    df
  }
  
  prep_res <- function(df, prefix) {
    df <- standardize_gene_col(df)
    
    req <- c("gene", "log2FoldChange", "padj")
    missing <- setdiff(req, colnames(df))
    if (length(missing) > 0) {
      stop("Missing required columns: ", paste(missing, collapse=", "),
           " in file with prefix ", prefix)
    }
    
    df %>%
      dplyr::select(gene, log2FoldChange, padj) %>%
      # enforce 1 row per gene (protect against prior merges/duplicates)
      dplyr::group_by(gene) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::rename(
        !!paste0("lfc_", prefix) := log2FoldChange,
        !!paste0("padj_", prefix) := padj
      )
  }
  
  results_list <- list()
  
  for (t in time_points) {
    message("Processing time point: ", t)
    
    f_phage <- file.path(in_dir, paste0("all_genes_", t, "_phage1-ab0.csv"))
    f_ab <- file.path(in_dir, paste0("all_genes_", t, "_phage0-ab4.csv"))
    f_combo <- file.path(in_dir, paste0("combo_", t, "_abconc4.phage1.csv"))
    f_int <- file.path(in_dir, paste0("all_genes_", t, "_abconc4.phage1.csv"))
    
    if (!all(file.exists(c(f_phage, f_ab, f_int)))) {
      warning("Skipping ", t, ": missing one or more files")
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
      full_join(ab_short, by = "gene") %>%
      full_join(int_short, by = "gene") %>%
      full_join(combo_short, by = "gene")
    
    # Diagnostics
    message("  N genes merged: ", nrow(res_all))
    message("  Duplicated genes: ", sum(duplicated(res_all$gene)))
    message("  NA padj rates (phage/ab/int): ",
            round(mean(is.na(res_all$padj_phage)), 3), " / ",
            round(mean(is.na(res_all$padj_ab)), 3), " / ",
            round(mean(is.na(res_all$padj_interaction)), 3))
    
    res_all <- res_all %>%
      mutate(
        sig_phage = !is.na(padj_phage) & padj_phage < 0.05,
        sig_ab    = !is.na(padj_ab)    & padj_ab    < 0.05,
        sig_interaction = !is.na(padj_interaction) & padj_interaction < 0.05,
        
        response_class = case_when(
          sig_interaction  ~ "Interaction-defined (Non-additive)",
          sig_ab & !sig_phage & !sig_interaction ~ "Antibiotic-driven (no interaction)",
          sig_phage & !sig_ab & !sig_interaction ~ "Phage-driven (no interaction)",
          TRUE                                   ~ "No significant change"
        )
      ) %>%
      left_join(annotation, by = c("gene" = "ID"))
    
    message("  Annotation hit rate: ", round(mean(res_all$gene %in% annotation$ID), 3))
    
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
  
  write.csv(df_results, file.path(out_dir, "combo_therapy_responses.csv"), row.names = FALSE)
  write.csv(summary_responses, file.path(out_dir, "combo_therapy_responses_summary.csv"), row.names = FALSE)
  
  invisible(list(data = df_results, summary = summary_responses))
}


plot_lfc_chicklets <- function(
    data,
    response_class,
    gene_category = NULL,
    filter_col = c("time_point","gene","lfc_ab","lfc_interaction","gene_name"),
    label_font_pt = 9,
    label_color = "#435663") {
  
  # --- your code, same flow ---
  df <- data %>%
    dplyr::filter(response_class == !!response_class) %>%
    dplyr::rename(gene_name = `gene.y`)
  
  # keep your select, but make it non-erroring if one column is missing
  # (this is the one change that prevents "full of errors" situations)
  
  lvls <- c("ab","phage","interaction")
  cols <-  c()
  filter_col <- intersect(filter_col, colnames(df))
  if("lfc_ab" %in% filter_col){
    cols <-  c(cols,"lfc_ab")
  }
  if("lfc_phage" %in% filter_col){
    
    cols <-  c(cols,"lfc_phage")
  }
  if("lfc_interaction" %in% filter_col){
    cols <-  c(cols,"lfc_interaction")
  }
  
  
  
  df <- df %>%
    dplyr::select(!!filter_col) %>%
    tidyr::pivot_longer(
      cols = !!cols,
      names_to = "group",
      values_to = "LFC",
      names_prefix = "lfc_"
    ) %>%
    dplyr::filter(!is.na(gene_name)) %>%
    dplyr::mutate(group = factor(group, levels = !!lvls))
  
  
  
  
  if (!is.null(gene_category)) {
    df <- df %>%
      merge(gene_category, by = "gene_name", all.x = TRUE) %>%
      arrange(functional_group) %>%
      mutate(
        gorder = dplyr::row_number(),
        label = paste0(
          "***", xfun::html_escape(gene_name), "***",
          "<span style='font-size:", label_font_pt, "pt; color:", label_color, ";'> (",
          xfun::html_escape(functional_group),
          ")</span>"
        )
      )
  } else {
    df <- df %>%
      mutate(
        gorder = dplyr::row_number(),
        label = paste0("***", xfun::html_escape(gene_name), "***")
      )
  }
  
  height_df <- df %>%
    dplyr::group_by(time_point) %>%
    dplyr::summarise(n_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
    dplyr::mutate(height = n_genes / min(n_genes))
  
  # dirty solution for legend fix
  dummy <- df[rep(1, length(lvls)), ]   # copy structure
  dummy[] <- NA                                   # set everything NA
  dummy$group <- factor(lvls, levels = lvls)
  dummy$time_point <- "4h"
  df <- rbind(df, dummy)
  
  p <- ggplot(data = df) +
    geom_chicklet(aes(x = reorder(label, gorder), y = LFC, color = group),
                  position = position_dodge2(reverse = TRUE, padding = 0.3),
                  size = .3, fill = "white", width = 0.6, radius = grid::unit(2, "pt"),na.rm = TRUE
    ) +
    geom_chicklet(aes(x = reorder(label, gorder), y = LFC, color = group, fill = group),
                  alpha = 0.6,
                  position = position_dodge2(reverse = TRUE, padding = 0.3),
                  width = 0.6, size = .3, radius = grid::unit(2, "pt"),na.rm = TRUE
    ) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.5, linetype = "dashed") +
    facet_wrap(~time_point, scales = "free_y", ncol = 1) +
    scale_fill_manual(
      values = c(phage = "#347928", ab = "#FF7100", interaction = "#912BBC"),
      breaks = c("ab", "phage", "interaction"),   # <- forces legend order + presence
      drop = FALSE,                               # <- keep missing levels in legend
      labels = c(ab = "Antibiotic Only", phage = "Phage Only", interaction = "Interactions"),
      na.translate = FALSE,
      name = "Response classes"
    ) +
    scale_color_manual(values = c("phage" = "#347928", "ab" = "#FF7100", "interaction" = "#912BBC"),
                       guide = "none") +
    labs(
      x = "Gene (functional group)",
      y = "Log2(Fold Change)",
      title = response_class
    ) +
    theme_minimal() +
    theme(
      plot.text = element_text(size = 11),
      axis.title = element_text(size = 11),
      
      # because coord_flip()
      axis.title.y = element_text(size=10),
      axis.title.x = element_text(size=10),
      axis.text.y = ggtext::element_markdown(size = 11),
      axis.text.y.left = ggtext::element_markdown(size = 11),
      axis.text.y.right = ggtext::element_markdown(size = 11),
      
      plot.title = element_text(hjust = 0.5),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.position = "top",
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_line(linetype = "dashed", size = .3)
    ) +
    
    facetted_pos_scales(x = lapply(height_df$height, function(h) scale_x_discrete())) +
    force_panelsizes(rows = height_df$height) +
    coord_flip()
  return(p)
}

