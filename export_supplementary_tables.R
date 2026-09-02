# Regenerates the 8 manuscript supplementary tables from the current repo
# and export/ state. Run after transcriptomics_pipeline.Rmd, time_kill_analysis.Rmd,
# and imc_heat_flow.Rmd have been (re)run, so their outputs are current.
#
# Requires ../time_kill_data/time_kill_imputed.csv (sibling directory, not
# part of this repo).

suppressMessages({library(tidyverse); library(dplyr); library(readr)})

# Adjust to your manuscript folder.
out_dir <- "/Users/mohammad.razavi/Library/CloudStorage/OneDrive-KarolinskaInstitutet/pas/supplementary_files"

# ── File 1: time-kill ─────────────────────────────────────────────────────
tk <- read.csv("../time_kill_data/time_kill_imputed.csv", stringsAsFactors = FALSE, check.names = FALSE)
write_tsv(tk, file.path(out_dir, "suppl.file_1_time_kill.tsv"))
cat("File 1 (time_kill): ", nrow(tk), "rows,", ncol(tk), "cols\n")

# ── File 2: IMC digitalization ───────────────────────────────────────────
imc <- read_tsv("output/heat_flow_metadata.tsv", show_col_types = FALSE) %>%
  select(well, conditions = condition, time_h, heat_flow_relative, antibiotic_concentration = ab_conc)
write_tsv(imc, file.path(out_dir, "suppl.file_2_IMC_digitalization.tsv"))
cat("File 2 (IMC): ", nrow(imc), "rows,", ncol(imc), "cols\n")

# ── File 3: PCoA VST euclidean distance matrix ───────────────────────────
pcoa <- read_tsv("export/pcoa_vst_mat_euclidean.tsv", show_col_types = FALSE)
write_tsv(pcoa, file.path(out_dir, "suppl.file_3_pcoa_vst_mat_euclidean.tsv"))
cat("File 3 (pcoa): ", nrow(pcoa), "rows,", ncol(pcoa), "cols\n")

# ── File 4: DESeq2 outputs (combo_therapy_responses) ─────────────────────
deseq <- read_csv("export/interaction/combo_therapy_responses.csv", show_col_types = FALSE)
write_tsv(deseq, file.path(out_dir, "suppl.file_4_deseq2_outputs.tsv"))
cat("File 4 (deseq2): ", nrow(deseq), "rows,", ncol(deseq), "cols\n")

# ── File 5: KEGG GSEA ─────────────────────────────────────────────────────
gsea <- read_csv("export/kegg_gsea/gsea_all_3_cat.csv", show_col_types = FALSE) %>%
  select(-1)  # drop leading unnamed row-index column from write.csv
write_tsv(gsea, file.path(out_dir, "suppl.file_5_kegg_gsea.tsv"))
cat("File 5 (kegg_gsea): ", nrow(gsea), "rows,", ncol(gsea), "cols\n")

# ── File 6: metadata (excludes 8 mg/L RNA-seq samples -- not part of this
#    study's transcriptomics design; keeps A38 since it was sequenced and
#    only excluded from the DESeq2 model as a QC decision, not from the
#    study design) ────────────────────────────────────────────────────────
meta <- read_tsv("export/metadata.tsv", show_col_types = FALSE) %>%
  filter(ab_conc != 8)
write_tsv(meta, file.path(out_dir, "suppl.file_6_metadata.tsv"))
cat("File 6 (metadata): ", nrow(meta), "rows,", ncol(meta), "cols\n")

# ── File 7: raw counts (filtered: post-8mg/L, post-A38 exclusion) ───────
rc <- read_tsv("export/raw_counts.tsv", show_col_types = FALSE)
write_tsv(rc, file.path(out_dir, "suppl.file_7_raw_counts.tsv"))
cat("File 7 (raw_counts): ", nrow(rc), "rows,", ncol(rc), "cols\n")

# ── File 8: EC45622 annotation (raw, unmodified -- includes known dup IDs;
#    the annotation-deduplication fix lives only in the derived DEG table,
#    not this raw reference table) ───────────────────────────────────────
ann <- read_csv("data/22ET500456__with-ecolik12.csv", show_col_types = FALSE)
write_tsv(ann, file.path(out_dir, "suppl.file_8_ec45622_annotation.tsv"))
cat("File 8 (annotation): ", nrow(ann), "rows,", ncol(ann), "cols\n")

message("\nDone. All 8 supplementary files written to ", out_dir)
