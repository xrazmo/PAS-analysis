# Transcriptomics Functions - Quick Reference Guide

## Overview

This guide provides quick examples for using the custom transcriptomics analysis functions.

## Installation and Setup

```r
# Required packages
install.packages(c("tidyverse", "readxl"))
BiocManager::install(c("DESeq2", "sva"))

# Load libraries
library(tidyverse)
library(readxl)
library(DESeq2)
library(sva)

# Source functions
source("transcriptomics_functions.R")
```

## Function Reference

### 1. `load_transcriptomics_data()`

**Purpose**: Load counts, metadata, and annotation files

**Basic Usage**:
```r
data <- load_transcriptomics_data(
  counts_path = "path/to/counts.csv",
  metadata_path = "path/to/metadata.xlsx"
)

counts <- data$counts
metadata <- data$metadata
annotation <- data$annotation  # NULL if not provided
```

**Key Parameters**:
- `counts_path` - Path to count matrix CSV (genes = rows, samples = columns)
- `metadata_path` - Path to metadata (Excel or CSV)
- `annotation_path` - (Optional) Path to gene annotations
- `metadata_filter_col` - Column to filter unwanted samples
- `metadata_filter_pattern` - Pattern to exclude (default: "missing")
- `fname_col` - Column name matching count column names (default: "fname")

**Advanced Example**:
```r
data <- load_transcriptomics_data(
  counts_path = "./counts.csv",
  metadata_path = "./metadata.xlsx",
  annotation_path = "./annotations.csv",
  metadata_filter_col = "qc_status",
  metadata_filter_pattern = "failed|missing",
  fname_col = "sample_name"
)
```

---

### 2. `create_deseq2_dataset()`

**Purpose**: Create DESeq2 object with optional SVA batch correction

**Basic Usage**:
```r
results <- create_deseq2_dataset(
  counts = counts,
  metadata = metadata,
  design_formula = "~ condition + time"
)

dds <- results$dds
vsd <- results$vsd
svobj <- results$sva
```

**Key Parameters**:
- `counts` - Count matrix from `load_transcriptomics_data()`
- `metadata` - Metadata from `load_transcriptomics_data()`
- `design_formula` - Design formula (string or formula object)
- `min_count` - Minimum total count per gene (default: 10)
- `use_sva` - Enable SVA batch correction (default: TRUE)
- `sva_formula` - Custom formula for SVA (default: uses design_formula)
- `n_sv` - Number of surrogate variables (NULL = auto-detect)
- `blind` - Blind dispersion for VST (default: TRUE)

**Without SVA**:
```r
results <- create_deseq2_dataset(
  counts = counts,
  metadata = metadata,
  design_formula = "~ treatment",
  use_sva = FALSE
)
```

**Custom SVA Formula**:
```r
results <- create_deseq2_dataset(
  counts = counts,
  metadata = metadata,
  design_formula = "~ batch + treatment + time",
  use_sva = TRUE,
  sva_formula = "~ treatment + time",  # Exclude batch from SVA model
  n_sv = 2  # Force 2 surrogate variables
)
```

---

### 3. `run_deseq2_pipeline()`

**Purpose**: One-step pipeline combining data loading and DESeq2 analysis

**Basic Usage**:
```r
results <- run_deseq2_pipeline(
  counts_path = "./counts.csv",
  metadata_path = "./metadata.xlsx",
  design_formula = "~ condition"
)

# Access everything
counts <- results$counts
metadata <- results$metadata
dds <- results$dds
vsd <- results$vsd
```

**Full Example**:
```r
results <- run_deseq2_pipeline(
  counts_path = './output/count_data_all.csv',
  metadata_path = "data/metadata.xlsx",
  annotation_path = "./output/annotations.csv",
  design_formula = "~ time_point + antibiotic + phage + ab_conc",
  min_count = 10,
  use_sva = TRUE,
  metadata_filter_col = "comment",
  metadata_filter_pattern = "missing",
  fname_col = "fname"
)
```

---

### 4. `summarize_deseq2()`

**Purpose**: Print quick summary of DESeq2 object

**Usage**:
```r
summarize_deseq2(dds)
```

**Output**:
```
=== DESeq2 Dataset Summary ===
Samples: 24
Genes: 4523
Design: ~SV1 + time_point + phage + ab_conc
...
```

---

### 5. `get_normalized_counts()`

**Purpose**: Extract normalized count matrices

**Usage**:
```r
# DESeq2 normalized counts
norm_counts <- get_normalized_counts(dds, method = "deseq")

# VST-transformed counts
vst_counts <- get_normalized_counts(dds, method = "vst", vsd = vsd)
```

---

## Common Workflows

### Workflow 1: Quick Analysis (One Function)

```r
# Load libraries
library(tidyverse)
library(DESeq2)
library(sva)
source("transcriptomics_functions.R")

# Run everything
results <- run_deseq2_pipeline(
  counts_path = "counts.csv",
  metadata_path = "metadata.xlsx",
  design_formula = "~ condition"
)

# Differential expression
res <- results(results$dds, contrast = c("condition", "treated", "control"))
summary(res)
```

### Workflow 2: Step-by-Step with Custom Options

```r
# Load libraries
library(tidyverse)
library(DESeq2)
library(sva)
source("transcriptomics_functions.R")

# Step 1: Load data
data <- load_transcriptomics_data(
  counts_path = "counts.csv",
  metadata_path = "metadata.xlsx"
)

# Step 2: Create DESeq2 object
deseq_res <- create_deseq2_dataset(
  counts = data$counts,
  metadata = data$metadata,
  design_formula = "~ batch + condition",
  min_count = 5,
  use_sva = TRUE
)

# Step 3: Differential expression
dds <- deseq_res$dds
res <- results(dds, contrast = c("condition", "treated", "control"))

# Step 4: Visualization
plotMA(res)
plotPCA(deseq_res$vsd, intgroup = "condition")
```

### Workflow 3: Complex Design with Multiple Factors

```r
source("transcriptomics_functions.R")

results <- run_deseq2_pipeline(
  counts_path = "counts.csv",
  metadata_path = "metadata.xlsx",
  design_formula = "~ batch + genotype + treatment + genotype:treatment",
  use_sva = TRUE,
  sva_formula = "~ genotype + treatment",  # Don't include batch in SVA
  min_count = 10
)

dds <- results$dds

# Different contrasts
res1 <- results(dds, contrast = c("treatment", "drug", "vehicle"))
res2 <- results(dds, contrast = c("genotype", "mutant", "wildtype"))
res3 <- results(dds, name = "genotypemutant.treatmentdrug")  # Interaction
```

---

## Tips and Best Practices

### Design Formula Tips

```r
# Simple design
"~ condition"

# Multiple factors (additive)
"~ time + treatment"

# With interaction term
"~ genotype + drug + genotype:drug"

# Nested design
"~ patient + treatment"

# Multiple covariates
"~ batch + sex + age + treatment"
```

### SVA Tips

1. **When to use SVA**: Use when you have unwanted variation (batch effects, unknown confounders)

2. **SVA formula**: Include only biological variables of interest, not batch variables:
   ```r
   design_formula = "~ batch + treatment"
   sva_formula = "~ treatment"  # Exclude batch
   ```

3. **Number of SVs**: Let it auto-detect (`n_sv = NULL`) unless you have good reason

### File Format Requirements

**Counts file (CSV)**:
```
Geneid,Sample1,Sample2,Sample3
Gene1,100,150,120
Gene2,50,60,55
```

**Metadata file (Excel/CSV)**:
```
sample_id,fname,condition,time_point
S1,Sample1,control,0h
S2,Sample2,treated,0h
```

---

## Troubleshooting

### Error: "Metadata and counts don't match!"

**Solution**: Ensure `fname` column in metadata matches count column names exactly

```r
# Check column names
colnames(counts)

# Check metadata fname column
metadata$fname

# They should match perfectly
```

### Error: "Some count samples not found in metadata"

**Solution**: Check for typos or missing samples in metadata

```r
# Find missing samples
missing <- setdiff(colnames(counts), metadata$fname)
print(missing)
```

### Low gene retention after filtering

```r
# Reduce min_count threshold
results <- create_deseq2_dataset(
  counts = counts,
  metadata = metadata,
  design_formula = "~ condition",
  min_count = 5  # Instead of 10
)
```

### SVA taking too long

```r
# Force fewer surrogate variables
results <- create_deseq2_dataset(
  counts = counts,
  metadata = metadata,
  design_formula = "~ condition",
  use_sva = TRUE,
  n_sv = 2  # Limit to 2 SVs
)

# Or disable SVA
results <- create_deseq2_dataset(
  counts = counts,
  metadata = metadata,
  design_formula = "~ condition",
  use_sva = FALSE
)
```

---

## Examples from Your Original Code

### Your Original Workflow

```r
# Original code translated
results <- run_deseq2_pipeline(
  counts_path = './output/count_data_all.csv',
  metadata_path = "data/metadata_ecphage.xlsx",
  annotation_path = "./output/22ET500456__with-ecolik12.csv",
  design_formula = "~ time_point + antibiotic + phage + ab_conc",
  min_count = 10,
  use_sva = TRUE,
  sva_formula = "~ time_point + phage + ab_conc",  # Your SVA model
  metadata_filter_col = "comment",
  metadata_filter_pattern = "missing",
  fname_col = "fname"
)

dds <- results$dds
vsd <- results$vsd
```

---

## Additional Resources

- DESeq2 vignette: `vignette("DESeq2")`
- SVA vignette: `vignette("sva")`
- Help: `?function_name`

---

**Created**: 2026-02-05
**Author**: Custom Transcriptomics Analysis Functions
