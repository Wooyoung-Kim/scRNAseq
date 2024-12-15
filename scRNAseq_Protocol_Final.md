
# Single-Cell RNA-seq Protocol

## Introduction
Single-cell RNA sequencing (scRNA-seq) is a powerful technique for analyzing the gene expression of individual cells. This protocol guides you through the major steps of scRNA-seq data analysis, including preprocessing, clustering, and annotation, using **R** and the **Seurat** package.

---

## Prerequisites

### Software Requirements
- R version ≥4.0.0
- RStudio
- Python (version ≥3.7)
- Cell Ranger (version ≥4.0)

### R Packages
Install the following packages:

```r
# Install Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# Install required packages
BiocManager::install(c("Seurat", "ggplot2", "dplyr", "Matrix", "SoupX", "DoubletFinder", "DropletUtils"))
```

### Python Libraries
Install the following Python libraries:

```bash
pip install pandas requests
```

---

## Workflow

### 1. Public Data Download
Use the following Python script to download publicly available datasets:

```python
import os
import requests

# Function to download files
def download_file(url, output_path):
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        with open(output_path, "wb") as file:
            for chunk in response.iter_content(1024):
                file.write(chunk)
        print(f"Downloaded: {output_path}")
    else:
        print(f"Failed to download: {url}")

# Example URLs and output paths
urls = {
    "Sample1": "https://example.com/sample1.tar.gz",
    "Sample2": "https://example.com/sample2.tar.gz"
}

output_dir = "data/"
os.makedirs(output_dir, exist_ok=True)

for sample, url in urls.items():
    output_path = os.path.join(output_dir, f"{sample}.tar.gz")
    download_file(url, output_path)
```

### 2. Cell Ranger Preprocessing
Prepare and process raw sequencing data using Cell Ranger:

#### Generate a Reference Genome
Use the following Python script to edit and prepare a custom reference genome if necessary:

```python
import os

def edit_gtf(input_gtf, output_gtf):
    with open(input_gtf, 'r') as infile, open(output_gtf, 'w') as outfile:
        for line in infile:
            if not line.startswith("#"):
                columns = line.strip().split("	")
                if "gene" in columns[2]:
                    outfile.write(line)
    print(f"Edited GTF saved to: {output_gtf}")

edit_gtf("path/to/input.gtf", "path/to/output.gtf")
```

#### Run Cell Ranger Count
Execute the `cellranger count` command to process raw sequencing data:

```bash
cellranger count   --id=sample1   --transcriptome=/path/to/reference   --fastqs=/path/to/fastqs   --sample=sample1   --localcores=8   --localmem=64
```

### 3. Load Required Libraries
```r
library(Seurat)
library(ggplot2)
library(dplyr)
library(SoupX)
library(DoubletFinder)
library(DropletUtils)
```

### 4. Import and Clean Data

#### Load raw and filtered matrices
```r
raw_matrix <- Read10X_h5("path/to/raw_feature_bc_matrix.h5")
filtered_matrix <- Read10X_h5("path/to/filtered_feature_bc_matrix.h5")

# Create SoupX object
soup_channel <- SoupChannel(raw_matrix, filtered_matrix)
soup_channel <- autoEstCont(soup_channel, forceAccept = TRUE)
adjusted_counts <- adjustCounts(soup_channel)

# Create Seurat object from cleaned data
seurat_object <- CreateSeuratObject(counts = adjusted_counts, min.cells = 3, min.features = 200)
```

#### Filter low-quality cells and genes
```r
# Add mitochondrial percentage as metadata
seurat_object["percent.mt"] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")

# QC filtering
seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```

### 5. Normalize Data
```r
seurat_object <- NormalizeData(
  seurat_object,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
```

### 6. Identify Highly Variable Features
```r
seurat_object <- FindVariableFeatures(
  seurat_object,
  selection.method = "vst",
  nfeatures = 2000
)

# Plot variable features
VariableFeaturePlot(seurat_object)
```

### 7. Scale Data
```r
seurat_object <- ScaleData(
  seurat_object,
  vars.to.regress = "percent.mt"
)
```

### 8. Perform PCA
```r
seurat_object <- RunPCA(seurat_object, npcs = 30)

# Visualize PCA
ElbowPlot(seurat_object)
```

### 9. Doublet Detection with DoubletFinder
```r
# Run DoubletFinder analysis
doublet_finder_analysis <- function(seurat_obj) {
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.1)

  sweep_res <- paramSweep(seurat_obj)
  stats <- summarizeSweep(sweep_res)
  pK <- find.pK(stats)
  optimal_pk <- as.numeric(as.character(pK$BCmetric[which.max(pK$BCmetric)]))

  nExp <- round(0.05 * nrow(seurat_obj@meta.data))
  seurat_obj <- doubletFinder_v3(seurat_obj, PCs = 1:10, pK = optimal_pk, nExp = nExp)

  seurat_obj <- subset(seurat_obj, subset = DF.classifications_0.25 == "Singlet")
  return(seurat_obj)
}

seurat_object <- doublet_finder_analysis(seurat_object)
```

### 10. Clustering and Dimensionality Reduction
```r
# Clustering
seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)

# UMAP
seurat_object <- RunUMAP(seurat_object, dims = 1:10)
DimPlot(seurat_object, reduction = "umap", label = TRUE)
```

### 11. Identify Marker Genes
```r
markers <- FindAllMarkers(
  seurat_object,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# View top markers
head(markers)
```

### 12. Annotate Cell Types
```r
# Example annotation
cell_type_annotation <- c(
  "0" = "Basal",
  "1" = "Ciliated",
  "2" = "Goblet",
  "3" = "Club"
)

seurat_object$cell_type <- cell_type_annotation[seurat_object$seurat_clusters]

# Visualize annotated cell types
DimPlot(seurat_object, group.by = "cell_type", label = TRUE)
```

### 13. Save Processed Data
```r
# Save Seurat object
saveRDS(seurat_object, file = "path/to/seurat_processed.RDS")
```

---

## Outputs
- **Quality Control Plots**: Histograms of feature counts, percentage of mitochondrial content, etc.
- **Dimensionality Reduction Plots**: PCA, UMAP, or t-SNE visualizations.
- **Cluster-Specific Marker Genes**: Table of differentially expressed genes for each cluster.
- **Cell Type Annotations**: UMAP with labeled cell types.

---

## References
- **Seurat Documentation**: [https://satijalab.org/seurat/](https://satijalab.org/seurat/)
- **Cell Ranger Documentation**: [https://support.10xgenomics.com/](https://support.10xgenomics.com/)
- **Single-Cell RNA-seq Analysis Workflow**: [https://www.bioconductor.org/](https://www.bioconductor.org/)

---

This protocol provides a general guide for scRNA-seq data analysis. Adjust parameters and methods based on your dataset and research objectives.
