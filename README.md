# Cell–Cell Interactions of T Cells in MASLD Across Tissue Compartments

This repository contains code from my research project studying **Building and characterising the T-cell communication network in multiple tissue compartments of metabolic dysfunction
associated steatotic liver disease (MASLD)** 
The project uses **single-cell RNA sequencing (scRNA-seq)** data and uses **NicheNet** and **MultiNicheNet** to investigate ligand–receptor signaling between T cells and other cell types (endothelial cells, hepatocytes) under *Fibrosis* and *No_fibrosis* conditions.  

⚠️ **Note:** This repository does **not** contain raw or processed scRNA-seq data.  
It only provides analysis scripts and representative plots for documentation and reproducibility.
---

## 🔬 Research Context

- **Goal:** Characterize how T cells communicate with endothelial cells and hepatocytes in different tissue compartments.  
- **Methods:**
  - Preprocessing and clustering using **Seurat**  
  - Cell–cell communication inference with **NicheNet** / **MultiNicheNet**  
  - Visualisation via UMAPs, circos plots, bubble plots, and regulatory networks and abundance plots  
- **Comparisons:** Fibrosis vs. No_fibrosis conditions across tissue compartments (Liver, SAT, VAT).  

---

## 📂 Repository Layout

This repository currently includes key analysis scripts and supporting files:

├── **T_cell_and_Endothelial_analysis.R**  
├── **T_cell_and_Hepatocyte_analysis.R**   
├── **Packages_requirement**  
└── **Curated Ligand-receptor database for HUMAN**

⚠️ Note: Input data (e.g. Seurat objects) and output plots are **not** included in this repository.  

---

## 📂 Recommended Project Structure (Example)

When running this analysis on a HPC or local environment, it is recommended that the following directory structure is used for clarity and reproducibility:
MASLD_research_project
├── **R_scripts # R scripts for analysis**
├── **Seurat_objects # Input Seurat objects (.rds, not included here)**
├── **T_cell_and_Endothelial_plots # Plots for T cell–endothelial analysis**
└── **T_cell_and_Hepatocyte_plots # Plots for T cell–hepatocyte analysis**



## 🖥️ Computational Environment and Resources

All bioinformatics analyses were conducted using **RStudio** and **R version 2024.09.0-375** with **R 4.4.1** in a **Rocky 9** environment, deployed on the **Apocrita High Performance Computing (HPC) cluster** provided by Queen Mary University of London. Interactive R sessions were launched via **OnDemand** with the following parameters:

- 1 CPU core  
- 128 GB RAM  
- 24-hour runtime per job  

Since most analyses in R and the **MultiNicheNet** pipeline were not explicitly parallelised, a single CPU core with high-memory allocation (128 GB) was chosen to ensure efficient resource usage during interactive sessions on the Apocrita HPC cluster.

The computational environment was further configured with the following modules:

- **GSL 2.7.1**  
- **HDF5 1.10.2**  
- **GDAL 2.3.1**  
- **PROJ 5.2.0**  
- **GEOS 3.7.1**  

These modules were essential for dependency resolution and performance optimisation of key R packages used throughout the pipeline.

## 📦 Following packages were installed:

- `Seurat`, `SeuratObject`  
- `SingleCellExperiment`, `SummarizedExperiment`  
- `dplyr`, `ggplot2`, `tidyverse`  
- `nichenetr`, `multinichenetr`  
- `igraph`, `ggraph`  
- `RColorBrewer`  

Install packages via:

```R
install.packages(c("Seurat", "dplyr", "ggplot2", "tidyverse", "igraph", "ggraph", "RColorBrewer"))
BiocManager::install(c("SingleCellExperiment", "SummarizedExperiment"))
remotes::install_github("saeyslab/nichenetr")
remotes::install_github("saeyslab/multinichenetr")


📖 Citations

If you use this repository, please cite:

NicheNet
Browaeys R, Saelens W, Saeys Y. NicheNet: modeling intercellular communication by linking ligands to target genes. Nat Methods. 2020.

MultiNicheNet (Preprint)
Browaeys R, Gilis J, Sang-aram C, De Bleser P, Hoste L, Tavernier S, Lambrechts D, Seurinck R, Saeys Y.
MultiNicheNet: a flexible framework for differential cell–cell communication analysis from multi-sample, multi-condition single-cell transcriptomics data. bioRxiv. 2023. DOI: 10.1101/2023.06.13.544751

MultiNicheNet Vignettes (GitHub)
For workflows and examples, see the vignette files in the vignettes/ directory of the MultiNicheNet GitHub repository
.

Resource ID (RRID)
MultiNicheNet (RRID: SCR_025903), registered in the dkNET SciCrunch registry.

This repository
[LaibaNaz]. Cell–cell interactions of T cells in MASLD across tissue compartments. GitHub repository: [MSc-Research-project-multinichenetr-analysis-T-cell-communication-in-MASLD](https://github.com/Laibafaisal26/MSc-Research-project-multinichenetr-analysis-T-cell-communication-in-MASLD)
