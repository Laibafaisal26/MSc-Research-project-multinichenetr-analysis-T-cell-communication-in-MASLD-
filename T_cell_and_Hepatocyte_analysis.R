# This file contains and well commented R script for running a multinichenetr analysis on T cells and Hepatocytes.

library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(nichenetr)
library(multinichenetr)
library(nichenetr) # Please update to v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(dplyr)

Seurat_Object <- readRDS("/data/scratch/bt24974/SeuObjx.rds")
DefaultAssay(Seurat_Object) <- "RNA"


#### cell type annotation again because c6 - endothelial were named as blood vessels ###


cluster_to_celltype <- c(
  "C0" = "CD8_Effector_Memory",
  "C1" = "CD4_Naïve_T_cells",
  "C2" = "NK_cells",
  "C3" = "CD14_Monocytes",
  "C4" = "CD8_Cytotoxic Tcells",
  "C5" = "Macrophage",
  "C6" = "Endothelial_cells",
  "C7" = "CD4_T_cells",
  "C8" = "Naïve_B_cells",
  "C9" = "ACKR1_Endothelial_cells",
  "C10" = "Myeloid_Dendritic_cells",
  "C11" = "CD16_Monocytes",
  "C12" = "CD56_dim_NK cells",
  "C13" = "Fibroblasts",
  "C14" = "Plasmacytoid_Dendritic_cells",
  "C15" = "Mast_cells",
  "C16" = "TAGLN_Endothelial_cells",
  "C17" = "Erythroid_cells",
  "C18" = "Liver_Ductal_cells",
  "C19" = "Plasma_Cells",
  "C20" = "Hepatocytes",
  "C21" = "CD4_Proliferating_T_cells",
  "C22" = "CD4_Dendritic_cells",
  "C23" = "TREM2_Dendritic_cells",
  "C24" = "Mesothelial_cells",
  "C25" = "Adipocytes",
  "C26" = "Adipocytes",
  "C27" = "B_cells",
  "C28" = "Platelets",
  "C29" = "Plasma cells",
  "C30" = "Plasma Cells"
)
Seurat_Object@meta.data[["cell_type"]] <- cluster_to_celltype[Seurat_Object@meta.data$cluster]


grouping <- c(
  # T cells cluster
  "CD8_Effector_Memory" = "T_cells",
  "CD8_Cytotoxic_Tcells" = "T_cells",
  "CD4_Naïve_T_cells" = "T_cells",
  "CD4_T_cells" = "T_cells",
  "CD4_Proliferating_T cells" = "T_cells",
  "NK_cells" = "T_cells",
  
  # Monocyte and macrophages cluster
  "CD14_Monocytes" = "Monocytes",
  "CD16_Monocytes" = "Monocytes",
  "Macrophage" = "Macrophages",
  
  # Dendritic cells cluster
  "Myeloid_Dendritic_cells" = "Dendritic_Cells",
  "Plasmacytoid_Dendritic cells" = "Dendritic_Cells",
  "CD4_Dendritic_cells" = "Dendritic_Cells",
  "TREM2_Dendritic_cells" = "Dendritic_Cells",
  
  # Hepatocytes
  "Hepatocytes" = "Parenchymal",
  
  # B cells
  "Naïve_B_cells" = "B_cells",
  "B_cells" = "B_cells",
  "Plasma_Cells" = "B_cells",
  
  "Endothelial_cells" = "Endothelial_cells", 
  "ACKR1_Endothelial_cells" = "Endothelial_cells",
  "TAGLN_Endothelial_cells" = "Endothelial_cells"
)


Seurat_Object@meta.data$celltype_grouping <- grouping[Seurat_Object@meta.data$cell_type]

Liver_T_cells_and_hepatocyte <- subset(Seurat_Object, subset = Tissue == "LIVER")
unique(Seurat_Object$cluster)



## subsetting to remove NK cells "C2"
T_and_Hepatocyte_object <- subset(Liver_T_cells_and_hepatocyte, subset = cluster %in% c("C0", "C1", "C20", "C21", "C7", "C4"))

## EXploratory Data Analysis ##

DimPlot(T_and_Hepatocyte_object, reduction = "pca", group.by = "cell_type")

VlnPlot(T_and_Hepatocyte_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "Stage")

#### MULTINICHENET analysis ###


organism = "human"
options(timeout = 120)

if(organism == "human"){
  
  lr_network_all = 
    readRDS(url(
      "https://zenodo.org/record/10229222/files/lr_network_human_allInfo_30112033.rds"
    )) %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  ligand_target_matrix = readRDS(url(
    "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"
  ))
  
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
  
}



sce = Seurat::as.SingleCellExperiment(T_and_Hepatocyte_object, assay = "RNA")
#sce = alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()
sce$cell_type <- gsub(" ", "_", sce$cell_type)
sample_id = "Sample" # identifies each sample and which cell came from which sample 
group_id = "Stage"   # for the biological conditions of interest, Fibrosis Vs No_fibrosis
celltype_id = "cell_type" # cell type annotations cluster for each cell type 


covariates = NA 
batches = NA 

#colData(sce)$Stage <- as.character(colData(sce)$Stage) # converts the Stage column to a Vector


contrasts_oi = c("'Fibrosis-No_fibrosis','No_fibrosis-Fibrosis'") # setting the contrast for conditions of interest

contrast_tbl = tibble(
  contrast = c("Fibrosis-No_fibrosis", "No_fibrosis-Fibrosis"),
  group = c("Fibrosis", "No_fibrosis")
) # creates a data frame that maps contrast to its corresponding group 


#celltypes = SummarizedExperiment::colData(sce)[, celltype_id] # extracts the cell type information from SCE

SummarizedExperiment::colData(sce)$Sample = SummarizedExperiment::colData(sce)$Sample %>% make.names()
#SummarizedExperiment::colData(sce)$Tissue = SummarizedExperiment::colData(sce)$Tissue %>% make.names()
SummarizedExperiment::colData(sce)$cell_type = SummarizedExperiment::colData(sce)$cell_type %>% make.names()
SummarizedExperiment::colData(sce)$Stage = SummarizedExperiment::colData(sce)$Stage %>% make.names()

senders_oi = c("Hepatocytes") # setting sender cell types and this line identifies which cells belong in that category 

receivers_oi = c(
  "CD8_Cytotoxic_Tcells", 
  "CD8_Effector_Memory",
  "CD4_Naïve_T_cells",
  "CD4_T_cells",
  "CD4_Proliferating_T_cells"
) # setting receiver cell types 

conditions_keep = c("Fibrosis","No_fibrosis") # Specifies conditions to keep for the analysis 


# filters the SCE object to only keep codnitions that were specified, Fibrosis Vs No_fibrosis 
sce = sce[, SummarizedExperiment::colData(sce)[,group_id] %in% 
            conditions_keep
] 

# Sets minimum number of cells per group per sample per condition 
# DEFAULT 10 but i chnaged it to 3 becuse one of the samples has less than 10
#Important for statistical power and validity—too few cells can lead to unreliable results.
min_cells = 3

## To make the names in sample_id syntactically correct. applied to the Sample Ids.

#SummarizedExperiment::colData(sce)[, sample_id] <- make.names(SummarizedExperiment::colData(sce)[, sample_id])

#Calculates the per sample cell abundance information for each cell type 
#Only includes cell that are previously specified in the "senders" and "receivers"
#Ensures that each sample has enough cells, specified by the "min_cells" ti be considered valid for downstream analysis 
abundance_info = get_abundance_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  batches = batches
)

#Plotting the abundance cell type information
abundance_info$abund_plot_sample

### removing absent cells ###

#Groups by condition and cell type, and counts in how many samples the cell type is present.
#Identifies cell types that are completely absent in at least one condition
abundance_df_summarized = abundance_info$abundance_data %>% 
  mutate(keep = as.logical(keep)) %>% 
  group_by(group_id, celltype_id) %>% 
  summarise(samples_present = sum((keep)))

celltypes_absent_one_condition = abundance_df_summarized %>% 
  filter(samples_present == 0) %>% pull(celltype_id) %>% unique() 

#Identifies cell types that are present in ≥2 samples in a condition
celltypes_present_one_condition = abundance_df_summarized %>% 
  filter(samples_present >= 2) %>% pull(celltype_id) %>% unique()

#Cell types only present in one condition, absent in the other. These are condition-specific.
condition_specific_celltypes = intersect(
  celltypes_absent_one_condition, 
  celltypes_present_one_condition)


#Determines how many unique conditions (e.g., "Fibrosis", "No_fibrosis") are in the dataset
total_nr_conditions = SummarizedExperiment::colData(sce)[,group_id] %>% 
  unique() %>% length()

#Finds cell types with too few cells (<2 samples) in every condition. These are unreliable for analysis.
absent_celltypes = abundance_df_summarized %>% 
  filter(samples_present < 2) %>% 
  group_by(celltype_id) %>% 
  summarise(n = n()) %>% 
  filter(n == total_nr_conditions) %>% 
  pull(celltype_id)

print("condition-specific celltypes:")
print(condition_specific_celltypes)
## character(0)

print("absent celltypes:")
## [1] "absent celltypes:"
print(absent_celltypes)

#Controls whether to include condition-specific cell types in the analysis.
# originally FALSE but chnaged it to TRUE as i want to see which cell types are condtion specific 
analyse_condition_specific_celltypes = FALSE

if(analyse_condition_specific_celltypes == TRUE){
  senders_oi = senders_oi %>% setdiff(absent_celltypes)
  receivers_oi = receivers_oi %>% setdiff(absent_celltypes)
} else {
  senders_oi = senders_oi %>% 
    setdiff(union(absent_celltypes, condition_specific_celltypes))
  receivers_oi = receivers_oi %>% 
    setdiff(union(absent_celltypes, condition_specific_celltypes))
}

sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)
]

min_sample_prop = 0.50
fraction_cutoff = 0.03 # default 0.05 but i chnaged as in some samples the number of cells is lower 

#Computes which genes are expressed in each cell type and sample based on frequency thresholds.
frq_list = get_frac_exprs(
  sce = sce, 
  sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id, 
  batches = batches, 
  min_cells = min_cells, 
  fraction_cutoff = fraction_cutoff, min_sample_prop = min_sample_prop)

#Filters the expression matrix to keep only the expressed genes.
genes_oi = frq_list$expressed_df %>% 
  filter(expressed == TRUE) %>% pull(gene) %>% unique() 
sce = sce[genes_oi, ]

#Combines abundance and expression data for downstream analysis.
#Uses ligand-receptor network info (lr_network) to model interactions.
abundance_expression_info = process_abundance_expression_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  lr_network = lr_network, 
  batches = batches, 
  frq_list = frq_list, 
  abundance_info = abundance_info)
abundance_expression_info$celltype_info$pb_df %>% head()
abundance_expression_info$celltype_info$pb_df_group %>% head()
abundance_expression_info$sender_receiver_info$pb_df %>% head()
abundance_expression_info$sender_receiver_info$pb_df_group %>% head()

#Computes the Differential expression genes for each cell type using the contrasts set 
DE_info = get_DE_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  batches = batches, covariates = covariates, 
  contrasts_oi = contrasts_oi, 
  min_cells = min_cells, 
  expressed_df = frq_list$expressed_df)

#Shows P values and DE results 
DE_info$celltype_de$de_output_tidy %>% head()
DE_info$hist_pvals


# computes Empirical p values but in this setting it is disabled beacuse P values calculated 
# in the above step are uniform mostly 
empirical_pval = FALSE

if(empirical_pval == TRUE){
  DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
  celltype_de = DE_info_emp$de_output_tidy_emp %>% select(-p_val, -p_adj) %>% 
    rename(p_val = p_emp, p_adj = p_adj_emp)
} else {
  celltype_de = DE_info$celltype_de$de_output_tidy
} 

# Matches DE genes with ligand-receptor pairs to create communication links.
sender_receiver_de = combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)

sender_receiver_de %>% head(20)

#Thresholds to define significantly DE genes.
logFC_threshold = 0.50
p_val_threshold = 0.05

p_val_adj = FALSE  # default in vignette is FALSE 

# scores gene set based on differential expression across different conditions 
geneset_assessment = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de, logFC_threshold, p_val_adj, p_val_threshold
  ) %>% 
  bind_rows() 
geneset_assessment

#geneset_assessment_adjustedPval = contrast_tbl$contrast %>% 
#  lapply(
#    process_geneset_data, 
#    celltype_de, logFC_threshold, p_val_adj = TRUE, p_val_threshold
#  ) %>% 
#  bind_rows() 
#geneset_assessment_adjustedPval

#Set parallelization and modeling parameters.
top_n_target = 200 #defalut is 250

verbose = TRUE
cores_system = 1
n.cores = min(cores_system, celltype_de$cluster_id %>% unique() %>% length()) 

#Predicts ligand activities by matching receiver DE signatures to target gene models.
#Uses a ligand-target matrix built from prior knowledge
ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(
  get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de,
    receivers_oi = intersect(receivers_oi, celltype_de$cluster_id %>% unique()),
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    p_val_adj = p_val_adj,
    top_n_target = top_n_target,
    verbose = verbose, 
    n.cores = n.cores
  )
))

ligand_activities_targets_DEgenes$ligand_activities %>% head(20)

# Builds a list of distinct sender-receiver cell type pairs for further analysis or plotting.
ligand_activity_down = FALSE
sender_receiver_tbl = sender_receiver_de %>% distinct(sender, receiver)

#Extracts cell metadata (e.g., sample ID, group/stage, cell type) from the (sce), 
#and converts it into a tibble for easier manipulation.
metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()

#Checking for batch variables and creates a smiple table if no batch vairables 
if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group")
}

#Calls MultiNicheNet’s core function to prioritize ligand-receptor interactions using multiple criteria.
#Inputs:

#  Abundance & DE info for sender-receiver pairs.
#Ligand activity info (based on expressed targets).
#Experimental contrasts (e.g., Fibrosis vs. No_fibrosis).
#Logical flags for whether to focus on ligand downregulation (ligand_activity_down).
#Returns a list of prioritization tables, ranking ligand-receptor pairs by relevance.
prioritization_tables = suppressMessages(generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  scenario = "regular", # all prioritization criteria will be weighted equally
  fraction_cutoff = fraction_cutoff, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender,
  ligand_activity_down = ligand_activity_down
))

prioritization_tables$group_prioritization_tbl %>% head(20)


#Computes how well the ligand-target interactions (from prior knowledge) match with observed differential expression 
#correlation in receiver cells.
lr_target_prior_cor = lr_target_prior_cor_inference(
  receivers_oi = prioritization_tables$group_prioritization_tbl$receiver %>% unique(), 
  abundance_expression_info = abundance_expression_info, 
  celltype_de = celltype_de, 
  grouping_tbl = grouping_tbl, 
  prioritization_tables = prioritization_tables, 
  ligand_target_matrix = ligand_target_matrix, 
  logFC_threshold = logFC_threshold, 
  p_val_threshold = p_val_threshold, 
  p_val_adj = p_val_adj
)

#Combines all the core results into a single list for easier access.
multinichenet_output = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor
) 
multinichenet_output = make_lite_output(multinichenet_output)

#Extracts top 50 overall ligand-receptor pairs, regardless of condition group.
prioritized_tbl_oi_all = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  top_n = 50, 
  rank_per_group = FALSE
)

#Filters the group-level prioritization table to include only those top 50 pairs.
#Joins with full ranking table to keep prioritization scores.
prioritized_tbl_oi = 
  multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% 
  left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)

circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)


group_oi = "Fibrosis"
prioritized_tbl_oi_M_50 = prioritized_tbl_oi_all %>% 
  filter(group == group_oi)
plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi_M_50 %>% inner_join(lr_network_all)
)
plot_oi

ggsave("Fibrosis_bubbleplot_FINAL.svg", plot = plot_oi, width = 23, height = 12, dpi = 300)




group_oi = "No_fibrosis"
prioritized_tbl_oi_S_50 = prioritized_tbl_oi_all %>% 
  filter(group == group_oi) 
plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi_S_50 %>% inner_join(lr_network_all)
)
plot_oi

ggsave("No_Fibrosis_bubbleplot.JPEG", plot = plot_oi, width = 23, height = 12, dpi = 300)

ggsave(
  filename = "No_Fibrosis_bubbleplot.svg",
  plot = plot_oi,
  width = 23,   # increase width
  height = 12,   # adjust height if needed
  units = "in"  # inches
)
#ggsave("bubble_plot_fibrosis_H_to_T.png", plot = plot_oi, width = 18, height = 12, dpi = 300)

group_oi = "Fibrosis"
prioritized_tbl_oi_A_50 = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  50, 
  groups_oi = group_oi) 

plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi_A_50 %>% inner_join(lr_network_all)
)
plot_oi
