
#### load all the necessary packages for multinichenetr, nichenet and single cell data analysis ###

library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(nichenetr)
library(multinichenetr)
library(nichenetr) 
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(dplyr)


# Loading the curated Omnipath ligand-receptor reference database for human for analysis 
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

# loading T cells and Endothelial object into R environment through Seurats readRDS function followed by the path to file. 
T_cell_and_Endothelial_subset_new <- readRDS("Path/to/Seurat/object")

#  Some cell type names had different structures, so we used gsub to standardize them
# by replacing spaces with underscores, reducing redundancy in downstream analysis.
T_cell_and_Endothelial_subset_new$cell_type <- gsub(" ", "_", T_cell_and_Endothelial_subset_new$cell_type) 

# The Seurat object contains multiple assays (e.g., RNA, ADT, etc.), 
# but since the analysis focuses on scRNA-seq data, we set the default assay to "RNA".
DefaultAssay(T_cell_and_Endothelial_subset_new) <- "RNA"

### because of space in one of the cell type names####
#### The reclutered T cell and endo has all the cell types so the next steps are:
# 1- subset for just Endo and T cells 
# 2- subset to take out "PBMCs" from it 

# filtering the "T_cell_and_Endothelial_subset_new" object to only keep the desired cell 
cell_types_retained <- c(
  "CD4 Naïve T cells", "CD4 Proliferating T cells", "CD4 T cells",
  "CD8 Cytotoxic Tcells", "CD8 Effector Memory",
  "Blood vessels",
  "ACKR1 Endothelial cells", "TAGLN Endothelial cells",
  "Hepatocytes"
)

# Subsetting the object
# Use Seurat's subset function
T_cell_and_Endothelial_filtered <- subset(
  T_cell_and_Endothelial_subset_new,
  subset = cell_type %in% c("CD4_T_cells", "CD4_Naïve_T_cells", "CD8_Effector_Memory", "Hepatocytes", "CD8_Cytotoxic_Tcells", "CD4_Proliferating_T_cells", "Blood_vessels", "TAGLN_Endothelial_cells", "ACKR1_Endothelial_cells")
)


# second subsetting to take out data of PBMCs

T_cell_and_Endothelial_final <- subset(
  T_cell_and_Endothelial_filtered,
  subset = Tissue != "PBMC"
)

# Save your object
saveRDS(T_cell_and_Endothelial_final, file = "/data/home/bt24974/T_cell_and_Endothelial_final.rds")

#### Exploratory Data Anlysis on Suerat Object ######

DimPlot(T_cell_and_Endothelial_final, group.by = "cell_type", label = TRUE)

# Normalize
T_cell_and_Endothelial_final <- NormalizeData(T_cell_and_Endothelial_final)

# Variable features
T_cell_and_Endothelial_final <- FindVariableFeatures(T_cell_and_Endothelial_final)

# Scaling
T_cell_and_Endothelial_final <- ScaleData(T_cell_and_Endothelial_final)

# PCA
T_cell_and_Endothelial_final <- RunPCA(T_cell_and_Endothelial_final)

# UMAP
T_cell_and_Endothelial_final <- RunUMAP(T_cell_and_Endothelial_final, dims = 1:20)

# Clustering
T_cell_and_Endothelial_final <- FindNeighbors(T_cell_and_Endothelial_final, dims = 1:20)
T_cell_and_Endothelial_final <- FindClusters(T_cell_and_Endothelial_final, resolution = 0.5)

# Visualise
DimPlot(T_cell_and_Endothelial_final, group.by = "cell_type", label = TRUE)


# Create a new metadata column "Stage.tissue" by concatenating the values 
# from the "Stage" and "Tissue" columns, separated by an underscore. 
# This provides a combined identifier (e.g., "Fibrosis_Liver") that can be useful 
# for grouping, stratifying, or visualising cells by both Stage and Tissue together
T_cell_and_Endothelial_final$Stage.tissue <- paste0(T_cell_and_Endothelial_final$Stage,"_", T_cell_and_Endothelial_final$Tissue)

# Convert Suerat Object into a Single Cell Experiment object because multinichenetr 
# requires data in SCE format to properly extract expression and meta data for predicted 
# ligand-receptor inference.
# Here "RNA" assay was explicitly used to ensure single cell RNAseq was analysed
sce = Seurat::as.SingleCellExperiment(T_cell_and_Endothelial_final, assay = "RNA")

# Create a new metadata column "patient_id_tissue_stage" by concatenating the values 
# from the "Patient_ID", Stage" and "Tissue" columns, separated by an underscore. 
# This provides a combined identifier (e.g., "x_SAT_fibrosis") that can be useful 
# for grouping, stratifying, or visualising cells by and patient_ID Stage and Tissue together
colData(sce)$patient_id_tissue_stage <- paste0(
  colData(sce)$Patient_ID, "_",
  colData(sce)$Tissue, "_",
  colData(sce)$Stage
)

# setting sample, group and celltype IDs for downstream analysis
sample_id = "patient_id_tissue_stage" # identifies each sample and which cell came from which sample, 
# which tissue type and what stage 
group_id = "Stage"   # for the biological conditions of interest, Fibrosis Vs No_fibrosis
celltype_id = "cell_type" # cell type annotations cluster for each cell type 


covariates = NA 
batches = NA 



# Clean up metadata column values by converting them into syntactically valid names.
# `make.names()` ensures that characters such as spaces, dashes, or special symbols 
# are replaced with underscores or adjusted so the values become valid R variable names. 
# This step prevents errors in downstream MultinicheNet functions (and other analyses) 
# that expect metadata categories like patient_id_tissue_stage, tissue_stage, and cell_type 
# to be free of spaces or special characters
SummarizedExperiment::colData(sce)$patient_id_tissue_stage = SummarizedExperiment::colData(sce)$patient_id_tissue_stage  %>% make.names()
SummarizedExperiment::colData(sce)$tissue_stage = SummarizedExperiment::colData(sce)$tissue_stage %>% make.names()
SummarizedExperiment::colData(sce)$cell_type = SummarizedExperiment::colData(sce)$cell_type %>% make.names()


# setting the contrast for conditions of interest
# Contrats setting must follow Muscat guidelines. Multiple contrasts are placed inside a single string.
# Contrasts are separated by commas (no spaces).
# The entire string is enclosed in double quotes 
contrasts_oi = c("'Fibrosis-No_fibrosis','No_fibrosis-Fibrosis'") 

# creates a data frame that maps contrasts set above to its corresponding group 
contrast_tbl = tibble(
  contrast = c("Fibrosis-No_fibrosis", "No_fibrosis-Fibrosis"),
  group = c("Fibrosis", "No_fibrosis")
) 

# knew biologically which cells we wanted to set as senders and which as recievrs based on disease progression
# hence these are set specifically 
# once this was investigated the senders and recievers were also switched around

senders_oi = c("ACKR1_Endothelial_cells", "Blood_vessels", "TAGLN_Endothelial_cells")

receivers_oi = c(
  "CD8_Cytotoxic_Tcells", 
  "CD8_Effector_Memory",
  "CD4_Naïve_T_cells",
  "CD4_T_cells",
  "CD4_Proliferating_T_cells",
  "Hepatocytes"
  
  
)

conditions_keep = c("Fibrosis","No_fibrosis") # Specifies conditions to keep for the analysis 


# filters the SCE object to only keep conditions that were specified, Fibrosis Vs No_fibrosis 
sce = sce[, SummarizedExperiment::colData(sce)[,group_id] %in% 
            conditions_keep
] 


# Sets minimum number of cells per group per sample per condition 
# DEFAULT 10 but changed it to 1 because one of the samples has less than 10
#Important for statistical power and validity—too few cells can lead to unreliable results.
min_cells = 1

# To make the names in sample_id syntactically correct. applied to the Sample Ids.
SummarizedExperiment::colData(sce)[, sample_id] <- make.names(SummarizedExperiment::colData(sce)[, sample_id])


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

abundance_info$abund_plot_sample
ggsave("abundance_endofinal.png", plot = abundance_info$abund_plot_sample, width = 28, height = 12, dpi = 600)

######## Visualisations for cell Abundance #########

abundance_info$abund_plot_sample +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 10) # for facet titles
  )



###   Number of cells per niche EDA ###
df_tissue <- as.data.frame(table(T_cell_and_Endothelial_final$Tissue))
colnames(df_tissue) <- c("Tissue", "CellCount")

ggplot(df_tissue, aes(x = Tissue, y = CellCount, fill = Tissue)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Number of cells per Tissue")

SampleStructurePlot(
  seurat_object = T_cell_and_Endothelial_final,
  group.by = "cell_type",
  split.by = "patient_id_tissue_stage"
)

cell_counts <- T_cell_and_Endothelial_final@meta.data %>%
  group_by(cell_type, Stage) %>%  # or replace Stage with contrast
  summarise(n_cells = n(), .groups = "drop")

# Visualise 
ggplot(cell_counts, aes(x = cell_type, y = n_cells, fill = Stage)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Number of cells per cell type per Stage",
       x = "Cell Type",
       y = "Number of cells")


### Filtering based on condition and cell type specificity for downstream aalysis 
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

print("absent celltypes:")

#Controls whether to include condition-specific cell types in the analysis.
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


min_sample_prop = 0.25
fraction_cutoff = 0.03 # default 0.05 but parameter changed as in some samples the number of cells is lower 

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
logFC_threshold = 0.25
p_val_threshold = 0.05


p_val_adj = FALSE  # default in vignette is FALSE 

# scores gene set based on differential expression across different conditions as set in the 
# contrasts_oi 
geneset_assessment = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de, logFC_threshold, p_val_adj, p_val_threshold
  ) %>% 
  bind_rows() 
geneset_assessment


############ Visualisation for upregulated and downregulated 
############ genes across cell type clusters and stages 

ggplot(geneset_assessment, aes(x = reorder(cluster_id, -prop_geneset_up), y = prop_geneset_up, fill = contrast)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Proportion of Geneset Upregulated Across Clusters",
    x = "Cell Type Cluster",
    y = "Proportion of Geneset Upregulated"
  ) +
  theme_minimal()


geneset_assessment_long <- geneset_assessment %>%
  select(cluster_id, contrast, prop_geneset_up, prop_geneset_down) %>%
  pivot_longer(cols = starts_with("prop_geneset"), names_to = "regulation", values_to = "proportion") %>%
  mutate(regulation = ifelse(regulation == "prop_geneset_up", "Upregulated", "Downregulated"))

ggplot(geneset_assessment_long, aes(x = reorder(cluster_id, -proportion), y = proportion, fill = regulation)) +
  geom_col(position = "dodge") +
  coord_flip() +
  facet_wrap(~ contrast) +
  labs(
    title = "Proportion of Geneset Up/Downregulated Across Clusters",
    x = "Cell Type Cluster",
    y = "Proportion"
  ) +
  theme_minimal()



# data preparation for a diverging bar plot
highlight_data <- geneset_assessment %>%
  mutate(
    proportion = prop_geneset_up - prop_geneset_down,  # net direction
    label = case_when(
      in_range_up ~ "Significant Up",
      in_range_down ~ "Significant Down",
      TRUE ~ "Not Significant"
    )
  )

# Plot
ggplot(highlight_data, aes(x = reorder(cluster_id, proportion), y = proportion, fill = label)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c(
    "Significant Up" = "firebrick",
    "Significant Down" = "steelblue",
    "Not Significant" = "gray70"
  )) +
  labs(
    title = "Net Regulation of Gene Set Expression Across Clusters",
    x = "Cell Type Cluster",
    y = "Net Proportion (Up - Down)",
    fill = "Significance"
  ) +
  theme_minimal()




top_n_target = 350 #default is 250 but wanted to investiagte mmore cell to cell communication 

# setting parameters for ligand-receptor interaction prediction 
#  Sets the number of cores used for parallelisation.
# It takes the minimum of available system cores and the number 
# of unique cell clusters to prevent over-allocation of resources.
verbose = TRUE
cores_system = 1
n.cores = min(cores_system, celltype_de$cluster_id %>% unique() %>% length()) 

#Predicts ligand activities by matching receiver DE signatures to target gene models.
#Uses a ligand-target matrix built from prior knowledge
ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(
  get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de,
    receivers_oi = intersect(receivers_oi, celltype_de$cluster_id %>% unique()),
    ligand_target_matrix = ligand_target_matrix, # prior knowledge reference 
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    p_val_adj = p_val_adj,
    top_n_target = top_n_target,
    verbose = verbose, 
    n.cores = n.cores
  )
))

# Displays the top 20 ligands with the highest predicted activity.
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

# Extract all unique senders and receivers from the prioritized 
# ligand–receptor table (prioritized_tbl_oi). 
# The union ensures every cell type that is either a sender 
# or receiver is included once. Sorting creates a consistent order
senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

# -----------------------------------------------------------
# Assign colors to sender and receiver cell types
# -----------------------------------------------------------

# Generate a color palette (Spectral) with as many distinct colors 
# as there are unique senders/receivers. 
# Colors are then named by the corresponding cell types 
# to ensure consistent mapping across visualisations.
colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)

# -----------------------------------------------------------
# Visualise sender–receiver interactions in a circos plot
# -----------------------------------------------------------

# Create a circos plot comparing cell–cell communication 
# between conditions of interest. This visualization highlights:
#   - Which cell types act as senders and receivers
#   - The strength and directionality of communication
#   - Condition-specific differences in interactions
circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)


# ---------------------------------------------------------------------------
# Visualise sender–receiver interactions and gene expression in a bubble plot 
# ---------------------------------------------------------------------------

group_oi = "Fibrosis"
prioritized_tbl_oi_M_50 = prioritized_tbl_oi_all %>% 
  filter(group == group_oi)
plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi_M_50 %>% inner_join(lr_network_all)
)
plot_oi
ggsave("bubble_plot_fibrosis_ENDO_LAST_new.png", plot = plot_oi, width = 26, height = 12, dpi = 600)


group_oi = "No_fibrosis"
prioritized_tbl_oi_S_50 = prioritized_tbl_oi_all %>% 
  filter(group == group_oi) 
plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi_S_50 %>% inner_join(lr_network_all)
)
plot_oi
ggsave("bubble_plot_NO_fib_new_ENDO_LAST.png", plot = plot_oi, width = 26, height = 12, dpi = 300)


# -----------------------------------------------------------
#                   Intercellular network
# -----------------------------------------------------------


lr_target_prior = prioritized_tbl_oi_all %>% inner_join(
  multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>%
    distinct(ligand, target, direction_regulation, contrast) %>% inner_join(contrast_tbl) %>% ungroup() 
) 
lr_target_df = lr_target_prior %>% distinct(group, sender, receiver, ligand, receptor, id, direction_regulation) 

lr_target_df %>% filter(receptor %in% union(lr_network$ligand, lr_network$receptor))

lr_target_df <- lr_target_df %>%
  dplyr::rename("target" = "receptor")

# Undo the rename
colnames(lr_target_df)
# If it's currently "target", rename it back:
lr_target_df <- lr_target_df %>%
  dplyr::rename("receptor" = "target")

network = infer_intercellular_regulatory_network(lr_target_df, prioritized_tbl_oi_all)
lr_target_df <- lr_target_df %>%
  dplyr::mutate(target = receptor)

network = infer_intercellular_regulatory_network(lr_target_df, prioritized_tbl_oi_all)
network$links %>% head()

network$nodes %>% head()

#colors_sender["ACKR1_Endothelial_cells", "Endothelial_cells", "TAGLN_Endothelial_cells"] = "pink" # 
colors_sender["ACKR1_Endothelial_cells"] <- "pink"
colors_sender["Endothelial_cells"] <- "lightblue"
colors_sender["TAGLN_Endothelial_cells"] <- "lightgreen"
#the  original yellow background with white font is not very readable
network_graph = visualize_network(network, colors_sender)
network_graph$plot

lr_target_df %>%
  filter(group == "Fibrosis") %>%
  semi_join(prioritized_tbl_oi_all, by = c("ligand", "target")) %>%
  nrow()


lr_network_df <- prioritized_tbl_oi %>%
  filter(prioritization_score > quantile(prioritization_score, 0.75)) %>%  # keep top 25%
  select(sender, receiver, ligand, receptor, group, prioritization_score)
edges <- lr_network_df %>%
  mutate(
    from = sender,
    to = receiver,
    interaction = paste0(ligand, "-", receptor)
  ) %>%
  select(from, to, interaction, prioritization_score)

# Create igraph object
lr_graph <- graph_from_data_frame(d = edges, directed = TRUE)

# Optional: assign vertex types (sender vs receiver)
V(lr_graph)$type <- ifelse(V(lr_graph)$name %in% edges$from, "Sender", "Receiver")


ggraph(lr_graph, layout = 'fr') +  # You can try 'kk', 'circle', etc.
  geom_edge_link(aes(width = prioritization_score), arrow = arrow(length = unit(3, 'mm')), end_cap = circle(3, 'mm')) +
  geom_node_point(aes(color = type), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  theme_void() +
  labs(title = "Ligand-Receptor Interaction Network")



# ----------------------------------------------
# Further EDA of the meta data in Seurat Object 
# ----------------------------------------------



# Compute total cells per group (fibrosis/no-fibrosis)
ggplot(abundance_info$abundance_data, aes(x = celltype_id, y = n, fill = group_id)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Number of Cells") +
  xlab("Cell Type") +
  ggtitle("Cell Counts by Cell Type and Condition") +
  geom_text(data = percent_data, 
            aes(label = paste0(round(percent, 1), "%"), y = n + 10), 
            position = position_dodge(width = 0.9), 
            vjust = 0)




# Get the data into a new variable 
df <- abundance_info$abundance_data

# Step 1: Identify the most abundant cell type per sample
top_celltypes <- df %>%
  group_by(sample_id) %>%
  top_n(1, n) %>%         # or slice_max(n, n = 1)
  ungroup()

# Custom distinct colors per cell type
custom_colors <- c(
  "ACKR1_Endothelial_cells" = "#E64B35",      # Red
  "CD4_Naïve_T_cells"       = "#F39B7F",      # Orange
  "CD4_T_cells"             = "#3C5488",      # Navy
  "CD8_Cytotoxic.Tcells"    = "#00A087",      # Teal
  "CD8_Effector_Memory"     = "#4DBBD5",      # Cyan
  "Endothelial_cells"       = "#8491B4"       # Purple-blue
)

# Plot with manual fill scale
ggplot(top_celltypes, aes(x = sample_id, y = n, fill = celltype_id)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ group_id, scales = "free_x") +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, size = 6),
    strip.text = element_text(size = 10)
  ) +
  labs(
    title = "Most Abundant Cell Type in Each Sample",
    subtitle = "Grouped by Condition (Fibrosis vs No-Fibrosis)",
    x = "Sample",
    y = "Cell Count",
    fill = "Top Cell Type"
  )


#Plotting the abundance cell type information
abundance_info$abund_plot_sample

abundance_info$abundance_data # to check for how many cells of each cell type in each sample are present. 

# saving as csv file 
write.csv(abundance_info$abundance_data, "celltype_counts_per_sample.csv", row.names = FALSE)

# visualising 

ggplot(abundance_info$abundance_data, aes(x = celltype_id, y = n)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Number of Cells") +
  xlab("Cell Type") +
  ggtitle("Cell Counts per Cell Type")


ggplot(abundance_info$abundance_data, aes(x = sample_id, y = n, fill = celltype_id)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Number of Cells") +
  xlab("Sample") +
  ggtitle("Cell Counts per Sample by Cell Type")


#### options to make th per sample less crowded and more understand easily

ggplot(abundance_info$abundance_data, aes(x = sample_id, y = n)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ celltype_id, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, size = 6),
        strip.text = element_text(size = 8)) +
  labs(title = "Cell Counts per Sample", x = "Sample", y = "Number of Cells")

### function for calling all  cell types ###

plot_patient_tissues <- function(data, patient_id) {
  df <- subset(data, grepl(paste0("^", patient_id), sample_id))
  if (nrow(df) == 0) {
    message("Patient ID not found.")
    return(NULL)
  }
  df$tissue <- sub(".*\\.", "", df$sample_id)
  
  lapply(unique(df$tissue), function(tissue_name) {
    ggplot(subset(df, tissue == tissue_name), aes(x = celltype_id, y = n, fill = celltype_id)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste(patient_id, "-", tissue_name),
           x = "Cell Type", y = "Cell Count") +
      guides(fill = "none")
  })
}


plots <- plot_patient_tissues(abundance_info$abundance_data, "X10203.SAT")
for (p in plots) print(p)

library(ggplot2)
library(dplyr)

df <- abundance_info$abundance_data

# Reorder group_id so No_fibrosis is first
df$group_id <- factor(df$group_id, levels = c("No_fibrosis", "Fibrosis"))

ggplot(df, aes(x = group_id, y = n, color = group_id)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_jitter(aes(fill = group_id), width = 0.2, shape = 21, size = 2, alpha = 0.8) +
  facet_wrap(~ celltype_id, scales = "free_y") +
  scale_color_manual(values = c("No_fibrosis" = "cyan3", "Fibrosis" = "tomato")) +
  scale_fill_manual(values = c("No_fibrosis" = "cyan3", "Fibrosis" = "tomato")) +
  theme_minimal(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey80", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Cell type abundances per group",
    x = "Group",
    y = "Cells per sample-celltype combination",
    color = "Condition",
    fill = "Condition"
  )




