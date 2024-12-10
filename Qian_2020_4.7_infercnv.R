source("Qian_2020_0_packages_functions.R")
library(infercnv)
# usethis::edit_r_environ()

blueprint_CRC_cleaned <- readRDS("R_objects/Qian_2020/blueprint_CRC_cleaned.rds")
blueprint_CRC_cleaned <- JoinLayers(blueprint_CRC_cleaned)
# blueprint_CRC_cleaned_Epi <- subset(blueprint_CRC_cleaned, annotation %in% c(paste0("Epithelial_", 1:5), "Goblet"))


annotation_group <- blueprint_CRC_cleaned@meta.data %>% select(PatientNumber_TableS1, annotation) %>%
  mutate(group = case_when(
    annotation %in% c("CD4+ T_naïve", "Treg", "CD4+ T_mem/eff", "CD4+ T_CXCL13+", 
                      "CD8+ T_GZMK+", "CD8+ T_cyto_ex", "CD8+ T_prolif",
                      "gd-like T", "NK") ~ "T/NK",
    annotation %in% c("B_mature_naïve", "B_mem GC-dependent IgM-", "B_mem GC-dependent IgM+", "B_mem GC-independent", "B_prolif") ~ "B",
    annotation %in% c("PlasmaB_IgA+_mature_1", "PlasmaB_IgA+_mature_2", "PlasmaB_IgG+_immature", "PlasmaB_IgG+_mature") ~ "PlasmaB",
    annotation %in% c("SPP1+ TAM", "C1QC+MRC1- Mph", "C1QC+ Mph", "LYVE1+ Mph", 
                      "CD14+ Mono", "CD16+ Mono", "CD14+CD16+ Mono", 
                      "cDC1", "cDC2", "Migratory cDC", "pDC", "Mast",
                      "Other Myeloid") ~ "Myeloid",
    annotation %in% c(paste0("Fibroblast_", 1:3), "CAF_1", "CAF_2", "Pericytes", "Telocytes", "Myofibroblast") ~ "Fibroblast",
    annotation %in% c(paste0("Epithelial_", 1:5), "Goblet") ~ "Epithelial",
    annotation %in% c("Endothelial_1", "Endothelial_2") ~ "Endothelial",
    TRUE ~ "Glial")) %>% 
  mutate(group = ifelse(group %in% "Epithelial", paste0("Epithelial_", PatientNumber_TableS1), as.character(group))) %>%
  mutate(group = factor(group, c(paste0("Epithelial_CRC_", 1:7), "T/NK", "B", "PlasmaB", "Myeloid", "Fibroblast", "Endothelial", "Glial")))

blueprint_CRC_cleaned_Epi <- subset(blueprint_CRC_cleaned, annotation %in% c(paste0("Epithelial_", 1:5), "Goblet"))

## Using epithelial cells from adjacent normal samples as reference; run infercnv by patient  
infercnv_obj_Epi_run_by_patient <- sapply(unique(blueprint_CRC_cleaned$PatientNumber_TableS1), function(patient) {
  annot_file <- blueprint_CRC_cleaned_Epi@meta.data %>% 
    filter(PatientNumber_TableS1 == patient) %>%
    mutate(group = paste0("Epithelial_", PatientNumber_TableS1, "_", TumorSite)) %>%
    select(group) 
  
  count_mat <- blueprint_CRC_cleaned_Epi@assays$RNA$counts[, rownames(annot_file)]
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = count_mat,
    annotations_file = annot_file, 
    delim = "\t",
    gene_order_file = "R_objects/Qian_2020/gencode_v19_gene_pos.txt",
    ref_group_names = paste0("Epithelial_", patient, "_N"))
  
  infercnv::run(
    infercnv_obj,
    cutoff = 0, 
    out_dir = paste0("infercnv/", patient, "_cutoff.0.1_sd_denoise"),
    cluster_by_groups = TRUE, 
    HMM = FALSE,
    save_rds = TRUE,
    denoise = TRUE,
    sd_amplifier=3,  # sets midpoint for logistic
    noise_logistic=TRUE, # turns gradient filtering on
    png_res = 300)
  
  }, simplify = FALSE, USE.NAMES = TRUE)


# names(infercnv_obj_default) <- unique(blueprint_CRC_cleaned$orig.ident)
annot_pal_group <- c(
  "skyblue2", "plum2", "#9c88b8",
  "darkolivegreen3", "orange",
  "#FED9A6", "#8AD9B1FF",  "#E0C0AA", "#DDD08C"
)


infercnv_obj_Epi_run_by_patient <- sapply(
  unique(blueprint_CRC_cleaned$PatientNumber_TableS1),
  function(x) {readRDS(paste0("infercnv/", x, "_cutoff.0.1_sd_denoise/run.final.infercnv_obj")) },
  USE.NAMES = TRUE, simplify = FALSE
)

epithelial_B_C_barcodes <- rownames(subset(blueprint_CRC_cleaned_Epi@meta.data, TumorSite %in% c("B", "C")))

cnv_df_Epi <- lapply(infercnv_obj_Epi_run_by_patient, function(infercnv_obj) {
  ## calculate CNV score and CNV correlation
  expr.data_Epi <- log2(infercnv_obj@expr.data)
  expr.data_Epi_B_C <- expr.data_Epi[, intersect(epithelial_B_C_barcodes, colnames(expr.data_Epi))]
  
  ## CNV score, epithelial
  cnv_scores <- apply(expr.data_Epi_B_C, 2, function(x) {x %*% x})
  cnv_scores <- cnv_scores / dim(expr.data_Epi_B_C)[1]
  cnv_scores <- sort(cnv_scores, decreasing = T)
  
  ## CNV score, all cells
  cnv_scores.all <- apply(expr.data_Epi, 2, function(x) {x %*% x})
  cnv_scores.all <- cnv_scores.all / dim(expr.data_Epi)[1]
  cnv_scores.all <- sort(cnv_scores.all, decreasing = T)
  
  ## select the epithelial cells with top 2% (you may change to other numbers) CNV scores. They are very likely to be tumor cells. Calculate the average CNV profile
  top_2_percent_cells <- names(cnv_scores[1:as.integer(length(cnv_scores) * 0.02)])
  cnv_correlations_mean <- rowMeans(expr.data_Epi_B_C[, top_2_percent_cells])
  
  ## calculate CNV correlations between the CNV profile of each cell and the average CNV profile 
  cnv_correlations = apply(expr.data_Epi_B_C, 2, function(x) {cor(x, cnv_correlations_mean)})
  cnv_correlations <- cnv_correlations[names(cnv_scores)]
  
  cnv_correlations.all = apply(expr.data_Epi, 2, function(x) {cor(x, cnv_correlations_mean)})
  cnv_correlations.all <- cnv_correlations.all[names(cnv_scores.all)]
  
  
  cnv_df <- data.frame(cnv.score=cnv_scores.all, cnv.correlation=cnv_correlations.all) %>%
    rownames_to_column("barcodes") %>%
    left_join(blueprint_CRC_cleaned_Epi@meta.data %>% rownames_to_column("barcodes") %>%
                select(barcodes, annotation, TumorSite))
  return(cnv_df)
})


cnv_df_Epi %>%
  data.table::rbindlist(idcol = "PatientNumber") %>%
  ggplot(aes(x= cnv.score, y = cnv.correlation)) + 
  geom_point(aes(color= annotation), size = 0.2) + 
  # scale_color_manual(values = pal_Epi) +
  # geom_vline(xintercept = 0.003, color="blue", linetype = 2) + 
  geom_hline(yintercept = 0.3, color="blue", linetype = 2) + 
  theme_cowplot(font_size = 12) + 
  labs(x="CNV Score", y="CNV Correlation") +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  coord_cartesian(xlim = c(0, 0.005)) +
  facet_wrap(~PatientNumber)  
