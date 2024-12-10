source("Qian_2020_0_packages_functions.R")
blueprint_CRC_cleaned <- readRDS("R_objects/Qian_2020/blueprint_CRC_cleaned.rds")

# blueprint_CRC <- readRDS(file = "R_objects/Qian_2020/blueprint_CRC_harmony.rds")
# blueprint_CRC_immune_cleaned <- readRDS("R_objects/Qian_2020/blueprint_CRC_immune_cleaned.rds")


fibroblast_clusters <- c(7, 11, 12, 18, 22, 23, 25, 30, 31)
epithelial_clusters <- c(0, 1, 9, 15, 27, 29)
immune_clusters <- c(
  8, 16,    # B cells,
  4, 14, 20,# Plasma B
  2, 3, 10, # T cells
  26,       # NK
  21,       # Mast
  6, 13, 19,# Myeloid
  24        # Proliferating
)

annotation_wrap <- function(x) {
  #str_replace(x, "_prolif", "\n_prolif") %>%
  str_replace(x, "B_mem ", "B_mem\n") %>%
    str_replace("(?<=\\+)_(?=mature|immature)", "\n_") %>%
    str_replace("(?<=^CD[48]\\+ T)_(?!ex)", "\n_") %>%
    str_replace(" Mono", "\nMono")
}

annot_all <- blueprint_CRC@meta.data %>%
  rownames_to_column("barcodes") %>%
  left_join(blueprint_CRC_immune_cleaned@meta.data %>% rownames_to_column("barcodes") %>% 
              select(barcodes, annotation_immune = annotation, annotation_immune_label = annotation_label),
            by = "barcodes") %>%
  mutate(annotation = case_when(
    !is.na(annotation_immune) ~ as.character(annotation_immune),
    harmony_clusters %in% immune_clusters ~ "low-quality",
    harmony_clusters == 7 ~ "Fibroblast_1",
    harmony_clusters == 22 ~ "Fibroblast_2",
    harmony_clusters == 23 ~ "Fibroblast_3",
    harmony_clusters == 12 ~ "CAF_1",
    harmony_clusters == 30 ~ "CAF_2",
    harmony_clusters == 11 ~ "Pericytes",
    harmony_clusters == 25 ~ "Telocytes",
    harmony_clusters == 31 ~ "Myofibroblast",
    harmony_clusters == 18 ~ "low-quality",
    harmony_clusters == 0 ~ "Epithelial_1",
    harmony_clusters == 1 ~ "Epithelial_2",
    harmony_clusters == 9 ~ "Epithelial_3",
    harmony_clusters == 15 ~ "Epithelial_4",
    harmony_clusters == 29 ~ "Epithelial_5",
    harmony_clusters == 27 ~ "Goblet",
    harmony_clusters == 17 ~ "Glial",
    harmony_clusters == 5 ~ "Endothelial_1",
    harmony_clusters == 28 ~ "Endothelial_2"
    )) %>%
  mutate(annotation = factor(
    annotation, c(levels(blueprint_CRC_immune_cleaned$annotation),
    paste0("Fibroblast_", 1:3), "CAF_1", "CAF_2", "Pericytes", "Telocytes", "Myofibroblast",
    paste0("Epithelial_", 1:5), "Goblet", "Endothelial_1", "Endothelial_2", "Glial", "low-quality"))) %>%
  mutate(annotation_label = factor(annotation_wrap(annotation), annotation_wrap(levels(annotation))))

blueprint_CRC$annotation <- annot_all$annotation
blueprint_CRC$annotation_label <- annot_all$annotation_label
blueprint_CRC$annotation_label_number <- paste0("c", str_pad(as.numeric(blueprint_CRC$annotation), 2, "left", pad = "0"))
blueprint_CRC$annotation_ordered <- paste0(blueprint_CRC$annotation_label_number,
                                           "_", as.character(blueprint_CRC$annotation))
  

blueprint_CRC_cleaned <- subset(blueprint_CRC, annotation != "low-quality")
blueprint_CRC_cleaned <- NormalizeData(blueprint_CRC_cleaned, normalization.method = "LogNormalize", scale.factor = 10000) # try sctransformation
blueprint_CRC_cleaned <- FindVariableFeatures(blueprint_CRC_cleaned, selection.method = "vst", nfeatures = 2000)
blueprint_CRC_cleaned <- ScaleData(blueprint_CRC_cleaned)
blueprint_CRC_cleaned <- RunPCA(blueprint_CRC_cleaned, npcs = 100)

ElbowPlot(blueprint_CRC_cleaned, 100)

## Perform Harmony integration
blueprint_CRC_cleaned <- IntegrateLayers(
  object = blueprint_CRC_cleaned, method = HarmonyIntegration,
  new.reduction = "harmony",
  verbose = TRUE, npcs = 70
)

blueprint_CRC_cleaned <- RunUMAP(
  blueprint_CRC_cleaned, reduction = "harmony", dims = 1:70,
  reduction.name = "umap.harmony", n.neighbors = 40)

blueprint_CRC_cleaned$annotation <- droplevels(blueprint_CRC_cleaned$annotation)
blueprint_CRC_cleaned$annotation_label <- droplevels(blueprint_CRC_cleaned$annotation_label)
blueprint_CRC_cleaned@meta.data <- blueprint_CRC_cleaned@meta.data %>%
  rownames_to_column("barcode") %>%
  left_join(PatientID_mapping, by = c("PatientNumber")) %>%
  left_join(TableS1, by = c("PatientNumber_TableS1" = "Patient_number")) %>%
  column_to_rownames("barcode")

blueprint_CRC_cleaned$TumorLocation <- str_extract(blueprint_CRC_cleaned$Pathological_subtype, "Left|Right")
blueprint_CRC_cleaned$PatientNumber <- factor(blueprint_CRC_cleaned$PatientNumber)
blueprint_CRC_cleaned$SampleID_Site <- paste(
  blueprint_CRC_cleaned$orig.ident,
  blueprint_CRC_cleaned$TumorSite, sep = "_")
blueprint_CRC_cleaned$TumorSite <- factor(
  blueprint_CRC_cleaned$TumorSite, c("C", "B", "N")
)

saveRDS(blueprint_CRC_cleaned, file = "R_objects/Qian_2020/blueprint_CRC_cleaned.rds")

# Visualization ------------------------------------------------------------------
annot_pal_immune <- c(viridis::turbo(9)[c(6:9, 4:2)], "grey80", "grey25",
                      "plum1", "plum3", "#EE3377", "#AA3377", "#420A68",
                      "#DADAEB", "#9c88b8", "#a062cc", "#0d026e", 
                      "goldenrod1", "#88D1EE", "darkolivegreen3", "darkgreen",
                      "rosybrown1", "tomato3", "saddlebrown",
                      "darkseagreen3", "#E78AC3", "lightblue3", "#FFD92F", "#BB8760", "#B2B8A3", "darkslategrey")
annot_pal_fibroblast <- c("#FED9A6", "#FBB4AE", "#C5B69B", "#BBDAB4", "#A2BCD2", "#FFDADA", "#ECC9DB", "#CDBAD3")
annot_pal_epi <- c(viridis::mako(15)[c(15, 13, 11, 9, 7)], "#CBD5E8")
annot_pal_endo_glial <- c("#E0C0AA", "#AAAAAA", "#DDD08C")
annot_pal <- c(annot_pal_immune, annot_pal_fibroblast, annot_pal_epi, annot_pal_endo_glial)

annotation_group <- data.frame(annotation = levels(blueprint_CRC_cleaned$annotation)) %>%
  mutate(group.l3 = case_when(
    annotation %in% c("CD4+ T_naïve", "Treg", "CD4+ T_mem/eff", "CD4+ T_CXCL13+") ~ "TCD4",
    annotation %in% c("CD8+ T_GZMK+", "CD8+ T_cyto_ex", "CD8+ T_prolif") ~ "TCD8", 
    annotation == "gd-like T" ~ "gdT", 
    annotation == "NK" ~ "NK",
    annotation %in% c("B_mature_naïve", "B_mem GC-dependent IgM-", "B_mem GC-dependent IgM+", "B_mem GC-independent", "B_prolif") ~ "B",
    annotation %in% c("PlasmaB_IgA+_mature_1", "PlasmaB_IgA+_mature_2", "PlasmaB_IgG+_immature", "PlasmaB_IgG+_mature") ~ "PlasmaB",
    annotation %in% c("SPP1+ TAM", "C1QC+MRC1- Mph", "C1QC+ Mph", "LYVE1+ Mph") ~ "Macrophage",
    annotation %in% c("CD14+ Mono", "CD16+ Mono", "CD14+CD16+ Mono") ~ "Monocyte", 
    annotation %in% c("cDC1", "cDC2", "Migratory cDC", "pDC") ~ "DC",
    annotation == "Mast" ~ "Mast",
    annotation %in% c("Neutrophil", "Myeloid_prolif") ~ "Other Myeloid",
    annotation %in% c(paste0("Fibroblast_", 1:3), "CAF_1", "CAF_2", "Pericytes", "Telocytes", "Myofibroblast") ~ "Fibroblast",
    annotation %in% c(paste0("Epithelial_", 1:5), "Goblet") ~ "Epithelial",
    annotation %in% c("Endothelial_1", "Endothelial_2") ~ "Endothelial",
    TRUE ~ "Glial")) %>%
  mutate(group.l2 = case_when(
    group.l3 %in% c("TCD4", "TCD8", "gdT", "NK") ~ "T/NK",
    group.l3 %in% c("B", "PlasmaB") ~ "B/PlasmaB",
    group.l3 %in% c("Macrophage", "Monocyte", "DC", "Other Myeloid") ~ "Myeloid",
    TRUE ~ group.l3)) %>%
  mutate(group.l1 = ifelse(group.l2 %in% c("T/NK", "B/PlasmaB", "Myeloid"), "Immune", group.l2)) %>%
  mutate(group.l1 = factor(group.l1, c("Immune", "Fibroblast", "Epithelial", "Endothelial", "Glial")))  %>%
  mutate(group.l2 = factor(group.l2, c("T/NK", "B/PlasmaB", "Myeloid", "Fibroblast", "Epithelial", "Endothelial", "Glial"))) %>%
  mutate(group.l3 = factor(group.l3, c("TCD4", "TCD8", "gdT", "NK", "B", "PlasmaB", "Macrophage", "Monocyte", "DC", "Mast", "Other Myeloid", "Fibroblast",  "Epithelial", "Endothelial", "Glial")))

annot_group <- blueprint_CRC_cleaned@meta.data %>%
  select(annotation) %>%
  left_join(annotation_group, by = "annotation")
blueprint_CRC_cleaned$group.l1 <- annot_group$group.l1
blueprint_CRC_cleaned$group.l2 <- annot_group$group.l2
blueprint_CRC_cleaned$group.l3 <- annot_group$group.l3

annot_pal_group.l3 <- c(
  "orangered2", "dodgerblue", "grey80", "grey25", "plum2", "#9c88b8",
  "darkolivegreen3", "orange",
  "lightblue3", 
  "#B2B8A3", "#FFD92F",
  "#FED9A6", "#8AD9B1FF",  "#E0C0AA", "#DDD08C"
)

p <- DimPlot(blueprint_CRC_cleaned, group.by = "annotation_label_number") +
  scale_color_manual(values = annot_pal,
                     breaks = sort(unique(blueprint_CRC_cleaned$annotation_label_number)), 
                     labels = sort(unique(blueprint_CRC_cleaned$annotation_ordered))) +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
LabelClusters(p, id = "annotation_label_number", fontface = 2, box.padding = 0.5,
              max.overlaps = Inf, min.segment.length = 1, segment.alpha = 0.7)
ggsave("figures/Qian_2020/all_cells/DimPlot_annotation.png", width = 12, height = 9)


for(position in c('stack', 'fill')) {
  perc_bar(blueprint_CRC_cleaned@meta.data, by = "TumorSite", cluster = "annotation",
           folder = "figures/Qian_2020/all_cells", suffix = "_immune_by_annotation", margin_l = 15,
           pal = ggpubfigs::friendly_pals$contrast_three[c(3, 2, 1)], width = 9, height = 3.5, position = position)
  perc_bar(blueprint_CRC_cleaned@meta.data, by = "orig.ident", cluster = "annotation",
           folder = "figures/Qian_2020/all_cells", suffix = "_immune_by_annotation",
           pal = c25, width = 9, height = 4, position = position, margin_l = 15)
  perc_bar(blueprint_CRC_cleaned@meta.data, by = "PatientNumber", cluster = "annotation",
           folder = "figures/Qian_2020/all_cells", suffix = "_immune_by_annotation",
           pal = c25, width = 9, height = 4, position = position, margin_l = 15)
  perc_bar(blueprint_CRC_cleaned@meta.data, by = "Gender", cluster = "annotation",
           folder = "figures/Qian_2020/all_cells", suffix = "_immune_by_annotation",
           pal = c25[2:1], width = 9, height = 4, position = position, margin_l = 15)
  perc_bar(blueprint_CRC_cleaned@meta.data, by = "TumorLocation", cluster = "annotation",
           folder = "figures/Qian_2020/all_cells", suffix = "_immune_by_annotation",
           pal = c25, width = 9, height = 4, position = position, margin_l = 15)
  
  perc_bar(blueprint_CRC_cleaned@meta.data, by = "annotation", cluster = "PatientNumber", legend_ncol = 2,
           folder = "figures/Qian_2020/all_cells", suffix = "celltype_by_Patient", position = position,
           width = 5, height = 5, use_default_pal = FALSE, pal = annot_pal)
  
  perc_bar(blueprint_CRC_cleaned@meta.data, by = "annotation", cluster = "orig.ident", legend_ncol = 2,
           folder = "figures/Qian_2020/all_cells", suffix = "celltype_by_Sample", position = position,
           width = 9, height = 5, facet = TRUE, facet_var = "TumorSite", pal = annot_pal)
  
  perc_bar(blueprint_CRC_cleaned@meta.data, by = "annotation", cluster = "SampleID_Site", legend_ncol = 2,
           folder = "figures/Qian_2020/all_cells", suffix = "celltype_by_Sample_group_by_Patient", position = position,
           width = 9, height = 5, facet = TRUE, facet_var = "PatientNumber", pal = annot_pal)
  perc_bar(blueprint_CRC_cleaned@meta.data, by = "group.l3", cluster = "SampleID_Site", legend_ncol = 2,
           folder = "figures/Qian_2020/all_cells", suffix = "celltype_by_Sample_group_by_Patient", position = position,
           width = 9, height = 5, facet = TRUE, facet_var = "PatientNumber", pal = annot_pal_group.l3)
}




# Proportion test ----------------------------------------------------------------
get_prop_ttest_res <- function(obj, design, design_names) {
  props <- getTransformedProps(obj$annotation, as.character(obj$orig.ident), transform="logit") 
  
  sample_info <- data.frame(orig.ident = colnames(props$Proportions)) %>%
    left_join(obj@meta.data %>% select(orig.ident, TumorSite, Gender, Molecular_status, TumorLocation) %>% distinct,
              by = "orig.ident")
  
  designAS <- model.matrix(as.formula(design), data = sample_info)
  colnames(designAS) <- c(design_names)
  
  ncell <- ncol(obj)
  nsample <- length(unique(obj$orig.ident))
  prop.list <- speckle::convertDataToList(
    props$Proportions, data.type="proportions",
    transform="logit", scale.fac=ncell/nsample)
  
  ttest_res <- list(
    Core_vs_Border = propeller.ttest(
      prop.list = prop.list, design = designAS, contrasts = makeContrasts("C-B", levels = designAS),
      robust = TRUE, trend = FALSE, sort = TRUE),
    Core_vs_Normal = propeller.ttest(
      prop.list = prop.list, design = designAS, contrasts = makeContrasts("C-N", levels = designAS),
      robust = TRUE, trend = FALSE, sort = TRUE),
    Border_vs_Normal = propeller.ttest(
      prop.list = prop.list, design = designAS, contrasts = makeContrasts("B-N", levels = designAS),
      robust = TRUE, trend = FALSE, sort = TRUE)) 
  
  return(list(ttest_res = ttest_res, props = props, sample_info = sample_info))
}

celltype_prop_ttest_allcells <- get_prop_ttest_res(
  blueprint_CRC_cleaned, design = "~ 0 + TumorSite + Gender + Molecular_status + TumorLocation",
  design_names = c("C", "B", "N", "MalevsFemale", "MSIhigh_vs_MSS", "LeftvsRight"))

celltype_prop_ttest_allcells_rm_Patient31 <- get_prop_ttest_res(
  subset(blueprint_CRC_cleaned, PatientNumber != 31),
  design = "~ 0 + TumorSite + Gender + TumorLocation",
  design_names = c("C", "B", "N", "MalevsFemale", "LeftvsRight"))


lapply(list(celltype_prop_ttest_allcells, celltype_prop_ttest_allcells_rm_Patient31), function(res) {
  res$ttest_res %>%
    lapply(as.data.frame) %>%
    lapply(rownames_to_column, "annotation") %>%
    lapply(select, -starts_with("PropMean")) %>%
    data.table::rbindlist(idcol = "Contrast") %>%
    mutate(Contrast = factor(Contrast, c("Core_vs_Border", "Core_vs_Normal", "Border_vs_Normal"))) %>%
    left_join(annotation_group, by = "annotation") %>%
    mutate(Significance = case_when(FDR < 0.05 ~ "FDR < 0.05", 
                                    FDR < 0.1 ~ "0.05 ≤ FDR < 0.1", 
                                    TRUE ~ "n.s.")) %>%
    mutate(Significance = factor(Significance, c("FDR < 0.05", "0.05 ≤ FDR < 0.1", "n.s."))) %>%
    ggplot(aes(x = PropRatio, y = fct_reorder(annotation, PropRatio))) +
    geom_vline(xintercept = 1, color = "grey50", linetype = 2) +
    geom_point(aes(color = Significance)) +
    facet_grid(group ~ Contrast, scales = "free_y", space = "free_y") +
    scale_x_continuous(transform = "log2", breaks = c(0.1, 0.5, 1, 2, 10, 25),
                       labels = function(x) {round(x, 1)}) +
    scale_color_manual(values = c("red", "pink", "black")) +
    labs(y = NULL, x = "Proportion Ratio") +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid = element_line(linetype = 3, color = "grey75", linewidth = 0.3))
  ggsave(paste0("figures/Qian_2020/all_cells/propeller_ttest",
                ifelse(nrow(res$sample_info) == 18, "_rm_Patient31", ""), ".png"),
         width = 9, height = 8)
  
  ttest_res.05 <- lapply(res$ttest_res, subset, FDR < 0.05)
  res$props$Proportions %>%
    as.data.frame %>%
    filter(clusters %in% unlist(sapply(ttest_res.05, rownames))) %>%
    left_join(res$sample_info, by = c("sample" = "orig.ident")) %>% 
    ggplot(aes(x = TumorSite, y = Freq)) +
    geom_point(aes(color = Molecular_status, shape = Gender, fill = Gender, group = "Gender"), 
               position = position_dodge(width = 0.8), stroke = 0.8, size = 0.8) +
    geom_line(aes(colour = Molecular_status, group = sample)) +
    scale_shape_manual(values = c(3, 19)) +
    scale_color_manual(values = c("red", "black"), breaks = c("MSI-high", "MSS")) +
    facet_wrap(~ clusters, scales = "free_y", ncol = 4) +
    cowplot::theme_cowplot(font_size = 12) +
    labs(y = "Proportion")
  
  ggsave(paste0("figures/Qian_2020/all_cells/stripchart_propeller_ttest_FDR.05",
                ifelse(nrow(res$sample_info) == 18, "_rm_Patient31", ""), ".png"),
         width = 10, height = 8)
})


DimPlot(blueprint_CRC_cleaned, group.by = "group.l3", cols = annot_pal_group.l3, label = TRUE, repel = TRUE, label.size = 3.5) 
ggsave('figures/Qian_2020/all_cells/DimPlot_group.l3.png', width = 7, height = 6)
