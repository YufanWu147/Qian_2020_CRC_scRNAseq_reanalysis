source("Qian_2020_0_packages_functions.R")

# blueprint_CRC_immune <- readRDS(file = "R_objects/Qian_2020/blueprint_CRC_immune_res1_updated_meta.rds")
# blueprint_CRC_immune_sub5_12 <- readRDS("R_objects/Qian_2020/blueprint_CRC_immune_sub5_12.rds")
# 
# blueprint_CRC_T_NK <- readRDS("R_objects/Qian_2020/blueprint_CRC_T_NK.rds")
# blueprint_CRC_Mono_Mac <- readRDS("R_objects/Qian_2020/blueprint_CRC_Mono_Mac_DC.rds")
# blueprint_CRC_B_PlasmaB <- readRDS("R_objects/Qian_2020/blueprint_CRC_B_PlasmaB.rds")
# 

blueprint_CRC_immune_cleaned <- readRDS("R_objects/Qian_2020/blueprint_CRC_immune_cleaned.rds")
blueprint_CRC_immune_cleaned$TumorLocation <- str_extract(blueprint_CRC_immune_cleaned$Pathological_subtype, "Left|Right")
blueprint_CRC_immune_cleaned$PatientNumber <- factor(blueprint_CRC_immune_cleaned$PatientNumber)
blueprint_CRC_immune_cleaned$SampleID_Site <- paste(
  blueprint_CRC_immune_cleaned$orig.ident,
  blueprint_CRC_immune_cleaned$TumorSite, sep = "_")
blueprint_CRC_immune_cleaned$TumorSite <- factor(
  blueprint_CRC_immune_cleaned$TumorSite, c("C", "B", "N")
)

# Assign annotation ----------------------------------------
sub5_12_annot <- blueprint_CRC_immune_sub5_12@meta.data %>%
  rownames_to_column("barcodes") %>%
  mutate(annot_sub5_12 = case_when(
    harmony_clusters_sub5_12 %in% c(0, 1, 4) ~ "",
    harmony_clusters_sub5_12 == 2 ~ "CD8+ T_prolif",
    harmony_clusters_sub5_12 == 3 ~ "B_prolif",
    harmony_clusters_sub5_12 == 5 ~ "Myeloid_prolif",
    harmony_clusters_sub5_12 == 6 ~ "PlasmaB_prolif"
  )) %>%
  select(barcodes, harmony_clusters_sub5_12, annot_sub5_12)

T_NK_annot <- blueprint_CRC_T_NK@meta.data %>%
  rownames_to_column("barcodes") %>%
  mutate(annot_T_NK = case_when(
    harmony_clusters_T_NK == 0 ~ "CD4+ T_naïve",
    harmony_clusters_T_NK == 1 ~ "Treg",
    harmony_clusters_T_NK == 5 ~ "CD4+ T_mem/eff",
    harmony_clusters_T_NK == 7 ~ "CD4+ T_CXCL13+",
    harmony_clusters_T_NK == 2 ~ "CD8+ T_GZMK+",
    harmony_clusters_T_NK == 4 ~ "CD8+ T_cyto_ex",
    harmony_clusters_T_NK == 6 ~ "CD8+ T_cyto",
    harmony_clusters_T_NK == 8 ~ "CD8+ T_prolif",
    harmony_clusters_T_NK == 10 ~ "gd-like T",
    harmony_clusters_T_NK == 9 ~ "NK",
    harmony_clusters_T_NK == 3 ~ ""               # low quality
  )) %>%
  select(barcodes, annot_T_NK)
  

Mono_Mac_annot <- blueprint_CRC_Mono_Mac@meta.data %>%
  rownames_to_column("barcodes") %>%
  mutate(annot_Mono_Mac = case_when(
    harmony_clusters_Mono_Mac == 0 ~ "C1QC+MRC1- Mph",
    harmony_clusters_Mono_Mac == 1 ~ "C1QC+ Mph",
    harmony_clusters_Mono_Mac == 2 ~ "SPP1+ TAM",
    harmony_clusters_Mono_Mac == 3 ~ "LYVE1+ Mph",
    harmony_clusters_Mono_Mac == 4 ~ "CD14+ Mono",
    harmony_clusters_Mono_Mac == 5 ~ "CD14+CD16+ Mono",
    harmony_clusters_Mono_Mac == 6 ~ "CD16+ Mono",
    harmony_clusters_Mono_Mac == 7 ~ ""
  )) %>%
  select(barcodes, annot_Mono_Mac)

B_PlasmaB_annot <- blueprint_CRC_B_PlasmaB@meta.data %>%
  rownames_to_column("barcodes") %>%
  mutate(annot_B_PlasmaB = case_when(
    harmony_clusters_B_PlasmaB == 5 ~ "B_mature_naïve",
    harmony_clusters_B_PlasmaB == 0 ~ "B_mem GC-dependent IgM+",
    harmony_clusters_B_PlasmaB == 9 ~ "B_mem GC-dependent IgM-",
    harmony_clusters_B_PlasmaB == 11 ~ "B_mem GC-independent",
    harmony_clusters_B_PlasmaB == 7 ~ "B_prolif",
    harmony_clusters_B_PlasmaB == 6  ~ "PlasmaB_IgG+_mature",
    harmony_clusters_B_PlasmaB == 8  ~ "PlasmaB_IgG+_immature",
    harmony_clusters_B_PlasmaB %in% c(1, 2) ~ "PlasmaB_IgA+_mature_1",
    harmony_clusters_B_PlasmaB %in% c(3, 4) ~ "PlasmaB_IgA+_mature_2",
    harmony_clusters_B_PlasmaB %in% c(12, 13) ~ "" # low-quality
  )) %>%
  select(barcodes, annot_B_PlasmaB) 


annotation_wrap <- function(x) {
  #str_replace(x, "_prolif", "\n_prolif") %>%
  str_replace(x, "B_mem ", "B_mem\n") %>%
    str_replace("(?<=\\+)_(?=mature|immature)", "\n_") %>%
    str_replace("(?<=^CD[48]\\+ T)_(?!ex)", "\n_") %>%
    str_replace(" Mono", "\nMono")
  }

annot <- blueprint_CRC_immune@reductions$umap.harmony@cell.embeddings %>%
  as.data.frame() %>%
  cbind(harmony_clusters_immune = blueprint_CRC_immune$harmony_clusters_immune) %>%
  rownames_to_column("barcodes") %>%
  left_join(sub5_12_annot, by = "barcodes") %>%
  left_join(T_NK_annot, by = "barcodes") %>%
  left_join(Mono_Mac_annot, by = "barcodes") %>%
  left_join(B_PlasmaB_annot, by = "barcodes") %>%
  mutate(annotation = ifelse(harmony_clusters_immune %in% c(5, 12), annot_sub5_12, NA)) %>%
  mutate(annotation = case_when(
    !is.na(annot_T_NK) ~ annot_T_NK,
    !is.na(annot_B_PlasmaB) ~ annot_B_PlasmaB,
    !is.na(annot_Mono_Mac) ~ annot_Mono_Mac,
    
    harmony_clusters_immune == 21 ~ "Neutrophil",
    harmony_clusters_immune == 24 ~ "cDC1",
    harmony_clusters_immune == 19 ~ "cDC2",
    harmony_clusters_immune == 22 ~ "Migratory cDC",
    harmony_clusters_immune == 25 ~ "pDC",
    harmony_clusters_immune == 13 ~ "Mast",
    
    TRUE ~ annotation
  )) %>%
  mutate(annotation = factor(annotation, levels = c(
    "CD4+ T_naïve", "Treg", "CD4+ T_mem/eff", "CD4+ T_CXCL13+",
    "CD8+ T_GZMK+", "CD8+ T_cyto_ex", "CD8+ T_prolif", 
    "gd-like T", "NK",
    "B_mature_naïve", "B_mem GC-dependent IgM-", "B_mem GC-dependent IgM+", "B_mem GC-independent", "B_prolif",
    "PlasmaB_IgA+_mature_1", "PlasmaB_IgA+_mature_2", "PlasmaB_IgG+_immature", "PlasmaB_IgG+_mature",
    "SPP1+ TAM", "C1QC+MRC1- Mph", "C1QC+ Mph", "LYVE1+ Mph",
    "CD14+ Mono", "CD16+ Mono", "CD14+CD16+ Mono", "Neutrophil", 
    "cDC1", "cDC2", "Migratory cDC", "pDC",
    "Mast", "Myeloid_prolif"
  ))) 
blueprint_CRC_immune$annotation <- annot$annotation



# Rerun UMAP
blueprint_CRC_immune_cleaned <- subset(blueprint_CRC_immune, annotation != "")
blueprint_CRC_immune_cleaned <- NormalizeData(blueprint_CRC_immune_cleaned, normalization.method = "LogNormalize", scale.factor = 10000) # try sctransformation
blueprint_CRC_immune_cleaned <- FindVariableFeatures(blueprint_CRC_immune_cleaned, selection.method = "vst", nfeatures = 2000)
blueprint_CRC_immune_cleaned <- ScaleData(blueprint_CRC_immune_cleaned)
blueprint_CRC_immune_cleaned <- RunPCA(blueprint_CRC_immune_cleaned, npcs = 100)

ElbowPlot(blueprint_CRC_immune_cleaned, 100)

## Perform Harmony integration
blueprint_CRC_immune_cleaned <- IntegrateLayers(
  object = blueprint_CRC_immune_cleaned, method = HarmonyIntegration,
  new.reduction = "harmony",
  verbose = TRUE, npcs = 60
)

blueprint_CRC_immune_cleaned <- RunUMAP(
  blueprint_CRC_immune_cleaned, reduction = "harmony", dims = 1:60,
  reduction.name = "umap.harmony", n.neighbors = 40)

saveRDS(blueprint_CRC_immune_cleaned, file = "R_objects/Qian_2020/blueprint_CRC_immune_cleaned.rds")

# Visualization ---------------------------------------------------------------
blueprint_CRC_immune_cleaned$annotation_label <- factor(
  annotation_wrap(blueprint_CRC_immune_cleaned$annotation), 
  annotation_wrap(levels(blueprint_CRC_immune_cleaned$annotation)))

annot_pal <- c(viridis::turbo(9)[c(6:9, 4:2)], "grey80", "grey25",
               "plum1", "plum3", "#EE3377", "#AA3377", "#420A68",
               "#DADAEB", "#9c88b8", "#a062cc", "#0d026e", 
               "goldenrod1", "#88D1EE", "darkolivegreen3", "darkgreen",
               "rosybrown1", "tomato3", "saddlebrown",
               "darkseagreen3", "#E78AC3", "lightblue3", "#FFD92F", "#BB8760", "#B2B8A3", "darkslategrey")

p <- DimPlot(blueprint_CRC_immune_cleaned, group.by = "annotation_label", 
             pt.size = 0.4, alpha = 0.45) +
  scale_color_manual(breaks = levels(blueprint_CRC_immune_cleaned$annotation_label), 
                     labels = levels(blueprint_CRC_immune_cleaned$annotation), 
                     values = annot_pal) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3, alpha = 1))) + labs(title = NULL)
LabelClusters(p, id = "annotation_label", size = 4, max.overlaps = Inf, lineheight = 0.9, hjust = 0.2,
              box = FALSE, fontface = 2, color = "black", box.padding = 0.6, min.segment.length = 1.5)
ggsave("figures/Qian_2020/final_annotation/DimPlot_final_annotation.png", 
       width = 11, height = 9)

p <- DimPlot(blueprint_CRC_immune_cleaned, group.by = "annotation", 
             pt.size = 0.1, alpha = 1) +
  scale_color_manual(values = annot_pal) +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3, alpha = 1))) + labs(title = NULL)
p$data$TumorSite <- blueprint_CRC_immune_cleaned$TumorSite
p + facet_grid(~ TumorSite) + theme(panel.border = element_rect(color = "black"))
ggsave("figures/Qian_2020/final_annotation/DimPlot_final_annotation_facet.png", 
       width = 15, height = 5)

for(position in c('stack', 'fill')) {
  perc_bar(blueprint_CRC_immune_cleaned@meta.data, by = "TumorSite", cluster = "annotation",
           folder = "figures/Qian_2020/final_annotation", suffix = "_immune_by_annotation", margin_l = 15,
           pal = ggpubfigs::friendly_pals$contrast_three[c(3, 2, 1)], width = 7, height = 3.5, position = position)
  perc_bar(blueprint_CRC_immune_cleaned@meta.data, by = "orig.ident", cluster = "annotation",
           folder = "figures/Qian_2020/final_annotation", suffix = "_immune_by_annotation",
           pal = c25, width = 8, height = 4, position = position, margin_l = 15)
  perc_bar(blueprint_CRC_immune_cleaned@meta.data, by = "PatientNumber", cluster = "annotation",
           folder = "figures/Qian_2020/final_annotation", suffix = "_immune_by_annotation",
           pal = c25, width = 8, height = 4, position = position, margin_l = 15)
  perc_bar(blueprint_CRC_immune_cleaned@meta.data, by = "Gender", cluster = "annotation",
           folder = "figures/Qian_2020/final_annotation", suffix = "_immune_by_annotation",
           pal = c25[2:1], width = 8, height = 4, position = position, margin_l = 15)
  perc_bar(blueprint_CRC_immune_cleaned@meta.data, by = "TumorLocation", cluster = "annotation",
           folder = "figures/Qian_2020/final_annotation", suffix = "_immune_by_annotation",
           pal = c25, width = 8, height = 4, position = position, margin_l = 15)
  
  perc_bar(blueprint_CRC_immune_cleaned@meta.data, by = "annotation", cluster = "PatientNumber",
           folder = "figures/Qian_2020/final_annotation", suffix = "celltype_by_Patient", position = position,
           width = 5, height = 5, use_default_pal = FALSE, pal = annot_pal)
  
  perc_bar(blueprint_CRC_immune_cleaned@meta.data, by = "annotation", cluster = "orig.ident",
           folder = "figures/Qian_2020/final_annotation", suffix = "celltype_by_Sample", position = position,
           width = 9, height = 5, facet = TRUE, facet_var = "TumorSite", pal = annot_pal)
  
  perc_bar(blueprint_CRC_immune_cleaned@meta.data, by = "annotation", cluster = "SampleID_Site",
           folder = "figures/Qian_2020/final_annotation", suffix = "celltype_by_Sample_group_by_Patient", position = position,
           width = 9, height = 5, facet = TRUE, facet_var = "PatientNumber", pal = annot_pal)
}




draw_dotplot(blueprint_CRC_immune_cleaned, Markers_Zhang_2023_FigS1, 
             group = "annotation", suffix = "_immune_Final_Annotation",
             folder = "Qian_2020/Final_annotation",
             pal = rev(brewer.pal(11, "RdBu")), fig_width = 12, fig_height = 11)


# Cell type proportion test -------------------------------------
# BiocManager::install("speckle")
# browseVignettes("speckle")
library(speckle)
library(limma)


# props$Proportions %>%
#   as.data.frame %>%
#   left_join(blueprint_CRC_immune_cleaned@meta.data %>% 
#               select(orig.ident, TumorSite) %>% distinct(), by = c("sample" = "orig.ident")) %>%
#   ggplot(aes(x = clusters, y = Freq, color = TumorSite)) +
#   geom_point(aes(shape = TumorSite), position = position_dodge(width = 0.7), size = 1.2) +
#   geom_boxplot(width = 0.7, outlier.alpha = 0) +
#   cowplot::theme_cowplot() + labs(x = NULL) +
#   scale_color_manual(values = ggpubfigs::friendly_pals$contrast_three[c(2, 3, 1)]) +
#   theme(axis.text.x = element_text(hjust = 1, angle = 45))
# 
# props$Proportions %>%
#   as.data.frame %>%
#   left_join(blueprint_CRC_immune_cleaned@meta.data %>% 
#               select(orig.ident, TumorSite) %>% distinct(), by = c("sample" = "orig.ident")) %>%
#   ggplot(aes(x = sample, y = Freq, color = clusters)) +
#   geom_point(aes(shape = TumorSite), position = position_dodge(width = 0.7), size = 1.2) +
#   cowplot::theme_cowplot() + labs(x = NULL) +
#   scale_color_manual(values = annot_pal) +
#   facet_grid(~TumorSite) +
#   theme(axis.text.x = element_text(hjust = 1, angle = 45))

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

celltype_prop_ttest <- get_prop_ttest_res(blueprint_CRC_immune_cleaned, design = "~ 0 + TumorSite + Gender + Molecular_status + TumorLocation",
                                          design_names = c("C", "B", "N", "MalevsFemale", "MSIhigh_vs_MSS", "LeftvsRight"))
celltype_prop_ttest_rm_Patient31 <- get_prop_ttest_res(subset(blueprint_CRC_immune_cleaned, PatientNumber != 31),
                                                       design = "~ 0 + TumorSite + Gender + TumorLocation",
                                                       design_names = c("C", "B", "N", "MalevsFemale", "LeftvsRight"))
celltype_prop_ttest.FDR.05 <- lapply(celltype_prop_ttest$ttest_res, subset, FDR < 0.05) 
celltype_prop_ttest_rm_Patient31.FDR.05 <- lapply(celltype_prop_ttest_rm_Patient31$ttest_res, subset, FDR < 0.05) 

celltype_prop_ttest$ttest_res %>%
  lapply(as.data.frame) %>%
  lapply(rownames_to_column, "annotation") %>%
  lapply(select, -starts_with("PropMean")) %>%
  data.table::rbindlist(idcol = "Contrast") %>%
  mutate(Contrast = factor(Contrast, c("Core_vs_Normal", "Border_vs_Normal", "Core_vs_Border"))) %>%
  ggplot(aes(y = fct_reorder2(annotation, desc(Contrast), desc(PropRatio)), x = PropRatio)) +
  geom_point() +
  facet_grid( ~ Contrast) +
  labs(y = NULL) +
  theme_bw()

annotation_group <- data.frame(annotation = levels(blueprint_CRC_immune_cleaned$annotation)) %>%
  mutate(group = case_when(
    annotation %in% c(
      "CD4+ T_naïve", "Treg", "CD4+ T_mem/eff", "CD4+ T_CXCL13+",
      "CD8+ T_GZMK+", "CD8+ T_cyto_ex", "CD8+ T_prolif", 
      "gd-like T", "NK") ~ "T/NK",
    annotation %in% c(
      "B_mature_naïve", "B_mem GC-dependent IgM-", "B_mem GC-dependent IgM+", "B_mem GC-independent", "B_prolif",
      "PlasmaB_IgA+_mature_1", "PlasmaB_IgA+_mature_2", "PlasmaB_IgG+_immature", "PlasmaB_IgG+_mature"
    ) ~ "B/PlasmaB",
    TRUE ~ "Myeloid")) %>%
  mutate(group = factor(group, c("T/NK", "B/PlasmaB", "Myeloid")))

lapply(list(celltype_prop_ttest, celltype_prop_ttest_rm_Patient31), function(res) {
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
    ggplot(aes(y = fct_reorder(annotation, PropRatio), 
               # y = factor(annotation, rev(blueprint_CRC_immune_cleaned$annotation %>% levels)),
               x = PropRatio)) +
    geom_vline(xintercept = 1, color = "grey50", linetype = 2) +
    geom_point(aes(color = Significance)) +
    facet_grid(group ~ Contrast, scales = "free_y", space = "free_y") +
    scale_x_continuous(transform = "log2", breaks = c(0.2, 0.5, 1, 2, 5, 25),
                       labels = function(x) {round(x, 1)}) +
    scale_color_manual(values = c("red", "pink", "black")) +
    labs(y = NULL, x = "Proportion Ratio") +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid = element_line(linetype = 3, color = "grey75", linewidth = 0.3))
  ggsave(paste0("figures/Qian_2020/final_annotation/propeller_ttest",
                ifelse(nrow(res$sample_info) == 18, "_rm_Patient31", ""), ".png"),
         width = 9, height = 6)
  
  ttest_res.05 <- lapply(res$ttest_res, subset, FDR < 0.05)
  res$props$Proportions %>%
    as.data.frame %>%
    filter(clusters %in% unlist(sapply(ttest_res.05, rownames))) %>%
    left_join(res$sample_info, by = c("sample" = "orig.ident")) %>%
    ggplot(aes(x = TumorSite, y = Freq)) +
    geom_point(aes(color = Molecular_status, shape = Gender, fill = Gender, group = "Gender"), 
               position = position_dodge(width = 0.8), stroke = 0.8, size = 0.8) +
    scale_shape_manual(values = c(3, 19)) +
    scale_color_manual(values = c("red", "black"), breaks = c("MSI-high", "MSS")) +
    facet_wrap(~clusters, scales = "free_y", ncol = 2) +
    cowplot::theme_cowplot(font_size = 12) +
    labs(y = "Proportion")
  
  ggsave(paste0("figures/Qian_2020/final_annotation/stripchart_propeller_ttest_FDR.05",
                ifelse(nrow(res$sample_info) == 18, "_rm_Patient31", ""), ".png"),
         width = 6, height = 5.5)
})





ncell_by_sample_annot <- blueprint_CRC_immune_cleaned@meta.data %>%
  count(annotation, orig.ident, PatientNumber, TumorSite) %>% ungroup %>%
  left_join(annotation_group, by = "annotation")

ncell_by_sample_annot_outliers <- ncell_by_sample_annot %>%
  group_by(annotation, group, TumorSite) %>%
  mutate(q1 = quantile(n, 0.25), q3 = quantile(n, 0.75)) %>%
  filter(n > q3 + 1.5 * (q3 - q1) | n < q1 - 1.5 * (q3 - q1)) %>%
  full_join(distinct(ncell_by_sample_annot[, c("annotation", "group", "TumorSite")]))

ncell_by_sample_annot %>%
  ggplot(aes(x = annotation, y = n, color = TumorSite)) +
  geom_boxplot(width = 0.8, alpha = 0.5, outlier.alpha = 0) +
  geom_point(aes(shape = PatientNumber, group = TumorSite), stroke = 0.8,
             position = position_dodge(width = 0.8), size = 1.3) +
  # geom_text(data = ncell_by_sample_annot_outliers, aes(label = PatientNumber), show.legend = FALSE,
  #           position = position_dodge(width = 0.7), size = 3.5, vjust = -0.3, fontface = "bold") +
  scale_color_manual(values = ggpubfigs::friendly_pals$contrast_three[c(3, 2, 1)]) +
  # scale_fill_manual(values = ggpubfigs::friendly_pals$contrast_three[c(3, 2, 1)]) +
  scale_shape_manual(values = c(16:17, 3:4, 7:9)) +
  cowplot::theme_cowplot(font_size = 12) + labs(x = NULL, y = "# Cells") +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  theme(axis.text.x = element_text(hjust = 1, angle = 45),
        plot.margin = margin(5, 5, 5, 15))
ggsave("figures/Qian_2020/final_annotation/boxplot_ncells_by_sample_annot.png", width = 14, height = 4.5)


prop_by_sample_annot <- celltype_prop_ttest$props$Proportions %>%
  as.data.frame() %>%
  left_join(celltype_prop_ttest$sample_info, by = c("sample" = "orig.ident")) %>%
  left_join(blueprint_CRC_immune_cleaned@meta.data[, c("orig.ident", "PatientNumber")] %>% distinct(),
            by = c("sample" = "orig.ident")) %>%
  left_join(annotation_group, by = c("clusters" = "annotation"))

prop_by_sample_annot_outliers <- prop_by_sample_annot %>%
  group_by(clusters, group, TumorSite) %>%
  mutate(q1 = quantile(Freq, 0.25), q3 = quantile(Freq, 0.75)) %>%
  filter(Freq > q3 + 1.5 * (q3 - q1) | Freq < q1 - 1.5 * (q3 - q1)) %>%
  full_join(distinct(prop_by_sample_annot[, c("clusters", "group", "TumorSite")]))


prop_by_sample_annot %>%
  ggplot(aes(x = clusters, y = Freq, color = TumorSite)) +
  # geom_point(aes(shape = TumorSite), position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.8, alpha = 0.5, outlier.alpha = 0) +
  geom_point(aes(shape = PatientNumber, group = TumorSite), stroke = 0.8,
             position = position_dodge(width = 0.8), size = 1.3) +
  # geom_text(data = prop_by_sample_annot_outliers, aes(label = PatientNumber), show.legend = FALSE,
  #           position = position_dodge(width = 0.7), size = 3.5, vjust = -0.3, fontface = "bold") +
  scale_color_manual(values = ggpubfigs::friendly_pals$contrast_three[c(2, 3, 1)]) +
  # scale_fill_manual(values = ggpubfigs::friendly_pals$contrast_three[c(2, 3, 1)]) +
  scale_shape_manual(values = c(16:17, 3:4, 7:9)) +
  cowplot::theme_cowplot(font_size = 12) + labs(x = NULL, y = "Proportion") +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  theme(axis.text.x = element_text(hjust = 1, angle = 45),
        plot.margin = margin(5, 5, 5, 15))
ggsave("figures/Qian_2020/final_annotation/boxplot_prop_by_sample_annot.png", width = 14, height = 4.5)


  # scale_color_manual(values = rev(ggpubfigs::friendly_pals$contrast_three[3:1]))
  # NoLegend()
# celltype_prop_test_ANOVA <- propeller( 
#   clusters = blueprint_CRC_immune_cleaned$annotation, 
#   sample = blueprint_CRC_immune_cleaned$orig.ident, 
#   group = blueprint_CRC_immune_cleaned$TumorSite)
# 
# celltype_prop_test_ANOVA_rm_Patient31 <- propeller(
#   clusters = subset(blueprint_CRC_immune_cleaned, PatientNumber != 31)$annotation, 
#   sample = subset(blueprint_CRC_immune_cleaned, PatientNumber != 31)$orig.ident, 
#   group = subset(blueprint_CRC_immune_cleaned, PatientNumber != 31)$TumorSite)
# 
# celltype_prop_test <- mapply(function(site1, site2) {
#   obj <- subset(blueprint_CRC_immune_cleaned, TumorSite %in% c(site1, site2))
#   obj$TumorSite <- factor(obj$TumorSite, c(site2, site1))
#   
#   propeller(clusters = obj$annotation, sample = obj$orig.ident, 
#             group = obj$TumorSite)}, 
#   site1 = c("C", "C", "B"), site2 = c("B", "N", "N"),
#   SIMPLIFY = FALSE) %>%
#   `names<-`(c("Core_vs_Border", "Core_vs_AdjacentNormal", 
#               "Border_vs_AdjacentNormal"))
# 
# celltype_prop_test_rm_Patient31 <- mapply(function(site1, site2) {
#   obj <- subset(blueprint_CRC_immune_cleaned, TumorSite %in% c(site1, site2))
#   obj <- subset(obj, PatientNumber != 31)
#   obj$TumorSite <- factor(obj$TumorSite, c(site2, site1))
#   
#   propeller(clusters = obj$annotation, sample = obj$orig.ident, 
#             group = obj$TumorSite)}, 
#   site1 = c("C", "C", "B"), site2 = c("B", "N", "N"),
#   SIMPLIFY = FALSE) %>%
#   `names<-`(c("Core_vs_Border", "Core_vs_AdjacentNormal", 
#               "Border_vs_AdjacentNormal"))
# 
# celltype_prop_test_rm_Patient31$Core_vs_Border
# 
