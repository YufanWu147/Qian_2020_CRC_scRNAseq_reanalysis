source("Qian_2020_0_packages_functions.R")
load("R_objects/Qian_2020/blueprint_CRC_immune_rm_sub5_Patient31_DE_by_TumorSite.RData")
blueprint_CRC_immune_rm_sub5 <- readRDS("R_objects/Qian_2020/blueprint_CRC_immune_rm_sub5.rds")
blueprint_CRC_immune_rm_sub5[["joined"]] <- JoinLayers(blueprint_CRC_immune_rm_sub5[["RNA"]])
DefaultAssay(blueprint_CRC_immune_rm_sub5) <- "joined"

blueprint_CRC_immune_rm_sub5_Patient31 <- subset(blueprint_CRC_immune_rm_sub5, PatientNumber != 31)

# no. of immune cells by cell type, excluding Patient 31 (MSI-high)
blueprint_CRC_immune_rm_sub5_Patient31_ncells_by_TumorSite <- table(
  blueprint_CRC_immune_rm_sub5_Patient31$annotation, 
  blueprint_CRC_immune_rm_sub5_Patient31$TumorSite)

as.data.frame(blueprint_CRC_immune_rm_sub5_Patient31_ncells_by_TumorSite) %>% 
  pivot_wider(names_from = "Var2", values_from = "Freq") %>%
  cbind(Total = rowSums(blueprint_CRC_immune_rm_sub5_Patient31_ncells_by_TumorSite)) %>% select(-Var1)


blueprint_CRC_immune_rm_sub5_Patient31_DE_by_TumorSite <- mapply(
  function(site1, site2) {
    sapply(levels(blueprint_CRC_immune_rm_sub5_Patient31$annotation), function(celltype) {
      obj <- subset(blueprint_CRC_immune_rm_sub5_Patient31, annotation == celltype)
      ncell_site1 <- sum(obj$TumorSite == site1)
      ncell_site2 <- sum(obj$TumorSite == site2)
      if(ncell_site1 <= 10 | ncell_site2 <= 10) {return(NULL)}
      else {
        FindMarkers(obj, only.pos = FALSE, ident.1 = site1, ident.2 = site2, 
                    test.use = "wilcox", group.by = "TumorSite", min.pct = 0.1, logfc.threshold = 0) %>%
          rownames_to_column("gene")
      }
    }, USE.NAMES = TRUE, simplify = FALSE)
  },
  site1 = c("C", "C", "B"),
  site2 = c("B", "N", "N"),
  SIMPLIFY = FALSE) %>%
  `names<-`(c("Core_vs_Border", "Core_vs_AdjacentNormal", "Border_vs_AdjacentNormal"))

blueprint_CRC_immune_rm_sub5_Patient31_DE_by_TumorSite_negbinom <- mapply(
  function(site1, site2) {
    sapply(levels(blueprint_CRC_immune_rm_sub5_Patient31$annotation), function(celltype) {
      obj <- subset(blueprint_CRC_immune_rm_sub5_Patient31, annotation == celltype)
      ncell_site1 <- sum(obj$TumorSite == site1)
      ncell_site2 <- sum(obj$TumorSite == site2)
      if(ncell_site1 <= 10 | ncell_site2 <= 10) {return(NULL)}
      else {
        FindMarkers(obj, only.pos = FALSE, ident.1 = site1, ident.2 = site2, 
                    test.use = "negbinom", group.by = "TumorSite", 
                    min.pct = 0.1, logfc.threshold = 0) %>%
          rownames_to_column("gene")
      }
    }, USE.NAMES = TRUE, simplify = FALSE)
  },
  site1 = c("C", "C", "B"),
  site2 = c("B", "N", "N"),
  SIMPLIFY = FALSE) %>%
  `names<-`(c("Core_vs_Border", "Core_vs_AdjacentNormal", "Border_vs_AdjacentNormal"))


save(blueprint_CRC_immune_rm_sub5_Patient31_DE_by_TumorSite,
     blueprint_CRC_immune_rm_sub5_Patient31_DE_by_TumorSite_negbinom,
     file = "R_objects/Qian_2020/blueprint_CRC_immune_rm_sub5_Patient31_DE_by_TumorSite.RData")

DEG_negbinom_FDR.05 <- lapply(blueprint_CRC_immune_rm_sub5_Patient31_DE_by_TumorSite, function(ls) {
  lapply(ls, function(df) {
    if(!is.null(df)) {df$FDR <- p.adjust(df$p_val, method = "BH")}
    return(df)}) %>%
    data.table::rbindlist(idcol = "celltype") %>%
    mutate(celltype = factor(celltype, levels(blueprint_CRC_immune_rm_sub5_Patient31$annotation))) %>%
    filter(p_val_adj < 0.05)
})
  
DEG_wilcox_FDR.05 <- lapply(blueprint_CRC_immune_rm_sub5_Patient31_DE_by_TumorSite, function(ls) {
  lapply(ls, function(df) {
    if(!is.null(df)) {df$FDR <- p.adjust(df$p_val, method = "BH")}
    return(df)}) %>%
    data.table::rbindlist(idcol = "celltype") %>%
    mutate(celltype = factor(celltype, levels(blueprint_CRC_immune_rm_sub5_Patient31$annotation))) %>%
    filter(p_val_adj < 0.05)
})
  

# df <- FetchData(
#   blueprint_CRC_immune_rm_sub5_Patient31, unique(CorevsBorder_negbinom_FDR.05_top10$gene)) %>%
#   t()
# library(ComplexHeatmap)
# blueprint_CRC_immune_rm_sub5_Patient31$annot_site <- paste(
#   blueprint_CRC_immune_rm_sub5_Patient31$annotation,
#   blueprint_CRC_immune_rm_sub5_Patient31$TumorSite, sep = "_"
# ) 
# 
# p <- DotPlot(blueprint_CRC_immune_rm_sub5_Patient31,
#         unique(CorevsBorder_negbinom_FDR.05_top10$gene),
#         group.by = "annot_site") +
#   coord_flip() +
#   cowplot::theme_cowplot(font_size = 12) +
#   labs(x = NULL, y = NULL)
# p$data$annotation <- str_remove(p$data$id, "_(B|C|N)$") %>%
#   factor(levels = levels(blueprint_CRC_immune_rm_sub5_Patient31$annotation))
# p$data$id <- str_extract(p$data$id, "(B|C|N)$")
# p + facet_grid( ~ annotation) +
#   scale_color_gradientn(colors = rev(brewer.pal(11, "RdBu")))
# 
# df <- p$data %>%
#   filter(!str_detect(id, "_N$")) %>%
#   select(features.plot, id, avg.exp) %>%
#   pivot_wider(names_from = "id", values_from = "avg.exp") %>%
#   column_to_rownames("features.plot")
# 
# Heatmap(df %>% t() %>% scale() %>% t(),
#         show_column_dend = FALSE,
#         show_row_dend = FALSE,
#         row_names_gp = gpar(fontsize = 10, fontface = 3),
#         column_names_gp = gpar(fontsize = 10, fontface = 3),
#         column_names_rot = 0,
#         column_split = factor(str_remove(colnames(df), "_(B|C|N)$"),
#                               levels(blueprint_CRC_immune_rm_sub5_Patient31$annotation)),
#         column_labels = str_extract(colnames(df), "(B|C|N)$"),
#         column_title_gp = gpar(fontsize = 10),
#         column_title_rot = 45,
#         cluster_columns = FALSE)


df <- DEG_negbinom_FDR.05$Core_vs_Border %>%
  group_by(gene) %>% filter(n() > 1) %>%
  filter(!str_detect(gene, "^RP[SL]")) %>%
  select(celltype, gene, avg_log2FC) %>%
  pivot_wider(names_from = "celltype", values_from = "avg_log2FC") %>%
  column_to_rownames("gene") 

Heatmap(df,
        show_column_dend = FALSE,
        show_row_dend = FALSE,
        row_names_gp = gpar(fontsize = 9, fontface = 3),
        column_names_gp = gpar(fontsize = 10, fontface = 3),
        column_names_rot = 45,
        cluster_rows = FALSE,
        column_split = c(rep("T_NK", 7), rep("B_PlasmaB", 5), rep("Myeloid", 8)) %>%
          factor(c("T_NK", "B_PlasmaB", "Myeloid")),
        column_title_gp = gpar(fontsize = 11, fontface = 2),
        cluster_columns = FALSE,
        name = "avg_log2FC")


df <- DEG_wilcox_FDR.05$Core_vs_Border %>%
  group_by(gene) %>% filter(n() > 1) %>%
  filter(!str_detect(gene, "^RP[SL]")) %>%
  select(celltype, gene, avg_log2FC) %>%
  pivot_wider(names_from = "celltype", values_from = "avg_log2FC") %>%
  column_to_rownames("gene") 

Heatmap(df,
        show_column_dend = FALSE,
        show_row_dend = FALSE,
        row_names_gp = gpar(fontsize = 9, fontface = 3),
        column_names_gp = gpar(fontsize = 10, fontface = 3),
        column_names_rot = 45,
        cluster_rows = FALSE,
        column_split = c(rep("T_NK", 7), rep("B_PlasmaB", 6), rep("Myeloid", 3)) %>%
          factor(c("T_NK", "B_PlasmaB", "Myeloid")),
        column_title_gp = gpar(fontsize = 11, fontface = 2),
        cluster_columns = FALSE,
        name = "avg_log2FC")
