source("Qian_2020_0_packages_functions.R")

blueprint_CRC <- readRDS(file = "R_objects/Qian_2020/blueprint_CRC_harmony.rds")
blueprint_CRC_Epithelial_markers <- readRDS("R_objects/Qian_2020/blueprint_CRC_Epithelial_markers.rds")

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

# 1. Fibroblast -----------------------------------------------------------
blueprint_CRC_Fibroblast <- subset(blueprint_CRC, harmony_clusters %in% fibroblast_clusters) # 33694 genes X 8205 cells
blueprint_CRC_Fibroblast[["joined"]] <- JoinLayers(blueprint_CRC_Fibroblast[["RNA"]])
DefaultAssay(blueprint_CRC_Fibroblast) <- "joined"

DimPlot(blueprint_CRC_Fibroblast, cols = brewer.pal(9, "Set1"), label = TRUE, repel = TRUE)

blueprint_CRC_Fibroblast_markers <- FindAllMarkers(
  blueprint_CRC_Fibroblast, logfc.threshold = 0.25, only.pos = TRUE, min.pct = 0.25)
saveRDS(blueprint_CRC_Fibroblast_markers, "R_objects/Qian_2020/blueprint_CRC_Fibroblast_markers.rds")

blueprint_CRC_Fibroblast_markers_top5 <- blueprint_CRC_Fibroblast_markers %>%
  group_by(cluster) %>% filter(p_val_adj < 0.01) %>%
  arrange(desc(avg_log2FC)) %>%
  slice(1:5)

blueprint_Fig3f_Fibroblast_markers <- list(
  C1_KCNN3 = c("KCNN3", "P2RY1", "THBS4", "SPRY2", "LY6H"),
  C2_ADAMDEC1 = c("ADAMDEC1", "APOE", "CCL8", "FABP5", "HAPLN1"),
  C3_SOX6 = c("SOX6", "BMP4", "BMP5", "WNT5A", "FRZB"),
  C7_MYH11 = c("MYH11", "RERGL", "ACTG2", "SORBS2"),
  C8_RGS5 = c("RGS5", "PDGFRB", "NDUFA4L2", "NOTCH3"),
  C9_CFD = c("CFD", "APOD", "MFAP5", "PI16"),
  C10_COMP = c("COMP", "CTHRC1", "FN1", "COL10A1", "COL11A1"),
  C11_SERPINE1 = c("SERPINE1", "CLDN1", "PTGIS", "RGS4", "CADM3", "WT1", "IGF1"))  %>%
  lapply(as.data.frame) %>%
  data.table::rbindlist(idcol = "cluster") %>%
  mutate(cluster = str_replace(cluster, "_", "_\n")) %>%
  `colnames<-`(c("cluster", "gene")) %>%
  mutate(cluster = factor(cluster, unique(.$cluster))) %>%
  mutate(gene = factor(gene, unique(.$gene)))

p <- DotPlot(blueprint_CRC_Fibroblast, blueprint_Fig3f_Fibroblast_markers$gene,
        cols = c("grey90", "brown2")) + 
  coord_flip() + scale_x_discrete(limits = rev) + 
  # scale_color_gradientn(colors = brewer.pal(9, "YlOrBr")) +
  theme(axis.text.y = element_text(face = 3),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = NULL, y = NULL)
p$data <- p$data %>% 
  left_join(blueprint_Fig3f_Fibroblast_markers, by= c("features.plot" = "gene"))
p + facet_grid(cluster ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme(strip.text.y.left = element_text(size = 11),
        panel.spacing.y = unit(1, "mm"))
ggsave("figures/Qian_2020/Dotplot_Fibroblast_markers_Qian_2020_Fig3f.png",
       width = 6, height = 7)




DotPlot(blueprint_CRC_Fibroblast, unique(blueprint_CRC_Fibroblast_markers_top5$gene),
        cols = c("grey90", "brown2")) + 
  coord_flip() + scale_x_discrete(limits = rev) + 
  # scale_color_gradientn(colors = brewer.pal(9, "YlOrBr")) +
  theme(axis.text.y = element_text(face = 3),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = NULL, y = NULL)
ggsave("figures/Qian_2020/Dotplot_blueprint_CRC_Fibroblast_top5Markers.png",
       width = 5, height = 7)


Fibroblast_markers <- c(
  "COL1A1", "COL1A2",        # Fibroblast
  "INHBA",                   # CAF (Pelka et al 2021)
  "FAP", "COL10A1", "MMP11", "CTHRC1", "SLC12A8", "F2R", "COL12A1", # CAF (Kang et al 2024)
  "ACTA2", "TAGLN", "PDGFA", # CAF-A
  "MMP2", "DCN",   # CAF-B
  "S100A4", "DPT", "MFAP5", "SFRP1", "SFRP2")    # NMF

DotPlot(blueprint_CRC_Fibroblast, Fibroblast_markers,
        cols = c("grey90", "brown2")) + 
  coord_flip() + scale_x_discrete(limits = rev) + 
  theme(axis.text.y = element_text(face = 3),
        axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = NULL, y = NULL)
ggsave("figures/Qian_2020/Dotplot_Fibroblast_markers_simplified.png",
       width = 5, height = 5.5)

Fibroblast_markers_Qi_2022 <- c(
  "COL1A1", "COL3A1", 
  "ACTG2", "MYH11", "MFAP5", "DES", # Myofibroblasts
  "F3", "FOXL1", "SOX6", "ICAM1", # Telocytes
  "CD24", "RSPO3",    # CD24+ fibroblasts
  "NT5E",             # NT5E+ fibroblast
  "FAP",              # FAP+ fibroblast (canonical CAF activation marker)
  "FGFR2", "DPP4",    # FGFR2+ fibroblast
  "MCAM", 
  "MKI67")            # proliferating fibroblast
  
Fibroblast_top_markers_Qi_2022 <- readxl::read_xlsx("data/signatures/Qi_2022_Suppl_Data_2.xlsx", sheet = "MSC")
Fibroblast_top5_markers_Qi_2022 <- Fibroblast_top_marker %>%
  filter(p_val_adj < 0.01) %>%
  # arrange(desc(avg_logFC)) %>%
  group_by(cluster) %>%
  slice(1:5) %>% 
  group_by(gene) %>% slice(1) 

DotPlot(blueprint_CRC_Fibroblast, fibroblast_markers_Qi_2022,
        cols = c("grey90", "brown2")) + 
  coord_flip() + scale_x_discrete(limits = rev) + 
  # scale_color_gradientn(colors = brewer.pal(9, "YlOrBr")) +
  theme(axis.text.y = element_text(face = 3),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = NULL, y = NULL)
ggsave("figures/Qian_2020/Dotplot_Fibroblast_markers_Qi_2022.png",
       width = 5, height = 3.5)

p <- DotPlot(blueprint_CRC_Fibroblast, Fibroblast_top5_markers_Qi_2022$gene,
        cols = c("grey90", "brown2")) + 
  coord_flip() + scale_x_discrete(limits = rev) + 
  # scale_color_gradientn(colors = brewer.pal(9, "YlOrBr")) +
  theme(axis.text.y = element_text(face = 3),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = NULL, y = NULL)
p$data <- p$data %>% 
  left_join(Fibroblast_top5_markers_Qi_2022, by= c("features.plot" = "gene")) %>%
  arrange(cluster) %>%
  mutate(cluster = str_replace(cluster, " (?=(fibroblast|myofibroblast|telocytes))", "\n")) %>%
  mutate(cluster = factor(cluster, levels = c(
    "Pericytes", "DES+\nmyofibroblasts", "MFAP5+\nmyofibroblasts",
    "ICAM1+\ntelocytes", "ICAM1-\ntelocytes",
    "CD24+\nfibroblasts", "CD73+\nfibroblasts", "FAP+\nfibroblasts",
    "FGFR2+\nfibroblasts",
    "Proliferating\nfibroblasts"
  )))
p + facet_grid(cluster ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme(strip.text.y.left = element_text(angle = 0, size = 11),
        panel.spacing.y = unit(2, "mm"))
ggsave("figures/Qian_2020/Dotplot_Fibroblast_markers_Qi_2022_top10DEGs.png",
       width = 6, height = 8)

perc_bar(blueprint_CRC_Fibroblast@meta.data, by = "TumorSite",
         folder = "figures/Qian_2020", suffix = "_fibroblast",
         pal = ggpubfigs::friendly_pals$contrast_three[c(2, 3, 1)], width = 3, height = 2)


# 2. Epithelial -----------------------------------------------------------
blueprint_CRC_Epithelial <- subset(blueprint_CRC, harmony_clusters %in% epithelial_clusters) # 33694 genes X 8205 cells
perc_bar(blueprint_CRC_Epithelial@meta.data, by = "TumorSite",
         folder = "figures/Qian_2020", suffix = "_epithelial",
         pal = ggpubfigs::friendly_pals$contrast_three[c(2, 3, 1)], width = 3, height = 2)


Epi_markers_Becker_2022_FigS6 <- list(
  tuft = c("DCLK1", "HTR3C", "HTR3E", "B4GALNT4"),       # tuft cells
  immature_goblet = c("RETNLB", "KLK1", "ITLN1", "WFDC2", "CLCA1"), # immature goblet cells
  goblet = c("MUC2", "TFF1", "FCGBP", "TBX10"),     # goblet cells
  immature_enterocytes = c("CA1"),   # immature enterocytes
  enterocytes = c("RAB6B"), # enterocytes
  "BEST4+_enterocytes" = c("BEST4", "CA7", "OTOP2", "OTOP3"), # BEST4+ enterocytes
  enteroendocrine = c("CRYBA2", "SCGN", "FEV"),    # enteroendocrine cells
  stem = c("SOX9", "SMOC2", "RGMB", "LGR5", "ASCL2") # stem cells
  ) %>%
  lapply(as.data.frame) %>%
  data.table::rbindlist(idcol = "cluster") %>%
  `colnames<-`(c("cluster", "gene")) %>%
  mutate(cluster = str_replace(cluster, "_", "\n")) %>%
  mutate(cluster = str_replace(cluster, "entero(?=endo)", "entero-\n")) %>%
  mutate(cluster = factor(cluster, unique(.$cluster)))

p <- DotPlot(blueprint_CRC_Epithelial, Epi_markers_Becker_2022_FigS6$gene,
        cols = c("grey90", "brown2")) + 
  coord_flip() + scale_x_discrete(limits = rev) + 
  theme(axis.text.y = element_text(face = 3),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = NULL, y = NULL)
p$data <- p$data %>% 
  left_join(Epi_markers_Becker_2022_FigS6, by= c("features.plot" = "gene"))
p + facet_grid(cluster ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme(strip.text.y.left = element_text(size = 11, angle = 0, lineheight = 0.8),
        panel.spacing.y = unit(1, "mm"))
ggsave("figures/Qian_2020/Dotplot_Epithelial_markers_Becker_2022.png",
       width = 6, height = 6)

# c("LRIG1", "SOX9", "RGMB", "SMOC2", "ASCL2",
#   "LGR5", "EPHB2", "ALCAM", "CD44", "EPCAM", "PROM1")


# 
Epi_markers_Shen_2022 <- c(
  "GUCA2B", "SLC26A3", "CA1",            # Mature colonocytes
  "KRT20",
  "TFF3", "SPDEF", "SPINK4", "REG4", "AGR2", # Goblet cell
  "CD74",
  "HLA-DPA1", "HLA-DPB1", "HLA-DRA", "HLA-DRB1"
  )
 

DotPlot(blueprint_CRC_Epithelial, Epi_markers_Shen_2022,
        cols = c("grey90", "brown2")) + 
  coord_flip() + scale_x_discrete(limits = rev) + 
  theme(axis.text.y = element_text(face = 3),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = NULL, y = NULL)


blueprint_CRC_Epithelial[["joined"]] <- JoinLayers(blueprint_CRC_Epithelial[["RNA"]])
DefaultAssay(blueprint_CRC_Epithelial) <- "joined"

blueprint_CRC_Epithelial_markers <- FindAllMarkers(
  blueprint_CRC_Epithelial, logfc.threshold = 0.25, only.pos = TRUE, min.pct = 0.25)
# saveRDS(blueprint_CRC_Epithelial_markers, "R_objects/Qian_2020/blueprint_CRC_Epithelial_markers.rds")

blueprint_CRC_Epithelial_markers_top5 <- blueprint_CRC_Epithelial_markers %>%
  group_by(cluster) %>% filter(p_val_adj < 0.01) %>%
  arrange(desc(avg_log2FC)) %>%
  slice(1:5)

blueprint_FigS3p_Epi_markers <- list(
  Epi_CRC = c("MT1E", "MT1G", "ITLN1", "ZG16"),
  C1_OLFM = c("OLFM4", "LEFTY1", "MLEC"),
  C2_CA1 = c("CA1", "FTH1", "FABP1"),
  C3_CEACAM7 = c("CEACAM7", "SLC26A3", "AQP8"),
  C4_SPINK4 = c("SPINK4", "RETNLB"),
  C5_TFF1 = c("TFF1", "IL3RA"),
  C6_BEST4 = c("BEST4", "OTOP2"),
  C7_PLCG2 = c("PLCG2", "RNU12", "RNU11"),
  C8_LY6G6D = c("CCDC170", "DPEP1", "C2") )%>%
  lapply(as.data.frame) %>%
  data.table::rbindlist(idcol = "cluster") %>%
  mutate(cluster = str_replace(cluster, "_", "_\n")) %>%
  `colnames<-`(c("cluster", "gene")) %>%
  mutate(cluster = factor(cluster, unique(.$cluster))) %>%
  mutate(gene = factor(gene, unique(.$gene)))

p <- DotPlot(blueprint_CRC_Epithelial, blueprint_FigS3p_Epi_markers$gene,
             cols = c("grey90", "brown2")) + 
  coord_flip() + scale_x_discrete(limits = rev) + 
  # scale_color_gradientn(colors = brewer.pal(9, "YlOrBr")) +
  theme(axis.text.y = element_text(face = 3),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = NULL, y = NULL)
p$data <- p$data %>% 
  left_join(blueprint_FigS3p_Epi_markers, by= c("features.plot" = "gene"))
p + facet_grid(cluster ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme(strip.text.y.left = element_text(size = 11),
        panel.spacing.y = unit(1, "mm"))
ggsave("figures/Qian_2020/Dotplot_Epithelial_markers_Qian_2020_FigS3p.png",
       width = 6, height = 7)


DotPlot(blueprint_CRC_Epithelial, unique(blueprint_CRC_Epithelial_markers_top5$gene),
        cols = c("grey90", "brown2"), group.by = "harmony_clusters") + 
  coord_flip() + scale_x_discrete(limits = rev) + 
  # scale_color_gradientn(colors = brewer.pal(9, "YlOrBr")) +
  theme(axis.text.y = element_text(face = 3),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = NULL, y = NULL)
ggsave("figures/Qian_2020/Dotplot_blueprint_CRC_Fibroblast_top5Markers.png",
       width = 5, height = 7)


