library(Seurat)
library(sctransform)
library(patchwork)
library(cowplot)
library(Matrix)
library(gdata)
library(reshape)
library(tidyverse)
library(parallel)
library(RColorBrewer)
library(speckle)
library(limma)

c25 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00","black", "gold1", 
         "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70", "khaki2", 
         "maroon", "orchid1","deeppink1","blue1","steelblue4","darkturquoise", "green3", 
         "yellow4", "yellow2", "#FEC44F", "brown4","darkcyan", "darkgoldenrod2", 
         "darkolivegreen4", "cornsilk", "aquamarine1", "darkseagreen1", "honeydew3", "salmon", 
         "slateblue3", "thistle2", "slategray2", "springgreen", "rosybrown1", "plum2", 
         "chocolate", "blanchedalmond", "brown1", "aliceblue", "cadetblue1", "coral1", 
         "firebrick1", "darkslategrey", "darkslateblue", "#745745", "#a01b68", "#ffc372",
         "#fff0f0", "#ebd4d4", "#eeeeee", "#de4463", "#e8ffc1", "#19d3da", "#ee6f57", "#ff9642",
         "#646464", "#0072ce", "#dc9cbf", "#ED8E7C", "#F1ECC3", "#D9DD6B", "#7C83FD", "#FDD2BF",
         "#E98580", "#C6B4CE", "#FFDADA", "#B5EAEA", "#BB8760", "#C9E4C5", "#FAEBE0", "#F7DBF0",
         "#66DE93", "#FFC074", "#B2B8A3", "#F2F4C3", "#D1D9D9", "#A7D0CD", "#114E60", "#28B5B5")


### Functions ---------------------------------------------------------------
perc_bar <- function(meta, by = "orig.ident", cluster = "seurat_clusters",
                     suffix = "", width = 7, height = 3, use_default_pal = FALSE,
                     pal = NULL, folder = "CRC_figures", margin_l = 5, position = "fill",
                     facet = FALSE, facet_var = "TumorSite", legend_ncol = 1, legend.pos = "right") {
  if(use_default_pal) {pal <- c(brewer.pal(3, "Paired")[-3], brewer.pal(8, "Set2"))}
  else {pal <- pal}
  p <- meta %>%
    mutate_at(.vars = by, .funs = as.factor) %>%
    group_by_at(.vars = c(cluster, by)) %>%
    summarise(n = n()) %>%
    ggplot(aes(x = .data[[cluster]], y = n, fill = .data[[by]]))+ 
    geom_bar(position=position, stat="identity", alpha = 1) +
    scale_fill_manual(values = pal, name = NULL)+
    coord_cartesian(expand = FALSE) +
    labs(y = NULL, x = NULL) +
    theme_bw(base_size = 11) +
    guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.7, ncol = legend_ncol)) +
    theme(panel.grid = element_blank(), legend.position = legend.pos,
          plot.margin = margin(5, 5, 5, margin_l),
          axis.text.x = element_text(angle = 45, hjust = 1))
  if(facet) {
    p$data <- p$data %>% left_join(distinct(meta[, c(cluster, facet_var)]))
    p <- p + facet_grid(reformulate(facet_var, "."), scales = "free_x")
  }
  ggsave(paste0(folder, "/percbar_by_", by, "_", suffix,
                ifelse(position == "stack", "_stacked", ""),
                ifelse(facet, "_facet", ""), ".png"), 
         p, width = width, height = height, dpi = 400)
  # ggsave(paste0(folder, "/percbar_by_", by, "_", suffix, ".pdf"), 
  #        width = width, height = height, dpi = 400)
  
}

draw_dotplot <- function(obj, marker_df, fig_width = 5, fig_height = 7, suffix = "", facet = TRUE,
                         pal = c("grey90", "brown2"), scale = TRUE, group = "seurat_clusters", folder = "Qian_2020",
                         panel_spacing = 2) {
  marker_df$gene <- factor(marker_df$gene, unique(marker_df$gene))
  p <- DotPlot(obj, unique(marker_df$gene), scale = scale, group.by = group) + 
    coord_flip() + scale_x_discrete(limits = rev) + 
    scale_color_gradientn(colors = pal) + 
    theme(axis.text.y = element_text(face = 3),
          axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(x = NULL, y = NULL)
  
  p$data <- p$data %>%  
    left_join(marker_df, by= c("features.plot" = "gene")) %>% 
    arrange(cluster)
  if(facet) {
    p <- p + facet_grid(cluster ~ ., scales = "free_y", space = "free_y", switch = "y") +
      theme(strip.text.y.left = element_text(angle = 0, size = 11),
            panel.spacing.y = unit(panel_spacing, "mm"))
  }
  
  ggsave(paste0("figures/", folder, "/Dotplot_", suffix, ".png"), p,
         width = fig_width, height = fig_height)
}

stack_vlnplot <- function(p) {
  for(i in 1:length(p)) {
    p[[i]] <- p[[i]] + labs(y = colnames(p[[i]]$data)[1], title = NULL, x = NULL) +
      theme(axis.title.y = element_text(angle = 0, face = 4, vjust = 0.5, hjust = 1, margin = margin(r = 10)),
            plot.margin = margin(0, 5, 0, 5), axis.text = element_text(size = 9))
    if(i < length(p)) {
      p[[i]] <- p[[i]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }}
  p
}

### Immune clusters ---------------------------------------------------------------
T_NK_clusters <- c(1, 2, 3, 8, 12, 15, 16, 18)
Myeloid_clusters <- c(6, 9, 10, 13, 14, 19, 21, 22, 24, 25)
B_PlasmaB_clusters <- c(0, 7, 11, 23, 4, 5, 17, 20, 12)

### Marker genes ---------------------------------------------------------------
Markers_Chu_2024_FigS1 <- c(
  "PTPRC",                   # Immune
  "MS4A1", "CD79A", "CD79B", # B cells
  "JCHAIN", "MZB1",          # Plasma B (added)
  "CD3D", "CD3E", "CD3G",    # T cells
  "CD8A", "CD8B",             # CD8+ T cells
  "CD4",                      # CD4+ T cells
  "FOXP3", "CTLA4", "IL2RA",  # Tregs
  "KLRF1",                    # NK cells
  "TPSAB1",                   # Mast cells
  "ITGAX", "ITGAM",           # Macrophage/Monocyte/DCs
  "CD68", "ADGRE1", "MRC1"
)

Markers_Zhang_2023_FigS1 <- list(
  immune = c("PTPRC"),
  Tcells = c("CD3D", "CD3E", "CD3G"),                                # T cells
  T_CD4 = c("CD4", "ICOS", "TNFRSF4", "KLRB1", "IL7R", "SPOCK2"),   # CD4+ T cells
  T_CD8 = c("CD8A", "CD8B", "NKG7", "GZMK", "GZMA", "CCL5"),        # CD8+ T cells
  NK = c("FGFBP2", "GNLY", "KLRF1", "GZMH", "KLRD1"),            # NK cells
  Bcells = c("MS4A1", "CD79A", "CD79B", "BANK1", "VPREB3", "SMIM14", "TNFRSF13C"),
  Plasma_B = c("JCHAIN", "MZB1", "DERL3", "TNFRSF17", "PRDX4"),        # Plasma B
  Macrophage = c("C1QA", "C1QB", "C1QC", "APOE", "APOC1"),     # Macrophage
  Monocyte = c("S100A8", "S100A9", "IL1B", "IL1RN", "G0S2"), # Monocyte
  DC = c("LILRA4", "TSPAN13", "PLD4", "SOX4", "IL3RA"), # DC
  Mast = c("TPSAB1", "TPSB2", "CPA3", "HPGDS")             # Mast cells
) %>%
  lapply(as.data.frame) %>%
  data.table::rbindlist(idcol = "cluster") %>%
  `colnames<-`(c("cluster", "gene")) %>%
  mutate(cluster = factor(cluster, unique(.$cluster)))

blueprint_immune_markers <- list(
  B_PlasmaB = list(Bcells = c("MS4A1", "CD19"), # CD27-: naive, CD27+: memory
                   "B_mature-naive" = c("IGHD", "IGHM"),
                   "B memory" = c("CD27"),
                   B_GC_memory = c("CCR7", "GPR183"),
                   B_GC_independent = c("CD38", "TCL1A", "IL4R", "FCER2", "RGS13", "NEIL1", "RFTN1", "CDCA7", "AICDA"),
                   PlasmaB = c("MZB1", "SDC1"),
                   PlasmaB_IgG = c("IGHG1", "IGHG2"),
                   PlasmaB_IgA = c("IGHA1", "IGHA2"),
                   PlasmaB_mature = c("PRDM1", "XBP1"),
                   Other = c("CD3D", "KRT8")),
  Myeloid = list(Macrophage = c("CD68", "MSR1", "MRC1"),
                 Macrophage_M1like = c("IL1B", "CXCL9", "CXCL10", "SOCS3"),
                 "Early Stage_Macrophage CCR2+" = c("CCR2", "PPA1", "GPR183"), 
                 "Early Stage_Macrophage CCL2+" = c("CCL2", "ANKRD28", "EMP1", "MMP19"), 
                 "TAM_CCL18+" = c("CCL18", "GPNMB", "NUPR1", "SEPP1", "STAB1", "CCL13"), 
                 "TAM_MMP9+" = c("MMP9", "MMP7", "CA2", "TM4SF19", "CCL22", "IL1RN", "CHI3L1"), 
                 "TAM_CX3CR1+" = c("CX3CR1", "C3", "FCGBP", "PLD4"), 
                 Perivascular_Macrophage = c("LYVE1", "EGFL7", "CD209", "CH25H"), 
                 Neutrophils = c("FCGR3B", "CXCL1", "CXCL8"),
                 "Monocyte_CD14+" = c("CD14", "SELL", "S100A8", "S100A9", "VCAN"), 
                 "Monocyte_CD16+" = c("FCGR3A", "CDKN1C", "MTSS1", "LYPD2"),
                 cDC1 = c("CLEC9A", "XCR1"),
                 cDC2 = c("CD1C", "CLEC10A", "SIRPA"),
                 migratory_cDC = c("CCR7", "CCL17", "CCL19"),
                 pDC = c("LILRA4", "CXCR3", "IRF7"),
                 "Langerhans_like DC"= c("CD1A", "CD207")),
  T_NK = list(Tcells = c("CD3D", "CD3E", "CD3G"),
              TCD4 = c("CD4"),
              TCD8 = c("CD8A", "CD8B"),
              Tregs = c("FOXP3", "IL2RA", "IKZF2"),
              Naive = c("CCR7", "LEF1", "SELL", "TCF7"),
              Cytotoxic_Effector = c("GNLY", "IFNG", "NKG7", "PRF1", "GZMA", "GZMB", "GZMH", "GZMK"),
              Inhibitory = c("HAVCR2", "LAG3", "PDCD1", "CTLA4", "TIGIT", "BTLA", "KLRC1"),
              Memory_Effector = c("ANXA1", "ANKRD28", "IL7R", "CD69", "CD40LG"),
              NK = c("NCR1", "NCAM1", "TYROBP", "FGFBP2", "KLRD1", "KLRF1", "KLRB1",
                     "CX3CR1", "FCGR3A", "XCL1", "XCL2"))
) %>%
  lapply(function(df) {
    df %>% lapply(as.data.frame) %>%
      data.table::rbindlist(idcol = "cluster") %>%
      `colnames<-`(c("cluster", "gene")) %>%
      mutate(cluster = str_replace(cluster, "B_GC_", "B GC-\n")) %>%
      mutate(cluster = str_replace(cluster, "_", "\n")) %>%
      mutate(cluster = factor(cluster, unique(.$cluster))) %>%
      mutate(gene = factor(gene, unique(.$gene)))
  })

Markers_Wang_2021 <- list(
  Tcells = c("CD3G", "CD8A", "CD4"),
  gdT = c("TRDC", "TRGC1"),
  Cytotxic = c("TNFSF10", "GZMA", "PRF1", "IFNG", "GZMB", "NKG7", "CST7"),
  Naive = c("SELL", "LEF1", "TCF7", "CCR7"),
  Regulatory = c("TGFB1", "IL7", "TGFBR1", "IL2RA", "TGFB3", "IL4R", "TGFB1", "FOXP3"),
  Exhausted = c("TIGIT", "HAVCR2", "LAG3", "KLRC1"),
  Co_stimulatory = c("TNFRSF14", "TNFRSF25", "TNFRSF9", "SLAMF1", "ICOS", "CD226"),
  "HSP+ Tcells" = c("HSPB1", "HSPA1A"),
  Cycling = c("CSK1B"),
  NK = c("NCAM1")
) %>%
  lapply(as.data.frame) %>%
  data.table::rbindlist(idcol = "cluster") %>%
  `colnames<-`(c("cluster", "gene")) %>%
  mutate(cluster = factor(cluster, unique(.$cluster)))

Markers_Mei_2021 <- list(
  "C1QA+ Macrophage" = c("LYZ", "PSAP", "HLA-DRB1", "APOC1", "CD74", "CTSD", "SEPP1", "FTL", "C1QC", "APOE", "C1QB", "C1QA"),
  "CD55+ Macrophage" = c("REL", "FLNA", "JUND", "EREG", "NR4A1", "METRNL", "VCAN", "RSRP1", "FCN1", "LINC00936", "KDM6B", "CD55", "SERPINB1", "BRI3", "SON"),
  "CXCL5+ TAM" = c("G0S2", "C15orf48", "TIMP1", "MALAT1", "HSP90AA1", "DNAJB1", "CCL3", "PLIN2", "HSPA1A", "HSPA1B"),
  "SPP1+ TAM" = c("ABCA1", "FMNL2", "ACP2", "SLCO2B1", "SLAMF8", "IFIT3", "TNS3", "SERPING1", "MMP14", "LHFPL2", "NR1H3", "NRP2"),
  Monocyte = c("UQCRH", "PPDPF", "COX4", "ATP6V1G1", "IGLL5", "IGJ", "PPIB", "YWHAZ", "XSWI"),
  Neutrophils = c("IL8", "MNDA", "SOD2", "HBA2", "FTH1", "IFITM2", "S100A9", "S100A8"), 
  pDC = c("NUB1", "BIRC3", "EEF1A1", "ACTB", "IDO1", "CRIP1", "TMSB10", "TXN", "MARCKSL1", "CCL19", "FSCN1"),
  cDC = c("HLA-DPB1", "HLA-DPA1", "DNASE1L3", "LTB", "LGALS2", "HLA-DQB1", "S100B")
) %>%
  lapply(as.data.frame) %>%
  data.table::rbindlist(idcol = "cluster") %>%
  `colnames<-`(c("cluster", "gene")) %>%
  mutate(cluster = factor(cluster, unique(.$cluster)))


Markers_Zhang_2020_Myeloid_TableS2 <- readxl::read_xlsx("data/signatures/Zhang_2020_TableS2_Myeloid_Markers.xlsx", skip = 1) %>%
  select(gene = `Gene symbol`, cluster = `Myeloid cluster`) %>%
  filter(!is.na(gene) & !is.na(cluster) ) %>%
  group_by(cluster) %>% slice(1:5) %>%
  group_by(gene) %>% slice(1) %>% ungroup

Markers_Zhang_2020_Myeloid_Fig2 <- list(
  "hM01_Mast-TPSAB1" = c("TPSAB1", "CPA3", "TPSB2", "TPSD1", "HDC", "SLC18A2", "GATA2", "MS4A2", "KIT", "IL1RL1", "HPGDS", "VWA5A"),
  "hM02_pDC-LILRA4" = c("GZMB", "TCF4", "JCHAIN", "MZB1", "IRF7", "IL3RA", "SLC15A4", "IRF4", "SPIB", "SERPINF1", "ITM2C", "PLAC8", "IRF8", "LILRB4"),
  "hM03_cDC2-CD1C" = c("AXL", "FCER1A", "CD1C", "FCGR2B", "CD1E", "CD1A", "PKIB", "CLEC10A"),
  "hM04_cDC1-BATF3" = c("CST7", "DAPP1", "FLT3", "C1orf54", "CCR7", "IDO1", "CD40", "BATF3", "XCR1", "ANPEP"),
  "hM05_Mono-CD14" = c("S100A8", "VCAN", "S100A9", "CD36", "S100A12", "CD14"),
  "hM06_Mono-CD16" = c("TCF7L2", "IFITM2", "MTSS1", "FCGR3A", "SERPINA1", "SIGLEC10", "RHOC", "CX3CR1", "MS4A7", "LILRB1", "LILRA1"),
  "hM08_Macro-NLRP3" = c("NLRP3", "EREG", "ITGAX", "IL1B", "VEGFA", "CCL3"),
  "hM12_TAM-C1QC" = c("APOE", "C1QB", "C1QA", "C1QC", "CD81", "TREM2", "ACP5", "APOC1", "PRDM1", "CSF1R", "SLAMF8"),
  "hM13_TAM-SPP1" = c("CSTB", "SPP1", "FN1", "ABL2", "PPARG", "SDC4", "RGCC", "ADM", "MARCO")) %>%
  lapply(as.data.frame) %>%
  data.table::rbindlist(idcol = "cluster") %>%
  `colnames<-`(c("cluster", "gene")) 


### Qi et al. 2022
Markers_Qi_2022 <- sapply(readxl::excel_sheets("data/signatures/Qi_2022_Suppl_Data_2.xlsx"),
                          function(sheet) {readxl::read_xlsx("data/signatures/Qi_2022_Suppl_Data_2.xlsx", sheet = sheet)},
                          USE.NAMES = TRUE, simplify = FALSE) 

Markers_Qi_2022_Tcell_top5 <- Markers_Qi_2022$`T Cells` %>%
  arrange(desc(avg_logFC)) %>%
  group_by(cluster) %>%
  slice(1:5) %>%
  # mutate(cluster = str_replace(cluster, " Memory", "\nMemory")) %>%
  mutate(cluster = factor(cluster, c(
    "Naive CD4+ T", "Effector CD4+T", "Memory CD4+ T", "Treg", "TFH",
    "Effector CD8+ T", "Exhausted CD8+ T", "gdT", "MAIT", "NK"
  )))

Markers_Qi_2022_Bcell_top5 <- Markers_Qi_2022$`B Cells` %>%
  arrange(desc(avg_logFC)) %>%
  group_by(cluster) %>%
  slice(1:5) %>%
  mutate(cluster = str_replace(cluster, " Memory", "\nMemory")) %>%
  mutate(cluster = factor(cluster, unique(.$cluster)))

Markers_Qi_2022_macrophage <- Markers_Qi_2022$`Myeloid Cells` %>%
  mutate(cluster = ifelse(cluster == "VCAN+ Monocyte", "VCAN+ Macrophage", cluster)) %>%
  mutate(cluster = ifelse(cluster == "MARCO+ Macrophage", "SPP1+ Macrophage", cluster)) %>%
  mutate(cluster = str_replace(cluster, " (?=Macrophage|Myeloid)", "\n")) 

Markers_Qi_2022_macrophage_top5 <- Markers_Qi_2022_macrophage %>%
  # filter(p_val_adj < 0.01) %>%
  filter(pct.1 > 0.25 | pct.2 > 0.25) %>%
  arrange(desc(avg_logFC)) %>%
  group_by(cluster) %>% slice(1:10) %>% 
  select(cluster, gene) %>%
  rbind(list("Activated DC" = c("CCR7", "FSCN1"),
             cDC1 = c("XCR1", "CLEC9A"),
             cDC2 = c("FCER1A", "CD1C"),
             "C1QC+ MRC1-\nMacrophage" = c("C1QC", "MRC1"),
             "THBS1+\nMacrophage" = "THBS1",
             "SPP1+\nMacrophage" = c("MARCO", "SPP1"),
             Neutrophil = "CSF3R") %>%
          lapply(as.data.frame) %>%
          data.table::rbindlist(idcol = "cluster") %>%
          `colnames<-`(c("cluster", "gene")), .) %>%
  distinct() %>%
  group_by(gene) %>% slice(1) %>%
  mutate(cluster = factor(cluster, unique(Markers_Qi_2022_macrophage$cluster)))


Chen_2024_signatures <- readxl::read_xlsx("data/signatures/Chen_2024_TableS2.xlsx", skip = 1) %>%
  select(cluster = `cell subtype`, gene, everything())

Markers_Ji_2024 <- readxl::read_xlsx("data/signatures/Markers_Ji_2024_TableS5.xlsx") %>%
  mutate(`Marker gene` = str_split(`Marker gene`, ", ")) %>%
  unnest(cols = c(`Marker gene`)) %>%
  rename("gene" = "Marker gene", "cluster" = "abbreviation") %>%
  mutate(cluster = factor(cluster, unique(.$cluster))) %>%
  mutate(gene = factor(gene, unique(.$gene))) %>%
  group_by(gene) %>% slice(1) %>% arrange(cluster, gene)
