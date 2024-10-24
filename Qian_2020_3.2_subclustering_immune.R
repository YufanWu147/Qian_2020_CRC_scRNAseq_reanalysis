source("Qian_2020_0_packages_functions.R")

# Subclustering immune clusters -----------------------------------------------------------
blueprint_CRC_immune <- readRDS(file = "R_objects/Qian_2020/blueprint_CRC_immune_res1_updated_meta.rds")
blueprint_CRC_immune_sub5_12 <- readRDS("R_objects/Qian_2020/blueprint_CRC_immune_sub5_12.rds")


immune_clusters <- c(
  8, 16,    # B cells,
  4, 14, 20,# Plasma B
  2, 3, 10, # T cells
  26,       # NK
  21,       # Mast
  6, 13, 19,# Myeloid
  24        # Proliferating
)


blueprint_CRC_immune <- subset(blueprint_CRC, harmony_clusters %in% immune_clusters) # 33694 genes X 21068 cells

blueprint_CRC_immune <- NormalizeData(blueprint_CRC_immune, normalization.method = "LogNormalize", scale.factor = 10000) # try sctransformation
blueprint_CRC_immune <- FindVariableFeatures(blueprint_CRC_immune, selection.method = "vst", nfeatures = 2000)
blueprint_CRC_immune <- ScaleData(blueprint_CRC_immune)
blueprint_CRC_immune <- RunPCA(blueprint_CRC_immune, npcs = 70)

## Perform Harmony integration
blueprint_CRC_immune <- IntegrateLayers(
  object = blueprint_CRC_immune, method = HarmonyIntegration,
  new.reduction = "integrated.harmony",
  verbose = FALSE
)

blueprint_CRC_immune <- RunUMAP(
  blueprint_CRC_immune, reduction = "integrated.harmony", dims = 1:50,
  reduction.name = "umap.harmony")

blueprint_CRC_immune <- FindNeighbors(blueprint_CRC_immune, reduction = "integrated.harmony", dims = 1:50)
blueprint_CRC_immune <- FindClusters(blueprint_CRC_immune, resolution = 1,
                                     cluster.name = "harmony_clusters_immune")

saveRDS(blueprint_CRC_immune, file = "R_objects/Qian_2020/blueprint_CRC_immune_res1.rds")

blueprint_CRC_immune[["joined"]] <- JoinLayers(blueprint_CRC_immune[["RNA"]])
DefaultAssay(blueprint_CRC_immune) <- "joined"

blueprint_CRC_immune_markers <- FindAllMarkers(
  blueprint_CRC_immune, logfc.threshold = 0.25, only.pos = TRUE, min.pct = 0.25)
saveRDS(blueprint_CRC_immune_markers, file = "R_objects/Qian_2020/blueprint_CRC_immune_markers.rds")
saveRDS(blueprint_CRC_immune@assays$joined$counts, 
        file = "R_objects/Qian_2020/blueprint_CRC_immune_counts.rds")


### FindAllMarkers -----------------------------------------------------------------
blueprint_CRC_immune_markers <- readRDS("R_objects/Qian_2020/blueprint_CRC_immune_markers.rds")
blueprint_CRC_immune_markers_top5 <- blueprint_CRC_immune_markers %>%
  group_by(cluster) %>% slice(1:5)
blueprint_CRC_immune_markers_top10 <- blueprint_CRC_immune_markers %>%
  group_by(cluster) %>% slice(1:10)

DotPlot(blueprint_CRC_immune, unique(blueprint_CRC_immune_markers_top5$gene)) +
  coord_flip() + 
  scale_x_discrete(limits = rev) + 
  scale_color_gradientn(colors = rev(brewer.pal(11, "RdBu"))) +
  theme(axis.text.y = element_text(face = 3),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = NULL, y = NULL)
ggsave("figures/Qian_2020/DotPlot_blueprint_CRC_immune_markers_top5.png", 
       width = 8, height = 15)

### QC metrics -----------------------------------------------------------------
blueprint_CRC_immune$log10.nCount_RNA <- log10(blueprint_CRC_immune$nCount_RNA)
blueprint_CRC_immune$log10.nFeature_RNA <- log10(blueprint_CRC_immune$nFeature_RNA)

FeaturePlot(blueprint_CRC_immune, c("log10.nCount_RNA", "log10.nFeature_RNA", "percent.mt", "percent.rb")) &
  scale_color_gradientn(colors = (brewer.pal(11, "YlOrRd")))
ggsave("figures/Qian_2020/FeaturePlot_subcluster_immune_QC_metrics.png", width = 9, height = 7)

VlnPlot(blueprint_CRC_immune, c("log10.nCount_RNA", "log10.nFeature_RNA", "percent.mt", "percent.rb"), 
        cols = c25, ncol = 1, pt.size = 0.1) & labs(x = NULL)
ggsave("figures/Qian_2020/VlnPlot_subcluster_immune_QC_metrics.png", width = 7, height = 8)

## Visualization -----------------------------------------------------
DimPlot(blueprint_CRC_immune, reduction = "umap.harmony", ncol = 2,
        group.by = c("harmony_clusters_immune", "harmony_clusters",
                     "PatientNumber", "TumorSite"), cols = c25,
        label = TRUE, repel = TRUE)
ggsave("figures/Qian_2020/DimPlot_subcluster_immune_test.png", width = 12, height = 9)

perc_bar(blueprint_CRC_immune@meta.data, by = "TumorSite",
         folder = "figures/Qian_2020", suffix = "_immune",
         pal = ggpubfigs::friendly_pals$contrast_three[c(2, 3, 1)], width = 5.4, height = 2)

perc_bar(blueprint_CRC@meta.data, by = "Molecular_status",
         folder = "figures/Qian_2020", 
         pal = ggpubfigs::friendly_pals$contrast_three[c(2, 3, 1)])

perc_bar(blueprint_CRC_immune@meta.data, by = "Molecular_status",
         folder = "figures/Qian_2020", suffix = "_immune", 
         pal = ggpubfigs::friendly_pals$contrast_three[c(2, 3, 1)])


p <- DimPlot(blueprint_CRC_immune, reduction = "umap.harmony", ncol = 1,
        group.by = "CellType", cols = c25)
p$data$cluster <- blueprint_CRC_immune$harmony_clusters_immune
LabelClusters(p, id = "cluster", size = 4, repel = TRUE, max.overlaps = Inf)
ggsave("figures/Qian_2020/DimPlot_subcluster_immune_author_annot.png", width = 5, height = 3.5)


p <- DimPlot(blueprint_CRC_immune, reduction = "umap.harmony", ncol = 1,
             group.by = "Molecular_status", cols = c("grey90", "brown2"))
p$data$cluster <- blueprint_CRC_immune$harmony_clusters_immune
LabelClusters(p, id = "cluster", size = 4, repel = TRUE, max.overlaps = Inf)

# Annotating imunne clusters --------------------------------------------------------------------

draw_dotplot(subset(blueprint_CRC_immune, harmony_clusters_immune %in% Myeloid_clusters), 
             Markers_Mei_2021, pal = rev(brewer.pal(11, "RdBu")),
             suffix = "subcluster_immune_myeloid_marker_Mei_2021.png", fig_width = 7, fig_height = 13)

perc_bar(subset(blueprint_CRC_immune@meta.data, harmony_clusters_immune %in% Myeloid_clusters),
         by = "TumorSite", 
         folder = "figures/Qian_2020", suffix = "_immune_subcluster_myeloid",
         pal = ggpubfigs::friendly_pals$contrast_three[c(2, 3, 1)], width = 3.5, height = 2)

FeaturePlot(blueprint_CRC_immune, c("TPSAB1", "LILRA4", "CD1C", "BATF3", "CD14", "FCGR3A",
                                    "CD163", "IL1B", "NLRP3", "PLTP", "C1QC", "SPP1"), ncol = 4, order = FALSE, pt.size = 0.1) #&
  #scale_color_gradientn(colors = brewer.pal(11, "Reds"))
ggsave('figures/Qian_2020/FeaturePlot_immune_subcluster_myeloid_markers.png', width = 15, height = 10)

FeaturePlot(blueprint_CRC_myeloid, c("TPSAB1", "LILRA4", "CD1C", "BATF3", "CD14", "FCGR3A",
                                    "CD163", "IL1B", "NLRP3", "PLTP", "C1QC", "SPP1"), ncol = 4, order = FALSE, pt.size = 0.1) #&
  #scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral")))
ggsave('figures/Qian_2020/FeaturePlot_immune_subcluster_myeloid_markers.png', width = 12, height = 7)

VlnPlot(subset(blueprint_CRC_immune, harmony_clusters_immune %in% Myeloid_clusters), 
        c("TPSAB1", "LILRA4", "CD1C", "BATF3", "CD14", "FCGR3A",
          "CD163", "IL1B", "NLRP3", "PLTP", "C1QC", "SPP1"), ncol = 4) &
  labs(x = NULL, y = NULL)
ggsave('figures/Qian_2020/VlnPlot_immune_subcluster_myeloid_markers.png', width = 12, height = 5)

VlnPlot(blueprint_CRC_myeloid, 
        c("TPSAB1", "LILRA4", "CD1C", "BATF3", "CD14", "FCGR3A",
          "CD163", "IL1B", "NLRP3", "PLTP", "C1QC", "SPP1"), ncol = 4, cols = c25) &
  labs(x = NULL, y = NULL)
ggsave('figures/Qian_2020/VlnPlot_immune_subcluster_myeloid_markers.png', width = 10, height = 5)

perc_bar(blueprint_CRC_myeloid@meta.data, by = "TumorSite",
         folder = "figures/Qian_2020", suffix = "_myeloid_subcluster", 
         pal = ggpubfigs::friendly_pals$contrast_three[c(2, 3, 1)], width = 3.5, height = 2)

draw_dotplot(blueprint_CRC_myeloid, 
             Markers_Zhang_2020_Myeloid_Fig2%>% filter(str_detect(cluster, "Mono|Macro|TAM")), 
             pal = rev(brewer.pal(11, "RdBu")),
             suffix = "subcluster_immune_macro_mono_marker_Zhang_2020_Fig2", fig_width = 7, fig_height = 7)

## T & NK cells -----------------------------------------------------------------
FeaturePlot(blueprint_CRC_immune, c("CD3D", "CD4", "CD8A", "CD8B", "NCR1", "NCAM1", "S.Score", "G2M.Score"),
            ncol = 4, pt.size = 0.1, order = TRUE)
ggsave("figures/Qian_2020/FeaturePlot_subcluster_Tcell.png", width = 12, height = 5)

VlnPlot(blueprint_CRC_immune, c("CD3D", "NCR1", "CD4", "NCAM1", "CD8A", "S.Score", "CD8B", "G2M.Score"), 
        ncol = 2, cols = c25) & labs(x = NULL, y = NULL)
ggsave("figures/Qian_2020/VlnPlot_subcluster_immune_Bcell_markers_CellCycle_Score.png", 
       width = 12, height = 5)

draw_dotplot(subset(blueprint_CRC_immune, harmony_clusters_immune %in% T_NK_clusters), Markers_Wang_2021,
             suffix = "Markers_Wang_2021", fig_width = 6, pal = rev(brewer.pal(11, "RdBu")))
draw_dotplot(subset(blueprint_CRC_immune, harmony_clusters_immune %in% T_NK_clusters), Pelka_2021_FigS2,
             suffix = "Markers_Pelka_2021", fig_width = 7, pal = rev(brewer.pal(11, "RdBu")), fig_height = 15)
draw_dotplot(subset(blueprint_CRC_immune, harmony_clusters_immune %in% T_NK_clusters), 
             data.frame(cluster= "TFH", gene = factor(c("CD4", "CXCR5", "ICOS", "PDCD1", "BCL6"),
                                                      levels = c("CD4", "CXCR5", "ICOS", "PDCD1", "BCL6"))),
             suffix = "TFH", fig_width = 7, pal = rev(brewer.pal(11, "RdBu")), fig_height = 2)
FeaturePlot(blueprint_CRC_immune, c("CD4", "CXCR5", "ICOS", "PDCD1", "BCL6"), order = TRUE)
VlnPlot(blueprint_CRC_immune, c("CD4", "CXCR5", "ICOS", "PDCD1", "BCL6"), stack = TRUE, flip = TRUE)


a = list(TCD4 = c("CCR7", "FAS", "CD28"),
         Treg = c("FOXP3", "ID2"),
         Tcm_Tfh = c("CXCR5", "PDCD1"), 
         Tcm_Th2 = c("GATA3", "CCR4"),
         Tcm_Th17 = c("RORC", "CCR6"),
         Tem_Th1_17 = c("TBX21", "RORC")) %>%
  lapply(as.data.frame) %>%
  data.table::rbindlist(idcol = "cluster") %>%
  `colnames<-`(c("cluster", "gene")) %>%
  mutate(cluster = factor(cluster, unique(.$cluster)))

draw_dotplot(subset(blueprint_CRC_immune, harmony_clusters_immune %in% T_NK_clusters), a,
             suffix = "test", fig_width = 7, pal = rev(brewer.pal(11, "RdBu")), fig_height = 3)

draw_dotplot(subset(blueprint_CRC_immune, harmony_clusters_immune %in% T_NK_clusters), 
             data.frame(cluster = "TCD4", 
                        gene = c("GNLY", "NKG7", "TCF7", "BCL6", "CXCR5", "CXCL13", "PDCD1",
                                 "HAVCR2", "IFNG", "BHLHE40", "TBX21", "CXCR3", "GZMK", 
                                 "IL17A", "RORC", "IL10", "FOXP3", "ICOS")),
             suffix = "test", fig_width = 7, pal = rev(brewer.pal(11, "RdBu")), fig_height = 5, facet = FALSE)

VlnPlot(
  subset(blueprint_CRC_immune, harmony_clusters_immune %in% T_NK_clusters), 
  features = c("CD4", "GNLY", "NKG7", "TCF7", "BCL6", "CXCR5", "CXCL13", "PDCD1",
           "HAVCR2", "IFNG", "BHLHE40", "TBX21", "CXCR3", "GZMK", 
           "IL17A", "RORC", "IL10", "FOXP3", "ICOS"),
  stack = TRUE, flip = TRUE
)

## B cells ---------------------------------------------------------------------
FeaturePlot(blueprint_CRC_immune, c("MS4A1", "CD79A", "CD79B", "JCHAIN", "MZB1", "S.Score", "G2M.Score"),
            ncol = 4, pt.size = 0.1, order = TRUE)
ggsave("figures/Qian_2020/FeaturePlot_subcluster_Bcell.png", width = 12, height = 5)

VlnPlot(blueprint_CRC_immune, c("MS4A1", "JCHAIN", "CD79A", "MZB1", "CD79B", "S.Score", "G2M.Score"), 
        ncol = 2, cols = c25) & labs(x = NULL, y = NULL)
ggsave("figures/Qian_2020/VlnPlot_subcluster_immune_B_markers.png", 
       width = 10, height = 7)

FeaturePlot(blueprint_CRC_immune, c("MS4A1", "MZB1", "CD27", "IGHD", "IGHM", 
                                    "CD38", "RGS13", "IGHG1", "IGHA2", "PRDM1"),
            ncol = 4, pt.size = 0.1, order = FALSE)
ggsave("figures/Qian_2020/FeaturePlot_subcluster_Bcell_markers.png", width = 12, height = 7)


VlnPlot(subset(blueprint_CRC_immune, harmony_clusters_immune %in% B_PlasmaB_clusters),
        c("MS4A1", "MZB1", "CD27", "IGHD", "IGHM", "CD38", "RGS13", "IGHG1", "IGHA2", "PRDM1", "CD3D", "KLRF1"),
        ncol = 4, cols = c25[B_PlasmaB_clusters + 1]) & labs(x = NULL, y = NULL)
ggsave("figures/Qian_2020/VlnPlot_subcluster_Bcell_markers.png", width = 13, height = 6)

### Cell cycle scoring ---------------------------------------------------------
blueprint_CRC_immune <- CellCycleScoring(
  blueprint_CRC_immune, s.features = cc.genes.updated.2019$s.genes,
  g2m.features = cc.genes.updated.2019$g2m.genes, assay = "joined", set.ident = FALSE)

FeaturePlot(blueprint_CRC_immune, reduction = "umap.harmony", 
            features = c("S.Score", "G2M.Score"), ncol = 1)
ggsave("figures/Qian_2020/Featureplot_subcluster_immune_CellCycle_Score.png", 
       width = 5, height = 7)

VlnPlot(blueprint_CRC_immune, features = c("S.Score", "G2M.Score"), ncol = 1, cols = c25)
ggsave("figures/Qian_2020/VlnPlot_subcluster_immune_CellCycle_Score.png", 
       width = 8, height = 6)

### Myeloid ---------------------------------------------------------
blueprint_CRC_immune_myeloid_markers <- FindAllMarkers(
  subset(blueprint_CRC_immune, harmony_clusters_immune %in% Myeloid_clusters), 
  logfc.threshold = 0.25, only.pos = TRUE, min.pct = 0.25)

blueprint_CRC_immune_macro_mono_markers <- FindAllMarkers(
  subset(blueprint_CRC_immune, harmony_clusters_immune %in% c(6, 9, 10, 14)), 
  logfc.threshold = 0.25, only.pos = TRUE, min.pct = 0.25)

blueprint_CRC_immune_myeloid_markers_top5 <- blueprint_CRC_immune_myeloid_markers %>%
  group_by(cluster) %>% #arrange(desc(avg_log2FC)) %>% 
  slice(1:10) %>% mutate(gene = factor(gene, unique(.$gene)))

blueprint_CRC_immune_macro_mono_markers_top5 <- blueprint_CRC_immune_macro_mono_markers %>%
  group_by(cluster) %>% #arrange(desc(avg_log2FC)) %>% 
  slice(1:10) %>% mutate(gene = factor(gene, unique(.$gene)))


draw_dotplot(subset(blueprint_CRC_immune, harmony_clusters_immune %in% Myeloid_clusters), 
             marker_df = blueprint_CRC_immune_myeloid_markers_top5, facet = FALSE, 
             pal = rev(brewer.pal(11, "RdBu")),
             fig_width = 6, fig_height = 13, suffix = "_immune_myeloid_markers_top5")

draw_dotplot(subset(blueprint_CRC_immune, harmony_clusters_immune %in% c(6, 9, 10, 14)), 
             marker_df = blueprint_CRC_immune_macro_mono_markers_top5, facet = FALSE, 
             pal = rev(brewer.pal(11, "RdBu")),
             fig_width = 6, fig_height = 8, suffix = "_immune_macro_mono_markers_top5")

p <- VlnPlot(subset(blueprint_CRC_immune, harmony_clusters_immune %in% Myeloid_clusters), 
        c("CD68", "MSR1", "MRC1", "CD14", "S100A8", "S100A9", "SELL", 
          "FCGR3A", "CDKN1C", "MTSS1"), pt.size = 0, ncol = 1) & labs(x = NULL) &
  theme(plot.title = element_text(size = 11), 
        axis.text.y = element_text(size = 10),
        plot.margin = margin(0, 5, 0, 5))
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + labs(y = colnames(p[[i]]$data)[1], title = NULL) +
    theme(axis.title.y = element_text(angle = 0, face = 4, vjust = 0.5, hjust = 1))
  if(i < length(p)) {
    p[[i]] <- p[[i]] + theme(axis.text.x = element_blank(), 
                             axis.ticks.x = element_blank())
  }
}
p  
ggsave("figures/Qian_2020/VlnPlot_myeloid_markers.png", width = 5, height= 7)

Markers_myeloid_pan_cancer_Cheng_2021 <- list(
  Macro_C1QC = c("C1QC", "C1QA", "APOE"),
  Macro_LYVE1 = c("LYVE1", "PLTP", "SEPP1"),
  Macro_NLRP3 = c("NLRP3", "EREG", "IL1B"),
  Macro_INHBA = c("INHBA", "IL1RN", "CCL4")
) %>%
  lapply(as.data.frame) %>%
  data.table::rbindlist(idcol = "cluster") %>%
  `colnames<-`(c("cluster", "gene")) %>%
  mutate(cluster = factor(cluster, unique(.$cluster)))


Markers_myeloid <- list(
  Macrophage = c("C1QA", "APOE", "TREM2"),
  Monocyte = c("FN1", "VCAN", "FCGR3A"),
  Patrolling_Monos = c("FCN1", "CDKN1C"),
  Classical_Monos = c("S100A8", "S100A9"),
  NLRP3_TAMs = c("NLRP3", "IL1B"),
  SPP1_TAMs = c("SPP1", "CHI3L1", "MT1G"),
  ISG15_TAMs = c("CXCL10", "CXCL11", "ISG15"),
  CXCL9_TAMs = c("CXCL9", "IL4l1"),
  "FOLR2+LYVE-TRMs" = c("FOLR2", "LYVE1", "MARCO"),
  C3_TRMs = c("ASP", "CX3CR1"),
  Proliferating_TAMs = "MKI67"
) %>%
  lapply(as.data.frame) %>%
  data.table::rbindlist(idcol = "cluster") %>%
  `colnames<-`(c("cluster", "gene")) %>%
  mutate(cluster = factor(cluster, unique(.$cluster)))

a = list(C1QC_TAM = c("C1QA", "C1QB", "ITM2B", "C1QC", "HLA-DMB", "MS4A6A", "CTSC", 
                      "TBXAS1", "TMEM176B", "SYNGR2", "ARHGDIB", "TMEM176A", "UCP2", 
                      "CAPZB", "MAF", "TREM2", "MSR1"),
         SPP1_TAMs = c("SPP1", "PCSK5", "SLC11A1", "VCAN", "SLC25A37", "FLNA", "UPP1", 
                       "BCL6", "AQP9", "TIMP1", "VEGFA", "ADM", "MARCO", "FN1", "IL1RN")) %>%
  lapply(as.data.frame) %>%
  data.table::rbindlist(idcol = "cluster") %>%
  `colnames<-`(c("cluster", "gene")) %>%
  mutate(cluster = factor(cluster, unique(.$cluster)))
draw_dotplot(subset(blueprint_CRC_immune, harmony_clusters_immune %in% Myeloid_clusters), 
             a, pal = rev(brewer.pal(11, "RdBu")),
             suffix = "subcluster_immune_Markers_myeloid.png", fig_width = 7, fig_height = 5)


## Immune marker genes from papers ----------------------------------------------
### Chu et al. 2024
FeaturePlot(blueprint_CRC_immune, Markers_Chu_2024_FigS1, ncol = 5, order = TRUE)
ggsave("figures/Qian_2020/FeaturePlot_subcluster_immune.png", width = 15, height = 12)

### Zhang et al. 2023
draw_dotplot(blueprint_CRC_immune, marker_df = Markers_Zhang_2023_FigS1,
             fig_width = 9, fig_height = 9, suffix = "subcluster_immune_markers_Zhang_2023_FigS1")


### Qian et al. 2020
mapply(function(celltype, fig_height) {
  marker_df <- blueprint_immune_markers[[celltype]]
  for(scale in c(TRUE, FALSE)) {
    if(scale) {pal <- rev(brewer.pal(11, "RdBu"))}
    else { pal <- brewer.pal(8, "Reds")}
    draw_dotplot(blueprint_CRC_immune, marker_df, fig_width = 9, fig_height = fig_height, pal = pal, scale = scale,
                 suffix = paste0("subcluster_blueprint_immune_markers_", celltype, ifelse(scale, "", "_unscaled")))}
  }, 
  celltype = names(blueprint_immune_markers), 
  fig_height = c(5, 9, 8)
)


mapply(function(celltype, fig_width, fig_height, clusters) {
  marker_df <- blueprint_immune_markers[[celltype]]
  for(scale in c(TRUE, FALSE)) {
    if(scale) { pal <- rev(brewer.pal(11, "RdBu"))}
    else { pal <- brewer.pal(8, "Reds")}
    draw_dotplot(subset(blueprint_CRC_immune, harmony_clusters_immune %in% clusters), 
                 marker_df, fig_width = fig_width, fig_height = fig_height, pal = pal, scale = scale,
                 suffix = paste0("subcluster_blueprint_immune_markers_sub_", celltype, ifelse(scale, "", "_unscaled")))
    }
  }, 
  celltype = names(blueprint_immune_markers), 
  clusters = list(B_PlasmaB_clusters, Myeloid_clusters, T_NK_clusters),
  fig_width = c(6, 7, 6),
  fig_height = c(5, 11, 8)
)



draw_dotplot(subset(blueprint_CRC_immune, harmony_clusters_immune %in% Myeloid_clusters), 
             Markers_Qi_2022_macrophage_top5, pal = rev(brewer.pal(8, "RdBu")),
             suffix = "subcluster_immune_macrophage_marker_Qi_2022.png", fig_width = 7, fig_height = 9)

draw_dotplot(subset(blueprint_CRC_immune, harmony_clusters_immune %in% B_PlasmaB_clusters), 
             Markers_Qi_2022_Bcell_top5, pal = rev(brewer.pal(8, "RdBu")),
             suffix = "subcluster_immune_Bcell_marker_Qi_2022.png", fig_width = 7, fig_height = 9)

draw_dotplot(subset(blueprint_CRC_immune, harmony_clusters_immune %in% T_NK_clusters), 
             Markers_Qi_2022_Tcell_top5, pal = brewer.pal(8, "Reds"),#rev(brewer.pal(11, "RdBu")),
             suffix = "subcluster_immune_T_NK_marker_Qi_2022.png", fig_width = 7, fig_height = 12)

### Zhang et al. 2020 ---------------------------------------------------------

draw_dotplot(subset(blueprint_CRC_immune, harmony_clusters_immune %in% Myeloid_clusters), 
             Markers_Zhang_2020_Myeloid_TableS2 #%>% filter(str_detect(cluster, "Mono|Macro|TAM"))
             , pal = rev(brewer.pal(11, "RdBu")),
             suffix = "subcluster_immune_myeloid_marker_Zhang_2020_TableS2", fig_width = 7, fig_height = 13)

draw_dotplot(subset(blueprint_CRC_immune, harmony_clusters_immune %in% Myeloid_clusters), 
             Markers_Zhang_2020_Myeloid_Fig2, #%>% filter(str_detect(cluster, "Mono|Macro|TAM")), 
             pal = rev(brewer.pal(11, "RdBu")),
             suffix = "subcluster_immune_myeloid_marker_Zhang_2020_Fig2", fig_width = 7, fig_height = 13)

# Subclustering cluster 5 & 12 together -------------------------------
## Cluster 5 & 12 are mix of T & B cells
FeatureScatter(subset(blueprint_CRC_immune, harmony_clusters_immune %in% c(5, 12)), "MS4A1", "CD3D") +
  facet_wrap(~colors) + scale_color_discrete(name = "cluster")

blueprint_CRC_immune_sub5_12 <- subset(blueprint_CRC_immune, harmony_clusters_immune %in% c(5, 12)) # 33694 genes X 1349 cells
blueprint_CRC_immune_sub5_12 <- FindVariableFeatures(blueprint_CRC_immune_sub5_12, selection.method = "vst", nfeatures = 2000)
blueprint_CRC_immune_sub5_12 <- ScaleData(blueprint_CRC_immune_sub5_12)
blueprint_CRC_immune_sub5_12 <- RunPCA(blueprint_CRC_immune_sub5_12, npcs = 70)
 
ElbowPlot(blueprint_CRC_immune_sub5_12, 50)

# ## Perform Harmony integration
blueprint_CRC_immune_sub5_12 <- IntegrateLayers(
  object = blueprint_CRC_immune_sub5_12, method = HarmonyIntegration,
  new.reduction = "harmony",
  verbose = TRUE, npcs = 30
)

blueprint_CRC_immune_sub5_12 <- RunUMAP(
  blueprint_CRC_immune_sub5_12, reduction = "harmony", dims = 1:30,
  reduction.name = "umap.harmony")

blueprint_CRC_immune_sub5_12 <- FindNeighbors(blueprint_CRC_immune_sub5_12, reduction = "harmony", dims = 1:30)
blueprint_CRC_immune_sub5_12 <- FindClusters(blueprint_CRC_immune_sub5_12, resolution = 0.2,
                                          cluster.name = "harmony_clusters_sub5_12")
saveRDS(blueprint_CRC_immune_sub5_12, file = "R_objects/Qian_2020/blueprint_CRC_immune_sub5_12.rds")

blueprint_CRC_immune_sub5_12$orig.UMAP1 <- subset(blueprint_CRC_immune, harmony_clusters_immune %in% c(5, 12))@reductions$umap.harmony@cell.embeddings[, 1]
blueprint_CRC_immune_sub5_12$orig.UMAP2 <- subset(blueprint_CRC_immune, harmony_clusters_immune %in% c(5, 12))@reductions$umap.harmony@cell.embeddings[, 2]

blueprint_CRC_immune_sub5_12@meta.data %>%
  ggplot(aes(x = orig.UMAP1, y = orig.UMAP2)) +
  geom_point(aes(color = harmony_clusters_sub5_12), size = 0.1) +
  scale_color_manual(values = c25, name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme_classic()
ggsave("figures/Qian_2020/DimPlot_immune_orig.UMAP_subcluster5_12.png", 
       width = 5, height = 3.5)

blueprint_CRC_immune_sub5_12[["joined"]] <- JoinLayers(blueprint_CRC_immune_sub5_12[["RNA"]])
DefaultAssay(blueprint_CRC_immune_sub5_12) <- "joined"

blueprint_CRC_immune_sub5_12_markers <- FindAllMarkers(
  blueprint_CRC_immune_sub5_12, logfc.threshold = 0.25, only.pos = TRUE, min.pct = 0.25)
blueprint_CRC_immune_sub5_12_markers_top5 <- blueprint_CRC_immune_sub5_12_markers %>%
  group_by(cluster) %>% arrange(desc(avg_log2FC)) %>% slice(1:5)
draw_dotplot(blueprint_CRC_immune_sub5_12, blueprint_CRC_immune_sub5_12_markers_top5,
             suffix = "blueprint_CRC_immune_sub5_12_markers_top5")

blueprint_CRC_immune_sub5_12 <- JoinLayers(blueprint_CRC_immune_sub5_12)
blueprint_CRC_immune_sub5_12 <- CellCycleScoring(
  blueprint_CRC_immune_sub5_12, s.features = cc.genes.updated.2019$s.genes,
  g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = FALSE)


DimPlot(blueprint_CRC_immune_sub5_12, 
        group.by = c("seurat_clusters", "harmony_clusters_immune",
                     "TumorSite", "PatientNumber"), order = TRUE,
        cols = c25, label = TRUE, repel = TRUE, label.size = 4, pt.size = 0.3)
ggsave("figures/Qian_2020/DimPlot_immune_subcluster5_12.png", width = 7, height = 6)

FeaturePlot(blueprint_CRC_immune_sub5_12, c(
  Markers_Chu_2024_FigS1, "S.Score", "G2M.Score", "percent.mt", "percent.rb"),
  ncol = 5, pt.size = 0.1, order = TRUE)
ggsave("figures/Qian_2020/FeauterPlot_immune_subcluster5_12.png", width = 14, height = 10)

blueprint_CRC_immune_sub5_12$log10.nCount_RNA <- log10(blueprint_CRC_immune_sub5_12$nCount_RNA)
blueprint_CRC_immune_sub5_12$log10.nFeature_RNA <- log10(blueprint_CRC_immune_sub5_12$nFeature_RNA)

p <- VlnPlot(blueprint_CRC_immune_sub5_12, c(Markers_Chu_2024_FigS1,
                                           "log10.nCount_RNA", "log10.nFeature_RNA",
                                           "percent.mt", "percent.rb", "S.Score", "G2M.Score"),
             pt.size = 0.1, cols = c25, ncol = 5) & labs(x = NULL, y = NULL)
p$data <- p$data %>% filter(!is.na(ident))
p
ggsave("figures/Qian_2020/VlnPlot_immune_subcluster5_12.png", width = 9, height = 8)


draw_dotplot(blueprint_CRC_immune_sub5_12, Markers_Zhang_2023_FigS1, pal = brewer.pal(8, "Reds"),
             suffix = "subcluster5_Markers_Zhang_2023_FigS1", fig_width = 6, fig_height = 9)

draw_dotplot(blueprint_CRC_immune_sub5_12, Markers_Qi_2022_macrophage_top5, pal = brewer.pal(8, "Reds"),
             suffix = "subcluster5_Markers_Qi_2022_macrophage_top5", fig_width = 6, fig_height = 9)

draw_dotplot(blueprint_CRC_immune_sub5_12, blueprint_immune_markers$B_PlasmaB, pal = brewer.pal(8, "Reds"),
             suffix = "subcluster5_Markers_blueprint_B_PlasmaB", fig_width = 6, fig_height = 5)


VlnPlot(blueprint_CRC_immune_sub5_12,
        c("MS4A1", "MZB1", "CD27", "IGHD", "IGHM", "CD38", "RGS13", "IGHG1", "IGHA2", "PRDM1", "CD3D", "KLRF1"),
        ncol = 4, cols = c25[B_PlasmaB_clusters + 1]) & labs(x = NULL, y = NULL)
ggsave("figures/Qian_2020/VlnPlot_subcluster_sub5_12_Bcell_markers.png", width = 13, height = 6)

# Subset myeloids ---------------------------------------------------------------
blueprint_CRC_myeloid <- subset(blueprint_CRC_immune, harmony_clusters_immune %in% c(6, 9, 10, 14, 19, 21)) # 33694 genes X 1349 cells
blueprint_CRC_myeloid <- FindVariableFeatures(blueprint_CRC_myeloid, selection.method = "vst", nfeatures = 2000)
blueprint_CRC_myeloid <- ScaleData(blueprint_CRC_myeloid)
blueprint_CRC_myeloid <- RunPCA(blueprint_CRC_myeloid, npcs = 70)

ElbowPlot(blueprint_CRC_myeloid, 70)

# ## Perform Harmony integration
blueprint_CRC_myeloid <- IntegrateLayers(
  object = blueprint_CRC_myeloid, method = HarmonyIntegration,
  new.reduction = "harmony",
  verbose = TRUE, npcs = 50
)

blueprint_CRC_myeloid <- RunUMAP(
  blueprint_CRC_myeloid, reduction = "harmony", dims = 1:50,
  reduction.name = "umap.harmony")

blueprint_CRC_myeloid <- FindNeighbors(blueprint_CRC_myeloid, reduction = "harmony", dims = 1:50)
blueprint_CRC_myeloid <- FindClusters(blueprint_CRC_myeloid, resolution = 0.4,
                                      cluster.name = "harmony_clusters_myeloid")

DimPlot(blueprint_CRC_myeloid, group.by = c("harmony_clusters_myeloid", "harmony_clusters_immune"), cols = c25,
        label = TRUE, label.size = 4, repel = TRUE)
ggsave("figures/Qian_2020/DimPlot_subcluster_myeloid.png", width = 8, height = 4)

draw_dotplot(blueprint_CRC_myeloid, 
             Markers_Zhang_2020_Myeloid_Fig2, pal = rev(brewer.pal(11, "RdBu")),
             suffix = "subcluster_myeloid_res.0.4_marker_Zhang_2020_Fig2", fig_width = 7, fig_height = 13)

draw_dotplot(blueprint_CRC_myeloid, 
             blueprint_immune_markers$Myeloid, pal = rev(brewer.pal(11, "RdBu")),
             suffix = "subcluster_myeloid_res.0.4_blueprint_immune_markers", fig_width = 7, fig_height = 13)

VlnPlot(blueprint_CRC_myeloid, 
             c("CD68", "MSR1", "MRC1", "CD14", "S100A8", "S100A9", "SELL", 
               "FCGR3A", "CDKN1C", "MTSS1"), pt.size = 0, stack = TRUE, flip = TRUE,
        split.by = "harmony_clusters_myeloid") & labs(x = NULL, y= NULL) &
  scale_fill_manual(values = c25)
ggsave("figures/Qian_2020/VlnPlot_myeloid_markers.png", width = 5, height= 5)

# Final annotation --------------------------------------------------------------
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

annotation_wrap <- function(x) {
  #str_replace(x, "_prolif", "\n_prolif") %>%
    str_replace(x, "B_mem ", "B_mem\n") %>%
    str_replace("_(?=mature|immature)", "\n_") %>%
    str_replace("(?<=^CD[48]\\+ T)_(?!ex)", "\n_")}

annot <- blueprint_CRC_immune@reductions$umap.harmony@cell.embeddings %>%
  as.data.frame() %>%
  cbind(harmony_clusters_immune = blueprint_CRC_immune$harmony_clusters_immune) %>%
  rownames_to_column("barcodes") %>%
  mutate(annotation = case_when(
    # CD8+ T 
    harmony_clusters_immune == 2 ~ "CD8+ T_GZMK+",
    harmony_clusters_immune == 8 ~ "CD8+ T_cyto_ex",
    
    # CD4+ T 
    harmony_clusters_immune == 1 ~ "CD4+ T_naïve",
    harmony_clusters_immune == 3 ~ "Tregs",
    harmony_clusters_immune == 15 ~ "CD4+ T_mem/eff",
    harmony_clusters_immune == 16 ~ "CD4+ T_CXCL13+",
    
    # NK
    harmony_clusters_immune == 18 ~ "NK",
    
    # B
    harmony_clusters_immune == 17 ~ "B_mature_naïve",
    harmony_clusters_immune == 4 ~ "B_mem GC-dependent",
    harmony_clusters_immune == 20 ~ "B_mem GC-independent",
    
    # Plasma B
    harmony_clusters_immune == 0 ~ "PlasmaB_IgA+_mature_1",
    harmony_clusters_immune == 7 ~ "PlasmaB_IgA+_mature_2",
    harmony_clusters_immune == 11 ~ "PlasmaB_IgG+_mature",
    harmony_clusters_immune == 23 ~ "PlasmaB_IgA+_immature",
    
    # Macrophage/Monocyte/Neutrophil
    harmony_clusters_immune == 6 ~ "SPP1+ TAM",
    harmony_clusters_immune == 10 ~ "LYVE1+ Macrophage",
    harmony_clusters_immune == 14 ~ "C1QC+ TAM",
    harmony_clusters_immune == 9 ~ "Monocyte",
    harmony_clusters_immune == 21 ~ "Neutrophil",
    
    # DC
    harmony_clusters_immune == 24 ~ "cDC1",
    harmony_clusters_immune == 19 ~ "cDC2",
    harmony_clusters_immune == 22 ~ "Migratory cDC",
    harmony_clusters_immune == 25 ~ "pDC",
    
    # Other
    harmony_clusters_immune == 13 ~ "Mast",
    TRUE ~ as.character(harmony_clusters_immune)
  )) %>%
  left_join(sub5_12_annot, by = "barcodes") %>%
  mutate(annotation = ifelse(harmony_clusters_immune %in% c(5, 12), 
                             annot_sub5_12, annotation)) %>%
  mutate(annotation = factor(annotation, levels = c(
    "CD4+ T_naïve", "Tregs", "CD4+ T_mem/eff", "CD4+ T_CXCL13+",
    "CD8+ T_GZMK+", "CD8+ T_cyto_ex", "CD8+ T_prolif",
    "NK",
    "B_mature_naïve", "B_mem GC-dependent", "B_mem GC-independent", "B_prolif",
    "PlasmaB_IgA+_immature", "PlasmaB_IgA+_mature_1", "PlasmaB_IgA+_mature_2", "PlasmaB_IgG+_mature", "PlasmaB_prolif",
    "SPP1+ TAM", "C1QC+ TAM", "LYVE1+ Macrophage", 
    "Monocyte", "Neutrophil", "Myeloid_prolif",
    "cDC1", "cDC2", "Migratory cDC", "pDC",
    "Mast"
  ))) %>%
  mutate(annotation_label = factor(annotation_wrap(annotation), annotation_wrap(levels(annotation))))

blueprint_CRC_immune$annotation <- annot$annotation
blueprint_CRC_immune$annotation_label <- annot$annotation_label

p <- DimPlot(subset(blueprint_CRC_immune, annotation != ""), group.by = "annotation_label", 
             pt.size = 0.4, alpha = 0.2) +
  scale_color_manual(breaks = levels(blueprint_CRC_immune$annotation_label), 
                     labels = levels(blueprint_CRC_immune$annotation), values = c25[-c(6, 29)]) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3, alpha = 0.8)))# + NoLegend()
             #cols = c(brewer.pal(8, "Set1"), brewer.pal(8, "Set2"), brewer.pal(8, "Set3"), c25))
LabelClusters(p, id = "annotation_label", size = 4, max.overlaps = Inf, lineheight = 0.9,
              box = FALSE, fontface = 2, color = "black", box.padding = 0.5, min.segment.length = 1.5)
ggsave("figures/Qian_2020/DimPlot_immune_final_annot.png", width = 11, height = 8)

p <- VlnPlot(subset(blueprint_CRC_immune, annotation != ""), c("S.Score", "G2M.Score"),
             group.by = "annotation", ncol= 1, cols = c25[-c(6, 29)]) & labs(x = NULL, y = NULL) & 
  theme(plot.margin = margin(5, 5, 5, 12))
p[[1]] <- p[[1]] + theme(axis.text.x = element_blank(),
                         axis.ticks.x = element_blank())
p
ggsave("figures/Qian_2020/VlnPlot_immune_final_annot.png", width = 8, height = 5)

draw_dotplot(subset(blueprint_CRC_immune, annotation != ""), 
             Markers_Zhang_2023_FigS1, group = "annotation", suffix = "_immune_Final_Annotation",
             pal = rev(brewer.pal(11, "RdBu")), fig_width = 11, fig_height = 11)
draw_dotplot(subset(blueprint_CRC_immune, annotation != ""), 
             rbind(blueprint_immune_markers$T_NK, blueprint_immune_markers$B_PlasmaB, blueprint_immune_markers$Myeloid), 
             group = "annotation", suffix = "_immune_Final_Annotation_blueprint_immune",
             pal = c("grey90", "brown2"), fig_width = 11, fig_height = 17)

draw_dotplot(subset(blueprint_CRC_immune, annotation != ""), 
             rbind(Markers_Qi_2022_Tcell_top5, Markers_Qi_2022_Bcell_top5, Markers_Qi_2022_macrophage_top5), 
             group = "annotation", suffix = "_immune_Final_Annotation_markers_Qi_2022",
             pal = c("grey90", "brown2"), fig_width = 11, fig_height = 17)

blueprint_CRC_immune_rm_sub5 <- subset(blueprint_CRC_immune, annotation != "")
blueprint_CRC_immune_rm_sub5 <- SetIdent(blueprint_CRC_immune_rm_sub5, value = "annotation")
saveRDS(blueprint_CRC_immune_rm_sub5, "R_objects/Qian_2020/blueprint_CRC_immune_rm_sub5.rds")

blueprint_CRC_immune_rm_sub5_markers <- FindAllMarkers(
  blueprint_CRC_immune_rm_sub5, logfc.threshold = 0.25, only.pos = TRUE, min.pct = 0.25)
saveRDS(blueprint_CRC_immune_rm_sub5_markers, "R_objects/Qian_2020/blueprint_CRC_immune_rm_sub5_markers.rds")
blueprint_CRC_immune_rm_sub5_markers_top5 <- blueprint_CRC_immune_rm_sub5_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>% arrange(desc(avg_log2FC)) %>% 
  slice(1:5) %>%
  group_by(gene) %>% slice(1) %>%
  mutate(gene = factor(gene, unique(.$gene)))

draw_dotplot(blueprint_CRC_immune_rm_sub5, blueprint_CRC_immune_rm_sub5_markers_top5,
             group = "annotation", suffix = "_immune_Final_Annotation_top5markers", panel_spacing = 1,
             pal = rev(brewer.pal(11, "RdBu")), fig_width = 12, fig_height = 23, facet = TRUE)


blueprint_CRC_immune$PatientNumber <- factor(blueprint_CRC_immune$PatientNumber)
blueprint_CRC_immune$Pathological_subtype_simple <- str_extract(blueprint_CRC_immune$Pathological_subtype, "left|right")
for(position in c('stack', 'fill')) {
  perc_bar(subset(blueprint_CRC_immune@meta.data, annotation != ""), by = "TumorSite", cluster = "annotation",
           folder = "figures/Qian_2020", suffix = "_immune_by_annotation", margin_l = 15,
           pal = ggpubfigs::friendly_pals$contrast_three[c(2, 3, 1)], width = 7, height = 3.5, position = position)
  perc_bar(subset(blueprint_CRC_immune@meta.data, annotation != ""), by = "orig.ident", cluster = "annotation",
           folder = "figures/Qian_2020", suffix = "_immune_by_annotation",
           pal = c25, width = 8, height = 4, position = position, margin_l = 15)
  perc_bar(subset(blueprint_CRC_immune@meta.data, annotation != ""), by = "PatientNumber", cluster = "annotation",
           folder = "figures/Qian_2020", suffix = "_immune_by_annotation",
           pal = c25, width = 8, height = 4, position = position, margin_l = 15)
  perc_bar(subset(blueprint_CRC_immune@meta.data, annotation != ""), by = "Gender", cluster = "annotation",
           folder = "figures/Qian_2020", suffix = "_immune_by_annotation",
           pal = c25[2:1], width = 8, height = 4, position = position, margin_l = 15)
  perc_bar(subset(blueprint_CRC_immune@meta.data, annotation != ""), by = "Pathological_subtype_simple", cluster = "annotation",
           folder = "figures/Qian_2020", suffix = "_immune_by_annotation",
           pal = c25, width = 8, height = 4, position = position, margin_l = 15)
  
  perc_bar(subset(blueprint_CRC_immune@meta.data, annotation != ""), by = "annotation", cluster = "PatientNumber",
           folder = "figures/Qian_2020", suffix = "celltype_by_Patient", position = position,
           width = 5, height = 4.5, use_default_pal = FALSE, pal = c25[-c(6, 29)])
  
  perc_bar(subset(blueprint_CRC_immune@meta.data, annotation != ""), by = "annotation", cluster = "orig.ident",
           folder = "figures/Qian_2020", suffix = "celltype_by_Sample", position = position,
           width = 9, height = 4.5, facet = TRUE, facet_var = "TumorSite", pal = c25[-c(6, 29)])
}

ncell_by_sample_annot <- subset(blueprint_CRC_immune@meta.data, annotation != "") %>%
  count(annotation, orig.ident, PatientNumber, TumorSite) %>% ungroup
ncell_by_sample_annot_outliers <- ncell_by_sample_annot %>%
  group_by(annotation, TumorSite) %>%
  mutate(q1 = quantile(n, 0.25), q3 = quantile(n, 0.75)) %>%
  filter(n > q3 + 1.5 * (q3 - q1) | n < q1 - 1.5 * (q3 - q1)) %>%
  full_join(distinct(ncell_by_sample_annot[, c("annotation", "TumorSite")]))
ncell_by_sample_annot %>%
  ggplot(aes(x = annotation, y = n, fill = TumorSite, color = TumorSite)) +
  geom_point(aes(shape = TumorSite), position = position_dodge(width = 0.7)) +
  geom_boxplot(width = 0.7, alpha = 0.5) +
  geom_text(data = ncell_by_sample_annot_outliers, aes(label = PatientNumber), show.legend = FALSE,
            position = position_dodge(width = 0.7), size = 3.5, vjust = -0.3, fontface = "bold") +
  scale_color_manual(values = ggpubfigs::friendly_pals$contrast_three[c(2, 3, 1)]) +
  scale_fill_manual(values = ggpubfigs::friendly_pals$contrast_three[c(2, 3, 1)]) +
  cowplot::theme_cowplot() + labs(x = NULL, y = "# Cells") +
  theme(axis.text.x = element_text(hjust = 1, angle = 45))
ggsave("figures/Qian_2020/boxplot_ncells_by_sample_annot.png", width = 12, height = 4)


VlnPlot(blueprint_CRC, c("MSH2", "MSH6", "MLH1", "PMS2"),
        group.by = "orig.ident", ncol = 1) & labs(x = NULL, y = NULL)

blueprint_CRC$SampleID_TumorSite <- paste(blueprint_CRC$orig.ident, blueprint_CRC$TumorSite, sep = "_")
p <- DotPlot(blueprint_CRC, c("MSH2", "MSH6", "MLH1", "PMS2"),
        group.by = "SampleID_TumorSite") & coord_flip()
p$data <- p$data %>% left_join(distinct(blueprint_CRC@meta.data[, c("PatientNumber", "SampleID_TumorSite")]),
                               by = c("id" = "SampleID_TumorSite"))
p + facet_grid( ~ PatientNumber, scales = "free_x") +
  scale_color_gradientn(colors = brewer.pal(8, "Reds")) + labs(x = NULL, y= NULL) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45))
ggsave("figures/Qian_2020/DotPlot_MLH1_gene_expression.png", width = 8, height = 3.2)

Chen_2024_signatures <- readxl::read_xlsx("data/signatures/Chen_2024_TableS2.xlsx", skip = 1) %>%
  select(cluster = `cell subtype`, gene, everything()) %>%
  group_by(cluster) %>% slice(1:10)
Chen_2024_signatures_TCD8 <- list(
  Tcell = c("CD3D", "CD3E", "CD3G"),
  CD8 = c("CD8A", "CD8B"),
  c16_CD8_Tn_SELL = c("CCR7", "TCF7", "SELL"),
  c18_CD8_Tcm_ANXA1 = c("ANXA1"),
  c20_CD8_Tem_GZMK = c("GZMK", "CCL4L2", "CCL4"),
  c21_CD8_Trm_XCL1 = c("CCR6", "XCL1", "XCL2"),
  c22_CD8_Trm_HSPA1B = c("HSPA1A", "HSPA1B"),
  c23_CD8_Tex_LAYN = c("CXCL13", "HAVCR2", "PDCD1", "LAG3", "ITGAE", "ENTPD1", "LAYN"),
  c27_CD8_MAIT_SLC4A10 = c("SLC4A10", "RORC", "RORA", "KLRB1"),
  c28_CD8_IEL_CD160 = c("CD160"),
  c29_T_MKI67 = c("MKI67", "STMN1"),
  cytotoxic = c("PRF1", "GZMB", "GZMA", "GZMH", "NKG7", "GNLY")
) %>%
  lapply(as.data.frame) %>%
  data.table::rbindlist(idcol = "cluster") %>%
  `colnames<-`(c("cluster", "gene")) %>%
  mutate(cluster = factor(cluster, unique(.$cluster)))

draw_dotplot(subset(blueprint_CRC_immune, harmony_clusters_immune %in% T_NK_clusters), subset(Chen_2024_signatures, cluster %in% unique(Chen_2024_signatures$cluster)[1:14]), 
             pal = rev(brewer.pal(11, "RdBu")), folder = "Qian_2020/subcluster_Mono_Mac",
             suffix = "subcluster_Chen_2024_TCD4", fig_width = 8, fig_height = 12)
draw_dotplot(subset(blueprint_CRC_immune, harmony_clusters_immune %in% T_NK_clusters), 
             Chen_2024_signatures_TCD8, 
             pal = rev(brewer.pal(11, "RdBu")), folder = "Qian_2020/subcluster_Mono_Mac",
             suffix = "subcluster_Chen_2024_TCD8", fig_width = 8, fig_height = 8)

