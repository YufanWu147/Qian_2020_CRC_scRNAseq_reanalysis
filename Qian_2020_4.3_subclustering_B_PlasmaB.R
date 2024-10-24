source("Qian_2020_0_packages_functions.R")
blueprint_CRC_B_PlasmaB <- readRDS("R_objects/Qian_2020/blueprint_CRC_B_PlasmaB.rds")

blueprint_CRC_B_PlasmaB <- subset(blueprint_CRC_immune, harmony_clusters_immune %in% c(0, 7, 11, 23, 4, 17, 20) | annotation %in% c("B_prolif", "PlasmaB_prolif")) # 33694 genes X 1349 cells
blueprint_CRC_B_PlasmaB <- NormalizeData(blueprint_CRC_B_PlasmaB, normalization.method = "LogNormalize", scale.factor = 10000) # try sctransformation
blueprint_CRC_B_PlasmaB <- FindVariableFeatures(blueprint_CRC_B_PlasmaB, selection.method = "vst", nfeatures = 2000)
blueprint_CRC_B_PlasmaB <- ScaleData(blueprint_CRC_B_PlasmaB, features = rownames(blueprint_CRC_B_PlasmaB))
blueprint_CRC_B_PlasmaB <- RunPCA(blueprint_CRC_B_PlasmaB, npcs = 70)

ElbowPlot(blueprint_CRC_B_PlasmaB, 70)

## Perform Harmony integration
blueprint_CRC_B_PlasmaB <- IntegrateLayers(
  object = blueprint_CRC_B_PlasmaB, method = HarmonyIntegration,
  new.reduction = "harmony",
  verbose = TRUE, npcs = 40
)

blueprint_CRC_B_PlasmaB <- RunUMAP(
  blueprint_CRC_B_PlasmaB, reduction = "harmony", dims = 1:40,
  reduction.name = "umap.harmony")

blueprint_CRC_B_PlasmaB <- FindNeighbors(blueprint_CRC_B_PlasmaB, reduction = "harmony", dims = 1:40)
blueprint_CRC_B_PlasmaB <- FindClusters(blueprint_CRC_B_PlasmaB, resolution = 0.4,
                                        cluster.name = "harmony_clusters_B_PlasmaB")
saveRDS(blueprint_CRC_B_PlasmaB, file = "R_objects/Qian_2020/blueprint_CRC_B_PlasmaB.rds")

blueprint_CRC_B_PlasmaB[["joined"]] <- JoinLayers(blueprint_CRC_B_PlasmaB[["RNA"]])
DefaultAssay(blueprint_CRC_B_PlasmaB) <- "joined"

# Visualization --------------------------------------------------------------------
p <- DimPlot(blueprint_CRC_B_PlasmaB, reduction = "umap.harmony", ncol = 2,
             group.by = c("harmony_clusters_B_PlasmaB", "annotation",
                          "PatientNumber", "TumorSite"), cols = c25)
p[[1]] <- LabelClusters(p[[1]], id = "harmony_clusters_B_PlasmaB", size = 4, repel = TRUE, fontface = 2)
p[[2]] <- p[[2]] + scale_color_manual(values = c(brewer.pal(4, "Reds"), c25[-2]))
p[[4]] <- p[[4]] + scale_color_manual(values = ggpubfigs::friendly_pals$contrast_three[c(2, 3, 1)])
ggsave("figures/Qian_2020/subcluster_B_PlasmaB/DimPlot_subcluster_B_PlasmaB.png", p, width = 11, height = 8)

## QC metrics & cell-cycle score
blueprint_CRC_B_PlasmaB$log10.nCount_RNA <- log10(blueprint_CRC_B_PlasmaB$nCount_RNA)
blueprint_CRC_B_PlasmaB$log10.nFeature_RNA <- log10(blueprint_CRC_B_PlasmaB$nFeature_RNA)

blueprint_CRC_B_PlasmaB <- CellCycleScoring(
  blueprint_CRC_B_PlasmaB, s.features = cc.genes.updated.2019$s.genes,
  g2m.features = cc.genes.updated.2019$g2m.genes, assay = "joined", set.ident = FALSE)

FeaturePlot(blueprint_CRC_B_PlasmaB, c("log10.nCount_RNA", "log10.nFeature_RNA", "percent.mt", "percent.rb", 
                                      "S.Score", "G2M.Score"), label = TRUE) &
  scale_color_gradientn(colors = (brewer.pal(9, "YlOrRd")))
ggsave("figures/Qian_2020/subcluster_B_PlasmaB/FeaturePlot_subcluster_B_PlasmaB_QC_metrics.png", width = 8, height = 10)


VlnPlot(blueprint_CRC_B_PlasmaB, c("log10.nCount_RNA", "log10.nFeature_RNA", "percent.mt", "percent.rb", 
                                  "S.Score", "G2M.Score"),  cols = c25, ncol = 1, pt.size = 0.1) & labs(x = NULL)
ggsave("figures/Qian_2020/subcluster_B_PlasmaB/VlnPlot_subcluster_B_PlasmaB_QC_metrics.png", width = 4, height = 9)


for(position in c("fill", "stack")) {
  perc_bar(blueprint_CRC_B_PlasmaB@meta.data, by = "TumorSite", cluster = "seurat_clusters",
           folder = "figures/Qian_2020/subcluster_B_PlasmaB", suffix = "B_PlasmaB", position = position,
           width = 4, height = 2, use_default_pal = FALSE, pal = ggpubfigs::friendly_pals$contrast_three[c(2, 3, 1)])
  perc_bar(blueprint_CRC_B_PlasmaB@meta.data, by = "annotation", cluster = "seurat_clusters", 
           folder = "figures/Qian_2020/subcluster_B_PlasmaB", suffix = "B_PlasmaB", position = position,
           width = 4, height = 2, use_default_pal = FALSE, pal = c(brewer.pal(4, "Reds"), c25[-2]))
  perc_bar(blueprint_CRC_B_PlasmaB@meta.data, by = "PatientNumber", cluster = "seurat_clusters",
           folder = "figures/Qian_2020/subcluster_B_PlasmaB", suffix = "B_PlasmaB", position = position,
           width = 4, height = 2, use_default_pal = FALSE, pal = c25)
  perc_bar(blueprint_CRC_B_PlasmaB@meta.data, by = "PatientNumber", cluster = "seurat_clusters",
           folder = "figures/Qian_2020/subcluster_B_PlasmaB", suffix = "B_PlasmaB", position = position, legend.pos = "none",
           width = 9, height = 4.5, facet = TRUE, facet_var = "TumorSite", pal = c25)
}


Markers_B_PlasmaB = c("MS4A1", "MZB1", "CD27", "IGHD", "IGHM", 
                      "CD38", "RGS13", "IGHG1", "IGHA2", "PRDM1")

FeaturePlot(
  blueprint_CRC_B_PlasmaB, c(Markers_B_PlasmaB), pt.size = 0.1, 
  ncol = 3, label = TRUE, repel = TRUE, label.size = 4, order = FALSE) +
  plot_layout(axis_titles = "collect") & guides(color = guide_colorbar(barwidth = 0.5)) &
  scale_color_gradientn(colors = viridis::plasma(11))

ggsave('figures/Qian_2020/subcluster_B_PlasmaB/FeaturePlot_B_PlasmaB_markers.png',
       width = 11, height = 13)

blueprint_CRC_B_PlasmaB$seurat_clusters <- factor(blueprint_CRC_B_PlasmaB$harmony_clusters_B_PlasmaB,
                                                  c(0, 5, 7, 9, 11, 1:4, 10, 6, 8, 12, 13))
VlnPlot(
  blueprint_CRC_B_PlasmaB, Markers_B_PlasmaB, group.by = "seurat_clusters",
  cols = c25[as.numeric(levels(blueprint_CRC_B_PlasmaB$seurat_clusters)) + 1], ncol = 1, pt.size = 0.1, alpha = 0.3) %>%
  stack_vlnplot()
ggsave("figures/Qian_2020/subcluster_B_PlasmaB/VlnPlot_subcluster_B_PlasmaB.png", width = 4, height = 6)

draw_dotplot(blueprint_CRC_B_PlasmaB, data.frame(cluster = "", gene = Markers_B_PlasmaB),
             pal = rev(brewer.pal(11, "RdBu")), folder = "Qian_2020/subcluster_B_PlasmaB", facet = FALSE,
             suffix = "subcluster_B_PlasmaB_marker_scaled", fig_width = 6, fig_height = 3)
draw_dotplot(blueprint_CRC_B_PlasmaB, data.frame(cluster = "", gene = Markers_B_PlasmaB),
             pal = brewer.pal(8, "YlOrRd"), folder = "Qian_2020/subcluster_B_PlasmaB", facet = FALSE, scale = FALSE,
             suffix = "subcluster_B_PlasmaB_marker", fig_width = 6, fig_height = 3)

# Annotation ------------------------------------------------------------------
draw_dotplot(blueprint_CRC_B_PlasmaB, data.frame(cluster = "", gene = Markers_Chu_2024_FigS1), 
             pal = rev(brewer.pal(11, "RdBu")), folder = "Qian_2020/subcluster_B_PlasmaB", facet = FALSE,
             suffix = "subcluster_B_PlasmaB_marker_Chu_2024.png", fig_width = 7, fig_height = 4)

draw_dotplot(blueprint_CRC_B_PlasmaB, blueprint_immune_markers$B_PlasmaB, group = "seurat_clusters",
             pal = rev(brewer.pal(11, "RdBu")), folder = "Qian_2020/subcluster_B_PlasmaB",
             suffix = "subcluster_B_PlasmaB_blueprint_markers.png", fig_width = 8, fig_height = 5.5)

draw_dotplot(blueprint_CRC_B_PlasmaB, Markers_Qi_2022_Bcell_top5, group = "seurat_clusters",
             pal = rev(brewer.pal(11, "RdBu")), folder = "Qian_2020/subcluster_B_PlasmaB",
             suffix = "subcluster_B_PlasmaB_marker_Qi_2022_top5.png", fig_width = 8, fig_height = 6.5)

draw_dotplot(blueprint_CRC_B_PlasmaB, Markers_Zhang_2023_FigS1, 
             pal = rev(brewer.pal(11, "RdBu")), folder = "Qian_2020/subcluster_B_PlasmaB",
             suffix = "subcluster_B_PlasmaB_marker_Zhang_2023_FigS1", fig_width = 8, fig_height = 9)

blueprint_CRC_B_PlasmaB_markers <- FindAllMarkers(
  blueprint_CRC_B_PlasmaB, logfc.threshold = 0.25, only.pos = TRUE, min.pct = 0.25)

saveRDS(blueprint_CRC_B_PlasmaB_markers, file = "R_objects/Qian_2020/blueprint_CRC_B_PlasmaB_markers.rds")

blueprint_CRC_B_PlasmaB_markers <- readRDS("R_objects/Qian_2020/blueprint_CRC_B_PlasmaB_markers.rds")
blueprint_CRC_B_PlasmaB_markers_top5 <- blueprint_CRC_B_PlasmaB_markers %>%
  mutate(cluster = factor(cluster, levels(blueprint_CRC_B_PlasmaB$seurat_clusters))) %>%
  group_by(cluster) %>% arrange(desc(avg_log2FC)) %>% 
  slice(1:5) %>% mutate(gene = factor(gene, unique(.$gene))) 

draw_dotplot(blueprint_CRC_B_PlasmaB, blueprint_CRC_B_PlasmaB_markers_top5, folder = "Qian_2020/subcluster_B_PlasmaB",
             group = "seurat_clusters", suffix = "blueprint_CRC_B_PlasmaB_markers_top5",
             fig_height = 9, fig_width = 6, facet = FALSE, pal = rev(brewer.pal(11, "RdBu")))


Chen_2024_signatures <- readxl::read_xlsx("data/signatures/Chen_2024_TableS2.xlsx", skip = 1) %>%
  select(cluster = `cell subtype`, gene, everything())
Chen_2024_signatures_B_PlasmaB_top10 <- subset(Chen_2024_signatures, cluster %in% unique(Chen_2024_signatures$cluster)[39:48]) %>%
  group_by(cluster) %>% slice(1:20) %>%
  group_by(gene) %>% slice(1)
Chen_2024_signatures_B_PlasmaB <- subset(Chen_2024_signatures, cluster %in% unique(Chen_2024_signatures$cluster)[39:48])  %>%
  split(f = .$cluster) %>% lapply(pull, "gene")
blueprint_CRC_B_PlasmaB <- AddModuleScore(
  blueprint_CRC_B_PlasmaB, features = Chen_2024_signatures_B_PlasmaB,
  name = names(Chen_2024_signatures_B_PlasmaB))
for(i in 1:length(Chen_2024_signatures_B_PlasmaB)) {
  sig_name <- names(Chen_2024_signatures_B_PlasmaB)[i]
  blueprint_CRC_B_PlasmaB[[sig_name]] <- blueprint_CRC_B_PlasmaB[[paste0(sig_name, i)]]
  blueprint_CRC_B_PlasmaB[[paste0(sig_name, i)]] <- NULL
}

VlnPlot(blueprint_CRC_B_PlasmaB, features = names(Chen_2024_signatures_B_PlasmaB),
        stack = TRUE, flip = TRUE, split.by = "harmony_clusters_B_PlasmaB", cols = c25) + NoLegend()

draw_dotplot(blueprint_CRC_B_PlasmaB, data.frame(gene = names(Chen_2024_signatures_B_PlasmaB), cluster = "")
             , pal = rev(brewer.pal(11, "RdBu")), folder = "Qian_2020/subcluster_B_PlasmaB",
             suffix = "subcluster_Chen_2024_B_PlasmaB", fig_width = 8, fig_height = 3, facet = FALSE)

draw_dotplot(blueprint_CRC_B_PlasmaB, Chen_2024_signatures_B_PlasmaB_top10, 
             , pal = rev(brewer.pal(11, "RdBu")), folder = "Qian_2020/subcluster_B_PlasmaB",
             suffix = "subcluster_Chen_2024_B_PlasmaB_top10", fig_width = 8, fig_height = 10)

draw_dotplot(blueprint_CRC_B_PlasmaB, data.frame(cluster = str_extract(names(Chen_2024_signatures_B_PlasmaB), "Mono|Mast|Neu|DC|Mph"), 
                                                gene = names(Chen_2024_signatures_B_PlasmaB)),
             pal = rev(brewer.pal(11, "RdBu")), folder = "Qian_2020/subcluster_B_PlasmaB",
             suffix = "subcluster_Chen_2024_myeloid", fig_width = 8, fig_height = 4, facet = TRUE)

DoHeatmap(blueprint_CRC_B_PlasmaB, features = Markers_B_PlasmaB, group.colors = c25,
          size = 4) 


Markers_B_PlasmaB_Ji_2024 <- list(
  B = c("MS4A1", "CD79A", "CD79B"), 
  Naive_B	= c("TCL1A", "IGHD"),
  Memory_B = c("CD27", "AIM2"),
  PlasmaB = c("JCHAIN", "MZB1"),
  IgG_Plasma = c("IGHG1"),
  IgA_Plasma = c("IGHA1"),
  Proli_Plasma = c("TOP2A"),
  GCBs = c("RGS13", "AICDA", "P2RY8", "IRF8")) %>%
  lapply(as.data.frame) %>%
  data.table::rbindlist(idcol = "cluster") %>%
  `colnames<-`(c("cluster", "gene")) %>%
  mutate(cluster = factor(cluster, unique(.$cluster)))


draw_dotplot(blueprint_CRC_B_PlasmaB, Markers_B_PlasmaB_Ji_2024, pal = rev(brewer.pal(11, "RdBu")), folder = "Qian_2020/subcluster_B_PlasmaB",
             suffix = "Markers_B_PlasmaB_Ji_2024", fig_width = 8, fig_height = 4)

