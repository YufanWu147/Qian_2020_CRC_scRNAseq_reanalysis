source("Qian_2020_0_packages_functions.R")
blueprint_CRC_Mono_Mac <- readRDS("R_objects/Qian_2020/blueprint_CRC_Mono_Mac_DC.rds")

blueprint_CRC_Mono_Mac <- subset(blueprint_CRC_immune, harmony_clusters_immune %in% c(6, 9, 10, 14#, 19, 21, 22, 24, 25, 13
                                                                                      )# | annotation == "Myeloid_prolif"
                                 ) # 33694 genes X 1349 cells
blueprint_CRC_Mono_Mac <- NormalizeData(blueprint_CRC_Mono_Mac, normalization.method = "LogNormalize", scale.factor = 10000) # try sctransformation
blueprint_CRC_Mono_Mac <- FindVariableFeatures(blueprint_CRC_Mono_Mac, selection.method = "vst", nfeatures = 2000)
blueprint_CRC_Mono_Mac <- ScaleData(blueprint_CRC_Mono_Mac, features = rownames(blueprint_CRC_Mono_Mac))
blueprint_CRC_Mono_Mac <- RunPCA(blueprint_CRC_Mono_Mac, npcs = 100)

ElbowPlot(blueprint_CRC_Mono_Mac, 100)

## Perform Harmony integration
blueprint_CRC_Mono_Mac <- IntegrateLayers(
  object = blueprint_CRC_Mono_Mac, method = HarmonyIntegration,
  new.reduction = "harmony",
  verbose = TRUE, npcs = 50
)

blueprint_CRC_Mono_Mac <- RunUMAP(
  blueprint_CRC_Mono_Mac, reduction = "harmony", dims = 1:50,
  reduction.name = "umap.harmony")

blueprint_CRC_Mono_Mac <- FindNeighbors(blueprint_CRC_Mono_Mac, reduction = "harmony", dims = 1:50)
blueprint_CRC_Mono_Mac <- FindClusters(blueprint_CRC_Mono_Mac, resolution = 0.4,
                                   cluster.name = "harmony_clusters_Mono_Mac")
saveRDS(blueprint_CRC_Mono_Mac, file = "R_objects/Qian_2020/blueprint_CRC_Mono_Mac_DC.rds")

# blueprint_CRC_Mono_Mac$annotation_temp <- as.character(blueprint_CRC_Mono_Mac$annotation) %>%
#   str_replace("LYVE", "LYVE1") %>%
#   factor(str_replace(levels(blueprint_CRC_Mono_Mac$annotation), "LYVE", "LYVE1"))
blueprint_CRC_Mono_Mac[["joined"]] <- JoinLayers(blueprint_CRC_Mono_Mac[["RNA"]])
DefaultAssay(blueprint_CRC_Mono_Mac) <- "joined"

# Visualization --------------------------------------------------------------------
p <- DimPlot(blueprint_CRC_Mono_Mac, reduction = "umap.harmony", ncol = 2,
             group.by = c("harmony_clusters_Mono_Mac", "annotation",
                          "PatientNumber", "TumorSite"), cols = c25)
p[[1]] <- LabelClusters(p[[1]], id = "harmony_clusters_Mono_Mac", size = 4, repel = TRUE, fontface = 2)
p[[2]] <- p[[2]] + scale_color_manual(values = c(brewer.pal(3, "Reds"), c25[-c(2, 9)]))
p[[4]] <- p[[4]] + scale_color_manual(values = ggpubfigs::friendly_pals$contrast_three[c(2, 3, 1)])
ggsave("figures/Qian_2020/subcluster_Mono_Mac/DimPlot_subcluster_Mono_Mac.png", p, width = 10, height = 8)

blueprint_CRC_Mono_Mac <- BuildClusterTree(blueprint_CRC_Mono_Mac)
PlotClusterTree(blueprint_CRC_Mono_Mac)

## QC metrics & cell-cycle score
blueprint_CRC_Mono_Mac$log10.nCount_RNA <- log10(blueprint_CRC_Mono_Mac$nCount_RNA)
blueprint_CRC_Mono_Mac$log10.nFeature_RNA <- log10(blueprint_CRC_Mono_Mac$nFeature_RNA)

blueprint_CRC_Mono_Mac <- CellCycleScoring(
  blueprint_CRC_Mono_Mac, s.features = cc.genes.updated.2019$s.genes,
  g2m.features = cc.genes.updated.2019$g2m.genes, assay = "joined", set.ident = FALSE)

FeaturePlot(blueprint_CRC_Mono_Mac, c("log10.nCount_RNA", "log10.nFeature_RNA", "percent.mt", "percent.rb", 
                                  "S.Score", "G2M.Score"), label = TRUE) &
  scale_color_gradientn(colors = (brewer.pal(9, "YlOrRd")))
ggsave("figures/Qian_2020/subcluster_Mono_Mac/FeaturePlot_subcluster_Mono_Mac_QC_metrics.png", width = 8, height = 9)


VlnPlot(blueprint_CRC_Mono_Mac, c("log10.nCount_RNA", "log10.nFeature_RNA", "percent.mt", "percent.rb", 
                              "S.Score", "G2M.Score"),  cols = c25, ncol = 1, pt.size = 0.1) & labs(x = NULL)
ggsave("figures/Qian_2020/subcluster_Mono_Mac/VlnPlot_subcluster_Mono_Mac_QC_metrics.png", width = 4, height = 9)


for(position in c("fill", "stack")) {
  perc_bar(blueprint_CRC_Mono_Mac@meta.data, by = "TumorSite", cluster = "seurat_clusters",
           folder = "figures/Qian_2020/subcluster_Mono_Mac", suffix = "Mono_Mac", position = position,
           width = 4, height = 2, use_default_pal = FALSE, pal = ggpubfigs::friendly_pals$contrast_three[c(2, 3, 1)])
  perc_bar(blueprint_CRC_Mono_Mac@meta.data, by = "PatientNumber", cluster = "seurat_clusters",
           folder = "figures/Qian_2020/subcluster_Mono_Mac", suffix = "Mono_Mac", position = position,
           width = 4, height = 2, use_default_pal = FALSE, pal = c25)
}


Markers_Mono_Mac = c(  
  "ITGAX", "ITGAM",           # Macrophage/Monocyte/DCs
  "CD68", "CD163", "MRC1", # "ADGRE1", "CSFR1", "MAFB",
  "S100A8", "S100A9", "S100A12",
  "CD14", "FCGR3A", "CD3D"
  #"CD1C", "FCGR3B", "XCR1", "TPSAB1"
)

FeaturePlot(
  blueprint_CRC_Mono_Mac, c(Markers_Mono_Mac), 
  pt.size = 0.1, ncol = 3, label = TRUE, repel = TRUE, label.size = 4, order = FALSE) +
  plot_layout(axis_titles = "collect") & guides(color = guide_colorbar(barwidth = 0.5)) #&
  # scale_color_gradientn(colors = (brewer.pal(9, "YlOrRd")))

ggsave('figures/Qian_2020/subcluster_Mono_Mac/FeaturePlot_Mono_Mac_markers.png',
       width = 8, height = 9)

VlnPlot(blueprint_CRC_Mono_Mac, Markers_Mono_Mac, cols = c25, ncol = 1, pt.size = 0.1) %>% stack_vlnplot
ggsave("figures/Qian_2020/subcluster_Mono_Mac/VlnPlot_subcluster_Mono_Mac.png", width = 4, height = 6)

draw_dotplot(blueprint_CRC_Mono_Mac, data.frame(cluster = "", gene = Markers_Mono_Mac),
             pal = rev(brewer.pal(11, "RdBu")), folder = "Qian_2020/subcluster_Mono_Mac", facet = FALSE,
             suffix = "subcluster_Mono_Mac_marker_scaled", fig_width = 6, fig_height = 3)
draw_dotplot(blueprint_CRC_Mono_Mac, data.frame(cluster = "", gene = Markers_Mono_Mac),
             pal = brewer.pal(8, "YlOrRd"), folder = "Qian_2020/subcluster_Mono_Mac", facet = FALSE, scale = FALSE,
             suffix = "subcluster_Mono_Mac_marker", fig_width = 6, fig_height = 3)

# Annotation ------------------------------------------------------------------
draw_dotplot(blueprint_CRC_Mono_Mac, data.frame(cluster = "", gene =Markers_Chu_2024_FigS1), 
             pal = rev(brewer.pal(11, "RdBu")), folder = "Qian_2020/subcluster_Mono_Mac",
             suffix = "subcluster_Mono_Mac_marker_Chu_2024.png", fig_width = 8, fig_height = 4)

draw_dotplot(blueprint_CRC_Mono_Mac, Markers_Mei_2021, 
             pal = rev(brewer.pal(11, "RdBu")), folder = "Qian_2020/subcluster_Mono_Mac",
             suffix = "subcluster_Mono_Mac_marker_Mei_2021.png", fig_width = 8, fig_height = 13)

draw_dotplot(blueprint_CRC_Mono_Mac, blueprint_immune_markers$Myeloid, 
             pal = rev(brewer.pal(11, "RdBu")), folder = "Qian_2020/subcluster_Mono_Mac",
             suffix = "subcluster_Mono_Mac_blueprint_markers.png", fig_width = 8, fig_height = 13)

draw_dotplot(blueprint_CRC_Mono_Mac, subset(Markers_Qi_2022_macrophage_top5, str_detect(cluster, "Macrophage")), 
             pal = rev(brewer.pal(11, "RdBu")), folder = "Qian_2020/subcluster_Mono_Mac",
             suffix = "subcluster_Mono_Mac_marker_Qi_2022_top5.png", fig_width = 8, fig_height = 4)

draw_dotplot(blueprint_CRC_Mono_Mac, subset(Markers_Zhang_2020_Myeloid_Fig2, str_detect(cluster, "Mono|Macro|TAM")), 
             pal = rev(brewer.pal(11, "RdBu")), folder = "Qian_2020/subcluster_Mono_Mac",
             suffix = "subcluster_Mono_Mac_marker_Zhang_2020_Myeloid_Fig2.png", fig_width = 8, fig_height = 7)

draw_dotplot(blueprint_CRC_Mono_Mac, Markers_Zhang_2023_FigS1, 
             pal = rev(brewer.pal(11, "RdBu")), folder = "Qian_2020/subcluster_Mono_Mac",
             suffix = "subcluster_Mono_Mac_marker_Zhang_2023_FigS1", fig_width = 8, fig_height = 9)

draw_dotplot(blueprint_CRC_Mono_Mac, data.frame(cluster = "", gene = c("SPP1", "C5AR1", "MMP3", "TIMP1", "ADAM8")), 
             pal = rev(brewer.pal(11, "RdBu")), folder = "Qian_2020/subcluster_Mono_Mac",
             suffix = "subcluster_Mono_Mac_test", fig_width = 8, fig_height = 3)

VlnPlot(blueprint_CRC_Mono_Mac, c("SPP1", "C1QC", "MRC1", "LYVE1", "CD14", "FCGR3A"), ncol = 1, cols = c25, pt.size = 0.1, alpha = 0.3) %>% stack_vlnplot()
ggsave("figures/Qian_2020/subcluster_Mono_Mac/VlnPlot_subcluster_Mono_Mac_SPP1.png", width = 4, height = 4)
FeaturePlot(blueprint_CRC_Mono_Mac, c("SPP1", "C1QC", "MRC1", "LYVE1", "CD14", "FCGR3A"), ncol = 3, 
            label = TRUE, repel = TRUE, label.size = 4, order = FALSE, pt.size = 0.1)
ggsave("figures/Qian_2020/subcluster_Mono_Mac/FeaturePlot_subcluster_Mono_Mac_SPP1.png", width = 9, height = 5.5)

VlnPlot(blueprint_CRC_Mono_Mac, c("CD68", "CD163", "MRC1", "CD14", "FCGR3A", "CD1C", "C1QC", "THBS1", 
                                  "ITGAM", "CD209", "CMKLR1", "VCAN", "MKI67", "CD3D"), cols = c25, ncol = 1)  %>%
  stack_vlnplot()
ggsave("figures/Qian_2020/subcluster_Mono_Mac/VlnPlot_subcluster_Mono_Mac_Markers_Qi_2022.png", width = 4, height = 7)

blueprint_CRC_Mono_Mac[["joined"]] <- JoinLayers(blueprint_CRC_Mono_Mac[["RNA"]])
DefaultAssay(blueprint_CRC_Mono_Mac) <- "joined"

blueprint_CRC_Mono_Mac_markers <- FindAllMarkers(
  blueprint_CRC_Mono_Mac, logfc.threshold = 0.25, only.pos = TRUE, min.pct = 0.25)

saveRDS(blueprint_CRC_Mono_Mac_markers, file = "R_objects/Qian_2020/blueprint_CRC_Mono_Mac_markers_new.rds")

blueprint_CRC_Mono_Mac_markers <- readRDS("R_objects/Qian_2020/blueprint_CRC_Mono_Mac_markers_new.rds")
blueprint_CRC_Mono_Mac_markers_top5 <- blueprint_CRC_Mono_Mac_markers %>%
  mutate(cluster = factor(cluster, levels(blueprint_CRC_Mono_Mac$seurat_clusters))) %>%
  group_by(cluster) %>% arrange(desc(avg_log2FC)) %>% 
  slice(1:5) %>% mutate(gene = factor(gene, unique(.$gene))) 

draw_dotplot(blueprint_CRC_Mono_Mac, blueprint_CRC_Mono_Mac_markers_top5, folder = "Qian_2020/subcluster_Mono_Mac",
             group = "seurat_clusters", suffix = "blueprint_CRC_Mono_Mac_markers_top5",
             fig_height = 6, fig_width = 6, facet = FALSE, pal = rev(brewer.pal(11, "RdBu")))


Chen_2024_signatures_myeloid <- subset(Chen_2024_signatures, cluster %in% unique(Chen_2024_signatures$cluster)[49:65])  %>%
  split(f = .$cluster) %>% lapply(pull, "gene")
Chen_2024_signatures_myeloid_top10 <- subset(Chen_2024_signatures, cluster %in% unique(Chen_2024_signatures$cluster)[49:65])  %>%
  group_by(gene) %>% slice(1) %>%
  group_by(cluster) %>% arrange(desc(avg_log2FC)) %>% slice(1:10)
blueprint_CRC_Mono_Mac <- AddModuleScore(
  blueprint_CRC_Mono_Mac, features = Chen_2024_signatures_myeloid,
  name = names(Chen_2024_signatures_myeloid))
for(i in 1:length(Chen_2024_signatures_myeloid)) {
  sig_name <- names(Chen_2024_signatures_myeloid)[i]
  blueprint_CRC_Mono_Mac[[sig_name]] <- blueprint_CRC_Mono_Mac[[paste0(sig_name, i)]]
  blueprint_CRC_Mono_Mac[[paste0(sig_name, i)]] <- NULL
}

VlnPlot(blueprint_CRC_Mono_Mac, features = names(Chen_2024_signatures_myeloid),
        stack = TRUE, flip = TRUE, split.by = "harmony_clusters_Mono_Mac", cols = c25) + NoLegend()

draw_dotplot(blueprint_CRC_Mono_Mac, Chen_2024_signatures_myeloid_top10
             , #%>% group_by(gene) %>% slice(1), 
             pal = rev(brewer.pal(11, "RdBu")), folder = "Qian_2020/subcluster_Mono_Mac",
             suffix = "subcluster_Chen_2024_myeloid_top10", fig_width = 8, fig_height = 15)


draw_dotplot(blueprint_CRC_Mono_Mac, data.frame(cluster = str_extract(names(Chen_2024_signatures_myeloid), "Mono|Mast|Neu|DC|Mph"), 
                                                gene = names(Chen_2024_signatures_myeloid)),
             pal = rev(brewer.pal(11, "RdBu")), folder = "Qian_2020/subcluster_Mono_Mac",
             suffix = "subcluster_Chen_2024_myeloid", fig_width = 8, fig_height = 4, facet = TRUE)

DoHeatmap(blueprint_CRC_Mono_Mac, features = Markers_Mono_Mac, group.colors = c25,
          size = 4) 

Markers_BloodCell <- read.csv("data/signatures/BloodCell.csv") %>%
  as.data.frame %>%
  `colnames<-`(c("cluster", "gene")) %>%
  mutate(cluster = str_remove(cluster, " Markers"))

draw_dotplot(blueprint_CRC_Mono_Mac, Markers_BloodCell,
             pal = rev(brewer.pal(11, "RdBu")), folder = "Qian_2020/subcluster_Mono_Mac", 
             suffix = "subcluster_Mono_Mac_marker_BloodCell", fig_width = 6, fig_height = 3)

Cell_marker_Seq <- readxl::read_xlsx("~/Downloads/Cell_marker_Seq.xlsx")
Cell_marker_Seq_CRC <- Cell_marker_Seq %>%
  filter(species == "Human" & cancer_type == "Colorectal Cancer") %>%
  filter(tissue_type %in% c("Colorectum", "Colon"))
