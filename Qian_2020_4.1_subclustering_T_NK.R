source("Qian_2020_0_packages_functions.R")
blueprint_CRC_T_NK <- readRDS("R_objects/Qian_2020/blueprint_CRC_T_NK.rds")

# Subclustering T & NK cells --------------------------------------------------
blueprint_CRC_T_NK <- subset(blueprint_CRC_immune, harmony_clusters_immune %in% setdiff(T_NK_clusters, 12) | annotation == "CD8+ T_prolif") # 33694 genes X 1349 cells
blueprint_CRC_T_NK <- NormalizeData(blueprint_CRC_T_NK, normalization.method = "LogNormalize", scale.factor = 10000) # try sctransformation
blueprint_CRC_T_NK <- FindVariableFeatures(blueprint_CRC_T_NK, selection.method = "vst", nfeatures = 2000)
blueprint_CRC_T_NK <- ScaleData(blueprint_CRC_T_NK)
blueprint_CRC_T_NK <- RunPCA(blueprint_CRC_T_NK, npcs = 70)

ElbowPlot(blueprint_CRC_T_NK, 70)

## Perform Harmony integration
blueprint_CRC_T_NK <- IntegrateLayers(
  object = blueprint_CRC_T_NK, method = HarmonyIntegration,
  new.reduction = "harmony",
  verbose = TRUE, npcs = 50
)

blueprint_CRC_T_NK <- RunUMAP(
  blueprint_CRC_T_NK, reduction = "harmony", dims = 1:50,
  reduction.name = "umap.harmony")

blueprint_CRC_T_NK <- FindNeighbors(blueprint_CRC_T_NK, reduction = "harmony", dims = 1:50)
blueprint_CRC_T_NK <- FindClusters(blueprint_CRC_T_NK, resolution = 0.,
                                      cluster.name = "harmony_clusters_T_NK")


saveRDS(blueprint_CRC_T_NK, file = "R_objects/Qian_2020/blueprint_CRC_T_NK.rds")

## Top markers -----------------------------------------------------------------
blueprint_CRC_T_NK[["joined"]] <- JoinLayers(blueprint_CRC_T_NK[["RNA"]])
DefaultAssay(blueprint_CRC_T_NK) <- "joined"

blueprint_CRC_T_NK_markers <- FindAllMarkers(
  blueprint_CRC_T_NK, logfc.threshold = 0.25, only.pos = TRUE, min.pct = 0.25)

saveRDS(blueprint_CRC_T_NK_markers, file = "R_objects/Qian_2020/blueprint_CRC_T_NK_markers.rds")


blueprint_CRC_T_NK_markers_top5 <- blueprint_CRC_T_NK_markers %>%
  mutate(cluster = factor(cluster, levels(blueprint_CRC_T_NK$seurat_clusters))) %>%
  group_by(cluster) %>% arrange(desc(avg_log2FC)) %>% 
  slice(1:5) %>% mutate(gene = factor(gene, unique(.$gene))) 

draw_dotplot(blueprint_CRC_T_NK, blueprint_CRC_T_NK_markers_top5, folder = "Qian_2020/subcluster_T_NK",
             group = "seurat_clusters", suffix = "blueprint_CRC_T_NK_markers_top5",
             fig_height = 9, fig_width = 6, facet = FALSE, pal = brewer.pal(8, "Reds"))

## Visualization ---------------------------------------------------------------

## DimPlot
p <- DimPlot(blueprint_CRC_T_NK, reduction = "umap.harmony", ncol = 2,
        group.by = c("harmony_clusters_T_NK", "annotation",
                     "PatientNumber", "TumorSite"), cols = c25)
p[[1]] <- LabelClusters(p[[1]], id = "harmony_clusters_T_NK", size = 4, repel = TRUE, fontface = 2)
p[[2]] <- p[[2]] + scale_color_manual(values = c(brewer.pal(4, "Reds"), brewer.pal(3, "Blues"), "#6A3D9A"))
p[[4]] <- p[[4]] + scale_color_manual(values = ggpubfigs::friendly_pals$contrast_three[c(2, 3, 1)])
ggsave("figures/Qian_2020/subcluster_T_NK/DimPlot_subcluster_T_NK.png", p, width = 10, height = 7)

## QC metrics & cell-cycle score
blueprint_CRC_T_NK$log10.nCount_RNA <- log10(blueprint_CRC_T_NK$nCount_RNA)
blueprint_CRC_T_NK$log10.nFeature_RNA <- log10(blueprint_CRC_T_NK$nFeature_RNA)

blueprint_CRC_T_NK <- CellCycleScoring(
  blueprint_CRC_T_NK, s.features = cc.genes.updated.2019$s.genes,
  g2m.features = cc.genes.updated.2019$g2m.genes, assay = "joined", set.ident = FALSE)

FeaturePlot(blueprint_CRC_T_NK, c("log10.nCount_RNA", "log10.nFeature_RNA", "percent.mt", "percent.rb", 
                                  "S.Score", "G2M.Score"), label = TRUE) &
  scale_color_gradientn(colors = (brewer.pal(9, "YlOrRd")))
ggsave("figures/Qian_2020/subcluster_T_NK/FeaturePlot_subcluster_T_NK_QC_metrics.png", width = 8, height = 9)


VlnPlot(blueprint_CRC_T_NK, c("log10.nCount_RNA", "log10.nFeature_RNA", "percent.mt", "percent.rb", 
                              "S.Score", "G2M.Score"),  cols = c25, ncol = 1, pt.size = 0.1) & labs(x = NULL)
ggsave("figures/Qian_2020/subcluster_T_NK/VlnPlot_subcluster_T_NK_QC_metrics.png", width = 5, height = 8)

## Percentage barplot

perc_bar(blueprint_CRC_T_NK@meta.data, by = "TumorSite", cluster = "seurat_clusters",
         folder = "figures/Qian_2020/subcluster_T_NK", suffix = "T_NK", position = position,
         width = 5, height = 3, use_default_pal = FALSE, pal = ggpubfigs::friendly_pals$contrast_three[c(2, 3, 1)])
perc_bar(blueprint_CRC_T_NK@meta.data, by = "PatientNumber", cluster = "seurat_clusters",
         folder = "figures/Qian_2020/subcluster_T_NK", suffix = "T_NK", position = position,
         width = 5, height = 3, use_default_pal = FALSE, pal = c25)

## Annotation ------------------------------------------------------------------
blueprint_CRC_T_NK$seurat_clusters <- factor(blueprint_CRC_T_NK$seurat_clusters, c(0, 1, 3, 5, 7, 2, 4, 6, 8, 9, 10))

FeaturePlot(blueprint_CRC_T_NK, c("CD3D", "CD4", "CD8A", "CD8B", "KLRF1", "TRGC1", "TRGC2", "TRDC", "ZBTB16", "SLC4A10"),
            label = TRUE, repel = TRUE, label.size = 3, pt.size = 0.1, order = TRUE)
ggsave("figures/Qian_2020/subcluster_T_NK/FeaturePlot_subcluster_T_NK.png", width = 11, height = 7)

p <- VlnPlot(blueprint_CRC_T_NK, c("CD3D", "CD4", "CD8A", "CD8B", "KLRF1", "TRGC1", "TRGC2", "TRDC", "ZBTB16", "SLC4A10"),
        group.by = "seurat_clusters", cols = c25, ncol = 1) & labs(x = NULL, y = NULL)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + labs(y = colnames(p[[i]]$data)[1], title = NULL) +
    theme(axis.title.y = element_text(angle = 0, face = 4, vjust = 0.5, hjust = 1, margin = margin(r = 10)),
          plot.margin = margin(0, 5, 0, 5), axis.text = element_text(size = 9))
  if(i < length(p)) {
    p[[i]] <- p[[i]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }}

ggsave("figures/Qian_2020/subcluster_T_NK/VlnPlot_subcluster_T_NK.png", p, width = 5, height = 7)


Pelka_2021_annotation <- readxl::read_xlsx("data/signatures/Pelka_2021_TableS2.xlsx", 
                                           sheet = "C. Compositional differences") %>%
  select(Cluster, Annotation) %>% distinct

Pelka_2021_TableS2 <- readxl::read_xlsx("data/signatures/Pelka_2021_TableS2.xlsx",
                                        sheet = "D. Program top genes") %>%
  slice(1:5) %>% select(contains("TNI")) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "label", values_to = "gene") %>%
  mutate(Cluster = paste0("c", str_extract(label, "TNI[0-9]+"))) %>%
  left_join(Pelka_2021_annotation, by = "Cluster") %>%
  mutate(cluster = paste(Cluster, Annotation, sep = "\n"))

draw_dotplot(blueprint_CRC_T_NK, Markers_Zhang_2023_FigS1, 
             fig_width = 7, fig_height = 8, folder = "Qian_2020/subcluster_T_NK",
             suffix = "subcluster_T_NK_Markers_Zhang_2023_FigS1", pal = rev(brewer.pal(11, "RdBu")))


draw_dotplot(blueprint_CRC_T_NK, blueprint_immune_markers$T_NK, 
             fig_width = 7, fig_height = 8, folder = "Qian_2020/subcluster_T_NK",
             suffix = "subcluster_T_NK", pal = rev(brewer.pal(11, "RdBu")))

draw_dotplot(blueprint_CRC_T_NK, Markers_Wang_2021, 
             fig_width = 7, fig_height = 7, folder = "Qian_2020/subcluster_T_NK",
             suffix = "subcluster_T_NK_Markers_Wang_2021", pal = rev(brewer.pal(11, "RdBu")))

draw_dotplot(blueprint_CRC_T_NK, Pelka_2021_TableS2, 
             fig_width = 7, fig_height = 15, folder = "Qian_2020/subcluster_T_NK",
             suffix = "subcluster_T_NK_Pelka_2021_TableS2", pal = rev(brewer.pal(11, "RdBu")))

Pelka_2021_FigS2 <- list(
  "T" = c("CD3D", "CD3E", "CD3G"),
  "CD4+ T" = c("CD4"),
  "CD4+ T IL7R+" = c("IL7R", "LEF1", "SELL", "CCR7", "TCF7"), 
  "CD4+ T IL7R+CCL5+" = c("CCL5"),
  "CD4+ TFH" = c("ICA1", "IL6ST", "PASK", "BCL6", "CXCR5"), 
  "CD4+ CXCL13+" = c("ZBED2", "CXCL13"), 
  "CD4+ Treg" = c("TNFRSF18", "CTLA4", "TNFRSF4", "BATF", "FOXP3", "IL2RA", "IL1RL1", "SOX4"), 
  "CD8+ T" = c("CD8A", "CD8B"),
  "CD8+ T IL7R+" = c("IL7R"), 
  "CD8+ T GZMK+" = c("TCF7", "GZMK"), 
  "CD8+ T CXCL13+" = c("ITGB7", "ITGAE", "RBPJ", "ETV1", "PDCD1", "LAG3", "HAVCR2", "TIGIT", "IFNG", "PRF1", "GNLY", "GZMH", "GZMB", "CXCL13", "TNFSF4", "ENTPD1", "ZBED2"), 
  "gd-like T" = c("TRGC1", "TRGC2", "KLRC2"), 
  "gd-like T PDCD1+" = c("LAG3", "HAVCR2", "PDCD1", "GZMB", "GNLY", "PRF1"), 
  "PLZF+_T" = c("ZBTB16", "TYROBP", "FCER1G"), 
  NK = c("KLF2", "GZMH", "LYAR", "KLRG1", "B3GAT1", "TBX21", "FGFBP2", "CX3CR1", "FCGR3A", "KLRF1", "CMC1", "KLRC1", "NCAM1", "XCL1", "XCL2"), 
  # ILC3 = c("SOX4", "KLRB1", "BCL6", "LST1", "RORC"),
  "HSP+" = c("HSPA1A", "HSPA1B"),
  "IL17+" = c("IL17A", "RORA", "KLRB1", "CD40LG", "TOX"), 
  "prolif" = c("STMN1", "HMGB2", "MKI67")
) %>%
  lapply(as.data.frame) %>%
  data.table::rbindlist(idcol = "cluster") %>%
  `colnames<-`(c("cluster", "gene")) %>%
  mutate(cluster = factor(cluster, unique(.$cluster))) %>%
  mutate(gene = factor(gene, unique(.$gene)))


draw_dotplot(blueprint_CRC_T_NK, Pelka_2021_FigS2, 
             fig_width = 8, fig_height = 14, folder = "Qian_2020/subcluster_T_NK",
             suffix = "subcluster_T_NK_Pelka_2021_FigS2", pal = rev(brewer.pal(11, "RdBu")))


Markers_Tcells_Zhang_2018 <- sapply(
  readxl::excel_sheets("data/signatures/Zhang_2018_TableS5.xlsx"),
  function(x) {
    readxl::read_xlsx("data/signatures/Zhang_2018_TableS5.xlsx", sheet = x, skip = 1) %>%
    filter(!is.na(`Gene ID`)) %>% arrange(F.adjusted.) %>% slice(1:5) %>%
      `colnames<-`(str_to_title(colnames(.))) %>%
      select(gene = `Gene Symbol`)}, 
  USE.NAMES = TRUE, simplify = FALSE) %>%
  data.table::rbindlist(idcol = "cluster") %>%
  select(cluster, gene)


test_marker <- list(
  naive = c("TCF7", "CCR7"),
  cytotoxic = c("PDCD1", "TNFRSF9", "CD137"),
  effector_memory = c("GZMK"),
  exhausted = c("PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2"),
  CD4_Th1_like = c("IFNGR1", "FASL25"),
  CD4_TFH = c("CXCR5", "TOX", "SLAMF6"),
  Treg = c("FOXP3")
) %>%
  lapply(as.data.frame) %>%
  data.table::rbindlist(idcol = "cluster") %>%
  `colnames<-`(c("cluster", "gene")) %>%
  mutate(cluster = factor(cluster, unique(.$cluster))) %>%
  mutate(gene = factor(gene, unique(.$gene)))
draw_dotplot(blueprint_CRC_T_NK, test_marker, 
             fig_width = 8, fig_height = 7, folder = "Qian_2020/subcluster_T_NK",
             suffix = "subcluster_T_NK_test", pal = rev(brewer.pal(11, "RdBu")))

Marker_Lee_2020 <- list(
  TCD4_TCD8 = c("CD4", "IL7R", "CD8A", "CD8B"),
  Naive = c("TCF7", "SELL", "LEF1", "CCR7"),
  Exhausted = c("LAG3", "TIGIT", "PDCD1", "HAVCR2", "CTLA4"),
  Cytotoxic = c("IL2", "GZMA", "GNLY", "PRF1", "GZMB", "GZMK", "IFNG", "NKG7"),
  Co_stimulatory = c("CD28", "TNFRSF14", "ICOS", "TNFRSF9"),
  Regulatory = c("IL2RA", "FOXP3", "IKZF2", "IL4R", "TGFB1", "TGFB3", "TGFBI", "TGFBR1"),
  gdT = c("TRGC1", "TRGC2", "TRDC"),
  Th17 = c("IL17A", "IL17F", "IL22", "CCR6", "KLRB1", "RORA"),
  Th1_like = c("STAT4", "IL12RB2", "IFNG"),
  Th2_like = c("GATA3", "STAT6", "IL4"),
  Tfh = c("MAF", "CXCL13", "CXCR5", "PDCD1"),
  NK = c("XCL1", "FCGR3A", "KLRD1", "KLRF1"),
  Proliferation = c("MKI67", "PCNA")
) %>%
  lapply(as.data.frame) %>%
  data.table::rbindlist(idcol = "cluster") %>%
  `colnames<-`(c("cluster", "gene")) %>%
  mutate(cluster = factor(cluster, unique(.$cluster))) %>%
  mutate(gene = factor(gene, unique(.$gene)))

draw_dotplot(blueprint_CRC_T_NK, Marker_Lee_2020, folder = "Qian_2020/subcluster_T_NK",
             group = "seurat_clusters", suffix = "blueprint_CRC_T_NK_markers_Lee_2020",
             fig_height = 10, fig_width = 7, facet = TRUE, pal = rev(brewer.pal(11, "RdBu")))


