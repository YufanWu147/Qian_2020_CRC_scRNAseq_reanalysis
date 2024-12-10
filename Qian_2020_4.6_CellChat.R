source("Qian_2020_0_packages_functions.R")
library(CellChat)
library(ggrepel)
options(future.globals.maxSize = 8000 * 1024^2)

blueprint_CRC_cleaned <- readRDS("R_objects/Qian_2020/blueprint_CRC_cleaned.rds")
blueprint_CRC_cleaned <- JoinLayers(blueprint_CRC_cleaned)
blueprint_CRC_cleaned$samples <- blueprint_CRC_cleaned$orig.ident


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
    group.l3 %in% c("Macrophage", "Monocyte", "DC", "Mast", "Other Myeloid") ~ "Myeloid",
    TRUE ~ group.l3)) %>%
  mutate(group.l1 = ifelse(group.l2 %in% c("T/NK", "B/PlasmaB", "Myeloid"), "Immune", group.l2)) %>%
  mutate(group.l1 = factor(group.l1, c("Immune", "Fibroblast", "Epithelial", "Endothelial", "Glial")))  %>%
  mutate(group.l2 = factor(group.l2, c("T/NK", "B/PlasmaB", "Myeloid", "Fibroblast", "Epithelial", "Endothelial", "Glial"))) %>%
  mutate(group.l3 = factor(group.l3, c("TCD4", "TCD8", "gdT", "NK", "B", "PlasmaB", "Macrophage", "Monocyte", "DC", "Mast", "Other Myeloid", "Fibroblast",  "Epithelial", "Endothelial", "Glial")))


# CellChat -----------------------------------------------------------------------
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)


blueprint_CRC_cleaned_cellChat <- sapply(c("C", "B", "N"), function(site) {
  cellchat <- createCellChat(object = subset(blueprint_CRC_cleaned, TumorSite == site),
                             group.by = "annotation", assay = "RNA")
  # set the used database in the object
  cellchat@DB <- CellChatDB.use 
  
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  return(cellchat)
}, simplify = FALSE, USE.NAMES = TRUE)



blueprint_CRC_cleaned_cellChat <- lapply(blueprint_CRC_cleaned_cellChat, function(x) {
  x@meta <- x@meta %>% 
    rownames_to_column("barcodes") %>%
    # select(-group.l1, -group.l2, -group.l3) %>%
    left_join(annotation_group, by = "annotation") %>%
    column_to_rownames("barcodes")
  return(x)
})

saveRDS(blueprint_CRC_cleaned_cellChat, file = "R_objects/Qian_2020/blueprint_CRC_cleaned_cellChat.rds")

blueprint_CRC_cleaned_cellChat <- readRDS("R_objects/Qian_2020/blueprint_CRC_cleaned_cellChat.rds")
blueprint_CRC_cleaned_cellChat_merged <- mergeCellChat(blueprint_CRC_cleaned_cellChat, add.names = names(blueprint_CRC_cleaned_cellChat))


# Part I: Identify altered interactions and cell populations -----------------------------------------------
## 1. Compare the total number of interactions and interaction strength ------------------------------------
compareInteractions(blueprint_CRC_cleaned_cellChat_merged, group = c("C", "B", "N"), show.legend = F) +
compareInteractions(blueprint_CRC_cleaned_cellChat_merged, group = c("C", "B", "N"), show.legend = F, measure = "weight")

## 2. Compare the number of interactions and interaction strength among different cell populations --------
### (A) Circle plot showing differential number of interactions or interaction strength among different cell populations across two datasets ---------
### Hard to read
# netVisual_diffInteraction(blueprint_CRC_cleaned_cellChat_merged, weight.scale = T, comparison = c(1, 2),
#                           sources.use = levels(blueprint_CRC_cleaned$annotation)[1:9])
# netVisual_diffInteraction(blueprint_CRC_cleaned_cellChat_merged, weight.scale = T, measure = "weight")

## Compare the number of interactions and interaction strength among different cell populations
### (B) Heatmap showing differential number of interactions or interaction strength among different cell populations across two datasets
png('figures/Qian_2020/cellchat/htmap_in_out_signal_strength_by_celltype_CvsB.png', 
    , width = 8, height = 8, units = "in", res = 400)
netVisual_heatmap(blueprint_CRC_cleaned_cellChat_merged, measure = "weight", 
                  comparison = c("C", "B"), color.use = annot_pal, 
                  title.name = "Differential Interaction Strength, Core vs. Border")
dev.off()
png('figures/Qian_2020/cellchat/htmap_in_out_signal_strength_by_celltype_CvsN.png', 
    , width = 8, height = 8, units = "in", res = 400)
netVisual_heatmap(blueprint_CRC_cleaned_cellChat_merged, measure = "weight",
                  comparison = c("C", "N"), color.use = annot_pal,
                  title.name = "Differential Interaction Strength, Core vs. Normal")
dev.off()
png('figures/Qian_2020/cellchat/htmap_in_out_signal_strength_by_celltype_BvsN.png', 
    , width = 8, height = 8, units = "in", res = 400)
netVisual_heatmap(blueprint_CRC_cleaned_cellChat_merged, measure = "weight", 
                  comparison = c("B", "N"), color.use = annot_pal,
                  title.name = "Differential Interaction Strength, Border vs. Normal")
dev.off()

### (C) Circle plot showing the number of interactions or interaction strength among different cell populations across multiple datasets -----------------

annot_pal_named <- annot_pal %>% `names<-`(levels(blueprint_CRC_cleaned$annotation))
png("figures/Qian_2020/cellchat/circos_strength.png", width = 15, height = 15, units = "in", res = 400)
par(mfrow = c(1, 3))
netVisual_circle(blueprint_CRC_cleaned_cellChat$C@net$weight,
                 weight.scale = T, label.edge= F,
                 sources.use = "SPP1+ TAM", color.use = annot_pal_named,
                 title.name = "Interaction weights/strength")
netVisual_circle(blueprint_CRC_cleaned_cellChat$B@net$weight,
                 weight.scale = T, label.edge= F,
                 sources.use = "SPP1+ TAM", color.use = annot_pal_named,
                 title.name = "Interaction weights/strength")
netVisual_circle(blueprint_CRC_cleaned_cellChat$N@net$weight,
                 weight.scale = T, label.edge= F,
                 sources.use = "SPP1+ TAM", color.use = annot_pal_named,
                 title.name = "Interaction weights/strength")
dev.off()

## (D) Circle plot showing the differential number of interactions or interaction strength among coarse cell types -------------
## Plots are hard to read
# blueprint_CRC_cleaned_cellChat_coarse <- lapply(blueprint_CRC_cleaned_cellChat, function(x) {mergeInteractions(x, annotation_group$group.l3)})
# blueprint_CRC_cleaned_cellChat_coarse_merged <- mergeCellChat(blueprint_CRC_cleaned_cellChat_coarse, add.names = names(blueprint_CRC_cleaned_cellChat_coarse))
# 
# weight.max <- getMaxWeight(blueprint_CRC_cleaned_cellChat_coarse, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
# names(annot_pal_group.l3) <- levels(annotation_group$group.l3)
# par(mfrow = c(1,3), xpd=TRUE)
# for (i in 1:length(blueprint_CRC_cleaned_cellChat_coarse)) {
#   netVisual_circle(blueprint_CRC_cleaned_cellChat_coarse[[i]]@net$count.merged, weight.scale = T, label.edge= T, 
#                    edge.weight.max = weight.max[3], edge.width.max = 12, color.use = annot_pal_group.l3,
#                    title.name = paste0("Number of interactions - ", names(blueprint_CRC_cleaned_cellChat_coarse)[i]))
# }

## 2. Compare the major sources and targets in a 2D space -----------------------------------------
### (A) Identify cell populations with significant changes in sending or receiving signals --------
num.link <- sapply(blueprint_CRC_cleaned_cellChat, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets

sapply(names(blueprint_CRC_cleaned_cellChat), function(i) {
  p <- netAnalysis_signalingRole_scatter(
    blueprint_CRC_cleaned_cellChat[[i]], dot.size = c(1, 6),
    title = i, weight.MinMax = weight.MinMax, group = annotation_group$group.l1,
    color.use = annot_pal) +
    coord_cartesian(xlim = c(0, 50))
    # scale_x_continuous(limits = c(0, 50)) +
    # scale_y_continuous(limits = c(0, 55))
  p$data <- p$data %>% left_join(annotation_group, by = c("labels" = "annotation")) %>%
    mutate(group = ifelse(group.l1 %in% c("Glial", "Endothelial"), "Other", as.character(group.l1))) %>%
    mutate(group = factor(group, c("Immune", "Fibroblast", "Epithelial", "Other")))
  
  # p <- p + facet_grid(group ~ ., scales = "free_y", space = "free_y")
}, simplify = FALSE, USE.NAMES = TRUE) %>%  
  data.table::rbindlist(idcol = "TumorSite") %>%
  mutate(labels = factor(labels, levels(blueprint_CRC_cleaned$annotation))) %>%
  mutate(TumorSite = factor(TumorSite, c("C", "B", "N"))) %>%
  left_join(blueprint_CRC_cleaned@meta.data[, c("annotation", "annotation_label")] %>% distinct(),
            by = c("labels" = "annotation")) %>%
  ggplot(aes(x = x, y = y, size = Count, color = labels, fill = labels)) +
  geom_point(alpha = 0.7) +
  geom_text_repel(aes(label = labels), size = 2.5, max.overlaps = Inf, fontface = 2,
                  lineheight = 0.8, box.padding = 0.2) +
  scale_color_manual(values = annot_pal, guide = "none") +
  scale_fill_manual(values = annot_pal, guide = "none") +
  scale_y_continuous(breaks = seq(0, 40, 10)) +
  facet_grid(group ~ TumorSite, scales = "free_y", space = "free_y") +
  labs(x = "Outgoing Interaction Strength", y = "Incoming Interaction Strength") +
  cowplot::theme_cowplot(font_size = 12)
  # wrap_plots() +
  # plot_layout(axis_titles = "collect", guides = "collect")
ggsave('figures/Qian_2020/cellchat/scatterplot_in_out_signal_strength_by_celltype.png', 
       width = 15, height = 9)

### (B) Identify the signaling changes of specific cell populations -------------------------
sig_change_CvsB_celltype <- blueprint_CRC_cleaned@meta.data %>%
  filter(TumorSite %in% c("C", "B")) %>%
  count(annotation, TumorSite) %>%
  pivot_wider(names_from = "TumorSite", values_from = "n", values_fill = 0) %>%
  filter(C > 20 & B > 20)

sig_change_CvsB <- sapply(levels(blueprint_CRC_cleaned$annotation), function(celltype) {
  p <- netAnalysis_signalingChanges_scatter(
    blueprint_CRC_cleaned_cellChat_merged, comparison = c(1, 2),
    idents.use = celltype)
  p$data
}, USE.NAMES = TRUE, simplify = FALSE) %>%
  data.table::rbindlist(idcol = "annotation") %>%
  filter(annotation %in% sig_change_CvsB_celltype$annotation) %>%
  mutate(annotation = factor(annotation, levels(blueprint_CRC_cleaned$annotation)))
sig_change_CvsB_labels <- sig_change_CvsB %>%
  group_by(annotation) %>%
  arrange()
sig_change_CvsB %>%
  ggplot(aes(x = outgoing, y = incoming, shape = specificity.out.in, color = specificity)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) + 
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) + 
  geom_point() +
  geom_text_repel(aes(label = labels), show.legend = FALSE, size = 3, max.overlaps = Inf) +
  facet_wrap(~annotation, scales = "free") +
  scale_color_manual(values = c("grey10", "#F8766D", "#00BFC4")) +
  labs(x = "Differential outgoing interaction strength",
       y = "Differential incoming interaction strength",
       title = "Signaling Changes, C vs. B")+
  cowplot::theme_cowplot(font_size = 11) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave('figures/Qian_2020/cellchat/scatterplot_in_out_signal_changes_by_celltype.png', 
       width = 15, height = 9)

# Part II: Identify altered signaling with distinct network architecture and interaction strength ---------------
##  Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2
#> 