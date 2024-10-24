blueprint_CRC <- readRDS(file = "R_objects/Qian_2020/blueprint_CRC_harmony.rds")

TableS1 <- readxl::read_xlsx("data/Qian_2020/Qian_2020_clinical_info.xlsx") %>%
  `colnames<-`(str_replace_all(colnames(.), " ", "_") %>% str_replace(
    "\\(years\\)", "years") ) 
  
TableS2 <- readxl::read_excel("data/Qian_2020/Qian_2020_sample_info.xlsx") %>%
  filter(`Cancer type` == "CRC") %>%
  `colnames<-`(str_replace(colnames(.), " ", "_")) 


# Matching Patient IDs
PatientID_mapping <- blueprint_CRC@meta.data %>%
  group_by(orig.ident, PatientNumber, TumorSite) %>%
  summarise(n = n()) %>% ungroup %>%
  left_join(TableS2, by = c("n" = "Cells")) %>%
  select(PatientNumber, PatientNumber_TableS1 = `Patient_number`) %>%
  distinct() 

blueprint_CRC@meta.data <- blueprint_CRC@meta.data %>%
  rownames_to_column("barcode") %>%
  left_join(PatientID_mapping, by = c("PatientNumber")) %>%
  left_join(TableS1, by = c("PatientNumber_TableS1" = "Patient_number")) %>%
  column_to_rownames("barcode")

blueprint_CRC_immune@meta.data <- blueprint_CRC_immune@meta.data %>%
  rownames_to_column("barcode") %>%
  left_join(PatientID_mapping, by = c("PatientNumber")) %>%
  left_join(TableS1, by = c("PatientNumber_TableS1" = "Patient_number")) %>%
  column_to_rownames("barcode")

saveRDS(blueprint_CRC, file = "R_objects/Qian_2020/blueprint_CRC_harmony_updated_meta.rds")
saveRDS(blueprint_CRC_immune, file = "R_objects/Qian_2020/blueprint_CRC_immune_res1_updated_meta.rds")

blueprint_CRC$Molecular_status <- factor(blueprint_CRC$Molecular_status, c("MSS", "MSI-high"))

p <- DimPlot(blueprint_CRC, group.by = "Molecular_status", pt.size = 0.1, 
        cols = c("grey80", "brown3"), label = FALSE)
p$data <- p$data %>% mutate(cluster = blueprint_CRC$harmony_clusters) %>% arrange(Molecular_status)
LabelClusters(p, id = "cluster", size = 4, repel = TRUE, max.overlaps = Inf)
ggsave("figures/Qian_2020/DimPlot_harmony_Molecular_status.png", width = 9, height = 7)


blueprint_CRC_immune$Molecular_status <- factor(blueprint_CRC_immune$Molecular_status, c("MSS", "MSI-high"))
DimPlot(blueprint_CRC_immune, group.by = "Molecular_status", pt.size = 0.1, order = TRUE)


perc_bar(blueprint_CRC@meta.data, by = "Molecular_status",
         folder = "figures/Qian_2020", suffix = "",
         pal = c("grey80", "brown3"), width = 7, height = 2.5)
perc_bar(blueprint_CRC_immune@meta.data, by = "Molecular_status",
         folder = "figures/Qian_2020", suffix = "_immune",
         pal = c("grey80", "brown3"), width = 6, height = 2.5)

p <- DimPlot(blueprint_CRC_immune, group.by = "Molecular_status", pt.size = 0.1, 
             cols = c("grey80", "brown3"), label = FALSE)
p$data <- p$data %>% mutate(cluster = blueprint_CRC_immune$harmony_clusters_immune) %>% arrange(Molecular_status)
LabelClusters(p, id = "cluster", size = 4, repel = TRUE, max.overlaps = Inf)
ggsave("figures/Qian_2020/DimPlot_harmony_immune_Molecular_status.png", width = 7, height = 5)

