selected <- "RNA_snn_res.0.5"
# Broad marker-supported labels on the existing clustering.
subtypes <- c(
  "0" = "Macrophage", "1" = "Macrophage", "2" = "Macrophage",
  "3" = "Monocyte", "4" = "Mast/basophil", "5" = "Macrophage",
  "6" = "Macrophage", "7" = "DC", "8" = "Macrophage",
  "9" = "Monocyte", "10" = "Mixed lineage", "11" = "Neutrophil",
  "12" = "Macrophage", "13" = "Macrophage", "14" = "Macrophage",
  "15" = "Epithelial/myeloid", "16" = "DC", "17" = "Mast/basophil"
)
annotation_markers <- list(
  Macrophage = c("CD68", "C1QA", "C1QC", "APOE", "FABP4", "MARCO", "FOLR2"),
  Monocyte = c("FCN1", "VCAN", "CD14", "FCGR3A"),
  DC = c("CD1C", "FCER1A", "CLEC10A", "CLEC4C", "IL3RA", "TCF4"),
  Mast_basophil = c("TPSAB1", "TPSB2", "CPA3", "KIT", "HDC", "IL1RL1"),
  Neutrophil = c("FCGR3B", "CSF3R", "CXCR2", "S100A8")
)
