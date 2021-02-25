## ----setup0, include=FALSE----------------------------------------------------
all_times <- list()  # store the time for each chunk
knitr::knit_hooks$set(time_it = local({
  now <- NULL
  function(before, options) {
    if (before) {
      now <<- Sys.time()
    } else {
      res <- difftime(Sys.time(), now, units = "secs")
      all_times[[options$label]] <<- res
    }
  }
}))
knitr::opts_chunk$set(
  tidy = TRUE,
  tidy.opts = list(width.cutoff = 95),
  message = FALSE,
  warning = FALSE,
  time_it = TRUE
)
celltypes_fast = readRDS("./fls/celltypes_fast_citeseq.rds")
celltypes = readRDS("./fls/celltypes_citeseq.rds")
# pbmc = readRDS("fls/pbmcs_signac_citeseq.rds")

## ----setupSeurat, message = F, eval = F---------------------------------------
#  library(Seurat)

## ----setup, message = F, eval = F---------------------------------------------
#  dir.create("fls")
#  download.file("https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_protein_v3/pbmc_10k_protein_v3_filtered_feature_bc_matrix.h5", destfile = "fls/pbmc_10k_protein_v3_filtered_feature_bc_matrix.h5")

## ----Seurat, message = T, eval = F--------------------------------------------
#  # load dataset
#  E = Read10X_h5(filename = "fls/pbmc_10k_protein_v3_filtered_feature_bc_matrix.h5")
#  pbmc <- CreateSeuratObject(counts = E$`Gene Expression`, project = "pbmc")
#  
#  # run sctransform
#  pbmc <- SCTransform(pbmc)

## ----Seurat 2, message = T, eval = F------------------------------------------
#  # These are now standard steps in the Seurat workflow for visualization and clustering
#  pbmc <- RunPCA(pbmc, verbose = FALSE)
#  pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)
#  pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)

## ----setup2, message = F, eval = F--------------------------------------------
#  install.packages("SignacX")

## ----Signac setup folder, message = T, eval = F-------------------------------
#  # load library
#  library(SignacX)

## ----Signac, message = T, eval = F--------------------------------------------
#  # Run Signac
#  labels <- Signac(pbmc, num.cores = 4)
#  celltypes = GenerateLabels(labels, E = pbmc)

## ----SignacFast, message = T, eval = F----------------------------------------
#  # Run Signac
#  labels_fast <- SignacFast(pbmc)
#  celltypes_fast = GenerateLabels(labels_fast, E = pbmc)

## ---- echo=FALSE, eval = T----------------------------------------------------
knitr::kable(table(Signac = celltypes$CellTypes, SignacFast = celltypes_fast$CellTypes), format = "html")

## ----Seurat Visualization 0, message = F, eval = F----------------------------
#  pbmc <- AddMetaData(pbmc, metadata=celltypes_fast$Immune, col.name = "immmune")
#  pbmc <- SetIdent(pbmc, value='immmune')
#  png(filename="fls/plot1_citeseq.png")
#  DimPlot(pbmc)
#  dev.off()

## ----Seurat Visualization 1, message = F, eval = F----------------------------
#  pbmc <- AddMetaData(pbmc, metadata=celltypes$L2, col.name = "celltypes")
#  pbmc <- SetIdent(pbmc, value='celltypes')
#  png(filename="fls/plot2_citeseq.png")
#  DimPlot(pbmc)
#  dev.off()

## ----Seurat Visualization 2, message = F, eval = F----------------------------
#  pbmc <- AddMetaData(pbmc, metadata=celltypes$CellTypes, col.name = "celltypes")
#  pbmc <- SetIdent(pbmc, value='celltypes')
#  png(filename="./fls/plot3_citeseq.png")
#  DimPlot(pbmc)
#  dev.off()

## ----Seurat Visualization 3, message = F, eval = F----------------------------
#  pbmc <- AddMetaData(pbmc, metadata=celltypes$CellTypes_novel, col.name = "celltypes_novel")
#  pbmc <- SetIdent(pbmc, value='celltypes_novel')
#  png(filename="./fls/plot4_citeseq.png")
#  DimPlot(pbmc)
#  dev.off()

## ----Seurat Visualization 4, message = F, eval = F----------------------------
#  pbmc <- AddMetaData(pbmc, metadata=celltypes$CellStates, col.name = "cellstates")
#  pbmc <- SetIdent(pbmc, value='cellstates')
#  png(filename="./fls/plot5_citeseq.png")
#  DimPlot(pbmc)
#  dev.off()

## ----Seurat visualize protein 33, message = F, eval = F-----------------------
#  # Downsample the clusters to a maximum of 500 cells each (makes the heatmap easier to see for
#  # small clusters)
#  pbmc <- SetIdent(pbmc, value='cellstates')
#  
#  # Find markers for all clusters, and draw a heatmap
#  markers <- FindAllMarkers(pbmc, only.pos = TRUE, verbose = F, logfc.threshold = 1)
#  library(dplyr)
#  top5 <- markers %>%  group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#  png(filename="./fls/plot6_citeseq.png")
#  DoHeatmap(pbmc, features = unique(top5$gene), angle = 90)
#  dev.off()

## ----Seurat add protein, eval = F---------------------------------------------
#  pbmc[["ADT"]] <- CreateAssayObject(counts = E$`Antibody Capture`[,colnames(E$`Antibody Capture`) %in% colnames(pbmc)])
#  pbmc <- NormalizeData(pbmc, assay = "ADT", normalization.method = "CLR")
#  pbmc <- ScaleData(pbmc, assay = "ADT")

## ----Seurat visualize protein 3, eval = F-------------------------------------
#  DefaultAssay(pbmc) <- "ADT"
#  # Find protein markers for all clusters, and draw a heatmap
#  adt.markers <- FindAllMarkers(pbmc.small, assay = "ADT", only.pos = TRUE, verbose = F)
#  png(filename="./fls/plot7_citeseq.png")
#  DoHeatmap(pbmc, features = unique(adt.markers$gene), angle = 90)
#  dev.off()

## ----save results, message = F, eval = F--------------------------------------
#  saveRDS(pbmc, file = "fls/pbmcs_signac_citeseq.rds")
#  saveRDS(celltypes, file = "fls/celltypes_citeseq.rds")
#  saveRDS(celltypes_fast, file = "fls/celltypes_fast_citeseq.rds")

## ----save.times, include = FALSE, eval = F------------------------------------
#  write.csv(x = t(as.data.frame(all_times)), file = "fls/tutorial_times_signac-Seurat_citeseq.csv")

## ---- echo=FALSE--------------------------------------------------------------
sessionInfo()

