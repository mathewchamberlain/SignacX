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
celltypes_fast = readRDS("./fls/celltypes_fast.rds")
celltypes = readRDS("./fls/celltypes.rds")
# pbmc = readRDS("fls/pbmcs_signac.rds")

## ----setupSeurat, message = F, eval = F---------------------------------------
#  library(Seurat)

## ----setup, message = F, eval = F---------------------------------------------
#  dir.create("fls")
#  download.file("https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.h5", destfile = "fls/pbmc_1k_v3_filtered_feature_bc_matrix.h5")

## ----Seurat, message = T, eval = F--------------------------------------------
#  # load data
#  E = Read10X_h5(filename = "fls/pbmc_1k_v3_filtered_feature_bc_matrix.h5")
#  pbmc <- CreateSeuratObject(counts = E, project = "pbmc")
#  
#  # run sctransform
#  pbmc <- SCTransform(pbmc, verbose = FALSE)

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
#  pbmc <- AddMetaData(pbmc, metadata=celltypes$Immune, col.name = "immmune")
#  pbmc <- SetIdent(pbmc, value='immmune')
#  png(filename="fls/plot1.png")
#  DimPlot(pbmc)
#  dev.off()

## ----Seurat Visualization 1, message = F, eval = F----------------------------
#  pbmc <- AddMetaData(pbmc, metadata=celltypes$L2, col.name = "L2")
#  pbmc <- SetIdent(pbmc, value='L2')
#  png(filename="fls/plot2.png")
#  DimPlot(pbmc)
#  dev.off()

## ----Seurat Visualization 2, message = F, eval = F----------------------------
#  lbls = factor(celltypes$CellTypes)
#  levels(lbls) <- sort(unique(lbls))
#  pbmc <- AddMetaData(pbmc, metadata=lbls, col.name = "celltypes")
#  pbmc <- SetIdent(pbmc, value='celltypes')
#  png(filename="./fls/plot3.png")
#  DimPlot(pbmc)
#  dev.off()

## ----Seurat Visualization 3, message = F, eval = F----------------------------
#  pbmc <- AddMetaData(pbmc, metadata=celltypes$CellTypes_novel, col.name = "celltypes_novel")
#  pbmc <- SetIdent(pbmc, value='celltypes_novel')
#  png(filename="./fls/plot4.png")
#  DimPlot(pbmc)
#  dev.off()

## ----Seurat Visualization 4, message = F, eval = F----------------------------
#  pbmc <- AddMetaData(pbmc, metadata=celltypes$CellStates, col.name = "cellstates")
#  pbmc <- SetIdent(pbmc, value='cellstates')
#  png(filename="./fls/plot5.png")
#  DimPlot(pbmc)
#  dev.off()

## ----Seurat visualize protein 3, message = F, eval = F------------------------
#  pbmc <- SetIdent(pbmc, value='celltypes_novel')
#  
#  # Find protein markers for all clusters, and draw a heatmap
#  markers <- FindAllMarkers(pbmc, only.pos = TRUE, verbose = F, logfc.threshold = 1)
#  library(dplyr)
#  top5 <- markers %>%  group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#  png(filename="./fls/plot6.png")
#  DoHeatmap(pbmc, features = unique(top5$gene), angle = 90)
#  dev.off()

## ----save results, message = F, eval = F--------------------------------------
#  saveRDS(pbmc, file = "fls/pbmcs_signac.rds")
#  saveRDS(celltypes, file = "fls/celltypes.rds")
#  saveRDS(celltypes_fast, file = "fls/celltypes_fast.rds")

## ----save.times, include = FALSE, eval = F------------------------------------
#  write.csv(x = t(as.data.frame(all_times)), file = "fls/tutorial_times_signac-Seurat_pbmcs.csv")

## ---- echo=FALSE--------------------------------------------------------------
sessionInfo()

