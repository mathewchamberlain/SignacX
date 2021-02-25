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
# kidney = readRDS(kidney, file = "fls/amp_kidney_signac.rds")
Q = readRDS(file = "fls/amp_kidney_MASC_result.rds")
celltypes = readRDS(file = "fls/amp_kidney_celltypes.rds")

## ----read CELSeq2, message = F, eval = F--------------------------------------
#  ReadCelseq <- function (counts.file, meta.file)
#  {
#    E = suppressWarnings(readr::read_tsv(counts.file));
#    gns <- E$gene;
#    E = E[,-1]
#    M = suppressWarnings(readr::read_tsv(meta.file))
#    S = lapply(unique(M$sample), function(x) {
#      logik = colnames(E) %in% M$cell_name[M$sample == x]
#      Matrix::Matrix(as.matrix(E[,logik]), sparse = TRUE)}
#      )
#    names(S) <- unique(M$sample)
#    S = lapply(S, function(x){
#      rownames(x) <- gns
#      x
#    })
#    S
#  }
#  
#  counts.file = "./SDY997_EXP15176_celseq_matrix_ru10_molecules.tsv.gz"
#  meta.file = "./SDY997_EXP15176_celseq_meta.tsv"
#  
#  E = ReadCelseq(counts.file = counts.file, meta.file = meta.file)
#  M = suppressWarnings(readr::read_tsv(meta.file))

## ----filter celseq, message = F, eval = F-------------------------------------
#  # keep only kidney cells
#  E = lapply(E, function(x){
#    x[,grepl("K", colnames(x))]
#  })
#  
#  # remove any sample with no cells
#  E = E[sapply(E, ncol) > 0]
#  
#  # merge into a single matrix
#  E = do.call(cbind, E)

## ----setupSeurat, message = F, eval = F---------------------------------------
#  library(Seurat)

## ----Seurat, message = T, eval = F--------------------------------------------
#  # load data
#  kidney <- CreateSeuratObject(counts = E, project = "celseq")
#  
#  # run sctransform
#  kidney <- SCTransform(kidney)

## ----Seurat 2, message = T, eval = F------------------------------------------
#  # These are now standard steps in the Seurat workflow for visualization and clustering
#  kidney <- RunPCA(kidney, verbose = FALSE)
#  kidney <- RunUMAP(kidney, dims = 1:30, verbose = FALSE)
#  kidney <- FindNeighbors(kidney, dims = 1:30, verbose = FALSE)

## ----setup signacX, message = F, eval = F-------------------------------------
#  install.packages("SignacX")

## ----Signac, message = T, eval = F--------------------------------------------
#  # Run Signac
#  library(SignacX)
#  labels <- Signac(kidney, num.cores = 4)
#  celltypes = GenerateLabels(labels, E = kidney)

## ----Seurat Visualization 0, eval = F-----------------------------------------
#  kidney <- AddMetaData(kidney, metadata=celltypes$Immune, col.name = "immmune")
#  kidney <- SetIdent(kidney, value='immmune')
#  png(filename="fls/plot1_amp.png")
#  DimPlot(kidney)
#  dev.off()

## ----Seurat Visualization 2, eval = F-----------------------------------------
#  kidney <- AddMetaData(kidney, metadata=celltypes$CellTypes, col.name = "celltypes")
#  kidney <- SetIdent(kidney, value='celltypes')
#  png(filename="fls/plot2_amp.png")
#  DimPlot(kidney)
#  dev.off()

## ----Seurat Visualization 3, eval = F-----------------------------------------
#  kidney <- AddMetaData(kidney, metadata=celltypes$CellStates, col.name = "cellstates")
#  kidney <- SetIdent(kidney, value='cellstates')
#  png(filename="fls/plot3_amp.png")
#  DimPlot(kidney)
#  dev.off()

## ----Seurat get IMAGES, message = F, eval = F---------------------------------
#  # Downsample just the immune cells
#  kidney.small <- kidney[, !celltypes$CellStates %in% c("NonImmune", "Fibroblasts", "Unclassified", "Endothelial", "Epithelial")]
#  
#  # Find protein markers for all clusters, and draw a heatmap
#  markers <- FindAllMarkers(kidney.small, only.pos = TRUE, verbose = F, logfc.threshold = 1)
#  require(dplyr)
#  top3 <- markers %>%  group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
#  png(filename="fls/plot4_amp.png")
#  DoHeatmap(kidney.small, features = unique(top3$gene), angle = 90)
#  dev.off()

## ----MASC, message = F, eval = F----------------------------------------------
#  Meta_mapped = M[match(colnames(kidney), M$cell_name), ]
#  Meta_mapped$CellStates = celltypes$CellStates
#  Meta_mapped$disease = factor(Meta_mapped$disease)
#  Q = MASC(dataset = Meta_mapped, cluster = Meta_mapped$CellStates, contrast = 'disease', random_effects = c("plate", "lane", "sample"))

## ---- results='asis'----------------------------------------------------------
writeLines("td, th { padding : 6px } th { background-color : brown ; color : white; border : 1px solid white; } td { color : brown ; border : 1px solid brown }", con = "mystyle.css")
knitr::kable(Q, format = "html")

## ----save results, message = F, eval = F--------------------------------------
#  saveRDS(kidney, file = "fls/amp_kidney_signac.rds")
#  saveRDS(Q, file = "fls/amp_kidney_MASC_result.rds")
#  saveRDS(celltypes, file = "fls/amp_kidney_celltypes.rds")

## ----save.times, include = FALSE, eval = F------------------------------------
#  write.csv(x = t(as.data.frame(all_times)), file = "fls/tutorial_times_signac-Seurat_AMP.csv")

## ---- echo=FALSE--------------------------------------------------------------
sessionInfo()

