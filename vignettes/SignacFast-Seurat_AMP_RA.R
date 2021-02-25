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
#celltypes_fast = readRDS("./fls/celltypes_fast_citeseq.rds")
#celltypes = readRDS("./fls/celltypes_citeseq.rds")
# pbmc = readRDS("fls/pbmcs_signac_citeseq.rds")
celltypes = readRDS(file = "fls/celltypes_amp_synovium.rds")
celltypes_fast = readRDS(file = "fls/celltypes_fast_synovium_celltypes.rds")

## ----read CELSeq2, message = F, eval = F--------------------------------------
#  ReadCelseq <- function (counts.file, meta.file)
#  {
#    E = suppressWarnings(readr::read_tsv(counts.file));
#    gns <- E$gene;
#    E = E[,-1]
#    E = Matrix::Matrix(as.matrix(E), sparse = TRUE)
#    rownames(E) <- gns
#    E
#  }
#  
#  counts.file = "./fls/celseq_matrix_ru10_molecules.tsv.gz"
#  meta.file = "./fls/celseq_meta.immport.723957.tsv"
#  
#  E = ReadCelseq(counts.file = counts.file, meta.file = meta.file)
#  M = suppressWarnings(readr::read_tsv(meta.file))
#  
#  # filter data based on depth and number of genes detected
#  kmu = Matrix::colSums(E != 0)
#  kmu2 = Matrix::colSums(E)
#  E = E[,kmu > 200 & kmu2 > 500]
#  
#  # filter by mitochondrial percentage
#  logik = grepl("^MT-", rownames(E))
#  MitoFrac = Matrix::colSums(E[logik,]) / Matrix::colSums(E) * 100
#  E = E[,MitoFrac < 20]

## ----setupSeurat, message = F, eval = F---------------------------------------
#  library(Seurat)

## ----Seurat, message = T, eval = F--------------------------------------------
#  # load data
#  synovium <- CreateSeuratObject(counts = E, project = "FACs")
#  
#  # run sctransform
#  synovium <- SCTransform(synovium)

## ----Seurat 2, message = T, eval = F------------------------------------------
#  # These are now standard steps in the Seurat workflow for visualization and clustering
#  synovium <- RunPCA(synovium, verbose = FALSE)
#  synovium <- RunUMAP(synovium, dims = 1:30, verbose = FALSE)
#  synovium <- FindNeighbors(synovium, dims = 1:30, verbose = FALSE)

## ----setup2, message = F, eval = F--------------------------------------------
#  library(SignacX)

## ----Signac, message = T, eval = F--------------------------------------------
#  labels <- Signac(synovium, num.cores = 4)
#  celltypes = GenerateLabels(labels, E = synovium)

## ----SignacFast, message = T, eval = F----------------------------------------
#  # Run SignacFast
#  labels_fast <- SignacFast(synovium, num.cores = 4)
#  celltypes_fast = GenerateLabels(labels_fast, E = synovium)

## ----Compare 1, echo = F, eval = T--------------------------------------------
knitr::kable(table(Signac = celltypes$CellTypes, SignacFast = celltypes_fast$CellTypes), format = "html")

## ----Compare 2, echo = F, eval = T--------------------------------------------
knitr::kable(table(Signac = celltypes$CellStates, SignacFast = celltypes_fast$CellStates), format = "html")

## ----save results, message = F, eval = F--------------------------------------
#  saveRDS(synovium, file = "fls/seurat_obj_amp_synovium.rds")
#  saveRDS(celltypes, file = "fls/celltypes_amp_synovium.rds")
#  saveRDS(celltypes_fast, file = "fls/celltypes_fast_amp_synovium_celltypes.rds")

## ----save.times, include = FALSE, eval = F------------------------------------
#  write.csv(x = t(as.data.frame(all_times)), file = "fls/tutorial_times_signac-Seurat_AMP.csv")

## ---- echo=FALSE--------------------------------------------------------------
sessionInfo()

