#' Fast classification of cellular phenotypes
#' 
#' \code{SignacFast} uses pre-computed neural network models to classify cellular phenotypes in single cell data:
#' these models were pre-trained with the HPCA training data. Any features that are
#' present in the training data and absent in the single cell data are set to zero. 
#' This is a factor of ~5-10 speed improvement over \code{\link{Signac}}.
#'
#' @param E a gene (rows) by cell (column) matrix, sparse matrix or a Seurat object. Rows are HUGO symbols.
#' @param Models if 'default', as returned by \code{\link{GetModels_HPCA}}. An ensemble of 1,800 neural network models.
#' @param threshold Probability threshold for assigning cells to "Unclassified." Default is 0.
#' @param smooth if TRUE, smooths the cell type classifications. Default is TRUE.
#' @param impute if TRUE, gene expression values are imputed prior to cell type classification (see \code{\link{KSoftImpute}}). Default is TRUE.
#' @param verbose if TRUE, code will report outputs. Default is TRUE.
#' @param do.normalize if TRUE, cells are normalized to the mean library size. Default is TRUE.
#' @param num.cores number of cores to use for parallel computation. Default is 1. 
#' @param return.probability if TRUE, returns the probability associated with each cell type label. Default is TRUE.
#' @param spring.dir If using SPRING, directory to categorical_coloring_data.json. Default is NULL.
#' @seealso \code{\link{Signac}} for another classification function.
#' @return A list of character vectors: cell type annotations (L1, L2, ...) at each level of the hierarchy
#' as well as 'clusters' for the Louvain clustering results.
#' @export
#' @seealso \code{\link{Signac}}
#' @examples
#' \dontrun{
#' # download single cell data for classification
#' file.dir = "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/"
#' file = "pbmc_1k_v3_filtered_feature_bc_matrix.h5"
#' download.file(paste0(file.dir, file), "Ex.h5")
#' 
#' # load data, process with Seurat
#' library(Seurat)
#' E = Read10X_h5(filename = "Ex.h5")
#' pbmc <- CreateSeuratObject(counts = E, project = "pbmc")
#' pbmc <- SCTransform(pbmc)
#' pbmc <- RunPCA(pbmc, verbose = FALSE)
#' pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)
#' pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
#' 
#' # classify cells
#' labels = SignacFast(E = pbmc)
#' celltypes = GenerateLabels(labels, E = pbmc)
#' 
#' # add labels to Seurat object, visualize
#' lbls <- factor(celltypes$CellStates)
#' levels(lbls) <- sort(unique(lbls))
#' pbmc <- AddMetaData(pbmc, metadata=celltypes$CellStates, col.name = "celltypes")
#' pbmc <- SetIdent(pbmc, value='celltypes')
#' DimPlot(pbmc, label = T)
#' 
#' # save results
#' saveRDS(pbmc, "pbmcs.rds")
#' saveRDS(celltypes, "celltypes.rds")
#' }
SignacFast <- function(E, Models = 'default', spring.dir = NULL, num.cores = 1, threshold = 0, smooth = TRUE, impute = TRUE, verbose = TRUE, do.normalize = TRUE, return.probability = FALSE)
{
  if (Models == 'default')
    Models = GetModels_HPCA()
  
  flag = class(E) == "Seurat"
  
  if (flag){
    default.assay <- Seurat::DefaultAssay(E)
    edges = E@graphs[[which(grepl(paste0(default.assay, "_nn"), names(E@graphs)))]]
  }
  
  if (verbose)
  {
    cat(" ..........  Entry in SignacFast \n");
    ta = proc.time()[3];
    
    # main function
    if (!flag)
    {
      cat(" ..........  Running SignacFast on input data matrix :\n");
    } else {
      cat(" ..........  Running SignacFast on Seurat object :\n");
    }
    cat("             nrow = ", nrow(E), "\n", sep = "");
    cat("             ncol = ", ncol(E), "\n", sep = "");
  }
  
  # keep only unique row names
  logik = CID.IsUnique(rownames(E))
  E = E[logik,]
  
  # intersect genes with reference set
  gns = sort(unique(unlist(sapply(Models, function(x){ x$genes }))))
  V = E[rownames(E) %in% gns, ]
  
  if (class(V) %in% "data.frame")
    V = Matrix::Matrix(as.matrix(V), sparse = TRUE)
  
  # normalize to the mean library size
  if (do.normalize)
  {
    if (!flag)
    {
      V = CID.Normalize(V)
    } else {
      V = CID.Normalize(V@assays[[default.assay]]@counts)
    }
  }
  
  # normalization function for gene expression scaling
  normalize <- function(x) {
    return ((x - min(x)) / (max(x) - min(x)))
  }
  
  # normalize V
  V = t(apply(V, 1, function(x){
    normalize(x)
  }))
  logik = apply(V, 1, function(x) {any(is.na(x))})
  V = V[!logik,]
  
  # set up imputation matrices
  if (flag) {
    dM = CID.GetDistMat(edges)
    louvain = CID.Louvain(edges = edges)
  } else {
    edges = CID.LoadEdges(data.dir = spring.dir)
    dM = CID.GetDistMat(edges)
    louvain = CID.Louvain(edges = edges)
  }
  res = pbmcapply::pbmclapply(Models, FUN = function(x){
    
    Z = V[rownames(V) %in% x$genes,]
    
    # run imputation (if desired)
    if (impute){
      Z = KSoftImpute(E = Z, dM = dM, verbose = FALSE)
      Z = t(apply(Z, 1, function(x){
        normalize(x)
      }))
    }
    
    # zero impute in any missing genes
    dummy = Matrix::Matrix(0, nrow = length(x$genes) - nrow(Z), ncol = ncol(Z))
    rownames(dummy) <- x$genes[!x$genes %in% rownames(Z)]
    Z = rbind(Z, dummy)
    Z = Z[order(rownames(Z)),]
    
    # generate predictions
    res = lapply(x$classifiers, function(y) {
      Predict = stats::predict(y, Matrix::t(Z))
      colnames(Predict) <- sort(y$model.list$response)
      return(Predict)
    })
    res = res[sapply(res, function(x) !is.null(x))]
    res.squared.mean <- Reduce("+", lapply(res, "^", 2)) / length(res)
    res = Reduce(res, f = '+') / length(res)
    res.variance <- res.squared.mean - res^2
    res.sd <- sqrt(res.variance)
    xx = apply(res, 1, which.max)
    celltypes = colnames(res)[xx]
    kmax = apply(res, 1, max)
    celltypes[kmax < threshold] = "Unclassified"
    errors = round(sapply(1:length(xx), function(x){res.sd[x, xx[x]]}), digits = 4)
    df = data.frame(celltypes = celltypes, probability = round(kmax, digits = 3), sd = errors, percent_features_detected = round(Matrix::colSums(Z != 0) / nrow(Z), digits = 3) * 100)
    
    # smooth the output classifications
    if (smooth & any(as.character(unique(x$celltypes)) %in% c("Immune", "Myeloid", "NonImmune", "Lymphocytes", "Monocytes.Neutrophils", "Monocytes", "Fibroblasts", "Epithelial")))
      df$celltypes = CID.smooth(df$celltypes, dM[[1]])
    
    # return probabilities and cell type classifications
    if (return.probability){
      return(df)
    } else {
      return(df$celltypes)
    }
  }, mc.cores = num.cores)
  
  res$louvain = louvain
  
  if (verbose) {
    tb = proc.time()[3] - ta;
    cat("\n ..........  Exit SignacFast.\n");
    cat("             Execution time = ", tb, " s.\n", sep = "");
  }
  return(res)
}
