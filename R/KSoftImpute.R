#' KNN-based imputation
#' 
#' \code{KSoftImpute} is an ultra-fast method for imputing missing gene expression values in single cell data.
#' \code{KSoftImpute} uses k-nearest neighbors to impute the expression of each gene by the weighted average of itself
#' and it's first-degree neighbors. Weights for imputation are determined by the number of detected genes. This method
#' works for large data sets (>100,000 cells) in under a minute.
#'
#' @param E A gene-by-sample count matrix (sparse matrix or matrix) with genes identified by their HUGO symbols.
#' @param dM see ?CID.GetDistMat
#' @param genes.to.use a character vector of genes to impute. Default is NULL.
#' @param verbose If TRUE, code reports outputs. Default is FALSE.
#' @return An expression matrix (sparse matrix) with imputed values.
#' @seealso \code{\link{Signac}} and \code{\link{SignacFast}}
#' @export
#' See \url{https://doi.org/10.1101/2021.02.01.429207} for more details.
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
#' 
#' # run Seurat pipeline
#' pbmc <- SCTransform(pbmc, verbose = FALSE)
#' pbmc <- RunPCA(pbmc, verbose = FALSE)
#' pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)
#' pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
#' 
#' # get edges from default assay from Seurat object
#' default.assay <- Seurat::DefaultAssay(pbmc)
#' edges = pbmc@graphs[[which(grepl(paste0(default.assay, "_nn"), names(pbmc@graphs)))]]
#' 
#' # get distance matrix
#' dM = CID.GetDistMat(edges)
#' 
#' # run imputation
#' Z = KSoftImpute(E = E, dM = dM, verbose = TRUE)
#' }
KSoftImpute <- function(E,  dM = NULL, genes.to.use = NULL, verbose = FALSE)
{
  # check inputs
  if (verbose)
  {
    cat(" ..........  Entry in KSoftImpute \n");
    ta = proc.time()[3];
    
    # main function
    cat(" ..........  Running K soft imputation on input data matrix :\n");
    cat("             nrow = ", nrow(E), "\n", sep = "");
    cat("             ncol = ", ncol(E), "\n", sep = "");
  }
  
  if (!is.null(genes.to.use))
    E = E[rownames(E) %in% genes.to.use,]
  
  # define the weight matrix Wjj, which is the number of genes detected in each cell
  weights = Matrix::colSums(E != 0)
  bas = Matrix::Matrix(0, nrow = length(weights), ncol = length(weights))
  diag(bas) <- weights
  
  g = dM[[1]] %*% bas 
  dd = g / (Matrix::rowSums(g) + 1)
  diag(dd) <- 1
  E_new = E %*% dd;
  
  if (verbose) {
    tb = proc.time()[3] - ta;
    cat("\n ..........  Exit KSoftImpute.\n");
    cat("             Execution time = ", tb, " s.\n", sep = "");
  }
  
  return(E_new)
}