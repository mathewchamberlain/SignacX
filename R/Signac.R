#' Classification of cellular phenotypes in single cell data
#' 
#' \code{Signac} trains and then uses an ensemble of neural networks to classify cellular phenotypes using an expression matrix or Seurat object.
#' The neural networks are trained with the HPCA training data using only features that are present
#' in both the single cell and HPCA training data set. \code{Signac} returns annotations at each level of the classification
#' hierarchy, which are then converted into cell type labels using \code{\link{GenerateLabels}}. For a faster alternative,
#' try \code{\link{SignacFast}}, which uses pre-computed neural network models.
#'
#' @param E a sparse gene (rows) by cell (column) matrix, or a Seurat object. Rows are HUGO symbols.
#' @param R Reference data. If 'default', R is set to GetTrainingData_HPCA().
#' @param spring.dir If using SPRING, directory to categorical_coloring_data.json. Default is NULL.
#' @param N Number of machine learning models to train (for nn and svm). Default is 100.
#' @param num.cores Number of cores to use. Default is 1.
#' @param threshold Probability threshold for assigning cells to "Unclassified." Default is 0.
#' @param smooth if TRUE, smooths the cell type classifications. Default is TRUE.
#' @param impute if TRUE, gene expression values are imputed prior to cell type classification (see \code{\link{KSoftImpute}}). Default is TRUE.
#' @param verbose if TRUE, code will report outputs. Default is TRUE.
#' @param do.normalize if TRUE, cells are normalized to the mean library size. Default is TRUE.
#' @param return.probability if TRUE, returns the probability associated with each cell type label. Default is TRUE.
#' @param hidden Number of hidden layers in the neural network. Default is 1.
#' @param set.seed If true, seed is set to ensure reproducibility of these results. Default is TRUE.
#' @param seed if set.seed is TRUE, seed is set to 42.
#' @return A list of character vectors: cell type annotations (L1, L2, ...) at each level of the hierarchy
#' as well as 'clusters' for the Louvain clustering results.
#' @seealso \code{\link{SignacFast}}, a faster alternative that only differs from \code{\link{Signac}} in nuanced T cell phenotypes.
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
#' # classify cells
#' labels = Signac(E = pbmc)
#' celltypes = GenerateLabels(labels, E = pbmc)
#' 
#' # add labels to Seurat object, visualize
#' pbmc <- Seurat::AddMetaData(pbmc, metadata=celltypes$CellTypes_novel, col.name = "immmune")
#' pbmc <- Seurat::SetIdent(pbmc, value='immmune')
#' DimPlot(pbmc)
#' 
#' # save results
#' saveRDS(pbmc, "example_pbmcs.rds")
#' }
Signac <- function(E, R = 'default', spring.dir = NULL, N = 100, num.cores = 1, threshold = 0, smooth = TRUE, impute = TRUE, verbose = TRUE, do.normalize = TRUE, return.probability = FALSE, hidden = 1, set.seed = TRUE, seed = '42')
{
  if (!is.null(spring.dir))
    spring.dir = gsub("\\/$", "", spring.dir, perl = TRUE)
  
  if (class(R) == 'character')
    R = GetTrainingData_HPCA()
  
  flag = class(E) == "Seurat"
  
  if (flag){
    default.assay <- Seurat::DefaultAssay(E)
    edges = E@graphs[[which(grepl(paste0(default.assay, "_nn"), names(E@graphs)))]]
  }
  
  if (verbose)
  {
    cat(" ..........  Entry in Signac \n");
    ta = proc.time()[3];
    
    # main function
    if (!flag)
    {
      cat(" ..........  Running Signac on input data matrix :\n");
    } else {
      cat(" ..........  Running Signac on Seurat object :\n");
    }
    cat("             nrow = ", nrow(E), "\n", sep = "");
    cat("             ncol = ", ncol(E), "\n", sep = "");
  }
  
  # keep only unique row names
  logik = CID.IsUnique(rownames(E))
  E = E[logik,]
  
  # intersect genes with reference set
  gns = intersect(rownames(E), R$genes)
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
  res = pbmcapply::pbmclapply(R$Reference, FUN = function(x){
    # keep same gene names
    gns = sort(intersect(rownames(V), colnames(x)))
    Z = V[rownames(V) %in% gns, ]
    dat = x[,colnames(x) %in% gns]
    Z = Z[order(rownames(Z)), ]
    dat = dat[, order(colnames(dat))]
    
    # remove any low variance genes
    kmu = apply(Z, 1, function(x){sum(x != 0)})
    logik = kmu > 0;
    Z = Z[logik,]
    dat = dat[,logik]
    
    # run imputation (if desired)
    if (impute){
      Z = KSoftImpute(E = Z, dM = dM, verbose = FALSE)
      Z = t(apply(Z, 1, function(x){
        normalize(x)
      }))
    }
    
    # build training set
    df = data.frame(dat, celltypes = x$celltypes)
    
    # train a neural network (N times)
    if (set.seed)
    {
      RNGkind("L'Ecuyer-CMRG")
      set.seed(seed = seed)
    }
    res = suppressWarnings(lapply(1:N, function(x) {
      nn=neuralnet::neuralnet(celltypes~.,hidden=hidden,data=df, act.fct = 'logistic', linear.output = FALSE)
      Predict = stats::predict(nn, Matrix::t(Z))
      colnames(Predict) <- sort(nn$model.list$response)
      return(Predict)
    }))
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
    if (smooth & any(as.character(unique(x$celltypes)) %in% c("Immune", "Myeloid", "NonImmune", "Lymphocytes", "Monocytes.Neutrophils", "Monocytes", "Fibroblasts", "Epithelial", "T", "NK", "T.CD8", "T.CD4")))
      df$celltypes = CID.smooth(df$celltypes, dM[[1]])
    
    # return probabilities and cell type classifications
    if (return.probability){
      return(df)
    } else {
      return(df$celltypes)
    }
  }, mc.cores = num.cores)
  
  res$louvain = louvain
  
  if (length(R$Reference) == 1)
    return(as.character(res[[1]]))
  
  if (verbose) {
    tb = proc.time()[3] - ta;
    cat("\n ..........  Exit Signac.\n");
    cat("             Execution time = ", tb, " s.\n", sep = "");
  }
  return(res)
}
