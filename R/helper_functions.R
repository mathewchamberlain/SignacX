#' Loads neural network models from GitHub
#'
#' \code{GetModels_HPCA} returns a list of neural network models that were trained with the HPCA training data.
#' The HPCA training data has 18 classification tasks with bootstrapped training data.
#' The training data are split into two distinct cellular phenotypes (i.e., immune and nonimmune).
#' \code{GetModels_HPCA} downloads pre-computed neural networks trained to classify cells for each task (n = 100 models trained for each tasks).
#'
#' @return list of 1,800 neural network models stacked to solve 18 classification problems.
#' @export
#' @seealso \code{\link{SignacFast}}, \code{\link{ModelGenerator}}
#' @examples
#' \dontrun{
#' P = GetModels()
#' }
GetModels_HPCA <- function(){
  return(readRDS(url("https://github.com/mathewchamberlain/SignacX/blob/master/assets/Models_HPCA.rds?raw=TRUE","rb")))
}

#' Loads bootstrapped HPCA training data from GitHub
#'
#' \code{GetTrainingData_HPCA} returns a list of bootstrapped HPCA training data.
#' The HPCA training data has 18 classification tasks, each split into two distinct cellular phenotypes (i.e., immune and nonimmune).
#' \code{GetTrainingData_HPCA} downloads bootstrapped training data to use for model-building.
#' 
#' @return A list with two elements; 'Reference' are the training data, and 'genes' -- the union of all features present in the training data.
#' @export
#' @examples
#' \dontrun{
#' P = GetTrainingData_HPCA()
#' }
GetTrainingData_HPCA <- function(){
  return(readRDS(url("https://github.com/mathewchamberlain/SignacX/blob/master/assets/training_HPCA.rds?raw=TRUE","rb")))
}

#' Generates cellular phenotype labels
#'
#' \code{GenerateLabels} returns a list of cell type and cell state labels, as well as novel cellular phenotypes and unclassified cells.
#'
#' @param cr list returned by \code{\link{Signac}} or by \code{\link{SignacFast}}.
#' @param E a sparse gene (rows) by cell (column) matrix, or a Seurat object. Rows are HUGO symbols.
#' @param spring.dir If using SPRING, directory to categorical_coloring_data.json. Default is NULL.
#' @param smooth if TRUE, smooths the cell type classifications. Default is TRUE.
#' @param new_populations Character vector specifying any new cell types that were learned by Signac. Default is NULL.
#' @param new_categories If new_populations are set to a cell type, new_category is a corresponding character vector indicating the population that the new population belongs to. Default is NULL.
#' @param min.cells If desired, any cell population with equal to or less than N cells is set to "Unclassified." Default is 10 cells.
#' @return A list of cell type labels for cell types, cell states and novel populations.
#' @export
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
#' # download bootstrapped reference data for training models
#' file.dir = "https://github.com/mathewchamberlain/Signac/blob/master/data/"
#' file = "training_HPCA.rda"
#' download.file(paste0(file.dir, file, "?raw=true"), destfile = "training_HPCA.rda")
#' load("training_HPCA.rda")
#' 
#' # classify cells
#' labels = SignacFast(E = pbmc, R = training_HPCA)
#' celltypes = GenerateLabels(labels, E = pbmc)
#' }
GenerateLabels = function(cr, E = NULL, smooth = TRUE, new_populations = NULL, new_categories = NULL, min.cells = 10,  spring.dir = NULL)
{
  
  # if using SPRING, load data
  if (!is.null(spring.dir)){
    edges = CID.LoadEdges(data.dir = spring.dir)
    dM = CID.GetDistMat(edges)
  }
  
  # check for Seurat object
  flag = class(E) == "Seurat"
  if (flag) {
    default.assay <- Seurat::DefaultAssay(E)
    edges = E@graphs[[which(grepl(paste0(default.assay, "_nn"), names(E@graphs)))]]
    dM = CID.GetDistMat(edges)
  }
  
  # remove louvain clusters
  res = list("")
  louvain = cr$louvain
  cr = cr[-which(names(cr) == "louvain")]
  
  # remove all fields from input
  if ("probability" %in% names(cr[[1]]))
    cr = lapply(cr, function(x) {x$celltypes})
  
  # change format to character
  cr = lapply(cr, function(x) {as.character(x)})
  
  # label categories according to cell type hierarchy
  for (j in 1:length(cr))
  {
    if (names(cr)[j] == "All"){
      res[[j]] = cr[[j]]
    } else {
      qq = res[[j - 1]]
      logik = qq == names(cr)[j]
      qq[logik] = cr[[j]][logik]
      res[[j]] = qq
    }
  }
  
  # define cell states and cell types 
  cellstates = res[[length(res)]]
  celltypes = cellstates
  celltypes[celltypes %in% c("B.memory", "B.naive")] = "B"
  celltypes[celltypes %in% c("DC", "Mon.Classical", "Mon.NonClassical", "Neutrophils", "Monocytes", "Macrophages")] = "MPh"
  celltypes[celltypes %in% c("NK", "T.CD4.naive", "T.CD4.memory", "T.regs", "T.CD8.naive", "T.CD8.memory", "T.CD8.cm","T.CD8.em")] = "TNK"
  celltypes[celltypes %in% c("Endothelial", "Fibroblasts", "Epithelial")] = "NonImmune"
  immune = res[[1]]  

  # allow for new cell populations
  if (!is.null(new_populations))
  {
    for (j in 1:length(new_populations))
      celltypes[celltypes %in% new_populations[j]] = new_categories[j]
  }
  
  # assign Unclassifieds
  if (!is.null(spring.dir) | flag){
  celltypes = CID.entropy(celltypes, dM)
  immune = CID.entropy(immune, dM)
  # smooth 
  if (smooth) {
    celltypes= CID.smooth(celltypes, dM[[1]])
    immune = CID.smooth(immune, dM[[1]])
  }
  }
  # set consistent unclassified cells
  logik = immune == "Unclassified" | celltypes == "Unclassified" | cellstates == "Unclassified"
  cellstates[logik] = "Unclassified"
  celltypes[logik] = "Unclassified"
  immune[logik] = "Unclassified"
  
  res$Immune = immune
  if (!is.null(spring.dir) | flag)
  {
  do = data.frame(table(louvain[cellstates == "Unclassified"]))
  df = data.frame(table(louvain[louvain %in% do$Var1]))
  logik = (1 - phyper(do$Freq, df$Freq , length(cellstates) - do$Freq, sum(cellstates == "Unclassified"))) < 0.01;
  if (any(logik)){
    do = do[logik,]
    logik = do$Freq > 3; # require at least 3 cell communities
    if (any(logik)){
      celltypes_novel  = celltypes;
      cellstates_novel = cellstates;
      cat("             Signac found", sum(logik), "novel celltypes!\n");
      lbls = rep("All", ncol(E))
      logik = louvain %in% do$Var1[logik] & cellstates == "Unclassified";
      lbls[logik] = louvain[logik]
      cellstates_novel[lbls != "All"] = paste0("+ Novel cluster", lbls[lbls != "All"])
      celltypes_novel[lbls != "All"] = paste0("+ Novel cluster", lbls[lbls != "All"])
      res$CellTypes_novel = celltypes_novel
      res$CellStates_novel = cellstates_novel
      }
    } 
  }

  df = data.frame(table(celltypes))
  logik = df$Freq < (min.cells + 1)
  if (any(logik))
    celltypes[celltypes %in% as.character(df[,1][logik])] <- "Unclassified"
  
  df = data.frame(table(cellstates))
  logik = df$Freq < (min.cells + 1)
  if (any(logik))
    cellstates[cellstates %in% as.character(df[,1][logik])] <- "Unclassified"
  
  res$CellTypes = celltypes
  res$CellStates = cellstates

  names(res)[names(res) == ""] = paste0("L", 1:sum(names(res) == ""))
  
  res$clusters = louvain
  
  # factor encode for Seurat visualization
  if (flag){
    res = lapply(res, function(x){
      xx = factor(x)
      levels(xx) <- sort(unique(xx))
      return(xx)
    })
  }

  return(res)
}

#' Generates bootstrapped single cell data
#'
#' \code{SignacBoot} uses a Seurat object or an expression matrix and performs feature selection, normalization and bootstrapping
#' to generate a training data set to be used for cell type or cluster classification.
#'
#' @param E a gene (rows) by cell (column) matrix, sparse matrix or a Seurat object. Rows are HUGO symbols.
#' @param L cell type categories for learning.
#' @param labels cell type labels corresponding to the columns of E.
#' @param size Number of bootstrapped samples for machine learning. Default is 1,000.
#' @param logfc.threshold Cutoff for feature selection. Default is 0.25.
#' @param p.val.adj Cutoff for feature selection. Default is 0.05.
#' @param impute if TRUE, performs imputation prior to bootstrapping (see \code{\link{KSoftImpute}}). Default is TRUE.
#' @param verbose if TRUE, code speaks. Default is TRUE.
#' @param spring.dir if using SPRING, directory to categorical_coloring_data.json. Default is NULL.
#' @return Training data set (data.frame) to be used for building new models=.
#' @seealso \code{\link{ModelGenerator}}
#' @export
#' @examples
#' \dontrun{
#' # load Seurat object from SignacFast example
#' P <- readRDS("pbmcs.rds")
#' 
#' # run feature selection + bootstrapping to generate 2,000 bootstrapped cells
#' x = P@meta.data$celltypes
#' R_learned = SignacBoot(P, L = c("B.naive", "B.memory"), labels = x)
#' }
SignacBoot <- function (E, L, labels, size = 1000, impute = TRUE, spring.dir = NULL, logfc.threshold = 0.25, p.val.adj = 0.05, verbose = TRUE)
{
  # check for Seurat object
  flag = class(E) == "Seurat"
  if (flag) {
    default.assay <- Seurat::DefaultAssay(E)
    edges = E@graphs[[which(grepl(paste0(default.assay, "_nn"), names(E@graphs)))]]
    dM = CID.GetDistMat(edges)
  }
  
  if (verbose)
  {
    cat(" ..........  Entry in SignacBoot \n");
    ta = proc.time()[3];
    
    # main function
    cat(ifelse(!flag, c(" ..........  Running SignacBoot on input data matrix :\n"), c(" ..........  Running SignacBoot on input Seurat object :\n")))
    cat("             nrow = ", nrow(E), "\n", sep = "");
    cat("             ncol = ", ncol(E), "\n", sep = "");
  }
  
  # set up imputation matrices
  if (impute & !flag){
    edges = CID.LoadEdges(data.dir = spring.dir)
    dM = CID.GetDistMat(edges)
    louvain = CID.Louvain(edges = edges)
  } else {
    louvain = CID.Louvain(edges = edges)
  }
  
  # keep only unique row names
  logik = CID.IsUnique(rownames(E))
  if (!flag) {
    V = E[logik,]
    colnames(V) <- 1:ncol(V)
  } else {
    V = E@assays[[default.assay]]@counts[logik,]
  }

  # run feature selection
  ctrl <- Seurat::CreateSeuratObject(counts = V[,labels %in% L], project = "Signac", min.cells = 0)
  ctrl <- Seurat::NormalizeData(ctrl)
  ctrl <- Seurat::AddMetaData(ctrl, metadata=labels[labels %in% L], col.name = "celltypes")
  ctrl <- Seurat::SetIdent(ctrl, value='celltypes')
  mrks = Seurat::FindMarkers(ctrl, ident.1 = L[1], ident.2 = L[2], only.pos = FALSE, logfc.threshold = logfc.threshold)
  mrks = mrks[mrks$p_val_adj < p.val.adj,]
  # bootstrap data
  dat = V[rownames(V) %in% rownames(mrks),]
  mrks$cluster = L[1]
  mrks$cluster[mrks$avg_logFC < 0] = L[2]
  mrks$gene = rownames(mrks)
  
  xx = labels
  
  # run imputation (if desired)
  if (impute)
    Z = KSoftImpute(E = dat, dM = dM, verbose = FALSE)
  
  cts = split.data.frame(mrks, f = mrks$cluster)
  N = lapply(cts, function(x){
      # first sample from cells in cluster 1, size cells
      logik = rownames(dat) %in% x$gene
      dummy = dat[logik,]
      logik = as.character(x$cluster[1]) == labels
      dummy = dummy[,logik]
      dd = t(apply(dummy, 1, function(z) {
        sample(z, size = size, replace = TRUE)}))
      dd = t(apply(dd, 1, function(z){
        stats::rnorm(n = length(z), mean = mean(z), sd = sd(z))
      }))
      # now sample from cells in cluster 2, size cells with True Negative Expr.
      logik = !rownames(dat) %in% x$gene
      dummy = dat[logik,]
      logik = as.character(x$cluster[1]) == labels
      dummy = dummy[,logik]
      dd2 = t(apply(dummy, 1, function(z) {
        sample(z, size = size, replace = TRUE)}))
      dd2 = t(apply(dd2, 1, function(z){
        stats::rnorm(n = length(z), mean = mean(z), sd = sd(z))
      }))
      rbind(dd, dd2)
    })
    N2 = merge(N[[1]],N[[2]],by="row.names")
    rownames(N2) <- N2$Row.names
    N2 = t(N2[,-1])
    # normalize
    normalize <- function(x) {
      return ((x - min(x)) / (max(x) - min(x)))
    }
    N2 = apply(N2, 2, function(x){
      normalize(x)
    })
  boot = data.frame(N2, celltypes = c(rep(names(cts)[1], size), rep(names(cts)[2], size)))
  colnames(boot) <- gsub(pattern = "\\.", replacement = "-", x = colnames(boot))
  
  boot = list(genes = colnames(boot)[-ncol(boot)],
              Reference = list(boot))
  
  if (verbose) {
    tb = proc.time()[3] - ta;
    cat("\n ..........  Exit SignacBoot.\n");
    cat("             Execution time = ", tb, " s.\n", sep = "");
  }
  return(boot)
}

#' Load edges from edge list for single cell network
#' 
#' \code{CID.LoadEdges} loads edges, typically after running the SPRING pipeline.
#'
#' @param data.dir A directory where "edges.csv" file is located
#' @return The edge list in data frame format
#' @export
#' @examples
#' \dontrun{
#' # Loads edges
#' file.dir = "https://kleintools.hms.harvard.edu/tools/client_datasets/"
#' file = "CITESEQ_EXPLORATORY_CITESEQ_5K_PBMCS/FullDataset_v1_protein/edges.csv"
#' download.file(paste0(file.dir, file, "?raw=true"), destfile = "edges.csv")
#' 
#' # data.dir is your path to the "edges.csv" file
#' edges = CID.LoadEdges(data.dir = ".")
#' }
CID.LoadEdges <- function(data.dir)
{
  edges = paste(data.dir, "edges.csv", sep = "/")
  file.exists(edges)
  flag = file.exists(edges);
  if (!flag) {
    cat("ERROR: from CID.LoadEdges:\n");
    cat("edges = ", edges, " does not exist.\n", sep = "");
    stop()
  }
  # read in knn graph edges
  cat(" ..........  Reading in graph edges from edgelist: \n")
  cat("             ", edges, "\n")
  edges <- read.delim(edges,sep = ";" , stringsAsFactors = FALSE, header = FALSE)
  if(min(edges) == 0)
    edges = edges + 1
  colnames(edges) <- c("V1", "V2")
  edges
}

#' Library size normalize
#' 
#' \code{CID.Normalize} normalizes the expression matrix to the mean library size.
#'
#' @param E Expression matrix
#' @return Normalized expression matrix where each cell sums to the mean total counts
#' @export
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
#' 
#' # observe average total counts
#' mean(Matrix::colSums(E))
#' 
#' # run normalization
#' E_norm = CID.Normalize(E)
#' 
#' # check normalization
#' kmu = mean(Matrix::colSums(E_norm))
#' head(kmu)
#' }
CID.Normalize <- function(E)
{
  xx = NULL
  if (!is.null(colnames(E)))
    xx <- colnames(E)
  if (class(E) == "matrix")
  {
    m = Matrix::Matrix(0, ncol(E), ncol(E))
    tots_use = Matrix::colSums(E)
    target_mean = mean(tots_use)
    diag(m) <- target_mean / tots_use
    E = as.matrix(E %*% m)
    if (!is.null(xx))
      colnames(E) <- xx
  } else {
    m = Matrix::Matrix(0, ncol(E), ncol(E))
    tots_use = Matrix::colSums(E)
    target_mean = mean(tots_use)
    diag(m) <- target_mean / tots_use
    E = Matrix::Matrix(E %*% m, sparse = TRUE)
    if (!is.null(xx))
      colnames(E) <- xx
  }
  return(E)
}

#' Smoothing function
#' 
#' \code{CID.smooth} uses k-nearest neighbors to identify cells which
#' correspond to a different label than the majority of their first-degree neighbors. If so,
#' those annotations are "smoothed."
#'
#' @param ac list containing a character vector where each element is a cell type or cell state assignment.
#' @param dM distance matrix (see ?CID.GetDistMat).
#' @return A character vector with smoothed labels
#' @export
#' @examples
#' \dontrun{
#' # load data classified previously (see SignacFast)
#' P <- readRDS("celltypes.rds")
#' S <- readRDS("pbmcs.rds")
#' 
#' # get edges from default assay from Seurat object
#' default.assay <- Seurat::DefaultAssay(S)
#' edges = S@graphs[[which(grepl(paste0(default.assay, "_nn"), names(S@graphs)))]]
#' 
#' # get distance matrix
#' D = CID.GetDistMat(edges)
#' 
#' # smooth labels
#' smoothed = CID.smooth(ac = P$CellTypes, dM = D[[1]])
#' }
CID.smooth <- function(ac, dM)
{
  Y = ac;
  # remove any unconnected cells from consideration
  logik = Matrix::rowSums(dM) != 0
  if(any(!logik)){
    dM = dM[logik,logik]
    Y = Y[logik]
  }
  if (sum(dM) == 0)
    return(ac)
  m = Matrix::Matrix(0, nrow = length(Y), ncol = length(unique(Y)), sparse = TRUE)
  m[cbind(1:nrow(m), as.numeric(factor(Y)))] <- 1
  res = dM %*% m
  res = res / Matrix::rowSums(res)
  mx.idx = apply(res, 1, function(x) which.max(x))
  mx = apply(res, 1, max)
  Y[mx > 0.5] = levels(factor(Y))[mx.idx[mx > 0.5]]
  ac[logik] = Y
  return(ac)
}

#' Normalized Shannon entropy-based "unclassified" assignment
#' 
#' \code{CID.entropy} calculates the normalized Shannon entropy of labels for each cell
#' among k-nearest neighbors less than four-degrees apart, and then sets cells with statistically
#' significant large Shannon entropy to be "Unclassified."
#' 
#' @param ac a character vector of cell type labels
#' @param distM the distance matrix, see ?CID.GetDistMat
#' @return A character vector like 'ac' but with cells type labels set to "Unclassified" if there was high normalized Shannon entropy.
#' @export
#' @examples
#' \dontrun{
#' # load data classified previously (see \code{SignacFast})
#' P <- readRDS("celltypes.rds")
#' S <- readRDS("pbmcs.rds")
#' 
#' # get edges from default assay from Seurat object
#' default.assay <- Seurat::DefaultAssay(S)
#' edges = S@graphs[[which(grepl(paste0(default.assay, "_nn"), names(S@graphs)))]]
#' 
#' # get distance matrix
#' D = CID.GetDistMat(edges)
#' 
#' # entropy-based unclassified labels labels
#' entropy = CID.entropy(ac = P$L2, distM = D)
#' }
CID.entropy <- function(ac,distM)
{
  # "Y" labels will be ac but with some "Unclassified" cells
  Y   = ac
  
  # Calculate normalized shannon entropy for each cell j for all connections with shortest path < N
  shannon = rep(0, length(Y))
  if (class(distM) == "list")
  {
    dM = Reduce('+', distM) > 0;
    # create identity matrix to subtract off self
    I <- methods::new("ngTMatrix", 
                      i = as.integer(1:nrow(dM) - 1L), 
                      j = as.integer(1:nrow(dM) - 1L),
                      Dim = as.integer(c(nrow(dM), nrow(dM))))
    dM = dM - I
  } else {
    dM = distM
  }
  
  # create matrix for cell labels 
  m = Matrix::Matrix(0, nrow = length(ac), ncol = length(unique(ac)), sparse = TRUE)
  m[cbind(1:nrow(m), as.numeric(factor(ac)))] <- 1
  
  res = dM %*% m
  res = res / Matrix::rowSums(res)
  N_unique = length(unique(Y))
  shannon = apply(res, 1, function(freqs) {-sum(freqs[freqs != 0] * log(freqs[freqs != 0])) / log(2) / log2(N_unique)})
  q = c(shannon, -shannon)
  Y[shannon > 2 * sd(q)] = "Unclassified"
  
  return(Y)
}

#' Writes JSON file for SPRING integration
#' 
#' \code{CID.writeJSON} is a SPRING-integrated function for writing a JSON file.
#'
#' @param cr output from \code{\link{GenerateLabels}}
#' @param json_new Filename where annotations will be saved for new SPRING color tracks. Default is "categorical_coloring_data_new.json".
#' @param spring.dir Directory where file 'categorical_coloring_data.json' is located. If set, will add Signac annotations to the file.
#' @param new_populations Character vector specifying any new cell types that Signac has learned. Default is NULL.
#' @param new_colors Character vector specifying the HEX color codes for new cell types. Default is NULL.
#' @return A categorical_coloring_data.json file with Signac annotations and Louvain clusters added.
#' @export
#' See \url{https://github.com/AllonKleinLab/SPRING_dev} for more details.
CID.writeJSON <- function(cr, json_new = "categorical_coloring_data.json", spring.dir, new_populations = NULL, new_colors = NULL)
{
  if (!is.null(spring.dir))
  {
    spring.dir = gsub("\\/$", "", spring.dir, perl = TRUE);
    if (file.exists(paste(spring.dir, 'categorical_coloring_data.json', sep = "/")))
    {
      json_file = 'categorical_coloring_data.json'
      gJ <- paste(spring.dir,json_file,sep = "/")
      json_data <- RJSONIO::fromJSON(gJ)
      json_out = jsonlite::toJSON(json_data, auto_unbox = TRUE)
      write(json_out,paste(spring.dir,'categorical_coloring_data_legacy.json',sep="/"))
    } else {
      json_data = list("")
    }
  }
  if ("clusters" %in% names(cr))
  {
    Q = as.character(cr$clusters)
    json_data$Clusters_Louvain$label_list = Q
    Ntypes = length(unique(Q))
    qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    col_palette <- as.list(col_vector[1:Ntypes]);
    names(col_palette) <- unique(Q)
    json_data$Clusters_Louvain$label_colors = col_palette
  }
  if ("CellTypes" %in% names(cr))
  {
    Q = cr$CellTypes
    json_data$CellTypes$label_list = Q
    C = get_colors(Q)
    json_data$CellTypes$label_colors = as.list(C[[1]])
  }
  if ("CellStates" %in% names(cr))
  {
    Q = cr$CellStates
    json_data$CellStates$label_list = Q
    C = get_colors(Q)
    json_data$CellStates$label_colors = as.list(C[[1]])
  }
  if ("CellTypes_novel" %in% names(cr))
  {
    Q = cr$CellTypes_novel
    json_data$CellTypes_novel$label_list = Q
    C = get_colors(Q)
    Ntypes = sum(C[[1]] == "")
    qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    C[[1]][C[[1]] == ""] <- col_vector[1:Ntypes];
    json_data$CellTypes_novel$label_colors = as.list(C[[1]])
  }
  if ("CellStates_novel" %in% names(cr))
  {
    Q = cr$CellStates_novel
    json_data$CellStates_novel$label_list = Q
    C = get_colors(Q)
    Ntypes = sum(C[[1]] == "")
    qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    C[[1]][C[[1]] == ""] <- col_vector[1:Ntypes]; # or sample if you wish
    json_data$CellStates_novel$label_colors = as.list(C[[1]])
  }
  if ("Immune" %in% names(cr))
  {
    Q = cr$Immune
    json_data$Immune$label_list = Q
    C = get_colors(Q)
    json_data$Immune$label_colors = as.list(C[[1]])
  }
  json_data = json_data[order(names(json_data))]
  json_data_backup = json_data;
  if (!is.null(new_populations)) {
    json_data = lapply(json_data, function(x) {
        for (j in 1:length(new_populations))
          x$label_colors[names(x$label_colors) == new_populations[j]] = new_colors[j]
        x
    })
  }
  fn = json_new
  json_out = jsonlite::toJSON(json_data, auto_unbox = TRUE)
  write(json_out,paste(spring.dir,fn,sep="/"))
  cat(paste(spring.dir,fn,sep="/"), "has been written to directory! \n")
  return(json_data)
}

#' Detects community substructure by Louvain community detection
#'
#' \code{CID.Louvain} determines Louvain clusters from an edge list.
#'
#' @param edges A data frame or matrix with edges
#' @return A character vector with the community substructures of the graph corresponding to Louvain clusters.
#' @export
#' @examples
#' \dontrun{
#' # Loads edges
#' file.dir = "https://kleintools.hms.harvard.edu/tools/client_datasets/"
#' file = "CITESEQ_EXPLORATORY_CITESEQ_5K_PBMCS/FullDataset_v1_protein/edges.csv"
#' download.file(paste0(file.dir, file, "?raw=true"), destfile = "edges.csv")
#' 
#' # data.dir is your path to the "edges.csv" file
#' edges = CID.LoadEdges(data.dir = ".")
#' 
#' # get louvain clusters from edge list
#' clusters = CID.Louvain(edges)
#' }
CID.Louvain <- function(edges)
{
  # convert edgelist to adjacency matrix
  if (nrow(edges) != ncol(edges))
  {
    adjmatrix <- methods::new("ngTMatrix", 
                              i = c(as.integer(edges$V1)-1L, as.integer(edges$V2)-1L), 
                              j = c(as.integer(edges$V2)-1L, as.integer(edges$V1)-1L),
                              Dim = as.integer(c(max(edges), max(edges))))
  } else {
    adjmatrix = edges; # just load straight adjacency matrix 
  }
  g = igraph::graph_from_adjacency_matrix(adjmatrix * 1, mode = c("undirected"), weighted = NULL, diag = TRUE, add.colnames = NULL, add.rownames = NA)
  as.character(igraph::cluster_louvain(g)$membership)
}

#' Computes distance matrix from edge list
#' 
#' \code{CID.GetDistMat} returns the distance matrix (i.e., adjacency matrix, second-degree adjacency matrix, ..., etc.)
#' from an edge list. Here, edges is the edge list, and n is the order of the connetions.
#'
#' @param edges a data frame with two columns; V1 and V2, specifying the cell-cell edge list for the network
#' @param n maximum network distance to subtend (n neighbors)
#' @return adjacency matrices for distances < n
#' @export
#' @examples
#' \dontrun{
#' # Loads edges
#' file.dir = "https://kleintools.hms.harvard.edu/tools/client_datasets/"
#' file = "CITESEQ_EXPLORATORY_CITESEQ_5K_PBMCS/FullDataset_v1_protein/edges.csv"
#' download.file(paste0(file.dir, file, "?raw=true"), destfile = "edges.csv")
#' 
#' # data.dir is your path to the "edges.csv" file
#' edges = CID.LoadEdges(data.dir = ".")
#' 
#' # get distance matrix (adjacency matrices with up to fourth-order connections) from edge list
#' distance_matrices = CID.GetDistMat(edges)
#' }
CID.GetDistMat <- function(edges, n = 4)
{
  "%^%" <- function(A, n) {if(n == 1) A else A %*% (A %^% (n-1)) }
  # create adjacency matrix
  if (nrow(edges) != ncol(edges))
  {
    m <- methods::new("ngTMatrix", 
                      i = c(as.integer(edges$V1)-1L, as.integer(edges$V2)-1L), 
                      j = c(as.integer(edges$V2)-1L, as.integer(edges$V1)-1L),
                      Dim = as.integer(c(max(edges), max(edges))))
  } else {
    m = edges
  }
  
  dm = list("") # initialize distance matrix
  for (j in 1:n)
    if(j == 1) dm[[j]] = m else dm[[j]] = m %^% j
  return(dm)
}

#' Extracts unique elements
#'
#' \code{CID.IsUnique} returns a Boolean for the unique elements of a vector.
#'
#' @param x A character vector 
#' @return boolean, unique elements are TRUE
#' @export
#' @examples
#' # generate a dummy variable
#' dummy = c("A", "A", "B", "C")
#' 
#' # get only unique elements
#' logik = CID.IsUnique(dummy) 
#' dummy[logik]
CID.IsUnique <- function (x) 
{
  rv = rep(TRUE, length(x))
  if (length(x) >= 2) {
    ord = order(x)
    ox = x[ord]
    neq = (ox[-length(ox)] != ox[-1])
    rv[ord] = c(neq, TRUE) & c(TRUE, neq)
  }
  return(rv)
}

#' Save count_matrix.h5 files for SPRING integration
#' 
#' To integrate with SPRING, \code{SaveCountsToH5} saves expression matrices
#' in a sparse ".h5" format to be read with SPRING notebooks in Jupyter.
#'
#' @param D A list of count matrices
#' @param data.dir directory (will be created if it does not exist) where results are saved
#' @param genome default is GRCh38.
#' @importFrom methods as
#' @importClassesFrom Matrix dgCMatrix
#' @return matrix.h5 file, where each is background corrected
#' @export
#' @importFrom rhdf5 h5createFile h5createGroup h5write h5write.default
#' @examples
#' \dontrun{
#' # download single cell data for classification
#' file.dir = "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/"
#' file = "pbmc_1k_v3_filtered_feature_bc_matrix.h5"
#' download.file(paste0(file.dir, file), "Ex.h5")
#' 
#' # load data
#' library(Seurat)
#' E = Read10X_h5(filename = "Ex.h5")
#' 
#' # save counts to h5 files
#' SaveCountsToH5(E, data.dir = "counts_h5")
#' 
#' }
SaveCountsToH5 <- function(D, data.dir, genome = "GRCh38")
{
  data.dir = gsub("\\/$", "", data.dir, perl = TRUE);
  
  if (!dir.exists(data.dir))
    dir.create(data.dir)
  
  if (!"list" %in% class(D))
    D = list(D)
  
  cat(" ..........  Entry in SaveCountsToH5: \n");
  ta = proc.time()[3];
  cat("             Number of cells:", sum(sapply(D, ncol)), "\n");
  cat("             Number of genes:", nrow(D[[1]]), "\n");
  cat("             Number of samples:", length(D), "\n");
  
  if (is.null(names(D)))
    names(D) <- paste0("Sample", seq_along(1:length(D)))
  
  flag = rownames(D[[1]])[1] != gsub( "_.*$", "", rownames(D[[1]])[1])
  
  if (flag)
  {
    genes = strsplit(x = rownames(D[[1]]), split = "_._")
    gene_names = sapply(genes, `[[`, 1)
    genes = sapply(genes, `[[`, 2)
  }
  
  D = lapply(D, function(x)
  {
    x@Dimnames[[1]] = rownames(D[[1]])
    x
  })
  
  data.dirs = paste(data.dir, names(D), sep = "/")
  
  q = lapply(data.dirs, function(x){
    if (!dir.exists(x))
      dir.create(x)})
  
  d = suppressWarnings(
    mapply(function(x,y){
      fn = "count_matrix.h5"
      fn = paste(y, fn, sep = "/")
      rhdf5::h5createFile(fn)
      rhdf5::h5createGroup(fn, genome)
      rhdf5::h5write.default(x@Dimnames[[2]], file = fn, name = paste0(genome, "/barcodes"))
      if (flag)
      {
        rhdf5::h5write.default(genes, file = fn, name=paste0(genome, "/genes"))
        rhdf5::h5write.default(gene_names, file = fn, name=paste0(genome, "/gene_names"))
      } else {
        rhdf5::h5write.default(x@Dimnames[[1]], file = fn, name=paste0(genome, "/genes"))
        rhdf5::h5write.default(x@Dimnames[[1]], file = fn, name=paste0(genome, "/gene_names"))
      }
      rhdf5::h5write.default(x@x, file = fn, name=paste0(genome, "/data"))
      rhdf5::h5write.default(dim(x), file = fn, name=paste0(genome, "/shape"))
      rhdf5::h5write.default(x@i, file = fn, name=paste0(genome, "/indices")) # already zero-indexed.
      rhdf5::h5write.default(x@p, file = fn, name=paste0(genome, "/indptr"))
    }, x = D, y = data.dirs)
  )
  rhdf5::h5closeAll()
  tb = proc.time()[3] - ta;
  cat(" ..........  Exit SaveCountsToH5 \n");
  cat("             Execution time = ", tb, " s.\n", sep = "");
}

#' Get HEX colors
#' 
#' @param P A character vector of cell type labels
#' @return A vector of HEX colors for each category of cell type labels
#' @keywords internal
get_colors <- function(P)
{
  P <- sort(unique(as.character(P)));
  colorcount = length(unique(P))
  cs = character(colorcount)
  # main cell type categories will be consistently labeled:
  main_types = c("B"     ,     "Epithelial", "Fibroblasts"        ,
                 "MPh"         , "Plasma.cells"  , "Endothelial",
                 "TNK"         ,     "Unclassified"     , "NonImmune"          ,
                 "Immune")
  main_colors = c("#aa3596"     ,     "#387f50"   , "#a2ba37" ,
                  "#f1ff51"     ,     "#d778de"   , "#73de97" ,
                  "#90c5f4"     ,     "#c0c0c0"   , "#7fc97f" ,
                  "#bb7fc7" )
  
  # sub cell types will be consistently labeled:
  sub_types  = c("B.memory", "B.naive", "DC"    , "Macrophages", "Mon.Classical", "Mon.NonClassical", "Monocytes", "Neutrophils", "NK", "T.CD4.memory", "T.CD4.naive", "T.CD8.cm", "T.CD8.em", "T.CD8.naive", "T.gd", "T.regs" )
  sub_colors  = c("#aa3596", "#edc5e6", "#f9a702" , "#f97501", "#d6b171", "#9e4c05", "#f1ff51", "#f7f2b2", "#ac9bf2", "#bb7fc7", "#8a5e83", "#c9c9ff", "#8c9fe1", "#90c5f4", "#64058a", "#2038b0" )
  
  for (i in 1:length(sub_types))
  {
    cs[P == sub_types[i]] = sub_colors[i]
    cs[P == main_types[i]] = main_colors[i]
  }
  
  cs[P == "None"] = "#808080"
  
  outs = list(cs)
  
  names(outs[[1]]) <- P
  
  return(outs)
}

#' Load data file from directory
#' 
#' Loads 'matrix.mtx' and 'genes.txt' files from a directory.
#'
#' @param data.dir Directory containing matrix.mtx and genes.txt.
#' @param mfn file name; default is 'matrix.mtx'
#' @return A sparse matrix with rownames equivalent to the names in genes.txt
#' @export
#' @examples
#' \dontrun{
#' # Loads data from SPRING
#' 
#' # dir is your path to the "categorical_coloring_data.json" file
#' dir = "./FullDataset_v1"
#' 
#' # load expression data
#' E = CID.LoadData(data.dir = dir)
#' }
CID.LoadData <- function(data.dir, mfn = "matrix.mtx")
{
  data.dir = gsub("\\/$", "", data.dir, perl = TRUE);
  if (! (file.exists(paste(data.dir, mfn, sep = "/")) & file.exists(paste(data.dir, "genes.txt", sep = "/"))))
    data.dir = dirname(data.dir)
  gE <- paste(data.dir,mfn,sep="/")
  flag = file.exists(gE);
  if (!flag) {
    cat("ERROR: from CID.LoadData:\n");
    cat("file = ", gE, " does not exist.\n", sep = "");
    stop()
  }
  E <- Matrix::readMM(gE)
  # read genes
  fn <-"genes.txt"
  gG <- paste(data.dir,fn, sep = "/")
  flag = file.exists(gG);
  if (!flag) {
    cat("ERROR: from CID.LoadData:\n");
    cat("file = ", gG, " does not exist.\n", sep = "");
    stop()
  }
  genes <- read.delim(gG, stringsAsFactors = FALSE, header = FALSE)$V1
  if (genes[1] != gsub( "_.*$", "", genes[1] ))
  {
    genes[!grepl("Total", genes)] = gsub( "_.*$", "", genes [!grepl("Total", genes)])
    genes[grepl("Total", genes)]  = gsub( "_", "-",   genes [grepl("Total", genes)])
  }
  if (any(grepl(pattern = "-ENSG00", x = genes)))
  {
    genes <-gsub(pattern = "-ENSG00", replacement = "_ENSG00",x = genes, fixed = TRUE)
    genes[!grepl("Total", genes)] = gsub( "_.*$", "", genes [!grepl("Total", genes)])
    genes[grepl("Total", genes)]  = gsub( "_", "-",   genes [grepl("Total", genes)])
  }
  
  if (grepl("^1", genes[1]))
    genes = do.call(rbind, strsplit(genes, " "))[,2]
  flag = length(genes) %in% c(nrow(E), ncol(E));
  if (!flag) {
    cat("ERROR: from CID.LoadData:\n");
    cat("length of genes in genes.txt = ", length(genes), " is not equal to nrow(E) = ", nrow(E), ", or ncol(E) = ", ncol(E), "\n", sep = "");
    stop()
  }
  if (nrow(E) != length(genes))
    E = Matrix::t(E)
  rownames(E) <- genes
  E
}