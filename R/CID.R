#' Load 10X count matrix from an h5 file
#'
#' @param filename Directory and filename of the h5 file
#' @param use.names Boolean TRUE, see ?Seurat::Read10X_h5
#' @return Count matrix with genes and barcodes
#' @export
CID.Read10Xh5 <- function (filename, use.names = TRUE)
{
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Please install hdf5r to read HDF5 files")
  }
  if (!file.exists(filename)) {
    stop("File not found")
  }
  
  infile <- hdf5r::H5File$new(filename)
  genomes <- names(infile)
  
  if ("matrix" %in% genomes)
  {
    if (!infile$attr_exists("PYTABLES_FORMAT_VERSION")) {
      if (use.names) {
        feature_slot <- "features/name"
      }
      else {
        feature_slot <- "features/id"
      }
    } else {
      if (use.names) {
        feature_slot <- "name"
        feature_slot2 = "feature_type"
      }
      else {
        feature_slot = "name"
      }
    }
  } else {
    if (!infile$attr_exists("PYTABLES_FORMAT_VERSION")) {
      if (use.names) {
        feature_slot <- "features/name"
      }
      else {
        feature_slot <- "features/id"
      }
      if (basename(filename) == "count_matrix.h5")
      {
        feature_slot <- "genes"
      }
      
    } else {
      if (use.names) {
        feature_slot <- "gene_names"
        feature_slot2 = "genes"
      }
      else {
        feature_slot = "genes"
      }
    }
  }
  
  output <- list()
  
  for (genome in genomes) {
    counts <- infile[[paste0(genome, "/data")]]
    indices <- infile[[paste0(genome, "/indices")]]
    indptr <- infile[[paste0(genome, "/indptr")]]
    shp <- infile[[paste0(genome, "/shape")]]
    features <- infile[[paste0(genome, "/", feature_slot)]][]
    barcodes <- infile[[paste0(genome, "/barcodes")]]
    sparse.mat <- Matrix::sparseMatrix(i = indices[] + 1, p = indptr[],
                                       x = as.numeric(counts[]), dims = shp[], giveCsparse = FALSE)
    rownames(sparse.mat) <- features
    colnames(sparse.mat) <- barcodes[]
    sparse.mat <- as(object = sparse.mat, Class = "dgCMatrix")
    output[[genome]] <- sparse.mat
  }
  infile$close_all()
  if (length(output) == 1) {
    return(output[[genome]])
  }
  else {
    return(output)
  }
}

#' Load data file from directory
#'
#' @param data.dir Directory containing matrix.mtx and genes.txt.
#' @return A sparse matrix with rownames equivalent to the names in genes.txt
#' @export
CID.LoadData <- function(data.dir)
{
  data.dir = gsub("\\/$", "", data.dir, perl = TRUE);
  if (! (file.exists(paste(data.dir, "matrix.mtx", sep = "/")) & file.exists(paste(data.dir, "genes.txt", sep = "/"))))
    data.dir = dirname(data.dir)
  gE <- paste(data.dir,"matrix.mtx",sep="/")
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
  genes <- read.delim(gG, stringsAsFactors = F, header = F)$V1
  if (genes[1] != gsub( "_.*$", "", genes[1] ))
    genes = gsub( "_.*$", "", genes )
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

#' Load imputed data file from directory
#'
#' @param imputed.data.dir Directory containing matrix.mtx and genes.txt.
#' @param type indicates the imputation method used. Default is "ksoft", the other option is "saver".
#' @return A sparse matrix with rownames equivalent to the names in genes.txt
#' @export
CID.LoadImputedData <- function(imputed.data.dir, type = "ksoft")
{
  imputed.data.dir = gsub("\\/$", "", imputed.data.dir, perl = TRUE);
  gE <- paste(imputed.data.dir,paste0("matrix_",type,"_imputed.mtx"), sep="/")
  E <- Matrix::readMM(gE)
  # read genes
  fn <-"genes_saver_imputed.txt"
  fn <- paste0("genes_",type,"_imputed.txt")
  gG <- paste(imputed.data.dir,fn, sep = "/")
  flag = file.exists(gG);
  if (!flag) {
    cat("ERROR: from CID.LoadData:\n");
    cat("file = ", gG, " does not exist.\n", sep = "");
    stop()
  }
  genes <- read.delim(gG, stringsAsFactors = F, header = T)$x
  if (genes[1] != gsub( "_.*$", "", genes[1] ))
    genes = gsub( "_.*$", "", genes )
  if (grepl("^1", genes[1]))
    genes = do.call(rbind, strsplit(genes, " "))[,2]
  flag = length(genes) %in% c(nrow(E), ncol(E));
  if (!flag) {
    cat("ERROR: from CID.LoadImputedData:\n");
    cat("length of genes = ", length(genes), " is not equal to nrow(E) = ", nrow(E), ", or ncol(E) = ", ncol(E), "\n", sep = "");
    stop()
  }
  if (nrow(E) != length(genes))
    E = Matrix::t(E)
  rownames(E) <- genes
  E
}

#' Load chunked imputed data files from directory
#'
#' @param chunk.dir As defined in ?CID.chunk
#' @return A sparse matrix with rownames equivalent to the names in genes.txt
#' @export
CID.LoadImputedDataFolder <- function(chunk.dir)
{
  chunk.dir = gsub("\\/$", "", chunk.dir, perl = TRUE)
  fls = list.files(chunk.dir, full.names = TRUE)
  idx = as.integer(do.call(rbind, strsplit(fls, split = "Chunk"))[,2])
  fls = fls[order(idx)]
  E = lapply(fls, CID.LoadImputedData)
  E = do.call(cbind,E)
  return(E)
}

#' Load edges from edge list
#'
#' @param data.dir A directory where "edges.csv" is located
#' @return The edgelist in data frame format
#' @export
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
  edges <- read.delim(edges,sep = ";" , stringsAsFactors = F, header = F)
  if(min(edges) == 0)
    edges = edges + 1
  colnames(edges) <- c("V1", "V2")
  edges
}

#' Library size normalize
#'
#' @param E Expression matrix
#' @return Normalized expression matrix to mean of total counts
#' @export
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

#' Chunk a dataset
#'
#' @param E A gene-by-sample count matrix (sparse matrix, matrix, or data.frame).
#' @param chunk.dir A directory where the chunked matrices will be stored. If the directory does not exist, it will be created.
#' @param number_of_chunks the number of chunks (e.g., the number of sub-sampled matrices after chunking)
#' @return A directory with chunked matrix
#' @export
#'
CID.Chunk <- function(E, chunk.dir, number_of_chunks = 30)
{
  chunk.dir = gsub("\\/$", "", chunk.dir, perl = TRUE);
  E = array_split(E, number_of_chunks = number_of_chunks)
  if (is.null(chunk.dir))
  {
    cat("Please specify a chunk.dir. \n")
    stop()
  }
  if (!dir.exists(chunk.dir))
    dir.create(chunk.dir)
  names(E) <- paste0("Chunk", seq_along(1:length(E)))
  large.dirs = paste(chunk.dir, names(E), sep = "/")
  q = lapply(large.dirs, function(x){
    if (!dir.exists(x))
      dir.create(x)})
  mapply(function(x,y){
    Matrix::writeMM(Matrix::Matrix(x, sparse = TRUE), file = paste(y, "matrix.mtx", sep = "/"))
    write.table(rownames(x), file = paste(y, "genes.txt", sep = "/"), col.names = F)
  }, x = E, y = large.dirs)
}

#' Imputation wrapper
#'
#' @param E A gene-by-sample count matrix (sparse matrix, matrix, or data.frame) with genes identified by their HUGO symbols.
#' @param data.dir directory for saving "matrix_saver_imputed.mtx" for future loading.
#' @param do.par Boolean. If true, imputation is performed in parallel on half of the machines available cores. Default = FALSE.
#' @return imputed expression matrix with only marker genes in rows.
#' @export
#'
CID.Impute <- function(E, data.dir = NULL, do.par = TRUE)
{
  if (ncol(E) > 80000)
    return("Chunk first, and then run on the chunks. See ?CID.chunk. \n")
  markers = Signac::markers
  cellstate_markers = Signac::cellstate_markers
  load("/site/ne/data/bh-results/C/CHAMBERLAIN.Mat/pipelines/Signac/data/immune_markers.rda")
  genes = do.call(rbind, cellstate_markers)
  genes.ind <- which(rownames(E) %in% unique(c(as.character(markers$`HUGO symbols`), as.character(genes$`HUGO symbols`), as.character(immune_markers$`HUGO symbols`))))
  if (do.par)
  { 
    numCores = parallel::detectCores()
  } else {
    numCores = 2;
  }
  I = SAVER::saver(E, pred.genes = genes.ind, pred.genes.only = T, ncores = numCores - 1, estimates.only = T)
  if (!is.null(data.dir))
  {
    Matrix::writeMM(Matrix::Matrix(I, sparse = TRUE), file = paste(data.dir, "matrix_saver_imputed.mtx", sep = "/"))
    write.table(rownames(E)[genes.ind], file = paste(data.dir, "genes_saver_imputed.txt", sep = "/"))
  }
  return(I)
}

#' K soft imputation
#'
#' @param E A gene-by-sample count matrix (sparse matrix, matrix, or data.frame) with genes identified by their HUGO symbols (see ?CID.geneconversion), or a list of such matrices, see ?CID.BatchMode.
#' @param data.dir directory of SPRING files "edges.csv" and "categorical_coloring_data.json"
#' @param genes a character vector of genes to impute
#' @return Imputed values for every gene
#' @export
KSoftImpute <- function(E,  dM = NULL, genes.to.use = NULL, do.save = F)
{
  # check inputs
  stopifnot(class(E) %in% c("dgCMatrix","dgTMatrix", "matrix", "data.frame"))
  stopifnot(!is.null(rownames(E)));
  if (class(E) %in% c("matrix", "data.frame"))
    E = Matrix::Matrix(as.matrix(E), sparse = T)
  
  cat(" ..........  Entry in KSoftImpute \n");
  ta = proc.time()[3];
  
  # main function
  cat(" ..........  Running K soft imputation on input data matrix :\n");
  cat("             nrow = ", nrow(E), "\n", sep = "");
  cat("             ncol = ", ncol(E), "\n", sep = "");
  
  if (!is.null(genes.to.use))
    E = E[rownames(E) %in% genes.to.use,]
  
  # define the weight matrix Wjj, which is the number of genes detected in each cell
  weights = Matrix::colSums(E != 0)
  bas = Matrix::Matrix(0, nrow = length(weights), ncol = length(weights))
  diag(bas) <- weights
  
  # if is null data.dir, run PCA + KNN
  if (is.null(dM))
  {
    dM = CID.GetNeighbors(E, normalize = normalize, min_counts = 3, min_cells = 3, min_vscore_pctl = 85, num_pc = 30, k_neigh = 3)
    dM = CID.GetDistMat(dM)
  }
  
  # g = dM[[1]] %*% bas + 1/2 * dM[[2]] %*% bas
  g = dM[[1]] %*% bas
  dd = g / Matrix::rowSums(g)
  diag(dd) <- 1
  E_new = E %*% dd;
  E_new = CID.Normalize(E_new)
  
  if (do.save)
  {
    data.dir = dirname(data.dir)
    Matrix::writeMM(Matrix::Matrix(E_new, sparse = TRUE), file = paste(data.dir, "matrix_ksoft_imputed.mtx", sep = "/"))
    write.table(rownames(E_new), file = paste(data.dir, "genes_ksoft_imputed.txt", sep = "/"))
  }
  
  tb = proc.time()[3] - ta;
  cat("\n ..........  Exit KSoftImpute.\n");
  cat("             Execution time = ", tb, " s.\n", sep = "");
  
  return(E_new)
}

#' Check input classes
#'
#' @return Checks for function calls
#' @export
RunChecks <- function(...)
{
  stopifnot(class(E) %in% c("dgCMatrix","dgTMatrix", "matrix", "data.frame"))
  stopifnot(!is.null(rownames(E)))
}

#' Get edges that are either pre-computed, or generate new edges
#'
#' @return edges for cell-cell similarity network
#' @export
CID.GetEdges <- function(...)
{
  if (is.null(data.dir))
  {
    edges = CID.GetNeighbors(E, normalize = T, min_counts = 3, min_cells = 3, min_vscore_pctl = 85, num_pc = 30, k_neigh = 4)
  } else {
    edges = CID.LoadEdges(data.dir)
  }
  return(edges)
}

#' Main function
#'
#' @param E A gene-by-sample count matrix (sparse matrix, matrix, or data.frame) with genes identified by their HUGO symbols (see ?CID.geneconversion), or a list of such matrices, see ?CID.BatchMode.
#' @param reference if default, uses the standard Signac markers (see ?CID.SeeMarkers).
#' @param full.dataset If E was subsetted or imputed, full.dataset is the full expression matrix, which is loaded for the detection of novel cell types. Default is NULL.
#' @param pval p-value cutoff for feature selection, as described in the manuscript + markdown file. Default is pval = 0.05.
#' @param data.dir directory of SPRING files "edges.csv" and "categorical_coloring_data.json"
#' @param entropy cells amended to high entropy labels with respect to their neighbors in the KNN graph are appended "Other" if entropy = TRUE. Default is TRUE.
#' @param omit Force remove specific cell types / states. Default is NULL.
#' @param force.keep Force keep specific cell types / states. Default is NULL.
#' @param sorted If cells are expected to be pure or mostly homogeneous (e.g., by FACs sorting), set sorted = TRUE. Default is FALSE.
#' @param do.impute Run SAVER or Ksoft imputation. Default is TRUE.
#' @param do.log log2 transform counts + 1. Default is TRUE.
#' @param min.cells minimum number of cells to call a discrete population. Default is 30.
#' @param nonimmune flag for whether we expect to see nonimmune cells in the tissue. Default is FALSE.
#' @return Filtered markers where each marker must have at least ncells that express at least ncounts
#' @export
CID.CellID <- function(E, reference = 'default', pval = 0.05, data.dir = NULL, entropy = T, omit = NULL, force.keep = NULL, sorted = F, do.impute = T, do.log = T, min.cells = 30, nonimmune = F, thold)
{
  # load markers
  if (class(reference) == 'character')
  {
    load("/site/ne/data/bh-results/C/CHAMBERLAIN.Mat/pipelines/Signac/data/markers.rda")
    load("/site/ne/data/bh-results/C/CHAMBERLAIN.Mat/pipelines/Signac/data/cellstate_markers.rda")
    load("/site/ne/data/bh-results/C/CHAMBERLAIN.Mat/pipelines/Signac/data/immune_markers.rda")
    all_markers = rbind(markers, do.call(rbind, cellstate_markers))
    nonimmune_markers = NULL
    reference = list("")
  } else {
    immune_markers = reference[[1]];
    markers = reference[[2]];
    cellstate_markers = reference[[3]];
    nonimmune_markers = reference[[4]];
    all_markers = rbind(markers, do.call(rbind, cellstate_markers))
  }
  
  # check inputs
  inputs = match.call()
  RunChecks(inputs)
  if (class(E) %in% c("matrix", "data.frame"))
    E = Matrix::Matrix(as.matrix(E), sparse = T)
  if (!is.null(data.dir))
    data.dir = gsub("\\/$", "", data.dir, perl = TRUE);
  
  # main function
  cat(" ..........  Entry in CID.CellID \n");
  ta = proc.time()[3];
  cat(" ..........  Computing Signac scores for cell types on input data matrix :\n");
  cat("             nrow = ", nrow(E), "\n", sep = "");
  cat("             ncol = ", ncol(E), "\n", sep = "");
  
  # normalize to the mean library size
  E = CID.Normalize(E)
  
  # user can omit any cell type or cell state
  if (!is.null(omit))
  {
    cat(" ..........  Forcibly omitting features :\n");
    all_markers = all_markers[! all_markers$`Cell population` %in% omit, ]
    markers = markers[! markers$`Cell population` %in% omit, ]
    cellstate_markers = cellstate_markers[! names(cellstate_markers) %in% omit]
    cellstate_markers = lapply(cellstate_markers, function(x){ x[! x$`Cell population` %in% omit,] })
    cat("             Omitted = ", paste0(omit, collapse = ", "), "\n");
  }
  
  # get edges for cell-cell similarity network
  edges = CID.GetEdges(inputs)
  
  # compute distance matrix
  distMat = CID.GetDistMat(edges = edges)
  
  # get Louvain clusters
  louvain = CID.Louvain(edges = edges)
  
  # we keep the full.dataset and segment the rest for efficiency
  full.dataset = E
  
  # user can let Signac auto-detect the presence of nonimmune cells
  if (nonimmune == "auto")
  {
    filtered_features=subset(immune_markers,get("HUGO symbols")%in%rownames(E))
    filtered_features=split(filtered_features[,c("HUGO symbols", "Cell population" ,"Polarity")], filtered_features[,"Cell population"])
    # compute cell type scores data.frame
    scores = CID.append(E,filtered_features, sorted = F)
    nonimmune = diptest::dip.test(c(as.numeric(scores[1,]), as.numeric(scores[2,])))$p.value < 0.2
  }

  # user can assert the presence of nonimmune cells
  if (nonimmune){
    all_markers = rbind(all_markers, nonimmune_markers)
    E = E[intersect(row.names(E), union(immune_markers$`HUGO symbols`, all_markers$`HUGO symbols`)),,drop=F]
  } else {
    E = E[intersect(row.names(E), all_markers$`HUGO symbols`),,drop=F]
  }
  
  # log2 + 1 transform
  if (do.log)
    E = Matrix::Matrix(log2(E + 1), sparse = T)
  
  # run KSoft imputation
  if (do.impute)
    E = KSoftImpute(E, dM = distMat)
  
  if (nonimmune)
  {
    filtered_features=subset(immune_markers,get("HUGO symbols")%in%rownames(E))
    filtered_features=split(filtered_features[,c("HUGO symbols", "Cell population" ,"Polarity")], filtered_features[,"Cell population"])
    # compute cell type scores data.frame
    scores = CID.append(E,filtered_features, sorted = F)
    indexMax = apply(scores, 2, which.max);
    immunetypes = rownames(scores)[indexMax];
    immunetypes = CID.smooth(immunetypes, distMat[[1]])
    df = table(immunetypes, louvain)
    indexMax = apply(df, 2, which.max)
    df = data.frame(lbls = colnames(df), type = rownames(df)[indexMax])
    immunetypes = as.character(df$type[match(louvain, df$lbls)])
    celltypes = immunetypes
    if ("nonimmune_markers" %in% names(reference)) {
      # get features
      filtered_features = CID.filter2(l = E, m  = nonimmune_markers, n = pval, o = FALSE)
      
      # compute cell type scores data.frame
      scores = CID.append(E, filtered_features, sorted = T)
      indexMax = apply(scores, 2, which.max);
      celltypes = immunetypes
      celltypes[immunetypes == "NonImmune"] <- rownames(scores)[indexMax][immunetypes == "NonImmune"];
    }
  } else {
    immunetypes = rep("Immune", ncol(E))
    immunestates = rep("Immune", ncol(E))
    celltypes = immunetypes
  }
  
  # get features
  filtered_features = CID.filter2(l = full.dataset, m = markers, n = pval, o = FALSE)
  
  if (!is.null(force.keep))
  {
    flag = ! force.keep %in% names(filtered_features)
    if (flag){
      q=subset(markers,get("HUGO symbols")%in%rownames(E))
      q=split(q[,c("HUGO symbols", "Cell population" ,"Polarity")], q[,"Cell population"])
      filtered_features = c(filtered_features, q[names(q) %in% force.keep])
      cat(" ..........  Forcibly keeping features :\n");
      cat("             Kept = ", paste0(names(filtered_features), collapse = ", "), "\n");
    }
  }
  
  # compute cell type scores data.frame
  scores = CID.append(E, filtered_features, sorted = T)
  
  # assign output classifications
  cat(" ..........  Assigning output classifications \n", sep ="");
  indexMax = apply(scores, 2, which.max);
  celltypes[immunetypes == "Immune"] <- rownames(scores)[indexMax][immunetypes == "Immune"];
  
  # amend low scoring states to "Other"
  kmax = apply(scores, 2 ,max)
  diff = kmax - apply(scores, 2, min)
  diff = diff[immunetypes == "Immune"]
  kmax = kmax[immunetypes == "Immune"]
  celltypes[immunetypes == "Immune"][diff < (mean(diff) - 2 * sd(diff))] = "Other"
  celltypes[immunetypes == "Immune"][kmax < 0] = "Other"
  
  # smooth the output classifications
  celltypes = CID.smooth(celltypes, distMat[[1]])
  
  # assign Others
  if (entropy)
    celltypes = CID.entropy(celltypes, distMat)
  
  # assign any populations of less than min.cells to "Other"
  q = data.frame(table(celltypes))
  celltypes[celltypes %in% q$celltypes[q$Freq <= min.cells]] = "Other"
  
  # cell state deep dive classifications
  cat(" ..........  Computing CID scores for cell states! \n");
  cellstates = celltypes
  cellstates_merged = celltypes
  
  # get only cell states with cell types in the data
  logik = names(cellstate_markers) %in% cellstates;
  
  # get output list
  k = 0
  outs  = list("")
  outs2 = list("")
  
  while(any(logik))
  {
    # get filtered features for cell states
    state_features = mapply(function(x, y){
      CID.filter2(l = full.dataset, m = x, n = pval, o = F)
    }, x = cellstate_markers[logik], y = names(cellstate_markers)[logik], SIMPLIFY = F)
    
    # it is possible that all of a given cell type has no discernable cell states
    flag = sapply(state_features, is.null)
    if (any(flag)){
      cellstates[cellstates == names(state_features)[flag]] = "Other";
      state_features = state_features[!flag]
    }
    
    # if a cell state has only one present type, we default to that label for all cells
    flag = sapply(state_features, function(x) length(x)) == 1
    if (any(flag)){
      cellstates[cellstates == names(state_features)[flag]] = as.character(state_features[flag][[1]][[1]][1,2]);
      state_features = state_features[!flag]
    }
    
    # get scores for cell states
    scores = mapply(function(x,y){
      q = CID.append(E[,cellstates == y], x, sorted = F)
      if (cor(t(q[1,]), t(q[2,])) < -0.95)
        q = CID.append(E[,cellstates == y], x, sorted = T)
      q
    }, x = state_features, y = names(state_features), SIMPLIFY = F)
    
    for (j in 1:length(scores))
    {
      logik = cellstates == names(scores)[j]
      if (any(logik))
      {
        indexMax   = apply(scores[[j]], 2, which.max);
        actual.max = apply(scores[[j]], 2, max);
        actual.min = apply(scores[[j]], 2, min);
        cellstates[logik] = rownames(scores[[j]])[indexMax];
        # amend low scoring states to "Other"
        diff = actual.max - actual.min;
        cellstates[logik][diff < (mean(diff) - 2 * sd(diff))] = "Other"
        qq = lapply(distMat, function(x) x[logik,logik])
        cellstates[logik] = CID.smooth(rownames(scores[[j]])[indexMax], qq[[1]])  
        cellstates[logik] = CID.smooth(cellstates[logik], qq[[1]])
        cellstates[logik] = CID.smooth(cellstates[logik], qq[[1]])
        if (entropy)
          cellstates[logik] = CID.entropy(cellstates[logik], qq)
        
        if (length(unique(cellstates[logik])) > 1 & names(scores)[j] %in% cellstates_merged)
        {
          # do DEG testing
          dummy = full.dataset[,logik]
          colnames(dummy) <- cellstates[logik]
          pos = CID.PosMarkers3(dummy, threshold = thold)
          q = lapply(pos, function(p){do.call(rbind, p)})
          pre = unique(colnames(dummy))
          
          # any categories with no differentially expressed genes are merged
          logik2 = sapply(q, function(x) {x$ident.1[1] == "NONE"})
          
          if (any(logik2)){
            cellstates_merged[cellstates %in% names(logik2)[logik2]] = names(scores)[j]
          } else {
            cellstates_merged[logik] = cellstates[logik]
          }
        }
        
        k = k + 1
        outs[[k]] = cellstates[logik]
        outs2[[k]] = cellstates
      }
    }
    logik = names(cellstate_markers) %in% cellstates;
  }
  
  do = data.frame(table(louvain[cellstates == "Other"]))
  df = data.frame(table(louvain[louvain %in% do$Var1]))
  logik = (1 - phyper(do$Freq, df$Freq , length(cellstates) - do$Freq, sum(cellstates == "Other"))) < 0.01;
  # logik = F
  if (any(logik)){
    do = do[logik,]
    logik = do$Freq > 3; # require at least 3 cell communities
    if (any(logik)){
      celltypes_novel  = celltypes;
      cellstates_novel = cellstates;
      cat("             Signac found", sum(logik), "novel celltypes!\n");
      lbls = rep("All", ncol(E))
      logik = louvain %in% do$Var1[logik] & cellstates == "Other";
      lbls[logik] = louvain[logik]
      colnames(full.dataset) <- lbls
      # FF = apply(full.dataset, 1, sd)^2 / Matrix::rowMeans(full.dataset)
      # full.dataset = full.dataset[FF > median(FF), ]
      new_lbls = CID.PosMarkers2(full.dataset, cellstates)
      cellstates_novel = new_lbls$lbls
      celltypes_novel[grepl("^[+]", cellstates_novel)] = new_lbls$lbls[grepl("^[+]", cellstates_novel)]
    } else {
      celltypes_novel  = celltypes
      cellstates_novel = cellstates
    }
  } else {
    celltypes_novel  = celltypes
    cellstates_novel = cellstates
  }
  
  # write nonimmune cell states
  if ("NonImmune" %in% immunetypes)
  { 
    immunestates = cellstates_novel
    immunestates[immunestates == "NonImmune"] = paste0(immunestates[immunestates == "NonImmune"], "_", louvain[immunestates == "NonImmune"])
  } else {
    immunestates = cellstates_novel
  }
  
  # Package output
  cr = list(celltypes = celltypes,
            celltypes_novel = celltypes_novel,
            cellstates = cellstates,
            cellstates_novel = cellstates_novel,
            cellstates_merged = cellstates_merged,
            immune = immunetypes,
            immunestates = immunestates,
            louvain = louvain)
  
  tb = proc.time()[3] - ta;
  cat("\n ..........  Exit CID.CellID.\n");
  cat("             Execution time = ", tb, " s.\n", sep = "");
  return (cr);
}

#' Learn new cell type markers
#'
#' @param E A gene-by-sample count matrix (sparse matrix, matrix, or data.frame) with genes identified by their HUGO symbols (see ?CID.geneconversion), or a list of such matrices, see ?CID.BatchMode.
#' @param cells.ind a named list of the indices of the cells, and the name of the putatitive label
#' @param reference if default, uses the standard Signac markers (see ?CID.SeeMarkers).
#' @return an updated reference set with the new cell type
#' @export
CID.Learn <- function(E, cells.ind, refererence = 'default')
{
  if (class(reference) == 'character')
    load("./data/LM22.rda")
  
  logik = CID.IsUnique(rownames(E))
  E = E[logik,]
  
  L = lapply(cells.ind,function(x){na.omit(E[, x])})
  M = do.call(cbind, L)
  
  # do cell type markers first
  logik = grepl(":", colnames(M))
  if (any(logik)){
    nms = do.call(rbind, strsplit(colnames(M)[logik], split = ":"))
    colnames(M)[logik] <- nms[,1]
  }
  
  S = LM22$data
  colnames(S) <- LM22$main_types
  
  gns = intersect(rownames(S), rownames(E))
  
  New_Ref = cbind(M[rownames(M) %in% gns,], S[rownames(S) %in% gns,])
  New_Ref = CID.Normalize(as.matrix(New_Ref))
  
  master_types = CID.GetLearningMarkers(E = New_Ref, full.dataset = E, logfc.threshold = 1)
  
  new_markers = rbind(markers, master_types)
  new_markers = unique(new_markers)
  new_markers = new_markers[order(new_markers$`Cell population`),]
  
  master_list = list( B.cells       =   list(B.cells.naive  = c("B.cells:B.cells.naive", "B_cell:Naive"),
                                             B.cells.memory = c("B.cells:B.cells.memory", "B_cell:Memory")) ,
                      MPh           =   list(Monocytes      = c("MPh:Monocytes",unique(colnames(S)[grepl("Monocyte", colnames(S))]), "Monocytes"),
                                             Macrophages    = c("MPh:Macrophages.M0"),
                                             Dendritic      = c("MPh:Dendritic.cells.activated","MPh:Dendritic.cells.resting", "DC:monocyte-derived",
                                                                unique(colnames(S)[grepl("DC:monocyte-derived:", colnames(S))]))),
                      Monocytes     =   list(Mon.Classical  = c("Monocyte:CD14+"),
                                             Mon.NonClass   = c("Monocyte:CD16+")),
                      TNK           =   list(T.cells        = c(unique(colnames(S)[grepl("^TNK:T", colnames(S))]),
                                                                unique(colnames(S)[grepl("T_cell", colnames(S))]),
                                                                "CD4+ T-cells", "CD8+ T-cells"),
                                             NK             = c("TNK:NK.cells.activated", "TNK:NK.cells.resting", "NK cells", "NK_cell:IL2")),
                      T.cells       =   list(T.cells.CD8    = c("TNK:T.cells.CD8", "CD8+ T-cells", "T_cell:CD8+_naive", "T_cell:CD8+"),
                                             T.CD4.FH.regs  = c("TNK:T.cells.CD4.memory.activated",
                                                                "TNK:T.cells.CD4.memory.resting", "TNK:T.cells.CD4.naive",
                                                                "TNK:T.cells.follicular.helper", "TNK:T.cells.regulatory..Tregs.",
                                                                "CD4+ T-cells", "T_cell:CD4+_effector_memory", "T_cell:CD4+_Naive",
                                                                "T_cell:CD4+", "T_cell:Treg:Naive", "T_cell:CD4+_central_memory")),
                      T.CD4.FH.regs =   list(T.regs         = c("TNK:T.cells.regulatory..Tregs.", "T_cell:Treg:Naive"),
                                             T.CD4.FH       = c("TNK:T.cells.CD4.memory.activated",
                                                                "TNK:T.cells.CD4.memory.resting", "TNK:T.cells.CD4.naive",
                                                                "TNK:T.cells.follicular.helper", "CD4+ T-cells", "T_cell:CD4+_central_memory")),
                      T.CD4.FH      =   list(T.cells.FH     = c("TNK:T.cells.follicular.helper"),
                                             T.cells.CD4    = c("TNK:T.cells.CD4.naive","TNK:T.cells.CD4.memory.activated",
                                                                "TNK:T.cells.CD4.memory.resting", "CD4+ T-cells", "T_cell:CD4+_central_memory")),
                      T.cells.CD4  =    list(T.cells.CD4n   = "TNK:T.cells.CD4.naive",
                                             T.cells.CD4m   = c("TNK:T.cells.CD4.memory.activated","TNK:T.cells.CD4.memory.resting", "T_cell:CD4+_central_memory")),
                      Granulocytes  =    list(Neutrophil = "Granulocytes:Neutrophils",
                                              Eosinophils    = "Granulocytes:Eosinophils")
  )
  
  logik = names(master_list) %in% nms[,1];
  logik2 = names(master_list[[which(logik)]]) %in% nms[,2]
  
  S = LM22$data
  colnames(S) <- LM22$types
  
  gns = intersect(rownames(S), rownames(E))
  
  New_Ref = cbind(M[rownames(M) %in% gns,], S[rownames(S) %in% gns,])
  New_Ref = CID.Normalize(as.matrix(New_Ref))
  
}

#' Learn new cell type markers
#'
#' @param E A gene-by-sample count matrix (sparse matrix, matrix, or data.frame) with genes identified by their HUGO symbols (see ?CID.geneconversion), or a list of such matrices, see ?CID.BatchMode.
#' @param cells.ind a named list of the indices of the cells, and the name of the putatitive label
#' @return an updated reference set with the new cell type
#' @export
CID.Learn2 <- function(E, cells.ind, markers = NULL)
{
  cat(" ..........  Learning and refining cell type markers on training subset: \n");
  cat("             cells used = ", length(unlist(cells.ind)), "\n", sep = "");
  
  # load default markers
  if (is.null(markers))
    load("./data/markers.rda")
  
  # get unique gene names only
  logik = CID.IsUnique(rownames(E))
  E = E[logik,]
  
  # get training set
  L = lapply(cells.ind,function(x){na.omit(E[, x])})
  M = do.call(cbind, L)
  
  # define markers based on Z-score transformation of log2 + 1 nUMI
  M = log2(M + 1)
  M = Matrix::t(scale(Matrix::t(M)))
  
  # remove any NA values
  logik = apply(M, 1, function(x) any(!is.na(x)))
  M = M[logik,]
  
  # split into train set into test + dev sets
  set.seed(42)
  train <- sample(ncol(M), 0.7*ncol(M), replace = FALSE)
  TrainSet <- M[,train]
  TestSet  <- M[,-train]
  
  # convert to numeric
  q = t(aggregate(t(M), list(colnames(M)), sum))
  colnames(q) <- q[1,]
  q = q[-1,]
  q = apply(q, 2, as.numeric)
  gns = rownames(M)
  
  # convert to marker format
  master_types = apply(q, 2, function(x) gns[order(x, decreasing = T)])
  
  # define function for optimizing K
  optim_K <- function(K, x = master_types, y = markers[markers$`Cell population` %in% colnames(master_types),], z = TrainSet)
  {
    cls = colnames(x)
    x = reshape2::melt(x[1:K,])
    x = data.frame(genes = x$value, celltype = x$Var2, polarity = "+")
    colnames(x) <- colnames(y)
    
    x = rbind(x, y)
    x = unique(x)
    x = x[sample(nrow(x), nrow(x), replace = T), ]
    x = unique(x)
    
    filtered_features=subset(x,get("HUGO symbols")%in%rownames(E))
    filtered_features=split(filtered_features[,c("HUGO symbols", "Cell population" ,"Polarity")], filtered_features[,"Cell population"])
    scores = CID.append(z, filtered_features, sorted = T)
    
    indexMax = apply(scores, 2, which.max);
    celltypes_train <- rownames(scores)[indexMax];
    acc = sum(celltypes_train == colnames(z)) / ncol(z)
    
    return(list(acc = acc, markers_used = x))
  }
  
  cat(" ..........  Training Signac classifier on training set: \n");
  cat("             nrow = ", nrow(TrainSet), "\n", sep = "");
  cat("             ncol = ", ncol(TrainSet), "\n", sep = "");
  
  q = sapply(2:30, function(x){
    N_markers = rep(x, 100)
    sapply(N_markers, optim_K)
  })
  
  accs = unlist(q[c(TRUE, FALSE)])
  mrks = q[c(FALSE, TRUE)]
  new_markers = mrks[which(accs == max(accs))]
  
  cat(" ..........  Done! Model accuracy on training set = ", max(accs), "\n");
  
  q_test = lapply(new_markers, function(x){
    filtered_features=subset(x,get("HUGO symbols")%in%rownames(TestSet))
    filtered_features=split(filtered_features[,c("HUGO symbols", "Cell population" ,"Polarity")], filtered_features[,"Cell population"])
    scores = CID.append(TestSet, filtered_features, sorted = T)
    
    indexMax = apply(scores, 2, which.max);
    celltypes_test <- rownames(scores)[indexMax];
    acc = sum(celltypes_test == colnames(TestSet)) / ncol(TestSet)
  })
  
  accs_test = unlist(q[c(TRUE, FALSE)])
  mrks_test = q[c(FALSE, TRUE)]
  
  cat(" ..........  Model accuracy on test set = ", max(accs_test), "\n");
  
  new_markers = mrks_test[which(accs_test == max(accs_test))]
  
  new_markers = lapply(new_markers, function(x) {
    x[order(x[,2]),]
  })
  
  scores_test = lapply(new_markers, function(x){
    x$`Cell population` = as.character(x$`Cell population`)
    filtered_features=subset(x,get("HUGO symbols")%in%rownames(TestSet))
    filtered_features=split(filtered_features[,c("HUGO symbols", "Cell population" ,"Polarity")], filtered_features[,"Cell population"])
    CID.append(TestSet, filtered_features, sorted = T)
  })
  
  cor_test = sapply(scores_test, function(x){
    d = cor(t(x))
    diag(d) <- 0
    sum(d)
  })
  
  cat(" ..........  Number of models with optimum performance = ", length(new_markers), "\n");
  
  new_markers = new_markers[[which.min(cor_test)]]
  
  return(new_markers)
  
}

#' Learn new cell type markers
#'
#' @param E A gene-by-sample count matrix (sparse matrix, matrix, or data.frame) with genes identified by their HUGO symbols (see ?CID.geneconversion), or a list of such matrices, see ?CID.BatchMode.
#' @return an updated reference set with the new cell type
#' @export
CID.Learn3 <- function(E, markers = NULL)
{
  cat(" ..........  Learning and refining cell type markers on training subset: \n");
  cat("             cells used = ", ncol(E), "\n", sep = "");
  
  # load default markers
  if (is.null(markers))
    load("./data/markers.rda")
  
  # get unique gene names only
  logik = CID.IsUnique(rownames(E))
  M = E[logik,]
  
  # split into train set into balanced test + dev sets, balanced by the minimum population
  if (length(unique(colnames(M))) != 2)
  {
    df = data.frame(table(colnames(M)))
    min.pop = min(df$Freq)
    dummy = M[,colnames(M) %in% df$Var1[df$Freq < 5 * min.pop]]
    dummy2 = M[,colnames(M) %in% df$Var1[df$Freq > 5 * min.pop]]
    balmat = list("")
    for (j in 1:length(unique(colnames(dummy2))))
      balmat[[j]] = dummy2[,which(colnames(dummy2) == unique(colnames(dummy2))[j])[sample((df$Freq[df$Var1 == unique(colnames(dummy2))[j]] / min.pop), replace = FALSE)]]
    balmat = do.call(cbind, balmat)
    M = cbind(dummy, balmat)
  }
  
  # log2 + 1 transform
  M = log2(M + 1)
  M = Matrix::Matrix(Matrix::t(scale(Matrix::t(M))), sparse = T)
  
  # remove any NA values
  logik = apply(M, 1, function(x) any(!is.na(x)))
  M = M[logik,]  
  
  set.seed(42)
  train <- sample(ncol(M), 0.7*ncol(M), replace = FALSE)
  TrainSet <- M[,train]
  TestSet  <- M[,-train]
  
  # convert to numeric
  q = Matrix::t(Matrix.utils::aggregate.Matrix(Matrix::t(TestSet), factor(colnames(TestSet)), fun = 'mean'))
  gns = rownames(M)
  
  # convert to marker format
  master_types = apply(q, 2, function(x) gns[order(x, decreasing = T)])
  
  # define function for optimizing K
  optim_K <- function(K, x = master_types, y = markers[markers$`Cell population` %in% colnames(master_types),], z = TrainSet)
  {
    cls = colnames(x)
    x = reshape2::melt(x[1:K,])
    x = data.frame(genes = x$value, celltype = x$Var2, polarity = "+")
    colnames(x) <- colnames(y)
    
    x = rbind(x, y)
    x = unique(x)
    x = x[sample(nrow(x), nrow(x), replace = T), ]
    x = unique(x)
    
    filtered_features=subset(x,get("HUGO symbols")%in%rownames(E))
    filtered_features=split(filtered_features[,c("HUGO symbols", "Cell population" ,"Polarity")], filtered_features[,"Cell population"])
    scores = CID.append(z, filtered_features, sorted = T)
    
    indexMax = apply(scores, 2, which.max);
    celltypes_train <- rownames(scores)[indexMax];
    #acc = sum(celltypes_train == colnames(z)) / ncol(z)
    
    acc = sapply(unique(celltypes_train), function(x){
      sum(celltypes_train == x & colnames(z) == x) / sum(celltypes_train == x)
    })
    
    acc = mean(acc)
    
    return(list(acc = acc, markers_used = x))
  }
  
  cat(" ..........  Training Signac classifier on training set: \n");
  cat("             nrow = ", nrow(TrainSet), "\n", sep = "");
  cat("             ncol = ", ncol(TrainSet), "\n", sep = "");
  
  q = sapply(100, function(x){
    N_markers = rep(x, 1000)
    numCores = parallel::detectCores()
    parallel::mclapply(N_markers, optim_K, mc.cores =  numCores - 1)
  })
  
  #accs = unlist(q[c(TRUE, FALSE)])
  #mrks = q[c(FALSE, TRUE)]
  #new_markers = mrks[which(accs == max(accs))]
  accs = sapply(q, function(x) x[[1]])
  mrks = lapply(q, function(x) x[[2]])
  new_markers = mrks[which(accs > mean(accs))]
  
  cat(" ..........  Done! Model accuracy on training set = ", max(accs), "\n");
  
  accs_test = sapply(new_markers, function(x){
    filtered_features=subset(x,get("HUGO symbols")%in%rownames(TestSet))
    filtered_features=split(filtered_features[,c("HUGO symbols", "Cell population" ,"Polarity")], filtered_features[,"Cell population"])
    scores = CID.append(TestSet, filtered_features, sorted = T)
    
    indexMax = apply(scores, 2, which.max);
    celltypes_test <- rownames(scores)[indexMax];
    # acc = sum(celltypes_test == colnames(TestSet)) / ncol(TestSet)
    
    acc = sapply(unique(celltypes_train), function(x){
      sum(celltypes_train == x & colnames(z) == x) / sum(celltypes_train == x)
    })
    
    mean(acc)
  })
  
  cat(" ..........  Model accuracy on test set = ", max(accs_test), "\n");
  
  new_markers = mrks[which(accs_test == max(accs_test))]
  
  new_markers = lapply(new_markers, function(x) {
    x[order(x[,2]),]
  })
  
  scores_test = lapply(new_markers, function(x){
    x$`Cell population` = as.character(x$`Cell population`)
    filtered_features=subset(x,get("HUGO symbols")%in%rownames(TestSet))
    filtered_features=split(filtered_features[,c("HUGO symbols", "Cell population" ,"Polarity")], filtered_features[,"Cell population"])
    CID.append(TestSet, filtered_features, sorted = T)
  })
  
  cor_test = sapply(scores_test, function(x){
    d = cor(t(x))
    diag(d) <- 0
    sum(d)
  })
  
  cat(" ..........  Number of models with optimum performance = ", length(which.min(cor_test)), "\n");
  
  new_markers = new_markers[[which.min(cor_test)]]
  
  return(new_markers)
  
}

#' Learn new cell type markers
#'
#' @param E A gene-by-sample count matrix (sparse matrix, matrix, or data.frame) with genes identified by their HUGO symbols (see ?CID.geneconversion), or a list of such matrices, see ?CID.BatchMode.
#' @return an updated reference set with the new cell type
#' @export
CID.Learn4 <- function(E)
{
  cat(" ..........  Learning and refining cell type markers on training subset: \n");
  cat("             cells used = ", ncol(E), "\n", sep = "");
  
  # get unique gene names only
  logik = CID.IsUnique(rownames(E))
  M = E[logik,]
  
  # split into train set into balanced test + dev sets, balanced by the minimum population
  if (length(unique(colnames(M))) != 2)
  {
    df = data.frame(table(colnames(M)))
    min.pop = min(df$Freq)
    dummy = M[,colnames(M) %in% df$Var1[df$Freq == min.pop]]
    dummy2 = M[,colnames(M) %in% df$Var1[df$Freq > min.pop]]
    balmat = list("")
    for (j in 1:length(unique(colnames(dummy2))))
      balmat[[j]] = dummy2[,which(colnames(dummy2) == unique(colnames(dummy2))[j])[sample(1:sum(colnames(dummy2) == unique(colnames(dummy2))[j]), min.pop , replace = FALSE)]]
    balmat = do.call(cbind, balmat)
    M = cbind(dummy, balmat)
  }
  
  
  M = log2(M + 1)
  
  # remove any NA values
  logik = apply(M, 1, function(x) any(!is.na(x)))
  M = M[logik,]
  
  set.seed(42)
  train <- sample(ncol(M), 0.7*ncol(M), replace = FALSE)
  TrainSet <- M[,train]
  TestSet  <- M[,-train]
  
  # run random forest
  lbls = as.factor(colnames(TrainSet))
  TrainSet = Matrix::t(TrainSet)
  xx <- colnames(TrainSet)
  xx = make.names(xx)
  colnames(TrainSet) <- xx
  model <- randomForest::randomForest(lbls ~ ., data = TrainSet)
  
  # test on the TestSet
  TestSet = Matrix::t(TestSet)
  colnames(TestSet) <- xx
  pred <- predict(model, newdata = TestSet)
  table(pred, rownames(TestSet))
  
  # convert to numeric
  q = Matrix::t(Matrix.utils::aggregate.Matrix(Matrix::t(TestSet), factor(colnames(TestSet)), fun = 'mean'))
  gns = rownames(M)
  
  # convert to marker format
  master_types = apply(q, 2, function(x) gns[order(x, decreasing = T)])
  
  # define function for optimizing K
  optim_K <- function(K, x = master_types, y = markers[markers$`Cell population` %in% colnames(master_types),], z = TrainSet)
  {
    cls = colnames(x)
    x = reshape2::melt(x[1:K,])
    x = data.frame(genes = x$value, celltype = x$Var2, polarity = "+")
    colnames(x) <- colnames(y)
    
    x = rbind(x, y)
    x = unique(x)
    x = x[sample(nrow(x), nrow(x), replace = T), ]
    x = unique(x)
    
    filtered_features=subset(x,get("HUGO symbols")%in%rownames(E))
    filtered_features=split(filtered_features[,c("HUGO symbols", "Cell population" ,"Polarity")], filtered_features[,"Cell population"])
    scores = CID.append(z, filtered_features, sorted = T)
    
    indexMax = apply(scores, 2, which.max);
    celltypes_train <- rownames(scores)[indexMax];
    #acc = sum(celltypes_train == colnames(z)) / ncol(z)
    
    acc = sapply(unique(celltypes_train), function(x){
      sum(celltypes_train == x & colnames(z) == x) / sum(celltypes_train == x)
    })
    
    acc = mean(acc)
    
    return(list(acc = acc, markers_used = x))
  }
  
  cat(" ..........  Training Signac classifier on training set: \n");
  cat("             nrow = ", nrow(TrainSet), "\n", sep = "");
  cat("             ncol = ", ncol(TrainSet), "\n", sep = "");
  
  q = sapply(100, function(x){
    N_markers = rep(x, 1000)
    numCores = parallel::detectCores()
    parallel::mclapply(N_markers, optim_K, mc.cores =  numCores - 1)
  })
  
  #accs = unlist(q[c(TRUE, FALSE)])
  #mrks = q[c(FALSE, TRUE)]
  #new_markers = mrks[which(accs == max(accs))]
  accs = sapply(q, function(x) x[[1]])
  mrks = lapply(q, function(x) x[[2]])
  new_markers = mrks[which(accs > mean(accs))]
  
  cat(" ..........  Done! Model accuracy on training set = ", max(accs), "\n");
  
  accs_test = sapply(new_markers, function(x){
    filtered_features=subset(x,get("HUGO symbols")%in%rownames(TestSet))
    filtered_features=split(filtered_features[,c("HUGO symbols", "Cell population" ,"Polarity")], filtered_features[,"Cell population"])
    scores = CID.append(TestSet, filtered_features, sorted = T)
    
    indexMax = apply(scores, 2, which.max);
    celltypes_test <- rownames(scores)[indexMax];
    # acc = sum(celltypes_test == colnames(TestSet)) / ncol(TestSet)
    
    acc = sapply(unique(celltypes_train), function(x){
      sum(celltypes_train == x & colnames(z) == x) / sum(celltypes_train == x)
    })
    
    mean(acc)
  })
  
  cat(" ..........  Model accuracy on test set = ", max(accs_test), "\n");
  
  new_markers = mrks[which(accs_test == max(accs_test))]
  
  new_markers = lapply(new_markers, function(x) {
    x[order(x[,2]),]
  })
  
  scores_test = lapply(new_markers, function(x){
    x$`Cell population` = as.character(x$`Cell population`)
    filtered_features=subset(x,get("HUGO symbols")%in%rownames(TestSet))
    filtered_features=split(filtered_features[,c("HUGO symbols", "Cell population" ,"Polarity")], filtered_features[,"Cell population"])
    CID.append(TestSet, filtered_features, sorted = T)
  })
  
  cor_test = sapply(scores_test, function(x){
    d = cor(t(x))
    diag(d) <- 0
    sum(d)
  })
  
  cat(" ..........  Number of models with optimum performance = ", length(which.min(cor_test)), "\n");
  
  new_markers = new_markers[[which.min(cor_test)]]
  
  return(new_markers)
  
}

#' Learn new cell type markers
#'
#' @param E A gene-by-sample count matrix (sparse matrix, matrix, or data.frame) with genes identified by their HUGO symbols (see ?CID.geneconversion), or a list of such matrices, see ?CID.BatchMode.
#' @return an updated reference set with the new cell type
#' @export
CID.Learn5 <- function(E, markers = NULL)
{
  cat(" ..........  Learning and refining cell type markers on training subset: \n");
  cat("             cells used = ", ncol(E), "\n", sep = "");
  
  # load default markers
  if (is.null(markers))
    load("./data/markers.rda")
  
  # get unique gene names only
  logik = CID.IsUnique(rownames(E))
  M = E[logik,]
  
  # log2 + 1 transform
  M = log2(M + 1)
  M = Matrix::Matrix(Matrix::t(scale(Matrix::t(M))), sparse = T)
  
  # remove any NA values
  logik = apply(M, 1, function(x) any(!is.na(x)))
  M = M[logik,]  
  
  # balance dataset
  min.pop = min(table(colnames(M)))
  min.grp = names(which.min(table(colnames(M))))
  
  dummy = M[,colnames(M) != min.grp]
  dummy = dummy[,sample(1:ncol(dummy), size = 2 * min.pop)]
  M = cbind(dummy, M[, colnames(M) == min.grp])
  
  set.seed(42)
  train <- sample(ncol(M), 0.7*ncol(M), replace = FALSE)
  TrainSet <- M[,train]
  TestSet  <- M[,-train]
  
  # convert to numeric
  q = Matrix::t(Matrix.utils::aggregate.Matrix(Matrix::t(TestSet), factor(colnames(TestSet)), fun = 'mean'))
  gns = rownames(M)
  
  # convert to marker format
  master_types = apply(q, 2, function(x) gns[order(x, decreasing = T)])
  
  # define function for optimizing K
  optim_K <- function(K, x = master_types, y = markers[markers$`Cell population` %in% colnames(master_types),], z = TrainSet)
  {
    cls = colnames(x)
    x = reshape2::melt(x[1:K,])
    x = data.frame(genes = x$value, celltype = x$Var2, polarity = "+")
    colnames(x) <- colnames(y)
    
    x = rbind(x, y)
    x = unique(x)
    
    filtered_features=subset(x,get("HUGO symbols")%in%rownames(E))
    filtered_features=split(filtered_features[,c("HUGO symbols", "Cell population" ,"Polarity")], filtered_features[,"Cell population"])
    
    # randomly sample the filtered_features
    filtered_features = lapply(filtered_features, function(x){unique(x[sample(nrow(x), nrow(x), replace = T), ])})
    
    scores = CID.append(z, filtered_features, sorted = T)
    
    indexMax = apply(scores, 2, which.max);
    celltypes_train <- rownames(scores)[indexMax];
    #acc = sum(celltypes_train == colnames(z)) / ncol(z)
    
    acc = sapply(unique(celltypes_train), function(x){
      sum(celltypes_train == x & colnames(z) == x) / sum(celltypes_train == x)
    })
    
    acc = mean(acc)
    
    p = do.call(rbind, filtered_features)
    rownames(p) <- NULL
    
    return(list(acc = acc, markers_used = p))
  }
  
  cat(" ..........  Training Signac classifier on training set: \n");
  cat("             nrow = ", nrow(TrainSet), "\n", sep = "");
  cat("             ncol = ", ncol(TrainSet), "\n", sep = "");
  
  q = sapply(100, function(x){
    N_markers = rep(x, 1000)
    numCores = parallel::detectCores()
    parallel::mclapply(N_markers, optim_K, mc.cores =  numCores - 1)
  })
  
  #accs = unlist(q[c(TRUE, FALSE)])
  #mrks = q[c(FALSE, TRUE)]
  #new_markers = mrks[which(accs == max(accs))]
  accs = sapply(q, function(x) x[[1]])
  mrks = lapply(q, function(x) x[[2]])
  new_markers = mrks[which(accs > mean(accs))]
  
  cat(" ..........  Done! Model accuracy on training set = ", max(accs), "\n");
  
  accs_test = sapply(new_markers, function(x){
    filtered_features=subset(x,get("HUGO symbols")%in%rownames(TestSet))
    filtered_features=split(filtered_features[,c("HUGO symbols", "Cell population" ,"Polarity")], filtered_features[,"Cell population"])
    scores = CID.append(TestSet, filtered_features, sorted = T)
    
    indexMax = apply(scores, 2, which.max);
    celltypes_test <- rownames(scores)[indexMax];
    # acc = sum(celltypes_test == colnames(TestSet)) / ncol(TestSet)
    
    acc = sapply(unique(celltypes_train), function(x){
      sum(celltypes_train == x & colnames(z) == x) / sum(celltypes_train == x)
    })
    
    mean(acc)
  })
  
  cat(" ..........  Model accuracy on test set = ", max(accs_test), "\n");
  
  new_markers = mrks[which(accs_test == max(accs_test))]
  
  new_markers = lapply(new_markers, function(x) {
    x[order(x[,2]),]
  })
  
  scores_test = lapply(new_markers, function(x){
    x$`Cell population` = as.character(x$`Cell population`)
    filtered_features=subset(x,get("HUGO symbols")%in%rownames(TestSet))
    filtered_features=split(filtered_features[,c("HUGO symbols", "Cell population" ,"Polarity")], filtered_features[,"Cell population"])
    CID.append(TestSet, filtered_features, sorted = T)
  })
  
  cor_test = sapply(scores_test, function(x){
    d = cor(t(x))
    diag(d) <- 0
    sum(d)
  })
  
  cat(" ..........  Number of models with optimum performance = ", length(which.min(cor_test)), "\n");
  
  new_markers = new_markers[[which.min(cor_test)]]
  
  return(new_markers)
  
}

#' Learn new cell type markers
#'
#' @param E A gene-by-sample count matrix (sparse matrix, matrix, or data.frame) with genes identified by their HUGO symbols (see ?CID.geneconversion), or a list of such matrices, see ?CID.BatchMode.
#' @return an updated reference set with the new cell type
#' @export
CID.Learn6 <- function(E, markers = NULL)
{
  cat(" ..........  Learning and refining cell type markers on training subset: \n")
  cat("             cells used = ", ncol(E), "\n", sep = "")
  
  # load default markers
  if (is.null(markers))
    load("./data/markers.rda")
  
  # get unique gene names only
  logik = CID.IsUnique(rownames(E))
  M = E[logik,]
  
  # log2 + 1 transform
  M = log2(M + 1)
  M = Matrix::Matrix(Matrix::t(scale(Matrix::t(M))), sparse = T)
  
  # remove any NA values
  logik = apply(M, 1, function(x) any(!is.na(x)))
  M = M[logik,]
  
  # balance dataset
  min.pop = min(table(colnames(M)))
  min.grp = names(which.min(table(colnames(M))))
  
  dummy = M[,colnames(M) != min.grp]
  dummy = dummy[,sample(1:ncol(dummy), size = 2 * min.pop)]
  M = cbind(dummy, M[, colnames(M) == min.grp])
  
  set.seed(42)
  train <- sample(ncol(M), 0.7*ncol(M), replace = FALSE)
  TrainSet <- M[,train]
  TestSet  <- M[,-train]
  
  # convert to numeric
  q = Matrix::t(Matrix.utils::aggregate.Matrix(Matrix::t(TestSet), factor(colnames(TestSet)), fun = 'mean'))
  gns = rownames(M)
  
  # convert to marker format
  master_types = apply(q, 2, function(x) gns[order(x, decreasing = T)])
  
  # define function for optimizing K
  optim_K <- function(K, x = master_types, y = markers[markers$`Cell population` %in% colnames(master_types),], z = TrainSet)
  {
    cls = colnames(x)
    x = reshape2::melt(x[1:K,])
    x = data.frame(genes = x$value, celltype = x$Var2, polarity = "+")
    colnames(x) <- colnames(y)
    
    x = rbind(x, y)
    x = unique(x)
    
    filtered_features=subset(x,get("HUGO symbols")%in%rownames(E))
    filtered_features=split(filtered_features[,c("HUGO symbols", "Cell population" ,"Polarity")], filtered_features[,"Cell population"])
    
    # randomly sample the filtered_features
    filtered_features = lapply(filtered_features, function(x){unique(x[sample(nrow(x), nrow(x), replace = T), ])})
    
    scores = CID.append(z, filtered_features, sorted = T)
    
    indexMax = apply(scores, 2, which.max);
    celltypes_train <- rownames(scores)[indexMax];
    #acc = sum(celltypes_train == colnames(z)) / ncol(z)
    
    acc = sapply(unique(celltypes_train), function(x){
      sum(celltypes_train == x & colnames(z) == x) / sum(celltypes_train == x)
    })
    
    acc = mean(acc)
    
    p = do.call(rbind, filtered_features)
    rownames(p) <- NULL
    
    return(list(acc = acc, markers_used = p))
  }
  
  cat(" ..........  Training Signac classifier on training set: \n");
  cat("             nrow = ", nrow(TrainSet), "\n", sep = "");
  cat("             ncol = ", ncol(TrainSet), "\n", sep = "");
  
  q = sapply(100, function(x){
    N_markers = rep(x, 1000)
    numCores = parallel::detectCores()
    parallel::mclapply(N_markers, optim_K, mc.cores =  numCores - 1)
  })
  
  #accs = unlist(q[c(TRUE, FALSE)])
  #mrks = q[c(FALSE, TRUE)]
  #new_markers = mrks[which(accs == max(accs))]
  accs = sapply(q, function(x) x[[1]])
  mrks = lapply(q, function(x) x[[2]])
  new_markers = mrks[which(accs > mean(accs))]
  
  cat(" ..........  Done! Model accuracy on training set = ", max(accs), "\n");
  
  accs_test = sapply(new_markers, function(x){
    filtered_features=subset(x,get("HUGO symbols")%in%rownames(TestSet))
    filtered_features=split(filtered_features[,c("HUGO symbols", "Cell population" ,"Polarity")], filtered_features[,"Cell population"])
    scores = CID.append(TestSet, filtered_features, sorted = T)
    
    indexMax = apply(scores, 2, which.max);
    celltypes_test <- rownames(scores)[indexMax];
    # acc = sum(celltypes_test == colnames(TestSet)) / ncol(TestSet)
    
    acc = sapply(unique(celltypes_train), function(x){
      sum(celltypes_train == x & colnames(z) == x) / sum(celltypes_train == x)
    })
    
    mean(acc)
  })
  
  cat(" ..........  Model accuracy on test set = ", max(accs_test), "\n");
  
  new_markers = mrks[which(accs_test == max(accs_test))]
  
  new_markers = lapply(new_markers, function(x) {
    x[order(x[,2]),]
  })
  
  scores_test = lapply(new_markers, function(x){
    x$`Cell population` = as.character(x$`Cell population`)
    filtered_features=subset(x,get("HUGO symbols")%in%rownames(TestSet))
    filtered_features=split(filtered_features[,c("HUGO symbols", "Cell population" ,"Polarity")], filtered_features[,"Cell population"])
    CID.append(TestSet, filtered_features, sorted = T)
  })
  
  cor_test = sapply(scores_test, function(x){
    d = cor(t(x))
    diag(d) <- 0
    sum(d)
  })
  
  cat(" ..........  Number of models with optimum performance = ", length(which.min(cor_test)), "\n");
  
  new_markers = new_markers[[which.min(cor_test)]]
  
  return(new_markers)
  
}

#' Get indices of training markers
#'
#' @param E A gene-by-sample count matrix (sparse matrix, matrix, or data.frame) with genes identified by their HUGO symbols (see ?CID.geneconversion), or a list of such matrices, see ?CID.BatchMode.
#' @param data.dir if default, uses the standard Signac markers (see ?CID.SeeMarkers).
#' @return a list of cell indices and names to use for training set
#' @export
CID.GetTrainingSet <- function(E, data.dir = NULL, method.use = "max.genes.detected")
{
  inputs = match.call()
  
  if (method.use == "min.entropy")
  {
    # get edges for cell-cell similarity network
    edges = CID.GetEdges(inputs)
    
    # compute distance matrix
    distM = CID.GetDistMat(edges = edges)
    
    # Calculate normalized shannon entropy for each cell j for all connections with shortest path < N
    shannon = rep(0, ncol(E))
    ac = colnames(E)
    
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
    m = Matrix::Matrix(0, nrow = length(ac), ncol = length(unique(ac)), sparse = T)
    m[cbind(1:nrow(m), as.numeric(factor(ac)))] <- 1
    
    res = dM %*% m
    res = res / Matrix::rowSums(res)
    N_unique = length(unique(ac))
    shannon = apply(res, 1, function(freqs) {-sum(freqs[freqs != 0] * log(freqs[freqs != 0])) / log(2) / log2(N_unique)})
    
    df = data.frame(cells = ac, entropy = shannon)
    
    q = lapply(unique(df$cells), function(x){
      which(df$cells == x)[order(df$entropy[df$cells == x])[1:round(sum(df$cells == x) * 0.1)]]
    })
    names(q) <- unique(df$cells)
    
    idx = which(names(q) == "Other")
    q[-idx]
  }
  
  if (method.use == "max.genes.detected") {
    ac = colnames(E)
    
    total_counts = Matrix::colSums(E != 0)
    
    df = data.frame(cells = ac, entropy = total_counts)
    
    q = lapply(unique(df$cells), function(x){
      which(df$cells == x)[order(df$entropy[df$cells == x], decreasing = TRUE)[1:round(0.1 * sum(df$cells == x))]]
    })
    names(q) <- unique(df$cells)
    
  }
  return(q)
  
}

#' Get indices of training markers
#'
#' @param E A gene-by-sample count matrix (sparse matrix, matrix, or data.frame) with genes identified by their HUGO symbols (see ?CID.geneconversion), or a list of such matrices, see ?CID.BatchMode.
#' @param data.dir if default, uses the standard Signac markers (see ?CID.SeeMarkers).
#' @return a list of cell indices and names to use for training set
#' @export
CID.GetTrainingSet2 <- function(E, cr,  data.dir = NULL)
{
  inputs = match.call()
  
  # get edges for cell-cell similarity network
  edges = CID.GetEdges(inputs)
  
  # compute distance matrix
  distM = CID.GetDistMat(edges = edges)
  
  # Calculate normalized shannon entropy for each cell j for all connections with shortest path < N
  shannon = rep(0, ncol(E))
  ac = colnames(E)
  
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
  m = Matrix::Matrix(0, nrow = length(ac), ncol = length(unique(ac)), sparse = T)
  m[cbind(1:nrow(m), as.numeric(factor(ac)))] <- 1
  
  res = dM %*% m
  res = res / Matrix::rowSums(res)
  N_unique = length(unique(ac))
  shannon = apply(res, 1, function(freqs) {-sum(freqs[freqs != 0] * log(freqs[freqs != 0])) / log(2) / log2(N_unique)})
  
  df = data.frame(cells = ac, entropy = shannon)
  
  q = lapply(unique(df$cells), function(x){
    which(df$cells == x)[order(df$entropy[df$cells == x])]
  })
  names(q) <- unique(df$cells)
  idx = which(names(q) == "Other")
  q[-idx]
  
  total_counts = Matrix::colSums(E != 0)
  
  df = data.frame(cells = ac, entropy = total_counts)
  
  q = lapply(unique(df$cells), function(x){
    which(df$cells == x)[order(df$entropy[df$cells == x], decreasing = TRUE)[1:round(0.1 * sum(df$cells == x))]]
  })
  names(q) <- unique(df$cells)
  
  
  return(q)
  
}

#' Unsupervised main function
#'
#' @param E A gene-by-sample count matrix (sparse matrix, matrix, or data.frame) with genes identified by their HUGO symbols (see ?CID.geneconversion), or a list of such matrices, see ?CID.BatchMode.
#' @param data.dir directory of SPRING files "edges.csv" and "categorical_coloring_data.json"
#' @return Filtered markers where each marker must have at least ncells that express at least ncounts
#' @export
CID.UnsupervisedLearning <- function(E, data.dir = NULL)
{
  # check input
  inputs = match.call()
  if (class(E) %in% c("matrix", "data.frame"))
    E = Matrix::Matrix(as.matrix(E), sparse = T)
  if (!is.null(data.dir))
    data.dir = gsub("\\/$", "", data.dir, perl = TRUE);
  cat(" ..........  Entry in CID.CellID \n");
  ta = proc.time()[3];
  
  # main function
  cat(" ..........  Running unsupervised classification for cell types on input data matrix :\n");
  cat("             nrow = ", nrow(E), "\n", sep = "");
  cat("             ncol = ", ncol(E), "\n", sep = "");
  
  # get edges for cell-cell similarity network
  edges = CID.GetEdges(inputs)
  
  # compute distance matrix
  distMat = CID.GetDistMat(edges = edges)
  
  # get Louvain clusters
  louvain = CID.Louvain(edges = edges)
  
  # get overlap between all Louvain clusters in KNN network
  m = Matrix::Matrix(0, nrow = length(louvain), ncol = length(unique(louvain)), sparse = T)
  m[cbind(1:nrow(m), as.numeric(factor(louvain)))] <- 1
  res = distMat[[1]] %*% m
  colnames(res) <- unique(louvain)
  q = aggregate(as.matrix(res), list(louvain), sum)
  colnames(q) <- c("Groups", levels(factor(louvain)))
  q2 = as.matrix(q[,-1])
  rownames(q2) <- colnames(q2)
  q2 = log2(q2 + 1)
  q2 = q2 / rowSums(q2)
  
  # get any Louvain clusters with >f %  edge connections
  idx = apply(q2, 1, function(x) colnames(q2)[which (x > 0.02)])
  names(idx) = rownames(q2)
  
  # merge them into groups
  celltypes = louvain
  
  for (j in 1:length(idx))
  {
    logik = which(sapply(idx, function(x) sum(idx[[j]] %in% x) > 0))
    celltypes[celltypes %in% unlist(idx[logik])] = paste("Group", j)
  }
  
  # get markers for these cell types
  cat("             Signac found", length(unique(celltypes)), "celltypes and ", length(unique(louvain)), "cellstates!\n");
  colnames(E) <- celltypes
  new_lbls = CID.PosMarkers2(E, celltypes)
  celltypes = new_lbls$lbls
  cat("             Establishing markers for cellstates...\n");
  
  # break these markers down into cell states
  cellstates = celltypes
  pb = txtProgressBar(min = 0, max = length(unique(celltypes)), initial = 0) 
  for (j in 1:length(unique(celltypes))){
    logik = celltypes == unique(celltypes)[j]
    for (k in 1:length(unique(louvain[logik])))
    {
      Q = E[,logik]
      dummy = louvain[logik]
      lbls = rep("All", ncol(Q))
      lbls[dummy == unique(louvain[logik])[k]] = dummy[dummy == unique(louvain[logik])[k]]
      colnames(Q) <- lbls
      new_lbls = CID.PosMarkers2(Q, lbls)
      cellstates[louvain == unique(louvain[logik])[k]] = new_lbls$lbls[new_lbls$lbls != "All"]
    }
    setTxtProgressBar(pb,j)
  }
  
  cr = list(celltypes_unsup  = celltypes,
            cellstates_unsup = cellstates)
  
  tb = proc.time()[3] - ta;
  cat("\n ..........  Exit CID.CellID.\n");
  cat("             Execution time = ", tb, " s.\n", sep = "");
  return (cr);
}



#' filters the geneset markers
#'
#' @param D An expression matrix with features (genes) in rows and samples (cells) in columns.
#' @param markersG A data frame with four columns ('HUGO Symbol', 'Cell Population', 'ENTREZ ID', 'Polarity'). Default is internally set to data(markers_v4).
#' @param pval p-value cutoff, as described in the manuscript. Default is 0.01.
#' @param sorted see ?CID.CellID
#' @return A list of markers and features.
#' @export
CID.filter2 <- function(l, m, n, o)
{
  if (o)
    m = m[m$Polarity == "+", ]
  l = l[intersect(row.names(l),m$"HUGO symbols"),,drop=F]
  l = l[apply(l, 1, sd) != 0,]
  markers.names = unique(m$`Cell population`)
  cat(" ..........  Filtering markers for features: \n","           ",paste(markers.names,collapse = ", ", sep = ""))
  features=subset(m,get("HUGO symbols")%in%rownames(l))
  features=split(features[,c("HUGO symbols", "Cell population" ,"Polarity")], features[,"Cell population"])
  features=features[intersect(markers.names,names(features))]
  features=features[sapply(features,function(x)sum(x$Polarity == "+") > 1)]
  missing.populations=setdiff(markers.names,names(features))
  features23=features[sapply(features,function(x)nrow(x)<3)]
  features=features[sapply(features,function(x)nrow(x)>=3)]
  if (NROW(features) >= 3)
  {
    L = lapply(features,function(x){na.omit(l[intersect(row.names(l),x$"HUGO symbols"[x$Polarity == "+"]),,drop=F])})
    # set.seed("42")
    
    # get pairwise combinations of markers with the top five markers
    # downsampling does not seem to affect results; done for speed.
    cols = lapply(L, function(x){
      t( combn(1:nrow(x), 2))
    })
    
    # run pairwise correlation test
    Q = mapply(function(x, y){
      apply( x , 1 , function(z) cor.test( y[z[1], ] , y[  z[2], ] )$p.value )
    }, x = cols, y = L)
    
    # remove features with weak overall correlation; median had the best sensitivity
    QC = sapply(Q, median)
    logik = QC > n & sapply(L, nrow) > 2;
    features = features[!logik]; L = L[!logik]
    if(!length(features))
      return(NULL)
    
    # combine features 
    features = c(features, features23)
    
  } else {
    for (j in levels(m$Polarity))
    {
      L = lapply(features,function(x){na.omit(l[intersect(row.names(l),x$"HUGO symbols"[x$Polarity == j]),,drop=F])})
      
      # remove markers that do not cluster with maximal group
      M = do.call(rbind, L)
      
      # generate correlation matrix, cluster and cut into two groups
      cor_mat = qlcMatrix::cosSparse(Matrix::t(M))
      ix = hclust(dist(cor_mat))
      grps = cutree(ix, k = (length(L)))
      vec = as.character(unlist(sapply(features, function(x) x$`Cell population`[x$Polarity == j])))
      q = table(data.frame(from = vec, to = grps))
      
      # find the maximal group
      idx = apply(q,1,which.max)
      genes_keep = names(grps)[grps %in% idx]
      
      # remove unstable markers from features
      features = lapply(features, function(x){
        x1 = x[x$`HUGO symbols` %in% genes_keep & x$Polarity == j,]
        rbind(x[x$Polarity != j,], x1)
      })
    }
    features = c(features, features23)
  }
  filtered.populations=setdiff(markers.names,setdiff(names(features),missing.populations))
  if(length(filtered.populations) > 0)
    cat(paste("\nNo markers exist for this feature(s) due to filtering:",paste(filtered.populations,collapse=", ")), "\n")
  if(length(filtered.populations)  == 0 & length(missing.populations) == 0)
    cat("\n ..........  Markers exist for all features! \n")
  features
}

#' filters the geneset markers
#'
#' @param expression An expression matrix with features (genes) in rows and samples (cells) in columns.
#' @param markersG A data frame with four columns ('HUGO Symbol', 'Cell Population', 'ENTREZ ID', 'Polarity'). Default is internally set to data(markers_v4).
#' @param pval p-value cutoff, as described in the manuscript. Default is 0.1.
#' @return A list of markers and features.
#' @export
CID.filter <- function(expression, markersG, pval = pval)
{
  expression = expression[intersect(row.names(expression),markersG$"HUGO symbols"),,drop=F]
  expression = log(expression + 1, 2)
  markers.names = unique(markersG$`Cell population`)
  cat(" ..........  Filtering markers for features: \n","           ",paste(markers.names,collapse = ", ", sep = ""))
  features=subset(markersG,get("HUGO symbols")%in%rownames(expression))
  features=split(features[,c("HUGO symbols", "Polarity")], features[,"Cell population"])
  features=features[intersect(markers.names,names(features))]
  features=features[sapply(features,function(x)nrow(x)>0)]
  missing.populations=setdiff(markers.names,names(features))
  if (NROW(features) > 2)
  {
    L = lapply(features,function(x){na.omit(expression[intersect(row.names(expression),x$"HUGO symbols"[x$Polarity == "+"]),,drop=F])})
    Q = lapply(L, function(x) {
      q = qlcMatrix::cosSparse(Matrix::t(x))
      diag(q) <- 0
      apply(q,1,function(x) max(abs(x)))})
    filter = summary(unlist(Q))[2]
    QC = lapply(Q, function(x) sum(x < filter))
    QC = mapply(function(x,y) { 1 - phyper(x - 1, sum(unlist(Q) < filter), length(unlist(Q)) - sum(unlist(Q) < filter), length(y)) }, x = QC, y = Q)
    logik = QC < pval;
    features = features[!logik]
    #filter = summary(unlist(Q))[1]
    filter = -1;
    nms = unlist(lapply(Q, function(x) names(x)[x > filter]))
    features = lapply(features, function(x) rbind(x[x[,1] %in% nms,], x[x$Polarity == "-",]))
  }
  features = features[sapply(features, function(x) sum(x$Polarity == "+") != 0)]
  filtered.populations=setdiff(markers.names,setdiff(names(features),missing.populations))
  if(length(missing.populations)>0)
    warning(paste("No markers exist for this feature(s):",paste(missing.populations,collapse=", ")), ".\n")
  if(length(filtered.populations) > 0)
    cat(paste("\nNo markers exist for this feature(s) due to filtering:",paste(filtered.populations,collapse=", ")), "\n")
  if(length(filtered.populations)  == 0 & length(missing.populations) == 0)
    cat("\n ..........  Markers exist for all features! \n")
  features
}

#' Calculates Signac scores
#'
#' @param expression An expression matrix with features (genes) in rows and samples (cells) in columns.
#' @param featuresG A list of features, where each list entry is a character vector of markers.
#' @param sorted see ?CID.CellID
#' @return A matrix of CellID scores; features by cells.
#' @export
CID.append=function(expression,featuresG, sorted)
{
  if (!sorted)
  {
    # subset expression matrix
    Z = expression[row.names(expression) %in% as.character(unlist(lapply(featuresG, function(x) x = x[,1]))), ]
    # Z-score transform
    Z = Matrix::t(scale(Matrix::t(Z)))
    # Get CellID scores
    res = as.data.frame(Matrix::t(do.call(cbind,
                                          lapply(featuresG,function(x){
                                            if (any(x$Polarity == "-")) {
                                              apply(Z[intersect(row.names(Z),x$`HUGO symbols`[x$Polarity == "+"]),,drop=F],2,function(x) mean(x, na.rm = T)) -
                                                apply(Z[intersect(row.names(Z),x$`HUGO symbols`[x$Polarity == "-"]),,drop=F],2,function(x) mean(x, na.rm = T))
                                            } else if (!any(x$Polarity == "-")) {
                                              apply(Z[intersect(row.names(Z),x$`HUGO symbols`[x$Polarity == "+"]),,drop=F],2,function(x) mean(x, na.rm = T))
                                            } else if (!any(x$Polarity == "+")) {
                                              apply(Z[intersect(row.names(Z),x$`HUGO symbols`[x$Polarity == "-"]),,drop=F],2,function(x) mean(-1 * x, na.rm = T))
                                            }
                                          }))))
    res = res[order(rownames(res)),]
    res
  } else {
    # subset expression matrix
    Z = expression[row.names(expression) %in% as.character(unlist(lapply(featuresG, function(x) x = x[,1]))), ]
    # Get CellID scores
    res = as.data.frame(Matrix::t(do.call(cbind,
                                          lapply(featuresG,function(x){
                                            apply(Z[intersect(row.names(Z),x$`HUGO symbols`[x$Polarity == "+"]),,drop=F],2,function(x) mean(x, na.rm = T))
                                          }))))
    res = res[order(rownames(res)),]
    res
  }
}

#' Smoothing function
#'
#' @param ac List containing a character vector where each element is a cell type or cell state assignment
#' @param dM Distance matrix (see ?CID.GetDistMat)
#' @return Smoothed cell type or cell state assignments
#' @export
CID.smooth <- function(ac,dM)
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
  m = Matrix::Matrix(0, nrow = length(Y), ncol = length(unique(Y)), sparse = T)
  m[cbind(1:nrow(m), as.numeric(factor(Y)))] <- 1
  res = dM %*% m
  res = res / Matrix::rowSums(res)
  mx.idx = apply(res, 1, function(x) which.max(x))
  mx = apply(res, 1, max)
  Y[mx > 0.5] = levels(factor(Y))[mx.idx[mx > 0.5]]
  ac[logik] = Y
  return(ac)
}

#' Entropy
#' 
#' @param ac A character vector of cell type labels
#' @param distM The distance matrix, see ?CID.GetDistMat
#' @export
CID.entropy <- function(ac,distM)
{
  # "Y" labels will be ac but with some "Others"
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
  m = Matrix::Matrix(0, nrow = length(ac), ncol = length(unique(ac)), sparse = T)
  m[cbind(1:nrow(m), as.numeric(factor(ac)))] <- 1
  
  res = dM %*% m
  res = res / Matrix::rowSums(res)
  N_unique = length(unique(Y))
  shannon = apply(res, 1, function(freqs) {-sum(freqs[freqs != 0] * log(freqs[freqs != 0])) / log(2) / log2(N_unique)})
  
  #df = data.frame(cells = Y, entropy = shannon)
  #ggplot(df, aes(x=cells, y=entropy, fill = cells)) + geom_boxplot()
  q = c(shannon, -shannon)
  # note symmetry
  # hist(q)
  # use Gaussian
  
  Y[shannon > 2 * sd(q)] = "Other"
  
  #df = data.frame(cells = Y, entropy = shannon)
  #ggplot(df, aes(x=cells, y=entropy, fill = cells)) + geom_boxplot()
  
  return(Y)
}

#' Write JSON file for SPRING visualization
#'
#' @param cr Output from CID.CellID. See ?CID.CellID
#' @param json_new Filename for new SPRING visualization. Default is json_new = "categorical_coloring_data_new.json".
#' @param data.dir Directory where file 'categorical_coloring_data.json' is located. If supplied, it will append this file to contain tracks for cell type / state classifications.
#' @return Smoothed cell type or cell state assignments
#' @export
CID.writeJSON <- function(cr, json_new = "categorical_coloring_data.json", data.dir = NULL)
{
  if (!is.null(data.dir))
  {
    data.dir = gsub("\\/$", "", data.dir, perl = TRUE);
    if (file.exists(paste(data.dir, 'categorical_coloring_data.json', sep = "/")))
    {
      json_file = 'categorical_coloring_data.json'
      gJ <- paste(data.dir,json_file,sep = "/")
      json_data <- rjson::fromJSON(file=gJ)
      json_out = jsonlite::toJSON(json_data, auto_unbox = TRUE)
      write(json_out,paste(data.dir,'categorical_coloring_data_legacy.json',sep="/"))
    } else {
      json_file = 'categorical_coloring_data_old.json'
      gJ <- paste(data.dir,json_file,sep = "/")
      json_data <- rjson::fromJSON(file=gJ)
      json_out = jsonlite::toJSON(json_data, auto_unbox = TRUE)
      write(json_out,paste(data.dir,'categorical_coloring_data_legacy.json',sep="/"))
    }
  }
  if ("louvain" %in% names(cr))
  {
    Q = as.character(cr$louvain)
    json_data$ClustersWT$label_list = Q
    Ntypes = length(unique(Q))
    qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) #len = 74
    #pie(rep(1,num_col), col=(col_vector[1:num_col]))
    col_palette <- as.list(col_vector[1:Ntypes]); # or sample if you wish
    names(col_palette) <- unique(Q)
    json_data$ClustersWT$label_colors = col_palette
  }
  if ("celltypes" %in% names(cr))
  {
    Q = cr$celltypes
    json_data$CellTypesID$label_list = Q
    C = get_colors(Q)
    json_data$CellTypesID$label_colors = as.list(C[[1]])
  }
  if ("cellstates" %in% names(cr))
  {
    Q = cr$cellstates
    json_data$CellStatesID$label_list = Q
    C = get_colors(Q)
    json_data$CellStatesID$label_colors = as.list(C[[1]])
  }
  if ("celltypes_novel" %in% names(cr))
  {
    Q = cr$celltypes_novel
    json_data$CellTypesID_Novel$label_list = Q
    C = get_colors(Q)
    Ntypes = sum(C[[1]] == "")
    qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) #len = 74
    #pie(rep(1,num_col), col=(col_vector[1:num_col]))
    C[[1]][C[[1]] == ""] <- col_vector[1:Ntypes]; # or sample if you wish
    json_data$CellTypesID_Novel$label_colors = as.list(C[[1]])
  }
  if ("cellstates_novel" %in% names(cr))
  {
    Q = cr$cellstates_novel
    json_data$CellStatesID_Novel$label_list = Q
    C = get_colors(Q)
    Ntypes = sum(C[[1]] == "")
    qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) #len = 74
    #pie(rep(1,num_col), col=(col_vector[1:num_col]))
    C[[1]][C[[1]] == ""] <- col_vector[1:Ntypes]; # or sample if you wish
    json_data$CellStatesID_Novel$label_colors = as.list(C[[1]])
  }
  if ("cellstates_merged" %in% names(cr))
  {
    Q = cr$cellstates_merged
    json_data$CellStatesID_merged$label_list = Q
    C = get_colors(Q)
    json_data$CellStatesID_merged$label_colors = as.list(C[[1]])
  }
  if ("cellstates_merged_novel" %in% names(cr))
  {
    Q = cr$cellstates_merged_novel
    json_data$CellStatesID_novel_merged$label_list = Q
    C = get_colors(Q)
    Ntypes = sum(C[[1]] == "")
    qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) #len = 74
    #pie(rep(1,num_col), col=(col_vector[1:num_col]))
    C[[1]][C[[1]] == ""] <- col_vector[1:Ntypes]; # or sample if you wish
    json_data$CellStatesID_novel_merged$label_colors = as.list(C[[1]])
  }
  if ("Immune" %in% names(cr))
  {
    Q = cr$Immune
    json_data$Immune$label_list = Q
    C = get_colors(Q)
    Ntypes = sum(C[[1]] == "")
    qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) #len = 74
    #pie(rep(1,num_col), col=(col_vector[1:num_col]))
    C[[1]][C[[1]] == ""] <- col_vector[1:Ntypes]; # or sample if you wish
    json_data$Immune$label_colors = as.list(C[[1]])
  }
  if ("immunestates" %in% names(cr))
  {
    Q = cr$immunestates
    json_data$ImmuneStates$label_list = Q
    C = get_colors(Q)
    Ntypes = sum(C[[1]] == "")
    qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) #len = 74
    #pie(rep(1,num_col), col=(col_vector[1:num_col]))
    C[[1]][C[[1]] == ""] <- col_vector[1:Ntypes]; # or sample if you wish
    json_data$ImmuneStates$label_colors = as.list(C[[1]])
  }
  json_data$CellTypesID$label_colors = replace(json_data$CellTypesID$label_colors, json_data$CellTypesID$label_colors == "", json_data$Immune_L0$label_colors[match(names(json_data$CellTypesID$label_colors), names(json_data$Immune_L0$label_colors))][json_data$CellTypesID$label_colors == ""])
  json_data$CellStatesID$label_colors = replace(json_data$CellStatesID$label_colors, json_data$CellStatesID$label_colors == "", json_data$Immune_L0$label_colors[match(names(json_data$CellStatesID$label_colors), names(json_data$Immune_L0$label_colors))][json_data$CellStatesID$label_colors == ""])
  json_data = json_data[order(names(json_data))]
  json_data_backup = json_data;
  fn = json_new
  if (is.null(data.dir))
  {
    fn = json_new;
    data.dir = getwd()
    data.dir = gsub("\\/$", "", data.dir, perl = TRUE);
  }
  json_out = jsonlite::toJSON(json_data, auto_unbox = TRUE)
  write(json_out,paste(data.dir,fn,sep="/"))
  cat(paste(data.dir,fn,sep="/"), "has been written to directory! \n")
  new.dirs = list.dirs(dirname(data.dir))
  if (length(new.dirs) != 1)
  {
    for (j in 1:length(new.dirs))
    {
      txt = list.files(new.dirs[j])
      flag = "cell_filter.txt" %in% list.files(new.dirs);
      new.dirs = new.dirs[flag]
      if ("cell_filter.txt" %in% txt)
      {
        idx = read.table(paste(new.dirs[j], "cell_filter.txt", sep = "/"), quote="\"", comment.char="", stringsAsFactors=FALSE)$V1 + 1;
        for (k in 1:length(json_data))
        {
          json_data[[k]]$label_list = json_data[[k]]$label_list[idx]
          json_data[[k]]$label_colors=json_data[[k]]$label_colors[names(json_data[[k]]$label_colors) %in% unique(json_data[[k]]$label_list)]
        }
        json_out = jsonlite::toJSON(json_data, auto_unbox = TRUE)
        write(json_out,paste(new.dirs[j],fn,sep="/"))
        #cat(paste(new.dirs[j],fn,sep="/"), "has been written to directory! \n")
        json_data = json_data_backup;
      }        
    }
    json_out = jsonlite::toJSON(json_data, auto_unbox = TRUE)
    write(json_out,paste(new.dirs[j],fn,sep="/"))
    #cat(paste(new.dirs[j],fn,sep="/"), "has been written to directory! \n")
    json_data = json_data_backup;
  }
  return(json_data)
}

#' Get HEX colors
#' 
#' @param P A character vector of cell type labels
get_colors <- function(P)
{
  P <- sort(unique(as.character(P)));
  colorcount = length(unique(P))
  cs = character(colorcount)
  # main cell type categories will be consistently labeled:
  main_types = c("B.cells"     ,     "Epithelial", "Fibroblasts"        ,
                 "Granulocytes",     "MPh"       , "Plasma.cells"       ,
                 "Secretory"   ,     "TNK"       , "Epithelial.Urinary" ,
                 "Other"       ,     "HSC"       , "pDCs"               ,
                 "Erythro"     ,     "Platelets" , "NonImmune"          ,
                 "Immune")
  main_colors = c("#aa3596", "#387f50", "#a3ba37",
                  "#D3D3D3", "#f0ff51", "#d878de",
                  "#8c933e", "#90c5f4", "#90b771",
                  "#C0C0C0", "#dbd0d0", "#e1f7d5",
                  "#ffbdbd", "#c9c9ff", "#7fc97f",
                  "#bb7fc7")
  
  # sub cell types will be consistently labeled:
  sub_types = c( "Dendritic.cells.activated"  ,  "Dendritic.cells.resting"     , "Eosinophils"       ,
                 "Macrophages.M0"             ,  "Macrophages.M1"              , "Macrophages.M2"    ,
                 "Mast.Progenitors"           ,  "Monocytes"                   , "Neutrophils"       ,
                 "NK"                         ,  "T.cells.CD4.memory.activated", "T.cells.CD4.naive" ,
                 "T.cells.CD8"                ,  "T.cells.FH"                  , "T.regs"            ,
                 "Mast.cells.activated"       ,  "T.cells.CD4.memory.resting"  , "B.cells.memory"    ,
                 "B.cells.naive"              ,  "T.cells"                     , "T.cells.CD4"       ,
                 "T.cells.CD8"                ,  "T.cells.CD4.FH.regs"         , "Not.Mast"          , 
                 "Mast"                       ,  "Mast.cells.activated"        , "Mast.cells.resting",
                 "Macrophages"                ,  "Dendritic"                   , "T.cells.follicular.helper",
                 "NK.cells"                   ,  "Mon.NonClass"                , "Mon.Classical"    ,
                 "DCs.Resting"                ,  "DCs.Activated"               , "T.cells.CD4n"     ,
                 "T.cells.CD4m"               ,  "Mon.NonClass"                , "Neutro.inflam"    ,
                 "T.cells"                    ,  "T.CD4.FH"                    , "T.CD4.FH.regs"   ,
                 "Neutrophil")
  sub_colors =c( "#a3b300"                    ,  "#f9d801"                     , "#d0e4e5"          ,
                 "#F9A602"                    ,  "#f97501"                     , "#d6b171"          ,
                 "#d4e881"                    ,  "#f0ff51"                     , "#f7f2b2"          ,
                 "#ad9bf2"                    ,  "#4ebdbd"                     , "#8c9ee1"          ,
                 "#4e59bd"                    ,  "#d4e881"                     , "#2038b0"          ,
                 "#e2e8c9"                    ,  "#def9f9"                     , "#aa3596"          ,
                 "#edc5e6"                    ,  "#4e59bd"                     , "#8c9ee1"          ,
                 "#90c5f4"                    ,  "#90c5f4"                     , "#D3D3D3"          ,
                 "#d4e881"                    ,  "#d4e881"                     , "#f7f2b2"          ,
                 "#F9A602"                    ,  "#f93a01"                     , "#d4e881"          ,
                 "#ad9bf2"                    , "#7fc97f"                      , "#bb7fc7"          ,
                 "#d4e881"                    , "#ad9bf2"                      , "#db1bbe"          ,
                 "#bb7fc7"                    , "#e68181"                      , "#ebbf8a"          ,
                 "#f7f2b2"                    , "#c75e4c"                      , "#92c74c"          , 
                 "#f7f2b2")
  
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

#' Community substructure by Louvain community detection
#'
#' @param edges A data frame or matrix with edges
#' @return The community substructures of the graph
#' @export
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

#' Get distance matrix
#'
#' @param edges a data frame with two columns; V1 and V2, specifying the cell-cell edge list for the network
#' @param n maximum network distance to subtend (n neighbors)
#' @return adjacency matrices for distances < n
#' @export
CID.GetDistMat <- function(edges, n = 4)
{
  "%^%" <- function(A, n) {if(n == 1) A else A %*% (A %^% (n-1)) }
  # create adjacency matrix
  m <- methods::new("ngTMatrix", 
                    i = c(as.integer(edges$V1)-1L, as.integer(edges$V2)-1L), 
                    j = c(as.integer(edges$V2)-1L, as.integer(edges$V1)-1L),
                    Dim = as.integer(c(max(edges), max(edges))))
  dm = list("") # initialize distance matrix
  for (j in 1:n)
    if(j == 1) dm[[j]] = m else dm[[j]] = m %^% j
  return(dm)
}

#' Split large data into smaller chunks
#'
#' @param E Count matrix
#' @param number_of_chunks Number of chunks to divide the data
#' @return list where each element is a chunk of the original matrix
#' @export
array_split <- function(E, number_of_chunks) {
  if (number_of_chunks > ncol(E))
    number_of_chunks = ncol(E)
  colIdx <- seq_len(ncol(E))
  lapply(split(colIdx, cut(colIdx, pretty(colIdx, number_of_chunks))), function(x) E[, x ])
}

#' Extract unique elements
#'
#' @param x A character vector 
#' @return boolean, unique elements are TRUE
#' @export
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

#' Get positive genes by z-score
#'
#' @param D Expression matrix with genes for rows, samples for columns
#' @param acn Vector of cell type labels 
#' @return List where each element contains the z score for the one vs. all comparison
#' @export
CID.PosMarkers2 <- function(D, acn)
{
  cts = unique(colnames(D));
  cts = cts[cts != "All"]
  mrks = list("")
  new_lbls = list("")
  
  # remove MT- and RP genes
  logik = grepl("^MT-", rownames(D)) | grepl("^RP", rownames(D))
  D = D[!logik,]
  
  # CV filter; keep top 2000 genes
  #kmu = Matrix::rowMeans(D)
  #ksd = apply(D,1,sd)
  #cv = ksd / (kmu + 10)
  #idx = sort(cv, decreasing = T)
  #D = D[rownames(D) %in% names(idx)[1:2000],]
  
  
  # establish all means
  if ("All" %in% colnames(D))
  {
    all_means = Matrix::rowMeans(D[,colnames(D) == "All"])
    
    mrks = lapply(cts, function(x){
      logik = colnames(D) == x
      if (sum(logik) > 200)
        logik = 1:ncol(D) %in% sample(which(colnames(D) == x), 200)
      kmu = Matrix::rowMeans(D[,logik])
      s  = apply(D[, logik], 1, sd)
      ff = (kmu - all_means) / (s + 0.02)
      ff = sort(ff, decreasing = T)
      return(ff)
    })
    names(mrks) <- cts
    
    new_lbls = lapply(mrks, function(x){
      gns = names(x)[1:2]
      paste("+", gns, collapse = " ", sep = "")
    })
    
    names(new_lbls) <- cts
    names(mrks) <- cts
    
    Y = acn
    
    for (j in 1:length(new_lbls))
      Y[colnames(D) == names(new_lbls)[j]] = new_lbls[[j]]
  } else {
    mrks = lapply(cts, function(x){
      logik = colnames(D) == x
      if (sum(logik) > 200)
        logik = 1:ncol(D) %in% sample(which(colnames(D) == x), 200)
      kmu = Matrix::rowMeans(D[,logik])
      all_means = Matrix::rowMeans(D[,!logik])
      s  = apply(D[, logik], 1, sd)
      ff = (kmu - all_means) / (s + 0.02)
      ff = sort(ff, decreasing = T)
      return(ff)
    })
    names(mrks) <- cts
    
    new_lbls = lapply(mrks, function(x){
      gns = names(x)[1:2]
      paste("+", gns, collapse = " ", sep = "")
    })
    
    names(new_lbls) <- cts
    names(mrks) <- cts
    
    Y = acn
    
    for (j in 1:length(new_lbls))
      Y[colnames(D) == names(new_lbls)[j]] = new_lbls[[j]]
  }
  
  
  
  return(list(lbls = Y, mrks = mrks))
}

#' Run DEG analysis with Seurat wrapper
#'
#' @param D Expression matrix with genes for rows, samples for columns
#' @param acn Vector of cell type labels 
#' @return List where each element contains the DEG tables for the one vs. all comparison
#' @export
CID.PosMarkers3 <- function(D, threshold)
{
  # get lbls
  lbls = colnames(D)
  # lbls = paste("Cluster", lbls)
  # set colnames
  colnames(D) <- seq(1, ncol(D))
  # make sure row names are not redundant
  logik = CID.IsUnique(rownames(D))
  D = D[logik,]
  # Set up object
  ctrl <- suppressWarnings(Seurat::CreateSeuratObject(counts = D))
  ctrl <- Seurat::NormalizeData(object = ctrl, verbose = F)
  ctrl <- Seurat::AddMetaData(ctrl, metadata=lbls, col.name = "celltypes")
  ctrl <- Seurat::SetIdent(ctrl, value='celltypes')
  cts = unique(lbls);
  mrks = list("")
  dd = list("")
  for (j in 1:length(cts))
  {
    cts2 = cts[-which(cts == cts[j])]
    for (k in 1:length(cts2))
    {
      dd[[k]] = tryCatch(Seurat::FindMarkers(ctrl, ident.1 = cts[j], ident.2 = cts2[k], logfc.threshold = threshold), error = function(e) {FALSE})
      if (class(dd[[k]]) != "logical") {
        dd[[k]]$ident.1 = cts[j]
        dd[[k]]$ident.2 = cts2[k]
        dd[[k]]$gene = rownames(dd)
        
        flag = any(dd[[k]]$avg_logFC > 0) & any(dd[[k]]$avg_logFC < 0)
        
        if (!flag)
          dd[[k]] = data.frame(ident.1 = "NONE", ident.2 = "NONE")
        
      } else {
        dd[[k]] = data.frame(ident.1 = "NONE", ident.2 = "NONE")
      }
    }
    mrks[[j]] = dd
  }
  names(mrks) <- cts
  
  return(mrks)
}

#' Get KNN edges from single cell data
#'
#' @param E Expression matrix with genes for rows, samples for columns
#' @param normalize Normalize expression matrix to mean counts per cell. Default is FALSE.
#' @param min_counts minimum number of counts per cell. Default is 3.
#' @param min_cells minimum bumber of cells expressing at least min_counts. Default is 3.
#' @param min_vscore_pctl Minimum v score percentile for genes to run with PCA. Default is 90.
#' @param num_pc Number of PCs to build the KNN graph. Default is 50. 
#' @param k_neigh k parameter in KNN. Default is 4.
#' @param genes_use if desired, manually set the genes for PCA
#' @return List where each element contains the DEG tables for the one vs. all comparison
#' @export
CID.GetNeighbors <- function(E, normalize = F, min_counts = 3, min_cells = 3, min_vscore_pctl = 85, num_pc = 30, k_neigh = 4, genes_use = NULL)
{
  cat(" ..........  Entry in CID.GetNeighbors \n");
  ta = proc.time()[3];
  
  if (normalize)
  {
    cat('             Normalizing \n')
    E = CID.Normalize(E)
  }
  
  if (is.null(genes_use))
  {
    cat('             Filtering genes \n')
    
    # Get gene stats (above Poisson noise, i.e. V-scores)
    outs = get_vscores_sparse(E)
    gene_ix = outs$gene_ix
    Vscores = outs$v_scores
    ix2 = Vscores>0
    gene_ix = gene_ix[ix2]
    
    # Filter genes: minimum V-score percentile and at least min_counts in at least min_cells
    min_log_vscore = quantile(log(Vscores), min_vscore_pctl/100)
    
    ix = (Matrix::rowSums(E[gene_ix,] >= min_counts) >= min_cells) & (log(Vscores) >= min_log_vscore)
    cat('             Using ', sum(ix), "genes \n")
    
    genes_use = rownames(E)[ix]
  }
  
  outs = get_knn_graph2(E, k=k_neigh, np=num_pc,genes_to_use=genes_use)
  
  outs = igraph::graph.adjacency(outs, diag = F)
  outs = igraph::get.data.frame(outs)
  outs = data.frame(apply(outs, 2, as.numeric))
  colnames(outs) <- c("V1", "V2")
  
  tb = proc.time()[3] - ta;
  cat("\n ..........  Exit CID.GetNeighbors.\n");
  cat("             Execution time = ", tb, " s.\n", sep = "");
  
  return(outs)
}

#' Get V scores
#'
#' @param E Expression matrix with genes for rows, samples for columns
#' @param min_mean Minimum mean gene expression. Default is zero.
#' @param nBins Number of bins for histogram binning. Default is 50.
#' @param fit_percentile Percentile for fitting. Default is 0.1.
#' @param error_wt Error for convergence of optimization. Default is 1.
#' @return V scores for each gene in the expression matrix E.
get_vscores_sparse <- function(E, min_mean=0, nBins=50, fit_percentile=0.1, error_wt=1)
{ 
  ncell = ncol(E)
  
  mu_gene = Matrix::rowMeans(E)
  gene_ix = mu_gene > min_mean
  mu_gene = mu_gene[gene_ix]
  tmp = E[gene_ix,]
  var_gene = Matrix::rowMeans(tmp ^ 2) - mu_gene ^ 2
  FF_gene = var_gene / mu_gene
  
  data_x = log(mu_gene)
  data_y = log(FF_gene / mu_gene)
  
  outs = runningquantile(data_x, data_y, fit_percentile, nBins)
  x = outs$xOut
  y = outs$yOut
  x = x[!is.na(x)]
  y = y[!is.na(y)]
  
  gLog <- function(input) {
    log(input[[2]] * exp(-input[[1]]) + input[[3]])
  }
  d = hist(log(FF_gene[mu_gene>0]), breaks=seq(min(log(FF_gene[mu_gene>0])), max(log(FF_gene[mu_gene>0])),length.out = 201), plot = F, warn.unused = F)
  h = d$counts
  b = d$breaks
  b = b[-length(b)] + diff(b)/2
  max_ix = which.max(h)
  c = max(exp(b[max_ix]))
  errFun <- function(b2) {
    sum(abs(gLog(list(x,c,b2))-y) * error_wt)
  }
  b0 = 0.1
  b = optimize(f = errFun, interval = c(0,1))$minimum
  a = c / (1 + b) - 1
  v_scores = FF_gene / ((1+a)*(1+b) + b * mu_gene);
  outs = list(v_scores = v_scores, gene_ix = gene_ix)
  return(outs)
}

#' Compute variance for sparse matrices.
#'
#' @param E Expression matrix with genes for rows, samples for columns
sparse_var <- function(E)
{
  mean_gene = Matrix::rowMeans(E)
  tmp = E;
  tmp = tmp ^ 2;
  return(Matrix::rowMeans(tmp) - mean_gene ^ 2)
}

#' Compute running quantile
#'
#' @param x log of mean gene expression.
#' @param y log of fano factor divided by mean gene expression.
#' @param p fit_percentile as defined in get_vscores_sparse
#' @param nBins number of bins for histogram binning
runningquantile <- function(x, y, p, nBins)
{
  ind = order(x)
  x = x[ind]
  y = y[ind]
  
  dx = (x[length(x)] - x[1]) / nBins
  xOut = seq(x[1]+dx/2, x[length(x)]-dx/2, length.out = nBins)
  
  yOut = matrix(0, 1, length(xOut))
  
  for (i in 1:length(xOut))
  {
    ind = (x >= xOut[i]-dx/2) & (x < xOut[i]+dx/2)
    if (sum(ind) > 0)
    {
      yOut[i] = quantile(y[ind], p)
    } else {
      if (i > 0)
      {
        yOut[i] = yOut[i-1]
      } else {
        yOut[i] = NA
      }
    }
  }
  return(list(xOut = xOut, yOut = yOut))
}

#' Get KNN graph
#'
#' @param X Expression matrix with genes for rows, samples for columns, after V score filtering.
#' @param k KNN parameter. Default is 5.
#' @param np Number of PCs to build the KNN graph. Default is 50.
#' @param run_force Runs ForceAtlas2. Default is FALSE.
#' @param genes_to_use gene list for KNN graph building
get_knn_graph2 <- function(X, k=4, np, run_force = F, genes_to_use)
{
  logik = CID.IsUnique(rownames(X))
  X = X[logik,]
  colnames(X) <- 1:ncol(X)
  ctrl <- suppressWarnings(Seurat::CreateSeuratObject(X))
  ctrl <- Seurat::ScaleData(ctrl, verbose = F)
  ctrl <- Seurat::RunPCA(ctrl, features = genes_to_use, pcs.compute = np, do.print = F)
  ctrl <- Seurat::FindNeighbors(object = ctrl, reduction = "pca", dims = 1:min(c(np, 50)), k.param = k)
  if (run_force)
    P <- RunForceAtlas(ctrl)
  return(ctrl@graphs$RNA_nn)
}

#' Run ForceAtlas2
#'
#' @param Q Seurat object created by get_knn_graph2
RunForceAtlas <- function(Q)
{
  cat('Running ForceAtlas2 \n')
  outs = Q@graphs$RNA_nn
  outs = igraph::graph.adjacency(outs, diag = F)
  outs = igraph::get.data.frame(outs)
  outs = data.frame(apply(outs, 2, as.numeric))
  outs$weights = 1;
  colnames(outs) <- c("from", "to", "weights")
  p <- ForceAtlas2::layout.forceatlas2(outs, directed = FALSE, iterations = 1000, plotstep = 100)
  return(p)
}

#' Run GSVA
#'
#' @param y expression matrix
#' @param geneSets list of genesets
RunGSVA <- function(y, geneSets)
{
  wt = colnames(y)
  
  y = y[rownames(y) %in% geneSets[[1]],]
  
  p <- nrow(y) ## number of genes
  n <- ncol(y) ## number of samples
  
  cts = unique(wt);
  outs = list("")
  
  for (j in 1:length(cts))
  {
    nGrp1 <- sum(wt == cts[j]) ## number of samples in group 1
    nGrp2 <- n - nGrp1 ## number of samples in group 2
    
    ## build design matrix
    design <- cbind(sampleGroup1=1, sampleGroup2vs1=rep(0, n))
    design[,2][wt == cts[j]] <- 1
    
    ## fit linear model
    fit <- limma::lmFit(y, design)
    
    ## estimate moderated t-statistics
    fit <- limma::eBayes(fit)
    
    ## genes in set1 are differentially expressed
    # topTable(fit, coef="sampleGroup2vs1")
    outs[[j]] = limma::topTable(fit, coef="sampleGroup2vs1", number = 50)
  }
  
  names(outs) <- cts
  
  outs
  
  
  ## estimate GSVA enrichment scores for the three sets
  # gsva_es <- gsva(as.matrix(y), geneSets, mx.diff=1,parallel.sz = 16)
  
  ## fit the same linear model now to the GSVA enrichment scores
  # colnames(design) <- colnames(y)
  # fit <- lmFit(gsva_es, design)
  
  ## estimate moderated t-statistics
  # fit <- eBayes(fit)
  
  ## set1 is differentially expressed
  # topTable(fit, coef="sampleGroup2vs1")
}

#' Run DEG analysis with Seurat wrapper
#'
#' @param E Expression matrix with genes for rows, samples for columns
#' @return List where each element contains the DEG tables for the one vs. all comparison
#' @export
CID.GetLearningMarkers <- function(E, full.dataset, logfc.threshold = 1, pseudocount.use = 1)
{
  library(scater)
  # get lbls
  lbls = colnames(E)
  # set colnames
  colnames(E) <- seq(1, ncol(E))
  # Set up object
  ctrl <- Seurat::CreateSeuratObject(counts = E, project = "CID", min.cells = 0)
  ctrl <- Seurat::NormalizeData(ctrl)
  ctrl <- Seurat::AddMetaData(ctrl, metadata=lbls, col.name = "celltypes")
  ctrl <- Seurat::SetIdent(ctrl, value='celltypes')
  outs = list("")
  cts = unique(lbls)
  cts = cts[cts != "Other"]
  n = length(unique(colnames(full.dataset)))
  for (j in 1:length(cts))
  {
    dd = tryCatch(Seurat::FindMarkers(ctrl, ident.1 = cts[j], ident.2 = NULL, test.use = "MAST", min.cells.group = 0, max.cells.per.ident = 200, logfc.threshold = logfc.threshold, pseudocount.use = pseudocount.use), error = function(e) {FALSE})
    if (class(dd) == "logical")
      return(NULL)
    dd$GeneSymbol = rownames(dd)
    dd$celltype = cts[j]
    outs[[j]] = dd
  }
  names(outs) <- cts
  # DD = do.call(rbind, outs)
  
  # first remove anything with p_val_adj > 0.01
  res = lapply(outs, function(x) x[x$p_val_adj < 0.01,])
  
  # second remove anything with log FC < 0
  res = lapply(res, function(x) x[x$avg_logFC > 0,])
  
  # Third we order each data frame by log FC
  res = lapply(res, function(x) x[order(x$avg_logFC, decreasing = T),])
  
  # Now we remove duplicated markers
  q = do.call(rbind, res)
  dups = q$GeneSymbol[!CID.IsUnique(q$GeneSymbol)]
  res = lapply(res, function(x) x[!x$GeneSymbol %in% dups,])
  
  # remove any cell type with zero markers
  res = res[sapply(res, function(x) nrow(x) != 0)]
  
  # now determine sample size
  res = lapply(res, function(x){
    D = full.dataset
    xx = colnames(D)
    D = CID.Normalize(D)
    logik = rownames(D) %in% x$GeneSymbol
    D = D[logik,]
    if (!is.null(nrow(D)))
    {
      idx = match(x$GeneSymbol, rownames(D))
      D = D[idx,]
      ksd = apply(D, 1, sd)
    } else {
      ksd = sd(D)
    }
    # get power estimate
    q = t(scale(Matrix::t(D)))
    q = Matrix::t(Matrix.utils::aggregate.Matrix(Matrix::t(q), factor(xx), FUN = 'mean'))
    P = apply(q, 1, max)
    s = 1 + 2 * 14.88 * n * (ksd / P) ^ 2
    x$samplesize = s
    s = cumsum(1 / x$samplesize)
    idx = which(s < 1);
    idx = c(idx, length(idx) + 1)
    na.omit(x[idx,])
  })
  
  q = do.call(rbind, res)
  
  geneset = data.frame(genes = q$GeneSymbol, celltype = q$celltype, Polarity = "+")
  colnames(geneset) <- c("HUGO symbols", "Cell population", "Polarity")
  
  # first remove anything with p_val_adj > 0.01
  res = lapply(outs, function(x) x[x$p_val_adj < 0.01,])
  
  # second remove anything with log FC > 0
  res = lapply(res, function(x) x[x$avg_logFC < 0,])
  
  # Third we order each data frame by log FC
  res = lapply(res, function(x) x[order(x$avg_logFC),])
  
  # Now we remove duplicated markers
  q = do.call(rbind, res)
  dups = q$GeneSymbol[!CID.IsUnique(q$GeneSymbol)]
  res = lapply(res, function(x) x[!x$GeneSymbol %in% dups,])
  
  # remove any cell type with zero markers
  res = res[sapply(res, function(x) nrow(x) != 0)]
  
  # now determine sample size
  res = lapply(res, function(x){
    if (nrow(x) > 1)
    {
      D = full.dataset
      xx = colnames(D)
      D = CID.Normalize(D)
      logik = rownames(D) %in% x$GeneSymbol
      D = D[logik,]
      if (!is.null(nrow(D)))
      {
        idx = match(x$GeneSymbol, rownames(D))
        D = D[idx,]
        ksd = apply(D, 1, sd)
      } else {
        ksd = sd(D)
      }
      # get power estimate
      q = t(scale(Matrix::t(D)))
      q = Matrix::t(Matrix.utils::aggregate.Matrix(Matrix::t(q), factor(xx), FUN = 'mean'))
      P = apply(q, 1, min)
      s = 1 + 2 * 14.88 * n * (ksd / P) ^ 2
      x$samplesize = s
      s = cumsum(1 / x$samplesize)
      idx = which(s < 1);
      idx = c(idx, length(idx) + 1)
      na.omit(x[idx,])
    }
    else{
      x$samplesize = 1
      x
    }
  })
  
  q = do.call(rbind, res)
  
  neg = data.frame(genes = q$GeneSymbol, celltype = q$celltype, Polarity = "-")
  colnames(neg) <- c("HUGO symbols", "Cell population", "Polarity")
  
  geneset = rbind(geneset, neg)
  geneset = geneset[order(geneset$`Cell population`), ]
  
  return(geneset)
}

grab_one_gene_dev <- function(data.dir, gene)
{
  require(hdf5r)
  data.dir = gsub("\\/$", "", data.dir, perl = TRUE);
  filename = paste(dirname(data.dir), "counts_norm_sparse_genes.hdf5", sep = "/")
  if (!file.exists(filename)) {
    stop("File not found")
  }
  infile = hdf5r::H5File$new(filename)
  E = Matrix::Matrix(0, ncol = 1, nrow = hdf5r::h5attributes(infile)$ncells, sparse = T)
  flag = tryCatch(test_data <- infile[[paste("/cell_ix", gene, sep = "/")]][], error = function(e) { FALSE})
  if (!class(flag) %in% "logical")
  {
    E[flag + 1] = infile[[paste("/counts", gene, sep = "/")]][]
    return(E)
  } else {
    genes = read.table(paste(dirname(data.dir), "/genes.txt", sep = ""), quote="\"", comment.char="", stringsAsFactors=FALSE)$V1
    logik = grepl(paste0("^", gene, "_"), genes)
    gene = genes[logik]
    flag2 = tryCatch(test_data <- infile[[paste("/cell_ix", gene, sep = "/")]][], error = function(e) { FALSE})
    if (! class(flag2) %in% "logical") {
      E[flag2 + 1] = infile[[paste("/counts", gene, sep = "/")]][]
      return(E)
    } else {
      cat("ERROR: Unable to find this gene. Can you check the spelling, and make sure that it matches the gene name in these data? \n")
    }
  }
}

grab_one_cell_dev <- function(data.dir, cell_idx)
{
  require(hdf5r)
  cell_idx = cell_idx - 1;
  data.dir = gsub("\\/$", "", data.dir, perl = TRUE);
  filename = paste(dirname(data.dir), "counts_norm_sparse_cells.hdf5", sep = "/")
  if (!file.exists(filename)) {
    stop("File not found")
  }
  infile = hdf5r::H5File$new(filename)
  E = Matrix::Matrix(0, ncol = 1, nrow = hdf5r::h5attributes(infile)$ngenes + 1, sparse = T)
  flag = tryCatch(test_data <- infile[[paste("/gene_ix", cell_idx, sep = "/")]][], error = function(e) { FALSE})
  if (!class(flag) %in% "logical")
  {
    E[flag + 1] = infile[[paste("/counts", cell_idx, sep = "/")]][]
    E = E[-length(E)]
    return(E)
  } else {
    return("Error: could not find cell idx in file")
  }
}

cor.test.p <- function(x){
  FUN <- function(x, y) cor.test(x, y)[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}