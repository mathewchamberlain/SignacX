#' Load data file from directory
#'
#' @param data.dir Directory containing matrix.mtx and genes.txt.
#' @return A sparse matrix with rownames equivalent to the names in genes.txt
#' @export
CID.LoadData <- function(spring.dir, fn = "matrix.mtx")
{
  data.dir = gsub("\\/$", "", spring.dir, perl = TRUE);
  if (! (file.exists(paste(spring.dir, "matrix.mtx", sep = "/")) & file.exists(paste(spring.dir, "genes.txt", sep = "/"))))
  data.dir = dirname(spring.dir)
  gE <- paste(data.dir,fn,sep="/")
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

#' Load count matrix from an h5 file
#'
#' @param filename directory and filename of the h5 file
#' @return count matrix with genes and barcodes
#' @export
CID.LoadH5 <- function(filename) 
{
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Please install hdf5r to read HDF5 files")
  }
  if (!file.exists(filename)) {
    stop("File not found")
  }
  infile <- hdf5r::H5File$new(filename)
  data.names <- names(infile)
  
  counts <- infile[["data"]]
  indices<- infile[["indices"]]
  indptr <- infile[["indptr"]]
  shp    <- infile[["shape"]]
  genes  <- infile[["genes"]]
  sparse.mat <- Matrix::t(Matrix::sparseMatrix(i = indices[] + 1, p = indptr[], 
                                               x = as.numeric(counts[]), dims = shp[], giveCsparse = FALSE))
  sparse.mat <- as(object = sparse.mat, Class = "dgCMatrix")
  rownames(sparse.mat) <- genes[]
  infile$close_all()
  return(sparse.mat)
}

#' Load edges from edge list
#'
#' @param spring.dir A directory where "edges.csv" is located
#' @return The edgelist in data frame format
#' @export
CID.LoadEdges <- function(spring.dir)
{
  spring.dir = gsub("\\/$", "", spring.dir, perl = TRUE);
  edges = paste(spring.dir, "edges.csv", sep = "/")
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
ne <- function(E)
{
m = Matrix::Matrix(0, ncol(E), ncol(E))
tots_use = Matrix::colSums(E)
target_mean = mean(tots_use)
diag(m) <- target_mean / tots_use
return(E %*% m)
}

#' Chunk a dataset
#'
#' @param E A gene-by-sample count matrix (sparse matrix, matrix, or data.frame).
#' @param chunk.dir A directory where the chunked matrices will be stored. If the directory does not exist, it will be created.
#' @return A directory with chunked matrix
#' @export
#'
CID.Chunk <- function(E, chunk.dir, number_of_chunks = 10)
{
  chunk.dir = gsub("\\/$", "", chunk.dir, perl = TRUE);
  data(markers)
  data(cellstate_markers)
  temp = markers[markers$`Cell population` == "Platelets",]
  temp = temp[1:3,]
  markers = markers[markers$`Cell population` != "Platelets",]
  markers = rbind(markers, temp)
  genes = do.call(rbind, cellstate_markers)
  genes.ind <- which(rownames(E) %in% unique(c(as.character(markers$`HUGO symbols`), as.character(genes$`HUGO symbols`))))
  E = E[genes.ind,]
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
#' @param do.par Boolean. If true, imputation is performed in parallel on half of the machines available cores. Default = FALSE.
#' @return imputed expression matrix with only marker genes in rows.
#' @export
#'
CID.Impute <- function(E, data.dir = NULL, do.par = TRUE)
{
  # SAVER wrapper
  # If imputation was performed already and we want to use it, load imputed matrix
  if (!is.null(data.dir))
  {
    data.dir = gsub("\\/$", "", data.dir, perl = TRUE);
    fn = paste(data.dir, "matrix_saver_imputed.mtx", sep = "/")
    flag = file.exists(fn)
  }
  if (flag) {
    I = Matrix::readMM(fn)
    fn <-"genes_saver_imputed.txt"
    gG <- paste(data.dir,fn, sep = "/")
    genes <- read.csv(gG, stringsAsFactors = F, sep = "")$x
    flag = length(genes) %in% c(nrow(I), ncol(I));
    if (!flag) {
      cat("ERROR: from CID.Impute:\n");
      cat("length of genes in genes.txt = ", length(genes), " is not equal to nrow(E) = ", nrow(I), "or ncol(E) = ", ncol(I), "\n", sep = "");
      stop()
    }
    if (nrow(I) != length(genes))
      I = Matrix::t(I)
    rownames(I) <- genes
    return(I)
  }  else {
    data(markers)
    data(cellstate_markers)
    genes = do.call(rbind, cellstate_markers)
    genes.ind <- which(rownames(E) %in% unique(c(as.character(markers$`HUGO symbols`), as.character(genes$`HUGO symbols`))))
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
  }
    return(I)
}

#' runs CID.CellID in batch mode
#'
#' @param E A list of expression matrices with features (genes) in rows and samples (cells) in columns.
#' @return A list where each element is a CID.CellID result. (See ?CID.CellID)
#' @export
CID.BatchMode <- function(E,f,pval,deep_dive,edges,entropy,sorted,walktrap)
{
  cat(" ..........  Running CID.CellID in batch mode on input matrices: \n");
  cat("             Detected N = ", NROW(E), " matrices \n", sep = "" );
  if (!is.null(edges))
    cr = mapply(function(x,y) CID.CellID(E = x, edges = y, deep_dive = deep_dive, f = f, pval = pval), x = E, y = edges, SIMPLIFY = FALSE)
  else
    cr = lapply(E, function(x) CID.CellID(E = x, edges = edges, deep_dive = deep_dive, f = f, pval = pval))
  if (is.null(names(E)))
    names(cr) <- paste("x", seq_along(E), sep = "")
  return(cr)
}

#' Main function
#'
#' @param E A gene-by-sample count matrix (sparse matrix, matrix, or data.frame) with genes identified by their HUGO symbols (see ?CID.geneconversion), or a list of such matrices, see ?CID.BatchMode.
#' @param pval p-value cutoff for feature selection, as described in the manuscript + markdown file. Default is pval = 0.1.
#' @param deep_dive Boolean, T will assign cell types and then cell states, F will only assign cell types. Default is deep_dive = T.
#' @return Filtered markers where each marker must have at least ncells that express at least ncounts
#' @export
CID.CellID <- function(E,pval = 0.1,deep_dive = TRUE,spring.dir = NULL, entropy = FALSE, walktrap = FALSE, omit = NULL, sorted = FALSE)
{
  # load markers
  data(markers)
  data(cellstate_markers)
  if (!length(markers) > 0) {
    cat("ERROR: from CID:\n");
    cat("required markers failed to load.\n", sep = "");
    stop()
  }
  
  # if list, run batch mode
  if(class(E) == "list")
    return(CID.BatchMode(E = E,pval = pval,deep_dive = deep_dive,spring.dir = spring.dir, entropy = entropy, sorted = sorted, walktrap = walktrap, omit = omit))

  # check inputs
  cat(" ..........  Entry in CID.CellID \n");
  ta = proc.time()[3];
  stopifnot(class(E) %in% c("dgCMatrix","dgTMatrix", "matrix", "data.frame"))
  stopifnot(!is.null(rownames(E)));

  # main function
  cat(" ..........  Computing Signac scores for cell types on input data matrix :\n");
  cat("             nrow = ", nrow(E), "\n", sep = "");
  cat("             ncol = ", ncol(E), "\n", sep = "");
  filtered_features = CID.filter(E, markersG  = markers, pval = pval)
  
  if (!is.null(omit))
  {
    if (sum(omit %in% names(filtered_features)) > 0)
    {
      omit = omit[omit %in% names(filtered_features)]
      filtered_features = filtered_features[-which(names(filtered_features) %in% omit)]
    }
  }
  
  dfY = CID.append(E,filtered_features, sorted)
  
  # assign output classifications
  cat(" ..........  Assigning output classifications \n", sep ="");
  indexMax = apply(dfY, 2, which.max);
  ac = rownames(dfY)[indexMax];
  
  # compute distance matrix
  if (!is.null(spring.dir) | entropy | walktrap)
    distMat = CID.GetDistMat(spring.dir)

  # smooth the output classifications
  if (!is.null(spring.dir))
  {
    cat(" ..........  Smoothing \n");
    # acOut_knn_smooth = CID.smooth(ac, distMat[[1]]) # smooth based on edges in knn graph
    acOut_knn_smooth = CID.smooth(ac, distMat[[1]])
    cat("\n ..........  Smoothing completed! \n");
  }

  if (entropy)
  {
    cat(" ..........  Assigning Others \n");
    acOut_knn_smooth = CID.entropy(acOut_knn_smooth, distMat)
    cat("\n ..........  Done! \n");
  }

  # cell state deep dive classifications
  if (deep_dive)
  {
    dummy = list("")
    cat(" ..........  Computing CID scores for cell states! \n");
    if (!is.null(spring.dir))
    {  ac_dd = acOut_knn_smooth
    } else {
      ac_dd = ac
    }
    for (j in names(cellstate_markers)){
      logik = ac_dd == j; sum(logik)
      if (sum(logik) != 0){
        filtered_features = CID.filter(E, markersG  = cellstate_markers[[which(names(cellstate_markers) == j)]])
      } else {
        filtered_features = NULL
      }
      if (!is.null(filtered_features))
      {
        ac_dd = CID.deepdive(filtered_features, j, acOut3 = ac_dd, expression = E, f = f, sorted = sorted)
        dummy[[j]] = ac_dd;
      }
    }
    dummy = dummy[sapply(dummy, function(x) length(x) != 1)]
    names(dummy) <- paste("L", seq_along(1:length(dummy)), sep = "")
    cat(" ..........  Deep dive completed!\n");
    if (!is.null(spring.dir))
    {
      cat(" ..........  Smoothing \n");
      ac_dd_knn = CID.smooth(ac_dd, distMat[[1]])
    }
  }

  if (walktrap)
  {
    cat("\n")
    wt = CID.Louvain(spring.dir)
    do = data.frame(table(wt[acOut_knn_smooth == "Other"]))
    df = data.frame(table(wt[wt %in% do$Var1]))
    logik = (1 - phyper(do$Freq, sum(acOut_knn_smooth == "Other") , length(ac) - sum(acOut_knn_smooth == "Other"), df$Freq)) < 0.01;
    if (sum(logik) > 0)
    {
      do = do[logik,]
      logik = do$Freq > 20; # require at least 20 cell communities
      if (sum(logik) > 0) 
      {
        lbls = rep("All", ncol(E))
        logik = wt %in% do$Var1[logik];
        lbls[logik] = wt[logik]
        lbls[acOut_knn_smooth != "Other"] = "All"
        colnames(E) <- lbls
        acOut_knn_smooth_wt = CID.PosMarkers(E, acOut_knn_smooth)
        ac_dd_knn_wt = ac_dd_knn
        logik = grepl("+", acOut_knn_smooth_wt)
        ac_dd_knn_wt[logik] = acOut_knn_smooth_wt[logik]
      }
    }
  }
  
  # Package output
  if (deep_dive & !is.null(spring.dir) & walktrap & sum(logik) != 0) {
    cr = list(scores = dfY,
              ctypes = ac,
              ctypessmoothed = acOut_knn_smooth,
              ctypessmoothed_dev = acOut_knn_smooth_wt,
              ddtypes = ac_dd,
              ddtypessmoothed = ac_dd_knn,
              ddtypessmoothed_dev = ac_dd_knn_wt,
              hierarchylabels = dummy,
              walktrap = wt)
  } else if (deep_dive & !is.null(spring.dir) & walktrap)
  {
    cr = list(scores = dfY,
              ctypes = ac,
              ctypessmoothed = acOut_knn_smooth,
              ddtypes = ac_dd,
              ddtypessmoothed = ac_dd_knn,
              hierarchylabels = dummy,
              walktrap = wt)
  }
  else if (deep_dive & !is.null(spring.dir)) {
    cr = list(scores = dfY,
              ctypes = ac,
              ctypessmoothed = acOut_knn_smooth,
              ddtypes = ac_dd,
              ddtypessmoothed = ac_dd_knn,
              hierarchylabels = dummy)
  }
  else if (!is.null(spring.dir)) {
    cr = list(scores = dfY,
              ctypes = ac,
              ctypessmoothed = acOut_knn_smooth);}
  else if (deep_dive & is.null(spring.dir)) {
    cr = list(scores = dfY,
              ctypes = ac ,
              ddtypes = ac_dd ,
              hierarchylabels = dummy)}
  else {
    cr = list(scores = dfY,
              ctypes = ac)
  }
  tb = proc.time()[3] - ta;
  cat("\n ..........  Exit CID.CellID.\n");
  cat("             Execution time = ", tb, " s.\n", sep = "");
  return (cr);
}

#' filters the geneset markers
#'
#' @param expression An expression matrix with features (genes) in rows and samples (cells) in columns.
#' @param markersG A data frame with four columns ('HUGO Symbol', 'Cell Population', 'ENTREZ ID', 'Polarity'). Default is internally set to data(markers_v4).
#' @param pval p-value cutoff, as described in the manuscript. Default is 0.1.
#' @return A list of markers and features.
#' @export
CID.filter <- function(expression, markersG, pval = 0.1)
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
    filter = summary(unlist(Q))[2]
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
#' @return A matrix of CellID scores; features by cells.
#' @export
CID.append=function(expression,featuresG, sorted)
{
  if (!sorted)
  {
    # subset expression matrix
    Z = expression[row.names(expression) %in% as.character(unlist(lapply(featuresG, function(x) x = x[,1]))), ]
    # transform expression -> log base 2 (expression + 1), and then Z-score transform
    Z = log(Z + 1, 2)
    Z = Matrix::t(scale(Matrix::t(Z)))
    # Get CellID scores
    res = as.data.frame(Matrix::t(do.call(cbind,
                                          lapply(featuresG,function(x){
                                            if (sum(x$Polarity == "-") > 0)
                                            {
                                              apply(Z[intersect(row.names(Z),x$`HUGO symbols`[x$Polarity == "+"]),,drop=F],2,function(x) mean(x, na.rm = T)) -
                                                apply(Z[intersect(row.names(Z),x$`HUGO symbols`[x$Polarity == "-"]),,drop=F],2,function(x) mean(x, na.rm = T))
                                            } else
                                            {
                                              apply(Z[intersect(row.names(Z),x$`HUGO symbols`[x$Polarity == "+"]),,drop=F],2,function(x) mean(x, na.rm = T))
                                            }
                                          }))))
    res = res[order(rownames(res)),]
    res
  } else {
    # subset expression matrix
    Z = expression[row.names(expression) %in% as.character(unlist(lapply(featuresG, function(x) x = x[,1]))), ]
    # transform expression -> log base 2 (expression + 1)
    Z = log(Z + 1, 2)
    # Get CellID scores
    res = as.data.frame(Matrix::t(do.call(cbind,
                                          lapply(featuresG,function(x){
                                              apply(Z[intersect(row.names(Z),x$`HUGO symbols`[x$Polarity == "+"]),,drop=F],2,function(x) mean(x, na.rm = T))
                                          }))))
    res = res[order(rownames(res)),]
    res
  }
}

CID.deepdive <- function(XX, YY, acOut3 = ac_dd, expression = E, f = f, sorted = sorted)
{
  cat(" ..........  Assigning cell states :\n");
  logik = acOut3 == YY; sum (logik)
  if (sum(logik) != 0 ){
    cat("             Hierarchy assignment: ", sum(logik), " ", YY, "\n", sep ="");
    ## get CID scores for cell types
    dfYe = CID.append(expression,XX, sorted)
    ar = rownames(dfYe)
    dfYe = dfYe[,logik]
    dfYe = na.omit(dfYe)
    if (!is.null(dfYe) & length(dfYe) > 1){
      indexMax = apply(as.matrix(dfYe), 2, which.max);
      dummy = ar[indexMax];
      acOut3[logik] = dummy;}
  } else
  {
    return(acOut3)
  }
  acOut3
}

#' Smoothing function
#'
#' @param ac List containing a character vector where each element is a cell type or cell state assignment
#' @param f Smoothing parameter
#' @param edges Location where knn edges are stored
#' @return Smoothed cell type or cell state assignments
#' @export
CID.smooth <- function(ac,dM)
{
  # pre-allocate
  Y   = ac # cells we may switch labels
  m = Matrix::Matrix(0, nrow = length(ac), ncol = length(unique(ac)), sparse = T)
  m[cbind(1:nrow(m), as.numeric(factor(ac)))] <- 1
  res = dM %*% m
  res = res / Matrix::rowSums(res)
  mx.idx = apply(res, 1, which.max)
  mx = apply(res, 1, max)
  Y[mx > 0.5] = levels(factor(ac))[mx.idx[mx > 0.5]]
  return(Y)
}

#' Distance matrix
#'
#' @param ac A character vector of cell type labels
#' @param edges A data frame or matrix with edges
#' @return The entropy for each node in the network
#' @export
CID.DistMatrix <- function(ac, spring.dir)
{
  # load edges
  edges    = CID.LoadEdges(spring.dir)
  # get distance matrix
  g = igraph::graph_from_edgelist(as.matrix(edges), directed = FALSE)
  igraph::shortest.paths(g, v=igraph::V(g), to=igraph::V(g))
}

#' Entropy
#' @export
CID.entropy <- function(ac,dM)
{
  # Pre-allocate cells for shannon calculation
  Y   = ac
  N_unique = length(unique(Y))
  
  # Calculate normalized shannon entropy for each cell j for all connections with shortest path < N
  shannon = rep(0, length(Y))
  dM = Reduce('+', dM) > 0;

  m = Matrix::Matrix(0, nrow = length(ac), ncol = length(unique(ac)), sparse = T)
  m[cbind(1:nrow(m), as.numeric(factor(ac)))] <- 1
  
  res = dM %*% m
  res = res / Matrix::rowSums(res)
  shannon = apply(res, 1, function(freqs) {-sum(freqs[freqs != 0] * log(freqs[freqs != 0])) / log(2) / log2(N_unique)})

  Y[shannon > (mean(shannon) + 3 * sd(shannon))] = "Other"
  
  return(Y)
}

#' Write JSON file for SPRING visualization
#'
#' @param cr Output from CID.CellID. See ?CID.CellID
#' @param json_new Filename for new SPRING visualization. Default is json_new = "categorical_coloring_data_new.json".
#' @param spring.dir Directory where file 'categorical_coloring_data.json' is located. If supplied, it will append this file to contain tracks for cell type / state classifications.
#' @return Smoothed cell type or cell state assignments
#' @export
CID.writeJSON <- function(cr, json_new = "categorical_coloring_data.json", spring.dir = NULL)
{
  if (!is.null(spring.dir))
  {
    spring.dir = gsub("\\/$", "", spring.dir, perl = TRUE);
    if (file.exists(paste(spring.dir, 'categorical_coloring_data_old.json', sep = "/")))
    {
      json_file = 'categorical_coloring_data_old.json'
      gJ <- paste(spring.dir,json_file,sep = "/")
      json_data <- rjson::fromJSON(file=gJ)
    } else {
      json_file = 'categorical_coloring_data.json'
      gJ <- paste(spring.dir,json_file,sep = "/")
      json_data <- rjson::fromJSON(file=gJ)
      json_out = jsonlite::toJSON(json_data, auto_unbox = TRUE)
      write(json_out,paste(spring.dir,'categorical_coloring_data_old.json',sep="/"))
    }
  }
  if ("walktrap" %in% names(cr))
  {
    Q = as.character(cr$walktrap)
    json_data$ClustersWT$label_list = Q
    Ntypes = length(unique(Q))
    qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) #len = 74
    #pie(rep(1,num_col), col=(col_vector[1:num_col]))
    col_palette <- as.list(col_vector[1:Ntypes]); # or sample if you wish
    names(col_palette) <- unique(Q)
    json_data$ClustersWT$label_colors = col_palette
  }
  if ("ctypes" %in% names(cr))
  {
    Q = cr$ctypes
    json_data$CellTypesID$label_list = Q
    C = get_colors(Q)
    json_data$CellTypesID$label_colors = as.list(C[[1]])
  }
  if ("ctypessmoothed" %in% names(cr))
  {
    Q = cr$ctypessmoothed
    json_data$CellTypesID$label_list = Q
    C = get_colors(Q)
    json_data$CellTypesID$label_colors = as.list(C[[1]])
  }
  if ("ddtypes" %in% names(cr))
  {
    Q = cr$ddtypes
    json_data$CellStatesID$label_list = Q
    C = get_colors(Q)
    json_data$CellStatesID$label_colors = as.list(C[[1]])
  }
  if ("ddtypessmoothed" %in% names(cr))
  {
    Q = cr$ddtypessmoothed
    json_data$CellStatesID$label_list = Q
    C = get_colors(Q)
    json_data$CellStatesID$label_colors = as.list(C[[1]])
  }
  if ("ctypessmoothed_dev" %in% names(cr))
  {
    Q = cr$ctypessmoothed_dev
    json_data$CellTypesID_dev$label_list = Q
    C = get_colors(Q)
    Ntypes = sum(C[[1]] == "")
    qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) #len = 74
    #pie(rep(1,num_col), col=(col_vector[1:num_col]))
    C[[1]][C[[1]] == ""] <- col_vector[1:Ntypes]; # or sample if you wish
    json_data$CellTypesID_dev$label_colors = as.list(C[[1]])
  }
  if ("ddtypessmoothed_dev" %in% names(cr))
  {
    Q = cr$ddtypessmoothed_dev
    json_data$CellStatesID_dev$label_list = Q
    C = get_colors(Q)
    Ntypes = sum(C[[1]] == "")
    qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) #len = 74
    #pie(rep(1,num_col), col=(col_vector[1:num_col]))
    C[[1]][C[[1]] == ""] <- col_vector[1:Ntypes]; # or sample if you wish
    json_data$CellStatesID_dev$label_colors = as.list(C[[1]])
  }
  json_data = json_data[order(names(json_data))]
  json_data_backup = json_data;
  fn = json_new
  if (is.null(spring.dir))
  {
    fn = json_new;
    spring.dir = getwd()
    spring.dir = gsub("\\/$", "", spring.dir, perl = TRUE);
  }
  json_out = jsonlite::toJSON(json_data, auto_unbox = TRUE)
  write(json_out,paste(spring.dir,fn,sep="/"))
  cat(paste(spring.dir,fn,sep="/"), "has been written to directory! \n")
  new.dirs = list.dirs(dirname(spring.dir))
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
                 "Erythro"     ,     "Platelets")
  main_colors = c("#aa3596", "#387f50", "#a3ba37",
                  "#D3D3D3", "#f0ff51", "#d878de",
                  "#8c933e", "#90c5f4", "#90b771",
                  "#C0C0C0", "#dbd0d0", "#e1f7d5",
                  "#ffbdbd", "#c9c9ff")

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
                 "NK.cells")
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
                 "#ad9bf2")

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
#' @param ac A character vector of cell type labels
#' @param edges A data frame or matrix with edges
#' @return The community substructures of the graph
#' @export
CID.Louvain <- function(edges)
{
  edges = CID.LoadEdges(edges)
  # convert edgelist to adjacency matrix
  adjmatrix <- methods::new("ngTMatrix", 
                    i = c(as.integer(edges$V1)-1L, as.integer(edges$V2)-1L), 
                    j = c(as.integer(edges$V2)-1L, as.integer(edges$V1)-1L),
                    Dim = as.integer(c(max(edges), max(edges))))
  g = igraph::graph_from_adjacency_matrix(adjmatrix * 1, mode = c("undirected"), weighted = NULL, diag = TRUE, add.colnames = NULL, add.rownames = NA)
  as.character(igraph::cluster_louvain(g)$membership)
}

#' Get distance matrix
#'
#' @param spring.dir directory where edges.csv is located
#' @param n maximum network distance to subtend (n neighbors)
#' @return adjacency matrices for distances < n
#' @export
CID.GetDistMat <- function(spring.dir, n = 4)
{
  "%^%" <- function(A, n) {if(n == 1) A else A %*% (A %^% (n-1)) }
  if(class(spring.dir) == "character") edges = CID.LoadEdges(spring.dir) else edges = spring.dir
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
#' @param E Count matrix (see ?Read_h5)
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

#' Run DEG analysis with Seurat wrapper
#'
#' @param E Expression matrix with genes for rows, samples for columns
#' @return List where each element contains the DEG tables for the one vs. all comparison
#' @export
CID.PosMarkers <- function(E, acn)
{
  # get lbls
  lbls = colnames(E)
  # set colnames
  colnames(E) <- seq(1, ncol(E))
  # make sure row names are not redundant
  logik = CID.IsUnique(rownames(E))
  E = E[logik,]
  # Set up object
  ctrl <- Seurat::CreateSeuratObject(counts = E)
  ctrl <- Seurat::NormalizeData(object = ctrl)
  ctrl <- Seurat::AddMetaData(ctrl, metadata=lbls, col.name = "celltypes")
  ctrl <- Seurat::SetIdent(ctrl, value='celltypes')
  #outs = list("")
  cts = unique(lbls);
  cts = cts[cts != "All"]
  for (j in 1:length(cts))
  {
    dd = Seurat::FindMarkers(ctrl, ident.1 = cts[j], ident.2 = NULL, min.cells.group = 0, max.cells.per.ident = 200, logfc.threshold = 1, pseudocount.use = 1, min.pct = 0)
    dd = dd[dd$p_val_adj < 0.01,]
    if (sum(dd$p_val_adj < 0.01) > 0)
    {
      gns = rownames(dd)[order(dd$avg_logFC, decreasing = TRUE)][1:2]
      gns = gns[!is.na(gns)]
      if (length(gns == 2))
      {
        acn[lbls == cts[j]] = rep( paste ("+", gns, collapse = " ", sep = ""), sum(lbls == cts[j]))
      } else {
        acn[lbls == cts[j]] = rep( paste ("+", gns, sep = ""), sum(lbls == cts[j]))
      }
    }
  }
    
    #dd = Seurat::FindMarkers(ctrl, ident.1 = cts[j], ident.2 = NULL, min.cells.group = 0, max.cells.per.ident = 200, logfc.threshold = 0, pseudocount.use = 1, min.pct = 0)
    #dd = Seurat::FindMarkers(ctrl, ident.1 = cts2$ident.1[j], ident.2 = cts2$ident.2[j], min.cells.group = 0, max.cells.per.ident = 200, logfc.threshold = 0, pseudocount.use = 0, min.pct = 0)
    #  dd$GeneSymbol = rownames(dd)
    #  dd$celltype = gs  " .*$", "", cts2$ident.1 )[j]
    #  outs[[j]] = dd
    #dd = Seurat::Fin  arkers(ctrl, ident.1 = "wt", ident.2 = "All", min.cells.group = 0, max.cells.per.ident = 200, logfc.threshold = 0, pseudocount.use = 1, min.pct = 0)
  #dd = dd[dd$p_val_adj < 0.01,]
  #acn[lbls == "wt"] = rep(paste( paste ("+", rownames(dd)[order(dd$avg_logFC, decreasing = TRUE)][1:2], collapse = " ", sep = ""), "cells", sep = " "), sum(lbls == "wt"))
  acn
  #}
  #names(outs) <- gsub( " .*$", "", cts2$ident.1 )
  #outs = do.call(rbind, outs)
  #outs = do.call(rbind, lapply(outs, function(x) x[x$p_val_adj < 0.05 & x$avg_logFC > 0.2, ]))
  #geneset = data.frame(genes = outs$GeneSymbol, celltype = outs$celltype, Polarity = "+")
  #colnames(geneset) <- c("HUGO symbols", "Cell population", "Polarity")
  #return(geneset)
}
  