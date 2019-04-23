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

#' Load edges from edge list
#'
#' @param data.dir A directory where "edges.csv" is located
#' @return The edgelist in data frame format
#' @export
CID.LoadEdges <- function(data.dir)
{
  data.dir = gsub("\\/$", "", data.dir, perl = TRUE);
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
#' @param number_of_chunks the number of chunks (e.g., the number of sub-sampled matrices after chunking)
#' @return A directory with chunked matrix
#' @export
#'
CID.Chunk <- function(E, chunk.dir, number_of_chunks = 10)
{
  chunk.dir = gsub("\\/$", "", chunk.dir, perl = TRUE);
  data('markers')
  data('cellstate_markers')
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
#' @param data.dir directory for saving "matrix_saver_imputed.mtx" for future loading.
#' @param do.par Boolean. If true, imputation is performed in parallel on half of the machines available cores. Default = FALSE.
#' @return imputed expression matrix with only marker genes in rows.
#' @export
#'
CID.Impute <- function(E = NULL, data.dir = NULL, do.par = TRUE)
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
    data('markers')
    data('cellstate_markers')
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

#' Main function
#'
#' @param E A gene-by-sample count matrix (sparse matrix, matrix, or data.frame) with genes identified by their HUGO symbols (see ?CID.geneconversion), or a list of such matrices, see ?CID.BatchMode.
#' @param pval p-value cutoff for feature selection, as described in the manuscript + markdown file. Default is pval = 0.1.
#' @param data.dir directory of SPRING files "edges.csv" and "categorical_coloring_data.json"
#' @param entropy cells amended to high entropy labels with respect to their neighbors in the KNN graph are appended "Other" if entropy = TRUE. Default is TRUE.
#' @param louvain Louvain community detection is performed, and then used together with Shannon entropy to detect potential novel cell types / states. Default is TRUE.
#' @param omit Force remove specific cell types / states with omit. Default is NULL.
#' @param sorted If cells are expected to be pure or mostly homogeneous (e.g., by FACs sorting), set sorted = TRUE. Default is FALSE.
#' @return Filtered markers where each marker must have at least ncells that express at least ncounts
#' @export
CID.CellID <- function(E,pval = 0.1,data.dir = NULL, entropy = T, louvain = T, omit = NULL, sorted = FALSE)
{
  # load markers
  data('markers')
  data('cellstate_markers')
  if (!length(markers) > 0) {
    cat("ERROR: from Signac Data:\n");
    cat("Required markers failed to load.\n", sep = "");
    stop()
  }
  
  # check inputs
  stopifnot(class(E) %in% c("dgCMatrix","dgTMatrix", "matrix", "data.frame"))
  stopifnot(!is.null(rownames(E)));
  cat(" ..........  Entry in CID.CellID \n");
  ta = proc.time()[3];

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
      cat(" ..........  Forcibly omitting features :\n");
      cat("             Omitted = ", omit, "\n", sep = "");
    }
  }
  
  # compute cell type score data.frame
  dfY = CID.append(E,filtered_features, sorted)
  
  # assign output classifications
  cat(" ..........  Assigning output classifications \n", sep ="");
  indexMax = apply(dfY, 2, which.max);
  ac = rownames(dfY)[indexMax];
  
  # if is null data.dir, run PCA + KNN
  if (is.null(data.dir))
    data.dir = CID.GetNeighbors(E, normalize = F, min_counts = 3, min_cells = 3, min_vscore_pctl = 90, num_pc = 50, k_neigh = 4)
   
  # compute distance matrix
  if (!is.null(data.dir) | entropy | louvain)
    distMat = CID.GetDistMat(data.dir)

  # smooth the output classifications
  if (!is.null(data.dir))
  {
    acOut_knn_smooth = CID.smooth(ac, distMat[[1]])
    cat("\n ..........  Smoothing completed! \n");
  }

  if (entropy)
  {
    cat(" ..........  Assigning Others \n");
    acOut_knn_smooth = CID.entropy(acOut_knn_smooth, distMat)
  }

  # cell state deep dive classifications
  dummy = list("")
  cat(" ..........  Computing CID scores for cell states! \n");
    if (!is.null(data.dir))
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
        ac_dd = CID.deepdive(filtered_features, j, acOut3 = ac_dd, expression = E, sorted = sorted)
        dummy[[j]] = ac_dd;
      }
    }
    dummy = dummy[sapply(dummy, function(x) length(x) != 1)]
    names(dummy) <- paste("L", seq_along(1:length(dummy)), sep = "")
    cat(" ..........  Deep dive completed!\n");
    if (!is.null(data.dir))
    {
      cat(" ..........  Smoothing \n");
      ac_dd_knn = CID.smooth(ac_dd, distMat[[1]])
    }

  if (louvain)
  {
    cat(" ..........  Running Louvain clustering\n");
    wt = CID.Louvain(edges = data.dir)
    do = data.frame(table(wt[acOut_knn_smooth == "Other"]))
    df = data.frame(table(wt[wt %in% do$Var1]))
    logik = (1 - phyper(do$Freq, sum(acOut_knn_smooth == "Other") , length(ac) - sum(acOut_knn_smooth == "Other"), df$Freq)) < 0.01;
    if (sum(logik) > 0)
    {
      do = do[logik,]
      logik = do$Freq > 20; # require at least 20 cell communities
      if (sum(logik) > 0) 
      {
        cat("             Signac detected n = ", sum(logik), "novel cell type populations!\n");
        lbls = rep("All", ncol(E))
        logik = wt %in% do$Var1[logik];
        lbls[logik] = wt[logik]
        lbls[acOut_knn_smooth != "Other"] = "All"
        colnames(E) <- lbls
        cat(" ..........  Identifying markers for the populations \n");
        acOut_knn_smooth_wt = CID.PosMarkers(E, acOut_knn_smooth)
        ac_dd_knn_wt = ac_dd_knn
        logik = grepl("^[+]", acOut_knn_smooth_wt)
        ac_dd_knn_wt[logik] = acOut_knn_smooth_wt[logik]
      }
    }
  }
  
  # Package output
  if (!is.null(data.dir) & louvain & sum(logik) != 0) {
    cr = list(scores = dfY,
              ctypes = ac,
              ctypessmoothed = acOut_knn_smooth,
              ctypessmoothed_dev = acOut_knn_smooth_wt,
              ddtypes = ac_dd,
              ddtypessmoothed = ac_dd_knn,
              ddtypessmoothed_dev = ac_dd_knn_wt,
              hierarchylabels = dummy,
              louvain = wt)
  } else if (!is.null(data.dir) & louvain)
  {
    cr = list(scores = dfY,
              ctypes = ac,
              ctypessmoothed = acOut_knn_smooth,
              ddtypes = ac_dd,
              ddtypessmoothed = ac_dd_knn,
              hierarchylabels = dummy,
              louvain = wt)
  }
  else if (!is.null(data.dir)) {
    cr = list(scores = dfY,
              ctypes = ac,
              ctypessmoothed = acOut_knn_smooth,
              ddtypes = ac_dd,
              ddtypessmoothed = ac_dd_knn,
              hierarchylabels = dummy)
  }
  else if (!is.null(data.dir)) {
    cr = list(scores = dfY,
              ctypes = ac,
              ctypessmoothed = acOut_knn_smooth);}
  else if (is.null(data.dir)) {
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
#' @param sorted see ?CID.CellID
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
#' Deep cell states
#'
#' @param XX filtered features for CID.append
#' @param YY cell type where we are going deeper
#' @param acOut3 Vector of cell type labels
#' @param expression the count matrix
#' @param sorted see ?CID.CellID
#' @return Vector for the deep cell state labels
CID.deepdive <- function(XX, YY, acOut3, expression, sorted = sorted)
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
#' @param dM Distance matrix (see ?CID.GetDistMat)
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

#' Entropy
#' 
#' @param ac A character vector of cell type labels
#' @param dM The distance matrix, see ?CID.GetDistMat
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
#' @param data.dir Directory where file 'categorical_coloring_data.json' is located. If supplied, it will append this file to contain tracks for cell type / state classifications.
#' @return Smoothed cell type or cell state assignments
#' @export
CID.writeJSON <- function(cr, json_new = "categorical_coloring_data.json", data.dir = NULL)
{
  if (!is.null(data.dir))
  {
    data.dir = gsub("\\/$", "", data.dir, perl = TRUE);
    if (file.exists(paste(data.dir, 'categorical_coloring_data_old.json', sep = "/")))
    {
      json_file = 'categorical_coloring_data_old.json'
      gJ <- paste(data.dir,json_file,sep = "/")
      json_data <- rjson::fromJSON(file=gJ)
    } else {
      json_file = 'categorical_coloring_data.json'
      gJ <- paste(data.dir,json_file,sep = "/")
      json_data <- rjson::fromJSON(file=gJ)
      json_out = jsonlite::toJSON(json_data, auto_unbox = TRUE)
      write(json_out,paste(data.dir,'categorical_coloring_data_old.json',sep="/"))
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
#' @param edges A data frame or matrix with edges
#' @return The community substructures of the graph
#' @export
CID.Louvain <- function(edges)
{
  if (class(edges) == "character"){
    edges = CID.LoadEdges(data.dir)
  } else{
    edges = data.dir
  }
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
#' @param data.dir directory where edges.csv is located
#' @param n maximum network distance to subtend (n neighbors)
#' @return adjacency matrices for distances < n
#' @export
CID.GetDistMat <- function(data.dir, n = 4)
{
  "%^%" <- function(A, n) {if(n == 1) A else A %*% (A %^% (n-1)) }
  if (class(data.dir) == "character"){
    edges = CID.LoadEdges(data.dir)
  } else{
    edges = data.dir
  }

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

#' Run DEG analysis with Seurat wrapper
#'
#' @param E Expression matrix with genes for rows, samples for columns
#' @param acn Vector of cell type labels 
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
  acn
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
#' @return List where each element contains the DEG tables for the one vs. all comparison
#' @export
CID.GetNeighbors <- function(E, normalize = F, min_counts = 3, min_cells = 3, min_vscore_pctl = 90, num_pc = 50, k_neigh = 4)
{
cat(" ..........  Entry in CID.GetKNNEdges \n");
ta = proc.time()[3];

  if (normalize)
  {
    cat('             Normalizing \n')
    E = CID.Normalize(E)
  }
  
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
  
  outs = get_knn_graph2(E[ix,], k=k_neigh, np = num_pc)
  
  outs = igraph::graph.adjacency(outs, diag = F)
  outs = igraph::get.data.frame(outs)
  outs = data.frame(apply(outs, 2, as.numeric))
  colnames(outs) <- c("V1", "V2")

  tb = proc.time()[3] - ta;
  cat("\n ..........  Exit CID.GetKNNEdges.\n");
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
get_knn_graph2 <- function(X, k=5, np, run_force = F)
{
  k = k + 1;
  logik = CID.IsUnique(rownames(X))
  X = X[logik,]
  colnames(X) <- 1:ncol(X)
  ctrl <- Seurat::CreateSeuratObject(X)
  ctrl <- Seurat::ScaleData(ctrl)
  ctrl <- Seurat::RunPCA(ctrl, features = rownames(X), pcs.compute = np, do.print = F)
  ctrl <- Seurat::FindNeighbors(object = ctrl, reduction = "pca", dims = 1:min(c(np, 50)), k.param = k)
  if (force)
    RunForceAtlas(ctrl)
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
  
  p <- layout.forceatlas2(outs, directed = FALSE, iterations = num_force_iter, plotstep = 100)
  
}