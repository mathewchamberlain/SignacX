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

#' Imputation wrapper
#'
#' @param E A gene-by-sample count matrix (sparse matrix, matrix, or data.frame) with genes identified by their HUGO symbols.
#' @param do.par Boolean. If true, imputation is performed in parallel on half of the machines available cores. Default = FALSE.
#' @return imputed expression matrix with only marker genes in rows.
#' @export
#'
CID.Impute <- function(E, spring.dir = NULL, do.par = TRUE, large = FALSE)
{
  # SAVER wrapper
  # If imputation was performed already and we want to use it, load imputed matrix
  flag = FALSE
  if (!is.null(spring.dir))
  {
    data.dir = dirname(spring.dir)
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
    temp = markers[markers$`Cell population` == "Platelets",]
    temp = temp[1:3,]
    markers = markers[markers$`Cell population` != "Platelets",]
    markers = rbind(markers, temp)
    genes = do.call(rbind, cellstate_markers)
    genes.ind <- which(rownames(E) %in% unique(c(as.character(markers$`HUGO symbols`), as.character(genes$`HUGO symbols`))))
  if (do.par)
  { 
    numCores = parallel::detectCores()
    if (!large) {
      I = SAVER::saver(E, pred.genes = genes.ind, pred.genes.only = T, ncores = numCores / 2, estimates.only = T)
    } else {
      q = split(1:nrow(E), ceiling(seq_along(1:nrow(E))/1500))
      L = lapply(q, function(x) {
        SAVER::saver(as.matrix(E), pred.genes = x, pred.genes.only = TRUE, ncores = numCores / 2, do.fast = FALSE, estimates.only = T)
      }
      )
      I <- SAVER::combine.saver(L)
       }
    if (!is.null(spring.dir))
    {
      data.dir = dirname(spring.dir)
      Matrix::writeMM(Matrix::Matrix(I, sparse = TRUE), file = paste(data.dir, "matrix_saver_imputed.mtx", sep = "/"))
      write.table(rownames(E)[genes.ind], file = paste(data.dir, "genes_saver_imputed.txt", sep = "/"))
    }  
  } else {
    I = SAVER::saver(E, pred.genes = genes.ind, pred.genes.only = T, estimates.only = T)
    if (!is.null(spring.dir))
    {
      data.dir = dirname(spring.dir)
      Matrix::writeMM(Matrix::Matrix(I, sparse = TRUE), file = paste(data.dir, "matrix_saver_imputed.mtx", sep = "/"))
      write.table(rownames(E)[genes.ind], file = paste(data.dir, "genes_saver_imputed.txt", sep = "/"))
    }  
  }
    return(I)
  }
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
CID.CellID <- function(E,pval = 0.1,deep_dive = TRUE,spring.dir = NULL, entropy = FALSE, walktrap = FALSE, omit = NULL)
{
  # load markers
  data(markers)
  data(cellstate_markers)
  if (!length(markers) > 0) {
    cat("ERROR: from CID:\n");
    cat("required markers failed to load.\n", sep = "");
    stop()
  }
  temp = markers[markers$`Cell population` == "Platelets",]
  temp = temp[1:3,]
  markers = markers[markers$`Cell population` != "Platelets",]
  markers = rbind(markers, temp)

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
    filtered_features = filtered_features[-which(names(filtered_features) %in% omit)]
  }
  
  dfY = CID.append(E,filtered_features)
  
  # assign output classifications
  cat(" ..........  Assigning output classifications \n", sep ="");
  indexMax = apply(dfY, 2, which.max);
  ac = rownames(dfY)[indexMax];
  
  # compute distance matrix
  if (!is.null(spring.dir) | entropy | walktrap)
    distMat = CID.DistMatrix(ac, spring.dir)

  # smooth the output classifications
  if (!is.null(spring.dir))
  {
    f = 0.5
    cat(" ..........  Smoothing with smoothing parameter f = ", f, "\n", sep ="");
    acOut_knn_smooth = CID.smooth(ac, distMat, f = f) # smooth based on edges in knn graph
    cat("\n ..........  Smoothing completed! \n");
  }

  if (entropy)
  {
    acOut_knn_smooth = CID.entropy(acOut_knn_smooth, distMat)
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
        ac_dd = CID.deepdive(filtered_features, j, acOut3 = ac_dd, expression = E, f = f)
        dummy[[j]] = ac_dd;
      }
    }
    dummy = dummy[sapply(dummy, function(x) length(x) != 1)]
    names(dummy) <- paste("L", seq_along(1:length(dummy)), sep = "")
    cat(" ..........  Deep dive completed!\n");
    if (!is.null(spring.dir))
    {
      if (is.null(f))
        f = 0.5
      cat(" ..........  Smoothing with smoothing parameter f = ", f, "\n", sep ="");
      ac_dd_knn = CID.smooth(ac_dd, distMat, f)
    }
  }

  if (walktrap)
  {
    cat("\n")
    wt = CID.WalkTrap(acOut_knn_smooth, spring.dir)
    do = data.frame(table(wt[acOut_knn_smooth == "Other"]))
    df = data.frame(table(wt[wt %in% do$Var1]))
    logik = 1 - (df$Freq - do$Freq) / df$Freq > 0.9; # logical select all walktrap communites > 90 % labeled "other"
    if (sum(logik) > 0)
    {
      do = do[logik,]
      logik = do$Freq > 20; # require at least 20 cell communities
      if (sum(logik) > 0) 
      {
        lbls = rep("All", ncol(E))
        logik = wt %in% do$Var1[logik]
        lbls[logik] = wt[logik]
        colnames(E) <- lbls
        acOut_knn_smooth_wt = CID.PosMarkers(E, acOut_knn_smooth)
        ac_dd_knn_wt = ac_dd_knn
        logik = grepl("+", acOut_knn_smooth_wt)
        ac_dd_knn_wt[logik] = acOut_knn_smooth_wt[logik]
      }
    }
  }
  
  # Package output
  if (deep_dive & !is.null(spring.dir) & walktrap) {
    cr = list(scores = dfY,
              ctypes = ac,
              ctypessmoothed = acOut_knn_smooth,
              ctypessmoothed_dev = acOut_knn_smooth_wt,
              ddtypes = ac_dd,
              ddtypessmoothed = ac_dd_knn,
              ddtypessmoothed_dev = ac_dd_knn_wt,
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
CID.append=function(expression,featuresG)
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
}

CID.deepdive <- function(XX, YY, acOut3 = ac_dd, expression = E, f = f)
{
  cat(" ..........  Assigning cell states :\n");
  logik = acOut3 == YY; sum (logik)
  if (sum(logik) != 0 ){
    cat("             Hierarchy assignment: ", sum(logik), " ", YY, "\n", sep ="");
    ## get CID scores for cell types
    dfYe = CID.append(expression,XX)
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
CID.smooth <- function(ac,dM,f = 0.5)
{
  # pre-allocate
  Y   = ac # cells we may switch labels

  # smooth based on majority vote of direct connections
  pb = txtProgressBar(min = 0, max = length(Y), initial = 0, style = 3
                      ) 
  
  for (i in 1:length(Y))
  {
    logik = dM[i,] == 1;
    x = data.frame(table(as.character(ac[logik])))
    logik =  (x$Freq/sum(x$Freq) > f)
    y = as.character(x$Var1)
    dummy = y[logik]
    if (length(dummy) > 0)
      Y[i] = dummy[1]
    setTxtProgressBar(pb,i)
  }
  
  # assign the smoothed labels
  ac = Y
  
  return(ac)
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
  
  # Calculate shannon entropy for each cell j for all connections with shortest path < N
  shannon = rep(0, length(Y))
  N = 4
  for (j in 1:length(Y))  {
    logik = dM[j,] <= N; sum(logik)
    freqs <- table(ac[logik])/sum(logik)
    shannon[j] = -sum(freqs * log2(freqs))
  }
  #df = data.frame(cells = Y, shannon = shannon)
  #ggplot2::ggplot(df, ggplot2::aes(x=cells, y=shannon, color = cells)) + ggplot2::geom_boxplot()
  logik = shannon > 0.75; sum(logik)
  Y[logik] = "Other"
  
  # assign Other labels
  ac = Y
  
  return(ac)
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
  if ("ctypessmoothed" %in% names(cr))
  {
    Q = cr$ctypessmoothed
    json_data$CellTypesID$label_list = Q
    C = get_colors(Q)
    json_data$CellTypesID$label_colors = as.list(C[[1]])
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
  if (is.null(spring.dir))
  {
    fn = json_new;
    spring.dir = getwd()
    spring.dir = gsub("\\/$", "", spring.dir, perl = TRUE);
  }
  fn = json_new
  new.dirs = list.dirs(dirname(spring.dir))
  if (length(new.dirs) != 1)
  {
    for (j in 1:length(new.dirs))
    {
      txt = list.files(new.dirs[j])
      flag = "cell_filter.txt" %in% list.files(new.dirs);
      if (!flag) {
        cat("ERROR: from CID.writeJSON:\n");
        cat("cell_filter.txt does not exist anywhere inside",new.dirs, ".\n", sep = "");
        stop()
      }
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
        cat(paste(new.dirs[j],fn,sep="/"), "has been written to directory! \n")
        json_data = json_data_backup;
      }        
    }
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

#' Makes probability distribution plots for genes across celltypes
#'
#' @param E expression matrix (genes by rows) where genes are rownames and columns are celltype labels (see ?CID.CellID)
#' @param genes a character vector of genes to plot.
#' @return A probability distribution plot for the queried genes.
#' @export

CID.VinPlot <- function(E,genes = c("MS4A1", "CD40", "NKG7"), lbls = NULL, do.print = FALSE, omit = NULL, filename = NULL)
{
  if (is.null(colnames(E)))
    warning("Make sure that the column names of E are set properly! (See ?CID.VinPlot) \n")

  missing.genes = genes[! genes %in% rownames(E)]

  if (length(missing.genes) != 0)
    cat("WARNING: Genes ", paste(missing.genes, collapse = ", "), "are not present in E. Consider checking for gene aliases. \n")

  genes = genes[!genes %in% missing.genes]
  # D = log(as.matrix(E[rownames(E) %in% genes,])+1,2)
  D = as.matrix(E[rownames(E) %in% genes,])
  names(D) <- colnames(E)
  if (!is.null(omit))
  {
    logik = !colnames(D) %in% omit
    D = D[,logik]
    lbls = lbls[logik]
  }
  df = reshape2::melt(D)
  if (!is.null(lbls))
  {
    colnames(D) = lbls;
    df2 = reshape2::melt(D)
    df$lbls = as.character(df2$Var2)
  }
  
  if (is.null(lbls)) {
    df$lbls = colnames(D)
  }
  df$Var1<- as.character(df$Var1)
  df = df[order(df$Var2),]
  df$Var2 <- factor(df$Var2, levels=unique(df$Var2))
  colors = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7",
             "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
             "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D",
             "#8A7C64", "#599861")
  p = list("")
  for (j in 1:length(genes))
  {
    p[[j]]<-ggplot(df[df$Var1 == genes[j],], aes(x=Var2, y=value, fill=Var2)) +
      geom_boxplot(aes(fill = lbls)) +
      labs(title=genes[j],x="Celltypes") +
      scale_fill_manual(values=colors[1:length(unique(df$Var2))]) +
      scale_x_discrete(labels = abbreviate) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.05)) + labs(fill = "Celltypes") +
      ylab(expression(atop("Expression level", paste("log" [2], "(nUMI + 1)"))))
  }
  if (do.print)
  {
    pdf(filename)
    for (i in p){
      print(i)
    }
    dev.off()
    invisible(NULL)
  }
  ###### getting distinguishable colours for clusters #####
  ##ref: http://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
  ##

  #qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  #col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) #len = 74
  #num_col <- 40
  #pie(rep(1,num_col), col=(col_vector[1:num_col]))
  #col_palette <- col_vector[1:num_col]; # or sample if you wish

  return(p)
}

#' Community substructure by random walk
#'
#' @param ac A character vector of cell type labels
#' @param edges A data frame or matrix with edges
#' @return The community substructures of the graph
#' @export
CID.WalkTrap <- function(ac,edges)
{
  edges = CID.LoadEdges(edges)
  # convert edgelist to adjacency matrix
  adjmatrix <- Matrix::Matrix(0, length(ac), length(ac), sparse = T)
  adjmatrix[cbind(edges$V1, edges$V2)] <- 1
  # get distance matrix
  g = igraph::graph_from_adjacency_matrix(adjmatrix, mode = c("directed", "undirected",
                                                              "max", "min", "upper", "lower", "plus"), weighted = NULL, diag = TRUE,
                                          add.colnames = NULL, add.rownames = NA)
  as.character(igraph::cluster_walktrap(g)$membership)
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
  #cts2 = data.frame(ident.1 = sort(cts[grepl("Samples1-2", cts)]), ident.2 = sort(cts[grepl("Samples3-12", cts)]))
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

    #dd = Seurat::FindMarkers(ctrl, ident.1 = cts[j], ident.2 = NULL, min.cells.group = 0, max.cells.per.ident = 200, logfc.threshold = 0, pseudocount.use = 1, min.pct = 0)
    #dd = Seurat::FindMarkers(ctrl, ident.1 = cts2$ident.1[j], ident.2 = cts2$ident.2[j], min.cells.group = 0, max.cells.per.ident = 200, logfc.threshold = 0, pseudocount.use = 0, min.pct = 0)
  #  dd$GeneSymbol = rownames(dd)
  #  dd$celltype = gsub( " .*$", "", cts2$ident.1 )[j]
  #  outs[[j]] = dd
  }
  #dd = Seurat::FindMarkers(ctrl, ident.1 = "wt", ident.2 = "All", min.cells.group = 0, max.cells.per.ident = 200, logfc.threshold = 0, pseudocount.use = 1, min.pct = 0)
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

#' Visualizations with Seurat wrapper
#'
#' @param spring.dir directory where SPRING files are located
#' @return Visualizations from Seurat
#' @export
CID.VisFeatures <- function(spring.dir, features.plot = c("LILRB2"))
{
  E = CID.LoadData(spring.dir)
  # make sure row names are not redundant
  logik = CID.IsUnique(rownames(E))
  E = E[logik,]
  # get lbls
  json_data <- rjson::fromJSON(file=paste(spring.dir, "categorical_coloring_data.json", sep = ""))
  # set unique colnames
  colnames(E) <- seq(1, ncol(E))
  # Set up object
  ctrl <- Seurat::CreateSeuratObject(counts = E)
  p = list("")
  q = list("")
  r = list("")
  k = 1;
  if ("CellTypesID" %in% names(json_data) & "Disease" %in% names(json_data))
  {
    lbls = paste(json_data$CellTypesID$label_list, json_data$Disease$label_list)
    ctrl <- Seurat::AddMetaData(ctrl, metadata=lbls, col.name = 'CellTypesID')
    ctrl <- Seurat::SetIdent(ctrl, value='CellTypesID')
    p[[k]]<-Seurat::RidgePlot(object = ctrl, features = features.plot, ncol = length(features.plot))
    q[[k]]<-Seurat::VlnPlot(object = ctrl, features = features.plot)
    r[[k]]<-Seurat::DotPlot(object = ctrl, features = features.plot)
    k = k + 1
  }
  if ("CellStatesID" %in% names(json_data) & "Disease" %in% names(json_data))
  {
    lbls = paste(json_data$CellStatesID$label_list, json_data$Disease$label_list)
    ctrl <- Seurat::AddMetaData(ctrl, metadata=lbls, col.name = 'CellStatesID')
    ctrl <- Seurat::SetIdent(ctrl, value='CellStatesID')
    p[[k]]<-Seurat::RidgePlot(object = ctrl, features = features.plot, ncol = length(features.plot))
    q[[k]]<-Seurat::VlnPlot(object = ctrl, features = features.plot)
    r[[k]]<-Seurat::DotPlot(object = ctrl, features = features.plot)
    k = k + 1
  }
  if ("CellTypesID" %in% names(json_data))
  {
    lbls = json_data$CellTypesID$label_list
    ctrl <- Seurat::AddMetaData(ctrl, metadata=lbls, col.name = 'CellTypesID')
    ctrl <- Seurat::SetIdent(ctrl, value='CellTypesID')
    p[[k]]<-Seurat::RidgePlot(object = ctrl, features = features.plot, ncol = length(features.plot))
    q[[k]]<-Seurat::VlnPlot(object = ctrl, features = features.plot)
    r[[k]]<-Seurat::DotPlot(object = ctrl, features = features.plot)
    k = k + 1
  }
  if ("CellStatesID" %in% names(json_data))
  {
    lbls = json_data$CellStatesID$label_list
    ctrl <- Seurat::AddMetaData(ctrl, metadata=lbls, col.name = 'CellStatesID')
    ctrl <- Seurat::SetIdent(ctrl, value='CellStatesID')
    p[[k]]<-Seurat::RidgePlot(object = ctrl, features = features.plot, ncol = length(features.plot))
    q[[k]]<-Seurat::VlnPlot(object = ctrl, features = features.plot)
    r[[k]]<-Seurat::DotPlot(object = ctrl, features = features.plot)
    k = k + 1
  }
  if ("Disease" %in% names(json_data))
  {
    lbls = json_data$Disease$label_list
    ctrl <- Seurat::AddMetaData(ctrl, metadata=lbls, col.name = 'Disease')
    ctrl <- Seurat::SetIdent(ctrl, value='Disease')
    p[[k]]<-Seurat::RidgePlot(object = ctrl, features = features.plot, ncol = length(features.plot))
    q[[k]]<-Seurat::VlnPlot(object = ctrl, features = features.plot)
    r[[k]]<-Seurat::DotPlot(object = ctrl, features = features.plot)
    k = k + 1
  }
  if ("Tissue" %in% names(json_data))
  {
    lbls = json_data$Tissue$label_list
    ctrl <- Seurat::AddMetaData(ctrl, metadata=lbls, col.name = 'Tissue')
    ctrl <- Seurat::SetIdent(ctrl, value='Tissue')
    p[[k]]<-Seurat::RidgePlot(object = ctrl, features = features.plot, ncol = length(features.plot))
    q[[k]]<-Seurat::VlnPlot(object = ctrl, features = features.plot)
    r[[k]]<-Seurat::DotPlot(object = ctrl, features = features.plot)
    k = k + 1
  }
  if ("CellType" %in% names(json_data))
  {
    lbls = json_data$CellType$label_list
    ctrl <- Seurat::AddMetaData(ctrl, metadata=lbls, col.name = 'Celltype')
    ctrl <- Seurat::SetIdent(ctrl, value='Celltype')
    p[[k]]<-Seurat::RidgePlot(object = ctrl, features = features.plot, ncol = length(features.plot))
    q[[k]]<-Seurat::VlnPlot(object = ctrl, features = features.plot)
    r[[k]]<-Seurat::DotPlot(object = ctrl, features = features.plot)
    k = k + 1
  }
  
  filename = paste(spring.dir, "geneviews_ridge.pdf", sep = "")
  pdf(filename)
    for (i in p){
      print(i)
    }
    dev.off()
    invisible(NULL)
    
    filename = paste(spring.dir, "geneviews_violin.pdf", sep = "")
    pdf(filename)
    for (i in q){
      print(i)
    }
    dev.off()
    invisible(NULL)
    
    filename = paste(spring.dir, "geneviews_dotplot.pdf", sep = "")
    pdf(filename)
    for (i in r){
      print(i)
    }
    dev.off()
    invisible(NULL)
}

#' Community substructure by random walk
#'
#' @param ac A character vector of cell type labels
#' @param edges A data frame or matrix with edges
#' @return The community substructures of the graph
#' @export
CID.Scatter <- function(spring.dir, features.plot)
{
  E = CID.LoadData(spring.dir)
  json_data <- rjson::fromJSON(file=paste(spring.dir, "categorical_coloring_data.json", sep = ""))
  genes.missing = features.plot[!features.plot %in% rownames(E)]
  if (length(genes.missing) > 0)
  {
    cat("WARNING: Gene(s) missing: ", genes.missing, "\n")
    stop()
  }
  # reshape for ggplot
  df = data.frame(Matrix::t(E[rownames(E) %in% features.plot,]))
  colnames(df) <- features.plot
  df$Disease = json_data$Disease$label_list
  df$CellTypes = json_data$CellTypesID$label_list
  # Change point shapes, colors and sizes
  ggplot2::ggplot(df, ggplot2::aes(x=PTPRC, y=ALDOB, color = CellTypes, shape = Disease)) + ggplot2::geom_point()
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

