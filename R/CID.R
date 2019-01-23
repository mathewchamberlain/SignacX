#' Load data file from directory
#'
#' @param data.dir Directory containing matrix.mtx and genes.txt.
#' @return A sparse matrix with rownames equivalent to the names in genes.txt
#' @export
CID.LoadData <- function(spring.dir)
{
  data.dir = dirname(spring.dir)
  data.dir = gsub("\\/$", "", data.dir, perl = TRUE);
  fn <- "matrix.mtx"
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
    cat("length of genes in genes.txt = ", length(genes), " is not equal to nrow(E) = ", nrow(E), "or ncol(E) - ", ncol(E), "\n", sep = "");
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
    cat("ERROR: from CID.entropy:\n");
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
CID.Impute <- function(E, do.par = FALSE)
{
  # SAVER wrapper
  data(markers)
  data(cellstate_markers)
  genes = do.call(rbind, cellstate_markers)
  genes.ind <- which(rownames(E) %in% unique(c(as.character(markers$`HUGO symbols`), as.character(genes$`HUGO symbols`))))
  if (do.par)
  {
    numCores = parallel::detectCores()
    SAVER::saver(E, pred.genes = genes.ind, pred.genes.only = T, ncores = numCores / 4, estimates.only = T)
  } else {
    SAVER::saver(E, pred.genes = genes.ind, pred.genes.only = T, estimates.only = T)
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
#' @param f Parameter for classification smoothing. Default is f = NULL if edges = NULL, else default is f = 0.5 automatically.
#' @param pval p-value cutoff for feature selection, as described in the manuscript + markdown file. Default is pval = 0.1.
#' @param deep_dive Boolean, T will assign cell types and then cell states, F will only assign cell types. Default is deep_dive = T.
#' @param edges Data-frame containing the edgelist for graph edges for smoothing, OR a directory where "edges.csv" is located, OR a character vector / list for batch mode. Default is edges = NULL.
#' @return Filtered markers where each marker must have at least ncells that express at least ncounts
#' @export
CID.CellID <- function(E,f = NULL,pval = 0.1,deep_dive = TRUE,edges = NULL, entropy = FALSE, sorted = FALSE, walktrap = FALSE)
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
    return(CID.BatchMode(E = E,f = f,pval = pval,deep_dive = deep_dive,edges = edges, entropy = entropy, sorted = sorted, walktrap = walktrap))

  # check inputs
  cat(" ..........  Entry in CID.CellID \n");
  ta = proc.time()[3];
  stopifnot(class(E) %in% c("dgCMatrix","dgTMatrix", "matrix", "data.frame"))
  stopifnot(!is.null(rownames(E)));

  if (sorted)  {
    markers = markers[markers$Polarity == "+",]
    cellstate_markers = lapply(cellstate_markers, function(x) {x[x$Polarity == "+",]})
  }

  # main function
  cat(" ..........  Computing Signac scores for cell types on input data matrix :\n");
  cat("             nrow = ", nrow(E), "\n", sep = "");
  cat("             ncol = ", ncol(E), "\n", sep = "");
  filtered_features = CID.filter(E, markersG  = markers, pval = pval)
  dfY = CID.append(E,filtered_features)

  # assign output classifications
  cat(" ..........  Assigning output classifications \n", sep ="");
  indexMax = apply(dfY, 2, which.max);
  ac = rownames(dfY)[indexMax];

  # compute distance matrix
  if (!is.null(edges) | entropy | walktrap)
  distMat = CID.DistMatrix(ac, edges)

  # smooth the output classifications
  if (!is.null(edges))
  {
    if (is.null(f))
      f = 0.5
    cat(" ..........  Smoothing with smoothing parameter f = ", f, "\n", sep ="");
    acOut_knn_smooth = CID.smooth(ac, distMat, f = f) # smooth based on edges in knn graph
    cat(" ..........  Smoothing completed! \n");
  }

  if (entropy)
  {
    ac = CID.entropy(ac, distMat)
    acOut_knn_smooth = CID.entropy(acOut_knn_smooth, distMat)
  }

  # cell state deep dive classifications
  if (deep_dive)
  {
    dummy = list("")
    cat(" ..........  Computing CID scores for cell states! \n");
    if (!is.null(edges))
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
    if (!is.null(edges))
    {
      if (is.null(f))
        f = 0.5
      cat(" ..........  Smoothing with smoothing parameter f = ", f, "\n", sep ="");
      ac_dd_knn = CID.smooth(ac_dd, f = f, distMat)
    }
  }

  if (walktrap)
    wt = CID.WalkTrap(acOut_knn_smooth, edges)
  
  # Package output
  if (deep_dive & !is.null(edges) & walktrap) {
    cr = list(scores = dfY,
              ctypes = ac,
              ctypessmoothed = acOut_knn_smooth,
              ddtypes = ac_dd,
              ddtypessmoothed = ac_dd_knn,
              hierarchylabels = dummy,
              walktrap = wt)
  }
  else if (deep_dive & !is.null(edges)) {
    cr = list(scores = dfY,
              ctypes = ac,
              ctypessmoothed = acOut_knn_smooth,
              ddtypes = ac_dd,
              ddtypessmoothed = ac_dd_knn,
              hierarchylabels = dummy)
  }
  else if (!is.null(edges)) {
    cr = list(scores = dfY,
              ctypes = ac,
              ctypessmoothed = acOut_knn_smooth);}
  else if (deep_dive & is.null(edges)) {
    cr = list(scores = dfY,
              ctypes = ac ,
              ddtypes = ac_dd ,
              hierarchylabels = dummy)}
  else {
    cr = list(scores = dfY,
              ctypes = ac)
  }
  tb = proc.time()[3] - ta;
  cat(" ..........  Exit CID.CellID.\n");
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
CID.smooth <- function(ac,distMat,f = 0.5)
{
  # pre-allocate
  Y = ac
  # smooth based on majority vote of direct connections
  for (i in 1:length(ac))
  {
    logik = distMat[,i] == 1;
    x = data.frame(table(as.character(ac[logik])))
    logik =  (x$Freq/sum(x$Freq) > f)
    y = as.character(x$Var1)
    dummy = y[logik]
    if (length(dummy) > 0)
      Y[i] = dummy[1]
  }

  return(Y)
}

#' Distance matrix
#'
#' @param ac A character vector of cell type labels
#' @param edges A data frame or matrix with edges
#' @return The entropy for each node in the network
#' @export
CID.DistMatrix <- function(ac,edges)
{
  # load edges
  edges = CID.LoadEdges(edges)
  # convert edgelist to adjacency matrix
  adjmatrix <- Matrix::Matrix(0, length(ac), length(ac), sparse = T)
  adjmatrix[cbind(edges$V1, edges$V2)] <- 1
  # get distance matrix
  g = igraph::graph_from_adjacency_matrix(adjmatrix, mode = c("directed", "undirected",
                                                              "max", "min", "upper", "lower", "plus"), weighted = NULL, diag = TRUE,
                                          add.colnames = NULL, add.rownames = NA)
  igraph::shortest.paths(g, v=igraph::V(g), to=igraph::V(g))
}
#' Entropy
#'
#' @param edges A data frame or matrix with edges
#' @return The entropy for each node in the network
#' @export
CID.entropy <- function(ac,distMat)
{
  # Calculate shannon entropy for each cell j for all connections with shortest path < N
  shannon = rep(0, length(ac))
  N = 4
  for (j in 1:length(ac))  {
    logik = distMat[,j] <= N; sum(logik)
    freqs <- table(ac[logik])/sum(logik)
    shannon[j] = -sum(freqs * log2(freqs))
  }
  #df = data.frame(cells = ac, shannon = shannon)
  #ggplot(df, aes(x=cells, y=shannon, color = cells)) + geom_boxplot()
  logik = shannon > 0.5; sum(logik)
  ac[logik] = "Other"
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
    json_data$Clusters$label_list = Q
    Ntypes = length(unique(Q))
    qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) #len = 74
    #pie(rep(1,num_col), col=(col_vector[1:num_col]))
    col_palette <- as.list(col_vector[1:Ntypes]); # or sample if you wish
    names(col_palette) <- unique(Q)
    json_data$Clusters$label_colors = col_palette
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
  for (i in 1:length(main_types))
  {
    cs[P == main_types[i]] = main_colors[i]
    cs[P == paste(main_types[i],"-like",sep ="")] = main_colors[i]
  }

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
                 "#4e59bd"                    ,  "#d4e881"                     , "#90c5f4"          ,
                 "#e2e8c9"                    ,  "#def9f9"                     , "#aa3596"          ,
                 "#edc5e6"                    ,  "#4e59bd"                     , "#8c9ee1"          ,
                 "#90c5f4"                    ,  "#90c5f4"                     , "#D3D3D3"          ,
                 "#d4e881"                    ,  "#d4e881"                     , "#f7f2b2"          ,
                 "#F9A602"                    ,  "#f93a01"                     , "#d4e881"          ,
                 "#ad9bf2")

  for (i in 1:length(sub_types))
  {
    cs[P == sub_types[i]] = sub_colors[i]
  }

  cs[P == "None"] = "#808080"

  colfunc <- colorRampPalette(sub_colors)

  coms = apply(expand.grid(sub_types, sub_types), 1, paste, collapse="-")
  col=colfunc(length(coms))
  for (i in 1:length(coms))
  {
    cs[P == coms[i]] = col[i]
  }
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

#' Community substructure by random walk
#'
#' @param ac A character vector of cell type labels
#' @param edges A data frame or matrix with edges
#' @return The community substructures of the graph
#' @export
CID.ViewHierarchy <- function()
{

df = data.frame(root = c(rep("Non-immune", 4), rep("Immune", 23)),
                leaf1 = c("Epithelial", "Epithelial.Urinary", "Fibroblasts", "Secretory", rep("Lymphocytes",10), rep("Myeloid",13)),
                leaf2 = c("Epithelial", "Epithelial.Urinary", "Fibroblasts", "Secretory", c(rep("TNK",7), rep("B.cells",2), "Plasma.Cells", rep("MPH",6), rep("GN",4), "pDCs", "Platelets","Erythro")),
                leaf3 = c("Epithelial", "Epithelial.Urinary", "Fibroblasts", "Secretory", c(rep("T",6),          "NK", "Memory", "Naive", "Plasma.Cells", rep("DC",2),            rep("Monocytes",  4),              rep("Mast",2),                              rep("Not.Mast",2)             , "pDCs", "Platelets","Erythro")),
                leaf4 = c("Epithelial", "Epithelial.Urinary", "Fibroblasts", "Secretory", c(rep("CD4",5),         "CD8", "NK", "Memory", "Naive", "Plasma.Cells", "activated", "resting", rep("Macrophages",3), "Monocytes", "Mast.cells.resting", "Mast.cells.activated", "Neutrophils", "Eosinophils", "pDCs", "Platelets","Erythro")),
                leaf5 = c("Epithelial", "Epithelial.Urinary", "Fibroblasts", "Secretory", c(rep("CD4",4), "Regs", "CD8", "NK", "Memory", "Naive", "Plasma.Cells", "activated", "resting", "M0", "M1", "M2"    , "Monocytes", "Mast.cells.resting", "Mast.cells.activated", "Neutrophils", "Eosinophils", "pDCs", "Platelets","Erythro")),
                leaf6 = c("Epithelial", "Epithelial.Urinary", "Fibroblasts", "Secretory", "naive", "memory.resting", "memory.activated", "follicular.helpere", "Regs", "CD8", "NK", "Memory", "Naive", "Plasma.Cells", "activated", "resting", "M0", "M1", "M2"    , "Monocytes", "Mast.cells.resting", "Mast.cells.activated", "Neutrophils", "Eosinophils", "pDCs", "Platelets","Erythro")
                )

setDT(df)

## Loading packages
library("data.table")
library("D3partitionR")

##Agregating data to have unique sequence for the 4 variables
var_names=c(" ","leaf1", "leaf2", "leaf3", "leaf4", "leaf5", "leaf6")
data_plot=df[,.N,by=var_names]
data_plot[,(var_names):=lapply(var_names,function(x){data_plot[[x]]=paste0(x,' ',data_plot[[x]])
})]

## Plotting the chart
library("magrittr")
D3partitionR() %>%
  add_data(data_plot,count = 'N',steps=c("root","leaf1", "leaf2", "leaf3", "leaf4", "leaf5", "leaf6")) %>%
  add_title('Titanic') %>%
  plot()
}
