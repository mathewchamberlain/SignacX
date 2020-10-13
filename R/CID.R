#' Generate a network with Graphical LASSO
#'
#' @param E Sparse expression matrix with rows genes and columns cells (see ?Signac)
#' @param gns_of_int <10,000 gene list of genes for network reconstruction
#' @param rho sparsity parameter for graphical LASSO
#' @param metadata data frame with metadata for each cell after signac classification
#' @return network object
#' @export
GenerateNetwork = function(E, gns_of_int, rho = c(0.3, 0.4, 0.5, 0.6), metadata)
{
  # subset expression matrix
  E = E[rownames(E) %in% gns_of_int,]
  
  # build networks around rho
  res = lapply(rho, function(z){
    
    ## seperate two matrices by cell state
    Q = lapply(unique(metadata$CellTypes), function(x) {
      Matrix::t(E[,metadata$CellTypes == x])
    })
    
    ## remove any genes with zero expression
    Q = lapply(Q, function(x) {x[,Matrix::colSums(x) != 0]})
    
    ## remove any cells with zero expression
    Q = lapply(Q, function(x) {x[Matrix::rowSums(x) != 0,]})
    
    # Graphical LASSO using QUIC
    Covs = lapply(Q, function(x){
      cov(as.matrix(x))
    })
    
    Nets = lapply(Covs, function(x){
      glasso::glasso(x, rho = z)
    })
    
    # apply rownames
    Nets = lapply(Nets, function(x) {x$wi})
    Nets = mapply(function(x,y) {
      rownames(x) <- colnames(y)
      colnames(x) <- colnames(y)
      x
    }, x = Nets, y = Q)
    
    names(Nets) <-unique(metadata$CellTypes)
    Nets
  })
  
  # keep stable networks
  N = lapply(res, function(x) {
    lapply(x, function(z){
      1*(z != 0)})
    })

  outs = list("")
  for (j in 1:length(res[[1]])){
    outs[[j]] = Reduce("+", lapply(N, function(x){x[[j]]}))
    logik = outs[[j]] == length(rho)
    outs[[j]][logik] = 1
    outs[[j]][!logik] = 0
    outs[[j]] = Matrix::Matrix(outs[[j]], sparse = T)
  }
  
  names(outs) <- names(N[[1]])
  
  # return networks
  return(outs)
}

#' Generate a probability plot
#'
#' @param df one element of the list returned by Signac function call (see ?Signac)
#' @return ggplot object
#' @export
ProbabilityPlot = function(df)
{
  ggplot(df, aes(x=genes_detected, y=probability, color=celltypes)) + geom_point() + geom_errorbar(aes(ymin=probability-error, ymax=probability+error), width=.2) + xlab("Genes detected") + ylab("Probability")
}

#' Get genes from xCell publication
#'
#' @param dataset default is NULL
#' @return a list of gene signatures
#' @export
getxCellGenes <- function(dataset = NULL)
{
  ## load signatures from xCell
  # devtools::install_github('dviraran/xCell')
  data("xCell.data", package = 'xCell')
  # get gene IDs
  signatures = lapply(xCell.data$signatures@.Data, function(x) {x@geneIds})
  
  # get celltypes
  celltypes = lapply(xCell.data$signatures@.Data, function(x) {sub("\\%.*", "", x@setName)})
  
  # get technologies
  datasets = lapply(xCell.data$signatures@.Data, function(x) {gsub(".*[%]([^.]+)[%].*", "\\1", x@setName)})
  
  if (!is.null(dataset)){
    logik = sapply(datasets, function(x) {x %in% dataset})
    celltypes = unlist(celltypes[logik])
    signatures = signatures[logik]
  } else {
    celltypes = unlist(celltypes)
  }
  
  names(signatures) = celltypes 
  
  return(signatures)   
}

#' Differential Expression Analysis for reference dataset
#'
#' @param R A reference list as returned by data("Reference_sim")
#' @return a list of DEG tables
#' @export
GetMarkers <- function(R)
{
  E = R$data
  # set colnames
  colnames(E) <- seq(1, ncol(E))
  # Set up object
  outs = lapply(1:ncol(R$celltypes),function(x){
    lbls = R$celltypes[,x]
    if (length(unique(lbls)) == 1)
    {
      return(NULL)
    } else {
      if (x > 1)
      {
        logik = R$celltypes[,x] != R$celltypes[,x-1]
        dat = E[,logik]
        lbls = lbls[logik]
      } else {
        dat = E
        lbls = lbls
      }
      ctrl <- Seurat::CreateSeuratObject(counts = dat, project = "CID", min.cells = 0)
      ctrl <- Seurat::NormalizeData(ctrl, normalization.method = "RC")
      ctrl <- Seurat::AddMetaData(ctrl, metadata=lbls, col.name = "celltypes")
      ctrl <- Seurat::SetIdent(ctrl, value='celltypes')
      Seurat::FindAllMarkers(ctrl, only.pos = T, min.cells.group = 0)
    }
  })
  
  return(outs)
}

#' Jackknife Differential Expression Analysis
#'
#' @param E A gene by cell expression matrix
#' @param Samples A character vector of samples.
#' @param Identities A character vector of cell type or cluster labels
#' @param Disease A character vector of disease labels
#' @param num.cores optionally, user can hard set the number of cores to use
#' @param omit Any samples to omit. Default is "none" and "empty"
#' @return a DEG table Jackknifed
#' @export
JackD <- function(E, Samples = NULL, Disease = NULL, Identities, num.cores = 1, omit = c("none", "Empty"))
{
  cat(" ..........  Entry in JackD: \n");
  ta = proc.time()[3];
  cat("             Number of cells:", ncol(E), "\n");
  cat("             Number of genes:", nrow(E), "\n");
  
  colnames(E) <- 1:ncol(E)
  
  if (! is.null(Disease))
  {
  logik = !Samples %in% omit
  E = E[,logik]
  Disease = Disease[logik]
  Identities = Identities[logik]
  Samples = Samples[logik]
  u = as.character(unique(Samples))
  cat("             Number of samples:", length(u), "\n");
  
  outs = parallel::mclapply(u, function(x){
    # dropout one sample
    if (length(u) > 1)
    {
      logik = Samples != x;
    } else {
      logik = Samples == x;
    }
    E_dummy          = E[,logik]
    Disease_dummy    = Disease[logik]
    Identities_dummy = Identities[logik]
    # create Seurat object for each cell type; run differential expression
    y = as.character(unique(Identities_dummy))
    qq = lapply(y, function(z){
      logik = Identities_dummy == z;
      E.          = E_dummy[,logik]
      Disease.    = Disease_dummy[logik]
      Identities. = Identities_dummy[logik]
      ctrl <- Seurat::CreateSeuratObject(counts = E.)
      ctrl <- Seurat::NormalizeData(object = ctrl)
      ctrl <- Seurat::AddMetaData(ctrl, metadata=Disease., col.name = "Disease")
      ctrl <- Seurat::SetIdent(ctrl, value='Disease')
      Seurat::FindAllMarkers(ctrl, only.pos = T)
    })
    names(qq) <- y
    qq
  }, mc.cores =  num.cores)
  
  ## pool together DE results
  res = list("")
  for (j in 1:length(unique(Identities)))
  {
    xx = unlist(lapply(outs, function(x) {
      x[[j]]$gene
    }))
    df = data.frame(table(xx))
    logik = df$Freq == max(df$Freq)
    gns_stable = sort(unique(as.character(df$xx)[logik]))
    xx = lapply(outs,function(x){
      x[[j]][x[[j]]$gene %in% gns_stable,]
    })
    xx = xx[sapply(xx, function(x) {nrow(x) == length(gns_stable)})]
    xx = lapply(xx, function(x) {x[order(x$gene),]})
  ## merge the numeric
    df = data.frame(
      p_val =  Reduce(lapply(xx, function(x) {x$p_val}), f = '+') / length(xx),
      avg_logFC =  Reduce(lapply(xx, function(x) {x$avg_logFC}), f = '+') / length(xx),
      pct.1 = Reduce(lapply(xx, function(x) {x$pct.1}), f = '+') / length(xx),
      pct.2 = Reduce(lapply(xx, function(x) {x$pct.2}), f = '+') / length(xx),
      p_val_adj = Reduce(lapply(xx, function(x) {x$p_val_adj}), f = '+') / length(xx),
      disease = xx[[1]]$cluster,
      gene = xx[[1]]$gene
    )
    df = df[order(df$disease),]
    df$celltype = unique(Identities)[j]
    df$pct.1 = round(df$pct.1, digits = 2) * 100
    df$pct.2 = round(df$pct.2, digits = 2) * 100
    df$p_val = round(df$p_val, digits = 4)
    df$avg_logFC = round(df$avg_logFC, digits = 2)
    df$p_val_adj = round(df$p_val_adj, digits = 4)
    res[[j]] = df
  }
  res = do.call(rbind, res)
  return(res)
  } else {
    logik = !Samples %in% omit
    E = E[,logik]
    Identities = Identities[logik]
    Samples = Samples[logik]
    u = as.character(unique(Samples))
    cat("             Number of samples:", length(u), "\n");
    
    outs = parallel::mclapply(u, function(x){
      # dropout one sample
      if (length(u) > 1)
      {
        logik = Samples != x;
      } else {
        logik = Samples == x;
      }
      E_dummy          = E[,logik]
      Identities_dummy = Identities[logik]
      ctrl <- Seurat::CreateSeuratObject(counts = E_dummy)
      ctrl <- Seurat::NormalizeData(object = ctrl)
      ctrl <- Seurat::AddMetaData(ctrl, metadata=Identities_dummy, col.name = "Identities")
      ctrl <- Seurat::SetIdent(ctrl, value='Identities')
      Seurat::FindAllMarkers(ctrl, only.pos = T)
    }, mc.cores =  num.cores)
    
    outs = lapply(outs, function(x) {split.data.frame(x, f = x$cluster )})
    outs = lapply(outs, function(x) {x[order(names(x))]})
    
    ## pool together DE results
    res = list("")
    for (j in 1:length(unique(Identities)))
    {
      xx = unlist(lapply(outs, function(x) {
        x[[j]]$gene
      }))
      df = data.frame(table(xx))
      logik = df$Freq == max(df$Freq)
      gns_stable = sort(unique(as.character(df$xx)[logik]))
      xx = lapply(outs,function(x){
        x[[j]][x[[j]]$gene %in% gns_stable,]
      })
      xx = xx[sapply(xx, function(x) {nrow(x) == length(gns_stable)})]
      xx = lapply(xx, function(x) {x[order(x$gene),]})
      ## merge the numeric
      df = data.frame(
        p_val =  Reduce(lapply(xx, function(x) {x$p_val}), f = '+') / length(xx),
        avg_logFC =  Reduce(lapply(xx, function(x) {x$avg_logFC}), f = '+') / length(xx),
        pct.1 = Reduce(lapply(xx, function(x) {x$pct.1}), f = '+') / length(xx),
        pct.2 = Reduce(lapply(xx, function(x) {x$pct.2}), f = '+') / length(xx),
        p_val_adj = Reduce(lapply(xx, function(x) {x$p_val_adj}), f = '+') / length(xx),
        celltype = xx[[1]]$cluster,
        gene = xx[[1]]$gene
      )
      df$pct.1 = round(df$pct.1, digits = 2) * 100
      df$pct.2 = round(df$pct.2, digits = 2) * 100
      df$p_val = round(df$p_val, digits = 6)
      df$avg_logFC = round(df$avg_logFC, digits = 2)
      df$p_val_adj = round(df$p_val_adj, digits = 2)
      df = df[order(df$avg_logFC, decreasing = T),]
      res[[j]] = df
    }
    res = do.call(rbind, res)
    return(res)
  }
}

#' Differential Expression Analysis
#'
#' @param E A gene by cell expression matrix
#' @param Samples A character vector of samples.
#' @param Identities A character vector of cell type or cluster labels
#' @param Disease A character vector of disease labels
#' @param omit Any samples to omit. Default is "none" and "Empty"
#' @return a DEG table
#' @export
GetAllMarkers <- function(E, Samples = NULL, Disease = NULL, Identities, omit = c("none", "Empty", "Double200-0109&200-0611"))
{
  cat(" ..........  Entry in GetAllMarkers: \n");
  ta = proc.time()[3];
  cat("             Number of cells:", ncol(E), "\n");
  cat("             Number of genes:", nrow(E), "\n");
  
  colnames(E) <- 1:ncol(E)
  
  if (! is.null(Disease))
  {
    if (!is.null(Samples)) {
      logik = !Samples %in% omit
      E = E[,logik]
      Disease = Disease[logik]
      Identities = Identities[logik]
      Samples = Samples[logik]
    }

      # create Seurat object for each cell type; run differential expression
      y = as.character(sort(unique(Identities)))
      qq = lapply(y, function(z){
        logik = Identities == z;
        E.          = E[,logik]
        Disease.    = Disease[logik]
        Identities. = Identities[logik]
        if (length(unique(Disease.))!=1) { 
          ctrl <- suppressWarnings(Seurat::CreateSeuratObject(E.))
          ctrl <- Seurat::NormalizeData(object = ctrl, verbose = F)
          ctrl <- Seurat::AddMetaData(ctrl, metadata=Disease., col.name = "Disease")
          ctrl <- Seurat::SetIdent(ctrl, value='Disease')
          dummy = Seurat::FindAllMarkers(ctrl, only.pos = T, verbose = F, min.cells.group = 0)
          if (nrow(dummy) != 0)
          {
            dummy$celltype = z
            dummy
          } else {
            NULL
          }
        } else {
          NULL
          }
      })
      names(qq) <- y
      
      qq = qq[sapply(qq, function(x){!is.null(x)})]
    
    ## pool together DE results
    df = do.call(rbind, qq)
    names(df) <- c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "disease", "gene", "identity")
    df = df[order(df$disease),]
    df$pct.1 = round(df$pct.1, digits = 2) * 100
    df$pct.2 = round(df$pct.2, digits = 2) * 100
    df$avg_logFC = round(df$avg_logFC, digits = 2)
    df$p_val = formatC(df$p_val, format = "e", digits = 2)
    df$p_val_adj = formatC(df$p_val_adj, format = "e", digits = 2)
    rownames(df) <- NULL
    df = data.frame(
      gene = df$gene,
      identity = df$identity,
      disease = df$disease,
      avg_logFC = df$avg_logFC,
      pct.1 = df$pct.1,
      pct.2 = df$pct.2,
      p_val = df$p_val,
      p_val_adj = df$p_val_adj
    )
    df = split.data.frame(df, f = df$identity)
    df = lapply(df, function(x){x[order(x$avg_logFC, decreasing = T),]})
    df = do.call(rbind, df)
    rownames(df) <- NULL
    return(df) 
  } else {
    if (!is.null(Samples)) {
      logik = !Samples %in% omit
      E = E[,logik]
      Identities = Identities[logik]
    }
    
    # create Seurat object for each cell type; run differential expression
    ctrl <- suppressWarnings(Seurat::CreateSeuratObject(E))
    ctrl <- Seurat::NormalizeData(object = ctrl, verbose = F)
    ctrl <- Seurat::AddMetaData(ctrl, metadata=Identities, col.name = "Identities")
    ctrl <- Seurat::SetIdent(ctrl, value='Identities')
    df = Seurat::FindAllMarkers(ctrl, only.pos = T, verbose = F, min.cells.group = 0)
    
    ## pool together DE results
    names(df) <- c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "identity", "gene")
    df = df[order(df$identity),]
    df$pct.1 = round(df$pct.1, digits = 2) * 100
    df$pct.2 = round(df$pct.2, digits = 2) * 100
    df$avg_logFC = round(df$avg_logFC, digits = 2)
    df$p_val = formatC(df$p_val, format = "e", digits = 2)
    df$p_val_adj = formatC(df$p_val_adj, format = "e", digits = 2)
    rownames(df) <- NULL
    df = data.frame(
      gene = df$gene,
      identity = df$identity,
      avg_logFC = df$avg_logFC,
      pct.1 = df$pct.1,
      pct.2 = df$pct.2,
      p_val = df$p_val,
      p_val_adj = df$p_val_adj
    )
    df = split.data.frame(df, f = df$identity)
    df = lapply(df, function(x){x[order(x$avg_logFC, decreasing = T),]})
    df = do.call(rbind, df)
    rownames(df) <- NULL
    return(df) 
  }
  
  tb = proc.time()[3] - ta;
  cat("\n ..........  Exit GetAllMarkers.\n");
  cat("             Execution time = ", tb, " s.\n", sep = "");
}


#' Main function for mixed effect modeling
#'
#' @param dataset data frame of covariate, cell type, clustering or disease information
#' @param cluster celltypes returned by Signac or cluster identities
#' @param contrast Typically disease
#' @param random_effects User specified random effect variables in dataset
#' @param fixed_effects User specific fixed effects in dataset
#' @param verbose If TRUE, algorithm reports outputs
#' @return mixed effect model results
#' @export
MASC <- function(dataset, cluster, contrast, random_effects = NULL, fixed_effects = NULL,
                 verbose = FALSE) {
  # Check inputs
  if (is.factor(dataset[[contrast]]) == FALSE) {
    stop("Specified contrast term is not coded as a factor in dataset")
  }
  
  cat(" ..........  Entry in MASC \n");
  ta = proc.time()[3];
  
  # Generate design matrix from cluster assignments
  cluster <- as.character(cluster)
  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset <- cbind(designmat, dataset)
  
  # Convert cluster assignments to string
  cluster <- as.character(cluster)
  # Prepend design matrix generated from cluster assignments
  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset <- cbind(designmat, dataset)
  # Create output list to hold results
  res <- vector(mode = "list", length = length(unique(cluster)))
  names(res) <- attributes(designmat)$dimnames[[2]]
  
  # Create model formulas
  if (!is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0(c(paste0(fixed_effects, collapse = " + "),
                          paste0("(1|", random_effects, ")", collapse = " + ")),
                        collapse = " + ")
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
    }
  } else if (!is.null(fixed_effects) && is.null(random_effects)) {
    model_rhs <- paste0(fixed_effects, collapse = " + ")
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
      # For now, do not allow models without mixed effects terms
      stop("No random effects specified")
    }
  } else if (is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0("(1|", random_effects, ")", collapse = " + ")
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
    }
  } else {
    model_rhs <- "1" # only includes intercept
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
      stop("No random or fixed effects specified")
    }
  }
  
  # Initialize list to store model objects for each cluster
  cluster_models <- vector(mode = "list",
                           length = length(attributes(designmat)$dimnames[[2]]))
  names(cluster_models) <- attributes(designmat)$dimnames[[2]]
  
  # Run nested mixed-effects models for each cluster
  for (i in seq_along(attributes(designmat)$dimnames[[2]])) {
    test_cluster <- attributes(designmat)$dimnames[[2]][i]
    if (verbose == TRUE) {
      message(paste("Creating logistic mixed models for", test_cluster))
    }
    null_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ 1 + "),
                                   model_rhs), collapse = ""))
    full_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ ", contrast, " + "),
                                   model_rhs), collapse = ""))
    # Run null and full mixed-effects models
    null_model <- lme4::glmer(formula = null_fm, data = dataset,
                              family = binomial, nAGQ = 1, verbose = 0,
                              control = lme4::glmerControl(optimizer = "bobyqa"))
    full_model <- lme4::glmer(formula = full_fm, data = dataset,
                              family = binomial, nAGQ = 1, verbose = 0,
                              control = lme4::glmerControl(optimizer = "bobyqa"))
    model_lrt <- anova(null_model, full_model)
    # calculate confidence intervals for contrast term beta
    contrast_lvl2 <- paste0(contrast, levels(dataset[[contrast]])[2])
    contrast_ci <- lme4::confint.merMod(full_model, method = "Wald",
                                  parm = contrast_lvl2)
    # Save model objects to list
    cluster_models[[i]]$null_model <- null_model
    cluster_models[[i]]$full_model <- full_model
    cluster_models[[i]]$model_lrt <- model_lrt
    cluster_models[[i]]$confint <- contrast_ci
  }
  
  # Organize results into output dataframe
  output <- data.frame(cluster = attributes(designmat)$dimnames[[2]],
                       size = colSums(designmat))
  output$model.pvalue <- sapply(cluster_models, function(x) x$model_lrt[["Pr(>Chisq)"]][2])
  output[[paste(contrast_lvl2, "OR", sep = ".")]] <- sapply(cluster_models, function(x) exp(fixef(x$full)[[contrast_lvl2]]))
  output[[paste(contrast_lvl2, "OR", "95pct.ci.lower", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "2.5 %"]))
  output[[paste(contrast_lvl2, "OR", "95pct.ci.upper", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "97.5 %"]))
  
  tb = proc.time()[3] - ta;
  cat("\n ..........  Exit MASC.\n");
  cat("             Execution time = ", tb, " s.\n", sep = "");
  
  # Return MASC results
    return(output)
}

#' Main function for solo-decision classification
#'
#' @param E a sparse gene (rows) by cell (column) matrix, or a Seurat object. Rows are HUGO symbols.
#' @param R Reference dataset; user should do data("Reference_sim") and then set R to Refernence_sim.
#' @param spring.dir If using SPRING, directory to categorical_coloring_data.json. Default is NULL.
#' @param model.use Machine learning model to use. Default option is neural network. Can also be set to 'svm' or 'rf'.
#' @param N Number of machine learning models to train (for nn and svm). Default is 25.
#' @param num.cores Number of cores to use. Default is 1.
#' @param threshold Probability threshold for assigning cells to "Unclassified." Default is 0.5.
#' @param smooth if TRUE, smooths the cell type classifications. Default is TRUE.
#' @param impute if TRUE, gene expression values are imputed prior to cell type classification. Default is TRUE.
#' @param verbose if TRUE, code will report outputs. Default is TRUE.
#' @param do.normalize if TRUE, cells are normalized to the mean library size. Default is TRUE.
#' @param probability if TRUE, returns the probability associated with each cell type label. Default is TRUE.
#' @param hidden Number of hidden layers in the neural network. Default is 1.
#' @return annotations
#' @export
Signac_Solo <- function(E, R , spring.dir = NULL, model.use = "nn", N = 25, num.cores = 1, threshold = 0.5, smooth = T, impute = T, verbose = T, do.normalize = T, probability = F, hidden = 1)
{
  
  flag = class(E) == "Seurat"
  
  if (flag & impute)
    edges = E@graphs$RNA_nn
  
  if (verbose)
  {
    cat(" ..........  Entry in Signac_Solo \n");
    ta = proc.time()[3];
    
    # main function
    if (!flag)
    {
      cat(" ..........  Running Signac_Solo on input data matrix :\n");
    } else {
      cat(" ..........  Running Signac_Solo on Seurat object :\n");
    }
    cat("             nrow = ", nrow(E), "\n", sep = "");
    cat("             ncol = ", ncol(E), "\n", sep = "");
  }
  
  # keep only unique row names
  logik = CID.IsUnique(rownames(E))
  E = E[logik,]
  
  # intersect genes with reference set
  gns = intersect(rownames(E), colnames(R))
  V = E[rownames(E) %in% gns, ]
  
  # make sure data are in the same order
  V = V[order(rownames(V)),]
  
  if (class(V) %in% "data.frame")
    V = Matrix::Matrix(as.matrix(V), sparse = T)
  
  # normalize to the mean library size
  if (do.normalize)
  {
    if (!flag)
    {
      V = CID.Normalize(V)
    } else {
      V = CID.Normalize(V@assays$RNA@counts)
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
  
    # keep same gene names
    gns = sort(intersect(rownames(V), colnames(R)))
    Z = V[rownames(V) %in% gns, ]
    dat = R[,colnames(R) %in% gns]
    Z = Z[order(rownames(Z)), ]
    dat = dat[, order(colnames(dat))]
    
    # remove any low variance genes
    kmu = apply(Z, 1, function(x){sum(x != 0)})
    logik = kmu > 0;
    Z = Z[logik,]
    dat = dat[,logik]
    
    # run imputation (if desired)
    if (impute){
      Z = KSoftImpute(E = Z, dM = dM, verbose = F)
      Z = t(apply(Z, 1, function(x){
        normalize(x)
      }))
    }
    
    # build training set
    df = data.frame(dat, celltypes = R$celltypes)
    
    # train a neural network (N times)
    if (model.use == "nn"){
      res = parallel::mclapply(1:N, function(x) {
        nn=neuralnet::neuralnet(celltypes~.,hidden=hidden,data=df, act.fct = 'logistic', linear.output = F)
        Predict = stats::predict(nn, Matrix::t(Z))
        colnames(Predict) <- sort(nn$model.list$response)
        return(Predict)
      }, mc.cores = num.cores)
      
      # compute standard deviation of predictions
      err = do.call(cbind, lapply(res, function(x) {x[,1]}))
      res.sd <- apply(err, 1, sd)
      err = do.call(cbind, lapply(res, function(x) {x[,2]}))
      res.sd <- cbind(res.sd, apply(err, 1, sd))
      
      ## compute average probability
      res = Reduce(res, f = '+') / N
      xx = apply(res, 1, which.max)
      celltypes = colnames(res)[xx]
      kmax = apply(res, 1, max)
      celltypes[kmax < threshold] = "Other"
      df = data.frame(celltypes = celltypes, probability = kmax, error = res.sd[xx], percent_features_detected = round(Matrix::colSums(Z != 0) / nrow(Z), digits = 3) * 100, genes_detected = round(Matrix::colSums(E != 0), digits = 3) * 100)
      df$celltypes = as.character(df$celltypes)
    }
    # if desired, run svm
    if (model.use == "svm") {
      model = e1071::svm(celltypes ~ ., data = df, probability=TRUE)
      xx = stats::predict(model, t(Z), probability=TRUE)
      res = attr(xx, "probabilities")
      xx = apply(res, 1, which.max)
      celltypes = colnames(res)[xx]
      kmax = apply(res, 1, max)
      celltypes[kmax < threshold] = "Other"
      df = data.frame(celltypes = celltypes, probability = kmax)
      df$celltypes = as.character(df$celltypes)
    }
    # if desired, run RF
    if (model.use == "rf") {
      res = parallel::mclapply(1:N, function(x) {
        model = randomForest::randomForest(celltypes ~ ., data=df, keep.forest = TRUE)
        stats::predict(model,newdata=Matrix::t(Z),type="prob")
      }, mc.cores = num.cores)
      res = Reduce(res, f = '+') / N
      xx = apply(res, 1, which.max)
      celltypes = colnames(res)[xx]
      kmax = apply(res, 1, max)
      celltypes[kmax < threshold] = "Other"
      df = data.frame(celltypes = celltypes, probability = kmax)
      df$celltypes = as.character(df$celltypes)
    }
    # smooth the output classifications
    if (smooth)
      df$celltypes = CID.smooth(df$celltypes, dM[[1]])

  
  if (verbose) {
    tb = proc.time()[3] - ta;
    cat("\n ..........  Exit Signac_Solo. \n");
    cat("             Execution time = ", tb, " s.\n", sep = "");
  }
    
    # return probabilities and cell type classifications
    if (probability){
      return(df)
    } else {
      return(df$celltypes)
    }
    
}
#' Main function for classification
#'
#' @param E a sparse gene (rows) by cell (column) matrix, or a Seurat object. Rows are HUGO symbols.
#' @param R Reference dataset; user should do data("Reference_sim") and then set R to Refernence_sim.
#' @param spring.dir If using SPRING, directory to categorical_coloring_data.json. Default is NULL.
#' @param model.use Machine learning model to use. Default option is neural network. Can also be set to 'svm' or 'rf'.
#' @param N Number of machine learning models to train (for nn and svm). Default is 25.
#' @param num.cores Number of cores to use. Default is 1.
#' @param threshold Probability threshold for assigning cells to "Unclassified." Default is 0.5.
#' @param smooth if TRUE, smooths the cell type classifications. Default is TRUE.
#' @param impute if TRUE, gene expression values are imputed prior to cell type classification. Default is TRUE.
#' @param verbose if TRUE, code will report outputs. Default is TRUE.
#' @param do.normalize if TRUE, cells are normalized to the mean library size. Default is TRUE.
#' @param probability if TRUE, returns the probability associated with each cell type label. Default is TRUE.
#' @param hidden Number of hidden layers in the neural network. Default is 1.
#' @return annotations
#' @export
Signac <- function(E, R , spring.dir = NULL, model.use = "nn", N = 25, num.cores = 1, threshold = 0.5, smooth = T, impute = T, verbose = T, do.normalize = T, probability = F, hidden = 1)
{

  flag = class(E) == "Seurat"
  
  if (flag & impute)
    edges = E@graphs$RNA_nn
  
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
  
  # make sure data are in the same order
  V = V[order(rownames(V)),]
  
  if (class(V) %in% "data.frame")
    V = Matrix::Matrix(as.matrix(V), sparse = T)
  
  # normalize to the mean library size
  if (do.normalize)
  {
    if (!flag)
    {
      V = CID.Normalize(V)
    } else {
      V = CID.Normalize(V@assays$RNA@counts)
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
  
  res = lapply(R$Reference, function(x){
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
      Z = KSoftImpute(E = Z, dM = dM, verbose = F)
      Z = t(apply(Z, 1, function(x){
        normalize(x)
      }))
    }
    
    # build training set
    df = data.frame(dat, celltypes = x$celltypes)
    
    # train a neural network (N times)
    if (model.use == "nn"){
      res = parallel::mclapply(1:N, function(x) {
        nn=neuralnet::neuralnet(celltypes~.,hidden=hidden,data=df, act.fct = 'logistic', linear.output = F)
        Predict = stats::predict(nn, Matrix::t(Z))
        colnames(Predict) <- sort(nn$model.list$response)
        return(Predict)
      }, mc.cores = num.cores)
      res.squared.mean <- Reduce("+", lapply(res, "^", 2)) / N
      res = Reduce(res, f = '+') / N
      res.variance <- res.squared.mean - res^2
      res.sd <- sqrt(res.variance)
      xx = apply(res, 1, which.max)
      celltypes = colnames(res)[xx]
      kmax = apply(res, 1, max)
      celltypes[kmax < threshold] = "Other"
      df = data.frame(celltypes = celltypes, probability = kmax, error = res.sd[xx], percent_features_detected = round(Matrix::colSums(Z != 0) / nrow(Z), digits = 3) * 100, genes_detected = round(Matrix::colSums(E != 0), digits = 3) * 100)
      df$celltypes = as.character(df$celltypes)
    }
    # if desired, run svm
    if (model.use == "svm") {
      model = e1071::svm(celltypes ~ ., data = df, probability=TRUE)
      xx = stats::predict(model, t(Z), probability=TRUE)
      res = attr(xx, "probabilities")
      xx = apply(res, 1, which.max)
      celltypes = colnames(res)[xx]
      kmax = apply(res, 1, max)
      celltypes[kmax < threshold] = "Other"
      df = data.frame(celltypes = celltypes, probability = kmax)
      df$celltypes = as.character(df$celltypes)
    }
    # if desired, run RF
    if (model.use == "rf") {
      res = parallel::mclapply(1:N, function(x) {
        model = randomForest::randomForest(celltypes ~ ., data=df, keep.forest = TRUE)
        stats::predict(model,newdata=Matrix::t(Z),type="prob")
      }, mc.cores = num.cores)
      res = Reduce(res, f = '+') / N
      xx = apply(res, 1, which.max)
      celltypes = colnames(res)[xx]
      kmax = apply(res, 1, max)
      celltypes[kmax < threshold] = "Other"
      df = data.frame(celltypes = celltypes, probability = kmax)
      df$celltypes = as.character(df$celltypes)
    }
      # smooth the output classifications
      if (smooth)
        df$celltypes = CID.smooth(df$celltypes, dM[[1]])
    
    # return probabilities and cell type classifications
    if (probability){
      return(df)
    } else {
      return(df$celltypes)
    }
    })
  
  res$louvain = louvain
  
  if (verbose) {
    tb = proc.time()[3] - ta;
    cat("\n ..........  Exit Signac.\n");
    cat("             Execution time = ", tb, " s.\n", sep = "");
  }
    return(res)
}

#' Generate labels from classifications
#'
#' @param cr list returned by Signac function call.
#' @param E a sparse gene (rows) by cell (column) matrix, or a Seurat object. Rows are HUGO symbols.
#' @param spring.dir If using SPRING, directory to categorical_coloring_data.json. Default is NULL.
#' @param smooth if TRUE, smooths the cell type classifications. Default is TRUE.
#' @return cell type labels (list) for each level of the hierarchy.
#' @export
Generate_lbls = function(cr, spring.dir = NULL, E = NULL, smooth = T)
{
  
  if (!is.null(spring.dir)){
    edges = CID.LoadEdges(data.dir = spring.dir)
    dM = CID.GetDistMat(edges)
  }
  
  flag = class(E) == "Seurat"
  if (flag) {
    edges = E@graphs$RNA_nn
    dM = CID.GetDistMat(edges)
  }
  
  res = list("")
  louvain = cr$louvain
  cr = cr[-which(names(cr) == "louvain")]
  if ("probability" %in% names(cr[[1]])){
    cr = lapply(cr, function(x) {x$celltypes})
  }
  
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
  
  cellstates = res[[length(res)]]
  celltypes = cellstates
  celltypes[celltypes %in% c("B.memory", "B.naive")] = "B"
  celltypes[celltypes %in% c("DC", "Mon.Classical", "Mon.NonClassical", "Neutrophils", "Monocytes", "Macrophages")] = "MPh"
  celltypes[celltypes %in% c("NK", "T.CD4.naive", "T.CD4.memory", "T.CD8", "T.regs", "T.CD8.naive", "T.CD8.memory", "T.CD8.cm","T.CD8.em" , "T.cyto")] = "TNK"
  celltypes[celltypes %in% c("Endothelial", "Fibroblasts", "HSC", "Epithelial")] = "NonImmune"
  immune = res[[1]]  

  # assign Others
  if (!is.null(spring.dir)){
  celltypes = CID.entropy(celltypes, dM)
  immune = CID.entropy(immune, dM)
  # smooth 
  if (smooth) {
    celltypes= CID.smooth(celltypes, dM[[1]])
    cellstates = CID.smooth(cellstates, dM[[1]])
    immune = CID.smooth(immune, dM[[1]])
  }
  }
  logik = immune == "Other" | celltypes == "Other" | cellstates == "Other"
  cellstates[logik] = "Other"
  celltypes[logik] = "Other"
  immune[logik] = "Other"
  
  res$Immune = immune
  if (!is.null(spring.dir))
  {
  do = data.frame(table(louvain[cellstates == "Other"]))
  df = data.frame(table(louvain[louvain %in% do$Var1]))
  logik = (1 - phyper(do$Freq, df$Freq , length(cellstates) - do$Freq, sum(cellstates == "Other"))) < 0.01;
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
      colnames(E) <- lbls
      new_lbls = CID.PosMarkers2(E, cellstates)
      cellstates_novel = new_lbls$lbls
      celltypes_novel[grepl("^[+]", cellstates_novel)] = new_lbls$lbls[grepl("^[+]", cellstates_novel)]
      res$CellTypes_novel = celltypes_novel
      res$CellStates_novel = cellstates_novel
    } 
  }
  }

  res$CellTypes = celltypes
  res$CellStates = cellstates

  names(res)[names(res) == ""] = paste0("L", 1:sum(names(res) == ""))
  
  res$clusters = louvain

  return(res)
}

#' Main function for learning
#'
#' @param E a sparse gene (rows) by cell (column) matrix. Rows are HUGO symbols.
#' @param learned_types Types for learning
#' @param labels cell type labels for the columns of E.
#' @param size Number of bootstrapped samples for machine learning. Default is 1,000.
#' @param logfc.threshold Cutoff for feature selection. Default is 0.25.
#' @param spring.dir if using SPRING, directory to categorical_coloring_data.json. Default is NULL.
#' @param impute if TRUE, performs imputation prior to bootstrapping. Default is TRUE.
#' @return training data frame to be added to a new reference training set
#' @export
SignacLearn <- function (E, learned_types, labels, size = 1000, impute = T, spring.dir = NULL, logfc.threshold = 0.25)
{
  # keep only unique row names
  logik = CID.IsUnique(rownames(E))
  E = E[logik,]
  colnames(E) <- 1:ncol(E)
  
  # run differential expression to find features
  logik = grepl("TotalSeq", rownames(E))
  E = E[!logik,]

  # run feature selection  
  ctrl <- Seurat::CreateSeuratObject(counts = E, project = "CID", min.cells = 0)
  ctrl <- Seurat::NormalizeData(ctrl)
  ctrl <- Seurat::AddMetaData(ctrl, metadata=labels, col.name = "celltypes")
  ctrl <- Seurat::SetIdent(ctrl, value='celltypes')
  mrks = Seurat::FindMarkers(ctrl, ident.1 = learned_types[1], ident.2 = learned_types[2], only.pos = F, logfc.threshold = logfc.threshold)
  # E = E[rownames(E) %in% rownames(mrks),]
  
  # down sample the cells
  #subs = round(sqrt(table(labels)))
  #subsample = unlist(sapply(unique(labels), function(x){
  #      idx = match(x, names(subs))
  #      paste(labels[labels == x], sample(letters[1:subs[idx]], size = sum(labels == x), replace = T), sep = "_")
  #    }))
  #tt = aggregate(Matrix::t(E),by=list(subsample),mean)
  #xx = tt[,1]
  #tt = as.matrix(tt[,-1])
  xx = labels
  
  # set up imputation matrices
  if (impute){
      edges = CID.LoadEdges(data.dir = spring.dir)
      dM = CID.GetDistMat(edges)
      louvain = CID.Louvain(edges = edges)
  }
  
  # bootstrap data
  #mrks2 = mrks[mrks$pct.1 > 0.5 | mrks$pct.2 > 0.5,]
  dat = E[rownames(E) %in% rownames(mrks),]
  mrks$cluster = learned_types[1]
  mrks$cluster[mrks$avg_logFC < 0] = learned_types[2]
  mrks$gene = rownames(mrks)
  
  # run imputation (if desired)
  if (impute)
    Z = KSoftImpute(E = dat, dM = dM, verbose = F)
  
  cts = split.data.frame(mrks, f = mrks$cluster)
    N = lapply(cts, function(x){
      # first sample from cells in cluster 1, size cells
      logik = rownames(dat) %in% x$gene
      dummy = dat[logik,]
      logik = grepl(as.character(x$cluster[1]), labels)
      dummy = dummy[,logik]
      dd = t(apply(dummy, 1, function(z) {
        sample(z, size = size, replace = T)}))
      dd = t(apply(dd, 1, function(z){
        stats::rnorm(n = length(z), mean = mean(z), sd = sd(z))
      }))
      # now sample from cells in cluster 2, size cells with True Negative Expr.
      logik = !rownames(dat) %in% x$gene
      dummy = dat[logik,]
      logik = grepl(as.character(x$cluster[1]), labels)
      dummy = dummy[,logik]
      dd2 = t(apply(dummy, 1, function(z) {
        sample(z, size = size, replace = T)}))
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
   # pca <- prcomp(x = boot[,-ncol(boot)], center = T, scale. = T) 
   # autoplot(pca, data = boot, colour = 'celltypes')
  #library(ggplot2)
  #ggplot(boot, aes(x=TRGC1, y=TRDC, color=celltypes)) + geom_point()
  return(boot)
}

#' Show hierarchical structure
#'
#' @param R Reference dataset
#' @return annotations by hierarchy
#' @export
ShowHierarchy <- function (R)
{
  df = lapply(1:ncol(R$celltypes), function(x){
    if (x == 1)
    {
      unique(R$celltypes[,x])
    } else {
      logik = R$celltypes[,x] != R$celltypes[,x-1]
      dat = R$celltypes[logik,]
      unique(dat[,x])
    }
  })
  names(df) <- colnames(R$celltypes)
  return(df)
}

#' Load data file from directory
#'
#' @param data.dir Directory containing matrix.mtx and genes.txt.
#' @param mfn file name; default is 'matrix.mtx'
#' @return A sparse matrix with rownames equivalent to the names in genes.txt
#' @export
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
  genes <- read.delim(gG, stringsAsFactors = F, header = F)$V1
  if (genes[1] != gsub( "_.*$", "", genes[1] ))
  {
    genes[!grepl("Total", genes)] = gsub( "_.*$", "", genes [!grepl("Total", genes)])
    genes[grepl("Total", genes)]  = gsub( "_", "-",   genes [grepl("Total", genes)])
  }
  if (any(grepl(pattern = "-ENSG00", x = genes)))
  {
    genes <-gsub(pattern = "-ENSG00", replacement = "_ENSG00",x = genes, fixed = T)
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

#' K soft imputation
#'
#' @param E A gene-by-sample count matrix (sparse matrix, matrix, or data.frame) with genes identified by their HUGO symbols (see ?CID.geneconversion), or a list of such matrices, see ?CID.BatchMode.
#' @param dM see ?CID.GetDistMat
#' @param genes.to.use a character vector of genes to impute. Default is NULL.
#' @param do.save If TRUE, imputed matrix is saved in directory. Default is FALSE.
#' @param verbose If TRUE, code reports outputs.
#' @return Imputed values for every gene
#' @export
KSoftImpute <- function(E,  dM = NULL, genes.to.use = NULL, do.save = F, verbose = T)
{
  # check inputs
  stopifnot(class(E) %in% c("dgCMatrix","dgTMatrix", "matrix", "data.frame"))
  stopifnot(!is.null(rownames(E)));
  if (class(E) %in% c("matrix", "data.frame"))
    E = Matrix::Matrix(as.matrix(E), sparse = T)
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
  
  # if is null data.dir, run PCA + KNN
  if (is.null(dM))
  {
    dM = CID.GetNeighbors(E, normalize = normalize, min_counts = 3, min_cells = 3, min_vscore_pctl = 85, num_pc = 30, k_neigh = 3)
    dM = CID.GetDistMat(dM)
  }
  
  # g = dM[[1]] %*% bas + 1/2 * dM[[2]] %*% bas
  g = dM[[1]] %*% bas
  dd = g / (Matrix::rowSums(g) + 1)
  diag(dd) <- 1
  E_new = E %*% dd;
  E_new = CID.Normalize(E_new)
  
  if (do.save)
  {
    data.dir = dirname(data.dir)
    Matrix::writeMM(Matrix::Matrix(E_new, sparse = TRUE), file = paste(data.dir, "matrix_ksoft_imputed.mtx", sep = "/"))
    write.table(rownames(E_new), file = paste(data.dir, "genes_ksoft_imputed.txt", sep = "/"))
  }
  
  if (verbose) {
    tb = proc.time()[3] - ta;
    cat("\n ..........  Exit KSoftImpute.\n");
    cat("             Execution time = ", tb, " s.\n", sep = "");
  }
  
  return(E_new)
}

#' Get edges that are either pre-computed, or generate new edges
#'
#' @param data.dir directory
#' @param E see Signac
#' @return edges for cell-cell similarity network
#' @export
CID.GetEdges <- function(E = NULL, data.dir = NULL)
{
  if (is.null(data.dir))
  {
    edges = CID.GetNeighbors(E, normalize = T, min_counts = 3, min_cells = 3, min_vscore_pctl = 85, num_pc = 30, k_neigh = 4)
  } else {
    edges = CID.LoadEdges(data.dir)
  }
  return(edges)
}

#' Get indices of training markers
#'
#' @param E A gene-by-sample count matrix (sparse matrix, matrix, or data.frame) with genes identified by their HUGO symbols (see ?CID.geneconversion), or a list of such matrices, see ?CID.BatchMode.
#' @param data.dir if default, uses the standard Signac markers (see ?CID.SeeMarkers).
#' @param method.use either 'max.genes.detected' or 'min.entropy'
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
      which(df$cells == x)[order(df$entropy[df$cells == x])[1:round(sum(df$cells == x) * 0.5)]]
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
      which(df$cells == x)[order(df$entropy[df$cells == x], decreasing = TRUE)[1:round(0.5 * sum(df$cells == x))]]
    })
    names(q) <- unique(df$cells)
    
  }
  
  return(q)
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
  if ("clusters" %in% names(cr))
  {
    Q = as.character(cr$clusters)
    json_data$Clusters_Louvain$label_list = Q
    Ntypes = length(unique(Q))
    qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) #len = 74
    #pie(rep(1,num_col), col=(col_vector[1:num_col]))
    col_palette <- as.list(col_vector[1:Ntypes]); # or sample if you wish
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
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) #len = 74
    #pie(rep(1,num_col), col=(col_vector[1:num_col]))
    C[[1]][C[[1]] == ""] <- col_vector[1:Ntypes]; # or sample if you wish
    json_data$CellTypes_novel$label_colors = as.list(C[[1]])
  }
  if ("CellStates_novel" %in% names(cr))
  {
    Q = cr$CellStates_novel
    json_data$CellStates_novel$label_list = Q
    C = get_colors(Q)
    Ntypes = sum(C[[1]] == "")
    qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) #len = 74
    #pie(rep(1,num_col), col=(col_vector[1:num_col]))
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
        json_data = json_data_backup;
      }        
    }
    json_out = jsonlite::toJSON(json_data, auto_unbox = TRUE)
    write(json_out,paste(new.dirs[j],fn,sep="/"))
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
  main_types = c("B"     ,     "Epithelial", "Fibroblasts"        ,
                 "MPh"         , "Plasma.cells"  , "Endothelial",
                 "TNK"         ,     "Other"     , "NonImmune"          ,
                 "Immune")
  main_colors = c("#aa3596"     ,     "#387f50"   , "#a2ba37" ,
                  "#f1ff51"     ,     "#d778de"   , "#73de97" ,
                  "#90c5f4"     ,     "#c0c0c0"   , "#bb7fc7" ,
                  "#7fc97f")
  
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

#' Get main cell types from hierarchy
#'
#' @param R reference matrix 
#' @return a vector of cell type identities given the hierarchy
#' @export
get_celltypes <- function(R) {
  
  # we need to populate a data frame with labels at each level of the hierarchy
  celltypes = unique(unlist(R$immune_hierarchy, use.names = F)) # bottom level
  qq = do.call(c, unlist(R$immune_hierarchy, recursive=FALSE))
  
  # construct a data frame from the list
  R$immune_hierarchy[names(R$immune_hierarchy)]

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
#' @param genes_to_use gene list for KNN graph building
get_knn_graph2 <- function(X, k=4, np, genes_to_use)
{
  logik = CID.IsUnique(rownames(X))
  X = X[logik,]
  colnames(X) <- 1:ncol(X)
  ctrl <- suppressWarnings(Seurat::CreateSeuratObject(X))
  ctrl <- Seurat::ScaleData(ctrl, verbose = F)
  ctrl <- Seurat::RunPCA(ctrl, features = genes_to_use, pcs.compute = np, do.print = F)
  ctrl <- Seurat::FindNeighbors(object = ctrl, reduction = "pca", dims = 1:min(c(np, 50)), k.param = k)
  return(ctrl@graphs$RNA_nn)
}

#' Save network h5 files
#'
#' @param A an adjacency matrix with row names and columns names genes (gene by gene)
#' @param data.dir directory where the networks are saved
#' @return a saved network file
#' @export
SaveNetworksH5 <- function(A, data.dir)
{
    
data.dir = gsub("\\/$", "", data.dir, perl = TRUE);
    
if (!dir.exists(data.dir))
  dir.create(data.dir)
    
if (!"list" %in% class(A))
  A = list(A)
    
    cat(" ..........  Entry in SaveNetworkH5: \n");
    ta = proc.time()[3];
    cat("             Number of genes:", nrow(A[[1]]), "\n");
    cat("             Number of networks:", length(A), "\n");
    
    if (is.null(names(A)))
      names(A) <- paste0("Network", seq_along(1:length(D)))
    
    data.dirs = paste(data.dir, names(A), sep = "/")
    
    q = lapply(data.dirs, function(x){
      if (!dir.exists(x))
        dir.create(x)})
    
    d = suppressWarnings(
      mapply(function(x,y){
        fn = "network_sparse_genes.h5"
        fn = paste(y, fn, sep = "/")
        rhdf5::h5createFile(fn)
        rhdf5::h5createGroup(fn, "edges")
        lapply(rownames(x), function(y){
          edges = which(x[rownames(x) == y,] != 0)
          rhdf5::h5write.default(edges, file = fn, name = paste0("edges/", y))
        })
        rhdf5::h5write.default(dim(x), file = fn, name="shape")
      }, x = A, y = data.dirs)
    )
    
    d = suppressWarnings(
      mapply(function(x,y){
        fn = "network_total.h5"
        fn = paste(y, fn, sep = "/")
        rhdf5::h5createFile(fn)
        rhdf5::h5write.default(x@Dimnames[[2]], file = fn, name = "genes")
        rhdf5::h5write.default(x@x, file = fn, name = "data")
        rhdf5::h5write.default(dim(x), file = fn, name="shape")
        rhdf5::h5write.default(x@i, file = fn, name="indices") # already zero-indexed.
        rhdf5::h5write.default(x@p, file = fn, name="indptr")
      }, x = A, y = data.dirs))
    
    rhdf5::h5closeAll()
    tb = proc.time()[3] - ta;
    cat(" ..........  Exit SaveNetworksH5 \n");
    cat("             Execution time = ", tb, " s.\n", sep = "");
}

#' Load gene from network h5 file
#'
#' @param gene A gene to query from the network
#' @param filename file name of network. Default is "network.hdf5"
#' @return a saved network file
#' @export
GetGeneFromNetwork <- function(gene, filename = "network_sparse_genes.h5")
{
  if (!requireNamespace("hdf5r", quietly = TRUE))
    stop("Please install hdf5r to read HDF5 files")
  if (!file.exists(filename))
    stop("File not found")
  infile = hdf5r::H5File$new(filename)
  E = Matrix::Matrix(0, ncol = 1, nrow = infile[["shape"]][][1], sparse = T)
  nodes <- infile[[paste("edges/", gene, sep = "/")]][]
  E[nodes] = 1
  #E = as.matrix(E)
  rownames(E) <- names(infile[["edges"]])
  return(E)
}

#' Load gene from network h5 file
#'
#' @param filename file name of network. Default is "network.hdf5"
#' @return a saved network file
#' @export
GetGeneListFromNetwork <- function(filename = "network.h5")
{
  if (!requireNamespace("hdf5r", quietly = TRUE))
    stop("Please install hdf5r to read HDF5 files")
  if (!file.exists(filename))
    stop("File not found")
  infile = hdf5r::H5File$new(filename)
  return(names(infile[["edges"]]))
}

#' Load gene from network h5 file
#'
#' @param filename file name of network. Default is "network.hdf5"
#' @return a saved network file
#' @export
LoadNetwork <- function(filename = "network_total.h5")
{
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Please install hdf5r to read HDF5 files")
  }
  if (!file.exists(filename)) {
    stop("File not found")
  }
  infile <- hdf5r::H5File$new(filename)
  
  counts <- infile[["data"]]
  indices <- infile[["indices"]]
  indptr <- infile[["indptr"]]
  shp <- infile[["shape"]]
  features <- infile[["genes"]][]
  sparse.mat <- Matrix::sparseMatrix(i = indices[] + 1, p = indptr[],
                                     x = as.numeric(counts[]), dims = shp[], giveCsparse = FALSE)
  rownames(sparse.mat) <- features
  colnames(sparse.mat) <- features
  sparse.mat <- as(object = sparse.mat, Class = "dgCMatrix")
 
  return(sparse.mat)
  
}
