genesig_wrapper_v2 <- function()

{
  ## differential expression analysis with edgeR
  #library(edgeR)
  #library(Biobase)
  #library(biomaRt)
  # Establish master cell types

  data("S")

  master_list = list( B.cells       =   list(B.cells.naive  = c("B.cells.naive", "naive B-cells"),
                                             B.cells.memory = c("B.cells.memory","Memory B-cells","Class-switched memory B-cells")) ,
                      Granulocytes  =   list(Not.Mast       = c("Eosinophils","Neutrophils"),
                                             Mast           = c("Mast.cells.resting","Mast.cells.activated", "Mast.Progenitors")),
                      Mast          =   c("Mast.cells.resting","Mast.cells.activated", "Mast.Progenitors"),
                      Not.Mast      =   c("Eosinophils","Neutrophils"),
                      MPh           =   list(Monocytes      = c("Monocytes"),
                                             Macrophages    = c("Macrophages.M0","Macrophages.M1","Macrophages M1","Macrophages.M2","Macrophages M2"),
                                             Dendritic      = c("Dendritic.cells.resting","Dendritic.cells.activated", "DC")),
                      Macrophages   =   list(Macrophages.M0 = c("Macrophages.M0"),
                                             Macrophages.M1 = c("Macrophages.M1","Macrophages M1"),
                                             Macrophages.M2 = c("Macrophages.M2","Macrophages M2")),
                      TNK           =   list(T.cells        = c("T.cells.CD8", "T.cells.CD4.naive","T.cells.CD4.memory.activated","T.cells.CD4.memory.resting",
                                             "T.cells.follicular.helper","T.cells.regs", "NK.cells.resting", "Cyto.T.cells",
                                             "Tregs","CD4+ T-cells","CD4+ Tcm", "CD4+ Tem", "CD8+ T-cells", "CD8+ Tcm", "CD8+ Tem"),
                                             NK             = c("NK cells", "NK.cells", "NK.cells.resting")),
                      T.cells       =   list(T.cells.CD8    = c("T.cells.CD8", "Cyto.T.cells", "CD8+ T-cells", "CD8+ Tcm", "CD8+ Tem"),
                                             T.CD4.FH.regs  = c("T.cells.follicular.helper","T.cells.regs", "T.cells.CD4.naive","T.cells.CD4.memory.activated",
                                                                "T.cells.CD4.memory.resting","Tregs","CD4+ T-cells", "CD4+ Tcm", "CD4+ Tem")),
                      T.CD4.FH.regs =   list(T.regs         = c("T.cells.regs","Tregs"),
                                             T.CD4.FH       = c("T.cells.follicular.helper","T.cells.CD4.naive","T.cells.CD4.memory.activated",
                                                                "T.cells.CD4.memory.resting","CD4+ T-cells", "CD4+ Tcm", "CD4+ Tem")),
                      T.CD4.FH      =   list(T.cells.FH     = c("T.cells.follicular.helper"),
                                             T.cells.CD4    = c("T.cells.CD4.naive","T.cells.CD4.memory.activated",
                                                                "T.cells.CD4.memory.resting","CD4+ T-cells", "CD4+ Tcm", "CD4+ Tem"))
                      )
                           
  outs = list("")
  
  for (k in 1:length(master_list)){
    Q = master_list[[k]]
    D = S[, colnames(S) %in% unlist(Q)] # initialize
    # rename columns
    if (class(Q) == "list")
    {
      for (j in 1:length(Q))
        colnames(D)[colnames(D) %in% Q[[j]]] = names(Q)[j]  
    }
    logfc = 1
    pos = NULL
    M = 0
    N = length(unique(colnames(D)))
    L = 0;
    P = 3;
    pseudo = 1
    while (is.null(pos) || M != N)
    {
      pos = CID.PosMarkers(D, pseudocount.use = pseudo, logfc.threshold = logfc)
      logfc = logfc - 0.25;
      if (logfc == 0)
      {
        pseudo = pseudo - 0.1;
        logfc = 0.25;
      }
      if (!is.null(pos))
      {
        M = length(unique(pos$`Cell population`))
        if (M == N)
        {
        L = min(table(pos$`Cell population`))
        if (L < P)
          M = 0;
        }
      }
    }
    neg = CID.NegMarkers(D)
    markers = rbind(pos, neg)
    markers = markers[order(markers$`Cell population`), ]
    outs[[k]] = markers
    cat("Step complete:", k, "\n")
  }
  names(outs) <- names(master_list)
  cellstate_markers = outs
  usethis::use_data(cellstate_markers, overwrite = TRUE)
  usethis::use_data(markers, cellstate_markers, overwrite = TRUE, internal = TRUE)
}

#' Run DEG analysis with Seurat wrapper
#'
#' @param E Expression matrix with genes for rows, samples for columns
#' @return List where each element contains the DEG tables for the one vs. all comparison
#' @export
CID.PosMarkers <- function(E, logfc.threshold = 1, pseudocount.use = 1)
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
  for (j in 1:length(cts))
  {
    dd = Seurat::FindMarkers(ctrl, ident.1 = cts[j], ident.2 = NULL, test.use = "MAST", min.cells.group = 0, max.cells.per.ident = 200, logfc.threshold = logfc.threshold, pseudocount.use = pseudocount.use)
    dd$GeneSymbol = rownames(dd)
    dd$celltype = cts[j]
    outs[[j]] = dd
  }
  names(outs) <- cts
  outs = do.call(rbind, lapply(outs, function(x) x[x$p_val < 0.05 & x$avg_logFC > logfc.threshold, ]))
  if (nrow(outs) > 0)
  {
    geneset = data.frame(genes = outs$GeneSymbol, celltype = outs$celltype, Polarity = "+")
    colnames(geneset) <- c("HUGO symbols", "Cell population", "Polarity")
  } else {
    geneset = NULL
  }

  return(geneset)

}

#' Run DEG analysis with Seurat wrapper
#'
#' @param E Expression matrix with genes for rows, samples for columns
#' @return List where each element contains the DEG tables for the one vs. all comparison
#' @export
CID.NegMarkers <- function(E, logfc.threshold = -1, pseudocount.use = 1)
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
  for (j in 1:length(cts))
  {
    dd = Seurat::FindMarkers(ctrl, ident.1 = cts[j], ident.2 = NULL, test.use = "MAST", min.cells.group = 0, max.cells.per.ident = 200, logfc.threshold = logfc.threshold, pseudocount.use = pseudocount.use)
    dd$GeneSymbol = rownames(dd)
    dd$celltype = cts[j]
    outs[[j]] = dd
  }
  names(outs) <- cts
  outs = do.call(rbind, lapply(outs, function(x) x[x$p_val < 0.05 & x$avg_logFC < logfc.threshold, ]))
  if (nrow(outs) > 0)
  {
    geneset = data.frame(genes = outs$GeneSymbol, celltype = outs$celltype, Polarity = "-")
    colnames(geneset) <- c("HUGO symbols", "Cell population", "Polarity")
  } else {
    geneset = NULL
  }
  
  return(geneset)
  
}