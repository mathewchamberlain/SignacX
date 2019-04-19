genesig_wrapper <- function()

{
  ## differential expression analysis with MAST
  data("S")

  # Establish master cell types
  master_list = list( B.cells   =   c("B.cells.naive","B.cells.memory", "naive B-cells", "Memory B-cells", "Class-switched memory B-cells") ,
                      Epithelial = "Epithelial",
                      Epithelial.Urinary = "Epithelial.Urinary",
                      Erythro = "Erythro",
                      Fibroblasts = "Fibroblasts",
                      Granulocytes  =   c("Eosinophils","Neutrophils","Mast.cells.resting","Mast.cells.activated", "Mast.Progenitors"),
                      MPh =   c("Monocytes","Macrophages.M0","Macrophages.M1","Macrophages.M2","Dendritic.cells.resting","Dendritic.cells.activated",
                                "Macrophages", "Macrophages M1", "Macrophages M2", "Monocytes") ,
                      pDCs = "pDCs",
                      Plasma.cells = "Plasma.cells",
                      Platelets = "Platelets",
                      Secretory = "Secretory",
                      TNK =   c("T.cells.CD8", "T.cells.CD4.naive","T.cells.CD4.memory.activated","T.cells.CD4.memory.resting",
                                "T.cells.follicular.helper","T.cells.regs", "NK.cells.resting", "Cyto.T.cells", "NK.cells",
                                "Tregs","CD4+ T-cells","CD4+ Tcm", "CD4+ Tem", "CD8+ T-cells", "CD8+ Tcm", "CD8+ Tem", "NK cells"),
                      HSC = c("CLP", "CMP", "GMP", "HSC", "Megakaryocytes", "MEP", "MPP")
                      )
  # rename columns
  for (j in 1:length(master_list))
    colnames(S)[colnames(S) %in% master_list[[j]]] = names(master_list)[j]

  # define positive marker sets
  pos = CID.PosMarkers(S)

  # define negative marker sets
  neg = CID.NegMarkers(S)

  # bind them, these are top of hierachy
  markers = rbind(pos, neg)
  markers = markers[order(markers$`Cell population`), ]

  usethis::use_data(markers, overwrite = TRUE)

}

#' Run DEG analysis with Seurat wrapper
#'
#' @param E Expression matrix with genes for rows, samples for columns
#' @return List where each element contains the DEG tables for the one vs. all comparison
#' @export
CID.PosMarkers <- function(E)
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
  cts = unique(lbls);
  #cts2 = data.frame(ident.1 = sort(cts[grepl("Samples 1-2", cts)]), ident.2 = sort(cts[grepl("Samples 3-12", cts)]))
  for (j in 1:length(cts))
  {
    dd = Seurat::FindMarkers(ctrl, ident.1 = cts[j], ident.2 = NULL, test.use = "MAST", min.cells.group = 0, max.cells.per.ident = 200, logfc.threshold = 1, pseudocount.use = 1)
    dd$GeneSymbol = rownames(dd)
    dd$celltype = cts[j]
    outs[[j]] = dd
  }
  names(outs) <- cts

  outs = do.call(rbind, lapply(outs, function(x) x[x$p_val_adj < 0.05 & x$avg_logFC > 0.75, ]))
  geneset = data.frame(genes = outs$GeneSymbol, celltype = outs$celltype, Polarity = "+")
  colnames(geneset) <- c("HUGO symbols", "Cell population", "Polarity")

  return(geneset)

}

#' Run DEG analysis with Seurat wrapper
#'
#' @param E Expression matrix with genes for rows, samples for columns
#' @return List where each element contains the DEG tables for the one vs. all comparison
#' @export
CID.NegMarkers <- function(E)
{
  library(scater)
  # get lbls
  lbls = colnames(S)
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
    dd = Seurat::FindMarkers(ctrl, ident.1 = cts[j], ident.2 = NULL, test.use = "MAST", min.cells.group = 0, max.cells.per.ident = 200, logfc.threshold = 2, pseudocount.use = 0)
    dd$GeneSymbol = rownames(dd)
    dd$celltype = cts[j]
    outs[[j]] = dd
  }
  names(outs) <- cts

  outs = do.call(rbind, lapply(outs, function(x) x[x$p_val_adj < 0.05 & x$avg_logFC < -2, ]))
  geneset = data.frame(genes = outs$GeneSymbol, celltype = outs$celltype, Polarity = "-")
  colnames(geneset) <- c("HUGO symbols", "Cell population", "Polarity")

  return(geneset)

}
