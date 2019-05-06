genesig_wrapper <- function()

{
  ## differential expression analysis with MAST
  data("S")

  # Establish master cell types
  master_list = list( B.cells   =   c("B.cells.naive","B.cells.memory", "naive B-cells", "Memory B-cells", "Class-switched memory B-cells") ,
                      Epithelial = "Epithelial",
                      Epithelial.Urinary = "Epithelial.Urinary",
                      Erythro = c("Erythro", "Erythrocytes"),
                      Fibroblasts = "Fibroblasts",
                      Granulocytes  =   c("Eosinophils","Neutrophils","Mast.cells.resting","Mast.cells.activated", "Mast.Progenitors"),
                      MPh =   c("Monocytes","Macrophages.M0","Macrophages.M1","Macrophages.M2","Dendritic.cells.resting","Dendritic.cells.activated",
                                "Macrophages", "Macrophages M1", "Macrophages M2", "Monocytes", "DC") ,
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
  pos = GetPosMarkers(S)

  # define negative marker sets
  neg = GetNegMarkers(S)

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
GetPosMarkers <- function(E)
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
  for (j in 1:length(cts))
  {
    dd = Seurat::FindMarkers(ctrl, ident.1 = cts[j], ident.2 = NULL, min.cells.group = 0, logfc.threshold = 1, test.use = "MAST")
    dd$GeneSymbol = rownames(dd)
    dd$celltype = cts[j]
    outs[[j]] = dd
  }
  names(outs) <- cts
  DD = do.call(rbind, outs)
  
  # first remove anything with p_val_adj > 0.05
  res = lapply(outs, function(x) x[x$p_val_adj < 0.05,])
  
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
  
  # Fourth we check correlation structure of the remaining markers.
  # To do this, I group all markers for all cell types into groups of two.
  # I then run hierarchical clustering on the distance matrix, and then 
  # cut the clustering result into two groups. Markers which belong to the 
  # wrong group are removed.
  tt = t(combn(1:length(res), 2))
  Q = apply(tt, 1, function(x) rbind(res[[x[1]]], res[[x[2]]]))
  
  # order Q by gene name
  Q = lapply(Q, function(x){
    x[order(x$GeneSymbol),]
  })
  
  # hierarchical clustering on distance matrix
  QC = lapply(Q, function(x){
    expr_mat = E[rownames(E) %in% x$GeneSymbol,]
    cor_mat = qlcMatrix::cosSparse(Matrix::t(expr_mat))
    ix = hclust(dist(cor_mat))
    grps = cutree(ix, k = 2)
    x$groups_by_two = grps;
    x
  })
  
  # inspect the results
  TB = lapply(QC, function(x){
    table(data.frame(from = x$celltype, to = x$groups_by_two))
  })
  
  # remove any marker assigned to the wrong cluster
  QC_idx = mapply(function(x,y){
    idx.max = apply(x , 1, which.max);
    logik = idx.max[match(y$celltype, names(idx.max))] == y$groups_by_two;
    y = y[logik,]
    y
  }, x = TB, y = QC, SIMPLIFY = FALSE)
  
  # study how many markers are stable predictors
  q = do.call(rbind, QC_idx)
  df = data.frame(table(q$GeneSymbol))
  logik = df$Freq == 12;
  gns_stable = as.character(df$Var1)[logik]
  
  q = do.call(rbind, res)
  # table(q$celltype[q$GeneSymbol %in% gns_stable])
  q = q[q$GeneSymbol %in% gns_stable,]
  
  geneset = data.frame(genes = q$GeneSymbol, celltype = q$celltype, Polarity = "+")
  colnames(geneset) <- c("HUGO symbols", "Cell population", "Polarity")

  return(geneset)

}

#' Run DEG analysis with Seurat wrapper
#'
#' @param E Expression matrix with genes for rows, samples for columns
#' @return List where each element contains the DEG tables for the one vs. all comparison
#' @export
GetNegMarkers <- function(E)
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
    dd = Seurat::FindMarkers(ctrl, ident.1 = cts[j], ident.2 = NULL, test.use = "MAST", min.cells.group = 0, logfc.threshold = 1)
    dd$GeneSymbol = rownames(dd)
    dd$celltype = cts[j]
    outs[[j]] = dd
  }
  names(outs) <- cts
  
  # first remove anything with p_val_adj > 0.05
  res = lapply(outs, function(x) x[x$p_val_adj < 0.05,])
  
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
  
  # Fourth we check correlation structure of the remaining markers.
  # To do this, I group all markers for all cell types into groups of two.
  # I then run hierarchical clustering on the distance matrix, and then 
  # cut the clustering result into two groups. Markers which belong to the 
  # wrong group are removed.
  tt = t(combn(1:length(res), 2))
  Q = apply(tt, 1, function(x) rbind(res[[x[1]]], res[[x[2]]]))
  
  # order Q by gene name
  Q = lapply(Q, function(x){
    x[order(x$GeneSymbol),]
  })
  
  # hierarchical clustering on distance matrix
  QC = lapply(Q, function(x){
    expr_mat = E[rownames(E) %in% x$GeneSymbol,]
    cor_mat = qlcMatrix::cosSparse(Matrix::t(expr_mat))
    ix = hclust(dist(cor_mat))
    grps = cutree(ix, k = 2)
    x$groups_by_two = grps;
    x
  })
  
  # inspect the results
  TB = lapply(QC, function(x){
    table(data.frame(from = x$celltype, to = x$groups_by_two))
  })
  
  # remove any marker assigned to the wrong cluster
  QC_idx = mapply(function(x,y){
    idx.max = apply(x , 1, which.max);
    logik = idx.max[match(y$celltype, names(idx.max))] == y$groups_by_two;
    y = y[logik,]
    y
  }, x = TB, y = QC, SIMPLIFY = FALSE)
  
  # study how many markers are stable predictors
  q = do.call(rbind, QC_idx)
  df = data.frame(table(q$GeneSymbol))
  logik = df$Freq == length(res); sum(logik)
  if (sum(logik) == 0)
  {
    return(NULL)
  } else {
    gns_stable = as.character(df$Var1)[logik]
    
    q = do.call(rbind, res)
    # table(q$celltype[q$GeneSymbol %in% gns_stable])
    q = q[q$GeneSymbol %in% gns_stable,]
    
    geneset = data.frame(genes = q$GeneSymbol, celltype = q$celltype, Polarity = "+")
    colnames(geneset) <- c("HUGO symbols", "Cell population", "Polarity")
    
    return(geneset)
  }

}
