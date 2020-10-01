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