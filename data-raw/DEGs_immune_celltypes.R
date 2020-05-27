genesig_wrapper <- function()
{
  source('./data-raw/get_celltype_markers.R')
  ## differential expression analysis with MAST
  load("./data/blueprint_encode.rda")
  load("./data/hpca.rda")
  load("./data/LM22.rda")
  
  Refs = merge(LM22$data, hpca$data, by = 'row.names', all = T)
  rownames(Refs) <- Refs$Row.names
  Refs = Refs[,-1]
  Refs = merge(Refs, blueprint_encode$data, by = 'row.names', all = T)
  xx <- Refs$Row.names
  Refs = Refs[,-1]
  Refs[is.na(Refs)] <- 0
  #Refs = Matrix::Matrix(as.matrix(Refs), sparse = T)
  rownames(Refs) <- xx
  
  # define markers
  colnames(Refs) <- c(LM22$immune_types, hpca$immune_types, blueprint_encode$immune_types)
  lfc = 0.25;
  pos = GetMarkers(Refs, lfc = lfc, pval = 0.01, polarity = "+")
  pos$`Cell population` = as.character(pos$`Cell population`)
  # create negative markers too
  neg = pos
  neg$Polarity = "-"
  logik = neg$`Cell population` == "Immune";
  neg$`Cell population`[logik] <- "NonImmune"
  neg$`Cell population`[!logik] <- "Immune"
  
  immune_markers = rbind(pos, neg)
  
  usethis::use_data(immune_markers, overwrite = TRUE)
}
