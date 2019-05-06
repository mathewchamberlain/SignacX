define_geneset <- function()
{
  ## Load LM22 expression data
  data_dir = "./data-raw/"
  LM22_full <- read.delim(paste(data_dir,"LM22-ref-sample.txt",sep=""), row.names=1, stringsAsFactors=FALSE) # expression matrix
  LM22 <- read.delim(paste(data_dir,"LM22.txt",sep=""), row.names=1, stringsAsFactors=FALSE) # expression matrix
  
  ## load blueprint_encode.rda data from singleR
  load(paste(data_dir, "blueprint_encode.rda", sep = ""))
  
  # remove RP + MT genes
  logik = grepl( "^MT-",rownames(LM22)) | grepl("^RP[0-9]", rownames(LM22)) | grepl("^RPS", rownames(LM22)) | grepl("^RPL", rownames(LM22))
  LM22 = LM22[!logik,]
  
  # load genes from single cell datasets
  data.dirs = list.files("../Data_Signac", full.names = T)
  gns = lapply(data.dirs, function(x) {read.delim(paste(x, "genes.txt", sep ="/"), stringsAsFactors = F, header = F)$V1})
  gns = lapply(gns, function(x){
    if (x[1] != gsub( "_.*$", "", x[1] )){
      x = gsub( "_.*$", "", x )
    } else {
      x
    }
  })
  
  gns = c(gns, list(rownames(LM22_full)), list(rownames(blueprint_encode$data)))
  
  # define geneset
  q = Reduce(intersect, gns)
  
  # remaining are 9,606 genes as candidates for cell type markers
  write.table(q, file = paste0(data_dir, "candidate_genes.txt"))
}

geneset_preprocessing <- function()
{
  library(edgeR)
  library(Matrix)
  
  # load candidate genes
  candidate_genes <- read.csv("/site/ne/data/bh-results/C/CHAMBERLAIN.Mat/pipelines/Signac/data-raw/candidate_genes.txt", sep="", stringsAsFactors=FALSE)$x

  ## Load LM22 expression data
  data_dir = "./data-raw/"
  LM22_full <- read.delim(paste(data_dir,"LM22-ref-sample.txt",sep=""), row.names=1, stringsAsFactors=FALSE) # expression matrix
  LM22 <- read.delim(paste(data_dir,"LM22.txt",sep=""), row.names=1, stringsAsFactors=FALSE) # expression matrix

  ## load singleR data
  load(paste(data_dir, "blueprint_encode.rda", sep = ""))
  
  ## the following columns were removed during outlier analysis:
  outliers = c("Monocyte.Day1.7..HG.U133A...IRIS_GSE22886.GSM565348."         ,"A_TS_RN_gdTcells_U133A..Chtanova_immune.A_TS_RN_gdTcells_U133A.",
               "A_TS_RN_gdTcellsREP_A..Chtanova_immune.A_TS_RN_gdTcellsREP_A.",
               "Bcell.naive.2..HG.U133A...IRIS_GSE22886.GSM565309."           ,"Bcell.naive.5..HG.U133A...IRIS_GSE22886.GSM565312.", "Bcell.naive.7..HG.U133A...IRIS_GSE22886.GSM565314.",
               "Bcell.Memory_IgM.1..HG.U133A...IRIS_GSE22886.GSM565319."      ,"Bcell.Memory_IgG_IgA.2..HG.U133A...IRIS_GSE22886.GSM565316.",
               "Bcell.Memory_IgM.2..HG.U133A...IRIS_GSE22886.GSM565320."      ,"Bcell.Memory_IgM.3..HG.U133A...IRIS_GSE22886.GSM565321." ,
               "Bcell.Memory_IgG_IgA.3..HG.U133A...IRIS_GSE22886.GSM565317."  ,"Bcell.Memory_IgG_IgA.4..HG.U133A...IRIS_GSE22886.GSM565318.",
               "NKcell.control.3..HG.U133A...IRIS_GSE22886.GSM565295."        ,"NKcell.IL2stimulated.4..HG.U133A...IRIS_GSE22886.GSM565300." ,
               "NKcell.IL15stimulated.4..HG.U133A...IRIS_GSE22886.GSM565305." ,"NKcell.IL2stimulated.5..HG.U133A...IRIS_GSE22886.GSM565301." ,
               "NKcell.IL15stimulated.5..HG.U133A...IRIS_GSE22886.GSM565306." ,"NKcell.IL2stimulated.1..HG.U133A...IRIS_GSE22886.GSM565297." ,
               "NKcell.IL15stimulated.1..HG.U133A...IRIS_GSE22886.GSM565302." ,"NKcell.IL15stimulated.2..HG.U133A...IRIS_GSE22886.GSM565303.",
               "NKcell.IL15stimulated.6..HG.U133A...IRIS_GSE22886.GSM565307." ,"NKcell.IL2stimulated.2..HG.U133A...IRIS_GSE22886.GSM565298." ,
               "NKcell.IL2stimulated.3..HG.U133A...IRIS_GSE22886.GSM565299."  ,"NKcell.IL15stimulated.3..HG.U133A...IRIS_GSE22886.GSM565304.",
               "A_LW_macroctrl_U133A_130503..Chtanova_immune.A_LW_macroctrl_U133A_130503.", "Monocyte.Day0.6..HG.U133A...IRIS_GSE22886.GSM565335.")

  LM22_full <- LM22_full[,!names(LM22_full) %in% outliers]

  ## clean cell labels; map to LM22 cell labels
  y = colnames(LM22)
  x = colnames(LM22_full)
  y = colnames(LM22) # B.cells.naive, B.cells.memory, Plasma.cells, T.cells.CD8, T.cells.CD4.naive, T.cells.CD4.memory.resting, T.cells.CD4.memory.activated, T.cells.follicular.helper
  # T.cells.regulatory..Tregs., T.cells.gamma.delta, NK.cells.resting, NK.cells.activated, Monocytes, Macrophages.M0, Macrophages.M1, Macrophages.M2, Dendritic.cells.resting,
  # Dendritic.cells.activated, Mast.cells.resting, Mast.cells.activated, Eosinophils, Neutrophils
  y = y[-10] # removed during outlier cleanup : T cells gamma delta
  y = y[-11] # removed during outlier cleanup : NK cells activated

  rm(LM22); gc() # not needed anymore

  # rename one mast cell column for grepl
  x[x == "A_LW_mastcellctrl_U133A..Chtanova_immune.A_LW_mastcellctrl_U133A." ] = "A_LW_ControlMASTCELL_U133A..Chtanova_immune.A_MF_ControlMASTCELL_U133A."

  ## get simplified expr matrix
  flag = c("naive", "Bcell.Memory", "PlasmaCell", "CD8", "TN_U133A", "MemoryTcell.RO.unactivated" ,"MemoryTcell.RO.activated", "CXCR5hiICOShi", "Treg_" , "NKcell.control" ,
           "Monocyte.Day0", "Monocyte.Day7", "classical.or.M1.", "Alternative.or.M2", "DendriticCell.Control", "DendriticCell.LPSstimulated", "ControlMASTCELL", "IgE", "sinophil", "eutr")

  logik2 = grepl(flag[1], x) # initialize logik2
  dummy = x; # initialize dummy

  for (i in 1:length(flag))
  {
    logik = grepl(flag[i], x, ignore.case = T)
    dummy[logik] <- y[i]
    logik2 = logik2 | logik
  }

  # remove anything in LM22 full that we are not using
  LM22 = LM22_full[,logik2]; rm(LM22_full); gc()
  dummy = dummy[logik2];

  # remove all non-candidate genes
  LM22 = LM22[rownames(LM22) %in% candidate_genes,]

  ## load expression matrix for epithelial, epithelial urinary, secretory cell types
  data.dir = "../Data_Signac/AMP_Phase1_SLE/"
  spring.dir = "../Data_Signac/AMP_Phase1_SLE/FullDataset_v1/"
  E <- CID.LoadData(data.dir = data.dir)

  # remove all non-candidate genes
  E = E[rownames(E) %in% candidate_genes,]
  
  # normalize to these genes
  E = CID.Normalize(E)
  
  # read the JSON from SPRING
  json_file = 'categorical_coloring_data.json'
  gJ <- paste(spring.dir,json_file,sep = "")
  json_data <- rjson::fromJSON(file=gJ)

  # define cell indices of interest
  EpUr_idx <- read.delim(paste(data_dir,"selected_cells_EpithelialUrinary.txt",sep=""), sep = "," , stringsAsFactors=FALSE, header = F)$V1 + 1# cell idxs
  SecC_idx <- read.delim(paste(data_dir,"selected_cells_Secretory.txt",sep=""), sep = ",", stringsAsFactors=FALSE, header = F)$V1 + 1# cell idxs
  logik = json_data$CellType$label_list == "Epithelial";
  
  EpUr = Matrix::rowMeans(E[,EpUr_idx])
  SecC = Matrix::rowMeans(E[,SecC_idx])
  Epis = Matrix::rowMeans(E[,logik])
  
  ## load expression matrix for Fibroblast cell types
  data.dir = "../Data_Signac/AMP_Phase1_RA/"
  spring.dir = "../Data_Signac/AMP_Phase1_RA/FullDataset_v1/"
  E <- CID.LoadData(data.dir = data.dir)

  # read the JSON from SPRING
  json_file = 'categorical_coloring_data.json'
  gJ <- paste(spring.dir,json_file,sep = "")
  json_data <- rjson::fromJSON(file=gJ)
  
  # remove all non-candidate genes
  E = E[rownames(E) %in% candidate_genes,]

  # normalize to these genes
  E = CID.Normalize(E)
  
  # create a subsetted expression matrix for imputation
  logik = json_data$CellTypeFACS$label_list == "Fibroblast"; sum(logik)
  
  Fibs = Matrix::rowMeans(E[,logik])

  LM22 = cbind(LM22, Epis[match(rownames(LM22), names(Epis))])
  LM22 = cbind(LM22, Fibs[match(rownames(LM22), names(Fibs))])
  LM22 = cbind(LM22, EpUr[match(rownames(LM22), names(EpUr))])
  LM22 = cbind(LM22, SecC[match(rownames(LM22), names(SecC))])

  LM22 = na.omit(LM22)

  # get genesets from sc data
  pDC_idx <- read.csv(paste(data_dir,"selected_cells_pDC.txt",sep=""), stringsAsFactors=FALSE, header = F)$V1 + 1 # cell idxs
  Erythro_idx <- read.csv(paste(data_dir,"selected_cells_Erythro.txt",sep=""), stringsAsFactors=FALSE, header = F)$V1 + 1 # cell idxs

  ## load expression matrix for T1D
  data.dir = "../Data_Signac/TS_T1D_Pilot/"
  E <- CID.LoadData(data.dir = data.dir)
  
  # remove all non-candidate genes
  E = E[rownames(E) %in% candidate_genes,]
  
  # normalize to these genes
  E = CID.Normalize(E)

  pDCs = Matrix::rowMeans(E[,pDC_idx])
  Erythro = Matrix::rowMeans(E[,Erythro_idx])

  LM22 = cbind(LM22, pDCs[match(rownames(LM22), names(pDCs))])
  LM22 = cbind(LM22, Erythro[match(rownames(LM22), names(Erythro))])

  LM22 = na.omit(LM22)

  ## load expression matrix for PBMC fresh frozen
  CytoT_indx  <-  read.csv(paste(data_dir,"selected_cells_CytotoxicTcellsFreshFrozen.txt",sep=""), stringsAsFactors=FALSE, header = F)$V1 + 1 # cell idxs
  NKcells_indx <- read.csv(paste(data_dir,"selected_cells_NKcellsFreshFrozen.txt"        ,sep=""), stringsAsFactors=FALSE, header = F)$V1 + 1 # cell idxs
  Platelets_indx <- read.csv(paste(data_dir,"selected_cells_PlateletsFreshFrozen.txt"    ,sep=""), stringsAsFactors=FALSE, header = F)$V1 + 1 # cell idxs
  MCP_indx <-read.csv(paste(data_dir,"selected_cells_MastCellProgenitors.txt"    ,sep=""), stringsAsFactors=FALSE, header = F)$V1 + 1 # cell idxs
  
  data.dir = "../Data_Signac/TS_Fresh_v_FrozenPBMC/"
  E <- CID.LoadData(data.dir = data.dir)
  
  # remove all non-candidate genes
  E = E[rownames(E) %in% candidate_genes,]
  
  # normalize to these genes
  E = CID.Normalize(E)
  
  CytoT = Matrix::rowMeans(E[,CytoT_indx])
  NKcells = Matrix::rowMeans(E[,NKcells_indx])
  Platelets = Matrix::rowMeans(E[,Platelets_indx])
  MCP = Matrix::rowMeans(E[,MCP_indx])
  
  LM22 = cbind(LM22, NKcells[match(rownames(LM22), names(NKcells))])
  LM22 = cbind(LM22, CytoT[match(rownames(LM22), names(CytoT))])
  LM22 = cbind(LM22, Platelets[match(rownames(LM22), names(Platelets))])
  LM22 = cbind(LM22, MCP[match(rownames(LM22), names(MCP))])
  LM22 = na.omit(LM22)

  colnames(LM22) <- c(dummy, "Epithelial", "Fibroblasts","Epithelial.Urinary","Secretory", "pDCs", "Erythro", "NK.cells", "Cyto.T.cells", "Platelets", "Mast.Progenitors");
  colnames(LM22)[colnames(LM22) == "T.cells.regulatory..Tregs."] = "T.cells.regs"

  logik = blueprint_encode$main_types %in% c("CD4+ T-cells", "CD8+ T-cells", "B-cells", "Fibroblasts", "HSC", "NK cells", "Monocytes", "Macrophages", "Erythrocytes","DC") & blueprint_encode$types != "Plasma cells";
  dd = blueprint_encode$data[,logik]
  colnames(dd) <- blueprint_encode$types[logik]
  LM22 = cbind(LM22, dd[match(rownames(LM22), rownames(dd)),])
  LM22 = na.omit(LM22)
  
  LM22 = as.matrix(LM22)

  # quantile normalize within biological replicates
  #for (i in 1:length(unique(colnames(LM22))))
  #{
  #  logik = colnames(LM22) == colnames(LM22)[i];
  #  LM22[,logik] = normalizeQuantiles(LM22[,logik])
  #}

  S = CID.Normalize(LM22)
  
  colnames(S) <- colnames(LM22)
  
  usethis::use_data(S, overwrite = TRUE)

}