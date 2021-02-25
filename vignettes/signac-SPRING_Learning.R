## ----setupSeurat, message = F, eval = F---------------------------------------
#  library(Seurat)
#  library(SignacX)

## ----setup, message = F, eval = F---------------------------------------------
#  # load CITE-seq data
#  data.dir = './CITESEQ_EXPLORATORY_CITESEQ_5K_PBMCS/FullDataset_v1_protein'
#  E = CID.LoadData(data.dir = data.dir)
#  
#  # Load labels
#  json_data = rjson::fromJSON(file=paste0(data.dir,'/categorical_coloring_data.json'))

## ----Seurat, eval = F---------------------------------------------------------
#  # separate protein and gene expression data
#  logik = grepl("Total", rownames(E))
#  P = E[logik,]
#  E = E[!logik,]
#  
#  # CLR normalization in Seurat
#  colnames(P) <- 1:ncol(P)
#  colnames(E) <- 1:ncol(E)
#  reference <- CreateSeuratObject(E)
#  reference[["ADT"]] <- CreateAssayObject(counts = P)
#  reference <- NormalizeData(reference, assay = "ADT", normalization.method = "CLR")

## ----Seurat 2, eval = F-------------------------------------------------------
#  # generate labels
#  lbls = json_data$CellStates$label_list
#  lbls[lbls != "NK"] = "Unclassified"
#  CD16 = reference@assays$ADT@counts[rownames(reference@assays$ADT@counts) == "CD16-TotalSeqB-CD16",]
#  CD56 = reference@assays$ADT@counts[rownames(reference@assays$ADT@counts) == "CD56-TotalSeqB-CD56",]
#  logik = log2(CD56) > 10 & log2(CD16) < 7.5 & lbls == "NK"; sum(logik)
#  lbls[logik] = "NK.CD56bright"

## ----Signac, message = T, eval = F--------------------------------------------
#  # generate bootstrapped single cell data
#  R_learned = SignacBoot(E = E, spring.dir = data.dir, L = c("NK", "NK.CD56bright"), labels = lbls, logfc.threshold = 1)
#  
#  # save the training data
#  save(R_learned, file = "training_NKBright_v207.rda")

## ----Seurat Visualization 0, message = F, eval = F----------------------------
#  # Classify another data set with new model
#  # load new data
#  new.data.dir = "./PBMCs_5k_10X/FullDataset_v1"
#  E = CID.LoadData(data.dir = new.data.dir)
#  # load cell types identified with Signac
#  json_data = rjson::fromJSON(file=paste0(new.data.dir,'/categorical_coloring_data.json'))

## ----Seurat Visualization 1, message = F, eval = F----------------------------
#  # generate new labels
#  cr_learned = Signac(E = E, R = R_learned, spring.dir = new.data.dir)

## ----Seurat Visualization 2, message = F, eval = F----------------------------
#  # modify the existing labels
#  cr = lapply(json_data, function(x) x$label_list)
#  logik = cr$CellStates == 'NK'
#  cr$CellStates[logik] = cr_learned[logik]
#  logik = cr$CellStates_novel == 'NK'
#  cr$CellStates_novel[logik] = cr_learned[logik]
#  new.data.dir = paste0(new.data.dir, "_Learned")

## ----Seurat Visualization 3, message = F, eval = F----------------------------
#  # save
#  dat = CID.writeJSON(cr, spring.dir = new.data.dir, new_colors = c('red'), new_populations = c( 'NK.CD56bright'))

## ---- echo=FALSE--------------------------------------------------------------
sessionInfo()

