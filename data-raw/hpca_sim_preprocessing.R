source('/site/ne/data/bh-results/C/CHAMBERLAIN.Mat/manuscripts/Vidi/Vidi/data-raw/GetMarkers.R')
source('/site/ne/data/bh-results/C/CHAMBERLAIN.Mat/manuscripts/Vidi/code/xCell.R')
## load hpca.rda data from singleR, LM22 data from cibersort
data_dir = "./data-raw/"
load(paste(data_dir, "hpca.rda", sep = ""))

# remove RP + MT genes
logik = grepl( "^MT-",rownames(hpca$data)) | grepl("^RP[0-9]", rownames(hpca$data)) | grepl("^RPS", rownames(hpca$data)) | grepl("^RPL", rownames(hpca$data))
hpca$data = hpca$data[!logik,]

# remove all progenitor cell types in bone marrow
logik = hpca$main_types %in% c("HSC_CD34+","HSC_-G-CSF","BM", "BM & Prog.", "CMP", "Embryonic_stem_cells", "Erythroblast", "GMP", "iPS_cells", "MEP", "MSC", "Myelocyte", "Pre-B_cell_CD34-", "Pro-B_cell_CD34+", "Pro-Myelocyte", "Tissue_stem_cells")
hpca$data = hpca$data[,!logik]
hpca$main_types = hpca$main_types[!logik]
hpca$types = hpca$types[!logik]

# switch effector and memory CD8 cell labels
logik = hpca$types %in% "T_cell:CD8+_effector_memory"
logik2 = hpca$types == "T_cell:CD8+_Central_memory"
hpca$types[logik] =  "T_cell:CD8+_Central_memory"
hpca$types[logik2] =  "T_cell:CD8+_effector_memory"

## top level; immune vs. non-immune
immune = hpca$types
immune[immune %in% c("B_cell", "B_cell:CXCR4-_centrocyte", "B_cell:CXCR4+_centroblast",                 
                           "B_cell:Germinal_center", "B_cell:immature", "B_cell:Memory",
                           "B_cell:Naive", "B_cell:Plasma_cell")] = "Immune"
immune[grepl("T_cell", immune)] = "Immune"
immune[immune %in% "Platelets"] = "Immune"
immune[immune %in%       c("NK_cell", "NK_cell:CD56hiCD62L+","NK_cell:IL2")] = "Immune"
immune[immune %in%       c("DC:monocyte-derived", "DC:monocyte-derived:A._fumigatus_germ_tubes_6h",        
                           "DC:monocyte-derived:AEC-conditioned", "DC:monocyte-derived:AM580",                            
                           "DC:monocyte-derived:anti-DC-SIGN_2h", "DC:monocyte-derived:antiCD40/VAF347",                  
                           "DC:monocyte-derived:CD40L", "DC:monocyte-derived:Galectin-1",                  
                           "DC:monocyte-derived:immature", "DC:monocyte-derived:LPS",                       
                           "DC:monocyte-derived:mature" ,  "DC:monocyte-derived:Poly(IC)",                      
                           "DC:monocyte-derived:rosiglitazone", "DC:monocyte-derived:rosiglitazone/AGN193109" ,          
                           "DC:monocyte-derived:Schuler_treatment","DC:monocyte-derived", "DC:monocyte-derived:A._fumigatus_germ_tubes_6h", 
                           "DC:monocyte-derived:AEC-conditioned","DC:monocyte-derived:AM580" ,"DC:monocyte-derived:anti-DC-SIGN_2h",
                           "DC:monocyte-derived:antiCD40/VAF347", "DC:monocyte-derived:CD40L","DC:monocyte-derived:Galectin-1",
                           "DC:monocyte-derived:immature","DC:monocyte-derived:LPS", "DC:monocyte-derived:mature","DC:monocyte-derived:Poly(IC)", "DC:monocyte-derived:rosiglitazone",            
                           "DC:monocyte-derived:rosiglitazone/AGN193109", "DC:monocyte-derived:Schuler_treatment", "Macrophage:Alveolar",
                           "Macrophage:Alveolar:B._anthacis_spores", "Macrophage:monocyte-derived", "Macrophage:monocyte-derived:IFNa",
                           "Macrophage:monocyte-derived:IL-4/cntrl", "Macrophage:monocyte-derived:IL-4/Dex/cntrl", "Macrophage:monocyte-derived:IL-4/Dex/TGFb",
                           "Macrophage:monocyte-derived:IL-4/TGFb","Macrophage:monocyte-derived:M-CSF" ,"Macrophage:monocyte-derived:M-CSF/IFNg","Macrophage:monocyte-derived:M-CSF/IFNg/Pam3Cys",
                           "Macrophage:monocyte-derived:M-CSF/Pam3Cys" , "Macrophage:monocyte-derived:S._aureus", "Monocyte", "Monocyte:anti-FcgRIIB", "Monocyte:CD14+", "Monocyte:CD16-",                            
                           "Monocyte:CD16+", "Monocyte:CXCL4", "Monocyte:F._tularensis_novicida" , "Monocyte:leukotriene_D4" ,  "Monocyte:MCSF" ,"Monocyte:S._typhimurium_flagellin")] = "Immune"
immune[immune %in% c("Neutrophil", "Neutrophil:commensal_E._coli_MG1655", "Neutrophil:GM-CSF_IFNg", "Neutrophil:inflam", "Neutrophil:LPS", "Neutrophil:uropathogenic_E._coli_UTI89" )] = "Immune"
immune[immune %in% c("Erythroblast", "MEP","HSC_-G-CSF","BM", "BM & Prog.", "CMP", "Embryonic_stem_cells", "Erythroblast", "GMP", "iPS_cells", "MEP", "MSC", "Myelocyte", "Pre-B_cell_CD34-", "Pro-B_cell_CD34+", "Pro-Myelocyte", "Tissue_stem_cells","HSC_CD34+")] = "Immune"
immune[! immune %in% c("Immune")] = "NonImmune"

## penultimate level; HSC and Not.HSC
#prog = immune
#prog[hpca$types %in% c("HSC_CD34+")] = "HSC"
#prog[prog == "Immune"] = "Not.HSC"

## lymphocytes and myeloid
hema = immune
hema[hpca$types %in% c("B_cell", "B_cell:CXCR4-_centrocyte", "B_cell:CXCR4+_centroblast",                 
                     "B_cell:Germinal_center", "B_cell:immature", "B_cell:Memory",
                     "B_cell:Naive", "B_cell:Plasma_cell")] = "Lymphocytes"
hema[grepl("T_cell", hpca$types)] = "Lymphocytes"
hema[hpca$types %in% "Platelets"] = "Lymphocytes"
hema[hpca$types %in%       c("NK_cell", "NK_cell:CD56hiCD62L+","NK_cell:IL2")] = "Lymphocytes"
hema[hpca$types %in%       c("DC:monocyte-derived", "DC:monocyte-derived:A._fumigatus_germ_tubes_6h",        
                           "DC:monocyte-derived:AEC-conditioned", "DC:monocyte-derived:AM580",                            
                           "DC:monocyte-derived:anti-DC-SIGN_2h", "DC:monocyte-derived:antiCD40/VAF347",                  
                           "DC:monocyte-derived:CD40L", "DC:monocyte-derived:Galectin-1",                  
                           "DC:monocyte-derived:immature", "DC:monocyte-derived:LPS",                       
                           "DC:monocyte-derived:mature" ,  "DC:monocyte-derived:Poly(IC)",                      
                           "DC:monocyte-derived:rosiglitazone", "DC:monocyte-derived:rosiglitazone/AGN193109" ,          
                           "DC:monocyte-derived:Schuler_treatment","DC:monocyte-derived", "DC:monocyte-derived:A._fumigatus_germ_tubes_6h", 
                           "DC:monocyte-derived:AEC-conditioned","DC:monocyte-derived:AM580" ,"DC:monocyte-derived:anti-DC-SIGN_2h",
                           "DC:monocyte-derived:antiCD40/VAF347", "DC:monocyte-derived:CD40L","DC:monocyte-derived:Galectin-1",
                           "DC:monocyte-derived:immature","DC:monocyte-derived:LPS", "DC:monocyte-derived:mature","DC:monocyte-derived:Poly(IC)", "DC:monocyte-derived:rosiglitazone",            
                           "DC:monocyte-derived:rosiglitazone/AGN193109", "DC:monocyte-derived:Schuler_treatment", "Macrophage:Alveolar",
                           "Macrophage:Alveolar:B._anthacis_spores", "Macrophage:monocyte-derived", "Macrophage:monocyte-derived:IFNa",
                           "Macrophage:monocyte-derived:IL-4/cntrl", "Macrophage:monocyte-derived:IL-4/Dex/cntrl", "Macrophage:monocyte-derived:IL-4/Dex/TGFb",
                           "Macrophage:monocyte-derived:IL-4/TGFb","Macrophage:monocyte-derived:M-CSF" ,"Macrophage:monocyte-derived:M-CSF/IFNg","Macrophage:monocyte-derived:M-CSF/IFNg/Pam3Cys",
                           "Macrophage:monocyte-derived:M-CSF/Pam3Cys" , "Macrophage:monocyte-derived:S._aureus", "Monocyte", "Monocyte:anti-FcgRIIB", "Monocyte:CD14+", "Monocyte:CD16-",                            
                           "Monocyte:CD16+", "Monocyte:CXCL4", "Monocyte:F._tularensis_novicida" , "Monocyte:leukotriene_D4" ,  "Monocyte:MCSF" ,"Monocyte:S._typhimurium_flagellin")] = "Myeloid"
hema[hpca$types %in% c("Neutrophil", "Neutrophil:commensal_E._coli_MG1655", "Neutrophil:GM-CSF_IFNg", "Neutrophil:inflam", "Neutrophil:LPS", "Neutrophil:uropathogenic_E._coli_UTI89" )] = "Myeloid"

## Layer separates TNK and B cells
celltypes_L0 = hema
celltypes_L0[hpca$types %in% c("NK_cell", "NK_cell:CD56hiCD62L+","NK_cell:IL2")] = "TNK"
celltypes_L0[grepl("T_cell", hpca$types)] = "TNK"
celltypes_L0[hpca$types %in% c("B_cell", "B_cell:CXCR4-_centrocyte", "B_cell:CXCR4+_centroblast",                 
                     "B_cell:Germinal_center", "B_cell:immature", "B_cell:Memory",
                     "B_cell:Naive", "B_cell:Plasma_cell")] = "B.cells"

## Layer separates MPh and DC
celltypes_L1 = celltypes_L0
celltypes_L1[hpca$types %in%       c("DC:monocyte-derived", "DC:monocyte-derived:A._fumigatus_germ_tubes_6h",        
                             "DC:monocyte-derived:AEC-conditioned", "DC:monocyte-derived:AM580",                            
                             "DC:monocyte-derived:anti-DC-SIGN_2h", "DC:monocyte-derived:antiCD40/VAF347",                  
                             "DC:monocyte-derived:CD40L", "DC:monocyte-derived:Galectin-1",                  
                             "DC:monocyte-derived:immature", "DC:monocyte-derived:LPS",                       
                             "DC:monocyte-derived:mature" ,  "DC:monocyte-derived:Poly(IC)",                      
                             "DC:monocyte-derived:rosiglitazone", "DC:monocyte-derived:rosiglitazone/AGN193109" ,          
                             "DC:monocyte-derived:Schuler_treatment","DC:monocyte-derived", "DC:monocyte-derived:A._fumigatus_germ_tubes_6h", 
                             "DC:monocyte-derived:AEC-conditioned","DC:monocyte-derived:AM580" ,"DC:monocyte-derived:anti-DC-SIGN_2h",
                             "DC:monocyte-derived:antiCD40/VAF347", "DC:monocyte-derived:CD40L","DC:monocyte-derived:Galectin-1",
                             "DC:monocyte-derived:immature","DC:monocyte-derived:LPS", "DC:monocyte-derived:mature","DC:monocyte-derived:Poly(IC)", "DC:monocyte-derived:rosiglitazone",            
                             "DC:monocyte-derived:rosiglitazone/AGN193109", "DC:monocyte-derived:Schuler_treatment", "Macrophage:Alveolar",
                             "Macrophage:Alveolar:B._anthacis_spores", "Macrophage:monocyte-derived", "Macrophage:monocyte-derived:IFNa",
                             "Macrophage:monocyte-derived:IL-4/cntrl", "Macrophage:monocyte-derived:IL-4/Dex/cntrl", "Macrophage:monocyte-derived:IL-4/Dex/TGFb",
                             "Macrophage:monocyte-derived:IL-4/TGFb","Macrophage:monocyte-derived:M-CSF" ,"Macrophage:monocyte-derived:M-CSF/IFNg","Macrophage:monocyte-derived:M-CSF/IFNg/Pam3Cys",
                             "Macrophage:monocyte-derived:M-CSF/Pam3Cys" , "Macrophage:monocyte-derived:S._aureus", "Monocyte", "Monocyte:anti-FcgRIIB", "Monocyte:CD14+", "Monocyte:CD16-",                            
                             "Monocyte:CD16+", "Monocyte:CXCL4", "Monocyte:F._tularensis_novicida" , "Monocyte:leukotriene_D4" ,  "Monocyte:MCSF" ,"Monocyte:S._typhimurium_flagellin")] = "MPh"
celltypes_L1[hpca$types %in% c("Neutrophil", "Neutrophil:commensal_E._coli_MG1655", "Neutrophil:GM-CSF_IFNg", "Neutrophil:inflam", "Neutrophil:LPS", "Neutrophil:uropathogenic_E._coli_UTI89" )] = "MPh"
celltypes_L1[hpca$types %in%       c("DC:monocyte-derived", "DC:monocyte-derived:A._fumigatus_germ_tubes_6h",        
                                     "DC:monocyte-derived:AEC-conditioned", "DC:monocyte-derived:AM580",                            
                                     "DC:monocyte-derived:anti-DC-SIGN_2h", "DC:monocyte-derived:antiCD40/VAF347",                  
                                     "DC:monocyte-derived:CD40L", "DC:monocyte-derived:Galectin-1",                  
                                     "DC:monocyte-derived:immature", "DC:monocyte-derived:LPS",                       
                                     "DC:monocyte-derived:mature" ,  "DC:monocyte-derived:Poly(IC)",                      
                                     "DC:monocyte-derived:rosiglitazone", "DC:monocyte-derived:rosiglitazone/AGN193109" ,          
                                     "DC:monocyte-derived:Schuler_treatment","DC:monocyte-derived", "DC:monocyte-derived:A._fumigatus_germ_tubes_6h", 
                                     "DC:monocyte-derived:AEC-conditioned","DC:monocyte-derived:AM580" ,"DC:monocyte-derived:anti-DC-SIGN_2h",
                                     "DC:monocyte-derived:antiCD40/VAF347", "DC:monocyte-derived:CD40L","DC:monocyte-derived:Galectin-1",
                                     "DC:monocyte-derived:immature","DC:monocyte-derived:LPS", "DC:monocyte-derived:mature","DC:monocyte-derived:Poly(IC)", "DC:monocyte-derived:rosiglitazone",            
                                     "DC:monocyte-derived:rosiglitazone/AGN193109", "DC:monocyte-derived:Schuler_treatment")] = "DC"

## Layer separates B cells from Plasma cells
celltypes_L2 = celltypes_L1
celltypes_L2[hpca$types %in% c("B_cell:Plasma_cell")] = "Plasma.cells"
celltypes_L2[hpca$types %in% c("B_cell", "B_cell:Memory", "B_cell:Naive")] = "B.cells.NoPlasma"

## Layer separates B memory and B naive
celltypes_L3 = celltypes_L2
celltypes_L3[hpca$types %in% c("B_cell:Naive")] = "B.cells.naive"
celltypes_L3[hpca$types %in% c("B_cell:Memory")] = "B.cells.memory"

# Layer separates NK from T cells
celltypes_L4 = celltypes_L3
celltypes_L4[grepl("T_cell", hpca$types)] = "T.cells"
celltypes_L4[hpca$types %in% c("NK_cell", "NK_cell:CD56hiCD62L+","NK_cell:IL2")] = "NK"

# Layer separates CD4 from CD8 T cells
celltypes_L5 = celltypes_L4
celltypes_L5[hpca$types %in% c("T_cell:CD4+", "T_cell:CD4+_central_memory", "T_cell:CD4+_effector_memory", "T_cell:CD4+_Naive",
                                           "T_cell:Treg:Naive")] = "T.cells.CD4"
celltypes_L5[hpca$types %in% c("T_cell:gamma-delta","T_cell:CD8+", "T_cell:CD8+_Central_memory", "T_cell:CD8+_effector_memory", "T_cell:CD8+_effector_memory_RA", "T_cell:CD8+_naive")] = "T.cells.CD8"

# CD4 memory from CD4 naive
celltypes_L6 = celltypes_L5
celltypes_L6[hpca$types %in% c("T_cell:CD4+_central_memory", "T_cell:CD4+_effector_memory", "T_cell:Treg:Naive")] = "T.cells.CD4.memory.regs"
celltypes_L6[hpca$types %in% c("T_cell:CD4+_Naive")] = "T.cells.CD4.naive"

# T regs from CD4 memory
celltypes_L7    = celltypes_L6
celltypes_L7[hpca$types %in% c("T_cell:CD4+_effector_memory")] = "T.cells.CD4.memory"
celltypes_L7[hpca$types %in% c("T_cell:Treg:Naive")] = "T.regs"

# CD8 naive from memory
celltypes_L8    = celltypes_L7
celltypes_L8[hpca$types %in% c("T_cell:CD8+_naive")] = "T.cells.CD8.naive"
celltypes_L8[hpca$types %in% c("T_cell:CD8+_effector_memory", "T_cell:CD8+_effector_memory_RA")] = "T.cells.CD8.memory"

# Monocytes from macrophages
celltypes_L9    = celltypes_L8
celltypes_L9[hpca$types %in% c( "Monocyte", "Monocyte:anti-FcgRIIB", "Monocyte:CD14+", "Monocyte:CD16-",                            
                                "Monocyte:CD16+", "Monocyte:CXCL4", "Monocyte:F._tularensis_novicida" , "Monocyte:leukotriene_D4" ,  "Monocyte:MCSF" ,"Monocyte:S._typhimurium_flagellin")] = "Monocytes.Neutrophils"
celltypes_L9[hpca$types %in% c("Macrophage:Alveolar",
                               "Macrophage:Alveolar:B._anthacis_spores", "Macrophage:monocyte-derived", "Macrophage:monocyte-derived:IFNa",
                               "Macrophage:monocyte-derived:IL-4/cntrl", "Macrophage:monocyte-derived:IL-4/Dex/cntrl", "Macrophage:monocyte-derived:IL-4/Dex/TGFb",
                               "Macrophage:monocyte-derived:IL-4/TGFb","Macrophage:monocyte-derived:M-CSF" ,"Macrophage:monocyte-derived:M-CSF/IFNg","Macrophage:monocyte-derived:M-CSF/IFNg/Pam3Cys",
                               "Macrophage:monocyte-derived:M-CSF/Pam3Cys" , "Macrophage:monocyte-derived:S._aureus")] = "Macrophages"
celltypes_L9[hpca$types %in% c("Neutrophil", "Neutrophil:commensal_E._coli_MG1655", "Neutrophil:GM-CSF_IFNg", "Neutrophil:inflam", "Neutrophil:LPS", "Neutrophil:uropathogenic_E._coli_UTI89" )] = "Monocytes.Neutrophils"

# Neutrophils from monocytes
celltypes_L10    = celltypes_L9
celltypes_L10[hpca$types %in% c("Monocyte", "Monocyte:anti-FcgRIIB", "Monocyte:CD14+", "Monocyte:CD16-",                            
                               "Monocyte:CD16+", "Monocyte:CXCL4", "Monocyte:F._tularensis_novicida" , "Monocyte:leukotriene_D4" ,  "Monocyte:MCSF" ,"Monocyte:S._typhimurium_flagellin")] = "Monocytes"
celltypes_L10[hpca$types %in% c("Neutrophil", "Neutrophil:commensal_E._coli_MG1655", "Neutrophil:GM-CSF_IFNg", "Neutrophil:inflam", "Neutrophil:LPS", "Neutrophil:uropathogenic_E._coli_UTI89" )] = "Neutrophils"

# classical and non-classical monocytes
celltypes_L11    = celltypes_L10
celltypes_L11[hpca$types %in% "Monocyte:CD16-" ] = "Mon.Classical"
celltypes_L11[hpca$types %in% "Monocyte:CD16+" ] = "Mon.NonClassical"

# NonImmune types
celltypes_L12    = celltypes_L11
celltypes_L12[celltypes_L11 %in% c("NonImmune")] = "Non.Fibroblasts"
celltypes_L12[hpca$main_types %in% c("Fibroblasts")] = "Fibroblasts"

celltypes_L13    = celltypes_L12
celltypes_L13[celltypes_L12 %in% c("Non.Fibroblasts")] = "Non.Epithelial"
celltypes_L13[hpca$main_types %in% c("Epithelial_cells")] = "Epithelial.cells"

celltypes_L14    = celltypes_L13
celltypes_L14[celltypes_L13 %in% c("Non.Epithelial")] = "NonImmune"
celltypes_L14[hpca$main_types %in% c("Endothelial_cells")] = "Endothelial.cells"

celltypes_L15 = celltypes_L14
celltypes_L15[hpca$types %in% c("T_cell:CD8+_effector_memory", "T_cell:CD8+_effector_memory_RA")] = "T.cells.CD8.em"
celltypes_L15[hpca$types %in% "T_cell:CD8+_Central_memory"] = "T.cells.CD8.cm"

# get xCell genes
#gns = getxCellGenes(dataset = "HPCA")
#plasma_markers = sort(unique(c(gns$`Plasma cells`)))

## create reference annotation set
gns = getxCellGenes()
Reference_hpca = list(data         =hpca$data[rownames(hpca$data) %in% xCell.data$genes,],
                                    celltypes    = data.frame(
                                    immune = as.character(immune),
                                   # prog = as.character(prog),
                                    hema   = as.character(hema),
                                    L0 = as.character(celltypes_L0),
                                    L1 = as.character(celltypes_L1),
                                    L2 = as.character(celltypes_L2),
                                    L3 = as.character(celltypes_L3),
                                    L4 = as.character(celltypes_L4),
                                    L5 = as.character(celltypes_L5),
                                    L6 = as.character(celltypes_L6),
                                    L7 = as.character(celltypes_L7),
                                    L8 = as.character(celltypes_L8),
                                    L9 = as.character(celltypes_L9),
                                    L10 = as.character(celltypes_L10),
                                    L11 = as.character(celltypes_L11),
                                    L12 = as.character(celltypes_L12),
                                    L13 = as.character(celltypes_L13),
                                    L14 = as.character(celltypes_L14),
                                    L15 = as.character(celltypes_L15)
                                  ))
Reference_hpca$celltypes = apply(Reference_hpca$celltypes, 2, as.character)

#mrks[[1]] = rbind(mrks[[1]], df)
# get xCell genes

#mrks$prog = list("")
#mrks$prog$genes = gns$HSC
#mrks$prog$cluster = "HSC"
#Reference_hpca$features = mrks
# run correlation heatmap code

# use xCell marker genes
#Reference_hpca$data = Reference_hpca$data[rownames(Reference_hpca$data) %in% gns,]

# get markers based on annotation for every level
mrks = GetMarkers(Reference_hpca)
names(mrks) <- colnames(Reference_hpca$celltypes)
#mrks = lapply(mrks, function(x) {x[x$avg_logFC > 0.5, ]})
#logik = sapply(mrks, function(x) {nrow(x) > 0})
#mrks = mrks[logik]
#Reference_hpca$celltypes = Reference_hpca$celltypes[,logik]
#df = data.frame(p_val = NA, avg_logFC = NA, pct.1 = 1, pct.2 = 1, p_val_adj = NA, cluster = "T.cells.CD4", gene = "IL6R")
#mrks$L5 = rbind(mrks$L5, df)
#mrks = lapply(mrks, function(x) {x[x$avg_logFC > 0.5, ]})

# Reference_hpca$data = hpca$data

gns = getxCellGenes()
Bcellmemory = sort(toupper(hpca$de.genes$`B_cell:Memory`$`B_cell:Naive`))
Bcellnaive = sort(toupper(hpca$de.genes$`B_cell:Naive`$`B_cell:Memory`))
Bcellmemory = Bcellmemory[Bcellmemory %in% rownames(Reference_hpca$data)]
Bcellnaive = Bcellnaive[Bcellnaive %in% rownames(Reference_hpca$data)]
df = data.frame(p_val = NA, avg_logFC = NA, pct.1 = 1, pct.2 = 1, p_val_adj = NA, cluster = c(rep("B.cells.memory", length(Bcellmemory)), rep("B.cells.naive", length(Bcellnaive))), gene = c(Bcellmemory, Bcellnaive))
mrks$L3 = df
#pca <- prcomp(x = t(Reference_hpca$data[rownames(Reference_hpca$data) %in% df$gene,]), center = T, scale. = T)
#df = data.frame(pca$x, celltypes = Reference_hpca$celltypes[,7])
#autoplot(pca, data = df, colour = 'celltypes')

plasma = sort(toupper(hpca$de.genes$`B_cell:Plasma_cell`$B_cell))
non.plasma = sort(toupper(hpca$de.genes$B_cell$`B_cell:Plasma_cell`))
plasma = plasma[plasma %in% rownames(Reference_hpca$data)]
non.plasma = non.plasma[non.plasma %in% rownames(Reference_hpca$data)]
df = data.frame(p_val = NA, avg_logFC = NA, pct.1 = 1, pct.2 = 1, p_val_adj = NA, cluster = c(rep("Plasma.cells", length(plasma)), rep("B.cells.NoPlasma", length(non.plasma))), gene = c(plasma, non.plasma))
mrks$L2 = df

cd4memory = sort(toupper(hpca$de.genes$`T_cell:CD4+_central_memory`$`T_cell:CD4+_Naive`)[1:20])
cd4naive = sort(toupper(hpca$de.genes$`T_cell:CD4+_Naive`$`T_cell:CD4+_central_memory`)[1:20])
cd4memory = cd4memory[cd4memory %in% rownames(Reference_hpca$data)]
cd4naive = cd4naive[cd4naive %in% rownames(Reference_hpca$data)]
df = data.frame(p_val = NA, avg_logFC = NA, pct.1 = 1, pct.2 = 1, p_val_adj = NA, cluster = c(rep("T.cells.CD4.memory.regs", length(cd4memory)), rep("T.cells.CD4.naive", length(cd4naive))), gene = c(cd4memory, cd4naive))
mrks$L6 = rbind(mrks$L6[mrks$L6$avg_logFC > 0.5,], df)

#Tregs = toupper(hpca$de.genes$`T_cell:Treg:Naive`$`T_cell:CD4+_effector_memory`)
Tregs = unique(unlist(gns[names(gns) == "Tregs"]))
cd4memory = toupper(hpca$de.genes$`T_cell:CD4+_effector_memory`$`T_cell:Treg:Naive`)[1:40]
cd4memory = cd4memory[cd4memory %in% rownames(Reference_hpca$data)]
Tregs = Tregs[Tregs %in% rownames(Reference_hpca$data)]
df = data.frame(p_val = NA, avg_logFC = NA, pct.1 = 1, pct.2 = 1, p_val_adj = NA, cluster = c(rep("T.regs", length(Tregs)), rep("T.cells.CD4.memory", length(cd4memory))), gene = c(Tregs, cd4memory))
mrks$L7 = df

#logik = names(gns) %in% c("CD8+ Tem")
#T.cells.CD8.memory =  unique(unlist(gns[which(logik)]))
T.cells.CD8.memory = unique(toupper(c ( hpca$de.genes$`T_cell:CD8+_naive`$`T_cell:CD8+_effector_memory`[1:30], hpca$de.genes$`T_cell:CD8+_effector_memory_RA`$`T_cell:CD8+_naive`[1:30])))
T.cells.CD8.naive = unique(toupper(c ( hpca$de.genes$`T_cell:CD8+_naive`$`T_cell:CD8+_Central_memory`[1:80])))
T.cells.CD8.naive = T.cells.CD8.naive[T.cells.CD8.naive %in% rownames(hpca$data)]
T.cells.CD8.memory = T.cells.CD8.memory[T.cells.CD8.memory %in% rownames(hpca$data)]
df = data.frame(p_val = NA, avg_logFC = NA, pct.1 = 1, pct.2 = 1, p_val_adj = NA, cluster = c(rep("T.cells.CD8.memory", length(T.cells.CD8.memory)), rep("T.cells.CD8.naive", length(T.cells.CD8.naive))), gene = c(T.cells.CD8.memory,T.cells.CD8.naive))
mrks$L8 = df

T.cells.CD8.em = unique(toupper( c( hpca$de.genes$`T_cell:CD8+_Central_memory`$`T_cell:CD8+_effector_memory`)))
T.cells.CD8.cm = unique(toupper( c( hpca$de.genes$`T_cell:CD8+_effector_memory`$`T_cell:CD8+_Central_memory`)))
T.cells.CD8.cm = T.cells.CD8.cm[T.cells.CD8.cm %in% rownames(hpca$data)]
T.cells.CD8.em = T.cells.CD8.em[T.cells.CD8.em %in% rownames(hpca$data)]
df = data.frame(p_val = NA, avg_logFC = NA, pct.1 = 1, pct.2 = 1, p_val_adj = NA, cluster = c(rep("T.cells.CD8.cm", length(T.cells.CD8.cm)), rep("T.cells.CD8.em", length(T.cells.CD8.em))), gene = c(T.cells.CD8.cm,T.cells.CD8.em))
mrks$L15 = df

# bootstrap data
#xx = rownames(Reference_hpca$data)
#Reference_hpca$data = preprocessCore::normalize.quantiles(Reference_hpca$data)
#rownames(Reference_hpca$data) <- xx
boots = list("")
size = 1000
set.seed('42')
 for (i in 1:length(mrks)) {
    x = mrks[[i]]
    gns = x$gene
    dat = Reference_hpca$data[rownames(Reference_hpca$data) %in% gns,]
    cts = split.data.frame(x, f = x$cluster)
    N = lapply(cts, function(x){
      # first sample from cells in cluster 1, size cells
      logik = rownames(dat) %in% x$gene
      dummy = dat[logik,]
      logik = Reference_hpca$celltypes[,i] == as.character(x$cluster[1])
      dummy = dummy[,logik]
      dd = t(apply(dummy, 1, function(z) {
        sample(z, size = size, replace = T)}))
      dd = t(apply(dd, 1, function(z){
        rnorm(n = length(z), mean = mean(z), sd = sd(z))
      }))
      # now sample from cells in cluster 2, size cells with True Negative Expr.
      logik = !rownames(dat) %in% x$gene
      dummy = dat[logik,]
      logik = Reference_hpca$celltypes[,i] == as.character(x$cluster[1])
      dummy = dummy[,logik]
      dd2 = t(apply(dummy, 1, function(z) {
        sample(z, size = size, replace = T)}))
      dd2 = t(apply(dd2, 1, function(z){
        rnorm(n = length(z), mean = mean(z), sd = sd(z))
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
    boots[[i]] = data.frame(N2, celltypes = c(rep(names(cts)[1], size), rep(names(cts)[2], size)))
 }

# i = 12
#pca <- prcomp(x = boots[[i]][,-ncol(boots[[i]])], center = T, scale. = T) 
#autoplot(pca, data = boots[[i]], colour = 'celltypes')

#i = 1
#E = t(boots[[i]][,-ncol(boots[[i]])])
#lbls = boots[[i]]$celltypes 
#logik = celltypes_L14 %in% c("T.cells.CD8.memory", "T.cells.CD8.naive")
#E = hpca$data[,logik]
#E = hpca$data
#lbls = celltypes_L14[logik]
#pbmc <- Seurat::CreateSeuratObject(counts = E[rownames(E) %in% mrks[[i]]$gene,], project = "HPCA", min.cells = 0)
#pbmc <- Seurat::CreateSeuratObject(counts = E, project = "HPCA", min.cells = 0)
#pbmc <- Seurat::NormalizeData(object = pbmc, verbose = F)
#pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
#pbmc <- ScaleData(pbmc)
#pbmc <- RunPCA(pbmc)
#pbmc <- RunUMAP(pbmc, dims = 1:10)
#lbls = Reference_hpca$celltypes[,i]
#pbmc <- Seurat::AddMetaData(pbmc, metadata=lbls, col.name = "celltypes")
#pbmc <- Seurat::SetIdent(pbmc, value='celltypes')
#DimPlot(pbmc, reduction = "umap")

#pca <- prcomp(x = t(Reference_hpca$data[rownames(Reference_hpca$data) %in% mrks$immune$gene,]), center = T, scale. = T)
#df = data.frame(pca$x, celltypes = Reference_hpca$celltypes[,1])
#autoplot(pca, data = df, colour = 'celltypes')

Reference_sim = list("")
Reference_sim$Reference = boots
sapply(mrks, function(x) {unique(x$cluster)})
names(Reference_sim$Reference) = c("All", "Immune", "Lymphocytes", "Myeloid", "B.cells", "B.cells.NoPlasma",
                                   "TNK", "T.cells", "T.cells.CD4", "T.cells.CD4.memory.regs", "T.cells.CD8", "MPh", 
                                   "Monocytes.Neutrophils", "Monocytes", "NonImmune", "Non.Fibroblasts", "Non.Epithelial", "T.cells.CD8.memory")
gns = sort(unique(unlist(sapply(mrks, function(x) {as.character(x$gene)}))))
Reference_sim$genes = gns
Reference_sim = Reference_sim[-1]
Reference_sim_v6 = Reference_sim
## save reference for internal loading
usethis::use_data(Reference_sim_v6, overwrite = TRUE)

## run logistic regressions
#D = lapply(Reference_sim_v5$Reference, function(x){
#  fits = apply(x[,-ncol(x)], 2, function(y){
#    fit = lm(y ~ x$celltypes)
#    as.numeric(summary(fit)$coefficients[,4][2])
#  })
#})

### define references by each level
All = Reference_sim_v6$Reference[[1]]
