## load blueprint_encode.rda data from singleR
data_dir = "./data-raw/"
load(paste(data_dir, "blueprint_encode.rda", sep = ""))

## remove NA values
blueprint_encode$data = na.omit(blueprint_encode$data)

# remove RP + MT genes
logik = grepl( "^MT-",rownames(blueprint_encode$data)) | grepl("^RP[0-9]", rownames(blueprint_encode$data)) | grepl("^RPS", rownames(blueprint_encode$data)) | grepl("^RPL", rownames(blueprint_encode$data))
blueprint_encode$data = blueprint_encode$data[!logik,]

# remove all progenitor cell types in bone marrow
logik = blueprint_encode$types %in%  c("MPP", "GMP", "CMP", "CLP", "MEP", "Megakaryocytes")
blueprint_encode$data = blueprint_encode$data[,!logik]
blueprint_encode$main_types = blueprint_encode$main_types[!logik]
blueprint_encode$types = blueprint_encode$types[!logik]

## create hierarchy;
## top level, immune vs. nonimmune
immune    = blueprint_encode$types
immune[immune %in% c("naive B-cells", "Memory B-cells", "Class-switched memory B-cells", "Plasma cells")] = "Immune"
immune[grepl("T", immune)] = "Immune"
immune[immune == "NK cells"] = "Immune"
immune[immune %in% c("Monocytes", "Macrophages", "Macrophages M1", "Macrophages M2", "DC")] = "Immune"
immune[immune %in% c("Neutrophils", "Eosinophils")] = "Immune"
immune[immune %in% c("Erythrocytes")] = "Immune"
immune[immune %in% c("HSC", "MPP", "GMP", "CMP", "CLP", "MEP", "Megakaryocytes")] = "Immune"
immune[immune %in% c("Adipocytes", "Astrocytes", "Chondrocytes",                  
                     "Endothelial cells", "Epithelial cells", "Fibroblasts",             
                     "Keratinocytes"  ,  "Megakaryocytes", "Melanocytes",
                     "Mesangial cells", "mv Endothelial cells", 
                     "Myocytes", "Neurons", "Pericytes", "Preadipocytes", "Skeletal muscle", "Smooth muscle")] = "NonImmune"

prog = immune
prog[blueprint_encode$types %in% c("HSC", "MPP", "GMP", "CMP", "CLP", "MEP", "Megakaryocytes")] = "HSC"

## penultimate level; lymphocytes vs. myeloid
hema = blueprint_encode$types
hema[hema %in% c("naive B-cells", "Memory B-cells", "Class-switched memory B-cells", "Plasma cells")] = "Lymphocytes"
hema[grepl("T", hema)] = "Lymphocytes"
hema[hema == "NK cells"] = "Lymphocytes"
hema[hema %in% c("Monocytes", "Macrophages", "Macrophages M1", "Macrophages M2", "DC")] = "Myeloid"
hema[hema %in% c("Neutrophils", "Eosinophils")] = "Myeloid"
hema[hema %in% c("Erythrocytes")] = "Myeloid"
hema[hema %in% c("HSC")] = "Myeloid"
hema[hema %in% c("Adipocytes", "Astrocytes", "Chondrocytes", "CLP", "CMP",                  
                     "Endothelial cells", "Epithelial cells", "Fibroblasts", "GMP",             
                     "Keratinocytes"  ,  "Megakaryocytes", "Melanocytes",
                     "MEP"              ,  "Mesangial cells",  "MPP", "mv Endothelial cells", 
                     "Myocytes", "Neurons", "Pericytes", "Preadipocytes", "Skeletal muscle", "Smooth muscle")] = "NonImmune"

## Layer separates TNK and B cells in the lymphocytes
celltypes = hema
celltypes[blueprint_encode$types %in% c("naive B-cells", "Memory B-cells", "Class-switched memory B-cells", "Plasma cells")] = "B.cells"
celltypes[grepl("T", blueprint_encode$types)] = "TNK"
celltypes[blueprint_encode$types == "NK cells"] = "TNK"

## Layer separates MPh and Granulocytes in myeloid
celltypes_L0 = celltypes
celltypes_L0[blueprint_encode$types %in% c("Monocytes", "Macrophages", "Macrophages M1", "Macrophages M2", "DC")] = "MPh"
celltypes_L0[blueprint_encode$types %in% c("Neutrophils", "Eosinophils")] = "Granulocytes"

## Layer separates B cells from Plasma cells
celltypes_L1 = celltypes_L0
celltypes_L1[blueprint_encode$types %in% c("naive B-cells", "Memory B-cells", "Class-switched memory B-cells")] = "B.cells.mem.naive"
celltypes_L1[blueprint_encode$types %in% c("Plasma cells")] = "Plasma.cells"

## Layer separates B memory and B naive
celltypes_L2 = celltypes_L1
celltypes_L2[blueprint_encode$types %in% c("naive B-cells")] = "B.cells.naive"
celltypes_L2[blueprint_encode$types %in% c("Memory B-cells", "Class-switched memory B-cells")] = "B.cells.memory"

# Layer separates NK from T cells
celltypes_L3 = celltypes_L2
celltypes_L3[grepl("T", blueprint_encode$types)] = "T.cells"
celltypes_L3[blueprint_encode$types == "NK cells"] = "NK"

# Layer separates CD4 from CD8 T cells
celltypes_L4 = celltypes_L3
celltypes_L4[blueprint_encode$types %in% c("CD4+ T-cells", "CD4+ Tcm", "CD4+ Tem", "Tregs")] = "T.cells.CD4"
celltypes_L4[blueprint_encode$types %in% c("CD8+ Tcm", "CD8+ Tem", "CD8+ T-cells")] = "T.cells.CD8"

# Tregs from CD4 (note: in hpca this level does gamma delta from T CD4s)
celltypes_L5 = celltypes_L4
celltypes_L5[blueprint_encode$types %in% c("CD4+ T-cells", "CD4+ Tcm", "CD4+ Tem")] = "T.cells.CD4.memory"
celltypes_L5[blueprint_encode$types %in% c("Tregs")] = "T.regs"

# CD4 memory from CD4 
celltypes_L6 = celltypes_L5
celltypes_L6[blueprint_encode$types %in% c("CD4+ Tem")] = "T.cells.CD4.effector.memory"
celltypes_L6[blueprint_encode$types %in% c("CD4+ Tcm")] = "T.cells.CD4.central.memory"

# Monocytes from DCs
celltypes_L7    = celltypes_L6
celltypes_L7[blueprint_encode$types %in% c("DC")] = "DC"
celltypes_L7[blueprint_encode$types %in% c("Monocytes", "Macrophages", "Macrophages M1", "Macrophages M2")] = "Monocytes"

# Granulocytes
celltypes_L8    = celltypes_L7
celltypes_L8[blueprint_encode$types %in% c("Neutrophils", "Eosinophils")] = blueprint_encode$types[blueprint_encode$types %in% c("Neutrophils", "Eosinophils")]

# NonImmune types
celltypes_L9    = immune
celltypes_L9[celltypes_L9 %in% c("NonImmune")] = blueprint_encode$main_types[celltypes_L9 %in% c("NonImmune")]

## create reference annotation set
Reference_blueprint_encode = list(data           = blueprint_encode$data,
                                    celltypes    = data.frame(
                                    immune = as.character(immune),
                                    prog = as.character(prog)
                                    #hema   = as.character(hema),
                                    #L0 = as.character(celltypes)#,
                                   # L1 = as.character(celltypes_L0),
                                  #  L2 = as.character(celltypes_L1),
                                  #  L3 = as.character(celltypes_L2),
                                  #  L4 = as.character(celltypes_L3),
                                  #  L5 = as.character(celltypes_L4),
                                  #  L6 = as.character(celltypes_L5),
                                  #  L7 = as.character(celltypes_L6),
                                  #  L8 = as.character(celltypes_L7),
                                  #  L9 = as.character(celltypes_L8),
                                  #  L10 = as.character(celltypes_L9)
                                  )
                      )
Reference_blueprint_encode$celltypes = apply(Reference_blueprint_encode$celltypes, 2, as.character)

## save reference for internal loading
usethis::use_data(Reference_blueprint_encode, overwrite = TRUE)
