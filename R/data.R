#' Reference expression matrix bootstrapped for cell type / cell state signatures.
#'
#' A reference dataset with 1,000 pure cell types
#'
#' @format A list with two entries
#' \describe{
#'   \item{genes}{genes, identified by HUGO symbols}
#'   \item{Reference}{different training datasets for binary cell type classification}
#'   ...
#' }
"training_HPCA"

#' Genes of interest for drug discovery / disease biology research
#'
#' 3,304 genes curated from external sources
#'
#' @format A data frame with five columns
#' \describe{
#'   \item{Genes}{genes, identified by HUGO symbols}
#'   \item{CellPhoneDB}{0 if not in CellPhoneDB, 1 if gene is listed as a receptor in CellPhoneDB}
#'   \item{GWAS}{0 if not GWAS, 1 if GWAS in GWAS catalog}
#'   \item{PI_drugs}{0 if not in Priority Index paper, 1 if gene is listed as a drug target in priority index paper}
#'   \item{PI_GWAS}{0 if not in Priority Index GWAS list, 1 if gene is listed as GWAS in priority index paper}
#'   ...
#' }
"Genes_Of_Interest"

#' A list of machine learning models to use with SignacFast
#'
#' List of neural network models (1,800) trained using the ModelGenerator function.
#'
#' @format A list with 18 entries, each with 100 neural networks.
#' \describe{
#'   \item{List_element}{Each list element is an ensemble of neural networks}
#'   \item{genes}{genes used in the neural network model}
#'   ...
#' }
"Models_HPCA"

#' A single reference data set used for identifying CD56+ NK cells
#'
#' For use with Signac_Solo
#'
#' @format A data frame with 32 columns (features + cell type labels) and 2,000 bootstrapped cells
#' \describe{
#'   \item{celltypes}{Labels for classification (celltypes)}
#'   ...
#' }
"R_learned"