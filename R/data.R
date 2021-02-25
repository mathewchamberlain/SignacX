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