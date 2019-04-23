#' Reference expression matrix for cell type / cell state signatures.
#'
#' A reference dataset with 9555 genes and 224 experiments (mostly sorted bulk RNA seq)
#'
#' @format A matrix with 9555 rows and 224 columns
#' \describe{
#'   \item{genes}{genes, identified by HUGO symbols}
#'   \item{sorted cells}{identified by common nomenclature}
#'   ...
#' }
"S"

#' The gene signatures used for cell state identification.
#'
#' A group of 3425 markers used for cell state identification.
#'
#' @format A named list (length 10), each containing a dataframe with three columns and N rows.
#' \describe{
#'   \item{HUGO symbols}{N genes, each identified by a HUGO symbol}
#'   \item{Cell population}{identified by common nomenclature}
#'   \item{Polarity}{+ for positive markers, - for negative markers}
#'   ...
#' }
"cellstate_markers"

#' The gene signatures used for cell type identification.
#'
#' A group of 2715 markers used for cell type identification.
#'
#' @format A data.frame with 2715 rows and three columns.
#' \describe{
#'   \item{HUGO symbols}{genes, identified by HUGO symbols}
#'   \item{Cell population}{identified by common nomenclature}
#'   \item{Polarity}{+ for positive markers, - for negative markers}
#'   ...
#' }
"markers"