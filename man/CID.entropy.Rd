% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{CID.entropy}
\alias{CID.entropy}
\title{Normalized Shannon entropy-based "unclassified" assignment}
\usage{
CID.entropy(ac, distM)
}
\arguments{
\item{ac}{a character vector of cell type labels}

\item{distM}{the distance matrix, see ?CID.GetDistMat}
}
\value{
A character vector like 'ac' but with cells type labels set to "Unclassified" if there was high normalized Shannon entropy.
}
\description{
\code{CID.entropy} calculates the normalized Shannon entropy of labels for each cell
among k-nearest neighbors less than four-degrees apart, and then sets cells with statistically
significant large Shannon entropy to be "Unclassified."
}
\examples{
\dontrun{
# load data classified previously (see \code{SignacFast})
P <- readRDS("celltypes.rds")
S <- readRDS("pbmcs.rds")

# get edges from default assay from Seurat object
default.assay <- Seurat::DefaultAssay(S)
edges = S@graphs[[which(grepl(paste0(default.assay, "_nn"), names(S@graphs)))]]

# get distance matrix
D = CID.GetDistMat(edges)

# entropy-based unclassified labels labels
entropy = CID.entropy(ac = P$L2, distM = D)
}
}
