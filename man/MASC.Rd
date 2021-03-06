% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MASC.R
\name{MASC}
\alias{MASC}
\title{Mixed effect modeling}
\usage{
MASC(
  dataset,
  cluster,
  contrast,
  random_effects = NULL,
  fixed_effects = NULL,
  verbose = FALSE
)
}
\arguments{
\item{dataset}{data frame of covariate, cell type, clustering or disease information}

\item{cluster}{celltypes returned by Signac or cluster identities}

\item{contrast}{Typically disease}

\item{random_effects}{User specified random effect variables in dataset}

\item{fixed_effects}{User specific fixed effects in dataset}

\item{verbose}{If TRUE, algorithm reports outputs}
}
\value{
mixed effect model results
}
\description{
\code{\link{MASC}} was imported from \url{https://github.com/immunogenomics/masc}.
Performs mixed-effect modeling.
}
\examples{
\dontrun{
# Load metadata
file.dir = "https://kleintools.hms.harvard.edu/tools/client_datasets/"
file = "AMP_Phase1_SLE_Apr2019/FullDataset_v1/categorical_coloring_data.json"
download.file(paste0(file.dir, file, "?raw=true"), destfile = "categorical_coloring_data.json")
d = rjson::fromJSON(file='categorical_coloring_data.json')
d = data.frame(sapply(d, function(x) x$label_list))

# run MASC
x = json_data$CellStates # optionally use clusters or cell types
Q = MASC(d, cluster = x, contrast = 'Disease', random_effects = c( "Tissue", "Plate", "Sample"))
}
}
