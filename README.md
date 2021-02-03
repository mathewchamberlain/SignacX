[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![GPL License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]

<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://www.sanofi.com/">
    <img src="images/GitHubFigure.svg" alt="Logo" width="600" height="300">
  </a>

  <h3 align="center">Signac</h3>

  <p align="center">
    Get the most out of your single cell data.
    <br />
    <a href="https://github.com/mathewchamberlain/Signac"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://htmlpreview.github.io/?https://github.com/mathewchamberlain/Signac/master/vignettes/signac-Seurat_CITEseq.html">View Demo</a>
    ·
    <a href="https://github.com/mathewchamberlain/Signac/issues">Report Bug</a>
    ·
    <a href="https://github.com/mathewchamberlain/Signac/issues">Request Feature</a>
  </p>
</p>



<!-- TABLE OF CONTENTS -->
## Table of Contents

* [What is Signac?](#about-the-project)
* [Getting Started](#getting-started)
  * [Installation](#installation)
  * [Quick start](#quickstart)
* [Usage](#usage)
* [Roadmap](#roadmap)
* [Contributing](#contributing)
* [License](#license)
* [Contact](#contact)
* [Acknowledgements](#acknowledgements)



<!-- ABOUT THE PROJECT -->
## What is Signac?

Signac helps solve the cell type classification problem in single cell RNA sequencing: We sequenced the RNA for each individual cell, but we do not know the identity of each cellular phenotype. Signac classifies each cell in scRNA-seq data using neural networks trained with sorted bulk gene expression data from the [Human Primary Cell Atlas](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-632). Check out our pre-print here: LINK.

To make life easier, Signac is integrated with [Seurat](https://satijalab.org/seurat/) (versions 3 and 4), [MASC](https://pubmed.ncbi.nlm.nih.gov/30333237/) and [SPRING](https://pubmed.ncbi.nlm.nih.gov/29228172/), allowing for easy downstream analysis of single cell data.

<!-- GETTING STARTED -->
## Getting Started

To install Signac in R, simply do:

### Installation

```r
devtools::install_github("mathewchamberlain/Signac")
```

### Quick start

The main functions in Signac are:

```r
# load the library
library(Signac)

# Load the reference training data set
data("training_HPCA")

# Generate initial labels
labels = Signac(your_data_here, R = training_HPCA)

# Get cell type labels
celltypes = Generate_lbls(labels, your_data_here)
```

<!-- USAGE EXAMPLES -->
## Usage

Running Signac is simple. We provide a few vignettes:

* The easiest way to use Signac is with Seurat; it is convenient to use Seurat functions (like differential expression, clustering and visualization) with cellular phenotypes identified with Signac. [Here, we perform multi-modal analysis of CITE-seq PBMCs from 10X Genomics in R using Seurat and Signac.](https://htmlpreview.github.io/?https://github.com/mathewchamberlain/Signac/master/vignettes/signac-Seurat_CITEseq.html)
* A different way to use Signac is to inetgrate it with [SPRING](https://pubmed.ncbi.nlm.nih.gov/29228172/). Using SPRING with Signac, we can easily explore cellular phenotypes interactively. Here, we provide a Jupyter notebook for [downloading scRNA-seq data from 10X Genomics and processing it with SPRING in Jupyter](https://github.com/mathewchamberlain/Signac/blob/master/vignettes/spring_notebook_10X.ipynb). Once this notebook has been run, Signac is integrated seamlessly with the output files. Like this:

```r
# load the Signac library
library(Signac)

# dir points to the "FullDataset_v1" directory generated by the Jupyter notebook
dir = "Your_Path_to_FullDataset_v1" 

# load the expression data
E = CID.LoadData(dir)

# load the reference data
data(training_HPCA)

# generate cellular phenotype labels
labels = Signac(E, R = training_HPCA, spring.dir = dir)
celltypes = Generate_lbls(labels, E = E)

# write cell types and Louvain clusters to SPRING
dat <- CID.writeJSON(celltypes, data.dir = dir)
```

After running these functions, cellular phenotypes and Louvain clusters are ready to be visualized with the SPRING interface, which can be setup locally as described [here](https://github.com/AllonKleinLab/SPRING_dev). 

<!-- ROADMAP -->
## Roadmap

See the [open issues](https://github.com/mathewchamberlain/Signac/issues) for a list of proposed features (and known issues).

<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to be learn, inspire, and create. Any contributions you make are **greatly appreciated**.

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<!-- LICENSE -->
## License

Distributed under the GPL v3.0 License. See `LICENSE` for more information.

<!-- CONTACT -->
## Contact

Mathew Chamberlain - [linkedin](https://linkedin.com/in/chamberlainmathew) - mathew.chamberlain@sanofi.com

Project Link: [https://github.com/mathewchamberlain/Signac](https://github.com/mathewchamberlain/Signac)

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/mathewchamberlain/Signac.svg?style=flat-square
[contributors-url]: https://github.com/mathewchamberlain/Signac/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/mathewchamberlain/Signac.svg?style=flat-square
[forks-url]: https://github.com/mathewchamberlain/Signac/network/members
[stars-shield]: https://img.shields.io/github/stars/mathewchamberlain/Signac.svg?style=flat-square
[stars-url]: https://github.com/mathewchamberlain/Signac/stargazers
[issues-shield]: https://img.shields.io/github/issues/mathewchamberlain/Signac.svg?style=flat-square
[issues-url]: https://github.com/mathewchamberlain/Signac/issues
[license-shield]: https://img.shields.io/github/license/mathewchamberlain/Signac.svg?style=flat-square
[license-url]: https://choosealicense.com/licenses/gpl-3.0/
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=flat-square&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/chamberlainmathew
[product-screenshot]: images/screenshot.png
