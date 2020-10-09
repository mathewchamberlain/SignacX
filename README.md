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
    <img src="images/sanofi_logo.png" alt="Logo" width="600" height="300">
  </a>

  <h3 align="center">Signac</h3>

  <p align="center">
    Get the most out of your single cell data.
    <br />
    <a href="https://htmlpreview.github.io/?https://github.com/mathewchamberlain/Signac/master/vignettes/signac-Seurat_CITEseq.html"><strong>Explore the docs »</strong></a>
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
  * [Built With](#built-with)
* [Getting Started](#getting-started)
  * [Installation](#installation)
* [Usage](#usage)
* [Roadmap](#roadmap)
* [Contributing](#contributing)
* [License](#license)
* [Contact](#contact)
* [Acknowledgements](#acknowledgements)



<!-- ABOUT THE PROJECT -->
## What is Signac?

Although Seurat is useful for performing basic visualization and quantification of single cell RNA-seq data, it does not perform statistical analyses that are crucial for clinical applications. We wrote Signac to help solve that problem. Here are the benefits of Signac:

* Signac implements regression-based ambient RNA correction, which corrects for a known confounding factor in droplet-based technologies.
* Often, a single sample in a dataset contributes >10-fold more cells than any other sample in the dataset. This one sample can dominate standard downstream analysis steps such as clustering, differential expression, cell type classification, differential abundance and gene co-expression network analysis. Signac implements jackknifing to solve this problem, which removes single-sample artifacts from clustering, differential gene expression analysis and gene co-expression analysis.
* Unlike Seurat, Signac is memory optimized to handle >100 datasets in R simultaneously. Signac utilizes sparse hdf5 files to load individual genes and cells from datasets rather than the entire expression matrix.

This package currently has publicly available only some of these features. I'll be adding more in the near future. You may also suggest changes by forking this repo and creating a pull request or opening an issue.

### Built With

* [Seurat](https://satijalab.org/seurat/)
* [rhdf5](https://www.bioconductor.org/packages/release/bioc/html/rhdf5.html)
* [parallel](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf)

<!-- GETTING STARTED -->
## Getting Started

To install Signac in R, simply do:

### Installation

```r
devtools::install_github("mathewchamberlain/Signac")
```

<!-- USAGE EXAMPLES -->
## Usage

Running Signac is simple. [Here is an example vignette](https://htmlpreview.github.io/?https://github.com/mathewchamberlain/Signac/blob/master/vignettes/Signac-Seurat_CITEseq.html) for processing a 1:1 mixture of human and mouse cells hosted on 10X. More to come.

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
