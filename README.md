# RegioMSI v1.1 (20241212)

Processing and Analysis of Mass Spectrometry Imaging Datasets

## Description

Provides functions for region-specific statistical analysis of biological replicates of mass spectrometry imaging datasets. Introduces improved image normalization methods to reduce intra- and intersample technical artifacts, and offers additional spatial segmentation strategies for unsupervised regional classification.    RegioMSI integrates the strengths of the Cardinal v3 package with a comprehensive framework for parallel data processing, normalization, and analysis of mass spectrometry imaging datasets containing multiple biological replicate samples. The package introduces the sparse LOESS normalization method for mass spectrometry imaging data, which addresses intrasample signal intensity drift and batch effects resulting from data acquisition. Spatial segmentation is performed by adapting the KNN algorithm used by the Seurat package to mass spectrometry imaging data.    Methods for statistical analysis are provided for comparing identified image regions both within individual samples and across treatment groups. Enrichment analysis based on compound class similarity is also supported. Tools are included for plotting ion images and statistical results, enabling seamless, detailed visualization of either single imaging runs or entire mass spectrometry imaging experiments.    RegioMSI is compatible with Windows, Linux, or WSL2. However, analyses conducted in Windows default to sequential processing.

* Important: Processing methods from Seurat (>= 5.1.0) and Cardinal (>= 3.6.2) are now the default. Upgrade these packages for compatibility. Additionally, updated methods run in Linux or WSL2 may require installation of additional dependencies within a conda environment.

## Getting Started

### Dependencies
* Windows 10-11, WSL Ubuntu 22.04 or higher, Linux Ubuntu 22.04 or higher, or macOS 12.7.1 or higher
* R version 4.4.1 or higher (https://cran.r-project.org/)
* (Optional) RStudio version 2023.06.2 or higher (https://posit.co/download/rstudio-desktop/)
* R-packages (downloaded from CRAN unless otherwise specified):
    * Suggests: 
        * knitr
        * rmarkdown
        * BiocManager
    * Imports: 
        * ggsci
        * viridis
        * ggplot2
        * Cardinal (>= 3.6.2)
        * BiocParallel
        * parallel
        * BiocGenerics
        * dplyr
        * fuzzyjoin
        * magrittr
        * ggnewscale
        * Seurat (>= 5.1.0)
        * harmony
        * ggpubr
        * gtools
        * SeuratObject (>= 5.0.2)
        * circlize
        * ComplexHeatmap
        * grid
        * reshape2
        * EnhancedVolcano

### Installation
* Run the following in a new R session on the command line or within R-Studio:

```
devtools::install_github(
  "cschasestevens/RegioMSI", 
  ref = "master", 
  build_vignettes = TRUE
)
```

## Help
* Browse vignettes by running the following:

```
browseVignettes("RegioMSI")
```

* Access function documentation by running the following:

```
# Type function name after package name
?RegioMSI::msi_stat_anova
```

## Authors

* Nathanial Chase Stevens, PhD, University of North Carolina at Chapel Hill
* Email: Nathanial_Stevens@med.unc.edu
* Alternate email: cschasestevens@gmail.com
* LinkedIn: https://www.linkedin.com/in/nathanial-chase-stevens-phd-08775180/

## Version History
* 1.1
    * Peak processing, including detection, alignment, and binning is now performed using a single function call
    * Added filters for automatic removal of artifacts from a processed peak list
    * Added root-mean-square normalization for benchmarking sparse LOESS normalization performance
    * Feature normalization can now be performed on both annotated and unknown features from a processed peak list
    * All ion image plotting for visualizing total pixel ion counts, feature intensity, colocalization, and segmented clusters is now supported by a single function
    * Image segmentation now supports options for either a single sample or entire dataset
    * Chemical similarity enrichment analysis now supports inputs directly generated from msi_stat_anova()
    * Simplified syntax for all package functions
    * Added default function parameters for greater cross-platform compatibility between Windows and Linux
* 1.0
    * Initial Release

## License

This project is licensed under the MIT License - see the LICENSE.md file for details.

## Acknowledgments

* Seurat package: Hao, Y., Stuart, T., Kowalski, M.H. et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nat Biotechnol 42, 293–304 (2024). https://doi.org/10.1038/s41587-023-01767-y
* ComplexHeatmap package: Gu Z, Eils R, Schlesner M (2016). “Complex heatmaps reveal patterns and correlations in multidimensional genomic data.” Bioinformatics. <doi:10.1093/bioinformatics/btw313>
* Circlize package: Gu, Z. circlize implements and enhances circular visualization in R. Bioinformatics 2014.
* Cardinal package: <doi:10.1093/bioinformatics/btv146>
