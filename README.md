# RegioMSI v1.1 (20241210)

Processing and Analysis of Mass Spectrometry Imaging Datasets

## Description

Provides functions for region-specific statistical analysis of biological replicates of mass spectrometry imaging datasets. Introduces improved image normalization methods to reduce intra- and intersample technical artifacts, and offers additional spatial segmentation strategies for unsupervised regional classification.    RegioMSI integrates the strengths of the Cardinal v3 package with a comprehensive framework for parallel data processing, normalization, and analysis of mass spectrometry imaging datasets containing multiple biological replicate samples. The package introduces the sparse LOESS normalization method for mass spectrometry imaging data, which addresses intrasample signal intensity drift and batch effects resulting from data acquisition. Spatial segmentation is performed by adapting the KNN algorithm used by the Seurat package to mass spectrometry imaging data.    Methods for statistical analysis are provided for comparing identified image regions both within individual samples and across treatment groups. Enrichment analysis based on compound class similarity is also supported. Tools are included for plotting ion images and statistical results, enabling seamless, detailed visualization of either single imaging runs or entire mass spectrometry imaging experiments.    RegioMSI is compatible with Windows, Linux, or WSL2. However, analyses conducted in Windows default to sequential processing.

* Important: Processing methods from Seurat (>= v5.0) and Cardinal (>= v3.6) are now the default. Upgrade these packages for compatibility. Additionally, updated methods run in Linux or WSL2 may require installation of additional dependencies within a conda environment.

## Getting Started

### Dependencies
* Windows 10-11, WSL Ubuntu 22.04 or higher, Linux Ubuntu 22.04 or higher, or macOS 12.7.1 or higher
* R version 4.4.1 or higher (https://cran.r-project.org/)
* (Optional) RStudio version 2023.06.2 or higher (https://posit.co/download/rstudio-desktop/)
* R-packages (downloaded from CRAN unless otherwise specified):
    * Suggests: 
        * knitr,
        * rmarkdown,
        * reticulate
    * Imports: 
        * ggplot2

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
    * Update here
* 1.0
    * Initial Release

## License

This project is licensed under the GNU General Public License Version 3 - see the LICENSE.md file for details

## Acknowledgments

* Seurat package: Hao, Y., Stuart, T., Kowalski, M.H. et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nat Biotechnol 42, 293–304 (2024). https://doi.org/10.1038/s41587-023-01767-y
* ComplexHeatmap package: Gu Z, Eils R, Schlesner M (2016). “Complex heatmaps reveal patterns and correlations in multidimensional genomic data.” Bioinformatics. <doi:10.1093/bioinformatics/btw313>
* Circlize package: Gu, Z. circlize implements and enhances circular visualization in R. Bioinformatics 2014.
* Cardinal: <doi:10.1093/bioinformatics/btv146>
