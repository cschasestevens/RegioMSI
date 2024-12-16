#' Sample Annotation Data
#'
#' A subset of data from Stevens et. al. 2024 containing LC-MS/MS based
#' annotations for validating peak IDs from a mass spectrometry imaging
#' dataset.
#'
#' @format ## `an1`
#' A data frame with 644 rows and 9 columns:
#' \describe{
#'   \item{Criteria}{matching criteria for annotation}
#'   \item{Level}{confidence level}
#'   \item{Alignment.ID}{peak ID from LC-MS/MS}
#'   \item{Average.Mz}{peak m/z}
#'   \item{Metabolite.name}{compound name}
#'   \item{Adduct.type}{adduct}
#'   \item{Reference.mz}{reference library m/z}
#'   \item{S_N.average}{signal/noise ratio}
#'   \item{Source}{tissue source}
#'   ...
#' }
#' @source Stevens et. al. 2024
"an1"

#' Sample Metadata Data
#'
#' A subset of metadata from Stevens et. al. 2024 for two
#' samples acquired by MALDI-TOF MS.
#'
#' @format ## `md1`
#' A data frame with 2 rows and 4 columns:
#' \describe{
#'   \item{ID}{sample ID}
#'   \item{file.name}{file name and extension}
#'   \item{sample.name}{alternate sample name}
#'   \item{Group}{treatment group}
#'   ...
#' }
#' @source Stevens et. al. 2024
"md1"

#' Demo Data
#'
#' A subset of data from Stevens et. al. 2024 for two
#' samples acquired by MALDI-TOF MS processed by msi_peak_proc().
#'
#' @format ## `d_norm`
#' A list of length 2:
#' \describe{
#'   \item{Sample.Info}{
#' A data frame containing input params. from msi_input_param()}
#'   \item{Data.files}{Two separate MSImagingExperiment objects
#' containing test data}
#'   ...
#' }
#' @source Stevens et. al. 2024
"d_norm"