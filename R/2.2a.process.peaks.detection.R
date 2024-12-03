#' Combined Peak Processing for MSI Datasets
#'
#' Performs peak picking and alignment of MSI data files.
#'
#' @param ld A list of Cardinal MSImagingExperiment objects.
#' @param lp A list containing input parameters and file names.
#' @param parl Should processing be run in parallel?
#' Note: only supported by WSL2 or Linux;
#' defaults to sequential processing on Windows OS.
#' @param core_perc Percentage of cores to be used during processing.
#' @return A list containing a data frame of processed peaks and QC metrics.
#' @examples
#' # d2 <- msi_peak_proc(
#' #   ld = d[["Data.files"]],
#' #   lp = lp1,
#' #   parl = TRUE,
#' #   core_perc = 0.2
#' # )
#'
#' @export
msi_peak_proc <- function(
  ld,
  lp,
  parl,
  core_perc
) {
  l1 <- ld
  l2 <- lp
  # Set parallelization parameters
  if(Sys.info()[["sysname"]] != "Windows" && # nolint
      parl == TRUE
  ) {
    Cardinal::setCardinalBPPARAM(
      BiocParallel::MulticoreParam(
        workers = floor(
          parallel::detectCores() * core_perc
        ),
        progressbar = TRUE
      )
    )
  }
  if(Sys.info()[["sysname"]] != "Windows" && # nolint
      parl == FALSE
  ) {
    Cardinal::setCardinalBPPARAM(
      BiocParallel::SerialParam()
    )
  }
  if(Sys.info()[["sysname"]] == "Windows") { # nolint
    Cardinal::setCardinalBPPARAM(
      BiocParallel::SerialParam()
    )
    print("Detected OS is Windows; Defaulting to serial processing...")
  }

  ### Combine MSI data and perform peak processing
  msi_d <- lapply(
    seq.int(1, length(l1), 1),
    function(x) l1[[x]]
  )
  d <- do.call(BiocGenerics::cbind, msi_d)

  if(unlist(packageVersion("Cardinal"))[2] < 6) { # nolint
    print(
      "An unsupported version of Cardinal was detected;
      Please update to the latest version..."
    )
  }

  if(unlist(packageVersion("Cardinal"))[2] >= 6) { # nolint
    d2 <- Cardinal::process(
      Cardinal::peakPick(
        d,
        method = "filter",
        SNR = 3
      )
    )
    d3 <- Cardinal::peakAlign(
      d2,
      tolerance = l2[["resolution"]],
      units = "ppm"
    )
  }
  return(
    list(
      "Processed.Peaks" = d3,
      "num.features" = nrow(Cardinal::spectra(d3))
    )
  )
}
