#' Combined Peak Detection for MSI Datasets
#'
#' Performs peak detection based on all data files included in MSI.LoadAll() using Cardinal.
#'
#' @param msi.list A list of Cardinal MSImagingExperiment objects created from MSI.LoadAll().
#' @param p.list A list containing input parameters and file names created from MSI.InputParams()
#' @param parl Should processing be run in parallel? Note: only supported by WSL2 or Linux; will default to sequential processing if toggled 'TRUE' on Windows OS.
#' @param core.perc Percentage of cores to be used during processing.
#' @return A list containing a data frame of processed peaks and processing details.
#' @examples
#' MSI.detectpeaks(list.d[["Data.files"]],list.p,TRUE,0.75)
#'
#' @export
MSI.detectpeaks <- function(
    msi.list,
    p.list,
    parl,
    core.perc){
  list.d <- msi.list
  list.p <- p.list

  # Parallelization parameters
  if(Sys.info()[["sysname"]] != "Windows" &
     parl == TRUE) {
    Cardinal::setCardinalBPPARAM(
      BiocParallel::MulticoreParam(
        workers = floor(
          parallel::detectCores()*core.perc
          ),
        progressbar = T
      )
    )
  }
  if(Sys.info()[["sysname"]] == "Windows" &
     parl == TRUE) {
    Cardinal::setCardinalBPPARAM(
      BiocParallel::SerialParam()
    )
    print("Detected OS is Windows; Defaulting to serial processing...")
  }
  if(Sys.info()[["sysname"]] != "Windows" &
     parl == FALSE) {
    Cardinal::setCardinalBPPARAM(
      BiocParallel::SerialParam()
      )
  }
  if(Sys.info()[["sysname"]] == "Windows" &
     parl == FALSE) {
    Cardinal::setCardinalBPPARAM(
      BiocParallel::SerialParam()
    )
  }

  # Processing parameters

  ifelse(
    sum(
      unlist(
        lapply(
          list.d,
          function(x) nrow(
            Cardinal::coord(x)
          )
        )
      )
    )/
      100 <
      2000,
    n.chunk <- 100,
    ifelse(
      sum(
        unlist(
          lapply(
            list.d,
            function(x)
              nrow(
                Cardinal::coord(x)
              )
          )
        )
      )/
        1000 <
        2000,
      n.chunk <- 1000,
      n.chunk <- 10000
    )
  )

if(unlist(packageVersion("Cardinal"))[2] < 6) {
    Cardinal::setCardinalNumBlocks(
      n = n.chunk
    )
  }

if(unlist(packageVersion("Cardinal"))[2] >= 6) {
  Cardinal::setCardinalNChunks(
    n = n.chunk
  )
}

  ### Align spectra (for every pixel), filter missing or low intensity peaks
  ### (<0.5% detected across all pixels; roughly equivalent to a peak that is detected in 10% of pixels in an individual image for
  ### a dataset that includes 12 separate images)
  ### Generates list of aligned m/z peaks detected for all samples
  ### Use this output as a reference peak list for binning profile spectra to a single centroided spectra suitable for comparing across samples
  ### Then normalize the binned spectra by sparse LOESS (sLOESS) correct for differences in intensity across sample run

  ### Combine MSI data and perform peak detection
  n.freq <- (((sum(
    unlist(
      lapply(
        list.d,
        function(x) nrow(
          Cardinal::coord(x)
        )
      )
    )
  )/
    length(
      list.d
    )))*0.1)/(
      sum(
        unlist(
          lapply(
            list.d,
            function(x) nrow(
              Cardinal::coord(x)
            )
          )
        )
      )
    )

  msi.d <- lapply(
    seq.int(
      1,
      length(list.d),
      1
      ),
    function(x)
      list.d[[x]]
    )

  d <- do.call(
    BiocGenerics::cbind,
    msi.d
    )

  if(unlist(packageVersion("Cardinal"))[1] < 4) {
    d.peaks <- Cardinal::process(
      Cardinal::peakFilter(
        Cardinal::process(
          Cardinal::peakAlign(
            d,
            tolerance = list.p[["resolution"]],
            units = "ppm"
          )
        ),
        freq.min = n.freq
      )
    )
  }

  if(unlist(packageVersion("Cardinal"))[1] > 4) {

    d.peaks <- Cardinal::peakAlign(
      d,
      tolerance = list.p[["resolution"]],
      units = "ppm"
    )

    d.peaks <- Cardinal::subsetFeatures(
      d.peaks,
      freq > n.freq
    )
  }

  return(
    list(
      "Processed.Peaks" = d.peaks,
      "num.features" = length(d.peaks)
      ))

  }




