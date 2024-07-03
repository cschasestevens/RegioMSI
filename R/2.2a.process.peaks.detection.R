#' Combined Peak Detection for MSI Datasets
#'
#' Formats data and returns processing parameters for peak detection based on all data files included in MSI.LoadAll() using Cardinal.
#'
#' @param msi.list A list of Cardinal MSImagingExperiment objects created from MSI.LoadAll().
#' @return A list containing a formatted 'MSImagingExperiment' object and processing parameters for Cardinal.
#' @examples
#' MSI.detectpeaks(list.d[["Data.files"]])
#'
#' @export
MSI.detectpeaks <- function(
    msi.list){
  list.d <- msi.list

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

  return(
    list(
      "data" = d,
      "chunks" = n.chunk,
      "min.freq" = freq
      )
  )

  }




