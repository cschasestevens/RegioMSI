#' Input Parameters for Loading MSI Datasets in .imzML Format
#'
#' Specifies input parameters to be read into R by the Cardinal package.
#' The user must provide the path to a folder containing .imzml files,
#' for creating a single list of MSImagingExperiment objects using
#' Cardinal. See doi:10.1007/978-1-60761-987-1_12 for more information
#' about .imzml format and doi:10.1093/bioinformatics/btv146 for more details
#' regarding the Cardinal R package. Note that for .imzml files exported from
#' commercial software, an additional .ibd file may accompany each .imzml,
#' which should both be located in the same folder.
#'
#' @param dp Data path relative to the "data/" folder. Data files
#' should be separated by subfolders according to ionization mode.
#' @param pol Ionization mode (either "pos" or "neg").
#' @param mz_range A vector specifying the limits of the m/z range.
#' @param mz_res Mass resolution, in ppm.
#' @return A data frame containing .imzML input parameters and file names.
#' @examples
#'
#' # lp1 <- msi_input_param(
#' #   dp = "pos/",
#' #   pol = "pos",
#' #   mz_range = c(300, 1300),
#' #   mz_res = 5
#' # )
#'
#' @export
msi_input_param <- function(
  dp,
  pol,
  mz_range,
  mz_res
) {
  d <- data.frame(
    "path.d" = paste("data/", dp, sep = ""),
    "polarity" = pol,
    "mz.range.low" = mz_range[[1]],
    "mz.range.high" = mz_range[[2]],
    "resolution" = mz_res
  )
  return(d)
}

#' Load .imzml Datasets
#'
#' Uses specified input parameters to read files into R using Cardinal.
#'
#' @param lp List containing input parameters from msi_input_param().
#' @return A list of Cardinal MSImagingExperiment objects.
#' @examples
#'
#' # d <- msi_load_data(lp1)
#'
#' @export
msi_load_data <- function(lp) {
  # List samples for processing
  list_s <- data.frame(
    "ID" = seq.int(
      1,
      length(
        list.files(
          path = lp$path.d
        )[grepl(
          ".imzML",
          list.files(
            path = lp$path.d
          )
        )]
      ),
      1
    ),
    "Sample.imzML" = list.files(
      path = lp$path.d
    )[grepl(
      ".imzML",
      list.files(
        path = lp$path.d
      )
    )],
    "Sample.ibd" = list.files(
      path = lp$path.d
    )[grepl(
      ".ibd",
      list.files(
        path = lp$path.d
      )
    )]
  )

  ## Import data
  list_d <- setNames(
    lapply(
      list_s[["ID"]],
      function(x) {
        d <- Cardinal::readMSIData(
          # data file name
          paste(
            lp$path.d,
            list_s[x, "Sample.imzML"],
            sep = ""
          ),
          # load dataset as sparse matrix
          memory = FALSE,
          # mass range
          mass.range = c(
            lp$mz.range.low,
            lp$mz.range.high
          ),
          # resolution
          resolution = lp$resolution,
          units = "ppm"
        )
        if(class(d) == "MSImagingExperiment") { # nolint
          Cardinal::centroided(d) <- FALSE
        }
        if(class(d) == "MSProcessedImagingExperiment") { # nolint
          Cardinal::centroided(d) <- FALSE
        }
        return(d)
      }
    ),
    c(
      list_s[["Sample.imzML"]]
    )
  )
  list_d <- list_d[lengths(list_d) > 0]
  list_d <- setNames(
    lapply(
      seq.int(1, length(list_d), 1),
      function(x) {
        Cardinal::run(list_d[[x]]) <- list_s[x, "ID"]
        return(list_d[[x]])
      }
    ),
    c(
      list_s[["Sample.imzML"]]
    )
  )
  return(
    list(
      "Sample.Info" = list_s,
      "Data.files" = list_d
    )
  )
}
