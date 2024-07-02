
#' Input Parameters for Loading MSI Datasets in .imzML Format
#'
#' Specifies input parameters to be read into R by the Cardinal package.
#' The user must provide the path to paired .imzML and .ibd files, if applicable,
#' for subsequent loading as a single list of MSImagingExperiment objects created
#' by Cardinal.
#'
#' @param dp data folder location, relative to R project directory
#' @param pol ionization mode, either 'pos' or 'neg'. The ionization mode must be present in the corresponding file names.
#' @param mz.low lower limit of mass range in units of m/z
#' @param mz.high upper limit of mass range in units of m/z
#' @param res mass resolution, in ppm
#' @return A list containing .imzML input parameters and file names.
#' @examples
#' MSI.InputParams("data/","neg",300,900,5)
#'
#' @export
MSI.InputParams <- function(
    dp,
    pol,
    mz.low,
    mz.high,
    res){

  data.frame(
    # data folder location (relative to R project directory)
    "path.d" = dp,
    # ionization mode (i.e. 'pos' or 'neg')
    ## Must be present in filename and both an .ibd and .imzML file from SCiLS
    ## Export must be located in same folder
    "polarity" = pol,
    # mass range in units of m/z
    "mz.range.low" = mz.low,
    "mz.range.high" = mz.high,
    # resolution (in ppm)
    "resolution" = res
  )

}

#' Load MSI Datasets in .imzML Format
#'
#' Uses input parameters from MSI.InputParams() to read files into R using the Cardinal function readMSIData().
#'
#' @param list.p list containing input parameters and file names created from MSI.InputParams()
#' @return A list of Cardinal MSImagingExperiment objects created for each file included in the specified data folder.
#' @examples
#' MSI.LoadAll(list.p)
#'
#' @export
MSI.LoadAll <- function(list.p){

  # List samples for processing
  list.s <- data.frame(
    "ID" = seq(
      1:length(
        list.files(
          path = list.p$path.d
        )[grepl(
          list.p$polarity,
          list.files(
            path = list.p$path.d
          )
        ) &
          grepl(
            ".imzML",
            list.files(
              path = list.p$path.d
            )
          )
        ])
    ),
    "Sample.imzML" = list.files(
      path = list.p$path.d
    )[grepl(
      list.p$polarity,
      list.files(
        path = list.p$path.d
      )
    ) &
      grepl(
        ".imzML",
        list.files(
          path = list.p$path.d
        )
      )
    ],
    "Sample.ibd" = list.files(
      path = list.p$path.d
    )[grepl(
      list.p$polarity,
      list.files(
        path = list.p$path.d
      )
    ) &
      grepl(
        ".ibd",
        list.files(
          path = list.p$path.d
        )
      )
    ]
  )


  ## Import data
  list.d <- setNames(
    lapply(
      list.s[["ID"]],
      function(x)
      {

        d <- Cardinal::readMSIData(
          # data file name
          paste(
            list.p$path.d,
            list.s[x,"Sample.imzML"],
            sep = ""
          ),
          # load dataset as sparse matrix
          attach.only = T,
          # mass range
          mass.range = c(
            list.p$mz.range.low,
            list.p$mz.range.high
          ),
          # resolution
          resolution = list.p$resolution,
          units = "ppm"
        )

        if(class(d) == "MSProcessedImagingExperiment") {

          Cardinal::centroided(d) <- F

        }

        return(d)

      }
    ),
    c(
      list.s[["Sample.imzML"]]
    )
  )

  list.d <- list.d[lengths(list.d) > 0]

  list.d <- setNames(
    lapply(
      seq(1:length(list.d)),
      function(x) {
        Cardinal::run(list.d[[x]]) <- list.s[x,"ID"]

        return(list.d[[x]])
      }
    ),
    c(
      list.s[["Sample.imzML"]]
    )
  )

  return(
    list(
      "Sample.Info" = list.s,
      "Data.files" = list.d
      )
    )

  }




