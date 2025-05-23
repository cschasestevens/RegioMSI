#' Normalization Formatting
#'
#' Annotates MSImagingExperiment objects by matching m/z values
#' in a provided reference list. Defaults to matching m/z values
#' within a 10 mDa threshold.
#'
#' @param ldp A processed MSImagingExperiment.
#' @param la A data frame containing reference-matched MSI peaks.
#' @param md A data frame containing treatment and group information.
#' @param ftype Data to use for normalization (either "annotated" or "all").
#' "Annotated" filters the MSI data based on features that are matched
#' to the reference m/z list, whereas "all" retains all features in the
#' dataset.
#' @return Formatted pixel, feature, and image data for downstream analysis.
#' @examples
#'
#' # d_norm <- msi_norm_form(
#' #   ldp = d2a[["Data.filtered"]],
#' #   la = d2a[["Annotated"]],
#' #   md = md1,
#' #   ftype = "annotated"
#' # )
#'
#' @export
msi_norm_form <- function(ldp, la, md, ftype) {
  d <- ldp
  d_anno <- la
  d_group <- md
  ## annotation list
  if(ftype == "annotated") { # nolint
    d_anno <- d_anno[!duplicated(d_anno[["mz.timsTOF"]]), ]
    d_anno[["mz.join"]] <- as.character(
      round(
        d_anno[["mz.timsTOF"]],
        digits = 4
      )
    )
    ## feature metadata
    md_feat <- dplyr::left_join(
      data.frame(
        "mz.join" = as.character(
          round(
            Cardinal::mz(d),
            digits = 4
          )
        ),
        "num.feat" = Cardinal::features(d)
      ),
      d_anno,
      by = "mz.join"
    )
  }
  if(ftype == "all") { # nolint
    md_feat <- d_anno
    md_feat[["Name"]] <- paste(
      rep("X", nrow(md_feat)),
      md_feat[["Feature"]],
      round(md_feat[["mz"]], digits = 4),
      sep = "."
    )
  }
  d_group[["ID"]] <- factor(
    as.character(d_group[["ID"]]),
    levels = c(seq.int(1, nrow(d_group), 1))
  )
  ## sample metadata
  md_samp <- data.frame(
    "ID" = Cardinal::run(d),
    "pixel" = Cardinal::pixels(d),
    "X" = d@elementMetadata[["x"]],
    "Y" = d@elementMetadata[["y"]]
  )
  md_samp <- dplyr::full_join(
    data.frame(
      "ID" = Cardinal::run(d),
      "pixel" = Cardinal::pixels(d),
      "X" = d@elementMetadata[["x"]],
      "Y" = d@elementMetadata[["y"]]
    ),
    d_group,
    by = "ID"
  )
  md_samp <- md_samp[!is.na(md_samp[["pixel"]]), ]

  ## bring spectra into memory
  d_mat <- t(Cardinal::as.matrix(Cardinal::spectra(d, "intensity")))

  return(
    list(
      "Data" = d_mat,
      "Feature" = md_feat,
      "Pixel" = md_samp
    )
  )
}
