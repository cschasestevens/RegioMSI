#' Peak Annotation of MSI Data
#'
#' Annotates processed peaks by matching m/z values
#' from a reference annotation list, such as a complimentary
#' LC-MS/MS dataset acquired independently from an MSI experiment.
#' Defaults to matching m/z values between MSI and reference data within 10 mDa.
#'
#' @param ldp A data frame of processed peaks.
#' @param la A reference m/z list, provided as a data frame.
#' @return A reference-matched data frame containing annotated MSI peaks.
#' @examples
#'
#' # d2a <- msi_peak_anno(
#' #   ldp = d2,
#' #   la = an1
#' # )
#'
#' @export
msi_peak_anno <- function(
  ldp,
  la
) {
  # Load peak and LC-MS/MS annotation list
  d_peaks <- ldp[["Processed.Peaks"]]
  list_anno <- la

  # Extract detected peaks and peak intensities
  list_f <- data.frame(
    "Feature" = Cardinal::features(d_peaks),
    "mz" = Cardinal::mz(d_peaks),
    "Count" = d_peaks@featureData[["count"]],
    "Frequency" = d_peaks@featureData[["freq"]]
  )

  list_anno_filt <- setNames(
    list_anno[!grepl(
      "iSTD",
      list_anno[["Metabolite.name"]]
    ),
    c("Metabolite.name",
      "Adduct.type",
      "Average.Mz",
      "Source"
    )],
    c("Name",
      "Adduct",
      "mz",
      "Source"
    )
  )

  ## Assign IDs based on LC-MS/MS data
  list_f_anno <- setNames(
    dplyr::select(
      fuzzyjoin::fuzzy_join(
        list_f,
        list_anno_filt,
        by = "mz",
        match_fun = ~ abs(.x - .y) < 0.01
      ),
      c("Feature",
        "Source", "Name",
        "Adduct",
        "mz.x", "mz.y",
        "Count", "Frequency"
      )
    ),
    c(
      "Feature",
      "Source", "Name", "Adduct",
      "mz.timsTOF", "mz.LCMS",
      "Count", "Frequency"
    )
  )
  list_f_anno <- list_f_anno[!duplicated(list_f_anno), ]
  ## Return unknowns
  list_f_unk <- list_f[list_f[-c(list_f_anno$mz.timsTOF), "mz"], ]
  ## Return filtered data
  list_sub <- unique(list_f_anno[["Feature"]])
  list_d_anno <- d_peaks[c(list_sub), ]

  return(
    list(
      "Data.filtered" = list_d_anno,
      "Data.original" = d_peaks,
      "Annotated" = list_f_anno,
      "Unassigned" = list_f_unk
    )
  )
}
