#' MSI Processing QC
#'
#' Calculates descriptive statistics for each feature in a formatted
#' MSI dataset. Can be used for subsequent artifact detection and
#' removal.
#'
#' @param dm A MSI data matrix.
#' @param feat A data frame containing metadata for each annotated compound.
#' @param pix A data frame containing sample metadata for each pixel
#' @param mtd Normalization method (either "none", "tic", "rms", or "sLOESS").
#' @param parl (optional) Logical indicating if normalization
#' should be performed in parallel.
#' @param core_perc (optional) Percentage of cores to use for parallel
#' normalization.
#' @return A list containing normalized pixel intensities and metadata.
#' @examples
#'
#' # data_qc <- msi_data_qc(
#' #   dm = d_norm[["Data"]],
#' #   feat = d_norm[["Feature"]],
#' #   pix = d_norm[["Pixel"]],
#' #   parl = FALSE,
#' #   core_perc = 0.25
#' # )
#'
#' @export
msi_data_qc <- function(
  dm,
  feat,
  pix,
  parl = FALSE,
  core_perc = NULL
) {
  d <- dm
  mft <- feat
  mpx <- pix
  if(length(mft[["Name"]]) != ncol(d)) { # nolint
    print("Error: Number of features does not match number of data columns!")
  }
  if(length(mft[["Name"]]) == ncol(d)) { # nolint
    if(parl == TRUE) { # nolint
      d_chk <- data.frame(
        setNames(
          as.data.frame(d), mft[["Name"]]
        ),
        "ID" = mpx[["ID"]],
        "X" = mpx[["X"]],
        "Y" = mpx[["Y"]]
      )
      # Flag potential artifacts based on descriptive stats
      d_qc <- magrittr::set_rownames(dplyr::bind_rows(setNames(
        parallel::mclapply(
          mc.cores = ceiling(parallel::detectCores() * core_perc),
          seq.int(1, ncol(d), 1),
          function(x) {
            d <- data.frame(
              "Name" = names(d_chk)[[x]],
              "int.mean" = round(mean(d_chk[[x]]), digits = 2),
              "int.min" = round(min(d_chk[[x]]), digits = 2),
              "int.med" = median(d_chk[[x]]),
              "int.max" = round(max(d_chk[[x]]), digits = 2),
              "sd" = round(sd(d_chk[[x]]), digits = 2),
              "freq" = round(dplyr::count(d_chk, d_chk[[x]] > 0)[2, 2] /
                  length(d_chk[[x]]), digits = 2
              ),
              "q1" = round(quantile(d_chk[[x]], 0.25), digits = 2),
              "q2" = round(quantile(d_chk[[x]], 0.5), digits = 2),
              "q3" = round(quantile(d_chk[[x]], 0.75), digits = 2),
              "q4" = round(quantile(d_chk[[x]], 0.99), digits = 2)
            )
            return(d)
          }
        ),
        mft[["Name"]]
      )), seq.int(1, ncol(d), 1))
      d_qc[["artifact.flag"]] <- ifelse(
        d_qc[["freq"]] == 0 | d_qc[["q4"]] == 0,
        TRUE,
        FALSE
      )
    }
    if(parl == FALSE | Sys.info()[["sysname"]] == "Windows") { # nolint
      d_chk <- data.frame(
        setNames(
          as.data.frame(d), mft[["Name"]]
        ),
        "ID" = mpx[["ID"]],
        "X" = mpx[["X"]],
        "Y" = mpx[["Y"]]
      )
      # Flag potential artifacts based on descriptive stats
      d_qc <- magrittr::set_rownames(dplyr::bind_rows(setNames(
        lapply(
          seq.int(1, ncol(d), 1),
          function(x) {
            d <- data.frame(
              "Name" = names(d_chk)[[x]],
              "int.mean" = round(mean(d_chk[[x]]), digits = 2),
              "int.min" = round(min(d_chk[[x]]), digits = 2),
              "int.med" = median(d_chk[[x]]),
              "int.max" = round(max(d_chk[[x]]), digits = 2),
              "sd" = round(sd(d_chk[[x]]), digits = 2),
              "freq" = round(dplyr::count(d_chk, d_chk[[x]] > 0)[2, 2] /
                  length(d_chk[[x]]), digits = 2
              ),
              "q1" = round(quantile(d_chk[[x]], 0.25), digits = 2),
              "q2" = round(quantile(d_chk[[x]], 0.5), digits = 2),
              "q3" = round(quantile(d_chk[[x]], 0.75), digits = 2),
              "q4" = round(quantile(d_chk[[x]], 0.99), digits = 2)
            )
            return(d)
          }
        ),
        mft[["Name"]]
      )), seq.int(1, ncol(d), 1))
      d_qc[["artifact.flag"]] <- ifelse(
        d_qc[["freq"]] == 0 | d_qc[["q4"]] == 0,
        TRUE,
        FALSE
      )
    }
  }
  return(d_qc)
}

#' MSI Data Normalization
#'
#' Includes multiple methods for MSI data normalization. RegioMSI introduces
#' sparse locally estimated scatterplot smoothing (sLOESS) normalization, which
#' excludes pixels with zero signal intensity to preserve biological differences
#' in compound spatial distribution while removing systematic artifacts
#' incurred during data acquisition.
#'
#' @param dm A MSI data matrix.
#' @param feat A data frame containing metadata for each annotated compound.
#' @param pix A data frame containing sample metadata for each pixel
#' @param mtd Normalization method (either "none", "tic", "rms", or "sLOESS").
#' @param parl (optional) Logical indicating if normalization
#' should be performed in parallel.
#' @param core_perc (optional) Percentage of cores to use for parallel
#' normalization.
#' @return A list containing normalized pixel intensities and metadata.
#' @examples
#'
#' # d_norm <- msi_data_norm(
#' #   dm = d_norm[["Data"]],
#' #   feat = d_norm[["Feature"]],
#' #   pix = d_norm[["Pixel"]],
#' #   mtd = "sLOESS",
#' #   parl = TRUE,
#' #   core_perc = 0.2
#' # )
#'
#' # d_norm <- msi_data_norm(
#' #   dm = d_norm[["Data"]],
#' #   feat = d_norm[["Feature"]],
#' #   pix = d_norm[["Pixel"]],
#' #   mtd = "tic"
#' # )
#'
#' @export
msi_data_norm <- function( # nolint
  dm,
  feat,
  pix,
  mtd,
  parl = FALSE,
  core_perc = NULL
) {
  # Load objects
  d <- dm
  mft <- feat
  mpx <- pix

  # No normalization
  if(mtd == "none") { # nolint
    print("Performing no normalization...")
    if(Sys.info()[["sysname"]] == "Windows" | parl == FALSE) { # nolint
      d <- as.data.frame(d)
      mpx[["TIC.norm.none"]] <- unlist(lapply(
        seq.int(1, ncol(t(as.matrix(d))), 1),
        function(x) sum(d[x, ])
      ))
    }
    if(parl == TRUE) { # nolint
      d <- as.data.frame(d)
      mpx[["TIC.norm.none"]] <- unlist(parallel::mclapply(
        mc.cores = ceiling(parallel::detectCores() * core_perc),
        seq.int(1, ncol(t(as.matrix(d))), 1),
        function(x) sum(d[x, ])
      ))
    }
    d1 <- list(
      "data" = d,
      "features" = mft,
      "pixels" = mpx,
      "norm.method" = "none"
    )
  }
  # TIC normalization
  if(mtd == "tic") { # nolint
    if(Sys.info()[["sysname"]] == "Windows" | parl == FALSE) { # nolint
      print(
        "Detected OS is Windows or parl is FALSE;
        Defaulting to sequential processing..."
      )
      d <- as.data.frame(d)
      dtic <- data.frame(
        "pixel.TIC" = unlist(
          lapply(
            seq.int(1, ncol(t(as.matrix(d))), 1),
            function(x) sum(d[x, ])
          )
        )
      )
      dtic <- dplyr::select(
        dplyr::left_join(
          data.frame("Group" = mpx[["Group"]], dtic),
          setNames(
            aggregate(
              dtic[["pixel.TIC"]],
              list(
                mpx[["Group"]]
              ),
              function(x) sum(x)
            ),
            c("Group", "TIC.sum")
          ),
          by = "Group"
        ),
        c("TIC.sum", dplyr::everything())
      )
      d <- as.data.frame(lapply(
        seq.int(1, ncol(d), 1),
        function(x) {
          (d[, x] / dtic[["TIC.sum"]]) *
            mean(dtic[["TIC.sum"]])
        }
      ))
      mpx[["TIC.norm.tic"]] <- unlist(lapply(
        seq.int(1, ncol(t(as.matrix(d))), 1),
        function(x) sum(d[x, ])
      ))
    }
    if(Sys.info()[["sysname"]] != "Windows" && # nolint
        parl == TRUE
    ) {
      d <- as.data.frame(d)
      dtic <- data.frame(
        "pixel.TIC" = unlist(
          parallel::mclapply(
            mc.cores = parallel::detectCores() * core_perc,
            seq.int(1, ncol(t(as.matrix(d))), 1),
            function(x) sum(d[x, ])
          )
        )
      )
      dtic <- dplyr::select(
        dplyr::left_join(
          data.frame("Group" = mpx[["Group"]], dtic),
          setNames(
            aggregate(
              dtic[["pixel.TIC"]],
              list(
                mpx[["Group"]]
              ),
              function(x) sum(x)
            ),
            c("Group", "TIC.sum")
          ),
          by = "Group"
        ),
        c("TIC.sum", dplyr::everything())
      )
      d <- as.data.frame(parallel::mclapply(
        mc.cores = parallel::detectCores() * core_perc,
        seq.int(1, ncol(d), 1),
        function(x) {
          (d[, x] / dtic[["TIC.sum"]]) *
            mean(dtic[["TIC.sum"]])
        }
      ))
      mpx[["TIC.norm.tic"]] <- unlist(parallel::mclapply(
        mc.cores = parallel::detectCores() * core_perc,
        seq.int(1, ncol(t(as.matrix(d))), 1),
        function(x) sum(d[x, ])
      ))
    }
    d1 <- list(
      "data" = d,
      "features" = mft,
      "pixels" = mpx,
      "norm.method" = "tic"
    )
  }
  # RMS normalization
  if(mtd == "rms") { # nolint
    if(Sys.info()[["sysname"]] == "Windows" | parl == FALSE) { # nolint
      print(
        "Detected OS is Windows or parl is FALSE;
        Defaulting to sequential processing..."
      )
      d <- as.data.frame(d)
      dtic <- data.frame(
        "pixel.RMS" = unlist(
          lapply(
            seq.int(1, ncol(t(d)), 1),
            function(x) sqrt(mean(sum((d[x, ]) ^ 2)))
          )
        )
      )
      d <- as.data.frame(lapply(
        seq.int(1, ncol(d), 1),
        function(x) {
          (d[, x] / dtic[["pixel.RMS"]]) *
            mean(dtic[["pixel.RMS"]])
        }
      ))
      mpx[["TIC.norm.rms"]] <- unlist(lapply(
        seq.int(1, ncol(t(as.matrix(d))), 1),
        function(x) sum(d[x, ])
      ))
      mpx[is.nan(mpx[["TIC.norm.rms"]]), "TIC.norm.rms"] <- 0
    }
    if(Sys.info()[["sysname"]] != "Windows" && # nolint
        parl == TRUE
    ) {
      d <- as.data.frame(d)
      dtic <- data.frame(
        "pixel.RMS" = unlist(
          parallel::mclapply(
            mc.cores = parallel::detectCores() * core_perc,
            seq.int(1, ncol(t(as.matrix(d))), 1),
            function(x) sqrt(mean(sum((d[x, ]) ^ 2)))
          )
        )
      )
      d <- as.data.frame(parallel::mclapply(
        mc.cores = parallel::detectCores() * core_perc,
        seq.int(1, ncol(d), 1),
        function(x) {
          (d[, x] / dtic[["pixel.RMS"]]) *
            mean(dtic[["pixel.RMS"]])
        }
      ))
      mpx[["TIC.norm.rms"]] <- unlist(parallel::mclapply(
        mc.cores = parallel::detectCores() * core_perc,
        seq.int(1, ncol(t(as.matrix(d))), 1),
        function(x) sum(d[x, ])
      ))
      mpx[is.nan(mpx[["TIC.norm.rms"]]), "TIC.norm.rms"] <- 0
    }
    d1 <- list(
      "data" = d,
      "features" = mft,
      "pixels" = mpx,
      "norm.method" = "rms"
    )
  }
  # sLOESS normalization (can divide into batches)
  if(mtd == "sLOESS") { # nolint
    if(Sys.info()[["sysname"]] != "Windows" && # nolint
        parl == TRUE
    ) {
      d <- as.data.frame(parallel::mclapply(
        mc.cores = ceiling(parallel::detectCores() * core_perc),
        seq.int(1, ncol(as.data.frame(d)), 1),
        function(x) {
          d <- data.frame(
            "ID" = mpx[["ID"]],
            "pixel" = mpx[["pixel"]],
            as.data.frame(d)[[x]]
          )
          # Filter 0 values to ignore in LOESS calculation
          dl2 <- d[d[[3]] > 0, ]
          # Run LOESS
          dl2_fit <- as.data.frame(
            lowess(
              x = dl2[["pixel"]],
              y = dl2[[3]],
              # use 10% of the average pixel number per sample for span
              f = ((length(dl2[["pixel"]]) / length(
                unique(dl2[["ID"]])
              )
              ) * 0.1) /
                length(dl2[["pixel"]]),
              iter = 5,
              delta = 0
            )
          )
          # Normalize based on fitted model
          dl3 <- setNames(
            aggregate(
              dl2[[3]],
              list(dl2[["ID"]]),
              function(x) mean(x)
            ),
            c("ID", "pixel.mean")
          )
          dl2[["norm"]] <- (dl2[[3]] / dl2_fit[["y"]]) *
            median(dl3[["pixel.mean"]])
          dl4 <- dplyr::left_join(
            d,
            dl2,
            by = "pixel"
          )[, "norm"]
          dl4[is.na(dl4)] <- 0
          dl4 <- setNames(as.data.frame(dl4), paste("X", x, sep = "."))
          return(dl4)
        }
      ))
      mpx[["TIC.norm.sloess"]] <- unlist(parallel::mclapply(
        mc.cores = parallel::detectCores() * core_perc,
        seq.int(1, ncol(t(as.matrix(d))), 1),
        function(x) sum(d[x, ])
      ))
    }
    if(Sys.info()[["sysname"]] == "Windows" | # nolint
        parl == FALSE
    ) {
      d <- as.data.frame(lapply(
        seq.int(1, ncol(as.data.frame(d)), 1),
        function(x) {
          d <- data.frame(
            "ID" = mpx[["ID"]],
            "pixel" = mpx[["pixel"]],
            as.data.frame(d)[[x]]
          )
          # Filter 0 values to ignore in LOESS calculation
          dl2 <- d[d[[3]] > 0, ]
          # Run LOESS
          dl2_fit <- as.data.frame(
            lowess(
              x = dl2[["pixel"]],
              y = dl2[[3]],
              # use 10% of the average pixel number per sample for span
              f = ((length(dl2[["pixel"]]) / length(
                unique(dl2[["ID"]])
              )
              ) * 0.1) /
                length(dl2[["pixel"]]),
              iter = 5,
              delta = 0
            )
          )
          # Normalize based on fitted model
          dl3 <- setNames(
            aggregate(
              dl2[[3]],
              list(dl2[["ID"]]),
              function(x) mean(x)
            ),
            c("ID", "pixel.mean")
          )
          dl2[["norm"]] <- (dl2[[3]] / dl2_fit[["y"]]) *
            median(dl3[["pixel.mean"]])
          dl4 <- dplyr::left_join(
            d,
            dl2,
            by = "pixel"
          )[, "norm"]
          dl4[is.na(dl4)] <- 0
          dl4 <- setNames(as.data.frame(dl4), paste("X", x, sep = "."))
          return(dl4)
        }
      ))
      mpx[["TIC.norm.sloess"]] <- unlist(lapply(
        seq.int(1, ncol(t(as.matrix(d))), 1),
        function(x) sum(d[x, ])
      ))
    }
    d1 <- list(
      "data" = as.matrix(d),
      "features" = mft,
      "pixels" = mpx,
      "norm.method" = "sloess"
    )
  }
  return(d1)
}

#' Normalization QC
#'
#' Scatter plot visualization to evaluate normalization performance.
#'
#' @param df A data frame containing pixel metadata.
#' @param var_x X-axis variable (usually the acquisition order of each pixel).
#' @param var_y Y-axis (usually the TIC intensity for each pixel).
#' @param var_g Grouping variable (usually the imaging run ID for each sample).
#' @param y_lim Upper limit of the y-axis.
#' @return A scatter plot visualizing pixel intensities over time.
#' @examples
#'
#' # ptic <- msi_plot_tic(
#' #   df = d_norm[["pixels"]],
#' #   var_x = "pixel",
#' #   var_y = "TIC.norm.tic",
#' #   var_g = "ID",
#' #   y_lim = 40000
#' # )
#'
#' @export
msi_plot_tic <- function(
  df,
  var_x,
  var_y,
  var_g,
  y_lim
) {
  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[[var_x]], # nolint
      y = .data[[var_y]]
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        color = as.factor(.data[[var_g]])
      ),
      shape = 16,
      size = 1,
      alpha = 0.5
    ) +
    ggplot2::geom_smooth(color = "firebrick1") +
    ggplot2::labs(y = "Intensity",
      x = "Pixel"
    ) +
    msi_theme1() + # nolint
    ggplot2::scale_color_manual(values = col_univ()) + # nolint
    ggplot2::scale_y_continuous(limits = c(0, y_lim))
}

#' Ion Image QC
#'
#' Plots the normalized pixel intensities as an ion image for each sample.
#'
#' @param df A normalized data matrix.
#' @param var_y Variable corresponding to the pixel intensity for each sample.
#' @param var_g Sample ID variable.
#' @param perc_int Quantile to be used for the intensity color scale
#' to improve contrast in ion images.
#' @param spl_samp Plot a specific sample
#' (default is to plot all sample images).
#' @param samp_no If spl_samp is TRUE, which sample number should be plotted?
#' @return An ion image visualizing total normalized pixel intensities.
#' @examples
#'
#' # ptic_img <- msi_plot_img(
#' #   df = d_norm[["pixels"]],
#' #   var_y = "TIC.norm.tic",
#' #   var_g = "ID",
#' #   perc_int = 0.95
#' # )
#'
#' @export
msi_plot_img <- function(
  df,
  var_y,
  var_g,
  perc_int,
  spl_samp = FALSE,
  samp_no = NULL
) {
  if(spl_samp == FALSE) { # nolint
    d1 <- df
    d1[is.na(d1[[var_y]])] <- 0
    p1 <- ggplot2::ggplot(
      d1,
      ggplot2::aes(
        x = X, # nolint
        y = Y # nolint
      )
    ) +
      ggplot2::geom_raster(
        ggplot2::aes(
          fill = .data[[var_y]] # nolint
        ),
        interpolate = TRUE
      ) +
      msi_theme2() + # nolint
      ggplot2::labs(fill = "Intensity") +
      ggplot2::scale_y_reverse() +
      ggplot2::scale_fill_gradientn(
        colors = col_grad(), # nolint
        limits = c(
          0,
          quantile(
            d1[[var_y]],
            perc_int
          )
        ),
        na.value = col_grad()[[12]]
      ) +
      ggplot2::facet_wrap(
        . ~ .data[[var_g]],
        ncol = ceiling(length(unique(d1[[var_g]])) / 2)
      )
  }
  if(spl_samp == TRUE) { # nolint
    d1 <- df[df[[var_g]] == samp_no, ]
    d1[is.na(d1[[var_y]])] <- 0
    p1 <- ggplot2::ggplot(
      d1,
      ggplot2::aes(
        x = X, # nolint
        y = Y # nolint
      )
    ) +
      ggplot2::geom_raster(
        ggplot2::aes(
          fill = .data[[var_y]] # nolint
        ),
        interpolate = TRUE
      ) +
      msi_theme2() + # nolint
      ggplot2::labs(fill = "Intensity") +
      ggplot2::scale_y_reverse() +
      ggplot2::scale_fill_gradientn(
        colors = col_grad(), # nolint
        limits = c(
          0,
          quantile(
            d1[[var_y]],
            perc_int
          )
        ),
        na.value = col_grad()[[12]]
      )
  }
  return(p1)
}
