#' Count Pixels
#'
#' Counts the pixel number per each specified variable in a MSI
#' dataset. Can select up to 3 different variables for counting.
#'
#' @param mdf A data frame containing MSI experiment metadata.
#' @param var_sel Names of metadata variables to count.
#' @return A data frame containing the number of pixels per each
#' group specified in var_sel.
#' @examples
#'
#' # cnt_pix <- msi_cnt_pix(
#' #   mdf = d_seg[["pixels"]],
#' #   var_sel = c("Group", "ID", "Cluster")
#' # )
#'
#' @export
msi_cnt_pix <- function(mdf, var_sel) {
  # Load data
  d <- mdf
  v1 <- var_sel
  if(length(v1) == 3) { # nolint
    d1 <- dplyr::count(
      d,
      .data[[v1[[1]]]], # nolint
      .data[[v1[[2]]]],
      .data[[v1[[3]]]]
    )
    d1[["prop"]] <- round(
      d1[["n"]] / nrow(d),
      digits = 3
    )
  }
  if(length(v1) == 2) { # nolint
    d1 <- dplyr::count(
      d,
      .data[[v1[[1]]]], # nolint
      .data[[v1[[2]]]]
    )
    d1[["prop"]] <- round(
      d1[["n"]] / nrow(d),
      digits = 3
    )
  }
  if(length(v1) == 1) { # nolint
    d1 <- dplyr::count(
      d,
      .data[[v1[[1]]]] # nolint
    )
    d1[["prop"]] <- round(
      d1[["n"]] / nrow(d),
      digits = 3
    )
  }
  return(d1)
}

#' Data Transformation and Quality Check
#'
#' Performs data imputation, log2-transformation, and
#' scaling of a MSI dataset. The default option is to use
#' segmented MSI data, stored as a variable named "d_seg."
#' The defaults can be overridden by specifying additional
#' function parameters. See ?RegioMSI::msi_data_check for
#' more details.
#'
#' @param ld A MSI data frame.
#' @return List containing formatted input data and metadata
#' for downstream analysis.
#' @examples
#'
#' # d_seg[["scale.data"]] <- msi_data_check()
#'
#' ## Specify data to use (replace with formatted data frames)
#' # d_seg[["scale.data"]] <- msi_data_check(ld = d_seg[["data"]])
#'
#' @export
msi_data_check <- function(ld = d_seg[["data"]]) { # nolint
  # input data from ms_input()
  ld1 <- list("data" = ld)
  # data imputation
  ## replace zeroes or NA with 1/10th of
  ## the lowest non-zero value present
  ld1[["imputed"]] <- ld1[["data"]]
  ld1[["imputed"]] <- setNames(
    as.data.frame(
      lapply(
        seq.int(1, ncol(ld1[["imputed"]]), 1),
        function(x) {
          ld1[["imputed"]][[x]][is.na(ld1[["imputed"]][[x]])] <- 0
          ld1[["imputed"]][[x]][is.nan(ld1[["imputed"]][[x]])] <- 0
          d1 <- ld1[["imputed"]][[x]]
          d1 <- ifelse(
            d1 == 0,
            round(0.1 * min(d1[d1 > 0]), digits = 2),
            round(d1, digits = 2)
          )
          return(d1)
        }
      )
    ),
    c(names(ld1[["imputed"]]))
  )
  # Log2-transformation
  ld1[["data.log2"]] <- setNames(
    as.data.frame(
      lapply(
        seq.int(
          1,
          ncol(ld1[["imputed"]]),
          1
        ),
        function(x) {
          log2(ld1[["imputed"]][[x]])
        }
      )
    ),
    c(names(ld1[["imputed"]]))
  )
  # Pareto scaling
  ld1[["data.pareto"]] <- setNames(
    as.data.frame(
      apply(
        apply(
          as.matrix(ld1[["data.log2"]]),
          2,
          function(x) x - mean(x)
        ),
        2,
        function(y) round(y / sqrt(sd(y)), digits = 2)
      )
    ),
    c(names(ld1[["data.log2"]]))
  )
  # Data distributions
  dst <- function(df, n1) {
    dsty <- as.data.frame(colMeans(df))
    names(dsty) <- c("Value")
    p <- ggplot2::ggplot(
      dsty,
      ggplot2::aes(x = Value) # nolint 
    ) +
      # Density Plot
      ggplot2::geom_density(
        color = "darkslategrey",
        fill = col_univ()[[1]] # nolint
      ) +
      # Plot Theme
      ggplot2::labs(
        title = paste(n1, "Data Distribution"),
        y = "Density"
      ) +
      msi_theme1() + # nolint
      ggplot2::theme(
        plot.margin = ggplot2::unit(c(rep(0.5, 4)), "cm")
      )
    return(p)
  }
  lp1 <- list(
    "plot.dist" = ggpubr::ggarrange(
      dst(ld1[["data"]], "Input"),
      dst(ld1[["data.log2"]], "Log2"),
      dst(ld1[["data.pareto"]], "Pareto Scaled"),
      nrow = 1,
      ncol = 3
    )
  )
  return(
    list(
      "data" = ld1,
      "plots" = lp1
    )
  )
}

#' Group-wise Fold Change Function
#'
#' Calculates the log2 fold changes for all
#' possible combinations of selected metadata variables.
#'
#' @param mat1 MSI data frame.
#' @param md MSI metadata, formatted as a data frame.
#' @param md_var Vector of metadata variables for calculating fold change.
#' @param an MSI feature metadata.
#' @param an_name Variable name containing compound names
#' (default value is "Name").
#' @param fc_class Logical indicating whether class-wide
#' fold changes should be calculated. Requires MSI feature metadata if TRUE.
#' @param grp_class Name of the class annotation column
#' to use. Ignored when fc_class is FALSE.
#' @param core_perc Not used on Windows; Proportion of cores to use for
#' conducting statistical analysis in parallel.
#' @return A data frame containing the group-wise
#' fold changes for each compound.
#' @examples
#'
#' ## Calculate individual fold changes
#' # msi_stat_fc(
#' #   mat1 = d_seg[["data"]],
#' #   md = d_seg[["pixels"]],
#' #   md_var = c("Group"),
#' #   an = d_seg[["features"]]
#' # )
#'
#' ## Calculate class-based fold changes
#' # msi_stat_fc(
#' #   mat1 = d_seg[["data"]],
#' #   md = d_seg[["pixels"]],
#' #   md_var = c("Group"),
#' #   an = d_seg[["features"]],
#' #   fc_class = TRUE,
#' #   grp_class = "label.saturation"
#' # )
#'
#' @export
msi_stat_fc <- function(
  mat1,
  md,
  md_var,
  an,
  an_name = "Name",
  fc_class = FALSE,
  grp_class = NULL,
  core_perc = 0.75
) {
  if(fc_class == FALSE) { # nolint
    d1 <- mat1
    md1 <- md
    an1 <- an
    an2 <- an_name
    # Group Combinations
    var_comb <- dplyr::bind_rows(lapply(
      seq.int(0, length(md_var) - 1, 1),
      function(x) {
        cb1 <- combn(md_var, x + 1)
        cb1 <- as.data.frame(t(
          dplyr::bind_cols(
            lapply(
              as.data.frame(cb1),
              function(y) paste(y, collapse = ":")
            )
          )
        ))
        return(cb1)
      }
    ))

    fold_comb <- dplyr::bind_rows(
      lapply(
        seq.int(1, nrow(var_comb), 1),
        function(x) {
          cb1 <- paste(unlist(strsplit(var_comb[x, 1], ":")), sep = ", ")
          # Change variable to factor
          cb2 <- setNames(
            as.data.frame(
              lapply(
                cb1,
                function(y) as.character(md1[, y])
              )
            ),
            c(cb1)
          )
          cb2 <- unique(cb2)
          # Combine columns
          cb2[["comb"]] <- unlist(
            lapply(
              as.data.frame(t(cb2)),
              function(z) paste(z, collapse = ":")
            )
          )
          # Determine combinations
          cb2 <- as.data.frame(
            t(
              unique(
                combn(
                  sort(cb2[["comb"]], decreasing = TRUE),
                  2
                )
              )
            )
          )
          return(cb2)
        }
      )
    )
    # Metabolite Group Means (excluding missing samples)
    fold_mean <- dplyr::bind_rows(
      setNames(
        parallel::mclapply(
          mc.cores = ceiling(parallel::detectCores() * core_perc),
          seq.int(1, nrow(var_comb), 1),
          function(x) {
            cm1 <- paste(unlist(strsplit(var_comb[x, 1], ":")), sep = ", ")
            cm2 <- aggregate(
              d1,
              lapply(cm1, function(z) as.character(md1[[z]])),
              function(y) mean(y)
            )
            cm2 <- setNames(
              data.frame(
                "Group" = unlist(
                  lapply(
                    as.data.frame(t(cm2[, 1:length(cm1)])), # nolint
                    function(z) paste(z, collapse = ":")
                  )
                ),
                cm2[, (length(cm1) + 1):ncol(cm2)]
              ),
              c("Group", names(cm2[, (length(cm1) + 1):ncol(cm2)]))
            )
            return(cm2)
          }
        ),
        var_comb[[1]]
      )
    )

    # Group Fold Changes
    fun.comb <- function(x, y) {x / y} # nolint
    fold_all <- setNames(
      reshape2::melt(
        dplyr::mutate(dplyr::bind_cols(
          lapply(
            seq.int(1, nrow(fold_comb), 1),
            function(x) {
              fm1 <- t(
                fold_mean[fold_mean[[1]] == fold_comb[x, 1], 2:ncol(fold_mean)]
              )
              fm2 <- t(
                fold_mean[fold_mean[[1]] == fold_comb[x, 2], 2:ncol(fold_mean)]
              )
              fc1 <- setNames(
                as.data.frame(fun.comb(fm1, fm2)),
                paste(fold_comb[x, 1], fold_comb[x, 2], sep = "-")
              )
              fc1 <- log2(fc1)
              return(fc1)
            }
          )
        ), "Name" = an1[[an2]]), id.vars = an2
      ),
      c("Name", "Comparison.fc", "Log2FC")
    )
    fold_all[["Comparison.fc"]] <- as.character(fold_all[["Comparison.fc"]])
  }



  if(fc_class == TRUE) { # nolint
    d1 <- mat1
    an1 <- an
    d1 <- aggregate(
      t(as.matrix(d1)),
      list(an1[[grp_class]]),
      function(x) mean(x)
    )
    d1 <- setNames(as.data.frame(t(d1[, 2:ncol(d1)])), d1[[1]])
    md1 <- md

    # Group Combinations
    var_comb <- dplyr::bind_rows(lapply(
      seq.int(0, length(md_var) - 1, 1),
      function(x) {
        cb1 <- combn(md_var, x + 1)
        cb1 <- as.data.frame(t(
          dplyr::bind_cols(
            lapply(
              as.data.frame(cb1),
              function(y) paste(y, collapse = ":")
            )
          )
        ))
        return(cb1)
      }
    ))

    fold_comb <- dplyr::bind_rows(
      lapply(
        seq.int(1, nrow(var_comb), 1),
        function(x) {
          cb1 <- paste(unlist(strsplit(var_comb[x, 1], ":")), sep = ", ")
          # Change variable to factor
          cb2 <- setNames(
            as.data.frame(
              lapply(
                cb1,
                function(y) as.character(md1[, y])
              )
            ),
            c(cb1)
          )
          cb2 <- unique(cb2)
          # Combine columns
          cb2[["comb"]] <- unlist(
            lapply(
              as.data.frame(t(cb2)),
              function(z) paste(z, collapse = ":")
            )
          )
          # Determine combinations
          cb2 <- as.data.frame(
            t(
              unique(
                combn(
                  sort(cb2[["comb"]], decreasing = TRUE),
                  2
                )
              )
            )
          )
          return(cb2)
        }
      )
    )
    # Metabolite Group Means (excluding missing samples)
    fold_mean <- dplyr::bind_rows(
      setNames(
        parallel::mclapply(
          mc.cores = ceiling(parallel::detectCores() * core_perc),
          seq.int(1, nrow(var_comb), 1),
          function(x) {
            cm1 <- paste(unlist(strsplit(var_comb[x, 1], ":")), sep = ", ")
            cm2 <- aggregate(
              d1,
              lapply(cm1, function(z) as.character(md1[[z]])),
              function(y) mean(y)
            )
            cm2 <- setNames(
              data.frame(
                "Group" = unlist(
                  lapply(
                    as.data.frame(t(cm2[, 1:length(cm1)])), # nolint
                    function(z) paste(z, collapse = ":")
                  )
                ),
                cm2[, (length(cm1) + 1):ncol(cm2)]
              ),
              c("Group", names(cm2[, (length(cm1) + 1):ncol(cm2)]))
            )
            return(cm2)
          }
        ),
        var_comb[[1]]
      )
    )

    # Group Fold Changes
    fun.comb <- function(x, y) {x / y} # nolint
    fold_all <- setNames(
      reshape2::melt(
        dplyr::mutate(dplyr::bind_cols(
          lapply(
            seq.int(1, nrow(fold_comb), 1),
            function(x) {
              fm1 <- t(
                fold_mean[fold_mean[[1]] == fold_comb[x, 1], 2:ncol(fold_mean)]
              )
              fm2 <- t(
                fold_mean[fold_mean[[1]] == fold_comb[x, 2], 2:ncol(fold_mean)]
              )
              fc1 <- setNames(
                as.data.frame(fun.comb(fm1, fm2)),
                paste(fold_comb[x, 1], fold_comb[x, 2], sep = "-")
              )
              fc1 <- log2(fc1)
              return(fc1)
            }
          )
        ), "Name" = names(d1)), id.vars = an2
      ),
      c("Name", "Comparison.fc", "Log2FC")
    )
    fold_all[["Comparison.fc"]] <- as.character(fold_all[["Comparison.fc"]])
  }

  return(fold_all)
}

#' Multivariate ANOVA
#'
#' Conducts an ANOVA for selected variables with Tukey's post-hoc
#' analysis for each compound and combines with fold changes
#' for each group comparison. The default option is to use
#' segmented MSI data, stored as a variable named "d_seg."
#' The defaults can be overridden by specifying additional
#' function parameters. See ?RegioMSI::msi_data_check for
#' more details.
#'
#' @param matpv Data matrix for conducting ANOVA.
#' @param matfc Data matrix for calculating fold change.
#' @param md Metadata from ms_input().
#' @param md_var Metadata variables to test by ANOVA.
#' @param an Reference annotation data frame.
#' @param an_name Variable name containing compound names
#' (default value is "Name").
#' @param fc_class1 Logical indicating if compound classes rather
#' than individual compounds should be tested for significant differences.
#' @param grp_class1 If fc_class1 is TRUE, which variable should be used
#' to group individual compounds into separate compound classes?
#' @param core_perc Not used on Windows; Proportion of cores to use for
#' conducting statistical analysis in parallel.
#' @return A data frame containing the combined ANOVA and fold change results
#' for the specified variable(s).
#' @examples
#'
#' ## Specify data
#' # d_stat <- msi_stat_anova(
#' #   matpv = d[["data"]][["data.pareto"]],
#' #   matfc = d[["data"]][["imputed"]],
#' #   md = d[["data"]][["meta"]],
#' #   md_var = c("Group", "Cluster"),
#' #   an = d[["data"]][["anno"]],
#' #   fc_class1 = FALSE,
#' #   grp_class1 = "Subclass"
#' # )
#'
#' ## or use default values
#' # d_stat <- msi_stat_anova(md_var = c("Group", "Cluster"))
#'
#' @export
msi_stat_anova <- function( # nolint
  matpv = d_seg[["scale.data"]][["data"]][["data.pareto"]], # nolint
  matfc = d_seg[["data"]],
  md = d_seg[["pixels"]],
  md_var,
  an = d_seg[["features"]],
  an_name = "Name",
  fc_class1 = FALSE,
  grp_class1 = NULL,
  core_perc = 0.75
) { # nolint
  if(Sys.info()[["sysname"]] != "Windows") { # nolint
    if(fc_class1 == FALSE){ # nolint
      # Determine interactions between selected variables
      # Load data
      d1 <- matpv
      md1 <- md
      an1 <- an
      an2 <- an_name
      ## ANOVA
      if(length(md_var) == 3) { # nolint
        d_mult <- dplyr::bind_rows(
          setNames(
            parallel::mclapply(
              mc.cores = ceiling(parallel::detectCores() * core_perc),
              seq.int(1, ncol(d1), 1),
              function(y) {
                d2 <- dplyr::bind_rows(
                  lapply(
                    seq.int(1, 7, 1),
                    function(x) {
                      as.data.frame(
                        as.matrix(
                          TukeyHSD(
                            aov(
                              d1[[y]] ~
                                as.factor(md1[[md_var[[1]]]]) *
                                  as.factor(md1[[md_var[[2]]]]) *
                                  as.factor(md1[[md_var[[3]]]])
                            )
                          )[[x]]
                        )
                      )
                    }
                  )
                )
                d2 <- dplyr::mutate(
                  d2,
                  "Name" = an1[[an2]][[y]],
                  "Comparison" = rownames(d2)
                )
                d2 <- d2[, c("Name", "Comparison", "p adj")]
                d2[["FDR"]] <- p.adjust(d2[["p adj"]], method = "fdr")
                return(d2)
              }
            ),
            an1[[an2]]
          )
        )
      }
      if(length(md_var) == 2) { # nolint
        d_mult <- dplyr::bind_rows(
          setNames(
            parallel::mclapply(
              mc.cores = ceiling(parallel::detectCores() * core_perc),
              seq.int(1, ncol(d1), 1),
              function(y) {
                d2 <- dplyr::bind_rows(
                  lapply(
                    seq.int(1, 3, 1),
                    function(x) {
                      as.data.frame(
                        as.matrix(
                          TukeyHSD(
                            aov(
                              d1[[y]] ~
                                as.factor(md1[[md_var[[1]]]]) *
                                  as.factor(md1[[md_var[[2]]]])
                            )
                          )[[x]]
                        )
                      )
                    }
                  )
                )
                d2 <- dplyr::mutate(
                  d2,
                  "Name" = an1[[an2]][[y]],
                  "Comparison" = rownames(d2)
                )
                d2 <- d2[, c("Name", "Comparison", "p adj")]
                d2[["FDR"]] <- p.adjust(d2[["p adj"]], method = "fdr")
                return(d2)
              }
            ),
            an1[[an2]]
          )
        )
      }
      if(length(md_var) == 1) { # nolint
        d_mult <- dplyr::bind_rows(
          setNames(
            parallel::mclapply(
              mc.cores = ceiling(parallel::detectCores() * core_perc),
              seq.int(1, ncol(d1), 1),
              function(y) {
                d2 <- dplyr::bind_rows(
                  lapply(
                    seq.int(1, 1, 1),
                    function(x) {
                      as.data.frame(
                        as.matrix(
                          TukeyHSD(
                            aov(
                              d1[[y]] ~
                                as.factor(md1[[md_var[[1]]]])
                            )
                          )[[x]]
                        )
                      )
                    }
                  )
                )
                d2 <- dplyr::mutate(
                  d2,
                  "Name" = an1[[an2]][[y]],
                  "Comparison" = rownames(d2)
                )
                d2 <- d2[, c("Name", "Comparison", "p adj")]
                d2[["FDR"]] <- p.adjust(d2[["p adj"]], method = "fdr")
                return(d2)
              }
            ),
            an1[[an2]]
          )
        )
      }

      # Add fold change
      fold_all <- msi_stat_fc( # nolint
        mat1 = matfc,
        md = md1,
        md_var = md_var,
        an = an1,
        an_name = an2,
        fc_class = fc_class1,
        grp_class = grp_class1
      )

      # Fix anova comparisons to match fold change comparisons
      dm_fix <- data.frame(
        "Comparison" = unique(d_mult[["Comparison"]]),
        "Comparison.fc" = unlist(
          lapply(
            seq.int(1, length(unique(d_mult[["Comparison"]])), 1),
            function(x) {
              fx1 <- unique(d_mult[["Comparison"]])[[x]]
              fx2 <- ifelse(
                fx1 %in% unique(fold_all[["Comparison.fc"]]) == FALSE,
                paste(
                  unlist(strsplit(fx1, "-"))[[2]],
                  unlist(strsplit(fx1, "-"))[[1]],
                  sep = "-"
                ),
                fx1
              )
              return(fx2)
            }
          )
        )
      )

      d_out <- dplyr::left_join(
        dplyr::left_join(
          d_mult,
          dm_fix,
          by = "Comparison"
        ),
        fold_all,
        by = c("Comparison.fc", "Name")
      )
      write.table(
        d_out,
        "analysis/table.stats.txt",
        col.names = TRUE,
        row.names = FALSE,
        sep = "\t"
      )
    }

    if(fc_class1 == TRUE){ # nolint
      # Determine interactions between selected variables
      # Load data
      d1 <- matpv
      md1 <- md
      an1 <- an
      an2 <- an_name
      d1 <- aggregate(
        t(as.matrix(d1)),
        list(an1[[grp_class1]]),
        function(x) mean(x)
      )
      d1 <- setNames(as.data.frame(t(d1[, 2:ncol(d1)])), d1[[1]])
      ## ANOVA
      if(length(md_var) == 3) { # nolint
        d_mult <- dplyr::bind_rows(
          setNames(
            parallel::mclapply(
              mc.cores = ceiling(parallel::detectCores() * core_perc),
              seq.int(1, ncol(d1), 1),
              function(y) {
                d2 <- dplyr::bind_rows(
                  lapply(
                    seq.int(1, 7, 1),
                    function(x) {
                      as.data.frame(
                        as.matrix(
                          TukeyHSD(
                            aov(
                              d1[[y]] ~
                                as.factor(md1[[md_var[[1]]]]) *
                                  as.factor(md1[[md_var[[2]]]]) *
                                  as.factor(md1[[md_var[[3]]]])
                            )
                          )[[x]]
                        )
                      )
                    }
                  )
                )
                d2 <- dplyr::mutate(
                  d2,
                  "Name" = names(d1)[[y]],
                  "Comparison" = rownames(d2)
                )
                d2 <- d2[, c("Name", "Comparison", "p adj")]
                d2[["FDR"]] <- p.adjust(d2[["p adj"]], method = "fdr")
                return(d2)
              }
            ),
            names(d1)
          )
        )
      }
      if(length(md_var) == 2) { # nolint
        d_mult <- dplyr::bind_rows(
          setNames(
            parallel::mclapply(
              mc.cores = ceiling(parallel::detectCores() * core_perc),
              seq.int(1, ncol(d1), 1),
              function(y) {
                d2 <- dplyr::bind_rows(
                  lapply(
                    seq.int(1, 3, 1),
                    function(x) {
                      as.data.frame(
                        as.matrix(
                          TukeyHSD(
                            aov(
                              d1[[y]] ~
                                as.factor(md1[[md_var[[1]]]]) *
                                  as.factor(md1[[md_var[[2]]]])
                            )
                          )[[x]]
                        )
                      )
                    }
                  )
                )
                d2 <- dplyr::mutate(
                  d2,
                  "Name" = names(d1)[[y]],
                  "Comparison" = rownames(d2)
                )
                d2 <- d2[, c("Name", "Comparison", "p adj")]
                d2[["FDR"]] <- p.adjust(d2[["p adj"]], method = "fdr")
                return(d2)
              }
            ),
            names(d1)
          )
        )
      }
      if(length(md_var) == 1) { # nolint
        d_mult <- dplyr::bind_rows(
          setNames(
            parallel::mclapply(
              mc.cores = ceiling(parallel::detectCores() * core_perc),
              seq.int(1, ncol(d1), 1),
              function(y) {
                d2 <- dplyr::bind_rows(
                  lapply(
                    seq.int(1, 1, 1),
                    function(x) {
                      as.data.frame(
                        as.matrix(
                          TukeyHSD(
                            aov(
                              d1[[y]] ~
                                as.factor(md1[[md_var[[1]]]])
                            )
                          )[[x]]
                        )
                      )
                    }
                  )
                )
                d2 <- dplyr::mutate(
                  d2,
                  "Name" = names(d1)[[y]],
                  "Comparison" = rownames(d2)
                )
                d2 <- d2[, c("Name", "Comparison", "p adj")]
                d2[["FDR"]] <- p.adjust(d2[["p adj"]], method = "fdr")
                return(d2)
              }
            ),
            names(d1)
          )
        )
      }

      # Add fold change
      fold_all <- msi_stat_fc( # nolint
        mat1 = matfc,
        md = md1,
        md_var = md_var,
        an = an1,
        an_name = an2,
        fc_class = fc_class1,
        grp_class = grp_class1
      )

      # Fix anova comparisons to match fold change comparisons
      dm_fix <- data.frame(
        "Comparison" = unique(d_mult[["Comparison"]]),
        "Comparison.fc" = unlist(
          lapply(
            seq.int(1, length(unique(d_mult[["Comparison"]])), 1),
            function(x) {
              fx1 <- unique(d_mult[["Comparison"]])[[x]]
              fx2 <- ifelse(
                fx1 %in% unique(fold_all[["Comparison.fc"]]) == FALSE,
                paste(
                  unlist(strsplit(fx1, "-"))[[2]],
                  unlist(strsplit(fx1, "-"))[[1]],
                  sep = "-"
                ),
                fx1
              )
              return(fx2)
            }
          )
        )
      )

      d_out <- dplyr::left_join(
        dplyr::left_join(
          d_mult,
          dm_fix,
          by = "Comparison"
        ),
        fold_all,
        by = c("Comparison.fc", "Name")
      )
      write.table(
        d_out,
        "analysis/table.stats.class.txt",
        col.names = TRUE,
        row.names = FALSE,
        sep = "\t"
      )
    }
  }
  if(Sys.info()[["sysname"]] == "Windows") { # nolint
    if(fc_class1 == FALSE){ # nolint
      # Determine interactions between selected variables
      # Load data
      d1 <- matpv
      md1 <- md
      an1 <- an
      an2 <- an_name
      ## ANOVA
      if(length(md_var) == 3) { # nolint
        d_mult <- dplyr::bind_rows(
          setNames(
            lapply(
              seq.int(1, ncol(d1), 1),
              function(y) {
                d2 <- dplyr::bind_rows(
                  lapply(
                    seq.int(1, 7, 1),
                    function(x) {
                      as.data.frame(
                        as.matrix(
                          TukeyHSD(
                            aov(
                              d1[[y]] ~
                                as.factor(md1[[md_var[[1]]]]) *
                                  as.factor(md1[[md_var[[2]]]]) *
                                  as.factor(md1[[md_var[[3]]]])
                            )
                          )[[x]]
                        )
                      )
                    }
                  )
                )
                d2 <- dplyr::mutate(
                  d2,
                  "Name" = an1[[an2]][[y]],
                  "Comparison" = rownames(d2)
                )
                d2 <- d2[, c("Name", "Comparison", "p adj")]
                d2[["FDR"]] <- p.adjust(d2[["p adj"]], method = "fdr")
                return(d2)
              }
            ),
            an1[[an2]]
          )
        )
      }
      if(length(md_var) == 2) { # nolint
        d_mult <- dplyr::bind_rows(
          setNames(
            lapply(
              seq.int(1, ncol(d1), 1),
              function(y) {
                d2 <- dplyr::bind_rows(
                  lapply(
                    seq.int(1, 3, 1),
                    function(x) {
                      as.data.frame(
                        as.matrix(
                          TukeyHSD(
                            aov(
                              d1[[y]] ~
                                as.factor(md1[[md_var[[1]]]]) *
                                  as.factor(md1[[md_var[[2]]]])
                            )
                          )[[x]]
                        )
                      )
                    }
                  )
                )
                d2 <- dplyr::mutate(
                  d2,
                  "Name" = an1[[an2]][[y]],
                  "Comparison" = rownames(d2)
                )
                d2 <- d2[, c("Name", "Comparison", "p adj")]
                d2[["FDR"]] <- p.adjust(d2[["p adj"]], method = "fdr")
                return(d2)
              }
            ),
            an1[[an2]]
          )
        )
      }
      if(length(md_var) == 1) { # nolint
        d_mult <- dplyr::bind_rows(
          setNames(
            lapply(
              seq.int(1, ncol(d1), 1),
              function(y) {
                d2 <- dplyr::bind_rows(
                  lapply(
                    seq.int(1, 1, 1),
                    function(x) {
                      as.data.frame(
                        as.matrix(
                          TukeyHSD(
                            aov(
                              d1[[y]] ~
                                as.factor(md1[[md_var[[1]]]])
                            )
                          )[[x]]
                        )
                      )
                    }
                  )
                )
                d2 <- dplyr::mutate(
                  d2,
                  "Name" = an1[[an2]][[y]],
                  "Comparison" = rownames(d2)
                )
                d2 <- d2[, c("Name", "Comparison", "p adj")]
                d2[["FDR"]] <- p.adjust(d2[["p adj"]], method = "fdr")
                return(d2)
              }
            ),
            an1[[an2]]
          )
        )
      }

      # Add fold change
      fold_all <- msi_stat_fc( # nolint
        mat1 = matfc,
        md = md1,
        md_var = md_var,
        an = an1,
        an_name = an2,
        fc_class = fc_class1,
        grp_class = grp_class1
      )

      # Fix anova comparisons to match fold change comparisons
      dm_fix <- data.frame(
        "Comparison" = unique(d_mult[["Comparison"]]),
        "Comparison.fc" = unlist(
          lapply(
            seq.int(1, length(unique(d_mult[["Comparison"]])), 1),
            function(x) {
              fx1 <- unique(d_mult[["Comparison"]])[[x]]
              fx2 <- ifelse(
                fx1 %in% unique(fold_all[["Comparison.fc"]]) == FALSE,
                paste(
                  unlist(strsplit(fx1, "-"))[[2]],
                  unlist(strsplit(fx1, "-"))[[1]],
                  sep = "-"
                ),
                fx1
              )
              return(fx2)
            }
          )
        )
      )

      d_out <- dplyr::left_join(
        dplyr::left_join(
          d_mult,
          dm_fix,
          by = "Comparison"
        ),
        fold_all,
        by = c("Comparison.fc", "Name")
      )
      head(d_out)
      write.table(
        d_out,
        "analysis/table.stats.txt",
        col.names = TRUE,
        row.names = FALSE,
        sep = "\t"
      )
    }

    if(fc_class1 == TRUE){ # nolint
      # Determine interactions between selected variables
      # Load data
      d1 <- matpv
      md1 <- md
      an1 <- an
      an2 <- an_name
      d1 <- aggregate(
        t(as.matrix(d1)),
        list(an1[[grp_class1]]),
        function(x) mean(x)
      )
      d1 <- setNames(as.data.frame(t(d1[, 2:ncol(d1)])), d1[[1]])
      ## ANOVA
      if(length(md_var) == 3) { # nolint
        d_mult <- dplyr::bind_rows(
          setNames(
            lapply(
              seq.int(1, ncol(d1), 1),
              function(y) {
                d2 <- dplyr::bind_rows(
                  lapply(
                    seq.int(1, 7, 1),
                    function(x) {
                      as.data.frame(
                        as.matrix(
                          TukeyHSD(
                            aov(
                              d1[[y]] ~
                                as.factor(md1[[md_var[[1]]]]) *
                                  as.factor(md1[[md_var[[2]]]]) *
                                  as.factor(md1[[md_var[[3]]]])
                            )
                          )[[x]]
                        )
                      )
                    }
                  )
                )
                d2 <- dplyr::mutate(
                  d2,
                  "Name" = names(d1)[[y]],
                  "Comparison" = rownames(d2)
                )
                d2 <- d2[, c("Name", "Comparison", "p adj")]
                d2[["FDR"]] <- p.adjust(d2[["p adj"]], method = "fdr")
                return(d2)
              }
            ),
            names(d1)
          )
        )
      }
      if(length(md_var) == 2) { # nolint
        d_mult <- dplyr::bind_rows(
          setNames(
            lapply(
              seq.int(1, ncol(d1), 1),
              function(y) {
                d2 <- dplyr::bind_rows(
                  lapply(
                    seq.int(1, 3, 1),
                    function(x) {
                      as.data.frame(
                        as.matrix(
                          TukeyHSD(
                            aov(
                              d1[[y]] ~
                                as.factor(md1[[md_var[[1]]]]) *
                                  as.factor(md1[[md_var[[2]]]])
                            )
                          )[[x]]
                        )
                      )
                    }
                  )
                )
                d2 <- dplyr::mutate(
                  d2,
                  "Name" = names(d1)[[y]],
                  "Comparison" = rownames(d2)
                )
                d2 <- d2[, c("Name", "Comparison", "p adj")]
                d2[["FDR"]] <- p.adjust(d2[["p adj"]], method = "fdr")
                return(d2)
              }
            ),
            names(d1)
          )
        )
      }
      if(length(md_var) == 1) { # nolint
        d_mult <- dplyr::bind_rows(
          setNames(
            lapply(
              seq.int(1, ncol(d1), 1),
              function(y) {
                d2 <- dplyr::bind_rows(
                  lapply(
                    seq.int(1, 1, 1),
                    function(x) {
                      as.data.frame(
                        as.matrix(
                          TukeyHSD(
                            aov(
                              d1[[y]] ~
                                as.factor(md1[[md_var[[1]]]])
                            )
                          )[[x]]
                        )
                      )
                    }
                  )
                )
                d2 <- dplyr::mutate(
                  d2,
                  "Name" = names(d1)[[y]],
                  "Comparison" = rownames(d2)
                )
                d2 <- d2[, c("Name", "Comparison", "p adj")]
                d2[["FDR"]] <- p.adjust(d2[["p adj"]], method = "fdr")
                return(d2)
              }
            ),
            names(d1)
          )
        )
      }

      # Add fold change
      fold_all <- msi_stat_fc( # nolint
        mat1 = matfc,
        md = md1,
        md_var = md_var,
        an = an1,
        an_name = an2,
        fc_class = fc_class1,
        grp_class = grp_class1
      )

      # Fix anova comparisons to match fold change comparisons
      dm_fix <- data.frame(
        "Comparison" = unique(d_mult[["Comparison"]]),
        "Comparison.fc" = unlist(
          lapply(
            seq.int(1, length(unique(d_mult[["Comparison"]])), 1),
            function(x) {
              fx1 <- unique(d_mult[["Comparison"]])[[x]]
              fx2 <- ifelse(
                fx1 %in% unique(fold_all[["Comparison.fc"]]) == FALSE,
                paste(
                  unlist(strsplit(fx1, "-"))[[2]],
                  unlist(strsplit(fx1, "-"))[[1]],
                  sep = "-"
                ),
                fx1
              )
              return(fx2)
            }
          )
        )
      )
      d_out <- dplyr::left_join(
        dplyr::left_join(
          d_mult,
          dm_fix,
          by = "Comparison"
        ),
        fold_all,
        by = c("Comparison.fc", "Name")
      )
      write.table(
        d_out,
        "analysis/table.stats.class.txt",
        col.names = TRUE,
        row.names = FALSE,
        sep = "\t"
      )
    }
  }
  return(d_out)
}

#' Student's T-test of MSI Data
#'
#' Conducts a Student's T-test for each compound.
#'
#' @param df Normalized data matrix containing cluster, sample ID,
#' pixel coordinate, and pixel intensity information.
#' @param md Number of metadata columns present in input data.
#' @param var_g Treatment group column name.
#' @param g1 Treatment group 1 name.
#' @param g2 Treatment group 2 name.
#' @param l_names Vector of compound names equal to the number
#' of compounds present in the given data matrix.
#' @return A data frame containing T-test results for each compound.
#' @examples
#'
#' # d_stat_t <- msi_stat_t(
#' #   df = d.p2,
#' #   md = 2,
#' #   var_g = "Group",
#' #   g1 = "F.HO3",
#' #   g2 = "F.SFA",
#' #   l_names = an1[["Name"]]
#' # )
#'
#' @export
msi_stat_t <- function(
  df,
  md,
  var_g,
  g1,
  g2,
  l_names
) {
  setNames(
    as.data.frame(
      lapply(
        names(
          df[
            ,
            (md + 1):ncol(df)
          ]
        ),
        function(x) {
          t.test(
            df[
              df[[var_g]] == g1 |
                df[[var_g]] == g2,
              x
            ] ~
              df[df[[var_g]] == g1 |
                  df[[var_g]] == g2,
                var_g
              ]
          )[["p.value"]]
        }
      )
    ),
    c(l_names)
  )
}
