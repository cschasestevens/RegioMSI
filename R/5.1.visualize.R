#' Volcano Plot
#'
#' Generates a volcano plot from a statistical results object for
#' a chosen comparison.
#'
#' @param d_stat Univariate results list generated by ms_stat_anova().
#' @param comp_name Comparison name, provided as a character string.
#' @param diff_col Column name containing effect size.
#' @param p_col Column name containing adjusted p-values.
#' @param p_cut Numeric value for p-value cutoff to indicate significance.
#' @param f_cut Fold change cutoff indicating high effect size.
#' @param f_lim Numeric value for fold change limits on the x-axis.
#' @param y_limit Upper-bound -log10 p-value to display on the y-axis.
#' @param x_title X-axis legend title, provided as a character string.
#' @return A volcano plot for the chosen treatment comparisons.
#' @examples
#'
#' # p_vol <- msi_plot_volcano(
#' #   d_stat = d[["stats"]],
#' #   comp_name = "KO vs. Control",
#' #   f_cut = 0.25,
#' #   f_lim = 6,
#' #   y_limit = 10,
#' #   x_title = "x-axis title"
#' # )
#'
#' @export
msi_plot_volcano <- function(
  d_stat,
  comp_name,
  diff_col = "Log2FC",
  p_col = "p adj",
  p_cut = 0.05,
  f_cut,
  f_lim,
  y_limit,
  x_title
) {
  # Subset Input data
  ld <- d_stat
  ld <- ld[ld[["Comparison"]] == comp_name, ]
  # Plot function
  v <- EnhancedVolcano::EnhancedVolcano(
    ld,
    lab = ld[["Name"]],
    title = ggplot2::element_blank(),
    subtitle = ggplot2::element_blank(),
    caption = ggplot2::element_blank(),
    x = diff_col,
    y = p_col,
    pCutoff = p_cut,
    FCcutoff = f_cut,
    cutoffLineType = "twodash",
    legendLabels = c("NS", "Fold Change",
                     "p-value", "FC+P"),
    legendLabSize = 12,
    labFace = "bold",
    col = ggsci::pal_npg("nrc")(10)[c(4, 3, 5, 8)],
    colAlpha = 0.7,
    legendIconSize = 4,
    pointSize = 2,
    border = "full",
    borderWidth = 1.5,
    legendPosition = "'right'",
    labSize = 3,
    drawConnectors = TRUE,
    max.overlaps = 10,
    typeConnectors = "open",
    min.segment.length = ggplot2::unit(
      1,
      "mm"
    )
  ) +
    msi_theme1() + # nolint
    ggplot2::labs(color = "Key") +
    ggplot2::labs(
      color = "Key",
      x = x_title,
      y = paste("-log10 P-value (FDR)")
    ) +
    ggplot2::theme(plot.margin = ggplot2::unit(c(.2, .2, .2, .2), "cm")) +
    ggplot2::coord_cartesian(xlim = c(-f_lim, f_lim),
                             ylim = c(0, y_limit)) +
    ggplot2::scale_x_continuous(breaks = seq(-f_lim, f_lim, (f_lim / 4)))
  return(v)
}

#' Heatmap 1
#'
#' Generates a heatmap from a statistical results object for
#' a vector of comparisons.
#'
#' @param d_stat Univariate results list generated by ms_stat_anova().
#' @param d_ref Reference annotation list generated by ms_input().
#' @param c_list Comparison list, provided as a vector.
#' @param diff_col Column name containing effect size.
#' @param p_col Column name containing adjusted p-values.
#' @param an1 Variable names used for annotating the heatmap.
#' @param hm_w Heatmap width.
#' @param hm_h Heatmap height.
#' @param fs_r Row fontsize.
#' @param cl_c Logical indicating if the columns should be clustered.
#' @param cl_r Logical indicating if the rows should be clustered.
#' @return A heatmap for the selected comparisons.
#' @examples
#'
#' # msi_plot_heat(
#' #   d_stat = d[["stats"]],
#' #   d_ref = d[["data"]][["anno"]],
#' #   c_list = c("KO vs. Control"),
#' #   diff_col = "Log2FC",
#' #   p_col = "FDR",
#' #   an1 = c("Major.Class", "Saturation"),
#' #   hm_w = 36,
#' #   hm_h = 10,
#' #   fs_r = 8,
#' #   cl_c = TRUE,
#' #   cl_r = TRUE
#' # )
#'
#' @export
msi_plot_heat <- function(
  d_stat, d_ref, c_list, diff_col,
  p_col, an1, hm_w, hm_h, fs_r,
  cl_c, cl_r
) {
  # Load and subset data
  hd1 <- d_stat
  href <- d_ref
  hm_list <- c_list
  an2 <- an1
  hm_in <- hd1[hd1[["Comparison.fc"]] %in% hm_list, ]
  # Format and prepare FC and P matrices
  hd <- reshape2::dcast(
    hm_in[, c(diff_col, "Comparison.fc", "Name")],
    Comparison.fc ~ Name,
    value.var = diff_col
  )
  hp <- reshape2::dcast(
    hm_in[, c(p_col, "Comparison.fc", "Name")],
    Comparison.fc ~ Name,
    value.var = p_col
  )
  hd_in <- as.matrix(hd[, 2:ncol(hd)])
  rownames(hd_in) <- hd[[1]]
  hd_in[is.na(hd_in)] <- max(hd_in[!is.na(hd_in)])
  hp_in <- as.matrix(hp[, 2:ncol(hp)])
  rownames(hp_in) <- hp[[1]]
  hp_in[is.na(hp_in)] <- min(hp_in[!is.na(hp_in)])
  # Figure annotations, colors, and scales
  qs <- quantile(
    hd_in,
    probs = c(
      0.05,
      0.95
    ),
    na.rm = TRUE
  )
  fun_hm_col <- circlize::colorRamp2(
    c(
      qs[[1]],
      (qs[[1]]) / 2,
      (qs[[2]]) / 2,
      qs[[2]]
    ),
    colors = c(
      col_grad()[[1]], # nolint
      col_grad()[[3]],
      "grey",
      col_grad()[[12]]
    )
  )
  if(missing(an2) == FALSE) { # nolint
    fun_hm_bar <- list(
      "Class" = setNames(
        col_univ()[1:length(unique(href[[an2[[1]]]]))], # nolint
        as.character(unique(href[[an2[[1]]]]))
      ),
      "Saturation" = setNames(
        col_univ()[1:length(unique(href[[an2[[2]]]]))], # nolint
        as.character(unique(unique(href[[an2[[2]]]])))
      )
    )
    ### Annotations
    hm_anno_list <- list(
      "column1" = ComplexHeatmap::HeatmapAnnotation(
        `Class` = (href[[an2[[1]]]]),
        `Saturation` = (href[[an2[[2]]]]),
        col = fun_hm_bar,
        show_annotation_name = FALSE
      )
    )
    # Create Plot
    h_out <- ComplexHeatmap::Heatmap(
      hd_in,
      # p-value filter
      layer_fun = function(j, i, x, y, w, h, fill) {
        # map positions of fold change matrix cells onto pvalue heatmap
        p_pos <- ComplexHeatmap::pindex(hp_in, i, j)
        # filter non-significant p values to fill in blanks
        p_blank <- p_pos > 0.05
        # create matrix for subsetting the heatmap
        # to display only significant cells
        p <- ComplexHeatmap::restore_matrix(j, i, x, y) # nolint
        # mask non-significant genes
        grid::grid.rect(
          x[p_blank],
          y[p_blank],
          w,
          h,
          gp = grid::gpar(col = "white")
        )
      },
      col = fun_hm_col,
      top_annotation = hm_anno_list[[1]],
      name = "Log2-Fold Change",
      show_column_names = FALSE,
      show_row_names = TRUE,
      heatmap_width = ggplot2::unit(hm_w, "cm"),
      heatmap_height = ggplot2::unit(hm_h, "cm"),
      row_names_side = "left",
      row_names_gp = grid::gpar(fontsize = fs_r),
      cluster_columns = cl_c,
      cluster_rows = cl_r
    )
  }
  if(missing(an2) == TRUE) { # nolint
    # Create Plot
    h_out <- ComplexHeatmap::Heatmap(
      hd_in,
      # p-value filter
      layer_fun = function(j, i, x, y, w, h, fill) {
        # map positions of fold change matrix cells onto pvalue heatmap
        p_pos <- ComplexHeatmap::pindex(hp_in, i, j)
        # filter non-significant p values to fill in blanks
        p_blank <- p_pos > 0.05
        # create matrix for subsetting the heatmap
        # to display only significant cells
        p <- ComplexHeatmap::restore_matrix(j, i, x, y) # nolint
        # mask non-significant genes
        grid::grid.rect(
          x[p_blank],
          y[p_blank],
          w,
          h,
          gp = grid::gpar(col = "white")
        )
      },
      col = fun_hm_col,
      name = "Log2-Fold Change",
      show_column_names = FALSE,
      show_row_names = TRUE,
      heatmap_width = ggplot2::unit(hm_w, "cm"),
      heatmap_height = ggplot2::unit(hm_h, "cm"),
      row_names_side = "left",
      row_names_gp = grid::gpar(fontsize = fs_r),
      cluster_columns = cl_c,
      cluster_rows = cl_r
    )
  }
  return(h_out)
}
