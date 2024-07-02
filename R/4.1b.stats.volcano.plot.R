#' Volcano Plot Function for Univariate Statistical Analysis
#'
#' Uses the EnhancedVolcano package to create a volcano plot from a selection of p-values, fold changes, and compound names.
#'
#' @param df Input data frame containing p-values and group fold changes for each compound.
#' @param group_p P-value variable name (provided as a character string).
#' @param group_fc Fold change variable name (provided as a character string).
#' @param name.lab Vector containing compound names.
#' @param fc.lim Numeric limit of the x-axis.
#' @param title Plot title.
#' @return A volcano plot displaying significantly altered compounds comparing a specific treatment.
#' @examples
#' ## Volcano plot
#' ggsave(
#'   "plotname.png",
#'   MSI.plot.vol(
#'     data.frame(
#'       "Name" = names(d.p2.result[,2:ncol(d.p2.result)]),
#'       "FC" = t(d.fold[1,2:ncol(d.fold)])[,1],
#'       "P" = t(d.p2.result[2,2:ncol(d.p2.result)])[,1]
#'     ),
#'     "FC","P","Name",4,
#'     "Group name"
#'   ),
#'   width = 8,
#'   height = 6,
#'   dpi = 800
#' )
#'
#' @export
  MSI.plot.vol <- function (
    df,
    group_fc,
    group_p,
    name.lab,
    fc.lim,
    title
    ) {
    # Plot function
    v <- EnhancedVolcano(
      df,
      lab = df[[name.lab]],
      title = element_blank(),
      subtitle = element_blank(),
      caption = element_blank(),
      x= group_fc,
      y= group_p,
      pCutoff = 0.05,
      FCcutoff = 0.5,
      cutoffLineType = 'twodash',
      legendLabels = c('NS','Fold Change',
                       'p-value','FC+P'),
      legendLabSize = 12,
      labFace = 'bold',
      col = ggsci::pal_npg("nrc")(10)[c(4,3,5,8)],
      colAlpha = 0.7,
      legendIconSize = 4,
      pointSize = 2,
      border = 'full',
      borderWidth = 1.5,
      legendPosition = 'right',
      labSize = 3,
      drawConnectors = T,
      typeConnectors = "open",
      min.segment.length = grid::unit(1,
                                "mm")
      ) +

      Regio.theme1() +
      ggplot2::labs(color='Key') +
      ggplot2::ggtitle(title) +
      ggplot2::theme(plot.margin = grid::unit(c(.2,.2,.2,.2),"cm")) +
      ggplot2::coord_cartesian(xlim=c(-fc.lim,fc.lim),
                               ylim = c(0,4)) +
      ggplot2::scale_x_continuous(breaks=seq(-4,4,1))

    return(v)

    }


