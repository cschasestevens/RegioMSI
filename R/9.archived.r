# #---- Calculate and visualize group RSDs ----
#
# MSI.calc.RSD <- function(
#     df,
#     md,
#     var.g,
#     df.feat,
#     var.n
#     ) {
#
#   ## calculate group RSDs for each compound
#   q.rsd <- setNames(
#     purrr::reduce(
#       lapply(
#         seq.int(
#           (md + 1),
#           ncol(df),
#           1
#           ),
#         function(x)
#           aggregate(
#             df[
#               df[[x]] > 0,
#               names(df[x])
#               ],
#             list(
#               df[
#                 df[[x]] > 0,
#                 var.g
#                 ]
#               ),
#             function(y)
#               (sd(y)/mean(y))*100)),
#       dplyr::left_join,
#       by = "Group.1"
#       ),
#     c(
#       "Group.RSD",
#       df.feat[[var.n]]
#       )
#     )
#
#   ## format output for plotting
#   p.rsd <- reshape2::melt(
#     q.rsd,
#     id.vars = "Group.RSD"
#     )
#
#   p.rsd[["value"]] <- round(
#     as.numeric(
#       p.rsd[["value"]]
#       ),
#     digits = 2
#     )
#
#   p.rsd <- p.rsd[
#     order(p.rsd[["Group.RSD"]]),
#     ]
#   p.rsd[["ID"]] <- seq.int(
#     1,
#     nrow(p.rsd),
#     1
#     )
#
#   return(
#     list(
#       "all.group.rsd" = p.rsd,
#       "median.rsd" = aggregate(
#         p.rsd[["value"]],
#         list(
#           p.rsd[["Group.RSD"]]
#           ),
#         function(x)
#           median(x)
#           )
#         )
#       )
#   }
#
#
# MSI.plot.RSD <- function(
#     df,
#     var.x,
#     var.rsd,
#     var.g
#     ) {
#
#   ggplot2::ggplot(
#     df,
#     ggplot2::aes(
#       x = .data[[var.x]],
#       y = .data[[var.rsd]]
#     )
#   ) +
#     ggplot2::geom_point(
#       ggplot2::aes(
#         color = .data[[var.g]]),
#       shape=16,
#       size = 1,
#       alpha = 0.5
#     ) +
#     ggplot2::geom_hline(
#       yintercept = max(
#         aggregate(
#           df[["value"]],
#           list(
#             df[[var.g]]
#           ),
#           function(x)
#             median(x)
#         )[["x"]]
#       ),
#       linetype = "dashed",
#       color = "firebrick1") +
#     ggplot2::labs(
#       y = "Group % RSD",
#       x = "Index"
#     ) +
#     Regio.theme1() +
#     ggplot2::scale_color_manual(values = Regio.col.univ())
#
#   }

MSI.image.check <- function(
    df,
    md.df,
    var.int,
    var.g,
    perc.int
  ) {

  d1 <- df
  md1 <- md.df
  d1[is.na(d1[[var.int]])] <- 0

  d1 <- cbind(
    md.df[,c("X","Y")],
    d1
    )

  d2 <- aggregate(
    d1[["X"]],
    list(d1[[var.g]],
         d1[["Y"]]),
    function(z)
      min(z)
  )
  d3 <- aggregate(
    d1[["X"]],
    list(d1[[var.g]],
         d1[["Y"]]),
    function(z)
      max(z)
  )

  d4 <- setNames(
    rbind(
      d2,
      d3[
        sort(d3[["Group.2"]],decreasing = T),
      ],
      d2[1:length(unique(d1[[var.g]])),]
    ),
    c(var.g,"Y","X")
  )

  p1 <- ggplot2::ggplot(
    d1,
    ggplot2::aes(
      x = .data[["X"]],
      y = .data[["Y"]]
    )
  ) +
    ggplot2::geom_raster(
      ggplot2::aes(
        fill = d1[[var.int]]
      ),
      interpolate = T
    ) +
    Regio.theme2() +
    ggplot2::labs(fill = "Intensity") +
    ggplot2::scale_y_reverse(
      limits = c(
        max(d1[["Y"]]
        ),
        1
      )
    ) +
    ggplot2::scale_x_continuous(
      limits = c(
        1,
        max(d1[["Y"]]
        )
      )
    ) +
    # geom_path(
    #   data = d4,
    #   aes(
    #     x = d4[["X"]],
    #     y = d4[["Y"]]
    #     ),
    #   linewidth = 10,
    #   lineend = "round",
    #   linejoin = "round",
    #   color = "white"
    #   ) +
  ggplot2::scale_fill_gradientn(
    colors = Regio.col.grad(),
    limits = c(
      0,
      quantile(d1[[var.int]],
               perc.int
      )
    ),
    na.value = Regio.col.grad()[[12]]
  ) +
    ggplot2::facet_wrap(
      . ~ .data[[var.g]],
      ncol = ceiling(
        length(
          unique(
            d1[[var.g]]
          )
        )/
          2
      )
    )

  return(p1)

}

#### Remove artifact images from dataset ####

## Apply noise/artifact filter
fun.art.filter <- function(md,p.md) {
  if(class(md.samp[["ID"]]) != "factor") {
    md.samp[["ID"]] <- as.factor(md.samp[["ID"]])
    }
  
  d2 <- dplyr::left_join(
    md.samp,
    data.frame(
      d[c(1:p.md)],
      d[,names(d) %in% md.feat[["Input.name"]]]
    ),
    by = c("pixel","ID","Group")
  )
  
  ### Sort colnames and assign final names
  d2 <- data.frame(
    d2[,c(1:md)],
    d2[,sort(names(d2[,(md + 1):ncol(d2)]))]
  )
  
  d2.names <- dplyr::left_join(
    data.frame("Input.name" = names(d2[,(md + 1):ncol(d2)])),
    md.feat,
    by = "Input.name"
  )
  
  names(d2) <- c(names(d2[,c(1:md)]),d2.names[["Name.adduct"]])
  
  return(d2)
  
}

#' Top-5 Cluster Heatmap
#'
#' Plots a heatmap of the top-5 compounds represented within each segmented cluster in a MSI Seurat object.
#'
#' @param so Seurat object containing segmented MSI data.
#' @param marks List of cluster marker compounds to use for subsetting Seurat object.
#' @param var.clus Cluster variable.
#' @param hm.w Heatmap width.
#' @param hm.h Heatmap height.
#' @param fs.c Heatmap column font size.
#' @param fs.r Heatmap row font size.
#' @param sample.no Sample ID number.
#' @return A heatmap displaying the top-5 markers for each cluster.
#' @examples
#' ## Change to lapply if not running on Linux or WSL2
#' d.seg.mark <- setNames(
#' parallel::mclapply(
#'   mc.cores = ceiling(
#'     parallel::detectCores()*
#'       0.5
#'   ),
#'   seq.int(1,length(d.seg),1),
#'   function(x){
#'     # Find marker compounds for each calculated cluster
#'     d.mark <- Seurat::FindAllMarkers(
#'       d.seg[[x]],
#'       test.use = "wilcox",
#'       densify = T,
#'       verbose = T
#'     )
#'     names(d.mark) <- c(
#'       names(
#'         d.mark[
#'           ,
#'           c(1:6)]
#'       ),
#'       "name"
#'     )
#'     d.mark[["ID"]] <- names(d.seg[x])
#'
#'     # Create top-5 heatmap
#'     h.out <- MSI.markers.top5(
#'       # Seurat object
#'       d.seg[[x]],
#'       # Marker gene list
#'       d.mark,
#'       # Cluster column
#'       "Cluster",
#'       # Heatmap width
#'       24,
#'       # Heatmap height
#'       12,
#'       # Column fontsize
#'       6,
#'       # Row fontsize
#'       8,
#'       # Sample number
#'       x
#'     )
#'
#'     return(
#'       list(
#'         "Markers" = d.mark,
#'         "Plot" = h.out
#'       )
#'     )
#'   }
#' ),
#' c(names(d.seg))
#' )
#'
#' @export
MSI.markers.top5 <- function(
    so,
    marks,
    var.clus,
    h.w,
    h.h,
    fs.c,
    fs.r,
    samp.no
) {
  ### Load data
  d.mark <- marks
  d.mark[["cluster.no"]] <- d.mark[["cluster"]]

  ### Top 5 genes per cluster
  d.mark <- dplyr::slice_max(
    dplyr::group_by(
      d.mark,
      .data[["cluster.no"]]),
    order_by = .data[["avg_log2FC"]],
    n = 5
  )[,c(
    "name",
    "cluster"
  )]


  ### Subset seurat and scale
  Seurat::DefaultAssay(so) <- "MSI"

  h <- Seurat::FetchData(
    so,
    vars = c(
      var.clus,
      unique(
        d.mark[["name"]]
      )
    )
  )

  ### Heatmap annotation (average expression)
  h.anno <- as.data.frame(
    lapply(
      h[,2:ncol(
        h
      )],
      function(x)
        mean(x)
    )
  )

  ### Scale and plot average expression per cell type
  h.in <- scale(
    as.matrix(
      magrittr::set_rownames(
        setNames(
          as.data.frame(
            lapply(
              h[,2:ncol(
                h
              )],
              function(x)
                dplyr::select(
                  aggregate(
                    x,
                    list(
                      h[,1]
                    ),
                    FUN = mean
                  ),
                  c(
                    2
                  )
                )
            )
          ),
          names(
            h[,2:ncol(
              h
            )]
          )
        ),
        levels(
          h[,1]
        )
      )
    ),
    center = T
  )

  qs <- quantile(
    h.in,
    probs = c(
      0.05,
      0.95
    )
  )



  fun.hm.col <- circlize::colorRamp2(
    c(
      qs[[1]],
      (qs[[1]])/2,
      (qs[[2]])/2,
      qs[[2]]
    ),
    colors = Regio.col.grad()[c(
      1,3,
      6,12
    )]
  )

  h.out <- ComplexHeatmap::Heatmap(
    h.in,
    col = fun.hm.col,
    name = "Scaled Intensity",
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      `Average.Intensity` = ComplexHeatmap::anno_barplot(
        as.vector(
          t(
            h.anno
          )
        ),
        gp = grid::gpar(
          fill = Regio.col.grad()
        )
      ),
      annotation_name_gp = grid::gpar(
        fontsize = 10
      )
    ),
    show_column_names = T,
    show_row_names = T,
    heatmap_width = grid::unit(
      h.w,"cm"
    ),
    heatmap_height = grid::unit(
      h.h,"cm"
    ),
    column_title = paste("Marker Annotations (Top 5)","Sample",samp.no,sep = " "),
    column_names_rot = 90,
    column_names_gp = grid::gpar(
      fontsize = fs.c
    ),
    row_names_side = "left",
    row_names_gp = grid::gpar(
      fontsize = fs.r
    ),
    cluster_columns = F,
    cluster_rows = F
  )

  return(
    h.out
  )

}