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