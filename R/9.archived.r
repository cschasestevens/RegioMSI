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