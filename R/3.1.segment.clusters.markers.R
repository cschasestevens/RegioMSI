#' Seurat-based KNN Clustering and Segmentation
#'
#' Converts MSI data to a Seurat object and performs spatial
#' segmentation using the default implementation of KNN
#' used by Seurat.
#'
#' @param df A normalized data matrix.
#' @param sample_no Sample ID number.
#' @param md Number of metadata columns present in data matrix.
#' @param span loess span parameter
#' (proportion of pixels to fit model to;
#' use lower proportion for tissues with small regions).
#' @param clus_res Resolution to use for clustering
#' (use ?Seurat::FindClusters() documentation for more details;
#' larger number = more clusters).
#' @return A Seurat object containing segmented clusters for a single sample.
#' @examples
#'
#' # test <- msi_segment()
#'
#' @export
msi_segment <- function(
  df,
  sample_no,
  md,
  span = 0.1,
  clus_res = 0.9
) {

  d.seg <- Seurat::CreateSeuratObject(
    counts = t(
      as.matrix(
        d[
          d[["ID"]] == sample.no,
          (md + 1):ncol(d)
        ]
      )
    ),
    meta.data = d[
      d[["ID"]] == sample.no,
      1:md
    ],
    assay = "MSI"
  )

  d.seg <- Seurat::NormalizeData(d.seg)

  d.seg <- Seurat::FindVariableFeatures(
    d.seg,
    loess.span = span
  )

  Seurat::DefaultAssay(
    object = d.seg
  ) <- "MSI"

  d <- Seurat::RunPCA(
    object = Seurat::ScaleData(
      object = d.seg,
      verbose = T
    ),
    verbose = T
  )

  d.seg <- MSI.PCA(d.seg)

  d.seg <- Seurat::FindNeighbors(
    d.seg,
    dims = 1:10
  )
  d.seg <- Seurat::FindClusters(
    d.seg,
    resolution = clus.res,
    n.iter = 25,
    random.seed = 1234
  )

  d.seg <- Seurat::AddMetaData(
    d.seg,
    metadata = as.factor(
      as.numeric(
        d.seg@meta.data$seurat_clusters
      )
    ),
    col.name = "Cluster"
  )

  return(d.seg)

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


#' Composite cluster images for all samples
#'
#' Visualizes combined cluster images for all samples. The input of this function requires a pre-existing normalized data matrix from
#' the previous steps, which is merged with the corresponding cluster column from the segmented Seurat object and indexed by sample ID and
#' pixel number.
#'
#' @param no Numeric sample ID.
#' @param clus1 Name of clustering variable (most often this is 'Cluster').
#' @return A panel of ion images stratified by sample ID for establishing distinct morphological regions in each sample.
#' @examples
#' # Add Seurat clustering column to image data frame and overwrite
#' d <- dplyr::left_join(
#'   dplyr::bind_rows(
#'     lapply(
#'       seq.int(1,length(d.seg),1),
#'       function(x)
#'         data.frame(
#'           "pixel" = d[
#'             d[["ID"]] == as.numeric(names(d.seg[x])),
#'             "pixel"
#'           ],
#'           "Cluster" = as.factor(as.numeric(
#'             d.seg[[x]]@meta.data[["seurat_clusters"]]
#'           )
#'           ))
#'     )
#'   ),
#'   d,
#'   by = "pixel"
#' )
#'
#' d <- dplyr::left_join(
#'   msi.pmd[
#'     ,
#'     c("pixel","X","Y")
#'   ],
#'   d,
#'   by = "pixel"
#' )
#'
#' # Plot cluster images for each sample
#'
#' d.seg.img <- setNames(
#'   lapply(
#'     seq.int(
#'       1,
#'       length(
#'         unique(
#'           d[,c("ID")]
#'         )
#'       ),
#'       1
#'     ),
#'     function(x){
#'       MSI.plot.AllClusters(
#'         # input data
#'         x,
#'         # cluster column
#'         "Cluster"
#'       )
#'     }
#'   ),
#'   c(
#'     unique(
#'       d[
#'         ,
#'         c("ID")
#'       ]
#'     )
#'   )
#' )
#'
#' @export
MSI.plot.AllClusters <- function(no,clus1) {

  ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = d[d[["ID"]] == no,],
      ggplot2::aes(
        x = X,
        y = Y,
        fill = .data[[clus1]]
      ),
      interpolate = T
    ) +
    Regio.theme2() +
    ggplot2::scale_y_reverse(
      limits = c(
        max(d[["Y"]]
        ),
        1
      )
    ) +
    ggplot2::scale_x_continuous(
      limits = c(
        1,
        max(d[["Y"]]
        )
      )
    ) +
    ggplot2::scale_fill_manual(
      name = "Cluster",
      values = Regio.col.univ()
    )
  }


#' Individual cluster images for each sample
#'
#' Visualizes split cluster images for an individual sample. The input of this function requires a pre-existing normalized data matrix from
#' the previous steps, which is merged with the corresponding cluster column from the segmented Seurat object and indexed by sample ID and
#' pixel number.
#'
#' @param df Normalized data matrix containing cluster, sample ID, and pixel coordinate information.
#' @param no.s Numeric sample ID.
#' @param no.c Numeric cluster ID.
#' @return A panel of ion images stratified by cluster ID for establishing distinct morphological regions in each sample.
#' @examples
#' # Split composite cluster images into individual clusters for validation
#' d.seg.img.split <- setNames(
#'   lapply(
#'     seq.int(
#'       1,
#'       length(unique(d[,c("ID")])),
#'       1
#'     ),
#'     function(x) {
#'       setNames(
#'         lapply(
#'           seq.int(
#'             1,
#'             length(
#'               unique(
#'                 d[
#'                   d[["ID"]] == x,
#'                 ]
#'                 [["Cluster"]]
#'               )
#'             ),
#'             1
#'           ),
#'           function(y)
#'             MSI.plot.IndClusters(
#'               # Input data
#'               d,
#'               # Sample number
#'               x,
#'               # Cluster number
#'               y
#'             )
#'         ),
#'         c(
#'           sort(
#'             unique(
#'               d[
#'                 d[["ID"]] == x,
#'               ]
#'               [["Cluster"]]
#'             )
#'           )
#'         )
#'       )
#'     }
#'   ),
#'   c(unique(
#'     d[
#'       ,
#'       c("ID")
#'     ]
#'   )
#'   )
#' )
#'
#' @export
MSI.plot.IndClusters <- function(df,no.s,no.c) {

  ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = df[df[["ID"]] == no.s,],
      ggplot2::aes(
        x = X,
        y = Y,
        fill = factor(
          ifelse(
            df[
              df[["ID"]] == no.s,
            ]
            [["Cluster"]] == no.c,
            no.c,
            (0)
          ),
          levels = c(
            no.c,
            (0)
          )
        )
      ),
      interpolate = T
    ) +
    Regio.theme2() +
    ggplot2::scale_y_reverse(
      limits = c(
        max(df[["Y"]]
        ),
        1
      )
    ) +
    ggplot2::scale_x_continuous(
      limits = c(
        1,
        max(df[["Y"]]
        )
      )
    ) +
    ggplot2::scale_fill_manual(
      name = "Cluster",
      values = c(Regio.col.univ()[[no.c]],"white")
    )

}






