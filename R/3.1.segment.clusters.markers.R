#' Seurat-based KNN Clustering and Segmentation
#'
#' Converts MSI data to a Seurat object and performs spatial
#' segmentation using the default implementation of KNN
#' used by Seurat.
#'
#' @param df A normalized data matrix.
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
#' # d_seg <- msi_segment(
#' #   df = t(as.matrix(d_norm2[["data"]])),
#' #   md = d_norm2[["pixels"]]
#' # )
#'
#' @export
msi_segment <- function(
  df,
  md,
  span = 0.1,
  clus_res = 0.9
) {
  # Load data (ensure that features are rows and pixels are columns)
  d <- df
  mdf <- md
  # Convert matrix into Seurat
  dser <- Seurat::CreateSeuratObject(
    counts = d,
    meta.data = mdf,
    assay = "MSI"
  )
  # Log-transform, scale data, and perform batch corrections using Harmony
  dser <- Seurat::NormalizeData(dser)
  dser <- Seurat::FindVariableFeatures(dser, loess.span = 0.1)
  dser <- Seurat::RunPCA(
    object = Seurat::ScaleData(
      object = dser,
      verbose = TRUE
    ),
    verbose = TRUE
  )
  Seurat::VizDimLoadings(
    object = dser,
    dims = 1:2,
    reduction = "pca"
  )
  Seurat::ElbowPlot(
    dser,
    reduction = "pca",
    ndims = 50
  )
  ## Run Harmony
  dser <- harmony::RunHarmony(
    dser,
    assay.use = "MSI",
    group.by.vars = "ID",
    reduction.use = "pca",
    reduction.save = "MSI.cor",
    project.dim = FALSE
  )
  ## Run UMAP
  dser <- Seurat::RunUMAP(
    dser,
    reduction = "MSI.cor",
    reduction.name = "umap.cor",
    reduction.key = "umapcor_",
    dims = 1:30,
    n.components = 3
  )
  dser <- Seurat::FindNeighbors(dser, dims = 1:30)
  dser <- Seurat::FindClusters(
    dser, resolution = 0.9, n.iter = 25, random.seed = 1234
  )
  dser <- Seurat::AddMetaData(
    dser,
    metadata = as.factor(as.numeric(dser@meta.data$seurat_clusters)),
    col.name = "cluster"
  )
  return(dser)
}

#' MSI UMAP
#'
#' Generates a panel of UMAPs given a
#' Seurat object containing PCA and UMAP results.
#'
#' @param so A Seurat object.
#' @param md_list A vector of character strings indicating
#' metadata columns for overlaying on a loadings plot.
#' @param slot1 A character string corresponding to the umap slot name to plot.
#' @return A series of UMAPs with specified metadata overlays.
#' @examples
#'
#' # p_umap <- msi_umap_panel(d_integrated,c("col1","col2","col3"),"wnn.umap")
#'
#' @export
msi_umap_panel <- function(
  so,
  md_list,
  slot1
) {
  d <- so
  if(ncol(d@reductions[[slot1]]@cell.embeddings) == 3) { #nolint
    d2 <- data.frame(
      d@meta.data,
      `UMAP.1` = d@reductions[[slot1]]@cell.embeddings[, 1],
      `UMAP.2` = d@reductions[[slot1]]@cell.embeddings[, 2],
      `UMAP.3` = d@reductions[[slot1]]@cell.embeddings[, 3]
    )
    d2_list <- setNames(
      lapply(
        c(md_list),
        function(x) {
          d2[, c(
            x,
            "UMAP.1",
            "UMAP.2",
            "UMAP.3"
          )
          ]
        }
      ),
      c(md_list)
    )

  }
  if(ncol(d@reductions[[slot1]]@cell.embeddings) == 2) { #nolint
    d2 <- data.frame(
      d@meta.data,
      `UMAP.1` = d@reductions[[slot1]]@cell.embeddings[, 1],
      `UMAP.2` = d@reductions[[slot1]]@cell.embeddings[, 2]
    )
    d2_list <- setNames(
      lapply(
        c(md_list),
        function(x) {
          d2[, c(
            x,
            "UMAP.1",
            "UMAP.2"
          )
          ]
        }
      ),
      c(md_list)
    )
  }
  # Generate plots
  d2_plot <- lapply(
    c(md_list),
    function(x) {
      p <-  ggplot2::ggplot(
        d2_list[[x]],
        ggplot2::aes(
          x=`UMAP.1`, # nolint
          y=`UMAP.2`, # nolint
          color = .data[[x]], # nolint
          label = .data[[x]] # nolint
        )
      ) +
        ggplot2::scale_color_manual(
          paste(""),
          values = col_univ() # nolint
        ) +
        # Add points
        ggplot2::geom_point(
          shape = 16,
          size = 1,
          alpha = 0.6
        ) +
        msi_theme1() + # nolint
        ggplot2::theme(
          panel.grid.major.y = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          plot.margin = ggplot2::unit(
            c(0.1, 0.1, 0.1, 0.1),
            "cm"
          ),
          legend.position = c(
            0.9,
            0.85
          )
        )
      return(p)
    }
  )
  # Combine output
  if(length(d2_plot) > 1) { # nolint
    d2_out <- ggpubr::ggarrange(
      plotlist = d2_plot,
      labels = c(
        names(
          d2_plot
        )
      ),
      ncol = ifelse(
        length(
          d2_plot
        ) <= 3,
        length(
          d2_plot
        ),
        3
      ),
      nrow = ifelse(
        length(
          d2_plot
        ) > 3,
        ceiling(
          length(
            d2_plot
          ) /
            3
        ),
        1
      )
    )
  }
  if(length(d2_plot) == 1) { # nolint
    d2_out <- d2_plot[[1]]
  }
  return(d2_out)
}

#' Top-10 Marker Gene Heatmap
#'
#' Generates a heatmap from a Seurat Object and
#' marker gene list based on the top-10 marker genes for each cluster.
#'
#' @param so An object of class Seurat.
#' @param asy Assay to use (ex. "RNA").
#' @param cl_var Character string containing the name
#' of the cluster variable for cell type predictions.
#' @param h_w Numeric value for heatmap width (passed to ComplexHeatmap).
#' @param h_h Numeric value for heatmap height (passed to ComplexHeatmap).
#' @param fs_c Numeric value for column fontsize (passed to ComplexHeatmap).
#' @param fs_r Numeric value for row fontsize (passed to ComplexHeatmap).
#' @param cl_c Cluster columns?
#' @param cl_r Cluster rows?
#' @param rot_c Rotation of column names.
#' @param col1 Gradient color scheme to use
#' (must be exactly 4 colors in length).
#' @return A ComplexHeatmap object containing a top-10 marker gene heatmap.
#' @examples
#'
#' # p_umap <- sc_top10_marker_heatmap(
#' #   d_annotated,
#' #   "RNA",
#' #   "seurat_clusters",
#' #   18,
#' #   24,
#' #   6,
#' #   8,
#' #   TRUE,
#' #   TRUE,
#' #   col_grad()[c(3, 6, 9, 12)]
#' # )
#'
#' @export
sc_top10_marker_heatmap <- function(
  so,
  asy,
  cl_var,
  h_w,
  h_h,
  fs_c,
  fs_r,
  cl_c,
  cl_r,
  rot_c,
  col1
) {
  d <- so
  if(!file.exists("analysis/table.marker.genes.txt") && asy == "RNA") { # nolint
    print(
      "No marker gene file has been created;
      calculating marker genes for each cluster..."
    )
    Seurat::DefaultAssay(d) <- "RNA"
    cl_mark <- Seurat::FindAllMarkers(d, verbose = TRUE)
    write.table(
      cl_mark,
      "analysis/table.marker.genes.txt",
      col.names = TRUE,
      row.names = FALSE,
      sep = "\t"
    )
  }

  if(!file.exists("analysis/table.marker.motifs.txt") && asy == "ATAC") { # nolint
    print(
      "No marker motif file has been created;
      calculating marker motifs for each cluster..."
    )
    Seurat::DefaultAssay(d) <- "ATAC"
    cl_mark <- Seurat::FindAllMarkers(
      d,
      min.pct = 0.05,
      verbose = TRUE
    )
    names_motif <- data.frame(
      "gene" = seq.int(1, nrow(d@assays$ATAC@meta.features), 1),
      "near.gene" = paste(
        d@assays$ATAC@meta.features[["nearestGene"]],
        seq.int(1, nrow(d@assays$ATAC@meta.features), 1),
        sep = "."
      ),
      "motif" = paste(
        d@assays$ATAC@meta.features[["seqnames"]],
        paste(
          d@assays$ATAC@meta.features[["start"]],
          d@assays$ATAC@meta.features[["end"]],
          sep = "-"
        ),
        sep = ":"
      )
    )
    cl_mark <- dplyr::left_join(
      cl_mark,
      names_motif,
      by = "gene"
    )
    write.table(
      cl_mark,
      "analysis/table.marker.motifs.txt",
      col.names = TRUE,
      row.names = FALSE,
      sep = "\t"
    )
  }

  if(asy == "RNA") { # nolint
    cl_mark <- read.table(
      "analysis/table.marker.genes.txt",
      sep = "\t",
      header = TRUE
    )
  }

  if(asy == "ATAC") { # nolint
    cl_mark <- read.table(
      "analysis/table.marker.motifs.txt",
      sep = "\t",
      header = TRUE
    )
  }
  ## Marker gene input matrix (top10 per cell type)
  if(class(cl_mark[["cluster"]]) == "character") { # nolint
    cl_mark <- cl_mark[gtools::mixedorder(cl_mark[["cluster"]]), ]
  }
  cl_mark[["CellType.no"]] <- cl_mark[["cluster"]]

  cl_mark <- dplyr::group_by(
    cl_mark,
    .data[["CellType.no"]] # nolint
  )
  ### Top 10 genes per cluster (by p value then by fold change)
  cl_mark <- dplyr::slice_max(
    cl_mark[cl_mark[["avg_log2FC"]] > 0, ],
    order_by = -.data[["p_val_adj"]], # nolint
    n = 25
  )[, c(
    "gene",
    "cluster",
    "avg_log2FC",
    "p_val_adj"
  )]

  cl_mark <- dplyr::group_by(
    cl_mark,
    .data[["cluster"]] # nolint
  )

  cl_mark <- dplyr::slice_max(
    cl_mark,
    order_by = .data[["avg_log2FC"]], # nolint
    n = 10
  )[, c(
    "gene",
    "cluster"
  )]

  #### Save table
  if(asy == "RNA") { # nolint
    write.table(
      cl_mark,
      "analysis/table.marker.genes.top10.txt",
      row.names = FALSE,
      col.names = TRUE,
      sep = "\t"
    )
    SeuratObject::DefaultAssay(d) <- "RNA"
  }
  if(asy == "ATAC") { # nolint
    write.table(
      cl_mark,
      "analysis/table.marker.motifs.top10.txt",
      row.names = FALSE,
      col.names = TRUE,
      sep = "\t"
    )
    SeuratObject::DefaultAssay(d) <- "ATAC"
  }
  ### Subset seurat and scale
  h <- SeuratObject::FetchData(
    d,
    vars = c(
      cl_var,
      unique(cl_mark[["gene"]])
    )
  )
  ### Heatmap annotation (average expression)
  h_anno <- as.data.frame(
    lapply(
      h[, 2:ncol(
        h
      )],
      function(x) {
        mean(x)
      }
    )
  )

  h_anno <- h_anno[, h_anno[1, ] > 0]
  ### Scale and plot average expression/accessibility per cell type
  h_in <- scale(
    as.matrix(
      magrittr::set_rownames(
        setNames(
          as.data.frame(
            lapply(
              h[, 2:ncol(
                h
              )],
              function(x) {
                dplyr::select(
                  aggregate(
                    x,
                    list(
                      h[, 1]
                    ),
                    FUN = mean
                  ),
                  c(2)
                )
              }
            )
          ),
          names(h[, 2:ncol(h)])
        ),
        levels(h[, 1])
      )
    ),
    center = TRUE
  )
  qs <- quantile(
    h_in,
    probs = c(
      0.05,
      0.95
    ),
    na.rm = TRUE
  )

  h_in <- as.matrix(
    as.data.frame(h_in)[, unlist(
      lapply(
        seq.int(1, ncol(as.data.frame(h_in)), 1),
        function(x) {
          !anyNA(as.data.frame(h_in)[x])
        }
      )
    )
    ]
  )

  fun_hm_col <- circlize::colorRamp2(
    c(
      qs[[1]],
      (qs[[1]]) / 2,
      (qs[[2]]) / 2,
      qs[[2]]
    ),
    colors = col1
  )
  # Create Plot
  h_out <- ComplexHeatmap::Heatmap(
    h_in,
    col = fun_hm_col,
    name = "Scaled Expression",
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      `Average.Expression` = ComplexHeatmap::anno_barplot(
        as.vector(t(h_anno)),
        gp = grid::gpar(fill = col1) # nolint
      ),
      annotation_name_gp = grid::gpar(
        fontsize = 10
      )
    ),
    show_column_names = TRUE,
    show_row_names = TRUE,
    heatmap_width = ggplot2::unit(h_w, "cm"),
    heatmap_height = ggplot2::unit(h_h, "cm"),
    column_title = "Top 10 Markers",
    column_names_rot = rot_c,
    column_names_gp = grid::gpar(fontsize = fs_c),
    row_names_side = "left",
    row_names_gp = grid::gpar(fontsize = fs_r),
    cluster_columns = cl_c,
    cluster_rows = cl_r,
  )
  return(h_out)
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






