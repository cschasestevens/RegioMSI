#' Seurat-based KNN Clustering and Segmentation
#'
#' Converts MSI data to a Seurat object and performs spatial
#' segmentation using the default implementation of KNN
#' used by Seurat.
#'
#' @param df A normalized data matrix.
#' @param md Number of metadata columns present in data matrix.
#' @param var_id Name of sample ID variable. If only one sample
#' is present, then batch effect correction based on Harmony is
#' excluded from the resulting image segmentation.
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
  var_id,
  span = 0.1,
  clus_res = 0.9
) {
  # Load data (ensure that features are rows and pixels are columns)
  d <- df
  mdf <- md
  if( # nolint
    !file.exists(
      paste(
        "analysis/data.", lp1[["polarity"]], ".segmented.seurat.rds", sep = "" # nolint
      )
    )
  ) {
    # Convert matrix into Seurat
    dser <- Seurat::CreateSeuratObject(
      counts = d,
      meta.data = mdf,
      assay = "MSI"
    )
    head(dser@meta.data)
    # Log-transform, scale data, and perform batch corrections using Harmony
    dser <- Seurat::NormalizeData(dser)
    dser <- Seurat::FindVariableFeatures(dser, loess.span = span)
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
    if(length(unique(as.character(mdf[[var_id]]))) == 1) { # nolint
      print("Only one sample is present; skipping batch correction step...")
      dser <- Seurat::RunUMAP(
        dser,
        reduction = "pca",
        reduction.name = "umap",
        reduction.key = "umap_",
        dims = 1:30,
        n.components = 3
      )
      dser <- Seurat::FindNeighbors(dser, dims = 1:30)
    }
    if(length(unique(as.character(mdf[[var_id]]))) > 1) { # nolint
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
    }
  }
  if( # nolint
    file.exists(
      paste(
        "analysis/data.", lp1[["polarity"]], ".segmented.seurat.rds", sep = "" # nolint
      )
    )
  ) {
    print("A Seurat object already exists for this dataset; loading existing object...") # nolint
    dser <- readRDS(
      paste(
        "analysis/data.", lp1[["polarity"]], ".segmented.seurat.rds", sep = "" # nolint
      )
    )
  }
  dser <- Seurat::FindClusters(
    dser, resolution = clus_res, n.iter = 25, random.seed = 1234
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

#' Top-10 Feature Heatmap
#'
#' Generates a heatmap from a Seurat Object containing
#' clustered MSI data. By default, the top-10 features per
#' cluster are returned.
#'
#' @param so A Seurat object.
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
#' # p_hmap_top <- msi_marker_heatmap(
#' #   so = d_seg,
#' #   h_w = 18,
#' #   h_h = 24,
#' #   fs_c = 6,
#' #   fs_r = 8,
#' #   cl_c = TRUE,
#' #   cl_r = TRUE,
#' #   rot_c = 45,
#' #   col1 = col_grad()[c(3, 6, 9, 12)]
#' # )
#'
#' @export
msi_marker_heatmap <- function(
  so,
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
  if(!file.exists("analysis/table.marker.features.txt")) { # nolint
    print(
      "No marker feature file has been created;
      calculating marker features for each cluster..."
    )
    Seurat::DefaultAssay(d) <- "MSI"
    cl_mark <- Seurat::FindAllMarkers(
      d,
      test.use = "wilcox",
      densify = TRUE,
      verbose = TRUE
    )
    write.table(
      cl_mark,
      "analysis/table.marker.features.txt",
      col.names = TRUE,
      row.names = FALSE,
      sep = "\t"
    )
  }
  if(file.exists("analysis/table.marker.features.txt")) { # nolint
    cl_mark <- read.table(
      "analysis/table.marker.features.txt",
      sep = "\t",
      header = TRUE
    )
  }
  ## Marker feature input matrix (top10 per region)
  if(class(cl_mark[["cluster"]]) == "character") { # nolint
    cl_mark <- cl_mark[gtools::mixedorder(cl_mark[["cluster"]]), ]
  }
  cl_mark[["region"]] <- cl_mark[["cluster"]]
  cl_mark[["feature"]] <- cl_mark[["gene"]]
  cl_mark <- dplyr::group_by(
    cl_mark,
    .data[["region"]] # nolint
  )
  ### Top 10 features per cluster (by p value then by fold change)
  cl_mark <- dplyr::slice_max(
    cl_mark[cl_mark[["avg_log2FC"]] > 0, ],
    order_by = -.data[["p_val_adj"]], # nolint
    n = 25
  )[, c(
    "feature",
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
    "feature",
    "cluster"
  )]
  #### Save table
  write.table(
    cl_mark,
    "analysis/table.marker.features.top10.txt",
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t"
  )
  ### Subset seurat and scale
  h <- SeuratObject::FetchData(
    d,
    vars = c(
      "cluster",
      unique(cl_mark[["feature"]])
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
  ### Scale and plot average intensity per region
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
    name = "Scaled Intensity",
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      `Average.Intensity` = ComplexHeatmap::anno_barplot(
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
