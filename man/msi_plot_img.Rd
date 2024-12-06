% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2.3b.normalize.sLOESS.R
\name{msi_plot_img}
\alias{msi_plot_img}
\title{Ion Image QC}
\usage{
msi_plot_img(df, var_y, var_g, perc_int, spl_samp = FALSE, samp_no = NULL)
}
\arguments{
\item{df}{A normalized data matrix.}

\item{var_y}{Variable corresponding to the pixel intensity for each sample.}

\item{var_g}{Sample ID variable.}

\item{perc_int}{Quantile to be used for the intensity color scale
to improve contrast in ion images.}

\item{spl_samp}{Plot a specific sample
(default is to plot all sample images).}

\item{samp_no}{If spl_samp is TRUE, which sample number should be plotted?}
}
\value{
An ion image visualizing total normalized pixel intensities.
}
\description{
Plots the normalized pixel intensities as an ion image for each sample.
}
\examples{

# ptic_img <- msi_plot_img(
#   df = d_norm[["pixels"]],
#   var_y = "TIC.norm.tic",
#   var_g = "ID",
#   perc_int = 0.95
# )

}
