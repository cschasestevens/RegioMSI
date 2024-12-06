% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4.1a.stats.functions.R
\name{MSI.stat.t}
\alias{MSI.stat.t}
\title{Student's T-test of MSI Data}
\usage{
MSI.stat.t(df, md, var.g, g1, g2, l.names)
}
\arguments{
\item{df}{Normalized data matrix containing cluster, sample ID, pixel coordinate, and pixel intensity information.}

\item{md}{Number of metadata columns present in input data.}

\item{var.g}{Treatment group column name.}

\item{g1}{Treatment group 1 name.}

\item{g2}{Treatment group 2 name.}

\item{l.names}{Vector of compound names equal to the number of compounds present in the given data matrix.}
}
\value{
A data frame containing T-test results for each compound.
}
\description{
Conducts a Student's T-test for each compound.
}
\examples{
## calculate sample means and conduct hypothesis test
d.p2 <- dplyr::left_join(
  data.frame(
    unique(
      d2[,c("ID","Group")]
    )
  ),
  setNames(
    purrr::reduce(
      lapply(
        seq.int(
          (8 + 1),
          ncol(d2),
          1
        ),
        function(x)
          aggregate(
            d2[,
               names(d2[x])
            ],
            list(
              d2[["ID"]]
            ),
            function(y)
              (mean(y)))),
      dplyr::left_join,
      by = "Group.1"
    ),
    c(
      "ID",
      names(d2[,(8 + 1):ncol(d2)])
    )
  ),
  by = "ID"
)

d.p2 <- data.frame(
  d.p2[,1:2],
  as.data.frame(
    lapply(
      d.p2[,3:ncol(d.p2)],
      function(x)
        log2(x)
    )
  )
)

d.p2.result <- cbind(
  data.frame(
    "Comparison" = c("M.HO3-M.SFA","F.HO3-F.SFA")
  ),
  rbind(
    setNames(
      as.data.frame(
        lapply(
          names(
            d.p2[,3:ncol(d.p2)]
          ),
          function(x)
            MSI.stat.t(d.p2,2,"Group","M.HO3","M.SFA",x))),
      c(names(d[,9:ncol(d)]))),
    setNames(
      as.data.frame(
        lapply(
          names(
            d.p2[,3:ncol(d.p2)]
          ),
          function(x)
            MSI.stat.t(d.p2,2,"Group","F.HO3","F.SFA",x))),
      c(names(d[,9:ncol(d)])))
  )
)

}
