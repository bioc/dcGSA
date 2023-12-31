% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readGMT.R
\name{readGMT}
\alias{readGMT}
\title{Read gene set file in gmt format}
\usage{
readGMT(file = NULL)
}
\arguments{
\item{file}{filename for the gmt file}
}
\value{
a list of gene sets with each element being a vector of gene names
}
\description{
Read gene set file in gmt format
}
\examples{
fpath <- system.file("extdata", "sample.gmt.txt", package="dcGSA")
GS <- readGMT(file=fpath)
}

