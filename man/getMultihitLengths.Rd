% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/intSiteRetriever.R
\name{getMultihitLengths}
\alias{getMultihitLengths}
\title{lengths distributions for multihits}
\usage{
getMultihitLengths(sample_ref, conn)
}
\arguments{
\item{sample_ref}{dataframe with 2 cols: sampleName, refGenome}

\item{conn}{connection to database}
}
\value{
df with cols: sampleName, refGenome, multihitID, length
}
\description{
lengths distributions for multihits
}

