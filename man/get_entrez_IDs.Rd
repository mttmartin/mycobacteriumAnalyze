% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mycobacteriumAnalyze.R
\name{get_entrez_IDs}
\alias{get_entrez_IDs}
\title{Retrieves entrez IDs for a given species that are compatible for
GO enrichment analysis.}
\usage{
get_entrez_IDs(data, species)
}
\arguments{
\item{data}{The orignal protein data composed of a 'gene' column and a 'protein' column}

\item{species}{The species you are analyzing (either 'avium' or 'abscessus')}
}
\value{
A vector of entrez IDs
}
\description{
Retrieves entrez IDs for a given species that are compatible for
GO enrichment analysis.
}
