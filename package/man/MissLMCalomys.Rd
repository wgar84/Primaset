% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MissLMCalomys.R
\name{MissLMCalomys}
\alias{MissLMCalomys}
\title{MissLMCalomys}
\usage{
MissLMCalomys(cleanup)
}
\arguments{
\item{cleanup}{output from CleanUpCalomys}
}
\value{
List with two elements:
filled landmark configuration (coord) and specimen information (info).
}
\description{
Fill up missing landmarks in Calomys database, grouping
individuals by age class to perform interpolations based on TPS.
}
\seealso{
CleanUpCalomys
}
\author{
Guilherme Garcia
}
