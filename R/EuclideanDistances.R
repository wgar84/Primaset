##' @title
##' EuclideanDistances
##'
##' @description
##' For a given landmark configuration, calculate euclidean distances between landmarks
##'
##' @param shape landmark configuration with named landmarks
##' @param dists vector of distances one wishes to calculate (as numeric indices for shape)
##'
##' @return vector containing calculated euclidean distances.
##'
##' @importFrom plyr aaply
##' 
##' @author Guilherme Garcia
##'
EuclideanDistances <- function(shape, dists)
{
    all.dists <- dist(shape, 'euclidean', upper = TRUE, diag = TRUE)
    all.dists <- as.matrix(all.dists)

    aaply(dists, 1, function(l) all.dists[l[1], l[2]])
}
