##' @title
##' GroupPrimates
##'
##' @description
##' This functions unifies Prosimian and Anthropoid databases in a single object.
##'
##' @param prosimians Prosimian DB, output of CleanUpProsimian
##' @param anthropoids Anthropoid DB, output of GroupAnthropoids
##'
##' @return some list
##' 
##' @author Guilherme Garcia
##'
##' @importFrom shapes procGPA
##' @importFrom geomorph estimate.missing
##' @importFrom plyr aaply alply
##' 
##' @seealso CleanUpCatarrhini CleanUpPlatyrrhini CleanUpHomo
GroupPrimates <- function(prosimians, anthropoids)
{
    ## find me the halflings
    halflings <-
        aaply(prosimians $ coord [, 1, 1, ], 2, function(c) sum(is.na(c)) > 21)

    ## mirror then symmetrize?
    
    
}
