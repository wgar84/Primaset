##' @title
##' ReadPlatyrrhini
##'
##' @description
##' This functions read Platyrrhini POLHEMUS data, stored in ascii files
##'
##' @param path path to Platyrrhini raw data
##'
##' @return
##' no idea whatsoever right now
##' 
##' @author Guilherme Garcia
##'
##' @importFrom gdata read.xls
##' @importFrom plyr aaply
##'
##' @examples
##'
##' ReadPlatyrrhini(path = '../Raw Data/Platyrrhini')

ReadPlatyrrhini <- function(path)
{
    id <- read.csv(paste0(path, '/semsp12.csv'))

    ## repeated individuals
    id <- id[- which(duplicated(id [, 1])), ]
    ## só 4

    num <- platyrrhini.id [, 1]


    
}
