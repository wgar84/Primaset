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
##' @importFrom plyr aaply alply
##'
##' @examples
##'
##' ReadPlatyrrhini(path = '../Raw Data/Platyrrhini')

ReadPlatyrrhini <- function(path)
{
    id.raw <- read.csv(paste0(path, '/semsp12.csv'))

    id.clean <- read.csv(paste0(path, '/GLMALL2_nodist.csv'))
    
    ## repeated individuals
    id.raw <- id.raw[- which(duplicated(id.raw $ ID.)), ]
    ## só 4

    id.clean <- id.clean[- which(duplicated(id.clean $ ID)), ]
    ## só 2 (que tb se repetiam no raw)

    num.raw <- as.character(id.raw $ ID.)
    num.clean <- as.character(id.clean $ ID)

    id.raw <- id.raw [num.raw %in% num.clean, ]
    num.raw <- num.raw [num.raw %in% num.clean]
    
    A <- dir(path, pattern = '\\.A', recursive = TRUE, ignore.case = TRUE)
    Z <- dir(path, pattern = '\\.Z', recursive = TRUE, ignore.case = TRUE)

    ### para cada elemento, achar os arquivos relevantes
    
    file.list <-
        alply(num.raw, 1,
              function(N)
              {
                  N.A <- paste0(N, '.A')
                  N.Z <- paste0(N, '.Z')
                  files <- c(dir(path, pattern = N.A, recursive = TRUE),
                             dir(path, pattern = N.Z, recursive = TRUE))
                      if(length(files) == 0)
                      {
                          N <- gsub("[^0-9]", "", N)
                          N.A <- paste0(N, '.A')
                          N.Z <- paste0(N, '.Z')
                          files <- c(dir(path, pattern = N.A, recursive = TRUE),
                                     dir(path, pattern = N.Z, recursive = TRUE))
                      }
                      files
              })
    
    table(laply(file.list, length))
    
    file.list[laply(file.list, length) == 519]
    
}
