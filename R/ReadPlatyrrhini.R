##' @title
##' ReadPlatyrrhini
##'
##' @description
##' This functions read Platyrrhini POLHEMUS data, stored in ascii files
##'
##' @param path path to Platyrrhini raw data
##'
##' @return
##' list with three elements
##' id: information
##' A: coordinate data for frontal view (most landmarks)
##' Z: coordinate data for back view (some neuro landmarks)
##'
##' @details
##' Although configurations are ordered, some individuals have been
##' measured more than once. Thus, elements in both A and Z exceed entries in id.
##' Individual configurations are labeled in both A and Z so that
##' replicates are identified, in the form (number in id).(replicate).
##' Also notice that the number of replicates for each A and Z may be different.
##' 
##' @author Guilherme Garcia
##'
##' @importFrom gdata read.xls
##' @importFrom plyr aaply alply
##' @importFrom shapes centroid.size
##' 
##' @examples
##'
##' ReadPlatyrrhini(path = '../Raw Data/Platyrrhini')

ReadPlatyrrhini <- function(path)
{
    readSingleFile <- function(filename)
    {
        shape <- scan(filename, what = '')
        shape <- shape [-1]
        shape <- matrix(as.numeric(shape), ncol = 3, byrow = TRUE)
        ## check for missing lms
        center.shape <- scale(shape, scale = FALSE)
        center.shape <- center.shape / centroid.size(center.shape)
        dist.to.centroid <- diag(center.shape %*% t(center.shape))
        miss.lm <- dist.to.centroid > 0.2 ## 10 times av. dist on A, 2.5 times av. on Z
        shape[miss.lm, ] <- NA
        if(nrow(shape) == 35)
            shape <- shape [1:33, ]
        shape
    }
    
    id.raw <- read.csv(paste0(path, '/semsp12.csv'))

    id.clean <- read.csv(paste0(path, '/GLMALL2_nodist.csv'))

    ## landmarks
    lms <- read.csv(paste0(path, 'vistas.csv'), header=FALSE)
    lms <- as.character(lms[, 1])
    lms <- gsub(' ', '', lms)
    lms.A <- lms[1:33]
    lms.Z <- lms[36:length(lms)]
    
    ## repeated individuals
    id.raw <- id.raw[- which(duplicated(id.raw $ ID.)), ]
    ## só 4

    id.clean <- id.clean[- which(duplicated(id.clean $ ID)), ]
    ## só 2 (que tb se repetiam no raw)

    id.raw $ GENUS. <- gsub(' ', '', id.raw $ GENUS.)
    
    num.raw <- as.character(id.raw $ ID.)
    num.clean <- as.character(id.clean $ ID)

    id.raw <- id.raw[num.raw %in% num.clean, ]
    num.raw <- num.raw[num.raw %in% num.clean]

    id.clean <- id.clean[num.clean %in% num.raw, ]
    num.clean <- num.clean[num.clean %in% num.raw]

    id.raw <- id.raw [order(num.raw), ]
    num.raw <- num.raw[order(num.raw)]

    id.clean <- id.clean [order(num.clean), ]
    num.clean <- num.clean[order(num.clean)]
    
    A <- dir(path, pattern = '\\.A.P', recursive = TRUE, ignore.case = TRUE)
    Z <- dir(path, pattern = '\\.Z.P', recursive = TRUE, ignore.case = TRUE)

    ## para cada elemento, achar os arquivos relevantes
    ## alguns não tinham a letra do número de tombo
    A.list <-
        alply(num.raw, 1,
              function(N)
              {
                  N.A <- paste0('/', N, '.A')
                  files <- A[grep(pattern = N.A, x = A)]
                  if(length(files) == 0)
                      {
                          N <- gsub("[^0-9]", "", N)
                          N.A <- paste0('/', N, '.A')
                          files <- A[grep(pattern = N.A, x = A)]
                      }
                  files
              })

    Z.list <-
        alply(num.raw, 1,
              function(N)
              {
                  N.Z <- paste0('/', N, '.Z')
                  files <- Z[grep(pattern = N.Z, x = Z)]
                  if(length(files) == 0)
                  {
                      N <- gsub("[^0-9]", "", N)
                      N.Z <- paste0('/', N, '.Z')
                      files <- Z[grep(pattern = N.Z, x = Z)]
                  }
                  files
              })
        
    names(A.list) <- names(Z.list) <- num.raw

    no.files <- which(laply(A.list, function(L) length(L) == 0) |
                      laply(Z.list, function(L) length(L) == 0))

    A.list <- A.list[-no.files]
    Z.list <- Z.list[-no.files]

    num.raw <- num.raw[-no.files]
    id.raw <- id.raw [-no.files, ]

    num.clean <- num.clean[-no.files]
    id.clean <- id.clean[-no.files, ]

    id.clean $ DIET <- id.raw $ DIET.
    
    id.clean $ REP.A <- laply(A.list, length)
    id.clean $ REP.Z <- laply(Z.list, length)
    
    A.num <- unlist(llply(1:length(A.list),
                          function(i)
                              paste(num.raw [i], 1:length(A.list [[i]]), sep = '.')))

    Z.num <- unlist(llply(1:length(Z.list),
                          function(i)
                              paste(num.raw [i], 1:length(Z.list [[i]]), sep = '.')))


    
    A.list <- unlist(A.list)

    A.misslm <- !(A.list == 'BBB/A2P/194355.A1P')
    
    A.list <- A.list[A.misslm]
    A.num <- A.num[A.misslm]
    A.sh <- aaply(A.list, 1, function(f) readSingleFile(paste0(path, f)))
   
    dimnames(A.sh) <- list(A.num, lms.A, c('X', 'Y', 'Z'))

    Z.list <- unlist(Z.list)
    Z.sh <- aaply(Z.list, 1, function(f) readSingleFile(paste0(path, f)))

    dimnames(Z.sh) <- list(Z.num, lms.Z, c('X', 'Y', 'Z'))
    
    A.sh <- aperm(A.sh, c(2, 3, 1))
    Z.sh <- aperm(Z.sh, c(2, 3, 1))
    
    list('id' = id.clean, 'A' = A.sh, 'Z' = Z.sh)
    
}
