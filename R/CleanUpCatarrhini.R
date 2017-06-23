##' @title
##' CleanUpCatarrhini
##' 
##' @description
##' This function cleans up Catarrhini data from excel files.
##' 
##' @param xls.list output obtained from ReadCatarrhini
##'
##' @return a list with two elements, info and coord
##'
##' @details will GlueSkulls
##' 
##' @author Guilherme Garcia
##'
##' @importFrom plyr llply aaply ldply
##'
##' @seealso ReadHomo

CleanUpCatarrhini <- function(xls.list)
{
    A <- llply(xls.list,
               function(L) suppressWarnings(aaply(L $ A, 1:4, as.numeric)))

    Z <- llply(xls.list,
               function(L) suppressWarnings(aaply(L $ Z, 1:4, as.numeric)))

        emptyA <-
        suppressWarnings(unlist(llply(A,
                                      function(L) aaply(L, 4,
                                                        function(m) all(is.na(m))))))

    emptyZ <-
        suppressWarnings(unlist(llply(Z,
                                      function(L) aaply(L, 4,
                                                        function(m) all(is.na(m))))))
    empty.any <- emptyA | emptyZ

    id <- ldply(xls.list, function(L) L $ id)
    
    names(empty.any) <- id $ numero

    Ns <- laply(A, function(L) dim(L) [4])
   
    coordA <- array(0, c(dim(A [[1]])[1:3], sum(Ns)))
    coordZ <- array(0, c(dim(Z [[1]])[1:3], sum(Ns)))

    dimnames(coordA)[1:3] <- dimnames(A [[1]]) [1:3]
    dimnames(coordZ)[1:3] <- dimnames(Z [[1]]) [1:3]
    dimnames(coordA)[[4]] <- dimnames(coordZ)[[4]] <- id $ numero

    counter <- 0
    for(i in 1:length(Ns))
    {
        coordA [, , , counter+(1:Ns[i])] <- A [[i]]
        coordZ [, , , counter+(1:Ns[i])] <- Z [[i]]

        counter <- counter + Ns[i]
    }

    ## tirar vazios
    id <- id [!empty.any, ]
    coordA <- coordA [, , , !empty.any]
    coordZ <- coordZ [, , , !empty.any]
