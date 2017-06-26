##' @title
##' GlueSkull
##'
##' @description
##' This function joins two hemiskulls
##'
##' @param Ai first set of landmark configurations
##' @param Zi second set of landmark configurations
##'
##' @return 
##' matrix with joined configurations
##'
##' @author Guilherme Garcia
##'
##' @importFrom shapes procOPA
##' 
##' @examples
##' 
GlueSkull <- function (Ai, Zi)
{
    mA <- dim (Ai) [1]
    mZ <- dim (Zi) [1]
    lmA <- dimnames (Ai) [[1]]
    lmZ <- dimnames (Zi) [[1]]
    cA <- which (lmA %in% lmZ)
    cZ <- which (lmZ %in% lmA)
    
    ## ordenar landmarks
    cA <- cA[order (lmA[cA])]
    cZ <- cZ[order (lmZ[cZ])]
    
    ## NA
    missing.entries <- is.na(Ai [cA, 1]) | is.na(Zi [cZ, 1])
    if(any(missing.entries)) {
        cA <- cA [!missing.entries]
        cZ <- cZ [!missing.entries] }
    
    ## centroides dos pontos em comum
    ccA <- t (array (colMeans (Ai[cA,]), c(3, mA)))
    ccZ <- t (array (colMeans (Zi[cZ,]), c(3, mZ)))
    
    ## centralizando no centroide dos pontos em comum
    Ai <- Ai - ccA
    Zi <- Zi - ccZ
    
    ## the easy way out
    R <- procOPA(Ai[cA, ], Zi[cZ, ], scale = FALSE, reflect = FALSE) $ R
    Zi <- Zi %*% R
    rbind(Ai, Zi[!(lmZ %in% lmA), ])    
}
