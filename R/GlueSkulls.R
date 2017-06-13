##' @title
##' GlueSkulls
##'
##' @description
##' This function joins two different sets of landmark configuration which correspond
##' to the same individuals
##'
##' @param A first set of landmark configurations
##' @param Z second set of landmark configurations
##'
##' @return 
##' array with joined configurations
##'
##' @author Guilherme Garcia
##'
##' @importFrom plyr aaply
##' @importFrom shapes procOPA
##' 
##' @examples
##' 
GlueSkulls <- function (A, Z)
  {
    glueSkull <- function (Ai, Zi)
      {
        mA <- dim (Ai) [1]
        mZ <- dim (Zi) [1]
        lmA <- dimnames (Ai) [[1]]
        lmZ <- dimnames (Zi) [[1]]
        cA <- which (lmA %in% lmZ)
        cZ <- which (lmZ %in% lmA)
        ### ordenar landmarks
        cA <- cA[order (lmA[cA])]
        cZ <- cZ[order (lmZ[cZ])]
        ### centroides dos pontos em comum
        ccA <- t (array (colMeans (Ai[cA,]), c(3, mA)))
        ccZ <- t (array (colMeans (Zi[cZ,]), c(3, mZ)))
        ### centralizando no centroide dos pontos em comum
        Ai <- Ai - ccA
        Zi <- Zi - ccZ
        ### the easy way out
        R <- procOPA(Ai[cA, ], Zi[cZ, ], scale = FALSE, reflect = FALSE) $ R
        Zi <- Zi %*% R
        rbind(Ai, Zi[!(lmZ %in% lmA), ])
      }
    if (dim (A) [3] != dim (Z) [3])
    {
        cat ('Número de vistas A não bate com número de vistas Z. Verifique!','\n')
        return (-1)
    }
    else
    {
        out <- aaply(1:dim (A)[3], 1, function(i) glueSkull (A[, , i], Z[, , i]))
        right.first = function (element)
          {
            return (ifelse (length (element) == 1, 0,
                            ifelse (element [1] == 'NLT', 3,
                                    ifelse (element [2] == 'D', 1, 2))))
          }
        ord = strsplit (dimnames (out) [[1]], split = '-')
        ord = sapply (ord, right.first)
        out = out [order (ord),,]
        return (list (out, dets))
      }
  }
