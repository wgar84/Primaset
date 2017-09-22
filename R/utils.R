##' @title
##' Auxiliary functions
##'
##' @description
##' Some internal functions for Primaset
##'
interpol.TPS <- function (x, y)
      {
        distance <- Norm (x - y)
        if (length (x) == 2)
          {
            if (distance < .Machine $ double.eps)
              return (0)
            else
              return (- log (distance ^ 2) * (distance ^ 2))
          }
        if (length (x) == 3)
          return (distance)
      }

interpol.L1.TPS <- function (y, X)
  return (aaply (X, 1, interpol.TPS, y = y))

interpol.L2.TPS <- function (X, Y)
  return (aaply (Y, 1, interpol.L1.TPS, X = X))

Rotate2MidlineMatrix <- function (X, midline)
  {
    ## returns the rotation matrix that aligns a specimen saggital line
    ## to plane y = 0 (2D) or z = 0 (3D)
    ncl <- ncol (X)
    Xm <- na.omit (X [midline, ])
    Mm <- matrix (apply (Xm, 2, mean), byrow = TRUE, nr = nrow (X), nc = ncl)
    Xc <- X - Mm
    W <- na.omit (Xc [midline, ])
    RM <-svd (var (W))$v
    return (RM)
  }
