##' @title
##' ThinPlateSpline
##'
##' @description
##' Estimates thin plate spline interpolation between reference and target shape
##'
##' @param reference.shape reference shape
##' @param target.shape target shape
##'
##' @return a list with some elements
##'
##' @export
##' @rdname ThinPlateSpline
##'
ThinPlateSpline <- function (target.shape, reference.shape)
  {
      ## calculates thin plate splines deformations between reference and target shapes
      theta.p <- function (p)
          return (b + (A %*% p) + t (W) %*% interpol.L1.TPS (p, Q))
      Q <- reference.shape
      P <- target.shape
      k <- nrow (Q)
      m <- ncol (Q)
      R <- interpol.L2.TPS (Q, Q)
      one.k <- array (1, c (k, 1))
      zero.11 <- array (0, c (1, 1))
      zero.m1 <- array (0, c (m, 1))
      zero.mm <- array (0, c (m, m))
      L <- rbind (R, t (one.k), t (Q))
      L <- cbind (L, rbind (one.k, zero.11, zero.m1))
      L <- cbind (L, rbind (Q, t (zero.m1), zero.mm))
      P.zero <- rbind (P, t (zero.m1), zero.mm)
      spline <- solve (L, P.zero)
      W <- spline [1:k, ]
      b <- t (t (spline [1+k, ]))
      A <- t (spline [(2:(m+1))+k, ])
      return (list ('W' = W, 'A' = A,
                    'b' = b, 'Q' = Q, 'P' = P, 'theta.p' = theta.p))
  }
