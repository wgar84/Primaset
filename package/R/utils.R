##' @title
##' Auxiliary functions
##'
##' @description
##' Some internal functions for Primaset
##'
##' @importFrom numDeriv jacobian
##' @importFrom stats var na.omit
##' @importFrom expm expm logm

Norm <- function(x) sqrt(sum(x * x))
Normalize <- function(x) x / Norm(x)

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
    Xm <- na.omit (X [midline, ])
    Mm <- matrix (apply (Xm, 2, mean), byrow = TRUE, nrow = nrow(X), ncol = ncol(X))
    Xc <- X - Mm
    W <- na.omit (Xc [midline, ])
    RM <-svd (var (W))$v
    return (RM)
  }

JacobianArray <- function (spline, tesselation, ...)
    {
        ## calculates jacobians for a given interpolation in a set of points
        ## determined from tesselation (as centroids of each tetrahedron defined, for now...)
        with (spline,
              {
                  Q.tetra <- Q [tesselation, ]
                  dim (Q.tetra) <- c (dim (tesselation), ncol (Q))
                  Q.centroids <- apply (Q.tetra, 1, colMeans)
                  jacobs <- apply (Q.centroids, 2, jacobian, func = theta.p, ...)
                  dim (jacobs) <- c (ncol (Q), ncol (Q), ncol (Q.centroids))
                  return (jacobs)
              })
    }

Center2MeanJacobian <- function (jacobArray, max.steps = 100)
    {
        ### calculates mean jacobian matrix for a set of jacobian matrices
        ### describing a local aspect of shape deformation for a given set of volumes,
        ### returning log determinants of deviations from mean jacobian (Woods, 2003).
        logm.single <- function (Ai, inv.Mk) return (logm (Ai %*% inv.Mk))
        A <- jacobArray
        N <- dim (A) [3]
        Mk <- diag (nrow (A))
        i <- 1
        repeat
        {
            inv.Mk <- solve (Mk)
            centered.now <- apply (A, 3, logm.single, inv.Mk = inv.Mk)
            o <- array (- rowMeans (centered.now), c(nrow (A), nrow (A)))
            if (all (abs (o) < .Machine$double.eps))
                break
            Mk <- expm (- o) %*% Mk
            i <- i + 1
            if (i == max.steps)
                stop ('Convergence has not been achieved in number of steps.')
        }
        centered.now <- array (centered.now, c(nrow (A), nrow (A), N))
        log.det <- apply (centered.now, 3, function (x) return (sum (diag (x))))
        return (log.det)
    }
