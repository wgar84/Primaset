##' @title
##' TwoBlockCommonFactor
##'
##' @description
##' Computes common factors for two sets of traits, as described in Mitteroecker and
##' Bookstein (2008).
##'
##' @param X n x p multivariate data set to project onto estimtated factors,
##'          with n being number of individuals and p the number of traits
##' @param cov.matrix covariance matrix to produce factors
##' @param modules vector distinguishing two modules, in binary form
##' @param permutations how many permutations we should use to evaluate how many factors
##'                     the analysis should retain; defaults to NULL, in which this
##'                     procedure is not done
##'
##' @importFrom plyr aaply
##' 
##' @references
##' Mitteroecker, P, and F L Bookstein. "The Evolutionary Role of Modularity and
##' Integration in the Hominoid Cranium." Evolution 62, no. 4 (2008): 943–958.
##' https://doi.org/10.1111/j.1558-5646.2008.00321.x.
##'
TwoBlockCommonFactor <- function(X, cov.matrix = NULL, modules, permutations = NULL)
{
    ## if cov.matrix is not provided, uses cov.matrix from data
    if(is.null(cov.matrix))
        cov.matrix <- var(X)

    ## gonna call subsets of F and N, just because
    
    f <- modules == 1

    Xf <- X[, f]
    Xn <- X[, !f]
    
    Sfn <- cov.matrix[f, !f]

    UDV <- svd(Sfn)

    U <- UDV $ u
    V <- UDV $ v
    D <- UDV $ d
    
    XfU <- Xf %*% U

    XnV <- Xn %*% V

    ## retention of k factors
    if(is.null(permutations))
        retain <- length(D)
    else
    {
        rand.sigma <- FALSE

        i <- 2

        while(!rand.sigma)
        {
            Zfi <- Xf - (XfU[, 1:(i - 1)] %*% t(U [, 1:(i - 1)]))
            
            Zni <- Xn - (XnV[, 1:(i - 1)] %*% t(V [, 1:(i - 1)]))

            di.dist <- c()
        
            for(j in 1:permutations)
            {
                Zni.perm <- Zni [sample(1:nrow(Zni)), ]
                covRes.perm <- var(Zfi, Zni.perm)
                di.dist[j] <- svd(covRes.perm) $ d [1]
            }

            p.perm <- sum(di.dist > D[i]) / permutations

            if(p.perm < 0.01) # good cutoff value is... ?
                i <- i + 1
            else
            {
                rand.sigma <- TRUE
                retain <- i - 1
                print(paste('Retained', retain, 'common factors'))
            }
        }
    }

    ## scaling

    F <- aaply(1:retain, 1, function(i)
        {
            Ci <- cbind(XfU[, i], XnV[, i])

            e1 <- eigen(t(Ci) %*% Ci) $ vectors [, 1]

            fi <- c(e1[1] * U[, i], e1[2] * V[, i])

            fi
        })

    F <- t(F) # n x k

    XF <- X %*% F

    R <- X - XF %*% t(F)

    list('integrated' = XF,
         'factors' = F, 
         'modular' = R)
    
}


