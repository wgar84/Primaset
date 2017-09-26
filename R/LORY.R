##' @title
##' LORY
##'
##' @description
##' R implementation of LORY routine from Marquez *et al.*, 2012.
##'
##' @importFrom plyr laply aaply
##' 
##' 
##' @export
##' @rdname LORY
##' 

LORY <- function (coords, tesselation)
    {
        gpa <- procGPA(coords)
        
        dimnames (mshape) <- dimnames (sym.mean)
        tps <- alply (gpa, 3, ThinPlateSpline,
                      reference.shape = mshape, .parallel = TRUE)
        print ('tps done')
        jacobs <- laply (tps, JacobianArray, tesselation = tesselation,
                         .parallel = TRUE)
        jacobs <- aperm (jacobs, c(2, 3, 1, 4), resize = TRUE)
        print ('jacobs done')
        local <- aaply (jacobs, 4, Center2MeanJacobian, .parallel = TRUE)
        local <- t (local)
        return (list ('tps' = tps,
                      'jacobians' = jacobs,
                      'local' = local,
                      'reference' = mshape,
                      'cs' = cs,
                      'info' = info))
    }
