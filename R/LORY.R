##' @title
##' LORY
##'
##' @description
##' R implementation of LORY routine from Marquez *et al.*, 2012.
##'
##' @param coords array of k landmarks, m dimensions, n individuals
##' @param tesselation matrix of landmark associations (one polygon/polyhedron per line)
##' @param parallel use registered parallel backend? (probably only works in linux now)
##' 
##' @importFrom plyr laply aaply alply
##' @importFrom shapes procGPA
##' 
##' @export
##' @rdname LORY
##' 

LORY <- function (coords, tesselation, parallel)
    {
        gpa <- procGPA(coords)
        mshape <- gpa $ mshape
        print('gpa done')
        dimnames (mshape) <- dimnames (coords) [1:2]
        tps <- alply (gpa $ rotated, 3, ThinPlateSpline,
                      reference.shape = mshape, .parallel = parallel)
        print ('tps done')
        jacobs <- laply (tps, JacobianArray, tesselation = tesselation,
                         .parallel = parallel)
        jacobs <- aperm (jacobs, c(2, 3, 1, 4), resize = parallel)
        print ('jacobs done')
        local <- aaply (jacobs, 4, Center2MeanJacobian, .parallel = parallel)
        local <- t (local)
        return (list ('tps' = tps,
                      'jacobians' = jacobs,
                      'local' = local,
                      'reference' = mshape,
                      'cs' = gpa $ size
                      ))
    }
