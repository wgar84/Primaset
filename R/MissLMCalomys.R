##' @title
##' MissLMCalomys
##'
##' @description
##' Fill up missing landmarks in Calomys database.
##'
##' @param cleanup output from CleanUpCalomys
##'
##' @return List with two elements:
##' filled landmark configuration (coord) and specimen information (info).
##'
##' @author Guilherme Garcia
##'
##' @importFrom geomorph estimate.missing
##' @importFrom plyr alply
##'
##' @export
##' @rdname MissLMCalomys
##'
##' @seealso CleanUpCalomys
##'

MissLMCalomys <- function(cleanup)
    {
        ## do it on replicates separately
        ## do it by age class

        ## create age categories
        info <- cleanup $ info
        coord <- cleanup $ coord

        ## duplicated names
        info $ SIB <-
            ifelse(duplicated(info $ ID),
                   paste0(as.character(info $ SIB), 'S'), as.character(info $ SIB))

        info $ ID <-
            as.factor(ifelse(duplicated(info $ ID),
                             paste0(as.character(info $ ID), 'S'),
                             as.character(info $ ID)))

        dimnames(coord) [[4]] <- info $ ID

        out <- array(0, dim(coord))
        dimnames(out) <- dimnames(coord)

        ## predefined based on colony (mis)design
        info $ AGECAT <-
            cut(info $ AGE, c(-Inf, 40, 75, 120, 239, 358, Inf))

        agecat.matrix <- model.matrix(~ 0 + AGECAT, info)

        alply(agecat.matrix, 2, function(agevec)
        {
            agevec <- as.logical(agevec)
            R1 <- coord[, , 'R1', agevec]
            R2 <- coord[, , 'R2', agevec]

            R1.est <- geomorph:::estimate.missing(R1, method = 'TPS')
            R2.est <- geomorph:::estimate.missing(R2, method = 'TPS')

            out [, , 'R1', agevec] <<- R1.est
            out [, , 'R2', agevec] <<- R2.est
        })

        cat('\n') ## or else prompt stays after last progress bar

        list('info' = info, 'coord' = out)

    }
