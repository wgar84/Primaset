##' @title
##' CleanUpPlatyrrhini
##'
##' @description
##' This function cleans up Platyrrhini data from excel files.
##'
##' @param raw.list raw output from reading ascii files
##'
##' @return a list with two elements, info and coord
##'
##' @details will GlueSkulls
##'
##' @author Guilherme Garcia
##'
##' @importFrom plyr llply aaply ldply
##'
##' @export
##' @rdname CleanUpPlatyrrhini
##'
##' @seealso ReadPlatyrrhini

CleanUpPlatyrrhini <- function(raw.list)
{
    ## build 1rep list
    A1 <- raw.list $ A [, , grepl('\\.1', dimnames(raw.list $ A) [[3]])]
    Z1 <- raw.list $ Z [, , grepl('\\.1', dimnames(raw.list $ Z) [[3]])]

    dimnames(A1) [[3]] <- gsub('\\.1', '', dimnames(A1) [[3]])
    dimnames(Z1) [[3]] <- gsub('\\.1', '', dimnames(Z1) [[3]])

    coord.1r <-
        aaply(1:nrow(raw.list $ id), 1,
              function(i) GlueSkull(A1 [, , i], Z1 [, , i]))

    coord.1r <- aperm(coord.1r, c(2, 3, 1))

    id <- raw.list $ id
    id $ REP <- raw.list $ id $ REP.A > 1 & raw.list $ id $ REP.Z > 1
    id <- id [, c(1:10, 13)]

    ## two rep

    two.rep.id <- id $ ID [id $ REP]

    coord.2r <- array(0, c(dim(coord.1r) [1:2], length(two.rep.id), 2))

    A.R1 <-
        raw.list $ A [, , match(paste0(two.rep.id, '.1'),
                                dimnames(raw.list $ A) [[3]])]

    A.R2 <-
        raw.list $ A [, , match(paste0(two.rep.id, '.2'),
                                dimnames(raw.list $ A) [[3]])]

    Z.R1 <-
        raw.list $ Z [, , match(paste0(two.rep.id, '.1'),
                                dimnames(raw.list $ Z) [[3]])]

    Z.R2 <-
        raw.list $ Z [, , match(paste0(two.rep.id, '.2'),
                                dimnames(raw.list $ Z) [[3]])]

    R1 <- aaply(1:length(two.rep.id), 1, function(i)
        GlueSkull(A.R1 [, , i], Z.R1 [, , i]))

    R2 <- aaply(1:length(two.rep.id), 1, function(i)
        GlueSkull(A.R2 [, , i], Z.R2 [, , i]))

    coord.2r [, , , 1] <- aperm(R1, c(2, 3, 1))
    coord.2r [, , , 2] <- aperm(R2, c(2, 3, 1))

    dimnames(coord.2r) [1:2] <- dimnames(coord.1r) [1:2]
    dimnames(coord.2r) [[3]] <- two.rep.id
    dimnames(coord.2r) [[4]] <- c('R1', 'R2')

    dimnames(coord.1r) [[3]] <- id $ ID

    miss.lm <- aaply(aaply(coord.1r, c(1, 2), is.na), 3, any)

    id $ MISS <- miss.lm

    list('id' = id, 'coord' = coord.1r, 'rep' = coord.2r)

}
