##' @title
##' CleanUpCalomys
##' 
##' @description
##' This function wraps up our work on cleaning up the Calomys data.
##' Repeated entries are removed, all individuals have the same number of landmarks,
##' that kind of stuff.
##' This is the second step towards a nice crispy Calomys DB.
##'
##' @param raw output obtained from ReadCalomys, which basically just reads excel files
##' @param id id from individuals (obtained from previous work with this DB)
##' 
##' @return
##' a list with two elements, id and coord, the cleaned up version of what came in
##'
##' @details
##' this function already glues left and right hemiskulls
##' 
##' @author Guilherme Garcia
##'
##' @importFrom plyr aaply
##'
##' @seealso ReadCalomys
##'
##' @examples
##' \dontrun{
##' }

CleanUpCalomys <- function(xls.list, id)
{
    Ns <- laply(xls.list, function(L) dim(L [[1]]) [4])
    ind <- unlist(laply(xls.list, function(L) dimnames(L [[1]]) [4]))
    
    lms.d <- dimnames(xls.list [[1]] [[1]]) [[1]]
    lms.e <- dimnames(xls.list [[1]] [[2]]) [[1]]

    coord.d <- array(0, c(length(lms.d), 3, 2, sum(Ns)))
    coord.e <- array(0, c(length(lms.e), 3, 2, sum(Ns)))

    dimnames(coord.d) <- list(lms.d, c('X', 'Y', 'Z'), c('R1', 'R2'), ind)
    dimnames(coord.e) <- list(lms.e, c('X', 'Y', 'Z'), c('R1', 'R2'), ind)
    
    counter <- 0
    for(i in 1:length(Ns))
    {
        coord.d [, , , counter+(1:Ns[i])] <- xls.list [[i]] $ D
        coord.e [, , , counter+(1:Ns[i])] <- xls.list [[i]] $ E

        counter <- counter + Ns[i]
    }

    ## empty entries
    coord.d <- coord.d [, , , !(ind == '')]
    coord.e <- coord.e [, , , !(ind == '')]
    ind <- ind[!(ind == '')]

    ## matching
    id $ num <- paste(id $ dam, id $ lit, id $ sib, sep = '-')

    ind[!ind %in% id $ num] [1] <- '69-6-D' ## star, no idea why

    coord.d <- coord.d[, , , ind %in% id $ num]
    coord.e <- coord.e[, , , ind %in% id $ num]
    ind <- ind[ind %in% id $ num]

    coord.d <- coord.d [, , , order(dimnames(coord.d) [[4]])]
    coord.e <- coord.e [, , , order(dimnames(coord.e) [[4]])]
    id <- id[order(id $ num), ]

    coord <- aaply(1:2, 1,
                   function(i) GlueSkulls(coord.d [, , i, ], coord.e [, , i, ]))

    coord <- aperm(coord, c(2, 3, 1, 4))

    dimnames(coord)[3:4] <- list(c('R1', 'R2'), id $ num)

    list('id' = id, 'coord' = coord)
    
}   
