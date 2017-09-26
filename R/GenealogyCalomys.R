##' @title
##' GenealogyCalomys
##'
##' @description
##' Creates inverse genealogy matrix for Calomys data.
##'
##' @param filled output from MissLMCalomys
##'
##' @return List with three elements:
##' coordinate data (coord), specimen information (info) and genealogy information (gen),
##' as produced by inverseA from MCMCglmm; thus, gen is itselg a list with three elements.
##'
##' @seealso inverseA
##'
##' @author Guilherme Garcia
##'
##' @importFrom MCMCglmm inverseA
##' @importFrom plyr aaply
##'
##' @export
##' @rdname GenealogyCalomys

GenealogyCalomys <- function(filled)
    {
        info <- filled $ info
        ped <- info[, c('ID', 'DAM', 'SIR')]
        p.gen <- c(unique(ped $ DAM), unique(ped $ SIR))
        p.gen <- cbind(p.gen, rep(NA, length(p.gen)), rep(NA, length(p.gen)))

        ped <- as.matrix(ped)
        ped <- aaply(ped, 1, gsub, pattern = ' ', replace = '')

        full.ped <- rbind(p.gen, ped)

        gen <- inverseA(full.ped)

        list('info' = info, 'coord' = filled $ coord, 'gen' = gen)
    }


