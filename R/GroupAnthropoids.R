##' @title
##' GroupAnthropoids
##'
##' @description
##' This functions unifies Catarrhine, Platyrrhine and Homo databases into a single object.
##'
##' @param catarrhine Catarrhini DB, output of CleanUpCatarrhini
##' @param platyrrhine Platyrrhini DB, output of CleanUpPlatyrrhini
##' @param homo Homo DB, output of CleanUpHomo
##' @param prev.info data frame, info on previous iteration of consolidated DB
##' (from G. Garcia's phD)
##'
##' @return some list
##' 
##' @details This function will solve taxonomic issues resolved during G. Garcia's phD; examples
##' include subspecies lifted to species (e.g. Chiropotes, Gorilla).
##' Since this package expands upon the original DB, specimens previously unused will have their
##' taxonomic information updated according to the original resolution.
##'
##' This function also averages replicates when available. However, it retains the replicates
##' thus allowing estimating repeatabilities.
##'
##' However, to do this, the function should also complete missing landmarks; ideally, it should
##' ##' do this by species.
##'
##' So, first we correct info
##' 
##' @author Guilherme Garcia
##'
##' @seealso CleanUpCatarrhini CleanUpPlatyrrhini CleanUpHomo
##'
GroupAnthropoids <- function(catarrhine, platyrrhine, homo, prev)
    {
        ## Homo sapiens
        ## Pan troglodytes
        ## Pan paniscus
        ## Gorilla gorilla
        ## Pongo abelii  (não tinha na DB original?)
        ## Nomascus leucogenys (idem ao Pongo)
        ## Chlorocebus sabaeus
        ## Macaca mulatta
        ## Papio anubis
        ## Callithrix jacchus

        ## IS PM NSL NA
        
        plat.info <- platyrrhine $ id
        cata.info <- catarrhine $ info
        homo.info <- homo $ info

        ## juntando ids
        
        cata.info $ ID <- laply(strsplit(as.character(cata.info $ ID), '_'),
                                function(L) paste(L[1], L[2], sep = '_'))

        plat.info $ ID <- paste(plat.info $ MUSEUM, plat.info $ ID, sep = '_')
                            
        rownames(cata.info) <- NULL
        
        full.cata.info <- rbind(cata.info, homo.info)

        common.data <- rbind.fill(plat.info, full.cata.info)
        
        ### fill this crap up
        
        common.coord <- array(0, c(36, 3, nrow(common.data)))

        common.rep <- array(0, c(36, 3, 2, nrow(common.data)))

        cata.mshape <- aaply(catarrhine $ coord, 4, function(shrep) procGPA(shrep) $ mshape)
        

    }

