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
##' do this by species.
##'
##' But, to complete landmarks properly, it should also deal with outliers... shit.
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

        colnames(plat.info) [c(2:3, 6)] <- c('GEN', 'SPE', 'MSM')
        
        ## juntando ids
        
        cata.info $ ID <- laply(strsplit(as.character(cata.info $ ID), '_'),
                                function(L) paste(L[1], L[2], sep = '_'))

        plat.info $ ID <- paste(plat.info $ MSM, plat.info $ ID, sep = '_')
                            
        rownames(cata.info) <- NULL
        
        full.cata.info <- rbind(cata.info, homo.info)

        common.data <- rbind.fill(plat.info, full.cata.info)

        ## corrections (legacy)

        common.data $ SEX [common.data $ SEX == '0'] = NA
        common.data $ SEX [common.data $ SEX == ''] = NA
        common.data $ SEX [common.data $ SEX == '?F'] = 'F'
        common.data $ SEX [common.data $ SEX == '?M'] = 'M'
        common.data $ SEX [common.data $ SEX == 'sexo'] = NA

        common.data $ SEX <- factor(as.character(common.data $ SEX))

        common.data $ SUB [common.data $ SUB == '0'] = NA
        common.data $ SUB [common.data $ SUB == ''] = NA

        common.data $ GEN <- as.character(common.data $ GEN)
        common.data $ SPE <- as.character(common.data $ SPE)
        common.data $ SUB <- as.character(common.data $ SUB)
                
        common.data [, 'SUB'] [common.data [, 'SUB'] == 'aequatoriali'] = 'aequatorialis'
        common.data [, 'GEN'] [common.data [, 'GEN'] == 'Leontopithec'] = 'Leontopithecus'
        common.data [, 'SUB'] [common.data [, 'SUB'] == 'aequatoriali'] = 'aequatorialis'
        common.data [, 'SPE'] [common.data [, 'SPE'] == 'senicula'] = 'seniculus'
        common.data [, 'SUB'] [common.data [, 'SUB'] == 'senicula'] = 'seniculus'

        common.data [, 'SUB'] [common.data [, 'SUB'] == 'palliates'] = 'palliatus'

        common.data [, 'SUB'] [common.data [, 'GEN'] == 'Cacajao' &
                               common.data [, 'SPE'] == 'melanocephal' &
                               common.data [, 'SUB'] == 'melanocephal'] = 'melanocephalus'

        common.data [, 'SPE'] [common.data [, 'GEN'] == 'Cacajao' &
                               common.data [, 'SPE'] == 'melanocephal'] = 'melanocephalus'

        common.data [, 'SPE'] [common.data [, 'GEN'] == 'Saimiri' &
                               common.data [, 'SPE'] == 'oerstedi'] = 'oerstedii'

        common.data [, 'SUB'] [common.data [, 'GEN'] == 'Saimiri' &
                               common.data [, 'SPE'] == 'oerstedii' &
                               common.data [, 'SUB'] == 'oerstedi'] = 'oerstedii'

        common.data [, 'SPE'] [common.data [, 'GEN'] == 'Chiropotes' &
                               common.data [, 'SPE'] == 'satanas' &
                               common.data [, 'SUB'] == 'chiropotes'] = 'chiropotes'

        common.data [, 'SPE'] [common.data [, 'GEN'] == 'Chiropotes' &
                               common.data [, 'SPE'] == 'satanas' &
                               common.data [, 'SUB'] == 'utahicki'] = 'utahickae'

        common.data [, 'GEN'] [common.data [, 'GEN'] == 'Cebuella'] = 'Callithrix'
        
        common.data [, 'GEN'] [common.data [, 'GEN'] == 'Bunopithecus'] = 'Hoolock'

        common.data [, 'SPE'] [common.data [, 'GEN'] == 'Pithecia' &
                               common.data [, 'SPE'] == 'monacha'] = 'monachus'

        common.data [, 'SPE'] [common.data [, 'GEN'] == 'Alouatta' &
                               common.data [, 'SPE'] == 'fusca'] = 'guariba'

        ## still don't have everything I did, maybe

        common.data [common.data $ 'GEN' == 'Kasi', 'GEN'] <- 'Trachypithecus'
        
        common.data $ MSM <-
            gsub("[[:digit:]]", "", laply(strsplit(common.data $ ID, '_'), function(L) L[1]))

        common.data $ MSM <-
            gsub(' ', '', common.data $ MSM)
        
        ## drop one replicate to deal with outliers

        common.coord <- array(0, c(36, 3, nrow(common.data)))

        dimnames(common.coord) [1:2] <- dimnames(platyrrhine $ coord) [1:2]
        dimnames(common.coord) [[3]] <- common.data $ ID
        
        common.coord[, , 1:nrow(plat.info)] <- platyrrhine $ coord

        common.coord[, , (1:nrow(cata.info)) + nrow(plat.info)] <- catarrhine $ coord [, , 1, ]

        common.coord[, , (1:nrow(homo.info)) + nrow(plat.info) + nrow(cata.info)] <-
            homo $ coord [, , 1, ]

        print('finding outliers')

        i <- 'Alouatta'
        
        outliers <- c()
        for (i in unique(common.data $ GEN))
            {
                
                sset <- grepl(i, common.data $ GEN) # 0.15

                sset.coord <- common.coord [, , sset]

                if(any(common.data $ MISS [sset]))
                {
                    sset.complete <- estimate.missing(sset.coord, method = 'TPS')
                    cat('\n')
                } else {sset.complete <- sset.coord}
                
                sset.gpa <- procGPA(sset.complete)

                outliers[sset] <- sset.gpa$rho > 2 * sset.gpa $ rmsrho
            }

        out.data <- common.data [!outliers, ]
        out.coord <- common.coord [, , !outliers]

        ## maybe do this twice, I guess (or use a while)

        pdf('procHist.pdf')
        for (i in unique(out.data $ GEN))
            {
               
                sset <- grepl(i, out.data $ GEN) # 0.15

                sset.coord <- out.coord [, , sset]

                if(any(out.data $ MISS [sset]))
                {
                    sset.complete <- estimate.missing(sset.coord, method = 'TPS')
                    cat('\n')
                } else { sset.complete <- sset.coord }
                
                sset.gpa <- procGPA(sset.complete)

                hist(sset.gpa$rho, main = i)
                abline(v = 2 * sset.gpa $ rmsrho, lty = 'dashed', col = 'red')
                
            }
        dev.off(dev.cur())

    }
