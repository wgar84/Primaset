##' @title
##' GroupAnthropoids
##'
##' @description
##' This functions unifies Catarrhine, Platyrrhine and Homo databases into a single object.
##'
##' @param platyrrhine Platyrrhini DB, output of CleanUpPlatyrrhini
##' @param catarrhine Catarrhini DB, output of CleanUpCatarrhini
##' @param homo Homo DB, output of CleanUpHomo
##'
##' @return some list
##' 
##' @details This function will solve taxonomic issues resolved during G. Garcia's phD; examples
##' include subspecies lifted to species (e.g. Chiropotes).
##' Since this package expands upon the original DB, specimens previously unused will have their
##' taxonomic information updated according to the original resolution.
##'
##' This function also averages replicates when available. However, it retains the replicates
##' thus allowing estimating repeatabilities.
##'
##' However, to do this, the function also completes missing landmarks; it does so by genus.
##'
##' It also deals with outliers, using the distribution of Procrustes distances for each genus.
##' 
##' @author Guilherme Garcia
##'
##' @importFrom shapes procGPA
##' @importFrom geomorph estimate.missing
##' @importFrom plyr aaply alply
##' 
##' @seealso CleanUpCatarrhini CleanUpPlatyrrhini CleanUpHomo
##'
GroupAnthropoids <- function(platyrrhine, catarrhine, homo)
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

        plat.info $ Group <- rep('Platyrrhini', nrow(plat.info))
        cata.info $ Group <- rep('Catarrhini', nrow(cata.info))
        homo.info $ Group <- rep('Homo', nrow(homo.info))
                       
        colnames(plat.info) [c(2:3, 6)] <- c('GEN', 'SPE', 'MSM')
        
        ## juntando ids
        
        cata.info $ ID <- laply(strsplit(as.character(cata.info $ ID), '_'),
                                function(L) paste(L[1], L[2], sep = '_'))

        plat.info $ ID <- paste(plat.info $ MSM, plat.info $ ID, sep = '_')
                            
        rownames(cata.info) <- NULL
        
        full.cata.info <- rbind(cata.info, homo.info)

        common.data <- rbind.fill(plat.info, full.cata.info)

        ## corrections (legacy)

        print('legacy corrections')
        
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

        ## now its cool, I think
        
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
        
        outliers <- c()

        for (i in unique(common.data $ GEN))
        {
            print(i)
            sset <- grepl(i, common.data $ GEN) # 0.15
            out.local <- rep(TRUE, sum(sset))
            repeat
            {
                sset.coord <- common.coord [, , sset]  [, , out.local]

                if(any(common.data $ MISS [sset]))
                {
                    sset.complete <- estimate.missing(sset.coord, method = 'TPS')
                    cat('\n')
                } else {sset.complete <- sset.coord}

                sset.gpa <- procGPA(sset.complete)
                
                out.it <- sset.gpa$rho > 1.96 * sset.gpa $ rmsrho

                if(any(out.it)) {
                    out.local [out.local] [which(out.it)] <- FALSE
                } else {
                    outliers[sset] <- !out.local
                    break
                }
            }
        }

        ## outliers in second replicates

        common.coord[, , (1:nrow(cata.info)) + nrow(plat.info)] <- catarrhine $ coord [, , 2, ]

        common.coord[, , (1:nrow(homo.info)) + nrow(plat.info) + nrow(cata.info)] <-
            homo $ coord [, , 2, ]

        outliers2 <- outliers

        print('outliers, second run')
        
        for (i in unique(subset(common.data, Group != 'Platyrrhini') $ GEN))
        {
            print(i)
            sset <- grepl(i, common.data $ GEN) # 0.15
            out.local <- rep(TRUE, sum(sset))
            repeat
            {
                sset.coord <- common.coord [, , sset]  [, , out.local]

                if(any(common.data $ MISS [sset]))
                {
                    sset.complete <- estimate.missing(sset.coord, method = 'TPS')
                    cat('\n')
                } else {sset.complete <- sset.coord}

                sset.gpa <- procGPA(sset.complete)
                
                out.it <- sset.gpa$rho > 1.96 * sset.gpa $ rmsrho

                if(any(out.it)) {
                    out.local [out.local] [which(out.it)] <- FALSE
                } else {
                    outliers2[sset] <- !out.local
                    break
                }
            }
        }

        common.data $ outliers <- outliers | outliers2
        
        ## aqueles callithrix que não tem sp
        common.data $ outliers[common.data $ SPE == 'RIO'] <- TRUE
        
        common.data $ REP[is.na(common.data $ REP)] <- TRUE

        ## replicates
        
        ## Platyrrhini

        rep.coord <- array(0, c(36, 3, 2, sum(common.data $ REP & !common.data $ outliers)))

        plat.rcoord.tmp <- aperm(platyrrhine $ rep, c(1, 2, 4, 3))

        rcoord.index <- which(subset(common.data, REP) $ Group == 'Platyrrhini' &
                              !subset(common.data, REP) $ outliers)
        
        rep.coord[, , , 1:length(rcoord.index)] <- plat.rcoord.tmp [, , , rcoord.index]

        ## Catarrhini

        rcoord.index.2 <- which(!subset(common.data, Group == 'Catarrhini') $ outliers)

        rep.coord[, , , length(rcoord.index) + (1:length(rcoord.index.2))] <-
            catarrhine $ coord [, , , rcoord.index.2]

        ## Homo

        rcoord.index.3 <- which(!subset(common.data, Group == 'Homo') $ outliers)
        
        rep.coord[, , , length(rcoord.index) +
                        length(rcoord.index.2) +
                        (1:length(rcoord.index.3))] <- homo $ coord [, , , rcoord.index.3]

        dimnames(rep.coord)[1:3] <- dimnames(catarrhine $ coord) [1:3]
        dimnames(rep.coord)[[4]] <- subset(common.data, REP & !outliers) $ ID
        
        ## full dataset
        
        common.coord <- array(0, c(36, 3, sum(!common.data $ outliers)))

        ## platyrrhini (complete but not average)

        print('estimating missing landmarks, averaging replicates')
        
        plat.sset <- which(common.data $ Group == 'Platyrrhini' &
                           !common.data $ outliers)
        
        plat.coord <- platyrrhine $ coord [, , plat.sset]

        plat.matrix <- model.matrix(~ 0 + GEN, data = common.data[plat.sset, ])

        alply(plat.matrix, 2, function(genus)
        {
            genus <- as.logical(genus)
            plat.coord [, , genus] <<- estimate.missing(plat.coord [, , genus], 'TPS')
            cat('\n')
        })

        common.coord[, , 1:dim(plat.coord)[3]] <- plat.coord

        ## catarrhini (complete & average)
       
        cata.sset <- which(!subset(common.data, Group == 'Catarrhini') $ outliers)
        
        cata.coord <- catarrhine $ coord [, , , cata.sset]

        cata.matrix <- model.matrix(~ 0 + GEN, data = subset(common.data,
                                                             Group == 'Catarrhini')[cata.sset, ])

        alply(cata.matrix, 2, function(genus)
        {
            genus <- as.logical(genus)
            if(any(subset(common.data, Group == 'Catarrhini')[cata.sset, 'MISS'] [genus]))
            {
                cata.coord [, , 1, genus] <<- estimate.missing(cata.coord [, , 1, genus], 'TPS')
                cata.coord [, , 2, genus] <<- estimate.missing(cata.coord [, , 2, genus], 'TPS')
            }
            cat('\n')
        })

        cata.mcoord <- aaply(cata.coord, 4, function(shp) procGPA(shp) $ mshape)

        cata.mcoord <- aperm(cata.mcoord, c(2, 3, 1))
        
        common.coord [, , dim(plat.coord)[3] + (1:dim(cata.mcoord)[3])] <- cata.mcoord
            
        ## homo (complete & average)
        homo.sset <- which(!subset(common.data, GEN == 'Homo') $ outliers)
        
        homo.coord <- homo $ coord [, , , homo.sset]

        homo.coord[, , 1, ] <- estimate.missing(homo.coord[, , 1, ], 'TPS')
        cat('\n')
        homo.coord[, , 2, ] <- estimate.missing(homo.coord[, , 2, ], 'TPS')
        cat('\n')

        homo.mcoord <- aaply(homo.coord, 4, function(shp) procGPA(shp) $ mshape)

        homo.mcoord <- aperm(homo.mcoord, c(2, 3, 1))
        
        common.coord [, , dim(plat.coord)[3] + dim(cata.mcoord)[3] + (1:dim(homo.mcoord)[3])] <-
            homo.mcoord

        dimnames(common.coord)[1:2] <- dimnames(rep.coord)[1:2]
        dimnames(common.coord)[[3]] <- subset(common.data, !outliers) $ ID

        return(
            list('info' = subset(common.data, !outliers)[, -15], # drop outliers col
                 'coord' = common.coord,
                 'rep' = rep.coord)
        )
    }
