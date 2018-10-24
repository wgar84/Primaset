##' @title
##' GroupPrimates
##'
##' @description
##' This functions unifies Prosimian and Anthropoid databases in a single object.
##'
##' @param prosimians Prosimian DB, output of CleanUpProsimian
##' @param anthropoids Anthropoid DB, output of GroupAnthropoids
##'
##' @return list with four elements:
##' extinct: fossil and subfossil prosimians (list with two elements, info and coord),
##' coord: array with individual coordinate data
##' rep: array with individuals with replicated measurements,
##' info: specimen information
##' 
##' @author Guilherme Garcia
##'
##' @importFrom shapes procGPA
##' @importFrom geomorph estimate.missing
##' @importFrom plyr aaply alply rbind.fill
##' 
##' @seealso CleanUpCatarrhini CleanUpPlatyrrhini CleanUpHomo
##' 
GroupPrimates <- function(prosimians, anthropoids)
{

    ## some empty entries out there
    empty <-
        aaply(prosimians $ coord [, 1, 1, ], 2, function(c) sum(is.na(c)) > 30)

    raw.coord <- prosimians $ coord [, , , !empty]
    raw.id <- prosimians $ id [!empty, ]
    
    ## warnings because of replicate landmark names
    
    ## find me the halflings
    halflings <-
        aaply(raw.coord [, 1, 1, ], 2, function(c) sum(is.na(c)) > 21)

    ##
    onerep <-
        aaply(raw.coord [, , 1, ] - raw.coord [, , 2, ], 3, sum) == 0 |
        aaply(raw.coord [1:8, 1, 2, ], 2, function(c) all(is.na(c)))

    onerep[is.na(onerep)] <- FALSE
    
    onerep <- onerep | halflings
    
    missLM <-
        aaply(raw.coord [, 1, 1, ], 2, function(c) any(is.na(c)) & !all(is.na(c))) &
        !halflings

    colnames(anthropoids $ info) [9] <- 'LOC'
    colnames(anthropoids $ info) [13] <- 'MAJOR'
    
    colnames(raw.id)[1:8] <-
        c('FILE', 'ID', 'REG', 'MSM', 'GEN', 'SPE', 'SEX', 'LOC')

    right.side.miss <-
        aaply(raw.coord [23:44, 1, 1, ], 2, function(c) all(is.na(c))) & halflings
    
    
    ## arrumar essa bagunça

    sp.solo <- raw.id $ SPE

    sp.solo <- gsub('Chero', 'Cheiro', sp.solo)
    sp.solo <- gsub('Vaecia', 'Varecia', sp.solo)
    sp.solo <- gsub('PROLEMUR', 'Prolemur', sp.solo)
    sp.solo <- gsub('SIMUS', 'simus', sp.solo)
    sp.solo <- gsub('coquerelli', 'coquereli', sp.solo)
    sp.solo <- gsub('\\?', '', sp.solo)
    sp.solo <- gsub('callabarensis', 'calabarensis', sp.solo)
    sp.solo <- gsub('madasgariensis', 'madagascariensis', sp.solo)
    
    for(i in unique(raw.id $ GEN))
    {
        sp.solo <- gsub(i, '', sp.solo)
    }

    sp.split <- strsplit(sp.solo, ' ')
    
    pri.sp <- laply(sp.split, function(L) L [1])
    
    sec.sp <- laply(sp.split, function(L) L [2])

    tir.sp <- laply(sp.split, function(L) L [3])

    ## sec na means sp is on pri
        
    sec.sp[is.na(sec.sp)] <- pri.sp[is.na(sec.sp)]

    tir.sp[pri.sp != '' & !is.na(sec.sp)] <- sec.sp [pri.sp != '' & !is.na(sec.sp)]

    sec.sp[pri.sp != '' & !is.na(sec.sp)] <- pri.sp [pri.sp != '' & !is.na(sec.sp)]
    
    sec.sp[sec.sp == 'sp.'] <- NA
    
    collections <-
        c('AMNH','FMNH', 'USNM', 'MCZ', 'FPDuke',
          'RMNH', 'AZM', 'MHNBV', 'ZMB', 'MNHN')

    coll.vec <- c()
    for(i in collections)
        coll.vec[grepl(i, raw.id $ ID)] <- i

    ## sex mess

    fixed.sex <- c()

    fixed.sex[grepl('F', raw.id $ SEX)] <- 'F'
    fixed.sex[grepl('M', raw.id $ SEX)] <- 'M'

    fossil.genera <-
        c('Archaeolemur', 'Mesopropithecus', 'Paleopropithecus',
          'Babakotia', 'Pachylemur', 'Megaladapis', 'Hadropithecus')

    who.fossil <- raw.id $ GEN %in% fossil.genera

    fixed.info <- data.frame('ID' = raw.id $ ID,
                             'GEN' = as.character(raw.id $ GEN),
                             'SPE' = sec.sp,
                             'SUB' = tir.sp,
                             'GROUP' = rep(NA, nrow(raw.id)),
                             'MSM' = coll.vec,
                             'IDORI' = raw.id $ SPE, ## pq é aquela bagunça
                             'SEX' = fixed.sex,
                             'LOC' = raw.id $ LOC,
                             'DIET' = rep(NA, nrow(raw.id)),
                             'REP' = !onerep,
                             'MISS' = missLM,
                             'HALF' = halflings,
                             'RIGHT' = right.side.miss,
                             'FOSSIL' = who.fossil,
                             'MAJOR' = rep('Prosimian', nrow(raw.id)),
                             'FILE' = raw.id $ FILE
                             )

    fixed.info <- fixed.info[!is.na(coll.vec), ] ## remover esses (microscribe test)
    raw.coord <- raw.coord [, , , !is.na(coll.vec)]

    pro.coord <- array(0, c(36, 3, 2, nrow(fixed.info)))

    fixed.info ['AMNH207949', ] $ REP <- FALSE ## miss 2 lms from R2, only one in sp.

    fixed.info $ REP [fixed.info $ GEN == 'Prolemur'] <- FALSE
    ## no idea why I can't use them
    
    ## halfling: duplicate sides
    
    ## destra
    half.dup.D <-
        aaply(raw.coord [23:44, , 1,
                         which(fixed.info $ HALF & !fixed.info $ RIGHT)], 3,
              function(S)
              {
                  E <- S %*% diag(c(-1, 1, 1))
                  dimnames(E) <- dimnames(S)
                  rownames(E) <- gsub('-D', '-E', rownames(E))
                  GlueSkull(E, S)
              })

    half.dup.D <- aperm(half.dup.D, c(2, 3, 1))
    pro.coord[, , 1, which(fixed.info $ HALF & !fixed.info $ RIGHT)] <- half.dup.D

    ## sinistra
    half.dup.E <-
        aaply(raw.coord [1:22, , 1,
                         which(fixed.info $ HALF & fixed.info $ RIGHT)], 3,
              function(S)
              {
                  D <- S %*% diag(c(-1, 1, 1))
                  dimnames(D) <- dimnames(S)
                  rownames(D) <- gsub('-E', '-D', rownames(D))
                  GlueSkull(S, D)
              })

    half.dup.E <- aperm(half.dup.E, c(2, 3, 1))
    pro.coord[, , 1, which(fixed.info $ HALF & fixed.info $ RIGHT)] <- half.dup.E

    dimnames(pro.coord) [1:2] <- dimnames(half.dup.E) [1:2]

    ## grudar o resto

    full.sk <- aaply(raw.coord[, , , !fixed.info $ HALF & fixed.info $ REP], 3:4,
                     function(S) GlueSkull(S[1:22, ], S[23:44, ]))

    full.sk <- aperm(full.sk, c(3, 4, 1, 2))

    pro.coord[, , , !fixed.info $ HALF & fixed.info $ REP] <- full.sk

    ## não replicados

    norep.sk <- aaply(raw.coord[, , 1, !fixed.info $ HALF & !fixed.info $ REP], 3,
                     function(S) GlueSkull(S[1:22, ], S[23:44, ]))

    norep.sk <- aperm(norep.sk, c(2, 3, 1))

    pro.coord[, , 1, !fixed.info $ HALF & !fixed.info $ REP] <- norep.sk

    ## extinct: keep separate, do nothing
    
    extinct <- list('coord' = pro.coord [, , , fixed.info $ FOSSIL],
                    'info' = subset(fixed.info, FOSSIL) [, -(13:14)])
    
    pro.coord <- pro.coord [, , , !fixed.info $ FOSSIL]
    fixed.info <- subset(fixed.info, !FOSSIL)
    
    ## missing: complete by genus, find outliers

    outliers <- c()
    
    for (i in unique(fixed.info $ GEN))
        {
            print(i)
            sset <- grepl(i, fixed.info $ GEN) # 0.15
            out.local <- rep(TRUE, sum(sset))
            if(sum(sset) > 1)
                {
                    repeat
                    {
                        sset.coord <- pro.coord [, , 1, sset]  [, , out.local]
                        
                        if(any(fixed.info $ MISS [sset]))
                        {
                            sset.complete <-
                                estimate.missing(sset.coord, method = 'TPS')
                            cat('\n')
                        } else {sset.complete <- sset.coord}
                        
                        sset.gpa <- procGPA(sset.complete)
                                            
                        out.it <- sset.gpa$rho > 1.96 * sset.gpa $ rmsrho
                        
                        if(any(out.it)) {
                            out.local [out.local] [which(out.it)] <- FALSE
                        } else {
                            pro.coord [, , 1, sset]  [, , out.local] <- sset.complete
                            outliers[sset] <- !out.local
                            break
                        }
                    }
                }
        }

    outliers[is.na(outliers)] <- FALSE

    ## replicates run

    rep.info <- subset(fixed.info, REP)
    rep.coord <- pro.coord[, , 2, fixed.info $ REP]
    
    outliers2 <- c()
    
    print('outliers, second run')
    
    for (i in unique(rep.info $ GEN))
        {
            print(i)
            sset <- grepl(i, rep.info $ GEN) 
            out.local <- rep(TRUE, sum(sset))
            if(sum(sset) > 1)
            {
                repeat
                {
                    sset.coord <- rep.coord [, , sset]  [, , out.local]

                    if(any(rep.info $ MISS [sset]))
                    {
                        sset.complete <- estimate.missing(sset.coord, method = 'TPS')
                        cat('\n')
                    } else {sset.complete <- sset.coord}
                    
                    sset.gpa <- procGPA(sset.complete)
                    
                    out.it <- sset.gpa$rho > 1.96 * sset.gpa $ rmsrho

                    if(any(out.it)) {
                        out.local [out.local] [which(out.it)] <- FALSE
                    } else {
                        rep.coord [, , sset]  [, , out.local] <- sset.complete
                        outliers2[sset] <- !out.local
                        break
                    }
                }
            }
        }

    outliers2 [is.na(outliers2)] <- FALSE
    
    outliers [fixed.info $ REP] <- outliers [fixed.info $ REP] | outliers2

    pro.coord[, , 2, fixed.info $ REP] <- rep.coord

    pro.coord <- pro.coord [, , , !outliers]
    fixed.info <- fixed.info [!outliers, ]

    ## averaging those we can
    av.coord <-
        aaply(pro.coord [, , , fixed.info $ REP], 4, function(S) procGPA(S) $ mshape)

    av.coord <- aperm(av.coord, c(2, 3, 1))

    final.coord <- array(0, c(36, 3, nrow(fixed.info)))

    final.coord [, , fixed.info $ REP] <- av.coord

    ## onerep: use one rep only
    final.coord [, , !fixed.info $ REP] <- pro.coord [, , 1, !fixed.info $ REP]

    dimnames(final.coord) [1:2] <- dimnames(pro.coord) [1:2]

    ## some individuals are reflected (how come I didn't pick this before?)
    refl <- c("RMNH16899", "RMNH9917ZMA",  "RMNH17202ZMA", "RMNH10728ZMA",
              "RMNH10754ZMA", "RMNH10727ZMA", "RMNH7787ZMA",  "RMNH7786ZMA",
              "RMNH046001a",  "RMNH046001b", "RMNH19181ZMA", "RMNH5544ZMA",
              "RMNH5484ZMA",  "RMNH10341ZMA", "RMNH21895ZMA", "RMNH045004e") 

    final.coord [, , rownames(fixed.info) %in% refl] <-
        aperm(aaply(final.coord [, , rownames(fixed.info) %in% refl], 3,
              function(S) S %*% diag(c(-1, 1, 1))), c(2, 3, 1))

    final.rep <- pro.coord [, , , fixed.info $ REP]
    
    ## test <- procGPA(final.coord)
    
    ## join all primates

    lm.order <- match(rownames(final.coord), rownames(anthropoids $ coord))

    antro.coord <- anthropoids $ coord [lm.order, , ]
    antro.rep <- anthropoids $ rep [lm.order, , , ]

    prima.coord <- array(0, c(36, 3, dim(antro.coord)[3] + dim(final.coord)[3]))
    prima.rep <- array(0, c(36, 3, 2, dim(antro.rep)[4] + dim(final.rep)[4]))

    prima.coord [, , 1:(dim(antro.coord)[3])] <- antro.coord
    prima.coord [, , 1:(dim(final.coord)[3]) + dim(antro.coord)[3]] <- final.coord

    prima.rep [, , , 1:(dim(antro.rep)[4])] <- antro.rep
    prima.rep [, , , 1:(dim(final.rep)[4]) + dim(antro.rep)[4]] <- final.rep

    
    prima.info <- rbind.fill(anthropoids $ info, fixed.info)

    dimnames(prima.coord) [1:2] <-
                              dimnames(prima.rep) [1:2] <-
                                                      dimnames(final.coord) [1:2]

    dimnames(prima.rep) [[3]] <- c('R1', 'R2')
    dimnames(prima.rep) [[4]] <- prima.info $ ID [prima.info $ REP]

    dimnames(prima.coord) [[3]] <- prima.info $ ID

    ## still some sp name corrections

    prima.info $ SPE <- gsub('macrocephalu', 'macrocephalus', prima.info $ SPE)
    prima.info $ SPE <- gsub('cassiquiaren', 'cassiquiarensis', prima.info $ SPE)
    prima.info $ SPE <- gsub('klossi', 'klossii', prima.info $ SPE)
    prima.info $ SPE <- gsub('thibethana', 'thibetana', prima.info $ SPE)
    prima.info $ SPE <- gsub('aequatoriali', 'aequatorialis', prima.info $ SPE)
    prima.info $ SPE <- gsub('xanthosterno', 'xanthosternos', prima.info $ SPE)
    prima.info $ SPE <- gsub('nigrivitattu', 'nigrivitattus', prima.info $ SPE)
    prima.info $ SPE <- gsub('sandfordi', 'sanfordi', prima.info $ SPE)
    prima.info $ SPE <- gsub('ruffifrons', 'rufifrons', prima.info $ SPE)

    prima.info $ SPE[prima.info $ SPE == 'klossi'] <- 'klossii' 
    prima.info $ SPE[prima.info $ SPE == ''] <- NA

    prima.info <-
        as.data.frame(lapply(prima.info, function(c) if(is.factor(c)) factor(c) else c))

    ## fucking POLHEMUS
    prima.coord [, , prima.info $ MAJOR == 'Platyrrhini'] <- 
        prima.coord [, , prima.info $ MAJOR == 'Platyrrhini'] * 10

    prima.rep [, , , subset(prima.info, REP) $ MAJOR == 'Platyrrhini'] <- 
        prima.rep [, , , subset(prima.info, REP) $ MAJOR == 'Platyrrhini'] * 10

    
    
    return(list('extinct' = extinct, 'info' = prima.info,
                'coord' = prima.coord, 'rep' = prima.rep))
}
