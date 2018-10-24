require(devtools)
require(shapes)
require(geomorph)
require(rstan)
require(plyr)
require(RColorBrewer)

require(expm)

require(evolqg)
require(plotrix)

require(doMC)
registerDoMC(cores = 32)

## stan uses different backend now
rstan_options(auto_write = TRUE)
options(mc.cores = 32)

load('Primates/06_grouped.RData')
load('../Raw Data/Aux.RData')

## remove 6337, 7157 & 7255 and it'll be ok (some traits with negative determinants)

negadet <- as.character(primates $ info $ ID [c(6337, 7157, 7255)])

primates $ info <- subset(primates $ info, !ID %in% negadet)

primates $ coord <-
    primates $ coord [, , !dimnames(primates $ coord) [[3]] %in% negadet]

primates $ rep <-
    primates $ rep [, , , !dimnames(primates $ rep) [[4]] %in% negadet]

## have to estimate symmetric components first, idiot
## no duplicated IDs!

primates $ info $ ID <- as.character(primates $ info $ ID)

primates $ info [duplicated(primates $ info $ ID), 'ID'] <-
    paste(primates $ info [duplicated(primates $ info $ ID), 'ID'],
          LETTERS[1:16], sep = '_')

dimnames(primates $ coord) [[3]] <- primates $ info $ ID

dimnames(primates $ rep) [[4]] <- primates $ info $ ID [primates $ info $ rep]

primates $ info $ ID <- factor(primates $ info $ ID)

prima.sym <- list()

prima.sym $ coord <- aaply(primates $ coord, 3, Symmetrize, .parallel = TRUE)

prima.sym $ coord <- aperm(prima.sym $ coord, c(2, 3, 1))

dimnames(prima.sym $ coord) [3]

## RENAME OBJECTS

prima.info <- primates $ info

prima.info $ SPE <- as.character(prima.info $ SPE)

prima.info $ SPE [prima.info $ SPE == 'lagothricha' & !is.na(prima.info $ SPE)] <-
    'lagotricha'

prima.info $ GSP <- paste(prima.info $ GEN, prima.info $ SPE, sep = '_')

save(prima.info, file = 'Primates/Info.RData')

prima.raw <- primates [-2]

save(prima.raw, file = 'Primates/Raw.RData')

save(prima.sym, file = 'Primates/Sym.RData')

rm(primates)

## tem que arrumar a ordem do single.tessel.38

prima.aux <- list()

prima.aux $ distance.matrix <-
    matrix(unlist(strsplit(rownames(Aux $ def.hyp) [-1], '\\.')), nrow = 38, byrow = TRUE)

prima.aux $ distance.numbers <- array(0, dim(prima.aux $ distance.matrix))

for(i in 1:38)
    for (j in 1:2)
        prima.aux $ distance.numbers[i, j] <-
            which(grepl(prima.aux $ distance.matrix [i, j],
                        rownames(prima.sym $ coord))) [1]

## ED

prima.ed <- list()

prima.ed $ raw <-
    aaply(prima.sym $ coord, 3, EuclideanDistances, dists = prima.aux $ distance.numbers,
          .parallel = TRUE)

colnames(prima.ed $ raw) <- aaply(prima.aux $ distance.matrix, 1, paste, collapse = '.')

rownames(prima.ed $ raw) <- prima.info $ ID

save(prima.ed, file = 'Primates/ed.RData')

## LORY

prima.lory <- LORY(prima.sym $ coord, prima.aux $ distance.numbers, TRUE)

## which(aaply(prima.lory $ jacobians, 3:4, det, .parallel = TRUE) < 0, arr.ind = TRUE)

prima.lory $ local <-
    aaply(prima.lory $ jacobians, 4, Center2MeanJacobian, .parallel = TRUE)

prima.lory $ local <- t (prima.lory $ local)

dimnames(prima.lory $ local) <- list(prima.info $ ID, colnames(prima.ed $ raw))

save(prima.lory, file = 'Primates/LORY.RData')

save(prima.aux, file = 'Primates/aux.RData')

