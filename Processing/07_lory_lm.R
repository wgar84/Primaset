require(devtools)
require(shapes)
require(geomorph)
require(doMC)

require(plyr)

registerDoMC(cores = 32)

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

primates $ sym.coord <- aaply(primates $ coord, 3, Symmetrize, .parallel = TRUE)

primates $ sym.coord <- aperm(primates $ sym.coord, c(2, 3, 1))

dimnames(primates $ sym.coord) [3]

## LORY

primates $ sym.lory <- LORY(primates $ sym.coord, Aux $ single.tessel.38, TRUE)

