require(devtools)
require(roxygen2)
require(shapes)
require(plyr)
require(rgl)
require(geomorph)
require(geometry)

load('Calomys/02_clean_up.RData')

calomys.complete <- MissLMCalomys(calomys.cleanup)

calomys.wgen <- GenealogyCalomys(calomys.complete)

str(calomys.wgen)

dimnames(calomys.wgen $ coord)[[1]] <-
    gsub('-D', '-R', dimnames(calomys.wgen $ coord)[[1]])

dimnames(calomys.wgen $ coord)[[1]] <-
    gsub('-E', '-L', dimnames(calomys.wgen $ coord)[[1]])

### got ourselves a reverted dude here
alledgedly.right <-
    calomys.wgen $ coord [grepl('-R', dimnames(calomys.wgen $ coord)[[1]]), , , 318]

calomys.wgen $ coord [grepl('-R', dimnames(calomys.wgen $ coord)[[1]]), , , 318] <-
    calomys.wgen $ coord [grepl('-L', dimnames(calomys.wgen $ coord)[[1]]), , , 318]

calomys.wgen $ coord [grepl('-L', dimnames(calomys.wgen $ coord)[[1]]), , , 318] <-
    alledgedly.right

### average replicates

calomys.wgen $ mean.coord <-
    aaply(calomys.wgen $ coord, 4, function(shrep) procGPA(shrep) $ mshape)

calomys.wgen $ mean.coord <- aperm(calomys.wgen $ mean.coord, c(2, 3, 1))

dimnames(calomys.wgen $ mean.coord)[1:2] <-
    dimnames(calomys.wgen $ coord)[1:2]

### symmetric component

dimnames(calomys.wgen $ sym.coord) [[1]]

calomys.wgen $ sym.coord <-
    aaply(calomys.wgen $ mean.coord, 3, Symmetrize)

calomys.wgen $ sym.coord <- aperm(calomys.wgen $ sym.coord, c(2, 3, 1))

calomys.wgen $ sym.cs <- aaply(calomys.wgen $ sym.coord, 3, centroid.size)

delaunayn(gpa.test $ mshape [1:22, ])
