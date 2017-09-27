require(devtools)
require(roxygen2)
require(shapes)
require(plyr)
require(rgl)
require(geomorph)
require(geometry)
require(viridis)
require(doMC)
registerDoMC(cores = 3)
require(plotrix)

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

### auxiliary stuff

calomys.wgen $ sym.gpa <- procGPA(calomys.wgen $ sym.coord)

calomys.wgen $ symmetry.corr <- cbind(9:22, 23:36)

calomys.wgen $ saggital <- rownames(calomys.wgen $ sym.coord) [1:8]

calomys.wgen $ bone.tri <-
    matrix(ncol = 3,
           data = c('MT', 'PNS', 'PM',     #  Oral      
                    'IS', 'PM', 'PNS',     #  Oral      
                    'PM', 'ZS', 'ZI',      #  Oral      
                    'PM', 'FM', 'NSL',     #  Oral      
                    'NSL', 'NA', 'FM',     #  Nasal     
                    'NSL', 'PM', 'IS',     #  Nasal     
                    'TSP', 'FM', 'PM',     #  Zygomatic 
                    'ZS', 'ZI', 'ZYGO',    #  Zygomatic 
                    'ZYGO', 'TS', 'TSP',   #  Zygomatic 
                    'FM', 'ZS', 'ZI',      #  Zygomatic 
                    'FM', 'PT', 'TSP',     #  Zygomatic 
                    'APET', 'BA', 'JP',    #  Base      
                    'PNS', 'APET', 'TS',   #  Base      
                    'APET', 'JP', 'PEAM',  #  Base      
                    'APET', 'EAM', 'TS',   #  Base      
                    'NA', 'BR', 'PT',      #  Vault     
                    'FM', 'NA', 'PT',      #  Vault     
                    'PT', 'AS', 'LD',      #  Vault     
                    'PT', 'BR', 'LD',      #  Vault     
                    'AS', 'TSP', 'PT',     #  Vault     
                    'AS', 'TSP', 'TS',     #  Vault     
                    'OPI', 'AS', 'JP',     #  Vault     
                    'LD', 'OPI', 'AS'      #  Vault     
                    ))

calomys.wgen $ bone.region <-
    c('Oral', 'Oral', 'Oral', 'Oral', 'Nasal', 'Nasal',
      'Zygomatic', 'Zygomatic', 'Zygomatic', 'Zygomatic', 'Zygomatic',
      'Base', 'Base', 'Base', 'Base',
      'Vault', 'Vault', 'Vault', 'Vault', 'Vault', 'Vault', 'Vault', 'Vault')

calomys.wgen $ bone.tri.sym <-
    rbind(
        matrix(ifelse(calomys.wgen $ bone.tri %in% calomys.wgen $ saggital,
                      calomys.wgen $ bone.tri, paste0(calomys.wgen $ bone.tri, '-R')),
               byrow = T, ncol = 3), 
        matrix(ifelse(calomys.wgen $ bone.tri %in% calomys.wgen $ saggital,
                      calomys.wgen $ bone.tri, paste0(calomys.wgen $ bone.tri, '-L')),
               byrow = T, ncol = 3))

calomys.wgen $ tri.sym.num <-
    matrix(match(calomys.wgen $ bone.tri.sym, rownames(calomys.wgen $ sym.coord)), ncol = 3)

sort.sum <-
    order(rowSums(calomys.wgen $ tri.sym.num [1:(nrow(calomys.wgen $ tri.sym.num)/2), ]))

sort.sum <- c(sort.sum, sort.sum + nrow(calomys.wgen $ tri.sym.num)/2)

calomys.wgen $ tri.sym.num <- calomys.wgen $ tri.sym.num [sort.sum, ]

calomys.wgen $ bone.tri.sym <- calomys.wgen $ bone.tri.sym [sort.sum, ]

#text3d(calomys.wgen $ sym.gpa $ mshape,
#       texts = rownames(calomys.wgen $ sym.coord))

coltest <- viridis(nrow(calomys.wgen $ tri.sym.num) / 2, option = 'C')

coltest <- c(coltest, coltest)

for(i in 1:nrow(calomys.wgen $ tri.sym.num))
    triangles3d(calomys.wgen $ sym.gpa $ mshape[calomys.wgen $ tri.sym.num[i, ], ],
                color = coltest[i])

### seems nice

calomys.wgen $ lory.out <- LORY(calomys.wgen $ sym.coord, calomys.wgen $ tri.sym.num, TRUE)

calomys.wgen $ local <- calomys.wgen $ lory.out $ local [, 1:23]

calomys.wgen $ age <- calomys.wgen $ info $ AGE
