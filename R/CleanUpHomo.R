##' @title
##' CleanUpHomo
##' 
##' @description
##' This function cleans up Homo data from excel files.
##' 
##' @param xls.list output obtained from ReadHomo
##'
##' @return a list with two elements, id and coord
##'
##' @details will GlueSkulls
##' 
##' @author Guilherme Garcia
##'
##' @importFrom plyr llply aaply ldply
##'
##' @seealso ReadHomo

CleanUpHomo <- function(xls.list)
{
    A <- llply(xls.list,
               function(L) suppressWarnings(aaply(L $ A, 1:4, as.numeric)))

    Z <- llply(xls.list,
               function(L) suppressWarnings(aaply(L $ Z, 1:4, as.numeric)))

    emptyA <-
        suppressWarnings(unlist(llply(A,
                                      function(L) aaply(L, 4,
                                                        function(m) all(is.na(m))))))

    emptyZ <-
        suppressWarnings(unlist(llply(Z,
                                      function(L) aaply(L, 4,
                                                        function(m) all(is.na(m))))))
    empty.any <- emptyA | emptyZ
    
    id <- ldply(xls.list, function(L) L $ id)
    
    names(empty.any) <- id $ numero

    Ns <- laply(A, function(L) dim(L) [4])
   
    coordA <- array(0, c(dim(A [[1]])[1:3], sum(Ns)))
    coordZ <- array(0, c(dim(Z [[1]])[1:3], sum(Ns)))

    dimnames(coordA)[[1]] <- gsub(' ', '', dimnames(A [[1]]) [[1]]) ### ops ops
    dimnames(coordA)[2:3] <- dimnames(A [[1]]) [2:3]
    dimnames(coordZ)[1:3] <- dimnames(Z [[1]]) [1:3]
    dimnames(coordA)[[4]] <- dimnames(coordZ)[[4]] <- id $ numero

    counter <- 0
    for(i in 1:length(Ns))
    {
        coordA [, , , counter+(1:Ns[i])] <- A [[i]]
        coordZ [, , , counter+(1:Ns[i])] <- Z [[i]]

        counter <- counter + Ns[i]
    }

    ## tirar vazios
    id <- id [!empty.any, ]
    coordA <- coordA [, , , !empty.any]
    coordZ <- coordZ [, , , !empty.any]

    ## ficar só com os que o fino usou
    lista.fino <- read.csv('Homo/lista_fino.csv')
    na.lista.do.fino <- id $ numero %in% lista.fino $ numero
    
    id <- id [na.lista.do.fino, ]
    coordA <- coordA [, , , na.lista.do.fino ]
    coordZ <- coordZ [, , , na.lista.do.fino ]

    lista.fino <- lista.fino [lista.fino $ numero %in% id $ numero, ]
    lista.fino <- lista.fino[order(lista.fino $ numero), ]
    lista.fino <- lista.fino[!duplicated(lista.fino $ numero), ]
    
    ord.id <- order(as.character(id $ numero))

    id <- id[ord.id, ]
    coordA <- coordA [, , , ord.id]
    coordZ <- coordZ [, , , ord.id]
    
    join.line.A <- rownames(coordA) %in% rownames(coordZ)
    join.line.Z <- rownames(coordZ) %in% rownames(coordA)

    screwed.join.A <- which(aaply(coordA[join.line.A, 1, 1, ], 1, is.na), arr.ind = T) [, 2]
    screwed.join.Z <- which(aaply(coordZ[join.line.Z, 1, 1, ], 1, is.na), arr.ind = T) [, 2]

    really.screwed <- c(screwed.join.A[!screwed.join.A %in% screwed.join.Z],
                        screwed.join.Z[!screwed.join.Z %in% screwed.join.A])

    id <- id[-really.screwed, ]
    lista.fino <- lista.fino [-really.screwed, ]
    coordA <- coordA[, , , -really.screwed]
    coordZ <- coordZ[, , , -really.screwed]
    
    ## glue    
    coord <-
        aaply(1:2, 1,
              function(i) aaply(1:nrow(id), 1,
                                function(j) GlueSkull (coordA [, , i, j], coordZ [, , i, j])))
              
    coord <- aperm(coord, c(3, 4, 1, 2))
    dimnames(coord)[3:4] <- list(c('R1', 'R2'), id $ numero)

    miss.lm <- aaply(aaply(coord[, , 1, ], c(1, 2), is.na), 3, any)
    
    info <- data.frame(
        'ID' = id $ numero,
        'GEN' = rep('Homo', times = nrow(id)),
        'SPE' = rep('sapiens', times = nrow(id)),
        'SEX' = lista.fino $ sex,
        'SUB' = lista.fino $ grupo,
        'MISS' = miss.lm,
        'FILE' = id [, 1])

    list('info' = info, 'coord' = coord)
    
}

