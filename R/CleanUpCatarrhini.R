##' @title
##' CleanUpCatarrhini
##' 
##' @description
##' This function cleans up Catarrhini data from excel files.
##' 
##' @param raw.list output from previous iteration of reading this DB
##' @param or.list consolidated DB from Oliveira et al (2009)
##' 
##' @return a list with two elements, info and coord
##'
##' @details will GlueSkulls
##' 
##' @author Guilherme Garcia
##'
##' @importFrom plyr llply aaply ldply
##'
##' @seealso 

CleanUpCatarrhini <- function(owm.raw, or.list)
{
    ### limpando raw (código do eu do passado)
    owm.raw$lmA = c("IS","PM-E","NSL","NA","BR","PT-E","FM-E",
                    "ZS-E","ZI-E","MT-E",
                    "PNS","APET-E","BA","OPI","EAM-E","PEAM-E",
                    "ZYGO-E","TSP-E","TS-E","JP-E",
                    "PM-D","PT-D","FM-D","ZS-D","ZI-D","MT-D","APET-D","EAM-D",
                    "PEAM-D","ZYGO-D",
                    "TSP-D","TS-D","JP-D")
    
    test = owm.raw$Z1[,1]
    dim (test) = c(3,12)
    test = t(test)
    rownames (test) = owm.raw$lmZ
    
    owm.raw$lmZ = c("LD","AS-D","JP-D","TS-D","PT-D","OPI",
                    "BA","BR","AS-E","JP-E","TS-E","PT-E")
    
    dim (owm.raw$A1) = c(3,length (owm.raw$lmA), ncol (owm.raw$A1))
    owm.raw$A1 = aperm (owm.raw$A1, c(2,1,3), resize = TRUE)
    
    dim (owm.raw$A2) = c(3,length (owm.raw$lmA), ncol (owm.raw$A2))
    owm.raw$A2 = aperm (owm.raw$A2, c(2,1,3), resize = TRUE)
    
    dim (owm.raw$Z1) = c(3,length (owm.raw$lmZ), ncol (owm.raw$Z1))
    owm.raw$Z1 = aperm (owm.raw$Z1, c(2,1,3), resize = TRUE)
    
    dim (owm.raw$Z2) = c(3,length (owm.raw$lmZ), ncol (owm.raw$Z2))
    owm.raw$Z2 = aperm (owm.raw$Z2, c(2,1,3), resize = TRUE)
    
    ## removendo vistas A sem vistas Z e vice-versa
    
    dataA = paste (owm.raw$dataA[,1], owm.raw$dataA[,5], sep = '')
    dataZ = paste (owm.raw$dataZ[,1], owm.raw$dataZ[,5], sep = '')
    
    a.in.z = dataA %in% dataZ
    z.in.a = dataZ %in% dataA
    
    owm.raw$A1 = owm.raw$A1[,,a.in.z]
    owm.raw$Z1 = owm.raw$Z1[,,z.in.a]
    owm.raw$A2 = owm.raw$A2[,,a.in.z]
    owm.raw$Z2 = owm.raw$Z2[,,z.in.a]
    
    owm.raw$dataA = owm.raw$dataA[a.in.z,]
    owm.raw$dataZ = owm.raw$dataZ[z.in.a,]
    
    dataA = paste (owm.raw$dataA[,1], owm.raw$dataA[,5], sep = '')
    dataZ = paste (owm.raw$dataZ[,1], owm.raw$dataZ[,5], sep = '')
    
    ## ordenando indivíduos
    
    order.A = order (dataA)
    order.Z = order (dataZ)
    
    owm.raw$A1 = owm.raw$A1[,,order.A]
    owm.raw$Z1 = owm.raw$Z1[,,order.Z]
    owm.raw$A2 = owm.raw$A2[,,order.A]
    owm.raw$Z2 = owm.raw$Z2[,,order.Z]
    
    owm.raw$dataA = owm.raw$dataA[order.A,]
    owm.raw$dataZ = owm.raw$dataZ[order.Z,]
    
    dataA = dataA[order.A]
    dataZ = dataZ[order.Z]

    ### dados do fino
    data.lm <- or.list

    specimen.lm = paste (data.lm$MUSEUM, data.lm$NUMBER, sep = '')
    
    sp.ord = order (specimen.lm)
    
    data.lm = data.lm [sp.ord,]
    specimen.lm = specimen.lm [sp.ord]
    
    sum (duplicated (specimen.lm))
    
    specimen.lm.plus = paste (data.lm$MUSEUM, data.lm$NUMBER,
                              data.lm$GENUS, data.lm$SPECIES, sep = '_')
   
    sp.dupl = !duplicated (specimen.lm.plus)
    
    data.lm = data.lm [sp.dupl,]
    specimen.lm = specimen.lm [sp.dupl]
    specimen.lm.plus = specimen.lm.plus [sp.dupl]

    dataA.plus = paste (owm.raw$dataA[,1], owm.raw$dataA[,5],
                        owm.raw$dataA[,2], owm.raw$dataA[,3], sep = '_')

    dataA.dupl <- !duplicated(dataA.plus)
    
    owm.raw$A1 = owm.raw$A1[, , dataA.dupl]
    owm.raw$Z1 = owm.raw$Z1[, , dataA.dupl]
    owm.raw$A2 = owm.raw$A2[, , dataA.dupl]
    owm.raw$Z2 = owm.raw$Z2[, , dataA.dupl]
    
    owm.raw$dataA = owm.raw$dataA[dataA.dupl, ]
    owm.raw$dataZ = owm.raw$dataZ[dataA.dupl,]
    
    dataA = dataA[dataA.dupl]
    dataZ = dataZ[dataA.dupl]

    dataA.plus <- dataA.plus[dataA.dupl]
    
    length(dataA)
    nrow(data.lm)

    no.sp <- owm.raw $ dataA [, 3] == 'sp' | is.na(owm.raw $ dataA [, 3])

    owm.raw$A1 = owm.raw$A1[, , !no.sp]
    owm.raw$Z1 = owm.raw$Z1[, , !no.sp]
    owm.raw$A2 = owm.raw$A2[, , !no.sp]
    owm.raw$Z2 = owm.raw$Z2[, , !no.sp]
    
    owm.raw$dataA = owm.raw$dataA[!no.sp, ]
    owm.raw$dataZ = owm.raw$dataZ[!no.sp,]
    
    dataA = dataA[!no.sp]
    dataZ = dataZ[!no.sp]

    dataA.plus <- dataA.plus [!no.sp]

    no.sp.lm <- as.character(data.lm $ SPECIES) == ''

    data.lm <- data.lm[!no.sp.lm, ]
    specimen.lm <- specimen.lm [!no.sp.lm]
    specimen.lm.plus <- specimen.lm.plus [!no.sp.lm]

    nrow(data.lm [duplicated(specimen.lm) | duplicated(specimen.lm, fromLast = TRUE),
                  c('MUSEUM', 'NUMBER', 'GENUS', 'SPECIES')])
 
    ###

    order.lm <- order(data.lm $ MUSEUM, data.lm $ NUMBER,
                      data.lm $ GENUS, data.lm $ SPECIES)
    
    data.lm <- data.lm [order.lm, ]
    specimen.lm <- specimen.lm [order.lm]
    specimen.lm.plus <- specimen.lm.plus [order.lm]

    order.data <- order(owm.raw $ dataA [, 1], owm.raw $ dataA [, 5],
                        owm.raw $ dataA [, 2], owm.raw $ dataA [, 3])

    owm.raw $ A1 <- owm.raw $ A1 [, , order.data]
    owm.raw $ A2 <- owm.raw $ A2 [, , order.data]
    owm.raw $ Z1 <- owm.raw $ Z1 [, , order.data]
    owm.raw $ Z2 <- owm.raw $ Z2 [, , order.data]
    
    owm.raw $ dataA <- owm.raw $ dataA [order.data, ]
    owm.raw $ dataZ <- owm.raw $ dataZ [order.data, ]

    dataA <- dataA [order.data]
    dataA.plus <- dataA.plus [order.data]
    
    specimen.lm.raw = paste (data.lm$MUSEUM, data.lm$NUMBER,
                              data.lm$GENUS1, data.lm$SPECIES1, sep = '_')
   
    no.lm.raw <- dataA.plus %in% specimen.lm.raw

    owm.raw $ A1 <- owm.raw $ A1 [, , no.lm.raw]
    owm.raw $ A2 <- owm.raw $ A2 [, , no.lm.raw]
    owm.raw $ Z1 <- owm.raw $ Z1 [, , no.lm.raw]
    owm.raw $ Z2 <- owm.raw $ Z2 [, , no.lm.raw]
    
    owm.raw $ dataA <- owm.raw $ dataA [no.lm.raw, ]
    owm.raw $ dataZ <- owm.raw $ dataZ [no.lm.raw, ]

    dataA <- dataA [no.lm.raw]
    dataA.plus <- dataA.plus [no.lm.raw]

    no.dataA <-  specimen.lm.raw %in% dataA.plus

    data.lm <- data.lm [no.dataA, ]
    specimen.lm <- specimen.lm [no.dataA]
    specimen.lm.plus <- specimen.lm.plus [no.dataA]
    specimen.lm.raw <- specimen.lm.raw [no.dataA]
    
    all(specimen.lm.raw == dataA.plus)

    rownames(owm.raw $ A1) <- rownames(owm.raw $ A2) <- owm.raw $ lmA
    rownames(owm.raw $ Z1) <- rownames(owm.raw $ Z2) <- owm.raw $ lmZ
    
    coord <-
        aaply(1:length(dataA.plus), 1,
              function(i)
              {
                  aaply(1:2, 1,
                        function(j)
                        {
                            Aij <- owm.raw [[j]] [, , i]
                            Zij <- owm.raw [[j + 2]] [, , i]
                            GlueSkull(Aij, Zij)
                        })})

    coord <- aperm(coord, c(3, 4, 2, 1))
    
    dimnames(coord)[2:4] <- list(c('X', 'Y', 'Z'), c('R1', 'R2'), specimen.lm.plus)
    miss.lm <- aaply(aaply(coord[, , 1, ], c(1, 2), is.na), 3, any)
    
    info <-
        data.frame('ID' = specimen.lm.plus,
                   'GEN' = data.lm $ GENUS,
                   'SPE' = data.lm $ SPECIES,
                   'SEX' = data.lm $ SEX,
                   'SUB' = data.lm $ SUBSPECIES,
                   'MISS' = miss.lm,
                   'FILE' = data.lm $ FILE)
    
    list('info' = info, 'coord' = coord)
    
}
