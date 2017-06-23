##' @title
##' CleanUpProsimian
##' 
##' @description
##' This function wraps up our work on cleaning up the Prosimian data.
##' Repeated entries are removed, all individuals have the same number of landmarks,
##' that kind of stuff.
##' This is the second step towards a nice crispy Prosimian DB.
##'
##' @param xls.list output obtained from ReadProsimian, which basically just reads
##' the Excel files
##'
##' @return
##' a list with two elements, id and coord, the cleaned up version of what came in
##'
##' @author Guilherme Garcia, Anna Penna
##'
##' @importFrom plyr ldply
##'
##' @seealso ReadProsimian
##' 
##' @examples
##' \dontrun{
##' prosimian.cleanup <- CleanUpProsimian(prosimian.raw)
##' }

CleanUpProsimian <- function(xls.list)
{
    ## landmarks que importam (PARA UMA BASE DE DADOS DO LEM)
    landmarks.that.matter <-
        c("IS", "NSL", "NA", "BR", "LD", "PNS", "BA", "OPI",
          "PT-E", "TSP-E", "FM-E", "ZS-E", "ZI-E", "ZYGO-E",
          "PM-E", "MT-E", "TS-E", "EAM-E", "PEAM-E", "AS-E", "APET-E", "JP-E",
          "IS", "NSL", "NA", "BR", "LD", "PNS", "BA", "OPI",
          "PT-D", "TSP-D", "FM-D", "ZS-D", "ZI-D", "ZYGO-D",
          "PM-D", "MT-D", "TS-D", "EAM-D" , "PEAM-D", "AS-D", "APET-D", "JP-D")

    id <- ldply(xls.list, function(L) data.frame(L $ id), .id = 'file')

    ## isso pesca todo mundo, acho
    all.missing.entries <-
        is.na(id $ Tombo) |
        id $ Tombo == '' |
        id $ Tombo == '0' |
        id $ Ind == '00' |
        id $ Ind == 'AMNH0'

    id <- id [!all.missing.entries, ]

    id $ Ind <- gsub('ZMBZMB', 'ZMB', id $ Ind)
    id $ Ind <- gsub('_Paris', '', id $ Ind)
    id $ Museu <- gsub('_Paris', '', id $ Museu)

    shapes <- array(0, c(length(landmarks.that.matter), 3, 2, 10 * length(xls.list)))

    for(i in 1:length(xls.list))
    {
        lms.to.keep <-
            which(dimnames(xls.list [[i]] $ shapes) [[1]] %in% landmarks.that.matter)
        
        shapes [, , , ((10*i)-9):(10*i)] <-
            xls.list [[i]] $ shapes [lms.to.keep, , , ]
    }

    shapes <- shapes[, , , !all.missing.entries]

    ## ok, not ok (doubled entries)
    ok.files <-
        id $ file [grepl('_ok', id $ file) &
                   (duplicated(id $ Ind) |
                    duplicated(id $ Ind, fromLast = TRUE))]

    not.ok.files <-
        id $ file [(!grepl('_ok', id $ file)) &
                   (duplicated(id $ Ind) |
                    duplicated(id $ Ind, fromLast = TRUE))]
    
    which.not.ok <- which((!grepl('_ok', id $ file)) &
                          (duplicated(id $ Ind) |
                           duplicated(id $ Ind, fromLast = TRUE)))

    id <- id [- which.not.ok [not.ok.files %in%
                              gsub('_ok', '', ok.files)], ]

    shapes <- shapes [, , , - which.not.ok [not.ok.files %in%
                                            gsub('_ok', '', ok.files)]]
    
    shapes <-
        shapes[, , ,
               - which(grepl('_ok', id $ file) &
                       (duplicated(id $ Ind) |
                        duplicated(id $ Ind, fromLast = TRUE))) [3:9]]

    id <-
        id [- which(grepl('_ok', id $ file) &
                    (duplicated(id $ Ind) |
                     duplicated(id $ Ind, fromLast = TRUE))) [3:9], ]

    id [duplicated(id $ Ind) |
        duplicated(id $ Ind, fromLast = TRUE),
        c('file', 'Ind', 'Especie')]
    
    still.double <-
        which(duplicated(id $ Ind) |
              duplicated(id $ Ind, fromLast = TRUE))
    
    id [id $ Ind == 'RMNH28535', 'Especie'] [2] <- 'Tarsius bancanus'
    
    still.double <- still.double [duplicated(id [still.double, 'Ind']) &
                                  duplicated(id [still.double, 'Especie'])]
    
    id <- id [- still.double, ]
    
    shapes <- shapes [, , , - still.double]
    
    id [id $ Ind == 'AMNH100589', 'Ind'] <-
        paste0(id [id $ Ind == 'AMNH100589', 'Ind'],
               c('E', 'L'))
    
    id [id $ Ind == 'MNHNMO-1910-101', 'Ind'] <-
        paste0(id [id $ Ind == 'MNHNMO-1910-101', 'Ind'],
               c('E', 'L'))

    sp <- as.character(id $ Especie)
    sp [is.na(sp)] <- paste(id $ Genero[is.na(sp)], 'sp.')

    id $ Especie <- factor(sp)

    dimnames(shapes) <- list(landmarks.that.matter,
                             LETTERS[24:26],
                             paste0("R", 1:2),
                             paste(id$Museu, id$Tombo, sep = "_"))
    
    cleanup <- list('id' = id, 'coord' = shapes)
    
    cleanup   
}
