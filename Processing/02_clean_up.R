require(devtools)
require(roxygen2)
require(shapes)
require(plyr)

### PROSIMIAN
load('Prosimians/01_from_files.RData')

### arquivos pra filtragem

load('prosimian_extant.RData')
load('strepBase.RData')

### landmarks que importam (PARA UMA BASE DE DADOS DO LEM)
landmarks.that.matter <-
    c("IS", "NSL", "NA", "BR", "LD", "PNS", "BA", "OPI",
      "PT-E", "TSP-E", "FM-E", "ZS-E", "ZI-E", "ZYGO-E",
      "PM-E", "MT-E", "TS-E", "EAM-E", "PEAM-E", "AS-E", "APET-E", "JP-E",
      "IS", "NSL", "NA", "BR", "LD", "PNS", "BA", "OPI",
      "PT-D", "TSP-D", "FM-D", "ZS-D", "ZI-D", "ZYGO-D",
      "PM-D", "MT-D", "TS-D", "EAM-D" , "PEAM-D", "AS-D", "APET-D", "JP-D")

prosimian.id <- ldply(prosimian.raw, function(L) data.frame(L $ id), .id = 'file')

### isso pesca todo mundo, acho
all.missing.entries <- 
    is.na(prosimian.id $ Tombo) |
    prosimian.id $ Tombo == '' |
    prosimian.id $ Tombo == '0' |
    prosimian.id $ Ind == '00' |
    prosimian.id $ Ind == 'AMNH0'

prosimian.id <- prosimian.id [!all.missing.entries, ]

names(table(prosimian.miss [, 1])) [table(prosimian.miss [, 1]) > 9]

prosimian.id $ Ind <- gsub('ZMBZMB', 'ZMB', prosimian.id $ Ind)
prosimian.id $ Ind <- gsub('_Paris', '', prosimian.id $ Ind)
prosimian.id $ Museu <- gsub('_Paris', '', prosimian.id $ Museu)

nrow(prosimian.id [, 2:4])

prosimian.shapes <-
    array(0, c(length(landmarks.that.matter), 3, 2, 10 * length(prosimian.raw)))

for(i in 1:length(prosimian.raw))
{
    lms.to.keep <-
        which(dimnames(prosimian.raw [[i]] $ shapes) [[1]] %in% landmarks.that.matter)
    
    prosimian.shapes [, , , ((10*i)-9):(10*i)] <-
        prosimian.raw [[i]] $ shapes [lms.to.keep, , , ]
}
               
prosimian.shapes <- prosimian.shapes[, , , !all.missing.entries]

### ok, not ok (doubled entries)

ok.files <-
    prosimian.id $ file [grepl('_ok', prosimian.id $ file) &
                         (duplicated(prosimian.id $ Ind) |
                          duplicated(prosimian.id $ Ind, fromLast = TRUE))]

unique(ok.files)

not.ok.files <-
    prosimian.id $ file [(!grepl('_ok', prosimian.id $ file)) &
                         (duplicated(prosimian.id $ Ind) |
                          duplicated(prosimian.id $ Ind, fromLast = TRUE))]

unique(not.ok.files)

not.ok.files %in% gsub('_ok', '', ok.files)

which.not.ok <- which((!grepl('_ok', prosimian.id $ file)) &
                      (duplicated(prosimian.id $ Ind) |
                       duplicated(prosimian.id $ Ind, fromLast = TRUE)))

prosimian.id <- prosimian.id [- which.not.ok [not.ok.files %in%
                                              gsub('_ok', '', ok.files)], ]

prosimian.shapes <- prosimian.shapes [, , , - which.not.ok [not.ok.files %in%
                                                            gsub('_ok', '', ok.files)]]

prosimian.shapes <-
    prosimian.shapes[, , ,
                     - which(grepl('_ok', prosimian.id $ file) &
                             (duplicated(prosimian.id $ Ind) |
                              duplicated(prosimian.id $ Ind, fromLast = TRUE))) [3:9]]

prosimian.id <-
    prosimian.id [- which(grepl('_ok', prosimian.id $ file) &
                        (duplicated(prosimian.id $ Ind) |
                         duplicated(prosimian.id $ Ind, fromLast = TRUE))) [3:9], ]

prosimian.id [duplicated(prosimian.id $ Ind) |
              duplicated(prosimian.id $ Ind, fromLast = TRUE),
              c('file', 'Ind', 'Especie')]

still.double <-
    which(duplicated(prosimian.id $ Ind) |
          duplicated(prosimian.id $ Ind, fromLast = TRUE))

prosimian.id [prosimian.id $ Ind == 'RMNH28535', 'Especie'] [2] <- 'Tarsius bancanus'

still.double <- still.double [duplicated(prosimian.id [still.double, 'Ind']) &
                              duplicated(prosimian.id [still.double, 'Especie'])]

prosimian.id <- prosimian.id [- still.double, ]

prosimian.shapes <- prosimian.shapes [, , , - still.double]

prosimian.id [prosimian.id $ Ind == 'AMNH100589', 'Ind'] <- 
    paste0(prosimian.id [prosimian.id $ Ind == 'AMNH100589', 'Ind'],
           c('E', 'L'))

prosimian.id [prosimian.id $ Ind == 'MNHNMO-1910-101', 'Ind'] <- 
    paste0(prosimian.id [prosimian.id $ Ind == 'MNHNMO-1910-101', 'Ind'],
           c('E', 'L'))

nrow(prosimian.id)

dim(prosimian.shapes)

sp <- as.character(prosimian.id $ Especie)
sp [is.na(sp)] <- paste(prosimian.id $ Genero[is.na(sp)], 'sp.')

prosimian.id $ Especie <- factor(sp)


### consolidados da papete
strep.base <- unique(strep.Base)

strep.base <- gsub('ZMBZMB', 'ZMB', strep.base)
strep.base <- gsub('_Paris', '', strep.base)

prosimian.id [strep.base %in% prosimian.id $ Ind, ]

sum(prosimian.id $ Especie == 'Indri indri')




prosimian.cleanup <- list(prosimian.id, prosimian.shapes)

save(prosimian.cleanup, file = 'Prosimians/02_clean_up.RData')
