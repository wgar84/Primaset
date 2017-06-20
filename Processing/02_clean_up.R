require(devtools)
require(roxygen2)
require(shapes)
require(plyr)

### PROSIMIAN
load('Prosimians/01_from_files.RData')

### arquivos pra filtragem

load('prosimian_extant.RData')
load('strepBase.RData')

## source('clean_up_backup.R')

## essa função contem o trampo de limpeza da semana anterior
prosimian.cleanup <- CleanUpProsimian(prosimian.raw)

### consolidados da papete
strep.base <- unique(strep.Base)

strep.base <- gsub('ZMBZMB', 'ZMB', strep.base)
strep.base <- gsub('_Paris', '', strep.base)

## dimnames(prosimian.shapes) <-
##     list(landmarks.that.matter,
##          LETTERS[24:26],
##          paste0("take", 1:2),
##          paste(prosimian.id$Museu, prosimian.id$Tombo, sep = "_"))

save(prosimian.cleanup, file = 'Prosimians/02_clean_up.RData')

### PLATYRRHINI

load('Platyrrhini/01_from_files.RData')
