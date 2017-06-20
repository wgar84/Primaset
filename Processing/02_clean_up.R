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

save(prosimian.cleanup, file = 'Prosimians/02_clean_up.RData')

### consolidados da papete
strep.base <- unique(strep.Base)

strep.base <- gsub('ZMBZMB', 'ZMB', strep.base)
strep.base <- gsub('_Paris', '', strep.base)


### CALOMYS

load('../Raw Data/Calomys/skull.RData')
load('../Raw Data/Calomys/skull2.RData')

calomys.id <- rbind(exp [[1]] [, 1:7], exp2 [[1]] [, 1:7])
rm(exp, exp2)

calomys.cleanup <- CleanUpCalomys(xls.list = calomys.raw, id = calomys.id)

save(calomys.cleanup, file = 'Calomys/02_clean_up.RData')

### PLATYRRHINI

## load('Platyrrhini/01_from_files.RData')

