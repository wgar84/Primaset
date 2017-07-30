require(devtools)
require(roxygen2)
require(shapes)
require(plyr)
require(rgl)

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

load('Calomys/01_from_files.RData')

load('../Raw Data/Calomys/skull.RData')
load('../Raw Data/Calomys/skull2.RData')

calomys.id <- rbind(exp [[1]] [, 1:7], exp2 [[1]] [, 1:7])
rm(exp, exp2)

calomys.cleanup <- CleanUpCalomys(xls.list = calomys.raw, id = calomys.id)

save(calomys.cleanup, file = 'Calomys/02_clean_up.RData')

## calomys.gpa <-
##     procGPA(calomys.cleanup $ coord [, , 1, !calomys.cleanup $ info $ MISS])

## for(i in 1:(dim(calomys.gpa $ rotated)[3]))
##     rgl::points3d(calomys.gpa $ rotated [, , i])

### HOMO

load('Homo/01_from_files.RData')

homo.cleanup <- CleanUpHomo(homo.raw)

save(homo.cleanup, file = 'Homo/02_clean_up.RData')

### CATARRHINI

load('../Raw Data/Catarrhini/raw.RData') ### raw from previous iteration
catarrhini.raw <- owm.raw
rm(owm.raw)

fino.list <- read.csv('../Raw Data/Catarrhini/catarrhini.csv')
fino.list <- subset(fino.list, GENUS != 'Homo')

catarrhini.cleanup <- CleanUpCatarrhini(catarrhini.raw, fino.list)

save(catarrhini.cleanup, file = 'Catarrhini/02_clean_up.RData')

### PLATYRRHINI

## load('Platyrrhini/01_from_files.RData')
platyrrhini.cleanup <- CleanUpPlatyrrhini(platyrrhini.raw)

save(platyrrhini.cleanup, file = 'Platyrrhini/02_clean_up.RData')
