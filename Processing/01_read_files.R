require(devtools)
require(roxygen2)
require(shapes)
require(plyr)

## prosimian

prosimian.list <-
    dir(path = '../Raw Data/Prosimians/Collections',
        pattern = '.xlsx', recursive = TRUE, include.dirs = TRUE, full.names = TRUE)

prosimian.raw <-
    alply(prosimian.list, 1, function(f)
        {
            print(f)
            ReadProsimian(f)
        })

prosimian.list <-
    prosimian.list [!is.na(prosimian.raw)]

prosimian.raw <-
    prosimian.raw [!is.na(prosimian.raw)]

prosimian.list <- gsub('../Raw Data/Prosimians/Collections/', '', prosimian.list)

names(prosimian.raw) <- prosimian.list

save(prosimian.raw, file = 'Prosimians/01_from_files.RData')

## attach('Prosimians/01_from_files.RData')

### catarrhine

catarrhine.list <- 
    dir(path = '../Raw Data/Catarrhini',
        pattern = '.xls', recursive = TRUE, include.dirs = TRUE, full.names = TRUE)

catarrhine.raw <-
    alply(catarrhine.list, 1, function(f)
        {
            print(f)
            ReadCatarrhine(f)
        })

catarrhine.list <-
    catarrhine.list [!is.na(catarrhine.raw)]

catarrhine.raw <-
    catarrhine.raw [!is.na(catarrhine.raw)]

catarrhine.list <- gsub('../Raw Data/Catarrhini/', '', catarrhine.list)

names(catarrhine.raw) <- catarrhine.list

save(catarrhine.raw, file = 'Catarrhini/01_from_files.RData')

## attach('Catarrhini/01_from_files.RData')

### Homo

homo.list <- 
    dir(path = '../Raw Data/Homo/dataset- homo sapiens',
        pattern = '.xls', recursive = TRUE, include.dirs = TRUE, full.names = TRUE)

homo.raw <-
    alply(homo.list, 1, function(f)
        {
            print(f)
            ReadHomo(f)
        })

homo.list <-
    homo.list [!is.na(homo.raw)]

homo.raw <-
    homo.raw [!is.na(homo.raw)]

homo.list <- gsub('../Raw Data/Homo/dataset- homo sapiens/', '', homo.list)

names(homo.raw) <- homo.list

save(homo.raw, file = 'Homo/01_from_files.RData')


### calomys

calomys.list <-
    dir(path = '../Raw Data/Calomys', pattern = 'input',
        recursive = TRUE, include.dirs = TRUE, full.names = TRUE)

calomys.raw <-
    alply(calomys.list, 1, function(f)
        {
            print(f)
            ReadCalomys(f)
        })

calomys.list <-
    calomys.list [!is.na(calomys.raw)]

calomys.raw <-
    calomys.raw [!is.na(calomys.raw)]

calomys.list <- gsub('../Raw Data/Calomys/', '', calomys.list)

names(calomys.raw) <- calomys.list

save(calomys.raw, file = 'Calomys/01_from_files.RData')

### platyrrhini
platyrrhini.raw <- ReadPlatyrrhini('../Raw Data/Platyrrhini/')

save(platyrrhini.raw, file = 'Platyrrhini/01_from_files.RData')
