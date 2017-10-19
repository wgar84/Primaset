require(devtools)
require(htmltools)
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
require(rmarkdown) 
require(knitr)

##
load('Anthropoidea/05_from_grouping.RData')
load('Prosimians/02_clean_up.RData')

prosimians <- prosimian.cleanup
rm(prosimian.cleanup)

dimnames(prosimians $ coord) [[4]] <- prosimians $ id $ Ind

## for half skulls, 

dimnames(prosimians $ coord) [[4]] [duplicated(dimnames(prosimians $ coord) [[4]])]



## for single replicates, we just discard one and 
