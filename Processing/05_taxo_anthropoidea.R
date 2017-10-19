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
load('Platyrrhini/02_clean_up.RData')
load('Catarrhini/02_clean_up.RData')
load('Homo/02_clean_up.RData')

## this takes a while to run
anthropoids <- GroupAnthropoids(platyrrhini.cleanup,
                                catarrhini.cleanup,
                                homo.cleanup)

save(anthropoids, file = 'Anthropoidea/05_from_grouping.RData')

