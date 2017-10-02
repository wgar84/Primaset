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

## consolidados da tese
load('../Raw Data/def.info.RData')


## batendo as planilhas

platyrrhine <- platyrrhini.cleanup
catarrhine <- catarrhini.cleanup
homo <- homo.cleanup
prev <- def.info

museum.ids <-
    c(paste(plat.info $ MUSEUM, plat.info $ ID, sep = '_'), 
      laply(strsplit(as.character(cata.info $ ID), '_'),
            function(L) paste(L[1], L[2], sep = '_')),
      paste(gsub("[[:digit:]]", "", homo.info $ ID),
            homo.info $ ID, sep = '_'))

prev.info [!paste(prev.info $ MSM, prev.info $ ID, sep = '_') %in% museum.ids, ]
