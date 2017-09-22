require(devtools)
require(roxygen2)
require(shapes)
require(plyr)
require(rgl)
require(geomorph)

load('Calomys/02_clean_up.RData')

calomys.complete <- MissLMCalomys(calomys.cleanup)

calomys.wgen <- GenealogyCalomys(calomys.complete)

str(calomys.wgen)

### now, for the real work, bring marquez to primaset
