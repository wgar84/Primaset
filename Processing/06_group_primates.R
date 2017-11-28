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
require(ggplot2)
require(RColorBrewer)

##
load('Anthropoidea/05_from_grouping.RData')
load('Prosimians/02_clean_up.RData')

prosimians <- prosimian.cleanup
rm(prosimian.cleanup)

## no idea why still using unchecked names
dimnames(prosimians $ coord) [[4]] <- prosimians $ id $ Ind

primates <- GroupPrimates(prosimians, anthropoids)

save(primates, file = 'Primates/06_grouped.RData')


## cause yeah, why not?
prima.gpa <- procGPA(primates $ coord)

prima.ss.df <- data.frame(prima.gpa $ scores [, 1:2],
                          'logCS' = log(prima.gpa $ size),
                          'GEN' = primates $ info $ GEN,
                          'MAJOR' = primates $ info $ MAJOR)

grad <- colorRampPalette(brewer.pal(8, "Spectral"), space="Lab")

ggplot(prima.ss.df) +
    geom_text(aes(x = PC1, y = PC2, color = GEN, label = substr(GEN, 1, 3)),
              size = 3, alpha = 0.8) +
    scale_color_manual('Genus', values = grad(length(unique(prima.ss.df $ GEN)))) +
    scale_shape_discrete('Major Group') +
    theme_bw() +
    xlab('Shape PC1') +
    ylab('Shape PC2') +
    guides(color = guide_legend (ncol = 3))

ggsave('primate_shape.pdf')

ggplot(prima.ss.df) +
    geom_text(aes(y = PC1, x = logCS, color = GEN, label = substr(GEN, 1, 3)),
              size = 3, alpha = 0.8) +
    scale_color_manual('Genus', values = grad(length(unique(prima.ss.df $ GEN)))) +
    scale_shape_discrete('Major Group') +
    theme_bw() +
    xlab('Shape PC1') +
    ylab('Shape PC2') +
    guides(color = guide_legend (ncol = 3))

ggsave('primate_allo.pdf')

## that pesky weird cercopithecus

text3d(primates $ coord [, , 7255], texts = rownames(primates $ coord))
