require(devtools)
require(shapes)
require(geomorph)
require(rstan)
require(plyr)
require(RColorBrewer)
require(evolqg)
require(plotrix)
require(doMC)
require(dplyr)
require(magrittr)
require(tidyr)
require(geiger)
require(ape)
require(ggplot2)

load('Primates/Info.RData')
load('Primates/LORY.RData')
load('Primates/newaux.RData')


## le screwed taxonomy, final part

## 10K TREES!!!!!

prima.info $ ARN <- prima.info $ GSP

## Alouatta_villosa -> Alouatta_pigra

prima.info $ ARN [which(prima.info $ GSP == 'Alouatta_villosa')] <- 'Alouatta_pigra'

## Aotus_nancymai -> Aotus_nancymaae

prima.info $ ARN [which(prima.info $ ARN == 'Aotus_nancymai')] <- 'Aotus_nancymaae'

## Ateles_belzebuth <- chamek e marginatus

prima.info $ ARN [which(prima.info $ ARN == 'Ateles_chamek')] <- 'Ateles_belzebuth'
prima.info $ ARN [which(prima.info $ ARN == 'Ateles_marginatus')] <- 'Ateles_belzebuth'

data.frame('sp' = prima.info $ SPE [prima.info $ GEN == 'Ateles'],
            'logCS' = log(prima.lory $ cs [prima.info $ GEN == 'Ateles']),
            prcomp(
                prima.lory $ local [prima.info $ GEN == 'Ateles', ]) $ x [, 1:2]) %>%
    ggplot(., aes(x = PC1, y = PC2, color = sp, group = sp)) +
    geom_point() +
    stat_ellipse()
    

## Hoolock_hoolock -> Bunopithecus_hoolock

prima.info $ ARN [which(
                 prima.info $ ARN == 'Hoolock_hoolock')] <- 'Bunopithecus_hoolock'

## Cebus

## data.frame('sp' = prima.info $ SPE [prima.info $ GEN == 'Cebus'],
##            'logCS' = log(prima.lory $ cs [prima.info $ GEN == 'Cebus']),
##            prcomp(
##                prima.lory $ local [prima.info $ GEN == 'Cebus', ]) $ x [, 1:2]) %>%
##     ggplot(.) +
##     geom_point(aes(x = PC1, y = PC2, color = sp))

prima.info $ ARN [which(
                 prima.info $ ARN == 'Cebus_nigritus')] <- 'Cebus_apella'

prima.info $ ARN [which(
                 prima.info $ ARN == 'Cebus_paraguayanus')] <- 'Cebus_apella'

## kuhliiiiiiii

prima.info $ ARN [which(
                 prima.info $ ARN == 'Callithrix_kuhlii')] <- 'Callithrix_kuhli'

## l8r

## Springer

ir.hb = read.tree ('../Springer2012/ir.hb.tre')
ir.sb = read.tree ('../Springer2012/ir.sb.tre')
ac.sb = read.tree ('../Springer2012/ac.sb.tre')
ac.hb = read.tree ('../Springer2012/ac.hb.tre')

Springer2012 = list ('ac.hb' = ac.hb, 'ac.sb' = ac.sb, 'ir.hb' = ir.hb, 'ir.sb' = ir.sb)

save (Springer2012, file = 'Primates/springer2012.RData')

prima.info $ SPR <- prima.info $ GSP

## Alouatta_villosa -> Alouatta_pigra

prima.info $ SPR [which(prima.info $ GSP == 'Alouatta_villosa')] <- 'Alouatta_pigra'

## Aotus_nancymai -> Aotus_nancymaae

prima.info $ SPR [which(prima.info $ SPR == 'Aotus_nancymai')] <- 'Aotus_nancymaae'

## Ateles_belzebuth <- chamek e marginatus

prima.info $ SPR [which(prima.info $ SPR == 'Ateles_chamek')] <- 'Ateles_belzebuth'
prima.info $ SPR [which(prima.info $ SPR == 'Ateles_marginatus')] <- 'Ateles_belzebuth'

## data.frame('sp' = prima.info $ SPR [prima.info $ GEN == 'Ateles'],
##             'logCS' = log(prima.lory $ cs [prima.info $ GEN == 'Ateles']),
##             prcomp(
##                 prima.lory $ local [prima.info $ GEN == 'Ateles', ]) $ x [, 1:2]) %>%
##     ggplot(., aes(x = PC1, y = PC2, color = sp, group = sp)) +
##     geom_point() +
##     stat_ellipse()

## Cebus

## data.frame('sp' = prima.info $ SPE [prima.info $ ARN == 'Cebus_apella'],
##            'logCS' = log(prima.lory $ cs [prima.info $ ARN == 'Cebus_apella']),
##            prcomp(
##                prima.lory $ local [prima.info $ ARN == 'Cebus_apella', ]) $
##            x [, 1:2]) %>%
##     ggplot(., aes(x = PC1, y = PC2, color = sp)) +
##     geom_point() +
##     stat_ellipse()

prima.info $ SPR [which(
                 prima.info $ SPR == 'Cebus_nigritus')] <- 'Cebus_apella'

prima.info $ SPR [which(
                 prima.info $ SPR == 'Cebus_paraguayanus')] <- 'Cebus_apella'


gsp.list <- table(prima.info $ SPR)

allo.tree <- treedata(Springer2012 [[1]], gsp.list)

pdf ('achb.allo.pdf', width = 8, height = 30)
plot (allo.tree [[1]], cex = 0.5, show.node.label = TRUE, direction = 'rightwards')
nodelabels ()
dev.off (dev.cur ())

table(prima.info $ SPR [!(prima.info $ SPR %in% rownames(allo.tree $ data))])

table(prima.info $ SPR [prima.info $ SPR %in% rownames(allo.tree $ data)])

save(allo.tree, file = 'Primates/allo.tree.RData')

### Loris outliers

data.frame('sp' = prima.info $ SPE [prima.info $ GEN == 'Loris'],
           'sex' = prima.info $ SEX [prima.info $ GEN == 'Loris'],
           'logCS' = log(prima.lory $ cs [prima.info $ GEN == 'Loris']),
           prcomp(
               prima.lory $ local [prima.info $ GEN == 'Loris', ]) $ x [, 1:3]) %>%
    ggplot(.) +
    geom_point(aes(x = PC1, y = PC3, color = sp, shape = sex))
