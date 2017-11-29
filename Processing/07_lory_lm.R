require(devtools)
require(shapes)
require(geomorph)
require(rstan)
require(plyr)
require(RColorBrewer)

require(evolqg)
require(plotrix)

require(doMC)
registerDoMC(cores = 32)

## stan uses different backend now
rstan_options(auto_write = TRUE)
options(mc.cores = 32)

## type III sum of squares
options(contrasts = c('contr.sum', 'contr.poly'))

load('Primates/06_grouped.RData')
load('../Raw Data/Aux.RData')

## remove 6337, 7157 & 7255 and it'll be ok (some traits with negative determinants)

negadet <- as.character(primates $ info $ ID [c(6337, 7157, 7255)])

primates $ info <- subset(primates $ info, !ID %in% negadet)

primates $ coord <-
    primates $ coord [, , !dimnames(primates $ coord) [[3]] %in% negadet]

primates $ rep <-
    primates $ rep [, , , !dimnames(primates $ rep) [[4]] %in% negadet]

## have to estimate symmetric components first, idiot
## no duplicated IDs!

primates $ info $ ID <- as.character(primates $ info $ ID)

primates $ info [duplicated(primates $ info $ ID), 'ID'] <-
    paste(primates $ info [duplicated(primates $ info $ ID), 'ID'],
          LETTERS[1:16], sep = '_')

dimnames(primates $ coord) [[3]] <- primates $ info $ ID

dimnames(primates $ rep) [[4]] <- primates $ info $ ID [primates $ info $ rep]

primates $ info $ ID <- factor(primates $ info $ ID)

primates $ sym.coord <- aaply(primates $ coord, 3, Symmetrize, .parallel = TRUE)

primates $ sym.coord <- aperm(primates $ sym.coord, c(2, 3, 1))

dimnames(primates $ sym.coord) [3]

## LORY

primates $ sym.lory <- LORY(primates $ sym.coord, Aux $ single.tessel.38, TRUE)

primates $ info $ GSP <- paste(primates $ info $ GEN, primates $ info $ SPE, sep = '_')

## RENAME OBJECTS

prima.lory <- primates $ sym.lory

prima.lory $ info <- primates $ info

prima.info <- prima.lory $ info

prima.lory <- prima.lory[-6]
    
save(prima.lory, file = 'Primates/LORY.RData')

save(prima.info, file = 'Primates/Info.RData')

prima.raw <- primates [-6]

rm(primates)

save(prima.raw, file = 'Primates/Raw.RData')

## start tests on stan (saguinus)

sag.local <-
    prima.lory $ local [prima.lory $ info $ GSP == 'Saguinus_fuscicollis', ]

sag.pca <- prcomp(sag.local, retx = TRUE)

sag.df <-
    data.frame('sub' = subset(primates $ info, GSP == 'Saguinus_fuscicollis') $ SUB,
               'cs' = log(primates $ sym.lory $ cs [primates $ info $ GSP ==
                                                'Saguinus_fuscicollis']),
               sag.pca $ x [, 1:5])

sag.df $ sub <- factor(sag.df $ sub)

table(sag.df $ sub)


sag.plot <-
    ggplot(sag.df, aes(x = PC1, y = PC2, group = sub, color = sub)) +
    geom_point() +
    stat_ellipse(linetype = 1) +
    scale_color_brewer(type = 'div', palette = 5) +
    theme_bw()

ggsave(file = 'sag.pdf', plot = sag.plot, width = 12, height = 10)

## generate stan input

sag.modmat <- model.matrix(~ sub, data = sag.df)

sag.data <-
    cbind(log(primates $ sym.lory $ cs [primates $ info $ GSP ==
                                        'Saguinus_fuscicollis']),
          sag.local)

sag.input <- list('N' = nrow(sag.data),
                  'J' = ncol(sag.modmat),
                  'K' = ncol(sag.data),
                  'Y' = sag.data,
                  'X' = sag.modmat)
sag.init <-
    function(i, traits, effects)
        list('Omega_P' = chol(RandomMatrix(traits)),
             'sigma_P' = rchisq(traits, 1),
             'beta' = matrix(rnorm(traits * effects, 0, 1), nrow = traits))

sag.test <- 
    stan('../Stan/pmatrix_lm.stan', data = sag.input, thin = 10, pars = c('beta', 'P'),
         init = alply(1:4, 1, sag.init, traits = sag.input $ K, effects = sag.input $ J))

## model diagnostics

sag.post <- extract(sag.test)

sag.lp <- extract(sag.test, pars = 'lp__', permuted = FALSE)

pdf('sag_lp.pdf', width = 10, height = 10)

par(mfcol = c(2, 2))
for(i in 1:4) ## chains
    plot(sag.lp [, i, 1], type = 'l', main = paste('Chain', i),
         xlab = 'iterate', ylab = 'log prob')

dev.off(dev.cur())


pdf('sampledMat.pdf', width = 10, height = 10)
color2D.matplot(cov2cor(sag.post $ P [2, , ]))
dev.off(dev.cur())

sag.lm <- lm(sag.data ~ sag.df $ sub)

sag.mlmat <- CalculateMatrix(sag.lm)

pdf('maxlikMat.pdf', width = 10, height = 10)
color2D.matplot(cov2cor(sag.mlmat))
dev.off(dev.cur())

## RS
sag.rs.dist <-
    aaply(sag.post $ P, 1, RandomSkewers, cov.y = sag.mlmat, .parallel = TRUE)

pdf('histRS.pdf', width = 10, height = 10)
hist(sag.rs.dist [, 1])
dev.off(dev.cur())
