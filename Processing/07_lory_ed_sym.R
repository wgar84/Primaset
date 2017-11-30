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

prima.sym <- list()

prima.sym $ coord <- aaply(primates $ coord, 3, Symmetrize, .parallel = TRUE)

prima.sym $ coord <- aperm(prima.sym $ coord, c(2, 3, 1))

dimnames(prima.sym $ coord) [3]

## RENAME OBJECTS

prima.info <- prima.lory $ info

prima.info $ GSP <- paste(prima.info $ GEN, prima.info $ SPE, sep = '_')

save(prima.info, file = 'Primates/Info.RData')

prima.raw <- primates [-2]

save(prima.raw, file = 'Primates/Raw.RData')

save(prima.sym, file = 'Primates/Sym.RData')

rm(primates)

## tem que arrumar a ordem do single.tessel.38

prima.aux <- list()

prima.aux $ distance.matrix <-
    matrix(unlist(strsplit(rownames(Aux $ def.hyp) [-1], '\\.')), nrow = 38, byrow = TRUE)

prima.aux $ distance.numbers <- array(0, dim(prima.aux $ distance.matrix))

for(i in 1:38)
    for (j in 1:2)
        prima.aux $ distance.numbers[i, j] <-
            which(grepl(prima.aux $ distance.matrix [i, j],
                        rownames(prima.sym $ coord))) [1]

## ED

prima.ed <- list()

prima.ed $ raw <-
    aaply(prima.sym $ coord, 3, EuclideanDistances, dists = prima.aux $ distance.numbers,
          .parallel = TRUE)

colnames(prima.ed $ raw) <- aaply(prima.aux $ distance.matrix, 1, paste, collapse = '.')

rownames(prima.ed $ raw) <- prima.info $ ID

save(prima.ed, file = 'Primates/ed.RData')

## LORY

## first no local to see if there's any negative determinants
prima.lory <- LORY(prima.sym $ coord, prima.aux $ distance.numbers, TRUE)


## down this point, I did some tests that can be ignored later

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

## generate stan input (saguinus)

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

## generate stan input (homo)

hom.local <- prima.lory $ local [prima.info $ GSP == 'Homo_sapiens', ]

hom.df <-
    data.frame('sub' = subset(prima.info, GSP == 'Homo_sapiens') $ SUB,
               'sex' = subset(prima.info, GSP == 'Homo_sapiens') $ SEX,
               'cs' = log(prima.lory $ cs [prima.info $ GSP == 'Homo_sapiens']))
               
hom.df $ sub <- factor(hom.df $ sub)
hom.df $ sex <- factor(hom.df $ sex)

hom.modmat <- model.matrix(~ sub * sex, data = hom.df)

hom.data <- cbind(hom.df $ cs, hom.local)

hom.input <- list('N' = nrow(hom.data),
                  'J' = ncol(hom.modmat),
                  'K' = ncol(hom.data),
                  'Y' = hom.data,
                  'X' = hom.modmat)

initialConditions <-
    function(i, traits, effects)
        list('Omega_P' = chol(RandomMatrix(traits)),
             'sigma_P' = rchisq(traits, 1),
             'beta' = matrix(rnorm(traits * effects, 0, 1), nrow = traits))

hom.test <- 
    stan('../Stan/pmatrix_lm.stan', data = hom.input, thin = 10, pars = c('beta', 'P'),
         init = alply(1:4, 1, initialConditions,
                      traits = hom.input $ K, effects = hom.input $ J))

## model diagnostics

hom.post <- extract(hom.test)

hom.lp <- extract(hom.test, pars = 'lp__', permuted = FALSE)

pdf('hom_lp.pdf', width = 10, height = 10)

par(mfcol = c(2, 2))
for(i in 1:4) ## chains
    plot(hom.lp [, i, 1], type = 'l', main = paste('Chain', i),
         xlab = 'iterate', ylab = 'log prob')

dev.off(dev.cur())


pdf('sampledMat.pdf', width = 10, height = 10)
color2D.matplot(cov2cor(hom.post $ P [2, , ]))
dev.off(dev.cur())

hom.lm <- lm(hom.data ~ hom.df $ sub * hom.df $ sex)

hom.mlmat <- CalculateMatrix(hom.lm)

pdf('maxlikMat.pdf', width = 10, height = 10)
color2D.matplot(cov2cor(hom.mlmat))
dev.off(dev.cur())

## RS
hom.rs.dist <-
    aaply(hom.post $ P, 1, RandomSkewers, cov.y = hom.mlmat, .parallel = TRUE)

pdf('histRS.pdf', width = 10, height = 10)
hist(hom.rs.dist [, 1])
dev.off(dev.cur())

