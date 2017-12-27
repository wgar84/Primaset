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

registerDoMC(cores = 4)

## stan uses different backend now
rstan_options(auto_write = TRUE)
options(mc.cores = 32)

load('Primates/Info.RData')
load('Primates/Sym.RData')
load('Primates/aux.RData')
load('Primates/ed.RData')
load('Primates/LORY.RData')
load('../Raw Data/Aux.RData')

devtools::load_all('/home/guilherme/lem_essentials/Primaset')

## start

prima.info [4475, 'SEX'] <- 'M' ## whatever, doesn't matter

save(prima.info, file = 'Primates/Info.RData')

## select sps for model fitting

sp2fit <- names(table(prima.info $ GSP)) [table(prima.info $ GSP) > 30]

sp2fit <- c(sp2fit, rownames(Aux $ data.man.sp) [!rownames(Aux $ data.man.sp) %in% sp2fit])

sp2fit [!sp2fit %in% rownames(Aux $ data.man.sp)]

## see if we have new anthropoid sps, find best models to fit
## add Prosimians to data.man.sp
## each sp inspected thru pca

add.data.man <- 
    matrix(byrow = TRUE, ncol = 3,
           data = c(
               'Avahi_laniger', 'NONE', 'NONE', 
               'Cebus_nigritus', 'SEX', 'SEX',
               'Cercopithecus_denti', 'SEX', 'SEX',
               'Cheirogaleus_major', 'NONE', 'NONE',
               'Cheirogaleus_medius', 'NONE', 'NONE',
               'Chiropotes_albinasus', 'SEX', 'SEX',
               'Daubentonia_madagascariensis', 'NONE', 'NONE',
               'Eulemur_albifrons', 'NONE', 'NONE',
               'Eulemur_collaris', 'NONE', 'NONE',
               'Eulemur_fulvus', 'NONE', 'NONE',
               'Eulemur_macaco', 'NONE', 'NONE',
               'Eulemur_mongoz', 'NONE', 'NONE',
               'Eulemur_rubriventer', 'NONE', 'NONE',
               'Eulemur_rufus', 'NONE', 'NONE',
               'Galago_senegalensis', 'NONE', 'NONE',
               'Hapalemur_griseus', 'NONE', 'NONE',
               'Indri_indri', 'NONE', 'NONE',
               'Lemur_catta', 'NONE', 'NONE',
               'Lepilemur_leucopus', 'NONE', 'NONE',
               'Lepilemur_ruficaudatus', 'NONE', 'NONE',
               'Macaca_assamensis', 'SEX', 'SEX',
               'Microcebus_griseorufus', 'NONE', 'NONE',
               'Microcebus_murinus', 'NONE', 'NONE',
               'Nycticebus_coucang', 'NONE', 'NONE',
               'Otolemur_crassicaudatus', 'NONE', 'NONE',
               'Perodicticus_potto', 'NONE', 'NONE',
               'Presbytis_melalophos', 'NONE', 'NONE', ## change presbytis genus to none
               'Propithecus_diadema', 'NONE', 'NONE',
               'Propithecus_verreauxi', 'NONE', 'NONE',
               'Saguinus_geoffroyi', 'NONE', 'NONE',
               'Saguinus_oedipus', 'NONE', 'NONE',
               'Saimiri_cassiquiarensis', 'NONE', 'NONE', ## change saimiri genus to none
               'Varecia_variegata', 'NONE', 'NONE'))

prima.aux $ models <- data.frame('GSP' = rownames(Aux $ data.man.sp),
                                 'MOD' = Aux $ data.man.sp [, 2], row.names = NULL)

prima.aux $ models[grepl('Saimiri', prima.aux $ models $ GSP), 'MOD'] <- 'NONE'
prima.aux $ models[grepl('Presbytis', prima.aux $ models $ GSP), 'MOD'] <- 'NONE'

add.df <- data.frame('GSP' = add.data.man [, 1],
                     'MOD' = add.data.man [, 3])

prima.aux $ models <- rbind.fill(prima.aux $ models, add.df)

prima.aux $ models [grepl('Trachypithecus',
                          prima.aux $ models $ GSP), 'MOD'] [3:4] <- 'SEX'

save(prima.aux, file = 'Primates/aux.RData')

## example code for inspection

cur.sp <- 'Trachypithecus_phayrei'

cur.data <-
    cbind(log(prima.lory $ cs), prima.lory $ local) [prima.info $ GSP == cur.sp, ]

cur.pca <- prcomp(cur.data, retx = TRUE)
cur.fix <- prima.info [prima.info $ GSP == cur.sp, c('SEX', 'SUB')]
cur.df <- data.frame(cur.fix, cur.pca $ x [, 1:3])

pdf('cur.pdf', width = 10, height = 7)

ggplot(cur.df, aes(x = PC1, y = PC2, color = SUB, group = SEX, shape = SEX)) +
    geom_point() +
    stat_ellipse() +
    theme_bw()

cur.data <- prima.ed $ raw [prima.info $ GSP == cur.sp, ]
cur.pca <- prcomp(cur.data, retx = TRUE)
cur.fix <- prima.info [prima.info $ GSP == cur.sp, c('SEX', 'SUB')]
cur.df <- data.frame(cur.fix, cur.pca $ x [, 1:3])

ggplot(cur.df, aes(x = PC1, y = PC2, color = SUB, group = SEX, shape = SEX)) +
    geom_point() +
    stat_ellipse() +
    theme_bw()

dev.off(dev.cur())

## now, some models

registerDoMC(cores = 53)

none.Pmat <- 
    aaply(as.character(prima.aux $ models $ GSP [prima.aux $ models $ MOD == 'NONE']), 1,
          function(cur.sp)
          {
              cur.data <- cbind(log(prima.lory $ cs),
                                prima.lory $ local) [prima.info $ GSP == cur.sp, ]
              
              cur.mod <- prima.aux $ models $ MOD [prima.aux $ models $ GSP == cur.sp]
              
              cur.info <- subset(prima.info, GSP == cur.sp)
              
              cur.Pmat <-
                  PmatrixPosterior(data = cur.data,
                                   model = cur.mod,
                                   info = cur.info,
                                   full.output = FALSE,
                                   thin = 10)
              
              cur.Pmat
          }, .parallel = TRUE)

save.image()

dimnames(none.Pmat) [[1]] <-
    as.character(prima.aux $ models $ GSP [prima.aux $ models $ MOD == 'NONE'])

registerDoMC(cores = 78)

sex.Pmat <- 
    aaply(as.character(prima.aux $ models $ GSP [prima.aux $ models $ MOD == 'SEX']), 1,
          function(cur.sp)
          {
              cur.data <- cbind(log(prima.lory $ cs),
                                prima.lory $ local) [prima.info $ GSP == cur.sp, ]
              
              cur.mod <- prima.aux $ models $ MOD [prima.aux $ models $ GSP == cur.sp]
              
              cur.info <- subset(prima.info, GSP == cur.sp)
              
              cur.Pmat <-
                  PmatrixPosterior(data = cur.data,
                                   model = cur.mod,
                                   info = cur.info,
                                   full.output = FALSE,
                                   thin = 10)
              
              cur.Pmat
          }, .parallel = TRUE); save.image()


dimnames(sex.Pmat) [[1]] <-
    as.character(prima.aux $ models $ GSP [prima.aux $ models $ MOD == 'SEX'])

save.image()

registerDoMC(cores = 11)

other.Pmat <-
    aaply(as.character(prima.aux $ models $ GSP [prima.aux $ models $ MOD != 'SEX' &
                                                 prima.aux $ models $ MOD != 'NONE']), 1,
          function(cur.sp)
          {
              cur.data <- cbind(log(prima.lory $ cs),
                                prima.lory $ local) [prima.info $ GSP == cur.sp, ]
              
              cur.mod <- prima.aux $ models $ MOD [prima.aux $ models $ GSP == cur.sp]
              
              cur.info <- subset(prima.info, GSP == cur.sp)
              
              cur.Pmat <-
                  PmatrixPosterior(data = cur.data,
                                   model = cur.mod,
                                   info = cur.info,
                                   full.output = FALSE,
                                   thin = 10)
              
              cur.Pmat
          }, .parallel = TRUE); save.image()

dimnames(other.Pmat)[[1]] <-
    as.character(prima.aux $ models $ GSP [prima.aux $ models $ MOD != 'SEX' &
                                           prima.aux $ models $ MOD != 'NONE'])

post.Pmat <- array(0, c(nrow(prima.aux $ models), 100, 39, 39))

post.Pmat [1:53, , , ] <- none.Pmat
post.Pmat [(1:78) + 53, , , ] <- sex.Pmat
post.Pmat [(1:11) + 53 + 78, , , ] <- other.Pmat

dimnames(post.Pmat) [[1]] <- c(dimnames(none.Pmat) [[1]],
                               dimnames(sex.Pmat) [[1]],
                               dimnames(other.Pmat) [[1]])

dimnames(post.Pmat)[[3]] <-
    dimnames(post.Pmat)[[4]] <- c('logCS', colnames(prima.lory $ local))

prima.lory $ post.Pmat <- aperm(post.Pmat, c(3, 4, 2, 1))

prima.lory $ post.Pmat <-
    prima.lory $ post.Pmat [, , , match(prima.aux $ models $ GSP,
                                        dimnames(prima.lory $ post.Pmat) [[4]])]

save(prima.lory, file = 'Primates/LORY.RData')

## compare with ML estimates

ml.Pmat <-
    aaply(as.matrix(prima.aux $ models), 1,
          function(spandmodel)
          {
              ## type III sum of squares
              options(contrasts = c('contr.sum', 'contr.poly'))
              
              sp <- spandmodel [1]
              model <- spandmodel [2]
              
              info <- subset(prima.info, GSP == sp)
              
              data <-
                  cbind(log(prima.lory $ cs),
                        prima.lory $ local) [prima.info $ GSP == sp, ]
              
              if(model == 'NONE')
                  form <- as.formula('~ 1')
              else
                  form <- as.formula(paste('~', model))
              
              modmat <- model.matrix(form, data = info)
              
              CalculateMatrix(lm(data ~ modmat))
          }, .parallel = TRUE)

prima.lory $ ml.Pmat <- aperm(ml.Pmat, c(2, 3, 1))

prima.lory $ post.ml.rs <-
    aaply(1:142, 1, function(i)
    {
        aaply(prima.lory $ post.Pmat [, , , i], 3, RandomSkewers,
              cov.y = prima.lory $ ml.Pmat [, , i]) [, 1]
    }, .parallel = TRUE)

prima.aux $ sample.sizes <- table(prima.info $ GSP)

ssize.filter <-
    prima.aux $ sample.sizes[names(prima.aux $ sample.sizes) %in% prima.aux $ models $ GSP]

prima.aux $ models $ sampleSize <-
    ssize.filter [match(prima.aux $ models $ GSP, names(ssize.filter))]

colnames(prima.aux $ models)[3] <- 'SSIZE'

post.ml.df <- data.frame('id' = prima.aux $ models $ GSP,
                         'ssize' = as.integer(prima.aux $ models [, 'SSIZE']), 
                         'it' = prima.lory $ post.ml.rs)

post.ml.df <- gather(post.ml.df, key = id, value = value, -id, -ssize) [, -3]

pdf('box_post_ml.pdf', width = 12, height = 7)

ggplot(post.ml.df) +
    geom_boxplot(aes(x = id, y = value, fill = ssize)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 3)) +
    scale_fill_continuous(trans = 'log')

dev.off(dev.cur())

save(prima.lory, file = 'Primates/LORY.RData')

save(prima.aux, file = 'Primates/aux.RData')

## try with sym, see if it works

prima.sym $ gpa <- procGPA(prima.sym $ coord)

prima.sym $ tan <- prima.sym $ gpa $ tan

## trim dimensions of tan

dim(prima.sym $ tan) <- c(36, 3, 10081)

dimnames(prima.sym $ tan) <- dimnames(prima.sym $ coord)

prima.sym $ tan <- prima.sym $ tan [1:22, , ]

dimnames(prima.sym $ tan) [[1]] <- gsub('-E', '', dimnames(prima.sym $ tan) [[1]])

dimnames(prima.sym $ tan) [[2]] <- c('X', 'Y', 'Z')

coord.names <- paste(rep(dimnames(prima.sym $ tan) [[1]], each = 3),
                     rep(dimnames(prima.sym $ tan) [[2]], times = 22), sep = '.')

prima.sym $ tan <- aperm(prima.sym $ tan, c(2, 1, 3))

dim(prima.sym $ tan) <- c(66, 10081)

dimnames(prima.sym $ tan) <- list(coord.names, prima.info $ ID)

prima.sym $ tan <- t(prima.sym $ tan)

save(prima.sym, file = 'Primates/Sym.RData')

cur.sp <- 'Pan_paniscus'

cur.data <-
    cbind(log(prima.sym $ gpa $ size),
          prima.sym $ tan / mean(prima.sym $ gpa $ size)) [prima.info $ GSP == cur.sp, ]

aaply(cur.data, 2, var)

cur.mod <- prima.aux $ models $ MOD [prima.aux $ models $ GSP == cur.sp]
              
cur.info <- subset(prima.info, GSP == cur.sp)
              
cur.Pmat <-
    PmatrixPosterior(data = cur.data,
                     model = cur.mod,
                     info = cur.info,
                     full.output = FALSE,
                     thin = 10)

pdf(file = 'bonobo_sym_post.pdf')
color2D.matplot(cov2cor(cur.Pmat [99, , ]), xlab = 'geonormED')
dev.off(dev.cur())

eigen(cur.Pmat [99, , ]) $ values

ml.sym.test <- aaply(cur.Pmat, 1, RandomSkewers, cov.y = var(cur.data)) [, 1]



## ED

prima.ed $ geo.norm <-
    aaply(prima.ed $ raw, 1, function(raw)
    {
        geomean <- mean(log(raw))
        c(geomean, raw/exp(geomean))
    })

colnames(prima.ed $ geo.norm) <- c('logGM', colnames(prima.ed $ raw))
rownames(prima.ed $ geo.norm) <- prima.info $ ID

save(prima.ed, file = 'Primates/ed.RData')

## raw first

registerDoMC(cores = 72)

none.Pmat <- 
    aaply(as.character(prima.aux $ models $ GSP [prima.aux $ models $ MOD == 'NONE']), 1,
          function(cur.sp)
          {
              cur.data <- prima.ed $ raw [prima.info $ GSP == cur.sp, ]
              
              cur.mod <- prima.aux $ models $ MOD [prima.aux $ models $ GSP == cur.sp]
              
              cur.info <- subset(prima.info, GSP == cur.sp)
              
              cur.Pmat <-
                  PmatrixPosterior(data = cur.data,
                                   model = cur.mod,
                                   info = cur.info,
                                   modelpath = '../Stan/pmatrix_lm_norm.stan',
                                   full.output = FALSE,
                                   thin = 10)
              
              cur.Pmat
          }, .parallel = TRUE)

sex.Pmat <- 
    aaply(as.character(prima.aux $ models $ GSP [prima.aux $ models $ MOD == 'SEX']), 1,
          function(cur.sp)
          {
              cur.data <- prima.ed $ raw [prima.info $ GSP == cur.sp, ]
              
              cur.mod <- prima.aux $ models $ MOD [prima.aux $ models $ GSP == cur.sp]
              
              cur.info <- subset(prima.info, GSP == cur.sp)
              
              cur.Pmat <-
                  PmatrixPosterior(data = cur.data,
                                   model = cur.mod,
                                   info = cur.info,
                                   modelpath = '../Stan/pmatrix_lm_norm.stan',
                                   full.output = FALSE,
                                   thin = 10)
              
              cur.Pmat
          }, .parallel = TRUE)

other.Pmat <-
    aaply(as.character(prima.aux $ models $ GSP [prima.aux $ models $ MOD != 'SEX' &
                                                 prima.aux $ models $ MOD != 'NONE']), 1,
          function(cur.sp)
          {
              cur.data <- prima.ed $ raw [prima.info $ GSP == cur.sp, ]
              
              cur.mod <- prima.aux $ models $ MOD [prima.aux $ models $ GSP == cur.sp]
              
              cur.info <- subset(prima.info, GSP == cur.sp)
              
              cur.Pmat <-
                  PmatrixPosterior(data = cur.data,
                                   model = cur.mod,
                                   info = cur.info,
                                   modelpath = '../Stan/pmatrix_lm_norm.stan',
                                   full.output = FALSE,
                                   thin = 10)
              
              cur.Pmat
          }, .parallel = TRUE)

post.Pmat <- array(0, c(nrow(prima.aux $ models), 100, 38, 38))

post.Pmat [1:53, , , ] <- none.Pmat
post.Pmat [(1:78) + 53, , , ] <- sex.Pmat
post.Pmat [(1:11) + 53 + 78, , , ] <- other.Pmat

dimnames(post.Pmat) [[1]] <- 
    c(as.character(prima.aux $ models $ GSP [prima.aux $ models $ MOD == 'NONE']), 
      as.character(prima.aux $ models $ GSP [prima.aux $ models $ MOD == 'SEX']),
      as.character(prima.aux $ models $ GSP [prima.aux $ models $ MOD != 'SEX' &
                                             prima.aux $ models $ MOD != 'NONE']))

dimnames(post.Pmat)[[3]] <-
    dimnames(post.Pmat)[[4]] <-
    colnames(prima.lory $ local)

prima.ed $ post.Pmat <- aperm(post.Pmat, c(3, 4, 2, 1))

prima.ed $ post.Pmat <-
    prima.ed $ post.Pmat [, , , match(prima.aux $ models $ GSP,
                                        dimnames(prima.ed $ post.Pmat) [[4]])]

save(prima.ed, file = 'Primates/ed.RData')

## compare with ML estimates

ml.Pmat <-
    aaply(as.matrix(prima.aux $ models), 1,
          function(spandmodel)
          {
              ## type III sum of squares
              options(contrasts = c('contr.sum', 'contr.poly'))
              
              sp <- spandmodel [1]
              model <- spandmodel [2]
              
              info <- subset(prima.info, GSP == sp)
              
              data <- prima.ed $ raw [prima.info $ GSP == sp, ]
              
              if(model == 'NONE')
                  form <- as.formula('~ 1')
              else
                  form <- as.formula(paste('~', model))
              
              modmat <- model.matrix(form, data = info)
              
              CalculateMatrix(lm(data ~ modmat))
          }, .parallel = TRUE)

prima.ed $ ml.Pmat <- aperm(ml.Pmat, c(2, 3, 1))

dimnames(prima.ed $ ml.Pmat)[[3]] <- as.character(prima.aux $ models [, 'GSP'])

prima.ed $ post.ml.rs <-
    aaply(1:142, 1, function(i)
    {
        aaply(prima.ed $ post.Pmat [, , , i], 3, RandomSkewers,
              cov.y = prima.ed $ ml.Pmat [, , i]) [, 1]
    }, .parallel = TRUE)

post.ml.df <- data.frame('id' = prima.aux $ models $ GSP,
                         'ssize' = as.integer(prima.aux $ models [, 'SSIZE']), 
                         'it' = prima.ed $ post.ml.rs)

post.ml.df <- gather(post.ml.df, key = id, value = value, -id, -ssize) [, -3]

pdf('box_post_ml.pdf', width = 12, height = 7)

ggplot(post.ml.df) +
    geom_boxplot(aes(x = id, y = value, fill = ssize)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 3)) +
    scale_fill_continuous(trans = 'log')

dev.off(dev.cur())

pdf('sag_mat.pdf', width = 14, height = 7)
par(mfrow = c(1, 2))
color2D.matplot(cov2cor(prima.ed $ post.Pmat [, , 74, 'Saguinus_midas']))
color2D.matplot(cov2cor(prima.ed $ ml.Pmat [, , 'Saguinus_midas']))
dev.off(dev.cur())

## test log

prima.log <- list()

prima.log $ raw <- log(prima.ed $ raw)

prima.log $ pca <- prcomp(prima.log $ raw, retx = TRUE)

grad <- colorRampPalette(brewer.pal(8, "Spectral"), space="Lab")

ggsave(
    'primate_log.pdf',
    width = 18, height = 10, 
    plot = data.frame(
        gen = prima.info $ GEN, 
        prima.log $ pca $ x [, 1:2]) %>%
        ggplot(., aes(x = PC1, y = PC2, color = gen, group = gen)) +
        geom_point() +
        stat_ellipse() +
        theme_bw() +
        scale_color_manual('Genus', values = grad(length(unique(prima.info $ GEN)))) +
        guides(color = guide_legend(ncol = 3))
    )

registerDoMC(cores = 32)

prima.log $ post.Pmat <- 
    aaply(as.character(prima.aux $ models $ GSP), 1,
          function(cur.sp)
          {
              cur.data <- prima.log $ raw [prima.info $ GSP == cur.sp, ]
              
              cur.mod <- prima.aux $ models $ MOD [prima.aux $ models $ GSP == cur.sp]
              
              cur.info <- subset(prima.info, GSP == cur.sp)
              
              cur.Pmat <-
                  PmatrixPosterior(data = cur.data,
                                   model = cur.mod,
                                   info = cur.info,
                                   modelpath = '../Stan/pmatrix_lm_norm.stan',
                                   full.output = FALSE,
                                   thin = 10)
              
              cur.Pmat
          }, .parallel = TRUE)

dimnames(post.Pmat) [[1]] <- 
    c(as.character(prima.aux $ models $ GSP [prima.aux $ models $ MOD == 'NONE']), 
      as.character(prima.aux $ models $ GSP [prima.aux $ models $ MOD == 'SEX']),
      as.character(prima.aux $ models $ GSP [prima.aux $ models $ MOD != 'SEX' &
                                             prima.aux $ models $ MOD != 'NONE']))

dimnames(post.Pmat)[[3]] <-
    dimnames(post.Pmat)[[4]] <-
    colnames(prima.lory $ local)

prima.log $ post.Pmat <- aperm(post.Pmat, c(3, 4, 2, 1))

prima.log $ post.Pmat <-
    prima.log $ post.Pmat [, , , match(prima.aux $ models $ GSP,
                                        dimnames(prima.log $ post.Pmat) [[4]])]

save(prima.log, file = 'Primates/ed.RData')

## compare with ML estimates

ml.Pmat <-
    aaply(as.matrix(prima.aux $ models), 1,
          function(spandmodel)
          {
              ## type III sum of squares
              options(contrasts = c('contr.sum', 'contr.poly'))
              
              sp <- spandmodel [1]
              model <- spandmodel [2]
              
              info <- subset(prima.info, GSP == sp)
              
              data <- prima.log $ raw [prima.info $ GSP == sp, ]
              
              if(model == 'NONE')
                  form <- as.formula('~ 1')
              else
                  form <- as.formula(paste('~', model))
              
              modmat <- model.matrix(form, data = info)
              
              CalculateMatrix(lm(data ~ modmat))
          }, .parallel = TRUE)

prima.log $ ml.Pmat <- aperm(ml.Pmat, c(2, 3, 1))

dimnames(prima.log $ ml.Pmat)[[3]] <- as.character(prima.aux $ models [, 'GSP'])

prima.log $ post.ml.rs <-
    aaply(1:142, 1, function(i)
    {
        aaply(prima.log $ post.Pmat [, , , i], 3, RandomSkewers,
              cov.y = prima.log $ ml.Pmat [, , i]) [, 1]
    }, .parallel = TRUE)

post.ml.df <- data.frame('id' = prima.aux $ models $ GSP,
                         'ssize' = as.integer(prima.aux $ models [, 'SSIZE']), 
                         'it' = prima.log $ post.ml.rs)

post.ml.df <- gather(post.ml.df, key = id, value = value, -id, -ssize) [, -3]

pdf('box_post_ml.pdf', width = 12, height = 7)

ggplot(post.ml.df) +
    geom_boxplot(aes(x = id, y = value, fill = ssize)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 3)) +
    scale_fill_continuous(trans = 'log')

dev.off(dev.cur())

pdf('sag_mat.pdf', width = 14, height = 7)
par(mfrow = c(1, 2))
color2D.matplot(cov2cor(prima.log $ post.Pmat [, , 74, 'Saguinus_midas']))
color2D.matplot(cov2cor(prima.log $ ml.Pmat [, , 'Saguinus_midas']))
dev.off(dev.cur())
