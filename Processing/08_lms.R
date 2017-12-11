require(devtools)
require(shapes)
require(geomorph)
require(rstan)
require(plyr)
require(RColorBrewer)
require(evolqg)
require(plotrix)
require(doMC)

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

prima.info [4475, 'SEX'] <- 'M' ## whatever

save(prima.info, file = 'Primates/Info.RData')

devtools::load_all('/home/guilherme/lem_essentials/Primaset')

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

dimnames(prima.lory $ post.ml.rs) [[1]] <- prima.aux $ models $ GSP

post.ml.df <- gather(data.frame('id' = prima.aux $ models $ GSP,
                                prima.lory $ post.ml.rs), key = id) [, c(1, 3)]

pdf('box_post_ml.pdf', width = 12, height = 7)

ggplot(post.ml.df) +
    geom_boxplot(aes(x = id, y = value)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 3))

dev.off(dev.cur())

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
