require(MCMCglmm)
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
require(viridis)
require(phytools)

## install.packages(
##     c(
##         'devtools', 
##         'shapes', 
##         'geomorph', 
##         'rstan', 
##         'plyr', 
##         'RColorBrewer', 
##         'evolqg', 
##         'plotrix', 
##         'doMC', 
##         'dplyr', 
##         'magrittr', 
##         'tidyr', 
##         'viridis', 
##         'phytools'
##      ), dependencies = TRUE)


registerDoMC(cores = 3)

load('Primates/Info.RData')
load('Primates/Sym.RData')
load('Primates/newaux.RData')
load('Primates/ed.RData')
load('Primates/LORY.RData')
load('Primates/newaux.RData')
load('../Raw Data/Aux.RData')
load('Primates/allo.tree.RData')

devtools::load_all('/home/guilherme/Dropbox/lem_essentials/Primaset')

allo.tree $ models <-
                prima.aux $ models [prima.aux $ models $ GSP %in% rownames(allo.tree $ data), ]

not.in.models <-
    allo.tree $ phy $ tip.label [!allo.tree $ phy $ tip.label %in% allo.tree $ models $ GSP]

add.to.models <-
    as.matrix(allo.tree $ data [rownames(allo.tree $ data) %in% not.in.models, ])

add.to.models <-
    data.frame(rownames(add.to.models),
               rep('NONE', nrow(add.to.models)),
               add.to.models[, 1],
               rep(FALSE, nrow(add.to.models)))

allo.tree $ models $ USE <- rep(TRUE, nrow(allo.tree $ models))

colnames(add.to.models) <- colnames(allo.tree $ models)

rownames(add.to.models) <- NULL

## belzebuth

add.to.models [add.to.models $ GSP == 'Ateles_belzebuth', 'USE'] <- TRUE

add.to.models <- subset(add.to.models, SSIZE > 1)

allo.tree $ models <- rbind(allo.tree $ models, add.to.models)

rownames(allo.tree $ models) <- allo.tree $ models $ GSP

allo.data <- treedata(allo.tree $ phy, allo.tree $ models)

allo.data $ models <- allo.tree $ models

allo.data <- allo.data [-2]

names(allo.data)

allo.data $ Pmat <-
                alply(as.matrix(allo.data $ models), 1,
                      function(spandmodel)
                      {
                          ## type III sum of squares
                          options(contrasts = c('contr.sum', 'contr.poly'))
                          
                          sp <- spandmodel [1]
                          model <- spandmodel [2]
                          
                          info <- subset(prima.info, SPR == sp)
              
                          data <- prima.lory $ local [prima.info $ SPR == sp, ]
              
                          if(model == 'NONE')
                              form <- as.formula('~ 1')
                          else
                              form <- as.formula(paste('~', model))
          
                          modmat <- model.matrix(form, data = info)
          
                          CalculateMatrix(lm(data ~ modmat))
                      }, .parallel = TRUE)

allo.data $ models $ DF <-
                         aaply(as.matrix(allo.data $ models), 1,
                               function(spandmodel)
                               {
                                   ## type III sum of squares
                                   options(contrasts = c('contr.sum', 'contr.poly'))
                                   
                                   sp <- spandmodel [1]
                                   model <- spandmodel [2]
                                   
                                   info <- subset(prima.info, SPR == sp)
                                   
                                   data <- prima.lory $ local [prima.info $ SPR == sp, ]
                                   
                                   if(model == 'NONE')
                                       form <- as.formula('~ 1')
                                   else
                                       form <- as.formula(paste('~', model))
                                   
                                   modmat <- model.matrix(form, data = info)
                                   
                                   lm(data ~ modmat) $ df.residual
                               }, .parallel = TRUE)


names(allo.data $ Pmat) <- allo.data $ models $ GSP

allo.data $ Wanc <- PhyloW(allo.data $ phy,
                           allo.data $ Pmat,
                           allo.data $ models $ DF)

allo.data $ le.Wanc  <- allo.data $ Wanc [[length(allo.data $ Wanc)]]

## plot W
## source('1A_plot_W.R')

allo.data $ eigenWanc <- eigen(allo.data $ le.Wanc)

## to project

spandmodel <- as.matrix(allo.data $ models) [1, ]

allo.data$to.project <-
    alply(as.matrix(allo.data $ models), 1,
          function(spandmodel)
          {
              ## type III sum of squares
              options(contrasts = c('contr.sum', 'contr.poly'))
              
              sp <- spandmodel [1]
              model <- spandmodel [2]
              
              info <- subset(prima.info, SPR == sp)
              
              data <- cbind(log(prima.lory $ cs),
                            prima.lory $ local)[prima.info $ SPR == sp, ]
              
              if(model == 'NONE')
                  form <- as.formula('~ 1')
              else
                  form <- as.formula(paste('~', model))
              
              modmat <- model.matrix(form, data = info)
              
              le.model <- lm(data ~ modmat)
              
              t(colMeans(fitted(le.model)) + t(residuals(le.model)))
          }, .parallel = TRUE)

allo.data $ to.project.mat <- NULL
for(i in 1:length(allo.data $ to.project))
    {
        rownames(allo.data $ to.project [[i]]) <-
            rep(allo.data $ models $ GSP [i], nrow(allo.data $ to.project [[i]]))

        allo.data $ to.project.mat <- rbind(allo.data $ to.project.mat,
                                            allo.data $ to.project [[i]])
    }

colnames(allo.data $ to.project.mat) [1] <- 'lnCS'

allo.data $ projected.Wanc <-
                t(allo.data $ eigenWanc $ vectors [, 1:4]) %*%
                t(allo.data $ to.project.mat [, -1])

allo.data$main.df <-
    data.frame('animal' = rownames(allo.data $ to.project.mat),
               'lnCS' = allo.data $ to.project.mat [, 1],
               'PC' = t(allo.data $ projected.W))

allo.data$main.df$genus <-
    laply(strsplit(as.character(allo.data $ main.df $ animal), '_'),
      function(L) L[1])

ggplot(allo.data $ main.df) +
    geom_point(aes(x = lnCS, y = PC.2, color = genus))

allo.data $ main.df $ lnCS <- scale(allo.data $ main.df $ lnCS, scale = FALSE)

ggplot(allo.data $ main.df) +
    geom_point(aes(x = lnCS, y = PC.1, color = genus))


save(allo.data, file = 'Primates/allo.RData')

## LE MODEL

pdf ('achb.allo.pdf', width = 8, height = 30)
plot (allo.data $ phy, cex = 0.5, show.node.label = TRUE, direction = 'rightwards')
#nodelabels ()
dev.off (dev.cur ())

diag(vcv.phylo(allo.data $ phy))

allo.data $ phy.force <- force.ultrametric(allo.data $ phy, method = 'nnls')

save(allo.data, file = 'Primates/allo.RData')

allo.data <-
  within (allo.data,
          {
            allo.mcmc.randomreg <-
              MCMCglmm(PC.1 ~ lnCS,
                       random = ~ us(1 + lnCS):animal,
                       rcov = ~ units,
                       data = main.df, pr = TRUE,
                       nitt = 100000, burnin = 50000, thin = 50,
                       pedigree = phy.force,
                       family = 'gaussian')
          })

allo.data $ main.df $ predPC1 <- predict(allo.data $ allo.mcmc.randomreg, marginal = NULL)

ggplot(allo.data $ main.df) +
    geom_line(aes(x = lnCS, y = predPC1, color = genus, group = animal))


allo.data <-
  within (allo.data,
          {
            multi.rreg <-
              MCMCglmm(cbind(PC.1, PC.2) ~ trait * lnCS,
                       random = ~ us(trait * lnCS):animal,
                       rcov = ~ us(trait):units,
                       data = main.df, pr = TRUE,
                       #nitt = 100000, burnin = 50000, thin = 50,
                       pedigree = phy.force,
                       family = rep('gaussian', 2))
          })

allo.data $ main.df[, c('mpPC1', 'mpPC2', 'mpPC3')] <-
                predict(allo.data $ multi.rreg, marginal = NULL)

ggplot(allo.data $ main.df) +
    geom_point(aes(x = lnCS, y = PC.1), alpha = 0.1) +
    geom_line(aes(x = lnCS, y = mpPC1, color = genus, group = animal))

ggplot(allo.data $ main.df) +
    geom_point(aes(x = lnCS, y = PC.2), alpha = 0.1) +
    geom_line(aes(x = lnCS, y = mpPC2, color = genus, group = animal))

ggplot(subset(allo.data $ main.df, genus == 'Loris')) +
##ggplot(allo.data $ main.df) +
    geom_text(aes(x = lnCS, y = PC.3, label = substr(genus, 1, 3)), alpha = 0.1) +
    geom_line(aes(x = lnCS, y = mpPC3, color = genus, group = animal))


allo.data $ ssizes <- table(allo.data $ main.df $ animal)

allo.data $ models $ SSIZE <-
                         allo.data $ ssizes [match(allo.data $ models $ GSP,
                                                   names(allo.data $ ssizes))]

## Stan

## rstan_options(auto_write = TRUE)
## options(mc.cores = 32)

## phylo (Ainv)
covphylo.to.stan <- inverseA(allo.data $ phy.force)

allo.data$main.df$animal <-
    factor(as.character(allo.data $ main.df $ animal),
           levels = allo.data $ phy.force $ tip.label)

pc.to.stan <-
    dlply(allo.data $ main.df, .(animal),
          function(sub.df)
          {
              as.matrix(sub.df [, c('lnCS', 'PC.1', 'PC.2', 'PC.3', 'PC.4')])
          })

max(allo.data $ models $ SSIZE)

         
pc.to.stan <-
    laply(pc.to.stan,
          function(block)
          {
              max.nrow <- max(allo.data $ models $ SSIZE)

              if(max.nrow == nrow(block))
                  return(block)

              rbind(block, 
                    matrix(0, nrow = max.nrow - nrow(block), ncol = ncol(block)))
          })

## montando lista

allo.data $ list4stan <- 
    list('k' = dim(pc.to.stan)[3] - 1,
         'm' = dim(pc.to.stan)[1],
         'cov_phylo' = as.matrix(covphylo.to.stan $ Ainv),
         'ni' = as.vector(table(allo.data$main.df$animal)),
         'ni_max' = max(as.vector(table(allo.data$main.df$animal))),
         'Y' = pc.to.stan[, , 2:dim(pc.to.stan)[3]],
         'X' = pc.to.stan[, , 1])

allo.data $ initialConditions4stan <-
    function(i = 1, k, m)
        list('as_term' = array(rnorm(k * m, 3), c(m, k)),
             'bs_term' = array(rnorm(k * m, 3), c(m, k)),
             'as_anc' = array(rnorm(k * (m - 2), 3), c(m - 2, k)),
             'bs_anc' = array(rnorm(k * (m - 2), 3), c(m - 2, k)),
             'as_root' = rnorm(k, 3),
             'bs_root' = rnorm(k, 3),
             'Sigma_e' = var(array(rnorm(k ^ 2), c(k, k))),
             'Sigma_a' = var(array(rnorm(k ^ 2), c(k, k))),
             'Sigma_b' = var(array(rnorm(k ^ 2), c(k, k))))

options(mc.cores = 2)

allo.data $ stan.out.4PC <-
                stan(file = '../Stan/random_reg_One.stan',
                     data = allo.data $ list4stan,
                     chains = 1)
                     ## init =
                     ##     list(allo.data $ initialConditions4stan(
                     ##                          k = allo.data $ list4stan $ k,
                     ##                          m = allo.data $ list4stan $ m)))

plot(allo.data $ stan.out.4PC)
