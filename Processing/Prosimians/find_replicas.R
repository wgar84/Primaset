if(!require(plyr)) {install.packages('plyr'); library(plyr)}
if(!require(dplyr)) {install.packages('dplyr'); library(dplyr)}
if(!require(magrittr)) {install.packages('magrittr'); library(magrittr)}
if(!require(lme4)) {install.packages('lme4'); library(lme4)}
if(!require(ggplot2)) {install.packages('ggplot2'); library(ggplot2)}
if(!require(tidyr)) {install.packages('tidyr'); library(tidyr)}
#if(!require(MCMCglmm)) {install.packages('MCMCglmm'); library(MCMCglmm)}
if(!require(reshape2)) {install.packages('reshape2'); library(reshape2)}
#if(!require(evolqg)) {install.packages('evolqg'); library(evolqg)}
if(!require(evolqg)) {devtools::install_github('lem-usp/evolqg'); library(evolqg)}
if(!require(readr)) {devtools::install_github('hadley/readr'); library(readr)}
if(!require(doParallel)) {install.packages('doParallel'); library(doParallel)}
#Registrando o numero de cores : 3 em casa, 7 no lab e até 10 no darwin
#para descobrir rodar no terminal: nproc
#abrir no terminal htop para ver os cores trabalhando
registerDoParallel(cores = 2)
#abrir no terminal htop para ver os cores trabalhando

load("~/Science_Code/LEM/Primaset/Primaset/Processing/.RData")

#Prosimian Shapes
str(prosimian.shapes)
#> num [1:44, 1:3, 1:2, 1:1867] # [landmarks, coordinates, replicas, individuals]
dimnames(prosimian.shapes) <- list(landmarks.that.matter,
                                LETTERS[24:26],
                                paste0("take", 1:2),
                                paste(prosimian.id$Ind))

#Exemplos de medidas
D2E2 = "AMNH170494"  # duas réplicas dos dois lados
D1E1 = "RMNH041010e" # uma réplica de ambos os lados
D0E1 = "USNM452953"  # uma réplica só lado esquerdo
D1E0 = "USNM452954"  # uma réplica só lado direito

teste = D2E2
teste = D1E1
teste = D0E1
teste = D1E0

#réplicas são diferentes?
diffR1R2 = prosimian.shapes[ , , 1,] - prosimian.shapes[ , , 2,]
dimnames(diffR1R2)[[1]] <- c(1:44)
diffR1R2 = diffR1R2 != 0 # tem uma réplica diferente da outra? T== replicas diferentes
diffR1R2 = aaply(diffR1R2, c(3, 2), function (x) sum(x == T, na.rm = T))
diffR1R2 = rowSums(diffR1R2)
table(diffR1R2 == 0)

dimnames(prosimian.shapes)[[4]] [diffR1R2 >67]

only.D = is.na(prosimian.shapes[1:22 , , 1, ]) # se só foi medido o lado direito, o lado esquerdo tem tudo NA
only.D = aaply(only.D, c(3, 2), function (x) all(x == T, na.rm = T))
only.D = rowSums(only.D)
table(only.D == 66)

only.D = is.na(prosimian.shapes[23:44 , , 1, ])
only.D = aaply(only.D, c(3, 2), function (x) sum(x == T, na.rm = T))
only.D = rowSums(only.D)
table(only.D == 66)

sum(diffR1R2 == T, na.rm = T) # == 132 para D2L2; == 0 para D0L1 ou D1L0; == 0 para D1L1.
sum(diffR1R2 == F, na.rm = T) # == 0 para D2L2; == 66 para D0L1 ou D1L0; == 132 para D1L1.

(is.na(diffR1R2) == T) # se for == 66 foi medido so um lado; se form >66 é lixo

table(prosimian.shapes[1:22 , , 1, ] == prosimian.shapes[23:44 , , 1, ]) # esquerdo diferente do direito?
sum(is.na(prosimian.shapes[1:22, , 1, 1234]) == TRUE) ==  66 | sum(is.na(prosimian.shapes[1:22, , 2, 1234]) == TRUE) ==  66  # tem só um dos lados medidos?
table(prosimian.shapes[1:22, , 1, "USNM_452954"] == prosimian.shapes[23:44 , , 1, "USNM_452954"])



names(prosimian.shapes)[1]
