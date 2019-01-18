require(geomorph)
require(shapes)
require(evolqg)
require(plotrix)

load('../Raw Data/Aux.RData')

load('Primates/Sym.RData')

load('Primates/Info.RData')

dim(prima.sym$coord)

sag.sym <- prima.sym $ coord [, , prima.info $ GSP == 'Saguinus_geoffroyi']

sag.sizeshape.gpa <- procGPA(sag.sym, scale = FALSE)

sag.ss.tan <- sag.sizeshape.gpa $ tan

dim(sag.ss.tan) <- c(36, 3, 109)

dimnames(sag.ss.tan)[1:2] <- dimnames(prima.sym $ coord)[1:2]

sag.ss.tan <- sag.ss.tan [1:22, , ]

dimnames(sag.ss.tan) [[1]] <- gsub('-E', '', dimnames(sag.ss.tan) [[1]])

dimnames(sag.ss.tan) [[2]] <- c('X', 'Y', 'Z')

coord.names <- paste(rep(dimnames(sag.ss.tan) [[1]], each = 3),
                     rep(dimnames(sag.ss.tan) [[2]], times = 22), sep = '.')

sag.ss.tan <- aperm(sag.ss.tan, c(2, 1, 3))

dim(sag.ss.tan) <- c(66, 109)

dimnames(sag.ss.tan) <- list(coord.names,
                             prima.info $ ID [prima.info $ GSP == 'Saguinus_geoffroyi'])

sag.ss.tan <- t(sag.ss.tan)

par(mfrow = c(1, 2))
color2D.matplot(cor(sag.ss.tan))
plot(eigen(var(sag.ss.tan)) $ values)

rownames(Aux $ sym.hyp [[1]])

sym.hyps <- Aux $ sym.hyp [[1]] [1:66, ]

rownames(sym.hyps) <- gsub('-D', '', rownames(sym.hyps))

sym.hyps <- sym.hyps [match(colnames(sag.ss.tan), rownames(sym.hyps)), ]

neuroface <-
    sym.hyps [, 'Neuro'] %*% t(sym.hyps [, 'Neuro']) +
    sym.hyps [, 'Face'] %*% t(sym.hyps [, 'Face'])

color2D.matplot(neuroface)

MantelModTest(neuroface, cor(sag.ss.tan), landmark.dim = 3,
              withinLandmark = FALSE, MHI = TRUE)

sag.cormat <- cor(sag.ss.tan)

sag.vcv <- var(sag.ss.tan)

sag.evec3 <- eigen(sag.vcv) $ vectors [, 1:3]

rownames(sag.evec3) <- colnames(sag.ss.tan)

hist(abs(sag.cormat) [which(sym.hyps [, 'Neuro'] == 1), which(sym.hyps [, 'Neuro'] == 1)])

hist(abs(sag.cormat) [which(sym.hyps [, 'Neuro'] == 1), which(sym.hyps [, 'Face'] == 1)])
