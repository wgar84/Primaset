load('Primates/Info.RData')
load('Primates/Raw.RData')

primaSym2D = prima.raw $ coord

dimnames(primaSym2D) [[2]] = c('X', 'Y', 'Z')

lms = rownames(primaSym2D)

dims = colnames(primaSym2D)

dim(primaSym2D) = c(36*3, 10081)

primaSym2D = t(primaSym2D)

colnames(primaSym2D) = paste(rep(lms, times = 3), rep(dims, each = 36), sep = '-')

whoNWM = which(prima.info $ MAJOR == 'Platyrrhini')

out = cbind(prima.info[whoNWM, ], primaSym2D [whoNWM, ])

write.csv(out, file = 'nwm.csv')
