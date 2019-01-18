require (shapes)

### banco de dados
data.dist = read.csv2 (file = 'GLMALL2.csv', header = TRUE, row.names = NULL)

### landmarks
vis = scan ("vistas.csv", what = "")
lmA = vis[1:35]
lmA[4] = "NA"
lmZ = vis[36:length (vis)]
rm (vis)

### AMERICAN MUSEUM OF NATURAL HISTORY

### get skulls
AMNH = get.skull (DIR = "AMNHPRO", lmA = lmA, lmZ = lmZ)
AMNH.rep = get.skull (DIR = "AMNHPRO", lmA = lmA, lmZ = lmZ, rep = TRUE)

### AMNH remove ind w/ missing lms
amnh.miss.A = apply (AMNH$A, 3, find.miss.lm)
amnh.miss.Z = apply (AMNH$Z, 3, find.miss.lm)

amnh.rem = amnh.miss.A & amnh.miss.Z

### 136215 zuado (falta MT-E)

amnh.rem [dimnames (AMNH$A) [[3]] == '136215'] = FALSE

### compare svd with ols
amnh.svd = glue.skulls (AMNH$A[,,amnh.rem], AMNH$Z[,,amnh.rem], soln = 'svd')
amnh.ols = glue.skulls (AMNH$A[,,amnh.rem], AMNH$Z[,,amnh.rem], soln = 'ols')

### distancias com as e ld

dists = t (array (dim = c (2,6),
  data = c (
    'BR', 'LD',
    'OPI', 'LD',
    'PT-D', 'AS-D',
    'JP-D', 'AS-D',
    'PT-E', 'AS-E',
    'JP-E', 'AS-E'
    )))

comp.dist = array (0, c(2 * dim (amnh.svd)[3], 6))

for (i in 1:6)
  {
    k = which (rownames (amnh.svd) %in% dists[i,])
    tmp1 = c (apply (amnh.svd [k,,], 3, dist), apply (amnh.ols [k,,], 3, dist))
    comp.dist [,i]  = tmp1
  }
rm (i, k)

comp.fac = factor (names (tmp1))

comp.lm = lm (comp.dist ~ comp.fac)
anova (comp.lm, test = 'Wilks') ### perfeito, tanto faz

### print.skull

print.skull (amnh.svd[,,1])

### carregar banco de dados de distâncias

dimnames (amnh.svd)[[3]] %in% as.character (data.dist [,1] [data.dist$MUSEUM == 'AMNH'])

amnh.raw = amnh.svd
amnh.data = data.dist [data.dist$MUSEUM == 'AMNH',1:10]

elem1 = function (element) {return (element[1])}
amnh.id = sapply (strsplit (as.character (amnh.data$ID), split = '.', fixed = TRUE), elem1)
amnh.data$ID = amnh.id

amnh.match = amnh.id %in% dimnames (amnh.raw) [[3]]

amnh.data = amnh.data [amnh.match,]

amnh.match = dimnames (amnh.raw) [[3]] %in% amnh.id

amnh.raw = amnh.raw [,,amnh.match]
amnh.order = match (amnh.data$ID, dimnames (amnh.raw) [[3]])

amnh.raw = amnh.raw [,,amnh.order]

amnh.export = list ('data' = amnh.data,
  'raw' = amnh.raw)

### checando formas

amnh.data$GENUS == 'Aotus'

vis.seq (amnh.raw[,, amnh.data$GENUS == 'Aotus'])

save (amnh.export, file = 'N.AMNH.RData')

### MUSEU DE ZOOLOGIA

### get skulls
MZUSP = get.skull (DIR = "spproc2", lmA = lmA, lmZ = lmZ)
dim (MZUSP$A)

### AMNH remove ind w/ missing lms
mzusp.miss.A = apply (MZUSP$A, 3, find.miss.lm)
mzusp.miss.Z = apply (MZUSP$Z, 3, find.miss.lm)

mzusp.rem = mzusp.miss.A & mzusp.miss.Z
sum (mzusp.rem)

### compare svd with ols
mzusp.raw = glue.skulls (MZUSP$A[,,mzusp.rem], MZUSP$Z[,,mzusp.rem], soln = 'svd')

### carregar banco de dados de distâncias

dimnames (mzusp.raw)[[3]] %in% as.character (data.dist [,1] [data.dist$MUSEUM == 'MZUSP'])

mzusp.data = data.dist [data.dist$MUSEUM == 'MZUSP',1:10]

mzusp.id = mzusp.data$ID

#elem1 = function (element) {return (element[1])}
#mzusp.id = sapply (strsplit (as.character (mzusp.data$ID), split = '.', fixed = TRUE), elem1)
#mzusp.data$ID = mzusp.id

mzusp.match = mzusp.id %in% dimnames (mzusp.raw) [[3]]

mzusp.data = mzusp.data [mzusp.match,]

mzusp.match = dimnames (mzusp.raw) [[3]] %in% mzusp.id

mzusp.raw = mzusp.raw [,,mzusp.match]
mzusp.order = match (mzusp.data$ID, dimnames (mzusp.raw) [[3]])

mzusp.raw = mzusp.raw [,,mzusp.order]
dim (mzusp.raw)

mzusp.export = list ('data' = mzusp.data,
  'raw' = mzusp.raw)

### checando formas

vis.seq (mzusp.raw, start = 1)

save (mzusp.export, file = 'N.MZUSP.RData')

### SMITHSONIAN

### get skulls
USNM = get.skull (DIR = "USNMPRO", lmA = lmA, lmZ = lmZ)
dim (USNM$A)

### AMNH remove ind w/ missing lms
usnm.miss.A = apply (USNM$A, 3, find.miss.lm)
usnm.miss.Z = apply (USNM$Z, 3, find.miss.lm)

usnm.rem = usnm.miss.A & usnm.miss.Z
sum (usnm.rem)

### compare svd with ols
usnm.raw = glue.skulls (USNM$A[,,usnm.rem], USNM$Z[,,usnm.rem], soln = 'svd')

### carregar banco de dados de distâncias

dimnames (usnm.raw)[[3]] %in% as.character (data.dist [,1] [data.dist$MUSEUM == 'USNM'])

usnm.data = data.dist [data.dist$MUSEUM == 'USNM',1:10]

elem1 = function (element) {return (element[1])}
usnm.id = sapply (strsplit (as.character (usnm.data$ID), split = '.', fixed = TRUE), elem1)
usnm.data$ID = usnm.id


usnm.match = usnm.id %in% dimnames (usnm.raw) [[3]]

usnm.data = usnm.data [usnm.match,]

usnm.match = dimnames (usnm.raw) [[3]] %in% usnm.id

usnm.raw = usnm.raw [,,usnm.match]
usnm.order = match (usnm.data$ID, dimnames (usnm.raw) [[3]])

usnm.raw = usnm.raw [,,usnm.order]
dim (usnm.raw)

usnm.export = list ('data' = usnm.data,
  'raw' = usnm.raw)

### checando formas

vis.seq (usnm.raw, start = 1)

save (usnm.export, file = 'N.USNM.RData')

### MUSEU NACIONAL

### get skulls
MNRJ = get.skull (DIR = "rjproc", lmA = lmA, lmZ = lmZ)
dim (MNRJ$A)

### AMNH remove ind w/ missing lms
mnrj.miss.A = apply (MNRJ$A, 3, find.miss.lm)
mnrj.miss.Z = apply (MNRJ$Z, 3, find.miss.lm)

mnrj.rem = mnrj.miss.A & mnrj.miss.Z
sum (mnrj.rem)

### compare svd with ols
mnrj.raw = glue.skulls (MNRJ$A[,,mnrj.rem], MNRJ$Z[,,mnrj.rem], soln = 'svd')

### carregar banco de dados de distâncias

dimnames (mnrj.raw)[[3]] %in% as.character (data.dist [,1] [data.dist$MUSEUM == 'MNRJ'])

mnrj.data = data.dist [data.dist$MUSEUM == 'MNRJ',1:10]

mnrj.id = mnrj.data$ID

#elem1 = function (element) {return (element[1])}
#mnrj.id = sapply (strsplit (as.character (mnrj.data$ID), split = '.', fixed = TRUE), elem1)
#mnrj.data$ID = mnrj.id

mnrj.match = mnrj.id %in% dimnames (mnrj.raw) [[3]]

mnrj.data = mnrj.data [mnrj.match,]

mnrj.match = dimnames (mnrj.raw) [[3]] %in% mnrj.id

mnrj.raw = mnrj.raw [,,mnrj.match]
mnrj.order = match (mnrj.data$ID, dimnames (mnrj.raw) [[3]])

mnrj.raw = mnrj.raw [,,mnrj.order]
dim (mnrj.raw)

mnrj.export = list ('data' = mnrj.data,
  'raw' = mnrj.raw)

### checando formas

vis.seq (mnrj.raw, start = 1)

save (mnrj.export, file = 'N.MNRJ.RData')

### FIELD MUSEUM OF CHICAGO

FMNH = get.skull (DIR = "FMPROC", lmA = lmA, lmZ = lmZ)

### AMNH remove ind w/ missing lms
fmnh.miss.A = apply (FMNH$A, 3, find.miss.lm)
fmnh.miss.Z = apply (FMNH$Z, 3, find.miss.lm)

fmnh.rem = fmnh.miss.A & fmnh.miss.Z

### compare svd with ols
fmnh.raw = glue.skulls (FMNH$A[,,fmnh.rem], FMNH$Z[,,fmnh.rem], soln = 'svd')

### carregar banco de dados de distâncias

dimnames (fmnh.raw)[[3]] %in% as.character (data.dist [,1] [data.dist$MUSEUM == 'FMNH'])

fmnh.data = data.dist [data.dist$MUSEUM == 'FMNH',1:10]

fmnh.id = fmnh.data$ID

fmnh.match = fmnh.id %in% dimnames (fmnh.raw) [[3]]

fmnh.data = fmnh.data [fmnh.match,]

fmnh.match = dimnames (fmnh.raw) [[3]] %in% fmnh.id

fmnh.raw = fmnh.raw [,,fmnh.match]
fmnh.order = match (fmnh.data$ID, dimnames (fmnh.raw) [[3]])

fmnh.raw = fmnh.raw [,,fmnh.order]

fmnh.export = list ('data' = fmnh.data,
  'raw' = fmnh.raw)

### checando formas

vis.seq (fmnh.raw, start = 33)

save (fmnh.export, file = 'N.FMNH.RData')

### MUSEU PARAENSE EMILIO GOELDI

MPEG = get.skull (DIR = "beproc", lmA = lmA, lmZ = lmZ)

### AMNH remove ind w/ missing lms
mpeg.miss.A = apply (MPEG$A, 3, find.miss.lm)
mpeg.miss.Z = apply (MPEG$Z, 3, find.miss.lm)

mpeg.rem = mpeg.miss.A & mpeg.miss.Z

### compare svd with ols
mpeg.raw = glue.skulls (MPEG$A[,,mpeg.rem], MPEG$Z[,,mpeg.rem], soln = 'svd')

### carregar banco de dados de distâncias

dimnames (mpeg.raw)[[3]] %in% as.character (data.dist [,1] [data.dist$MUSEUM == 'MPEG'])

mpeg.data = data.dist [data.dist$MUSEUM == 'MPEG',1:10]

mpeg.id = mpeg.data$ID

mpeg.match = mpeg.id %in% dimnames (mpeg.raw) [[3]]

mpeg.data = mpeg.data [mpeg.match,]

mpeg.match = dimnames (mpeg.raw) [[3]] %in% mpeg.id

mpeg.raw = mpeg.raw [,,mpeg.match]
mpeg.order = match (mpeg.data$ID, dimnames (mpeg.raw) [[3]])

mpeg.raw = mpeg.raw [,,mpeg.order]

mpeg.export = list ('data' = mpeg.data,
  'raw' = mpeg.raw)

### checando formas

vis.seq (mpeg.raw, start = 33)

save (mpeg.export, file = 'N.MPEG.RData')


### total
rowSums (sapply (list (mpeg.raw, mnrj.raw, fmnh.raw, amnh.raw, usnm.raw, mzusp.raw), dim))

### agora, e os caras sem nome de museu?

sem.museu = data.dist$ID [data.dist$MUSEUM == '']
length (sem.museu)

all.a1.files = dir (pattern = '.A1P', recursive = TRUE, full.names = FALSE, include.dirs = FALSE)

last.elem = function (obj) {return (obj[length(obj)])}
all.a1.files = sapply (strsplit (all.a1.files, split = '/'), last.elem)
all.a1.files = sapply (strsplit (all.a1.files, split = '.', fixed = TRUE), elem1)

sum (as.character (sem.museu) %in% all.a1.files)

all.a1.files [all.a1.files %in% as.character (sem.museu)]

dir (pattern = '.A1P', recursive = TRUE, full.names = FALSE, include.dirs = FALSE) [all.a1.files %in% as.character (sem.museu)]
### JIMGAB/P* e um do MZ

### JIM AMNH


AMNH2 = get.skull.vl (DIR = "JIMAMNH", lmA = lmA, lmZ = lmZ)

### AMNH remove ind w/ missing lms
amnh2.miss.A = apply (AMNH2$A, 3, find.miss.lm)
amnh2.miss.Z = apply (AMNH2$Z, 3, find.miss.lm)

amnh2.rem = amnh2.miss.A & amnh2.miss.Z

### compare svd with ols
amnh2.raw = glue.skulls (AMNH2$A[,,amnh2.rem], AMNH2$Z[,,amnh2.rem], soln = 'svd')

### carregar banco de dados de distâncias

dimnames (amnh2.raw)[[3]] %in% as.character (data.dist [,1] [data.dist$MUSEUM == 'AMNH2'])

amnh2.data = data.dist [data.dist$MUSEUM == 'AMNH2',1:10]

amnh2.id = amnh2.data$ID

amnh2.match = amnh2.id %in% dimnames (amnh2.raw) [[3]]

amnh2.data = amnh2.data [amnh2.match,]

amnh2.match = dimnames (amnh2.raw) [[3]] %in% amnh2.id

amnh2.raw = amnh2.raw [,,amnh2.match]
amnh2.order = match (amnh2.data$ID, dimnames (amnh2.raw) [[3]])

amnh2.raw = amnh2.raw [,,amnh2.order]

amnh2.export = list ('data' = amnh2.data,
  'raw' = amnh2.raw)

### checando formas

vis.seq (amnh2.raw, start = 33)

save (amnh2.export, file = 'N.AMNH2.RData')


ALL = get.skull.mf (lmA = lmA, lmZ = lmZ)
