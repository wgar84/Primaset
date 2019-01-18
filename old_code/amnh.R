require (shapes)

### banco de dados
data.dist = read.csv2 (file = 'GLMALL2.csv', header = TRUE, row.names = NULL)

### landmarks
vis = scan ("vistas.csv", what = "")
lmA = vis[1:35]
lmA[4] = "NA"
lmZ = vis[36:length (vis)]
rm (vis)

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

