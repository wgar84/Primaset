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

