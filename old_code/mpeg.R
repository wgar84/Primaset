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

