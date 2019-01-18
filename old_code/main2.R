require (shapes)

### banco de dados
data.dist = read.csv2 (file = 'GLMALL2.csv', header = TRUE, row.names = NULL)

### landmarks
vis = scan ("vistas.csv", what = "")
lmA = vis[1:35]
lmA[4] = "NA"
lmZ = vis[36:length (vis)]
rm (vis)

### pegar todo mundo
ALL = get.skull.mf (lmA = lmA, lmZ = lmZ)

all.miss.A = apply (ALL$A, 3, find.miss.lm)
all.miss.Z = apply (ALL$Z, 3, find.miss.lm)

all.rem = all.miss.A & all.miss.Z
sum (all.rem)

### 136215 zuado (falta MT-E)
all.rem [dimnames (ALL$A) [[3]] == '136215'] = FALSE

### identificações

all.svd = glue.skulls (ALL$A[,,all.rem], ALL$Z[,,all.rem], soln = 'svd')

### crânios esquisitos
all.svd = all.svd [,,-c(1626,1666)]

all.data = data.dist[,1:10]
all.id = all.data$ID

elem1 = function (element) {return (element[1])}
all.id = sapply (strsplit (as.character (all.id), split = '.', fixed = TRUE), elem1)

svd.id = dimnames (all.svd) [[3]]

sum (all.id %in% svd.id)
sum (svd.id %in% all.id)

all.data$ID = all.id

all.data = all.data [all.id %in% svd.id,]

all.data = all.data [!duplicated (all.data$ID),]
all.id = all.data$ID

all.raw = all.svd [,,svd.id %in% all.id]
all.id = dimnames (all.raw) [[3]]

all.raw = all.raw [,,order (all.id)]
all.data = all.data [order (all.data$ID),]

table (all.data$GENUS)

vis.seq (all.raw [,,all.data$GENUS == 'Brachyteles'])

data.id = sapply (strsplit (as.character (data.dist$ID), split = '.', fixed = TRUE), elem1)
ALL.id = dimnames (ALL$A) [[3]]

length (data.id) - sum (ALL.id %in% data.id)
### quase 1200 crânios que não estão aqui (ou não estão processados...)

miss.data = data.dist [!(data.id %in% all.data$ID),1:10]
miss.data$ID = sapply (strsplit (as.character (miss.data$ID), split = '.', fixed = TRUE), elem1)

a1p.my.ass = z1p.my.ass = c()
for (i in miss.data$ID)
  {
    az.file = dir (pattern = i, recursive = TRUE)
    a1p.my.ass = c(a1p.my.ass, az.file [grep (pattern = 'A1P', x = az.file)])
    z1p.my.ass = c(z1p.my.ass, az.file [grep (pattern = 'Z1P', x = az.file)])
  }

REST = get.skull.mf2 (a1p.my.ass, z1p.my.ass, lmA = lmA, lmZ = lmZ)

rest.miss.A = apply (REST$A, 3, find.miss.lm)
rest.miss.Z = apply (REST$Z, 3, find.miss.lm)

rest.rem = rest.miss.A & rest.miss.Z
sum (rest.rem)

### identificações

rest.svd = glue.skulls (REST$A[,,rest.rem], REST$Z[,,rest.rem], soln = 'svd')

### crânios esquisitos
rest.svd = rest.svd [,,-c(44,356)] ### 1a
rest.svd = rest.svd [,,-c(94)] ### 2a

rest.svd.id = dimnames (rest.svd) [[3]]

miss.id = miss.data$ID
miss.data = miss.data [miss.id %in% rest.svd.id,]

miss.data = miss.data [!duplicated (miss.data$ID),]
miss.id = miss.data$ID

rest.raw = rest.svd [,,rest.svd.id %in% miss.id]
rest.id = dimnames (rest.raw) [[3]]

rest.raw = rest.raw [,,order (rest.id)]
miss.data = miss.data [order (miss.data$ID),]

table (miss.data$GENUS)

all.data = rbind(all.data, miss.data)

rest.raw = rest.raw * 10

## 1a
all.raw.new2 = array (0, c (dim (all.raw.new)[1:2], dim (all.raw.new)[3] + dim (rest.raw)[3]))
all.raw.new2[,,1:dim (all.raw.new) [3]] = all.raw.new
all.raw.new2[,,(dim (all.raw.new) [3]+1):dim (all.raw.new2) [3]] = rest.raw

dimnames (all.raw.new2) [[1]] = dimnames (all.raw.new) [[1]]
dimnames (all.raw.new2) [[2]] = dimnames (all.raw.new) [[2]]
dimnames (all.raw.new2) [[3]] = c(dimnames (all.raw.new) [[3]], dimnames (rest.raw) [[3]])

table (all.data$GENUS)

#### cm to mm
# all.raw.new = all.raw.new * 10

two.sides.data = read.csv (file = 'semsp12.csv')

two.sides.full = !is.na (apply (two.sides.data [,61:130], 1, sum))

two.sides.id = as.character (two.sides.data[two.sides.full,1])
two.sides.data [!(two.sides.id %in% all.data$ID),]$SPECIES.

all.data = all.data [!duplicated (all.data$ID),]
all.raw.new2 = all.raw.new2 [,,!duplicated (dimnames (all.raw.new2)[[3]])]

data.id [!(data.id %in% all.data$ID)]


### tentar excluir rotina de verificação de nomes.
