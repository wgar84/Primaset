require (shapes)

base.glm = read.csv2 ('GLMALL2.csv', header = TRUE, row.names = NULL)
base.full = read.csv (file = 'semsp12.csv')

base.full.set = !is.na (apply (base.full[,61:130], 1, sum))
sum (base.full.set)

base.full.id = as.character (base.full$ID.[base.full.set])

elem1 = function (element) {return (element[1])}
base.glm.id = sapply (strsplit (as.character (base.glm$ID), split = '.', fixed = TRUE), elem1)

### elementos na base final que possuem os dois lados intactos
base.common.index = base.glm.id %in% base.full.id

base.common = base.glm.id [base.common.index]
main.data = base.glm [base.common.index, 1:10]

### removendo duplicatas
main.data = main.data [!duplicated (base.common),]
base.common = base.common [!duplicated (base.common)]

length (base.common)


### achando os arquivos
a1p.files = z1p.files = list ()
a1p.count = z1p.count = c()
for (i in 1:length (base.common))
  {
    tmp = dir (pattern = base.common[i], recursive = TRUE)
    tmp2 = grep (pattern = '.A1P', x = tmp)
    a1p.count[i] = length (tmp2) ### registro de arquivos encontrados
    a1p.files[[i]] = tmp [tmp2] ### nomes de arquivo
    tmp2 = grep (pattern = '.Z1P', x = tmp)
    z1p.count[i] = length (tmp2)
    z1p.files[[i]] = tmp [tmp2]
  }

main.data [!files.found,] ### não achou esses

a1p.files.final = sapply (a1p.files, elem1)
z1p.files.final = sapply (z1p.files, elem1)

files.found = a1p.count != 0 & z1p.count != 0

a1p.files.final = a1p.files.final [files.found]
z1p.files.final = z1p.files.final [files.found]
main.id = base.common [files.found]
main = main.data [files.found, ]

main$ID = main.id

### verificando
length (main.id)
length (a1p.files.final)
length (z1p.files.final)
dim (main)
### tudo igual

### generos
table (main$GENUS)

### carregando landmarks
vis = scan ("vistas.csv", what = "")
lmA = vis[1:35]
lmA[4] = "NA"
lmZ = vis[36:length (vis)]
rm (vis)

### pega tudo
main.all = get.skull.mf2 (a1p.files.final, z1p.files.final, lmA, lmZ)

dimnames (main.all$A) [[3]] == main.id ### true

test.glue = glue.skull (main.all$A[,,1], main.all$Z[,,1])

main.glued = glue.skulls (main.all$A, main.all$Z)

### todos com det(R) = 1
main.raw = main.glued[[1]]

main.raw = main.raw * 10

dimnames (main.raw) [[3]] == main.id

save (main.raw, main, file = 'NWM.RData')

vis.seq (main.raw[,,main$GENUS == 'Aotus'], start = 1)
