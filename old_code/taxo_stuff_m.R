require (geomorph)
require (shapes)
require (multicore)
load ('ss.RData')
load ('ds.RData')
load ('back2models.RData')
load ('back.again.RData')
source ('../Func/models2.R')
source ('../Func/input.R')
source ('~/Dropbox/Code/matrix.func.r')

options (contrasts = c('contr.sum', 'contr.poly'))

common.data = apply (common.data, 2, as.character)

common.data [, 6] [common.data [, 6] == '0'] = NA
common.data [, 6] [common.data [, 6] == ''] = NA
common.data [, 6] [common.data [, 6] == '?F'] = 'F'
common.data [, 6] [common.data [, 6] == '?M'] = 'M'
common.data [, 6] [common.data [, 6] == 'sexo'] = NA

common.data [, 4] [common.data [, 4] == '0'] = NA
common.data [, 4] [common.data [, 4] == ''] = NA

common.data [, 'SUB'] [common.data [, 'SUB'] == 'aequatoriali'] = 'aequatorialis'
common.data [, 'GEN'] [common.data [, 'GEN'] == 'Leontopithec'] = 'Leontopithecus'
common.data [, 'SUB'] [common.data [, 'SUB'] == 'aequatoriali'] = 'aequatorialis'
common.data [, 'SPE'] [common.data [, 'SPE'] == 'senicula'] = 'seniculus'
common.data [, 'SUB'] [common.data [, 'SUB'] == 'senicula'] = 'seniculus'

common.data

# 1. verificar grupos
# 2. outliers
# 3. ajustar modelos

main.data = list ()
## for (i in 1:length (unique (common.data [, 'GEN'])))
##   {
##     who = unique (common.data [, 'GEN']) [i]
##     main.data [[i]] = createOTU (who, common.data, common.ss, common.ds)
##   }

### empty sp

remove.empty.sp = which (common.data [, 'SPE'] == '')
common.data = common.data [- remove.empty.sp, ]
common.ss = common.ss [, , - remove.empty.sp]
common.ds = common.ds [, , - remove.empty.sp]

# taxonomic corrections (names spelled wrong)
# common.data [, 'GEN'] [common.data [, 'GEN'] == 'Simias'] = 'Nasalis'

common.data [, 'SUB'] [common.data [, 'SUB'] == 'palliates'] = 'palliatus'

common.data [, 'SUB'] [common.data [, 'GEN'] == 'Cacajao' &
                       common.data [, 'SPE'] == 'melanocephal' &
                       common.data [, 'SUB'] == 'melanocephal'] = 'melanocephalus'

common.data [, 'SPE'] [common.data [, 'GEN'] == 'Cacajao' &
                       common.data [, 'SPE'] == 'melanocephal'] = 'melanocephalus'

common.data [, 'SPE'] [common.data [, 'GEN'] == 'Saimiri' &
                       common.data [, 'SPE'] == 'oerstedi'] = 'oerstedii'

common.data [, 'SUB'] [common.data [, 'GEN'] == 'Saimiri' &
                       common.data [, 'SPE'] == 'oerstedii' &
                       common.data [, 'SUB'] == 'oerstedi'] = 'oerstedii'

common.data [, 'SPE'] [common.data [, 'GEN'] == 'Chiropotes' &
                       common.data [, 'SPE'] == 'satanas' &
                       common.data [, 'SUB'] == 'chiropotes'] = 'chiropotes'

common.data [, 'SPE'] [common.data [, 'GEN'] == 'Chiropotes' &
                       common.data [, 'SPE'] == 'satanas' &
                       common.data [, 'SUB'] == 'utahicki'] = 'utahickae'


common.data [, 'GEN'] [common.data [, 'GEN'] == 'Cebuella'] = 'Callithrix'

common.data [, 'GEN'] [common.data [, 'GEN'] == 'Bunopithecus'] = 'Hoolock'

common.data [, 'SPE'] [common.data [, 'GEN'] == 'Pithecia' &
                       common.data [, 'SPE'] == 'monacha'] = 'monachus'

common.data [, 'SPE'] [common.data [, 'GEN'] == 'Alouatta' &
                       common.data [, 'SPE'] == 'fusca'] = 'guariba'

which.to.remove = paste (common.data [, 'ID'],
       common.data [, 'MSM']) %in% paste (info.back [, 'ID'],
                                          info.back [, 'MSM'])

sum (which.to.remove)

which.to.remove2 = ! (paste (common.data [, 'ID'],
    common.data [, 'MSM']) %in% paste (info.back2 [, 'ID'],
                                       info.back2 [, 'MSM']))

sum (which.to.remove & which.to.remove2)

common.data = common.data [which.to.remove & which.to.remove2, ]
common.ss = common.ss [, , which.to.remove & which.to.remove2]
common.ds = common.ds [, , which.to.remove & which.to.remove2]

dim (common.data)
dim (common.ss)
dim (common.ds)

common.data = cbind (common.data, paste (common.data [, 'GEN'], common.data [, 'SPE']))
colnames (common.data) [7] = 'OTU'

main.data = list ()

for (i in 1:length (unique (common.data [, 'OTU'])))
  {
    who = unique (common.data [, 'OTU']) [i]
    main.data [[i]] = list ('info' = common.data [common.data [, 'OTU'] == who, ],
                'ss' = common.ss [, , common.data [, 'OTU'] == who],
                'ds' = common.ds [, , common.data [, 'OTU'] == who])
    main.data [[i]] $ n = nrow (main.data [[i]] $ info)
  }

names (main.data) = unique (common.data [, 'OTU'])
main.data = main.data [which (sapply (main.data, function (x) return (x $ n)) > 20)]

main.data = lapply (main.data, gpagenWrap.ds)
main.data = lapply (main.data, OSymm.wrap)
main.data = lapply (main.data, calcDistWrap.ds)
main.data = lapply (main.data, symPCA.wrap)

for (i in 1:length (main.data))
  {
    main.data [[i]] $ info = data.frame (main.data [[i]] $ info)
  }

names (main.data)

load ('sp.post.RData')
rm (main.sp)

### funcionar
rownames (data.man.sp) = names (main.data)
# data.man.sp = edit (data.man.sp)

for (i in 1:length (main.data))
  {
    main.data [[i]] = adjustModel (main.data [[i]],
                model = data.man.sp [i, 2])
  }

### ajustar modelos pra ed, tirar médias por grupo e arrumar objetos de output
### (assimetria?)

mshape.groups = strsplit (data.man.sp [, 1], split = ' ')

for (i in 1:length (main.data))
  {
    main.data [[i]] = meanShapesGroup (main.data [[i]],
                what = mshape.groups [[i]])
  }

for (i in 1:length (main.data))
  {
    main.data [[i]] = edModel (main.data [[i]],
                model = data.man.sp [i,2])
  }

main.sp = list ()
for (i in 1:length (main.data))
  {
    main.sp [[i]] = fixOTU.postmodel (main.data [[i]],
                model = data.man.sp [i,2])
  }
names (main.sp) = names (main.data)

save (main.sp, data.man.sp, file = 'sp.post.RData')

