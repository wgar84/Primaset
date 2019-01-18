require (ape)
require (phytools)
require (geiger)
require (plyr)

load ('springer2012.RData')

load ('out.of.this.mess.RData')
load ('outofhere.RData')
load ('who.RData')

matrices = c (nwm.mat, owm.mat)
sex = c (as.character (nwm.info$SEX), as.character (owm.info$SEX))
sex = factor (sex)
# levels (sex) = c ('F', 'M', 'X')
taxa = c(nwm.taxa, owm.taxa)

rm (nwm.MAT, owm.MAT) ### safely in outofhere.RData
length (matrices)

mfd <- array (0, c(39, 3, 61))
mfd [, , 1:29] <- nwm.vec 
mfd [, , 30:61] <- owm.vec

ed = rbind (nwm.ed, owm.ed)

dimnames (mfd) [[1]] <- dimnames (nwm.vec) [[1]]
dimnames (mfd) [[2]] <- dimnames (nwm.vec) [[2]]
dimnames (mfd) [[3]] <- c (dimnames (nwm.vec) [[3]], dimnames (owm.vec) [[3]])

info <- list ('nwm' = nwm.info, 'owm' = owm.info)

rm (list = ls (pat = 'nwm'))
rm (list = ls (pat = 'owm'))

### arrancar

### Homo_sapiens_negro
### Saguinus_fuscicollis_nigrifrons
### Saguinus_fuscicollis_weddelli
### Saguinus_midas_midas
### Chiropotes albinasus
### Cebus albifrons
### Cebus nigritus
### Cercopithecus denti

rip.it.off <- c('Homo sapiens negro', 'Saguinus fuscicollis nigrifrons',
                'Saguinus fuscicollis weddelli', 'Saguinus midas midas',
                'Chiropotes albinasus', 'Cebus albifrons', 'Cebus nigritus',
                'Cercopithecus denti')

rip.who <- taxa %in% rip.it.off
rip.uni <- names (matrices) %in% rip.it.off

ed <- ed [!rip.who, ]
mfd <- mfd [, , !rip.uni]
sex <- sex [!rip.who]
matrices <- matrices [!rip.uni]
taxa <- taxa [!rip.who]

### mudar nomes

### mudar nomes

### Alouatta_fusca -> Alouatta_guariba
### Alouatta_senicula -> Alouatta_seniculus
### Aotus_lemurinus_lemurinus -> Aotus_lemurinus
### Cebuella_pygmaea -> Callithrix_pygmaea
### Cercopithecus_mitis_stuhlmanni -> Cercopithecus_mitis
### Chiropotes_satanas_chiropotes -> Chiropotes_chiropotes
### Chlorocebus_pygerythrus_hilgerti -> Chlorocebus_pygerythrus
### Gorilla_beringei_graueri -> Gorilla_beringei
### Homo_sapiens_branco -> Homo_sapiens
### Macaca_fascicularis_fascicularis -> Macaca_fascicularis
### Pithecia_irrorata_irrorata -> Pithecia_irrorata
### Pithecia_monacha -> Pithecia_monachus
### Saguinus_fuscicollis_illigeri -> Saguinus_fuscicollis
### Saguinus_midas_niger -> Saguinus_midas
### Saguinus_mystax_mystax -> Saguinus_mystax

change.names <- c ('Alouatta fusca', 'Alouatta guariba',
                   'Alouatta senicula', 'Alouatta seniculus',
                   'Aotus lemurinus lemurinus', 'Aotus lemurinus',
                   'Cebuella pygmaea', 'Callithrix pygmaea',
                   'Cercopithecus mitis stuhlmanni', 'Cercopithecus mitis',
                   'Chiropotes satanas chiropotes', 'Chiropotes chiropotes',
                   'Chlorocebus pygerythrus hilgerti', 'Chlorocebus pygerythrus',
                   'Gorilla beringei graueri', 'Gorilla beringei',
                   'Homo sapiens branco', 'Homo sapiens',
                   'Macaca fascicularis fascicularis', 'Macaca fascicularis',
                   'Pithecia irrorata irrorata', 'Pithecia irrorata',
                   'Pithecia monacha', 'Pithecia monachus',
                   'Saguinus fuscicollis illigeri', 'Saguinus fuscicollis',
                   'Saguinus midas niger', 'Saguinus midas',
                   'Saguinus mystax mystax', 'Saguinus mystax')

dim (change.names) <- c(2, 15)
change.names <- t (change.names)

for (i in 1:dim (change.names) [1])
    {
        change.who <- taxa %in% change.names [i, 1]
        change.uni <- names (matrices) %in% change.names [i, 1]
        dimnames (mfd) [[3]] [change.uni] <- change.names [i, 2]
        names (matrices) [change.uni] <- change.names [i, 2]
        taxa [change.who] <- rep (change.names [i, 2], sum (change.who))
    }

names (matrices) = gsub (' ', '_', names (matrices))
dimnames (mfd) [[3]] = gsub (' ', '_', dimnames (mfd) [[3]])
taxa = gsub (' ', '_', taxa)

names (sex) = taxa
rownames (ed) = taxa

rm (list = ls(pat = 'rip'))
rm (list = ls(pat = 'change'))

trees = llply (Springer2012, function (x) return (treedata (x, matrices) $ phy))

matrices = treedata (trees [[1]], matrices, sort = TRUE) $ data
taxa.used = dimnames (matrices) [[1]]
dim (matrices) = NULL
names (matrices) = taxa.used

order (names (matrices))
mfd = aperm (mfd, c(3, 1, 2), resize = TRUE)
names.mfd = dimnames (mfd)

mfd [1, , ]

dim (mfd) = c (53, 39 * 3)

dimnames (mfd) [[1]] = names.mfd [[1]]

mfd = treedata (trees [[1]], mfd, sort = TRUE) $ data
names.mfd [[1]] = dimnames (mfd) [[1]]
dim (mfd) = c(53, 39, 3)
dimnames (mfd) = names.mfd
mfd = aperm (mfd, c(2, 3, 1), resize = TRUE)

dimnames (mfd)
