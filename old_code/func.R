get.skull.mf2 = function (A.file, Z.file, lmA, lmZ)
  {
    indA = length (A.file)
    indZ = length (Z.file)
    labA = labZ = NULL
    A = Z = NULL
    for (i in 1:indA)
      {
        ### pega vista A
        cat (A.file[i], '\n')
        tmpA = scan (A.file[i], what = "")
        labtA = tmpA[1]
        tmpA = as.numeric (tmpA[-1])
        if (length (tmpA) %% length (lmA) != 0)
          tmpA = c (tmpA, rep (NA, times = 6))
        labA = c (labA, labtA)
        A = cbind (A, tmpA)
      }
    for (i in 1:indZ)
      {
        ### pega vista Z
        cat (Z.file[i], '\n')
        tmpZ = scan (Z.file[i], what = "")
        labtZ = tmpZ[1]
        tmpZ = as.numeric (tmpZ[-1])
        if (length (tmpZ) %% length (lmZ) != 0)
          cat ("Não adicionou", Z.file[i],". Tamanho não conforma.\n")
        else
          {
            labZ = c(labZ, labtZ)
            Z = cbind (Z, tmpZ)
          }
      }
    dim (A) = c (3, length (lmA), length (labA))
    dim (Z) = c (3, length (lmZ), length (labZ))
    A = aperm (A, c(2,1,3), resize = TRUE)
    Z = aperm (Z, c(2,1,3), resize = TRUE)
    ### removendo duplicatas
    dupA = duplicated (labA)
    dupZ = duplicated (labZ)
    labA = labA [!dupA]
    labZ = labA [!dupZ]
    A = A [,,!dupA]
    Z = Z [,,!dupZ]
    dimnames (A) [[1]] = lmA
    dimnames (Z) [[1]] = lmZ
    dimnames (A) [[2]] = dimnames (Z) [[2]] = c("X","Y","Z")
    dimnames (A) [[3]] = labA
    dimnames (Z) [[3]] = labZ
    ### remover NLT
    A = A [-length (lmA) + c(0,1),,]
    return (list ("A" = A, "Z" = Z))
  }

glue.skulls = function (A, Z, soln = 'svd')
  {
    glue.skull = function (Ai, Zi, sol = "svd")
      {
        mA = dim (Ai) [1]
        mZ = dim (Zi) [1]
        lmA = dimnames (Ai) [[1]]
        lmZ = dimnames (Zi) [[1]]
        cA = which (lmA %in% lmZ)
        cZ = which (lmZ %in% lmA)
        ### ordenar landmarks
        cA = cA[order (lmA[cA])]
        cZ = cZ[order (lmZ[cZ])]
        ### centroides dos pontos em comum
        ccA = t (array (colMeans (Ai[cA,]), c(3, mA)))
        ccZ = t (array (colMeans (Zi[cZ,]), c(3, mZ)))
        ### centralizando no centroide dos pontos em comum
        Ai = Ai - ccA
        Zi = Zi - ccZ
        ### angulos entre pontos em comum
        M = t (Ai[cA,]) %*% Zi[cZ,]
        UDV = svd (M)
        ### R é ortonormal, mas faz reflexões quando M faz
        D = diag (ifelse (UDV$d > 0, 1, -1))
        R = UDV$v %*% D %*% t(UDV$u)
        ### segunda rotação (sem expansão)
        Zi = Zi %*% R
        out = list (rbind (Ai, Zi [!(lmZ %in% lmA),]), det (R))
        return (out)
      }
    if (dim (A) [3] != dim (Z) [3])
      {
        cat ('Número de vistas A não bate com número de vistas Z. Verifique!','\n')
        return (-1)
      }
    else
      {
        out = list ()
        dets = c()
        for (i in 1:dim (A)[3])
          {
            tmp = glue.skull (A[,,i], Z[,,i], sol = soln)
            if (i == 1)
              out = array (0, c(dim (tmp[[1]]), dim (A)[3]))
            out[,,i] = tmp[[1]]
            dets[i] = tmp[[2]] 
          }
        dimnames (out) = list (rownames (tmp[[1]]),
                   colnames (tmp[[1]]),
                   dimnames (A) [[3]])
        right.first = function (element)
          {
            return (ifelse (length (element) == 1, 0,
                            ifelse (element [1] == 'NLT', 3,
                                    ifelse (element [2] == 'D', 1, 2))))
          }
        ord = strsplit (dimnames (out) [[1]], split = '-')
        ord = sapply (ord, right.first)
        out = out [order (ord),,]
        return (list (out, dets))
      }
  }

print.skull = function (skull, double = TRUE, rgl.open = TRUE, spit.wire = FALSE)
  {
    require (shapes)
    x = 27
    lines = t (array (dim = c(2,x),
      data = c (
        'IS', 'PM',
        'NA', 'NSL',
        'NA', 'BR',
        'PM', 'MT',
        'MT', 'ZI',
        'ZS', 'ZI',
        'BR', 'PT',
        'FM', 'ZS',
        'FM', 'NA',
        'ZS', 'NSL',
        'PT', 'TSP',
        'IS', 'PNS',
        'ZI', 'ZYGO',
        'BR', 'LD',
        'FM', 'PT',
        'BA', 'PNS',
        'LD', 'AS',
        'TS', 'ZYGO',
        'AS', 'TSP',
        'EAM', 'PEAM',
        'AS', 'PEAM',
        'JP', 'AS',
        'JP', 'BA',
        'LD', 'OPI',
        'JP', 'APET',
        'AS', 'BR',
        'PM', 'NSL'
        )))
    if (double)
      {
        lines.d = lines.e = lines
        lines.d = ifelse (lines.d == 'IS' | lines.d == 'NA' | lines.d == 'NSL' | lines.d == 'PNS' |
          lines.d == 'BR' | lines.d == 'BA' | lines.d == 'OPI' | lines.d == 'LD',
          lines.d, paste (lines.d, '-D', sep = ''))
        lines.e = ifelse (lines.e == 'IS' | lines.e == 'NA' | lines.e == 'NSL' | lines.e == 'PNS' |
          lines.e == 'BR' | lines.e == 'BA' | lines.e == 'OPI' | lines.e == 'LD',
          lines.e, paste (lines.e, '-E', sep = ''))
        lines = rbind (lines.d, lines.e)
      }
    on.skull = rownames (skull)
    links = array (match (lines, on.skull), dim (lines))
    if (spit.wire)
      return (links)
    shapes3d (skull, rglopen = rgl.open)
    for (i in 1:dim (links)[1])
      {
        if (!any (is.na(links[i,])))
          shapes3d (skull[links[i,],], type = 'l', joinline = 1:2, rglopen = FALSE)
      }
    
  }

to.morphologika <- function (mtarray,wire,poly) {
    rp <- mtarray
    ind <- dim(rp)[3]
    if (length(dim(rp))==4)
      ind <- ind*dim(rp)[4]
    dm <- dim(rp)[2]
    lms <- dim(rp)[1]
    nom <- dimnames(rp)[[3]]
    if (length(dim(rp))==4)
      nom <- paste(rep(dimnames(rp)[[4]],each=dim(rp)[3]),rep(nom,times=dim(rp)[4]),sep="-")
    rp <- array(rp,c(lms,dm,ind))
    dimnames(rp)[[3]] <- paste("'#",1:ind,sep="")
    output <- list("{individuals}"=ind,
                   "{landmarks}"=lms,
                   "{dimensions}"=dm,
                   "{names}"=paste(nom,c("PL"),sep=""),
                   "{rawpoints}"=rp,
                   "{wireframe}"=wire,
                   "{polygons}"=poly)
    return(output) }

vis.seq = function (sk, start = 1)
  {
    ind = dim (sk) [3]
    for (i in start:ind)
      {
        print.skull (sk[,,i], rgl.open = ifelse (i == 1, TRUE, FALSE))
        tmp = readline(prompt = paste (i, ': ', dimnames (sk) [[3]] [i], ' (\'int\' to leave) > '))
        if (tmp == 'int')
          break
        rgl.clear ()
      }
  }

 (setq last-kbd-macro
    [?\M-< ?\C-  ?\M-> ?\C-\M-% ?\\ ?\[ ?. ?* ?\\ ?\] return return ?! ?\M-< ?\C-  ?\M-> ?\C-\M-% ?\C-o right ?  ?+ return ?\C-o return ?! ?\M-< ?\C-  ?\M-> ?\M-% ?, ?  ?, ?  return return ?! ?\M-< ?\C-  ?\M-> ?\M-% ?\C-o ?\C-o return ?\C-o return ?! ?\M-< ?\C-  ?\M-> ?\C-\M-% ?` return return ?! ?\M-< ?\C-  ?\M-> ?\M-% ?\{ return ?\[ return ?! ?\M-< ?\C-  ?\M-> ?\M-% ?\} return ?\] return ?! ?\M-< ?\C-  ?\M-> ?\M-% ?$ return return ?! ?\M-< ?\C-  ?\M-> ?\C-\M-% ?  ?  ?+ return ?  return ?! ?\M-< ?\C-  ?\M-> ?\M-% ?\" return return ?! ?\M-< ?\C-  ?\M-> ?\M-% ?P ?L ?\S-  return ?\C-o return ?! ?\C-  ?\M-< ?\C-  ?\M-> ?\M-% ?P ?L return return ?! ?\M-< ?\C-  ?\M-> ?\M-% ?\C-o ?\C-o return ?\C-o return ?!])
