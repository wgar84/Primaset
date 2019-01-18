GlueSkullsSVD = function (A, Z)
  {
    glue.skull = function (Ai, Zi)
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
            tmp = glue.skull (A[,,i], Z[,,i])
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

GlueSkullsOLS = function (A, Z)
  {
    glue.skull = function (Ai, Zi)
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
        ### solução OLS
        R = solve (solve (t (Ai[cA,]) %*% Ai[cA,]) %*% (t (Ai[cA,]) %*% Zi[cZ,]))
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
            tmp = glue.skull (A[,,i], Z[,,i])
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
