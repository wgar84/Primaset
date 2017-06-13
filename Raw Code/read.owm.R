read.owm.files = function ()
  {
    require (XLConnect)
    files = dir (recursive = TRUE, pattern = '.xls')
    files.n = length (files)
    dataA = dataZ = viewA1 = viewA2 = viewZ1 = viewZ2 = NULL
    for (i in 1:files.n)    
      {
        wb = loadWorkbook (files[i])
        if (any (!existsSheet (wb, c('a','z'))))
          {
            cat ('skipped', files[i], '\n')
            next
          }
        ### A
        ws.a = readWorksheet (wb, 'a', endRow = 34, header = FALSE)
        museum = ws.a [1,1]
        ws.a = ws.a [,-1]
        ind = ncol (ws.a) / 8
        if (i == 1)
          lmA = ws.a [-1,2]
        for (j in 1:ind)
          {
            m = 8 * (j - 1)
            if (is.na(ws.a[2,m+3]))
              next
            dataA = rbind (dataA, c(museum, ws.a [1:20, m + 1]))
            ind.now.1 = ws.a [-1, m + (3:5)]
            ind.now.1 = t (ind.now.1)
            dim (ind.now.1) = prod (dim (ind.now.1))
            viewA1 = cbind (viewA1, ind.now.1)
            ind.now.2 = ws.a [-1, m + (6:8)]
            ind.now.2 = t (ind.now.2)
            dim (ind.now.2) = prod (dim (ind.now.2))
            viewA2 = cbind (viewA2, ind.now.2)
          }
        ### Z
        ws.z = readWorksheet (wb, 'z', endRow = 12, header = FALSE)
        museum = ws.z [1,1]
        ws.z = ws.z [,-1]
        ind = ncol (ws.z) / 8
        if (i == 1)
          lmZ = ws.z [,2]
        for (j in 1:ind)
          {
            m = 8 * (j - 1)
            if (is.na(ws.z[1,m+3]))
              next
            dataZ = rbind (dataZ, c(museum, ws.z [1:5, m + 1]))
            ind.now.1 = ws.z [, m + (3:5)]
            ind.now.1 = t (ind.now.1)
            dim (ind.now.1) = prod (dim (ind.now.1))
            viewZ1 = cbind (viewZ1, ind.now.1)
            ind.now.2 = ws.z [, m + (6:8)]
            ind.now.2 = t (ind.now.2)
            dim (ind.now.2) = prod (dim (ind.now.2))
            viewZ2 = cbind (viewZ2, ind.now.2)
          }
        cat ('read', files[i], '\n')
        cat (dim (dataA) [1], dim(dataZ) [1], '\n')
      }
    ## dim (viewA1) = c(3, length (lmA), ncol (viewA1))
    ## viewA1 = aperm (viewA1, c(2,1,3), resize = TRUE)
    ## dim (viewA2) = c(3, length (lmA), ncol (viewA2))
    ## viewA2 = aperm (viewA2, c(2,1,3), resize = TRUE)
    ## dimnames (viewA1) = dimnames (viewA2) =
    ##   list (lmA, c('X','Y','Z'), paste (dataA[,1], dataA[,5], sep = '-'))
    ## dim (viewZ1) = c(3, length (lmZ), ncol (viewZ1))
    ## viewZ1 = aperm (viewZ1, c(2,1,3), resize = TRUE)
    ## dim (viewZ2) = c(3, length (lmZ), ncol (viewZ2))
    ## viewZ2 = aperm (viewZ2, c(2,1,3), resize = TRUE)
    ## dimnames (viewZ1) = dimnames (viewZ2) =
    ##   list (lmZ, c('X','Y','Z'), paste (dataA[,1], dataA[,5], sep = '-'))

    return (list ('A1' = viewA1, 'A2' = viewA2,
                  'Z1' = viewZ1, 'Z2' = viewZ2,
                  'dataA' = dataA, 'dataZ' = dataZ,
                  'lmA' = lmA, 'lmZ' = lmZ))
  }
