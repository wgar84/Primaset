##' @title
##' ReadCatarrhine
##'
##' @description
##' This function reads Excel files with catarrhine coordinate data.
##' 
##' @param filename catarrhine register file (xls format)
##' @return
##' list with three elements:
##' id: data.frame with specimen information
##' A: array with A view coordinates
##' Z: array with Z view coordinates
##'
##' @author Guilherme Garcia
##'
##' @importFrom gdata read.xls
##' @importFrom plyr aaply alply
##' @examples
##' \dontrun{
##' catarrhine.list <- 
##'     dir(path = '../Raw Data/Catarrhini',
##'         pattern = '.xls', recursive = TRUE, include.dirs = TRUE, full.names = TRUE)
##' 
##' catarrhine.raw <-
##'     alply(catarrhine.list, 1, function(f)
##'         {
##'             print(f)
##'             ReadCatarrhine(f)
##'         })
##' 
##' catarrhine.list <-
##'     catarrhine.list [!is.na(catarrhine.raw)]
##' 
##' catarrhine.raw <-
##'     catarrhine.raw [!is.na(catarrhine.raw)]
##' 
##' catarrhine.list <- gsub('../Raw Data/Catarrhini/', '', catarrhine.list)
##' 
##' names(catarrhine.raw) <- catarrhine.list
##' 
##' save(catarrhine.raw, file = 'Catarrhini/01_from_files.RData')
##' } 
ReadCatarrhine <- function(filename)
{
    error.found <- FALSE
    ## checking if filename is a valid measurement file
    tryCatch(expr = {a.view <- read.xls(filename, sheet = 'a', header = FALSE)},
             error = function(cond) {error.found <<- TRUE})
    tryCatch(expr = {z.view <- read.xls(filename, sheet = 'z', header = FALSE)},
             error = function(cond) {error.found <<- TRUE})
    if(error.found) return(NA)
    
    a.view <- a.view[1:34, ]
    z.view <- z.view[1:12, ]

    Museum <- a.view [1, 1]

    id.head <-
        c('Genero','Especie','Subespecie','numero','sexo','local','pais',
          'latitude','longitude','medidas','data','coletor','altitude',
          'obs','obs','obs','obs','obs','obs','obs')

    ### A
    a.start <- which(a.view == 'IS', arr.ind = TRUE)
    a.end <- which(a.view == 'JP', arr.ind = TRUE)
    a.end <- a.end[seq(2, nrow(a.end), 2), ]

    a.lm <- as.character(a.view[a.start[1, 'row']:a.end[1, 'row'], a.start[1, 'col']])
    a.lm[is.na(a.lm)] <- 'NA'
    a.left <- duplicated(a.lm)
    a.right <- duplicated(a.lm, fromLast = TRUE)
    a.lm[a.left] <- paste0(a.lm[a.left], '-E')
    a.lm[a.right] <- paste0(a.lm[a.right], '-D')

    ### Z
    z.start <- which(z.view == 'LD', arr.ind = TRUE)
    z.end <- which(z.view == 'PT', arr.ind = TRUE)
    z.end <- z.end[seq(2, nrow(z.end), 2), ]

    z.lm <- as.character(z.view[z.start[1, 'row']:z.end[1, 'row'], z.start[1, 'col']])
    z.left <- duplicated(z.lm)
    z.right <- duplicated(z.lm, fromLast = TRUE)
    z.lm[z.left] <- paste0(z.lm[z.left], '-E')
    z.lm[z.right] <- paste0(z.lm[z.right], '-D')

    ## ID
    a.id <- 
        aaply(1:nrow(a.start), 1, function(i)
        {
            id.start <- a.start [i, ] - c(1, 1) 
            id.end <- 20
            current.id <- as.character(a.view[id.start['row']:id.end, id.start['col']])
            current.id
        })

    z.id <- 
        aaply(1:nrow(z.start), 1, function(i)
        {
            id.start <- z.start [i, ] - c(1, 1) 
            id.end <- 5
            current.id <- as.character(z.view[id.start['row']:id.end, id.start['col']])
            current.id
        })
    
    ## WORK ON THAT ID
    incomplete <-
        suppressWarnings(aaply(a.id, 1, function(l) all(l == id.head)))

    if(any(incomplete))
    {
        a.id <- a.id [1:(which(incomplete) [1]), ]
        z.id <- z.id [1:(which(incomplete) [1]), ]
    }

    obs <- which(id.head == 'obs')
    id.head[obs] <- paste0('obs', 1:length(obs))
    
    colnames(a.id) <- id.head
    colnames(z.id) <- id.head[1:5]

    a.id[, 'numero'] <- paste0(Museum, a.id[, 'numero'])
    z.id[, 'numero'] <- paste0(Museum, z.id[, 'numero'])

    if(!all (a.id[, 'numero'] == z.id[, 'numero']))
        warning('a.id and z.id don\'t match')

    id <- a.id

    a.shapes <- aaply(1:nrow(id), 1, function(i)
    {
        row.frame <- a.start[1, 'row']:a.end[1, 'row']
        col.frame <- (1:6) + a.start[i, 'col']
        current.shape <- as.matrix(a.view[row.frame, col.frame])
        dim(current.shape) <- c(length(a.lm), 3, 2)
        current.shape
    })

    z.shapes <- aaply(1:nrow(id), 1, function(i)
    {
        row.frame <- z.start[1, 'row']:z.end[1, 'row']
        col.frame <- (1:6) + z.start[i, 'col']
        current.shape <- as.matrix(z.view[row.frame, col.frame])
        dim(current.shape) <- c(length(z.lm), 3, 2)
        current.shape
    })

    a.shapes <- aperm(a.shapes, c(2, 3, 4, 1))
    z.shapes <- aperm(z.shapes, c(2, 3, 4, 1))

    dimnames(a.shapes) <-
        list(a.lm, c('X', 'Y', 'Z'), c('R1', 'R2'), id [, 'numero'])

    dimnames(z.shapes) <-
        list(z.lm, c('X', 'Y', 'Z'), c('R1', 'R2'), id [, 'numero'])

    list('id' = as.data.frame(id), 'A' = a.shapes, 'Z' = z.shapes)
}
