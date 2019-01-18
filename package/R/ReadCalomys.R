##' @title
##' ReadCalomys
##'
##' @description
##' This function reads Excel files containing the Calomys colony data.
##'
##' @param filename excel register file (xls)
##'
##' @return
##' list with two elements:
##' D: right side
##' E: left side
##'
##' @author Guilherme Garcia
##'
##' @importFrom gdata read.xls
##' @importFrom plyr aaply
##'
##' @export
##' @rdname ReadCalomys
##'
##' @examples
##' \dontrun{
##' calomys.list <-
##'     dir(path = '../Raw Data/Calomys', pattern = 'input',
##'         recursive = TRUE, include.dirs = TRUE, full.names = TRUE)
##'
##' calomys.raw <-
##'     alply(calomys.list, 1, function(f)
##'         {
##'             print(f)
##'             ReadCalomys(f)
##'         })
##'
##' calomys.list <-
##'     calomys.list [!is.na(calomys.raw)]
##'
##' calomys.raw <-
##'     calomys.raw [!is.na(calomys.raw)]
##'
##' calomys.list <- gsub('../Raw Data/Calomys/', '', calomys.list)
##' }

ReadCalomys <- function(filename)
{
    error.found <- FALSE
    ## checking if file is a valid measurement file
    tryCatch(expr = {de.view <- read.xls(filename, sheet = 1, header = FALSE)},
             error = function(cond) {error.found <<- TRUE})

    de.view <- de.view[, 1:16]

    de.start <- which(de.view == 'IS', arr.ind = TRUE)
    de.end <- which(de.view == 'PM', arr.ind = TRUE)

    d.start <- de.start[which(de.start[, 'col'] == 2), ]
    e.start <- de.start[which(de.start[, 'col'] == 10), ]

    d.end <- de.end[which(de.end[, 'col'] == 2), ]
    e.end <- de.end[which(de.end[, 'col'] == 10), ]

    ids <- de.view [d.start - c(1, 1)]

    lms <- as.character(de.view[d.start[1, 'row']:d.end[1, 'row'], 2])
    lms[is.na(lms)] <- 'NA'

    midline <- c('IS', 'NSL', 'NA', 'BR', 'LD', 'BA', 'OPI', 'PNS')
    d.lms <- paste0(lms, ifelse(lms %in% midline, '', '-D'))
    e.lms <- paste0(lms, ifelse(lms %in% midline, '', '-E'))

    d.view <-
        aaply(1:nrow(d.start), 1, function(i)
        {
            start <- d.start[i, ]
            end <- d.end[i, ]
            current <- as.matrix(de.view[start[1]:end[1], start[2]+(1:6)])
            aaply(current, 1, as.numeric)
        })

    e.view <-
        aaply(1:nrow(e.start), 1, function(i)
        {
            start <- e.start[i, ]
            end <- e.end[i, ]
            current <- as.matrix(de.view[start[1]:end[1], start[2]+(1:6)])
            aaply(current, 1, as.numeric)
        })

    dim(d.view) <- dim(e.view) <-
        c(length(ids), length(lms), 3, 2)

    d.view <- aperm(d.view, c(2, 3, 4, 1))
    e.view <- aperm(e.view, c(2, 3, 4, 1))

    dimnames(d.view) <-list(d.lms, c('X', 'Y', 'Z'), c('R1', 'R2'), ids)
    dimnames(e.view) <-list(e.lms, c('X', 'Y', 'Z'), c('R1', 'R2'), ids)

    list('D' = d.view, 'E' = e.view)

}
