#' ReadProsimian
#'
#' This function reads Excel files with prosimian coordinate data.
#'
#' @param file prosimian register file (xlsx format)
#' @return
#' List with two elements:
#' id: data frame with information
#' shapes: array of coordinate data.
#'
#' @details
#' ZMB_Berlin/Propithecus/Propithecus_Berlin_D.xlsx has more individuals listed in info than
#' actual individuals measured; had to delete extra entries to read that file.
#'
#' FNHM_Chicago/Out_Madagascar/Loris/Loris_A.xlsx also had the same issue;
#' furthermore, localities were in the wrong column.
#'
#' FNHM_Chicago/Out_Madagascar/Tarsius/Tarsius_B.xlsx and
#' FNHM_Chicago/Out_Madagascar/Tarsius/Tarsius_C.xlsx also had this issue
#'
#' FPD_Duke Lemur Center/Madagascar/Extant/Prolemur/Prolemur_A.xlsx
#' had wrong tag in first lm of third specimen
#'
#' MNHN_Paris/Madagascar/Mirza/Mirza_A.xlsx had some weird
#' characters and extra columns in 'Info ID'; had to remove extra columns
#'
#' we glued a skull in AMNH_NewYork/Madagascar/Cheirogaleus/Cheirogaleus_A_ok.xlsx 
#' 
#' all xlsx files in MNHN_Paris/Madagascar/Propithecus/ also required a similar intervention;
#' delete extra columns in 'Info ID'
#'  
#' @author Guilherme Garcia
#' @rdname ReadProsimian
#' @export
#' @importFrom gdata read.xls
#' @importFrom plyr aaply alply
#' @examples
#' \dontrun{
#' ## you can change the path to collection here
#' prosimian.list <-
#'    dir(path = '../Raw Data/Prosimians/Collections',
#'        pattern = '.xlsx', recursive = TRUE, include.dirs = TRUE, full.names = TRUE)
#'
#' ## if you solve the issues listed in Details, this should read the entire collection
#' prosimian.raw <-
#'    alply(prosimian.list, 1, function(f)
#'        {
#'            print(f)
#'            ReadProsimian(f)
#'        })
#'
#' ## this removes whatever other xlsx files there were on the list (marked with NA)
#' prosimian.list <-
#'     prosimian.list [!is.na(prosimian.raw)]
#'
#' prosimian.raw <-
#'    prosimian.raw [!is.na(prosimian.raw)]
#'
#' ## names each element in the list with its file name, removing the path
#'
#' prosimian.list <- gsub('../Raw Data/Prosimians/Collections/', '', prosimian.list)
#'
#' names(prosimian.raw) <- prosimian.list
#' 
#' ## this saves the current state of the DB
#' save(prosimian.raw, file = '01_from_files.RData')
#' }

ReadProsimian <- function(file)
{
    error.found <- FALSE
    ## checking if file is a valid measurement file
    tryCatch(expr = {id <- read.xls(file, sheet = 'Info ID')},
             error = function(cond) {error.found <<- TRUE})
    tryCatch(expr = {raw.data <- read.xls(file, sheet = 'All-in-one')},
             error = function(cond) {error.found <<- TRUE})
    if(error.found) return(NA)

    ## catch info and trim
    id <- as.matrix(id[, 1:11])
    id <- id[!is.na(id [, 1]), ]
    id <- id [, -1]
    indiv <- paste0(id [, 2], id [, 1])
    indiv <- gsub(' ', '', indiv)
    
    id <- cbind(indiv, id)
    
    colnames(id) <- c('Ind', 'Tombo', 'Museu', 'Genero', 'Especie', 'Sexo',
                      'Localidade', 'DataColeta', 'DataDado', 'Obs1', 'Obs2')
    
    ## trim raw.data
    raw.data <- raw.data [, 1:9]
    
    ## define frame
    starts <- which(raw.data == 'IS', arr.ind = TRUE)
    starts <- starts[seq(1, nrow(starts), 2), ]
    ends <- which(raw.data == 'JPd', arr.ind = TRUE)

    lm.names <- as.character(raw.data[starts[1, 'row']:ends[1, 'row'], starts[1, 'col']])
    ## NA SUCKS BIGTIME
    lm.names[is.na(lm.names)] <- 'NA'
    lm.names <- gsub('e', '-E', lm.names)
    lm.names <- gsub('d', '-D', lm.names)
    lm.names <- toupper(lm.names)
    shapes <- aaply(1:nrow(id), 1, function(i)
    {
        row.frame <- starts[i, 'row']:ends[i, 'row']
        col.frame <- (1:6) + starts[1, 'col']
        current.shape <- as.matrix(raw.data[row.frame, col.frame])
        dim(current.shape) <- c(length(lm.names), 3, 2)
        current.shape
    })

    shapes <- aperm(shapes, c(2, 3, 4, 1))
    
    dimnames(shapes) <-
        list(lm.names, c('X', 'Y', 'Z'), c('R1', 'R2'), indiv)

    list('id' = id, 'shapes' = shapes)
}    
