% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ReadProsimian.R
\name{ReadProsimian}
\alias{ReadProsimian}
\title{ReadProsimian}
\usage{
ReadProsimian(filename)
}
\arguments{
\item{filename}{prosimian register file (xlsx format)}
}
\value{
List with two elements:
id: data frame with information
shapes: array of coordinate data.
}
\description{
This function reads Excel files with prosimian coordinate data.
}
\details{
ZMB_Berlin/Propithecus/Propithecus_Berlin_D.xlsx has more individuals listed in info than
actual individuals measured; had to delete extra entries to read that file.

FNHM_Chicago/Out_Madagascar/Loris/Loris_A.xlsx also had the same issue;
furthermore, localities were in the wrong column.

FNHM_Chicago/Out_Madagascar/Tarsius/Tarsius_B.xlsx and
FNHM_Chicago/Out_Madagascar/Tarsius/Tarsius_C.xlsx also had this issue

FPD_Duke Lemur Center/Madagascar/Extant/Prolemur/Prolemur_A.xlsx
had wrong tag in first lm of third specimen

MNHN_Paris/Madagascar/Mirza/Mirza_A.xlsx had some weird
characters and extra columns in 'Info ID'; had to remove extra columns

we glued a skull in AMNH_NewYork/Madagascar/Cheirogaleus/Cheirogaleus_A_ok.xlsx 

all xlsx files in MNHN_Paris/Madagascar/Propithecus/ also required a similar intervention;
delete extra columns in 'Info ID'
}
\examples{
\dontrun{
## you can change the path to collection here
prosimian.list <-
   dir(path = '../Raw Data/Prosimians/Collections',
       pattern = '.xlsx', recursive = TRUE, include.dirs = TRUE, full.names = TRUE)

## if you solve the issues listed in Details, this should read the entire collection
prosimian.raw <-
   alply(prosimian.list, 1, function(f)
       {
           print(f)
           ReadProsimian(f)
       })

## this removes whatever other xlsx files there were on the list (marked with NA)
prosimian.list <-
    prosimian.list [!is.na(prosimian.raw)]

prosimian.raw <-
   prosimian.raw [!is.na(prosimian.raw)]

## names each element in the list with its file name, removing the path

prosimian.list <- gsub('../Raw Data/Prosimians/Collections/', '', prosimian.list)

names(prosimian.raw) <- prosimian.list

## this saves the current state of the DB
save(prosimian.raw, file = '01_from_files.RData')
}
}
\author{
Guilherme Garcia
}
