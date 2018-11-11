#source("functions.R")
#require(here)

require(ecoli2.db)
Ecoli_data <- read.table("datos/Ecoli_data.csv", sep= ",", dec = ".",header = T)
#anno_data <- read.table("datos/Ecoli_ASv2.na36.annot.csv",sep=",",na.strings = "---")
correspondencia <- read.table("datos/correspondencia.csv",sep=",",na.strings = "---", header = T)
