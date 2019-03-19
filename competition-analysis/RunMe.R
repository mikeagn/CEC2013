# clean memory
rm(list=ls(all = TRUE)) 

# macro for loading libraries
library <- function(package, ..., character.only=FALSE){
  if (!character.only)
    package <- as.character(substitute(package))
  yesno <- require(package, ..., character.only = TRUE)
  if(!yesno){
    try(install.packages(package, dependencies=TRUE))
    yesno <- require(package, ..., character.only = TRUE)
  }
  invisible(yesno)
} 

suppressWarnings(library(data.table))
suppressWarnings(library(ggplot2))
suppressWarnings(library(reshape2))
suppressWarnings(library(progress))

source("Analysis1-classic.R")
source("Analysis2-staticF1.R")
source("Analysis3-dynF1.R")

print(results.classic)
print(results.staticF1)
print(results.dynF1)
