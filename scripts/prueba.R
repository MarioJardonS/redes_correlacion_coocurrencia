
args = commandArgs(trailingOnly=TRUE)

if (!require(vegan)) install.packages('vegan')
library(vegan)
if (!require(igraph)) install.packages('igraph')
library(igraph)
#library(ggplot2)
if (!require(apcluster)) install.packages('apcluster')
library(apcluster)
if (!require(plyr)) install.packages('plyr')
library(plyr)


#getwd()
setwd("..")

#####DATOS#####

data <- paste0("./data/tables/", args[1])
data <- read.table(data , row.names = 1, header = TRUE , sep= "" )
print(colnames(data))

metadata <- paste0("./data/metadata/", args[2] )
metadata <- read.csv(metadata , row.names = 1, colClasses = "character")#el 2 me inquieta todavia
r_n_metadata <- metadata[,1] 
for (i in 1:length(r_n_metadata)){
  r_n_metadata[i] <- make.names(r_n_metadata[i])
  }
metadata[,1] <- r_n_metadata
colnames(metadata) <- c("ID","Grupos")
print(metadata)

grupos <- unique(metadata[,"Grupos"])
grouping <- list()
for (i in 1:length(grupos)){
  grouping_i <- c()
  for (j in 1:dim(metadata)[1]){
    if (metadata[j,"Grupos"] == grupos[i]){
      grouping_i <- c(grouping_i , metadata[j,"ID"])
    }
  }
  grouping[[i]] <- grouping_i
}

print(grouping)