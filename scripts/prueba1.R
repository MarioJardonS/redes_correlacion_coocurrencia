#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------------------------------------------------


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
data <- read.table(data , row.names = 1, header = TRUE , sep = "" )


metadata <- paste0("./data/metadata/", args[2])
metadata <- read.csv(metadata , colClasses = "character")
r_n_metadata <- metadata[,1] 
for (i in 1:length(r_n_metadata)){
  r_n_metadata[i] <- make.names(r_n_metadata[i])
}
metadata[,1] <- r_n_metadata
colnames(metadata) <- c("ID","Grupos")

print(metadata)

n_grupos <- unique(metadata[,"Grupos"])
n_grupos <- setdiff( n_grupos , c(NA) )#parche


vector_no_na <- which(is.na(metadata[,"Grupos"]) == FALSE)

metadata <- metadata[vector_no_na,]

print(metadata)
grupos <- list()
for (i in 1:(length(n_grupos))){
  grupos_i <- c()
  for (j in 1:dim(metadata)[1]){
    if (metadata[j,"Grupos"] == n_grupos[i]){
      grupos_i <- c(grupos_i , metadata[j,"ID"])
    }
  }
  grupos[[i]] <- grupos_i
}

print(grupos)