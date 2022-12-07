#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------------------------------------------------


#####DATOS#####

#data <- paste0("./data/", args[1])
#data <- read.table(data , row.names = 1, header = FALSE , sep= "" )

#####ANALISIS Y AGRUPACION DE LAS MUESTRAS####

#Agrupación de las muestras según metadatos
#produccion <- c("V2", "V3" ,"V4","V5")
#llenado_de_fruto <- c("V6", "V7", "V8")
#plantacion <- c()
#por_transplantar <- c("V9")
#desarrollo <- c("V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17","V18","V19","V20","V21","V22","V23")


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


metadata <- paste0("./data/metadata/", args[2] )
metadata <- read.csv(metadata , row.names = 1, colClasses = "character")#el formato de estos metadatos me inquieta
r_n_metadata <- metadata[,1] 
for (i in 1:length(r_n_metadata)){
  r_n_metadata[i] <- make.names(r_n_metadata[i])
}
metadata[,1] <- r_n_metadata
colnames(metadata) <- c("ID","Grupos")



n_grupos <- unique(metadata[,"Grupos"])

print(n_grupos)
grupos <- list()
for (i in 1:(length(n_grupos))){
  grupos_i <- c()
  for (j in 1:(dim(metadata)[1]-1)){
    if (metadata[j,"Grupos"] == n_grupos[i]){
      grupos_i <- c(grupos_i , metadata[j,"ID"])
    }
    }
  print(grupos_i)
  grupos[[i]] <- grupos_i
}

print(grupos)