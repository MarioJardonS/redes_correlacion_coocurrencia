
if (!require(vegan)) install.packages('vegan')
library(vegan)
if (!require(igraph)) install.packages('igraph')
library(igraph)
#library(ggplot2)
if (!require(apcluster)) install.packages('apcluster')
library(apcluster)
if (!require(plyr)) install.packages('plyr')
library(plyr)
if (!require(stringr)) install.packages('stringr')
library(stringr)


setwd("/home/mario/proyecto_redes/redes_correlacion_coocurrencia/")

#####DATOS#####
.
data <- paste0("./data/tables/tomate_desarrollo.csv")

if (str_sub( data , -4 , 1) == ".csv" ){
  data <- read.csv(data , row.names = 1 , header = TRUE)
} else {
  data <- read.table(data , row.names = 1, header = TRUE , sep = "" )  
}



