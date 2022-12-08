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
for (i in 1:(length(n_grupos)-1)){
  grupos_i <- c()
  for (j in 1:dim(metadata)[1]){
    if (metadata[j,"Grupos"] == n_grupos[i]){
      grupos_i <- c(grupos_i , metadata[j,"ID"])
    }
  }
  grupos[[i]] <- grupos_i
}



#Análisis pcoa y/o clusterización para detectar outliers

bc_dist <- vegdist(t(data), method = "bray")
#PCoA <- cmdscale(bc_dist, eig = TRUE, k = 2)
s <- 1 - bc_dist
s <- as.matrix(s)
clustering <- apcluster(s)

clusters <- clustering@clusters
filtro_0 <- lapply(clusters, length)
clusters_no_outliers <- clusters[which(filtro_0 > 1)]
no_outliers <- unlist(clusters_no_outliers)
no_outliers <- names(no_outliers)
data <- data[,no_outliers]

#Agrupación sin outliers

for (i in 1:length(grupos)){
  grupos[[i]] <- intersect(grupos[[i]], no_outliers)
}

#Grupos necesarios para el análisis auc
len_list <- llply(grupos , length)
len_list <- which(len_list > 1)
grupos <- grupos[len_list]


#######ANÁLISIS DE OTUS######

data$nodos <- 0:(dim(data)[1]-1)

#Eliminación de otus según su aparición en muestras
filt <- c()
for (i in 1:dim(data)[1]) {
  
  v_i <- as.vector(data[i,1:(dim(data)[2]-1)])
  #el siguiente 1 es filtro
  if (length(v_i [ v_i > 0 ]) > 1 ) {
    filt <- c(filt, i)
  }
}

data <- data[filt,]


######CARGA DE RED Y AJUSTE A FILTRACIÓN DE OTUS######
red <- paste0("./data/networks/", args[3])
red <- read.csv(red)
red = red[,1:3]#se asume la forma del archivo de red

#Dado que se han filtrado otus, solo retendremos las aristas que se refieren a los otus conservados en nuestros datos
edges <- c()
for (i in 1:dim(red)[1]) {
  if (is.element(red[i,1], data$nodos) && is.element(red[i,2], data$nodos) && red[i,3] > 0  ){
    edges <- c(edges , i)
  }
}

red <- red[edges, 1:2]

#####AJUSTES PREVIOS AL ISO DE igraph######
#red <- red + 1

for (i in 1:dim(red)[1]){
  for (j in 1:dim(red)[2]){
    red[i,j] <- paste0("v_",as.character(red[i,j]))
  }
}


#data$nodos <- data$nodos + 1

for (i in 1:dim(data)[1]){
  data[i,"nodos"] <- paste0("v_" , as.character(data[i,"nodos"]))
}


########CARGA DE RED CON igraph Y ELECCION DE COMPONENTE CONEXA PRINCIPAL##### 


net_work <- graph_from_edgelist(as.matrix(red) , directed = FALSE )



##componente(s) conexa(s) principal(es)
compo_conexas <- components(net_work)
size_compo_conexas <- compo_conexas$csize
princ <- which(size_compo_conexas == max(size_compo_conexas))
pertenencia <- compo_conexas$membership
compo_princ <- which(pertenencia == princ )
compo_princ <- names(compo_princ)

##nuevos datos

filtro_componente <- c()
for (i in 1:dim(data)[1]){
  if(is.element(data[i,"nodos"],compo_princ)){
    filtro_componente <- c(filtro_componente, i)
  }
}


data <- data[filtro_componente,]
print(dim(data))
net_work <- induced_subgraph(net_work, compo_princ ,"auto")

degrees <- c()
for (i in 1:dim(data)[1]) {
  d_i <- degree(net_work, data[i,"nodos"])
  degrees <- c(degrees, d_i)
}
data$degrees <- degrees


## -----------------------------------------------------------------------------------------------------------------------
#closeness_cent <- c()
#for (i in 1:dim(data)[1]) {
#  c_i <- closeness(net_work, data[i,"nodos"])
#  closeness_cent <- c(closeness_cent, c_i)
#}
#data$closeness <- closeness_cent


#betweenness_cent <- c()
#for (i in 1:dim(data)[1]) {
#  b_i <- betweenness(net_work, data[i,"nodos"])
#  print(b_i)
#  betweenness_cent <- c(betweenness_cent, b_i)
#}
#data$betweenness <- betweenness_cent


data_deg <- data[order(data$degrees, decreasing = TRUE),]
#data_close <- data[order(data$closeness , decreasing = TRUE),]
#data_between <- data[order(data$betweenness, decreasing = TRUE),]





file <- args[4]

write.csv(data_deg , paste0("./results/",file,"_bydegree.csv") , row.names = TRUE)
#write.csv(data_close,paste0("./results/",file,"_bycloseness.csv") , row.names = TRUE)
#write.csv(data_between, paste0("./results/",file,"_bybetweenness.csv") , row.names = TRUE)
######CALCULO DE MEDIDAS DE CENTRALIDAD############


#file <- args[4]

#write.csv(data_deg , paste0("./results/",file,"_bydegree.csv") , row.names = TRUE)
#write.csv(data_close,paste0("./results/",file,"_bycloseness.csv") , row.names = TRUE)
#write.csv(data_between, paste0("./results/",file,"_bydegree.csv") , row.names = TRUE)