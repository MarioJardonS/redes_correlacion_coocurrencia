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


#Normalización
for (i in 1:dim(data)[1]){
  data[,i] <- data[,i]/sum(data[,i])
}



metadata <- paste0("./data/metadata/", args[2])
metadata <- read.csv(metadata , colClasses = "character")
r_n_metadata <- metadata[,1] 
for (i in 1:length(r_n_metadata)){
  r_n_metadata[i] <- make.names(r_n_metadata[i])
}
metadata[,1] <- r_n_metadata
colnames(metadata) <- c("ID","Grupos")



n_grupos <- unique(metadata[,"Grupos"])
n_grupos <- setdiff( n_grupos , c(NA) )#parche


vector_no_na <- which(is.na(metadata[,"Grupos"]) == FALSE)

metadata <- metadata[vector_no_na,]


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
closeness_cent <- c()
for (i in 1:dim(data)[1]) {
  c_i <- closeness(net_work, data[i,"nodos"])
  closeness_cent <- c(closeness_cent, c_i)
}
data$closeness <- closeness_cent


betweenness_cent <- c()
for (i in 1:dim(data)[1]) {
  b_i <- betweenness(net_work, data[i,"nodos"])
  
  betweenness_cent <- c(betweenness_cent, b_i)
}
data$betweenness <- betweenness_cent


data_deg <- data[order(data$degrees, decreasing = TRUE),]
data_close <- data[order(data$closeness , decreasing = TRUE),]
data_between <- data[order(data$betweenness, decreasing = TRUE),]





file <- args[4]

write.csv(data_deg , paste0("./results/",file,"_bydegree.csv") , row.names = TRUE)
write.csv(data_close,paste0("./results/",file,"_bycloseness.csv") , row.names = TRUE)
write.csv(data_between, paste0("./results/",file,"_bybetweenness.csv") , row.names = TRUE)

report_1 <- args[5]

hdeg <- which(data$degrees >= quantile(data$degrees , probs = seq(0, 1, 0.33))[3])
hclose <- which(data$closeness >= quantile(data$closeness , probs = seq(0, 1, 0.33))[3])
lbetween <- which(data$betweenness <= quantile(data$betweenness , probs = seq(0, 1, 0.33))[2])

results_1 <- intersect(hdeg,hclose)
results_1 <- intersect(results_1 , lbetween)

data_report_1 <- data[results_1,]

write.csv(data_report_1 , paste0("./results/",report_1) , row.names = TRUE)



if (length(grupos) < 2){
  stop("Solo hay un reporte de otus clave")}
## -----------------------------------------------------------------------------------------------------------------------

#Esta función calcula el área bajo la curva de una función representada como un dataframe, donde la columna 1 es el dominio, y sus imágenes están en la columna con nombre feature. 
auc <- function(df, centrality){
  sm_p <- 0
  for (i in 1:dim(df)[1]) {
  f_i = df[i, centrality]
  if (is.nan(f_i) == FALSE)
  {sm_p = sm_p + f_i}
}

return(sm_p)}


## -----------------------------------------------------------------------------------------------------------------------
#Con un dataframe que represente una función, y una cota superior (se espera que sea una fracción del área bajo la curva de dicha funcion). Nos devuelve la cantidad de "otus", como están ordenados según el dataframe, que alcanzan dicha cota.
auc_percent <- function(df, centrality, bound){
  #df es nuestro dataframe de muestras y medidas de centralidad, centrality un string con el nombre de nuestra centralidad, y bound un número positivo 
  
  i <- 1
  sum_par <- 0
  #c <- c() 
  while(i <= dim(df)[1] && sum_par < bound) {
    #c <- c(c,i)
    g_i = df[i, centrality]
    sum_par = sum_par + g_i
    i = i+1
  }
  return(i-1)
}


## -----------------------------------------------------------------------------------------------------------------------
library(plyr)
pseudo_F <- function(df , groups){
  #df es dataframe de muestras, groups, una lista de vectores de etiquetas por grupo , por grupos
  N <- dim(df)[2] #número de muestras
  a <- length(groups) #número de grupos
  
  df_groups <- list() #se crea la lista que incluirá los subdataframes por grupo
  
  for (i in 1:length(groups)){
    df_i <- df[,groups[[i]]]
    df_groups[[i]] <- df_i
  } 
  
  
  
  
  n_s <- llply( .data = df_groups , .fun = ncol )
  n_s <- unlist(n_s)
  n <- mean(n_s) #número promedio de muestras por grupo
  
  dist <- vegdist(t(df))
  dist <- as.vector(dist) #distancias bray_curtis entre todo par de muestras
  
 
  df_groups <- llply(.data = df_groups , .fun = t)
  in_dist <- llply(.data = df_groups , .fun = vegdist )
  in_dist <- llply(.data = in_dist , .fun = as.vector )
  in_dist <- unlist(in_dist) #distancias bray-curtis intra-grupo
  
  #calculo del estadistico pseudo-f
  
  ss_t <- sum(dist^2)/N 
  
  ss_w <- sum(in_dist^2)/n
  
  ss_a <- ss_t - ss_w

  f_stat <- (ss_a/(a-1))/(ss_w/(N-a))
  
  return(f_stat)
}



## -----------------------------------------------------------------------------------------------------------------------
#Se o



## -----------------------------------------------------------------------------------------------------------------------
area_deg <- auc(data_deg , "degrees")
auc5_percent_deg <- c()
for (x in 1:20){
  auc5_percent_deg = c(auc5_percent_deg , auc_percent(data_deg, "degrees" ,(area_deg/20)*x))
  
}


## -----------------------------------------------------------------------------------------------------------------------
area_close <- auc(data_close , "closeness")
auc5_percent_close <- c()
for (x in 1:20){
  auc5_percent_close = c(auc5_percent_close , auc_percent(data_close, "closeness" ,(area_close/20)*x))
  
}


## -----------------------------------------------------------------------------------------------------------------------
area_between <- auc(data_between , "betweenness")
auc5_percent_between <- c()
for (x in 1:20){
  auc5_percent_between = c(auc5_percent_between , auc_percent(data_between, "betweenness" ,(area_between/20)*x))
  
}

analisis_auc <- args[6]


## -----------------------------------------------------------------------------------------------------------------------
f_stat_deg <- c()
var_deg <- c()
for (i in auc5_percent_deg){
  df_i <- data_deg[1:i ,1:length(no_outliers)]
  f_i <- pseudo_F(df_i , grupos)
  f_stat_deg <- c(f_stat_deg , f_i) 
  
  var_i <- var(f_stat_deg)
  var_deg <- c(var_deg, var_i)
}



analisis_auc_deg <- data.frame(auc5_percent_deg, f_stat_deg , var_deg)
write.csv(analisis_auc_deg, paste0("/results/", analisis_auc ,"degree.csv"))

analisis_auc_deg <- analisis_auc_deg[2:dim(analisis_auc_deg)[1],]
n_auc_deg <- which((analisis_auc_deg[,2] - analisis_auc_deg[,3]) == max(analisis_auc_deg[,2] - analisis_auc_deg[,3]))
n_auc_deg <- n_auc_deg + 1
n_auc_deg <- auc5_percent_deg[n_auc_deg]


## -----------------------------------------------------------------------------------------------------------------------
f_stat_close <- c()
var_close <- c()
for (i in auc5_percent_close){
  df_i <- data_close[1:i ,1:length(no_outliers)]
  f_i <- pseudo_F(df_i , grupos)
  f_stat_close <- c(f_stat_close , f_i) 
  
  var_i <- var(f_stat_close)
  var_close <- c(var_close, var_i)
}


analisis_auc_close <- data.frame(auc5_percent_close, f_stat_close , var_close)
write.csv(analisis_auc_deg, paste0("/results/", analisis_auc ,"closeness.csv"))

analisis_auc_close <- analisis_auc_close[2:dim(analisis_auc_close)[1],]
n_auc_close <- which((analisis_auc_close[,2] - analisis_auc_close[,3]) == max(analisis_auc_close[,2] - analisis_auc_close[,3]))
n_auc_close <- n_auc_close + 1
n_auc_close <- auc5_percent_close[n_auc_close]


## -----------------------------------------------------------------------------------------------------------------------
f_stat_between <- c()
var_between <- c()
for (i in auc5_percent_between){
  df_i <- data_between[1:i ,1:length(no_outliers)]
  f_i <- pseudo_F(df_i , grupos)
  f_stat_between <- c(f_stat_between , f_i) 
  
  var_i <- var(f_stat_between)
  var_between <- c(var_between ,var_i)
}

analisis_auc_between <- data.frame(auc5_percent_between, f_stat_between , var_between)
write.csv(analisis_auc_deg, paste0("/results/", analisis_auc ,"betweenness.csv"))

analisis_auc_between <- analisis_auc_between[2:dim(analisis_auc_between)[1],]
n_auc_between <- which((analisis_auc_between[,2] - analisis_auc_between[,3]) == max(analisis_auc_between[,2] - analisis_auc_between[,3]))
n_auc_between <- n_auc_between + 1
n_auc_between <- auc5_percent_between[n_auc_between]


report_2 <- args[7]

results_2_deg <- row.names(data_deg[1:n_auc_deg,])
results_2_close <- row.names(data_close[1:n_auc_close])
results_2_between <- row.names(data_between[1:n_auc_between,])

results_2 <- union(results_2_deg , results_2_close)
results_2 <- union(results_2 , results_2_between)

