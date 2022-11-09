## -----------------------------------------------------------------------------------------------------------------------
getwd()
setwd("~/redes_correlacion_coocurrencia")

#Se cargan los datos de las muestras
data <- read.table("table.from_tomate.txt", row.names = 1, header = FALSE , sep= "" )

#Se excluyeron 3 outliers según un análisis pcoa con diversidad bray-curtis
data <- data[,c(1,3:6,9:21)]

#Se agrega una columna que a cada otu asigna la etiqueta correspondiente en la red
data$nodos <- 0:(dim(data)[1]-1)


#Dado que hay varios otus solo presentes en una muestra, y tienen por lo tanto grado artificialmente alto, dichos otus son descartados. Esta filtración puede modificarse para descartar otus presentes en a lo más otra cota de muestras

#Se crea el vector que escogerá los otus en más de una 
filt <- c()
for (i in 1:dim(data)[1]) {
  #nos concentramos en las columnas referentes a las muestras
  v_i <- as.vector(data[i,1:7])
  #el siguiente 1 es filtro
  if (length(v_i [ v_i > 0 ]) > 1 ) {
    filt <- c(filt, i)
  }
}

data <- data[filt,]

head(data)
dim(data)


## -----------------------------------------------------------------------------------------------------------------------
#La separación es por etapa fenológica; los siguientes vectores describen qué muestras corresponden a cada etapa
produccion <- c("V2","V4","V5")
llenado_de_fruto <- c("V6", "V7")
#plantacion <- c()
#por_transplantar <- c("V9")
desarrollo <- c("V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17","V18","V19","V20","V21","V22")

#no_des <- c()

grupos <- list()
grupos[[1]] <- produccion
grupos[[2]] <- llenado_de_fruto
grupos[[3]] <- desarrollo




## -----------------------------------------------------------------------------------------------------------------------
#Se carga la red
red <- read.csv("networks/tomate_species_raw_network.csv")
red = red[,1:3]

#Dado que se han filtrado otus, solo retendremos las aristas que se refieren a los otus conservados en nuestros datos
edges <- c()
for (i in 1:dim(red)[1]) {
  if (is.element(red[i,1], data$nodos) && is.element(red[i,2], data$nodos) && red[i,3] > 0  ){
    edges <- c(edges , i)
  }
}

red <- red[edges, 1:2]
red <- red + 1

data$nodos <- data$nodos + 1

head(red)
dim(red)


## -----------------------------------------------------------------------------------------------------------------------
library(igraph)
red <- graph_from_edgelist(as.matrix(red) , directed = FALSE )
red


## -----------------------------------------------------------------------------------------------------------------------
degrees <- c()
for (i in 1:dim(data)[1]) {
  d_i <- degree(red, data[i,"nodos"])
  degrees <- c(degrees, d_i)
}
data$degrees <- degrees


## -----------------------------------------------------------------------------------------------------------------------
closeness_cent <- c()
for (i in 1:dim(data)[1]) {
  c_i <- closeness(red, data[i,"nodos"])
  closeness_cent <- c(closeness_cent, c_i)
}
data$closeness <- closeness_cent


## -----------------------------------------------------------------------------------------------------------------------
betweenness_cent <- c()
for (i in 1:dim(data)[1]) {
  b_i <- betweenness(red, data[i,"nodos"])
  betweenness_cent <- c(betweenness_cent, b_i)
}
data$betweenness <- betweenness_cent


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
library(vegan)
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
data_deg <- data[order(data$degrees, decreasing = TRUE),]
data_close <- data[order(data$closeness , decreasing = TRUE),]
#data_between <- data[order(data$betweenness, decreasing = TRUE),]


## -----------------------------------------------------------------------------------------------------------------------
area_deg <- auc(data_deg , "degrees")
auc5_percent_deg <- c()
for (x in 1:20){
  auc5_percent_deg = c(auc5_percent_deg , auc_percent(data_deg, "degrees" ,(area_deg/20)*x))
  print(auc_percent(data_deg, "degrees", (area_deg/20)*x))
}


## -----------------------------------------------------------------------------------------------------------------------
area_close <- auc(data_close , "closeness")
auc5_percent_close <- c()
for (x in 1:20){
  auc5_percent_close = c(auc5_percent_close , auc_percent(data_close, "closeness" ,(area_close/20)*x))
  print(auc_percent(data_close, "closeness" , (area_close/20)*x))
}


## -----------------------------------------------------------------------------------------------------------------------
#area_between <- auc(data_between , "betweenness")
#auc5_percent_between <- c()
#for (x in 1:20){
 # auc5_percent_between = c(auc5_percent_between , auc_percent(data_between, "betweenness" ,(area/20)*x))
  #print(auc_percent(data_between, "betweenness" , (area/20)*x))
#}


## -----------------------------------------------------------------------------------------------------------------------
f_stat_deg <- c()
for (i in auc5_percent_deg){
  df_i <- data_deg[1:i ,1:18]
  f_i <- pseudo_F(df_i , grupos)
  f_stat_deg <- c(f_stat_deg , f_i) 
  print(c(f_i ,var(f_stat_deg)))
}


## -----------------------------------------------------------------------------------------------------------------------
f_stat_close <- c()
for (i in auc5_percent_close){
  df_i <- data_close[1:i ,1:18]
  f_i <- pseudo_F(df_i , grupos)
  f_stat_close <- c(f_stat_close , f_i) 
  print(c(f_i ,var(f_stat_close)))
}


## -----------------------------------------------------------------------------------------------------------------------
f_stat_between <- c()
for (i in auc5_percent_between){
  df_i <- data_between[1:i ,1:18]
  f_i <- pseudo_F(df_i , grupos)
  f_stat_between <- c(f_stat_between , f_i) 
  print(c(f_i ,var(f_stat_between)))
}

