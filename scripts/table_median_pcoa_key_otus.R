#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------------------------------------------------

#argumentos: tabla de otus clave, tabla de abundancias, taxonomía, tabla de metadatos, metadato a considerar

args = commandArgs(trailingOnly=TRUE)

#paquetes
library(phyloseq)
library(xtable)
library(ggplot2)

#construcción de phyloseqs desde argumentos 

key_otus <- paste0( "../results/central_otus/" , args[1] )
data <- paste0( "../data/tables/" , args[2] )
taxonomy <- paste0( "../data/taxonomy/" , args[3] )
metadata <- paste0("../data/metadata/" , args[4] ) #los argumentos deberían de estar en la carpeta 
#correspondiente
tipo <- args[5] #nombre de la columna de metadatos a considerar

##carga de los argumentos a dataframes de R
key_otus <- read.csv(key_otus , row.names = 1 ) #se asume que es una tabla de salida del script ./first_analysis.R

#arreglar nombres de tabla de otus clave

coln <- c()
for (j in 1:(dim(key_otus)[2]-4)) {
  col_j <- make.names(colnames(key_otus)[j])
  col_j <- substr(col_j  , 1 ,  nchar(col_j)-21)
  
  coln <- c(coln , col_j )
}
colnames(key_otus) <- c(coln , colnames(key_otus)[(length(colnames(key_otus))-3):length(colnames(key_otus))])




#if (substr( data ,  length(data) - 3 , length(data)  ) == ".csv"){
data <- read.csv( data , row.names = 1 , header = TRUE )
#} else {
  #data <- read.table( data , row.names = 1 , header = TRUE , sep = "" )
#}


#arreglar nombres de tabla de muestras
col <- c()
for (j in 1:dim(data)[2]) {
  colj <- make.names(colnames(data)[j])
  colj <- substr(colj  , 1 ,  nchar(colj)-21)
  #print(colj)
  col <- c(col , colj)
  
}

colnames(data) <- col



taxonomy <- read.csv(taxonomy , header = FALSE , sep = ";" , row.names = 1)

metadata <- read.csv(metadata ,  colClasses = "character")




##utilización de tipo para sacar samp_data desde metadata

###los nombres de las muestras están en una columna llamada "ID"
###los nombres de las muestras se corrigen según la sintaxis de R
for (i in 1:dim(metadata)[1]){
  metadata[i , "ID"] <- make.names(metadata[i , "ID"])
}

###se restringen los metadatos a las muestras usadas en el análisis de OTUs clave
metadata <- metadata[which(is.element(metadata[ , "ID"], colnames(key_otus) )) ,  ]

tipo <- data.frame( ID = metadata[ , "ID"], Type = metadata [ , tipo ],row.names = metadata[ , "ID"])

tipo <- sample_data(tipo)

##obtención de tax_table y otu_table desde key_otus, data y taxonomy
o_table_key <- otu_table(key_otus[intersect(row.names(key_otus) , row.names(taxonomy)) ,  intersect(row.names(tipo) , colnames(key_otus))  ] , taxa_are_rows = TRUE)  


o_table <- otu_table(data[ intersect(row.names(data) , row.names(taxonomy)) ,  intersect(row.names(tipo) ,colnames(key_otus)) ], taxa_are_rows = TRUE)


taxonomy <- taxonomy[ intersect(row.names(taxonomy) ,  row.names(key_otus) ) , ]
taxonomy_table <- tax_table(taxonomy)
row.names(taxonomy_table@.Data) <- row.names(taxonomy)

##creación y normalización de los phyloseqs
phy_key <- phyloseq(otu_table = o_table_key , tax_table = taxonomy_table , sample_data = tipo)
#print(phy_key@tax_table)
phy <- phyloseq(otu_table = o_table , sample_data = tipo)

phy_key <- transform_sample_counts(phy_key , function(x) x / sum(x) )
phy <- transform_sample_counts(phy , function(x) x / sum(x) )


#código latex para tabla de OTUs con columnas de especie/género, abundancia media, abundancia mediana y
#error standard

##vectores que serán columnas de la tablas
otu <- c()
medias <- c()
medianas <- c()
se <- c()

for (i in row.names(taxonomy)){
  #print(i)
  #de momento la taxonomía es "V9" = especie y "V8" = género
  if (taxonomy[i, "V9"] != ""){
    otu <- c(otu , as.character(taxonomy[i, "V9"]))
    #print(taxonomy[i, "V9"])
  } else {
    otu <- c(otu ,as.character(taxonomy[i, "V8"]))
   # print(as.character(taxonomy[i, "V8"]))
  }
  
  medias <- c(medias , mean(phy@otu_table@.Data[i , ] ))
  medianas <- c(medianas , median(phy@otu_table@.Data[i , ] ))
  se <- c(se , mean_se(phy@otu_table@.Data[i , ])[ ,"ymax" ]- mean_se(phy@otu_table@.Data[i , ])[ ,"y" ] )
}

tabla <- data.frame( "OTU" = otu ,  "MeanRA" = medias , "MedianRA" = medianas, "SE" = se  , row.names = row.names(taxonomy) )

print(xtable(tabla , caption = "Keystone OTUs of "  , digits = 8))

key_not_key <- c()
for (i in 1:dim(phy@otu_table@.Data)[1]){
  
  if (is.element( row.names(phy@otu_table@.Data)[i]   , row.names(phy_key@otu_table@.Data)  )   ){
    key_not_key <- c(key_not_key , "key")
  } else {
    key_not_key <- c(key_not_key , "not_key")
    
  }
}

#figuras que contrastan abundancia de otus clave con media y mediana através de las muestras

##código para una figura que contrasta distribuciones individuales con las medianas
keys_vs_median <- phy@otu_table@.Data[which(key_not_key == "key") , ]
medians <- c()
means <- c()
for (i in 1:dim(phy@otu_table@.Data)[2]){
  
  medians <- c(medians , median(phy@otu_table@.Data[, i]))
  
}

keys_vs_median <- rbind(medians , keys_vs_median)
#keys_vs_median <- rbind(keys_vs_median , means)
keys_vs_median <-  psmelt(otu_table(keys_vs_median  , taxa_are_rows = TRUE))

plot <- ggplot(keys_vs_median , aes(x = Sample , y = Abundance , col = OTU , group = OTU))  + geom_line( ) + facet_wrap(~OTU)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave( paste0("../results/analisis/abundance_" , args[2] , "_key_otus_medians.png" ), plot ,device = "png")




##codigo para obtener una figura que contrasta la distribución de todos los OTUs con la de los OTUs clave
key_vs_no_key <-  psmelt(phy@otu_table)

key_not_key <- c()
for (i in 1:dim(key_vs_no_key)[1]){
  
  if (is.element( key_vs_no_key[i, "OTU"]   , row.names(phy_key@otu_table@.Data)  )   ){
    key_not_key <- c(key_not_key , "Keystone")
  } else {
    key_not_key <- c(key_not_key , "Not keystone")
    
  }
}
key_vs_no_key <- cbind(key_vs_no_key , key_not_key)
colnames(key_vs_no_key)[4] <- c("Type")

plot <- ggplot(key_vs_no_key , aes(x = Sample , y= Abundance , color = Type  , group = Type ))+ #geom_dotplot(binaxis='y', stackdir='center') +
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), 
               geom="crossbar", width=0.5) +
  stat_summary(fun="median", fun.args = list(mult=1), 
               geom="line", width=0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0("../results/analisis/mean_median_key_vs_not_key_" , args[2] , ".png") , plot , device = "png")

#Análisis PCoA con distancia bray-curtis

meta_ord <- ordinate( phy , method = "PCoA", distance = "bray")
plot_pcoa_muestras <- plot_ordination(physeq = phy, ordination = meta_ord , color = "Type")
ggsave( paste0("../results/analisis/pcoa_muestras_" , args[2] , ".png"), plot_pcoa_muestras , device = "png")

meta_ord <- ordinate( phy_key , method = "PCoA", distance = "bray")
plot_pcoa_key <- plot_ordination(physeq = phy_key, ordination = meta_ord , color = "Type")
ggsave( paste0("../results/analisis/pcoa_key_otus_" , args[2] , ".png"), plot_pcoa_key , device = "png")