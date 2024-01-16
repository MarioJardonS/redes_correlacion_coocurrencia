#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------------------------------------------------

#argumentos: tabla de otus clave, tabla de abundancias, taxonomía, tabla de metadatos, metadato a considerar

args = commandArgs(trailingOnly=TRUE)

#paquetes
library(phyloseq)
library(ggplot2)
library(RColorBrewer)

#construcción de phyloseqs desde argumentos 

key_otus <- paste0( "../results/central_otus/" , args[1] )
data <- paste0( "../data/tables/" , args[2] )
taxonomy <- paste0( "../data/taxonomy/" , args[3] )
metadata <- paste0("../data/metadata/" , args[4] ) #los argumentos deberían de estar en la carpeta 
#correspondiente
tipo <- args[5] #nombre de la columna de metadatos a considerar
nivel <- args[6]

poco_abundante <- args[7]

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




#figuras para abundancia relativa de filo, familia y género

##transformación de phy_key en el dataframe con el que se construirán la figura de abundancias relativas de
##filo, aquí descrito como ta3

if (nivel == "Phylum"){
  nivel0 <- "ta3"
} else {
  if (nivel == "Family"){
    nivel0 <- "ta6"
  } else {
    nivel0 <- "ta7"
  }
}

df_key <- tax_glom(phy_key, taxrank = nivel0)
df_key <- psmelt(df_key)
colnames(df_key)[dim(df_key)[2]] <- nivel

##poco abundante
if (poco_abundante == "Si"){
  abundance <- c()
  for (i in 1:length(levels(df_key[ , nivel]))){
    which_i <- which(as.vector(df_key[ , nivel]) == levels(df_key[ , nivel])[i])
    #print(which_i)
    sum_i <- sum(as.vector(df_key[which_i,]$Abundance))
    abundance <- c(abundance , sum_i)
  }
  nivel_abundante <- which(abundance > summary(abundance)[args[8]] )
  nivel_abundante <- levels(df_key[ , nivel])[nivel_abundante]
  
  ab <- c()
  for (i in 1:dim(df_key)[1]){
    if (is.element(as.vector(df_key[ , nivel])[i] , nivel_abundante )){
      ab <- c(ab , as.vector(df_key[ , nivel])[i]   )
    } else {
      ab <- c(ab , "Other")
    }
  }
  
  df_key[ , nivel] <- ab
}

  

##construcción de la figura de filo 
df_key[ , nivel] <- as.factor(df_key[ , nivel])
colors_rel<- colorRampPalette(brewer.pal(8,"Dark2")) (length(levels(df_key[ , nivel])))


if (nivel == "Phylum"){
  relative_plot <- ggplot(data=df_key, aes(x=Sample, y=Abundance, fill=Phylum))+ 
    geom_bar(aes(), stat="identity", position="stack")+
    scale_fill_manual(values = colors_rel)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
} else {
  if (nivel == "Family"){
    relative_plot <- ggplot(data=df_key, aes(x=Sample, y=Abundance, fill=Family))+ 
      geom_bar(aes(), stat="identity", position="stack")+
      scale_fill_manual(values = colors_rel)+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
  } else {
    relative_plot <- ggplot(data=df_key, aes(x=Sample, y=Abundance, fill=Genus))+ 
      geom_bar(aes(), stat="identity", position="stack")+
      scale_fill_manual(values = colors_rel)+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
  }
}




ggsave( paste0("../results/analisis/" , substr(args[2] , 1 , nchar(args[2])-4 ), "_relative_abundance_" , nivel  , ".png") , relative_plot , device = 'png' )




