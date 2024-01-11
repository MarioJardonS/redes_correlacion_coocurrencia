#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------------------------------------------------

#argumentos: tabla de otus clave, tabla de abundancias, taxonomía, tabla de metadatos, metadato a considerar

args = commandArgs(trailingOnly=TRUE)

#paquetes
library(phyloseq)

#construcción de phyloseqs desde argumentos 

key_otus <- paste0( "../results/central_otus/" , args[1] )
data <- paste0( "../data/tables/" , args[2] )
taxonomy <- paste0( "../data/taxonomy/" , args[3] )
metadata <- paste0("../data/metadata/" , args[4] ) #los argumentos deberían de estar en la carpeta 
#correspondiente
tipo <- args[5] #nombre de la columna de metadatos a considerar

##carga de los argumentos a dataframes de R
key_otus <- read.csv(key_otus , row.names = 1 ) #se asume que es una tabla de salida del script ./first_analysis.R


if (substr( data ,  length(data) - 4 , length(data) - 1 ) == ".csv"){
  data <- read.csv( data , row.names = 1 , header = TRUE )
} else {
  data <- read.table( data , row.names = 1 , header = TRUE , sep = "" )
}

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
o_table <- otu_table(data[ intersect(row.names(data) , row.names(taxonomy)) ,  intersect(row.names(etapa) ,colnames(key_otus)) ], taxa_are_rows = TRUE)


phy_key <- phyloseq(otu_table = o_table_key , sample_data = tipo)
phy <- phyloseq(otu_table = o_table , sample_data = tipo)
#código latex para tabla de OTUs con columnas de especie/género, abundancia media, abundancia mediana y
#error standard
