#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------------------------------------------------

#argumentos: tabla de otus clave, tabla de abundancias, taxonomía, tabla de metadatos, metadato a considerar

args = commandArgs(trailingOnly=TRUE)

#paquetes
library(phyloseq)
library(ggplot2)
library(RColorBrewer)



#argumentos generales a todos los phyloseqs, comparten taxonomía, metadatos, y metadato elegido como tipo
taxonomy <- paste0( "../data/taxonomy/" , args[1] )
metadata <- paste0("../data/metadata/" , args[2] ) #los argumentos deberían de estar en la carpeta 
#correspondiente
tipo <- args[3] #nombre de la columna de metadatos a considerar

#nivel taxonómico al que se hará figura
nivel <- args[4]
#criterio para separar otus que se graficaran por separado y los demás
poco_abundante <- args[5]

#lista de tablas de otus a comparar´
lista <- list(args[6:length(args)])

##cada uno de los elementos de la lista se 
lista_data <- list()
for (i in 1:length(lista)){
  
  data <- read.csv( paste0("../data/tables/" , lista[[i]]) , row.names = 1 , header = TRUE )
  
  col <- c()
  for (j in 1:dim(data)[2]) {
    colj <- make.names(colnames(data)[j])
    colj <- substr(colj  , 1 ,  nchar(colj)-21)
    #print(colj)
    col <- c(col , colj)
    
  }
  
  colnames(data) <- col
  lista_data[[i]] <- data
  
}


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

##utilización de tipo para sacar samp_data desde metadata

###los nombres de las muestras están en una columna llamada "ID"
###los nombres de las muestras se corrigen según la sintaxis de R
for (i in 1:dim(metadata)[1]){
  metadata[i , "ID"] <- make.names(metadata[i , "ID"])
}

###se restringen los metadatos a las muestras usadas en el análisis de OTUs clave
metadata <- metadata[which(is.element(metadata[ , "ID"], colnames(tabla) )) ,  ]

tipo <- data.frame( ID = metadata[ , "ID"], Type = metadata [ , tipo ],row.names = metadata[ , "ID"])

tipo <- sample_data(tipo)

##lista en que se guardarán los phyloseq
lista_phy <- list()

for (i in 1:length(lista)){
  o_table_i <- otu_table(lista_data[[i]][ intersect(row.names(lista_data[[i]]) , row.names(taxonomy)) ,  intersect(row.names(tipo) ,colnames(lista_data[[i]])) ], taxa_are_rows = TRUE)
  
  
  taxonomy_i <- taxonomy[ intersect(row.names(taxonomy) ,  row.names(lista_data[[i]]) ) , ]
  taxonomy_table <- tax_table(taxonomy_i)
  row.names(taxonomy_table@.Data) <- row.names(taxonomy_i)
  
}

lista_df <- list()

for (i in 1:length(lista)){
  
  df_key <- tax_glom(lista_phy[[i]], taxrank = nivel0)
  df_key <- psmelt(df_key)
  colnames(df_key)[dim(df_key)[2]] <- nivel
  lista_df[[i]] <- df_key
}

df <- lista_df[[1]]
for (i in 1:(length(lista))){
  df <- merge(lista_df[[i]] , df ,  all = TRUE)
}


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




ggsave( paste0("../results/analisis/" , args[2], "_relative_abundance_" , nivel  , ".png") , relative_plot , device = 'png' )




