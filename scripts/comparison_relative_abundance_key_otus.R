#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------------------------------------------------

#argumentos: tabla de otus clave, tabla de abundancias, taxonomía, tabla de metadatos, metadato a considerar

args = commandArgs(trailingOnly=TRUE)

#paquetes
library(phyloseq)
library(ggplot2)
library(RColorBrewer)

#construcción de phyloseqs desde argumentos 

lista <- args[1]
taxonomy <- paste0( "../data/taxonomy/" , args[2] )
metadata <- paste0("../data/metadata/" , args[4] ) #los argumentos deberían de estar en la carpeta 
#correspondiente
tipo <- args[5] #nombre de la columna de metadatos a considerar
nivel <- args[6]

poco_abundante <- args[7]






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




ggsave( paste0("../results/analisis/" , args[2], "_relative_abundance_" , nivel  , ".png") , relative_plot , device = 'png' )




