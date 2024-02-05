#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------------------------------------------------

#argumentos: tabla de otus clave, tabla de abundancias, taxonomía, tabla de metadatos, metadato a considerar

args = commandArgs(trailingOnly=TRUE)

#paquetes
library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(patchwork)


#argumentos generales a todos los phyloseqs, comparten taxonomía, metadatos, y metadato elegido como tipo
taxonomy <- paste0( "../data/taxonomy/" , args[1] )
metadata <- paste0("../data/metadata/" , args[2] ) #los argumentos deberían de estar en la carpeta 
#correspondiente
tipo <- args[3] #nombre de la columna de metadatos a considerar

#nivel taxonómico al que se hará figura
nivel <- args[4]
#criterio para separar otus que se graficaran por separado y los demás




#lista de tablas de otus a comparar´

lista <- list()
for (i in 5:length(args)){
  lista[[i-4]] <- args[i]
}



##cada uno de los elementos de la lista se 
lista_data <- list()
for (i in 1:length(lista)){
  
  data <- read.csv( paste0("../results/central_otus/" , lista[[i]]) , row.names = 1 , header = TRUE )
  
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

taxonomy <- read.csv(taxonomy , header = FALSE , sep = ";" , row.names = 1)

metadata <- read.csv(metadata ,  colClasses = "character")


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
#metadata <- metadata[which(is.element(metadata[ , "ID"], colnames(tabla) )) ,  ]

#tipo <- data.frame( ID = metadata[ , "ID"], Type = metadata [ , tipo ],row.names = metadata[ , "ID"])

#tipo <- sample_data(tipo)

##lista en que se guardarán los phyloseq
lista_phy <- list()

for (i in 1:length(lista)){
  
  o_table_i <- otu_table(lista_data[[i]][ intersect(row.names(lista_data[[i]]) , row.names(taxonomy)) , intersect( colnames(lista_data[[i]]) , metadata[ , "ID"] )  ], taxa_are_rows = TRUE)
  
  
  taxonomy_i <- taxonomy[ intersect(row.names(taxonomy) ,  row.names(lista_data[[i]]) ) , ]
  taxonomy_table <- tax_table(taxonomy_i)
  row.names(taxonomy_table@.Data) <- row.names(taxonomy_i)
  phy_i <- phyloseq(otu_table = o_table_i , tax_table = taxonomy_table)
  lista_phy[[i]] <- transform_sample_counts(phy_i , function(x) x / sum(x) )
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

  abundance <- c()
  #tax_niv <- c()
  for (i in 1:length(levels(df[ , nivel]))){
    which_i <- which(as.vector(df[ , nivel]) == levels(df[ , nivel])[i])
    #print(which_i)
    sum_i <- sum(as.vector(df[which_i,]$Abundance))
    abundance <- c(abundance , sum_i)
    #tax_niv <- c(tax_niv , levels(df[ , nivel])[i])
  }
  
  if (nivel == "Phylum"){
    il <- 0
    #il <- summary(abundance)["1st Qu."] 
  } else {
    il <- 0
     # il <- summary(abundance)["1st Qu."] 
    
  }
  print(il)
  
  nivel_abundante <- which(abundance >= il )
  nivel_abundante <- levels(df_key[ , nivel])[nivel_abundante]
  #print(nivel_abundante)


  for (i in 1:length(lista)){
    ab_i <- c()
    for (j in 1:dim(lista_df[[i]])[1]){
      if (is.element(as.vector(lista_df[[i]][ , nivel])[j] , nivel_abundante )){
        ab_i <- c(ab_i , as.vector(lista_df[[i]][ , nivel])[j]   )
      } else {
        ab_i <- c(ab_i , "Other")
      }
    }
    
    lista_df[[i]][ , nivel] <- ab_i
    lista_df[[i]][ , nivel] <- as.factor(lista_df[[i]][ , nivel ])
  }
  
  ab <- c()
  for (j in 1:dim(df)[1]){
    if (is.element(as.vector(df[ , nivel])[j] , nivel_abundante )){
      ab <- c(ab , as.vector(df[ , nivel])[j]   )
    } else {
      ab <- c(ab , "Other")
    }
  }
  df[ , nivel] <- ab
    
  df[ , nivel] <- as.factor(df[ , nivel ])
      
  
    



##construcción de la figura de filo 
#df_key[ , nivel] <- as.factor(df_key[ , nivel])
colors_rel <- colorRampPalette(brewer.pal(8,"Dark2")) (length(levels(df[ , nivel])))
names(colors_rel) <- levels(df[ , nivel])
print(colors_rel)

lista_rp <- list()
nombres <- c()
if (nivel == "Phylum"){
  for (i in 1:length(lista)){
  rp_i <- paste0("rp_" , as.character(i))  
  nombres <- c(nombres , rp_i)
  
  relative_plot_i <- ggplot( lista_df[[i]] ,  aes(x=Sample, y=Abundance, fill=Phylum))  +
      geom_bar(aes(), stat="identity", position="stack")+
      scale_fill_manual(values = colors_rel[levels(df[ , nivel])] , drop = FALSE)+ 
      theme(axis.text.x = element_text(angle = 45, hjust = 1) )
  lista_rp[[i]] <- relative_plot_i

  }

  names(lista_rp) <- nombres  
  impar <- which(1:length(lista) %% 2 == 1)
  
  
  for (i in impar){
    if (is.element(i+1 , 1:length(lista))){
    plot <- ggplot_add(lista_rp[[i+1]]  , lista_rp[[i]] , names(lista_rp)[i+1])
    plot <- plot + plot_layout(ncol = 1) 
    } else {
      plot <- lista_rp[[i]]
    }
    
    ggsave( paste0("../results/analisis/",args[4+i] , "_" , args[4+i+1] , "_relative_abundance_" , nivel  , ".png") , plot , device = 'png' )
  }
    
  
  
  
} else {
  if (nivel == "Family"){
    for (i in 1:length(lista)){
      rp_i <- paste0("rp_" , as.character(i))  
      nombres <- c(nombres , rp_i)
      
      relative_plot_i <- ggplot( lista_df[[i]] ,  aes(x=Sample, y=Abundance, fill=Family))  +
        geom_bar(aes(), stat="identity", position="stack")+
        scale_fill_manual(values = colors_rel[levels(df[ , nivel])] , drop = FALSE)+ 
        theme(axis.text.x = element_text(angle = 45, hjust = 1) )
      lista_rp[[i]] <- relative_plot_i
      }
    
    
    names(lista_rp) <- nombres  
    
    impar <- which(1:length(lista) %% 2 == 1)
    
    
    for (i in impar){
      if (is.element(i+1 , 1:length(lista))){
        plot <- ggplot_add(lista_rp[[i+1]]  , lista_rp[[i]] , names(lista_rp)[i+1])
        plot <- plot + plot_layout(ncol = 1) 
      } else {
        plot <- lista_rp[[i]]
      }
      
      ggsave( paste0("../results/analisis/",args[4+i] , "_" , args[4+i+1] , "_relative_abundance_" , nivel  , ".png") , plot , device = 'png' )
    }
    
    
    
    
  } else {
    for (i in 1:length(lista)){
      rp_i <- paste0("rp_" , as.character(i))  
      nombres <- c(nombres , rp_i)
      
      relative_plot_i <- ggplot( lista_df[[i]] ,  aes(x=Sample, y=Abundance, fill=Genus))  +
        geom_bar(aes(), stat="identity", position="stack")+
        scale_fill_manual(values = colors_rel , drop = FALSE)+ 
        theme(axis.text.x = element_text(angle = 45, hjust = 1) )
      lista_rp[[i]] <- relative_plot_i
    }
    
    
    names(lista_rp) <- nombres  
    impar <- which(1:length(lista) %% 2 == 1)
    
    
    for (i in impar){
      if (is.element(i+1 , 1:length(lista))){
        plot <- ggplot_add(lista_rp[[i+1]]  , lista_rp[[i]] , names(lista_rp)[i+1])
        plot <- plot + plot_layout(ncol = 1) 
      } else {
        plot <- lista_rp[[i]]
      }
      
      ggsave( paste0("../results/analisis/",args[4+i] , "_" , args[4+i+1] , "_relative_abundance_" , nivel  , ".png") , plot , device = 'png' )
    }
    
  }
}


#ggsave( paste0("../results/analisis/", args[5], args[6], "_relative_abundance_" , nivel  , ".png") , plot , device = 'png' )

#ggsave( paste0("../results/analisis/" ,as.character(length(lista)), "_relative_abundance_" , nivel  , ".png") , relative_plot , device = 'png' )




