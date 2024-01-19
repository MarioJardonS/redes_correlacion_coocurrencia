#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------------------------------------------------

#construcción de phyloseqs 
#los argumentos son una tabla de muestras (totales o parciales) una taxonomía posiblemnte con más otus 
#que la tabla, unos matadatos posiblemente con más muestras, y un tipo que tiene que estar en los metadatos
#esta funcion está para aplicarse en scripts que comparen una tabla de muestras con su reporte de otus clave
#o que comparen varios de estos reportes

phy <- function(tabla , taxonomy , metadata , tipo){
  
  
  
  
  
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
  
  ##obtención de tax_table y otu_table desde key_otus, data y taxonomy
  
  
  
  o_table <- otu_table(tabla[ intersect(row.names(tabla) , row.names(taxonomy)) ,  intersect(row.names(tipo) ,colnames(data)) ], taxa_are_rows = TRUE)
  
  
  taxonomy <- taxonomy[ intersect(row.names(taxonomy) ,  row.names(tabla) ) , ]
  taxonomy_table <- tax_table(taxonomy)
  row.names(taxonomy_table@.Data) <- row.names(taxonomy)
  
  phy_data <- phyloseq(otu_table = o_table , tax_table = taxonomy_table , sample_data = tipo)
  
  phy_data <- transform_sample_counts(phy_data , function(x) x / sum(x) )} return( phy_data )

