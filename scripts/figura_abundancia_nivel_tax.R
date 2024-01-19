#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------------------------------------------------

#figuras que comparen distribuciones de varias tablas a un nivel taxonómico
#para scripts que comparen distintos reportes


figura <- function( lista , nivel ){
  lista_df <- list()
  for (i in 1:length(lista)){
    lista_df[[i]] <- tax_glom(lista[[i]] , nivel)
  }
  
  lista_df <- lapply(lista_df , psmelt)
  
  
}


