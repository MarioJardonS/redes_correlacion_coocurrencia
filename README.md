# Obtención de especies clave desde las redes de correlación o coocurrencia.

Una especie clave es aquella que tiene un efecto desproporcionado en la salud de su ecosistema. Dadas las dificultades en la determinación empírica de las especies clave de las comunidades microbianas se ha optado por buscarlas a través de las redes de correlación o de coocurrencia. En este repositorio se encuentra un análisis de centralidad de OTUs dentro de sus microbiomas con el que se pueden postular candidatos a especie clave. La composición de este repositorio es aproximadamente la siguiente:

* En la carpeta data se encuentran datos de muestras de microbiomas de las rizósferas de maíz, de tomate y de chile, así como sus correspondientes redes de correlación de Spearman, presentes en la subcarpeta networks.
* En la carpeta scripts se encuentra el script de R con el que desde los datos de microbioma de tomate se obtuvieron tablas que ordenan a los OTUs según tres distintas medidas de centralidad. 
* Dichas tablas están en la carpeta results. 
* En la carpeta docs se encuentra un archivo .Rmd que describe los pasos que sigue el script de R.
* En la carpeta example se encuentra un ejemplo de juguete, incluidos datos y archivos .Rmd, .html y .pdf sobre el funcionamiento de estos códigos.
