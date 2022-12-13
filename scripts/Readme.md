El script aquí presente nos da dos reportes de candidatos a OTUs clave desde muestras de microbioma. Sus primeros tres argumentos consisten de: 
1. Una tabla de abundancias (en este caso se han usado los resultados de la función biom.table desde archivos biom)
2. Un archivo .csv consistente de una columna consistente de las muestras de la tabla anterior y una columna de metadatos 
3. Un archivo .csv cuyas primeras dos entradas representen pares de OTUs (sin redundancias) y la tercera sea una medida de correlación entre los OTUs correspondientes.

Estos tres archivos deben estar respectivamente en las carpetas ./data/tables/ , ./data/metadata/ y ./data/networks/ , paralelas a la carpeta ./scripts/ donde se encuentre este script. Los siguientes cuatro argumentos son los nombres de los resultados:
1. El nombre de las tablas en las que la mayor componente conexa de la red estará ordenada según medidas de centralidad
2. El nombre del primer reporte (intersección de grado alto, centralidad alta e intermediación baja)
3. El nombre de las tablas con pseudo F-estadísticos según cortes de 5% de las medidas de centralidad
4. El nombre del segundo reporte (OTUs centrales que maximicen la diferencia entre las muestras).
Estos archivos .csv aparecerán en la carpeta ./results. 
