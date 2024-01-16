#!/bin/bash

Rscript table_key_otus.R  $1 $2 $3 $4 $5 > ../results/analisis/$1.tex
Rscript figures_relative_abundance_key_otus.R $1 $2 $3 $4 $5 Phylum No   >> ../results/analisis/$1.tex

echo "\begin{figure}" >> ../results/analisis/$1.tex
 echo "\centering" >> ../results/analisis/$1.tex

echo "\includegraphics[scale = 0.8]{"$2"_relative_abundance_Phylum.png}" >> ../results/analisis/$1.tex

echo "\caption{Relative abundance by phyla of keystone OTUs }" >> ../results/analisis/$1.tex

echo "\label{fig:"$2"_phyla}" >> ../results/analisis/$1.tex

echo "\end{figure}" >> ../results/analisis/$1.tex

Rscript figures_relative_abundance_key_otus.R $1 $2 $3 $4 $5 Family Si Median  >> ../results/analisis/$1.tex

echo "\begin{figure}" >> ../results/analisis/$1.tex
 echo "\centering" >> ../results/analisis/$1.tex

echo "\includegraphics[scale = 0.8]{"$2"_relative_abundance_Family.png}" >> ../results/analisis/$1.tex

echo "\caption{Relative abundance by families of keystone OTUs }" >> ../results/analisis/$1.tex

echo "\label{fig:"$2"_family}" >> ../results/analisis/$1.tex

echo "\end{figure}" >> ../results/analisis/$1.tex

Rscript figures_relative_abundance_key_otus.R $1 $2 $3 $4 $5 Genus Si Median  >> ../results/analisis/$1.tex

echo "\begin{figure}" >> ../results/analisis/$1.tex
 echo "\centering" >> ../results/analisis/$1.tex

echo "\includegraphics[scale = 0.8]{"$2"_relative_abundance_Genus.png}" >> ../results/analisis/$1.tex

echo "\caption{Relative abundance by genera of keystone OTUs }" >> ../results/analisis/$1.tex

echo "\label{fig:"$2"_genus}" >> ../results/analisis/$1.tex

echo "\end{figure}" >> ../results/analisis/$1.tex


Rscript figure_mean_median_key_otus.R $1 $2 $3 $4 $5 
