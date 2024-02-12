#!/bin/bash

Rscript table_median_pcoa_key_otus.R  $1 $2 $3 $4 $5 > ../results/analisis/$1.tex
Rscript figures_relative_abundance_key_otus.R $1 $2 $3 $4 $5 Phylum No   >> ../results/analisis/$1.tex
Rscript figures_relative_abundance_key_otus.R $1 $2 $3 $4 $5 Family Si Median  >> ../results/analisis/$1.tex
Rscript figures_relative_abundance_key_otus.R $1 $2 $3 $4 $5 Genus Si Median  >> ../results/analisis/$1.tex

#figura phylum
echo "\begin{figure}" >> ../results/analisis/$1.tex
echo "\centering" >> ../results/analisis/$1.tex
echo "\includegraphics[scale = 0.8]{"$2"_relative_abundance_Phylum.png}" >> ../results/analisis/$1.tex
echo "\caption{Relative abundance by phyla of keystone OTUs }" >> ../results/analisis/$1.tex
echo "\label{fig:"$2"_phyla}" >> ../results/analisis/$1.tex
echo "\end{figure}" >> ../results/analisis/$1.tex

#figura familia
echo "\begin{figure}" >> ../results/analisis/$1.tex
echo "\centering" >> ../results/analisis/$1.tex
echo "\includegraphics[scale = 0.8]{"$2"_relative_abundance_Family.png}" >> ../results/analisis/$1.tex
echo "\caption{Relative abundance by families of keystone OTUs }" >> ../results/analisis/$1.tex
echo "\label{fig:"$2"_family}" >> ../results/analisis/$1.tex
echo "\end{figure}" >> ../results/analisis/$1.tex

#figura genero
echo "\begin{figure}" >> ../results/analisis/$1.tex
echo "\centering" >> ../results/analisis/$1.tex
echo "\includegraphics[scale = 0.8]{"$2"_relative_abundance_Genus.png}" >> ../results/analisis/$1.tex
echo "\caption{Relative abundance by genera of keystone OTUs }" >> ../results/analisis/$1.tex
echo "\label{fig:"$2"_genus}" >> ../results/analisis/$1.tex
echo "\end{figure}" >> ../results/analisis/$1.tex

#figura key indivduales vs mediana
echo "\begin{figure}" >> ../results/analisis/$1.tex

echo "   \centering" >> ../results/analisis/$1.tex

echo "   \includegraphics[scale = 0.8]{abundance_"$2"_key_otus_medians.png}" >> ../results/analisis/$1.tex

echo "   \caption{Plots representing relative abundance of each keystone OTU and one representing the median relative abundance  across samples of rhizosphere of "$2". Most keystone OTUs have relative abundance bigger than the median across all samples.  }">> ../results/analisis/$1.tex

echo "   \label{key_otus_vs_medians_"$2"}">> ../results/analisis/$1.tex

echo "\end{figure}" >> ../results/analisis/$1.tex


#figura distribuciones todas vs key
echo "\begin{figure}">> ../results/analisis/$1.tex

echo   " \centering">> ../results/analisis/$1.tex

echo   " \includegraphics[scale = 0.75]{mean_median_key_vs_not_key_"$2".png}">> ../results/analisis/$1.tex

echo    "\caption{Boxes represent mean and standard error in the distribution of corresponding samples. Lines represent the corresponding medians. In these samples of rhizosphere of"$2"}">> ../results/analisis/$1.tex

echo   "\label{mean_median_"$2"}">> ../results/analisis/$1.tex

echo "\end{figure}">> ../results/analisis/$1.tex


#figura pcoa muestras
echo "\begin{figure}">> ../results/analisis/$1.tex

echo "   \centering">> ../results/analisis/$1.tex

echo "   \includegraphics[scale = 0.7]{pcoa_muestras_"$2".png}">> ../results/analisis/$1.tex

echo   " \caption{PCoA analysis with Bray-Curtis distance of rhizosphere samples of "$2".}">> ../results/analisis/$1.tex

echo   " \label{fig:"$2"_pcoa}">> ../results/analisis/$1.tex

echo "\end{figure}">> ../results/analisis/$1.tex


#figura pcoa key
echo "\begin{figure}">> ../results/analisis/$1.tex

echo  "  \centering">> ../results/analisis/$1.tex

echo  "  \includegraphics[scale = 0.7]{pcoa_key_otus_"$2".png}">> ../results/analisis/$1.tex

echo  "  \caption{PCoA analysis with Bray-Curtis distance of rhizosphere samples of "$2", restricted to keystone OTUs.}">> ../results/analisis/$1.tex

echo  "  \label{fig:"$2"_pcoa_key_otus}">> ../results/analisis/$1.tex

echo "\end{figure}">> ../results/analisis/$1.tex


