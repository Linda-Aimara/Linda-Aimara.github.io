
#########Paquetes para instalar##########
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("ensembldb")
BiocManager::install("biomaRt")
install.packages('gplots')
BiocManager::install("DESeq2")
install.packages('pheatmap')
install.packages('RColorBrewer')
install.packages('tidyverse')
BiocManager::install('EnhancedVolcano')
install.packages("here")
BiocManager::install("clusterProfiler")


#############Si se instalaron correctamente deberían de poder cargar 
#########la libreria del paquete
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ensembldb)
library(biomaRt)
library(gplots)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(EnhancedVolcano)
library(here)
library(clusterProfiler)