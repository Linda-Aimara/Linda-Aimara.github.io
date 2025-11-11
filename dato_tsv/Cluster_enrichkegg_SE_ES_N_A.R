######  Analisis de RNAseq todas las comparaciones 2020   #####
#Programa: Anotacion de datos de Expresion diferencial con Clusterprofiler
#Lugar: My home
#Desarrolladora: Linda Aimara Kempis Calanis
#########################################################
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
####Librerias necesarias####
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
packageDescription("clusterProfiler")
#1. Debemos elegir la direccion de los archivos
setwd("C:/Users/aimar/Documents/RNA seq kallisto/kallisto_155_dv_55_031019/abundance_tsv_trimmomatic/htseq_results/Tablas/")

#2 Leer las tablas necesarias para las comparaciones neonatos vs adultos
DE_US_neo_vs_adul<- read.table("DEG_LFC_1_NS_Neonate_CD4_+_T_Cells_vs_NS_Adult_CD4_+_T_Cells_.txt")
DE_CD3_CD28_neo_vs_adul<-read.table("DEG_LFC_1_CD3_CD28_Neonate_CD4_+_T_Cells_vs_CD3_CD28_Adult_CD4_+_T_Cells_.txt")
DE_CD3_Flag_neo_vs_adul<-read.table("DEG_LFC_1_CD3_Flag_Neonate_CD4_+_T_Cells_vs_CD3_Flag_Adult_CD4_+_T_Cells_.txt")

#3 Seleccionar el entrezgeid de los genes up and down de las tablas
padj.cutoff <- 0.05
lfc.cutoff <- 1
#US basales
US_sigN1 <- DE_US_neo_vs_adul %>% dplyr::filter(padj <= padj.cutoff & log2FoldChange >= lfc.cutoff)
dim(US_sigN1)
US_Neo_entrez<-US_sigN1$entrezgene_id
US_sigA1 <- DE_US_neo_vs_adul %>% dplyr::filter(padj <= padj.cutoff & log2FoldChange <= lfc.cutoff)
dim(US_sigA1)
US_Adul_entrez<-US_sigA1$entrezgene_id
#CD3_CD28
CD3_CD28_sigN1 <- DE_CD3_CD28_neo_vs_adul %>% dplyr::filter(padj <= padj.cutoff & log2FoldChange >= lfc.cutoff)
dim(CD3_CD28_sigN1)
CD3_CD28_Neo_entrez<-CD3_CD28_sigN1$entrezgene_id
CD3_CD28_sigA1 <- DE_CD3_CD28_neo_vs_adul %>% dplyr::filter(padj <= padj.cutoff & log2FoldChange <= lfc.cutoff)
dim(CD3_CD28_sigA1)
CD3_CD28_Adul_entrez<-CD3_CD28_sigA1$entrezgene_id
#CD3_Flag
CD3_Flag_sigN1 <- DE_CD3_Flag_neo_vs_adul %>% dplyr::filter(padj <= padj.cutoff & log2FoldChange >= lfc.cutoff)
dim(CD3_Flag_sigN1)
CD3_Flag_Neo_entrez<-CD3_Flag_sigN1$entrezgene_id
CD3_Flag_sigA1 <- DE_CD3_Flag_neo_vs_adul %>% dplyr::filter(padj <= padj.cutoff & log2FoldChange <= lfc.cutoff)
dim(CD3_Flag_sigA1)
CD3_Flag_Adul_entrez<-CD3_Flag_sigA1$entrezgene_id


#4 Crear una lista que tenga los cluster que quieres comparar
#mi_lista <- list("US_Neonatos" = US_Neo_entrez, "US_Adultos" = US_Adul_entrez, "CD3_C28_Neonatos" = CD3_CD28_Neo_entrez,"CD3_C28_Adultos" = CD3_CD28_Adul_entrez,"CD3_Flag_Neonatos" = CD3_Flag_Neo_entrez,"CD3_Flag_Adultos" = CD3_Flag_Adul_entrez )
#figura SE y ES
mi_lista <- list("US_Neonatos" = US_Neo_entrez, "US_Adultos" = US_Adul_entrez, "CD3_C28_Neonatos" = CD3_CD28_Neo_entrez,"CD3_C28_Adultos" = CD3_CD28_Adul_entrez)
mi_lista <- list("US_Neonatos" = US_Neo_entrez, "US_Adultos" = US_Adul_entrez, "CD3_C28_Neonatos" = CD3_CD28_Neo_entrez,"CD3_C28_Adultos" = CD3_CD28_Adul_entrez,"CD3_Flag_Neonatos" = CD3_Flag_Neo_entrez,"CD3_Flag_Adultos" = CD3_Flag_Adul_entrez )

head(ck1)
#view(mi_lista$US_neonatos)
#
if(!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("YuLab-SMU/clusterProfiler")
ck1 <- compareCluster(geneCluster = mi_lista, fun = enrichKEGG ,organism="", pvalueCutoff = 0.05, qvalueCutoff= 0.05)
pdf("cluster_enrichkegg_neo_vs_adultos_SE_ES.pdf")
head(as.data.frame(ck1))
x<-(as.data.frame(ck1))
write.table(x,file = "listas_cluster_SE_ES.txt", sep="\t", col.names = TRUE)
View(ck1)
dotplot(ck1,showCategory = 18,font.size=10.5)+ ggtitle("Functional enrichment analysis visualization")
dev.off()
