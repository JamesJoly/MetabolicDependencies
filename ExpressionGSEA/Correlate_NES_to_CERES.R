rm(list=ls(all=TRUE))
library(data.table)
library(dplyr)
library(tidyr)

setwd("/.../Data/")

##Read in DEPMAP data
setwd("C:/Users/James/Documents/GitHub/MetabolicDependencies/Data")
DEPMAP <- fread("Achilles_gene_effect.csv", header = T, sep =",")
##Read in cell line info / identifiers
cell_line_info <- read.csv("sample_info.csv", header = T, sep = ",")

##Read in GSEA data
GSEA.data <- fread("CCLE_GSEA_Metabolic_Pathways.csv", header = T)

# Reformat gene names
DEPMAP_genes <- colnames(DEPMAP)
DEPMAP_genes <- gsub(" .*$", "", DEPMAP_genes)
DEPMAP_genes[1] <- "DepMap_ID"
colnames(DEPMAP) <- DEPMAP_genes

#DEPMAP.Met <- as.data.frame(DEPMAP.Met)
DEPMAP <- as.data.frame(DEPMAP)

rownames(DEPMAP) <- DEPMAP$DepMap_ID
#remove Depmap id
DEPMAP$DepMap_ID <- NULL
identified.genes <- colnames(DEPMAP)

# Change CCLE name to DepMap ID in DPEA
identifiers <- select(cell_line_info, "CCLE_Name","DepMap_ID","culture_type")

GSEA <- select(GSEA.data, "Cell","Media","Pathway","KS_Normalized")
#GSEA <- GSEA[,1:5]
GSEA <- as.data.frame(GSEA)

GSEA <- merge(GSEA, identifiers, by.x = "Cell",by.y = "CCLE_Name")
GSEA$Cell <- NULL #remove Cell column

depmap.cells <- rownames(DEPMAP)
gsea.cells <- unique(GSEA$DepMap_ID)

GSEA.filtered <- GSEA[which(GSEA$DepMap_ID %in% depmap.cells),]
DepMap.filtered <- DEPMAP[which(rownames(DEPMAP) %in% gsea.cells),]

media.types <- factor(unique(GSEA$Media))
culture.types <- factor(unique(GSEA$culture_type))
Pathways <- unique(GSEA$Pathway)
Correlations.All <- matrix(data = NA, nrow = 0, ncol = 7)
colnames(Correlations.All) <- c("Media","Culture_type","Pathway","Gene","corr.coefficient","p_value","p_adj")
Correlations.All <- as.data.frame(Correlations.All)


#Correlations.All$Gene <- identified.genes

#GSEA.filtered$culture_type <- NULL

library(tidyverse)

    for (k in 1:length(media.types)){
    temp.df <- GSEA.filtered[GSEA.filtered$Media == media.types[k] & GSEA.filtered$culture_type == "Adherent",]
    
    temp.cells <- factor(unique(temp.df$DepMap_ID))
    if (length(temp.cells) > 10){
      
      
      temp.depmap <- DepMap.filtered[which(rownames(DepMap.filtered) %in% temp.cells),]
    
      temp.df$Media <- NULL
      temp.df$culture_type <- NULL
      temp.df.wide <- as.data.frame(pivot_wider(temp.df, names_from = "DepMap_ID", values_from = "KS_Normalized"))
      rownames(temp.df.wide) <- temp.df.wide$Pathway
      #Pathways <- as.character(temp.df.wide$Pathway)
      temp.df.wide$Pathway <- NULL

    
      for (i in 1:length(Pathways)){
        temp.pathway <- t(temp.df.wide[Pathways[i],])
        #temp.pathway <- temp.pathway[temp.pathway > 1.3 | temp.pathway < -1.3, , drop = F]
        
        temp.corr.mat <- merge(temp.pathway,temp.depmap, by = 0)
        
        Correlations <- data.frame(matrix(data = NA, nrow = length(identified.genes), ncol = 7))
        colnames(Correlations) <- c("Media","Culture_type","Pathway","Gene","corr.coefficient","p_value","p_adj")
        #rownames(Correlations) <- identified.genes
        Correlations[,"Media"] <- as.character(media.types[k])
        Correlations[,"Culture_type"] <- as.character(culture.types[1])
        Correlations[,"Pathway"] <- Pathways[i]
        Correlations[,"Gene"] <- identified.genes
        
        
        for (j in 1:length(identified.genes)){
          temp.cor <- cor.test(x = -1*temp.corr.mat[,identified.genes[j]], y = temp.corr.mat[,Pathways[i]], method = "spearman")
          
          Correlations[Correlations$Gene == identified.genes[j],]$corr.coefficient <- temp.cor$estimate
          Correlations[Correlations$Gene == identified.genes[j],]$p_value <- temp.cor$p.value
          
        }
        Correlations$p_adj <- p.adjust(Correlations$p_value, method = "BH")
        #print(paste("pathway complete:", i)) #uncomment for progress
      
      #Correlations <- as.data.frame(Correlations)
      #Correlations$Gene <- rownames(Correlations)
      #rownames(Correlations) <- NULL
      Correlations.All <- rbind(Correlations.All, Correlations)
      }
      print(paste("medium done:", media.types[k]))
    }
  }

write.csv(Correlations.final, file = "Corr_coefficients_NES_to_CERES.csv")


