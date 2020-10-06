rm(list=ls(all=TRUE))
setwd("/.../Data/")

# below code is by Brandon Chew
PRISM.data.full <- data.table::fread("secondary-screen-dose-response-curve-parameters.csv",header=T)

setwd("/.../ExpressionGSEA/")
GSEA.data.full <- data.table::fread("CCLE GSEA by culture medium AND culture type- new pathway names 04_07_2020.csv",header=T)

#output dir
setwd("/.../IndividualDrugs/")

Adherent.RPMI <- GSEA.data.full[GSEA.data.full$Media == "RPMI" & GSEA.data.full$Culture.Type == "Adherent",]
Adherent.DMEM <- GSEA.data.full[GSEA.data.full$Media == "DMEM" & GSEA.data.full$Culture.Type == "Adherent",]

#Proceed w/ Adherent RPMI first, change it here

#GSEA.data.subsetted <- Adherent.RPMI
#GSEA.data.subsetted <- Adherent.DMEM

# Data preparation --------------------------------------------------------
library(dplyr)
PRISM.data <- dplyr::select(PRISM.data.full, matches("ccle_name"), matches("auc"), matches("name"), -matches("row_name"), matches("direction"))
DMEM.data <- dplyr::select(Adherent.DMEM, matches("Cell"), matches("Pathway"), matches("KS_Normalized"))
RPMI.data <- dplyr::select(Adherent.RPMI, matches("Cell"), matches("Pathway"), matches("KS_Normalized"))

rm(PRISM.data.full)
rm(GSEA.data.full)
#rm(GSEA.data.subsetted)

PRISM.cells <- unique(PRISM.data$ccle_name) 
DMEM.cells <- unique(DMEM.data$Cell) 
RPMI.cells <- unique(RPMI.data$Cell) 

GSEA.pathways <- unique(DMEM.data$Pathway) 
PRISM.drugs <- unique(PRISM.data$name)

DMEM.PRISM.cells <- intersect(DMEM.cells,PRISM.cells) #initial array trimming by overall cell overlap
RPMI.PRISM.cells <- intersect(RPMI.cells,PRISM.cells) #initial array trimming by overall cell overlap

PRISM.DMEM.data <- PRISM.data[PRISM.data$ccle_name %in% DMEM.PRISM.cells,]
PRISM.RPMI.data <- PRISM.data[PRISM.data$ccle_name %in% RPMI.PRISM.cells,]

PRISM.DMEM.data$auc <- -PRISM.DMEM.data$auc
PRISM.RPMI.data$auc <- -PRISM.RPMI.data$auc

GSEA.DMEM.data <- DMEM.data[DMEM.data$Cell %in% DMEM.PRISM.cells,]
GSEA.RPMI.data <- RPMI.data[RPMI.data$Cell %in% RPMI.PRISM.cells,]


pathway.queried <- "KEGG_LINOLEIC_ACID_METABOLISM"
drug.queried <- "afatinib"

make.combined.scatter.plots <- function(pathway.queried, drug.queried){

  temp.DMEM.data <- GSEA.DMEM.data[GSEA.DMEM.data$Pathway == pathway.queried,]
  temp.RPMI.data <- GSEA.RPMI.data[GSEA.RPMI.data$Pathway == pathway.queried,]
  
  temp.PRISM.RPMI.data <- PRISM.RPMI.data[PRISM.RPMI.data$name == drug.queried,]
  temp.PRISM.DMEM.data <- PRISM.DMEM.data[PRISM.DMEM.data$name == drug.queried,]
  
  dupes.RPMI <- length(duplicated(temp.PRISM.RPMI.data$ccle_name))
  
  if (dupes.RPMI>0) { 
    PRISM.duped <- temp.PRISM.RPMI.data[(duplicated(temp.PRISM.RPMI.data$ccle_name)|duplicated(temp.PRISM.RPMI.data$ccle_name,fromLast = TRUE)),]
    PRISM.not.duped <- temp.PRISM.RPMI.data[!(duplicated(temp.PRISM.RPMI.data$ccle_name)|duplicated(temp.PRISM.RPMI.data$ccle_name,fromLast = TRUE)),]
    u.dupes <- unique(PRISM.duped$ccle_name)
    PRISM.merge <- PRISM.not.duped[c(1:length(u.dupes)),]
    
    for (k in 1:length(u.dupes)) {
      holder <- PRISM.duped[PRISM.duped$ccle_name==u.dupes[k],]
      PRISM.merge[k,] <- holder[1,]
      PRISM.merge[k,2] <- mean(holder$auc)
    }
    PRISM.RPMI.hold <- rbind(PRISM.not.duped,PRISM.merge)
  }
  
  dupes.DMEM <- length(duplicated(temp.PRISM.DMEM.data$ccle_name))
  
  if (dupes.DMEM>0) { 
    PRISM.duped <- temp.PRISM.DMEM.data[(duplicated(temp.PRISM.DMEM.data$ccle_name)|duplicated(temp.PRISM.DMEM.data$ccle_name,fromLast = TRUE)),]
    PRISM.not.duped <- temp.PRISM.DMEM.data[!(duplicated(temp.PRISM.DMEM.data$ccle_name)|duplicated(temp.PRISM.DMEM.data$ccle_name,fromLast = TRUE)),]
    u.dupes <- unique(PRISM.duped$ccle_name)
    PRISM.merge <- PRISM.not.duped[c(1:length(u.dupes)),]
    
    for (k in 1:length(u.dupes)) {
      holder <- PRISM.duped[PRISM.duped$ccle_name==u.dupes[k],]
      PRISM.merge[k,] <- holder[1,]
      PRISM.merge[k,2] <- mean(holder$auc)
    }
    PRISM.DMEM.hold <- rbind(PRISM.not.duped,PRISM.merge)
  }
  
  
  merged.DMEM <- merge(PRISM.DMEM.hold, temp.DMEM.data, by.x = "ccle_name", by.y = "Cell")
  merged.RPMI <- merge(PRISM.RPMI.hold, temp.RPMI.data, by.x = "ccle_name", by.y = "Cell")
  
  DMEM.corr.test <- cor.test(merged.DMEM$auc, merged.DMEM$KS_Normalized, method = "spearman")
  DMEM.spearman <- signif(DMEM.corr.test$estimate, digits = 3)
  DMEM.pval <- p.adjust(signif(DMEM.corr.test$p.value, digits = 3), method = "BH", n = 1448)
  DMEM.title <- paste("DMEM","\nrho =", DMEM.spearman, "\np_adj =")
  
  
  RPMI.corr.test <- cor.test(merged.RPMI$auc, merged.RPMI$KS_Normalized, method = "spearman")
  RPMI.spearman <- signif(RPMI.corr.test$estimate, digits = 3)
  RPMI.pval <- p.adjust(signif(RPMI.corr.test$p.value, digits = 3), method = "BH", n = 1448)
  RPMI.title <- paste("RPMI","\nrho =", RPMI.spearman, "\np_adj =")
  
  library(ggplot2)
  scatter.dmem <- ggplot(merged.DMEM, aes(x = rank(auc), y = rank(KS_Normalized))) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm") +
    labs(x = paste(drug.queried, "Rank"), y = "",
         title = DMEM.title) +
    cowplot::theme_cowplot(font_size = 8) +
    theme(plot.title = element_text(hjust = 0.5))
  #scatter.dmem
  
  scatter.rpmi <- ggplot(merged.RPMI, aes(x = rank(auc), y = rank(KS_Normalized))) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm") +
    labs(x = paste(drug.queried, "Rank"), y = paste(pathway.queried,"Rank"),
         title = RPMI.title) +
    cowplot::theme_cowplot(font_size = 8) +
    theme(plot.title = element_text(hjust = 0.5))
  #scatter.rpmi
  
  library(patchwork)
  compiled <- scatter.rpmi + scatter.dmem
  
  setwd("C:/Users/James/Documents/Research/Metabolic\ Dependencies/Individual_Drug_Response_Correlations/Combine_RPMI_DMEM")
  ggsave(compiled, file = paste(drug.queried, pathway.queried, "RPMI and DMEM scatter plots.pdf"), height = 4, width = 5)
}

make.combined.scatter.plots(pathway.queried = "KEGG_Core_Glycolysis.Gluconeogenesis_hsa_M00001",
                            drug.queried = "AZD8931")

make.combined.scatter.plots(pathway.queried = "KEGG_LINOLEIC_ACID_METABOLISM",
                            drug.queried = "afatinib")

make.combined.scatter.plots(pathway.queried = "KEGG_PURINE_METABOLISM",
                            drug.queried = "poziotinib")


make.combined.scatter.plots(pathway.queried = "KEGG_STEROID_HORMONE_BIOSYNTHESIS",
                            drug.queried = "poziotinib")

make.combined.scatter.plots(pathway.queried = "KEGG_O_GLYCAN_BIOSYNTHESIS",
                            drug.queried = "NMS-E973")

make.combined.scatter.plots(pathway.queried = "KEGG_PHENYLALANINE_METABOLISM",
                            drug.queried = "atorvastatin")

