rm(list=ls(all=TRUE))

### Generate correlation matrix for all drugs in drug set
setwd("/.../Data/")
PRISM.data.full <- data.table::fread("secondary-screen-dose-response-curve-parameters.csv",header=T)

setwd("/.../ExpressionGSEA/")
GSEA.data.full <- data.table::fread("CCLE_GSEA_Metabolic_Pathways.csv",header=T)

#output dir
setwd("/.../IndividualDrugs/")

Adherent.RPMI <- GSEA.data.full[GSEA.data.full$Media == "RPMI" & GSEA.data.full$Culture.Type == "Adherent",]
Adherent.DMEM <- GSEA.data.full[GSEA.data.full$Media == "DMEM" & GSEA.data.full$Culture.Type == "Adherent",]

#Proceed w/ Adherent RPMI first, change it here for DMEM
GSEA.data.subsetted <- Adherent.RPMI
output.filename <- "Individual_Drug_Correlation_Results_Adherent_RPMI.csv"

#GSEA.data.subsetted <- Adherent.DMEM
#output.filename <- "Individual_Drug_Correlation_Results_Adherent_DMEM.csv"

# Data preparation --------------------------------------------------------
library(dplyr)
PRISM.data <- dplyr::select(PRISM.data.full, matches("ccle_name"), matches("auc"), matches("name"), -matches("row_name"), matches("direction"))
GSEA.data <- dplyr::select(GSEA.data.subsetted, matches("Cell"), matches("Pathway"), matches("KS_Normalized"))

rm(PRISM.data.full)
rm(GSEA.data.full)
rm(GSEA.data.subsetted)

PRISM.cells <- unique(PRISM.data$ccle_name) 
GSEA.cells <- unique(GSEA.data$Cell) 
GSEA.pathways <- unique(GSEA.data$Pathway) 
PRISM.drugs <- unique(PRISM.data$name)

cells <- intersect(GSEA.cells,PRISM.cells) #initial array trimming by overall cell overlap

PRISM.data <- PRISM.data[PRISM.data$ccle_name %in% cells,]
PRISM.data$auc <- -PRISM.data$auc #multiply by -1 to keep positive direction as negative growth
GSEA.data <- GSEA.data[GSEA.data$Cell %in% cells,]

### The following code was written by BC and reviewed by JJ

# Creating and filling correlation matrix -------------------------------

corr.matrix <- as.data.frame(matrix(data = 0,nrow = length(PRISM.drugs),ncol = length(GSEA.pathways)))
rownames(corr.matrix) <- PRISM.drugs
colnames(corr.matrix) <- GSEA.pathways

corr.results <- as.data.frame(matrix(data = 0,nrow = 0,ncol = 5))
colnames(corr.results) <- c("Metabolic_Pathway", "Drug_Name","Spearman_Corr","Corr_p_val","p_adj")

for (i in 1:length(PRISM.drugs)) { 
  PRISM.hold.drug <- PRISM.data[PRISM.data$name==PRISM.drugs[i],]
  temp.corr.results <- as.data.frame(matrix(data = NA,nrow = length(GSEA.pathways),ncol = 5))
  colnames(temp.corr.results) <- c("Metabolic_Pathway", "Drug_Name","Spearman_Corr","Corr_p_val","p_adj")
  temp.corr.results$Metabolic_Pathway <- GSEA.pathways
  
  for (j in 1:length(GSEA.pathways)) {
    GSEA.hold <- GSEA.data[GSEA.data$Pathway==GSEA.pathways[j],]
    #GSEA.hold <- GSEA.hold[abs(GSEA.hold$KS_Normalized) > 1.3,] #Filter for |NES| > 1.3
    
    hold.cells <- intersect(GSEA.hold$Cell,PRISM.hold.drug$ccle_name)
    PRISM.hold <- PRISM.hold.drug[PRISM.hold.drug$ccle_name %in% hold.cells,] #there are some drugs with multiple aucs per cell line
    dupes <- sum(duplicated(PRISM.hold$ccle_name))
    
    if (dupes>0) { 
      PRISM.duped <- PRISM.hold[(duplicated(PRISM.hold$ccle_name)|duplicated(PRISM.hold$ccle_name,fromLast = TRUE)),]
      PRISM.not.duped <- PRISM.hold[!(duplicated(PRISM.hold$ccle_name)|duplicated(PRISM.hold$ccle_name,fromLast = TRUE)),]
      u.dupes <- unique(PRISM.duped$ccle_name)
      PRISM.merge <- PRISM.not.duped[c(1:length(u.dupes)),]
      
      for (k in 1:length(u.dupes)) {
        holder <- PRISM.duped[PRISM.duped$ccle_name==u.dupes[k],]
        PRISM.merge[k,] <- holder[1,]
        PRISM.merge[k,2] <- mean(holder$auc)
      }
      PRISM.hold <- rbind(PRISM.not.duped,PRISM.merge)
    }
    
    GSEA.hold <- GSEA.hold[GSEA.hold$Cell %in% hold.cells,]
    #PRISM.hold <- PRISM.hold[order(PRISM.hold$ccle_name),]
    
    merged.data <- merge(GSEA.hold, PRISM.hold, by.x = "Cell", by.y = "ccle_name")
    #sometimes theres a single value with na, let's just remove it
    merged.data <- na.omit(merged.data)
    
    #value.holder <- cor(merged.data$KS_Normalized,(-1 * merged.data$auc),method = "spearman")
    corr.test <- cor.test(merged.data$KS_Normalized,merged.data$auc,method = "spearman")
    
    temp.corr.results[temp.corr.results$Metabolic_Pathway == GSEA.pathways[j],]$Drug_Name <- PRISM.drugs[i]
    temp.corr.results[temp.corr.results$Metabolic_Pathway == GSEA.pathways[j],]$Spearman_Corr <- corr.test$estimate
    temp.corr.results[temp.corr.results$Metabolic_Pathway == GSEA.pathways[j],]$Corr_p_val <- corr.test$p.value
    
  }
  corr.results <- rbind(corr.results, temp.corr.results)
  if(i %% 100 == 0){
    print(paste("Done with drug:", i))
  }
}

for (i in 1:length(GSEA.pathways)){
  corr.results[corr.results$Metabolic_Pathway == GSEA.pathways[i],]$p_adj <- p.adjust(corr.results[corr.results$Metabolic_Pathway == GSEA.pathways[i],]$Corr_p_val,
                                                                                      method= "BH")
}

write.csv(corr.results, file = output.filename, row.names = FALSE)


