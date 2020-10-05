

setwd("/.../Data/")
Drug.Pathway.Targets <- read.csv("Metabolic_Pathways_drug_targets.csv", stringsAsFactors = F)

total.drugs <- unique(unlist(Drug.Pathway.Targets))

### Generate correlation matrix for all drugs in drug set
PRISM.data.full <- data.table::fread("secondary-screen-dose-response-curve-parameters.csv",header=T)

GSEA.data.full <- data.table::fread("CCLE_GSEA_Metabolic_Pathways.csv",header=T)

Adherent.RPMI <- GSEA.data.full[GSEA.data.full$Media == "RPMI" & GSEA.data.full$Culture.Type == "Adherent",]
Adherent.DMEM <- GSEA.data.full[GSEA.data.full$Media == "DMEM" & GSEA.data.full$Culture.Type == "Adherent",]

### Need to modify the PRISM data to account for mechanisms of action
mechanisms.of.action <-unique(PRISM.data.full$moa) #531 total

inhibitors <- stringr::str_detect(mechanisms.of.action,"inhibitor")
inhibitors <- mechanisms.of.action[inhibitors] #333

activators <- stringr::str_detect(mechanisms.of.action,"activator")
activators <- mechanisms.of.action[activators] #19

antagonist <- stringr::str_detect(mechanisms.of.action,"antagonist")
antagonist <- mechanisms.of.action[antagonist] #62

agonist <- stringr::str_detect(mechanisms.of.action," agonist") #need to include space to not get antagonists
agonist <- mechanisms.of.action[agonist] #51

blockers <- stringr::str_detect(mechanisms.of.action,"blocker")
blockers <- mechanisms.of.action[blockers] #12

positive.regulators <- c(activators, agonist)
negative.regulators <- c(inhibitors, antagonist, blockers)

### Add directionality for PRISM data
PRISM.data.full$direction <- 1
PRISM.data.full[PRISM.data.full$moa %in% negative.regulators,]$direction <- 1
PRISM.data.full[PRISM.data.full$moa %in% positive.regulators,]$direction <- -1

inhibitor.names <- PRISM.data.full[PRISM.data.full$moa %in% negative.regulators,]$name
activator.names <- PRISM.data.full[PRISM.data.full$moa %in% positive.regulators,]$name

#Proceed w/ Adherent RPMI first, change it here

GSEA.data.subsetted <- Adherent.RPMI
output.filename <- "Drug_corr_matrix_Adherent_RPMI.csv"
cell.num.threshold <- 150

#GSEA.data.subsetted <- Adherent.DMEM
#output.filename <- ""Drug_corr_matrix_Adherent_DMEM.csv"
#cell.num.threshold <- 90


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


## NOTE: For DSEA we need a background, so we don't want to do this step.
## JJ - lets filter for the PRISM.drugs that we already know affect metabolic pathways
#PRISM.drugs.filtered <- PRISM.drugs[PRISM.drugs %in% total.drugs]
PRISM.drugs.filtered <- PRISM.drugs #so we don't have to change below

# Creating and filling correlation matrix -------------------------------

corr.matrix <- as.data.frame(matrix(data = 0,nrow = length(PRISM.drugs.filtered),ncol = length(GSEA.pathways)))
rownames(corr.matrix) <- PRISM.drugs.filtered
colnames(corr.matrix) <- GSEA.pathways

for (i in 1:length(PRISM.drugs.filtered)) { 
  PRISM.hold.drug <- PRISM.data[PRISM.data$name==PRISM.drugs.filtered[i],]
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
    
    value.holder <- cor(merged.data$KS_Normalized,(-1 * merged.data$auc),method = "spearman")
    
    corr.matrix[PRISM.drugs.filtered[i],GSEA.pathways[j]] <- value.holder
  }
}

corr.matrix <- corr.matrix[,colSums(is.na(corr.matrix))<nrow(corr.matrix)]

## Something is weird with the drug aripiprazole
PRISM.hold.drug <- PRISM.data[PRISM.data$name=="aripiprazole",]
### It is becuase there are only 3 cell lines overlapping between the data sets!
### Lets see how pervasive the problem is

corr.matrix <- as.data.frame(corr.matrix)
corr.matrix$Drug.Name <- row.names(corr.matrix)
corr.matrix$num.cell.lines <- 0

for (i in 1:length(PRISM.drugs.filtered)){
  PRISM.hold.drug <- PRISM.data[PRISM.data$name==PRISM.drugs.filtered[i],]
  
  GSEA.hold <- GSEA.data[GSEA.data$Pathway==GSEA.pathways[1],] #for this step we can use just one pathway
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
  
  corr.matrix[corr.matrix$Drug.Name == PRISM.drugs.filtered[i],]$num.cell.lines <- length(unique(merged.data$Cell))
}

corr.matrix <- corr.matrix[,c(71:72,1:70)]

#filter for drugs with more than 150 for adherent RPMI and over 90 for adherent DMEM
corr.matrix <- corr.matrix[,"num.cell.lines" >= cell.num.threshold]

write.csv(corr.matrix, file = output.filename, row.names = FALSE)
