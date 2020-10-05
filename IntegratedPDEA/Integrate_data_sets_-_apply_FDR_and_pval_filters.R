### Integrating Dependency GSEA and Pathway Inhibition Score data

setwd("/.../DependencyGSEA/")

Dependency.data <- read.csv("Genetic_PDEA_Adherent_RPMI.csv")
Dependency.data$Pathway.Expression <- Dependency.data$Sample
Dependency.data$Sample <- NULL
Dependency.data <- Dependency.data[,c(7,1:6)]

setwd("/.../PharmaPDEA/")

Inhibition.data <- read.csv("PharmaPDEA_Adherent_RPMI.csv")
Inhibition.data$Pathway.Expression <- Inhibition.data$ï..Pathway.Expression
Inhibition.data$ï..Pathway.Expression <- NULL
Inhibition.data <- Inhibition.data[,c(6,1:5)]

#output dir
setwd("/.../IntegratedPDEA/")


Dependency.data$Pathway.Expression <- gsub(" Correlation","", Dependency.data$Pathway.Expression)
Dependency.data$Gene.Set <- gsub(" Dependency","", Dependency.data$Gene.Set)

colnames(Dependency.data) <- gsub("Gene.Set","Dependency.Set", colnames(Dependency.data))

Inhibition.data$Pathway.Expression <- gsub("KEGG_Core_Glycolysis.Gluconeogenesis_hsa_M00001", "KEGG_Core_Glycolysis",Inhibition.data$Pathway.Expression)
Inhibition.data$Pathway.Expression <- gsub("KEGG_OXIDATIVE_PHOSPHORYLATION", "Oxidative_phosphorylation",Inhibition.data$Pathway.Expression)


Dependency.data$Pathway.Expression <- gsub("KEGG_","", Dependency.data$Pathway.Expression)
Dependency.data$Pathway.Expression <- gsub("_"," ", Dependency.data$Pathway.Expression)
Dependency.data$Pathway.Expression <- stringr::str_to_title(Dependency.data$Pathway.Expression)

Dependency.data$Dependency.Set <- gsub("KEGG_","", Dependency.data$Dependency.Set)
Dependency.data$Dependency.Set <- gsub("_"," ", Dependency.data$Dependency.Set)
Dependency.data$Dependency.Set <- stringr::str_to_title(Dependency.data$Dependency.Set)

Inhibition.data$Pathway.Expression <- gsub("KEGG_","", Inhibition.data$Pathway.Expression)
Inhibition.data$Pathway.Expression  <- gsub("_"," ", Inhibition.data$Pathway.Expression)
Inhibition.data$Pathway.Expression  <- stringr::str_to_title(Inhibition.data$Pathway.Expression)

Inhibition.data$Pathway.Dependency <- gsub("KEGG_","", Inhibition.data$Pathway.Dependency)
Inhibition.data$Pathway.Dependency  <- gsub("_"," ", Inhibition.data$Pathway.Dependency)
Inhibition.data$Pathway.Dependency  <- stringr::str_to_title(Inhibition.data$Pathway.Dependency)

colnames(Inhibition.data)[1:6] <- paste(colnames(Inhibition.data)[1:6], "PRISM") 
colnames(Dependency.data)[1:7] <- paste(colnames(Dependency.data)[1:7], "DepMap") 


Inhibition.data$Identifier <- paste(Inhibition.data$Pathway.Expression, Inhibition.data$`Pathway.Dependency PRISM`, sep = "|")
Dependency.data$Identifier <- paste(Dependency.data$Pathway.Expression, Dependency.data$`Dependency.Set DepMap`, sep = "|")


Merged.data.frames <- merge(Dependency.data, Inhibition.data, by = "Identifier")

write.csv(Merged.data.frames, file = "Merged_Genetic_and_Pharma_PDEA_results.csv", row.names = F)


