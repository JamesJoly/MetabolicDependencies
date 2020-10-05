setwd("/.../IndividualDrugs/")

DMEM <- read.csv(file = "Individual_Drug_Correlation_Results_Adherent_DMEM.csv", header = TRUE)
RPMI <- read.csv(file = "Individual_Drug_Correlation_Results_Adherent_RPMI.csv", header = TRUE)

RPMI.gene.sets <- unique(RPMI$Metabolic_Pathway)
DMEM.gene.sets <- unique(DMEM$Metabolic_Pathway)

all.equal(RPMI.gene.sets, DMEM.gene.sets)

colnames(RPMI)[1] <- "Metabolic_Pathway"
colnames(RPMI) <- paste(colnames(RPMI), "RPMI")
colnames(DMEM) <- paste(colnames(DMEM), "DMEM")

RPMI$ID <- paste(RPMI$`Metabolic_Pathway RPMI`, RPMI$`Drug_Name RPMI`, sep = "|")
DMEM$ID <- paste(DMEM$`Metabolic_Pathway DMEM`, DMEM$`Drug_Name DMEM`, sep = "|")


merged.results <- merge(RPMI, DMEM, by = "ID")

merged.results$Signif.DMEM <- "No"
merged.results[merged.results$`p_adj DMEM` < 0.05,]$Signif.DMEM <- "Yes"

merged.results$Signif.RPMI <- "No"
merged.results[merged.results$`p_adj RPMI` < 0.05,]$Signif.RPMI <- "Yes"

merged.results$signif.both <- "No"
merged.results[merged.results$Signif.DMEM == "Yes" & merged.results$Signif.RPMI == "Yes",]$signif.both <- "Yes"

merged.results$same.sign <- "No"
merged.results[merged.results$`Spearman_Corr DMEM` < 0 & merged.results$`Spearman_Corr RPMI` < 0 |
                 merged.results$`Spearman_Corr DMEM` > 0 & merged.results$`Spearman_Corr RPMI` > 0,]$same.sign <- "Yes"

#count number successes
nrow(merged.results[merged.results$same.sign == "Yes" & merged.results$signif.both == "Yes",]) #66


write.csv(merged.results, file = "Combined_RPMI_and_DMEM_spearman_correlation_coefficients_individual_drugs.csv", row.names = F)

