rm(list=ls(all=TRUE))

setwd("/.../Data/")
data_in_DMEM <- read.csv(file = "Genetic_PDEA_Adherent_DMEM.csv", header = TRUE, sep = ",")
data_in_RPMI <- read.csv(file = "Genetic_PDEA_Adherent_RPMI.csv", header = TRUE, sep = ",")


data_in_DMEM$Weight <- -log10(data_in_DMEM$FDR_q_value + 1e-6)
data_in_DMEM$Weighted.KS_Normalized <- data_in_DMEM$KS_Normalized * data_in_DMEM$Weight

data_in_RPMI$Weight <- -log10(data_in_DMEM$FDR_q_value + 1e-6)
data_in_RPMI$Weighted.KS_Normalized <- data_in_RPMI$KS_Normalized * data_in_RPMI$Weight

library(dplyr)
library(tidyr)

data_in_RPMI.filtered <- select(data_in_RPMI, Weighted.KS_Normalized, Gene.Set,Sample)

Weighted.KS_Normalized.Mat.RPMI <- pivot_wider(data_in_RPMI.filtered, values_from = "Weighted.KS_Normalized", names_from = "Sample")

Weighted.KS_Normalized.Mat.RPMI$RowSums <- rowSums(Weighted.KS_Normalized.Mat.RPMI[,2:69])
Weighted.KS_Normalized.Mat.RPMI$Mean.KS_Normalized <- rowMeans(Weighted.KS_Normalized.Mat.RPMI[,2:69])

#repeat for DMEM
data_in_DMEM.filtered <- select(data_in_DMEM, Weighted.KS_Normalized, Gene.Set,Sample)

Weighted.KS_Normalized.Mat.DMEM <- pivot_wider(data_in_DMEM.filtered, values_from = "Weighted.KS_Normalized", names_from = "Sample")

Weighted.KS_Normalized.Mat.DMEM$RowSums <- rowSums(Weighted.KS_Normalized.Mat.DMEM[,2:69])
Weighted.KS_Normalized.Mat.DMEM$Mean.KS_Normalized <- rowMeans(Weighted.KS_Normalized.Mat.DMEM[,2:69])

Weighted.Mean.RPMI <- select(Weighted.KS_Normalized.Mat.RPMI, matches("Gene.Set"), matches("Mean.KS_Normalized"))
Weighted.Mean.DMEM <- select(Weighted.KS_Normalized.Mat.DMEM, matches("Gene.Set"), matches("Mean.KS_Normalized"))

merged.results <- merge(Weighted.Mean.DMEM, Weighted.Mean.RPMI, by = "Gene.Set")
colnames(merged.results)[2] <- "DMEM"
colnames(merged.results)[3] <- "RPMI"

merged.results$Gene.Set <- gsub(" Dependency","", merged.results$Gene.Set)
merged.results$Gene.Set <- gsub("KEGG_", "", merged.results$Gene.Set)
merged.results$Gene.Set <- gsub("_", " ", merged.results$Gene.Set)
merged.results$Gene.Set <- stringr::str_to_title(merged.results$Gene.Set)

#going to use Excel to VLOOKUP the recipe fold changes
setwd("/.../DependencyGSEA/")
write.csv(merged.results, file = "Weighted_Mean_NES_for_Adherent_RPMI_and_DMEM_-_log10fdr_weight.csv", row.names = F)

## read back in data
for.plotting <- read.csv("Weighted_Mean_NES_for_Adherent_RPMI_and_DMEM_-_log10fdr_weight.csv")

library(tidyverse)
library(ggplot2)
library(aplot)
library(cowplot)

for.plotting$Differece <- for.plotting$RPMI - for.plotting$DMEM

for.plotting$Gene.Set <- reorder(for.plotting$Gene.Set, -for.plotting$Differece)

for.plotting.long <- pivot_longer(for.plotting, cols = c("RPMI","DMEM"))


dmem.vs.rpmi <- ggplot(for.plotting.long) +
  geom_bar(aes(x = Gene.Set, y = value, fill = name), position = "dodge", stat = "identity", color = "black") +
  scale_fill_manual(values = c("red1","blue1")) + coord_flip() +
  labs(x = "Pathway", y = "Weighted Mean NES", fill = "", title = "Differences in Pathway Dependency based on Media\n Postive scores indicate increased dependency") +
  cowplot::theme_cowplot() +
  theme(legend.position = "bottom", legend.title.align = 0.5, legend.background = element_blank(),
        panel.border = element_blank(), plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.line.y = element_line(color = "black", size = 0.7),
        axis.line.x = element_line(color = "black", size = 0.7),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10))


tiles <- ggplot(for.plotting, aes(x = Gene.Set, y = 1)) + geom_tile(aes(fill = Log2.RPMI.DMEM..Recipe), size = 10) +
  scale_fill_gradient2(low = "#AF8DC3", high = "#7FBF7B", mid = "#F7F7F7", midpoint = 0,
                       limits = c(-2,2), oob = scales::squish)+
  cowplot::theme_nothing() + theme(legend.position = "bottom") + labs(fill = "Log2 (RPMI / DMEM) Recipe") +
  xlim2(dmem.vs.rpmi)  + ylim2(dmem.vs.rpmi) + coord_flip()
#tiles


compiled <- plot_grid(dmem.vs.rpmi,tiles, nrow = 1, rel_widths = c(10,10), align = "h")
ggsave(compiled, file = "DMEM_vs_RPMI_Weighted_value_NES_comparison_-_-log10_fdr_weight.pdf", height = 12, width = 18)
