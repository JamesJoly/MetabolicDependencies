rm(list=ls(all=TRUE))

setwd("/.../SimulationStudy/")
data_in <- read.csv(file = "Synthetic_gene_set_correlation_coefficients.csv", header = TRUE, sep = ",")

library(tidyverse)
data_in$Gene <- NULL

long.data <- data_in %>% pivot_longer(-Replicate, names_to = "Gradients_added", values_to = "correlation_coefficient")

long.data$exp.value.added <- as.character(long.data$Gradients_added)
long.data$exp.value.added <- substr(long.data$exp.value.added, 1, 20)
long.data$exp.value.added <- gsub("Exp_value_added_","",long.data$exp.value.added)
long.data$exp.value.added <- gsub("_","",long.data$exp.value.added)
long.data$exp.value.added <- gsub("Ce","",long.data$exp.value.added)


long.data$ceres.value.added <- as.character(long.data$Gradients_added)

long.data$ceres.value.added <- substr(long.data$ceres.value.added, 22, nchar(long.data$ceres.value.added))
long.data$ceres.value.added <- gsub("[^0-9.-]", "",long.data$ceres.value.added)

## Calculate t test p-value
long.data$t.stat <- long.data$correlation_coefficient * sqrt(300-2) / sqrt(1-(long.data$correlation_coefficient)^2)
long.data$pval <- 2*pt(-abs(long.data$t.stat), df = 299)


## apply BH correction
replicates <- unique(long.data$Replicate)
long.data$FDR <- NA
long.data <- as.data.frame(long.data)
for(i in 1:length(replicates)){
  long.data[long.data$Replicate == replicates[i],]$FDR <- p.adjust(long.data[long.data$Replicate == replicates[i],]$pval, method = "BH", n = 16643)
}

## We want to calculate the % of p-values < some threshold (0.01, 0.001, 1e-4)

long.data.grouped <- group_by(long.data, exp.value.added, ceres.value.added)
long.data.percent <- summarise(long.data.grouped,
                               n_replicates  = n(),
                               percent.significant.01 = sum(FDR < 0.01) / n_replicates,
                               percent.significant.001 = sum(FDR < 0.001) / n_replicates,
                               percent.significant.0001 = sum(FDR < 0.0001) / n_replicates,)

library(viridis)
tile.plot <- ggplot(long.data.percent, aes(x = exp.value.added, y = ceres.value.added, fill = percent.significant.01)) +
  geom_tile() +
  #geom_text(aes(label = max.NES), size = 3) +
  scale_fill_viridis() + 
  scale_x_discrete(breaks = seq(0,1,0.1)) + scale_y_discrete(breaks = seq(0,1,0.1)) +
  labs(x = "Expression value added", y = "Ceres value added", fill = "",
       title = "Percent of correlation coefficients wih FDR < 1e-2") +
  theme(axis.text = element_text(color = "black", size = 10))


ggsave(tile.plot, file = "Percent_of_corr_coefficinets_with_FDR_less_than_1e-2.pdf", height = 7, width = 7)

