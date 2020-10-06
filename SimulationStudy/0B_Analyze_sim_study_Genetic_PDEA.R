setwd("/.../SimulationStudy/")

data_in <- read.csv("Simulation_study_Genetic_PDEA.csv", header = TRUE)

replicates <- unique(data_in$Simulation_Number)

## keep synthetic gene set only
data_in <- data_in[data_in$Gene.Set == "Synthetic_Gene_Set",]

data_in$Exp.value.added <- as.character(data_in$Sample)

data_in$Exp.value.added <- substr(data_in$Exp.value.added, 1, 20)
data_in$Exp.value.added <- gsub("Exp_value_added_","",data_in$Exp.value.added)
data_in$Exp.value.added <- gsub("_","",data_in$Exp.value.added)
data_in$Exp.value.added <- gsub("Ce","",data_in$Exp.value.added)

data_in$Ceres.value.added <- as.character(data_in$Sample)

data_in$Ceres.value.added <- substr(data_in$Ceres.value.added, 22, nchar(data_in$Ceres.value.added))
data_in$Ceres.value.added <- gsub("[^0-9.-]", "",data_in$Ceres.value.added)


### We want to ask how many p values < 0.05 for each of the 50 replicates at each of the combinations of ceres and exp values

exp.values <- unique(data_in$Exp.value.added)
ceres.values <- unique(data_in$Ceres.value.added)

results <- data.frame(matrix(data= NA, nrow = 0, ncol = 11))
colnames(results) <- c("exp.value.added", "ceres.value.added","num.pval.05", "num.pval.01", "num.qval.25", "num.qval.05",
                       "num.pval.05.qval.25", "num.pval.05.qval.05", "num.pval.01.qval.05", "mean.NES", "median.NES")
for (i in 1:length(exp.values)){
  
  temp.results <- data.frame(matrix(data= NA, nrow = length(ceres.values), ncol = 11))
  colnames(temp.results) <- c("exp.value.added", "ceres.value.added","num.pval.05", "num.pval.01", "num.qval.25", "num.qval.05",
                              "num.pval.05.qval.25", "num.pval.05.qval.05", "num.pval.01.qval.05", "mean.NES", "median.NES")
  temp.results$exp.value.added <- exp.values[i]
  
  for (j in 1:length(ceres.values)){
    temp.results$ceres.value.added[j] <- ceres.values[j]
    
    temp.df <- data_in[data_in$Exp.value.added == exp.values[i] & data_in$Ceres.value.added == ceres.values[j],]
    
    num.pvals.05 <- sum(temp.df$p_value < 0.05)
    num.pvals.01 <- sum(temp.df$p_value < 0.01)
    
    temp.results$num.pval.05[j] <- sum(temp.df$p_value < 0.05)
    temp.results$num.pval.01[j] <- sum(temp.df$p_value < 0.01)
    
    temp.results$num.qval.25[j] <- sum(temp.df$FDR_q_value < 0.25)
    temp.results$num.qval.05[j] <- sum(temp.df$FDR_q_value < 0.05)
    
    temp.results$num.pval.05.qval.25[j] <- sum(temp.df$p_value < 0.05 & temp.df$FDR_q_value < 0.25)
    temp.results$num.pval.05.qval.05[j] <- sum(temp.df$p_value < 0.05 & temp.df$FDR_q_value < 0.05)
    temp.results$num.pval.01.qval.05[j] <- sum(temp.df$p_value < 0.01 & temp.df$FDR_q_value < 0.25)
    
    temp.results$mean.NES[j] <- signif(mean(temp.df$KS_Normalized),digits = 3)
    temp.results$median.NES[j] <- signif(median(temp.df$KS_Normalized), digits = 3)
  }
  results <- rbind(results,temp.results)
}


results$num.pval.05.percent <- results$num.pval.05 / 50 * 100
results$num.pval.01.percent <- results$num.pval.01 / 50 * 100

results$num.qval.25.percent <- results$num.qval.25 / 50 * 100
results$num.qval.05.percent <- results$num.qval.05 / 50 * 100

results$num.pval.05.qval.25.percent <- results$num.pval.05.qval.25 / 50 * 100
results$num.pval.05.qval.05.percent <- results$num.pval.05.qval.05 / 50 * 100
results$num.pval.01.qval.05.percent <- results$num.pval.01.qval.05 / 50 * 100

library(tidyr)
results.long <- pivot_longer(results, cols = c("num.pval.05.percent","num.pval.01.percent",
                                               "num.qval.25.percent","num.qval.05.percent",
                                               "num.pval.05.qval.25.percent", "num.pval.05.qval.05.percent",
                                               "num.pval.01.qval.05.percent"),
                             names_to = "percent")
  
plots <- list()  

results.long$percent <- gsub("num.","", results.long$percent)
results.long$percent <- gsub(".percent","", results.long$percent)
types.of.percents <- unique(results.long$percent)

i = 1
library(viridis)
library(ggplot2)
for(i in 1:length(types.of.percents)){
  temp.df <- results.long[results.long$percent == types.of.percents[i],]
  temp.plot <- ggplot(temp.df, aes(x = exp.value.added, y = ceres.value.added, fill = value)) +
    geom_tile() +
    #geom_text(aes(label = median.NES), size = 3) +
    scale_fill_viridis() + 
    labs(x = "Expression value added", y = "Ceres value added", fill = "Percent",
         title = paste("Filter applied:", types.of.percents[i]))
  
  plots[[i]] <- temp.plot
}


library(gridExtra)
compiled.plots <- marrangeGrob(plots, nrow = 2, ncol = 2, top = NULL)


ggsave(compiled.plots, file ="Heatmaps_of_simulation_study_percent_significant_-_median_NES_label.pdf", height = 15, width = 15)  
