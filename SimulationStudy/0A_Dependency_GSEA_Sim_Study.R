### We have 300 cell lines in Adherent RPMI
### We want to know the true positive / false negative rates of enrichment
### So we will create synthetic data with enrichments for both Gene Expression and CERES Score.
###
### How?
###
### We will vary gene expression data and CERES scores in the same directions, then calculate:
### 1) ssGSEA NES's (for 1 pathway, N = 25 genes)
### 2) Correlation coefficients for each gene (16k genes x 1 pathway)
### 3) Dependency GSEA - p values are what we want here
### 
### How will we visualize the data?
###
### A scatter plot where x-axis is CERES score perturbation and y-axis is gene expression pertubation
### Dots will be % of true significant results (p < 0.05)


### How to build this?

### Step 1: Generate gene expression matrix with a normal distribution for each cell line, then perturb 25 genes in gene set
###         Note: need to name the cell lines such that it's easy to track, we will want to perturb CERES in the same direction
### Step 2: Calculate GSEA NES for the synthetic gene set.
###
### Step 3: Generate CERES scores that have a normal distribution for each cell line, then perturb the same 25 genes
###         Note: since we will be doing spearman correlations, we don't have to worry about the range for CERES scores being not centered around 0.
### Step 4: Generate correlation coefficients
###
### Step 5: Run dependency GSEA for each different variant (+/- X gene expression, +/- X CERES)
###
### Step 6: Re-run for N replicates. Depending on how long it takes, 20? 50? 100?


# Let's use the correct # of genes and gene names
setwd("/.../ExpressionGSEA/")

real.data <- read.csv("Corr_coefficients_NES_to_CERES.csv", header = TRUE)
genes <- unique(real.data[,"Gene"])

setwd("/.../Data/")
KEGG.pathways <- GSA::GSA.read.gmt(file = "KEGG_metabolic_pathways_for_sim_study.gmt")

setwd("/.../SimulationStudy/")
synthetic.gene.set.gmt <- GSA::GSA.read.gmt(file = "synthetic_gene_set_for_sim_study.gmt")

num.genes <- length(genes)

# Make a synthetic gene set of size 25 genes
gene.set.size <- 25
#synthetic.gene.set <- sample(genes,gene.set.size) # we did this once and now we want to use the same gene set for each time, so lets save it separately

synthetic.gene.set <-  c("TACR3","PFDN6","ARHGAP10","STPG3","ADAMTS14","UACA","BTBD2","REM1", "MRC1", "RBM4B","RPN2",
                         "ARL8B","HCN1", "ZNF773","CEACAM6","ABCD3","KARS", "COQ5", "ARPP19","SOCS3","VLDLR","SAMD11","CMTM1","XKR3", "NUP85")

# Make synthetic cell line names
synthetic.cell.names <- seq(from = 1, to = 300, by = 1)
synthetic.cell.names <- paste("Cell_No_", synthetic.cell.names, sep = "")


GSEA.results.compiled.all.replicates <- data.frame(matrix(data = NA, nrow = 0, ncol = 8))
colnames(GSEA.results.compiled.all.replicates) <- c("Sample","Gene.Set", "KS", "KS_Normalized","p_value",
                                                    "Position_at_max", "FDR_q_value", "Simulation_Number")

GSEA_custom <- function(input.df, gmt.list,
                        num.permutations = 1000,
                        stat.type = "Weighted"){
  loop.time <- Sys.time()
  
  nperm = num.permutations #number of permutations
  if (stat.type == "Classic"){
    score.weight = 0
  }
  if (stat.type == "Weighted"){
    score.weight = 1
  }
  
  
  #Read in gene expression data
  #Genes should be first column, named "Gene"
  #Samples should be columns 2:N
  data_in <- input.df
  
  gmt.for.reformat <- gmt.list
  Gene.Sets <- t(plyr::ldply(gmt.for.reformat$genesets, rbind)) #reformat gmt list to desired format
  colnames(Gene.Sets) <- gmt.for.reformat$geneset.names
  
  Gene.Sets <- as.data.frame(Gene.Sets)
  
  testthat::expect_is(data_in, "data.frame")
  testthat::expect_is(Gene.Sets, "data.frame")
  
  GSEA.EnrichmentScore <- function(gene.list, gene.set, weighted.score.type = score.weight, correl.vector = NULL){
    tag.indicators <- sign(match(gene.list, gene.set, nomatch = 0))
    no.tag.indicator <- 1 - tag.indicators
    N <- length(gene.list)
    Nh <- numhits_pathway
    Nm <- N - Nh
    if (weighted.score.type == 0){
      correl.vector <- rep(1,N)
    }
    alpha <- weighted.score.type
    correl.vector <- abs(correl.vector**alpha)
    sum.correl.tag <- sum(correl.vector[tag.indicators == 1])
    norm.tag <- 1.0/sum.correl.tag
    norm.no.tag <- 1.0/Nm
    RES <- cumsum(tag.indicators * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)
    max.ES <- max(RES)
    min.ES <- min(RES)
    if (max.ES > - min.ES) {
      #      ES <- max.ES
      ES <- signif(max.ES, digits = 5)
      arg.ES <- which.max(RES)
    } else {
      #      ES <- min.ES
      ES <- signif(min.ES, digits=5)
      arg.ES <- which.min(RES)
    }
    return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicators))
  } #for real ES
  
  GSEA.EnrichmentScore2 <- function(gene.list, gene.set, weighted.score.type = score.weight, correl.vector = NULL) {
    
    N <- length(gene.list)
    Nh <- numhits_pathway
    Nm <-  N - Nh
    
    loc.vector <- vector(length=N, mode="numeric")
    peak.res.vector <- vector(length=Nh, mode="numeric")
    valley.res.vector <- vector(length=Nh, mode="numeric")
    tag.correl.vector <- vector(length=Nh, mode="numeric")
    tag.diff.vector <- vector(length=Nh, mode="numeric")
    tag.loc.vector <- vector(length=Nh, mode="numeric")
    
    loc.vector[gene.list] <- seq(1, N)
    tag.loc.vector <- loc.vector[gene.set]
    
    tag.loc.vector <- sort(tag.loc.vector, decreasing = F)
    
    if (weighted.score.type == 0) {
      tag.correl.vector <- rep(1, Nh)
    } else if (weighted.score.type == 1) {
      tag.correl.vector <- correl.vector[tag.loc.vector]
      tag.correl.vector <- abs(tag.correl.vector)
    } else if (weighted.score.type == 2) {
      tag.correl.vector <- correl.vector[tag.loc.vector]*correl.vector[tag.loc.vector]
      tag.correl.vector <- abs(tag.correl.vector)
    } else {
      tag.correl.vector <- correl.vector[tag.loc.vector]**weighted.score.type
      tag.correl.vector <- abs(tag.correl.vector)
    }
    
    norm.tag <- 1.0/sum(tag.correl.vector)
    tag.correl.vector <- tag.correl.vector * norm.tag
    norm.no.tag <- 1.0/Nm
    tag.diff.vector[1] <- (tag.loc.vector[1] - 1)
    tag.diff.vector[2:Nh] <- tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
    tag.diff.vector <- tag.diff.vector * norm.no.tag
    peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
    valley.res.vector <- peak.res.vector - tag.correl.vector
    max.ES <- max(peak.res.vector)
    min.ES <- min(valley.res.vector)
    ES <- signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits=5)
    
    return(ES)
    
  } #for permutation ES
  
  Samples <- colnames(data_in)
  if (Samples[1] != "Gene"){
    stop("Please ensure that your data frame is organized with the first column to be named 'Gene'")
  }
  Samples <- Samples[-1]
  
  Gene.Sets.All <- colnames(Gene.Sets)
  
  annotations <- matrix(data = 0, nrow = nrow(data_in), ncol = length(Gene.Sets.All))
  colnames(annotations) <- Gene.Sets.All
  
  annotations <- as.data.frame(annotations)
  
  annotations <- cbind(data_in$Gene,annotations)
  colnames(annotations) <- c("Gene", Gene.Sets.All)
  annotations <- as.matrix(annotations)
  
  num.hits.pathways <- list()
  
  ### Annotate gene sets
  for (j in 1:length(Gene.Sets.All)){
    temp.pathway <- Gene.Sets[,Gene.Sets.All[j]]
    for (i in 1:nrow(annotations)){
      if (annotations[i,"Gene"] %in% temp.pathway){
        annotations[i,j+1] = "X";
      }
    }
    num.hits.pathways[[Gene.Sets.All[j]]] <- sum(annotations[,Gene.Sets.All[j]] == "X")
  }
  
  num.hits.pathways.df <- matrix(unlist(num.hits.pathways))
  row.names(num.hits.pathways.df) = Gene.Sets.All
  num.gene.sets.under.5 <- which(num.hits.pathways.df < 5)
  if (length(num.gene.sets.under.5) > 1){
    print("Warning: Removing gene sets with less than 5 genes observed in data set.")
    gene.sets.to.remove <- Gene.Sets.All[num.gene.sets.under.5]
    annotations[,which(colnames(annotations) %in% gene.sets.to.remove)] <- NULL
  }
  annotations <- as.data.frame(annotations)
  data_in <- merge(data_in, annotations, by = "Gene")
  
  data_in <- stats::na.omit(data_in)
  
  GSEA.Results.All.Samples <- matrix(data = NA, nrow = 0, ncol = 7)
  colnames(GSEA.Results.All.Samples) <- c("Sample","Gene.Set","KS","KS_Normalized",
                                          "p-value","Position at Max",
                                          "FDR q-value")
  Mountain.Plot.Info.All.Samples <- list()
  rank_metric.All.Samples <- list()
  
  #Find out how many cores are available (if you don't already know)
  cores<-parallel::detectCores()
  #Create cluster with desired number of cores, leave one open for the machine
  #core processes
  cl <- snow::makeCluster(cores[1]-1)
  #Register cluster
  doSNOW::registerDoSNOW(cl)
  
  rm(annotations)
  data_in2 <- array(data = NA)
  for (u in 1:length(Samples)){
    
    data_in2 <- cbind(subset(data_in, select = Gene.Sets.All),
                      dplyr::select(data_in, Samples[u]))  #select one Sample type and the genes and Gene.Sets.A.and.B
    data_in2[,Samples[u]] <- as.numeric(as.character(data_in2[,Samples[u]]))
    data_in2 <- data_in2[order(-data_in2[,Samples[u]]),] #sort by descending order for the rank metric
    rownames(data_in2) <- 1:nrow(data_in2) #reorder row indices for counting in for loop below
    
    ## Assuming first two columns in data table are Genes and Rank Metric (e.g. Foldchange, SNR)
    
    GSEA.Results <- matrix(data = NA, nrow = length(Gene.Sets.All), ncol = 7)
    colnames(GSEA.Results) <- c("Sample","Gene.Set","KS","KS_Normalized",
                                "p_value","Position_at_max",
                                "FDR_q_value")
    GSEA.Results <- as.data.frame(GSEA.Results)
    GSEA.Results$Gene.Set <- Gene.Sets.All
    GSEA.Results$Sample <- Samples[u]
    
    ions <- nrow(data_in2)
    
    #for plotting
    ks_results_plot <- list()
    positions.of.hits <- list()
    
    #ks_results_plot <- as.data.frame(ks_results_plot)
    gene.list <- 1:ions
    rank_metric <- data_in2[,Samples[u]] #Save the rank metric
    
    pos_gene_set <- array(data = 0, dim = nrow(data_in2), dimnames = NULL);
    
    ## Calculate Real KS Statistic
    for (i in 1:length(Gene.Sets.All)){
      data_in3 <- data_in2[,Gene.Sets.All[i]]
      numhits_pathway <- sum(data_in3 == "X"); #check to see if there is anything in the column (e.g. X)
      if (numhits_pathway > 1){
        pos_gene_set <- which(data_in2[,Gene.Sets.All[i]] %in% c("X"))
        KS_real <- GSEA.EnrichmentScore(gene.list, pos_gene_set, weighted.score.type = score.weight, correl.vector = rank_metric)
        GSEA.Results[GSEA.Results$Gene.Set == Gene.Sets.All[i],]$KS <- KS_real$ES;
        GSEA.Results[GSEA.Results$Gene.Set == Gene.Sets.All[i],]$Position_at_max <- KS_real$arg.ES;
        ks_results_plot[[Gene.Sets.All[i]]] = KS_real$RES
        positions.of.hits[[Gene.Sets.All[i]]] = pos_gene_set
      }
    }
    
    Mountain.Plot.Info <- list(MountainPlot = ks_results_plot, Position.of.hits = positions.of.hits)
    rm(pos_gene_set)
    rm(numhits_pathway)
    rm(data_in3)
    rm(KS_real)
    
    #print("Calculating permutations...")
    
    #pb <- utils::txtProgressBar(max = num.permutations, style = 3)
    #progress <- function(n) utils::setTxtProgressBar(pb, n)
    #opts <- list(progress = progress)
    
    KSRandomArray <- matrix(data = NA, nrow = nperm, ncol = length(Gene.Sets.All))
    num.gene.sets.all <- length(Gene.Sets.All)
    `%dopar%` <- foreach::`%dopar%`
    
    #KSRandomArray <- foreach::foreach(L = 1:nperm, .combine = "rbind",.options.snow = opts) %dopar% {
    KSRandomArray <- foreach::foreach(L = 1:nperm, .combine = "rbind") %dopar% {
      temp.KSRandomArray <- matrix(data = NA, nrow = 1, ncol = num.gene.sets.all)
      for(i in 1:length(Gene.Sets.All)){
        numhits_pathway <- length(positions.of.hits[[Gene.Sets.All[i]]])
        pos_gene_set <- sample(1:ions,numhits_pathway)
        temp.KSRandomArray[,i] <- GSEA.EnrichmentScore2(gene.list, pos_gene_set, weighted.score.type = score.weight, correl.vector = rank_metric)
      }
      temp.KSRandomArray
    }
    colnames(KSRandomArray) <- Gene.Sets.All
    
    #rm(opts)
    #rm(pb)
    KSRandomArray <- data.frame(matrix(unlist(KSRandomArray), nrow = nperm, byrow = T))
    colnames(KSRandomArray) <- Gene.Sets.All
    KSRandomArray <- stats::na.omit(KSRandomArray)
    
    #print("Normalizing enrichment scores...")
    KSRandomArray <- as.data.frame(KSRandomArray)
    ###normalize the GSEA distribution
    KSRandomArray.Norm <- matrix(data = NA, nrow = nrow(KSRandomArray), ncol = ncol(KSRandomArray))
    colnames(KSRandomArray.Norm) <- colnames(KSRandomArray)
    avg <- 0
    KSRandomArray.temp <- 0
    for (i in 1:ncol(KSRandomArray.Norm)){
      avg <- 0
      KSRandomArray.temp <- KSRandomArray[,i]
      pos.temp <- KSRandomArray.temp[which(KSRandomArray.temp >= 0)]
      neg.temp <- KSRandomArray.temp[which(KSRandomArray.temp < 0)]
      
      avg.pos <- mean(pos.temp)
      avg.neg <- mean(neg.temp)
      
      norm.pos.temp <- pos.temp / avg.pos
      norm.neg.temp <- neg.temp / avg.neg * -1
      
      norm.perms <- c(norm.pos.temp,norm.neg.temp)
      
      KSRandomArray.Norm[,i] <- norm.perms
      
    }
    
    GSEA.NES.perms <- as.vector(KSRandomArray.Norm)
    rm(KSRandomArray.Norm)
    GSEA.NES.perms.pos <- GSEA.NES.perms[which(GSEA.NES.perms >= 0)]
    GSEA.NES.perms.neg <- GSEA.NES.perms[which(GSEA.NES.perms < 0)]
    rm(GSEA.NES.perms)
    percent.pos.GSEA <- sum(GSEA.Results$KS > 0) / length(GSEA.Results$KS)
    percent.neg.GSEA <- sum(GSEA.Results$KS < 0) / length(GSEA.Results$KS)
    
    # Calculate GSEA NES and p-value and FDR
    #print("Calculating GSEA FDR...")
    for (i in 1:length(Gene.Sets.All)){
      temp.gene.set <- Gene.Sets.All[i]
      temp.KS <- GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$KS
      if (temp.KS > 0){
        pos.perms <- KSRandomArray[,temp.gene.set]
        pos.perms <- pos.perms[which(pos.perms > 0)]
        #p-val
        GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$p_value = signif(sum(pos.perms > temp.KS) / length(pos.perms),digits = 3)
        #NES
        GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$KS_Normalized = signif(temp.KS / mean(pos.perms), digits = 3)
        #FDR
        percent.temp <- sum(GSEA.NES.perms.pos > GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$KS_Normalized) / length(GSEA.NES.perms.pos)
        GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$FDR_q_value = ifelse(signif(percent.temp / percent.pos.GSEA, digits = 3) < 1, signif(percent.temp / percent.pos.GSEA, digits = 3), 1)
      } else if (temp.KS < 0){
        neg.perms <- KSRandomArray[,temp.gene.set]
        neg.perms <- neg.perms[which(neg.perms < 0)]
        #p-val
        GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$p_value = signif(sum(neg.perms < temp.KS) / length(neg.perms),digits = 3)
        #NES
        GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$KS_Normalized = signif(temp.KS / mean(neg.perms) * -1, digits = 3)
        #FDR
        percent.temp <- sum(GSEA.NES.perms.neg < GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$KS_Normalized) / length(GSEA.NES.perms.neg)
        GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$FDR_q_value = ifelse(signif(percent.temp / percent.neg.GSEA, digits = 3) < 1, signif(percent.temp / percent.neg.GSEA, digits = 3), 1)
      }
    }
    
    GSEA.Results.All.Samples <- rbind(GSEA.Results.All.Samples,GSEA.Results)
    Mountain.Plot.Info.All.Samples <- c(Mountain.Plot.Info.All.Samples,Mountain.Plot.Info)
    rank_metric.All.Samples <- c(rank_metric.All.Samples,rank_metric)
    
    #print(paste("Sample #: ", u))
    
    #end.loop.time <- Sys.time()
    #total.loop.time <- signif(end.loop.time - loop.time, digits = 3)
    #print(paste("Time per Sample:" , total.loop.time))
  }
  
  snow::stopCluster(cl)
  rm(cl)
  
  return(list(GSEA.Results = GSEA.Results.All.Samples,
              Mountain.Plot.Info = Mountain.Plot.Info.All.Samples,
              ranking.metric = rank_metric.All.Samples))
  
}

# How far do we want to vary values here? Start with 0 to 1
values.to.vary <- seq(from = 0, to = 1, by = 0.05)


num.replicates <- 100
for (N in 1:num.replicates){
  loop.start.time <- Sys.time()

### Step 1: Generate gene expression matrix with a normal distribution for each cell line, then perturb 25 genes in gene set
###         Note: need to name the cell lines such that it's easy to track, we will want to perturb CERES in the same direction
#####    

  exp.matrix <- matrix(data = NA, nrow = length(genes), ncol = 301)
  row.names(exp.matrix) <- genes
  colnames(exp.matrix) <- c("Gene",synthetic.cell.names)
  exp.matrix[,"Gene"] <- genes
  
  # Fill out matrix with normal distribution of expression values
  # And replace synthetic gene set with a value of X
  
  for(i in 1:length(synthetic.cell.names)){
    temp.cell <- synthetic.cell.names[i]
    exp.matrix[,temp.cell] <- rnorm(num.genes, mean = 0, sd = 1)
  }
  
  expression.matrices <- list()
  
  ## Need to do 0 separately
  varied.value <- values.to.vary[1]
  varied.value.name <- paste("value_added_",varied.value, sep = "")
  for (k in 1:length(synthetic.cell.names)){
    exp.matrix[synthetic.gene.set, synthetic.cell.names[k]] <- rnorm(length(synthetic.gene.set),
                                                                     mean = 0, sd = 1)
  }
  expression.matrices[[varied.value.name]] <- exp.matrix
  
  for(j in 2:length(values.to.vary)){
    varied.value <- values.to.vary[j]
    varied.value.name <- paste("value_added_",varied.value, sep = "")
    
    range.to.vary <- seq(from = -varied.value, to = varied.value, by = varied.value / length(synthetic.cell.names) * 2)
    for (k in 1:length(synthetic.cell.names)){
      exp.matrix[synthetic.gene.set, synthetic.cell.names[k]] <- rnorm(length(synthetic.gene.set),
                                                                              mean = range.to.vary[k], sd = 1)
    }
    expression.matrices[[varied.value.name]] <- exp.matrix
  }



#####
### Step 2: Calculate GSEA NES for the synthetic gene set.
###

GSEA.matrices <- list()
  for (i in 1:length(names(expression.matrices))){
    temp.exp.matrix <- as.data.frame(expression.matrices[[i]])
    #temp.exp.matrix <- temp.exp.matrix[,1:21] #remove this later
    start <- Sys.time()
    temp.GSEA.results <- GSEA_custom(temp.exp.matrix, gmt.list = synthetic.gene.set.gmt)
    end <- Sys.time()
    total.time <- signif(end - start,digits = 3)
    print(paste("Time taken",total.time))
    NES.results <- temp.GSEA.results$GSEA.Results
    NES.results <- NES.results[NES.results$Gene.Set == "Synthetic_Gene_Set",]
    GSEA.matrices[[names(expression.matrices[i])]] <- NES.results
  }
  
#####
### Step 3: Generate CERES scores that have a normal distribution for each cell line, then perturb the same 25 genes
###         Note: since we will be doing spearman correlations, we don't have to worry about the range for CERES scores being not centered around 0.

#####
ceres.matrix <- matrix(data = NA, nrow = nrow(exp.matrix), ncol = (ncol(exp.matrix)-1))
row.names(ceres.matrix) <- genes
colnames(ceres.matrix) <- synthetic.cell.names

# Fill out matrix with normal distribution of expression values
# And replace synthetic gene set with a value of X

ceres.matrices <- list()

for(i in 1:length(synthetic.cell.names)){
  temp.cell <- synthetic.cell.names[i]
  ceres.matrix[,temp.cell] <- rnorm(num.genes, mean = 0, sd = 1)
}

values.to.vary <- seq(from = 0, to = 1, by = 0.05)

ceres.matrices <- list()


## Need to do 0 separately
varied.value <- values.to.vary[1]
varied.value.name <- paste("value_added_",varied.value, sep = "")
for (k in 1:length(synthetic.cell.names)){
  ceres.matrix[synthetic.gene.set, synthetic.cell.names[k]] <- rnorm(length(synthetic.gene.set),
                                                                   mean = 0, sd = 1)
}
ceres.matrices[[varied.value.name]] <- ceres.matrix

for(j in 2:length(values.to.vary)){
  varied.value <- values.to.vary[j]
  varied.value.name <- paste("value_added_",varied.value, sep = "")
  
  range.to.vary <- seq(from = -varied.value, to = varied.value, by = varied.value / length(synthetic.cell.names) * 2)
  for (k in 1:length(synthetic.cell.names)){
    ceres.matrix[synthetic.gene.set, synthetic.cell.names[k]] <- rnorm(length(synthetic.gene.set),
                                                                     mean = range.to.vary[k], sd = 1)
  }
  ceres.matrices[[varied.value.name]] <- ceres.matrix
}




#####

### Step 4: Generate correlation coefficients
###
#####
#exp.matrix.names <- names(expression.matrices)
ceres.matrix.names <- names(ceres.matrices)
GSEA.results.names <- names(GSEA.matrices)

correlation.matrices <- list()

for(i in 1:length(GSEA.results.names)){
  loop.start.time <- Sys.time()
  for (j in 1:length(ceres.matrix.names)){
    GSEA.value <- paste("Exp_",GSEA.results.names[i],sep = "")
    ceres.value <- paste("Ceres_",ceres.matrix.names[j], sep = "")
    
    corr.name <- paste(GSEA.value, ceres.value,"correlation", sep = "_")
    
    temp.GSEA.results <- (GSEA.matrices[[GSEA.results.names[i]]])
    temp.ceres.matrix <- t(ceres.matrices[[ceres.matrix.names[j]]])
    
    corr.matrix <- matrix(data = NA, nrow = num.genes, ncol = 2)
    colnames(corr.matrix) <- c("Gene", corr.name)
    row.names(corr.matrix) <- genes
    corr.matrix[,"Gene"] <- genes
    
    temp.GSEA.results.matrix <- as.matrix(temp.GSEA.results$KS_Normalized)
    
    for(k in 1:length(genes)){
      temp.ceres.gene <- temp.ceres.matrix[,genes[k], drop = FALSE]

      #merged <- merge(temp.exp.gene, temp.ceres.gene, by = 0) #don't need to merge bc these are already in order! nice
      corr.matrix[genes[k],corr.name] <- cor(temp.ceres.gene, temp.GSEA.results.matrix, method = "spearman")
    }
    correlation.matrices[[corr.name]] <- corr.matrix
  }
  # if(i == 1){
  #   loop.end.time <- Sys.time()
  #   loop.total.time <- loop.end.time - loop.start.time
  #   total.time.expected <- (length(GSEA.results.names) - 1) * loop.total.time
  #   current.time <- Sys.time()
  #   estimated.end.time <- current.time + total.time.expected
  #   print(paste("Estimated end time:", estimated.end.time))
  # }
}

#length(names(correlation.matrices))
# We have 441 dependency GSEAs to run

#####

##### Step 4: Run Dependency GSEA

corr.matrices.df <- data.frame(matrix(data = NA, nrow = num.genes, ncol = 1))
colnames(corr.matrices.df) <- c("Gene")
corr.matrices.df$Gene <- genes

for (i in 1:length(names(correlation.matrices))){
  temp.df <- correlation.matrices[[i]]
  
  corr.matrices.df <- merge(corr.matrices.df, temp.df, by = "Gene")
}


temp.df <- as.data.frame(correlation.matrices[[1]])

#start <- Sys.time()
GSEA.results <- GSEA_custom(temp.df, gmt.list = KEGG.pathways)

GSEA.results.compiled <- GSEA.results$GSEA.Results


for (i in 2:length(names(correlation.matrices))){
  temp.df <- as.data.frame(correlation.matrices[[i]])
  
  #start <- Sys.time()
  GSEA.results <- GSEA_custom(temp.df, gmt.list = KEGG.pathways)
  #end <- Sys.time()
  
  GSEA.results.compiled <- rbind(GSEA.results.compiled, GSEA.results$GSEA.Results)
}

GSEA.results.compiled$Simulation_Number <- N


GSEA.results.compiled.all.replicates <- rbind(GSEA.results.compiled.all.replicates,GSEA.results.compiled)


if(N == 1){
  loop.end.time <- Sys.time()
  loop.total.time <- loop.end.time - loop.start.time
  total.time.expected <- (length(exp.matrix.names) - 1) * loop.total.time
  current.time <- Sys.time()
  estimated.end.time <- current.time + total.time.expected
  print(paste("Time per loop:", loop.total.time))
  print(paste("Estimated end time:", estimated.end.time))
}

}

write.csv(GSEA.results.compiled.all.replicates, file = "Simulation_study_Genetic_PDEA.csv",row.names = F)


