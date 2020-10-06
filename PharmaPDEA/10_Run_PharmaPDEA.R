setwd("/.../Data/")
Drug.Pathway.Targets <- read.csv("Metabolic_Pathways_drug_targets.csv", stringsAsFactors = F)
PRISM.data.full <- data.table::fread("secondary-screen-dose-response-curve-parameters.csv",header=T)

setwd("/.../PharmaPDEA/")
corr.matrix <- read.csv(file = "Drug_corr_matrix_Adherent_RPMI.csv", header = TRUE)
output.filename <- "PharmaPDEA_Adherent_RPMI.csv"

# corr.matrix <- read.csv(file = "Drug_corr_matrix_Adherent_DMEM.csv", header = TRUE)
# output.filename <- "PharmaPDEA_Adherent_DMEM.csv"

ssGsea <- function(input.df, Gene.Sets,
                   num.permutations = 1000,
                   stat.type = "Weighted"){
  nperm = num.permutations #number of permutations
  
  if (stat.type == "Classic"){
    score.weight = 0
  }
  if (stat.type == "Weighted"){
    score.weight = 1
  }
  
  library(dplyr)
  library(foreach)
  library(doParallel)
  library(doSNOW)
  library(testthat); library(usethis)
  
  #Read in gene expression data
  #Genes should be first column, named "Gene"
  #Samples should be columns 2:N
  data_in <- input.df
  Gene.Sets <- Gene.Sets
  
  expect_is(data_in, "data.frame")
  expect_is(Gene.Sets, "data.frame")
  
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
  ### Annotate gene sets
  for (j in 1:length(Gene.Sets.All)){
    temp.pathway <- Gene.Sets[,Gene.Sets.All[j]]
    for (i in 1:nrow(annotations)){
      if (annotations[i,"Gene"] %in% temp.pathway){
        annotations[i,j+1] = "X";
      }
    }
  }
  annotations <- as.data.frame(annotations)
  data_in <- merge(data_in, annotations, by = "Gene")
  
  data_in <- na.omit(data_in)
  
  GSEA.Results.All.Samples <- matrix(data = NA, nrow = 0, ncol = 9)
  colnames(GSEA.Results.All.Samples) <- c("Sample","Gene.Set","KS","NES",
                              "p_value","Position_at_max", "Number_of_hits",
                              "FDR_q_value", "Leading_Edge_Genes")
  Mountain.Plot.Info.All.Samples <- list()
  rank_metric.All.Samples <- list()
  
  #Find out how many cores are available (if you don't already know)
  cores<-detectCores()
  #Create cluster with desired number of cores, leave one open for the machine         
  #core processes
  cl <- makeCluster(cores[1]-1)
  #Register cluster
  registerDoSNOW(cl)
  
  rm(annotations)
  data_in2 <- array(data = NA)
  for (u in 1:length(Samples)){
    loop.time <- Sys.time()
    
    data_in2 <- cbind(dplyr::select(data_in, Samples[u], "Gene"),subset(data_in, select = Gene.Sets.All)) 
    
    data_in2 <- data_in2[order(-data_in2[,Samples[u]]),] #sort by descending order for the rank metric
    rownames(data_in2) <- 1:nrow(data_in2) #reorder row indices for counting in for loop below
    
    ## Assuming first two columns in data table are Genes and Rank Metric (e.g. Foldchange, SNR)
    
    GSEA.Results <- matrix(data = NA, nrow = length(Gene.Sets.All), ncol = 9)
    colnames(GSEA.Results) <- c("Sample","Gene.Set","KS","NES",
                                "p_value","Position_at_max", "Number_of_hits",
                                "FDR_q_value", "Leading_Edge_Genes")
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
      if (numhits_pathway > 0){
        pos_gene_set <- which(data_in2[,Gene.Sets.All[i]] %in% c("X"))
        KS_real <- GSEA.EnrichmentScore(gene.list, pos_gene_set, weighted.score.type = score.weight, correl.vector = rank_metric)
        GSEA.Results[GSEA.Results$Gene.Set == Gene.Sets.All[i],]$KS <- KS_real$ES;
        GSEA.Results[GSEA.Results$Gene.Set == Gene.Sets.All[i],]$Position_at_max <- KS_real$arg.ES;
        GSEA.Results[GSEA.Results$Gene.Set == Gene.Sets.All[i],]$Number_of_hits <- numhits_pathway
        ks_results_plot[[Gene.Sets.All[i]]] = KS_real$RES
        positions.of.hits[[Gene.Sets.All[i]]] = pos_gene_set
        ### Find leading edge genes (e.g. genes before position at max)
        if (KS_real$ES > 0){
          leading.edge.positions <- pos_gene_set[which(pos_gene_set < KS_real$arg.ES)]
          leading.edge.genes <- data_in2[leading.edge.positions,"Gene"]
          GSEA.Results[GSEA.Results$Gene.Set == Gene.Sets.All[i],]$Leading_Edge_Genes <- paste(leading.edge.genes, collapse = ", ")
        } else if (KS_real$ES < 0){
          leading.edge.positions <- pos_gene_set[which(pos_gene_set > KS_real$arg.ES)]
          leading.edge.genes <- data_in2[leading.edge.positions,"Gene"]
          GSEA.Results[GSEA.Results$Gene.Set == Gene.Sets.All[i],]$Leading_Edge_Genes <- paste(leading.edge.genes, collapse = ", ")
        }
      }
    }
    
    Mountain.Plot.Info <- list(MountainPlot = ks_results_plot, Position.of.hits = positions.of.hits)
    rm(pos_gene_set)
    rm(numhits_pathway)
    rm(data_in3)
    rm(KS_real)
    
    rand_order_gene_set <- subset(data_in2, select = Gene.Sets.All)
    rand_order_gene_set <- as.matrix(rand_order_gene_set)
    ## Permutations using the same random shuffled gene order for each pathway
    
    
    print("Calculating permutations...")
    
    pb <- txtProgressBar(max = num.permutations, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    KSRandomArray <- matrix(data = NA, nrow = nperm, ncol = length(Gene.Sets.All))
    num.gene.sets.all <- length(Gene.Sets.All)
    KSRandomArray <- foreach(L = 1:nperm, .combine = "rbind",.options.snow = opts) %dopar% {
      temp.KSRandomArray <- matrix(data = NA, nrow = 1, ncol = num.gene.sets.all)
      for(i in 1:length(Gene.Sets.All)){
        numhits_pathway <- length(positions.of.hits[[Gene.Sets.All[i]]])
        pos_gene_set <- sample(1:ions,numhits_pathway)
        if (numhits_pathway > 0) {
          temp.KSRandomArray[,i] <- GSEA.EnrichmentScore2(gene.list, pos_gene_set, weighted.score.type = score.weight, correl.vector = rank_metric)
        }
      }
      temp.KSRandomArray
    }
    colnames(KSRandomArray) <- Gene.Sets.All
    
    rm(opts)
    rm(pb)
    KSRandomArray <- data.frame(matrix(unlist(KSRandomArray), nrow = nperm, byrow = T))
    colnames(KSRandomArray) <- Gene.Sets.All
    KSRandomArray <- na.omit(KSRandomArray)
    
    rm(rand_order_gene_set)
    print("Normalizing enrichment scores...")
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
    print("Calculating GSEA FDR...")
    for (i in 1:length(Gene.Sets.All)){
      temp.gene.set <- Gene.Sets.All[i]
      temp.KS <- GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$KS
      if (is.numeric(temp.KS)){
        if (temp.KS > 0){
          pos.perms <- KSRandomArray[,temp.gene.set]
          pos.perms <- pos.perms[which(pos.perms > 0)]
          #p-val
          GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$p_value = signif(sum(pos.perms > temp.KS) / length(pos.perms),digits = 3)
          #NES
          GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$NES = signif(temp.KS / mean(pos.perms), digits = 3)
          #FDR
          percent.temp <- sum(GSEA.NES.perms.pos > GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$NES) / length(GSEA.NES.perms.pos)
          if (is.numeric(percent.temp)){
            GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$FDR_q_value = ifelse(signif(percent.temp / percent.pos.GSEA, digits = 3) < 1, signif(percent.temp / percent.pos.GSEA, digits = 3), 1)
          } 
        }else if (temp.KS < 0){
          neg.perms <- KSRandomArray[,temp.gene.set]
          neg.perms <- neg.perms[which(neg.perms < 0)]
          #p-val
          GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$p_value = signif(sum(neg.perms < temp.KS) / length(neg.perms),digits = 3)
          #NES
          GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$NES = signif(temp.KS / mean(neg.perms) * -1, digits = 3)
          #FDR
          percent.temp <- sum(GSEA.NES.perms.neg < GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$NES) / length(GSEA.NES.perms.neg)
          if (is.numeric(percent.temp)){
            GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$FDR_q_value = ifelse(signif(percent.temp / percent.neg.GSEA, digits = 3) < 1, signif(percent.temp / percent.neg.GSEA, digits = 3), 1)
          }
        }
      }
    }
    rm(GSEA.NES.perms.pos)
    rm(GSEA.NES.perms.neg)
    
    GSEA.Results.All.Samples <- rbind(GSEA.Results.All.Samples,GSEA.Results)
    Mountain.Plot.Info.All.Samples <- c(Mountain.Plot.Info.All.Samples,Mountain.Plot.Info)
    rank_metric.All.Samples <- c(rank_metric.All.Samples,rank_metric)
    
    ####################################################
    ####### Below is only for data visualizaiton #######
    
    ## Create GSEA plot
    # Save default for resetting
    
    print(paste("Sample #: ", u))
    
    end.loop.time <- Sys.time()
    total.loop.time <- signif(end.loop.time - loop.time, digits = 3)
    print(paste("Time per Sample:" , total.loop.time))
  }
  
  stopCluster(cl)
  rm(cl)
  
  return(list( GSEA.Results = GSEA.Results.All.Samples, 
               Mountain.Plot.Info = Mountain.Plot.Info.All.Samples,
               ranking.metric = rank_metric.All.Samples))
  
}

### Need to adjust the correlation coefficients for activators


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

inhibitor.names <- PRISM.data.full[PRISM.data.full$moa %in% negative.regulators,]$name
activator.names <- PRISM.data.full[PRISM.data.full$moa %in% positive.regulators,]$name

rm(PRISM.data.full)


corr.matrix.adjusted <- corr.matrix
corr.matrix.adjusted[which(corr.matrix.adjusted$ï..Drug.Name %in% activator.names),3:72] <- corr.matrix.adjusted[which(corr.matrix.adjusted$ï..Drug.Name %in% activator.names),3:72]  * -1

colnames(corr.matrix.adjusted)[1] <- "Gene"
corr.matrix.adjusted$num.cell.lines <- NULL
colnames(corr.matrix.adjusted)[2:71] <- paste(colnames(corr.matrix.adjusted)[2:71], "Correlation", sep = "_")

DSEA.v1 <- ssGsea(corr.matrix.adjusted, Gene.Sets = Drug.Pathway.Targets)

DSEA.res <- DSEA.v1$GSEA.Results

DSEA.res.filtered <- DSEA.res[DSEA.res$Number_of_hits >= 4,]

write.csv(DSEA.res.filtered, file = output.filename ,row.names = FALSE)

