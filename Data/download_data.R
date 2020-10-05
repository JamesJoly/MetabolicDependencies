setwd("/.../Data/")

setwd("C:/Users/James/Documents/GitHub/MetabolicDependencies/Data")
#Expression
download.file(url = "https://ndownloader.figshare.com/files/20234346",
              destfile = "CCLE_expression_19Q4.csv",
              method = "libcurl", mode = "wb")

#Gene Dependency
download.file(url = "https://ndownloader.figshare.com/files/22543640", 
              destfile = "Achilles_gene_effect.csv",
              method = "libcurl", mode = "wb")

#PRISM drug response
download.file(url = "https://ndownloader.figshare.com/files/20237739",
              destfile = "secondary-screen-dose-response-curve-parameters.csv",
              method = "libcurl", mode = "wb")

#Cell line info
download.file(url = "https://ndownloader.figshare.com/files/20274744",
              destfile = "sample_info.csv", mode = "wb")
