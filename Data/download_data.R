setwd("/.../Data/")

#Gene Dependency
download.file(url = "https://ndownloader.figshare.com/files/22543640", destfile = "Achilles_gene_effect.csv")

#Expression
download.file(url = "https://depmap.org/portal/download/api/download?file_name=ccle%2Fccle_2019%2FCCLE_RNAseq_genes_rpkm_20180929.gct.gz&bucket=depmap-external-downloads",
              destfile = "CCLE_expression_19Q4.csv")

#PRISM drug response
download.file(url = "https://ndownloader.figshare.com/files/20237739",
              destfile = "secondary-screen-dose-response-curve-parameters.csv")

#Cell line info
download.file(url = "https://ndownloader.figshare.com/files/20274744",
              destfile = "sample_info.csv")
