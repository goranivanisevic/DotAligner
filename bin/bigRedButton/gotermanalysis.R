## gotermanalysis.R <INPUT FILE> <OUTPUT FILE>
#  GOterm enrichment test
args <- commandArgs(trailingOnly = TRUE)
input <- as.character(args[1])
output <- as.character(args[2])

#source("http://bioconductor.org/biocLite.R")
library(org.Hs.eg.db)
library(GO.db)
library("GOstats")
ens.ids <- readLines(input)
ens.ids.entrez <- unique(unlist(mget(ens.ids, org.Hs.egENSEMBL2EG, ifnotfound=NA)))
entrez_object <- org.Hs.egGO
universe <- mappedkeys(entrez_object)
params <- new('GOHyperGParams',geneIds=ens.ids.entrez,universeGeneIds=universe,ontology='BP',pvalueCutoff=0.001,conditional=F,testDirection='over',annotation="org.Hs.eg.db")
hgOver <- hyperGTest(params)
write.table(summary(hgOver), output, sep = "\t")


