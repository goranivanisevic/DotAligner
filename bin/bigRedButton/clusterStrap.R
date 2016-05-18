#!/share/ClusterShare/software/contrib/gi/R/3.0.0/bin/R
#library(parallel)
library(snow)
library(pvclust)
source(paste(Sys.getenv("PATH_TO_SGE_SCRIPTS"),"pvrect.bigRedButton.R",sep="/"))   ###NEW
source(paste(Sys.getenv("PATH_TO_SGE_SCRIPTS"),"pvclust.bigRedButton.R",sep="/"))   ###NEW
source(paste(Sys.getenv("PATH_TO_SGE_SCRIPTS"),"parPvclust.bigRedButton.R",sep="/"))   ###NEW
source(paste(Sys.getenv("PATH_TO_SGE_SCRIPTS"),"pvpick.bigRedButton.R",sep="/"))   ###NEW

############
args <- commandArgs(TRUE)
print("[ R ] processed args")
# Use fread() instead, much faster. Make sure the data is clean! (data.table, plyr and dplyr libraries)
Matrix<-read.table(args[1])
Matrix.df<-as.data.frame.matrix(Matrix)
print("[ R ] loaded Matrix data frame")
############
Ids<-read.table(args[2])[,1]
print("[ R ] loaded IDs")
############
colnames(Matrix.df)<-Ids
rownames(Matrix.df)<-Ids  ###NEW
#args[3] = # CPUs
############
print("[ R ] Loading cluster")
cluster <- makeSOCKcluster(rep("localhost", args[3]))
clusterExport(cluster, ls())
############
print("[ R ] Clustering and bootstrapping ")
# convert similarity matrix to distance matrix
distances <- as.dist(1 - Matrix.df)
# runs bootstrapping in parallel
boots <- as.numeric(args[4])
#result <- pvclust.bigRedButton( Matrix.df, distances, method.dist="cor", method.hclust="average", nboot=boots)
result <- parPvclust.bigRedButton( cluster, Matrix.df, distances, method.dist="cor", method.hclust="average", nboot=boots)   ###NEW
#result <- parPvclust( cluster, Matrix.df, method.dist="cor", method.hclust="average", nboot=boots)
stopCluster(cluster)
# print out dendrograms
outfile <- gsub(".matrix","",args[1])
############
# out95 <- paste(outfile,"_a95.pdf", sep="")
# pdf( out95 ,width=ncol(Matrix.df),height=15)
# plot(result)
# pvrect(result, alpha=0.95)
# dev.off()
############
alpha <- as.numeric(args[5])
beta <- as.numeric(args[6])
out <- paste0(outfile,"_a",alpha,".pdf")
pdf( out ,width=ncol(Matrix.df)/3,height=10)
plot(result)
pvrect.bigRedButton(result, alpha=alpha, beta=beta, type="geq", max.only=TRUE )   ###NEW
dev.off()
############
print("[ R ] Extracting Clusters")
list2ascii <- function(x,file=paste(deparse(substitute(x)),".txt",sep="")) {
	tmp.wid = getOption("width")  # save current width
	options(width=10000)          # increase output width
	sink(file)                    # redirect output to file
	print(x)                      # print the object
	sink()                        # cancel redirection
	options(width=tmp.wid)        # restore linewidth
	return(invisible(NULL))       # return (nothing) from function
}
#list2ascii(pvpick(result, alpha=0.95, pv="au", type="geq", max.only=TRUE)$clusters, file=paste0(outfile,"_a95_clusters"))
list2ascii(pvpick.bigRedButton(result, alpha=alpha, beta=beta, type="geq", max.only=TRUE)$clusters, file=paste0(outfile,"_a",alpha,"_clusters"))   ###NEW
#list2ascii(pvpick(result, alpha=alpha, pv="au", type="geq", max.only=TRUE)$clusters, file=paste0(outfile,"_a",alpha,"_clusters"))
############
print("[ R ] Creating newick trees from pvclust output")
library(ape)
result_phylo <- as.phylo(result$hclust)
write(write.tree(result_phylo), file = paste0(outfile,".newick") )
print("[ R ] Done! ")
q(save="no")
