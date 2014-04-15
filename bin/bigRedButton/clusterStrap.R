#!/share/ClusterShare/software/contrib/gi/R/3.0.0/bin/R
#library(parallel)
library(snow)
library(pvclust)
############
args <- commandArgs(TRUE)
print("[ R ] processed args")
Matrix<-read.table(args[1])
Matrix.df<-as.data.frame.matrix(Matrix)
print("[ R ] loaded MAtrix data frame")
############
Ids<-read.table(args[2])[,1]
print("[ R ] loaded IDs")
############
colnames(Matrix.df)<-Ids
#args[3] = # CPUs
############
print("[ R ] Loading cluster")
cluster <- makeSOCKcluster(rep("localhost", args[3]))
clusterExport(cluster, ls())
############
print("[ R ] Clustering and bootstrapping ")
# runs bootstrapping in parallel
result <- parPvclust( cluster, Matrix.df, method.dist="cor", method.hclust="average", nboot=10000)
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
out99 <- paste(outfile,"_a99.pdf", sep="")
pdf( out99 ,width=ncol(Matrix.df),height=15)
plot(result)
pvrect(result, alpha=0.99)
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
list2ascii(pvpick(result, alpha=0.99, pv="au", type="geq", max.only=TRUE)$clusters, file=paste0(outfile,"_a99_clusters"))
############
print("[ R ] Creating newick trees from pvclust output")
library(ape)
result_phylo <- as.phylo(result$hclust)
write(write.tree(result_phylo), file = paste(outfile,".newick", sep="") )
print("[ R ] Done! ")
q(save="no")