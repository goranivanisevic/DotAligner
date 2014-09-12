#!/share/ClusterShare/software/contrib/gi/R/3.0.0/bin/R
#library(parallel)
library(snow)
library(pvclust)
############
args <- commandArgs(TRUE)
# Use fread() instead, much faster. Make sure the data is clean! (data.table, plyr and dplyr libraries)
print("[ R ] ...loading Matrix as data.frame object")
Matrix<-read.table(args[1])
Matrix.df<-as.data.frame.matrix(Matrix)
############
print("[ R ] ...loading IDs")
Ids<-read.table(args[2])[,1]
############
colnames(Matrix.df)<-Ids
#args[3] = # CPUs
############
print("[ R ] ...creating parallel environemnt")
cluster <- makeSOCKcluster(rep("localhost", args[3]))
clusterExport(cluster, ls())
############
print("[ R ] ...performing clustering and bootstrapping ")
# runs bootstrapping in parallel
boots <- as.numeric(args[4])
result <- parPvclust( cluster, Matrix.df, method.dist="cor", method.hclust="average", nboot=boots)
stopCluster(cluster)
# print out dendrograms
print("[ R ] ...preparing ouput")
outfile <- gsub(".matrix","",args[1])
############
# out95 <- paste(outfile,"_a95.pdf", sep="")
# pdf( out95 ,width=ncol(Matrix.df),height=15)
# plot(result)
# pvrect(result, alpha=0.95)
# dev.off()
############
alpha <- as.numeric(args[5])
print("[ R ] ...creating pdf file")
out <- paste0(outfile,"_a",alpha,".pdf")
pdf( out ,width=ncol(Matrix.df)/3,height=10)
plot(result)
pvrect(result, alpha=alpha )
dev.off()
############
print("[ R ] ...extracting significant clusters")
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
list2ascii(pvpick(result, alpha=alpha, pv="au", type="geq", max.only=TRUE)$clusters, file=paste0(outfile,"_a",alpha,"_clusters"))
############
print("[ R ] ...creating newick trees from pvclust output")
library(ape)
result_phylo <- as.phylo(result$hclust)
write(write.tree(result_phylo), file = paste0(outfile,".newick") )
print("[ R ] Done! ")
q(save="no")
