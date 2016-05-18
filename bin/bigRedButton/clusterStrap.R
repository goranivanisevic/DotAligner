#!/share/ClusterShare/software/contrib/gi/R/3.0.0/bin/R
#library(parallel)
library(snow)
library(pvclust)
library(ape)
source(paste(Sys.getenv("PATH_TO_SGE_SCRIPTS"),"pvrect.bigRedButton.R",sep="/"))
source(paste(Sys.getenv("PATH_TO_SGE_SCRIPTS"),"pvclust.bigRedButton.R",sep="/"))   
source(paste(Sys.getenv("PATH_TO_SGE_SCRIPTS"),"parPvclust.bigRedButton.R",sep="/"))
source(paste(Sys.getenv("PATH_TO_SGE_SCRIPTS"),"pvpick.bigRedButton.R",sep="/"))   

############
# Args from postAlign.sge
# $R_BIN ${PATH_TO_SGE_SCRIPTS}/clusterStrap.R ${1}/${2}/srtd.matrix ${1}/${2}/srtd.ids ${NSLOTS} ${BOOTSTRAPS} ${ALPHA_STAT} ${BETA_STAT} ${MAX_PER_CLUSTER}
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
# convert similarity matrix to distance matrix
distances <- as.dist(1 - Matrix.df)
############
print("[ R ] ...creating parallel environemnt")
cluster <- makeSOCKcluster(rep("localhost", args[3]))
clusterExport(cluster, ls())
############
print("[ R ] ...performing clustering and bootstrapping ")
# runs bootstrapping in parallel
boots <- as.numeric(args[4])
result <- parPvclust.bigRedButton( cluster, Matrix.df, distances, method.dist="cor", method.hclust="average", nboot=boots ) 
#result <- parPvclust( cluster, Matrix.df, method.dist="cor", method.hclust="average", nboot=boots)
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
beta <- as.numeric(args[6])
print("[ R ] ...creating pdf file")
out <- paste0(outfile,"_a",alpha,".pdf")
pdf( out ,width=ncol(Matrix.df)/3,height=10)
plot(result)
pvrect.bigRedButton(result, alpha=alpha )
dev.off()
############
print("[ R ] ...parsing pvclust output")
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
clusters <- pvpick.bigRedButton(result, alpha=alpha, beta=beta, type="geq", max.only=TRUE)
list2ascii(clusters$clusters, file=paste0(outfile,"_a",alpha,"_clusters"))
############
print("[ R ] ...extracting significant clusters")
phy <- as.phylo( result$hclust )
x=as.numeric(args[7]) # max sequences per cluster
counter=1
for (i in clusters$clusters) { 
	# remove all leaves not in cluster form the tree	
	############
	cluster_tree <- drop.tip( phy, subset( phy$tip.label, !( phy$tip.label %in% i ) ) ) 
	###### IF cluster size > x (8 by default)
	if ( cluster_tree$Nnode > x ) {
		print("[ R ] ...cluster is too big, selecting most similar alignments" )
		# matrix of all distances between  terminal nodes
		tree_dists <- cophenetic.phylo( cluster_tree ) 
		# Gets row (first cell) of matrix with lowest non-null value
		top <- which( tree_dists ==  min( tree_dists[ tree_dists > 0 ] ) , arr.ind = TRUE )[1,1]
		# Get the (x+1)th least similar rownames
		most_similar <- order( tree_dists[top,] )[ 1 : x ]
		# Get their names
		to_keep <- rownames( tree_dists )[ most_similar ]
		#write the sub-cluster to a file
		write( to_keep, file = paste0("clusters/sub_cluster_",counter,".ids")
		# remove all the other tips from the cluster
		sub_cluster <- drop.tip( cluster_tree, subset( cluster_tree$tip.label, !( cluster_tree$tip.label %in% to_keep ) ) ) 
		write( write.tree( cluster_tree ), file = paste0("clusters/cluster_",counter,".newick") )
		write( write.tree( sub_cluster ), file = paste0("clusters/sub_cluster_",counter,".newick") )
	} else {
		print("[ R ] ...cluster is not too big, writing directly to file" )
		write( write.tree( cluster_tree ), file = paste0("clusters/cluster_",counter,".newick") )
	}
	counter=counter+1
}
#extract.clade #does the inverse operation: it keeps all the tips from a given node, and deletes all the other tips.
print("[ R ] Done! ")
q(save="no")