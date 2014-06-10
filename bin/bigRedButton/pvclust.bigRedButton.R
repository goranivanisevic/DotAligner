#library(pvclust)
#source("/home/stesee/DotAligner/example/rfam/temp.clustering/pvrect.bigRedButton.R")
#source("/home/stesee/DotAligner/example/rfam/temp.clustering/pvclust.bigRedButton.R")
#fix(pvclust)
#Matrix<-read.table("../dotaligner/srtd.matrix")
#Ids<-read.table("../dotaligner/srtd.ids")[,1]
#Matrix.df<-as.data.frame.matrix(Matrix)
#colnames(Matrix.df)<-Ids
#rownames(Matrix.df)<-Ids
#Matrix.dist<-as.dist(Matrix.df)
#getAnywhere("dist.pvclust")
#res <- as.dist(1 - cor(Matrix.df, method = "pearson", use = "pairwise.complete.obs"))
#res <- as.dist(1 - Matrix.df)
#result <- pvclust.bigRedButton( Matrix.df, res, method.dist="cor", method.hclust="average", nboot=100)
#pdf( "temp.cluster4.pdf" ,width=ncol(Matrix.df)/3,height=10)
#plot(result)
#pvrect.bigRedButton(result, alpha=0.8, beta=0.4, type="geq", max.only=TRUE )
#pvrect(result, alpha=0.99, pv="au", type="geq", max.only=TRUE )
#pvrect(result, alpha=0.5, pv="bp", type="geq", max.only=TRUE )
#dev.off()

pvclust.bigRedButton <- function (data, data.dist, method.hclust = "average", method.dist = "correlation", use.cor = "pairwise.complete.obs", nboot = 1000, r = seq(0.5, 1.4, by = 0.1), store = FALSE, weight = FALSE)
{
  n <- nrow(data)
  p <- ncol(data)
  METHODS <- c("ward", "single", "complete", "average", "mcquitty", "median", "centroid")
  method.hclust <- METHODS[pmatch(method.hclust, METHODS)]
  data.hclust <- hclust(data.dist, method = method.hclust)
  size <- floor(n * r)
  rl <- length(size)
  if (rl == 1) {
    if (r != 1)
      warning("Relative sample size r is set to 1.0. AU p-values are not calculated\n")
    r <- list(1)
  }
  else r <- as.list(size/n)
  mboot <- lapply(r, pvclust:::boot.hclust, data = data, object.hclust = data.hclust,
  nboot = nboot, method.dist = method.dist, use.cor = use.cor,
  method.hclust = method.hclust, store = store, weight = weight)
  result <- pvclust:::pvclust.merge(data = data, object.hclust = data.hclust,
  mboot = mboot)
  return(result)
}


