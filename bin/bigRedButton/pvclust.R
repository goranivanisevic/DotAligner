#!/usr/bin/R
library(pvclust)
seed_56_95.new.nw<-read.table("/home/stesee/DotAligner/analysis/nw/seed_56_95.new.nw.2x2")
seed_56_95.new.nw.df<-as.data.frame.matrix(seed_56_95.new.nw)
identifier<-read.table("/home/stesee/DotAligner/analysis/nw/seed_56_95.new.identifier")[,1]
colnames(seed_56_95.new.nw.df)<-identifier
result <- pvclust(seed_56_95.new.nw.df, method.dist="cor", method.hclust="average", nboot=10000)
pdf("/home/stesee/DotAligner/analysis/nw/seed_56_95.new.nw.bootstrap10000.pdf",width=180,height=10)
plot(result)
pvrect(result, alpha=0.95)
dev.off()
list2ascii <- function(x,file=paste(deparse(substitute(x)),".txt",sep="")) {
 tmp.wid = getOption("width")  # save current width
 options(width=10000)          # increase output width
 sink(file)                    # redirect output to file
 print(x)                      # print the object
 sink()                        # cancel redirection
 options(width=tmp.wid)        # restore linewidth
 return(invisible(NULL))       # return (nothing) from function
}
list2ascii(pvpick(result, alpha=0.95, pv="au", type="geq", max.only=TRUE)$clusters,file="/home/stesee/DotAligner/analysis/nw/seed_56_95.new.nw.clusters")

