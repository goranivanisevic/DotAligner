#Plot interactive 3d plot 
#
#Input example:
# ./DotAligner -d ../example/test4.pp -d ../example/test1.pp -k 0.5 -a -1 -b -0.1 --local1 &> testit
# awk '/^Local/{print $5, $7, $9}' testit > testit2
#
#how to launch in R:
# commandArgs <- function() "/media/data/phd/projects/DotAligner/bin/testit2"
# source("draw3d.R")

input<-commandArgs()

t<-read.table(input)
#library(scatterplot3d)
library(rgl)
plot3d(t[,1], t[,2], t[,3], col="red", size=3)
plot3d(t[,1], t[,2], t[,3], col="red", size=3, xlab="reference base of S_a", ylab="reference base of S_b", zlab="Local similarity")
