parPvclust.bigRedButton <- function (cl, data, data.dist, method.hclust = "average", method.dist = "correlation", use.cor = "pairwise.complete.obs", nboot = 1000, r = seq(0.5, 1.4, by = 0.1), store = FALSE, weight = FALSE, init.rand = TRUE, seed = NULL)
{
  if (!(require(snow)))
    stop("Package snow is required for parPvclust.")
  if ((ncl <- length(cl)) < 2 || ncl > nboot) {
    warning("Too small value for nboot: non-parallel version is executed.")
    return(pvclust.bigRedButton(data, data.dist, method.hclust, method.dist, use.cor, nboot, r, store))
  }
  if (init.rand) {
    if (is.null(seed))
      seed <- 1:length(cl)
    else if (length(seed) != length(cl))
      stop("seed and cl should have the same length.")
      parLapply(cl, as.list(seed), set.seed)
    }
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
    nbl <- as.list(rep(nboot%/%ncl, times = ncl))
    if ((rem <- nboot%%ncl) > 0)
      nbl[1:rem] <- lapply(nbl[1:rem], "+", 1)
    cat("Multiscale bootstrap... ")
    mlist <- parLapply(cl, nbl, pvclust:::pvclust.node, r = r, data = data, object.hclust = data.hclust, method.dist = method.dist, use.cor = use.cor, method.hclust = method.hclust, store = store, weight = weight)
    cat("Done.\n")
    mboot <- mlist[[1]]
    for (i in 2:ncl) {
      for (j in 1:rl) {
        mboot[[j]]$edges.cnt <- mboot[[j]]$edges.cnt + mlist[[i]][[j]]$edges.cnt
        mboot[[j]]$nboot <- mboot[[j]]$nboot + mlist[[i]][[j]]$nboot
        mboot[[j]]$store <- c(mboot[[j]]$store, mlist[[i]][[j]]$store)
      }
    }
    result <- pvclust:::pvclust.merge(data = data, object.hclust = data.hclust, mboot = mboot)
    return(result)
}
