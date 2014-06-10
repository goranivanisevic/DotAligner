pvpick.bigRedButton <- function (x, alpha = 0.95, beta = 0.80, type = "geq", max.only = TRUE)
{
  len <- nrow(x$edges)
  member <- pvclust:::hc2split(x$hclust)$member
  order <- x$hclust$order
  ht <- c()
  a <- list(clusters = list(), edges = c())
  j <- 1
  if (is.na(pm <- pmatch(type, c("geq", "leq", "gt", "lt")))) 
    stop("Invalid type argument: see help(pickup)")
  for (i in (len - 1):1) {
    if (pm == 1)
      wh <- (x$edges[i, "au"] >= alpha && x$edges[i, "bp"] >= beta)
    else if (pm == 2)
      wh <- (x$edges[i, "au"] <= alpha && x$edges[i, "bp"] <= beta)
    else if (pm == 3)
      wh <- (x$edges[i, "au"] > alpha && x$edges[i, "bp"] > beta)
    else if (pm == 4)
      wh <- (x$edges[i, "au"] < alpha && x$edges[i, "bp"] < beta)
    if (wh) {
      mi <- member[[i]]
      ma <- match(mi, order)
      if (max.only == FALSE || (max.only && sum(match(ma, ht, nomatch = 0)) == 0)) {
        a$clusters[[j]] <- x$hclust$labels[mi]
        a$edges <- c(a$edges, i)
        j <- j + 1
      }
      ht <- c(ht, ma)
    }
  }
  a$edges <- a$edges[length(a$edges):1]
  a$clusters <- a$clusters[length(a$edges):1]
  return(a)
}
  
