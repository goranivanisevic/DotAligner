pvrect.bigRedButton <- function (x, alpha = 0.95, beta = 0.80, type = "geq", max.only = TRUE, border = 2, ...)
{
  len <- nrow(x$edges)
  member <- pvclust:::hc2split(x$hclust)$member
  order <- x$hclust$order
  usr <- par("usr")
  xwd <- usr[2] - usr[1]
  ywd <- usr[4] - usr[3]
  cin <- par()$cin
  ht <- c()
  j <- 1
  if (is.na(pm <- pmatch(type, c("geq", "leq", "gt", "lt"))))
    stop("Invalid type argument: see help(pvrect)")
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
        xl <- min(ma)
        xr <- max(ma)
        yt <- x$hclust$height[i]
        yb <- usr[3]
        mx <- xwd/length(member)/3
        my <- ywd/200
        rect(xl - mx, yb + my, xr + mx, yt + my, border = border, shade = NULL, ...)
        j <- j + 1
      }
      ht <- c(ht, ma)
    }
  }
}

