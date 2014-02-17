
topics <- function(factor, threshold)
{
  c <- sprintf("%d", factor)
  t <- beta.with.titles[order(-beta.with.titles[[c]])[1:20],]
  s <- beta.with.titles.b[order(-beta.with.titles.b[[c]])[1:20],]
  if (all(beta.with.titles[[c]] == 0)) {
    return()
  } 
  if (t[1,c] < threshold) {
   return()
  }  
  buf <- sprintf("TOPIC %d", factor)
  print(buf)
  print(na.omit(s[c("title", factor)]), row.names=FALSE)
}

alltopics <- function(threshold)
{
  for (i in c(1:K)) {
    topics(i, threshold)
  }
}

show_top_user_factors <- function(u)
{
  df <- data.frame(cbind(c(1:K), as.matrix(t(subset(theta, theta$userid == u)))[3:(K+2)]))
  colnames(df) <- c("fac", "theta")
  print(df)
  betam <- as.matrix(beta[, -c(1,2)])
  facs <- as.matrix(df$fac[order(df$theta, decreasing=T)[1:100]])
  n <- 0
  for (c in facs[1:100]) {
      # skip unused factor
      q <- length(which(as.logical(beta[[c]])))
      if (q <= 1) {
        next;
      }
      topics(c, 0.01)
      n <- n + 1
      if (n > 4) {
        break;
      }
  }
}
