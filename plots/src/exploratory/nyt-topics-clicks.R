# set the directory with the model fit 
#dir <- "/n/fs/learning/run/nyt/000"
dir <- "../output/nyt/bpf.hier.bin"

# set the number of factors
K <- 100

#
# SET THE THRESHOLD TO APPLY ON TOPICS
#
threshold <- 0.01


# global objects
nusers <- 1615675
ntitles <- 107523
ndocs <- 107523

etitles <- data.frame(read.table("nyt-titles.tsv", sep="|"))
colnames(etitles) <- c("docid", "title")
x <- read.table("str2id.tsv")
colnames(x) <- c("docid", "seq")
titles <- merge(etitles, x, by="docid")

#
# CHANGE beta.tsv to hbeta.tsv for HIERARCHICAL MODEL
#
betaf <- sprintf ("%s/hbeta.tsv", dir)
beta <- data.frame(read.table(betaf))
colnames(beta) <- c("seq0", "seq")

theta.fname <- sprintf("%s/htheta.tsv", dir)
theta <- read.table(theta.fname)
colnames(theta) <- c("seq", "userid", c(1:K))

rf <- sprintf("%s/ranking.tsv", dir)
r <- data.frame(read.table(rf))
colnames(r) <- c("userid", "seq", "pred", "rating")

users.fname <- "../data/nyt/train.tsv"
users <- read.table(users.fname)
colnames(users) <- c("userid", "seq", "rating")
articles <- merge(users, titles, by="seq")

show_top_articles <- function(u) {
    articles.of.user <- subset(articles, articles$userid == u)
    a <- articles.of.user[order(-articles.of.user$rating)[1:100],]
    print(na.omit(a[c("title","rating")]), row.names=FALSE)
}


beta.with.titles <- merge(beta, titles, by="seq", all.x=TRUE)
beta.with.titles.b <- merge(beta, titles, by="seq")
colnames(beta.with.titles) <- c("seq", "seq0", c(1:K), "docid", "title")
colnames(beta.with.titles.b) <- c("seq", "seq0", c(1:K), "docid", "title")
#save(beta.with.titles, file="beta.with.titles.rda")

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
  df <- data.frame(cbind(c(1:K),
        as.matrix(t(subset(theta, theta$userid == u)))[3:(K+2)]))
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
      topics(c, 0)
      n <- n + 1
      if (n > 2) {
        break;
      }
  }
}

show_related_movies <- function(moviename)
{
  mname <- paste(moviename, "(,|\\s+).*", sep="")
  itemid <- grep(mname, titles$title, perl=TRUE)
  if (length(itemid) > 1) {
    itemid <- itemid[1]
  }
  type <- as.character(titles[itemid,"type"])
  itemid <- as.matrix(titles[itemid,]["itemid"])
  itemid <- as.integer(itemid[1][1])
  betam <- as.matrix(beta[, -c(1,2)])
  movseq <- which(beta$itemid == itemid)
  printf("item = %d, title = %s, type = %s", 
         itemid, moviename, type)
  df <- data.frame(cbind(c(1:100), betam[movseq,]))
  colnames(df) <- c("seq", "beta")
  top.mov.seq <- df$seq[order(df$beta, decreasing=T)] 
  for (c in top.mov.seq[1:3]) {
    top_movies_by_factor(c)
  }
  barplot(betam[movseq,])
}

top_movies_by_factor <- function(x)
{ 
  printf("FACTOR %d", x)
  factorname <- sprintf("K%d",x)
  a <- merge(beta, titles, by="itemid")
  a <- a[order(-a[[factorname]])[1:10],]
  print(na.omit(a[c("title", "type")]), row.names=FALSE)
}



show_top_articles(10)

quit()


a <- sprintf("%s/topics.txt", dir)
sink(a)
alltopics(0)
sink(NULL)

b <- sprintf("%s/subset-topics3.txt", dir)
sink(b)
alltopics(threshold)
sink(NULL)
