library(ggplot2)
library(plyr)

topdir <- "nyt"
K <- 100

theme_set(theme_bw())
datdir <- "../output/nyt/bpf.bias.hier.bin2";
titles <- data.frame(read.table("../data/nyt/nyt-titles.tsv", sep="|"))
colnames(titles) <- c("id", "title")

byusers.fname <- sprintf ("../data/%s/byusers.tsv", topdir)
byitems.fname <- sprintf ("../data/%s/byitems.tsv", topdir)
byusers <- read.table(byusers.fname)
colnames(byusers) <- c("seq", "id", "count", "tratings")
byitems <- read.table(byitems.fname)
colnames(byitems) <- c("seq", "id", "count", "tratings")

user.rate.fname <- sprintf("%s/thetarate.tsv", datdir)
item.rate.fname <- sprintf("%s/betarate.tsv", datdir)
user.rate <- read.table(user.rate.fname)
colnames(user.rate) <- c("seq", "id", "val")
item.rate <- read.table(item.rate.fname)
colnames(item.rate) <- c("seq", "id", "val")

user.rate <- merge(byusers, user.rate, by=c("seq", "id"))
item.rate <- merge(byitems, item.rate, by=c("seq", "id"))

beta.fname <- sprintf("%s/hbeta.tsv", datdir)
beta <- read.table(beta.fname)
colnames(beta) <- c("seq", "id", c(1:K))

theta.fname <- sprintf("%s/htheta.tsv", datdir)
theta <- read.table(theta.fname)
colnames(theta) <- c("seq", "id", c(1:K))

beta <- cbind(beta[1:2], beta[,-c(1,2)] / rowSums(beta[,-c(1,2)]))
beta <- merge(item.rate, beta, by="id")
beta <- merge(beta, titles, by="id")

q <- beta[order(1/beta$val, decreasing=TRUE),]
r <- beta[order(1/beta$val, decreasing=FALSE),]

y1 <- q[1:100,]
z1 <- cbind(y1, rowSums(y1[,s] > 0.1))
t1 <- z1[,-6]
colnames(t1) <- c("id", "seq", "count", "tratings", "val", c(1:K), "title", "memberships")
print(t1[,c("memberships", "tratings", "count", "val", "title")])

y2 <- r[1:100,]
z2 <- cbind(y2, rowSums(y2[,s] > 0.1))
t2 <- z2[,-6]
colnames(t2) <- c("id", "seq", "count", "tratings", "val", c(1:K), "title", "memberships")
print(t2[,c("memberships", "tratings", "count", "val", "title")])


topics <- function(factor, threshold)
{
    c <- sprintf("%d", factor)
    t <- beta[order(-beta[[c]])[1:20],]
    if (all(beta[[c]] == 0)) {
        return()
    }
    if (t[1,c] < threshold) {
        return()
    }
    buf <- sprintf("TOPIC %d", factor)
    print(buf)
    print(na.omit(t[c("title", factor)]), row.names=FALSE)    
}

alltopics <- function(threshold)
{
   for (i in c(1:K)) {
    topics(i, threshold)
   }
}

threshold <- 0.3
b <- sprintf("%s/subset-topics2.txt", datdir)
sink(b)
alltopics(threshold)
sink(NULL)
    

quit()

# -- scratch --
#theta <- merge(user.rate, theta, by="id")
#theta.fname <- sprintf("%s/htheta.tsv", datdir)
#theta <- read.table(theta.fname)
#colnames(theta) <- c("seq", "id", c(1:K))
#t(apply(beta,1,norm<-function(x){return (x/sum(x[c(7:K+7)]))}))
