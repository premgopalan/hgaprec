library(ggplot2)
library(plyr)

theme_set(theme_bw())

analyze <- function(topdir, bin, bsize) {
  if (bin == 1) {
      fname <- "bpf.hier.bin"
  } else {
      fname <- "bpf.hier"
  }
  xdir <- sprintf("../output/%s/%s", topdir, fname)

  titles <- data.frame(read.table("../data/nyt/nyt-titles.tsv", sep="|"))
  colnames(titles) <- c("id", "title")

  byusers.fname <- sprintf ("../data/%s/byusers.tsv", topdir)
  byitems.fname <- sprintf ("../data/%s/byitems.tsv", topdir)
  byusers <- read.table(byusers.fname)
  colnames(byusers) <- c("seq", "id", "count", "tratings")
  byitems <- read.table(byitems.fname)
  colnames(byitems) <- c("seq", "id", "count", "tratings")

  user.rate.fname <- sprintf("%s/thetarate.tsv", xdir)
  item.rate.fname <- sprintf("%s/betarate.tsv", xdir)
  user.rate <- read.table(user.rate.fname)
  colnames(user.rate) <- c("seq", "id", "val")
  
  item.rate <- read.table(item.rate.fname)
  colnames(item.rate) <- c("seq", "id", "val")
  
  user.rate <- merge(byusers, user.rate, by=c("seq", "id"))
  item.rate <- merge(byitems, item.rate, by=c("seq", "id"))

  i.rate.plot <- ggplot()
  u.rate.plot <- ggplot()
  i.rate.plot2 <- ggplot()
  u.rate.plot2 <- ggplot()
  
  if (bin == 1) {
      
      user.rate.sum <- ddply(user.rate, c("count"), summarize, mean.val=mean(1/val), stdev.val=sqrt(var(1/val)), freq=length(count))
      item.rate.sum <- ddply(item.rate, c("count"), summarize, mean.val=mean(1/val), stdev.val=sqrt(var(1/val)), freq=length(count))
      
      u.rate.plot <- u.rate.plot + geom_point(data=user.rate.sum, aes(x=count, y=mean.val), alpha=0.4, size=2)
      i.rate.plot <- i.rate.plot + geom_point(data=item.rate.sum, aes(x=count, y=mean.val), alpha=0.4, size=2)
      
  } else {
      
      user.rate.sum <- ddply(user.rate, c("tratings"), summarize, mean.val=mean(1/val), stdev.val=sqrt(var(1/val)), freq=length(tratings))
      user.rate.sum2 <- transform(user.rate, avg_clicks=as.integer(tratings/count))
      user.rate.sum2 <- ddply(user.rate.sum2, c("avg_clicks"), summarize, mean.val=mean(1/val), stdev.val=sqrt(var(1/val)), freq=length(avg_clicks), nitems=mean(count))

      item.rate.sum <- ddply(item.rate, c("tratings"), summarize, mean.val=mean(1/val), stdev.val=sqrt(var(1/val)), freq=length(tratings), nusers=mean(count))
      item.rate.sum2 <- transform(item.rate, avg_clicks=as.integer(tratings/count))
      
      #item.rate.sum2 <- subset(item.rate.sum2, count >= 5)
      item.rate.sum2 <- ddply(item.rate.sum2, c("avg_clicks"), summarize, max.val=max(1/val), min.val=min(1/val), freq=length(avg_clicks), nusers=mean(count))
      
      u.rate.plot <- u.rate.plot + geom_point(data=user.rate.sum, aes(x=tratings, y=mean.val), alpha=0.4, size=2)
      u.rate.plot2 <- u.rate.plot2 + geom_linerange(data=user.rate.sum2, aes(x=avg_clicks, ymin=(mean.val+stdev.val), ymax=(mean.val-stdev.val),
                                                    size=as.factor(as.integer(log(freq)+1)), color=as.factor(as.integer(log(nitems)+1))), alpha=0.4)
      i.rate.plot <- i.rate.plot + geom_point(data=item.rate.sum, aes(x=tratings, y=mean.val), alpha=0.4, size=2)
      i.rate.plot2 <- i.rate.plot2 + geom_linerange(data=item.rate.sum2, aes(x=avg_clicks, ymin=min.val, ymax=max.val,
                                                    size=as.factor(as.integer(log(freq)+1)), color=as.factor(as.integer(log(nusers)+1))), alpha=0.4)
      
  }
  i.rate.plot <- i.rate.plot + xlab("Total ratings for item") + ylab("Expected item popularity")
  i.rate.plot <- i.rate.plot + opts(axis.title.x = theme_text(size=24), axis.title.y=theme_text(size=24, angle=90))
  i.rate.plot <- i.rate.plot + scale_color_discrete()
  out <- sprintf("../output/figures/exploratory/%s/%s/itemsrate.pdf", topdir, fname)
  ggsave(i.rate.plot, filename=out)

  u.rate.plot <- u.rate.plot + xlab("Total ratings by user") + ylab("Expected user activity")
  u.rate.plot <- u.rate.plot + opts(axis.title.x = theme_text(size=24), axis.title.y=theme_text(size=24, angle=90))
  u.rate.plot <- u.rate.plot + scale_color_discrete()
  out <- sprintf("../output/figures/exploratory/%s/%s/userrate.pdf", topdir, fname)
  ggsave(u.rate.plot, filename=out)

  if (bin == 0) {
      u.rate.plot2 <- u.rate.plot2 + xlab("Mean user clicks") + ylab("Expected user activity")
      u.rate.plot2 <- u.rate.plot2 + opts(axis.title.x = theme_text(size=24), axis.title.y=theme_text(size=24, angle=90))
      u.rate.plot2 <- u.rate.plot2 + scale_color_discrete()
      u.rate.plot2 <- u.rate.plot2 + labs(size="Log(users in bucket)", color="Log(mean items per user)")
      out <- sprintf("../output/figures/exploratory/%s/%s/userrate2.pdf", topdir, fname)
      ggsave(u.rate.plot2, filename=out)      

      i.rate.plot2 <- i.rate.plot2 + xlab("Mean user clicks") + ylab("Expected item popularity")
      i.rate.plot2 <- i.rate.plot2 + opts(axis.title.x = theme_text(size=24), axis.title.y=theme_text(size=24, angle=90))
      i.rate.plot2 <- i.rate.plot2 + scale_color_discrete()
      i.rate.plot2 <- i.rate.plot2 + labs(size="Log(items in bucket)", color="Log(mean users per item)")
      out <- sprintf("../output/figures/exploratory/%s/%s/itemsrate2.pdf", topdir, fname)
      ggsave(i.rate.plot2, filename=out)
  }
}

analyze("nyt", 0, 1)
analyze("nyt", 1, 1)
analyze("netflix", 0, 1)
analyze("netflix", 1, 1)
analyze("echonest", 1, 5)
analyze("echonest", 0, 5)
#analyze("mendeley", 1, 1)
