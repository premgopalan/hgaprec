library(ggplot2)
library(plyr)

theme_set(theme_bw())

analyze <- function(topdir, bias, hier, bin) {
  if (bin == 1) {
     if (hier == 1 && bias == 1) {
         fname <- "bpf.bias.hier.bin"
     } else if (hier) {
         fname <- "bpf.hier.bin"
    }
  } else {
     if (hier == 1 && bias == 1) {
     	fname <- "bpf.bias.hier"
     } else if (hier) {
     	fname <- "bpf.hier"
     }
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

  if (hier == 1) {
      user.rate.fname <- sprintf("%s/thetarate.tsv", xdir)
      item.rate.fname <- sprintf("%s/betarate.tsv", xdir)
      user.rate <- read.table(user.rate.fname)
      colnames(user.rate) <- c("seq", "id", "val")
      
      item.rate <- read.table(item.rate.fname)
      colnames(item.rate) <- c("seq", "id", "val")

      #item.rate <- item.rate[order(item.rate$val), ]
      #user.rate <- user.rate[order(user.rate$val), ]
    
      user.rate <- merge(byusers, user.rate, by=c("seq", "id"))
      item.rate <- merge(byitems, item.rate, by=c("seq", "id"))

      i.rate.plot <- ggplot()
      u.rate.plot <- ggplot()
      i.rate.plot2 <- ggplot()      
    
      if (bin == 1) {
          #user.rate.sum <- transform(user.rate, bucket=as.integer(log(count, base=1.2)))
          #user.rate.sum <- transform(user.rate, bucket=as.integer(count / 1000))
          #user.rate.sum <- ddply(user.rate.sum, c("bucket"), summarize, mean.val=mean(val), stdev.val=sqrt(var(val)), freq=length(bucket))
          user.rate.sum <- ddply(user.rate, c("count"), summarize, mean.val=mean(val), stdev.val=sqrt(var(val)), freq=length(count))

          #item.rate.sum <- transform(item.rate, bucket=as.integer(count / 1000))
          #item.rate.sum <- ddply(item.rate.sum, c("bucket"), summarize, mean.val=mean(val), stdev.val=sqrt(var(val)), freq=length(bucket))
          item.rate.sum <- ddply(item.rate, c("count"), summarize, mean.val=mean(val), stdev.val=sqrt(var(val)), freq=length(count))

          #u.rate.plot <- u.rate.plot + geom_point(data=user.rate.sum, aes(x=count, y=1/mean.val, color=as.factor(as.integer(freq))), size=2)
          #i.rate.plot <- i.rate.plot + geom_point(data=item.rate.sum, aes(x=count, y=1/mean.val, color=as.factor(as.integer(freq))), size=2)          
          
          u.rate.plot <- u.rate.plot + geom_point(data=user.rate.sum, aes(x=count, y=1/mean.val, alpha=0.1), size=2)
          i.rate.plot <- i.rate.plot + geom_point(data=item.rate.sum, aes(x=count, y=1/mean.val, alpha=0.1), size=2)
          
      } else {
          
          #user.rate.sum <- transform(user.rate, bucket=as.integer(log(tratings, base=1.2)))
          #user.rate.sum <- transform(user.rate, bucket=as.integer(count / 1000))
          #user.rate.sum <- ddply(user.rate.sum, c("bucket"), summarize, mean.val=mean(val), stdev.val=sqrt(var(val)), freq=length(bucket))

          user.rate.sum <- ddply(user.rate, c("tratings"), summarize, mean.val=mean(val), stdev.val=sqrt(var(val)), freq=length(count))          

          #item.rate.sum <- transform(item.rate, bucket=as.integer(log(tratings, base=1.2)))
          #item.rate.sum <- transform(item.rate, bucket=as.integer(count / 1000))
          #item.rate.sum <- ddply(item.rate.sum, c("bucket"), summarize, mean.val=mean(val), stdev.val=sqrt(var(val)), freq=length(bucket))

          item.rate.sum <- ddply(item.rate, c("tratings"), summarize, mean.val=mean(val), stdev.val=sqrt(var(val)), freq=length(count))
          item.rate.sum2 <- transform(item.rate, avg_clicks=as.integer(tratings/count))
          item.rate.sum2 <- ddply(item.rate.sum2, c("avg_clicks"), summarize, mean.val=mean(val), stdev.val=sqrt(var(val)), freq=length(count))
          
          #u.rate.plot <- u.rate.plot + geom_point(data=user.rate.sum, aes(x=1000*bucket, y=1/mean.val, color=as.factor(as.integer(freq))), size=2)
          #i.rate.plot <- i.rate.plot + geom_point(data=item.rate.sum, aes(x=1000*bucket, y=1/mean.val, color=as.factor(as.integer(freq))), size=2)

          u.rate.plot <- u.rate.plot + geom_linerange(data=user.rate.sum, aes(x=tratings, ymin=1/(mean.val + stdev.val), ymax=1/(mean.val-stdev.val), alpha=0.1), size=2)
          #i.rate.plot <- i.rate.plot + geom_point(data=item.rate.sum, aes(x=tratings, y=1/mean.val, alpha=0.1), size=2)
          i.rate.plot <- i.rate.plot + geom_linerange(data=item.rate.sum, aes(x=tratings, ymin=1/(mean.val + stdev.val), ymax=1/(mean.val-stdev.val), alpha=0.1), size=2)
          #i.rate.plot2 <- i.rate.plot2 + geom_point(data=item.rate.sum2, aes(x=avg_clicks, y=1/mean.val, alpha=0.1, color=as.factor(as.integer(freq))), size=2)
          i.rate.plot2 <- i.rate.plot2 + geom_linerange(data=item.rate.sum2, aes(x=avg_clicks, ymin=1/(mean.val + stdev.val), ymax=1/(mean.val-stdev.val),
                                                            alpha=0.1, color=as.factor(as.integer(freq))), size=2)
      }
      i.rate.plot <- i.rate.plot + xlab("Total ratings for item") + ylab("Expected item scale")
      i.rate.plot <- i.rate.plot + opts(axis.title.x = theme_text(size=24), axis.title.y=theme_text(size=24, angle=90))
      i.rate.plot <- i.rate.plot + scale_color_discrete()
#      i.rate.plot <- i.rate.plot + scale_x_log10(breaks=10^(0:4))
      out <- sprintf("../output/figures/exploratory/%s/%s/itemsrate.pdf", topdir, fname)
      ggsave(i.rate.plot, filename=out)

      u.rate.plot <- u.rate.plot + xlab("Total ratings by user") + ylab("Expected user scale")
      u.rate.plot <- u.rate.plot + opts(axis.title.x = theme_text(size=24), axis.title.y=theme_text(size=24, angle=90))
      u.rate.plot <- u.rate.plot + scale_color_discrete()
 #     u.rate.plot <- u.rate.plot + scale_x_log10(breaks=10^(0:4))
      out <- sprintf("../output/figures/exploratory/%s/%s/userrate.pdf", topdir, fname)
      ggsave(u.rate.plot, filename=out)

      if (bin == 0) {
          i.rate.plot2 <- i.rate.plot2 + xlab("Avg. clicks") + ylab("Expected item scale")
          i.rate.plot2 <- i.rate.plot2 + opts(axis.title.x = theme_text(size=24), axis.title.y=theme_text(size=24, angle=90))
          i.rate.plot2 <- i.rate.plot2 + scale_color_discrete()
          out <- sprintf("../output/figures/exploratory/%s/%s/itemsrate2.pdf", topdir, fname)
          ggsave(i.rate.plot2, filename=out)
      }
      
  }
  
  if (bias == 1) {
    user.bias.fname <- sprintf ("%s/thetabias.tsv", xdir)
    item.bias.fname <- sprintf ("%s/betabias.tsv", xdir)
    user.bias <- read.table(user.bias.fname)
    colnames(user.bias) <- c("seq", "id", "val")
    
    item.bias <- read.table(item.bias.fname)
    colnames(item.bias) <- c("seq", "id", "val")
    
    #item.bias <- item.bias[order(item.bias$val), ]
    #user.bias <- user.bias[order(user.bias$val), ]
    
    user.bias <- merge(byusers, user.bias, by=c("seq", "id"))
    item.bias <- merge(byitems, item.bias, by=c("seq", "id"))

    i.bias.plot <- ggplot()
    u.bias.plot <- ggplot()
    if (bin == 1) {
        user.bias.sum <- ddply(user.bias, c("count"), summarize, mean.val=mean(val))
        item.bias.sum <- ddply(item.bias, c("count"), summarize, mean.val=mean(val))
        u.bias.plot <- u.bias.plot + geom_point(data=user.bias.sum, aes(x=count, y=mean.val), color="blue")                                                              
        i.bias.plot <- i.bias.plot + geom_point(data=item.bias.sum, aes(x=count, y=mean.val), color="blue")
    } else {
        user.bias.sum <- ddply(user.bias, c("tratings"), summarize, mean.val=mean(val))
        item.bias.sum <- ddply(item.bias, c("tratings"), summarize, mean.val=mean(val))
        u.bias.plot <- u.bias.plot + geom_point(data=user.bias.sum, aes(x=tratings, y=mean.val), color="blue")                                       
        i.bias.plot <- i.bias.plot + geom_point(data=item.bias.sum, aes(x=tratings, y=mean.val), color="blue")
    }
    i.bias.plot <- i.bias.plot + scale_x_continuous("Total ratings for item") + scale_y_continuous("Expected item bias")
    i.bias.plot <- i.bias.plot + opts(axis.title.x = theme_text(size=24), axis.title.y=theme_text(size=24, angle=90))   
    out <- sprintf("../output/figures/exploratory/%s/%s/itembias.pdf", topdir, fname)
    ggsave(i.bias.plot, filename=out)
    
    u.bias.plot <- u.bias.plot + xlab("Total ratings by user") + ylab("Expected user bias")
    u.bias.plot <- u.bias.plot + opts(axis.title.x = theme_text(size=24), axis.title.y=theme_text(size=24, angle=90))
    out <- sprintf("../output/figures/exploratory/%s/%s/userbias.pdf", topdir, fname)
    ggsave(u.bias.plot, filename=out)
    
    #theta.fname <- sprintf ("%s/theta.tsv", xdir)
    #beta.fname <- sprintf ("%s/beta.tsv", xdir)
    #theta <- read.table(theta.fname)
    #beta <- read.table(beta.fname)
    #colnames(theta) <- c("seq", "id", seq(0:24))
    #colnames(beta) <- c("seq", "id", seq(0:24))
    
  } 
}

#analyze("nyt", 1, 1, 1)
#analyze("nyt", 1, 1, 0)
#analyze("nyt", 0, 1, 1)
#analyze("nyt", 0, 1, 0)

analyze("nyt", 0, 1, 0)
analyze("netflix", 0, 1, 0)
analyze("nyt", 0, 1, 1)
analyze("echonest", 0, 1, 1)
analyze("echonest", 0, 1, 0)
analyze("netflix", 0, 1, 1)


#p <- ggplot(vb) + geom_line(data=vb, aes(x=seq(0:10), y=V1))
#p
