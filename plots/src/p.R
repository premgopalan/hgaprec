require(ggplot2)
require(data.table)
require(plyr)
require(scales)

theme_set(theme_bw())

precision.by.user <- data.frame()
recall.by.user <- data.frame()
coverage.by.item <- data.frame()
for (dataset in c("mendeley", "echonest", "nyt", "netflix", "netflix45")) {
#for (dataset in c("nyt")) {
  print(dataset)
  for (method in c("bpf.hier", "bpf", "nmf", "lda")) {
      tsv <- sprintf('../output/%s/%s/precision.txt', dataset, method)
      print(tsv)
      if (file.exists(tsv))
        precision.by.user <- rbind(precision.by.user,
                                   read.table(tsv, header=T))
      tsv <- sprintf('../output/%s/%s/recall.txt', dataset, method)
      print(tsv)
      if (file.exists(tsv))
       recall.by.user <- rbind(recall.by.user,
                               read.table(tsv, header=T))
    }
}

precision.by.user$dataset <- revalue(precision.by.user$dataset, c("mendeley"="mdy", "echonest"="ecn", "nyt"="nyt","netflix"="nfx", "netflix45"="n45"))
recall.by.user$dataset <- revalue(recall.by.user$dataset, c("mendeley"="mdy", "echonest"="ecn", "nyt"="nyt","netflix"="nfx", "netflix45"="n45"))

plot.data <- ddply(precision.by.user, c("dataset","method","K","num.recs"), summarize, mean.precision=mean(precision))
plot.data <- subset(plot.data, num.recs %in%  c(10,100))
p <- ggplot(plot.data, aes(x=as.factor(num.recs), y=mean.precision))
p <- p + geom_point(aes(color=as.factor(method)), size=3)
p <- p + geom_hline(aes(yintercept=mean.precision, colour=as.factor(method)), size=3)
p <- p + xlab("") + ylab('Normalized mean precision')
p <- p + scale_y_continuous(labels=percent)
p <- p + theme(legend.title=element_blank()) #, legend.position="none") #c(0.8,0.75))
p <- p + facet_wrap(dataset ~ num.recs, nrow=1, scale="free")
p <- p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank())
p <- p + theme(strip.text.x = element_text(size = 8, colour = "blue", face="bold"))
ggsave(p, filename='../output/figures/meanprecision2.pdf', width=10, height=2.5)
p

plot.data <- ddply(recall.by.user, c("dataset","method","K","num.recs"), summarize, mean.recall=mean(recall))
plot.data <- subset(plot.data, num.recs %in%  c(10,100))
p <- ggplot(plot.data, aes(x=as.factor(num.recs), y=mean.recall))
p <- p + geom_point(aes(color=as.factor(method)), size=3)
p <- p + geom_hline(aes(yintercept=mean.recall, colour=as.factor(method)), size=3)
p <- p + xlab("") + ylab('Mean recall')
p <- p + scale_y_continuous(labels=percent)
p <- p + theme(legend.title=element_blank())
p <- p + facet_wrap(dataset ~ num.recs, nrow=1, scale="free")
p <- p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank())
p <- p + theme(strip.text.x = element_text(size = 8, colour = "blue", face="bold"))
ggsave(p, filename='../output/figures/meanrecall2.pdf', width=10, height=2.5)
p

N <- 20
q <- subset(precision.by.user, num.recs==N)
plot.data1 <- ddply(q, .(), subset, activity < quantile(activity, 0.1))
plot.data1 <- ddply(plot.data1, c("dataset", "method", "K"), summarize, mean.precision=mean(precision), stdev=sqrt(var(precision)), freq=length(K))
plot.data1 <- transform(plot.data1, qtl=10)
plot.data2 <- ddply(q, .(), subset, activity < quantile(activity, 0.9))
plot.data2 <- ddply(plot.data2, c("dataset", "method", "K"), summarize, mean.precision=mean(precision), stdev=sqrt(var(precision)), freq=length(K))
plot.data2 <- transform(plot.data2, qtl=90)
plot.data <- rbind(plot.data1, plot.data2)
p <- ggplot(plot.data, aes(x=as.factor(qtl), y=mean.precision))
p <- p + geom_point(aes(color=as.factor(method)), size=3)
p <- p + geom_hline(aes(yintercept=mean.precision, colour=as.factor(method)), size=3)
p <- p + xlab("User activity") + ylab('Normalized mean precision')
p <- p + scale_y_continuous(labels=percent)
p <- p + theme(legend.title=element_blank()) #, legend.position="none") #c(0.8,0.75))
p <- p + facet_wrap(dataset ~ qtl, nrow=1, scale="free")
p <- p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank())
p <- p + theme(strip.text.x = element_text(size = 8, colour = "blue", face="bold"))
ggsave(p, filename='../output/figures/useractivity-meanprecision2.pdf', width=10, height=2.5)
p

N <- 20
q <- subset(recall.by.user, num.recs==N)
plot.data1 <- ddply(q, .(), subset, activity < quantile(activity, 0.1))
plot.data1 <- ddply(plot.data1, c("dataset", "method", "K"), summarize, mean.recall=mean(recall), stdev=sqrt(var(recall)), freq=length(K))
plot.data1 <- transform(plot.data1, qtl=10)
plot.data2 <- ddply(q, .(), subset, activity < quantile(activity, 0.9))
plot.data2 <- ddply(plot.data2, c("dataset", "method", "K"), summarize, mean.recall=mean(recall), stdev=sqrt(var(recall)), freq=length(K))
plot.data2 <- transform(plot.data2, qtl=90)
plot.data <- rbind(plot.data1, plot.data2)
p <- ggplot(plot.data, aes(x=as.factor(qtl), y=mean.recall))
p <- p + geom_point(aes(color=as.factor(method)), size=3)
p <- p + geom_hline(aes(yintercept=mean.recall, colour=as.factor(method)), size=3)
p <- p + xlab("User activity") + ylab('Mean recall')
p <- p + scale_y_continuous(labels=percent)
p <- p + theme(legend.title=element_blank()) #, legend.position="none") #c(0.8,0.75))
p <- p + facet_wrap(dataset ~ qtl, nrow=1, scale="free")
p <- p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank())
p <- p + theme(strip.text.x = element_text(size = 8, colour = "blue", face="bold"))
ggsave(p, filename='../output/figures/useractivity-meanrecall2.pdf', width=10, height=2.5)
p

