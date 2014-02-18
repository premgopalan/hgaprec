require(ggplot2)
require(data.table)
require(plyr)
require(scales)

theme_set(theme_bw())

precision.by.user <- data.frame()
recall.by.user <- data.frame()
coverage.by.item <- data.frame()
for (dataset in c("mendeley", "echonest", "nyt", "netflix", "netflix45")) {
  print(dataset)
  for (method in c("bpf.hier", "bpf", "nmf", "lda", "mfpop", "mfunif")) {
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


########################################
# preprocessing
########################################
method.levels <- c("BPF.HIER", "BPF", "LDA", "MF", "NMF")
dataset.levels <- c("mendeley"="Mendeley","nyt"="New York Times","echonest"="Echo Nest","netflix45"="Netflix (implicit)","netflix"="Netflix (explicit)")

# remove users with missing activity from the training file
precision.by.user <- subset(precision.by.user, !is.na(activity))
recall.by.user <- subset(recall.by.user, !is.na(activity))

# keep only mfpop, renaming to mf
# clean up dataset names
precision.by.user <- subset(precision.by.user, method != "MFUNIF")
precision.by.user <- transform(precision.by.user,
                               method=revalue(method, c("MFPOP"="MF")),
                               dataset=revalue(dataset, dataset.levels))
recall.by.user <- subset(recall.by.user, method != "MFUNIF")
recall.by.user <- transform(recall.by.user,
                            method=revalue(method, c("MFPOP"="MF")),
                            dataset=revalue(dataset, dataset.levels))
#recall.by.user <- transform(recall.by.user, method=as.factor(gsub('MFPOP','MF',method)))

# set order of methods and datasets for all plots
precision.by.user <- transform(precision.by.user,
                               dataset=factor(as.character(dataset), dataset.levels),
                               method=factor(as.character(method), method.levels))
recall.by.user <- transform(recall.by.user,
                            dataset=factor(as.character(dataset), dataset.levels),
                            method=factor(as.character(method), method.levels))


# make all evaluations for rank 100 and 20 recommendations
N <- 20
rank <- 100

########################################
# precision/recall at N
########################################

# plot mean precision at N recs for methods and datasets
plot.data <- subset(precision.by.user, num.recs==N & K==rank)
plot.data <- ddply(plot.data, c("dataset","method","K","num.recs"), summarize, mean.precision=mean(precision))
p <- ggplot(plot.data, aes(x=dataset, y=mean.precision))
p <- p + geom_point(aes(color=method), size=1)
p <- p + geom_hline(aes(yintercept=mean.precision, colour=method, linetype=method), size=1)
p <- p + facet_wrap(~ dataset, nrow=1, scale="free")
p <- p + xlab("") + ylab('Normalized mean precision')
p <- p + scale_y_continuous(labels=percent)
p <- p + theme(legend.title=element_blank())
p <- p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank())
ggsave(p, filename=sprintf('../../KDD-paper/figures/mean_precision_at_%d.pdf', N), width=10, height=2.5)
p


# plot mean recall at N recs for methods and datasets
plot.data <- subset(recall.by.user, num.recs==N & K==rank)
plot.data <- ddply(plot.data, c("dataset","method","K","num.recs"), summarize, mean.recall=mean(recall))
p <- ggplot(plot.data, aes(x=dataset, y=mean.recall))
p <- p + geom_point(aes(color=method), size=1)
p <- p + geom_hline(aes(yintercept=mean.recall, colour=method, linetype=method), size=1)
p <- p + facet_wrap(~ dataset, nrow=1, scale="free")
p <- p + xlab("") + ylab('Mean recall')
p <- p + scale_y_continuous(labels=percent)
p <- p + theme(legend.title=element_blank())
p <- p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank())
ggsave(p, filename=sprintf('../../KDD-paper/figures/mean_recall_at_%d.pdf', N), width=10, height=2.5)
p


########################################
# precision/recall by number of recs
########################################

# plot mean precision by number of recs for methods and datasets
plot.data <- subset(precision.by.user, K==rank)
plot.data <- ddply(plot.data, c("dataset","method","K","num.recs"), summarize, mean.precision=mean(precision))
p <- ggplot(plot.data, aes(x=num.recs, y=mean.precision))
p <- p + geom_line(aes(linetype=as.factor(method), colour=as.factor(method)))
p <- p + xlab('Number of recommendations') + ylab('Normalized mean precision')
p <- p + scale_x_continuous(breaks=c(10,50,100)) + scale_y_continuous(labels=percent)
p <- p + theme(legend.title=element_blank())
p <- p + facet_wrap(~ dataset, nrow=1, scale="free_y")
ggsave(p, filename='../../KDD-paper/figures/mean_precision_by_num_recs.pdf', width=10, height=2.5)
p


# plot mean recall by number of recs for methods and datasets
plot.data <- subset(recall.by.user, K==rank)
plot.data <- ddply(plot.data, c("dataset","method","K","num.recs"), summarize, mean.recall=mean(recall))
p <- ggplot(plot.data, aes(x=num.recs, y=mean.recall))
p <- p + geom_line(aes(linetype=as.factor(method), colour=as.factor(method)))
p <- p + xlab('Number of recommendations') + ylab('Mean recall')
p <- p + scale_x_continuous(breaks=c(10,50,100)) + scale_y_continuous(labels=percent)
p <- p + theme(legend.title=element_blank())
p <- p + facet_wrap(~ dataset, nrow=1, scale="free_y")
ggsave(p, filename='../../KDD-paper/figures/mean_recall_by_num_recs.pdf', width=10, height=2.5)
p


########################################
# precision/recall by user activity
########################################

# plot mean precision by user activity percentile
percentiles <- seq(0.05,1,0.05)
plot.data <- subset(precision.by.user, num.recs==N)
plot.data <- ddply(plot.data, c("dataset","method"), function(df) {
  adply(percentiles, 1, function(p) {
    with(subset(df, activity <= quantile(activity, p)), mean(precision, na.rm=T))
  })
})
plot.data$X1 <- percentiles[plot.data$X1]
names(plot.data) <- c("dataset","method","percentile","mean.precision")
p <- ggplot(plot.data, aes(x=percentile, y=mean.precision))
p <- p + geom_line(aes(color=method, linetype=method))
p <- p + facet_wrap(~ dataset, nrow=1, scale="free_y")
p <- p + scale_x_continuous(labels=percent, breaks=c(0.1, 0.5, 0.9))
p <- p + scale_y_continuous(labels=percent)
p <- p + xlab('User percentile by activity') + ylab('Normalized mean precision')
p <- p + theme(legend.title=element_blank())
ggsave(p, filename=sprintf('../../KDD-paper/figures/mean_precision_at_%d_by_user_percentile.pdf', N), width=10, height=2.5)
p


# plot mean recall by user activity percentile
percentiles <- seq(0.05,1,0.05)
plot.data <- subset(recall.by.user, num.recs==N)
plot.data <- ddply(plot.data, c("dataset","method"), function(df) {
  adply(percentiles, 1, function(p) {
    with(subset(df, activity <= quantile(activity, p)), mean(recall, na.rm=T))
  })
})
plot.data$X1 <- percentiles[plot.data$X1]
names(plot.data) <- c("dataset","method","percentile","mean.recall")
p <- ggplot(plot.data, aes(x=percentile, y=mean.recall))
p <- p + geom_line(aes(color=method, linetype=method))
p <- p + facet_wrap(~ dataset, nrow=1, scale="free_y")
p <- p + scale_x_continuous(labels=percent, breaks=c(0.1, 0.5, 0.9))
p <- p + scale_y_continuous(labels=percent)
p <- p + xlab('User percentile by activity') + ylab('Mean recall')
p <- p + theme(legend.title=element_blank())
ggsave(p, filename=sprintf('../../KDD-paper/figures/mean_recall_at_%d_by_user_percentile.pdf', N), width=10, height=2.5)
p



stop("scratch")
########################################
# SCRATCH
########################################

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
ggsave(p, filename='../../KDD-paper/figures/meanprecision2.pdf', width=10, height=2.5)
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
ggsave(p, filename='../../KDD-paper/figures/meanrecall2.pdf', width=10, height=2.5)
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
ggsave(p, filename='../../KDD-paper/figures/useractivity-meanprecision2.pdf', width=10, height=2.5)
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
ggsave(p, filename='../../KDD-paper/figures/useractivity-meanrecall2.pdf', width=10, height=2.5)
p
