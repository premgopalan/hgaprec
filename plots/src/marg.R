library(MASS)
library(ggplot2)
library(plyr)

f <- function(topdir) {

byusers.fname <- sprintf ("../data/%s/byusers.tsv", topdir)
byitems.fname <- sprintf ("../data/%s/byitems.tsv", topdir)
byusers <- read.table(byusers.fname)
colnames(byusers) <- c("seq", "id", "count","tratings")
byusers <- subset(byusers, count >= 5)

byitems <- read.table(byitems.fname)
colnames(byitems) <- c("seq", "id", "count","tratings")
byitems <- subset(byitems, count >= 5)

x <- sprintf("%s users (activity)", topdir)
plot.data <- ddply(byusers, c("count"), summarize, num_users=length(count))
p <- ggplot(data=plot.data, aes(x=count, y=num_users))
p <- p + geom_point()
p <- p + scale_x_log10(labels=comma, breaks=10^(0:4)) + scale_y_log10(labels=comma, breaks=10^(0:6))
p <- p + xlab('User activity') + ylab('Number of users')
p <- p + opts(title=x)
out <- sprintf("../output/figures/marginals5/%s-users.pdf", topdir)
ggsave(p, filename=out, width=12, height=4)

x <- sprintf("%s items (popularity)", topdir)
plot.data <- ddply(byitems, c("count"), summarize, num_items=length(count))
p <- ggplot(data=plot.data, aes(x=count, y=num_items))
p <- p + geom_point()
p <- p + scale_x_log10(labels=comma, breaks=10^(0:4)) + scale_y_log10(labels=comma, breaks=10^(0:6))
p <- p + xlab('Item popularity') + ylab('Number of items')
p <- p + opts(title=x)
out <- sprintf("../output/figures/marginals5/%s-items.pdf", topdir)
ggsave(p, filename=out, width=12, height=4)


nbfit <- fitdistr(byusers$count, "negative binomial")
gsfit <- fitdistr(byusers$count, "normal")

e <- as.data.frame(rnegbin(10000, mu = nbfit$estimate[["mu"]], theta=nbfit$estimate[["size"]]))
colnames(e) <- c("e")
g <- as.data.frame(rnorm(10000, mean=gsfit$estimate[["mean"]], sd = gsfit$estimate[["sd"]]))
colnames(g) <- c("g")

x <- sprintf("%s users (activity)", topdir)
tu <- ggplot()
tu <- tu + geom_density(data=g, aes(x=g, color="gauss"), size=1)
tu <- tu + geom_histogram(data=byitems, aes(x=count, y=..density.., fill=..count..), size=1, alpha=0.2) + scale_fill_gradient("Count", low = "green", high = "red")
tu <- tu + geom_density(data=e, aes(x=e, color="nbinom"), size=1)
#tu <- tu + stat_function(fun=dnorm, aes(color="gauss"), args = list(mean = gsfit$estimate[["mean"]], sd = gsfit$estimate[["sd"]]))
tu <- tu + scale_color_discrete() + opts(legend.position="right", legend.title=theme_blank())
tu <- tu + scale_y_continuous("density") + scale_x_log10("User activity",  breaks=10^(0:4))
tu <- tu + opts(title=x)

  
out <- sprintf("../output/figures/marginals5/%s-users-density.pdf", topdir)
ggsave(tu, filename=out)

nbfit <- fitdistr(byitems$count, "negative binomial")
gsfit <- fitdistr(byitems$count, "normal")

e <- as.data.frame(rnegbin(10000, mu = nbfit$estimate[["mu"]], theta=nbfit$estimate[["size"]]))
colnames(e) <- c("e")
g <- as.data.frame(rnorm(10000, mean=gsfit$estimate[["mean"]], sd = gsfit$estimate[["sd"]]))
colnames(g) <- c("g")

x <- sprintf("%s items (popularity)", topdir)
ti <- ggplot()
ti <- ti + geom_density(data=g, aes(x=g, color="gauss"), size=1)
ti <- ti + geom_histogram(data=byitems, aes(x=count, y=..density.., fill=..count..), size=1, alpha=0.2) + scale_fill_gradient("Count", low = "green", high = "red")
ti <- ti + geom_density(data=e, aes(x=e, color="nbinom"), size=1)
#ti <- ti + stat_function(fun=dnorm, aes(color="gauss"), args = list(mean = gsfit$estimate[["mean"]], sd = gsfit$estimate[["sd"]]))
ti <- ti + scale_color_discrete() + opts(legend.position="right", legend.title=theme_blank())
ti <- ti + scale_y_continuous("density") + scale_x_log10("Item popularity", breaks=10^(0:4))
ti <- ti + opts(title=x)

out <- sprintf("../output/figures/marginals/%s-items-density.pdf", topdir)
ggsave(ti, filename=out)

#q1 <- qqplot(e$e, byitems$count, main="QQ plot")
#q2 <- qqplot(g$g, byitems$count, main="QQ plot")
#save(q1, file="q1.pdf")
#save(q2, file="q2.pdf")
}

f("movielens")
f("echonest")
f("mendeley")
f("netflix")
f("nyt")

