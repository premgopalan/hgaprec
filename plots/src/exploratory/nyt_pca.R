require(ggplot2)
require(reshape)

theme_set(theme_bw())

# load latent factor weights for each article
load('../../data/nyt/dataframes/nyt-beta.rda')

# remove old title, which is partially corrupted
beta.with.titles.b <- beta.with.titles.b[, setdiff(names(beta.with.titles.b), "title")]

# load article metadata from output of parse_nyt_json.py
articles <- read.table('../../data/nyt/articles.tsv',
                       col.names=c('id','section','title','url'),
                       header=F, sep='\t', comment.char='', quote='')

# merge latent factors with metadata
# note: why are we dropping so many rows here?
articles <- merge(beta.with.titles.b, articles, by.x="docid", by.y="id")

# run pca
pca <- princomp(articles[, 2:101])


# add metadata to pca output
plot.data <- data.frame(pca$scores, articles[, c("seq","docid","title","section","url")])
plot.data <- transform(plot.data,
                       section=gsub('\\..*', '', section),
                       Comp.1=Comp.1 - min(Comp.1) + 1,
                       Comp.2=Comp.2 - min(Comp.2) + 1)

# 2-d scatter plot
p <- ggplot(plot.data, aes(x=Comp.1, y=Comp.2, label=title, color=section))
p <- p + geom_text(size=2)
p <- p + xlab('') + ylab('')
p <- p + scale_x_log10() + scale_y_log10()
ggsave(p, file='../../output/figures/exploratory/nyt_pca_2d.pdf', width=10, height=8)
p

# zoomed 2-d scatter plot
p <- p + ylim(c(1,1.005))
ggsave(p, file='../../output/figures/exploratory/nyt_pca_2d_zoom.pdf', width=10, height=8)
p

# 1-d densities
p <- ggplot(plot.data, aes(x=Comp.1, fill=section))
p <- p + geom_density()
p <- p + scale_x_log10()
ggsave(p, file='../../output/figures/exploratory/nyt_pca_1d.pdf', width=10, height=8)
p
