# SET DIRECTORY WITH THE MODEL FIT 
dir <- "/n/fs/learning/run/nyt/003"

# SET THE NUMBER OF FACTORS
K <- 100

# SET DATASET SETTINGS
nusers <- 1615675
ntitles <- 64919
ndocs <- 64919
nvocab <- 2655

  thetaf <- sprintf ("%s/theta.tsv", dir)
  theta <- data.frame(read.table(thetaf))
  colnames(theta) <- c("seq", "docid")

  betaf <- sprintf ("%s/beta.tsv", dir)
  beta <- data.frame(read.table(betaf))
  colnames(beta) <- c("seq", "wordid")

  words <- data.frame(read.table("/n/fs/learning/dat/nyt/content/vocab.dat"))
  words <- data.frame(cbind(c(0:(nvocab-1)), words))
  colnames(words) <- c("wordid", "word")

  beta.with.words <- merge(beta, words, by="wordid")
  colnames(beta.with.words) <- c("wordid", "seq", c(1:K), "word")
  #save(beta.with.words, file="beta.with.words.rda")

topics <- function(factor)
{
  c <- sprintf("%d", factor)
  t <- beta.with.words[order(-beta.with.words[[c]])[1:20],]
  #if (t[1,factor] < 1e-5) {
  #  return(0) 
  #}
  print(na.omit(t[c("word", factor)]), row.names=FALSE)
  return(1)
}

alltopics <- function()
{

  for (i in c(1:K)) {
    buf <- sprintf("TOPIC %d", i)
    print(buf)
    topics(i)
  }
}

a <- sprintf("%s/topics.txt", dir)
sink(a)
alltopics()
sink(NULL)
