
setwd("~/git/antismash-transat/antismash/specific_modules/nrpspks/nrpspksdomainalign/data/")

l <- read.table('RAxML_bestTree.647KS_RAxML.leaves', header=F, col.names=c('leaf_name'))
c <- read.table('Reference_data2607.csv', sep=',', header=T)
colnames(c)[1] <- 'leaf_name'
library(dplyr)
c <- c %>% filter(leaf_name != '?' & Clade != '?') %>% select(leaf_name, Clade) %>% distinct(leaf_name, Clade)
j <- l %>% left_join(c, by='leaf_name')

n <- j %>% filter(is.na(Clade)) %>% select(leaf_name)
write.csv(n, 'noclade.csv', row.names = F, quote = F)

l2c <- j %>% filter(!(is.na(Clade)))
write.table(j, 'clade_mc.csv', sep=',', row.names = F, quote = F)


j$spec <- sub("^.+_(.+)$", "\\1", j$leaf_name)

c2s <- j %>% select(Clade, spec) %>% distinct(Clade, spec)
cct <- na.omit(c2s %>% group_by(Clade) %>% summarize(count=n()))
write.table(cct, 'clade_counts.csv', sep=',', row.names = F, quote = F)

c2 <- read.table('KSinCompound_CladeInfo_v5b.txt', sep='\t', header=F, skip = 1, col.names=c('leaf_name', 'n', 'fullinfo', 'compound', 'KSindex', 'geneSeqID', 'Clade'))
c2 <- c2 %>% filter(leaf_name != '?' & Clade != '?') %>% select(leaf_name, Clade) %>% distinct(leaf_name, Clade)

c.comb <- rbind(c, c2) %>% distinct(leaf_name, Clade)
