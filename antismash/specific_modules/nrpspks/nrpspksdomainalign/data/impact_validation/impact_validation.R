setwd("~/git/antismash-transat/antismash/specific_modules/nrpspks/nrpspksdomainalign/data/impact_validation/")

v <- read.table("res.tsv", header=T, sep="\t")
library(dplyr)
split.num <- v %>% group_by(Split) %>% summarize(tot=n())
correct_clade <- v %>% filter(CladeMatch=='Y') %>% group_by(Split) %>% summarize(n_correct_clade=n()) %>% left_join(split.num, by='Split') %>% mutate(frac_correct_clade=n_correct_clade/tot)
correct_desc <- v %>% filter(DescMatch=='Y') %>% group_by(Split) %>% summarize(n_correct_desc=n()) %>% left_join(correct_clade, by='Split') %>% mutate(frac_correct_desc=n_correct_desc/tot)

called <- v %>% filter(!(Clade=='clade_not_conserved'))
called.split.num <- called %>% group_by(Split) %>% summarize(tot=n())
called_correct_clade <- called %>% filter(CladeMatch=='Y') %>% group_by(Split) %>% summarize(n_correct_clade=n()) %>% left_join(called.split.num, by='Split') %>% mutate(frac_correct_clade=n_correct_clade/tot)
called_correct_desc <- called %>% filter(DescMatch=='Y') %>% group_by(Split) %>% summarize(n_correct_desc=n()) %>% left_join(called_correct_clade, by='Split') %>% mutate(frac_correct_desc=n_correct_desc/tot)

write.table(correct_desc, "full.summary.tsv", row.names=F, quote=F, sep="\t")
write.table(called_correct_desc, "covered.summary.tsv", row.names=F, quote=F, sep="\t")

truthintrain <- v %>% filter(TruthInTrain > 1)
t.split.num <- truthintrain %>% group_by(Split) %>% summarize(tot=n())
t_correct_clade <- truthintrain %>% filter(CladeMatch=='Y') %>% group_by(Split) %>% summarize(n_correct_clade=n()) %>% left_join(t.split.num, by='Split') %>% mutate(frac_correct_clade=n_correct_clade/tot)
t_correct_desc <- truthintrain %>% filter(DescMatch=='Y') %>% group_by(Split) %>% summarize(n_correct_desc=n()) %>% left_join(t_correct_clade, by='Split') %>% mutate(frac_correct_desc=n_correct_desc/tot)

tcalled <- truthintrain %>% filter(!(Clade=='clade_not_conserved'))
tcalled.split.num <- tcalled %>% group_by(Split) %>% summarize(tot=n())
tcalled_correct_clade <- tcalled %>% filter(CladeMatch=='Y') %>% group_by(Split) %>% summarize(n_correct_clade=n()) %>% left_join(tcalled.split.num, by='Split') %>% mutate(frac_correct_clade=n_correct_clade/tot)
tcalled_correct_desc <- tcalled %>% filter(DescMatch=='Y') %>% group_by(Split) %>% summarize(n_correct_desc=n()) %>% left_join(tcalled_correct_clade, by='Split') %>% mutate(frac_correct_desc=n_correct_desc/tot)

write.table(t_correct_desc, "truthintrain.full.summary.tsv", row.names=F, quote=F, sep="\t")
write.table(tcalled_correct_desc, "truthintrain.covered.summary.tsv", row.names=F, quote=F, sep="\t")


res_by_sz.num.c <- called %>% group_by(Split, TruthInTrain) %>% summarize(tot=n()) 
res_by_sz.c <- called %>% filter(CladeMatch=='Y') %>% group_by(Split, TruthInTrain) %>% summarize(n_correct_clade=n()) %>% left_join(res_by_sz.num.c, by=c('Split', 'TruthInTrain')) %>% mutate(frac_correct_clade=n_correct_clade/tot)
sz.sum.c <- res_by_sz.c %>% group_by(TruthInTrain) %>% summarize(n=n(), m=mean(frac_correct_clade), se=sd(frac_correct_clade)/sqrt(n))
sz.sum.c[is.na(sz.sum.c$se),]$se <- 0
library(ggplot2)
ggplot(sz.sum.c, aes(x=TruthInTrain, y=m)) +
  geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.2)+
  geom_point()+geom_smooth()+
  theme_classic()+
  ylab("Fraction Correct (Covered)")+xlab("Members of Clade in Training Set")+
  scale_x_continuous(breaks=seq(2,18,2))+
  scale_y_continuous(breaks = seq(.8,1,.1), limits=c(0.78,1.06))
ggsave("covered.50-50.pdf", dpi=600, height=5, width=8)

res_by_sz.num <- v %>% group_by(Split, TruthInTrain) %>% summarize(tot=n()) 
res_by_sz <- v %>% filter(CladeMatch=='Y') %>% group_by(Split, TruthInTrain) %>% summarize(n_correct_clade=n()) %>% left_join(res_by_sz.num, by=c('Split', 'TruthInTrain')) %>% mutate(frac_correct_clade=n_correct_clade/tot)
sz.sum <- res_by_sz %>% group_by(TruthInTrain) %>% summarize(n=n(), m=mean(frac_correct_clade), se=sd(frac_correct_clade)/sqrt(n))
sz.sum[is.na(sz.sum$se),]$se <- 0
library(ggplot2)
ggplot(sz.sum, aes(x=TruthInTrain, y=m)) +
  geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.2)+
  geom_point()+geom_smooth()+
  theme_classic()+
  ylab("Fraction Correct (Covered)")+xlab("Members of Clade in Training Set")+
  scale_x_continuous(breaks=seq(2,22,2))+
  scale_y_continuous(breaks = seq(0,1,.2), limits=c(0.1,1.35))
ggsave("full.50-50.pdf", dpi=600, height=5, width=8)










setwd("~/git/antismash-transat/antismash/specific_modules/nrpspks/nrpspksdomainalign/data/impact_validation_90ten/")

v <- read.table("res.tsv", header=T, sep="\t")
library(dplyr)
split.num <- v %>% group_by(Split) %>% summarize(tot=n())
correct_clade <- v %>% filter(CladeMatch=='Y') %>% group_by(Split) %>% summarize(n_correct_clade=n()) %>% left_join(split.num, by='Split') %>% mutate(frac_correct_clade=n_correct_clade/tot)
correct_desc <- v %>% filter(DescMatch=='Y') %>% group_by(Split) %>% summarize(n_correct_desc=n()) %>% left_join(correct_clade, by='Split') %>% mutate(frac_correct_desc=n_correct_desc/tot)

called <- v %>% filter(!(Clade=='clade_not_conserved'))
called.split.num <- called %>% group_by(Split) %>% summarize(tot=n())
called_correct_clade <- called %>% filter(CladeMatch=='Y') %>% group_by(Split) %>% summarize(n_correct_clade=n()) %>% left_join(called.split.num, by='Split') %>% mutate(frac_correct_clade=n_correct_clade/tot)
called_correct_desc <- called %>% filter(DescMatch=='Y') %>% group_by(Split) %>% summarize(n_correct_desc=n()) %>% left_join(called_correct_clade, by='Split') %>% mutate(frac_correct_desc=n_correct_desc/tot)

write.table(correct_desc, "9010.full.summary.tsv", row.names=F, quote=F, sep="\t")
write.table(called_correct_desc, "9010.covered.summary.tsv", row.names=F, quote=F, sep="\t")

truthintrain <- v %>% filter(TruthInTrain > 1)
t.split.num <- truthintrain %>% group_by(Split) %>% summarize(tot=n())
t_correct_clade <- truthintrain %>% filter(CladeMatch=='Y') %>% group_by(Split) %>% summarize(n_correct_clade=n()) %>% left_join(t.split.num, by='Split') %>% mutate(frac_correct_clade=n_correct_clade/tot)
t_correct_desc <- truthintrain %>% filter(DescMatch=='Y') %>% group_by(Split) %>% summarize(n_correct_desc=n()) %>% left_join(t_correct_clade, by='Split') %>% mutate(frac_correct_desc=n_correct_desc/tot)

tcalled <- truthintrain %>% filter(!(Clade=='clade_not_conserved'))
tcalled.split.num <- tcalled %>% group_by(Split) %>% summarize(tot=n())
tcalled_correct_clade <- tcalled %>% filter(CladeMatch=='Y') %>% group_by(Split) %>% summarize(n_correct_clade=n()) %>% left_join(tcalled.split.num, by='Split') %>% mutate(frac_correct_clade=n_correct_clade/tot)
tcalled_correct_desc <- tcalled %>% filter(DescMatch=='Y') %>% group_by(Split) %>% summarize(n_correct_desc=n()) %>% left_join(tcalled_correct_clade, by='Split') %>% mutate(frac_correct_desc=n_correct_desc/tot)

write.table(t_correct_desc, "9010.truthintrain.full.summary.tsv", row.names=F, quote=F, sep="\t")
write.table(tcalled_correct_desc, "9010.truthintrain.covered.summary.tsv", row.names=F, quote=F, sep="\t")


res_by_sz.num.c <- called %>% group_by(Split, TruthInTrain) %>% summarize(tot=n()) 
res_by_sz.c <- called %>% filter(CladeMatch=='Y') %>% group_by(Split, TruthInTrain) %>% summarize(n_correct_clade=n()) %>% left_join(res_by_sz.num.c, by=c('Split', 'TruthInTrain')) %>% mutate(frac_correct_clade=n_correct_clade/tot)
sz.sum.c <- res_by_sz.c %>% group_by(TruthInTrain) %>% summarize(n=n(), m=mean(frac_correct_clade), se=sd(frac_correct_clade)/sqrt(n))
sz.sum.c[is.na(sz.sum.c$se),]$se <- 0
library(ggplot2)
ggplot(sz.sum.c, aes(x=TruthInTrain, y=m)) +
  geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.2)+
  geom_point()+geom_smooth()+
  theme_classic()+
  ylab("Fraction Correct (Covered)")+xlab("Members of Clade in Training Set")+
  scale_x_continuous(breaks=seq(2,28,2))+
  scale_y_continuous(breaks = seq(.8,1,.1), limits=c(0.78,1.06))
ggsave("covered.90-10.pdf", dpi=600, height=5, width=8)

res_by_sz.num <- v %>% group_by(Split, TruthInTrain) %>% summarize(tot=n()) 
res_by_sz <- v %>% filter(CladeMatch=='Y') %>% group_by(Split, TruthInTrain) %>% summarize(n_correct_clade=n()) %>% left_join(res_by_sz.num, by=c('Split', 'TruthInTrain')) %>% mutate(frac_correct_clade=n_correct_clade/tot)
sz.sum <- res_by_sz %>% group_by(TruthInTrain) %>% summarize(n=n(), m=mean(frac_correct_clade), se=sd(frac_correct_clade)/sqrt(n))
sz.sum[is.na(sz.sum$se),]$se <- 0
library(ggplot2)
ggplot(sz.sum, aes(x=TruthInTrain, y=m)) +
  geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.2)+
  geom_point()+geom_smooth()+
  theme_classic()+
  ylab("Fraction Correct (Covered)")+xlab("Members of Clade in Training Set")+
  scale_x_continuous(breaks=seq(2,38,2))+
  scale_y_continuous(breaks = seq(0,1,.2), limits=c(0.1,1.35))
ggsave("full.90-10.pdf", dpi=600, height=5, width=8)









