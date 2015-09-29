# Originally created on August 30th, 2015
# Written by: Nick Sard

#ABOUT: Written to create an allelic richness function

#loading libraries
library(adegenet)
library(hierfstat)
library(ggplot2)
library(reshape2)

#loading in functions
source("C:/Users/Nick.Sard/Documents/Research/R/Source scripts/share.R")

#setting working directory
setwd("C:/Users/Nick.Sard/Documents/Research/R/Library development/my.library")

#####################################
### Loading in original genotypes ###
#####################################

#loading in original adult genotypes - for this example I make some fake data
all.gts <-sim.pop.gts(n = 1000,loci=10,pops=2, k=20, k.sd =3)
my.pops <- unique(all.gts[,1])
head(all.gts)

#getting allelic richness calculations from hierfstat
#but first need to get it into the right format
agts2 <- six.char(population = all.gts, missing.data = "0",one.pop = F)
head(agts2)
str(agts2)

#now moving to genind
obj <- df2genind(agts2[,-1:-2], ploidy=2,pop = agts2[,1],ind.names = agts2[,2], loc.names = colnames(agts2)[-c(1,2)])
obj

#put to hierfstat
obj1 <- genind2hierfstat(obj)
head(obj1)

#calculating per pair fst
fsts <- pp.fst(obj1)
fsts <- data.frame(fsts$fst.pp)
str(fsts)
colnames(fsts) <- my.pops
row.names(fsts) <- NULL
fsts$pop <-my.pops
fsts <- fsts[,c(ncol(fsts),1:ncol(fsts)-1)]
fsts

################
### graphics ###
################

#percent genotyped
ggplot(df, aes(x=Locus, y=P.GT))+
  facet_wrap(~Pop, ncol=1)+
  geom_bar(stat="identity")+
  labs(x= "Loci",y="Percent individuals genotyped")+
  theme_bw()+
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 28),
        strip.text = element_text(size = 26))

#allelic richness
ggplot(df, aes(x=Locus, y=Ar))+
  facet_wrap(~Pop, ncol=1)+
  geom_bar(stat="identity")+
  labs(x= "Loci",y="Allelic richness")+
  theme_bw()+
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 28),
        strip.text = element_text(size = 26))

#Fis
ggplot(df, aes(x=Locus, y=Fis))+
  facet_wrap(~Pop, ncol=1)+
  geom_point(stat="identity", size=5)+
  labs(x= "Loci",y="Fis")+
  geom_hline(y=0, linetype=2, size=1)+
  theme_bw()+
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 28),
        strip.text = element_text(size = 26))
