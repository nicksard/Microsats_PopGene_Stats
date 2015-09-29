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


#this script relies on the assumption that your genotypic data are in this format:
# population (first column), sample ID (second column), genotypes (all following columns)
# genotypes are assumed to be in three character format (111) and be in two columns (right next to each other)
# finally no other columns should be in the data frame

#here I used simulated allele frequencies

all.gts <- sim.pop.gts(n=100, loci = 10, missing.data = 0, pops = 2, k = 20, k.sd = 3)
head(all.gts)

#saving this for later
my.pops <- unique(all.gts$Pops)

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

#getting some basic stats
bs <- basic.stats(obj1)
head(bs)

#Now extracting the information I want out of bs

#first counting the number of individuals in each population that were successfully genotyped
df1 <- data.frame(bs$n.ind.samp)
df1$locus <- as.character(obj$loc.names)
row.names(df1) <- NULL
df1 <- melt(df1)
df1$variable <- as.character(df1$variable)
for(i in 1:length(my.pops)){
  my.vars <- unique(df1$variable)
  df1$variable[grepl(my.vars[i],df1$variable)==T] <- my.pops[i]
}
names(df1)[3] <- "n.GTed"
df <- df1
head(df)

#now getting observed heterozygosity
df1 <- data.frame(bs$Ho)
df1$locus <- as.character(obj$loc.names)
row.names(df1) <- NULL
df1 <- melt(df1)
df1$variable <- as.character(df1$variable)
for(i in 1:length(my.pops)){
  my.vars <- unique(df1$variable)
  df1$variable[grepl(my.vars[i],df1$variable)==T] <- my.pops[i]
}
df$Ho <- df1$value
head(df)

#now getting expected heterozygosity
df1 <- data.frame(bs$Hs)
df1$locus <- as.character(obj$loc.names)
row.names(df1) <- NULL
df1 <- melt(df1)
df1$variable <- as.character(df1$variable)
for(i in 1:length(my.pops)){
  my.vars <- unique(df1$variable)
  df1$variable[grepl(my.vars[i],df1$variable)==T] <- my.pops[i]
}
df$He <- df1$value
head(df)

#now getting Fis
df1 <- data.frame(bs$Fis)
df1$locus <- as.character(obj$loc.names)
row.names(df1) <- NULL
df1 <- melt(df1)
df1$variable <- as.character(df1$variable)
for(i in 1:length(my.pops)){
  my.vars <- unique(df1$variable)
  df1$variable[grepl(my.vars[i],df1$variable)==T] <- my.pops[i]
}
df$Fis <- df1$value
head(df)

#now getting allelic richness
#getting allelic richness
Ar <- allelic.richness(obj1)
df1 <- data.frame(Ar$Ar)
df1$locus <- as.character(obj$loc.names)
row.names(df1) <- NULL
df1 <- melt(df1)
df1$variable <- as.character(df1$variable)
for(i in 1:length(my.pops)){
  my.vars <- unique(df1$variable)
  df1$variable[grepl(my.vars[i],df1$variable)==T] <- my.pops[i]
}
df$Ar <- df1$value
head(df)

#now getting actual number of alleles observed
A <- nb.alleles(obj1)
df1 <- data.frame(A)
df1$locus <- as.character(obj$loc.names)
row.names(df1) <- NULL
df1 <- melt(df1)
df1$variable <- as.character(df1$variable)
for(i in 1:length(my.pops)){
  my.vars <- unique(df1$variable)
  df1$variable[grepl(my.vars[i],df1$variable)==T] <- my.pops[i]
}
df$A <- df1$value
head(df)

#now using my own function to get the actual count of alleles, the total number of samples in each population,
#the number of missing individuals missing genotypes and calculating the percent individuals genotyped at each locus
y <- gt.summary(population=all.gts,pop.sum=F,one.pop=F)
head(y)

#just getting the columns I want and adding them to the end of df
df <- cbind(df,y[,c(1,2,7:9)])
head(df)

#making sure the loci and populations line up
table(df$locus == df$Locus)
table(df$variable == df$Pop)

#removing the first two columns because they are repeatitive
head(df)
head(df[,c(9,10,11,3,12,13,8,7,4:6)])
df <- df[,c(9,10,11,3,12,13,8,7,4:6)]
head(df)

#using a function I made to identify unique alleles in each population
ua <- uniq.alleles(alfs = alf.freq(population = all.gts,opt = 1))
head(ua)

#now summaryizing this information a format more akin to df
loci <- colnames(agts2)[-c(1,2)]
loci <- data.frame(loci)
pops <- unique(agts2$Pop)
pops

j <- NULL
i <- NULL
for(i in 1:length(pops)){
  ua1 <- ua[ua$Pops == pops[i],]
  for(j in 1:nrow(loci)){
    loci[j,i+1] <-  nrow(ua1[ua1$Locus == loci$loci[j],])
  }
}

loci <- melt(loci)
loci$variable <- as.character(loci$variable)
for(i in 1:length(my.pops)){
  my.vars <- unique(df1$variable)
  df1$variable[grepl(my.vars[i],df1$variable)==T] <- my.pops[i]
}
df$Au <- loci$value
df <- df[,c(1:7,ncol(df),8:11)]
head(df)

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

#calculating averages for all statistics for each population
head(df)


OUT <- NULL

#now for the summary
i <- 1
for(i in 1:length(my.pops)){
  
  #going to fill with information
  df1 <- df[1,]
  df1[1,] <- 0
  
  #getting just the pop of interest
  df2 <- df[df$Pop == my.pops[i],]
  
  #now assigning information
  df1$Locus <- NULL
  
  #adding the pops name
  df1$Pop <- my.pops[i]
  
  #means
  df1$N <- mean(df2$N)
  df1$n.GTed <- mean(df2$n.GTed)
  df1$N.Missing <- mean(df2$N.Missing)
  df1$P.GT <- mean(df2$P.GT)
  df1$A <- mean(df2$A)
  df1$Au <- mean(df2$Au)
  df1$Ar <- mean(df2$Ar)
  df1$Ho <- mean(df2$Ho)
  df1$He <- mean(df2$He)
  df1$Fis <- mean(df2$Fis)

  #standard error
  df1$N.se <- sd(df2$N)/sqrt(df1$N)
  df1$n.GTed.se <- sd(df2$n.GTed)/sqrt(df1$N)
  df1$N.Missing.se <- sd(df2$N.Missing)/sqrt(df1$N)
  df1$P.GT.se <- sd(df2$P.GT)/sqrt(df1$N)
  df1$A.se <- sd(df2$A)/sqrt(df1$N)
  df1$Au.se <- sd(df2$Au)/sqrt(df1$N)
  df1$Ar.se <- sd(df2$Ar)/sqrt(df1$N)
  df1$Ho.se <- sd(df2$Ho)/sqrt(df1$N)
  df1$He.se <- sd(df2$He)/sqrt(df1$N)
  df1$Fis.se <- sd(df2$Fis)/sqrt(df1$N)
  
  #adding to OUT
  OUT <- rbind(OUT,df1)
}

OUT

#writing df and OUT to file
write.table(x = df,file = "per.locus.stats.per.pop.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
write.table(x = OUT,file = "summary.stats.per.pop.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)

#fin!