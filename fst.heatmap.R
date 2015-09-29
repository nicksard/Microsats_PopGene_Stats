# Originally created on August 30th, 2015
# Written by: Nick Sard

#ABOUT: Written to create an allelic richness function

#loading libraries
library(adegenet)
library(hierfstat)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
#loading in functions
source("C:/Users/Nick.Sard/Documents/Research/R/Source scripts/share.R")

#setting working directory
setwd("C:/Users/Nick.Sard/Documents/Research/R/Library development/my.library")

#####################################
### Loading in original genotypes ###
#####################################

#loading in original adult genotypes
adults <- read.table("Input/04DEC13_cougar_adult_master.txt", header=T,stringsAsFactors=F)
head(adults)

#preping HOR/NOR ids
table(adults$type)
adults$type <- ifelse(adults$typ == "H", "HOR", "NOR")
table(adults$type)

head(adults[,c(1,6:27)])

#saving year info later
Pop <- paste0("Adult.",adults[,31])
adult.years <- adults[,31]

#alternative adult ids
Adult.Pop <- rep("Adult",times=nrow(adults))
Adult.Pop1 <- paste("Adult",adults[,2],adults[,31],sep=".")
head(Adult.Pop1)

#getting just the gts
adults <- adults[,c(1,6:27)]
adults <- cbind(Pop,adults)
head(adults)

adults$Pop <- as.character(adults$Pop)
adults$Sample.Name <- as.character(adults$Sample.Name)
all.gts <- adults
my.pops <- unique(all.gts$Pop)[c(1,4)]
my.pops
all.gts <- all.gts[all.gts$Pop  %in% my.pops,]
head(all.gts)

#getting allelic richness calculations from hierfstat
#but first need to get it into the right format
agts2 <- six.char(population = all.gts, missing.data = "0",one.pop = F)
head(agts2)
str(agts2)

#now moving to genind
obj <- df2genind(agts2[,-1:-2], ploidy=2,pop = agts2[,1],ind.names = agts2[,2], loc.names = colnames(agts2)[-c(1,2)])
head(obj)

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

#saving the order for later
my.sites <- colnames(fsts)[-1]
my.sites

#melting the df so its computer reable and fixing the names of the columns
fsts <- melt(fsts, id.vars=colnames(fsts)[1])
names(fsts)=c("popA","popB","value")
head(fsts)
str(fsts)

#getting rid of all the zeros from populations being compared to themselves
fsts <- fsts[(fsts$popA == fsts$popB)==F,]
fsts <- fsts[is.na(fsts$value)==F,]
head(fsts)

#making a column to identify the duplicated values and removing them
fsts$dups <- duplicated(paste(pmin(as.character(fsts$popA), as.character(fsts$popB)),
                             pmax(as.character(fsts$popA), as.character(fsts$popB)), sep="_"))
table(fsts$dups) # should be equal in size
fsts <- fsts[fsts$dups ==F,-ncol(fsts)]
head(fsts)

#getting ride of the significance
head(fsts)

#checking out the fsts values
unique(fsts$value)
summary(fsts$value)
fsts$value <- as.numeric(fsts$value)
summary(fsts$value)

#making some reordering base on factors
fsts[["popA"]]<-factor(fsts[["popA"]],levels= my.sites,ordered=T)
#fsts[["popB"]]<-factor(fsts[["popB"]],levels=rev(my.sites),ordered=T)

#making all values less than 0 just 0
head(fsts)
fsts$value2 <- fsts$value
fsts$value2[fsts$value <0]<-0

#picking pretty colors
my.colors <- brewer.pal(7,"OrRd")

################
### graphics ###
################

#works best if there are several populations to compare to
ggplot(fsts, aes(x = popA, y=popB, fill=value2))+
  geom_tile(color="black")+
  scale_fill_gradientn(colours = my.colors)+ 
  theme_bw()+
  labs(x = "",y="", fill="fsts")+
  theme(legend.title= element_text(size=28),
        axis.text = element_text(size=24),
        axis.text.x = element_text(angle = 45, vjust=1, hjust=1))

#the plot
ggplot(fsts, aes(x = popA, y=popB, fill=value2, label=round(value2,3)))+
  geom_tile(color="black")+
  geom_text()+
  scale_fill_gradientn(colours = my.colors)+ 
  theme_bw()+
  labs(x = "",y="", fill="fsts")+
  theme(legend.title= element_text(size=28),
        axis.text = element_text(size=24),
        axis.text.x = element_text(angle = 45, vjust=1, hjust=1))


#fin!

