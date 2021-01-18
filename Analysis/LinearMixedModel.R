# Residual matrix network analysis platform

# linear mixed effects model (categories and covariates)
# doesn't need a random effect, but can be included


library(phyloseq)
library(lme4)
library(ggplot2)
library(dplyr)
library(reshape2)
library(doParallel)
library(parallel)
library(foreach)
library(corrplot)
library(vegan)
library(igraph)

# functions for processing samples:
GAD.QScale<-function(ps, type){
  
  if(type=="B"){
    scale<-as.numeric(as.character(sample_data(ps)$Bac_QPCR))
    out<-GAD.Qscale(ps, val=1, scale)
  }
  if(type=="F"){
    scale<-as.numeric(as.character(sample_data(ps)$Fun_QPCR))
    out<-GAD.Qscale(ps, val=1, scale)
  }
  out
}
GAD.Qscale<-function(ps, val, scale){
  scaled<-data.frame(mapply(`*`, data.frame(t(as.matrix(otu_table(transform_sample_counts(ps, function(x) x/sum(x)))))), scale * val))# sample_data(ps)$val))
  names<-rownames(data.frame(t(as.matrix(otu_table(ps)))))
  rownames(scaled)<-names
  scaled<-round(scaled)
  
  p2<-ps
  otu_table(p2)<- otu_table(scaled, taxa_are_rows=T)
  p2
  
}
taxmodel<-function(d){ 
  # input is phyloseq object
  # already normalized
  d1<-as.data.frame(t(as.matrix(otu_table(d)))) # dim =
  d2<-as.data.frame(as.matrix(sample_data(d)))# dataframes must be samples as rows
  if(!identical(rownames(d1), rownames(d2))){stop("dataframe orientation does not match")}
  treatment<-as.factor(d2$Treatment)
  depth<-as.factor(d2$Depth)
  out<-NULL
  out$residuals<-matrix(ncol=ncol(d1), nrow=length(predict(glm(d1[,1] ~ treatment*depth,family="gaussian"))))
  out$predicted<-matrix(ncol=ncol(d1), nrow=length(predict(glm(d1[,1] ~ treatment*depth,family="gaussian"))))
  out$fit<-list(1:length(ncol(d1)))
  colnames(out$predicted)<-colnames(d1)
  for(i in c(1:ncol(d1))){
    fit<-NULL
    fit <- glm(d1[,i] ~ treatment*depth,family="gaussian")
    
    out$residuals[,i]<-residuals(fit)
    out$predicted[,i]<-predict(fit)
    out$fit[[i]]<-fit
    out$tab[[i]]<-cbind(summary(aov(fit))[[1]],"varExplained"=summary(aov(fit))[[1]]$`Sum Sq`/sum(summary(aov(fit))[[1]]$`Sum Sq`)*100)
  }
  out$predcor<-cor(out$predicted, use="pairwise.complete.obs")
  out$spcor<-cor(d1, use="pairwise.complete.obs")
  #out$hc<-hclust(out$cor, method="ward.D")
  # regress for combined value of original data
  # get residuals of the regression
  # correlate the residuals
  rownames(out$spcor)<-colnames(d1)
  colnames(out$spcor)<-colnames(d1)
  rownames(out$predcor)<-colnames(d1)
  colnames(out$predcor)<-colnames(d1)
  names(out$tab)<-colnames(d1)
  names(out$fit)<-colnames(d1)
  out
}

taxEnvmodel<-function(d){ 
  # input is phyloseq object
  # already normalized
  d1<-as.data.frame(t(as.matrix(otu_table(d)))) # dim =
  d2<-as.data.frame(as.matrix(sample_data(d)))# dataframes must be samples as rows
  if(!identical(rownames(d1), rownames(d2))){stop("dataframe orientation does not match")}
  pH<-as.numeric(as.character(d2$pH))
  Cpercent<-as.numeric(as.character(d2$C_percent))
  Ammonia<-as.numeric(as.character(d2$Nh4_ugPerg))
  Nitrate<-as.numeric(as.character(d2$No3_ugPerg))
  out<-NULL
  out$residuals<-matrix(ncol=ncol(d1), nrow=length(predict(glm(d1[,1] ~ pH+Cpercent+Ammonia+Nitrate,family="gaussian"))))
  out$predicted<-matrix(ncol=ncol(d1), nrow=length(predict(glm(d1[,1] ~ pH+Cpercent+Ammonia+Nitrate,family="gaussian"))))
  out$fit<-list(1:length(ncol(d1)))
  colnames(out$predicted)<-colnames(d1)
  for(i in c(1:ncol(d1))){
    fit<-NULL
    fit <- glm(d1[,i] ~ pH+Cpercent+Ammonia+Nitrate,family="gaussian")
    
    out$residuals[,i]<-residuals(fit)
    out$predicted[,i]<-predict(fit)
    out$fit[[i]]<-fit
    out$tab[[i]]<-cbind(summary(aov(fit))[[1]],"varExplained"=summary(aov(fit))[[1]]$`Sum Sq`/sum(summary(aov(fit))[[1]]$`Sum Sq`)*100)
  }
  out$predcor<-cor(out$predicted, use="pairwise.complete.obs")
  out$spcor<-cor(d1, use="pairwise.complete.obs")
  #out$hc<-hclust(out$cor, method="ward.D")
  # regress for combined value of original data
  # get residuals of the regression
  # correlate the residuals
  rownames(out$spcor)<-colnames(d1)
  colnames(out$spcor)<-colnames(d1)
  rownames(out$predcor)<-colnames(d1)
  colnames(out$predcor)<-colnames(d1)
  names(out$tab)<-colnames(d1)
  names(out$fit)<-colnames(d1)
  out
}
get.no<-function(v, d2){
  tax1<-v
  o<-c(rep(NA,length(d1)))
  for(i in 1:length(d1)){
    tax2<-NULL
    fit<-NULL
    tab<-NULL
    tax2<-d1[,i]
    treatment<-as.factor(d2$Treatment)
    depth<-as.factor(d2$Depth)
    pH<-as.numeric(as.character(d2$pH))
    Cpercent<-as.numeric(as.character(d2$C_percent))
    Ammonia<-as.numeric(as.character(d2$Nh4_ugPerg))
    Nitrate<-as.numeric(as.character(d2$No3_ugPerg))
    fit <- glm(tax1 ~ treatment*depth+pH+Cpercent+Ammonia+Nitrate+tax2,family="gaussian")
    tab<-cbind(summary(aov(fit))[[1]],"varExplained"=summary(aov(fit))[[1]]$`Sum Sq`/sum(summary(aov(fit))[[1]]$`Sum Sq`)*100)
    o[i]<-tab$varExplained[7]*(fit$coefficients[8]/abs(fit$coefficients[8]))/100
   
  }
  o
}
sppInt<-function(d){
  require(foreach)
  require(doParallel)
  require(phyloseq)
  # prepare data
  d1<-as.data.frame(t(as.matrix(otu_table(d)))) # dim =
  d2<-as.data.frame(as.matrix(sample_data(d)))# dataframes must be samples as rows
  if(!identical(rownames(d1), rownames(d2))){stop("dataframe orientation does not match")}
  #treatment<-as.factor(d2$Treatment)
  #depth<-as.factor(d2$Depth)
  #pH<-as.numeric(as.character(d2$pH))
  #Cpercent<-as.numeric(as.character(d2$C_percent))
  #Ammonia<-as.numeric(as.character(d2$Nh4_ugPerg))
  #Nitrate<-as.numeric(as.character(d2$No3_ugPerg))
  
  # prepare environment
  # mc.cores=5,
  out<-do.call(cbind, mclapply(d1, get.no, d2, d1, mc.preschedule = T, mc.cores=10, mc.cleanup = T))
  # return results
  #rownames(out)<-colnames(out)
  #colnames(out)<-colnames(d1)
  out
  
  
}

get.no2<-function(v, d2, d1){
  tax1<-v
  o<-c(rep(NA,length(d1)))
  for(i in 1:length(d1)){
    tax2<-NULL
    fit<-NULL
    tab<-NULL
    tax2<-d1[,i]
    #treatment<-as.factor(d2$Treatment)
    depth<-as.factor(d2$Depth)
    pH<-as.numeric(as.character(d2$pH))
    Cpercent<-as.numeric(as.character(d2$C_percent))
    Ammonia<-as.numeric(as.character(d2$Nh4_ugPerg))
    Nitrate<-as.numeric(as.character(d2$No3_ugPerg))
    fit <- glm(tax1 ~ depth+pH+Cpercent+Ammonia+Nitrate+tax2,family="gaussian")
    tab<-cbind(summary(aov(fit))[[1]],"varExplained"=summary(aov(fit))[[1]]$`Sum Sq`/sum(summary(aov(fit))[[1]]$`Sum Sq`)*100)
    o[i]<-tab$varExplained[6]*(fit$coefficients[8]/abs(fit$coefficients[8]))/100
    
  }
  o
}
sppInt2<-function(d){
  require(foreach)
  require(doParallel)
  require(phyloseq)
  # prepare data
  d1<-as.data.frame(t(as.matrix(otu_table(d)))) # dim =
  d2<-as.data.frame(as.matrix(sample_data(d)))# dataframes must be samples as rows
  if(!identical(rownames(d1), rownames(d2))){stop("dataframe orientation does not match")}
  #treatment<-as.factor(d2$Treatment)
  #depth<-as.factor(d2$Depth)
  #pH<-as.numeric(as.character(d2$pH))
  #Cpercent<-as.numeric(as.character(d2$C_percent))
  #Ammonia<-as.numeric(as.character(d2$Nh4_ugPerg))
  #Nitrate<-as.numeric(as.character(d2$No3_ugPerg))
  
  # prepare environment
  # mc.cores=5,
  out<-do.call(cbind, mclapply(d1, get.no2, d2, d1, mc.preschedule = T, mc.cores=10, mc.cleanup = T))
  # return results
  #rownames(out)<-colnames(out)
  #colnames(out)<-colnames(d1)
  out
  
  
}
# d<-readRDS("/Volumes/Seagate\ Expansion\ Drive/Illumina/Processing/DES/GAD/GADbac2020.rds")
# fd<-readRDS("/Volumes/Seagate\ Expansion\ Drive/Illumina/Processing/DES/GAD/GADfun2020.rds")

# import fungal-bacterial dataset
GAD.bac<-readRDS("/Users/dietrich/Documents/GitHub/Anchoring/Data/GAD/GADbac2020.rds")
GAD.Fun<-readRDS("/Users/dietrich/Documents/GitHub/Anchoring/Data/GAD/GADfun2020.rds")
sample_data(GAD.bac)$SeqDepth<-sample_sums(GAD.bac)
sample_data(GAD.Fun)$SeqDepth<-sample_sums(GAD.Fun)

# normalize by QPCR:
fGA.Q<-GAD.QScale(GAD.Fun,type="F")
bGA.Q<-GAD.QScale(GAD.bac,type="B")

# filter samples to top 4 levels
fGA.Q<-subset_samples(fGA.Q, Depth=="0_5"|Depth=="5_10"|Depth=="10Ap"|Depth=="AP30")
bGA.Q<-subset_samples(bGA.Q, Depth=="0_5"|Depth=="5_10"|Depth=="10Ap"|Depth=="AP30")

# alpha diversity
a.ffit<-lm(unlist(estimate_richness(fGA.Q, measures="Observed"))~sample_sums(fGA.Q))
a.bfit<-lm(unlist(estimate_richness(bGA.Q, measures="Observed"))~sample_sums(bGA.Q))
# sanity check
f.meta<-as.data.frame(as.matrix(sample_data(fGA.Q)))
identical(rownames(f.meta), rownames(estimate_richness(fGA.Q, measures="Observed")))
resids<-a.ffit$residuals
summary(aov(lm(resids~f.meta$Depth*f.meta$Treatment)))
boxplot(resids~f.meta$Depth+f.meta$Treatment, las=2, xlab=NULL,cex.names=0.1,main="Fungal Alpha Diversity")
abline(h = 0, col = 'black', lty=c(2)) 

b.meta<-as.data.frame(as.matrix(sample_data(bGA.Q)))
identical(rownames(b.meta), rownames(estimate_richness(bGA.Q, measures="Observed")))
resids<-a.bfit$residuals
summary(aov(lm(resids~b.meta$Depth*b.meta$Treatment)))
boxplot(resids~b.meta$Depth+b.meta$Treatment, las=2, xlab=NULL,cex.names=0.1, main="Bacterial Alpha Diversity")
abline(h = 0, col = 'black', lty=c(2)) 

# beta diversity / PERMANOVA
f.meta$pH<-as.numeric(f.meta$pH)
f.meta$C_percent<-as.numeric(f.meta$C_percent)
f.meta$No3_ugPerg<-as.numeric(f.meta$No3_ugPerg)
f.meta$Nh4_ugPerg<-as.numeric(f.meta$Nh4_ugPerg)
f.meta$C_N_ratio<-as.numeric(f.meta$C_N_ratio)
f.otu<-as.data.frame(t(as.matrix(otu_table(fGA.Q))))
f.dist<-vegdist(f.otu, method="bray")
adonis(f.dist~Depth+Treatment+C_N_ratio+pH+C_percent+No3_ugPerg+Nh4_ugPerg, data=f.meta, permutations = 999, method="bray")
adonis(f.dist~C_N_ratio+pH+C_percent+No3_ugPerg+Nh4_ugPerg, data=f.meta, permutations = 999, method="bray")

b.otu<-as.data.frame(as.matrix(otu_table(bGA.Q)))
b.beta<-adonis(b.otu~Depth+Treatment+pH+C_percent+No3_ugPerg+Nh4_ugPerg, data=b.meta, permutations = 999, method="bray")
# filter uncommon species:

bGA.Qf<-filter_taxa(bGA.Q, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
fGA.Qf<-filter_taxa(fGA.Q, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

# determined that order doesn't seem to make a difference
# gaussian seems to get highest Rsquared
# variance explained by each factor? 

# species interaction model:

f.int<-sppInt(fGA.Qf)
# not centroid or median or ward.D
corrplot::corrplot(f.int, method="color", 
                   col=colorRampPalette(c("red", "white", "blue"))(200), 
                   order="hclust",
                   hclust.method = "average",
                   tl.pos="n")

hist(sample(f.int, 1000, replace=F))

# full community assembly
t1<-Sys.time()
b.int<-sppInt(bGA.Qf)
Sys.time()-t1

# by treatment!!
# full community assembly
bQ.NT<-subset_samples(bGA.Qf, Treatment=="NT")
bQ.CT<-subset_samples(bGA.Qf, Treatment=="CT")
bQ.ORG<-subset_samples(bGA.Qf, Treatment=="Org3")
t1<-Sys.time()
b.NTnet<-sppInt2(bQ.NT)
b.CTnet<-sppInt2(bQ.CT)
b.Org3net<-sppInt2(bQ.ORG)
Sys.time()-t1

b.associationSummary<-matrix(nrow=3, ncol=2)
rownames(b.associationSummary)<-c("NT", "CT", "Org")
colnames(b.associationSummary)<-c("Positive", "Negative")

b.associationSummary[1,1]<-sum(b.NTnet>0.8 & b.NTnet<1)
b.associationSummary[2,1]<-sum(b.CTnet>0.8 & b.CTnet<1)
b.associationSummary[3,1]<-sum(b.ORGnet>0.8 & b.ORGnet<1)
b.associationSummary[1,2]<-sum(-0.8>b.NTnet)
b.associationSummary[2,2]<-sum(-0.8>b.CTnet)
b.associationSummary[3,2]<-sum(-0.8>b.ORGnet)
b.associationSummary

fQ.NT<-subset_samples(fGA.Qf, Treatment=="NT")
fQ.CT<-subset_samples(fGA.Qf, Treatment=="CT")
fQ.ORG<-subset_samples(fGA.Qf, Treatment=="Org3")
t1<-Sys.time()
f.NTnet<-sppInt2(fQ.NT)
f.CTnet<-sppInt2(fQ.CT)
f.Org3net<-sppInt2(fQ.ORG)
Sys.time()-t1

f.associationSummary<-matrix(nrow=3, ncol=2)
rownames(b.associationSummary)<-c("NT", "CT", "Org")
colnames(b.associationSummary)<-c("Positive", "Negative")

f.associationSummary[1,1]<-sum(f.NTnet>0.8 & f.NTnet<1)
f.associationSummary[2,1]<-sum(f.CTnet>0.8 & f.CTnet<1)
f.associationSummary[3,1]<-sum(f.ORGnet>0.8 & f.ORGnet<1)
f.associationSummary[1,2]<-sum(-0.8>f.NTnet)
f.associationSummary[2,2]<-sum(-0.8>f.CTnet)
f.associationSummary[3,2]<-sum(-0.8>f.ORGnet)
f.associationSummary

corrplot::corrplot(b.int, method="color", 
                   col=colorRampPalette(c("red", "white", "blue"))(200), 
                   order="hclust",
                   hclust.method = "average",
                   tl.pos="n")
hist(sample(b.int, 1000, replace=F))

# igraph network:

# model species by farming system.
t1<-Sys.time()
mtax.fun<-taxmodel(fGA.Qf)
Sys.time()-t1
trt<-rep(NA, length(mtax.fun$tab))
for(i in c(1:length(mtax.fun$tab))){
  trt[i]<-mtax.fun$tab[[i]]$varExplained[1]
}
hist(trt, main="Variance explained by farming system")
names(trt)<-names(mtax.fun$fit)
f.taxtab<-as.data.frame(as.matrix(tax_table(fGA.Qf)))
fs.tab<-subset(f.taxtab, rownames(f.taxtab) %in% names(trt[trt>30]))
f.otutab<-as.data.frame(as.matrix(otu_table(fGA.Qf)))
fs.otutab<-subset(f.otutab, rownames(f.otutab) %in% names(trt[trt>30]))
sdat<-as.data.frame(as.matrix(sample_data(fGA.Qf)))
identical(sdat$Sample, colnames(fs.otutab))
boxplot(unlist(fs.otutab[rownames(fs.otutab)=="GTAAAAGTCGTAACAAGGTCTCCGTTGGTGAACCAGCGGAGGGATCATTACAGGACTCGCAAGACTCCTTAAACCCCTGTGAACTTACTGTTTATACGTTGCTTCGGCGGGTGCTCCGGGGTCCGCCCCGGGGCGCTGCGCCCGCCGGCAGCCTACTTAATTCTGTTTCTCTGCGTTGGCATCTCGAGTAAGCAAAATAAGTTAAAACTTTCAACAACGGATCTCTTGGTTCTGG",])~sdat$Depth+sdat$Treatment,las=2)



# tax glom ...
# boxplot ...
View(fs.otutab)

fs.ps<-phyloseq(otu_table(fs.otutab, taxa_are_rows = T), tax_table(as.matrix(fs.tab)), sample_data(fGA.Qf))
p = plot_bar(fs.ps, x="Treatment", fill="Order", facet_grid="Depth")
p + geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack")

glom<-subset_taxa(fs.ps, Order=="o__Glomerales")
p2 = plot_bar(glom, x="Treatment", fill="Species", facet_grid="Depth")
p2 + geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack")

dep<-rep(NA, length(mtax.fun$tab))
for(i in c(1:length(mtax.fun$tab))){
  dep[i]<-mtax.fun$tab[[i]]$varExplained[2]
}
hist(dep, main="Variance explained by depth")


# bacteria ####
t1<-Sys.time()
mtax.bac<-taxmodel(bGA.Qf)
Sys.time()-t1

trt<-rep(NA, length(mtax.bac$tab))
for(i in c(1:length(mtax.bac$tab))){
 trt[i]<-mtax.bac$tab[[i]]$varExplained[1]
}
hist(trt, main="Percent var explained by Farming System")

dep<-rep(NA, length(mtax.bac$tab))
for(i in c(1:length(mtax.bac$tab))){
  dep[i]<-mtax.bac$tab[[i]]$varExplained[2]
}
hist(dep, main="Percent var explained by Depth")

deptrt<-rep(NA, length(mtax.bac$tab))
for(i in c(1:length(mtax.bac$tab))){
  deptrt[i]<-mtax.bac$tab[[i]]$varExplained[3]
}
hist(deptrt, main="Variance explained by interaction of farming system and depth")

model<-rep(NA, length(mtax.bac$tab))
for(i in c(1:length(mtax.bac$tab))){
  model[i]<-sum(mtax.bac$tab[[i]]$varExplained[1:3])
}
hist(model, main="Variance explained by all factors")


stack<-matrix(data=NA, ncol=length(mtax.bac$tab), nrow=8)
for(i in c(1:length(mtax.bac$tab))){
  stack[,i]<-mtax.bac$tab[[i]]$varExplained
}
stack<-t(stack)
colnames(stack)<-c("System", "Depth", "S-D interaction", "residual")
rownames(stack)<-paste0("spp", 1:nrow(stack))
stack<-data.frame(stack)
stack$spp<-paste0("spp", 1:nrow(stack))
stack[]
stack<-melt(stack, id.vars=c("spp"))#, "System", "Depth", "pH", "Carbon", "Ammonium", "Nitrate", "S.D.interaction"))
#lvl<-names(sort(tapply(stack$variable=="spp", stack$value, mean)))
ggplot(stack, aes(x=spp, y=value, fill=variable))+geom_bar(position="fill", stat="identity")

hist(sample(mtax.bac$predcor, 10000, replace=F), main="Exp. model correlation coefficients")
hist(sample(mtax.bac$spcor, 10000, replace=F), main="Exp. species correlation coefficients")

# environment model ####


t1<-Sys.time()
mtaxEnv.bac<-taxEnvmodel(bGA.Q)
Sys.time()-t1

p<-rep(NA, length(mtaxEnv.bac$tab))
for(i in c(1:length(mtaxEnv.bac$tab))){
  p[i]<-mtaxEnv.bac$tab[[i]]$varExplained[1]
}
hist(p, main="Percent var explained by pH")

Cp<-rep(NA, length(mtaxEnv.bac$tab))
for(i in c(1:length(mtaxEnv.bac$tab))){
  Cp[i]<-mtaxEnv.bac$tab[[i]]$varExplained[2]
}
hist(Cp, main="Variance explained by Carbon")

No<-rep(NA, length(mtaxEnv.bac$tab))
for(i in c(1:length(mtaxEnv.bac$tab))){
  No[i]<-mtaxEnv.bac$tab[[i]]$varExplained[4]
}
hist(No, main="Variance explained by Nitrate")

Nh<-rep(NA, length(mtaxEnv.bac$tab))
for(i in c(1:length(mtaxEnv.bac$tab))){
  Nh[i]<-mtaxEnv.bac$tab[[i]]$varExplained[3]
}
hist(Nh, main="Variance explained by Ammonium")

model<-rep(NA, length(mtaxEnv.bac$tab))
for(i in c(1:length(mtaxEnv.bac$tab))){
  model[i]<-sum(mtaxEnv.bac$tab[[i]]$varExplained[1:4])
}
hist(model, main="Variance explained by whole model")

sp.model1<-rep(NA, length(mtaxEnv.bac$tab))
#sp.model2<-rep(NA, length(mtaxEnv.bac$tab))
for(i in c(1:length(mtaxEnv.bac$tab))){
  sp.model1[i]<-sum(mtaxEnv.bac$tab[[i]]$varExplained[1:4])
  #sp.model1[i]<-sum(mtaxEnv.bac$tab[[i]]$varExplained[1:4])
}
hist(sample(mtaxEnv.bac$predcor[sp.model1>60,sp.model1>60], 10000, replace=F), main="best model env correlation coefficients")
hist(sample(mtaxEnv.bac$spcor[sp.model1>60,sp.model1>60], 10000, replace=F), main="best species env correlation coefficients")
#hist(sp.model2[sp.model1>60], main="Variance explained by whole model")

# subset to high correlation of models

stack<-matrix(data=NA, ncol=length(mtaxEnv.bac$tab), nrow=6)
for(i in c(1:length(mtaxEnv.bac$tab))){
  stack[,i]<-mtaxEnv.bac$tab[[i]]$varExplained
}
stack<-t(stack)
colnames(stack)<-c("System", "Depth", "pH", "Carbon", "Ammonium", "Nitrate", "S-D interaction", "residual")
rownames(stack)<-paste0("spp", 1:nrow(stack))
stack<-data.frame(stack)
stack$spp<-paste0("spp", 1:nrow(stack))
stack[]
stack<-melt(stack, id.vars=c("spp"))#, "System", "Depth", "pH", "Carbon", "Ammonium", "Nitrate", "S.D.interaction"))
#lvl<-names(sort(tapply(stack$variable=="spp", stack$value, mean)))
ggplot(stack, aes(x=spp, y=value, fill=variable))+geom_bar(position="fill", stat="identity")

hist(sample(mtaxEnv.bac$predcor[mtaxEnv.bac$predcor<1], 10000, replace=F), main="model env correlation coefficients")
hist(sample(mtaxEnv.bac$spcor[mtaxEnv.bac$spcor<1], 10000, replace=F), main="species env correlation coefficients")


# heatmap


#scratch:
######################################
d1<-as.data.frame(t(as.matrix(otu_table(fGA.Qf)))) 
d2<-as.data.frame(as.matrix(sample_data(fGA.Qf)))
treatment<-as.factor(d2$Treatment)
depth<-as.factor(d2$Depth)
pH<-as.numeric(as.character(d2$pH))
Cpercent<-as.numeric(as.character(d2$C_percent))
Ammonia<-as.numeric(as.character(d2$Nh4_ugPerg))
Nitrate<-as.numeric(as.character(d2$No3_ugPerg))
 # test mclapply
v<-c(1:ncol(d1))
names(v)<-colnames(d1)

t1<-Sys.time()
df2<-do.call(cbind, mclapply(d1, get.no, mc.preschedule = T, mc.cores=10, mc.cleanup = T))
Sys.time-t1

df2<-df2/100

sum(df2[df2<0.3])
hist(df2)
# test corplot output
plot.new()
dev.off()
corrplot::corrplot(df2, method="color", 
                   col=colorRampPalette(c("red", "white", "blue"))(200), 
                   order="hclust",
                   hclust.method = "ward.D2",
                   tl.pos="n")

###########################################################


sppInt<-function(d){
  require(foreach)
  require(doParallel)
  require(phyloseq)
  # prepare data
  d1<-as.data.frame(t(as.matrix(otu_table(d)))) # dim =
  d2<-as.data.frame(as.matrix(sample_data(d)))# dataframes must be samples as rows
  if(!identical(rownames(d1), rownames(d2))){stop("dataframe orientation does not match")}
  treatment<-as.factor(d2$Treatment)
  depth<-as.factor(d2$Depth)
  pH<-as.numeric(as.character(d2$pH))
  Cpercent<-as.numeric(as.character(d2$C_percent))
  Ammonia<-as.numeric(as.character(d2$Nh4_ugPerg))
  Nitrate<-as.numeric(as.character(d2$No3_ugPerg))
  
  # prepare environment
  cores<-detectCores(logical=F)
  cores<-cores-1
  cl<-makeCluster(cores)
  registerDoParallel(cl,cores=cores)
  chunk.size<-ncol(d1)/(cores)
  #chunk<-c(1,round(chunk.size*1), round(chunk.size*2), round(chunk.size*3), round(chunk.size*4), round(chunk.size*5))
  out<-matrix(ncol=nrow(d1), nrow=ncol(d1))
  print("data set up")
  # run regressions in parallel
  # .combine="cbind"
  foreach(i=1:cores,.combine="cbind", .inorder=T) %dopar% {
    out<-matrix(ncol=chunk.size, nrow=nrow(d1))
    for(x in ((i-1)*chunk.size+1):(i*chunk.size)){
      for(y in 1:ncol(d1)){
        tax<-d1[,y]
        fit<-NULL
        fit <- glm(df[,x] ~ treatment*depth+pH+Cpercent+Ammonia+Nitrate+tax,family="gaussian")
        tab<-cbind(summary(aov(fit))[[1]],"varExplained"=summary(aov(fit))[[1]]$`Sum Sq`/sum(summary(aov(fit))[[1]]$`Sum Sq`)*100)
        out[y,x]<-tab$varExplained[7]
    }}
  }
  # shut down cluster
  stopImplicitCluster()
  stopCluster(cl)
  
  # return results
  rownames(out)<-colnames(d1)
  colnames(out)<-colnames(d1)
  out
  

}
t1<-Sys.time()
f.taxassoc<-sppInt(fGA.Q)
Sys.time()-t1

# simple foreach
m<-matrix(nrow=10, ncol=20, rnorm(200,36, 10))
cores<-detectCores(logical=F)
cores<-cores-1
cl<-makeCluster(cores)
registerDoParallel(cl,cores=cores)
chunk.size<-ncol(d1)/(cores)
out<-foreach(i=1:cores,.combine="cbind", .inorder=T) %dopar% { 
  df<-d1[,c(((i-1)*chunk[i]+1):chunk[i])]
  o<-matrix(ncol=ncol(df), nrow=nrow(df))
  for(x in 1:ncol(df)){
    
    for(y in 1:ncol(d1)){
      tax<-d1[,y]
      fit<-NULL
      fit <- glm(df[,x] ~ treatment*depth+pH+Cpercent+Ammonia+Nitrate+tax,family="gaussian")
      tab<-cbind(summary(aov(fit))[[1]],"varExplained"=summary(aov(fit))[[1]]$`Sum Sq`/sum(summary(aov(fit))[[1]]$`Sum Sq`)*100)
      o[y,x]<-tab$varExplained[7]
    }
    
  }
  o
}


sppInt<-function(d){
  d1<-as.data.frame(t(as.matrix(otu_table(d)))) # dim =
  d2<-as.data.frame(as.matrix(sample_data(d)))# dataframes must be samples as rows
  if(!identical(rownames(d1), rownames(d2))){stop("dataframe orientation does not match")}
  treatment<-as.factor(d2$Treatment)
  depth<-as.factor(d2$Depth)
  pH<-as.numeric(as.character(d2$pH))
  Cpercent<-as.numeric(as.character(d2$C_percent))
  Ammonia<-as.numeric(as.character(d2$Nh4_ugPerg))
  Nitrate<-as.numeric(as.character(d2$No3_ugPerg))
  out<-NULL
  s<-c(1:ncol(d1))
  out$cor<-matrix(ncol=ncol(d1), nrow=ncol(d1))
  for(i in c(1:ncol(d1))){
    
    for(k in s){
      tax<-d1[,k]
      fit<-NULL
      fit <- glm(d1[,i] ~ treatment*depth+pH+Cpercent+Ammonia+Nitrate+tax,family="gaussian")
      tab<-cbind(summary(aov(fit))[[1]],"varExplained"=summary(aov(fit))[[1]]$`Sum Sq`/sum(summary(aov(fit))[[1]]$`Sum Sq`)*100)
      out$cor[i,k]<-tab$varExplained[7]
    }
    
  }
  #rownames(out$cor)<-colnames(d1)
  #colnames(out$cor)<-colnames(d1)
  out
}
t1<-Sys.time()
test<-sppInt(bGA.Q)
Sys.time()-t1


View(testcor$residuals[c(1:5), c(1:5)])
dim(testcor$cor)
dim(d1)
hist(testcor$cor)
heatmap(testcor$cor[testcor$cor>0.9 | -0.9>testcor$cor])

effectmodel<-function(d){ 
  # input is phyloseq object
  # already normalized
  d1<-as.data.frame(as.matrix(otu_table(d))) # dim =
  d2<-as.data.frame(as.matrix(sample_data(d)))# dataframes must be samples as rows
  if(!identical(rownames(d1), rownames(d2))){stop("dataframe orientation does not match")}
  out<-NULL
  out$residuals<-matrix(ncol=ncol(d1), nrow=nrow(d1))
  out$predicted<-matrix(ncol=ncol(d1), nrow=nrow(d1))
  out$fit<-list(1:length(ncol(d1)))
  colnames(out$predicted)<-colnames(d1)
  for(i in c(1:ncol(d1))){
  fit<-NULL
  fit<-glmer.nb(d1[,i]~d2$Treatment + d2$Depth + d2$Depth^2 +(1|)) # make sure formula is consistent with output
  out$residuals[,i]<-residuals(fit)
  out$predicted[,i]<-predict(fit)
  out$fit[[i]]<-fit
  }
  out$cor<-cor(out$predicted, use="pairwise.complete.obs")
  out$hc<-hclust(out$cor, method="ward.D")
  # regress for combined value of original data
  # get residuals of the regression
  # correlate the residuals
  rownames(out$cor)<-colnames(d1)
  colnames(out$cor)<-colnames(d1)
  out
}


fit <- glm(d1[,10] ~ 1+ as.factor(d2$Treatment)*as.factor(d2$Depth) + as.numeric(as.character(d2$pH)) + as.numeric(as.character(d2$C_percent)) + as.numeric(as.character(d2$Nh4_ugPerg))+as.numeric(as.character(d2$No3_ugPerg)),family="gaussian")

fit2 <- lm(d1[,10] ~ 1+ as.factor(d2$Treatment)*as.factor(d2$Depth) + as.numeric(as.character(d2$pH)) + as.numeric(as.character(d2$C_percent)) + as.numeric(as.character(d2$Nh4_ugPerg))+as.numeric(as.character(d2$No3_ugPerg)))
out<-cbind(summary(aov(fit))[[1]],"varExplained"=summary(aov(fit))[[1]]$`Sum Sq`/sum(summary(aov(fit))[[1]]$`Sum Sq`)*100)
length(predict(fit2))

plot(fit)
rsq(fit)
fit$aic
f<-summary(aov(fit))
f[[1]]$`Sum Sq`

predicted <- stats::predict(fit, type="response")   # Save the predicted values
residuals <- residuals(fit) # Save the residual values

pdf<-data.frame(predicted, predicted2, predicted3)

cor(pdf)
# cor to make a matrix of correlation values
# hclust to condense the data / define groups
#

for(i in c(1:ncol(d1))){
  fit<-NULL
  fit<-lm(d1[,i]~d2$Depth) # make sure formula is consistent with output
  out$residuals[,i]<-residuals(fit)
  print('resids2')
  out$predicted[,i]<-predict(fit)
  print('preds2')
  out$fit[[i]]<-fit
  print('fit3')
  }
