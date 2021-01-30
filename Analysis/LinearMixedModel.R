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
    Clay<-as.numeric(as.character(d2$Clay_percent))
    Silt<-as.numeric(as.character(d2$Silt_percent))
    fit <- glm(tax1 ~ depth+pH+Cpercent+Ammonia+Nitrate+Clay+Silt+tax2,family="gaussian")
    tab<-cbind(summary(aov(fit))[[1]],"varExplained"=summary(aov(fit))[[1]]$`Sum Sq`/sum(summary(aov(fit))[[1]]$`Sum Sq`)*100)
    o[i]<-tab$varExplained[8]*(fit$coefficients[8]/abs(fit$coefficients[8]))/100
    
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
  rownames(out)<-colnames(out)
  tm<-as.data.frame(as.matrix(tax_table(d)))
  tmOrder<-rownames(tm) #tm$Class, tm$Order, 
  if(!identical(colnames(d1), as.character(tmOrder))){stop("taxa names not match")} # sanity check
  tmOrder1<-tmOrder[order(tm$Phylum,tm$Class,tm$Order,tm$Family,tm$Genus,tm$Species)]
  out<-out[as.character(tmOrder1), as.character(tmOrder1)]
  out
  
  
}

get.no3<-function(v, d2, d1){
  tax1<-v
  o<-c(rep(NA,length(d1)))
  for(i in 1:length(d1)){
    tax2<-NULL
    fit<-NULL
    tab<-NULL
    tax2<-d1[,i]
    #treatment<-as.factor(d2$Treatment)
    trt<-as.factor(d2$Treatment)
    pH<-as.numeric(as.character(d2$pH))
    Cpercent<-as.numeric(as.character(d2$C_percent))
    Ammonia<-as.numeric(as.character(d2$Nh4_ugPerg))
    Nitrate<-as.numeric(as.character(d2$No3_ugPerg))
    Clay<-as.numeric(as.character(d2$Clay_percent))
    Silt<-as.numeric(as.character(d2$Silt_percent))
    fit <- glm(tax1 ~ trt+pH+Cpercent+Ammonia+Nitrate+Clay+Silt+tax2,family="gaussian")
    tab<-cbind(summary(aov(fit))[[1]],"varExplained"=summary(aov(fit))[[1]]$`Sum Sq`/sum(summary(aov(fit))[[1]]$`Sum Sq`)*100)
    o[i]<-tab$varExplained[8]*(fit$coefficients[8]/abs(fit$coefficients[8]))/100
    
  }
  o
}
sppInt3<-function(d){
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
  out<-do.call(cbind, mclapply(d1, get.no3, d2, d1, mc.preschedule = T, mc.cores=10, mc.cleanup = T))
  # return results
  
  rownames(out)<-colnames(out)
  tm<-as.data.frame(as.matrix(tax_table(d)))
  tmOrder<-rownames(tm) #tm$Class, tm$Order, 
  #print(rownames(d1)[1:5])
  #print(as.character(tmOrder)[1:5])
  if(!identical(colnames(d1), as.character(tmOrder))){stop("taxa names not match")} # sanity check
  tmOrder1<-tmOrder[order(tm$Phylum,tm$Class,tm$Order,tm$Family,tm$Genus,tm$Species)]
  out<-out[as.character(tmOrder1), as.character(tmOrder1)]
  #colnames(out)<-colnames(d1)
  out
} # object returned ordered by taxonomy

# d<-readRDS("/Volumes/Seagate\ Expansion\ Drive/Illumina/Processing/DES/GAD/GADbac2020.rds")
# fd<-readRDS("/Volumes/Seagate\ Expansion\ Drive/Illumina/Processing/DES/GAD/GADfun2020.rds")

# import fungal-bacterial dataset
GAD.bac<-readRDS("/Users/dietrich/Documents/GitHub/Anchoring/Data/GAD/GADbac2020.rds")
GAD.Fun<-readRDS("/Users/dietrich/Documents/GitHub/Anchoring/Data/GAD/GADfun2020.rds")
sample_data(GAD.bac)$SeqDepth<-sample_sums(GAD.bac)
sample_data(GAD.Fun)$SeqDepth<-sample_sums(GAD.Fun)
sample_data(GAD.Fun)$SampleDepth<-sample_data(GAD.Fun)$SeqDepth/sample_data(GAD.Fun)$Fun_QPCR
# annotate fungi by fungal traits database 

# normalize by QPCR:
fGA.Q<-GAD.QScale(GAD.Fun,type="F")
bGA.Q<-GAD.QScale(GAD.bac,type="B")

# filter samples to top 4 levels
fGA.Ql<-subset_samples(fGA.Q, Depth=="0_5"|Depth=="5_10"|Depth=="10Ap"|Depth=="Ap30") # rethink this!
bGA.Ql<-subset_samples(bGA.Q, Depth=="0_5"|Depth=="5_10"|Depth=="10Ap"|Depth=="Ap30")

# alpha diversity
a.ffit<-lm(unlist(estimate_richness(fGA.Ql, measures="Observed"))~sample_sums(fGA.Ql))
a.bfit<-lm(unlist(estimate_richness(bGA.Ql, measures="Observed"))~sample_sums(bGA.Ql))
# sanity check
f.meta<-as.data.frame(as.matrix(sample_data(fGA.Ql)))
f.meta$Depth<-factor(f.meta$Depth, levels=c("0_5", "5_10", "10Ap"))
f.meta$Treatment<-factor(f.meta$Treatment, levels=c("NT", "CT", "Org3"))
identical(rownames(f.meta), rownames(estimate_richness(fGA.Q, measures="Observed")))
resids<-a.ffit$residuals
summary(aov(lm(resids~f.meta$Depth*f.meta$Treatment)))
#boxplot
boxplot(resids~f.meta$Depth+f.meta$Treatment, las=2, xlab=NULL,cex.names=0.1,main="Fungal Alpha Diversity")
abline(h = 0, col = 'black', lty=c(2)) 

b.meta<-as.data.frame(as.matrix(sample_data(bGA.Q)))
b.meta
identical(rownames(b.meta), rownames(estimate_richness(bGA.Q, measures="Observed")))
resids<-a.bfit$residuals
summary(aov(lm(resids~b.meta$Depth*b.meta$Treatment)))
#boxplot:
boxplot(resids~b.meta$Depth+b.meta$Treatment, las=2, xlab=NULL,cex.names=0.1, main="Bacterial Alpha Diversity")
abline(h = 0, col = 'black', lty=c(2)) 

# beta diversity / PERMANOVA
f.meta$pH<-as.numeric(f.meta$pH)
f.meta$C_percent<-as.numeric(f.meta$C_percent)
f.meta$No3_ugPerg<-as.numeric(f.meta$No3_ugPerg)
f.meta$Nh4_ugPerg<-as.numeric(f.meta$Nh4_ugPerg)
f.meta$C_N_ratio<-as.numeric(f.meta$C_N_ratio)
f.otu<-as.data.frame(t(as.matrix(otu_table(fGA.Ql))))
f.dist<-vegdist(f.otu, method="bray")
adonis(f.dist~Depth+Treatment+C_N_ratio+pH+C_percent+No3_ugPerg+Nh4_ugPerg, data=f.meta, permutations = 999, method="bray")
adonis(f.dist~C_N_ratio+pH+C_percent+No3_ugPerg+Nh4_ugPerg, data=f.meta, permutations = 999, method="bray")

b.otu<-as.data.frame(as.matrix(otu_table(bGA.Ql)))
b.beta<-adonis(b.otu~Depth+Treatment+pH+C_percent+No3_ugPerg+Nh4_ugPerg, data=b.meta, permutations = 999, method="bray")

# ignore bacteria for now []
# subset taxa to depth and farming system does not work because it lacks power ...

fGA.Q.5<-subset_samples(fGA.Ql, Depth=="0_5")
fGA.Q.10<-subset_samples(fGA.Ql, Depth=="5_10")
fGA.Q.Ap<-subset_samples(fGA.Ql, Depth=="10Ap")
fGA.Q.30<-subset_samples(fGA.Ql, Depth=="Ap30")
fGA.Qf.5<-filter_taxa(fGA.Q.5, function(x) sum(x > 3) > 5, TRUE)
fGA.Qf.10<-filter_taxa(fGA.Q.10, function(x) sum(x > 3) > 5, TRUE)
fGA.Qf.Ap<-filter_taxa(fGA.Q.Ap, function(x) sum(x > 3) > 5, TRUE)
fGA.Qf.30<-filter_taxa(fGA.Q.30, function(x) sum(x > 3) > 5, TRUE)

#fGA.Q.5<-phyloseq::index_reorder(fGA.Q.5, index_type="both")

t1<-Sys.time()
f.5net<-sppInt3(fGA.Qf.5)
f.10net<-sppInt3(fGA.Qf.10)
f.Anet<-sppInt3(fGA.Qf.Ap)
f.30net<-sppInt3(fGA.Qf.30)
Sys.time()-t1

f.5net[is.na(f.5net)]<-0
f.10net[is.na(f.10net)]<-0
f.Anet[is.na(f.Anet)]<-0
f.30net[is.na(f.30net)]<-0

fD.associationSummary<-matrix(nrow=4, ncol=2)
rownames(fD.associationSummary)<-c("D5", "D10", "DAp", "D30")
colnames(fD.associationSummary)<-c("Positive", "Negative")

fD.associationSummary[1,1]<-sum(f.5net>0.5 & f.5net<1)
fD.associationSummary[2,1]<-sum(f.10net>0.5 & f.10net<1)
fD.associationSummary[3,1]<-sum(f.Anet>0.5 & f.Anet<1)
fD.associationSummary[4,1]<-sum(f.30net>0.5 & f.30net<1)
fD.associationSummary[1,2]<-sum(-0.5>f.5net)
fD.associationSummary[2,2]<-sum(-0.5>f.10net)
fD.associationSummary[3,2]<-sum(-0.5>f.Anet)
fD.associationSummary[4,2]<-sum(-0.5>f.30net)
fD.associationSummary
# test some hclust methods for dendrogram

corrplot::corrplot(f.5net, method="color", 
                   col=colorRampPalette(c("red", "white", "blue"))(200), 
                   order="original",
                   #hclust.method = "average",
                   tl.pos="n",
                   title="top 5 cm")
corrplot::corrplot(f.10net, method="color", 
                   col=colorRampPalette(c("red", "white", "blue"))(200), 
                   order="original",
                   #hclust.method = "average",
                   tl.pos="n",
                   title="5-10 cm")
corrplot::corrplot(f.Anet, method="color", 
                   col=colorRampPalette(c("red", "white", "blue"))(200), 
                   order="original",
                   #hclust.method = "average",
                   tl.pos="n",
                   title="10-Ap cm")
corrplot::corrplot(f.30net, method="color", 
                   col=colorRampPalette(c("red", "white", "blue"))(200), 
                   order="original",
                   #hclust.method = "average",
                   tl.pos="n",
                   title="Ap-30 cm")



fQ.NT<-subset_samples(fGA.Ql, Treatment=="NT")
fQ.CT<-subset_samples(fGA.Ql, Treatment=="CT")
fQ.ORG<-subset_samples(fGA.Ql, Treatment=="Org3")
fQ.NTf<-filter_taxa(fQ.NT, function(x) sum(x > 3) > 5, TRUE)
fQ.CTf<-filter_taxa(fQ.CT, function(x) sum(x > 3) > 5, TRUE)
fQ.ORGf<-filter_taxa(fQ.ORG, function(x) sum(x > 3) > 5, TRUE)
t1<-Sys.time()
f.NTnet<-sppInt2(fQ.NTf)
f.CTnet<-sppInt2(fQ.CTf)
f.Org3net<-sppInt2(fQ.ORGf)
Sys.time()-t1

f.NTnet[is.na(f.NTnet)]<-0
f.CTnet[is.na(f.CTnet)]<-0
f.Org3net[is.na(f.Org3net)]<-0

fFS.associationSummary<-matrix(nrow=3, ncol=2)
rownames(fFS.associationSummary)<-c("NT", "CT", "Org")
colnames(fFS.associationSummary)<-c("Positive", "Negative")

fFS.associationSummary[1,1]<-sum(f.NTnet>0.5 & f.NTnet<1)
fFS.associationSummary[2,1]<-sum(f.CTnet>0.5 & f.CTnet<1)
fFS.associationSummary[3,1]<-sum(f.Org3net>0.5 & f.Org3net<1)
fFS.associationSummary[1,2]<-sum(-0.5>f.NTnet)
fFS.associationSummary[2,2]<-sum(-0.5>f.CTnet)
fFS.associationSummary[3,2]<-sum(-0.5>f.Org3net)
fFS.associationSummary

#f.associationSummary2<-matrix(nrow=3, ncol=2)
#rownames(f.associationSummary2)<-c("NT", "CT", "Org")
#colnames(f.associationSummary2)<-c("Positive", "Negative")

corrplot::corrplot(f.NTnet, method="color", 
                   col=colorRampPalette(c("red", "white", "blue"))(200), 
                   order="original",
                   #hclust.method = "average",
                   tl.pos="n",
                   title="NT")
corrplot::corrplot(f.CTnet, method="color", 
                   col=colorRampPalette(c("red", "white", "blue"))(200), 
                   order="original",
                   #hclust.method = "average",
                   tl.pos="n",
                   title="CT")
corrplot::corrplot(f.Org3net, method="color", 
                   col=colorRampPalette(c("red", "white", "blue"))(200), 
                   order="original",
                   #hclust.method = "average",
                   tl.pos="n",
                   title="Org3")

# omit 30-60 depth because there are not enough samples
# do grouped analysis by depth and by treatment for adequate coverage.
# test space for glm
df1<-as.data.frame(t(as.matrix(otu_table(fGA.Q.O5))))
df2<-as.data.frame(as.matrix(sample_data(fGA.Q.O5)))
p<-df2$pH
c<-df2$C_percent
n1<-df2$No3_ugPerg
n2<-df2$Nh4_ugPerg
t2<-df1[,1]
t1<-df1[,5]

oy<-glm(t1~p+c+n1+n2+t2)
summary(oy)


f.int<-sppInt(fGA.Qf)

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

# by treatment and depth!!
# full community assembly
bQ.NT<-subset_samples(bGA.Ql, Treatment=="NT")
bQ.CT<-subset_samples(bGA.Ql, Treatment=="CT")
bQ.ORG<-subset_samples(bGA.Ql, Treatment=="Org3")


t1<-Sys.time()
b.NTnet<-sppInt2(bQ.NT)
Sys.time()-t1
t1<-Sys.time()
b.CTnet<-sppInt2(bQ.CT)
Sys.time()-t1
t1<-Sys.time()
b.Org3net<-sppInt2(bQ.ORG)
Sys.time()-t1

b.NTnet[is.na(b.NTnet)]<-0
b.CTnet[is.na(b.CTnet)]<-0
b.Org3net[is.na(b.Org3net)]<-0
b.associationSummary<-matrix(nrow=3, ncol=2)
rownames(b.associationSummary)<-c("NT", "CT", "Org")
colnames(b.associationSummary)<-c("Positive", "Negative")

b.associationSummary[1,1]<-sum(b.NTnet>0.5 & b.NTnet<1)
b.associationSummary[2,1]<-sum(b.CTnet>0.5 & b.CTnet<1)
b.associationSummary[3,1]<-sum(b.Org3net>0.5 & b.Org3net<1)
b.associationSummary[1,2]<-sum(-0.5>b.NTnet)
b.associationSummary[2,2]<-sum(-0.5>b.CTnet)
b.associationSummary[3,2]<-sum(-0.5>b.Org3net)
b.associationSummary



f.associationSummary2[1,1]<-sum(f.NTnet>0.5 & f.NTnet<1)
f.associationSummary2[2,1]<-sum(f.CTnet>0.5 & f.CTnet<1)
f.associationSummary2[3,1]<-sum(f.Org3net>0.5 & f.Org3net<1)
f.associationSummary2[1,2]<-sum(-0.5>f.NTnet)
f.associationSummary2[2,2]<-sum(-0.5>f.CTnet)
f.associationSummary2[3,2]<-sum(-0.5>f.Org3net)
f.associationSummary2


corrplot::corrplot(f.NTnet*((f.NTnet>0.5) + (-0.5 > f.NTnet)), method="color", 
                   col=colorRampPalette(c("red", "white", "blue"))(200), 
                   order="hclust",
                   hclust.method = "average",
                   tl.pos="n")
corrplot::corrplot(f.CTnet*((f.CTnet>0.5) + (-0.5 > f.CTnet)), method="color", 
                   col=colorRampPalette(c("red", "white", "blue"))(200), 
                   order="hclust",
                   hclust.method = "average",
                   tl.pos="n")
corrplot::corrplot(f.Org3net*((f.Org3net>0.5) + (-0.5 > f.Org3net)), method="color", 
                   col=colorRampPalette(c("red", "white", "blue"))(200), 
                   order="hclust",
                   hclust.method = "average",
                   tl.pos="n")

# igraph network

c<-f.NTnet[((f.NTnet>0.7) + (-0.7 > f.NTnet))]
sum(c)
n<-graph_from_incidence_matrix(c, directed=T, mode="out")
#cfg<-(n)

plot(n, layout=layout.fruchterman.reingold.grid(n),
     vertex.label=NA, main="NT", vertex.size=1)


corrplot::corrplot(b.int, method="color", 
                   col=colorRampPalette(c("red", "white", "blue"))(200), 
                   order="hclust",
                   hclust.method = "average",
                   tl.pos="n")
hist(sample(b.int, 1000, replace=F))


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

mort<-subset_taxa(fs.ps, Order=="o__Mortierellales")
p3 = plot_bar(mort, x="Treatment", fill="Genus", facet_grid="Depth")
p3 + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")

pip<-subset_taxa(fs.ps, Family=="f__Piptocephalidaceae")
p4 = plot_bar(pip, x="Treatment", fill="Species", facet_grid="Depth")
p4 + geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack")

hyp<-subset_taxa(fs.ps, Order=="o__Hypocreales")
p5 = plot_bar(hyp, x="Treatment", fill="Species", facet_grid="Depth")
p5 + geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack")

can<-subset_taxa(fs.ps, Order=="o__Cantharellales")
#can<-subset_samples(can, sample_sums(can)!=0)
p6 = plot_bar(can, x="Treatment", fill="Family", facet_grid="Depth")
p6 + geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack")

ple<-subset_taxa(fs.ps, Order=="o__Pleosporales")
p7 = plot_bar(ple, x="Treatment", fill="Family", facet_grid="Depth")
p7 + geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack")

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
