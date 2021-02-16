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

# update:
taxmodel<-function(d){ 
  # input is phyloseq object
  # already normalized
  d1<-as.data.frame(t(as.matrix(otu_table(d)))) # dim =
  d2<-as.data.frame(as.matrix(sample_data(d)))# dataframes must be samples as rows
  if(!identical(rownames(d1), rownames(d2))){stop("dataframe orientation does not match")}
  treatment<-as.factor(d2$Treatment)
  effort<-as.numeric(as.character(d2$SampleDepth))
  depth<-as.factor(d2$Depth) # has 4 levels (intercept + 3 levels, then count rest)
  pH<-as.numeric(as.character(d2$pH))
  Cpercent<-as.numeric(as.character(d2$C_percent))
  Ammonia<-as.numeric(as.character(d2$Nh4_ugPerg))
  Nitrate<-as.numeric(as.character(d2$No3_ugPerg))
  Clay<-as.numeric(as.character(d2$Clay_percent))
  Silt<-as.numeric(as.character(d2$Silt_percent))
  out<-NULL
  out$residuals<-matrix(ncol=ncol(d1), nrow=length(predict(glm(d1[,1] ~ effort+treatment*depth+pH+Cpercent+Ammonia+Nitrate+Clay+Silt,family="gaussian"))))
  out$predicted<-matrix(ncol=ncol(d1), nrow=length(predict(glm(d1[,1] ~ effort+treatment*depth+pH+Cpercent+Ammonia+Nitrate+Clay+Silt,family="gaussian"))))
  out$fit<-list(1:length(ncol(d1)))
  colnames(out$predicted)<-colnames(d1)
  for(i in c(1:ncol(d1))){
    fit<-NULL
    fit <- glm(d1[,i] ~ effort+treatment*depth+pH+Cpercent+Ammonia+Nitrate+Clay+Silt,family="gaussian")
    
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
# need to update:
taxEnvmodel<-function(d){ 
  # input is phyloseq object
  # already normalized
  d1<-as.data.frame(t(as.matrix(otu_table(d)))) # dim =
  d2<-as.data.frame(as.matrix(sample_data(d)))# dataframes must be samples as rows
  if(!identical(rownames(d1), rownames(d2))){stop("dataframe orientation does not match")}
  effort<-as.numeric(as.character(d2$SampleDepth))
  pH<-as.numeric(as.character(d2$pH))
  Cpercent<-as.numeric(as.character(d2$C_percent))
  Ammonia<-as.numeric(as.character(d2$Nh4_ugPerg))
  Nitrate<-as.numeric(as.character(d2$No3_ugPerg))
  Clay<-as.numeric(as.character(d2$Clay_percent))
  Silt<-as.numeric(as.character(d2$Silt_percent))
  out<-NULL
  out$residuals<-matrix(ncol=ncol(d1), nrow=length(predict(glm(d1[,1] ~ effort+pH+Cpercent+Ammonia+Nitrate+Clay+Silt,family="gaussian"))))
  out$predicted<-matrix(ncol=ncol(d1), nrow=length(predict(glm(d1[,1] ~ effort+pH+Cpercent+Ammonia+Nitrate+Clay+Silt,family="gaussian"))))
  out$fit<-list(1:length(ncol(d1)))
  colnames(out$predicted)<-colnames(d1)
  for(i in c(1:ncol(d1))){
    fit<-NULL
    fit <- glm(d1[,i] ~ effort+pH+Cpercent+Ammonia+Nitrate+Clay+Silt,family="gaussian")
    
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
    effort<-as.numeric(as.character(d2$SampleDepth))
    depth<-as.factor(d2$Depth) # has 4 levels (intercept + 3 levels, then count rest)
    pH<-as.numeric(as.character(d2$pH))
    Cpercent<-as.numeric(as.character(d2$C_percent))
    Ammonia<-as.numeric(as.character(d2$Nh4_ugPerg))
    Nitrate<-as.numeric(as.character(d2$No3_ugPerg))
    Clay<-as.numeric(as.character(d2$Clay_percent))
    Silt<-as.numeric(as.character(d2$Silt_percent))
    fit <- glm(tax1 ~ effort+depth+pH+Cpercent+Ammonia+Nitrate+Clay+Silt+tax2,family="gaussian")
    tab<-cbind(summary(aov(fit))[[1]],"varExplained"=summary(aov(fit))[[1]]$`Sum Sq`/sum(summary(aov(fit))[[1]]$`Sum Sq`)*100)
    o[i]<-tab$varExplained[9]*(fit$coefficients[12]/abs(fit$coefficients[12]))/100
    
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
    effort<-as.numeric(as.character(d2$SampleDepth))
    trt<-as.factor(d2$Treatment)
    pH<-as.numeric(as.character(d2$pH))
    Cpercent<-as.numeric(as.character(d2$C_percent))
    Ammonia<-as.numeric(as.character(d2$Nh4_ugPerg))
    Nitrate<-as.numeric(as.character(d2$No3_ugPerg))
    Clay<-as.numeric(as.character(d2$Clay_percent))
    Silt<-as.numeric(as.character(d2$Silt_percent))
    fit <- glm(tax1 ~ effort+trt+pH+Cpercent+Ammonia+Nitrate+Clay+Silt+tax2,family="gaussian")
    tab<-cbind(summary(aov(fit))[[1]],"varExplained"=summary(aov(fit))[[1]]$`Sum Sq`/sum(summary(aov(fit))[[1]]$`Sum Sq`)*100)
    o[i]<-tab$varExplained[9]*(fit$coefficients[11]/abs(fit$coefficients[11]))/100 # count the position of tax2 in the effects; the position of var explained = num of levels in factor + number of linear covariates ...
    
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


# m = matrix of effects
# d = phyloseq object of community
# type= on of pp, ps, sp, ss, tp, ts
# pp = pathogen - pathogen interaction
# sp = symbiont - pathogen interaction
# ps = pathogen - symbiont interaction
# ss = symbtiont - symbiont interaction
# tp = total pathogen interactions
# ts = total symbiont interaction
getinteraction<-function(m,d,type){
  out<-c(NA, NA)
  if(type=="pp"){
    tm<-as.data.frame(as.matrix(tax_table(d)))
    tmOrder<-tm$primary_lifestyle #tm$Class, tm$Order, 
    #print(rownames(d1)[1:5])
    #print(as.character(tmOrder)[1:5])
    #if(!identical(colnames(d1), as.character(tmOrder))){stop("taxa names not match")} # sanity check
    tmOrder1<-tmOrder[order(tm$Phylum,tm$Class,tm$Order,tm$Family,tm$Genus,tm$Species)]
    tmOrder1[is.na(tmOrder1)]<-"Unknown"
    m<-m[tmOrder1=="plant_pathogen", tmOrder1=="plant_pathogen"]
    out[1]<-sum(m>0.3& m<1)
    out[2]<-sum(-0.3>m)
  }
  if(type=="sp"){
    tm<-as.data.frame(as.matrix(tax_table(d)))
    tmOrder<-tm$primary_lifestyle #tm$Class, tm$Order, 
    #print(rownames(d1)[1:5])
    #print(as.character(tmOrder)[1:5])
    #if(!identical(colnames(d1), as.character(tmOrder))){stop("taxa names not match")} # sanity check
    tmOrder1<-tmOrder[order(tm$Phylum,tm$Class,tm$Order,tm$Family,tm$Genus,tm$Species)]
    tmOrder1[is.na(tmOrder1)]<-"Unknown"
    m<-m[tmOrder1=="ectomycorrhizal"|tmOrder1=="arbuscular_mycorrhizal"|tmOrder1=="unspecified_symbiotroph"|tmOrder1=="moss_symbiont", tmOrder1=="plant_pathogen"]
    out[1]<-sum(m>0.3& m<1)
    out[2]<-sum(-0.3>m)
  }
  if(type=="ps"){
    tm<-as.data.frame(as.matrix(tax_table(d)))
    tmOrder<-tm$primary_lifestyle #tm$Class, tm$Order, 
    #print(rownames(d1)[1:5])
    #print(as.character(tmOrder)[1:5])
    #if(!identical(colnames(d1), as.character(tmOrder))){stop("taxa names not match")} # sanity check
    tmOrder1<-tmOrder[order(tm$Phylum,tm$Class,tm$Order,tm$Family,tm$Genus,tm$Species)]
    tmOrder1[is.na(tmOrder1)]<-"Unknown"
    m<-m[tmOrder1=="plant_pathogen", tmOrder1=="ectomycorrhizal"|tmOrder1=="arbuscular_mycorrhizal"|tmOrder1=="unspecified_symbiotroph"|tmOrder1=="moss_symbiont"]
    out[1]<-sum(m>0.3& m<1)
    out[2]<-sum(-0.3>m)
  }
  if(type=="ss"){
    tm<-as.data.frame(as.matrix(tax_table(d)))
    tmOrder<-tm$primary_lifestyle #tm$Class, tm$Order, 
    #print(rownames(d1)[1:5])
    #print(as.character(tmOrder)[1:5])
    #if(!identical(colnames(d1), as.character(tmOrder))){stop("taxa names not match")} # sanity check
    tmOrder1<-tmOrder[order(tm$Phylum,tm$Class,tm$Order,tm$Family,tm$Genus,tm$Species)]
    tmOrder1[is.na(tmOrder1)]<-"Unknown"
    m<-m[tmOrder1=="ectomycorrhizal"|tmOrder1=="arbuscular_mycorrhizal"|tmOrder1=="unspecified_symbiotroph"|tmOrder1=="moss_symbiont", tmOrder1=="ectomycorrhizal"|tmOrder1=="arbuscular_mycorrhizal"|tmOrder1=="unspecified_symbiotroph"|tmOrder1=="moss_symbiont"]
    out[1]<-sum(m>0.3& m<1)
    out[2]<-sum(-0.3>m)
  }
  
  if(type=="ts"){
    tm<-as.data.frame(as.matrix(tax_table(d)))
    tmOrder<-tm$primary_lifestyle #tm$Class, tm$Order, 
    #print(rownames(d1)[1:5])
    #print(as.character(tmOrder)[1:5])
    #if(!identical(colnames(d1), as.character(tmOrder))){stop("taxa names not match")} # sanity check
    tmOrder1<-tmOrder[order(tm$Phylum,tm$Class,tm$Order,tm$Family,tm$Genus,tm$Species)]
    tmOrder1[is.na(tmOrder1)]<-"Unknown"
    m<-m[tmOrder1=="ectomycorrhizal"|tmOrder1=="arbuscular_mycorrhizal"|tmOrder1=="unspecified_symbiotroph"|tmOrder1=="moss_symbiont",]
    out[1]<-sum(m>0.3& m<1)
    out[2]<-sum(-0.3>m)
  }
  if(type=="tp"){
    tm<-as.data.frame(as.matrix(tax_table(d)))
    tmOrder<-tm$primary_lifestyle #tm$Class, tm$Order, 
    #print(rownames(d1)[1:5])
    #print(as.character(tmOrder)[1:5])
    #if(!identical(colnames(d1), as.character(tmOrder))){stop("taxa names not match")} # sanity check
    tmOrder1<-tmOrder[order(tm$Phylum,tm$Class,tm$Order,tm$Family,tm$Genus,tm$Species)]
    tmOrder1[is.na(tmOrder1)]<-"Unknown"
    m<-m[tmOrder1=="plant_pathogen",]
    out[1]<-sum(m>0.3& m<1)
    out[2]<-sum(-0.3>m)
  }
  out
}

interact.summary<-function(m, d){
  out<-matrix(nrow=2, ncol=6)
  rownames(out)<-c("positive", "negative")
  colnames(out)<-c("pp", "sp", "ps", "ss", "tp", "ts")
  out[,1]<-getinteraction(m,d,"pp")
  out[,2]<-getinteraction(m,d,"sp")
  out[,3]<-getinteraction(m,d,"ps")
  out[,4]<-getinteraction(m,d,"ss")
  out[,5]<-getinteraction(m,d,"tp")
  out[,6]<-getinteraction(m,d,"ts")
  out
  }

ismatrix<-function(m,d, direction){
  tm<-as.data.frame(as.matrix(tax_table(d)))
  tmOrder<-tm$primary_lifestyle
  tmOrder1<-tmOrder[order(tm$Phylum,tm$Class,tm$Order,tm$Family,tm$Genus,tm$Species)]
  tmOrder1[is.na(tmOrder1)]<-"Unknown"
  out<-matrix(ncol=length(unique(tmOrder1)), nrow=length(unique(tmOrder1)))
  rownames(out)<-unique(tmOrder1)
  colnames(out)<- unique(tmOrder1)
  
  if(direction=="P"){
    for(i in 1:length(unique(tmOrder1))){
      for(j in 1:length(unique(tmOrder1))){
        out[j,i]<-sum(m[tmOrder1==unique(tmOrder1)[j],tmOrder1==unique(tmOrder1)[i]]>0.3 & m[tmOrder1==unique(tmOrder1)[j],tmOrder1==unique(tmOrder1)[i]]<1)
      }
    }
  }
  if(direction=="N"){
    for(i in 1:length(unique(tmOrder1))){
      for(j in 1:length(unique(tmOrder1))){
        out[j,i]<-sum(-0.3>m[tmOrder1==unique(tmOrder1)[j],tmOrder1==unique(tmOrder1)[i]])
      }
    }
  }
  out
}


mematrix<-function(m,d, direction){
  tm<-as.data.frame(as.matrix(tax_table(d)))
  tmOrder<-tm$primary_lifestyle
  tmOrder1<-tmOrder[order(tm$Phylum,tm$Class,tm$Order,tm$Family,tm$Genus,tm$Species)]
  tmOrder1[is.na(tmOrder1)]<-"Unknown"
  out<-matrix(ncol=length(unique(tmOrder1)), nrow=length(unique(tmOrder1)))
  rownames(out)<-unique(tmOrder1)
  colnames(out)<- unique(tmOrder1)
  
  if(direction=="P"){
    for(i in 1:length(unique(tmOrder1))){
      for(j in 1:length(unique(tmOrder1))){
        out[j,i]<-sum(m[tmOrder1==unique(tmOrder1)[j],tmOrder1==unique(tmOrder1)[i]]>0.5 & m[tmOrder1==unique(tmOrder1)[j],tmOrder1==unique(tmOrder1)[i]]<1)
      }
    }
  }
  if(direction=="N"){
    for(i in 1:length(unique(tmOrder1))){
      for(j in 1:length(unique(tmOrder1))){
        out[j,i]<-sum(-0.5>m[tmOrder1==unique(tmOrder1)[j],tmOrder1==unique(tmOrder1)[i]])
      }
    }
  }
  out
}
  

# start analysis ####
# d<-readRDS("/Volumes/Seagate\ Expansion\ Drive/Illumina/Processing/DES/GAD/GADbac2020.rds")
# fd<-readRDS("/Volumes/Seagate\ Expansion\ Drive/Illumina/Processing/DES/GAD/GADfun2020.rds")

# import fungal-bacterial dataset
GAD.bac<-readRDS("/Users/dietrich/Documents/GitHub/Anchoring/Data/GAD/GADbac2020.rds")
GAD.Fun<-readRDS("/Users/dietrich/Documents/GitHub/Anchoring/Data/GAD/GADfun2020.rds")
sample_data(GAD.bac)$SeqDepth<-sample_sums(GAD.bac)
sample_data(GAD.Fun)$SeqDepth<-sample_sums(GAD.Fun)
sample_data(GAD.Fun)$SampleDepth<-sample_data(GAD.Fun)$SeqDepth/sample_data(GAD.Fun)$Fun_QPCR
tt1<-as.data.frame(as.matrix(tax_table(GAD.Fun)))
sum(is.na(tt1$Species))/nrow(tt1) # 75 % unknown to species
sum(is.na(tt1$Genus))/nrow(tt1) # 62 % unknown to genus
sum(is.na(tt1$Family))/nrow(tt1) # 56 % unknown to family
sum(is.na(tt1$Order))/nrow(tt1) # 50 % unknown to order
sum(is.na(tt1$Class))/nrow(tt1) # 46 % unknown to class

# if I keep unknown spp by pasting genus; aggregate to spp, I lose 6% compared to aggregating to family and 12%$ compared to aggregating to order (16% aggregated to class level). Class and above are more likely to be random amplification targets.
GAD.Fun.nf<-GAD.Fun
GAD.Fun.nf<-subset_samples(GAD.Fun.nf, Depth=="0_5"|Depth=="5_10"|Depth=="10Ap"|Depth=="Ap30") # dataset that is not taxa filtered (compare overarching patterns vs completely subset dataset)

GAD.Fun<-tax_glom(GAD.Fun, taxrank="Genus")
tt1<-as.data.frame(as.matrix(tax_table(GAD.Fun)))
tt1$Species[is.na(tt1$Species)]<-paste(tt1$Genus[is.na(tt1$Species)], "Undefined")
tax_table(GAD.Fun)<-tax_table(as.matrix(tt1))


GAD.Fun<-tax_glom(GAD.Fun, taxrank="Species")

# annotate fungi by fungal traits database 
fungalTraits<-read.csv("/Users/dietrich/Documents/GitHub/Plant-Health-Project/Analysis/Data/FungalTraits2021.csv")
View(fungalTraits)
fungalTraits<-as.data.frame(fungalTraits)
fungalTraits$GENUS<-paste0("g__", fungalTraits$GENUS, sep="")
tt<-as.data.frame(as.matrix(tax_table(GAD.Fun)))
tt2<-left_join(tt, fungalTraits, by=c("Genus"="GENUS"))
identical(tt$Family, tt2$Family.x) # sanity check
rownames(tt2)<-rownames(tt)
tt2<-tt2[,-c(8:13,16,21,25,29:31)] # remove columns not being used in the analysis
tax_table(GAD.Fun)<-tax_table(as.matrix(tt2))
# normalize by QPCR:
fGA.Q<-GAD.QScale(GAD.Fun,type="F")
fGA.Qnf<-GAD.QScale(GAD.Fun.nf, type="F")
bGA.Q<-GAD.QScale(GAD.bac,type="B")

# filter samples to top 4 levels
fGA.Ql<-subset_samples(fGA.Q, Depth=="0_5"|Depth=="5_10"|Depth=="10Ap"|Depth=="Ap30") 
bGA.Ql<-subset_samples(bGA.Q, Depth=="0_5"|Depth=="5_10"|Depth=="10Ap"|Depth=="Ap30")

# alpha diversity
a.ffit<-lm(unlist(estimate_richness(fGA.Ql, measures="Observed"))~sample_data(fGA.Ql)$SeqDepth)
#a.bfit<-lm(unlist(estimate_richness(bGA.Ql, measures="Observed"))~sample_sums(bGA.Ql))
summary(aov(a.ffit))  # by seq depth = 17.8% of variance explained at spp level



# sanity check
f.meta<-as.data.frame(as.matrix(sample_data(fGA.Ql)))
f.meta$Depth<-factor(f.meta$Depth, levels=c("0_5", "5_10", "10Ap", "Ap30"))
f.meta$Treatment<-factor(f.meta$Treatment, levels=c("NT", "CT", "Org3"))
f.meta$pH<-as.numeric(as.character(f.meta$pH))
f.meta$C_percent<-as.numeric(as.character(f.meta$C_percent))
f.meta$No3_ugPerg<-as.numeric(as.character(f.meta$No3_ugPerg))
f.meta$Nh4_ugPerg<-as.numeric(as.character(f.meta$Nh4_ugPerg))
f.meta$C_N_ratio<-as.numeric(as.character(f.meta$C_N_ratio))
f.meta$Clay_percent<-as.numeric(as.character(f.meta$Clay_percent))
f.meta$Silt_percent<-as.numeric(as.character(f.meta$Silt_percent))

identical(rownames(f.meta), rownames(estimate_richness(fGA.Ql, measures="Observed")))
resids<-a.ffit$residuals
alphadiveffect<-summary(aov(lm(resids~f.meta$Depth*f.meta$Treatment)))
#boxplot
boxplot(resids~f.meta$Depth+f.meta$Treatment, las=2, xlab=NULL,cex.names=0.1,main="Fungal Alpha Diversity", ylab="Sequencing Depth Residual", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8))
abline(h = 0, col = 'black', lty=c(2)) 

# just sequencing effort (for dissertation):
a.ffit2<-lm(unlist(estimate_richness(fGA.Ql, measures="Observed"))~sample_data(fGA.Ql)$SampleDepth) # by relative sampling effort
summary(aov(a.ffit2)) # by sampleing effort = 24.1% of variance explained
resids2<-a.ffit2$residuals
alphadiveffect2<-summary(aov(lm(resids2~f.meta$Depth*f.meta$Treatment)))
#boxplot
boxplot(resids2~f.meta$Depth+f.meta$Treatment, las=2, xlab=NULL,cex.names=0.1,main="Fungal Alpha Diversity", ylab="Sampling Effort Residual", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8))
abline(h = 0, col = 'black', lty=c(2)) 
# partitioning of residual variance:
# 17.8, 24.1 each comes from variance explained by the model
divVar<-data.frame("SeqDepth"=(100-17.8)*alphadiveffect[[1]]$`Sum Sq`/sum(alphadiveffect[[1]]$`Sum Sq`), "SampleEffort"=(100-24.1)*alphadiveffect2[[1]]$`Sum Sq`/sum(alphadiveffect2[[1]]$`Sum Sq`))
rownames(divVar)<-c("Depth", "Farming System", "Interaction", "Residual")
divVar

# add biodiversity not filtered


# test sources of alph diversity by funcitonal guild
# what does alpha diversity mean in an organic farming system?
# add unspecified            

p.path<-subset_taxa(fGA.Ql, primary_lifestyle=="plant_pathogen")
soil.sap<-subset_taxa(fGA.Ql, primary_lifestyle=="soil_saprotroph")
an.par<-subset_taxa(fGA.Ql, primary_lifestyle=="animal_parasite")
wood.sap<-subset_taxa(fGA.Ql, primary_lifestyle=="wood_saprotroph")
myco.par<-subset_taxa(fGA.Ql, primary_lifestyle=="mycoparasite")
unsp.sap<-subset_taxa(fGA.Ql, primary_lifestyle=="unspecified_saprotroph")
litter.sap<-subset_taxa(fGA.Ql, primary_lifestyle=="litter_saprotroph")
dung.sap<-subset_taxa(fGA.Ql, primary_lifestyle=="dung_saprotroph")
AMF<-subset_taxa(fGA.Ql, primary_lifestyle=="arbuscular_mycorrhizal")
fol.end<-subset_taxa(fGA.Ql, primary_lifestyle=="foliar_endophyte")
nect.sap<-subset_taxa(fGA.Ql, primary_lifestyle=="nectar/tap_saprotroph")
pollen.sap<-subset_taxa(fGA.Ql, primary_lifestyle=="pollen_saprotroph")
lich.par<-subset_taxa(fGA.Ql, primary_lifestyle=="lichen_parasite")
ECM<-subset_taxa(fGA.Ql, primary_lifestyle=="ectomycorrhizal")
unsp.path<-subset_taxa(fGA.Ql, primary_lifestyle=="unspecified_pathotroph")
epi<-subset_taxa(fGA.Ql, primary_lifestyle=="epiphyte")
alg.par<-subset_taxa(fGA.Ql, primary_lifestyle=="algal_parasite")
soot.mold<-subset_taxa(fGA.Ql, primary_lifestyle=="sooty_mold")

# richness of each group
p.path.r<-unlist(estimate_richness(p.path, measures="Observed"))
soil.sap.r<-unlist(estimate_richness(soil.sap, measures="Observed"))
an.par.r<-unlist(estimate_richness(an.par, measures="Observed"))
wood.sap.r<-unlist(estimate_richness(wood.sap, measures="Observed"))
#nknown.r<-unlist(estimate_richness(p.path, measures="Observed")
myco.par.r<-unlist(estimate_richness(myco.par, measures="Observed"))
unsp.sap.r<-unlist(estimate_richness(unsp.sap, measures="Observed"))
litter.sap.r<-unlist(estimate_richness(litter.sap, measures="Observed"))
dung.sap.r<-unlist(estimate_richness(dung.sap, measures="Observed"))
AMF.r<-unlist(estimate_richness(AMF, measures="Observed"))
fol.end.r<-unlist(estimate_richness(fol.end, measures="Observed"))
nect.sap.r<-unlist(estimate_richness(nect.sap, measures="Observed"))
pollen.sap.r<-unlist(estimate_richness(pollen.sap, measures="Observed"))
lich.par.r<-unlist(estimate_richness(lich.par, measures="Observed"))
ECM.r<-unlist(estimate_richness(ECM, measures="Observed"))
unsp.path.r<-unlist(estimate_richness(unsp.path, measures="Observed"))
epi.r<-unlist(estimate_richness(epi, measures="Observed"))
alg.par.r<-unlist(estimate_richness(alg.par, measures="Observed"))
soot.mold.r<-unlist(estimate_richness(soot.mold, measures="Observed"))

div.df<-data.frame(p.path.r,
                   soil.sap.r,
                   an.par.r,
                   wood.sap.r,
                   myco.par.r,
                   unsp.sap.r,
                   litter.sap.r,
                   dung.sap.r,
                   AMF.r,
                   fol.end.r,
                   nect.sap.r,
                   pollen.sap.r,
                   lich.par.r,
                   ECM.r,
                   unsp.path.r,
                   epi.r,
                   alg.par.r,
                   soot.mold.r,
                   "depth"=f.meta$Depth,
                   "treatment"= f.meta$Treatment,
                   "seqs"=as.numeric(as.character(f.meta$SeqDepth)))

with(div.df,boxplot(p.path.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Plant Pathogen Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(p.path.r~seqs+depth*treatment))))

with(div.df,boxplot(soil.sap.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Soil Saprotroph Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(soil.sap.r~seqs+depth*treatment))))

with(div.df,boxplot(an.par.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Animal Parasite Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(an.par.r~seqs+depth*treatment))))

with(div.df,boxplot(wood.sap.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Wood Saprotroph Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(wood.sap.r~seqs+depth*treatment))))


#unique(as.data.frame(as.matrix(tax_table(fGA.Ql)))$primary_lifestyle)
#estimate_richness(fGA.Ql, measures="Observed")

with(div.df,boxplot(myco.par.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Mycoparasite Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(myco.par.r~seqs+depth*treatment))))

with(div.df,boxplot(unsp.sap.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Unspecified Saprotroph Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(unsp.sap.r~seqs+depth*treatment))))

with(div.df,boxplot(litter.sap.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Litter Saprotroph Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(litter.sap.r~seqs+depth*treatment))))

with(div.df,boxplot(dung.sap.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Dung Saprotroph Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(dung.sap.r~seqs+depth*treatment))))

with(div.df,boxplot(AMF.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="AMF Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(AMF.r~seqs+depth*treatment))))

with(div.df,boxplot(fol.end.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Foliar Endophyte Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(fol.end.r~seqs+depth*treatment))))

with(div.df,boxplot(pollen.sap.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Pollen Saprotroph Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(pollen.sap.r~seqs+depth*treatment))))

with(div.df,boxplot(lich.par.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Lichen Parasite Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(lich.par.r~seqs+depth*treatment))))

with(div.df,boxplot(ECM.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="ECM Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(ECM.r~seqs+depth*treatment))))

with(div.df,boxplot(unsp.path.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Unspecified Pathogen Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(unsp.path.r~seqs+depth*treatment))))

with(div.df,boxplot(epi.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Epiphite Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(epi.r~seqs+depth*treatment))))

with(div.df,boxplot(alg.par.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Algal Parasite Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(alg.par.r~seqs+depth*treatment))))

with(div.df,boxplot(soot.mold.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Soot Mold Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(soot.mold.r~seqs+depth*treatment))))

# table of mean (+/- se) for diversity for groups
#library(qwraps2)
library(dplyr)
summary.divdf<-as_tibble(div.df) %>% group_by(depth,treatment) %>% 
  summarise(across(.cols = where(is.numeric), .fns = list(Mean = mean, SD = sd), na.rm = TRUE, .names = "{col}_{fn}"))
summary.divdf2<-as.data.frame(summary.divdf)
summary.divdf3<-data.frame("treatment"=summary.divdf2$treatment,
                           "depth"=summary.divdf2$depth,
                           "plant pathogen"=paste(summary.divdf2$p.path.r_Mean, " (",round(summary.divdf2$p.path.r_SD, digits=2), ")"),
                           "Soil saprotroph"=paste(summary.divdf2$soil.sap.r_Mean, " (",round(summary.divdf2$soil.sap.r_SD,digits=2), ")"),
                           "animal parasite"=paste(summary.divdf2$an.par.r_Mean, " (",round(summary.divdf2$an.par.r_SD,digits=2), ")"),
                           "wood saprotroph"=paste(summary.divdf2$wood.sap.r_Mean, " (",round(summary.divdf2$wood.sap.r_SD,digits=2), ")"),
                           "mycoparasite"=paste(summary.divdf2$myco.par.r_Mean, " (",round(summary.divdf2$myco.par.r_SD,digits=2), ")"),
                           "unsp. saprotroph"=paste(summary.divdf2$unsp.sap.r_Mean, " (",round(summary.divdf2$unsp.sap.r_SD,digits=2), ")"),
                           "litter saprotroph"=paste(summary.divdf2$litter.sap.r_Mean, " (",round(summary.divdf2$litter.sap.r_SD,digits=2), ")"),
                           "dung saprotroph"=paste(summary.divdf2$dung.sap.r_Mean, " (",round(summary.divdf2$dung.sap.r_SD,digits=2), ")"),
                           "arbuscular myc."=paste(summary.divdf2$AMF.r_Mean, " (",round(summary.divdf2$AMF.r_SD,digits=2), ")"),
                           "foliar endophyte"=paste(summary.divdf2$fol.end.r_Mean, " (",round(summary.divdf2$fol.end.r_SD,digits=2), ")"),
                           "nectar saprotroph"=paste(summary.divdf2$nect.sap.r_Mean, " (",round(summary.divdf2$nect.sap.r_SD,digits=2), ")"),
                           "pollen saprotroph"=paste(summary.divdf2$pollen.sap.r_Mean, " (",round(summary.divdf2$pollen.sap.r_SD,digits=2), ")"),
                           "lichen parasite"=paste(summary.divdf2$lich.par.r_Mean, " (",round(summary.divdf2$lich.par.r_SD,digits=2), ")"),
                           "Ectomycorrhiza"=paste(summary.divdf2$ECM.r_Mean, " (",round(summary.divdf2$ECM.r_SD,digits=2), ")"),
                           "unsp. pathogen"=paste(summary.divdf2$unsp.path.r_Mean, " (",round(summary.divdf2$unsp.path.r_SD,digits=2), ")"),
                           "epiphyte"=paste(summary.divdf2$epi.r_Mean, " (",round(summary.divdf2$epi.r_SD,digits=2), ")"),
                           "algal parasite"=paste(summary.divdf2$alg.par.r_Mean, " (",round(summary.divdf2$alg.par.r_SD,digits=2), ")"),
                           "soot mold"=paste(summary.divdf2$soot.mold.r_Mean, " (", round(summary.divdf2$soot.mold.r_SD,digits=2), ")"))

View(summary.divdf3)

# table of mean (+/- se) for abundance of groups ####


#b.meta<-as.data.frame(as.matrix(sample_data(bGA.Q)))
#b.meta
#identical(rownames(b.meta), rownames(estimate_richness(bGA.Q, measures="Observed")))
#resids<-a.bfit$residuals
#summary(aov(lm(resids~b.meta$Depth*b.meta$Treatment)))
#boxplot:
#boxplot(resids~b.meta$Depth+b.meta$Treatment, las=2, xlab=NULL,cex.names=0.1, main="Bacterial Alpha Diversity")
#abline(h = 0, col = 'black', lty=c(2)) 

# beta diversity / PERMANOVA

f.otu<-as.data.frame(t(as.matrix(otu_table(fGA.Ql))))
f.dist<-vegdist(f.otu, method="bray")
adonis(f.dist~Depth+Treatment+C_N_ratio+pH+C_percent+No3_ugPerg+Nh4_ugPerg+Clay_percent+Silt_percent, data=f.meta, permutations = 999, method="bray")
adonis(f.dist~C_N_ratio+pH+C_percent+No3_ugPerg+Nh4_ugPerg+Clay_percent+Silt_percent, data=f.meta, permutations = 999, method="bray")

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



interact.summary(f.5net, fGA.Qf.5)
interact.summary(f.10net, fGA.Qf.10)
interact.summary(f.Anet, fGA.Qf.Ap)
interact.summary(f.30net, fGA.Qf.30)

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

p.guild.int.NT<-ismatrix(f.NTnet, fQ.NTf, "P")
p.guild.int.CT<-ismatrix(f.CTnet, fQ.CTf, "P")
p.guild.int.O3<-ismatrix(f.Org3net, fQ.ORGf, "P")
n.guild.int.NT<-ismatrix(f.NTnet, fQ.NTf, "N")
n.guild.int.CT<-ismatrix(f.CTnet, fQ.CTf, "N")
n.guild.int.O3<-ismatrix(f.Org3net, fQ.ORGf, "N")

View(p.guild.int.NT)
View(p.guild.int.CT)
View(p.guild.int.O3)
View(n.guild.int.NT)
View(n.guild.int.CT)
View(n.guild.int.O3)

interact.summary(f.NTnet, fQ.NTf)
interact.summary(f.CTnet, fQ.CTf)
interact.summary(f.Org3net, fQ.ORGf)
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

# network analysis!!!n####
colrs<-c("#FF0000", "#FF9900", "#CCFF00", "#33FF00", "#00FF66", "#00FFFF", "#0066FF", "#3300FF", "#CC00FF", "#FF0099","#4DAF4A", "#7FC97F", "#7570B3", "#CCCCCC")
group<-as.data.frame(as.matrix(tax_table(fQ.NTf)))$primary_lifestyle
c<-f.NTnet*(f.NTnet>0.3 + (-0.3 > f.NTnet))
sum(c)
n<-graph_from_incidence_matrix(c, directed=T, mode="out")
#cfg<-(n)
V(n)$color<-colrs[as.factor(group)]

n<-delete.vertices(simplify(n), degree(n)==0)
l<-layout_with_dh(n)
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
plot(n,vertex.label=NA, main="NT", vertex.size=5, edge.arrow.size=0.5, 
     #rescale=F, 
     layout=l,
     vertex.color=colrs[as.factor(group)])


c2<-f.CTnet*(f.CTnet>0.3 + (-0.3 > f.CTnet))
sum(c2)
n2<-graph_from_incidence_matrix(c2, directed=T, mode="out")

#cfg<-(n)
n2<-delete.vertices(simplify(n2), degree(n2)==0)
l<-layout_with_dh(n2)
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
plot(n2,vertex.label=NA, main="CT", vertex.size=5, edge.arrow.size=0.5, 
     #rescale=F, 
     layout=l)

closeness(n2)
hist(degree(n2))
tc <- cluster_walktrap(n2)
modularity(n2, unlist(tc))

c3<-f.Org3net*(f.Org3net>0.3 + (-0.3 > f.Org3net))
sum(c3)
n3<-graph_from_incidence_matrix(c3, directed=T, mode="out")
table(V(n3)$type)

plot(n3, layout=layout.fruchterman.reingold.grid(n3),
     vertex.label=NA, main="ORG", vertex.size=5)
# guild questions:

# main effect of farming system ####
t1<-Sys.time()
f.5fs<-sppInt3(fGA.Qf.5)
Sys.time()-t1

f.5fs[is.na(f.5nfs)]<-0
f.10fs[is.na(f.10fs)]<-0
f.Afs[is.na(f.Afs)]<-0
f.30fs[is.na(f.30fs)]<-0

p.guild.int.NT<-ismatrix(f.NTfs, fQ.NTf, "P")
p.guild.int.CT<-ismatrix(f.CTfs, fQ.CTf, "P")
p.guild.int.O3<-ismatrix(f.Org3fs, fQ.ORGf, "P")
n.guild.int.NT<-ismatrix(f.NTfs, fQ.NTf, "N")
n.guild.int.CT<-ismatrix(f.CTnfs, fQ.CTf, "N")
n.guild.int.O3<-ismatrix(f.Org3fs, fQ.ORGf, "N")


fGA.Qlf<-filter_taxa(fGA.Ql, function(x) sum(x > 3) > 5, TRUE)

mp.path<-subset_taxa(fGA.Qlf, primary_lifestyle=="plant_pathogen")
msoil.sap<-subset_taxa(fGA.Qlf, primary_lifestyle=="soil_saprotroph")
man.par<-subset_taxa(fGA.Qlf, primary_lifestyle=="animal_parasite")
mwood.sap<-subset_taxa(fGA.Qlf, primary_lifestyle=="wood_saprotroph")
#munknown<-subset_taxa(fGA.Qlf, primary_lifestyle=="NA")
mmyco.par<-subset_taxa(fGA.Qlf, primary_lifestyle=="mycoparasite")
munsp.sap<-subset_taxa(fGA.Qlf, primary_lifestyle=="unspecified_saprotroph")
mlitter.sap<-subset_taxa(fGA.Qlf, primary_lifestyle=="litter_saprotroph")
mdung.sap<-subset_taxa(fGA.Qlf, primary_lifestyle=="dung_saprotroph")
mAMF<-subset_taxa(fGA.Qlf, primary_lifestyle=="arbuscular_mycorrhizal")
mfol.end<-subset_taxa(fGA.Qlf, primary_lifestyle=="foliar_endophyte")
mnect.sap<-subset_taxa(fGA.Qlf, primary_lifestyle=="nectar/tap_saprotroph")
mpollen.sap<-subset_taxa(fGA.Qlf, primary_lifestyle=="pollen_saprotroph")
mlich.par<-subset_taxa(fGA.Qlf, primary_lifestyle=="lichen_parasite")
mECM<-subset_taxa(fGA.Qlf, primary_lifestyle=="ectomycorrhizal")
munsp.path<-subset_taxa(fGA.Qlf, primary_lifestyle=="unspecified_pathotroph")
mepi<-subset_taxa(fGA.Qlf, primary_lifestyle=="epiphyte")
malg.par<-subset_taxa(fGA.Qlf, primary_lifestyle=="algal_parasite")
msoot.mold<-subset_taxa(fGA.Qlf, primary_lifestyle=="sooty_mold")

# sample sums by category
maineffects<-function(ps){
  f.meta<-as.data.frame(as.matrix(sample_data(ps)))
  f.meta$Depth<-factor(f.meta$Depth, levels=c("0_5", "5_10", "10Ap", "Ap30"))
  f.meta$Treatment<-factor(f.meta$Treatment, levels=c("NT", "CT", "Org3"))
  f.meta$pH<-as.numeric(as.character(f.meta$pH))
  f.meta$C_percent<-as.numeric(as.character(f.meta$C_percent))
  f.meta$No3_ugPerg<-as.numeric(as.character(f.meta$No3_ugPerg))
  f.meta$Nh4_ugPerg<-as.numeric(as.character(f.meta$Nh4_ugPerg))
  f.meta$C_N_ratio<-as.numeric(as.character(f.meta$C_N_ratio))
  f.meta$Clay_percent<-as.numeric(as.character(f.meta$Clay_percent))
  f.meta$Silt_percent<-as.numeric(as.character(f.meta$Silt_percent))
  f.meta$SeqDepth<-as.numeric(as.character(f.meta$SeqDepth))
  
  print(summary(aov(glm(unlist(sample_sums(ps))~f.meta$SeqDepth+f.meta$Depth*f.meta$Treatment+f.meta$pH+f.meta$pH+f.meta$C_percent+f.meta$No3_ugPerg+f.meta$Nh4_ugPerg+f.meta$C_N_ratio+f.meta$Clay_percent+f.meta$Silt_percent))))
  
  boxplot(unlist(sample_sums(ps))~f.meta$Depth*f.meta$Treatment,las=2, xlab=NULL,cex.names=0.1)
  
}

maineffects(mp.path)
maineffects(msoil.sap)
maineffects(man.par)
maineffects(mwood.sap)
maineffects(mmyco.par)
maineffects(munsp.sap)
maineffects(mlitter.sap)
maineffects(mdung.sap)
maineffects(mAMF)
maineffects(mfol.end)
maineffects(mnect.sap)
maineffects(mpollen.sap)
maineffects(mlich.par)
maineffects(mECM)
maineffects(munsp.path)
maineffects(malg.par)
maineffects(mepi)
maineffects(msoot.mold)
# is growth-form associated with farming system?
# abundance of filamentous_mycelium by farming system (boxplot) + aov
# abundance of biflagellate-rhizomycelial boxplot + aov

# are there more plant pathogens?
# plant pathogens by farming system (boxplot) + aov

# are there more plant symbionts
# plant symbiont by farming system (boxplot) + aov

# does farming system affect inter and intra guild competition?
# number of significant interactions 
# directed (+/- for each):
# symbiont -> pathotroph
# pathotroph -> symbiotroph
# soil saprotroph -> pathotroph
# pathotroph -> soil saprotroph
# soil saprotroph -> symbiotroph
# symbiotroph -> pathotroph



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

c<-f.NTnet[((f.NTnet>0.3) + (-0.3 > f.NTnet))]
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


# model species by farming system ####
fGA.Qlf<-filter_taxa(fGA.Ql, function(x) sum(x > 3) > 5, TRUE)
t1<-Sys.time()
mtax.fun<-taxmodel(fGA.Qlf)
Sys.time()-t1 # 10.7 seconds

trt<-rep(NA, length(mtax.fun$tab))
for(i in c(1:length(mtax.fun$tab))){
  trt[i]<-mtax.fun$tab[[i]]$varExplained[2]
}
hist(trt, main="Variance explained by farming system") # 16 taxa trt explained greater than 30%
names(trt)<-names(mtax.fun$fit)
f.taxtab<-as.data.frame(as.matrix(tax_table(fGA.Qlf)))
fs.tab<-subset(f.taxtab, rownames(f.taxtab) %in% names(trt[trt>30]))
f.otutab<-as.data.frame(as.matrix(otu_table(fGA.Qlf)))
fs.otutab<-subset(f.otutab, rownames(f.otutab) %in% names(trt[trt>30]))
sdat<-as.data.frame(as.matrix(sample_data(fGA.Qlf)))
identical(sdat$Sample, colnames(fs.otutab))
#boxplot(unlist(fs.otutab[rownames(fs.otutab)=="GTAAAAGTCGTAACAAGGTCTCCGTTGGTGAACCAGCGGAGGGATCATTACAGGACTCGCAAGACTCCTTAAACCCCTGTGAACTTACTGTTTATACGTTGCTTCGGCGGGTGCTCCGGGGTCCGCCCCGGGGCGCTGCGCCCGCCGGCAGCCTACTTAATTCTGTTTCTCTGCGTTGGCATCTCGAGTAAGCAAAATAAGTTAAAACTTTCAACAACGGATCTCTTGGTTCTGG",])~sdat$Depth+sdat$Treatment,las=2)

to<-as.data.frame(as.matrix(tax_table(fGA.Qlf)))
toOrder<-to$primary_lifestyle
tsOrder<-to$Species
tgOrder<-to$Genus
tfOrder<-to$Growth_form_template
#if(!identical(colnames(d1), as.character(tmOrder))){stop("taxa names not match")} # sanity check
#tmOrder1<-tmOrder[order(tm$Phylum,tm$Class,tm$Order,tm$Family,tm$Genus,tm$Species)]
toOrder[is.na(toOrder)]<-"Unknown"
tsOrder[is.na(tsOrder)]<-"Unknown"
tgOrder[is.na(tgOrder)]<-"Unknown"
tfOrder[is.na(tfOrder)]<-"Unknown"

trtvar<-rep(NA, length(mtax.fun$tab))
for(i in c(1:length(mtax.fun$tab))){
  trtvar[i]<-mtax.fun$tab[[i]]$varExplained[2]
}
table(toOrder[trtvar>20])
table(tsOrder[trtvar>20])
table(tgOrder[trtvar>20])
table(tfOrder[trtvar>20])

plot_bar(fGA.Qlf, x="Treatment", fill="Growth_form_template") + geom_bar(aes(color=Growth_form_template, fill=Growth_form_template), stat="identity", position="stack")
plot_bar(fGA.Qlf, x="Treatment", fill="primary_lifestyle")+geom_bar(aes(color=primary_lifestyle, fill=primary_lifestyle), stat="identity", position="stack")

# boxplot of filamentous_mycelium by farming system

mycelium<-subset_taxa(fGA.Qlf, Growth_form_template=="filamentous_mycelium")
mycss<-sample_sums(mycelium)
mycss.fs<-as.data.frame(as.matrix(sample_data(mycelium)))$Treatment
mycss.fs<-factor(mycss.fs, levels=c("NT", "CT", "Org3"))
mycss.d<-as.data.frame(as.matrix(sample_data(mycelium)))$Depth
mycss.d<-factor(mycss.d, levels=c("0_5", "5_10", "10Ap", "Ap30"))
mycss.sd<-as.data.frame(as.matrix(sample_data(mycelium)))$SampleDepth
mycss.sd<-as.numeric(as.character(mycss.sd))
boxplot(mycss~mycss.fs+mycss.d, las=2)
stripchart(mycss~mycss.fs+mycss.d, vertical = TRUE,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')
summary(aov(glm(mycss ~ mycss.sd + mycss.d * mycss.fs)))


mycelium.o<-psmelt(mycelium)
mycelium.o$Treatment<-factor(mycelium.o$Treatment, levels=c("NT", "CT", "Org3"))
mycelium.o$Depth<-factor(mycelium.o$Depth, levels=c("0_5", "5_10", "10Ap", "Ap30"))

plot_bar(mycelium, x="Treatment", fill="Growth_form_template") + geom_bar(aes(color=Growth_form_template, fill=Growth_form_template), stat="identity", position="stack")
plot_bar(mycelium, x="Treatment", fill="primary_lifestyle", facet_grid="Depth")+geom_bar(aes(color=primary_lifestyle, fill=primary_lifestyle), stat="identity", position="stack")

mycelium.o %>% ggplot( aes(x=Treatment, y=Abundance, col=primary_lifestyle)) +
  geom_boxplot() + facet_grid(. ~ Depth) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(size=0.4, alpha=0.9) +
  #theme_ipsum() +
  theme(legend.position="bottom",
    plot.title = element_text(size=11)) + ggtitle("mycelium by farming system boxplot")


mycelium.am<-subset_taxa(fGA.Qlf, Growth_form_template=="filamentous_mycelium"&primary_lifestyle=="arbuscular_mycorrhizal")
Soilam<-sample_sums(mycelium.am)
Soilam.fs<-as.data.frame(as.matrix(sample_data(mycelium.am)))$Treatment
Soilam.fs<-factor(Soilam.fs, levels=c("NT", "CT", "Org3"))
Soilam.d<-as.data.frame(as.matrix(sample_data(mycelium.am)))$Depth
Soilam.d<-factor(Soilam.d, levels=c("0_5", "5_10", "10Ap", "Ap30"))
Soilam.sd<-as.data.frame(as.matrix(sample_data(mycelium.am)))$SampleDepth
Soilam.sd<-as.numeric(as.character(Soilam.sd))
boxplot(Soilam~Soilam.fs+Soilam.d, las=2)
stripchart(Soilam~Soilam.fs+Soilam.d, vertical = TRUE,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')
summary(aov(glm(Soilam ~ Soilam.sd + Soilam.d * Soilam.fs))) # not sure if I want to include the sequencing effort ... it partially explains depth structure because sampling effort is calculated in part by total amount od DNA to capture

#mycelium.o.am<-psmelt(mycelium.am)
#mycelium.o.am$Treatment<-factor(mycelium.o.am$Treatment, levels=c("NT", "CT", "Org3"))
#mycelium.o.am$Depth<-factor(mycelium.o.am$Depth, levels=c("0_5", "5_10", "10Ap", "Ap30"))

mycelium.o.am %>% ggplot( aes(x=Treatment, y=Abundance, col=Family.x)) +
  geom_boxplot(outlier.shape=NA) + facet_grid(. ~ Depth) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(size=0.4, alpha=0.9) +
  #theme_ipsum() +
  theme(legend.position="bottom",
        plot.title = element_text(size=11)) + ggtitle("AM fungi by farming system boxplot")

mycelium.du<-subset_taxa(fGA.Qlf, Growth_form_template=="filamentous_mycelium"&primary_lifestyle=="dung_saprotroph")
Soildu<-sample_sums(mycelium.du)
Soildu.fs<-as.data.frame(as.matrix(sample_data(mycelium.du)))$Treatment
Soildu.fs<-factor(Soildu.fs, levels=c("NT", "CT", "Org3"))
Soildu.d<-as.data.frame(as.matrix(sample_data(mycelium.du)))$Depth
Soildu.d<-factor(Soildu.d, levels=c("0_5", "5_10", "10Ap", "Ap30"))
Soildu.sd<-as.data.frame(as.matrix(sample_data(mycelium.du)))$SampleDepth
Soildu.sd<-as.numeric(as.character(Soildu.sd))
boxplot(Soildu~Soildu.fs+Soildu.d, las=2)
stripchart(Soildu~Soildu.fs+Soildu.d, vertical = TRUE,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')
summary(aov(glm(Soildu ~ Soildu.sd + Soildu.d * Soildu.fs))) # not sure if I want to include the sequencing effort ... it partially explains depth structure because sampling effort is calculated in part by total amount od DNA to capture


mycelium.p<-subset_taxa(fGA.Qlf, Growth_form_template=="filamentous_mycelium"&primary_lifestyle=="animal_parasite")
Soilp<-sample_sums(mycelium.p)
Soilp.fs<-as.data.frame(as.matrix(sample_data(mycelium.p)))$Treatment
Soilp.fs<-factor(Soilp.fs, levels=c("NT", "CT", "Org3"))
Soilp.d<-as.data.frame(as.matrix(sample_data(mycelium.p)))$Depth
Soilp.d<-factor(Soilp.d, levels=c("0_5", "5_10", "10Ap", "Ap30"))
Soilp.sd<-as.data.frame(as.matrix(sample_data(mycelium.p)))$SampleDepth
Soilp.sd<-as.numeric(as.character(Soilp.sd))
boxplot(Soilp~Soilp.fs+Soilp.d, las=2)
stripchart(Soilp~Soilp.fs+Soilp.d, vertical = TRUE,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')

summary(aov(Soilp ~ Soilp.d * Soilp.fs))

mycelium.o.p<-psmelt(mycelium.p)
mycelium.o.p$Treatment<-factor(mycelium.o.p$Treatment, levels=c("NT", "CT", "Org3"))
mycelium.o.p$Depth<-factor(mycelium.o.p$Depth, levels=c("0_5", "5_10", "10Ap", "Ap30"))

mycelium.o.p %>% ggplot( aes(x=Treatment, y=Abundance, col=Genus)) +
  geom_boxplot(outlier.shape=NA) + facet_grid(. ~ Depth) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(size=0.4, alpha=0.9) +
  #theme_ipsum() +
  theme(legend.position="bottom",
        plot.title = element_text(size=11)) + ggtitle("Animal Parasite fungi by farming system boxplot")

mycelium.sap<-subset_taxa(fGA.Qlf, Growth_form_template=="filamentous_mycelium"&primary_lifestyle=="soil_saprotroph")
SoilSap<-sample_sums(mycelium.sap)
SoilSap.fs<-as.data.frame(as.matrix(sample_data(mycelium.sap)))$Treatment
SoilSap.fs<-factor(SoilSap.fs, levels=c("NT", "CT", "Org3"))
SoilSap.d<-as.data.frame(as.matrix(sample_data(mycelium.sap)))$Depth
SoilSap.d<-factor(SoilSap.d, levels=c("0_5", "5_10", "10Ap", "Ap30"))
boxplot(SoilSap~SoilSap.fs+SoilSap.d, las=2)

mycelium.o.sap<-psmelt(mycelium.sap)
mycelium.o.sap$Treatment<-factor(mycelium.o.sap$Treatment, levels=c("NT", "CT", "Org3"))
mycelium.o.sap$Depth<-factor(mycelium.o.sap$Depth, levels=c("0_5", "5_10", "10Ap", "Ap30"))
mycelium.o.sap$sample_sums<-sample_sums(mycelium.sap)

mycelium.o.sap %>% ggplot( aes(x=Treatment, y=sample_sums)) +
  geom_boxplot(outlier.shape=NA) + facet_grid(. ~ Depth) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #geom_jitter(size=0.4, alpha=0.9) +
  #theme_ipsum() +
  theme(legend.position="bottom",
        plot.title = element_text(size=11)) + ggtitle("Soil Saprotroph fungi by farming system boxplot")

mycelium.pp<-subset_taxa(fGA.Qlf, Growth_form_template=="filamentous_mycelium"&primary_lifestyle=="plant_pathogen")
pp<-sample_sums(mycelium.pp)
pp.fs<-as.data.frame(as.matrix(sample_data(mycelium.pp)))$Treatment
pp.fs<-factor(pp.fs, levels=c("NT", "CT", "Org3"))
pp.d<-as.data.frame(as.matrix(sample_data(mycelium.pp)))$Depth
pp.d<-factor(pp.d, levels=c("0_5", "5_10", "10Ap", "Ap30"))
pp.sd<-as.data.frame(as.matrix(sample_data(mycelium.pp)))$SampleDepth
pp.sd<-as.numeric(as.character(pp.sd))
boxplot(pp~pp.fs+pp.d, las=2)
stripchart(pp~pp.fs+pp.d, vertical = TRUE,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')

summary(aov(pp ~ pp.d * pp.fs))


# boxplot of yest by farming system
yeast<-subset_taxa(fGA.Qlf, Growth_form_template=="yeast")
yss<-sample_sums(yeast)
yss.fs<-as.data.frame(as.matrix(sample_data(yeast)))$Treatment
yss.fs<-factor(yss.fs, levels=c("NT", "CT", "Org3"))
yss.d<-as.data.frame(as.matrix(sample_data(yeast)))$Depth
yss.d<-factor(yss.d, levels=c("0_5", "5_10", "10Ap", "Ap30"))
yss.sd<-as.data.frame(as.matrix(sample_data(yeast)))$SampleDepth
yss.sd<-as.numeric(as.character(yss.sd))
boxplot(yss~yss.fs+yss.d, las=2)
stripchart(yss~yss.fs+yss.d, vertical = TRUE,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')
summary(aov(glm(yss ~ yss.sd + yss.d * yss.fs)))

unique(as.data.frame(as.matrix(tax_table(yeast)))$primary_lifestyle)

# soil saprotroph
yeast.ss<-subset_taxa(fGA.Qlf, Growth_form_template=="yeast"&primary_lifestyle=="soil_saprotroph")
yeastss<-sample_sums(yeast.ss)
yeastss.fs<-as.data.frame(as.matrix(sample_data(yeast.ss)))$Treatment
yeastss.fs<-factor(yeastss.fs, levels=c("NT", "CT", "Org3"))
yeastss.d<-as.data.frame(as.matrix(sample_data(yeast.ss)))$Depth
yeastss.d<-factor(yeastss.d, levels=c("0_5", "5_10", "10Ap", "Ap30"))
yeastss.sd<-as.data.frame(as.matrix(sample_data(yeast.ss)))$SampleDepth
yeastss.sd<-as.numeric(as.character(yeastss.sd))
boxplot(yeastss~yeastss.fs+yeastss.d, las=2)
stripchart(yeastss~yeastss.fs+yeastss.d, vertical = TRUE,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')

summary(aov(glm(yeastss~yeastss.fs*yeastss.d)))



# animal parasite

yeast.p<-subset_taxa(fGA.Qlf, Growth_form_template=="yeast"&primary_lifestyle=="animal_parasite")
yeastp<-sample_sums(yeast.p)
yeastp.fs<-as.data.frame(as.matrix(sample_data(yeast.p)))$Treatment
yeastp.fs<-factor(yeastp.fs, levels=c("NT", "CT", "Org3"))
yeastp.d<-as.data.frame(as.matrix(sample_data(yeast.p)))$Depth
yeastp.d<-factor(yeastp.d, levels=c("0_5", "5_10", "10Ap", "Ap30"))
yeastp.sd<-as.data.frame(as.matrix(sample_data(yeast.p)))$SampleDepth
yeastp.sd<-as.numeric(as.character(yeastp.sd))
boxplot(yeastp~yeastp.fs+yeastp.d, las=2)
stripchart(yeastp~yeastp.fs+yeastp.d, vertical = TRUE,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')

summary(aov(glm(yeastp~yeastp.fs*yeastp.d)))

# mycoparasite

yeast.mp<-subset_taxa(fGA.Qlf, Growth_form_template=="yeast"&primary_lifestyle=="mycoparasite")
yeastmp<-sample_sums(yeast.mp)
yeastmp.fs<-as.data.frame(as.matrix(sample_data(yeast.mp)))$Treatment
yeastmp.fs<-factor(yeastmp.fs, levels=c("NT", "CT", "Org3"))
yeastmp.d<-as.data.frame(as.matrix(sample_data(yeast.mp)))$Depth
yeastmp.d<-factor(yeastmp.d, levels=c("0_5", "5_10", "10Ap", "Ap30"))
yeastmp.sd<-as.data.frame(as.matrix(sample_data(yeast.mp)))$SampleDepth
yeastmp.sd<-as.numeric(as.character(yeastmp.sd))
boxplot(yeastmp~yeastmp.fs+yeastmp.d, las=2)
stripchart(yeastmp~yeastmp.fs+yeastmp.d, vertical = TRUE,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')

summary(aov(glm(yeastmp~yeastmp.fs*yeastmp.d)))

# litter saprotroph

yeast.ls<-subset_taxa(fGA.Qlf, Growth_form_template=="yeast"&primary_lifestyle=="litter_saprotroph")
yeastls<-sample_sums(yeast.ls)
yeastls.fs<-as.data.frame(as.matrix(sample_data(yeast.ls)))$Treatment
yeastls.fs<-factor(yeastmp.fs, levels=c("NT", "CT", "Org3"))
yeastls.d<-as.data.frame(as.matrix(sample_data(yeast.ls)))$Depth
yeastls.d<-factor(yeastmp.d, levels=c("0_5", "5_10", "10Ap", "Ap30"))
yeastls.sd<-as.data.frame(as.matrix(sample_data(yeast.ls)))$SampleDepth
yeastls.sd<-as.numeric(as.character(yeastmp.sd))
boxplot(yeastls~yeastls.fs+yeastls.d, las=2)
stripchart(yeastls~yeastls.fs+yeastls.d, vertical = TRUE,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')

summary(aov(glm(yeastls~yeastls.fs*yeastls.d)))

# epiphyte
yeast.e<-subset_taxa(fGA.Qlf, Growth_form_template=="yeast"&primary_lifestyle=="epiphyte")
yeaste<-sample_sums(yeast.e)
yeaste.fs<-as.data.frame(as.matrix(sample_data(yeast.e)))$Treatment
yeaste.fs<-factor(yeaste.fs, levels=c("NT", "CT", "Org3"))
yeaste.d<-as.data.frame(as.matrix(sample_data(yeast.e)))$Depth
yeaste.d<-factor(yeaste.d, levels=c("0_5", "5_10", "10Ap", "Ap30"))
yeaste.sd<-as.data.frame(as.matrix(sample_data(yeast.e)))$SampleDepth
yeaste.sd<-as.numeric(as.character(yeaste.sd))
boxplot(yeaste~yeaste.fs+yeaste.d, las=2)
stripchart(yeaste~yeaste.fs+yeaste.d, vertical = TRUE,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')

summary(aov(glm(yeaste~yeaste.fs*yeaste.d)))
# nectar ...



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
d1<-as.data.frame(t(as.matrix(otu_table(fGA.Qf.5)))) 
d2<-as.data.frame(as.matrix(sample_data(fGA.Qf.5)))
effort<-as.numeric(as.character(d2$SampleDepth))
treatment<-as.factor(d2$Treatment)
depth<-as.factor(d2$Depth)
pH<-as.numeric(as.character(d2$pH))
Cpercent<-as.numeric(as.character(d2$C_percent))
Ammonia<-as.numeric(as.character(d2$Nh4_ugPerg))
Nitrate<-as.numeric(as.character(d2$No3_ugPerg))
 
fit<-glm(d1[,1]~effort+treatment+pH+Cpercent+d1[,2])
fit$coefficients

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
get.model<-function(v, d2){
    fit<-NULL
    tab<-NULL
    treatment<-as.factor(d2$Treatment)
    effort<-as.numeric(as.character(d2$SampleDepth))
    depth<-as.factor(d2$Depth) # has 4 levels (intercept + 3 levels, then count rest)
    pH<-as.numeric(as.character(d2$pH))
    Cpercent<-as.numeric(as.character(d2$C_percent))
    Ammonia<-as.numeric(as.character(d2$Nh4_ugPerg))
    Nitrate<-as.numeric(as.character(d2$No3_ugPerg))
    Clay<-as.numeric(as.character(d2$Clay_percent))
    Silt<-as.numeric(as.character(d2$Silt_percent))
    fit <- glm(tax1 ~ effort+treatment+depth+pH+Cpercent+Ammonia+Nitrate+Clay+Silt,family="gaussian")
    fit <- glm(tax1 ~ treatment*depth+pH+Cpercent+Ammonia+Nitrate+tax2,family="gaussian")
    tab<-cbind(summary(aov(fit))[[1]],"varExplained"=summary(aov(fit))[[1]]$`Sum Sq`/sum(summary(aov(fit))[[1]]$`Sum Sq`)*100)
    
  o$fit<-fit
  o$tab<-tab
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