# FSP GAD fungal analysis 

# load libraries, functions:
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

# Models the relationship between taxon and environment with treatment factors
# outputs:  $fit = glm object
#           $tab = significance table
#           $residuals = residuals from the model
#           $predicted = model predicted values
#           $predcor = correlation matrix of the predicted values for each spp
#           $spcor = species correlation matrix
taxmodel<-function(d){ 
  # input is phyloseq object
  # already normalized
  d1<-as.data.frame(t(as.matrix(otu_table(d)))) # dim =
  d2<-as.data.frame(as.matrix(sample_data(d)))# dataframes must be samples as rows
  d3<-as.data.frame(as.matrix(tax_table(d)))
  if(!identical(rownames(d1), rownames(d2))){stop("dataframe orientation does not match")}
  if(!identical(colnames(d1), rownames(d3))){stop("tax table not in correct order")}
  treatment<-as.factor(d2$Treatment)
  effort<-as.numeric(as.character(d2$SeqDepth))
  depth<-as.factor(d2$Depth)
  depth<-factor(depth, levels=c("0_5", "5_10", "10Ap", "Ap30"))
  treatment<-factor(treatment, levels=c("NT", "CT", "Org3"))
  pH<-as.numeric(as.character(d2$pH))
  Cpercent<-as.numeric(as.character(d2$C_percent))
  Ammonia<-as.numeric(as.character(d2$Nh4_ugPerg))
  Nitrate<-as.numeric(as.character(d2$No3_ugPerg))
  Clay<-as.numeric(as.character(d2$Clay_percent))
  Silt<-as.numeric(as.character(d2$Silt_percent))
  CN_ratio<-as.numeric(as.character(d2$C_N_ratio))
  N_percent<-as.numeric(as.character(d2$N_percent))
  Sand<-as.numeric(as.character(d2$Sand_percent))
  B.Density<-as.numeric(as.character(d2$B.Density_gcm3))
  species<-as.character(d3$Species)
  out<-NULL
  out$residuals<-matrix(ncol=ncol(d1), nrow=length(predict(glm(d1[,1] ~ effort+treatment*depth+pH+Cpercent*N_percent+Ammonia+Nitrate+Clay+Silt+Sand+B.Density,family="gaussian"))))
  out$predicted<-matrix(ncol=ncol(d1), nrow=length(predict(glm(d1[,1] ~ effort+treatment*depth+pH+Cpercent*N_percent+Ammonia+Nitrate+Clay+Silt+Sand+B.Density,family="gaussian"))))
  out$fit<-list(1:length(ncol(d1)))
  out$plots<-list(1:length(ncol(d1)))
  colnames(out$predicted)<-colnames(d1)
  for(i in c(1:ncol(d1))){
    fit<-NULL
    fit <- glm(d1[,i] ~ effort+treatment*pH+Cpercent*N_percent+Ammonia+Nitrate+Clay+Silt+Sand+B.Density,family="gaussian")
    
    out$residuals[,i]<-residuals(fit)
    out$predicted[,i]<-predict(fit)
    out$fit[[i]]<-fit
    out$tab[[i]]<-cbind(summary(aov(fit))[[1]],"varExplained"=summary(aov(fit))[[1]]$`Sum Sq`/sum(summary(aov(fit))[[1]]$`Sum Sq`)*100)
    out$plots[[i]]<-boxplot(d1[,i]~depth + treatment, las=2, xlab=NULL,cex.names=0.1, ylab="Abundance",main=species[i], names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8))
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
  names(out$plots)<-as.character(d3$Species)
  
  out$effectDF<-matrix(nrow=ncol(d1), ncol=2)
  rownames(out$effectDF)<-as.character(d3$Species)
  colnames(out$effectDF)<-c("CT", "ORG")
  for(i in 1:nrow(out$effectDF)){
    out$effectDF[i,1]<-out$fit[[i]]$coefficients["treatmentCT"]
    out$effectDF[i,2]<-out$fit[[i]]$coefficients["treatmentOrg3"]
    
    # out$effectDF[i,1]<-out$fit[[i]]$coefficients["treatmentCT"]/out$fit[[i]]$coefficients["(Intercept)"]
    #  out$effectDF[i,2]<-out$fit[[i]]$coefficients["treatmentOrg3"]/out$fit[[i]]$coefficients["(Intercept)"]
  }
  
  
  out
}
# Models the relationship between taxon and environment without treatment factors
# outputs:  $fit = glm object
#           $tab = significance table
#           $residuals = residuals from the model
#           $predicted = model predicted values
#           $predcor = correlation matrix of the predicted values for each spp
#           $spcor = species correlation matrix
taxEnvmodel<-function(d){ 
  # input is phyloseq object
  # already normalized
  d1<-as.data.frame(t(as.matrix(otu_table(d)))) # dim =
  d2<-as.data.frame(as.matrix(sample_data(d)))# dataframes must be samples as rows
  if(!identical(rownames(d1), rownames(d2))){stop("dataframe orientation does not match")}
  effort<-as.numeric(as.character(d2$SeqDepth))
  pH<-as.numeric(as.character(d2$pH))
  Cpercent<-as.numeric(as.character(d2$C_percent))
  Ammonia<-as.numeric(as.character(d2$Nh4_ugPerg))
  Nitrate<-as.numeric(as.character(d2$No3_ugPerg))
  Clay<-as.numeric(as.character(d2$Clay_percent))
  Silt<-as.numeric(as.character(d2$Silt_percent))
  CN_ratio<-as.numeric(as.character(d2$C_N_ratio))
  N_percent<-as.numeric(as.character(d2$N_percent))
  Sand<-as.numeric(as.character(d2$Sand_percent))
  B.Density<-as.numeric(as.character(d2$B.Density_gcm3))
  out<-NULL
  out$residuals<-matrix(ncol=ncol(d1), nrow=length(predict(glm(d1[,1] ~ effort+pH+Cpercent*N_percent+Ammonia+Nitrate+Clay+Silt+Sand+B.Density,family="gaussian"))))
  out$predicted<-matrix(ncol=ncol(d1), nrow=length(predict(glm(d1[,1] ~ effort+pH+Cpercent*N_percent+Ammonia+Nitrate+Clay+Silt+Sand+B.Density,family="gaussian"))))
  out$fit<-list(1:length(ncol(d1)))
  colnames(out$predicted)<-colnames(d1)
  for(i in c(1:ncol(d1))){
    fit<-NULL
    fit <- glm(d1[,i] ~ effort+pH+Cpercent*N_percent+Ammonia+Nitrate+Clay+Silt+Sand+B.Density,family="gaussian")
    
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

# for interaction network, whole
# still needs effort in it!!
get.no<-function(v, d2){
  tax1<-v
  o<-c(rep(NA,length(d1)))
  for(i in 1:length(d1)){
    tax2<-NULL
    fit<-NULL
    tab<-NULL
    tax2<-d1[,i]
    effort<-as.numeric(as.character(d2$SeqDepth))
    treatment<-as.factor(d2$Treatment)
    depth<-as.factor(d2$Depth)
    pH<-as.numeric(as.character(d2$pH))
    Cpercent<-as.numeric(as.character(d2$C_percent))
    Ammonia<-as.numeric(as.character(d2$Nh4_ugPerg))
    Nitrate<-as.numeric(as.character(d2$No3_ugPerg))
    Clay<-as.numeric(as.character(d2$Clay_percent))
    Silt<-as.numeric(as.character(d2$Silt_percent))
    CN_ratio<-as.numeric(as.character(d2$C_N_ratio))
    N_percent<-as.numeric(as.character(d2$N_percent))
    Sand<-as.numeric(as.character(d2$Sand_percent))
    B.Density<-as.numeric(as.character(d2$B.Density_gcm3))
    fit <- glm(tax1 ~ treatment*depth+pH+Cpercent*N_percent+Ammonia+Nitrate+Clay+Silt+Sand+B.Density+tax2,family="gaussian")
    tab<-cbind(summary(aov(fit))[[1]],"varExplained"=summary(aov(fit))[[1]]$`Sum Sq`/sum(summary(aov(fit))[[1]]$`Sum Sq`)*100)
    o[i]<-tab$varExplained[13]*(fit$coefficients[17]/abs(fit$coefficients[17]))/100
    
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

# for interaction network data sliced by farming system
# get.no2 MIGHT NEED TO BE UPDATED FOR DEGREES OF FREEDOM (REDUCE FACTORS!!)
get.no2<-function(v, d2, d1){
  tax1<-v
  o<-c(rep(NA,length(d1)))
  for(i in 1:length(d1)){
    tax2<-NULL
    fit<-NULL
    tab<-NULL
    tax2<-d1[,i]
    #treatment<-as.factor(d2$Treatment)
    effort<-as.numeric(as.character(d2$SeqDepth))
    depth<-as.factor(d2$Depth) # has 4 levels (intercept + 3 levels, then count rest)
    pH<-as.numeric(as.character(d2$pH))
    Cpercent<-as.numeric(as.character(d2$C_percent))
    Ammonia<-as.numeric(as.character(d2$Nh4_ugPerg))
    Nitrate<-as.numeric(as.character(d2$No3_ugPerg))
    Clay<-as.numeric(as.character(d2$Clay_percent))
    Silt<-as.numeric(as.character(d2$Silt_percent))
    CN_ratio<-as.numeric(as.character(d2$C_N_ratio))
    N_percent<-as.numeric(as.character(d2$N_percent))
    Sand<-as.numeric(as.character(d2$Sand_percent))
    B.Density<-as.numeric(as.character(d2$B.Density_gcm3))
    fit <- glm(tax1 ~ effort+depth+pH+Cpercent*N_percent+Ammonia+Nitrate+Clay+Silt+Sand+B.Density+tax2,family="gaussian")
    tab<-cbind(summary(aov(fit))[[1]],"varExplained"=summary(aov(fit))[[1]]$`Sum Sq`/sum(summary(aov(fit))[[1]]$`Sum Sq`)*100)
    o[i]<-tab$varExplained[12]*(fit$coefficients[15]/abs(fit$coefficients[15]))/100
    
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
  
  out<-do.call(cbind, mclapply(d1, get.no2, d2, d1, mc.preschedule = T, mc.cores=10, mc.cleanup = T))
  rownames(out)<-colnames(out)
  tm<-as.data.frame(as.matrix(tax_table(d)))
  tmOrder<-rownames(tm) #tm$Class, tm$Order, 
  if(!identical(colnames(d1), as.character(tmOrder))){stop("taxa names not match")} # sanity check
  tmOrder1<-tmOrder[order(tm$Phylum,tm$Class,tm$Order,tm$Family,tm$Genus,tm$Species)]
  out<-out[as.character(tmOrder1), as.character(tmOrder1)]
  out
}

# for interaction network with data sliced by depth
# get.no3 MIGHT NEED TO BE UPDATED FOR DEGREES OF FREEDOM (REDUCE FACTORS!!)
get.no3<-function(v, d2, d1){
  tax1<-v
  o<-c(rep(NA,length(d1)))
  for(i in 1:length(d1)){
    tax2<-NULL
    fit<-NULL
    tab<-NULL
    tax2<-d1[,i]
    #treatment<-as.factor(d2$Treatment)
    effort<-as.numeric(as.character(d2$SeqDepth))
    trt<-as.factor(d2$Treatment)
    pH<-as.numeric(as.character(d2$pH))
    Cpercent<-as.numeric(as.character(d2$C_percent))
    Ammonia<-as.numeric(as.character(d2$Nh4_ugPerg))
    Nitrate<-as.numeric(as.character(d2$No3_ugPerg))
    Clay<-as.numeric(as.character(d2$Clay_percent))
    Silt<-as.numeric(as.character(d2$Silt_percent))
    CN_ratio<-as.numeric(as.character(d2$C_N_ratio))
    N_percent<-as.numeric(as.character(d2$N_percent))
    Sand<-as.numeric(as.character(d2$Sand_percent))
    B.Density<-as.numeric(as.character(d2$B.Density_gcm3))
    fit <- glm(tax1 ~ effort+trt+pH+Cpercent*N_percent+Ammonia+Nitrate+Clay+Silt+Sand+B.Density+tax2,family="gaussian")
    tab<-cbind(summary(aov(fit))[[1]],"varExplained"=summary(aov(fit))[[1]]$`Sum Sq`/sum(summary(aov(fit))[[1]]$`Sum Sq`)*100)
    o[i]<-tab$varExplained[12]*(fit$coefficients[14]/abs(fit$coefficients[14]))/100 # count the position of tax2 in the effects; the position of var explained = num of levels in factor + number of linear covariates ...
    
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

# import and format data
GAD.Fun<-readRDS("/Users/dietrich/Documents/GitHub/Anchoring/Data/GAD/GADfun2020.rds")
#sample_data(GAD.bac)$SeqDepth<-sample_sums(GAD.bac)
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

GAD.Fun<-subset_taxa(GAD.Fun, !is.na(Genus)) # remove anything not ID'd to Genus
# now give all unknown spp a same name:
tt1<-as.data.frame(as.matrix(tax_table(GAD.Fun)))
tt1$Species[is.na(tt1$Species)]<-paste0(tt1$Genus[is.na(tt1$Species)], "Undefined", sep="")
tax_table(GAD.Fun)<-tax_table(as.matrix(tt1))
# now aggregate to species level
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
#bGA.Q<-GAD.QScale(GAD.bac,type="B")

# filter samples to top 4 levels
fGA.Ql<-subset_samples(fGA.Q, Depth=="0_5"|Depth=="5_10"|Depth=="10Ap"|Depth=="Ap30") 
#bGA.Ql<-subset_samples(bGA.Q, Depth=="0_5"|Depth=="5_10"|Depth=="10Ap"|Depth=="Ap30")

# covariate analysis ####
f.meta<-as.data.frame(as.matrix(sample_data(fGA.Ql)))
f.meta$Fun_QPCR<-as.numeric(as.character(f.meta$Fun_QPCR))
f.meta$Depth<-factor(f.meta$Depth, levels=c("0_5", "5_10", "10Ap", "Ap30"))
f.meta$Treatment<-factor(f.meta$Treatment, levels=c("NT", "CT", "Org3"))
f.meta$pH<-as.numeric(as.character(f.meta$pH))
f.meta$C_percent<-as.numeric(as.character(f.meta$C_percent))
f.meta$No3_ugPerg<-as.numeric(as.character(f.meta$No3_ugPerg))
f.meta$Nh4_ugPerg<-as.numeric(as.character(f.meta$Nh4_ugPerg))
f.meta$C_N_ratio<-as.numeric(as.character(f.meta$C_N_ratio))
f.meta$N_percent<-as.numeric(as.character(f.meta$N_percent))
f.meta$Clay_percent<-as.numeric(as.character(f.meta$Clay_percent))
f.meta$Silt_percent<-as.numeric(as.character(f.meta$Silt_percent))
f.meta$Sand_percent<-as.numeric(as.character(f.meta$Sand_percent))
f.meta$B.Density_gcm3<-as.numeric(as.character(f.meta$B.Density_gcm3))
f.meta$SeqDepth<-as.numeric(as.character(f.meta$SeqDepth))
f.meta$SampleDepth<-as.numeric(as.character(f.meta$SampleDepth))
# analysis of correltion of seqdepth with other factors:

with(f.meta, summary(aov(glm(SampleDepth~Depth))))


cor(f.meta$Clay_percent, f.meta$Silt_percent)
cor(f.meta$Clay_percent, f.meta$Sand_percent)
cor(f.meta$Sand_percent, f.meta$Silt_percent)
cor(f.meta$Sand_percent, f.meta$B.Density_gcm3)
cor(f.meta$Clay_percent,f.meta$B.Density_gcm3)

plot(f.meta$Clay_percent, f.meta$Silt_percent)
plot(f.meta$Clay_percent, f.meta$Sand_percent)
plot(f.meta$Sand_percent, f.meta$Silt_percent)

plot(f.meta$Sand_percent, f.meta$B.Density_gcm3)
plot(f.meta$Clay_percent,f.meta$B.Density_gcm3)

# check effect of farming system and depth on total fungi
boxplot(f.meta$Fun_QPCR~f.meta$Depth+f.meta$Treatment,las=2, xlab=NULL,cex.names=0.1,main="Total Fungi QPCR", ylab="Gene Count", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8))
summary(aov(f.meta$Fun_QPCR~f.meta$Depth*f.meta$Treatment))
# this shows that there is significant effects of depth, trt;
# also supports interaction between depth and trt

# colinearity:
summary(aov(f.meta$pH~f.meta$Depth)) # 0.433
summary(aov(f.meta$C_percent~f.meta$Depth)) # 0.870
summary(aov(f.meta$No3_ugPerg~f.meta$Depth)) # 0.466
summary(aov(f.meta$Nh4_ugPerg~f.meta$Depth)) # NS: 0.11
summary(aov(f.meta$C_N_ratio~f.meta$Depth)) # 0.448
summary(aov(f.meta$Clay_percent~f.meta$Depth)) # 0.76
summary(aov(f.meta$Silt_percent~f.meta$Depth)) # 0.484
summary(aov(f.meta$Sand_percent~f.meta$Depth)) # 0.39
summary(aov(f.meta$N_percent~f.meta$Depth)) # 0.85007
summary(aov(f.meta$B.Density_gcm3~f.meta$Depth)) # 0.61

summary(aov(f.meta$pH~f.meta$Depth*f.meta$Treatment)) 
aov(f.meta$pH~f.meta$Depth*f.meta$Treatment)$coefficients
boxplot(f.meta$pH~f.meta$Depth+f.meta$Treatment, las=2, xlab=NULL,cex.names=0.1,main="Soil pH", ylab="Soil pH", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8))

summary(aov(f.meta$C_percent~f.meta$Depth*f.meta$Treatment)) 
aov(f.meta$C_percent~f.meta$Depth*f.meta$Treatment)$coefficients
boxplot(f.meta$C_percent~f.meta$Depth+f.meta$Treatment, las=2, xlab=NULL,cex.names=0.1,main="C Percent", ylab="C Percent", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8))

summary(aov(f.meta$No3_ugPerg~f.meta$Depth*f.meta$Treatment)) 
aov(f.meta$No3_ugPerg~f.meta$Depth*f.meta$Treatment)$coefficients
boxplot(f.meta$No3_ugPerg~f.meta$Depth+f.meta$Treatment, las=2, xlab=NULL,cex.names=0.1,main="Nitrate", ylab="Nitrate ug/g", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8))

summary(aov(f.meta$Nh4_ugPerg~f.meta$Depth*f.meta$Treatment)) 
aov(f.meta$Nh4_ugPerg~f.meta$Depth*f.meta$Treatment)$coefficients
boxplot(f.meta$Nh4_ugPerg~f.meta$Depth+f.meta$Treatment, las=2, xlab=NULL,cex.names=0.1,main="Ammonia", ylab="Ammonia ug/g", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8))

summary(aov(f.meta$C_N_ratio~f.meta$Depth*f.meta$Treatment)) 
aov(f.meta$C_N_ratio~f.meta$Depth*f.meta$Treatment)$coefficients
boxplot(f.meta$C_N_ratio~f.meta$Depth+f.meta$Treatment, las=2, xlab=NULL,cex.names=0.1,main="C:N Ratio", ylab="C:N ratio", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8))

summary(aov(f.meta$Clay_percent~f.meta$Depth*f.meta$Treatment)) 
aov(f.meta$Clay_percent~f.meta$Depth*f.meta$Treatment)$coefficients
boxplot(f.meta$Clay_percent~f.meta$Depth+f.meta$Treatment, las=2, xlab=NULL,cex.names=0.1,main="Clay Percent", ylab="Clay Percent", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8))

summary(aov(f.meta$Silt_percent~f.meta$Depth*f.meta$Treatment)) 
aov(f.meta$Silt_percent~f.meta$Depth*f.meta$Treatment)$coefficients
boxplot(f.meta$Silt_percent~f.meta$Depth+f.meta$Treatment, las=2, xlab=NULL,cex.names=0.1,main="Percent Silt", ylab="Percent Silt", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8))

summary(aov(f.meta$Sand_percent~f.meta$Depth*f.meta$Treatment))
aov(f.meta$Sand_percent~f.meta$Depth*f.meta$Treatment)$coefficients
boxplot(f.meta$Sand_percent~f.meta$Depth+f.meta$Treatment, las=2, xlab=NULL,cex.names=0.1,main="Percent Sand", ylab="Percent Sand", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8))

summary(aov(f.meta$N_percent~f.meta$Depth*f.meta$Treatment)) 
aov(f.meta$N_percent~f.meta$Depth*f.meta$Treatment)$coefficients
boxplot(f.meta$N_percent~f.meta$Depth+f.meta$Treatment, las=2, xlab=NULL,cex.names=0.1,main="N Percent", ylab="Percent Nitrogen", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8))

summary(aov(f.meta$B.Density_gcm3~f.meta$Depth*f.meta$Treatment)) 
aov(f.meta$B.Density_gcm3~f.meta$Depth*f.meta$Treatment)$coefficients
boxplot(f.meta$B.Density_gcm3~f.meta$Depth+f.meta$Treatment, las=2, xlab=NULL,cex.names=0.1,main="Bulk Density", ylab="Bulk Density", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8))

# alpha diversity analysis (all of it) ####
a.ffit<-lm(unlist(estimate_richness(fGA.Ql, measures="Observed"))~sample_data(fGA.Ql)$SeqDepth)
#a.bfit<-lm(unlist(estimate_richness(bGA.Ql, measures="Observed"))~sample_sums(bGA.Ql))
summary(aov(a.ffit))  # by seq depth = 17.8% of variance explained at spp level

identical(rownames(f.meta), rownames(estimate_richness(fGA.Ql, measures="Observed")))
resids<-a.ffit$residuals
alphadiveffect<-summary(aov(lm(resids~f.meta$Depth*f.meta$Treatment)))
plot(fitted(a.ffit), resid(a.ffit))
abline(0,0)
#boxplot
boxplot(resids~f.meta$Depth+f.meta$Treatment, las=2, xlab=NULL,cex.names=0.1,main="Fungal Alpha Diversity", ylab="Sequencing Depth Residual", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8))
abline(h = 0, col = 'black', lty=c(2)) 

# just sequencing effort (for dissertation):
# run as log because it reduces the structure of the residuals ...
a.ffit2<-lm(unlist(estimate_richness(fGA.Ql, measures="Observed"))~log(sample_data(fGA.Ql)$SampleDepth)) # by relative sampling effort
plot(fitted(a.ffit2), resid(a.ffit2)) # clear structure in residuals
abline(0,0)

summary(aov(a.ffit2)) # by sampleing effort = 24.1% of variance explained
# sampling effort explains 56 % of log transformed variance
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


# test sources of alph diversity by functional guild
# what does alpha diversity mean in an organic farming system?
# this data is not filtered as needed for abundance data       

# sanity check: list uniuqe primary liftestyles:
unique(as.data.frame(as.matrix(tax_table(fGA.Ql)))$primary_lifestyle)

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
root.endophyte<-subset_taxa(fGA.Ql, primary_lifestyle=="root_endophyte")
pollen.sap<-subset_taxa(fGA.Ql, primary_lifestyle=="pollen_saprotroph")
lich.par<-subset_taxa(fGA.Ql, primary_lifestyle=="lichen_parasite")
ECM<-subset_taxa(fGA.Ql, primary_lifestyle=="ectomycorrhizal")
unsp.path<-subset_taxa(fGA.Ql, primary_lifestyle=="unspecified_pathotroph")
epi<-subset_taxa(fGA.Ql, primary_lifestyle=="epiphyte")
alg.par<-subset_taxa(fGA.Ql, primary_lifestyle=="algal_parasite")
soot.mold<-subset_taxa(fGA.Ql, primary_lifestyle=="sooty_mold")
unknown<-subset_taxa(fGA.Ql, primary_lifestyle=="unspecified")
lichenized<-subset_taxa(fGA.Ql, primary_lifestyle=="lichenized")
epiphyte<-subset_taxa(fGA.Ql, primary_lifestyle=="epiphyte")
myc<-subset_taxa(fGA.Ql, Growth_form_template=="filamentous_mycelium")
yeast<-subset_taxa(fGA.Ql, Growth_form_template=="yeast")

# richness of each group
p.path.r<-unlist(estimate_richness(p.path, measures="Observed"))
soil.sap.r<-unlist(estimate_richness(soil.sap, measures="Observed"))
an.par.r<-unlist(estimate_richness(an.par, measures="Observed"))
wood.sap.r<-unlist(estimate_richness(wood.sap, measures="Observed"))
unknown.r<-unlist(estimate_richness(unknown, measures="Observed"))
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
root.endophyte.r<-unlist(estimate_richness(root.endophyte, measures="Observed"))
epiphyte.r<-unlist(estimate_richness(epiphyte, measures="Observed"))
lichenized.r<-unlist(estimate_richness(lichenized, measures="Observed"))
myc.r<-unlist(estimate_richness(myc, measures="Observed"))
yeast.r<-unlist(estimate_richness(yeast, measures="Observed"))
tot.r<-unlist(estimate_richness(fGA.Ql, measures="Observed"))


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
                   yeast.r,
                   myc.r,
                   tot.r,
                   lichenized.r,
                   unknown.r,
                   root.endophyte.r,
                   epiphyte.r,
                   "depth"=f.meta$Depth,
                   "treatment"= f.meta$Treatment,
                   "seqs"=as.numeric(as.character(f.meta$SeqDepth)))

# plot and statistics of abundance of each group by sequences, depth and treatment

with(div.df,boxplot(p.path.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Plant Pathogen Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(p.path.r~seqs+depth*treatment))))

with(div.df,boxplot(soil.sap.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Soil Saprotroph Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(soil.sap.r~seqs+depth*treatment))))
with(div.df,TukeyHSD(aov(lm(soil.sap.r~seqs+depth*treatment))))

with(div.df,boxplot(an.par.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Animal Parasite Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(an.par.r~seqs+depth*treatment))))

with(div.df,boxplot(wood.sap.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Wood Saprotroph Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(wood.sap.r~seqs+depth*treatment))))
with(div.df,TukeyHSD(aov(lm(wood.sap.r~seqs+depth*treatment))))

#unique(as.data.frame(as.matrix(tax_table(fGA.Ql)))$primary_lifestyle)
#estimate_richness(fGA.Ql, measures="Observed")

with(div.df,boxplot(myco.par.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Mycoparasite Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(myco.par.r~seqs+depth*treatment))))

with(div.df,boxplot(unsp.sap.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Unspecified Saprotroph Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(unsp.sap.r~seqs+depth*treatment))))
with(div.df,TukeyHSD(aov(lm(unsp.sap.r~seqs+depth*treatment))))
p<-plot_bar(unsp.sap, x="Treatment",facet_grid = "Depth", fill="Family.x")+ geom_bar(aes(color=Family.x, fill=Family.x), stat="identity", position="stack")
p

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

with(div.df,boxplot(unknown.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Unknown Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(unknown.r~seqs+depth*treatment))))

with(div.df,boxplot(root.endophyte.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Root Endophyte Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(root.endophyte.r~seqs+depth*treatment))))

with(div.df,boxplot(lichenized.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="lichenized Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(lichenized.r~seqs+depth*treatment))))

with(div.df,boxplot(myc.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="mycelium Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(myc.r~seqs+depth*treatment))))
with(div.df,TukeyHSD(aov(lm(myc.r~seqs+depth*treatment))))
with(div.df,boxplot(yeast.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Yeast Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(yeast.r~seqs+depth*treatment))))

with(div.df,boxplot(tot.r~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Total Diversity", ylab="Taxa", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8)))
with(div.df,summary(aov(lm(tot.r~seqs+depth*treatment))))

# organism abundnace analysis ####

fGA.Qlf<-filter_taxa(fGA.Ql, function(x) sum(x > 3) > 10, TRUE)

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
#munknown<-subset_taxa(fGA.Qlf, primary_lifestyle=="unspecified")
#mlichenized<-subset_taxa(fGA.Qlf, primary_lifestyle=="lichenized")
mepiphyte<-subset_taxa(fGA.Qlf, primary_lifestyle=="epiphyte")
# sample sums by category
# this function makes a plot of abundance by farming system and depth
# it also runs an anova to test difference in abundance by factor categories
# sequencing depth is included as a covariate

maineffects<-function(ps){
  f.meta<-as.data.frame(as.matrix(sample_data(ps)))
  f.meta$SeqDepth<-as.numeric(as.character(f.meta$SeqDepth))
  f.meta$Depth<-factor(f.meta$Depth, levels=c("0_5", "5_10", "10Ap", "Ap30"))
  f.meta$Treatment<-factor(f.meta$Treatment, levels=c("NT", "CT", "Org3"))
  f.meta$pH<-as.numeric(as.character(f.meta$pH))
  f.meta$C_percent<-as.numeric(as.character(f.meta$C_percent))
  f.meta$No3_ugPerg<-as.numeric(as.character(f.meta$No3_ugPerg))
  f.meta$Nh4_ugPerg<-as.numeric(as.character(f.meta$Nh4_ugPerg))
  f.meta$C_N_ratio<-as.numeric(as.character(f.meta$C_N_ratio))
  f.meta$N_percent<-as.numeric(as.character(f.meta$N_percent))
  f.meta$Clay_percent<-as.numeric(as.character(f.meta$Clay_percent))
  f.meta$Silt_percent<-as.numeric(as.character(f.meta$Silt_percent))
  f.meta$Sand_percent<-as.numeric(as.character(f.meta$Sand_percent))
  f.meta$B.Density_gcm3<-as.numeric(as.character(f.meta$B.Density_gcm3))
  
  
  print(summary(aov(glm(unlist(sample_sums(ps))~f.meta$SeqDepth+f.meta$Depth*f.meta$Treatment+f.meta$pH+f.meta$C_percent*f.meta$N_percent+f.meta$No3_ugPerg+f.meta$Nh4_ugPerg+f.meta$Clay_percent+f.meta$Sand_percent+f.meta$Silt_percent+f.meta$B.Density_gcm3))))
  
  boxplot(unlist(sample_sums(ps))~f.meta$Depth*f.meta$Treatment,las=2, xlab=NULL,cex.names=0.1)
  
  #  plot_bar(ps, x="Treatment", fill="Species", facet_grid="Depth")+geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack")
}


# determine significance of treatment...
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
maineffects(mepiphyte)
t1<-Sys.time()
emp.path<-taxmodel(mp.path)
emsoil.sap<-taxmodel(msoil.sap)
eman.par<-taxmodel(man.par)
emwood.sap<-taxmodel(mwood.sap)
emmyco.par<-taxmodel(mmyco.par)
emunsp.sap<-taxmodel(munsp.sap)
emlitter.sap<-taxmodel(mlitter.sap)
emdung.sap<-taxmodel(mdung.sap)
emAMF<-taxmodel(mAMF)
emfol.end<-taxmodel(mfol.end)
emnect.sap<-taxmodel(mnect.sap)
empollen.sap<-taxmodel(mpollen.sap)
emlich.par<-taxmodel(mlich.par)
emECM<-taxmodel(mECM)
emunsp.path<-taxmodel(munsp.path)
emalg.par<-taxmodel(malg.par)
emepi<-taxmodel(mepi)
emsoot.mold<-taxmodel(msoot.mold)
Sys.time()-t1

# plot the relative effect size for each of the significant taxa
# significant taxa by farming system
# ratio of coefficients:
# CT/Intercept or Org3/Intercept <- double check this! I changed to try effect size to determine which demonstrated the effects better
# alpha = 0.05
plotdiff<-function(x){
  n<-c(rep(0, length(x$tab)))
  for(i in 1:length(x$tab)){
    n[i]<-x$tab[[i]]$`Pr(>F)`[2]
  }
  par(mar=c(10, 4, 4, 2) + 0.1)
  barplot(t(x$effectDF[n<0.05,]), beside = T, las=2)
}

plotdiff(emp.path)
plotdiff(emsoil.sap)
plotdiff(eman.par)
plotdiff(emwood.sap)
plotdiff(emmyco.par)
plotdiff(emunsp.sap)
plotdiff(emlitter.sap)
plotdiff(emdung.sap)
plotdiff(emAMF)
plotdiff(emfol.end) 
plotdiff(emnect.sap)
plotdiff(empollen.sap)
plotdiff(emlich.par) #none 
plotdiff(emECM) # none
plotdiff(emunsp.path) # only one significant
plotdiff(emalg.par) # only one sigificant
plotdiff(emepi) # none
plotdiff(emsoot.mold) # none

# plot the relative effect size for each of the significant taxa
# significant taxa by farming system
# ratio of coefficients:
# CT/Intercept or Org3/Intercept
# alpha = 0.05
plotdiff<-function(x){
  n<-c(rep(0, length(x$tab)))
  for(i in 1:length(x$tab)){
    n[i]<-x$tab[[i]]$`Pr(>F)`[2]
  }
  par(mar=c(10, 4, 4, 2) + 0.1)
  barplot(t(x$effectDF[n<0.05,]), beside = T, las=2)
}



# beta diversity analysis ####

f.otu<-as.data.frame(t(as.matrix(otu_table(fGA.Ql))))
f.dist<-vegdist(f.otu, method="bray")
adonis(f.dist~Depth+Treatment+C_N_ratio+pH+C_percent+No3_ugPerg+Nh4_ugPerg+Clay_percent+Silt_percent, data=f.meta, permutations = 999, method="bray")
adonis(f.dist~C_N_ratio+pH+C_percent+No3_ugPerg+Nh4_ugPerg+Clay_percent+Silt_percent, data=f.meta, permutations = 999, method="bray")

plot_ordination(fGA.Ql,plot(ordinate(fGA.Ql, method="DCA", distance = "bray")), type="samples", color = "Treatment") + facet_wrap(~Depth) + theme_bw()
plot_ordination(fGA.Ql,plot(ordinate(fGA.Ql, method="DCA", distance = "bray")), type="taxa", color = "Phylum.x") + facet_wrap(~primary_lifestyle) + theme_bw()

# network analysis ####

get.no4<-function(v, d2, d1){
  tax1<-as.numeric(as.character(v))
  o<-c(rep(NA,length(d1)))
  for(i in 1:length(d1)){
    tax2<-NULL
    fit<-NULL
    tab<-NULL
    tax2<-as.numeric(as.character(d1[,i]))
    effort<-as.numeric(as.character(d2$SeqDepth))
    QPCR<-as.numeric(as.character(d2$Fun_QPCR))
    treatment<-as.factor(d2$Treatment)
    depth<-as.factor(d2$Depth)
    pH<-as.numeric(as.character(d2$pH))
    Cpercent<-as.numeric(as.character(d2$C_percent))
    Ammonia<-as.numeric(as.character(d2$Nh4_ugPerg))
    Nitrate<-as.numeric(as.character(d2$No3_ugPerg))
    Clay<-as.numeric(as.character(d2$Clay_percent))
    Silt<-as.numeric(as.character(d2$Silt_percent))
    CN_ratio<-as.numeric(as.character(d2$C_N_ratio))
    N_percent<-as.numeric(as.character(d2$N_percent))
    Sand<-as.numeric(as.character(d2$Sand_percent))
    B.Density<-as.numeric(as.character(d2$B.Density_gcm3))
    fit <- glm(tax1 ~ effort*QPCR+treatment*depth+pH+Cpercent*N_percent+Ammonia+Nitrate+Clay+Silt+Sand+B.Density+tax2*treatment,family="gaussian")
    tab<-cbind(summary(aov(fit))[[1]],"varExplained"=summary(aov(fit))[[1]]$`Sum Sq`/sum(summary(aov(fit))[[1]]$`Sum Sq`)*100)
    o[i]<-tab$varExplained[15]*(fit$coefficients[19]/abs(fit$coefficients[19]))/100
    
  }
  o
}
sppInt4<-function(d){
  require(foreach)
  require(doParallel)
  require(phyloseq)
  # prepare data
  d1<-as.data.frame(t(as.matrix(otu_table(d)))) # dim =
  d2<-as.data.frame(as.matrix(sample_data(d)))# dataframes must be samples as rows
  if(!identical(rownames(d1), rownames(d2))){stop("dataframe orientation does not match")}
  
  out<-do.call(cbind, mclapply(d1, get.no4, d2, d1, mc.preschedule = T, mc.cores=10, mc.cleanup = T))
  rownames(out)<-colnames(out)
  tm<-as.data.frame(as.matrix(tax_table(d)))
  tmOrder<-rownames(tm) #tm$Class, tm$Order, 
  if(!identical(colnames(d1), as.character(tmOrder))){stop("taxa names not match")} # sanity check
  tmOrder1<-tmOrder[order(tm$Phylum,tm$Class,tm$Order,tm$Family,tm$Genus,tm$Species)]
  out<-out[as.character(tmOrder1), as.character(tmOrder1)]
  out
  out
  
  
}

test<-sppInt4(fGA.Qlf)
test[is.na(test)]<-0

corrplot::corrplot(test, method="color", 
                   col=colorRampPalette(c("red", "white", "blue"))(200), 
                   order="original",
                   #hclust.method = "average",
                   tl.pos="n",
                   title="test")
max(test)
hist(test)
c<-test
c1<-c
c2<-c
c1[-.15<c1]<-0
c2[c2<0.15]<-0

#*(test>0.1 + (-0.1 > test))
n4<-graph_from_incidence_matrix(c2, directed=T, mode="out", weighted=T)
n4<-delete.vertices(simplify(n4), degree(n4)==0)
E(n4)$color <- ifelse(E(n4)$weight > 0,'navy','maroon')
n4<-plotIgraphNet(n4, fGA.Qlf, main="total network")
runplots1(fGA.Qlf, c)
d<-degree(n4, mode="in", loops=F)
degree_distribution(n4, mode="in", loops=F)
barplot(degree_distribution(n4, mode="in", loops=F),names=as.character(c(0:(length(degree_distribution(n4, mode="in", loops=F))-1))), main="In degree distribution")
barplot(degree_distribution(n4, mode="out", loops=F), names=as.character(c(0:(length(degree_distribution(n4, mode="out", loops=F))-1))), main="Out degree distribution")
V(n4)$primary_lifestyle[V(n4)$name==c(V(n4)$name[degree(n4, mode="in")==top_n(degree(n4, mode="in"), 15)]]

nb.cols<-length(unique(V(n4)$primary_lifestyle))
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)
coul<-mycolors[as.numeric(as.factor(V(n4)$primary_lifestyle))]

#legend("bottomleft", legend=levels(as.factor(V(updatedn2)$primary_lifestyle)), col = mycolors, bty = "n", text.col=mycolors, horiz = FALSE)


plot(jitter(degree(n4, mode= "in")),jitter(degree(n4, mode = "out")), col = coul, pch=16)
legend("topright", legend=levels(as.factor(V(n4)$primary_lifestyle)), col = mycolors, bty = "n", text.col=mycolors, horiz = F, cex=0.75)