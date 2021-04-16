# FSP Fungal Analysis ####

# load libraries ####

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
library(lawstat)
library(factoextra)

# functions for processing samples: ####

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

# Models the relationship between taxon and treatment factors
# d = ps object
# name = name of column from tax table to use for labeling output
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
  colnames(d1)<-paste(sub("g__", "", d3$Genus),sub("s__", "", d3$Species), sep="_")
  treatment<-as.factor(d2$Treatment)
  effort<-as.numeric(as.character(d2$SeqDepth))
  depth<-as.factor(d2$Depth)
  depth<-factor(depth, levels=c("0_5", "5_10", "10Ap"))#, "Ap30"))
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
  out$residuals<-matrix(ncol=ncol(d1), nrow=length(predict(glm(d1[,1] ~ depth*treatment,family="gaussian"))))#+pH+Cpercent*N_percent+Ammonia+Nitrate+Clay+Silt+Sand+B.Density,family="gaussian"))))
  out$predicted<-matrix(ncol=ncol(d1), nrow=length(predict(glm(d1[,1] ~ depth*treatment,family="gaussian"))))#+pH+Cpercent*N_percent+Ammonia+Nitrate+Clay+Silt+Sand+B.Density,family="gaussian"))))
  out$fit<-list(1:length(ncol(d1)))
  #out$plots<-list(1:length(ncol(d1)))
  colnames(out$predicted)<-colnames(d1)
  for(i in c(1:ncol(d1))){
    fit<-NULL
    fit <- glm(d1[,i] ~ depth*treatment, family="gaussian")#+pH+Cpercent*N_percent+Ammonia+Nitrate+Clay+Silt+Sand+B.Density,family="gaussian")
    
    out$residuals[,i]<-residuals(fit)
    out$predicted[,i]<-predict(fit)
    out$fit[[i]]<-fit
    out$tab[[i]]<-cbind(summary(aov(fit))[[1]],"varExplained"=summary(aov(fit))[[1]]$`Sum Sq`/sum(summary(aov(fit))[[1]]$`Sum Sq`)*100)
    #out$plots[[i]]<-boxplot(d1[,i]~depth + treatment, las=2, xlab=NULL,cex.names=0.1, ylab="Abundance",main=species[i], names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8))
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
  #names(out$plots)<-as.character(d3$Species)
  
  out$effectDF<-matrix(nrow=ncol(d1), ncol=4)
  rownames(out$effectDF)<-paste(sub("g__", "", d3$Genus),sub("s__", "", d3$Species), sep="_")
  colnames(out$effectDF)<-c("CT", "ORG", "Lifestyle", "SciName")
  for(i in 1:nrow(out$effectDF)){
    #out$effectDF[i,1]<-out$fit[[i]]$coefficients["(Intercept)"]
    out$effectDF[i,1]<-out$fit[[i]]$coefficients["treatmentCT"]
    out$effectDF[i,2]<-out$fit[[i]]$coefficients["treatmentOrg3"]
    
    # out$effectDF[i,1]<-out$fit[[i]]$coefficients["treatmentCT"]/out$fit[[i]]$coefficients["(Intercept)"]
    #  out$effectDF[i,2]<-out$fit[[i]]$coefficients["treatmentOrg3"]/out$fit[[i]]$coefficients["(Intercept)"]
  }
  out$effectDF[,3]<-d3$primary_lifestyle
  out$effectDF[,4]<-paste(d3$Genus, d3$Species)
  
  out
}
# Models the relationship between taxon and environment
# d = ps object
# name = name of column from tax table to use for labeling output
# outputs:  $fit = glm object
#           $tab = significance table
#           $residuals = residuals from the model
#           $predicted = model predicted values
#           $predcor = correlation matrix of the predicted values for each spp
#           $spcor = species correlation matrix
taxenvmodel<-function(d){ 
  # input is phyloseq object
  # already normalized
  d1<-as.data.frame(t(as.matrix(otu_table(d)))) # dim =
  d2<-as.data.frame(as.matrix(sample_data(d)))# dataframes must be samples as rows
  d3<-as.data.frame(as.matrix(tax_table(d)))
  if(!identical(rownames(d1), rownames(d2))){stop("dataframe orientation does not match")}
  if(!identical(colnames(d1), rownames(d3))){stop("tax table not in correct order")}
  colnames(d1)<-paste(sub("g__", "", d3$Genus),sub("s__", "", d3$Species), sep="_")
  
  #treatment<-as.factor(d2$Treatment)
  effort<-as.numeric(as.character(d2$SeqDepth))
  #depth<-as.factor(d2$Depth)
  #depth<-factor(depth, levels=c("0_5", "5_10", "10Ap"))#, "Ap30"))
  #treatment<-factor(treatment, levels=c("NT", "CT", "Org3"))
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
  out$residuals<-matrix(ncol=ncol(d1), nrow=length(predict(glm(d1[,1] ~ pH+Cpercent*N_percent+Ammonia+Nitrate+Clay+Silt+Sand+B.Density,family="gaussian"))))
  out$predicted<-matrix(ncol=ncol(d1), nrow=length(predict(glm(d1[,1] ~ pH+Cpercent*N_percent+Ammonia+Nitrate+Clay+Silt+Sand+B.Density,family="gaussian"))))
  out$fit<-list(1:length(ncol(d1)))
  #out$plots<-list(1:length(ncol(d1)))
  colnames(out$predicted)<-colnames(d1)
  for(i in c(1:ncol(d1))){
    fit<-NULL
    fit <- glm(d1[,i] ~ pH+Cpercent*N_percent+Ammonia+Nitrate+Clay+Silt+Sand+B.Density,family="gaussian")
    
    out$residuals[,i]<-residuals(fit)
    out$predicted[,i]<-predict(fit)
    out$fit[[i]]<-fit
    out$tab[[i]]<-cbind(summary(aov(fit))[[1]],"varExplained"=summary(aov(fit))[[1]]$`Sum Sq`/sum(summary(aov(fit))[[1]]$`Sum Sq`)*100)
    #out$plots[[i]]<-boxplot(d1[,i]~depth + treatment, las=2, xlab=NULL,cex.names=0.1, ylab="Abundance",main=species[i], names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","NT Ap-30cm","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","CT Ap-30cm","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap","Org Ap-30cm"), par(cex.axis=0.8))
  }
  #names(out$tab)<-paste(sub("g__", "", d3$Genus),sub("s__", "", d3$Species), sep="_")
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
  #names(out$plots)<-as.character(d3$Species)
  
  out$effectDF<-matrix(nrow=ncol(d1), ncol=4)
  rownames(out$effectDF)<-paste(sub("g__", "", d3$Genus),sub("s__", "", d3$Species), sep="_")
  colnames(out$effectDF)<-c("CT", "ORG", "Lifestyle", "SciName")
  for(i in 1:nrow(out$effectDF)){
    #out$effectDF[i,1]<-out$fit[[i]]$coefficients["(Intercept)"]
    out$effectDF[i,1]<-out$fit[[i]]$coefficients["treatmentCT"]
    out$effectDF[i,2]<-out$fit[[i]]$coefficients["treatmentOrg3"]
    
    # out$effectDF[i,1]<-out$fit[[i]]$coefficients["treatmentCT"]/out$fit[[i]]$coefficients["(Intercept)"]
    #  out$effectDF[i,2]<-out$fit[[i]]$coefficients["treatmentOrg3"]/out$fit[[i]]$coefficients["(Intercept)"]
  }
  out$effectDF[,3]<-d3$primary_lifestyle
  out$effectDF[,4]<-paste(d3$Genus, d3$Species)
  
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
    fit <- glm(tax1 ~ depth+pH+Cpercent+Ammonia+Nitrate+Clay+Silt+Sand+B.Density+tax2,family="gaussian")
    tab<-cbind(summary(aov(fit))[[1]],"varExplained"=summary(aov(fit))[[1]]$`Sum Sq`/sum(summary(aov(fit))[[1]]$`Sum Sq`)*100)
    o[i]<-tab$varExplained[10]*(fit$coefficients[12]/abs(fit$coefficients[12]))/100
    
  }
  o
}
sppInt2<-function(d){
  require(foreach)
  require(doParallel)
  require(phyloseq)
  # prepare data
  d1<-as.data.frame(as.matrix(otu_table(d))) # dim =
  d2<-as.data.frame(as.matrix(sample_data(d)))# dataframes must be samples as rows
  if(!identical(rownames(d1), rownames(d2))){
    d1<-as.data.frame(t(as.matrix(otu_table(d))))
    #stop("dataframe orientation does not match")
    }
  
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
# type = one of: bb, bp, bc, pp, pb, pc, cc, cb, cp
# bb = beneficial - beneficial interaction
# bp = beneficial - pathogen interaction
# bc = beneficial - commensal interaction
# pp = pathogen - pathogen interaction
# pb = pathogen - beneficial interaction
# pc = pathogen = commensal interaction
# cc = commensal - commensal interaction
# cb = commensal - beneficial interaction
# cp = commensal - pathogen interaction
# tp = total pathogen interactions
# tb = total beneficial interaction
# tc = total commensal interaction
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
  if(type=="bp"){
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
  if(type=="pb"){
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
  if(type=="pc"){
    tm<-as.data.frame(as.matrix(tax_table(d)))
    tmOrder<-tm$primary_lifestyle #tm$Class, tm$Order, 
    #print(rownames(d1)[1:5])
    #print(as.character(tmOrder)[1:5])
    #if(!identical(colnames(d1), as.character(tmOrder))){stop("taxa names not match")} # sanity check
    tmOrder1<-tmOrder[order(tm$Phylum,tm$Class,tm$Order,tm$Family,tm$Genus,tm$Species)]
    tmOrder1[is.na(tmOrder1)]<-"Unknown"
    m<-m[tmOrder1=="plant_pathogen", tmOrder1=="dung_saprotroph"|tmOrder1=="litter_saprotroph"|tmOrder1=="nectar/tap_saprotroph"|tmOrder1=="soil_saprotroph"|tmOrder1=="unspecified_saprotroph"|tmOrder1=="wood_saprotroph"]
    out[1]<-sum(m>0.3& m<1)
    out[2]<-sum(-0.3>m)
  }
  if(type=="bb"){
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
  if(type=="bc"){
    tm<-as.data.frame(as.matrix(tax_table(d)))
    tmOrder<-tm$primary_lifestyle #tm$Class, tm$Order, 
    #print(rownames(d1)[1:5])
    #print(as.character(tmOrder)[1:5])
    #if(!identical(colnames(d1), as.character(tmOrder))){stop("taxa names not match")} # sanity check
    tmOrder1<-tmOrder[order(tm$Phylum,tm$Class,tm$Order,tm$Family,tm$Genus,tm$Species)]
    tmOrder1[is.na(tmOrder1)]<-"Unknown"
    m<-m[tmOrder1=="ectomycorrhizal"|tmOrder1=="arbuscular_mycorrhizal"|tmOrder1=="unspecified_symbiotroph"|tmOrder1=="moss_symbiont", tmOrder1=="dung_saprotroph"|tmOrder1=="litter_saprotroph"|tmOrder1=="nectar/tap_saprotroph"|tmOrder1=="soil_saprotroph"|tmOrder1=="unspecified_saprotroph"|tmOrder1=="wood_saprotroph"]
    out[1]<-sum(m>0.3& m<1)
    out[2]<-sum(-0.3>m)
  }
  
  if(type=="cc"){
    tm<-as.data.frame(as.matrix(tax_table(d)))
    tmOrder<-tm$primary_lifestyle #tm$Class, tm$Order, 
    #print(rownames(d1)[1:5])
    #print(as.character(tmOrder)[1:5])
    #if(!identical(colnames(d1), as.character(tmOrder))){stop("taxa names not match")} # sanity check
    tmOrder1<-tmOrder[order(tm$Phylum,tm$Class,tm$Order,tm$Family,tm$Genus,tm$Species)]
    tmOrder1[is.na(tmOrder1)]<-"Unknown"
    m<-m[tmOrder1=="dung_saprotroph"|tmOrder1=="litter_saprotroph"|tmOrder1=="nectar/tap_saprotroph"|tmOrder1=="soil_saprotroph"|tmOrder1=="unspecified_saprotroph"|tmOrder1=="wood_saprotroph", tmOrder1=="dung_saprotroph"|tmOrder1=="litter_saprotroph"|tmOrder1=="nectar/tap_saprotroph"|tmOrder1=="soil_saprotroph"|tmOrder1=="unspecified_saprotroph"|tmOrder1=="wood_saprotroph"]
    out[1]<-sum(m>0.3& m<1)
    out[2]<-sum(-0.3>m)
  }
  if(type=="cb"){
    tm<-as.data.frame(as.matrix(tax_table(d)))
    tmOrder<-tm$primary_lifestyle #tm$Class, tm$Order, 
    #print(rownames(d1)[1:5])
    #print(as.character(tmOrder)[1:5])
    #if(!identical(colnames(d1), as.character(tmOrder))){stop("taxa names not match")} # sanity check
    tmOrder1<-tmOrder[order(tm$Phylum,tm$Class,tm$Order,tm$Family,tm$Genus,tm$Species)]
    tmOrder1[is.na(tmOrder1)]<-"Unknown"
    m<-m[tmOrder1=="dung_saprotroph"|tmOrder1=="litter_saprotroph"|tmOrder1=="nectar/tap_saprotroph"|tmOrder1=="soil_saprotroph"|tmOrder1=="unspecified_saprotroph"|tmOrder1=="wood_saprotroph", tmOrder1=="ectomycorrhizal"|tmOrder1=="arbuscular_mycorrhizal"|tmOrder1=="unspecified_symbiotroph"|tmOrder1=="moss_symbiont"]
    out[1]<-sum(m>0.3& m<1)
    out[2]<-sum(-0.3>m)
  }
  if(type=="cp"){
    tm<-as.data.frame(as.matrix(tax_table(d)))
    tmOrder<-tm$primary_lifestyle #tm$Class, tm$Order, 
    #print(rownames(d1)[1:5])
    #print(as.character(tmOrder)[1:5])
    #if(!identical(colnames(d1), as.character(tmOrder))){stop("taxa names not match")} # sanity check
    tmOrder1<-tmOrder[order(tm$Phylum,tm$Class,tm$Order,tm$Family,tm$Genus,tm$Species)]
    tmOrder1[is.na(tmOrder1)]<-"Unknown"
    m<-m[tmOrder1=="dung_saprotroph"|tmOrder1=="litter_saprotroph"|tmOrder1=="nectar/tap_saprotroph"|tmOrder1=="soil_saprotroph"|tmOrder1=="unspecified_saprotroph"|tmOrder1=="wood_saprotroph", tmOrder1=="plant_pathogen"]
    out[1]<-sum(m>0.3& m<1)
    out[2]<-sum(-0.3>m)
  }
  if(type=="tb"){
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
  if(type=="tc"){
    tm<-as.data.frame(as.matrix(tax_table(d)))
    tmOrder<-tm$primary_lifestyle #tm$Class, tm$Order, 
    #print(rownames(d1)[1:5])
    #print(as.character(tmOrder)[1:5])
    #if(!identical(colnames(d1), as.character(tmOrder))){stop("taxa names not match")} # sanity check
    tmOrder1<-tmOrder[order(tm$Phylum,tm$Class,tm$Order,tm$Family,tm$Genus,tm$Species)]
    tmOrder1[is.na(tmOrder1)]<-"Unknown"
    m<-m[tmOrder1=="dung_saprotroph"|tmOrder1=="litter_saprotroph"|tmOrder1=="nectar/tap_saprotroph"|tmOrder1=="soil_saprotroph"|tmOrder1=="unspecified_saprotroph"|tmOrder1=="wood_saprotroph",]
    out[1]<-sum(m>0.3& m<1)
    out[2]<-sum(-0.3>m)
  }
  out
}

interact.summary<-function(m, d){
  out<-matrix(nrow=2, ncol=12)
  rownames(out)<-c("positive", "negative")
  colnames(out)<-c("pathogen-pathogen", "pathogen-beneficial", "pathogen-commensal", "beneficial-beneficial", "beneficial-pathogen", "beneficial-commensal", "commensal-commensal", "commensal-beneficial", "commensal-pathogen", "total beneficial", "total pathogen", "total commensal")
  out[,1]<-getinteraction(m,d,"pp")
  out[,2]<-getinteraction(m,d,"pb")
  out[,3]<-getinteraction(m,d,"pc")
  out[,4]<-getinteraction(m,d,"bb")
  out[,5]<-getinteraction(m,d,"bp")
  out[,6]<-getinteraction(m,d,"bc")
  out[,7]<-getinteraction(m,d,"cc")
  out[,8]<-getinteraction(m,d,"cb")
  out[,9]<-getinteraction(m,d,"cp")
  out[,10]<-getinteraction(m,d,"tb")
  out[,11]<-getinteraction(m,d,"tp")
  out[,12]<-getinteraction(m,d,"tc")
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

maineffects<-function(ps){
  f.meta<-as.data.frame(as.matrix(sample_data(ps)))
  f.meta$SeqDepth<-as.numeric(as.character(f.meta$SeqDepth))
  f.meta$Depth<-factor(f.meta$Depth, levels=c("0_5", "5_10", "10Ap"))#, "Ap30"))
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
  
  print(summary(aov(glm(unlist(sample_sums(ps))~f.meta$Depth*f.meta$Treatment+f.meta$pH+f.meta$C_percent*f.meta$N_percent+f.meta$No3_ugPerg+f.meta$Nh4_ugPerg+f.meta$Clay_percent+f.meta$Sand_percent+f.meta$Silt_percent+f.meta$B.Density_gcm3))))
  
  boxplot(unlist(sample_sums(ps))~f.meta$Depth*f.meta$Treatment,las=2, xlab=NULL,cex.names=0.1)
  o<-glm(unlist(sample_sums(ps))~f.meta$Depth*f.meta$Treatment+f.meta$pH+f.meta$C_percent*f.meta$N_percent+f.meta$No3_ugPerg+f.meta$Nh4_ugPerg+f.meta$Clay_percent+f.meta$Sand_percent+f.meta$Silt_percent+f.meta$B.Density_gcm3)
  o
}
calc.maineffects<-function(list){
  o<-base::lapply(list, taxmodel)
  o
}
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
  #par(mar=c(10, 4, 4, 2) + 0.1)
  bp<-as.data.frame(x$effectDF[n<0.05,])
  CT<-as.numeric(as.character(bp$CT))
  ORG<-as.numeric(as.character(bp$ORG))
  names<-rownames(bp)
  bp<-data.frame(CT, ORG)
  rownames(bp)<-names
  par(mar=c(10,4,4,1)+.1)
  barplot(t(bp), beside = T, las=2, density=30, angle=50, legend=T)
  #legend("topright", legend=colnames(x$effectDF), fill=grey.colors(2))
  abline(0,0)
}

# plot 2.1 ####
# adding coloring to the border of the species plot by lifestyle 
# select top 8 effects
plotdiff2<-function(x){
  n<-c(rep(0, length(x$tab)))
  for(i in 1:length(x$tab)){
    n[i]<-x$tab[[i]]$`Pr(>F)`[2]
  }
  
  bp<-as.data.frame(x$effectDF[n<0.05,])
  CT<-as.numeric(as.character(bp$CT))
  ORG<-as.numeric(as.character(bp$ORG))
  names<-bp$Lifestyle
  bp<-data.frame(CT, ORG)
  rownames(bp)<-names
  nb.cols<-length(unique(names))
  mycolors <- colorRampPalette(brewer.pal(4, "PuOr"))(nb.cols)
  coul<-mycolors[as.factor(names)]
  #par(mar=c(10, 4, 4, 2) + 0.1)
  par(mar=c(8,4,4,1)+.1)
  barplot(t(as.matrix(bp)), beside = T, las=2, density=30, angle=50, legend=T) #border=rep(coul, each=2))
  #legend("topright", legend=unique(bp$Lifestyle), fill=coul)
  #legend("topright", legend=colnames(bp[,c(1,2)]), fill=grey.colors(2), border=coul)
  abline(0,0)
}

plotdiff3<-function(x){
  n<-c(rep(0, length(x$tab)))
  for(i in 1:length(x$tab)){
    n[i]<-x$tab[[i]]$`Pr(>F)`[2]
  }
  
  bp<-as.data.frame(x$effectDF[n<0.05,])
  #rownames(bp)<-bp$SciNames
  bp$CT<-as.numeric(as.character(bp$CT))
  bp$ORG<-as.numeric(as.character(bp$ORG))
  bp$Lifestyle<-as.character(bp$Lifestyle)
  bp$Lifestyle[bp$Lifestyle=="soil_saprotroph"|bp$Lifestyle=="wood_saprotroph"|bp$Lifestyle=="litter_saprotroph"|bp$Lifestyle=="dung_saprotroph"|bp$Lifestyle=="unspecified_saprotroph"|bp$Lifestyle=="nectar/tap_saprotroph"|bp$Lifestyle=="pollen_saprotroph"]<-"Saprotroph"
  bp$Lifestyle[bp$Lifestyle=="plant_pathogen"]<-"Plant Pathogen"
  bp$Lifestyle[bp$Lifestyle=="animal_parasite"]<-"Animal Parasite"
  bp$Lifestyle[bp$Lifestyle=="arbuscular_mycorrhizal"|bp$Lifestyle=="ectomycorrhizal"]<-"Mutualist"
  bp$Lifestyle[bp$Lifestyle=="algal_parasite"|bp$Lifestyle=="unspecified_pathotroph"|bp$Lifestyle=="foliar_endophyte"|bp$Lifestyle=="mycoparasite"]<-"Other"
  #bp<-data.frame(CT, ORG, names)
  #rownames(bp)<-names
  nb.cols<-length(unique(bp$Lifestyle))
  mycolors <- colorRampPalette(brewer.pal(4, "PuOr"))(nb.cols)
  coul<-mycolors[as.factor(bp$Lifestyle)]
  #par(mar=c(10, 4, 4, 2) + 0.1)
  par(mar=c(8,4,4,1)+.1)
  barplot(t(as.matrix(bp[,c(1,2)])), beside = T, las=2, density=30, angle=50, legend=T, border=rep(coul, each=2))
  legend("topright", legend=unique(bp$Lifestyle), fill=unique(coul))
  #legend("topright", legend=colnames(bp[,c(1,2)]), fill=grey.colors(2), border=coul)
  abline(0,0)
}
# import data ####

# import fungal-bacterial dataset
GAD.bac<-readRDS("/Users/dietrich/Documents/GitHub/Anchoring/Data/GAD/GADbac2020.rds")
GAD.Fun<-readRDS("/Users/dietrich/Documents/GitHub/Anchoring/Data/GAD/GADfun2020.rds")

# data preparation ####
# calculate % archaea in each sample
arc<-subset_taxa(GAD.bac, Kingdom=="Archaea")
bac<-subset_taxa(GAD.bac, Kingdom=="Bacteria")
ratio<-sample_sums(bac)/(sample_sums(arc) + sample_sums(bac))

# Using the Bacterial metadata; different sample order
bmeta<-as.data.frame(as.matrix(sample_data(GAD.bac)))
identical(rownames(bmeta), names(ratio)) # make sure the datasets are organized the same way
bmeta$bacteria<-as.numeric(as.character(bmeta$Bac_QPCR))*ratio # take ratio of bacteria from 16s
sample_data(GAD.bac)<-sample_data(bmeta) # return updated metadata to bacterial dataset
GAD.bac<-subset_samples(GAD.bac, Depth=="0_5"|Depth=="5_10"|Depth=="10Ap")#|Depth=="Ap30") # subset the dataset to depths with enough replication
GAD.Fun<-subset_samples(GAD.Fun, Depth=="0_5"|Depth=="5_10"|Depth=="10Ap")
# format data for visualization and modeling
b.meta<-as.data.frame(as.matrix(sample_data(GAD.bac))) # format data to insert
b.meta$Fun_QPCR<-as.numeric(as.character(b.meta$Fun_QPCR)) 
b.meta$Bac_QPCR<-as.numeric(as.character(b.meta$Bac_QPCR))
b.meta$bacteria<-as.numeric(as.character(b.meta$bacteria))
b.meta$Depth<-factor(b.meta$Depth, levels=c("0_5", "5_10", "10Ap"))#, "Ap30")) # put factor levels in the right order
b.meta$Treatment<-factor(b.meta$Treatment, levels=c("NT", "CT", "Org3")) # put factor levels in the right order
b.meta$FB_Ratio<-b.meta$Fun_QPCR/b.meta$bacteria # calculate the fungal:bacteria ratio
b.meta2<-b.meta[b.meta$Description!="Org30_54"&b.meta$Description!="Org35_104"&b.meta$Description!="Org310Ap4"&b.meta$Description!="Org3Ap304",] # remove outlier for ratio analysis
b.meta3<-b.meta[b.meta$Description=="Org30_54"|b.meta$Description=="Org35_104"|b.meta$Description=="Org310Ap4"|b.meta$Description=="Org3Ap304",] # look at the data for only the outlier core

# correlation of Bactera vs fungi across systems ####
# total:
cor(b.meta2$Bac_QPCR,b.meta2$Fun_QPCR)
# NT: r=0.93
cor(b.meta2$Bac_QPCR[b.meta2$Treatment=="NT"],b.meta2$Fun_QPCR[b.meta2$Treatment=="NT"])
# CT: r=0.697
cor(b.meta2$Bac_QPCR[b.meta2$Treatment=="CT"],b.meta2$Fun_QPCR[b.meta2$Treatment=="CT"])
# Org3: r=0.733
cor(b.meta2$Bac_QPCR[b.meta2$Treatment=="Org3"],b.meta2$Fun_QPCR[b.meta2$Treatment=="Org3"])

# correlation of Carbon and Nitrogen across systems
# NT: r=0.97
cor(as.numeric(b.meta2$C_percent[b.meta2$Treatment=="NT"]),as.numeric(b.meta2$N_percent[b.meta2$Treatment=="NT"]))
# CT: r=0.98
cor(as.numeric(b.meta2$C_percent[b.meta2$Treatment=="CT"]),as.numeric(b.meta2$N_percent[b.meta2$Treatment=="CT"]))
# Org3: r=0.83 or 0.778 if using b.meta2 (filter out the outliers)
cor(as.numeric(b.meta2$C_percent[b.meta2$Treatment=="Org3"]),as.numeric(b.meta2$N_percent[b.meta2$Treatment=="Org3"]))
# correlation of Carbon:Nitrogen and F:B across systems
# NT: r=-0.033
cor(as.numeric(b.meta2$FB_Ratio[b.meta2$Treatment=="NT"]),as.numeric(b.meta2$C_N_ratio[b.meta2$Treatment=="NT"]))
# CT: r=0.259
cor(as.numeric(b.meta2$FB_Ratio[b.meta2$Treatment=="CT"]),as.numeric(b.meta2$C_N_ratio[b.meta2$Treatment=="CT"]))
# Org3: r=0.172 or -0.0056 if using b.meta2
cor(as.numeric(b.meta2$FB_Ratio[b.meta2$Treatment=="Org3"]),as.numeric(b.meta2$C_N_ratio[b.meta2$Treatment=="Org3"]))
cor(as.numeric(b.meta2$FB_Ratio),as.numeric(b.meta2$C_N_ratio))

# total abundance ####
boxplot(b.meta$Fun_QPCR~b.meta$Depth+b.meta$Treatment,las=2, xlab=NULL,cex.names=0.1,main="Total Fungi QPCR", ylab="Gene Count", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap"), par(cex.axis=0.8))
stripchart(b.meta$Fun_QPCR~b.meta$Depth+b.meta$Treatment, vertical = TRUE,
           method = "jitter", add = TRUE, pch = 20, col = 'black')
summary(aov(b.meta$Fun_QPCR~b.meta$Depth*b.meta$Treatment))
TukeyHSD(aov(b.meta$Fun_QPCR~b.meta$Depth*b.meta$Treatment))
fungiQfit<-aov(b.meta$Fun_QPCR~b.meta$Depth*b.meta$Treatment)$coefficients
fungiQfit[5]/fungiQfit[1] # % change organic vs NT = 33
fungiQfit[4]/fungiQfit[1] # % CT vs NT = 19

# figure 1.2 ###
boxplot(b.meta$bacteria~b.meta$Depth+b.meta$Treatment,las=2, xlab=NULL,cex.names=0.1,main="Total Bacteria QPCR", ylab="Gene Count", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap"), par(cex.axis=0.8))
stripchart(b.meta$bacteria~b.meta$Depth+b.meta$Treatment, vertical = TRUE,
           method = "jitter", add = TRUE, pch = 20, col = 'black')
summary(aov(b.meta$bacteria~b.meta$Depth*b.meta$Treatment))
TukeyHSD(aov(b.meta$bacteria~b.meta$Depth*b.meta$Treatment))
# figure 1.3 ####
boxplot(b.meta$FB_Ratio~b.meta$Depth+b.meta$Treatment,las=2, xlab=NULL,cex.names=0.1,main="Fungal Bacterial Ratio", ylab="Ratio of Fungal to Bacterial Genes", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap"), par(cex.axis=0.8))
stripchart(b.meta$FB_Ratio~b.meta$Depth+b.meta$Treatment, vertical = TRUE,
           method = "jitter", add = TRUE, pch = 20, col = 'black')
summary(aov(b.meta$FB_Ratio~b.meta$Depth*b.meta$Treatment))

# plots/stats w/o outlier core
boxplot(b.meta2$Bac_QPCR~b.meta2$Depth+b.meta2$Treatment,las=2, xlab=NULL,cex.names=0.1,main="Bacterial Gene Count", ylab="Gene Count", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap"), par(cex.axis=0.8))
stripchart(b.meta2$Bac_QPCR~b.meta2$Depth+b.meta2$Treatment, vertical = TRUE,
           method = "jitter", add = TRUE, pch = 20, col = 'black')
summary(aov(b.meta2$Bac_QPCR~b.meta2$Depth*b.meta2$Treatment))
bacQfit<-aov(b.meta2$Bac_QPCR~b.meta2$Depth*b.meta2$Treatment)$coefficients
bacQfit[5]/bacQfit[1] # % change organic vs NT = -12
bacQfit[4]/bacQfit[1] # % CT vs NT = -11


boxplot(b.meta2$FB_Ratio~b.meta2$Depth+b.meta2$Treatment,las=2, xlab=NULL,cex.names=0.1,main="Fungal Bacterial Ratio", ylab="Gene Count", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap"), par(cex.axis=0.8))
stripchart(b.meta2$FB_Ratio~b.meta2$Depth+b.meta2$Treatment, vertical = TRUE,
           method = "jitter", add = TRUE, pch = 20, col = 'black')
summary(aov(b.meta2$FB_Ratio~b.meta2$Depth*b.meta2$Treatment))

ggplot(b.meta2, aes(x=Bac_QPCR, y=Fun_QPCR, color=Treatment, shape=Depth, size=2)) + 
  geom_point() + 
  ggtitle("Fungal vs Bacterial Gene Abundance")+ 
  labs(y="Fungal Gene Abundance", x="Bacterial Gene Abundance")+
  scale_color_grey() + 
  theme_bw()

summary(aov(b.meta2$FB_Ratio~b.meta2$Depth*b.meta2$Treatment))

sd<-with(b.meta2, summary(aov(glm(FB_Ratio~Depth+Treatment+as.numeric(as.character(pH))+as.numeric(as.character(C_percent))+as.numeric(as.character(N_percent))+as.numeric(as.character(C_N_ratio))+as.numeric(as.character(No3_ugPerg))+as.numeric(as.character(Nh4_ugPerg))+as.numeric(as.character(Clay_percent))+as.numeric(as.character(Silt_percent))+as.numeric(as.character(Sand_percent))+as.numeric(as.character(B.Density_gcm3))))))

v<-100*sd[[1]]$`Sum Sq`/sum(sd[[1]]$`Sum Sq`)
names(v)<-c("Depth", "Farming System", "pH", "C Percent", "N Percent", "CN Ratio", "Nitrate", "Ammonium", "Clay Percent", "Silt Percent", "Sand Percent", "Bulk Density", "Residuals")
par(mar=c(8,4,4,1)+.1)
barplot(v, main="Predictors of F:B", ylab="Percent Variation Explained", las=2)

# model fungal vs bacterial depth dependence
bdepthfit<-glm(bacteria~Depth*Treatment,data = b.meta)
fdepthfit<-glm(Fun_QPCR~Depth*Treatment,data = b.meta)
fdepthfit2<-glm(Fun_QPCR~Depth*Treatment,data = b.meta2)
rdepthfit<-glm(FB_Ratio~Depth*Treatment,data = b.meta)
#p1<-effect_plot(bdepthfit, pred=Treatment, plot.points = T)
#p2<-effect_plot(fdepthfit, pred=Treatment)

bdepthfit$coefficients
fdepthfit$coefficients # 33.46 %
fdepthfit2$coefficients # 32.87 %

# filtering and scaling ####
sample_data(GAD.Fun)$SeqDepth<-sample_sums(GAD.Fun)
sample_data(GAD.Fun)$SampleDepth<-sample_data(GAD.Fun)$SeqDepth/sample_data(GAD.Fun)$Fun_QPCR
tt1<-as.data.frame(as.matrix(tax_table(GAD.Fun)))
sum(is.na(tt1$Species))/nrow(tt1) # 75 % unknown to species
sum(is.na(tt1$Genus))/nrow(tt1) # 62 % unknown to genus
sum(is.na(tt1$Family))/nrow(tt1) # 56 % unknown to family
sum(is.na(tt1$Order))/nrow(tt1) # 50 % unknown to order
sum(is.na(tt1$Class))/nrow(tt1) # 46 % unknown to class

# if I keep unknown spp by pasting genus; aggregate to spp, I lose 6% compared to aggregating to family and 12%$ compared to aggregating to order (16% aggregated to class level). Class and above are more likely to be random amplification targets.
#GAD.Fun.nf<-GAD.Fun
#GAD.Fun.nf<-subset_samples(GAD.Fun.nf, Depth=="0_5"|Depth=="5_10"|Depth=="10Ap"|Depth=="Ap30") # dataset that is not taxa filtered (compare overarching patterns vs completely subset dataset)

# quantitative scaling step:
fGA.Q<-GAD.QScale(GAD.Fun,type="F")
# filtering steps:
# These steps prepare the dataset for annotation, which requires taxonomy at genus level

fGA.Q<-subset_taxa(fGA.Q, !is.na(Genus)) # remove anything not ID'd to Genus
fGA.Qa<-fGA.Q # use Qa for alpha diversity (set aside for alpha div analysis; next step is to aggregate taxa to spp level + unknowns grouped together)
fGA.Qnf<-subset_taxa(GAD.Fun, is.na(Genus)) # diversity of removed organisms

# now give all unknown spp a same name:
tt1<-as.data.frame(as.matrix(tax_table(fGA.Q)))
tt1$Species[is.na(tt1$Species)]<-paste0(tt1$Genus[is.na(tt1$Species)], "Undefined", sep="")
tax_table(fGA.Q)<-tax_table(as.matrix(tt1))
# now aggregate to species level
fGA.Q<-tax_glom(fGA.Q, taxrank="Species") # automatically removes taxa not assigned to spp level! shouldn't be a problem if previous step worked (unknowns should be "genusUndefined")

# annotate fungi by fungal traits database ####
fungalTraits<-read.csv("/Users/dietrich/Documents/GitHub/Plant-Health-Project/Analysis/Data/FungalTraits2021.csv")
#View(fungalTraits)
fungalTraits<-as.data.frame(fungalTraits) # convert reference data to dataframe 
fungalTraits$GENUS<-paste0("g__", fungalTraits$GENUS, sep="") # Add "g__" to reference database
tt<-as.data.frame(as.matrix(tax_table(fGA.Q))) # extracting the tax table from phyloseq = fGA.Q
tt2<-left_join(tt, fungalTraits, by=c("Genus"="GENUS")) # merge the dataframes by genus column
identical(tt$Family, tt2$Family.x) # sanity check: make sure family match after merging

# stop here and make sure output is TRUE !!

rownames(tt2)<-rownames(tt) # adding rownames 
tt2<-tt2[,-c(8:13,16,21,25,29:31)] # remove columns not being used in the analysis
tax_table(fGA.Q)<-tax_table(as.matrix(tt2)) # merge tax table with annotations back into phyloseq object

# edaphic factors ####
f.meta<-as.data.frame(as.matrix(sample_data(fGA.Q)))
f.meta$Fun_QPCR<-as.numeric(as.character(f.meta$Fun_QPCR))
f.meta$Bac_QPCR<-as.numeric(as.character(f.meta$Bac_QPCR))
f.meta$FB_Ratio<-f.meta$Fun_QPCR/f.meta$Bac_QPCR
f.meta$Depth<-factor(f.meta$Depth, levels=c("0_5", "5_10", "10Ap"))#, "Ap30"))
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

# env PCA
#library(factoextra)
dat<-f.meta[,c(8:11,13:18)]
group<-f.meta$Treatment
group2<-f.meta$Depth

res.pca <- prcomp(dat, scale = TRUE)

fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)



fviz_pca_ind(res.pca,
             col.ind = group, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Farming System",
             label="none",
             pointsize = 3,
             ellipse.level=0.99999,
             repel = TRUE
)

fviz_pca_ind(res.pca,
             col.ind = group2, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07", "#696969"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Depth",
             label="none",
             pointsize = 3,
             ellipse.level=0.99999,
             repel = TRUE
)

fviz_pca_biplot(res.pca, label="var",palette = c("#00AFBB", "#E7B800", "#FC4E07", "#696969"), habillage=f.meta$Treatment,addEllipses=TRUE,pointsize = 3,ellipse.level=0.99,ellipse.type = "confidence")

fviz_pca_biplot(res.pca, label="var",palette = c("#00AFBB", "#E7B800", "#FC4E07", "#696969"), habillage=f.meta$Depth,addEllipses=TRUE,pointsize = 3,ellipse.level=0.99,ellipse.type = "confidence")


# alpha diversity ####
# model total number of OTU's captured as a function of sequencing depth
# use several different data sources;
# unfiltered data and filtered data
# Check consistency
# unfiltered believed to be more accurate because it includes data for taxa that only occure in one sample
# e.g. less likely to be biased from filtering parameters, which includes databases

# alpha diversity of whole dataset, no filtering
a.n1ffit<-glm(unlist(estimate_richness(GAD.Fun, measures="Observed"))~sample_data(GAD.Fun)$SeqDepth)
#a.bfit<-lm(unlist(estimate_richness(bGA.Ql, measures="Observed"))~sample_sums(bGA.Ql))
summary(aov(a.n1ffit))  # by seq depth = 17.8% of variance explained at spp level
plot(fitted(a.n1ffit), resid(a.n1ffit)) # check residual plot
resids<-a.n1ffit$residuals
alphadiveffect<-summary(aov(glm(resids~f.meta$Depth*f.meta$Treatment))) # no effect of farming system on alpha diversity
alphadiveffect

#boxplot
boxplot(resids~f.meta$Depth+f.meta$Treatment, las=2, xlab=NULL,cex.names=0.1,main="Fungal Unfiltered Alpha Diversity", ylab="Sampling Effort Residual", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap"), par(cex.axis=0.8))
abline(h = 0, col = 'black', lty=c(2)) 

# partitioning of residual variance:
# 17.8, 24.1 each comes from variance explained by the model
#divVar<-data.frame("SeqDepth"=(100-17.8)*alphadiveffect[[1]]$`Sum Sq`/sum(alphadiveffect[[1]]$`Sum Sq`), "SampleEffort"=(100-24.1)*alphadiveffect2[[1]]$`Sum Sq`/sum(alphadiveffect2[[1]]$`Sum Sq`))
#rownames(divVar)<-c("Depth", "Farming System", "Interaction", "Residual")
#divVar


# alpha diversity of filtered taxa
identical(rownames(f.meta), rownames(estimate_richness(fGA.Qa, measures="Observed"))) # sanity check
a.ffit<-glm(unlist(estimate_richness(fGA.Qa, measures="Observed"))~sample_data(fGA.Qa)$SeqDepth)
nresids1<-a.ffit$residuals # extract residuals
n1alphadiveffect<-summary(aov(glm(nresids1~f.meta$Depth*f.meta$Treatment))) # anova of residuals against depth and treatment
n1alphadiveffect
plot(fitted(a.ffit), resid(a.ffit)) # plot residuals
abline(0,0)
#boxplot
boxplot(nresids1~f.meta$Depth+f.meta$Treatment, las=2, xlab=NULL,cex.names=0.1,main="Fungal filtered Alpha Diversity", ylab="Sequencing Depth Residual", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap"), par(cex.axis=0.8))
abline(h = 0, col = 'black', lty=c(2)) 

# abundance of functional groups & beta diversity ####
# filter to interpretable taxa only for the species-wise analysis
fGA.Qlf<-filter_taxa(fGA.Q, function(x) sum(x > 3) > 10, TRUE) # filter to interpretable taxa

#filter_taxa(GAD.fun, function(x) sum(x)>0, TRUE)
f.otu1<-as.data.frame(as.matrix(otu_table(filter_taxa(GAD.Fun, function(x) sum(x)>0, TRUE))))
f.meta2<-as.data.frame(as.matrix(sample_data(GAD.Fun)))
f.meta2$Fun_QPCR<-as.numeric(as.character(f.meta2$Fun_QPCR))
f.meta2$Bac_QPCR<-as.numeric(as.character(f.meta2$Bac_QPCR))
f.meta2$FB_Ratio<-f.meta2$Fun_QPCR/f.meta2$Bac_QPCR
f.meta2$Depth<-factor(f.meta2$Depth, levels=c("0_5", "5_10", "10Ap"))#, "Ap30"))
f.meta2$Treatment<-factor(f.meta2$Treatment, levels=c("NT", "CT", "Org3"))
f.meta2$pH<-as.numeric(as.character(f.meta2$pH))
f.meta2$C_percent<-as.numeric(as.character(f.meta2$C_percent))
f.meta2$No3_ugPerg<-as.numeric(as.character(f.meta2$No3_ugPerg))
f.meta2$Nh4_ugPerg<-as.numeric(as.character(f.meta2$Nh4_ugPerg))
f.meta2$C_N_ratio<-as.numeric(as.character(f.meta2$C_N_ratio))
f.meta2$N_percent<-as.numeric(as.character(f.meta2$N_percent))
f.meta2$Clay_percent<-as.numeric(as.character(f.meta2$Clay_percent))
f.meta2$Silt_percent<-as.numeric(as.character(f.meta2$Silt_percent))
f.meta2$Sand_percent<-as.numeric(as.character(f.meta2$Sand_percent))
f.meta2$B.Density_gcm3<-as.numeric(as.character(f.meta2$B.Density_gcm3))
f.meta2$SeqDepth<-as.numeric(as.character(f.meta2$SeqDepth))
f.meta2$SampleDepth<-as.numeric(as.character(f.meta2$SampleDepth))

# unfiltered dataset beta diversity:
f.dist<-vegdist(f.otu1, method="bray")
adonis(f.dist~Depth+Treatment+C_N_ratio+pH+C_percent+No3_ugPerg+Nh4_ugPerg+Clay_percent+Silt_percent, data=f.meta2, permutations = 999, method="bray")
adonis(f.dist~C_N_ratio+pH+C_percent+No3_ugPerg+Nh4_ugPerg+Clay_percent+Silt_percent, data=f.meta2, permutations = 999, method="bray")

# filtered dataset beta diversity:
f.otu2<-as.data.frame(t(as.matrix(otu_table(fGA.Qlf))))
f.dist2<-vegdist(f.otu2, method="bray")
adonis(f.dist2~Depth+Treatment+C_N_ratio+pH+C_percent+No3_ugPerg+Nh4_ugPerg+Clay_percent+Silt_percent, data=f.meta, permutations = 999, method="bray")
adonis(f.dist2~C_N_ratio+pH+C_percent+No3_ugPerg+Nh4_ugPerg+Clay_percent+Silt_percent, data=f.meta, permutations = 999, method="bray")

fGA.Qlf2<-fGA.Q # prepare data for aggregated to whole functional group
tax<-as.data.frame(as.matrix(tax_table(fGA.Qlf2))) # 
tax<-tax[,c(8,1:7,9:19)] # 
tax_table(fGA.Qlf2)<-tax_table(as.matrix(tax)) #

unique(tax$primary_lifestyle)
# for main effects include total taxa
mp.lifestyle<-tax_glom(fGA.Qlf2, taxrank = "primary_lifestyle")

# prepare data for the individual functional groups
total.sap<-subset_taxa(fGA.Qlf2,primary_lifestyle=="soil_saprotroph"|primary_lifestyle=="wood_saprotroph"|primary_lifestyle=="unspecified_saprotroph"|primary_lifestyle=="litter_saprotroph"|primary_lifestyle=="dung_saprotroph"|primary_lifestyle=="nectar/tap_saprotroph"|primary_lifestyle=="pollen_saprotroph")
mp.path<-subset_taxa(fGA.Qlf2, primary_lifestyle=="plant_pathogen")
msoil.sap<-subset_taxa(fGA.Qlf2, primary_lifestyle=="soil_saprotroph")
man.par<-subset_taxa(fGA.Qlf2, primary_lifestyle=="animal_parasite")
mwood.sap<-subset_taxa(fGA.Qlf2, primary_lifestyle=="wood_saprotroph")
#munknown<-subset_taxa(fGA.Qlf2, primary_lifestyle=="NA")
mmyco.par<-subset_taxa(fGA.Qlf2, primary_lifestyle=="mycoparasite")
munsp.sap<-subset_taxa(fGA.Qlf2, primary_lifestyle=="unspecified_saprotroph")
mlitter.sap<-subset_taxa(fGA.Qlf2, primary_lifestyle=="litter_saprotroph")
mdung.sap<-subset_taxa(fGA.Qlf2, primary_lifestyle=="dung_saprotroph")
mAMF<-subset_taxa(fGA.Qlf2, primary_lifestyle=="arbuscular_mycorrhizal")
mfol.end<-subset_taxa(fGA.Qlf2, primary_lifestyle=="foliar_endophyte")
mnect.sap<-subset_taxa(fGA.Qlf2, primary_lifestyle=="nectar/tap_saprotroph")
mpollen.sap<-subset_taxa(fGA.Qlf2, primary_lifestyle=="pollen_saprotroph")
#mlich.par<-subset_taxa(fGA.Qlf2, primary_lifestyle=="lichen_parasite")
mECM<-subset_taxa(fGA.Qlf2, primary_lifestyle=="ectomycorrhizal")
munsp.path<-subset_taxa(fGA.Qlf2, primary_lifestyle=="unspecified_pathotroph")
#mepi<-subset_taxa(fGA.Qlf2, primary_lifestyle=="epiphyte")
malg.par<-subset_taxa(fGA.Qlf2, primary_lifestyle=="algal_parasite")
#msoot.mold<-subset_taxa(fGA.Qlf2, primary_lifestyle=="sooty_mold")
#munknown<-subset_taxa(fGA.Qlf2, primary_lifestyle=="unspecified")
#mlichenized<-subset_taxa(fGA.Qlf2, primary_lifestyle=="lichenized")

# prepare data for boxplots
ECM.a<-unlist(sample_sums(mECM))
AMF.a<-unlist(sample_sums(mAMF))
p.path.a<-unlist(sample_sums(mp.path))
soilsap.a<-unlist(sample_sums(msoil.sap))
man.par.a<-unlist(sample_sums(man.par))
total.sap.a<-unlist(sample_sums(total.sap))
wood.sap.a<-unlist(sample_sums(mwood.sap))
litter.sap.a<-unlist(sample_sums(mlitter.sap))
dung.sap.a<-unlist(sample_sums(mdung.sap))

ab.df<-data.frame(ECM.a, 
                  AMF.a, 
                  p.path.a,
                  soilsap.a,
                  man.par.a,
                  total.sap.a,
                  wood.sap.a,
                  litter.sap.a,
                  dung.sap.a,
                  "depth"=as.data.frame(as.matrix(sample_data(fGA.Qlf2)))$Depth,
                  "treatment"= as.data.frame(as.matrix(sample_data(fGA.Qlf2)))$Treatment,
                  "C"= as.data.frame(as.matrix(sample_data(fGA.Qlf2)))$C_percent,
                  "N"= as.data.frame(as.matrix(sample_data(fGA.Qlf2)))$N_percent,
                  "nitrate"= as.data.frame(as.matrix(sample_data(fGA.Qlf2)))$No3_ugPerg,
                  "ammonium"= as.data.frame(as.matrix(sample_data(fGA.Qlf2)))$Nh4_ugPerg,
                  "Clay"= as.data.frame(as.matrix(sample_data(fGA.Qlf2)))$Clay_percent,
                  "Silt"= as.data.frame(as.matrix(sample_data(fGA.Qlf2)))$Silt_percent,
                  "Sand"= as.data.frame(as.matrix(sample_data(fGA.Qlf2)))$Sand_percent,
                  "pH"= as.data.frame(as.matrix(sample_data(fGA.Qlf2)))$pH,
                  "Density"= as.data.frame(as.matrix(sample_data(fGA.Qlf2)))$B.Density_gcm3)
ab.df$depth<-factor(ab.df$depth,levels=c("0_5", "5_10", "10Ap"))
ab.df$treatment<-factor(ab.df$treatment,levels=c("NT", "CT", "Org3"))
ab.df$pH<-as.numeric(as.character(ab.df$pH))
ab.df$C<-as.numeric(as.character(ab.df$C))
ab.df$N<-as.numeric(as.character(ab.df$N))
ab.df$ammonium<-as.numeric(as.character(ab.df$ammonium))
ab.df$nitrate<-as.numeric(as.character(ab.df$nitrate))
ab.df$Clay<-as.numeric(as.character(ab.df$Clay))
ab.df$Silt<-as.numeric(as.character(ab.df$Silt))
ab.df$Sand<-as.numeric(as.character(ab.df$Sand))
ab.df$Density<-as.numeric(as.character(ab.df$Density))

total.sap1<-subset_taxa(fGA.Qlf,primary_lifestyle=="soil_saprotroph"|primary_lifestyle=="wood_saprotroph"|primary_lifestyle=="unspecified_saprotroph"|primary_lifestyle=="litter_saprotroph"|primary_lifestyle=="dung_saprotroph"|primary_lifestyle=="nectar/tap_saprotroph"|primary_lifestyle=="pollen_saprotroph")
mp.path1<-subset_taxa(fGA.Qlf, primary_lifestyle=="plant_pathogen")
msoil.sap1<-subset_taxa(fGA.Qlf, primary_lifestyle=="soil_saprotroph")
man.par1<-subset_taxa(fGA.Qlf, primary_lifestyle=="animal_parasite")
mwood.sap1<-subset_taxa(fGA.Qlf, primary_lifestyle=="wood_saprotroph")
#munknown<-subset_taxa(fGA.Qlf2, primary_lifestyle=="NA")
mmyco.par1<-subset_taxa(fGA.Qlf, primary_lifestyle=="mycoparasite")
munsp.sap1<-subset_taxa(fGA.Qlf, primary_lifestyle=="unspecified_saprotroph")
mlitter.sap1<-subset_taxa(fGA.Qlf, primary_lifestyle=="litter_saprotroph")
mdung.sap1<-subset_taxa(fGA.Qlf, primary_lifestyle=="dung_saprotroph")
mAMF1<-subset_taxa(fGA.Qlf, primary_lifestyle=="arbuscular_mycorrhizal")
mfol.end1<-subset_taxa(fGA.Qlf, primary_lifestyle=="foliar_endophyte")
mnect.sap1<-subset_taxa(fGA.Qlf, primary_lifestyle=="nectar/tap_saprotroph")
mpollen.sap1<-subset_taxa(fGA.Qlf, primary_lifestyle=="pollen_saprotroph")
#mlich.par<-subset_taxa(fGA.Qlf2, primary_lifestyle=="lichen_parasite")
mECM1<-subset_taxa(fGA.Qlf, primary_lifestyle=="ectomycorrhizal")
munsp.path1<-subset_taxa(fGA.Qlf, primary_lifestyle=="unspecified_pathotroph")
#mepi<-subset_taxa(fGA.Qlf2, primary_lifestyle=="epiphyte")
malg.par1<-subset_taxa(fGA.Qlf, primary_lifestyle=="algal_parasite")
#msoot.mold<-subset_taxa(fGA.Qlf2, primary_lifestyle=="sooty_mold")
#munknown<-subset_taxa(fGA.Qlf2, primary_lifestyle=="unspecified")
#mlichenized<-subset_taxa(fGA.Qlf2, primary_lifestyle=="lichenized")

total.sp<-taxmodel(fGA.Qlf2)
total.lifestyle<-taxmodel(mp.lifestyle)
emp.path<-taxmodel(mp.path1)
emsoil.sap<-taxmodel(msoil.sap1)
eman.par<-taxmodel(man.par1)
emwood.sap<-taxmodel(mwood.sap1)
emmyco.par<-taxmodel(mmyco.par1)
emunsp.sap<-taxmodel(munsp.sap1)
emlitter.sap<-taxmodel(mlitter.sap1)
emdung.sap<-taxmodel(mdung.sap1)
emAMF<-taxmodel(mAMF1)
emfol.end<-taxmodel(mfol.end1)
emnect.sap<-taxmodel(mnect.sap1)
empollen.sap<-taxmodel(mpollen.sap1)
#emlich.par<-taxmodel(mlich.par)
emECM<-taxmodel(mECM1)
emunsp.path<-taxmodel(munsp.path1)
emalg.par<-taxmodel(malg.par1)
#emepi<-taxmodel(mepi)
#emsoot.mold<-taxmodel(msoot.mold)

# biplot based on lifestyle
# plots are based on filtered to 0.05 alpha
par(mar=c(8,4,4,1)+.1)
plotdiff2(total.lifestyle)

# analysis of decomposers ####

# soil saprotrophs
plotdiff(emsoil.sap)
with(ab.df,boxplot(soilsap.a~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Soil Saprotroph Abundance", ylab="Abundance", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap"), par(cex.axis=0.8)))
with(ab.df,summary(aov(glm(soilsap.a~depth+treatment+pH+C*N+ammonium+nitrate+Clay+Silt+Sand+Density))))
with(ab.df,TukeyHSD(aov(glm(soilsap.a~depth*treatment))))
with(ab.df,aov(glm(soilsap.a~depth*treatment))$coefficients)# 99.7% org vs NT

# wood saprotrophs
plotdiff(emwood.sap)
with(ab.df,boxplot(wood.sap.a~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Wood Saprotroph Abundance", ylab="Abundance", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap"), par(cex.axis=0.8)))
with(ab.df,summary(aov(glm(wood.sap.a~depth*treatment))))
with(ab.df,TukeyHSD(aov(glm(wood.sap.a~depth*treatment))))
with(ab.df,aov(glm(wood.sap.a~depth*treatment))$coefficients)

# litter saprotrophs
plotdiff(emlitter.sap)
with(ab.df,boxplot(litter.sap.a~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Litter Saprotroph Abundance", ylab="Abundance", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap"), par(cex.axis=0.8)))
with(ab.df,summary(aov(glm(litter.sap.a~depth*treatment))))
with(ab.df,TukeyHSD(aov(glm(litter.sap.a~depth*treatment))))
with(ab.df,aov(glm(litter.sap.a~depth*treatment))$coefficients)
with(ab.df,aov(glm(litter.sap.a~depth*treatment))$coefficients[5]/aov(glm(litter.sap.a~depth*treatment))$coefficients[1])*100
# dung saprotroph
plotdiff(emdung.sap)
with(ab.df,boxplot(dung.sap.a~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Dung Saprotroph Abundance", ylab="Abundance", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap"), par(cex.axis=0.8)))
with(ab.df,summary(aov(glm(dung.sap.a~depth*treatment))))
with(ab.df,TukeyHSD(aov(glm(dung.sap.a~depth*treatment))))
with(ab.df,aov(glm(dung.sap.a~depth*treatment))$coefficients)
with(ab.df,aov(glm(dung.sap.a~depth*treatment))$coefficients[5]/aov(glm(dung.sap.a~depth*treatment))$coefficients[1])*100 # percent increase in Organic vs NT
# analysis of other beneficials ####

# AMF
plotdiff(emAMF)
with(ab.df,boxplot(AMF.a~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="AMF Abundance", ylab="Abundance", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap"), par(cex.axis=0.8)))
with(ab.df,summary(aov(glm(AMF.a~depth*treatment))))
with(ab.df,TukeyHSD(aov(glm(AMF.a~depth*treatment))))
with(ab.df,aov(glm(AMF.a~depth*treatment))$coefficients)

amfGradient<-taxenvmodel(mAMF1)
emAMF$fit$Glomus_aggregatum$coefficients
emAMF$fit$Septoglomus_viscosum$coefficients
emAMF$fit$Rhizophagus_irregularis$coefficients
emAMF$fit$Claroideoglomus_g__ClaroideoglomusUndefined$coefficients

amfGradient$tab$Glomus_aggregatum
amfGradient$fit$Glomus_aggregatum$coefficients

amfGradient$tab$Septoglomus_viscosum
amfGradient$fit$Septoglomus_viscosum$coefficients
amfGradient$tab$Rhizophagus_irregularis
amfGradient$fit$Rhizophagus_irregularis$coefficients
amfGradient$tab$Claroideoglomus_g__ClaroideoglomusUndefined
amfGradient$fit$Claroideoglomus_g__ClaroideoglomusUndefined$coefficients
# ECM
plotdiff(emECM)
with(ab.df,boxplot(ECM.a~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="ECM Abundance", ylab="Abundance", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap"), par(cex.axis=0.8)))
with(ab.df,summary(aov(glm(ECM.a~depth*treatment))))

# Animal Parasites
plotdiff(eman.par)
with(ab.df[-c(30),],boxplot(man.par.a~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Animal Parasite Abundance", ylab="Abundance", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap"), par(cex.axis=0.8)))
with(ab.df[-c(30),],summary(aov(glm(man.par.a~depth*treatment))))

Anparenv<-taxenvmodel(man.par1)
Anparenv$fit$Metarhizium_marquandii
Anparenv$tab$Metarhizium_marquandii
eman.par$fit$Metarhizium_marquandii$coefficients

# Plant Pathogens
plotdiff(emp.path)
with(ab.df,boxplot(p.path.a~depth+treatment,las=2, xlab=NULL,cex.names=0.1, main="Plant Pathogen Abundance", ylab="Abundance", names=c("NT 0-5cm","NT 5-10cm","NT 10cm-Ap","CT 0-5cm","CT 5-10cm","CT 10cm-Ap","Org 0-5cm", "Org 5-10cm","Org 10cm-Ap"), par(cex.axis=0.8)))
with(ab.df,summary(aov(glm(p.path.a~depth*treatment))))
with(ab.df,aov(glm(p.path.a~depth*treatment))$coefficients)
with(ab.df,aov(glm(p.path.a~depth*treatment))$coefficients["treatmentOrg3"]/aov(glm(p.path.a~depth*treatment))$coefficients["(Intercept)"]) # calculate the percent increase in plant pathogens due to organic treatment relative to no-till

pathGradient<-taxenvmodel(mp.path1)
pathGradient$tab$Fusarium_acutatum
pathGradient$fit$Fusarium_acutatum
pathGradient$tab$Alternaria_tenuissima
pathGradient$fit$Alternaria_tenuissima

emp.path$fit$Fusarium_acutatum
emp.path$fit$Alternaria_tenuissima
emp.path$fit$Paraphoma_pye # 2800% higher in ORG than NT
emp.path$fit$Paraphoma_radicina # 258% higher in ORG than NT
emp.path$fit$Alternaria_dactylidicola # 471% higher in ORG than in NT
emp.path$fit$Paraphoma_rhaphiolepidis # 1926% higher in ORG than NT
pathGradient$fit$Paraphoma_pye
pathGradient$tab$Paraphoma_pye
pathGradient$fit$Paraphoma_radicina
pathGradient$tab$Paraphoma_radicina
pathGradient$fit$Paraphoma_rhaphiolepidis
pathGradient$tab$Paraphoma_rhaphiolepidis
pathGradient$fit$Alternaria_dactylidicola
pathGradient$tab$Alternaria_dactylidicola
# network analysis ####

fQ.NT<-subset_samples(fGA.Qlf, Treatment=="NT")
fQ.CT<-subset_samples(fGA.Qlf, Treatment=="CT")
fQ.ORG<-subset_samples(fGA.Qlf, Treatment=="Org3")
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

fFS.associationSummary[1,1]<-sum(f.NTnet>0.3 & f.NTnet<1)
fFS.associationSummary[2,1]<-sum(f.CTnet>0.3 & f.CTnet<1)
fFS.associationSummary[3,1]<-sum(f.Org3net>0.3 & f.Org3net<1)
fFS.associationSummary[1,2]<-sum(-0.3>f.NTnet)
fFS.associationSummary[2,2]<-sum(-0.3>f.CTnet)
fFS.associationSummary[3,2]<-sum(-0.3>f.Org3net)
fFS.associationSummary

p.guild.int.NT<-ismatrix(f.NTnet, fQ.NTf, "P")
p.guild.int.CT<-ismatrix(f.CTnet, fQ.CTf, "P")
p.guild.int.O3<-ismatrix(f.Org3net, fQ.ORGf, "P")
n.guild.int.NT<-ismatrix(f.NTnet, fQ.NTf, "N")
n.guild.int.CT<-ismatrix(f.CTnet, fQ.CTf, "N")
n.guild.int.O3<-ismatrix(f.Org3net, fQ.ORGf, "N")

#View(p.guild.int.NT)
#View(p.guild.int.CT)
#View(p.guild.int.O3)
#View(n.guild.int.NT)
#View(n.guild.int.CT)
#View(n.guild.int.O3)

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

plotIgraphNet<-function(net, ps, main){
  require(RColorBrewer)
  require(igraph)
  require(phyloseq)
  require(dplyr)
  mdf<-as.data.frame(as.matrix(tax_table(ps)))
  mdf$id<-rownames(mdf)
  mdf$id<-as.character(mdf$id)
  updatedn2<-igraph::set_vertex_attr(net, 
                                     name="primary_lifestyle",
                                     index=V(net),
                                     value = sapply(V(net)$name, function(x){
                                       mdf$primary_lifestyle[mdf$id==x]
                                     }))
  l<-layout_with_dh(updatedn2)
  nb.cols<-length(unique(V(updatedn2)$primary_lifestyle))
  mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)
  coul<-mycolors[as.numeric(as.factor(V(updatedn2)$primary_lifestyle))]
  plot(updatedn2, #layout=l,
       #laout=layout.bipartite,
       layout=layout.circle,
       vertex.label=NA, main=main, vertex.size=5,edge.arrow.size=0.5,vertex.color=coul)
  legend("bottomleft", legend=levels(as.factor(V(updatedn2)$primary_lifestyle)), col = mycolors, bty = "n", text.col=mycolors, horiz = FALSE)
  updatedn2
}


c<-f.NTnet#*(f.NTnet>0.2 + (-0.2 > f.NTnet))
c[c>(-.4)&c<0.4]<-0
#sum(c)
n<-graph_from_incidence_matrix(c, directed=T, mode="out", weighted=T)
#cfg<-(n)
n<-delete.vertices(simplify(n), degree(n)==0)
E(n)$color <- ifelse(E(n)$weight > 0,'navy','maroon')


c2<-f.CTnet
c2[c2>(-.3)&c2<0.3]<-0
#sum(c2)
n2<-graph_from_incidence_matrix(c2, directed=T, mode="out", weighted=T)
#cfg<-(n)
n2<-delete.vertices(simplify(n2), degree(n2)==0)
E(n2)$color <- ifelse(E(n2)$weight > 0,'navy','maroon')

#closeness(n2)
hist(degree(n2))
#tc <- cluster_walktrap(n2)
#modularity(n2, unlist(tc))

c3<-f.Org3net
c3[c3>(-.3)&c3<0.3&c3!=1]<-0
#sum(c3)0.1
n3<-graph_from_incidence_matrix(c3, directed=T, mode="out", weighted=T)
n3<-delete.vertices(simplify(n3), degree(n3)==0)
l<-layout_with_dh(n3)
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
E(n3)$color <- ifelse(E(n3)$weight > 0,'navy','maroon')

plotIgraphNet(n, fQ.NTf, main="NT")
plotIgraphNet(n2, fQ.CTf, main="CT")
plotIgraphNet(n3, fQ.ORGf, main="ORG")
