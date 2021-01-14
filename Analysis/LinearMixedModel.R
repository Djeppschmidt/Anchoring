# Residual matrix network analysis platform

# linear mixed effects model (categories and covariates)
# doesn't need a random effect, but can be included


library(phyloseq)
library(lme4)
library(ggplot2)
library(dplyr)
library(reshape2)
library(doParallel)
library(foreach)

# d<-readRDS("/Volumes/Seagate\ Expansion\ Drive/Illumina/Processing/DES/GAD/GADbac2020.rds")
# fd<-readRDS("/Volumes/Seagate\ Expansion\ Drive/Illumina/Processing/DES/GAD/GADfun2020.rds")

# import fungal-bacterial dataset
GAD.bac<-readRDS("/Users/dietrich/Documents/GitHub/Anchoring/Data/GAD/GADbac2020.rds")
GAD.Fun<-readRDS("/Users/dietrich/Documents/GitHub/Anchoring/Data/GAD/GADfun2020.rds")
sample_data(GAD.bac)$SeqDepth<-sample_sums(GAD.bac)
sample_data(GAD.Fun)$SeqDepth<-sample_sums(GAD.Fun)

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

# normalize by QPCR:
fGA.Q<-GAD.QScale(GAD.Fun,type="F")
bGA.Q<-GAD.QScale(GAD.bac,type="B")

# filter uncommon species:

#<-filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE)

bGA.Q<-filter_taxa(bGA.Q, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

fGA.Q<-filter_taxa(fGA.Q, function(x) sum(x > 3) > (0.2*length(x)), TRUE)




# steps for the anchoring study:
#fGA.RA<-transform_sample_counts(GAD.Fun, function(x) x/sum(x))
#bGA.RA<-transform_sample_counts(GAD.bac, function(x) x/sum(x))
#GAD<-list(bGA.Q, bGA.RA, GAD.bac, fGA.Q, fGA.RA, GAD.Fun)

# determined that order doesn't seem to make a difference
# gaussian seems to get highest Rsquared
# variance explained by each factor? 

# FSP species abundance model:
# experimental design
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

# bacteria ####
t1<-Sys.time()
mtax.bac<-taxmodel(bGA.Q)
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
  
  print("data set up")
  # prepare environment
  cores<-detectCores(logical=F)
  cores<-cores-1
  cl<-makeCluster(cores)
  registerDoParallel(cl,cores=cores)
  chunk.size<-ncol(d1)/(cores)
  
  # run regressions in parallel
  out<-foreach(i=1:cores, .combine="cbind", .inorder=T) %dopar% { 
    calc<-matrix(ncol=chunk.size, nrow=ncol(d1))
    for(x in ((i-1)*chunk.size+1):(i*chunk.size)){
      
      for(y in 1:ncol(d1)){
        tax<-d1[,y]
        fit<-NULL
        fit <- glm(d1[,x] ~ treatment*depth+pH+Cpercent+Ammonia+Nitrate+tax,family="gaussian")
        tab<-cbind(summary(aov(fit))[[1]],"varExplained"=summary(aov(fit))[[1]]$`Sum Sq`/sum(summary(aov(fit))[[1]]$`Sum Sq`)*100)
        calc[y,x]<-tab$varExplained[7]
      }
      calc
    }
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

sppInt(fGA.Q)



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
