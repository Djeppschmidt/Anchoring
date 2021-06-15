# using modeling to account for differences ...


# gluseen datasets ####


# FSP datasets ####
FSP.bac<-readRDS("/Users/dietrich/Documents/GitHub/Anchoring/Data/FSP/FSPbac2020.rds")
FSP.Fun<-readRDS("/Users/dietrich/Documents/GitHub/Anchoring/Data/FSP/FSPfun2020.rds")
# Using the Bacterial metadata; different sample order


#fsp.bmeta$bacteria<-as.numeric(as.character(bmeta$Bac_QPCR)) 
#sample_data(FSP.bac)<-sample_data(fsp.bmeta) # return updated metadata to bacterial dataset
FSP.bac<-subset_samples(FSP.bac, Depth=="0_5"|Depth=="5_10"|Depth=="10Ap")# subset the dataset to depths with enough replication
FSP.Fun<-subset_samples(FSP.Fun, Depth=="0_5"|Depth=="5_10"|Depth=="10Ap")
# format data for visualization and modeling
#b.meta<-as.data.frame(as.matrix(sample_data(FSP.bac))) # format data to insert
fsp.bmeta<-as.data.frame(as.matrix(sample_data(FSP.bac)))
fsp.bmeta$Fun_QPCR<-as.numeric(as.character(fsp.bmeta$Fun_QPCR)) 
fsp.bmeta$Bac_QPCR<-as.numeric(as.character(fsp.bmeta$Bac_QPCR))
fsp.bmeta$Depth<-factor(fsp.bmeta$Depth, levels=c("0_5", "5_10", "10Ap"))#, "Ap30")) # put factor levels in the right order
fsp.bmeta$Treatment<-factor(fsp.bmeta$Treatment, levels=c("NT", "CT", "Org3")) # put factor levels in the right order

fsp.fmeta<-as.data.frame(as.matrix(sample_data(FSP.Fun)))
fsp.fmeta$Fun_QPCR<-as.numeric(as.character(fsp.fmeta$Fun_QPCR)) 
fsp.fmeta$Bac_QPCR<-as.numeric(as.character(fsp.fmeta$Bac_QPCR))
fsp.fmeta$Depth<-factor(fsp.fmeta$Depth, levels=c("0_5", "5_10", "10Ap"))# put factor levels in the right order
fsp.fmeta$Treatment<-factor(fsp.fmeta$Treatment, levels=c("NT", "CT", "Org3")) # put factor levels in the right order


# relative abundance vs Qseq abundance
# compare variance explained by model
# compare variance explained by environment
# compare models that account for seq depth to seq effort to no seq info
# explore probit model structure