# GLUSEEN Network and facilitation study
library(Hmsc)
library(phyloseq) # need to check to make sure know how to install in scinet
library(parallel)


GLUFUN_RAW<-readRDS("/Users/dietrich/Documents/GitHub/Anchoring/Data/GLUSEEN/GLUSEENFungi_2020.RDS")
# Export fungi for funguilds
# f<-as.data.frame(as.matrix(tax_table(GLUFUN_RAW)))
# write.csv(f, "/Users/dietrich/Documents/GitHub/Anchoring/Data/F_taxonomy.csv")
GLUBac_RAW<-readRDS("/Users/dietrich/Documents/GitHub/Anchoring/Data/GLUSEEN/Bac16s_RAW.RDS")
meta<-as.data.frame(as.matrix(sample_data(GLUBac_RAW)))

GLU_meta<-as.data.frame(as.matrix(read.csv("/Users/dietrich/Documents/GitHub/Anchoring/Data/GLUSEEN/QPCR_Data_100615.csv")))
EW<-as.data.frame(as.matrix(read.csv("/Users/dietrich/Documents/GitHub/Anchoring/Data/GLUSEEN/data_EW_total_biomass.csv", sep=";")))
EW$site<-gsub("Lahti", "Lakti", EW$site)
EW<-EW[,-c(1,3)]
names(EW)[names(EW)=="site"]<-"Codes.on.samples"

#head(GLU_meta)
GLU_meta$Sample_ID<-paste0("GLU0", c(GLU_meta$A))
GLU_meta$Sample_ID
GLU_meta$Sample_ID[GLU_meta$Sample_ID=="GLU0100"]<-"GLU100"
GLU_meta$Sample_ID[GLU_meta$Sample_ID=="GLU01"]<-"GLU001"
GLU_meta$Sample_ID[GLU_meta$Sample_ID=="GLU02"]<-"GLU002"
GLU_meta$Sample_ID[GLU_meta$Sample_ID=="GLU03"]<-"GLU003"
GLU_meta$Sample_ID[GLU_meta$Sample_ID=="GLU04"]<-"GLU004"
GLU_meta$Sample_ID[GLU_meta$Sample_ID=="GLU05"]<-"GLU005"
GLU_meta$Sample_ID[GLU_meta$Sample_ID=="GLU06"]<-"GLU006"
GLU_meta$Sample_ID[GLU_meta$Sample_ID=="GLU07"]<-"GLU007"
GLU_meta$Sample_ID[GLU_meta$Sample_ID=="GLU08"]<-"GLU008"
GLU_meta$Sample_ID[GLU_meta$Sample_ID=="GLU09"]<-"GLU009"
rownames(GLU_meta)<-GLU_meta$Sample_ID
GLU_meta<-left_join(meta, GLU_meta, by="Sample_ID")
GLU_meta<-left_join(GLU_meta, EW, by="Codes.on.samples")
#View(GLU_meta)
rownames(GLU_meta)<-GLU_meta$Sample_ID
sample_data(GLUFUN_RAW)<-GLU_meta
sample_data(GLUBac_RAW)<-GLU_meta
sample_data(GLUFUN_RAW)$SeqDepth<-sample_sums(GLUFUN_RAW)
sample_data(GLUBac_RAW)$SeqDepth<-sample_sums(GLUBac_RAW)


# scale by qpcr
QScale<-function(ps, type){
  
  if(type=="B"){
    scale<-as.numeric(as.character(sample_data(ps)$X16s))
    out<-Qscale(ps, val=1, scale)
  }
  if(type=="F"){
    scale<-as.numeric(as.character(sample_data(ps)$its))
    out<-Qscale(ps, val=1, scale)
  }
  out
}
Qscale<-function(ps, val, scale){
  scaled<-data.frame(mapply(`*`, data.frame(t(as.matrix(otu_table(transform_sample_counts(ps, function(x) x/sum(x)))))), scale * val))# sample_data(ps)$val))
  names<-rownames(data.frame(t(as.matrix(otu_table(ps)))))
  rownames(scaled)<-names
  scaled<-round(scaled)
  
  p2<-ps
  otu_table(p2)<- otu_table(scaled, taxa_are_rows=T)
  p2
  
}
bGLU.Q<-QScale(GLUBac_RAW, type="B") # quantitatively scaled bacteria
fGLU.Q<-QScale(GLUFUN_RAW, type="F") # quantitatively scaled fungi
# merge fungi and bacteria
bfGLU<-merge_phyloseq(bGLU.Q, fGLU.Q)

# merge taxa to genus level
bfGLU<-tax_glom(bfGLU, "Genus")
# filter chlorophyll, archaea, mitochondria
bfGLU<-subset_taxa(bfGLU, Kingdom == "Bacteria" | Kingdom =="Fungi", TRUE)
bfGLU<-subset_taxa(bfGLU, Family != "Mitochondria", TRUE)

# remove samples without enough metadata
bfGLU<-prune_samples(!is.na(sample_data(bfGLU)$C_org), bfGLU)
bfGLU<-prune_samples(!is.na(sample_data(bfGLU)$Ni_avail), bfGLU)
# subset to city groups

Balt.ref<-subset_samples(bfGLU, CITY=="C1" & TRT=="T1")
Balt.rem<-subset_samples(bfGLU, CITY=="C1" & TRT=="T2")
Balt.turf<-subset_samples(bfGLU, CITY=="C1" & TRT=="T3")
Balt.rud<-subset_samples(bfGLU, CITY=="C1" & TRT=="T4")

Bud.ref<-subset_samples(bfGLU, CITY=="C4" & TRT=="T1")
Bud.rem<-subset_samples(bfGLU, CITY=="C4" & TRT=="T2")
Bud.turf<-subset_samples(bfGLU, CITY=="C4" & TRT=="T3")
Bud.rud<-subset_samples(bfGLU, CITY=="C4" & TRT=="T3")

Hels.ref<-subset_samples(bfGLU, CITY=="C2" & TRT=="T1")
Hels.rem<-subset_samples(bfGLU, CITY=="C2" & TRT=="T2")
Hels.turf<-subset_samples(bfGLU, CITY=="C2" & TRT=="T3")
Hels.rud<-subset_samples(bfGLU, CITY=="C2" & TRT=="T4")

Laht.ref<-subset_samples(bfGLU, CITY=="C3" & TRT=="T1")
Laht.rem<-subset_samples(bfGLU, CITY=="C3" & TRT=="T2")
Laht.turf<-subset_samples(bfGLU, CITY=="C3" & TRT=="T3")
Laht.rud<-subset_samples(bfGLU, CITY=="C3" & TRT=="T4")

Potch.ref<-subset_samples(bfGLU, CITY=="C5" & TRT=="T1")
Potch.rem<-subset_samples(bfGLU, CITY=="C5" & TRT=="T2")
Potch.turf<-subset_samples(bfGLU, CITY=="C5" & TRT=="T3")
Potch.rud<-subset_samples(bfGLU, CITY=="C5" & TRT=="T4")

# run HMSC on each
GLU<-list(Balt.ref,Balt.rem,Balt.turf, Balt.rud, Hels.ref,Hels.rem,Hels.turf,Hels.rud,Laht.ref, Laht.rem, Laht.turf, Laht.rud, Bud.ref, Bud.rem, Bud.turf, Bud.rud, Potch.ref, Potch.rem, Potch.turf, Potch.rud)
names(GLU)<-c("Balt.ref","Balt.rem","Balt.turf", "Balt.rud", "Hels.ref","Hels.rem","Hels.turf","Hels.rud","Laht.ref", "Laht.rem", "Laht.turf", "Laht.rud", "Bud.ref", "Bud.rem", "Bud.turf", "Bud.rud", "Potch.ref", "Potch.rem", "Potch.turf", "Potch.rud")
# remove taxa that don't exist in the subnets
GLU<-lapply(GLU, filter_taxa, function(x) sum(x)==0, TRUE)

# GLUSEEN HMSC study

applyGLHMSC<-function(ps){
  out<-NULL
  
  Ydat<-as.data.frame(t(as.matrix(otu_table(ps))))
  XData<-as.data.frame(as.matrix(sample_data(ps)), stringsAsFactors = TRUE)
  XDat1<-XData[,c(6:26)] # subset data to what I need for this run!!
  rownames(XDat1)<-c(1:nrow(XDat1))
  XDat1<-XDat1[,colSums(is.na(XDat1))==0]
  
  XFormula1= ~ pH.H2O + C_org + NO3_N + NH4_N + Ni_avail + Zn_avail 
  XDat1$Sample<-as.factor(c(1:nrow(XDat1)))
  studyDesign = data.frame("Sample"=XDat1$Sample)
  run<-function(Ydat, XDat, XFormula, studyDesign){
    out<-NULL
    out$probit<-NULL
    out$LogPoi<-NULL
    rL1 <- HmscRandomLevel(units=studyDesign$Sample) # reduce model complexity as much as possible for speed of calculations!
    bf<-Hmsc(Y=as.matrix(Ydat), XData=XDat1, XFormula=XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1), distr="probit")
    out$probit$bf<-sampleMcmc(bf, thin=1, samples=1000, transient=500, nChains=1, nParallel=1, verbose=100)
    
    # examine correlation matrix for probit model
    # get outputs for HMSC analysis
    out$probit$OmegaCor.probit=computeAssociations(out$probit$bf)
    out$probit$preds<-computePredictedValues(out$probit$bf)
    out$probit$MF<-evaluateModelFit(out$probit$bf, predY=out$probit$preds)
    out$probit$VP<-computeVariancePartitioning(out$probit$bf, group=c(1,1,2,3,3,4,4), groupnames=c("pH", "OM", "Nitrogen", "Metals")) 
    out$probit$postBeta<-getPostEstimate(out$probit$bf, parName="Beta") 
    
    Y2<-as.matrix(Ydat) # global OTU table
    Y2[Y2==0]<-NA # remove zero counts
    out$LogPoi$bf<-Hmsc(Y=Y2, XData=XDat1, XFormula=XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1), distr="poisson") # use same formulas as previous
    out$LogPoi$bf<-sampleMcmc(out$LogPoi$bf, thin=1, samples=1000, transient=500, nChains=1, nParallel=1, verbose=100)
    print("bf Completed") # troubleshooting
    out$LogPoi$OmegaCor.lp=computeAssociations(out$LogPoi$bf)
    print("associations Completed") # troubleshooting
    out$LogPoi$preds<-computePredictedValues(out$LogPoi$bf)
    print("preds Completed") # troubleshooting
    out$LogPoi$MF<-evaluateModelFit(out$LogPoi$bf, predY=out$LogPoi$preds)
    print("model fit Completed") # troubleshooting
    out$LogPoi$VP<-computeVariancePartitioning(out$LogPoi$bf, group=c(1,1,2,3,3,4,4), groupnames=c("pH", "OM", "Nitrogen", "Metals"))
    out$LogPoi$postBeta<-getPostEstimate(out$LogPoi$bf, parName="Beta") #beta is species abundance ; gamma is traits; rho is phylogenetic signal
    out
  }
  
  out<-run(Ydat=Ydat, XDat=XDat1, XFormula=XFormula1, studyDesign=studyDesign)
  out
}

runGLHMSC<-function(list){
  out<-mclapply(list, applyGLHMSC, simplify=F, USE.NAMES = T, mc.cores=16) # each computation takes a very long time
  out
}

GLU.out<-runGLHMSC(GLU)
saveRDS(GLU.out, "GLUSEEN_Network.RDS")