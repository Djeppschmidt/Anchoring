# GLUSEEN
# import data
setwd("")
GADFUN_RAW.rds<-readRDS("")
GADBac_RAW.rds<-readRDS("")
GAD_meta<-read.csv("")

# scale by qpcr

# scale by RA

# run HMSC on each


# Anchoring Study
tdat$PLFA.Gram.Negative<-as.numeric(as.character(tdat$PLFA.Gram.Negative))
tdat$PLFA.Gram.Positive<-as.numeric(as.character(tdat$PLFA.Gram.Positive))
tdat$NLFA.Gram.Negative<-as.numeric(as.character(tdat$NLFA.Gram.Negative))
tdat$NLFA.Gram.Positive<-as.numeric(as.character(tdat$NLFA.Gram.Positive))
tdat$NLFA.Tbac<-as.numeric(as.character(tdat$NLFA.Tbac))
tdat$PLFA.Tbac<-as.numeric(as.character(tdat$PLFA.Tbac))
tdat$Quantity..picograms.<-as.numeric(as.character(tdat$Quantity..picograms.))
tdat$NLFAT
plot(tdat$NLFA.Tbac~tdat$PLFA.Tbac)
plot(tdat$NLFA.Tbac~tdat$Quantity..picograms.)
plot(tdat$PLFA.Tbac~tdat$Quantity..picograms.)

# QC data:
# Filter to sites with QPCR, PH, OM, NFLA, PLFA
# Calculate the sequencing depth for each
# Run for each farming system [so there are replicates]
# 

# reconcile metadata
library(dplyr)
library(plyr)
PLFA<-as.data.frame(read.csv("/Users/dietrich/Desktop/PhD/PHP/Hmsc_reproducible/PLFA.csv")) # from excel sheets provided by Jude Maul
Chem<-as.data.frame(read.csv("/Users/dietrich/Desktop/PhD/PHP/key.csv")) # sheet labeled key
CFU<-as.data.frame(read.csv("/Users/dietrich/Desktop/PhD/PHP/CFU.csv"))
meta<-as.data.frame(as.matrix(sample_data(bac))) # from workflow output on scinet (prepared by Ryan?)

meta<-meta[meta$System.loc=="Org_6" | meta$System.loc=="Org_3" | meta$System.loc=="CT-MD" | meta$System.loc=="NT-MD",]

A.Belt<-subset_samples(bac, Site=="Beltsville" & Soil_Zone=="rhizosphere" & year == "2013")

# make data match formats between data tables
PLFA$TimeOfyear<-as.character(PLFA$TimeOfyear)
PLFA$TimeOfyear[PLFA$TimeOfyear=="Pre"]<-"pre"
PLFA$TimeOfyear[PLFA$TimeOfyear=="Post"]<-"post"
PLFA$genotype<-as.character(PLFA$genotype)
PLFA$genotype[PLFA$genotype=="RR2"]<-"RR"
PLFA$genotype[PLFA$genotype=="AF2"]<-"Non_RR"
PLFA$genotype[PLFA$genotype=="allen"]<-"RR"
PLFA$genotype[PLFA$genotype=="usg56"]<-"Non_RR"
PLFA$crop<-as.character(PLFA$crop)
PLFA$crop[PLFA$crop=="Corn"]<-"corn"
PLFA$crop[PLFA$crop=="Soy"]<-"soy"
PLFA$crop[PLFA$crop=="Soybean"]<-"soy"
PLFA$Glyphosphate.Treatment<-as.character(PLFA$Glyphosphate.Treatment)
PLFA$Glyphosphate.Treatment[PLFA$Glyphosphate.Treatment=="+"]<-"spray"
PLFA$Glyphosphate.Treatment[PLFA$Glyphosphate.Treatment=="-"]<-"no_spray"
PLFA$Bl.Rh<-as.character(PLFA$Bl.Rh)
PLFA$Bl.Rh[PLFA$Bl.Rh=="Bulk"]<-"bulk"
PLFA$Bl.Rh[PLFA$Bl.Rh=="Rhizo"]<-"rhizosphere"
PLFA$system.ID<-as.character(PLFA$system.ID)
PLFA$system.ID[PLFA$system.ID=="CT"]<-"CT-MD"
PLFA$system.ID[PLFA$system.ID=="NT"]<-"NT-MD"
PLFA$system.ID[PLFA$system.ID=="Org 6"]<-"Org_6"
PLFA$system.ID[PLFA$system.ID=="Org 3"]<-"Org_3"
PLFA$Bl.Rh[PLFA$Bl.Rh==""]<-NA

NLFA<-NLFA[NLFA$MarkerType=="nlfa",]
NLFA2013<-NLFA[NLFA$year==2013,]
NLFA2014<-NLFA[NLFA$year==2014,]

PLFA<-PLFA[PLFA$MarkerType=="plfa",]
PLFA2013<-PLFA[PLFA$year==2013,]
PLFA2014<-PLFA[PLFA$year==2014,]
PLFA$FSPplot<-as.factor(PLFA$FSPplot)

meta2013<-meta[meta$year==2013,]
meta2014<-meta[meta$year==2014,]

PLFA2013$Merge<-paste0(PLFA2013$FSPplot, PLFA2013$Glyphosphate.Treatment, PLFA2013$TimeOfyear, PLFA2013$system.ID, PLFA2013$genotype, PLFA2013$crop, sep="")
meta2013$Merge<-paste0(meta2013$Loc_plot_ID, meta2013$Glyphosphate_Treatment, meta2013$Sampling_date, meta2013$System.loc, meta2013$genotype, meta2013$crop, sep="")

PLFA2014$Merge<-paste0(PLFA2014$FSPplot, PLFA2014$Glyphosphate.Treatment, PLFA2014$TimeOfyear, PLFA2014$Bl.Rh, PLFA2014$system.ID, PLFA2014$genotype, PLFA2014$crop, sep="")
meta2014$Merge<-paste0(meta2014$Loc_plot_ID, meta2014$Glyphosphate_Treatment, meta2014$Sampling_date, meta2014$Soil_Zone, meta2014$System.loc, meta2014$genotype, meta2014$crop, sep="")

m.PLFA2013<-merge(meta2013, PLFA2013, by="Merge", all.x=TRUE)
m.PLFA2014<-merge(meta2014, PLFA2014, by="Merge", all.x=TRUE)

PLFA2<-rbind(m.PLFA2013,m.PLFA2014)
colnames(PLFA2)
PLFA2<-PLFA2[,c(2,5:15,20,39,42,52:60)]
names(PLFA2)[17]<-"PLFA.General.FAME"
names(PLFA2)[18]<-"PLFA.AM.Fungi"
names(PLFA2)[19]<-"PLFA.Gram.Negative"
names(PLFA2)[20]<-"PLFA.Fungi"
names(PLFA2)[21]<-"PLFA.Gram.Positive"
names(PLFA2)[22]<-"PLFA.Anaerobe"
names(PLFA2)[23]<-"PLFA.Actinomycetes"
names(PLFA2)[24]<-"PLFA.Protozoa"
colnames(PLFA2)
View(PLFA2) # PLFA2 is a reconciled dataset!
RH.PLFA2<-PLFA2[PLFA2$Soil_Zone=="rhizosphere",]

NLFA2013$Merge<-paste0(NLFA2013$FSPplot, NLFA2013$Glyphosphate.Treatment, NLFA2013$TimeOfyear, NLFA2013$system.ID, NLFA2013$genotype, NLFA2013$crop, sep="")

NLFA2014$Merge<-paste0(NLFA2014$FSPplot, NLFA2014$Glyphosphate.Treatment, NLFA2014$TimeOfyear, NLFA2014$Bl.Rh, NLFA2014$system.ID, NLFA2014$genotype, NLFA2014$crop, sep="")

m.NLFA2013<-merge(meta2013, NLFA2013, by="Merge", all.x=TRUE)
m.NLFA2014<-merge(meta2014, NLFA2014, by="Merge", all.x=TRUE)

NLFA2<-rbind(m.NLFA2013,m.NLFA2014)

# select columns and change rownames
colnames(NLFA2)
NLFA2<-NLFA2[,c(2,5:15,20,39,42,52:60)]
names(NLFA2)[17]<-"NLFA.General.FAME"
names(NLFA2)[18]<-"NLFA.AM.Fungi"
names(NLFA2)[19]<-"NLFA.Gram.Negative"
names(NLFA2)[20]<-"NLFA.Fungi"
names(NLFA2)[21]<-"NLFA.Gram.Positive"
names(NLFA2)[22]<-"NLFA.Anaerobe"
names(NLFA2)[23]<-"NLFA.Actinomycetes"
names(NLFA2)[24]<-"NLFA.Protozoa"
colnames(NLFA2)
NLFA2<-NLFA2[NLFA2$Soil_Zone=="rhizosphere",] # reconciled!!

View(Chem)
Chem<-Chem[-c(96:101),]
Chem$Glyphosphate_Treatment<-Chem$Microplot.treatment
Chem$Glyphosphate_Treatment<-as.character(Chem$Glyphosphate_Treatment)
Chem$Glyphosphate_Treatment[Chem$Glyphosphate_Treatment=="allen"]<-"no_spray"
Chem$Glyphosphate_Treatment[Chem$Glyphosphate_Treatment=="allen+gly"]<-"spray"
Chem$Glyphosphate_Treatment[Chem$Glyphosphate_Treatment=="usg56"]<-"no_spray"
Chem$Glyphosphate_Treatment[Chem$Glyphosphate_Treatment=="RR2"]<-"no_spray"
Chem$Glyphosphate_Treatment[Chem$Glyphosphate_Treatment=="RR2+gly"]<-"spray"
Chem$Glyphosphate_Treatment[Chem$Glyphosphate_Treatment=="AF2"]<-"no_spray"
Chem$Microplot.treatment<-as.character(Chem$Microplot.treatment)
Chem$Microplot.treatment[Chem$Microplot.treatment=="allen"]<-"RR"
Chem$Microplot.treatment[Chem$Microplot.treatment=="allen+gly"]<-"RR"
Chem$Microplot.treatment[Chem$Microplot.treatment=="usg56"]<-"Non_RR"
Chem$Microplot.treatment[Chem$Microplot.treatment=="RR2"]<-"RR"
Chem$Microplot.treatment[Chem$Microplot.treatment=="RR2+gly"]<-"RR"
Chem$Microplot.treatment[Chem$Microplot.treatment=="AF2"]<-"Non_RR"

Chem$crop<-as.character(Chem$crop)
Chem$crop[Chem$crop=="soybean"]<-"soy"
Chem$system<-as.character(Chem$system)
Chem$system[Chem$system=="ct"]<-"CT-MD"
Chem$system[Chem$system=="nt"]<-"NT-MD"
Chem$system[Chem$system=="Org6"]<-"Org_6"
Chem$system[Chem$system=="Org3"]<-"Org_3"
Chem$Plot.number<-as.character(Chem$Plot.number)

meta2013<-meta[meta$year==2013,]
meta2014<-meta[meta$year==2014,]
plyr::ddply(Chem, c("Plot.number"), summarise,
            OM... = mean(OM...),
            pH = mean(pH))->Chem2014

meta2014<-merge(meta2014, Chem2014, by.x="Loc_plot_ID", by.y="Plot.number")

Chem$MergeC<-paste0(Chem$Plot.number, Chem$Glyphosphate_Treatment, Chem$system, Chem$Microplot.treatment, Chem$crop,sep="")
#Chem$Merge2014<-paste0(Chem$Plot.number, Chem$Glyphosphate_Treatment, Chem$system, Chem$Microplot.treatment,sep="")
meta$MergeC<-paste0(meta$Loc_plot_ID, meta$Glyphosphate_Treatment, meta$System.loc, meta$genotype, meta$crop, sep="")

setdiff(Chem$MergeC, meta2013$MergeC) # should be = 0

Chem2<-merge(meta2013, Chem, by="MergeC", all.x=TRUE)
RH.Chem3<-bind_rows(RH.Chem2, meta2014)
RH.Chem3<-RH.Chem3[RH.Chem3$Soil_Zone=="rhizosphere",]

# merge into phyloseq
RH.meta<-left_join(NLFA2, RH.PLFA2)
RH.meta<-left_join(RH.meta, RH.Chem3)

rownames(RH.meta)<-RH.meta$X.SampleID
RH<-merge_phyloseq(RH, sample_data(RH.meta))

Belt<-merge_phyloseq(Belt, sample_data(B.meta))


b.RH<-readRDS("Data/BacterialRhizospherePhyloseq18112020.rds")

View(as.data.frame(as.matrix(sample_data(b.RH))))

ARH<-subset_samples(b.RH, year=="2013" & soil_zone=="rhizosphere")
# subset taxa to exclude mitochondria etc; and ...
zz<-
# measure total bacterial by NLFA
sample_data(b.RH)$NLFA.Tbac<-sample_data(b.RH)$NLFA.Gram.Positive+sample_data(b.RH)$NLFA.Gram.Negative
# measure fungal to bacterial proportion by PLFA
sample_data(b.RH)$PLFA.Tbac<-sample_data(b.RH)$PLFA.Gram.Positive+sample_data(b.RH)$PLFA.Gram.Negative
# define total counts per sample (as a function of anchor value * relative total PLFA)
hist(sample_sums(subset_taxa(b.RH, Kingdom=="Bacteria"))/(100*(sample_data(b.RH)$NLFA.Tbac)), breaks=30, main="bacterial sequence counts per unit NLFA", xlab="Sequence counts per unit NLFA")
hist(sample_sums(subset_taxa(b.RH, Kingdom=="Bacteria"))/(100*(sample_data(b.RH)$PLFA.Tbac)), breaks=30, main="bacterial sequence counts per unit NLFA", xlab="Sequence counts per unit NLFA")


trfm<-function(ps, val){
  scaled<-data.frame(mapply(`*`, data.frame(t(as.matrix(otu_table(transform_sample_counts(ps, function(x) x/sum(x)))))), val))# sample_data(ps)$val))
  
  rnames<-rownames(data.frame(t(as.matrix(otu_table(ps)))))
  cnames<-sample_names(ps)
  rownames(scaled)<-rnames
  colnames(scaled)<-cnames
  scaled<-round(scaled)
  
  p2<-ps
  otu_table(p2)<- otu_table(scaled, taxa_are_rows=T)
  p2
}

# transform data into the final version!!
RAW.RH<-b.RH
RA.RH<-transform_sample_counts(b.RH, function(x) x/sum(x))
NLFA.RH<-trfm(b.RH, 100*(sample_data(b.RH)$NLFA.Tbac))
PLFA.RH<-trfm(b.RH, 100*(sample_data(b.RH)$PLFA.Tbac))
QPCR.RH<-trfm(b.RH, sample_data(b.RH)$QPCR)

b2.RH<-prune_taxa(taxa_sums(b2.RH) > 0, b2.RH)  


RA.RH.CT<-subset_samples()
RA.RH.NT<-subset_samples()
RA.RH.O3<-subset_samples()
RA.RH.O6<-subset_samples()

RAW.RH.CT<-subset_samples()
RAW.RH.NT<-subset_samples()
RAW.RH.O3<-subset_samples()
RAW.RH.O6<-subset_samples()

NLFA.RH.CT<-subset_samples()
NLFA.RH.NT<-subset_samples()
NLFA.RH.O3<-subset_samples()
NLFA.RH.O6<-subset_samples()

PLFA.RH.CT<-subset_samples()
PLFA.RH.NT<-subset_samples()
PLFA.RH.O3<-subset_samples()
PLFA.RH.O6<-subset_samples()

QPCR.RH.CT<-subset_samples()
QPCR.RH.NT<-subset_samples()
QPCR.RH.O3<-subset_samples()
QPCR.RH.O6<-subset_samples()

l.Anchor<-list(RAW.RH.CT, RAW.RH.NT, RAW.RH.O3, RAW.RH.O6, RA.RH.CT, RA.RH.NT, RA.RH.O3, RA.RH.O6,NLFA.RH.CT, NLFA.RH.NT, NLFA.RH.O3, NLFA.RH.O6,PLFA.RH.CT, PLFA.RH.NT, PLFA.RH.O3, PLFA.RH.O6,QPCR.RH.CT, QPCR.RH.NT, QPCR.RH.O3, QPCR.RH.O6)


saveRDS(fb.RH, "Data/fb_Rhizosphere_2020.RDS")

# HMSC work !!!!

# PHP model
fb.RH<-tax_glom(fb.RH, taxrank="Genus")
  
  Ydat<-as.data.frame(as.matrix(otu_table(ps)))
  XData<-as.data.frame(as.matrix(sample_data(ps)), stringsAsFactors = TRUE)
  XDat1<-XData[,c(2,3,6,8:10,15,47,75,81)] # subset data to what I need for this run!!
  rownames(XDat1)<-c(1:nrow(XDat1))
  XFormula= ~Glyphosphate_Treatment + genotype + System.loc + crop + pH + OM... # consider genotype 
  XDat1$Sample<-as.factor(c(1:nrow(XDat1)))
  studyDesign = data.frame("Sample"=XDat1$Sample,"Sampling_date"=XDat1$Sampling_date, "Depth"=XDat1$depth)#
  
    PHP<-NULL
    PHP$probit<-NULL
    PHP$LogPoi<-NULL
    rL1 <- HmscRandomLevel(units=studyDesign$Sample)
    rL2 <- HmscRandomLevel(units=studyDesign$Sampling_date) # set random level to loc_plot_ID
    rL3 <- HmscRandomLevel(units=studyDesign$depth)
    bf<-Hmsc(Y=as.matrix(Ydat), XData=XDat1, XFormula=XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Sampling_date"=rL2), distr="probit")#, "year"=rL3), distr="probit") 
    PHP$probit$bf<-sampleMcmc(bf, thin=2, samples=1000, transient=500, nChains=1, nParallel=1, verbose=100)
    
    # examine correlation matrix for probit model
    # get PHPputs for HMSC analysis
    PHP$probit$OmegaCor.probit=computeAssociations(PHP$probit$bf)

    PHP$probit$preds<-computePredictedValues(PHP$probit$bf)

    PHP$probit$MF<-evaluateModelFit(PHP$probit$bf, predY=PHP$probit$preds)
  
    print(PHP$probit$bf$X)
    PHP$probit$VP<-computeVariancePartitioning(PHP$probit$bf, group=c(1,1), groupnames=c("Glyphosate"))
    PHP$probit$postBeta<-getPostEstimate(PHP$probit$bf, parName="Beta") #beta is species abundance ; gamma is traits; rho is phylogenetic signal
    
    Y2<-as.matrix(Ydat) # global OTU table
    Y2[Y2==0]<-NA # remove zero counts
    PHP$LogPoi$bf<-Hmsc(Y=Y2, XData=XDat1, XFormula=XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Sampling_date"=rL2), distr="poisson") # use same formulas as previous
    PHP$LogPoi$bf<-sampleMcmc(PHP$LogPoi$bf, thin=2, samples=1000, transient=500, nChains=1, nParallel=1, verbose=100)

    PHP$LogPoi$OmegaCor.lp=computeAssociations(PHP$LogPoi$bf)

    PHP$LogPoi$preds<-computePredictedValues(PHP$LogPoi$bf)

    PHP$LogPoi$MF<-evaluateModelFit(PHP$LogPoi$bf, predY=PHP$LogPoi$preds)

    if(AM==0){PHP$LogPoi$VP<-computeVariancePartitioning(PHP$LogPoi$bf, group=c(1,1), groupnames=c("Glyphosate"))}
    if(AM==1){PHP$LogPoi$VP<-computeVariancePartitioning(PHP$LogPoi$bf, group=c(1,1,2), groupnames=c("Glyphosate", "AMF"))}

    PHP$LogPoi$postBeta<-getPostEstimate(PHP$LogPoi$bf, parName="Beta") #beta is species abundance ; gamma is traits; rho is phylogenetic signal

    PHP
  }
  
  out$REF<-run(Ydat=Ydat, XDat=XDat1, XFormula=XFormula1, studyDesign=studyDesign, AM=0)
  out$AM<-run(Ydat=Ydat, XDat=XDat1, XFormula=XFormula2, studyDesign=studyDesign, AM=1)
  out$ps<-ps
  out
}