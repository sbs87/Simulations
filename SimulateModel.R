rm(list=ls())
source("/Users/stevensmith/bin/R_source_functions/ngs.R")
setwd("/Users/stevensmith/bin/Qualifiers/Simulations")
## each subject has her own subject-specfic effect. There are 50 subjects total.
N<-50

## There are 6 samples taken from each indigidual. 
n<-6

## Eavery 6th sample belongs to a different subject. 
(subject_table<-data.frame(subj=rep(1:N,each=n),subject.effect=rep(-99,times=n*N),BV.status=rep(rep(c(0,1),each=3),times=N)))
for(i in 1:N){
  subject_table[subject_table$subj==i,]$subject.effect<-rnorm(n,mean=runif(1,min=0,max=2),sd=runif(1))
}
(subject_table)

## There are 500 miRNAs that sequencing at once. Only 50 of them are truly differentially expressed
m<-50
k<-0.1*m ## Truly differential expressed

(miRNA<-data.frame(miRNA.id=seq(1:m),bio.effect=rnorm(m,mean=runif(1,min=0,max=1),sd=runif(1))))
miRNA$exp.pbv<-rpois(m,6)
miRNA$exp.nbv<-rpois(m,3)
miRNA$BV<-c(rep(1,times=k),rep(0,times=m-k))
(miRNA)
(counts<-data.frame(miRNA_id=rep(1:m,each=N*6),subject_id=rep(rep(1:N,each=6),times=m),sample_id=rep(1:n,times=m*N),miRNA.expression=-99))

for(miRNA_n in 1:m){
  #miRNA_n<-10
  print(str(miRNA_n))
  for(subj_n in 1:N){
    #subj_n<-3
    for(sample_n in 1:n){
      #sample_n<-1
      currentIter.miRNA<-miRNA[miRNA$miRNA.id==miRNA_n,]
      currentIter.subj.sample.row<-subject_table[subject_table$subj==subj_n,][sample_n,]
      miRNA_expression_sample_bit = miRNA_expression_sample<--1000
      miRNA_expression_sample_bit<-1*currentIter.subj.sample.row$BV.status+2*currentIter.miRNA$BV
      switch(miRNA_expression_sample_bit+1,miRNA_expression_sample<-currentIter.miRNA$exp.nbv,
             miRNA_expression_sample<-currentIter.miRNA$exp.nbv,
             miRNA_expression_sample<-currentIter.miRNA$exp.nbv,
             miRNA_expression_sample<-currentIter.miRNA$exp.pbv)
      counts[(counts$miRNA_id==miRNA_n) & (counts$subject_id==subj_n) & (counts$sample_id==sample_n),]$miRNA.expression <-
        abs(currentIter.subj.sample.row$subject.effect) + abs(currentIter.miRNA$bio.effect) + miRNA_expression_sample
    }
  }
}
(counts)

head(counts)
counts$BV<-rep(rep(c(0,1),each=3),times=N*m)
attach(counts)
counts$uniqid<-paste("mir",miRNA_id,"_SID_",subject_id,"_",sample_id,sep="")
counts$subj_rep<-paste("SID_",subject_id,"_",sample_id,sep="")
detach(counts)

simulated_miRNA_counts<-counts
save(list=c("N","n","m","k","simulated_miRNA_counts"),file="simulated_miRNA_counts_objects.R")

##------------------------------------------------------
## Start from here
rm(list=ls())
source("/Users/stevensmith/bin/R_source_functions/ngs.R")
setwd("/Users/stevensmith/bin/Qualifiers/Simulations")
load("simulated_miRNA_counts_objects.R")
counts<-simulated_miRNA_counts

(counts.h<-head(counts,n=10*N*n))

ggplot(counts.h,aes(x=miRNA_id,y=miRNA.expression,fill=factor(BV)))+geom_bar(stat="identity")+facet_wrap(~BV,ncol=1)

#attach(counts.h)
#CrossTable(miRNA_id,miRNA.expression,prop.r = F,prop.c=F,prop.t = F,prop.chisq = F)
#detach(counts.h)
# Histogram Colored (blue and red)
#hist(miRNA[miRNA$BV==0,]$exp.raw, col=rgb(1,0,0,0.5),xlim=c(0,1.2*max(miRNA$exp.raw)), main="BV vs Non BV Raw Expression", xlab="Expression")
#hist(miRNA[miRNA$BV==1,]$exp.raw, col=rgb(0,0,1,0.5), add=T)
#box()
library("reshape")

counts.reshape<-cast(counts,formula = miRNA_id~subj_rep,value="miRNA.expression")
head(counts.reshape)
## Next, test ability of model to correctly identify differentially expressed genes (top 50) given other sources of varibility

## As a start, use edgeR, which allows a GLM approach
library(edgeR)
group<-rep(c(rep(0,times=3),rep(1,times=3)),times=ncol(counts.reshape)/6)
dge<-DGEList(counts=counts.reshape,group=group)
dge<-estimateCommonDisp(dge)
dge<-estimateTagwiseDisp(dge)
et<-exactTest(dge)
topTags(et)
subject<-rep(1:50,each=6)
design<-model.matrix(~subject+group)
dge<-estimateGLMCommonDisp(dge)
fit<-glmFit(counts.reshape,design,dispersion = estimateGLMCommonDisp(counts.reshape,design))
topTags(glmLRT(fit,coef=2))
counts$BV<-as.numeric(counts$BV)
counts.mirna11<-counts[counts$miRNA_id==11,]
counts.mirna11$miRNA.expression<-round(counts.mirna11$miRNA.expression)
glm(miRNA.expression~BV,family=neg.bin,data=counts.mirna11)
warnings()

## ROC analysis
specificty<-tp/(tp+fn)
sensititivty<-fp(fp+tn)

load("/Users/stevensmith/Documents/School/Maryland/D4_CT_rep1.easyRNAseq.transcriptCounts")
head (la.counts)
la.counts[grep("ENST00000385060",row.names(la.counts)),]
test<-data.frame(cbind(la.counts,la.counts))

library(limma)
limmaUsersGuide()

## model
## the model as I;ve defined it is:
y_i=f(S_j,BVS_k,CS_k,BVD_j,P_k,B_k,e_k)

#y is modeled as log2 transformed miRNA read counts and initially assumed to follow a Negative Binomial distribution with parameters μ (per miRNA and sample normalized expected value) and dispersion α, but can also follow a Normal distribution (and will be modeled as such if found to fit this distribution)
#S_j is effect of subject j and is assumed to be N(0,σ_j^2). The σ_j^2 can come from a seperate distribiton as well. This will be the main crux of the simulation: what is the effect of changing this parameter on the degree of specificity and sensitivity
#BVS_k the effect of BV state (negative or positive) of sample k within subject j. There will be some variability to this as well, to allow to class misassignment
#CS_k is effect of community state type (I,II,III,IV-a,IV-b,V, see ref [16]) of sample k within subject j. Same as above, allow for class misassignment
#BVD_j is the effect of subject-specific chronological direction of BV transition (from NBV to PBV or PBV to NBV). This can be simulated as anothe random process depending on complexity. For now, leave it out. 
#P_k and B_k are proxies for sample k’s physiological (i.e., menstruation) or behavioral (i.e., sexual activity) effects, respectively
#and e_k is the error term. Since both BV and community state assignments may have inherent biases or inaccuracies, the 〖BVS〗_kand 〖CS〗_kterms will be modeled as weighted functions depending on the confidence of group assignment, i.e., 〖BVS〗_k=〖bw〗_k*〖[BVS==b]〗_kand 〖CS〗_k=〖cw〗_k*〖[BVS==c]〗_k, where 〖bw〗_k and 〖cw〗_k are a weighted measure of the confidence in BV and community state assignments b and c, respectively, and [X==x] is encoding for each respective state type. The weight will be a distance metric such as Euclidian distance from group assignment centroid. 

library(DESeq2)
#source("http://bioconductor.org/biocLite.R")
library("parathyroidSE")
data(parathyroidGenesSE)
se<-parathyroidGenesSE
colnames(se)<-colData(se)$run
ddsPara<-DESeqDataSet(se=se,design=~patient+treatment)
colData(ddsPara)$treatment <-factor(colData(ddsPara)$treatment,levels=c("Control","DPN","OHT"))
ddsPara

library(Biobase)
library(pasilla)
data(pasillaGenes)
countData<-counts(pasillaGenes)
colData<-pData(pasillaGenes)[,c("condition","type")]
dds<-DESeqDataSetFromMatrix(countData=countData,
                            colData=colData,
                            design=~condition
                            )
colData(dds)$condition<-factor(colData(dds)$condition,levels=c("untreated","treated"))
dds

#colData(dds)$condition<-relevel(colData(dds)$condition,"control")
dds<-DESeq(dds)
res<-results(dds)
res<-res[order(res$padj),]
head(res)
plotMA(dds,ylim=c(-2,2),main="DESeq2")
mcols(res,use.names=T)
colData(dds)
design(dds)<-formula(~type+condition)
dds<-DESeq(dds)
res<-results(dds)
head(res)
rld<-rlogTransformation(dds,blind=T)
print(plotPCA(rld,intgroup=c("condition","type")))
nbinomWaldTest

DESeq::fitNbinomGLMs  
DESeq::fitNbinomGLMsForMatrix
