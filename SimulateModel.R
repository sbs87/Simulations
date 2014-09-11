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