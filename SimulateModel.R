rm(list=ls())
#source("/Users/stevensmith/bin/R_source_functions/ngs.R")
setwd("/Users/stevensmith/bin/Qualifiers/Simulations")
library("reshape")
library(DESeq2)
library(ggplot2)
## each subject has her own subject-specfic effect. There are 50 subjects total.
N<-50

## There are 6 samples taken from each indigidual. 
n<-6


## There are 500 miRNAs that sequencing at once. Only 50 of them are truly differentially expressed
m<-50
ROC<-20
nsig<-seq(from=.5,to=5,length.out = 5)
spec.sens<-data.frame(sensitivity=rep(0,times=ROC),specificity=rep(0,times=ROC))
(miRNA<-data.frame(miRNA.id=seq(1:m),bio.effect=rnorm(m,mean=runif(1,min=0,max=1),sd=runif(1))))

for(i in 1:m){
  up_or_down<-rbinom(1,1,0.6)
  miRNA$exp.pbv[i]<-rnbinom(1000,mu = runif(1,min = 1,max=1000),size=1/0.001)*(2*up_or_down-1) ## x when up_or_down =1; -x when up_or_down=0
  miRNA$exp.nbv[i]<-rnorm(1,mean = 0,sd = 1)
}

k<-round(.25*m) ## Truly differential expressed
miRNA$BV<-c(rep(1,times=k),rep(0,times=m-k))

(counts<-data.frame(miRNA_id=rep(1:m,each=N*6),subject_id=rep(rep(1:N,each=6),times=m),sample_id=rep(1:n,times=m*N),miRNA.expression=-99,random.effect=rnorm(n*m*N,0,1)))

for (r in 1:ROC){
## Eavery 6th sample belongs to a different subject. 
(subject_table<-data.frame(subj=rep(1:N,each=n),subject.effect=rep(-99,times=n*N),BV.status=rep(rep(c(0,1),each=3),times=N)))
for(i in 1:N){
  #i<-1
  #r<-5
  subject_table[subject_table$subj==i,]$subject.effect<-rnorm(1,mean=0,sd=nsig[r])
}

for(miRNA_n in 1:m){
  print(str(miRNA_n))
  for(subj_n in 1:N){
    for(sample_n in 1:n){
      currentIter.miRNA<-miRNA[miRNA$miRNA.id==miRNA_n,]
      currentIter.subj.sample.row<-subject_table[subject_table$subj==subj_n,][sample_n,]
      miRNA_expression_sample_bit = miRNA_expression_sample<--1000
      miRNA_expression_sample_bit<-1*currentIter.subj.sample.row$BV.status+2*currentIter.miRNA$BV
      switch(miRNA_expression_sample_bit+1,miRNA_expression_sample<-currentIter.miRNA$exp.nbv,#+rnorm(1,0,1),
             miRNA_expression_sample<-currentIter.miRNA$exp.nbv,#+rnorm(1,0,1),
             miRNA_expression_sample<-currentIter.miRNA$exp.nbv,#+rnorm(1,0,1),
             miRNA_expression_sample<-currentIter.miRNA$exp.pbv)#+rnorm(1,0,1))
      u<-exp(miRNA_expression_sample )#+ currentIter.subj.sample.row$subject.effect)
      alpha<-0.01
      (neg_bin_sample<-rnbinom(1000,mu=u,size=1/alpha))
      #hist(neg_bin_sample,main=paste("mirna",miRNA_n,"subj",subj_n,"sample",sample_n,"mu",log(u,base=exp(1))),col='blue')
      #abline(v=neg_bin_sample[1],col='red',lwd=4)
      counts[(counts$miRNA_id==miRNA_n) & (counts$subject_id==subj_n) & (counts$sample_id==sample_n),]$miRNA.expression <-neg_bin_sample[1]
    }
  }
}
#(counts)
#head(counts)
counts$BV<-rep(rep(c(0,1),each=3),times=N*m)
attach(counts)
counts$uniqid<-paste("mir",miRNA_id,"_SID_",subject_id,"_",sample_id,sep="")
counts$subj_rep<-paste("SID_",subject_id,"_",sample_id,sep="")
detach(counts)
hist(log(counts$miRNA.expression))
#(counts.h<-head(counts,n=35*N*n))
#ggplot(counts.h,aes(x=miRNA_id,y=miRNA.expression,fill=factor(BV)))+geom_bar(stat="identity")+facet_wrap(~BV,ncol=1)
counts$miRNA.expression<-counts$miRNA.expression+1
counts.reshape<-cast(counts,formula = miRNA_id~subj_rep,value="miRNA.expression")
head(counts.reshape)

(colData<-data.frame(id=paste("SID_",subject_table$subj,"_",rep(1:6),sep=""),subject_table))
#counts.reshape[counts.reshape<0]<-0
#counts.reshape<-floor(counts.reshape)
deseq_se<-DESeqDataSetFromMatrix(countData=as.matrix(counts.reshape),colData=colData,design=~subj+BV.status)

dds<-DESeq(deseq_se)
res<-results(dds)
res[is.na(res$padj),]$padj<-1
spec.sens[r,]<-calculateSensitivitySpecificity(res,miRNA)
print(str(r))
}

res<-res[order(res$padj),]

res[res$padj<0.01,]
mcols(res,use.names=T)
colData(dds)
rld<-rlogTransformation(dds,blind=T)
print(plotPCA(rld,intgroup=c("BV.status"))) ## Double chec this- the BV points should be a subset of the non BV points. They look about equal

## NExt steps:
## Get rid of integer effect
## Wider range of mean expression values
## Clean up scrtipt
## MOre diff exp genes
## Define specificity and sentitivty script- realized an issue with the first run- the test doesnt pick up DE genes based on FC becasue the "true" expression isnt above the thresholds. 
## Re think the way in which fc values are sampled. Perhaps the "true" de need to be smapled on a FC al 1.5 basis- in other words, they are linked and 95% if the time, the fc values are al 2. 

## - ok this may need tobe rethought. Real de seq data has counts that range from 1, 10s to 100s and 500s. Need a way to simulate this
## Historgram of real count data shows a mixed/uniform that has representatin from 0 to 10^8.

## First, make counts plot look like ones from existing data. Then, figure out parameters needed to make that happen. 
## Generate counts for both sig and non sig based on what this dist should look like
plot(res$log2FoldChange,-log(res$padj,base=10))
abline(h=-log(0.01,base=10),col='red')
abline(v=c(-1.5,1.5),col='red')
hist(res$log2FoldChange)
miRNA$exp.pbv/miRNA$exp.nbv

calculateSensitivitySpecificity<-function(res,miRNA){
  diff.exp.test<-as.character(row.names(res[res$padj<0.01 ,]))#& abs(res$log2FoldChange)>1.5,]))
  ndiff.exp.test<-as.character(row.names(res[res$padj>=0.01 ,]))#| abs(res$log2FoldChange)<=1.5,]))
  diff.exp.true<-as.character(miRNA[miRNA$BV==1,]$miRNA.id )
  diff.exp.false<-as.character(miRNA[miRNA$BV==0,]$miRNA.id)

  true.positives<-sum(diff.exp.test %in% diff.exp.true)
  false.positives<-sum((diff.exp.test %in% diff.exp.false))
  true.negatives<-sum(ndiff.exp.test %in% diff.exp.false)
  false.negatives<-sum((ndiff.exp.test %in% diff.exp.true))
  sensitivity<-true.positives/(true.positives+false.negatives)
  specificity<-true.negatives/(true.negatives+false.positives)
  
  return(data.frame(sensitivity=sensitivity,specificity=specificity))
}
res.combined<-data.frame(miRNA[order(as.character(miRNA$miRNA.id)),],res[order(row.names(res)),])
data.frame(res.combined,de=(abs(res.combined$log2FoldChange)>1.5 & res.combined$padj<0.01))

qbeta(.2,2,4)
hist(rbeta(1000,2,7))


library("parathyroidSE")
data("parathyroidGenesSE")
se<-parathyroidGenesSE
colnames(se)<-colData(se)$run
ddsPara<-DESeqDataSet(se=se,design=~patient + treatment)
colData(ddsPara)$treatment <- factor(colData(ddsPara)$treatment,levels=c("Control","DPN","OHT"))
dds<-DESeq(ddsPara)
res<-results(dds)
x<-abs(res[abs(res$log2FoldChange)>1.5 & !is.na(res$log2FoldChange),]$log2FoldChange)
hist(x)
table(x)
library(MASS)
fitdistr(round(assays(ddsPara)$counts),"Negative Binomial")
warnings()
hist(log(assays(ddsPara)$counts))
hist(rnbinom(1000,mu=2.04134881,size=1/262.24898769))
plot(-10:10,exp(-10:10))
