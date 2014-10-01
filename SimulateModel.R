rm(list=ls())
#source("/Users/stevensmith/bin/R_source_functions/ngs.R")
setwd("/Users/stevensmith/bin/Qualifiers/Simulations")
library("reshape")
library(DESeq2)
library(ggplot2)


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


topDESeq<- function(results=res,p.thres=0.01,FC.thres=1.5){
  DE.pval<-results[!is.na(results$padj) & results$padj<p.thres,]
  DE.fc<-results[!is.na(results$log2FoldChange) & abs(results$log2FoldChange)>FC.thres,]
  DE<-results[rownames(results) %in% intersect(rownames(DE.pval),rownames(DE.fc)),]
  DE.counts<-data.frame(Pval= nrow(DE.pval),FC=nrow(DE.fc),Both=nrow(DE))
  DEList<-list(pval=DE.pval,fc=DE.fc,DE=DE,counts=DE.counts)
  return(DEList)
}

#library("parathyroidSE")
#data("parathyroidGenesSE")
#se<-parathyroidGenesSE
#colnames(se)<-colData(se)$run
#ddsPara<-DESeqDataSet(se=se,design=~patient + treatment)
#colData(ddsPara)$treatment <- factor(colData(ddsPara)$treatment,levels=c("Control","DPN","OHT"))
#dds<-DESeq(ddsPara)
#res<-results(dds)
#hist(res$log2FoldChange)

# res[abs(res$log2FoldChange)>1.5 & res$padj<0.01 & !is.na(res$padj),]
# plot(res$log2FoldChange,-log(res$padj,base=10))
# abline(v=c(-1.5,1.5))
# abline(h=-log(0.01,base=10))
# 
# mu.estimates<-assays(dds)$mu
# disp.estimates<-dispersions(dds)
# disp.estimates[disp.estimates>10]
# hist(abs(log(mu.estimates,base=exp(1))))
# plot(mu.estimates[,1]<10,disp.estimates)
# hist(log(disp.estimates,base=10))
# 
# hist(rnbinom(1000,mu =mu.estimates[2,1],size = 1/disp.estimates[2]))
# abline(v=mu.estimates[2,1],lwd=3,col='red')
# hist(rnbinom(1000,mu =mu.estimates[2,2],size = 1/disp.estimates[2]))
# abline(v=mu.estimates[2,2],lwd=3,col='blue')
# hist(rnbinom(1000,mu =mu.estimates[2,3],size = 1/disp.estimates[2]))
# abline(v=mu.estimates[2,3],lwd=3,col='blue')
# 
# library(MASS)
# x<-abs(log(mu.estimates,base=exp(1)))
# x<-x[!is.na(x)]
# x<-x[1:round(.10*length(x))]
# hist(x)
# fitdistr(x,"Negative Binomial")
# nbin(length(x))
# hist(rnbinom(length(x),mu = 4,size = 1/.05))
# 
 mu.intrcept<-exp((2*rbinom(500,1,0.5)-1)*rnbinom(500,mu=4,size=1/0.5))
 hist(log(mu.intrcept,base=exp(1)))
# hist(coef(dds)$treatmentDPN)
# 
 mu.BV<-exp((2*rbinom(500,1,0.5)-1)*rnbinom(500,mu=4,size=1/0.5))
 hist(log(mu.BV,base=exp(1)))
# hist(coef(dds)$treatmentDPN)
# 
# coef<-data.frame(genes=row.names(coef(dds)),b=coef(dds)$treatmentDPN)
# hist(coef$b)
# coef.siggenes<-coef[abs(coef$b>1.5) & !is.na(coef$b),]
# hist(coef.siggenes$b)


## Generate counts using NB distribution and seperate u/alpha pararmeters
M<-length(mu.intrcept)
nsubj<-25
numreps<-6
numsiggenes<-0.1*M
J<-nsubj*numreps
count<-matrix(nrow=M,ncol=J)
BV.effect<-matrix(0,nrow=M)
BV.effect[1:(numsiggenes)]<-1
Subj.state<-rep(rep(c(0,1),each=3),times=25)
fc.min<-1.5

d<-rgamma(M,shape = 0.87 ,rate = 1.36 )
direction<-2*rbinom(M,prob = .5,size = 1)-1
mu.DE<-exp(direction*(d+fc.min))
mu.NDE<-exp(rep(1,times=M))

for(j in 1:J){
for(miRNA in 1:M){
  count[miRNA,j]<-Subj.state[j]*BV.effect[miRNA]*rnbinom(1,mu=mu.DE[miRNA],size=1/0.01)+
    ((1-Subj.state[j])|(1-BV.effect[miRNA]))*rnbinom(1,mu=mu.NDE[miRNA],size=1/0.01)+
    round(abs(rnorm(1,0,1)))
}
}

##similar states should come from the same distributin. Think about this...
warnings()
mu.intrcept[3]
mu.BV[3]
hist(log(as.matrix(count)))
count[1:10,1:10]
head(count)
count.pl<-data.frame(count=log(count,base=exp(1)),BV.effect)
ggplot(count.pl)+geom_histogram(aes(x=count,fill=factor(BV.effect)))
count.pl.bv<-data.frame(rbind(count.pl[1:(0.1*M),],count.pl[(0.1*M+1):((0.1*M+0.1*M)),]))
ggplot(count.pl.bv)+geom_histogram(aes(x=count,fill=factor(BV.effect)))
hist(log(as.matrix(assays(se)$counts)))
count.pl[1,1]

head(count)

count<-data.frame(count)
names(count)<-paste(rep(1:25,each=6),Subj.state,sep=".")
row.names(count)<-seq(1:nrow(count))
count.dds<-DESeqDataSetFromMatrix(countData = count,design = ~BVState,
                       colData = data.frame(ID=names(count),Subj=rep(1:25,each=6),BVState=Subj.state))
count.dds<-DESeq(count.dds)
res<-results(count.dds)
(res<-res[order(res$padj),])
plot(res$log2FoldChange,-log(res$padj,base=10))
(top<-topDESeq(res))
sort(as.numeric(rownames(top$DE)))
sort(as.numeric(rownames(top$pval)))
res[rownames(res)==1,]
count[rownames(count)==1,]



