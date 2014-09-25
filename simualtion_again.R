library(DESeq2)
(counts.simulation<-matrix(floor(abs(1*rnorm(6*20,5,1))),nrow=20,ncol=6))
(counts.simulation[1:5,1:3]<-matrix(floor(abs(1*rnorm(5*3,50,1))),nrow=5,ncol=3))
(rownames(counts.simulation)<-1:20)
(colData.simulation<-data.frame(treatments=c(rep("untreated",times=3),rep("treated",times=3))))
(simulation<-DESeqDataSetFromMatrix(countData = counts.simulation,colData=colData.simulation,design=~treatments))
simulation<-estimateSizeFactors(simulation)
simulation <- estimateDispersionsGeneEst(simulation)
dispersions(simulation) <- mcols(simulation)$dispGeneEst

simulation<-DESeq(simulation)
