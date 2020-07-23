mypca<-function(x,zeromean=T,fdiag=T){
  ## do PCA on x %*% t(x) , normalising first if desired
  xN<-mypcanorm(x,zeromean,fdiag)
  tcov<-xN %*% t(xN)
  eigen(tcov)	
}


mypcanorm<-function(x,zeromean=T,fdiag=T){
  ## Normalise for PCA by setting the diagonal to be the average of the rows, then zero meaning by substracting the row mean (so the diagonal ends up as zero)
  if(fdiag){
    diag(x)<-rowSums(x)/(dim(x)[1]-1)
  }
  if(zeromean){
    return(x-rowMeans(x))
  }else{
    return(x)
  }
}

##########################################
#source("FinestructureLibrary.R") # read in the R functions, which also calls the needed packages
library(plyr)
library(RColorBrewer)
load('RData/metadata.RData')
sort(table(metadata$k13Class))

metadata$`Kelch` = metadata$k13Class
metadata$`Kelch`[!metadata$`Kelch` %in% c('R539T','Y493H','WT','C580Y','P553L','I543T')] = 'Other'
metadata$k13colors =
  mapvalues(metadata$Kelch,
            from = unique(metadata$Kelch),
            to=brewer.pal(name = 'Set1', n = length(unique(metadata$Kelch))))

metadata$Plasmepcolors =
  mapvalues(metadata$PLA1,
            from = unique(metadata$PLA1),
            to=brewer.pal(name = 'Dark2', n = 3)[1:2])

chunkfile<-"RData/output_Ne_zeroone.chunkcounts.out" 
dataraw<-as.matrix(read.table(chunkfile,row.names=1,header=T,skip=0)) # read in the pairwise coincidence 
pcares<-mypca(dataraw)
plot(pcares$vectors[,1],pcares$vectors[,2],xlab=paste("PC",1),
     ylab=paste("PC",2),main=paste("PC",1,"vs",2,": Ne= 0.1"), col = metadata$k13colors)


chunkfile<-"RData/output_Ne_one.chunkcounts.out" 
dataraw<-as.matrix(read.table(chunkfile,row.names=1,header=T,skip=0)) # read in the pairwise coincidence 
pcares<-mypca(dataraw)
plot(pcares$vectors[,1],pcares$vectors[,2],xlab=paste("PC",1),
     ylab=paste("PC",2),main=paste("PC",1,"vs",2,": Ne=1"), 
     col = metadata$k13colors)

