library(mclust)
library(ggplot2)

fishersTest_DFs <- function(cls_df1,cls_df2) {
  levels1  = levels(as.factor(cls_df1[,1]))
  levels2  = levels(as.factor(cls_df2[,1]))
  #print(levels1)
  #print(levels2)
  fisherMat = matrix(nrow=length(levels1),ncol=length(levels2))
  for (i in 1:length(levels1)) {
    cls_1 = rownames(subset(cls_df1, cls_df1[,1] == levels1[i]))
    #print(cls_1) 
    for (j in 1:length(levels2)) {
      cls_2 = rownames(subset(cls_df2,cls_df2[,1] == levels2[j]))
      #print(cls_2)
      ovl_ij = which(cls_1 %in% cls_2)
      fisherMat[i,j] = length(ovl_ij)
    }
  }
  print(fisherMat)
  fisher_test = fisher.test(fisherMat,alternative="two.sided",workspace=2e6)
  print(paste("P.value=",fisher_test$p.value))
}


plotPCA_syn <- function(mol_dta){
               clus_vec = unlist(lapply(rownames(mol_dta),FUN=function(x){unlist(strsplit(x,"_"))[1]}))
               mol_dta.pca = as.data.frame(prcomp(mol_dta,retx=T)$x)
               mol_dta.pca$cluster = clus_vec
               plt = ggplot(mol_dta.pca,aes(x=PC1,y=PC2,col=as.factor(cluster))) + geom_point(size=4)
               return(plt)
}


getKmeans <- function(molData,k=2){
          cls = kmeans(molData,k)
          cls.df = data.frame(cluster = as.numeric(cls$cluster),row.names=names(cls$cluster))
          return(cls.df)
  
}

getMclust <- function(molData){
          molData.mcl = Mclust(molData)$classification
          cls.df = data.frame(cluster=as.numeric(molData.mcl),row.names=names(molData.mcl))
          return(cls.df)
}


