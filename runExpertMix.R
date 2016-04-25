source("expertMix.R")

clFile="SyntheticResponse.txt"
molFile ="SyntheticMixture.txt"

molDF = read.table(molFile,check.names=F,stringsAsFactors=F)
clDF = read.table(clFile,check.names=F,stringsAsFactors=F)

print(molDF[1:5,1:5])
print(clDF[1:5,])

bic_clsw = parRuns(molDF,clDF)

print(bic_clsw)

