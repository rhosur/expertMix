#written by Raghavendra Hosur, Biogen Idec
#           October 2012

#source("/home/rhosur/Projects/src/sparseRegression.R")
library(snow)

#first read in the data and make it consistent



readData <- function(molFile,clVarFile){
#clDir <- "/home/rhosur/Projects/DAC/"
#clDir <- "/home/rhosur/Projects/src/"
allDta <- read.table(molFile,header=TRUE,check.names=F,stringsAsFactors=F,sep="\t")
#allDta <- read.table(paste(clDir,"SyntheticMixture.txt",sep=""),header=TRUE)
#allDta_allProbes <- read.table(paste(clDir,"ACP_SigProbes_t.txt",sep=""),header=TRUE)
#rownames(allDta) <- sub(".","-",rownames(allDta),fixed=TRUE)
#clVars <- read.table(paste(clDir,"SyntheticResponse.txt",sep=""),header=TRUE)
clVars <- read.table(clVarFile,header=TRUE) 
#row.names(clVars) <- clVars$SAMPLE_ID
print(clVars)
print(allDta)
cl_e <- clVars[which(rownames(clVars) %in% rownames(allDta)),]
allDta <- allDta[which(rownames(allDta) %in% rownames(clVars)),]
print(nrow(allDta))
print(nrow(clVars))
#allDta_allProbes <- allDta_allProbes[which(rownames(allDta_allProbes) %in% clVars$SAMPLE_ID),]
clVars_nums <- clVars[rownames(allDta),,drop=TRUE]
clinicVars <- clVars_nums #- mean(clVars$EDSS)
#allDta <- as.data.frame(scale(allDta,center=TRUE,scale=FALSE))
names(clinicVars) <- rownames(allDta)
print(paste("Final Data sets: ",paste(dim(allDta),collapse=","),paste(dim(clinicVars),collapse=",")))
print(clinicVars[1:5]) 
print(allDta[1:5,1:5])
dta_list <- list(allDta,clinicVars)

return(dta_list)
}


##### function to run on slave nodes -- all other functions are within this

mainFT <- function(rank_slave, args_list1) {

          #result_FT = tryCatch(main(rank_slave,args_list1),error=function(e) c(rep(NULL,args_list1[[3]])))
          #return(result_FT)


main <- function(rankSlave,arg_list) {
        source("sparseRegression.R")
        #library(glmnet)
        #log file
        #last value in the args_list is the debug var
        debug = arg_list[[length(arg_list)]]
        if (debug) {
           sink(paste("Par_expertMix.log.",rankSlave,sep=""),append=TRUE)
        }


L2DistDF <- function(dtaFrame){
            dta_mu <- mean(dtaFrame)
            #subtract mean
            dta_norm <- t(dtaFrame)-dta_mu
            dta_norm <- dta_norm*dta_norm
            dta_L2Dist <- apply(dta_norm,2,sum)
            dta_L2Dist <- sqrt(dta_L2Dist)
            return(dta_L2Dist)
}

probGauss <- function(pt,mu,SigmaInv) {
             #figure out if this is going to be a column or row matrix
             ptSc = pt-mu
             #print(pt)
             #print(t(ptSc))
             #print("*****prob gauss")
             #print(data.class(pt))
             #print(data.class(mu))
             #print(data.class(ptSc))
             prb = t(ptSc) %*% SigmaInv 
             #print(prb)
             prb = prb %*% (ptSc)
             #print(paste("Exponent error term:",prb))
             prb = exp(-0.5*prb)
             #print(paste("Gauss prob:",prb))
             #print(pt-mu) 
             #print(SigmaInv)
             #print(prb)
             #print("####") 
             return(prb)
} 

probLinRegression <- function(y,x,beta,sigma){
                     #treating regression as gaussian prob
                     #add an intercept term, because beta contains an intercept (the last coeff)
                     #x$Intercept = 1.0
                     #print("######prob lin regression")
                     numfeats = length(beta)  
                     mu = sum(x * beta[1:numfeats-1]) + 1.0*beta[numfeats]
                     sigma_inv = 1.0/sigma^2
                     #print(paste("x:",x,"beta:",beta))
                     #print(paste("mu: ",mu,"y:",y,"sigma_inv:",sigma_inv))
                     prb <- probGauss(y,mu,sigma_inv)
                     twopiterm = sqrt(2*3.142)
                     constF = 1.0/(twopiterm*sigma)
                     #print(paste("lin regression prb:",prb)) 
                     return(constF*prb)
}


ClusterSigmaInv <- function(clsPts) {
             #given a cluster, this function estimates inverse of the covariance matrix
             covMat <- cov(clsPts) 
             covInv <- solve(covMat)
             return(covInv) 
}   

clusterLL <- function(clsPts,clsMu,clsSigInv,clsSigDet) {
             #given a df with points belonging to a cluster, evaluate the log-likelihood -- this summed over all the clusters (weightd by cls wts) gives the lower bound to log-likelihood of the data
             #remember that clsPts has to be features in rows and samples in columns (same for Mu)
             ll = 0.0
             #print(colnames(clsPts))
             #print(clsMu)
             #print(clsSigInv)
             nFeat = nrow(clsPts)
             clsSize = length(colnames(clsPts))
             for (col in 1:clsSize){
                 #print(col)
                 ll = ll + log(probGauss(clsPts[,col],clsMu,clsSigInv))
             }
             constantTerm = (2*3.142)^(nFeat/2) * sqrt(clsSigDet) 
             return(ll-clsSize*log(constantTerm))
} 


predictY_newSampleX <- function(newPt,listMu,listSigInv,listDetSig,clsW,listBeta) {
                       #this function is to calculate the prediction based on cluster assignment calculated from posterior probability
                       print("Computing prediction for new point")
                       classign = ''
                       pp_max = ''
                       pp_sum = 0.0
                        #if the intercept is already added, dont add intercept
                       if (length(newPt) == length(listBeta[[1]])-1) {
                          newPt_int = as.numeric(c(newPt,1.0))
                       } else {
                          newPt_int = as.numeric(newPt)
                       }

                       yPred = 0.0

                       for (clIdx  in 1:length(listMu)) {
                           #print("Computing postProb: ")
                           ptsc = as.numeric(newPt-listMu[[clIdx]])
                           #print(ptsc)
                           #print(listSigInv[[clIdx]])
                           #print(length(ptsc))
                           #print(dim(listSigInv[[clIdx]]))
                           #print(data.class(ptsc))
                           #print(data.class(listSigInv[[clIdx]]))
                           prb = t(ptsc) %*% listSigInv[[clIdx]]
                           #print(ptsc)
                           #print(prb)
                           pp_curr = probGauss(as.numeric(newPt),listMu[[clIdx]],listSigInv[[clIdx]])
                           pp_curr = pp_curr*clsW[clIdx]/sqrt(listDetSig[[clIdx]])
                           #print(paste("Curr post prob:",pp_curr))
                           pp_sum = pp_sum + pp_curr
                           if (pp_curr > pp_max) {
                              classign = clIdx
                              pp_max = pp_curr
                           }
                           yPred = pp_curr*(newPt_int %*% listBeta[[clIdx]]) + yPred
                       }
                        
                       #if the intercept is already added, dont add intercept
                            
                       return(yPred/pp_sum)
                       #return(classign)
}



posteriorProb <- function(allPts,clsMu,clsSigInv, clsW,respVec,beta,sigma){
                 #same function as DF but multiplied by the cluster weights clsW
                 numData = length(colnames(allPts))
                 #print(paste("PostProb:beta",beta))
                 pp = c(rep(0.0,numData))
                 for (i in 1:numData){
                     pp[i] = probLinRegression(respVec[i],allPts[,i],beta,sigma) * probGauss(allPts[,i],clsMu,clsSigInv) * clsW
                 }
                 #print("Psoterior prob calculation")
                 #print(pp)
                 #return as a numeric vector maybe
                 names(pp) <- colnames(allPts)
                 return(pp)
}   

clusterExpectation <- function(ptFeat,postpr) {
                      #make sure the order of the samples and the post pr are same, I am not checking it here...
                      #this function calculates the expected value of a function over the cluster whose postpr is given
                      #both vectors have to be of size=number of data points
                      ept = 0.0
                      #print(ptFeat)
                      for (i in 1:length(postpr)) {
                          ept = ept + postpr[i]*ptFeat[i]
                      }
                      #normalize
                      ept = ept/sum(postpr) 
                      return(ept)
}


computeNewMean <- function(allPts,postpr) {

                  #allPts has features as rows and columns as samples. 
                  nFeat = nrow(allPts)
                  newMu = c(rep(0.0,nFeat))
                  #print("Computing new mean..")
                  #print(allPts)
                  #print(postpr)
                  for (i in 1:nFeat) {
                      newMu[i] =  clusterExpectation(allPts[i,],postpr) 
                  }
                  #newMu_matmult = allPts %*% postpr
                  #newMu_matmult = newMu_matmult/sum(postpr)
                  #print("NewMu:")
                  #print(newMu)
                  #print("NewMuMatMult:")
                  #print(newMu_matmult)
                  return(newMu)
}


computeNewSigma <- function(allPts,newMu,postpr,covpr,ess){
                   #given the new mu, compute the ML estimate for the cov variance matrix
                   nFeat = length(newMu)
                   nSam = length(colnames(allPts))
                   nu0= 2.0
                   nprime = ess 
                   #print(paste("Computing new sigma:",nFeat,nSam))
                   for (i in 1:nSam) {
                       if (i==1){
                          #print(allPts[,1])
                          #print(newMu)
                          scVec = allPts[,1]-newMu
                          #print(scVec)
                          newCov = (scVec %o% scVec)*postpr[i]
                          #print(newCov)
                       } else {
                         scVec = allPts[,i] - newMu
                         newCov = newCov + (scVec %o% scVec)*postpr[i]
                         #print(paste("Sample number: ",i))
                         #print(newCov)
                       }
                    }       
                    #add the prior (for Jeffreys prior, no initial covpr is added)
                    newCov = newCov + nprime*(covpr)
                    denM = sum(postpr) + nprime 
                    newCov = newCov*(1.0/denM)   
                    return(newCov)
} 

computeNewClusterWt <- function(postpr){
                        #new clstr weights are \sum postpr/\sum_over_all_clusters \sum_over_all_points
                        return(sum(postpr))          
}







computeRegressionMat_B <- function(allpts,postprob) {
                          #make sure the allpts are arranged so that columns are samples, rows are feats
                          # newBeta = B^-1 A
                          #where B(i,j) = <xi*xj>  
                          nFeat = nrow(allpts)
                          B = diag(nFeat)*0.0
                          for (i in 1:nFeat) {
                              for (j in i:nFeat) {
                                  xixj = allpts[i,]*allpts[j,]
                                  B[i,j] = clusterExpectation(xixj,postprob)
                                  B[j,i] = B[i,j]
                              }
                          }
                          return(B)
}     

computeRegressionMat_A <- function(allpts,respVec,postprob){
                          #newBeta = B^-1 A
                          #where A(i,j) = <yixj>
                          #if y is one dimensional, as in our case
                          #A(j) = <yxj>
                          nFeat = nrow(allpts)
                          A = c()
                          for (j in 1:nFeat) {
                              yxj = respVec * allpts[j,]
                              A = c(A,clusterExpectation(yxj,postprob))   
                          }
                          return(A)
}



computeNewRegressionCoeffs <- function(allpts,respVec,postprob) {
                              #newBeta <- B^-1 A
                              B = computeRegressionMat_B(allpts,postprob)
                              B = B + diag(dim(B)[1])*0.000001
                              Binv = solve(B)
                              A = computeRegressionMat_A(allpts,respVec,postprob) 
                              newBeta = Binv %*% A
                              #print(paste("newBeta:",newBeta))
                              return(newBeta)
}


computeNewRegressionVariance <- function(allpts,respVec,betaVec,postprob) {

                             # sigma^2 = \sum_{i} ppi*(yi-xi*beta)^2 / sum ppi
                             nsam = length(colnames(allpts))
                             sqResiduals = 0.0
                             for (i in 1:nsam) {
                                 yhat = sum(betaVec * allpts[,i])
                                 #print(paste("betaVec:",betaVec,"x:",allpts[,i]))
                                 #print(paste("yhat:",yhat))
                                 sqResiduals = sqResiduals + postprob[i] * (respVec[i]-yhat)^2
                             }      
                             regSig = sqrt(sqResiduals/sum(postprob))
                             #print(paste("newregsig:",regSig))
                             return(regSig)
} 




computeNewEffSampleSize <- function(lambda,detSig,SigInv,priorCov) {
               #lambda is the parameter defining the poisson prior on effective sample size
               #log(n') = -0.5*log(detSig) - 0.5*trace(S*SigInv) + log(lambda)
               expFac = sum(diag(priorCov %*% SigInv))
               print(paste("lambda:",lambda,"DetSig:",detSig," expFac:",expFac))
               lognprime = -0.5*log(detSig) - 0.5*expFac + log(lambda)
               return(exp(lognprime))
}

logprobPoisson <- function(counts,lambda) {
                  #use stirling's approximation to avoid overflow
                  lp = counts*log(lambda) - lambda - counts*log(counts) + counts
                  print(paste("poisson prob: ",lp))
                  return(lp)
}


Debug <- function(df,strToWrite) {
         str2w <- paste("Debug: ",strToWrite)
         print(paste(str2w,df))        
}

rowNormalize <- function(pp) {

                normpp = pp
                cols = length(pp)
                ndat = length(pp[[1]])
                for (i in 1:ndat) {
                    rowSum = 0.0
                    for (j in 1:cols){
                        rowSum = rowSum + pp[[j]][i]
                    }
                    for (j in 1:cols) {
                        normpp[[j]][i] = pp[[j]][i]/rowSum
                    }
                }
                #print("Normalizing post.prob")
                #print(pp)
                #print(normpp)
                return(normpp)
}





EM <- function(listOfMu,listOfSig,initClsWts,MAX_ITER,listcovPrior,effSampleSizes,responses,listOfBeta,initRegSig,initAlpha,writeFit=F) {
   
      currll = 0.0
      numClusters = length(listOfMu)
      nFeat = length(listOfMu[[1]])
      numData = nrow(allDta)
      print("Running EM algorithm..")
      Debug(nFeat," Number of features")
      Debug(numClusters, "Number of clusters")
      Debug(numData, "Number of data points")
      Debug(length(responses),"Number of responses") 


      ###testData
      #testDta <- read.table('/home/rhosur/Projects/DAC/Dacl_1801merged_SNPPCs_NoOutliers.txt',header=TRUE)
      #numTest <- nrow(testDta)  
      ####



      currListOfSigInv = list()
      currListOfDetSig = list()

      for (i in 1:numClusters) {
          #to prevent poor condition number
          sigInv = solve(listOfSig[[i]]+diag(nFeat)*0.000001)
          currListOfSigInv[[i]]= sigInv
          #print(currListOfSigInv[[i]])
          detSig = det(listOfSig[[i]])
          currListOfDetSig[[i]]= detSig
          #print(detSig)
      } 
      currListOfMu = listOfMu
      #currListOfSig = listOfSig
      currClsWts = initClsWts
      currEffSamp = effSampleSizes
      currlistOfBeta = listOfBeta
      currlistOfRegSig = initRegSig
      currlistOfAlpha = initAlpha
      #compute loglikelihood 
      currll = LLK(t(allDta),currListOfMu,currListOfSigInv,currListOfDetSig,currClsWts,listcovPrior,effSampleSizes,responses,currlistOfBeta,currlistOfRegSig,currlistOfAlpha)


      hasConverged = FALSE
      iter = 0
      while ( !hasConverged && (iter <= MAX_ITER) ) {
            iter = iter + 1
            ######Estep###### 
            currPostProb = list()
            for (clIdx in 1:numClusters) {
                pp = posteriorProb(t(allDta),currListOfMu[[clIdx]],currListOfSigInv[[clIdx]],currClsWts[clIdx],responses,currlistOfBeta[[clIdx]],currlistOfRegSig[[clIdx]])
                #store the pp
                 
                constFac = (2*3.142)^(nFeat/2) * sqrt(currListOfDetSig[[clIdx]])
                #print("PostProb before const:")
                #print(pp)
                #print(constFac)
                #post prob only makes sense if you have more than 1 cluster
                if (numClusters == 1) {
                   currPostProb[[clIdx]] = c(rep(1.0,numData))
                } else {
                  currPostProb[[clIdx]]= pp/constFac
                }
                #check to see if the post prob make sense, i.e. the cluster is not empty
                if (!is.finite(currPostProb[[clIdx]]) || max(currPostProb[[clIdx]]) < 1e-30){
                   print(paste("Warning!! posterior probabilities very low for a cluster -- Exiting! ",max(currPostProb[[clIdx]])))  
                   hasConverged = TRUE
                } 
            }

            if (hasConverged) {
               #give a high negative number so that this is not favoured
               newll = -10000
               break
            }
            #print(currPostProb)
            currPostProb <- rowNormalize(currPostProb)
            currPostProb_str <- paste(unlist(lapply(currPostProb,FUN=function(x){paste(x,collapse=",")})),collapse="\t")
            print(paste("CURRPOSTPROB:",currPostProb_str))




            #########Mstep######
            newll = 0.0
            for (clIdx in 1:numClusters) {
                currClsWts[clIdx] = computeNewClusterWt(currPostProb[[clIdx]])
            }
            clsWtSum = sum(currClsWts)
             for (clIdx in 1:numClusters) {
                currClsWts[clIdx] = currClsWts[clIdx]/clsWtSum
            } 
            
            currSBL_llk <- list()
            for (clIdx in 1:numClusters) {
                #first update clustering params
                currListOfMu[[clIdx]] = computeNewMean(t(allDta),currPostProb[[clIdx]])
                currSig = computeNewSigma(t(allDta),currListOfMu[[clIdx]],currPostProb[[clIdx]],listcovPrior[[clIdx]],currEffSamp[clIdx])
                #add stuff to the diagonal for numerical stability
                currSig = currSig + diag(nFeat)*0.00001
                #print(currSig)
                currListOfSigInv[[clIdx]] = solve(currSig)
                currListOfDetSig[[clIdx]] = det(currSig)

                #estimate regression params using sparse bayesian regression, alpha is the prior covariance on each regression coefficient

                optParams_SBL = mainSBL(allDta,responses,currPostProb[[clIdx]],currlistOfBeta[[clIdx]],currlistOfRegSig[[clIdx]],currlistOfAlpha[[clIdx]])
                currlistOfBeta[[clIdx]] = optParams_SBL[[1]]
                print(paste("CURRLISTOFBETA:",paste(currlistOfBeta[[clIdx]],collapse=",")))
                currlistOfRegSig[[clIdx]] = optParams_SBL[[2]]
                print(paste("CURRLISTOFREGSIG:",paste(currlistOfRegSig[[clIdx]],collapse=",")))
                currlistOfAlpha[[clIdx]] = optParams_SBL[[3]]
                print(paste("CURRLISTOFALPHA:",paste(currlistOfAlpha[[clIdx]],collapse=",")))
                currSBL_llk[[clIdx]] = optParams_SBL[[4]]
                #currlistOfBeta[[clIdx]] = computeNewRegressionCoeffs(t(allDta),responses,currPostProb[[clIdx]])
                #currlistOfRegSig[[clIdx]]  = computeNewRegressionVariance(t(allDta),responses,currlistOfBeta[[clIdx]],currPostProb[[clIdx]]) 
                #print(paste("CurrListOfBeta:",currlistOfBeta[[clIdx]]))
            }
            newll = LLK(t(allDta),currListOfMu,currListOfSigInv,currListOfDetSig,currClsWts,listcovPrior,currEffSamp,responses,currlistOfBeta,currlistOfRegSig,currlistOfAlpha) 
          
            predFile = paste("./ClusterOut/Tys_SNPs_DACPredictions.txt.",rankSlave,sep="") 
            #exit conditions
            if (sum(newll) < sum(currll)) {
               print("Warning: decrease in log likelihood!!")
               hasConverged = TRUE
               print("Mu:")
               print(currListOfMu)
               print("Beta:")
               print(currlistOfBeta)
               #print("SigmaInv:")
               #print(currListOfSigInv)

            } else {
            if ((sum(newll)-sum(currll)) < 0.0000001) {
               hasConverged = TRUE
               print("LLK Converged!")
               print("Mu:")
               print(currListOfMu)
               #print(currListOfSigInv)
               print(currListOfDetSig)
               #print(currClsWts)
               print("Beta:")
               print(currlistOfBeta)
               #print(currlistOfRegSig)  
               #print("Final Posterior probs: ")
               #print(currPostProb)
               ##print out predictions
               if (numClusters == 2) {
                  #for (i in 1:numTest){
                    #ptSum = 0.0
                    #ypred = predictY_newSampleX(testDta[i,],currListOfMu,currListOfSigInv,currListOfDetSig,currClsWts,currlistOfBeta)
                    #print(paste(rownames(allDta)[i],responses[i],ypred))
                    #write("\n",predFile,append=TRUE,ncolumns=2000)
                    #write(paste(rownames(testDta)[i],ypred),predFile,append=TRUE,ncolumns=2000)
                    #}
                 }
              }
            }
            # if weight of any one cluster becomes greater than 0.99, then also say converged 
            for (clw in currClsWts) {
                if (clw > 0.99) {
                   #hasConverged = TRUE
                   print("Warning! Cluster Weight very high!!")
                   #print(currListOfMu)
                   #print(currListOfSigInv)
                   #print(currListOfDetSig)
                   #print(currClsWts)
                }
            } 
            print(paste("Iter: ",iter))
            #print(paste("prevLL: ", sum(currll), " newLL: ", sum(newll)))
            currll <- newll
      }
      print(paste("Exiting EM -- ",iter,"iterations"))
      print(currListOfMu)
      mmbList <- getMembership(currPostProb)
      if (writeFit){
         print("Got here -- write1")
         fout="expertMix_out.txt"
         write("######",fout,append=T)
         write(paste(as.POSIXlt(Sys.time()),collapse="_"),fout,append=T)
         for (w_k in 1:numClusters) {
             mmbList.cls = mmbList[[w_k]]
             regVec.cls = currlistOfBeta[[w_k]]
             mmbList.cls.str = paste(mmbList.cls,collapse=",")
             regVec.cls.str = paste(regVec.cls,collapse=",")
             str_to_write1 = paste("CLUSTERMEMBERSHIP(",numClusters,":",sum(currll),"):",mmbList.cls.str)
             write(str_to_write1,fout,append=T)
             str_to_write2 = paste("REGCOEFFS(",numClusters,":", sum(currll),"):",regVec.cls.str)
             write(str_to_write2,fout,append=T)
         }
      }  
      #}
      
      return(list("LLK_tot"=sum(currll),"membership"=mmbList,"reg_coeffs"=currlistOfBeta))
}


logprobWishart <- function(precMatrix,priorInvCov,detSig,ess) {
                  #this is an alternative formulation, taken from Jaakkola's slides. n' is a fuzzy parameter -- "effective sample size" (number of training samples you would have to see for
                  #the prior and the empirical to have the same effect.  
                  #log prob = n'/2*log(det(precMatrix)) - 0.5*n'*tr(priorInvCov*precMatrix)
                  #nprime = nu0+p+1 = 2p+4
                  nF = dim(precMatrix)[1]
                  nu0 = nF+1
                  nprime = ess
                  lp = -1.0*(nprime)/2.0 * log(detSig) - 0.5*nprime*sum(diag(priorInvCov %*% precMatrix))           
                  return(lp)
}


logprobRegressionWtPrior <- function(regWeights,alpha) {
                           # this is to calculate the prior probability of the regression coefficients. We assume (in SBL) that the coefficients are normally distributed with mean 
                           # zero and inv variance given by alpha
                           probWt = 1.0
                           #iterate over all dimensions
                           #print("*******Regression weights prior")
                           for (i in 1:length(regWeights)) {
                               currDprob = sqrt(alpha[i]/(2*3.142)) * exp(-0.5*alpha[i]*regWeights[i]*regWeights[i])
                               probWt = probWt * currDprob
                               #print(paste("alpha_i:",alpha[i],"regWeight:",regWeights[i],"currDprob:",currDprob))
                           }
                           return(log(probWt))
}


logprobGamma <- function(x,sh,rt) {
                #return the prbability of gamma distirbution - Gamma(x|a,b) ~ gammaF(a)^-1 * b^a * x^(a-1) * exp(-b*x)  
                #print(paste("X:",x))
                lp= pgamma(x,shape=sh,rate=rt,log=TRUE)
                #print("LogP:")
                return(lp)
}



logprobJeffreys <- function(detSig,numFeat) {
                   lp = (-1.0*(numFeat+2)/2.0) * log(detSig)
                   return(lp)
} 


LLK <- function(allDta, listOfMu, listOfSigI, listOfDetSig, listOfClsWts,listOfCovprior,effSizes,respVec,listofbeta,listofregsig,listofAlpha) {
       #llk of each point is a weighted sum of the llk w.r.t each cluster
       llk = 0.0;
       nSam = length(colnames(allDta))
       numCl = length(listOfMu)
       numDims = length(listOfMu[[1]])
       #print(paste("LLK:",nSam,numCl))
       llk_clust = 0.0
       llk_reg = 0.0
       for (i in 1:nSam){
           ptSum = 0.0
           for (clIdx in 1:numCl){
               #print(listOfDetSig[[clIdx]])
               constF = sqrt(listOfDetSig[[clIdx]])
               twoPiTerm = (2*3.142)^(numDims/2) 
               constF = 1.0/(twoPiTerm*constF)
               ptClust = constF*probGauss(allDta[,i],listOfMu[[clIdx]],listOfSigI[[clIdx]])
               ptReg = probLinRegression(respVec[i],allDta[,i],listofbeta[[clIdx]],listofregsig[clIdx])
               #ptProb = probGauss(allDta[,i],listOfMu[[clIdx]],listOfSigI[[clIdx]]) * probLinRegression(respVec[i],allDta[,i],listofbeta[[clIdx]],listofregsig[clIdx]) 
               ptProb_const = ptClust * ptReg  
               ptSum = ptSum + listOfClsWts[clIdx] * ptProb_const
               #if (constF == 1.0){
               print(paste("NUMBEROFCLUSTERS:",numCl,"ptSum:",ptSum,"ptProb:",ptProb_const,"constF:",constF,"ptClust:",ptClust,"ptReg:",ptReg))
               #}
            }
           #jtllk is the
            
           llk = llk + log(ptSum)
       }
       #add logPrior
       logprior = 0.0
       for(clIdx in 1:numCl) {
           lp_SigmaInvClus =  logprobWishart(listOfSigI[[clIdx]],listOfCovprior[[clIdx]],listOfDetSig[[clIdx]],effSizes[clIdx]) 
           lp_RegWts = logprobRegressionWtPrior(listofbeta[[clIdx]],listofAlpha[[clIdx]])
           ###for the hyperpriors, a=b=c=d=1e-5 (if you change this, change it in SBL as well)
           lp_Alpha = sum(logprobGamma(listofAlpha[[clIdx]],1e-5,1e-5)) 
           ### regSig^(-2) is distributed as gamma
           lp_RegSig = logprobGamma(listofregsig[clIdx]^(-2),1e-5,1e-5) 
           #print(paste("RegWeightsPrior:",lp_RegWts,"RegSigPrior:",lp_RegSig,"AlphaPrior:",lp_Alpha,"SigmaInvPrior:",lp_SigmaInvClus)) 
           logprior = logprior + lp_SigmaInvClus + lp_RegWts+ lp_Alpha+ lp_RegSig
       }

       print(paste("LLK: ",llk," logprior: ",logprior))
       return(c(llk,logprior))
}

getMembership <- function(listOfPostProbs) {
                 numClusters = length(listOfPostProbs)
                 numDta = length(listOfPostProbs[[1]])
                 clustMembers <- list()
                 for (k in 1:numClusters){
                     clustMembers[[k]]=NA
                 }
                 
                 for (i in 1:numDta){
                     maxj = 1
                     maxpostProb = 0.0
                     for (j in 1:numClusters){
                         pp = listOfPostProbs[[j]][i]
                         if (is.finite(pp)){
                            if (pp > maxpostProb) {
                               maxpostProb = pp
                               maxj = j 
                            }
                         }else {
                            print(paste("PostProb not finite:",i))
                         }
                         
                     }
                     #print(maxj)
                     clustMembers[[maxj]]  = c(clustMembers[[maxj]],names(listOfPostProbs[[j]])[i])
                     #print(paste(names(listOfPostProbs[[j]])[i],maxj))
                 }
                 print(clustMembers) 
                 return(clustMembers)                         
}             



heuristicStableClusters <- function() {
clIdx = 0
par(mfrow=c(5,4))
for (clf in clFiles){
    dtaF <- read.table(paste(clDir,clf,sep=""),header=TRUE)
    currdtaF <- dtaF
    isStable <- FALSE
    iter <- 0
    print(paste("Cluster",clIdx))
    print(paste("Initial cluster size: ", nrow(dtaF)))
    while(!isStable) {
         dta_L2Dist <- L2DistDF(currdtaF)
         mean_d <- mean(dta_L2Dist)
         std_d <- sd(dta_L2Dist)
         robust <- dta_L2Dist[which(dta_L2Dist <= 2*std_d + mean_d)]
         hist(dta_L2Dist,main=paste("Cluster",clIdx,"Iter: ",iter))
         #print(robust)
         currdtaF <- currdtaF[names(robust),,drop=TRUE]
         iter = iter+1
         if ( (length(robust) < 0.70* nrow(dtaF)) || (length(robust) == length(dta_L2Dist))) {
            isStable <- TRUE
            print(paste("Final cluster size: ",length(dta_L2Dist)))
         } else {
                remSamples = length(dta_L2Dist) - length(robust)
                print(paste("Removed",remSamples,"iter:",iter))
         }
   }

   #once you get the stable clusters, get the likelihood
   currMu = mean(currdtaF)
   currSigI = ClusterSigmaInv(currdtaF)

   ##### posterior on all data
   llk = posteriorProb(t(allDta),currMu,currSigI,nrow(currdtaF)/110.0)
   newMu = computeNewMean(t(allDta),llk)
   print("CurrMu: ")
   print(currMu)
   print("NewMu:")
   print(newMu)
   print(newMu-currMu)
   newSigI = computeNewSigmaInv(t(allDta),newMu,llk)
   print("CurrSigI:")
   print(currSigI)
   print("NewSigI:")
   print(newSigI)
   #print(paste("Cluster",clIdx," LLK: ",llk)) 
   clIdx = clIdx + 1
}  
}






#main <- function(arg_list) {
#heuristicStableClusters()
     #will need to initialize the EM
     print(arg_list)
     print("Initializing EM..")
     allDta <- arg_list[[1]]
     clinicVars <- arg_list[[2]]
     print(allDta[1:5,1:2])
     print(clinicVars[1:5,]) 

     ######REMOVE THIS######
     ###randomize the msss values
     #oldOrder = names(clinicVars)
     #newOrder = sample(oldOrder,length(oldOrder),replace=FALSE)
     #names(clinicVars) = newOrder
     ###    
 
     numData <- nrow(allDta)
     #test for Kmax clusters
     Kmax = arg_list[[3]]
     print("Got here")
     #allDta$Intercept = c(rep(1.0,numData))
     #priorCov = solve(priorCov)
     BIC <- list()
     sparsity <- list()
     membership_all <- list()
     regression_coeffs <- list()
     BIC[[Kmax+1]] = 0.0
     listofBICs <- c()
     numFeats = length(colnames(allDta))
     if (numFeats==1){
        priorCov = cov(as.matrix(allDta))
     }else {
        print(numFeats)
        priorCov = 1.01*cov(allDta)
     }

     #repeat the whole procedure 10 times
     #for (rerun in 1:10) {
     #print(paste("Rerun: ", rerun))
      for (k in 1:Kmax) {


         #first randomize
         positions <- sample(numData) 
         listOfMu <- list()
         listOfSig <- list()
         listOfEffSam <- c()
         lambda = numData/k
         initClsWts = c()
         listofBeta = list()
         listofregsig = c()
         listofAlpha = list()
         everyThingOK = TRUE

         clIdx = 1
         sizeOfCls = round(numData/k)
         ##initialize via kmeans
         if (k>1){
                #allDta_initcls <- kmeans(allDta,centers=k,nstart=5)
                cluster.init=sample(seq(1,k,1),nrow(allDta),replace=T)
                names(cluster.init) = rownames(allDta)
                allDta_initcls = list("cluster"=cluster.init) 
                print(allDta_initcls)
         }
         ####initializing####
         for (clf in 1:k){
             lowerLim = (clf-1)*sizeOfCls+1
             upperLim = min(numData,(clf)*sizeOfCls)
             #print(paste(lowerLim,upperLim)) 
             #dtaF <- allDta[positions[lowerLim:upperLim],]
             if (k>1){
                  dtaF <- allDta[which(allDta_initcls$cluster==clf),]
             } else {
                  dtaF <- allDta
             }
             print(paste("Size of initial Cluster: ",nrow(dtaF)))
             currmu <- unlist(lapply(dtaF,mean))
             print(paste("Dim of dtaF:",dim(dtaF)))
             print(dtaF[1:5,1:2])
             print(paste("currmu:",currmu)) 
             currcov <- cov(dtaF)
             listOfMu[[clIdx]] = currmu
             if (nrow(dtaF) > numFeats){
                listOfSig[[clIdx]]= currcov
             } else {
                
                #lazy to deal with very small clusters...so just skip over that clustering
                #everyThingOK = FALSE
                #break
                listOfSig[[clIdx]] = currcov
             }
            
             listOfEffSam <- c(listOfEffSam,1.0) 
             initClsWts <- c(initClsWts,1.0/k)

             #initialize regression params
             #cl_e <- clVars[which(rownames(clVars) %in% rownames(dtaF)),]
             clFeat <- dtaF[which(rownames(dtaF) %in% rownames(clinicVars)),]
             cl_e <- clinicVars[rownames(clFeat),]
             #add intercept, so the number of features in greater by 1
             print(clFeat[1:5,1:2])
             #clFeat$Intercept = c(rep(1.0,nrow(dtaF)))
             print("Building initial model...")
             if (nrow(dtaF) > 2*(numFeats+1)){
                 #model <- lm(cl_e~0+.,data=clFeat)
                 #go for L2 regularized solution 
                 model <- glmnet(as.matrix(clFeat),cl_e,alpha=0)
                 # need a penalty -- using an arbitrary one for now
                 #need to put intercept at the back. 
                 model$coefficients <- coefficients(model,s=0.01)[2:(numFeats+1),1]
                 model$coefficients[numFeats+1] = coefficients(model,s=0.01)[1,1]
                 print("Model formula for general..")
                 print(model$call) 
                 listofBeta[[clIdx]] = model$coefficients
                 print("Initial Model:")
                 print(model$coefficients)
                 listofregsig = c(listofregsig,0.5*var(cl_e))
                 alpha_cls=c(rep(1000*max(abs(model$coefficients)),numFeats),1)
                 listofAlpha[[clIdx]] = alpha_cls
                 #print(alpha_cls)

             } else {
                    model <- glmnet(as.matrix(clFeat),cl_e,alpha=0)
                    # need a penalty -- using an arbitrary one for now
                    #need to put intercept at the back. 
                    model$coefficients <- coefficients(model,s=0.01)[2:(numFeats+1),1]
                    model$coefficients[numFeats+1] = coefficients(model,s=0.01)[1,1]
                    #fit only one feature, put all others to zero
                    #model_coeff = c(rep(1e-4,numFeats+1))
                    #select feature to include randomly
                    #feat_int <- sample(seq(1:numFeats),1)
                    #feat_incl = colnames(clFeat)[feat_int]
                    #print(feat_incl)
                    #fmla <- as.formula(paste("cl_e~0+",feat_incl,sep="")) 
                    #print("Model formula for fitting only one feature...")
                    #print(fmla)
                    #model <- lm(fmla,data=clFeat) 
                    #model_coeff[feat_int] = model$coefficients
                    #print(clFeat)
                    #print(cl_e) 
                    #model_coeff[numFeats+1] = 1.0
                    print("Initial Model:")
                    #print(model$coefficients)
                    listofBeta[[clIdx]] = model$coefficients
                    if (length(cl_e) > 1) {
                       listofregsig = c(listofregsig,0.5*var(cl_e))
                    } else {
                        listofregsig = c(listofregsig,0.1)
                    }
                    alpha_cls=c(rep(1000*max(abs(model$coefficients)),numFeats),1)
                    listofAlpha[[clIdx]] = alpha_cls
             } 

             clIdx = clIdx + 1
         }
         #print(listOfMu)

         if (!everyThingOK) {
            listofBICs <- c(listofBICs,'NA')
            break
         }
         #### list of prior cov matrices is the same as list of sigmas
         listOfPriorCov = listOfSig
         #calculate BIC based on log-likelihood
         
         em_results = EM(listOfMu,listOfSig,initClsWts,1000,listOfPriorCov,listOfEffSam,clinicVars[,1],listofBeta,listofregsig,listofAlpha,writeFit=F)
         loglik = em_results[["LLK_tot"]]
         #number of params = k*dim_of_mu + k*dim_of_mu^2(Sigma) + (k-1) cluster weghts + (dim_of_mu+1)*k (regression wts, std_dev) + (dim_of_mu+1)*k (alpha)
         penalty =  (k*numFeats+k*(numFeats*numFeats) + k-1 + 2*k*(numFeats+1))*log(numData)
         print(paste("LLK(data)",loglik," Penalty: ",penalty))   
         listofBICs <- c(listofBICs,-2*loglik + penalty)
         sp_beta = unlist(lapply(listofBeta,FUN=function(x){(length(which(x > 0.000005)))/numFeats}))
         sparsity[[k]] = mean(sp_beta)
         #print(sparsity)
         #if (rerun == 1){
         BIC[[k]] = -2*loglik + penalty
         if (is.infinite(BIC[[k]])){
            BIC[[k]] = NaN
         }
         membership_all[[k]] = em_results[["membership"]]
         regression_coeffs[[k]] = em_results[["reg_coeffs"]]
         #} else {
            #BIC[[k]] = c(BIC[[k]],-2*loglik + penalty) 
         #}
     
    }
    #print(listofBICs)
    #print(BIC[1:Kmax])
    #svg("ACP_BICs_ChaussInit.svg")
    #boxplot(BIC[1:Kmax],main="K-means Initialization",xlab="Number of Clusters",ylab="BIC")
    #dev.off()
    return(list("BIC"=BIC,"membership"=membership_all,"reg_coeffs"=regression_coeffs))
}#main

    print("Got here:Flag 1")
    result_FT = tryCatch(main(rank_slave,args_list1),error=function(e) print(e$message)) #c(rep(NULL,args_list1[[3]])))
    print("Got here:Flag2")
    return(result_FT)


}#mainFT
#run the main function


#############################################

parRuns <- function(molDF,clinDF,reruns=5,maxClust=5,plotFile="parEM_boxplot",debug=F) {

           numReruns  = reruns
           Kmax = maxClust
           dta_args <- list(molDF,clinDF,Kmax,debug)
           numSlaves = 5
           bic_clsw <- list()
           for (i in 1:Kmax){
               bic_clsw[[i]] = c(rep(NA,numSlaves*numReruns))
           }
           
           min_bic = 99999999999
           best_membership = list()
           best_reg_coeffs = list()
           
           for (rnIdx in 1:numReruns) {
               print("Got here -1")  
               jpClusters <- makeCluster(numSlaves,type="SOCK")
           
               #print(dta_args[[2]]) 
               print("Got here 0")
               parEM_all <- clusterApply(jpClusters,1:numSlaves,mainFT,dta_args) 
               #print(parEM_all)
               ##plotting function
               if (length(parEM_all) > 1) {
                  parEM = lapply(parEM_all,FUN=function(x){if (length(names(x))>1){return(x$BIC)}else{return(NULL)}})
                  parEM_mmb = lapply(parEM_all,FUN=function(x){if (length(names(x))>1){return(x$membership)}else{return(NULL)}})
                  parEM_reg = lapply(parEM_all,FUN=function(x){if (length(names(x))>1){return(x$reg_coeffs)}else{return(NULL)}})
                  #print(parEM_sparsity)
                  print(parEM)
               
                  for (j in 1:numSlaves){
                      if (length(parEM[[j]]) > 0) {
                        for (i in 1:Kmax) {
                           if (!is.nan(parEM[[j]][[i]])) {
                              if (!is.na(parEM[[j]][[i]])) {
                                 bic_clsw[[i]][numSlaves*(rnIdx-1)+j] = parEM[[j]][[i]]
                                 if (parEM[[j]][[i]] < min_bic){
                                    min_bic = parEM[[j]][[i]]
                                    best_reg_coeffs = parEM_reg[[j]][[i]]
                                    best_membership = parEM_mmb[[j]][[i]]
                                 }
                              }   
                           }   
                        }
                      }
                  }
               }else{
                  print("Something went wrong -- try re-running with Debug option")
               }    
               stopCluster(jpClusters)
           } 
           #print(bic_clsw)
           #remove NaNs
           bic_clsw_clean = lapply(bic_clsw,FUN=function(x){x[!is.na(x)]})
           #print("BIC:")
           #print(bic_clsw_clean)
           print(paste("minBIC:",min_bic))
           print("Membership:")
           print(na.omit(best_membership))
           print("RegCoeffs:")
           print(best_reg_coeffs)
           pdf(paste(plotFile,".pdf",sep=""))
           boxplot(bic_clsw_clean)
           #stopCluster(jpClusters) 
           dev.off()
           
           return(bic_clsw_clean)
           
}



#args <- commandArgs(trailingOnly=TRUE)
#if (length(args) <= 2) {
   #print(length(args))
   #print(paste(args,collapse=" "))
   #print("Usage: ")
   #print("wtFit.R  MolecularDataFile ClinicalDataFile plotFile")
#}else{ 
   #if (length(args)==3) {
   #parRuns(args[1],args[2],args[3])
   #} else {
   #parRuns(args[1],args[2])
   #}
#}
