#Written by: Raghavendra Hosur, Biogen Idec, October 2012
#sparse regression as a bayesian learning procedure (Tipping et al.,2001)
# A = diag(numFeats) * alpha  -- we assume equal variance (1/alpha) for all the regression coefficients 
# X = design matrix
# Y = response vector (assume scalar for the moment)

# In the E-step, compute the posterior probability of the regression weights as follows:
#               P(W|Y,alpha,regSig) = N(mu,Sig)

# where         mu = regSig^-2 * Sig %*% X^t %*% Y
#               SigInv = regSig^-2 X^t %*% X + A 
#               Sig = inv(SigInv)

# In the M-step we maximize the expected value (w.r.t to the above distribution) of the complete log likelihood, and find alpha and regSig
#               alpha_d = 1/(Sig_dd + mu_d*mu_d)
#               if using single alpha
#               alpha = d/\sum_d(mu_d^2 + Sigma_dd)    
#               regSig^2 = (||Y-XW||^2 + regSig_prev^2 \sum gamma_d) / N
   

library(glmnet)

computeNewSigma_w <- function(dtaPts,regSig,weightPrior) {
                  #we will assume an equal prior variance for all the regression weights
                  nmFeats = length(colnames(dtaPts))
                  XtX = t(as.matrix(dtaPts)) %*% as.matrix(dtaPts) 
                  A = diag(nmFeats) * weightPrior
                  #print("A:")
                  #print(A)
                  SigInv = regSig^(-2) * XtX + A
                  Sig = solve(SigInv)
                  return(Sig)
}


computeNewMean_w <- function(dtaPts,respVect,Sigma,regSig){
                  XtY = t(as.matrix(dtaPts)) %*% respVect
                  mu = regSig^(-2) * Sigma %*% XtY                  
                  return(mu) 
}

updatePriorVariance <- function(Sigma,mu){
                       # alpha[i] is distributed as Gamma(a,b) with a = b = 1e-5 (uninformative)
                       a = 1e-5
                       b = 1e-5
                       diagEl <- diag(Sigma)
                       nfeats <- length(diagEl)
                       newvar = c()
                       for (i in 1:nfeats){
                           newvar[i] =  (1+2*a)/(mu[i]*mu[i] + diagEl[i]+2*b)
                       }  
                       
                       return(newvar)
}


updatePriorRegSig <- function(dtaPts,respVect,mu,Sigma,regSig,alpha) {
                    

                     ####### 1/regSig^2 is distributed as Gamma(c,d) with c=d=1e-5   (uninformative)
                     c = 1e-5
                     d = 1e-5 
                     errorT = respVect - as.matrix(dtaPts) %*% mu 
                     errorSq = sum(errorT * errorT)
                     sum_gamma = sum(1-alpha*diag(Sigma))
                     new_regSig2 = (errorSq + (regSig^2)*sum_gamma + 2*d)/ (nrow(dtaPts)+2*c)
                     return(sqrt(new_regSig2))
}


LogLikelihood_type2 <- function(dtaPts,respVec,regSig,Sigma,mu,alpha){

                       #evaluating equation #36 in Tipping and Faul, 2003
                       beta = regSig^(-2)
                       A = diag(length(colnames(dtaPts))) * alpha
                       term1 = -1.0*log(det(Sigma)) - nrow(dtaPts) * log(beta) - log(det(A))
                       error2 = respVec - as.matrix(dtaPts) %*% mu
                       #take dot product with Y
                       term2 = beta * t(respVec) %*% error2

                       llk = -0.5*(term1+term2)  
                       return(llk)
}



sparseBayesianRegression <- function(dta,respVec,initW,initAlpha,initRegSig,MAX_ITER) {
    
                           iter = 0
                           nFeats <- length(colnames(dta))
                           currRegSig = initRegSig
                           currW = initW
                           currAlpha = initAlpha 
                           currA = diag(nFeats)*currAlpha
                           hasConverged = FALSE
                           #list to store the optimal params after convergence -- (currMu, currRegSig, currAlpha)
                           optParams = list()

                           while (!hasConverged && (iter <= MAX_ITER)) {
            
                           ## given initial values ,first compute the E-step, i.e. get new mus and sigmas
                           ### E-step 
                                    currSigma = computeNewSigma_w(dta,currRegSig,currAlpha)
                                    currMu  = computeNewMean_w(dta,respVec,currSigma,currRegSig)
                                    if (iter == 0){ 
                                        currllk = LogLikelihood_type2(dta,respVec,currRegSig,currSigma,currMu,currAlpha)
                                    }
                           ### M-step
                                    newAlpha = updatePriorVariance(currSigma,currMu)
                                    currRegSig = updatePriorRegSig(dta,respVec,currMu,currSigma,currRegSig,currAlpha)
                                    currAlpha = newAlpha
                                    newllk = LogLikelihood_type2(dta,respVec,currRegSig,currSigma,currMu,currAlpha)
                                    #check for convergence
                                    if ((newllk-currllk) < 0.00000001) { 
                                       print(paste("SBL EM has converged -- num_iter:",iter))
                                       hasConverged = TRUE
                                       #print(paste("SBL Converged Iter: ",iter," regSig: ", currRegSig, "prevLLK: ",currllk," newLLK: ",newllk))
                                       #print("SBLCurrAlpha:")
                                       #print(currAlpha)
                                       #print("SBLCurrMu:")
                                       #print(currMu)
                                       #print("SBLCurrSig:")
                                       #print(currSigma)
                                    }
                                     

                                    if (newllk < currllk) {
                                       print(paste("SBL Warning: decrease in log-likelihood!", currllk, newllk))
                                    } 

                                    #print(paste("SBL Iter: ",iter," regSig: ", currRegSig, "prevLLK: ",currllk," newLLK: ",newllk))
                                    iter = iter + 1
                                    currllk = newllk 
                           }
                           if (!hasConverged){
                              print("SBL Exiting after MAX_ITER")
                              #print("SBL CurrMu:")
                              #print(currMu)
                              #print("SBL regSig:")
                              #print(currRegSig)
                           }
                           optParams[[1]] = currMu
                           optParams[[2]] = currRegSig
                           optParams[[3]] = currAlpha
                           optParams[[4]] = currllk
                           return(optParams)
}



mainSBL <- function(dta,respVec,sampleWts,initW,initRegSig,initAlpha) {

        #####read the data

        #clDir <- "/home/rhosur/workspace/ExpressionAnalysis/"
        #clDir <- "/home/rhosur/Projects/src/"
        #allDta <- read.table(paste(clDir,"ACP_ChaussabelProjections_NOoutliers.txt",sep=""),header=TRUE)
        #allDta_2 <- read.table("SyntheticMixture.txt",header=TRUE,check.names=F,stringsAsFactors=F)
        #allDta_allProbes <- read.table(paste(clDir,"ACP_SigProbes_t.txt",sep=""),header=TRUE)

        #clVars <- read.table("SyntheticResponse.txt",header=TRUE,check.names=F,stringsAsFactors = F)
        #clVars <- read.table('/home/rhosur/Projects/ACP/ACP_ClinicalVars_MSSS.txt',header=TRUE)
        #row.names(clVars) <- clVars$SAMPLE_ID
        #cl_e <- clVars[which(clVars$SAMPLE_ID %in% rownames(allDta)),]
        #allDta <- allDta[which(rownames(allDta) %in% clVars$SAMPLE_ID),]
        #allDta_allProbes <- allDta_allProbes[which(rownames(allDta_allProbes) %in% clVars$SAMPLE_ID),]
        #clVars <- cl_e[rownames(allDta),,drop=FALSE]
        #clinicVars <- clVars$EDSS #- mean(clVars$EDSS)
        #sampleWts = 1.0 
        #add the intercept and scale the data and the responses with the weights 
        dta$Intercept = c(rep(1.0,nrow(dta)))
        allDta <- dta * sqrt(sampleWts)
        clinicVars <- respVec * sqrt(sampleWts)

        #print(sampleWts)

        #initial values
        nfeats = length(colnames(allDta))
        #print(nfeats)
        #if (nrow(allDta) > nfeats) {
            #initModel <- lm(clinicVars ~ 0+.,data=allDta)
            #W = initModel$coefficients
        #} else {
        #W = c(rep(1.0,nfeats))
        #}
        #print("Initial Weights:")
        W = initW
        #print(W)
        #don't want to penalize the intercept too much, so put a very small prior inv cov on that....
        #alpha=c(rep(0.1,nfeats-1),1e-5) 
        #print("Initial Alpha:")
        alpha = initAlpha
        #print(alpha)
        regSig = initRegSig
        #regSig = 0.1*var(clinicVars)

        #print("Initial regSig:")
        #print(regSig) 
        MAX_ITER = 10000
        optimalParams = sparseBayesianRegression(allDta,clinicVars,W,alpha,regSig,MAX_ITER)  
        return(optimalParams)
}

#mainSBL()


 
