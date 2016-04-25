import os,sys,numpy, random


numClusters = 2
numDims = 10
numPts =  30


def getSigma(d,delta):
    #generate a cov dxd matrix
    A = numpy.random.rand(d,d)
    B = A*A.transpose() 
    return B

def getMu(d,delta):
    A = numpy.random.rand(d)
    return A+delta


def getBeta(d,sc,clNum):
    Beta = numpy.zeros(d)
    #generate some sparsity
    #put the clNum'th interval of coefficients to non-zero
    k = int(numDims/numClusters)
    dmin = clNum*k 
    dmax =  clNum*k+k 
    for i in range(dmin,dmax):
        r = numpy.random.rand(1)
        Beta[i] = r
    return sc*Beta 


def replaceBrackets(string): 
    str1 = string.replace('[','')
    str1 = str1.replace(']','')
    return str1

def featureString(numFeats):
    featstr = ""
    for i in range(0,numFeats):
        if len(featstr) == 0:
           featstr += "Feat"+str(i)
        else:
           featstr += "\tFeat"+str(i)
    return featstr

def responseString():
    return "EDSS"    

if __name__ == "__main__":
   f = open('SyntheticMixture.txt','w')
   print >> f, featureString(numDims)
   f1 = open('SyntheticResponse.txt','w')
   print >> f1, responseString()
   numpy.set_printoptions(precision=4,threshold=20000)
   old_sig = getSigma(numDims,1)
   old_mu = getMu(numDims,0.01)

   for k in range(0,numClusters):
       sigmaRand = (k+1)*random.random()
       currS = getSigma(numDims,k+1)
       currS= old_sig
       #currS = numpy.identity(numDims)
       currMu = getMu(numDims,k*0.001)
       currMu = old_mu
       currBeta = getBeta(numDims,10.0000,k)
       print "Cluster: ",k
       print currMu
       print currBeta
       print currS
       #print numpy.linalg.inv(currS)
       dtps = numpy.random.multivariate_normal(currMu,currS,numPts)
       #print len(dtps)
       dtps_str = numpy.array_str(dtps,max_line_width=200000,precision=None,suppress_small=None)
       #print dtps_str
       dtps_str = dtps_str.replace('[','')
       dtps_str = dtps_str.replace(']','')
       dtpts = dtps_str.split('\n')
       for i,dtp in enumerate(dtpts):
           str_w = "Cls"+str(k)+"_"+str(i)+" " + dtp
           print >> f, str_w
           ###synthetic response
           pt = dtps[i]
           resp = numpy.mat(pt) * numpy.mat(currBeta).transpose()
           resp_str = numpy.array_str(resp)
           resp_str= replaceBrackets(resp_str)
           #add an intercept
           yInt = 2+1*numpy.random.rand(1)
           #print yInt
           yInt = replaceBrackets(numpy.array_str(yInt))
           print >> f1, "Cls" + str(k) + "_"+str(i) + "  "+str(eval(yInt)+eval(resp_str))



   f.close()
   f1.close()        
