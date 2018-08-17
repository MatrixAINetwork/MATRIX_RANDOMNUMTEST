from RandomTestModule import *
from RandomGenerator import *
import warnings
import numpy as np
from itertools import islice
from scipy.stats import chisquare
from scipy.stats import kstest
from scipy.stats import norm
import scipy.special as spc
import math
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")

def GeneratemdFile(columnName = [], rowname = [], content = [], outputpath = None):
    with open(outputpath,'w') as f:
        f.write('|'.join(columnName)+'\n')
        f.write('|'.join(['---' for i in range(len(columnName))])+'\n')
        for c,e in zip(rowname,content):
            temp = [c]+e
            f.write('|'.join(temp)+'\n')
    return


def FromData2Index(datalist=None,bins=10,params=None):
    rndtest = randomtest(datalist=datalist)
    index = rndtest.GeneratemdTable(bins=bins,params = params)
    return ["%.5f"%e if type(e) != str else e for e in index]

def GetRngDpFunc(func='mt19937',seed=0):
    rng = MT19937(0)
    if func == 'mt19937':
        rng = MT19937(seed)
    elif func == 'pythonsys':
        rng = pythonsytem(seed)
    elif func == 'well1024':
        rng = Well1024a(seed)
    elif func == 'xorshift':
        rng = xorshift1024(seed)
    elif func == 'ThreeFry':
        rng = MyThreeFry(seed)
    elif func == 'LCG64':
        rng = MyPCG64(seed)
    return rng

def GenerateRandom01Arr(groups = 5, N = 10000, seed=0, func='mt19937'):
    seedlist = PrimeLessThanN(100000)
    #seedlist = [i for i in range(groups)]
    #print(len(seedlist))
    seedlist = seedlist[-groups:]
    #print(len(seedlist))
    #seedlist = [13,17,19,23,29]
    res = []
    for i in range(groups):
        rng = GetRngDpFunc(func=func,seed=seedlist[i])
        res.append(rng.GenerateRandomU01(N=N))
    return res

def PrimeLessThanN(N):
    ls = [True for i in range(N+1)]
    for i in range(2,N,1):
        if ls[i] is not True:
            continue
        for j in range(2*i,N+1,i):
            ls[j] = False
    res = []
    for i in range(2,N+1,1):
        if ls[i] == True:
            res.append(i)
    return  res


def ChiSquare(datalist, bins=10, confidence = 0.8):
    def BinCountForData(data,bins = 10):
        bincount,temp= np.histogram(data, bins = bins, range = (0.0,1.0), density = False)
        return bincount

    def ExpectedCountForData(data,bins = 10):
        expectedcount = np.ones((bins), dtype=np.float) * len(data) / bins
        return expectedcount
    # uniformity
    # p-value is the probability of assert H0 is true
    # H0 hypethesis: bincount Obey expected distribution
    chisqlist = []
    res = []
    bincountSum = np.zeros((bins))
    for data in datalist:
        bincount = BinCountForData(data,bins=bins)
        expectedcount = ExpectedCountForData(data,bins=bins)
        chisq,pvalue = chisquare(bincount,expectedcount)
        chisqlist.append(chisq)
        res.append(pvalue)
    return sum(chisqlist)/len(datalist),sum(res)/len(datalist)

def KSTest(datalist):
    passcnt = 0
    chisqlist = []
    pvaluelist = []
    for data in datalist:
        chisq,pvalue = kstest(data,'uniform',args=(0,1))

        chisqlist.append(chisq)
        pvaluelist.append(pvalue)
        #if pvalue > alpha:
        #    passcnt += 1
    return sum(chisqlist)/len(datalist),sum(pvaluelist)/len(datalist)

def Serial2DTest(datalist, binsX = 10, binsY = 10, alpha = 0.8):
    def Bin2DCountForData(data,binsX = 10,binsY = 10):
        oddlist = data[::2]
        evenlist = data[1::2]
        H,xedges,yedges = np.histogram2d(oddlist,evenlist,bins=[binsX,binsY],range=[[0.0,1.0],[0.0,1.0]])
        return H.ravel()

    def Expected2DCountForData(data, binsX = 10, binsY = 10):
        bins = binsX * binsY
        expectedcount = np.ones((bins), dtype=np.float) * len(data) / 2 / bins
        return expectedcount

    passcnt = 0
    chisqlist = []
    pvaluelist = []
    for data in datalist:
        binCount = Bin2DCountForData(data, binsX = binsX, binsY = binsY)
        expectedCount = Expected2DCountForData(data, binsX = binsX, binsY = binsY)
        chisq,pvalue = chisquare(binCount,expectedCount)
        chisqlist.append(chisq)
        pvaluelist.append(pvalue)
    return sum(chisqlist)/len(chisqlist),sum(pvaluelist)/len(pvaluelist)

def kldistance(datalist,bins = 10, epsilon = 0.00001):
    def BinProbForData(data,bins = 10):
        binprob,temp= np.histogram(data, bins = bins, range = (0.0,1.0), density = True)
        return binprob

    def ExpectedProbForData(data, bins = 10):
        expectedProb = np.ones((bins), dtype=np.float) / bins
        return expectedProb
    kl_dissum = 0.0
    for data in datalist:
        binProb = BinProbForData(data,bins=bins) + epsilon
        expectedProb = ExpectedProbForData(data,bins=bins) + epsilon
        kl_dissum += np.sum(binProb*np.log(binProb/expectedProb))/bins
    return kl_dissum/bins/len(datalist)

def MonobitTest(datalist, alpha = 0.01):
    def get01Sequence(data):
        return [0 if e < 0.5 else 1 for e in data]

    passcnt = 0
    Sobslist = []
    pvaluelist = []
    for data in datalist:
        sequence01 = get01Sequence(data)
        Sn = sum([2*e-1 for e in sequence01])
        Sobs = math.fabs(Sn)/math.sqrt(len(datalist))
        p_value = spc.erfc(Sobs / math.sqrt(2.0))
        Sobslist.append(Sobs)
        pvaluelist.append(p_value)
        #if p_value > alpha:
        #    passcnt += 1
    return sum(Sobslist)/len(Sobslist),sum(pvaluelist)/len(pvaluelist)


def BlockFrequency(datalist, blocksize=1000, alpha=0.01):
    def get01Sequence(data):
        return [0 if e < 0.5 else 1 for e in data]

    def PropotionOneInList(data):
        return 1.0*sum(data)/len(data)
    passcnt = 0
    for data in datalist:
        sequence01 = get01Sequence(data)
        block_num = int(len(sequence01) / blocksize)
        # calc chisq for all blocks
        chisq_obs = 0.0
        for i in range(block_num):
            pi_i = PropotionOneInList(sequence01[(i*blocksize):(i+1)*blocksize])
            chisq_obs += 4.0*blocksize*(pi_i-0.5)**2
        pvalue = spc.gammaincc(block_num/2.0,chisq_obs/2.0)
    return chisq_obs,pvalue

def AutoCorTest(datalist, lagklist = [], confidence = 0.8):
    def AutoCovariance(data,lagk = 10):
        Rksum = 0.0
        for i in range(len(data)-lagk):
            Rksum += (data[i]-0.5)*(data[i+lagk]-0.5)
        Rk = Rksum/(len(data)-lagk)
        return Rk
    passcnt = 0
    assert confidence <= 1.0 and confidence > 0.5
    alpha = 1-confidence
    x0 = norm(0,1).ppf(1-alpha/2)
    for e in lagklist:
        for data in datalist:
            Rk = AutoCovariance(data,lagk=e)
            if Rk < x0:
                passcnt+=1
    return 1.0*passcnt/len(lagklist)/len(datalist)




if __name__=="__main__":
    # MT19937
    #res=MT19937(autocorrelation_tests0x12341).GenerateRandomU01(N=100000)
    # python
    N=10000
    bins=10
    seed=0x123451
    rowname = []
    data=[]
    res = []
    #columnName = ["N=%d,seed=%d"%(N,seed),"chisquare-pvalue","kstest-pvalue","autocor","serial2D","kl/bins","monobit","blockfreq"]
    columnName = ["N=%d,seed=%d"%(N,seed),"chisquare-pvalue"]
    rngtypes = ['pythonsys','well1024','xorshift','ThreeFry','LCG64']
    for rngtype in rngtypes:
        datalist = GenerateRandom01Arr(groups=2000,N=N,seed=seed,func=rngtype)
        passratio = AutoCorTest(datalist,lagklist=[5,7,11,13,17,23,25,50,100,200,500])
        print(rngtype)
        print(passratio)
        #Echisq,Epvlue = ChiSquare(datalist,bins=bins,confidence=0.8)
