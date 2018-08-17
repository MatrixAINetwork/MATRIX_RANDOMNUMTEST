from RandomTestModule import *
from RandomGenerator import *
from scipy.stats import chisquare
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
    res = []
    passcnt = 0
    bincountSum = np.zeros((bins))
    for data in datalist:
        bincount = BinCountForData(data,bins=bins)
        expectedcount = ExpectedCountForData(data,bins=bins)
        chisq,pvalue = chisquare(bincount,expectedcount)
        #print(bincount)
        #print(expectedcount)
        res.append(pvalue)
        if pvalue > confidence:
            passcnt += 1
            bincountSum = bincountSum + bincount
    #print(bincount)
    #print(expectedcount)
    # sum([(a-e)**2/e for a,e in zip(list(bincount), list(expectedcount))])
    return 1.0*passcnt/len(datalist),res,bincountSum

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
    rngtypes = ['pythonsys']#,'pythonsys','well1024','xorshift','ThreeFry','LCG64']
    for rngtype in rngtypes:
        datalist = GenerateRandom01Arr(groups=2000,N=N,seed=seed,func=rngtype)
        passratio, pvalueres, bincountsum = ChiSquare(datalist,bins=bins,confidence=0.8)
        print(passratio)
        print(pvalueres)
        print(bincountsum)
