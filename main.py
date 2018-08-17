from RandomTestModule import *
from RandomGenerator import *
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
        print(seedlist[i])
        rng = GetRngDpFunc(func=func,seed=seedlist[i])
        res.append(rng.GenerateRandomU01())
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

if __name__=="__main__":
    # MT19937
    #res=MT19937(autocorrelation_tests0x12341).GenerateRandomU01(N=100000)
    # python
    N=1000
    bins=10
    seed=0x123451
    rowname = []
    data=[]
    res = []
    columnName = ["randomtype","chisquare-pvalue","kstest-pvalue","autocor","serial2D","kl/bins","monobit","blockfreq"]
    #columnName = ["N=%d,seed=%d"%(N,seed),"chisquare-pvalue"]
    #rngtypes = ['mt19937']#,'pythonsys','well1024','xorshift','ThreeFry','LCG64']
    rngtypes = ['pythonsys','well1024','xorshift','ThreeFry','LCG64']
    for rngtype in rngtypes:
        datalist = GenerateRandom01Arr(groups=20,N=N,seed=seed,func=rngtype)
        rndtest = randomtest(datalist=datalist)
        passratio, res = rndtest.ChiSquare(bins=bins,confidence=0.8)
        print(passratio)
        #print(res)
        #print("%s chisq pv is %f"%(rngtype,rndtest.ChiSquare(bins=bins,confidence=0.8)))
        #index = rndtest.GeneratemdTable(bins=bins,params = columnName)

    #GeneratemdFile(columnName=columnName,rowname=rowname,content=data,outputpath='output.md')


