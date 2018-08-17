import numpy as np
from itertools import islice
from scipy.stats import chisquare
from scipy.stats import kstest
from scipy.stats import norm
import scipy.special as spc
import math
import warnings
warnings.filterwarnings("ignore",message="numpy.dtype size changed")

"""
all metric:

"""

class randomtest:
    """
    chi-square:
        Most commonly used test
        Can be used for any distribution
        It's none fo bin k
        prepare a hist for the obeserved data
        compare obeseved frequencies with theoretical one
        D = \sum_{i=1}^k\frac{(o_i-e_i)^2}{e_i}
        D = 0 => Exact fit
        D has a chi-square distribution with k-1 degrees of freedom

    """
    def __init__(self,datalist = None):
        self.datalist = datalist
        self.groups = len(datalist)
        # every datapoint must belong to [0.0,1.0)
        assert self.isUniform() is True

    def isUniform(self):
        for e in self.datalist:
            for f in e:
                if f < 0.0 or f >= 1.0:
                    return False
        return True

    """
    The below 3 function implemented for Chi-square Test
    Most commonly ,so ...
    """


    def ChiSquare(self, bins=10, confidence = 0.8):
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
        for data in self.datalist:
            bincount = BinCountForData(data,bins=bins)
            expectedcount = ExpectedCountForData(data,bins=bins)
            chisq,pvalue = chisquare(bincount,expectedcount)
            #print(bincount)
            #print(expectedcount)
            res.append(chisq)
            if pvalue > confidence:
                passcnt += 1
        #print(bincount)
        #print(expectedcount)
        # sum([(a-e)**2/e for a,e in zip(list(bincount), list(expectedcount))])
        return 1.0*passcnt/len(self.datalist),res

    """
    the below function implemented for Kolmogorov-Smirnov Test
    K^{+} is maximum ..
    transfer data to cdf compared with uniform/u(0,1)
    """
    def KSTest(self, alpha):
        passcnt = 0
        for data in self.datalist:
            chisq,pvalue = kstest(data,'uniform',args=(0,1))
            if pvalue > alpha:
                passcnt += 1
        return 1.0*passcnt/len(self.datalist)

    """
    the below 2 function implemented for Serial-Correlation Test
    Rk = Autocovariance at lag 
    R_k = \frac{1}{n-k}\sum_{i=1}^{n-k}(U_i-0.5)(U_{i+k}-0.5)
    R_k is normally distributed with a mean of 0 and a variance of 1/(144*(n-k))
    when all Rk less than z_{1-\alpha/2}/(12*\sqrt{n-k})
    """


    def AutoCorTest(self, lagklist = [], confidence = 0.8):
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
            for data in self.datalist:
                Rk = AutoCovariance(lagk=e)
                if Rk < x0:
                    passcnt+=1
        return 1.0/passcnt/len(lagklist)/len(self.datalist)

    """
    Two level Test,
    when teh sample size is too small, the test results may apply to locally, but not
    globally, to the complete cycle.
    Use Chi-square test on n samples of size k each and use a chi-square test on \\
    the set of n Chi-square statistics so obtained
    Chi-square on Chi-square test
    """
    def TwoLevelTest(self):
        return

    #def kdistributivity(self):
    #    return

    """
    the below 3 function implemented for serial test
    the odd element as X-Axis, the even element as Y-Axis
    In 2D dimensions, divide the space between 0 and 1 into K^2 cell of equal area
    obey to the chi-square distribution with the freedom of K^2-1
    """


    def Serial2DTest(self, binsX = 10, binsY = 10, alpha = 0.8):
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
        for data in self.datalist:
            binCount = Bin2DCountForData(data, binsX = binsX, binsY = binsY)
            expectedCount = Expected2DCountForData(data, binsX = binsX, binsY = binsY)
            chisq,pvalue = chisquare(binCount,expectedCount)
            if pvalue >= alpha:
                passcnt += 1
        return 1.0*passcnt/len(self.datalist)

    """
    the below function implemented for 
    maybe every lcgs  has some hyperplane.
    necessary?
    """
    #def SpectralTest(self):
    #    return

    """
    def autocorrelationForStm(self, st=0, m=5):
        # Independence
        # refer to http://www.aritzhaupt.com/resource/autocorrelation/
        # z_score and norm?
        N = len(self.datalist)
        M = int((N-1-st+0.5)/m)
        # Ri,Ri+m,...
        Rm1 = self.datalist[0::m]
        rho_im = sum([e*f for e,f in zip(Rm1[0:M],Rm1[1:M+1])])/(M+1) - 0.25
        sigma_rhoim = np.sqrt(13*M+7)/(12*(M+1))
        z_0 = rho_im / sigma_rhoim
        return z_0    
    def autocorrelation_test(self):
        maxzvalue = 0.0
        for m in range(1,50,1):
            thiszvalue = self.autocorrelationForStm(st=0,m=m)
            if math.fabs(thiszvalue) > maxzvalue:
                maxzvalue = math.fabs(thiszvalue)
        return maxzvalue
    """

    """
    The below function is implemented for kl distance divided by bins
    kl_distance(binprob,expectedprob)/bins
    """


    def kldistance(self,bins = 10, epsilon = 0.00001):
        def BinProbForData(data,bins = 10):
            binprob,temp= np.histogram(data, bins = bins, range = (0.0,1.0), density = True)
            return binprob

        def ExpectedProbForData(data, bins = 10):
            expectedProb = np.ones((bins), dtype=np.float) / bins
            return expectedProb
        kl_dissum = 0.0
        for data in self.datalist:
            binProb = BinProbForData(data,bins=bins) + epsilon
            expectedProb = ExpectedProbForData(data,bins=bins) + epsilon
            kl_dissum += np.sum(binProb*np.log(binProb/expectedProb))/bins
        return kl_dissum/bins/len(self.datalist)

    """
    The below code implemented nist test codes
    monobit Test
    """
    def MonobitTest(self, alpha = 0.01):
        def get01Sequence(data):
            return [0 if e < 0.5 else 1 for e in data]

        passcnt = 0
        for data in self.datalist:
            sequence01 = get01Sequence(data)
            Sn = sum([2*e-1 for e in sequence01])
            Sobs = math.fabs(Sn)/math.sqrt(len(self.datalist))
            p_value = spc.erfc(Sobs / math.sqrt(2.0))
            if p_value > alpha:
                passcnt += 1
        return 1.0*passcnt/len(self.datalist)

    """
    Block frequency 
    """
    def BlockFrequency(self, blocksize=1000, alpha=0.01):
        def get01Sequence(data):
            return [0 if e < 0.5 else 1 for e in data]

        def PropotionOneInList(data):
            return 1.0*sum(data)/len(data)
        passcnt = 0
        for data in self.datalist:
            sequence01 = get01Sequence(data)
            block_num = int(len(sequence01) / blocksize)
            # calc chisq for all blocks
            chisq_obs = 0.0
            for i in range(block_num):
                pi_i = PropotionOneInList(sequence01[(i*blocksize):(i+1)*blocksize])
                chisq_obs += 4.0*blocksize*(pi_i-0.5)**2
            pvalue = spc.gammaincc(block_num/2.0,chisq_obs/2.0)
            passcnt = 0
            if pvalue > alpha:
                passcnt += 1
        return 1.0*passcnt/len(self.datalist)

    """
    runtest
    """
    def RunsTest(self):
        def get01Sequence(data):
            return [0 if e < 0.5 else 1 for e in data]
        chisq_obslist = []
        pvaluelist = []
        for data in self.datalist:
            seq01 = get01Sequence(data)
            n = len(seq01)
            pi = 1.0*sum(seq01)/len(seq01)
            chisq_obs = sum([0 if e == f else 1 for e,f in zip(seq01[0:len(seq01)-1],seq01[1:])]) + 1
            pvalue = spc.erfc(math.fabs(chisq_obs-2*n*pi*(1-pi))/(2*pi*(1-pi)*np.sqrt(2*n)))
            chisq_obslist.append(chisq_obs)
            pvaluelist.append(pvalue)
        return sum(chisq_obslist)/len(chisq_obslist),sum(pvaluelist)/len(pvaluelist)

    """
    CumultiveSumsTest
    """
    """
    def CumulativeSumsTest(self):
        def GetNPSequence(data):
            return [-1 if e < 0.5 else 1 for e in data]
        
        def ForwardMaxCumlativeSum(data):
            seq01 = ForwardMaxCumlativeSum(data)
            lastmax = 0
            lastsum = 0
            for e in seq01:
                lastsum += e
                if math.fabs(lastsum) > lastmax:
                    lastmax = math.fabs(lastsum)
            return lastmax
        for data in self.datalist:
            Sobs = max(ForwardMaxCumlativeSum(data),ForwardMaxCumlativeSum(reversed(data)))
    """


    def GeneratemdTable(self,bins=10,params=[]):
        res = []
        chisq,pvalue = self.ChiSquare(bins=bins)
        #if "chisquare" in params:
        #    res.append(chisq)
        if "chisquare-pvalue" in params:
            res.append(pvalue)
        ks,pvalue = self.KSTest()
        #if "kstest" in params:
        #    res.append(ks)
        if "kstest-pvalue" in params:
            res.append(pvalue)
        passornot,passcnt = self.AutoCorTest(lagklist = [5,8,10,15,20,30,60], confidence = 0.9)
        if "autocor" in params:
            res.append(passornot)

        chisq,pvalue = self.Serial2DTest(binsX = 10, binsY = 10)
        if "serial2D" in params:
            res.append(pvalue)

        kldis = self.kldistance(bins=bins)
        if "kl/bins" in params:
            res.append(kldis)

        passornot = self.MonobitTest(alpha=0.01)
        if "monobit" in params:
            res.append('Pass' if passornot else 'Failed')

        pvalue = self.BlockFrequency(blocksize=100)
        if "blockfreq" in params:
            res.append(pvalue)
        return res
