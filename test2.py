from RandomTestModule import *
from RandomGenerator import *
import warnings
import matplotlib.pyplot as plt
from scipy.stats import chisquare
from functools import reduce

warnings.filterwarnings("ignore", message="numpy.dtype size changed")


if __name__=="__main__":

    a = [498,520,538,501,475,525,493,474,448,500,471,495,512,532,510,553,499,456
,523,477]
    b = [500 for i in range(20)]
    chisq = sum([(e-f)**2/f for e,f in zip(a,b)])
    print(chisq)
