#!/usr/bin/python3
#coding=utf-8
import random as rnd
from randomgen import RandomGenerator, ThreeFry, PCG64, Xorshift1024
from wellrng.well1024a import WELL1024a

def _int32(x):
    return int(0xFFFFFFFF & x)

class MT19937:
    def __init__(self, seed):
        self.mt = [0] * 624
        self.mt[0] = seed
        for i in range(1, 624):
            self.mt[i] = _int32(1812433253 * (self.mt[i - 1] ^ self.mt[i - 1] >> 30) + i)


    def extract_number(self):
        self.twist()
        y = self.mt[0]
        y = y ^ y >> 11
        y = y ^ y << 7 & 2636928640
        y = y ^ y << 15 & 4022730752
        y = y ^ y >> 18
        return _int32(y)

    def Uniform(self,low=0.0,high=1.0):
        return (high-low)*self.extract_number()/(2**32)


    def twist(self):
        for i in range(0, 624):
            y = _int32((self.mt[i] & 0x80000000) + (self.mt[(i + 1) % 624] & 0x7fffffff))
            self.mt[i] = y ^ self.mt[(i + 397) % 624] >> 1
            if y % 2 != 0:
                self.mt[i] = self.mt[i] ^ 0x9908b0df

    def GenerateRandomU01(self,N=1000):
        res = []
        for i in range(N):
            res.append(self.Uniform())
        return res

    def Generate01Sequence(self, N=1000):
        res = self.GenerateRandomU01(N=N)
        return [0 if e < 0.5 else 1 for e in res]

class pythonsytem:
    def __init__(self,seed = 0):
        self.rnd = rnd
        self.rnd.seed(seed)

    def GenerateRandomU01(self,N=1000):
        res = []
        for i in range(N):
            res.append(self.rnd.uniform(0.0,1.0))
        return res

    def Generate01Sequence(self, N=1000):
        res = self.GenerateRandomU01(N=N)
        return [0 if e < 0.5 else 1 for e in res]


class Well1024a:
    def __init__(self,seed=0):
        self.rng = WELL1024a(seed)

    def GenerateRandomU01(self,N=1000):
        res = []
        for i in range(N):
            res.append(self.rng.random())
        return res

    def Generate01Sequence(self, N=1000):
        res = self.GenerateRandomU01(N=N)
        return [0 if e < 0.5 else 1 for e in res]


class xorshift1024:
    def __init__(self,seed=0):
        self.rng = RandomGenerator(Xorshift1024(seed))
        #self.rng.seed(seed)

    def GenerateRandomU01(self,N=1000):
        res = []
        for i in range(N):
            res.append(self.rng.random_sample())
        return res

    def Generate01Sequence(self, N=1000):
        res = self.GenerateRandomU01(N=N)
        return [0 if e < 0.5 else 1 for e in res]

class MyPCG64:
    def __init__(self,seed=0):
        self.rng = RandomGenerator(PCG64(seed))
        #self.rng.seed(seed)

    def GenerateRandomU01(self,N=1000):
        res = []
        for i in range(N):
            res.append(self.rng.random_sample())
        return res

    def Generate01Sequence(self, N=1000):
        res = self.GenerateRandomU01(N=N)
        return [0 if e < 0.5 else 1 for e in res]

class MyThreeFry:
    def __init__(self,seed=0):
        self.rng = RandomGenerator(ThreeFry(seed))
        #self.rng.seed(seed)

    def GenerateRandomU01(self,N=1000):
        res = []
        for i in range(N):
            res.append(self.rng.random_sample())
        return res

    def Generate01Sequence(self, N=1000):
        res = self.GenerateRandomU01(N=N)
        return [0 if e < 0.5 else 1 for e in res]



class ReadFromPath:
    def __init__(self, path):
        self.path = path

    def readbyline(self):
        res = []
        with open(self.path,'r') as f:
            for line in f.readlines():
                res.append(float(line))
        return res

    def GenerateRandomU01(self,N=1000):
        res = self.readbyline()
        assert len(res) >= N
        res = res[0:N]
        return res

    def Generate01Sequence(self, N=1000):
        res = self.GenerateRandomU01()
        return [0 if e < 0.5 else 1 for e in res]

if __name__=="__main__":
    print('#'*20+' Test MT19937 '+'#'*20)
    #rand=MT19937(0x12341)
    #for i in range(30):
    #    print(rand.Uniform())

    rand = MyPCG64(1010104)
    print(rand.GenerateRandomU01(N=30))
