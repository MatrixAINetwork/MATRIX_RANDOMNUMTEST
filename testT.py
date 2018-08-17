from RandomTestModule import *
from RandomGenerator import *
import warnings
from contextlib import contextmanager
import time
warnings.filterwarnings("ignore", message="numpy.dtype size changed")

@contextmanager
def timer(title):
    t0 = time.time()
    yield
    print("{} - done in {:.3f}s".format(title, time.time() - t0))

groups = 200
seedlist = [i for i in range(groups)]
with timer('elapsed time for mt19937'):
    for i in range(groups):
        rng = MT19937(i)
        res = rng.GenerateRandomU01(N=10000)

with timer('elapsed time for pythonsys'):
    for i in range(groups):
        rng = pythonsytem(i)
        res = rng.GenerateRandomU01(N=10000)

with timer('elapsed time for Well1024a'):
    for i in range(groups):
        rng = Well1024a(i)
        res = rng.GenerateRandomU01(N=10000)

with timer('elapsed time for xorshift'):
    for i in range(groups):
        rng = xorshift1024(i)
        res = rng.GenerateRandomU01(N=10000)

with timer('elapsed time for MyThreeFry'):
    for i in range(groups):
        rng = MyThreeFry(i)
        res = rng.GenerateRandomU01(N=10000)

with timer('elapsed time for LCG64'):
    for i in range(groups):
        rng = MyPCG64(i)
        res = rng.GenerateRandomU01(N=10000)

