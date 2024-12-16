import math
import scipy.special
N = 32000
w = scipy.special.lambertw(N/math.e).real
k = math.log2(math.e)*(1+w)
print(k)
print(math.log2(N))
print(0.75*math.log2(N))