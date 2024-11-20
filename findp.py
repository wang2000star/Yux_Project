import sympy
from math import log2

def find_prime_numbers(p_min, p_max, n_min, N):
    primes = []  # 保存找到的素数
    n = n_min
    
    nmax = int(log2(p_max-1))+1  # n 的最大值
    print(f"nmax = {nmax}")
    while len(primes) < N: 
        print(f"n = {n}")
        kmin = (p_min-1)//(2**n)  # k 的最小值
        print(f"kmin = {kmin}")
        kmax = (p_max-1)//(2**n)  # k 的最大值
        print(f"kmax = {kmax}")
        for k in range(1, kmax+1, 2):
            # 计算 p = k * 2^n + 1
            if k == 78557:
                continue
            p = k * (2**n) + 1
            print(f"p = {p}")
            # 检查 p 是否大于 p_min 并且是素数
            if sympy.isprime(p):
                primes.append((p, n))
        n += 1
        if n > nmax:
            break
    
    return primes

# 示例输入
p_min = 2**16  # 最小 p 的值, p=k*2^n+1,k是奇数
p_max = 2**17  # 最大 p 的值
n_min = 15    # 最小 n 的值
N = 100       # 找到的素数个数

# 调用函数
primes = find_prime_numbers(p_min, p_max, n_min, N)

# 对primes元素进行排序，首先是p，然后是n
primes.sort(key=lambda x: (x[1], x[0]))
# 存入文件findp.txt
with open("findp.txt", "w") as f:
    for p, n in primes:
        f.write(f"p = {p}, n = {n}\n")