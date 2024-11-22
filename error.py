from estimator import Param, estimate_lwe

# 定义 LWE 参数
n = 256        # LWE 维数
q = 2**15      # 模数
sd = 3.2       # 噪声标准差

# 初始化参数
params = Param(n=n, q=q, secret_distribution="normal", m=float("inf"), sd=sd)

# 估算安全性
result = estimate_lwe(params)

# 输出估算结果
print("安全性估算结果:")
print(result)
print("Classical security level:", result["rop"])
print("Quantum security level:", result["red_cost"])
