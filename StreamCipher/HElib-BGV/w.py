import numpy as np

M = np.array([
    [-8192, 2048, -512, 128, -32, 8, -2, -32768],
    [-8192, -32768, 2, 8, 32, 128, 512, 2048],
    [-8192, -128, -2, 2048, 32, -32768, -512, -8],
    [-8192, -8, 512, -32768, -32, 2048, 2, -128],
    [-8192, -2048, -512, -128, -32, -8, -2, 32768],
    [-8192, 32768, 2, -8, 32, -128, 512, -2048],
    [-8192, 128, -2, -2048, 32, 32768, -512, 8],
    [-8192, 8, 512, 32768, -32, -2048, 2, 128]
])

scalar = 8
modulus = 65537

result = (M * scalar) % modulus
print(result)