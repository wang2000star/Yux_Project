import math
import re
import numpy as np
import scipy.special
# 统计"+="和"="的次数
def count_operators(code):
    plus_equal_count = code.count("+=")
    equal_count = code.count("=") - plus_equal_count
    return plus_equal_count, equal_count
def generate_code(num_cols, num_rows, M, L,idbool):
    code1 = ""
    code3 = ""
    def generate_recursive_code(prefix, indices, depth, max_depth, s1, s2):
        nonlocal code1, code3
        if depth == max_depth:
            j_indices = '_'.join(map(str, indices))
            if idbool:
                code1 += f"Ctxt temp{j_indices} = temp{prefix};\n"
                code1 += f"temp{j_indices} += temp[id{indices[-1]}];\n"
                code3 += f"temp{j_indices} = temp{prefix};\n"
                code3 += f"temp{j_indices} += temp[id{indices[-1]}];\n"
            else:
                code1 += f"Ctxt temp{j_indices} = temp{prefix};\n"
                code1 += f"temp{j_indices} += temp[{indices[-1]}];\n"
                code3 += f"temp{j_indices} = temp{prefix};\n"
                code3 += f"temp{j_indices} += temp[{indices[-1]}];\n"
        else:
            start = indices[-1] + 1 if indices else s1
            for j in range(start, s2 - (max_depth - depth - 1)):
                generate_recursive_code(f"{prefix}_{j}" if prefix else f"{j}", indices + [j], depth + 1, max_depth, s1, s2)

    for block in range(0, math.ceil(num_cols / L)):
        s1 = L * block
        s2 = min(L * (block + 1), num_cols)
        for depth in range(2, L + 1):
            generate_recursive_code("", [], 0, depth, s1, s2)   
    code2 = ""
    for i in range(num_rows):
        for jj in range(num_cols):
            if M[i][jj] == 1:
                min1 = jj
                break
        min1block = min1 // L
        index = []
        for k in range(min1, L * (min1block + 1)):
            if M[i][k] == 1 and k != i:
                index.append(k)
        if len(index) > 0:
            if M[i][i] == 0:
                if len(index) > 1:
                    index_str = "_".join(map(str, index))
                    if idbool:
                        expression = f"eData[id{i}] = temp{index_str};\n"
                    else:
                        expression = f"eData[{i}] = temp{index_str};\n"
                else:
                    if idbool:
                        expression = f"eData[id{i}] = temp[id{index[0]}];\n"
                    else:
                        expression = f"eData[{i}] = temp[{index[0]}];\n"
                code2 += expression
            else:
                if len(index) > 1:
                    index_str = "_".join(map(str, index))
                    if idbool:
                        expression = f"eData[id{i}] += temp{index_str};\n"
                    else:
                        expression = f"eData[{i}] += temp{index_str};\n"
                else:
                    if idbool:
                        expression = f"eData[id{i}] += temp[id{index[0]}];\n"
                    else:
                        expression = f"eData[{i}] += temp[{index[0]}];\n"
                code2 += expression
        for block in range(min1block + 1, num_cols // L):
            index = []
            for k in range(L * block, L * (block + 1)):
                if M[i][k] == 1 and k != i:
                    index.append(k)
            if len(index) > 0:
                if len(index) > 1:
                    index_str = "_".join(map(str, index))
                    if idbool:
                        expression = f"eData[id{i}] += temp{index_str};\n"
                    else:
                        expression = f"eData[{i}] += temp{index_str};\n"
                else:
                    if idbool:
                        expression = f"eData[id{i}] += temp[id{index[0]}];\n"
                    else:
                        expression = f"eData[{i}] += temp[{index[0]}];\n"
                code2 += expression
        if L * (num_cols // L) < num_cols:
            index = []
            for k in range(L * (num_cols // L), num_cols):
                if M[i][k] == 1 and k != i:
                    index.append(k)
            if len(index) > 0:
                if len(index) > 1:
                    index_str = "_".join(map(str, index))
                    if idbool:
                        expression = f"eData[id{i}] += temp{index_str};\n"
                    else:
                        expression = f"eData[{i}] += temp{index_str};\n"
                else:
                    if idbool:
                        expression = f"eData[id{i}] += temp[id{index[0]}];\n"
                    else:
                        expression = f"eData[{i}] += temp[id{index[0]}];\n"
                code2 += expression
    return code1, code2, code3
def optimize_code(data,num_cols, num_rows, M, L, idbool):
    code1, code2, code3 = generate_code(num_cols, num_rows, M, L, idbool)
    plus_equal_count1, equal_count1 = count_operators(code1)
    plus_equal_count2, equal_count2 = count_operators(code2)
    # print("初步优化结果：")
    # print("+=个数：", plus_equal_count1, " =个数：", equal_count1)
    all6 = set()
    all5 = set()
    all4 = set()
    all3 = set()
    all2 = set()
    S6 = set()
    S5 = set()
    S4 = set()
    S3 = set()
    S2 = set()
    Q2 = set()
    Q3 = set()
    Q4 = set()
    Q5 = set()
    Q6 = set()
    for block in range(0, math.ceil(num_cols / L)):
        s1 = L * block
        if L * (block + 1) <= num_cols:
            s2 = L * (block + 1)
        else:
            s2 = num_cols
        if L >= 6:
            for j1 in range(s1, s2 - 5):
                for j2 in range(j1 + 1, s2 - 4):
                    for j3 in range(j2 + 1, s2 - 3):
                        for j4 in range(j3 + 1, s2 - 2):
                            for j5 in range(j4 + 1, s2 - 1):
                                for j6 in range(j5 + 1, s2):
                                    all6.add(f"temp{j1}_{j2}_{j3}_{j4}_{j5}_{j6}")
                                    if re.search(rf'\btemp{j1}_{j2}_{j3}_{j4}_{j5}_{j6}\b', code2):
                                        S6.add(f"temp{j1}_{j2}_{j3}_{j4}_{j5}_{j6}")
                                        Q6.add(f"temp{j1}_{j2}_{j3}_{j4}_{j5}")
                                        Q5.add(f"temp{j1}_{j2}_{j3}_{j4}")
                                        Q4.add(f"temp{j1}_{j2}_{j3}")
                                        Q3.add(f"temp{j1}_{j2}")
        S5 = Q6
        if L >= 5:
            for j1 in range(s1, s2 - 4):
                for j2 in range(j1 + 1, s2 - 3):
                    for j3 in range(j2 + 1, s2 - 2):
                        for j4 in range(j3 + 1, s2 - 1):
                            for j5 in range(j4 + 1, s2):
                                all5.add(f"temp{j1}_{j2}_{j3}_{j4}_{j5}")
                                if re.search(rf'\btemp{j1}_{j2}_{j3}_{j4}_{j5}\b', code2):
                                    S5.add(f"temp{j1}_{j2}_{j3}_{j4}_{j5}")
                                    Q5.add(f"temp{j1}_{j2}_{j3}_{j4}")
                                    Q4.add(f"temp{j1}_{j2}_{j3}")
                                    Q3.add(f"temp{j1}_{j2}")
        S4 = Q5
        if L >= 4:
            for j1 in range(s1, s2 - 3):
                for j2 in range(j1 + 1, s2 - 2):
                    for j3 in range(j2 + 1, s2 - 1):
                        for j4 in range(j3 + 1, s2):
                            all4.add(f"temp{j1}_{j2}_{j3}_{j4}")
                            if re.search(rf'\btemp{j1}_{j2}_{j3}_{j4}\b', code2):
                                S4.add(f"temp{j1}_{j2}_{j3}_{j4}")
                                Q4.add(f"temp{j1}_{j2}_{j3}")
                                Q3.add(f"temp{j1}_{j2}")
        S3 = Q4
        if L >= 3:
            for j1 in range(s1, s2 - 2):
                for j2 in range(j1 + 1, s2 - 1):
                    for j3 in range(j2 + 1, s2):
                        all3.add(f"temp{j1}_{j2}_{j3}")
                        if re.search(rf'\btemp{j1}_{j2}_{j3}\b', code2):
                            S3.add(f"temp{j1}_{j2}_{j3}")
                            if re.search(rf'\btemp{j1}_{j2}\b', code2):
                                Q3.add(f"temp{j1}_{j2}")
        S2 = Q3
        if L >= 2:
            for j1 in range(s1, s2 - 1):
                for j2 in range(j1 + 1, s2):
                    all2.add(f"temp{j1}_{j2}")
                    if re.search(rf'\btemp{j1}_{j2}\b', code2):
                        S2.add(f"temp{j1}_{j2}")
    T = len(all2) + len(all3) + len(all4) + len(all5) + len(all6) - (len(S2) + len(S3) + len(S4) + len(S5) + len(S6))
    # print("最终优化结果：")
    # print("+=个数：", plus_equal_count1 + plus_equal_count2 - T, " =个数：", equal_count1 + equal_count2 - T)
    data.append([L, plus_equal_count1 + plus_equal_count2 - T, equal_count1 + equal_count2 - T])
    new_code1 = ""
    new_code3 = ""
    if L >= 2:
        for str2 in S2:
            parts = str2.split('_')
            j1 = int(parts[0][4:])  # 去掉'temp'前缀并转换为整数
            j2 = int(parts[1])  # 直接转换为整数
            if idbool:
                new_code1 += f"Ctxt temp{j1}_{j2} = temp[id{j1}];\n"
                new_code1 += f"temp{j1}_{j2} += temp[id{j2}];\n"
                new_code3 += f"temp{j1}_{j2} = temp[id{j1}];\n"
                new_code3 += f"temp{j1}_{j2} += temp[id{j2}];\n"
            else:
                new_code1 += f"Ctxt temp{j1}_{j2} = temp[{j1}];\n"
                new_code1 += f"temp{j1}_{j2} += temp[{j2}];\n"
                new_code3 += f"temp{j1}_{j2} = temp[{j1}];\n"
                new_code3 += f"temp{j1}_{j2} += temp[{j2}];\n"
    if L >= 3:
        for str3 in S3:
            parts = str3.split('_')
            j1 = int(parts[0][4:])
            j2 = int(parts[1])
            j3 = int(parts[2])
            if idbool:
                new_code1 += f"Ctxt temp{j1}_{j2}_{j3} = temp{j1}_{j2};\n"
                new_code1 += f"temp{j1}_{j2}_{j3} += temp[id{j3}];\n"
                new_code3 += f"temp{j1}_{j2}_{j3} = temp{j1}_{j2};\n"
                new_code3 += f"temp{j1}_{j2}_{j3} += temp[id{j3}];\n"
            else:
                new_code1 += f"Ctxt temp{j1}_{j2}_{j3} = temp{j1}_{j2};\n"
                new_code1 += f"temp{j1}_{j2}_{j3} += temp[{j3}];\n"
                new_code3 += f"temp{j1}_{j2}_{j3} = temp{j1}_{j2};\n"
                new_code3 += f"temp{j1}_{j2}_{j3} += temp[{j3}];\n"
    if L >= 4:
        for str4 in S4:
            parts = str4.split('_')
            j1 = int(parts[0][4:])
            j2 = int(parts[1])
            j3 = int(parts[2])
            j4 = int(parts[3])
            if idbool:
                new_code1 += f"Ctxt temp{j1}_{j2}_{j3}_{j4} = temp{j1}_{j2}_{j3};\n"
                new_code1 += f"temp{j1}_{j2}_{j3}_{j4} += temp[id{j4}];\n"
                new_code3 += f"temp{j1}_{j2}_{j3}_{j4} = temp{j1}_{j2}_{j3};\n"
                new_code3 += f"temp{j1}_{j2}_{j3}_{j4} += temp[id{j4}];\n"
            else:
                new_code1 += f"Ctxt temp{j1}_{j2}_{j3}_{j4} = temp{j1}_{j2}_{j3};\n"
                new_code1 += f"temp{j1}_{j2}_{j3}_{j4} += temp[{j4}];\n"
                new_code3 += f"temp{j1}_{j2}_{j3}_{j4} = temp{j1}_{j2}_{j3};\n"
                new_code3 += f"temp{j1}_{j2}_{j3}_{j4} += temp[{j4}];\n"
    return new_code1, code2, new_code3


def generate_random_matrix(n, m, sparsity):
    """
    生成一个随机的n x m的01矩阵，指定稀疏度（1的概率）
    
    参数:
    n (int): 矩阵的行数
    m (int): 矩阵的列数
    sparsity (float): 1的概率，范围在0到1之间
    
    返回:
    np.ndarray: 生成的随机01矩阵
    """
    return np.random.choice([0, 1], size=(n, m), p=[1-sparsity, sparsity])

if __name__ == "__main__":
    M1 = np.array(
        [[1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0],
        [1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1],
        [0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0],
        [0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1],
        [0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1],
        [1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0],
        [1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1],
        [0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1],
        [1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1],
        [0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1],
        [0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1],
        [1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0],
        [1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0],
        [1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1],
        [0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1],
        [0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1],
        [0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0],
        [1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0],
        [1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0],
        [1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1],
        [1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1],
        [0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1],
        [1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0],
        [0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1],
        [0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1],
        [0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1],
        [1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1],
        [1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0],
        [1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
        [0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1],
        [1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0],
        [0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1],
        [1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1],
        [1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1],
        [0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1],
        [0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0]]);
    M2 = np.array(
[[1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1],
[0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0],
[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0],
[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1],
[0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1],
[0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1],
[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0],
[1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0],
[0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
[0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0],
[0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0],
[1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1],
[0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
[0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0],
[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0],
[0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
[0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1],
[1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0],
[0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
[1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1],
[0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0]])
    # sparsity = 0.15 #1的比例
    # M_row = 32
    # M_col = 32
    # M = generate_random_matrix(M_row, M_col, sparsity)
    M = M2
    M_row = M.shape[0]
    M_col = M.shape[1]
    # 统计M中1的个数
    plus_equal_count = np.sum(M)-M_row
    print("不优化结果: ")
    print("+=个数: ", plus_equal_count)
    #统计M对角线上的1的个数
    equal_count = M_row - np.sum(np.diag(M))
    print("=个数: ", equal_count)
    # 定义转换数
    rate = 0.82
    print("等效转换+=个数：", plus_equal_count + rate*equal_count)  
    num_rows = M.shape[0]
    num_cols = M.shape[1] if num_rows > 0 else 0
    idbool = 0 # 是否含有id
    data = []
    N = M_col
    w = scipy.special.lambertw(N/math.e).real
    k = math.log2(math.e)*(1+w)
    print("k=",k)
    Lmax  = 6#23,4,5,6
    for L in range(4, 5):
        # code1_1, code1_2, code1_3 = optimize_code(data,num_cols, num_rows, M1, L, idbool)
        # code2_1, code2_2, code2_3 = optimize_code(data,num_cols, num_rows, M2, L, idbool)
        code1_1, code1_2, code1_3 = generate_code(num_cols, num_rows, M1, L, idbool)
    
    # 求code1_3与code2_3的差集
    #print(len(list(set(code2_3)-set(code1_3))))
    
    # # 求解最优解data第二个指标最小
    # data = sorted(data, key=lambda x: x[1])
    # print(data)
    # #print(f"{data[0][1] / plus_equal_count:.4f}")
    # # 进行等效转换
    # for i in range(len(data)):
    #     print(data[i][0], data[i][1]+data[i][2]*rate)
    
    # 写入文件
    with open("autoFR1.txt", "w") as file1:
        file1.write(code1_1)
    print("autoFR1.txt has been generated")

    with open("autoFR2.txt", "w") as file2:
        file2.write(code1_2)
    print("autoFR2.txt has been generated")

    with open("autoFR3.txt", "w") as file3:
        file3.write(code1_3)
    print("autoFR3.txt has been generated")
