import math
import re
# 统计"+="和"="的次数
def count_operators(code):
    plus_equal_count = code.count("+=")
    equal_count = code.count("=") - plus_equal_count
    return plus_equal_count, equal_count
def generate_code(num_cols, num_rows, M, L,idbool):
    code1 = ""
    code3 = ""
    for block in range(0, math.ceil(num_cols / L)):
        s1 = L * block
        if L * (block + 1) <= num_cols:
            s2 = L * (block + 1)
        else:
            s2 = num_cols
        if L >= 2:
            for j1 in range(s1, s2 - 1):
                for j2 in range(j1 + 1, s2):
                    if idbool:
                        code1 += f"Ctxt temp{j1}_{j2} = temp[id{j1}];\n"
                        code1 += f"temp{j1}_{j2} += temp[id{j2}];\n"
                        code3 += f"temp{j1}_{j2} = temp[id{j1}];\n"
                        code3 += f"temp{j1}_{j2} += temp[id{j2}];\n"
                    else:
                        code1 += f"Ctxt temp{j1}_{j2} = temp[{j1}];\n"
                        code1 += f"temp{j1}_{j2} += temp[{j2}];\n"
                        code3 += f"temp{j1}_{j2} = temp[{j1}];\n"
                        code3 += f"temp{j1}_{j2} += temp[{j2}];\n"
        if L >= 3:
            for j1 in range(s1, s2 - 2):
                for j2 in range(j1 + 1, s2 - 1):
                    for j3 in range(j2 + 1, s2):
                        if idbool:
                            code1 += f"Ctxt temp{j1}_{j2}_{j3} = temp{j1}_{j2};\n"
                            code1 += f"temp{j1}_{j2}_{j3} += temp[id{j3}];\n"
                            code3 += f"temp{j1}_{j2}_{j3} = temp{j1}_{j2};\n"
                            code3 += f"temp{j1}_{j2}_{j3} += temp[id{j3}];\n"
                        else:
                            code1 += f"Ctxt temp{j1}_{j2}_{j3} = temp{j1}_{j2};\n"
                            code1 += f"temp{j1}_{j2}_{j3} += temp[{j3}];\n"
                            code3 += f"temp{j1}_{j2}_{j3} = temp{j1}_{j2};\n"
                            code3 += f"temp{j1}_{j2}_{j3} += temp[{j3}];\n"
        if L >= 4:
            for j1 in range(s1, s2 - 3):
                for j2 in range(j1 + 1, s2 - 2):
                    for j3 in range(j2 + 1, s2 - 1):
                        for j4 in range(j3 + 1, s2):
                            if idbool:
                                code1 += f"Ctxt temp{j1}_{j2}_{j3}_{j4} = temp{j1}_{j2}_{j3};\n"
                                code1 += f"temp{j1}_{j2}_{j3}_{j4} += temp[id{j4}];\n"
                                code3 += f"temp{j1}_{j2}_{j3}_{j4} = temp{j1}_{j2}_{j3};\n"
                                code3 += f"temp{j1}_{j2}_{j3}_{j4} += temp[id{j4}];\n"
                            else:
                                code1 += f"Ctxt temp{j1}_{j2}_{j3}_{j4} = temp{j1}_{j2}_{j3};\n"
                                code1 += f"temp{j1}_{j2}_{j3}_{j4} += temp[{j4}];\n"
                                code3 += f"temp{j1}_{j2}_{j3}_{j4} = temp{j1}_{j2}_{j3};\n"
                                code3 += f"temp{j1}_{j2}_{j3}_{j4} += temp[{j4}];\n"
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
    # print("+=个数：", plus_equal_count1 + plus_equal_count2, " =个数：", equal_count1 + equal_count2)
    all4 = set()
    all3 = set()
    all2 = set()
    S4 = set()
    S3 = set()
    S2 = set()
    Q2 = set()
    Q3 = set()
    Q4 = set()
    for block in range(0, math.ceil(num_cols / L)):
        s1 = L * block
        if L * (block + 1) <= num_cols:
            s2 = L * (block + 1)
        else:
            s2 = num_cols
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
    T = len(all4) - len(S4) + len(all3) - len(S3) + len(all2) - len(S2)
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
import numpy as np

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
    M1 = [[ 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1 ],
        [ 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1 ],
        [ 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0 ],
        [ 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0 ],
        [ 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0 ],
        [ 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1 ],
        [ 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1 ],
        [ 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0 ],
        [ 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0 ],
        [ 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1 ],
        [ 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0 ],
        [ 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1 ],
        [ 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0 ],
        [ 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0 ],
        [ 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1 ],
        [ 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1 ],
        [ 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0 ],
        [ 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0 ],
        [ 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1 ],
        [ 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1 ],
        [ 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1 ],
        [ 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0 ],
        [ 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0 ],
        [ 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1 ],
        [ 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1 ],
        [ 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0 ],
        [ 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1 ],
        [ 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0 ],
        [ 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1 ],
        [ 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0 ],
        [ 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1 ],
        [ 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0 ],
        [ 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0 ],
        [ 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1 ],
        [ 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0 ],
        [ 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0 ]
    ];
    M2 = [[ 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1 ],
    [ 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0 ],
    [ 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0 ],
    [ 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1 ],
    [ 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0 ],
    [ 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1 ],
    [ 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0 ],
    [ 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1 ],
    [ 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1 ],
    [ 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1 ],
    [ 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1 ],
    [ 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0 ],
    [ 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1 ],
    [ 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1 ],
    [ 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0 ],
    [ 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0 ],
    [ 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1 ],
    [ 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0 ],
    [ 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1 ],
    [ 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0 ],
    [ 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0 ],
    [ 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0 ],
    [ 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1 ],
    [ 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1 ],
    [ 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0 ],
    [ 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0 ],
    [ 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1 ],
    [ 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1 ],
    [ 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0 ],
    [ 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0 ],
    [ 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1 ],
    [ 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0 ],
    [ 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0 ],
    [ 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1 ],
    [ 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0 ],
    [ 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1 ]]

    sparsity = 1 #1的比例
    M_row = 36
    M_col = 36
    M = generate_random_matrix(M_row, M_col, sparsity)

    # 统计M中1的个数
    plus_equal_count = np.sum(M)-M_row
    print("不优化结果: ")
    print("+=个数: ", plus_equal_count)
    #统计M对角线上的1的个数
    equal_count = M_row - np.sum(np.diag(M))
    print("=个数: ", equal_count)
    num_rows = M.shape[0]
    num_cols = M.shape[1] if num_rows > 0 else 0
    idbool = 0 # 是否含有id
    data = []
    Lmax  = 6#M_col
    for L in range(2, Lmax+1):
        new_code1, code2, new_code3 = optimize_code(data,num_cols, num_rows, M, L, idbool)
    print(data)
    # 写入文件
    with open("autoFR1.txt", "w") as file1:
        file1.write(new_code1)
    print("autoFR1.txt has been generated")

    with open("autoFR2.txt", "w") as file2:
        file2.write(code2)
    print("autoFR2.txt has been generated")

    with open("autoFR3.txt", "w") as file3:
        file3.write(new_code3)
    print("autoFR3.txt has been generated")
