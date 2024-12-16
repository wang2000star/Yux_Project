import itertools
import numpy as np
def IndexOfoneInMatrix(matrix, m, n):
    # 定义元胞数组
    Index = [[] for _ in range(m)]
    for i in range(m):
        for j in range(n):
            if matrix[i][j] == 1:# and j != i:
                Index[i].append(j)
    return Index

def find_max_submatrix(Index, N):
    max_index_len = 0
    max_combinations = []
    
    for rows in itertools.combinations(range(len(Index)), N):
        intersection = set(Index[rows[0]])
        for row in rows[1:]:
            intersection &= set(Index[row])
        index_len = (len(intersection) - 1) * (N - 1)
        if index_len > max_index_len:
            max_index_len = index_len
            max_combinations = [rows]
        elif index_len == max_index_len:
            max_combinations.append(rows)
    
    return max_index_len, max_combinations

def add_to_set(N, intersection, num, sets):
    if num >= len(sets[N]):
        sets[N].append([])
    for r in intersection:
        sets[N][num].append(r)

def check_intersections(sets,setsN,NN,colM,Max,sumNNMax):
    sum_max=0
    if NN == 2:
        for i1 in range(0,setsN):
            flag = 0
            for i2 in range(i1,setsN): 
                for x1 in sets[i1]:
                    for x2 in sets[i2]:
                        # if len(set(x1))+len(set(x2)) > colM:
                        #     continue
                        if len(set(x1) & set(x2)) >= 0:
                            sum = (Max[i1] + Max[i2])-len(set(x1) & set(x2))
                            if sum_max < sum:
                                sum_max = sum
                                flag = 1
                                print(sum)
                                print(i1,i2)
                                print(x1)
                                print(x2)
                                print("************************")
                                if i1<setsN-1:
                                    if sum_max > NN*Max[i1+1]: 
                                        return sum_max
                                    break
                                break
                    if flag == 1:
                        break
                if flag == 1:
                    break
    if NN == 3: 
        for i1 in range(0,setsN):
            if NN * Max[i1] < sumNNMax:
                return sum_max 
            for i2 in range(i1,setsN):
                if Max[i1] + (NN-1)*Max[i2] < sumNNMax:
                    break
                flag = 0
                for i3 in range(i2,setsN): 
                    if Max[i1] + Max[i2] + Max[i3] < sumNNMax:
                        break
                    for x1 in sets[i1]:
                        for x2 in sets[i2]:
                            for x3 in sets[i3]:
                                if len(set(x1))+len(set(x2))+len(set(x3)) > colM:
                                    continue
                                if len(set(x1) & set(x2)) == 0 and len(set(x1) & set(x3)) == 0 and len(set(x2) & set(x3)) == 0:     
                                    sum = Max[i1] + Max[i2] + Max[i3]
                                    if sum_max < sum:
                                        sum_max = sum
                                        flag = 1
                                        print(sum)
                                        print(i1,i2,i3)
                                        print(x1)
                                        print(x2)
                                        print(x3)
                                        print("************************")
                                        if i1 < setsN-1:
                                            if sum_max > NN*Max[i1+1]:
                                                return sum_max
                                            break
                                        break        
                            if flag == 1:
                                break
                        if flag == 1:
                            break
                    if flag == 1:
                        break
    if NN == 4: 
        for i1 in range(0,setsN):
            if NN * Max[i1] < sumNNMax:
                return sum_max 
            for i2 in range(i1,setsN):
                if Max[i1] + (NN-1)*Max[i2] < sumNNMax:
                    break
                for i3 in range(i2,setsN):
                    if Max[i1] + Max[i2] + (NN-2)*Max[i3] < sumNNMax:
                        break
                    flag = 0
                    for i4 in range(i3,setsN):
                        if Max[i1] + Max[i2] + Max[i3] + Max[i4] < sumNNMax:
                            break
                        for x1 in sets[i1]:
                            for x2 in sets[i2]:
                                for x3 in sets[i3]:
                                    for x4 in sets[i4]:
                                        if len(set(x1))+len(set(x2))+len(set(x3))+len(set(x4)) > colM:
                                            continue
                                        if len(set(x1) & set(x2)) == 0 and len(set(x1) & set(x3)) == 0 and len(set(x1) & set(x4)) == 0 and len(set(x2) & set(x3)) == 0 and len(set(x2) & set(x4)) == 0 and len(set(x3) & set(x4)) == 0:
                                            sum = Max[i1] + Max[i2] + Max[i3] + Max[i4]
                                            if sum_max < sum:
                                                sum_max = sum
                                                flag = 1
                                                print(sum)
                                                print(i1,i2,i3,i4)
                                                print(x1)
                                                print(x2)
                                                print(x3)
                                                print(x4)
                                                print("************************")
                                                if i1 < setsN-1:
                                                    if sum_max > NN*Max[i1+1]:
                                                        return sum_max
                                                    break 
                                                break
                                    if flag == 1:
                                        break        
                                if flag == 1:
                                    break
                            if flag == 1:
                                break
                        if flag == 1:
                            break
    if NN == 5:
        for i1 in range(0,setsN):
            if NN * Max[i1] < sumNNMax:
                return sum_max 
            for i2 in range(i1,setsN):
                if Max[i1] + (NN-1)*Max[i2] < sumNNMax:
                    break
                for i3 in range(i2,setsN):
                    if Max[i1] + Max[i2] + (NN-2)*Max[i3] < sumNNMax:
                        break
                    for i4 in range(i3,setsN):
                        if Max[i1] + Max[i2] + Max[i3] + (NN-3)*Max[i4] < sumNNMax:
                            break
                        flag = 0
                        for i5 in range(i4,setsN):
                            if Max[i1] + Max[i2] + Max[i3] + Max[i4] + Max[i5] < sumNNMax:
                                break
                            for x1 in sets[i1]:
                                for x2 in sets[i2]:
                                    for x3 in sets[i3]:
                                        for x4 in sets[i4]:
                                            for x5 in sets[i5]:
                                                if len(set(x1))+len(set(x2))+len(set(x3))+len(set(x4))+len(set(x5)) > colM:
                                                    continue
                                                if len(set(x1) & set(x2)) == 0 and len(set(x1) & set(x3)) == 0 and len(set(x1) & set(x4)) == 0 and len(set(x1) & set(x5)) == 0 and len(set(x2) & set(x3)) == 0 and len(set(x2) & set(x4)) == 0 and len(set(x2) & set(x5)) == 0 and len(set(x3) & set(x4)) == 0 and len(set(x3) & set(x5)) == 0 and len(set(x4) & set(x5)) == 0:
                                                    sum = Max[i1] + Max[i2] + Max[i3] + Max[i4] + Max[i5]
                                                    if sum_max < sum:
                                                        sum_max = sum
                                                        flag = 1
                                                        print(sum)
                                                        print(i1,i2,i3,i4,i5)
                                                        print(x1)
                                                        print(x2)
                                                        print(x3)
                                                        print(x4)
                                                        print(x5)
                                                        print("************************")
                                                        if i1 < setsN-1:
                                                            if sum_max > NN*Max[i1+1]:
                                                                return sum_max
                                                            break
                                                        break
                                            if flag == 1:
                                                    break
                                        if flag == 1:
                                            break
                                    if flag == 1:
                                        break
                                if flag == 1:
                                    break
                            if flag == 1:
                                break
    return sum_max
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
M1 =np.array([
    [ 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1 ],
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
]);
M = M1
M_row = M.shape[0]
M_col = M.shape[1]

# sparsity = 1 #1的比例
# M_row = 36
# M_col = 36
# M = generate_random_matrix(M_row, M_col, sparsity)

# 统计M中1的个数
plus_equal_count = np.sum(M)-M_row
print("不优化结果: ")
print("+=个数: ", plus_equal_count)
#统计M对角线上的1的个数
equal_count = M_row - np.sum(np.diag(M))
print("=个数: ", equal_count)
Index = IndexOfoneInMatrix(M, M_row, M_col)

# 初始化集合

sets = {2: [], 3: [], 4: [], 5: [], 6: []}#, 7: [], 8: [], 9:[], 10:[], 11:[]}
Max = []
setsN = 7
# 测试不同的 N 值
for N in range(2, setsN):
    max_index_len, max_combinations = find_max_submatrix(Index, N)
    if max_index_len == 0:
        break
    print(f"N = {N}")
    Max.append(max_index_len)
    print(f"最大子矩阵大小: {max_index_len}")
    num = 0
    for combination in max_combinations:
        #print(f"行组合: {combination}")
        intersection = set(Index[combination[0]])
        for row in combination[1:]:
            intersection &= set(Index[row])
        #print(f"列组合: {list(intersection)}")
        add_to_set(N, intersection, num, sets)
        num += 1

# 获取排序后的下标（递减排序）
index = sorted(range(len(Max)), key=lambda i: Max[i], reverse=True)
sortMax = sorted(Max, reverse=True)
sortsets = {}
num =0
for i in index:
    sortsets.update({num: sets[i+2]})
    num += 1
# 检查交集
sum2Max = sortMax[0] + sortMax[1]
print("sum2max",sum2Max)
NN = 2
sum2Max = check_intersections(sortsets,setsN-2,NN,M_col,sortMax,sum2Max)
print("sum2max",sum2Max)
NN = 3
sum3Max = check_intersections(sortsets,setsN-2,NN,M_col,sortMax,sum2Max)
print("sum3max",sum3Max)
NN = 4
sum4Max = check_intersections(sortsets,setsN-2,NN,M_col,sortMax,sum3Max)
print("sum4max",sum4Max)
NN = 5
sum5Max = check_intersections(sortsets,setsN-2,NN,M_col,sortMax,sum4Max)
print("sum5max",sum5Max)
#这个函数可以优化一下，设置终止条件，比如最大值一定小于某个值