import itertools

def IndexOfoneInMatrix(matrix, m, n):
    # 定义元胞数组
    Index = [[] for _ in range(m)]
    for i in range(m):
        for j in range(n):
            if matrix[i][j] == 1 and j != i:
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

def check_intersections(sets,NN,colM,Max):
    if NN == 2:
        sum_max=0
        for i1 in range(2,8-(NN-1)):
            for i2 in range(i1,8-(NN-2)):
                flag = 0
                for x1 in sets[i1]:
                    for x2 in sets[i2]:
                        if len(set(x1))+len(set(x2)) > colM:
                            continue
                        if len(set(x1) & set(x2)) == 0:
                            sum = Max[i1-2] + Max[i2-2]
                            if sum_max < sum:
                                sum_max = sum
                                flag = 1
                                print(sum)
                                print(i1,i2)
                                print(x1)
                                print(x2)
                                print("************************")
                                break
                    if flag == 1:
                        break
    if NN == 3: 
        sum_max = 0 
        for i1 in range(2,8-(NN-1)):
            for i2 in range(i1,8-(NN-2)):
                for i3 in range(i2,8-(NN-3)):
                    flag = 0
                    for x1 in sets[i1]:
                        for x2 in sets[i2]:
                            for x3 in sets[i3]:
                                if len(set(x1))+len(set(x2))+len(set(x3)) > colM:
                                    continue
                                if len(set(x1) & set(x2)) == 0 and len(set(x1) & set(x3)) == 0 and len(set(x2) & set(x3)) == 0:     
                                    sum = Max[i1-2] + Max[i2-2] + Max[i3-2]
                                    if sum_max < sum:
                                        sum_max = sum
                                        flag = 1
                                        print(sum)
                                        print(i1,i2,i3)
                                        print(x1)
                                        print(x2)
                                        print(x3)
                                        print("************************")
                                    break
                            if flag == 1:
                                break
                        if flag == 1:
                            break
    if NN == 4:
        for i1 in range(2,8-(NN-1)):
            for i2 in range(2,8-(NN-2)):
                for i3 in range(2,8-(NN-3)):
                    for i4 in range(2,8-(NN-4)):
                        for x1 in sets[i1]:
                            for x2 in sets[i2]:
                                for x3 in sets[i3]:
                                    for x4 in sets[i4]:
                                        if len(set(x1))+len(set(x2))+len(set(x3))+len(set(x4)) > colM:
                                            continue
                                        if len(set(x1) & set(x2)) == 0 and len(set(x1) & set(x3)) == 0 and len(set(x1) & set(x4)) == 0 and len(set(x2) & set(x3)) == 0 and len(set(x2) & set(x4)) == 0 and len(set(x3) & set(x4)) == 0:
                                            print(x1)
                                            print(x2)
                                            print(x3)
                                            print(x4)
                                            print("************************")
                                            return

M = [
    [ 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1 ],
    [ 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0 ],
    [ 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0 ],
    [ 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0 ],
    [ 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1 ],
    [ 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0 ],
    [ 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1 ],
    [ 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1 ],
    [ 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1 ],
    [ 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1 ],
    [ 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0 ],
    [ 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1 ]
]

RowM = 12
ColM = 12
Index = IndexOfoneInMatrix(M, RowM, ColM)

# 初始化集合
sets = {2: [], 3: [], 4: [], 5: [], 6: [], 7: []}
Max = []

# 测试不同的 N 值
for N in range(2, RowM + 1):
    max_index_len, max_combinations = find_max_submatrix(Index, N)
    if max_index_len == 0:
        break
    print(f"N = {N}")
    Max.append(max_index_len)
    print(f"最大子矩阵大小: {max_index_len}")
    num = 0
    for combination in max_combinations:
        print(f"行组合: {combination}")
        intersection = set(Index[combination[0]])
        for row in combination[1:]:
            intersection &= set(Index[row])
        print(f"列组合: {list(intersection)}")
        add_to_set(N, intersection, num, sets)
        num += 1

# 检查交集

NN = 3
flag = 0

check_intersections(sets,NN,ColM,Max)