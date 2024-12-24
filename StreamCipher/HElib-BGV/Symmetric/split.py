def split_matrix(matrix):
    matrix1 = []
    matrix2 = []
    
    for row in matrix:
        new_row1 = []
        new_row2 = []
        for element in row:
            if element == 2:
                new_row1.append(1)
                new_row2.append(1)
            else:
                new_row1.append(element)
                new_row2.append(0)
        matrix1.append(new_row1)
        matrix2.append(new_row2)
    
    return matrix1, matrix2
M = [
[ 0, 0, 1, 2, 1, 1, 2, 0, 2, 2, 0, 1, 0, 1, 2, 2, 1, 0 ],
[ 0, 2, 1, 1, 0, 2, 2, 2, 0, 1, 1, 2, 2, 2, 2, 1, 2, 0 ],
[ 2, 1, 1, 0, 2, 1, 2, 2, 1, 1, 0, 1, 1, 2, 0, 2, 0, 1 ],
[ 2, 1, 0, 0, 0, 1, 2, 1, 1, 2, 0, 2, 2, 0, 1, 0, 1, 2 ],
[ 1, 2, 0, 0, 2, 1, 1, 0, 2, 2, 2, 0, 1, 1, 2, 2, 2, 2 ],
[ 2, 0, 1, 2, 1, 1, 0, 2, 1, 2, 2, 1, 1, 0, 1, 1, 2, 0 ],
[ 0, 1, 2, 2, 1, 0, 0, 0, 1, 2, 1, 1, 2, 0, 2, 2, 0, 1 ],
[ 2, 2, 2, 1, 2, 0, 0, 2, 1, 1, 0, 2, 2, 2, 0, 1, 1, 2 ],
[ 1, 2, 0, 2, 0, 1, 2, 1, 1, 0, 2, 1, 2, 2, 1, 1, 0, 1 ],
[ 2, 0, 1, 0, 1, 2, 2, 1, 0, 0, 0, 1, 2, 1, 1, 2, 0, 2 ],
[ 1, 1, 2, 2, 2, 2, 1, 2, 0, 0, 2, 1, 1, 0, 2, 2, 2, 0 ],
[ 1, 0, 1, 1, 2, 0, 2, 0, 1, 2, 1, 1, 0, 2, 1, 2, 2, 1 ],
[ 2, 0, 2, 2, 0, 1, 0, 1, 2, 2, 1, 0, 0, 0, 1, 2, 1, 1 ],
[ 2, 2, 0, 1, 1, 2, 2, 2, 2, 1, 2, 0, 0, 2, 1, 1, 0, 2 ],
[ 2, 2, 1, 1, 0, 1, 1, 2, 0, 2, 0, 1, 2, 1, 1, 0, 2, 1 ],
[ 2, 1, 1, 2, 0, 2, 2, 0, 1, 0, 1, 2, 2, 1, 0, 0, 0, 1 ],
[ 1, 0, 2, 2, 2, 0, 1, 1, 2, 2, 2, 2, 1, 2, 0, 0, 2, 1 ],
[ 0, 2, 1, 2, 2, 1, 1, 0, 1, 1, 2, 0, 2, 0, 1, 2, 1, 1 ]
];
def write_matrices_to_file(matrix1, matrix2, filename):
    with open(filename, 'w') as f:
        f.write("Matrix 1:\n")
        for row in matrix1:
            f.write('[' + ', '.join(map(str, row)) + '],\n')
        f.write("\nMatrix 2:\n")
        for row in matrix2:
            f.write('[' + ', '.join(map(str, row)) + '],\n')

# 示例使用
M1, M2 = split_matrix(M)
# 将M1，M2输出到split.txt
write_matrices_to_file(M1, M2, 'split.txt')
