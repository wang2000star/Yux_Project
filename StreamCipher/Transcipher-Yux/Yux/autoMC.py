def generate_matrix_multiplication_expressions_withPlainMod(n, M):
    expressions = []
    for i in range(n):
        expression = f"A[id{i}] = ((("
        terms = []
        for j in range(n):
            if M[i][j] != 0:
                term = f"{M[i][j]}*temp[id{j}]"
                if M[i][j] == 1:
                    term = f"temp[id{j}]"
                elif M[i][j] < 0:
                    term = f"({M[i][j]})*temp[id{j}]"
                terms.append(term)
        expression += " + ".join(terms) + ") % PlainMod)+PlainMod)%PlainMod;"
        expressions.append(expression)
    return expressions

def generate_matrix_multiplication_expressions(n, M):
    expressions = []
    for i in range(n):
        expression = f"A[id{i}] = "
        terms = []
        for j in range(n):
            if M[i][j] != 0:
                term = f"{M[i][j]}*temp[id{j}]"
                if M[i][j] == 1:
                    term = f"temp[id{j}]"
                elif M[i][j] < 0:
                    term = f"({M[i][j]})*temp[id{j}]"
                terms.append(term)
        expression += " + ".join(terms) + ";"
        expressions.append(expression)
    return expressions
#0,1，-1矩阵
def generate_matrix_multiplication_expressions_withHElib(n, M):
    expressions = []
    for i in range(n):
        expression = f"eData[id{i}] = 0;"  
        expressions.append(expression)   
        for j in range(n):
            if M[i][j] == 1:
                expression = f"eData[id{i}] += temp[id{j}];"
                expressions.append(expression)
            if M[i][j] == -1:
                expression = f"eData[id{i}] -= temp[id{j}];"
                expressions.append(expression)
    return expressions
# 示例矩阵 M 和列向量 temp
n = 12
M = [[ 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1 ],
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
    [ 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1 ]]

expressions = generate_matrix_multiplication_expressions(n, M)
expressions2 = generate_matrix_multiplication_expressions_withPlainMod(n, M)
expressions3 = generate_matrix_multiplication_expressions_withHElib(n, M)

with open('autoMC.txt', 'w') as f:
    for expr in expressions:
        f.write(expr + '\n')
    f.write('\n')
    for expr in expressions2:
        f.write(expr + '\n')
    f.write('\n')
    for expr in expressions3:
        f.write(expr + '\n')
    
