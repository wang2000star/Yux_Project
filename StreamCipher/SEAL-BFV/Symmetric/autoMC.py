def generate_matrix_multiplication_expressions_withPlainMod(n, M,idbool):
    expressions = []
    for i in range(n):
        if idbool == 1:
            expression = f"A[id{i}] = ((("
        else:
            expression = f"A[{i}] = ((("
        terms = []
        for j in range(n):
            if M[i][j] != 0: 
                if M[i][j] == 1:
                    if idbool == 1:
                        term = f"temp[id{j}]"
                    else:
                        term = f"temp[{j}]"
                elif M[i][j] < 0:
                    if idbool == 1:
                        term = f"({M[i][j]})*temp[id{j}]"
                    else:
                        term = f"({M[i][j]})*temp[{j}]"
                else:
                    if idbool == 1:
                       term = f"{M[i][j]}*temp[id{j}]"
                    else:
                       term = f"{M[i][j]}*temp[{j}]"
                terms.append(term)
        expression += " + ".join(terms) + ") % PlainMod)+PlainMod)%PlainMod;"
        expressions.append(expression)
    return expressions

def generate_matrix_multiplication_expressions(n, M,idbool):
    expressions = []
    for i in range(n):
        if idbool == 1:
            expression = f"A[id{i}] = "
        else:
            expression = f"A[{i}] = "
        terms = []
        for j in range(n):
            if M[i][j] != 0:
                if M[i][j] == 1:
                    if idbool == 1:
                        term = f"temp[id{j}]"
                    else:
                        term = f"temp[{j}]"
                elif M[i][j] < 0:
                    if idbool == 1:
                        term = f"({M[i][j]})*temp[id{j}]"
                    else:
                        term = f"({M[i][j]})*temp[{j}]"
                else:
                    if idbool == 1:
                       term = f"{M[i][j]}*temp[id{j}]"
                    else:
                       term = f"{M[i][j]}*temp[{j}]"
                terms.append(term)
        expression += " + ".join(terms) + ";"
        expressions.append(expression)
    return expressions
#0,1，2矩阵
def generate_matrix_multiplication_expressions_withHElib(n, M,idbool):
    expressions = []
    for i in range(n):
        flag = 0
        for jj in range(n):
            if M[i][jj] > 0:
                min1 = jj
                break
        if M[i][i] == 0:
            flag = 1
            if M[i][min1] == 1:
                if idbool == 1:
                    expression = f"eData[id{i}] = temp[id{min1}];"
                else:
                    expression = f"eData[{i}] = temp[{min1}];"
            if M[i][min1] == 2:
                if idbool == 1:
                    expression = f"eData[id{i}] = temp[id{min1}];eData[id{i}] += temp[id{min1}];"
                else:
                    expression = f"eData[{i}] = temp[{min1}];eData[{i}] += temp[{min1}];"
        if M[i][i] == 2:
           if idbool == 1:
                expression = f"eData[id{i}] += temp[id{i}];"
           else:
                expression = f"eData[{i}] += temp[{i}];"
        if M[i][i] != 1:
           expressions.append(expression)   
        Set = list(range(n))
        Set.remove(i)
        if flag == 1:
            Set.remove(min1)
        for j in Set:
            if M[i][j] == 1:
                if idbool == 1:
                    expression = f"eData[id{i}] += temp[id{j}];"
                else:
                    expression = f"eData[{i}] += temp[{j}];"
                expressions.append(expression)
            if M[i][j] == 2:
                if idbool == 1:
                    expression = f"eData[id{i}] += temp[id{j}];eData[id{i}] += temp[id{j}];"
                else:
                    expression = f"eData[{i}] += temp[{j}];eData[{i}] += temp[{j}];"
                expressions.append(expression)
    return expressions
# 示例矩阵 M 和列向量 temp
n = 36
M =[
    [ 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1 ],
    [ 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0 ],
    [ 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1 ],
    [ 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0 ],
    [ 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1 ],
    [ 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1 ],
    [ 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1 ],
    [ 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0 ],
    [ 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1 ],
    [ 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0 ],
    [ 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1 ],
    [ 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1 ],
    [ 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0 ],
    [ 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1 ],
    [ 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1 ],
    [ 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1 ],
    [ 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1 ],
    [ 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1 ],
    [ 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1 ],
    [ 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1 ],
    [ 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1 ],
    [ 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1 ],
    [ 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0 ],
    [ 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1 ],
    [ 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1 ],
    [ 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0 ],
    [ 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0 ],
    [ 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1 ],
    [ 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1 ],
    [ 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1 ],
    [ 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1 ],
    [ 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0 ],
    [ 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1 ],
    [ 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0 ],
    [ 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1 ],
    [ 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0 ]
];
idbool = 0;
expressions1 = generate_matrix_multiplication_expressions(n, M,idbool)
expressions2 = generate_matrix_multiplication_expressions_withPlainMod(n, M,idbool)
expressions3 = generate_matrix_multiplication_expressions_withHElib(n, M,idbool)

with open('autoMC.txt', 'w') as f:
    for expr in expressions1:
        f.write(expr + '\n')
    f.write('\n')
    for expr in expressions2:
        f.write(expr + '\n')
    f.write('\n')
    for expr in expressions3:
        f.write(expr + '\n')
    
