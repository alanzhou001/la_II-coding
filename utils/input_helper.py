from sympy import Matrix, symbols, sympify

def input_matrix():
    print("请输入矩阵的行数:")
    rows = int(input("行数: "))
    print("请输入矩阵的列数:")
    cols = int(input("列数: "))

    print(f"请输入矩阵元素，按行输入，用空格隔开，每行回车:")
    data = []
    for i in range(rows):
        row = list(map(eval, input(f"第{i+1}行: ").strip().split()))
        if len(row) != cols:
            print(f"错误：第{i+1}行输入了{len(row)}列，应该是{cols}列。重新输入。")
            return input_matrix()
        data.append(row)
    return Matrix(data)

def input_polynomial(var_list):
    print(f"请输入以变量 {', '.join(str(v) for v in var_list)} 为变量的多项式表达式：")
    expr_str = input("多项式: ")
    expr = sympify(expr_str)
    return expr
