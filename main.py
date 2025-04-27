from polynomial.multivar_sym_poly import MultiVarSymPloy
from matrix.jordan_form import jordan_form
from utils.input_helper import input_matrix, input_polynomial
from sympy import pprint

def main():
    print("请选择要进行的操作：")
    print("1. 多项式操作（输入示例，乘号：*；幂：**;变量：x_1, x_2, y）")
    print("2. 矩阵Jordan分解")
    choice = input("请输入数字 (1 或 2): ").strip()

    if choice == "1":
        print("\n=== 多项式操作 ===")
        n = int(input("请输入变量个数 (n): "))
        order = int(input("请输入多项式最大阶次 (order): "))
        msp = MultiVarSymPloy(n, order)

        print("\n请选择多项式操作类型：")
        print("1. 获取多项式的系数和幂次，例如:2*x1**2*x2**5，返回系数2")
        print("2. 合并同类项")
        print("3. 初等对称多项式基")
        print("4. 初等对称多项式值")
        print("5. 幂和多项式表示")
        print("6. 自定义多项式转初等对称多项式表示")
        poly_choice = input("请输入数字 (1-6): ").strip()

        if poly_choice in {"1", "2", "6"}:
            expr = input_polynomial(msp.var_list)

        if poly_choice == "1":
            result = msp.getPolyCoffAndPow(expr)
            pprint(result)
        elif poly_choice == "2":
            result = msp.getPolyCollect(expr)
            pprint(result)
        elif poly_choice == "3":
            result = msp.all_primary_sym_ploy_base()
            pprint(result)
        elif poly_choice == "4":
            result = msp.primary_sym()
            pprint(result)
        elif poly_choice == "5":
            eq = msp.simplify_to_primary_sym()
            pprint(eq)
        elif poly_choice == "6":
            eq = msp.primary_sym_representation(expr)
            pprint(eq)
        else:
            print("无效选择。")
    elif choice == "2":
        print("\n=== 矩阵 Jordan 分解 ===")
        M = input_matrix()
        P, J = jordan_form(M)

        print("\n输入矩阵 M:")
        pprint(M)
        print("\nJordan 标准型 J:")
        pprint(J)
        print("\n变换矩阵 P:")
        pprint(P)
    else:
        print("无效选择。")

if __name__ == "__main__":
    main()
