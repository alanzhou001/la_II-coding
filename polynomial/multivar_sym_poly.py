from sympy import symbols, Poly, simplify, Mul, Eq
import itertools

class MultiVarSymPloy:
    def __init__(self, n, order):
        self.n = n
        self.order = order
        self.var_list = symbols(f'x_1:{n+1}')
        self.sigma = symbols(f'\sigma_1:{n+1}')

    def getPolyCoffAndPow(self, poly_expr):
        '''给定一个多元多项式，返回系数以及各元的幂次'''
        alist = []
        for term in poly_expr.args: # 遍历多项式的每一项
            # 对每一项调用 return_cof 方法，获取系数和幂次
            coff, power_list = self.return_cof(term)
            alist.append([coff, power_list])
        return alist

    def return_cof(self, expr):
        '''返回多元多项式单项的系数,系数不能包含变元'''
        # 获取每个变量的幂次
        power_list = [Poly(expr, v).degree() for v in self.var_list]
        term = [v**power for power, v in zip(power_list, self.var_list)]
        coff = simplify(expr / Mul(*term[:len(term)]))
        return coff, power_list

    def base_expr(self, power):
        '''返回各个变元幂次的乘积'''
        # 计算每个变量的幂次，并返回它们的乘积
        term = [v**p for p, v in zip(power, self.var_list)]
        return Mul(*term[:len(term)])

    def getPolyCollect(self, poly_expr):
        '''合并同类项'''
        # 将多项式展开并获取系数和幂次
        poly_expr = poly_expr.expand()
        plist = self.getPolyCoffAndPow(poly_expr)

        # 将系数和幂次转换为字典，并按幂次排序
        dict_list = [{"coff": term[0], "base": term[1]} for term in plist]
        dict_list.sort(key=lambda x: x['base'])
        grouped = itertools.groupby(dict_list, key=lambda x: x['base'])

        # 计算每组的系数和幂次的乘积
        qlist = []
        for base, group in grouped:
            sum_ = sum([item["coff"] for item in group])
            ppp = sum_ * self.base_expr(base)
            qlist.append(ppp)
        return sum(qlist)

    def all_primary_sym_ploy_base(self):
        '''返回n元order次初等对称多项式的基'''
        # 生成从0到order的整数列表
        arr = [i for i in range(0, self.order+1)]
        alist = []
        result = itertools.product(arr, repeat=self.n)

        # 遍历所有可能的组合
        for term in result:
            sum_ = sum([k*v for k, v in zip(term, range(1, self.n+1))])
            if sum_ <= self.order:
                expr = Mul(*[b**p for b, p in zip(self.sigma, term)])
                if expr != 1:
                    alist.append(expr)
        return alist

    def primary_sym(self):
        '''返回初等对称多项式的值'''
        # 计算初等对称多项式的值
        x = symbols("x")
        term = [(x + var_) for var_ in self.var_list]
        expr = Mul(*term[:len(term)]).expand()
        ploy = Poly(expr, x)
        return ploy.all_coeffs()[1:]

    def get_dicts(self):
        '''返回变量和初等对称多项式的映射字典'''
        # 获取变量和初等对称多项式的映射关系
        lhs = self.sigma
        rhs = self.primary_sym()
        return dict(zip(lhs, rhs))

    def simplify_to_primary_sym(self):
        '''将多项式简化为初等对称多项式的表示'''
        dicts = self.get_dicts()
        degree = self.order
        origin_expr = sum([x**degree for x in self.var_list])
        alist = self.all_primary_sym_ploy_base()
        ks = [symbols(f'k{i}') for i in range(len(alist))]
        expr = sum([s*k for s, k in zip(alist, ks)])
        expr1 = expr.subs(dicts) - origin_expr
        wlist = self.getPolyCollect(expr1)
        coff_list = []

        # 获取多项式的系数
        for expr_ in wlist.args:
            power_list = [Poly(expr_, v).degree() for v in self.var_list]
            term = [v**power for power, v in zip(power_list, self.var_list)]
            coff = simplify(expr_ / Mul(*term[:len(term)]))
            coff_list.append(coff)

        from sympy import linsolve
        rlist = list(set(coff_list))
        sols = linsolve(rlist, ks)

        result_subs = {}
        for sol in sols:
            for k, sol0 in zip(ks, sol):
                result_subs[k] = sol0

        lhs = origin_expr
        rhs = expr.subs(result_subs)
        return Eq(lhs, rhs, evaluate=False)

    def primary_sym_representation(self, origin_expr):
        '''将给定的多项式转换为初等对称多项式的表示'''
        dicts = self.get_dicts()
        alist = self.all_primary_sym_ploy_base()
        ks = [symbols(f'k{i}') for i in range(len(alist))]
        expr = sum([s*k for s, k in zip(alist, ks)])
        expr1 = expr.subs(dicts) - origin_expr
        wlist = self.getPolyCollect(expr1)
        coff_list = []

        for expr_ in wlist.args: # 遍历多项式的每一项
            # 对每一项调用 return_cof 方法，获取系数和幂次
            power_list = [Poly(expr_, v).degree() for v in self.var_list]
            term = [v**power for power, v in zip(power_list, self.var_list)]
            coff = simplify(expr_ / Mul(*term[:len(term)]))
            coff_list.append(coff)

        from sympy import linsolve
        rlist = list(set(coff_list))
        sols = linsolve(rlist, ks)

        result_subs = {}
        for sol in sols:
            for k, sol0 in zip(ks, sol):
                result_subs[k] = sol0

        lhs = origin_expr
        rhs = expr.subs(result_subs)
        return Eq(lhs, rhs, evaluate=False)
