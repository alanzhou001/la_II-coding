from sympy import Matrix

def jordan_form(M):
    n = M.rows
    eigs = M.eigenvals()

    def eig_mat(lam, k):
        """计算(M - lam*I)**k"""
        return (M - lam * Matrix.eye(n))**k

    def nullity_chain(lam, alg_mult):
        """计算零空间维数链"""
        nulls = [0]
        k = 1
        while True:
            nullity = n - eig_mat(lam, k).rank()
            nulls.append(nullity)
            if nullity == alg_mult or nullity == nulls[-2]:
                break
            k += 1
        return nulls

    def blocks_from_nullity_chain(d):
        """计算Jordan块的大小"""
        mids = [2*d[i] - d[i-1] - d[i+1] for i in range(1, len(d)-1)]
        end = [d[-1] - d[-2]] if len(d) > 1 else [d[0]]
        return mids + end

    def pick_vec(small, big):
        """选择一个向量，使得它与small线性无关且与big线性相关"""
        if not small:
            return big[0]
        for v in big:
            combined = Matrix.hstack(*(small + [v]))
            if combined.rank() > Matrix.hstack(*small).rank():
                return v

    # Jordan块结构
    block_structure = []
    for lam, alg_mult in eigs.items():
        chain = nullity_chain(lam, alg_mult)
        sizes = blocks_from_nullity_chain(chain)
        for size, count in enumerate(sizes, start=1):
            block_structure.extend([(lam, size)] * count)

    # 构建Jordan标准型J
    blocks = []
    for lam, size in block_structure:
        B = Matrix.zeros(size)
        for i in range(size):
            B[i, i] = lam
            if i < size-1:
                B[i, i+1] = 1
        blocks.append(B)
    J = Matrix.diag(*blocks)

    # 构建变换矩阵P
    P_cols = []
    for lam, size in block_structure:
        null_big = eig_mat(lam, size).nullspace()
        null_small = eig_mat(lam, size-1).nullspace() if size > 1 else []
        v0 = pick_vec(P_cols + null_small, null_big)
        chain_vecs = [v0]
        for k in range(1, size):
            chain_vecs.append((M - lam*Matrix.eye(n)) * chain_vecs[-1])
        P_cols.extend(chain_vecs[::-1])

    P = Matrix.hstack(*P_cols)
    return P, J
