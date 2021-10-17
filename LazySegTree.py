import numpy as np
from numba import jitclass, i8

spec = [
    ('N', i8),
    ('id_elem_X', i8),
    ('id_elem_A', i8),
    ('X', i8[:]),
    ('A', i8[:]),
]

@jitclass(spec)
class LazySegTree():
    """
    0-indexedの配列[x0, x1, ..., x(N-1)]に対して以下のクエリをそれぞれO(logx)で行う：
        1. xiをxに更新
        2. 半開区間[i, j)にモドイド間の作用f:A->Xをかける
        3. 半開区間[i, j)の積(xi * ... * x(j-1))を返す
    参考：https://maspypy.com/segment-tree-のお勉強2、ACL Beginner Contest E
    注意：セグ木自体の配列は1-indexed
    """
    # >>> SET YOURSELF >>>
    def op_X(self, x, y): return min(x, y)
    
    def op_A(self, a, b): return a + b
 
    def act(self, x, a): return x + a
    # <<< SET YOURSELF <<<

    def __init__(self, N):
        # >>> SET YOURSELF >>>
        self.id_elem_X = 1 << 60
        self.id_elem_A = 0
        # <<< SET YOURSELF <<<

        # Nは2べきでなくても全く同様に動作する
        self.N = N
        self.X = np.full((2 * N,), self.id_elem_X, dtype=np.int64)
        self.A = np.full((2 * N,), self.id_elem_A, dtype=np.int64)
    
    def _eval_at(self, i):
        return self.act(self.X[i], self.A[i])
    
    def _propagate_at(self, i):
        # 自ノードの溜め込んだ作用をかける
        self.X[i] = self._eval_at(i)
        # 子ノードに作用を伝搬
        self.A[i << 1] = self.op_A(self.A[i << 1], self.A[i])
        self.A[i << 1 | 1] = self.op_A(self.A[i << 1 | 1], self.A[i])
        # 自ノードの溜め込んだ作用を初期化
        self.A[i] = self.id_elem_A

    def _propagate_above(self, i):
        bitlen = 0
        while pow(2, bitlen) <= i: bitlen += 1
        H = bitlen - 1  # 根まで遡る回数
        for h in range(H, 0, -1):
            self._propagate_at(i >> h)  # 根から順に作用の伝搬を実行

    def _recalc_above(self, i):
        while i > 1:
            # 親ノードに移って子ノード2つの値を統合
            i >>= 1
            self.X[i] = self.op_X(self._eval_at(i << 1), self._eval_at(i << 1 | 1))

    def build(self, seq):
        """
        配列を入力としてセグ木の初期化を行う
        """
        for i, x in enumerate(seq, self.N):
            self.X[i] = x
        for i in range(self.N - 1, 0, -1):
            self.X[i] = self.op_X(self.X[i << 1], self.X[i << 1 | 1])

    def update(self, i, x):
        """
        i番目の要素をxに更新
        """
        i += self.N
        self._propagate_above(i)
        self.X[i] = x
        self.A[i] = self.id_elem_A
        self._recalc_above(i)

    def action(self, i, j, a):
        """
        半開区間[i, j)に作用をかける
        """
        i += self.N
        j += self.N
        i0 = i // (i & -i)  # 奇数になるまでiを2で割ったもの
        j0 = j // (j & -j)  # 奇数になるまでjを2で割ったもの - 1
        self._propagate_above(i0)
        self._propagate_above(j0)
        while i < j:
            if i & 1:
                self.A[i] = self.op_A(self.A[i], a)
                i += 1
            if j & 1:
                j -= 1
                self.A[j] = self.op_A(self.A[j], a)
            i >>= 1
            j >>= 1
        self._recalc_above(i0)
        self._recalc_above(j0)
    
    def get_elem(self, i):
        """
        i番目の要素を返す
        """
        return self.mul(i, i+1)

    def mul(self, i, j):
        """
        半開区間[i, j)の積を返す
        """
        i += self.N
        j += self.N
        self._propagate_above(i // (i & -i))
        self._propagate_above(j // (j & -j))
        vL = self.id_elem_X
        vR = self.id_elem_X
        while i < j:
            if i & 1:
                vL = self.op_X(vL, self._eval_at(i))
                i += 1
            if j & 1:
                j -= 1
                vR = self.op_X(self._eval_at(j), vR)
            i >>= 1
            j >>= 1
        return self.op_X(vL, vR)
    
    def debug(self):
        for i in range(1, self.N * 2):
            x, a = self.X[i], self.A[i]
            print("i = {:2}, (x={} / a={})".format(i, x, a))