import numpy as np
from numba import jitclass, i8


spec = [
    ('n', i8),
    ('mod', i8),
    ('data', i8[:]),
    ('el', i8[:])
]

@jitclass(spec)
class BIT:
    """
    https://tjkendev.github.io/procon-library/python/range_query/bit.html
    Binary index treeの実装
    1-indexedの配列[a1, a2,...,an]に対して以下のクエリをO(logn)で行う:
        1. aiにxを加える
        2. 区間和 ai + a(i+1) + ... + aj の和を求める
    isom法を使えば、
        1. aiの値を取得する
        2. 区間[i, j]の全ての数にxを加算する
    もO(logN)で行うことができる
    """

    def __init__(self, n: int, mod: int = -1):
        """
        添字は1スタート、modは必要であれば設定しないと内部でオーバーフローする。
        """
        self.n = n
        self.mod = mod
        self.data = np.zeros(n + 1, dtype=np.int64)
        self.el = np.zeros(n + 1, dtype=np.int64)

    def add(self, i: int, x: int):
        """
        i>0に対してaiにxを加算(x < 0でもOK)
        """
        if i <= 0 or self.n < i:
            raise ValueError("i should be within 1 to n")
        else:
            self.el[i] += x
            if self.mod != -1: self.el[i] %= self.mod
            while i <= self.n:
                self.data[i] += x
                if self.mod != -1:
                    self.data[i] %= self.mod
                i += i & -i

    def sum(self, i: int = -1):
        """
        a1+a2+...+aiを求める。
        i=Noneの場合はj=self.nとして計算される。
        """
        if i < 0: i = self.n
        s = 0
        while i > 0:
            s += self.data[i]
            if self.mod != -1: s %= self.mod
            i -= i & -i # i&(-i)でiの最下位ビットのみ立った値を得る
        return s

    def get(self, i: int, j: int = -1):
        """
        ai+a(i+1)+...+ajを求める。
        j=Noneの場合はj=self.nとして計算される。
        """
        if j < 0: j = self.n
        assert i <= j
        return self.sum(j) - self.sum(i - 1)

    def binary_search(self, x: int):
        """
        a1+a2+...+ai <= x < a1+a2+...+ai+a(i+1)となるiを求める。
        BITをmultisetの代替として使用する場合に使える。
        数列は正の数からのみ成ることが前提。必要であれば下2行のチェックを省いて高速化すること。
        """
        for i in range(1, self.n):
            assert self.get(i, i) >= 0

        if self.get(1, self.n) <= x:
            return self.n + 1
        left = 1  # a1+...+a(left) <= x
        right = self.n  # x < a1+...+a(right)
        while right - left > 1:
            mid = (left + right) // 2
            if self.get(1, mid) <= x:
                left = mid
            else:
                right = mid
        return left

    def debug(self):
        """
        BITが仮想的に見ている配列を返す
        """
        return [self.get(i) for i in range(1, self.n + 1)]