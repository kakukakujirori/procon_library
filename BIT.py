class BIT:

    """
    https://tjkendev.github.io/procon-library/python/range_query/bit.html
    Binary index treeの実装
    配列[a1, a2,...,an]に対して以下のクエリをO(logn)で行う:
        1. aiにxを加える
        2. 区間和 ai + a(i+1) + ... + aj の和を求める
    isom法を使えば、
        1. aiの値を取得する
        2. 区間[i, j]の全ての数にxを加算する
    もO(logN)で行うことができる
    """

    def __init__(self, n):
        """
        添字は1スタート
        """
        self.n = n
        self.data = [0] * (n + 1)
        self.el = [0] * (n + 1)

    def add(self, i, x):
        """
        i>0に対してaiにxを加算(x < 0でもOK)
        """
        if i <= 0 or self.n < i:
            print("i should be within 1 to n")
        else:
            self.el[i] += x
            while i <= self.n:
                self.data[i] += x
                i += i & -i

    def sum(self, i):
        """
        添字1からiまでの累積和を求める
        """
        s = 0
        while i > 0:
            s += self.data[i]
            i -= i & -i # i $ (-i)でiの最下位ビットのみ立った値を得る
        return s

    def get(self, i, j=None):
        """
        添字iからjまでの累積和を求める
        j=Noneの場合はaiの値を返す
        """
        if j is None:
            return self.el[i]
        return self.sum(j) - self.sum(i - 1)
