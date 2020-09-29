class SegTree():
    """
    0-indexedの配列[a0, a1, a2, ..., a(N-1)]に対して以下のクエリをそれぞれO(logx)で行う：
        1. i番目の要素にxを代入
        2. 半開区間[i, j)内の最小値を返す
    self.id_elemとself.funcを変えることでRmQ以外にも例えばRMQに対応可能。
    参考：https://juppy.hatenablog.com/entry/2019/05/02/蟻本_python_セグメント木_競技プログラミング_Atcoder
    """
    def __init__(self, N):
        #####identity element######
        self.id_elem = 0
        self.func = max

        #num_max: the smallest power of two over N
        self.num_max = 2 ** ((N - 1).bit_length())
        self.x = [self.id_elem] * (2 * self.num_max - 1)

    def get_elem(self, i):
        """
        i番目の要素を返す
        """
        return self.x[self.num_max + i - 1]

    def update(self, i, x):
        """
        i番目の要素にxを代入
        """
        i += self.num_max - 1
        self.x[i] = x
        while (i > 0):
            i = (i - 1) // 2
            self.x[i] = self.func(self.x[i * 2 + 1], self.x[i * 2 + 2])

    def query(self, i=0, j=-1):
        """
        半開区間[i, j)内の最小値を返す
        query()で配列全体に対してクエリを実行する
        """
        if j == -1: j = self.num_max
        i += self.num_max - 1
        j += self.num_max - 1
        res = self.id_elem
        
        while i < j:
            if i % 2 == 0:
                res = self.func(res, self.x[i])
                i += 1
            if j % 2 == 0:
                res = self.func(res, self.x[j - 1])
            i = (i - 1) // 2
            j = (j - 1) // 2
        return res