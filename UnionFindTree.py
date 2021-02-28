import numpy as np

class UnionFind():
    """
    https://note.nkmk.me/python-union-find/
    DFSの上位互換と考えて良い
    ２要素x, yがpath-connectedかどうかをほぼ定数オーダーで判定する（螺旋本の14.1節参照）
    さらに連結成分の要素数がO(1)で取得可能なように改造してある
    """
    def __init__(self, n):
        """
        要素数をnとして、各ノードを0,1,...,(n-1)の番号で管理する
        parentsは各ノードの属する木の根を表す
        ただし根ノードのparentには(その木のノード数)*(-1)を格納する
        """
        self.n = n
        self.parents = np.ones(n, dtype=np.int64) * -1

    def find(self, x):
        """
        xの属する木の根を返す
        このとき同時に経路圧縮して、探索途中のノードを全て根に繋ぎ直す
        """
        if self.parents[x] < 0:
            return x
        else:
            self.parents[x] = self.find(self.parents[x])
            return self.parents[x]

    def union(self, x, y):
        """
        x, yのそれぞれ属する木Tx, Tyの根同士を繋ぐ
        このとき木の要素数が小さい方を大きい方に繋ぐ（rankではなくsizeを用いる）
        """
        x = self.find(x)
        y = self.find(y)

        if x == y:
            return

        if self.parents[x] > self.parents[y]:
            x, y = y, x

        self.parents[x] += self.parents[y]
        self.parents[y] = x

    def size(self, x):
        """
        xの属する木の要素数を返す
        根の親を要素数の(-1)倍で定めておいたおかげでO(1)で取得可能
        """
        return -self.parents[self.find(x)]

    def same(self, x, y):
        """
        xとyがpath-connectedかを判定する
        """
        return self.find(x) == self.find(y)

    def members(self, x):
        """
        xの属する木の要素を列挙する
        """
        root = self.find(x)
        return [i for i in range(self.n) if self.find(i) == root]

    def roots(self):
        """
        連結成分の代表元のリストを返す
        """
        return [i for i, x in enumerate(self.parents) if x < 0]

    def group_count(self):
        """
        連結成分の個数を返す
        """
        return len(self.roots())

    def all_group_members(self):
         """
         連結成分およびそれぞれの代表元をまとめた辞書を返す
         代表元がキーになってる
         """
         return {r: self.members(r) for r in self.roots()}

    def __str__(self):
        """
        連結成分およびその代表元を出力
        """
        return '\n'.join('{}: {}'.format(r, self.members(r)) for r in self.roots())
