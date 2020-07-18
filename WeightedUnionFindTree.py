class WeightedUnionFind():
    """
    https://qiita.com/drken/items/cce6fc5c579051e64fab
    UnionFindTreeに親ノードへ向かう辺の重みを追加したもの。
    これによりノードx, yがpath-connectedな場合にその間の符号付き距離をO(logN)で取得できる。
    ただしこれは「有向辺」であることに注意。
    というのも距離の導出は各ノードから根ノードまでの距離を取ってその差分で求めているため、
    無向辺の場合はLCAを求めないといけない。そしてUnionFindでは経路圧縮のせいでLCAを求めるのは不可能。
    """
    def __init__(self, n):
        """
        要素数をnとして、各ノードを0,1,...,(n-1)の番号で管理する。
        parentsは各ノードの属する木の根を表す。
        ただし根ノードのparentには(その木のノード数)*(-1)を格納する。
        diff_weightは各ノードから親ノードへ向かう辺の重みを格納する。
        根ノードに格納する重みは0とする。
        """
        self.n = n
        self.parents = [-1] * n
        self.diff_weight = [0] * n

    def find(self, x):
        """
        xの属する木の根を返す。
        このとき同時に経路圧縮して、探索途中のノードを全て根に繋ぎ直す。
        """
        if self.parents[x] < 0:
            return x
        else:
            edge_weight = self.diff_weight[x]
            prev_parent = self.parents[x]
            self.parents[x] = self.find(self.parents[x])
            self.diff_weight[x] = self.diff_weight[prev_parent] + edge_weight
            return self.parents[x]

    def union(self, x, y, w):
        """
        xからyへ重みwの辺を繋ぐ。
        既にx, yがpath-connectedな場合はFalseを返す。
        そうでなければTrueを返す。
        """
        rx = self.find(x)
        ry = self.find(y)
        w += self.diff_weight[y] - self.diff_weight[x]
        w *= -1 # 内部処理的にはfind(y) -> find(x)へ辺を繋ぐ

        if rx == ry:
            return False

        if self.parents[rx] > self.parents[ry]:
            rx, ry = ry, rx # xの方が木のサイズが大きいようにする
            w *= -1 # それに応じて結ぶ辺の向きも反転させる

        self.parents[rx] += self.parents[ry]
        self.parents[ry] = rx
        self.diff_weight[ry] = w
        return True

    def weight(self, x):
        """
        xからxの属する木の根ノードまでの重みを求める
        """
        _ = self.find(x)
        return self.diff_weight[x]

    def diff(self, x, y):
        """
        xからyへ向かうのにかかる重みを求める
        """
        return self.weight(x) - self.weight(y)

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
