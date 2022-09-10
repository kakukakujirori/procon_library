class LCA():
    """
    木の最小共通祖先(Least Common Ancestor)をダブリングで求める。
    頂点番号が0, 1, ..., n-1で与えられた木のrepnを入力としてインスタンス化、2頂点を入力してLCAを返す。
    """
    def __init__(self, repn, root=0):
        import sys
        sys.setrecursionlimit(10**8)
        assert root >= 0
        self.n = len(repn)
        self.repn = repn
        self.log_depth = len(bin(self.n - 1)[2:])
        self.parent = None  # [[-1] * self.n for _ in range(self.log_depth)]  # parent[k][v] = (頂点vから2^k上の親)
        self.depth = None  # [-1] * self.n
        self.root = root
        self.reroot(root)
        
    def reroot(self, r=None):
        """
        Reroot the tree from r.
        If r is None, r is set to be an endpoint of a diameter of the tree,
        which can be useful when accounting for a diameter of the tree.
        """
        if r is None:
            self.root = self.farthest_vertex()
        
        # initialize
        self.depth = [-1] * self.n
        self.parent = [[-1] * self.n for _ in range(self.log_depth)]
        self.depth[self.root] = 0

        # dfs
        def dfs(v, p, d):
            for nv in self.repn[v]:
                if nv == p: continue
                self.parent[0][nv] = v
                self.depth[nv] = d + 1
                dfs(nv, v, d + 1)
        
        dfs(self.root, -1, 0)

        # doubling preprocess
        for k in range(self.log_depth - 1):
            for v in range(self.n):
                if self.parent[k][v] < 0:
                    self.parent[k + 1][v] = -1
                else:
                    self.parent[k + 1][v] = self.parent[k][self.parent[k][v]]

    def farthest_vertex(self):
        """
        Return a farthest vertex from self.root
        This can be useful when detecting a diameter (p, q) by calling
        self.reroot(), and then (p, q) = (lca.root, lca.farthest_vertex())
        """
        v = 0
        max_depth = -1
        for i, x in enumerate(self.depth):
            if x > max_depth:
                v = i
                max_depth = x
        assert max_depth >= 0
        return v

    def __call__(self, u, v):
        """
        Return LCA of u and v
        """
        if self.depth[u] < self.depth[v]: u, v = v, u
        for k in range(self.log_depth):
            if ((self.depth[u] - self.depth[v]) >> k & 1):
                u = self.parent[k][u]

        if u == v: return u
        for k in range(self.log_depth - 1, -1, -1):
            if self.parent[k][u] != self.parent[k][v]:
                u = self.parent[k][u]
                v = self.parent[k][v]

        return self.parent[0][u]