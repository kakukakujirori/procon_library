# Graph Algorithms
import collections, heapq


class Dijkstra():
    def __init__(self):
        self.e = collections.defaultdict(list)

    def add(self, u, v, d, directed=False):
        """
        add an edge from u to v with cost d
        """
        if directed == False:
            self.e[u].append([v, d])
            self.e[v].append([u, d])
        else:
            self.e[u].append([v, d])

    def delete(self, u, v):
        """
        delete all the edges between u and v
        """
        self.e[u] = [_ for _ in self.e[u] if _[0] != v]
        self.e[v] = [_ for _ in self.e[v] if _[0] != u]

    def Dijkstra_search(self, s):
        """
        s: start point
        return: d[u] = (minimal distance from s to u)
                prev[u] = (previous vertex in the directed path from s to u)
        """
        d = collections.defaultdict(lambda: floar('inf'))
        prev = collections.defaultdict(lambda: None)
        d[s] = 0
        q = []
        heapq.heappush(q, (0, s))
        visited = collections.defaultdict(bool)
        while q:
            k, u = heapq.heappop(q)
            if visited[u]:
                continue
            visited[u] = True

            for uv, ud in self.e[u]:
                if visited[uv]:
                    continue
                vd = k + ud
                if d[uv] > vd:
                    d[uv] = vd
                    prev[uv] = u
                    heapq.heappush(q, (vd, uv))
        return d, prev

    def getDijkstraShortestPath(self, start, goal):
        _, prev = self.Dijkstra_search(start)
        shortestPath = []
        node = goal
        while node is not None:
            shortestPath.append(node)
            node = prev[node]
        return shortestPath[::-1]



class BellmanFord():
    def __init__(self, N):
        self.N = N
        self.edges = []

    def add(self, u, v, d, directed=False):
        """
        add an edge from u to v with cost d
        """
        if directed == False:
            self.edges.append([u, v, d])
            self.edges.append([v, u, d])
        else:
            self.edges.append([u, v, d])

    def BellmanFord_search(self, s):
        """
        s: start point
        return: d[u] = (minimal distance from s to u)
        """
        d = [float('inf')] * self.N
        d[s] = 0
        numEdges = len(self.edges)
        while True:
            update = False
            for i in range(numEdges):
                u, v, cost = self.edges[i]
                if d[u] != float('inf') and d[v] > d[u] + cost:
                    d[v] = d[u] + cost
                    update = True
            if not update:
                break
        return d

    def BellmanFord_negative_bool(self, start, numNodes):
        """
        If True, there exists negative closed loops
        NO GUARANTEE TO WORK PROPERLY!!!
        """
        d = [float('inf')] * self.N
        d[start] = 0
        numEdges = len(self.edges)
        for i in range(numEdges):
            for j in range(numEdges):
                u, v, cost = self.edges[j]
                if d[v] > d[u] + cost:
                    d[v] = d[u] + cost
                    if i == numNodes - 1:
                        return True, d
        return False, d



class Prim():
    """
    NEED TO IMPLEMENT!!!
    Prim's algorithm works only for 'indirected' graphs!
    Also this code is valid only for '0-indexed' graphs!
    """
    def __init__(self, N):
        self.edges = [[] for _ in range(N)]
        self.N = N

    def add(self, u, v, d):
        """
        add an edge from u to v with cost d
        Be careful for the order [cost, to]
        """
        self.edges[u].append([d, v])
        self.edges[v].append([d, u])

    def delete(self, u, v):
        """
        delete all the edges between u and v
        """
        self.edges[u] = [_ for _ in self.edges[u] if _[1] != v]
        self.edges[v] = [_ for _ in self.edges[v] if _[1] != u]

    def Prim(self):
        """
        return = summation of the costs of Minimum Spanning Tree
        """
        not_used = [True] * self.N
        edgelist = []
        for e in self.edges[0]:
            heapq.heappush(edgelist, e)
        not_used[0] = False
        res = 0
        while edgelist:
            minedge = heapq.heappop(edgelist)
            if not not_used[minedge[1]]:
                continue
            v = minedge[1]
            not_used[v] = False
            for e in self.edges[v]:
                if not_used[e[1]]:
                    heapq.heappush(edgelist, e)
            res += minedge[0]
        return res



class WarshallFroyd():
    def __init__(self, N):
        self.N = N
        self.d = [[float('inf')] * N for _ in range(N)]

    def add(self, u, v, c, directed=False):
        """
        Valid only for '0-indexed' graphs!
        add an edge from u to v with cost d
        """
        if directed == False:
            self.d[u][v] = c
            self.d[v][u] = c
        else:
            self.d[u][v] = c

    def WarshallFloyd_search(self, DetectNegativeCycle=False):
        """
        Change self.d to self.d[i][j] = (minimal cost from i to j)
        If self.d[i][j] < 0, this graph has negative loops
        """
        for k in range(self.N):
            for i in range(self.N):
                for j in range(self.N):
                    self.d[i][j] = min(self.d[i][j], self.d[i][k] + self.d[k][j])
        if DetectNegativeCycle == True:
            hasNegativeCycle = False
            for i in range(self.N):
                if self.d[i][i] < 0:
                    hasNegativeCycle = True
                    break
        for i in range(self.N):
            self.d[i][i] = 0
        if DetectNegativeCycle == True:
            return hasNegativeCycle, self.d
        else:
            return self.d



class Kruskal_UnionFind():
    # 無向グラフであるという前提に注意
    def __init__(self, N):
        self.edges = []
        self.rank = [0] * N
        self.par = [i for i in range(N)]
        self.counter = [1] * N

    def add(self, u, v, d):
        """
        u = from, v = to, d = cost
        """
        self.edges.append([u, v, d])

    def find(self, x):
        if self.par[x] == x:
            return x
        else:
            self.par[x] = self.find(self.par[x])
            return self.par[x]

    def unite(self, x, y):
        x = self.find(x)
        y = self.find(y)
        if x != y:
            z = self.counter[x] + self.counter[y]
            self.counter[x], self.counter[y] = z, z
        if self.rank[x] < self.rank[y]:
            self.par[x] = y
        else:
            self.par[y] = x
            if self.rank[x] == self.rank[y]:
                self.rank[x] += 1

    def size(self, x):
        x = self.find(x)
        return self.counter[x]

    def same(self, x, y):
        return self.find(x) == self.find(y)

    def Kruskal(self):
        """
        return: 最小全域木のコストの和
        """
        edges = sorted(self.edges, key=lambda x: x[2])  # costでself.edgesをソートする
        res = 0
        for e in edges:
            if not self.same(e[0], e[1]):
                self.unite(e[0], e[1])
                res += e[2]
        return res
