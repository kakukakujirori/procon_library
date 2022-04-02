import numpy as np
from numba import njit

@njit("(i8,i8[:],i8[:])", cache=True)
def stronglyConnectedComponents(N: int, edges_from: np.ndarray, edges_to: np.ndarray):
    """
    グラフの強連結成分分解を行う。
    さらにこの時返されるgroupはトポロジカル順になる。
    頂点番号は0スタートであること。

    !!! まだ実装の中身が理解できてない !!!
    """
    M = len(edges_from)
    start = [0] * (N+1)
    elist = [0] * M
    for s in edges_from:
        start[s+1] += 1
    for i in range(1, N+1):
        start[i] += start[i-1]
    counter = start[:]
    for s, t in zip(edges_from, edges_to):
        elist[counter[s]] = t
        counter[s] += 1

    visited = [i for i in range(0)]  # dfsで訪れた順番に頂点を格納、帰りがけにpopする
    low = [0] * N  # ??????
    Ord = [-1] * N  # 最初に訪れた時の時刻を格納
    ids = [0] * N  # 属する連結成分番号を格納
    NG = [0, 0]  # [時刻、連結成分番号]

    for i in range(N):
        if Ord[i] == -1:
            # dfs
            stack = [(i, -1, 0),(i, -1, 1)]  # 順方向と逆方向のdfsをあらかじめ仕込んでおく
            while stack:
                v, p, t = stack.pop()
                if t:  # 順方向
                    if p != -1 and Ord[v] != -1:
                        low[p] = min(low[p], Ord[v])
                        stack.pop()
                        continue
                    low[v], Ord[v] = NG[0], NG[0]
                    NG[0] += 1
                    visited.append(v)
                    for i in range(start[v], start[v+1]):
                        to = elist[i]
                        if Ord[to] == -1:
                            stack.append((to, v, 0))
                            stack.append((to, v, 1))
                        else:
                            low[v] = min(low[v], Ord[to])
                else:  # 逆方向
                    if low[v] == Ord[v]:
                        while True:
                            u = visited.pop()
                            Ord[u] = N
                            ids[u] = NG[1]
                            if u == v:
                                break
                        NG[1] += 1
                    low[p] = min(low[p], low[v])
    
    # 連結成分の番号を昇順にしてgroupsがトポロジカル順になるよう変更
    ids = [NG[1] - 1 - x for x in ids]

    groups = [[i for i in range(0)] for _ in range(NG[1])]
    for i in range(N):
        groups[ids[i]].append(i)

    return groups


if __name__ == '__main__':
    """AtCoder Library Practice Contest G-SCC
    https://atcoder.jp/contests/practice2/tasks/practice2_g?lang=ja"""
    import sys

    N, M = np.fromstring(sys.stdin.readline(), dtype=int, sep=' ')
    edges_from, edges_to = np.fromstring(sys.stdin.read(), dtype=int, sep=' ').reshape(M, 2).T
    ans = stronglyConnectedComponents(N, edges_from, edges_to)
    print(len(ans))
    for line in ans:
        print(*([len(line)] + line))