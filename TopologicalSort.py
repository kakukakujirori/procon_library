from collections import deque

def topologicalSort(repn):
    # total vertex num
    n = len(repn)

    # count incoming degrees
    indegrees = [0] * n
    for i in range(n):
        for j in repn[i]:
            indegrees[j] += 1
    
    # search start
    sorted_vertices = []
    q = deque([i for i, deg in enumerate(indegrees) if deg == 0])
    while q:
        i = q.popleft()
        sorted_vertices.append(i)
        for j in repn[i]:
            indegrees[j] -= 1
            if indegrees[j] == 0:
                q.append(j)
    
    # success check
    if len(sorted_vertices) == n:
        return sorted_vertices
    else:
        return []

