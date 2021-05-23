def EulerTour(X, i0=0):
    n = len(X)
    """
    「部分木」＝「オイラーツアーの区間」は典型！
    頂点vを根とする部分木は、ET[ET1[v]:(ET2[v] + 1)]に対応する。
    応用面ではオイラーツアーで訪れた順に頂点番号をrenumberingするのが良い？
    その場合元の頂点番号との対応はETを使うことで可能。

    注意：以下はスタックを用いた非再帰処理で記述されてる。再帰だと４行ほどで書ける。
    参考：https://qiita.com/Kiri8128/items/2b0023bed9af642c751c

    例：
        ET, ET1, ET2 = EulerTour(X, 0)
        print("ET =", ET) # Pathのi番目の頂点番号
        print("ET1 =", ET1) # Start
        print("ET2 =", ET2) # End
        print("DE =", DE) # Depth
    """
    P = [-1] * n
    Q = [~i0, i0]
    cnt = -1
    ET = []
    ET1 = [0] * n
    ET2 = [0] * n
    DE = [0] * n
    depth = -1
    while Q:
        i = Q.pop()
        if i < 0:
            # ↓ 戻りも数字を足す場合はこれを使う
            # ct += 1
            # ↓ 戻りもETに入れる場合はこれを使う
            # ET.append(P[~i])
            ET2[~i] = cnt
            depth -= 1
            continue
        if i >= 0:
            ET.append(i)
            cnt += 1
            if ET1[i] == 0: ET1[i] = cnt
            depth += 1
            DE[i] = depth
        for a in X[i][::-1]:
            if a != P[i]:
                P[a] = i
                Q.append(~a)
                Q.append(a)
    return ET, ET1, ET2, DE, P
