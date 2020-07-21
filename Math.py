import math


class Combination:
    """
    O(n)の前計算を1回行うことで，O(1)でnCr mod mを求められる
    n_max = 10**6のとき前処理は約950ms (PyPyなら約340ms, 10**7で約1800ms)
    使用例：
    comb = Combination(1000000)
    print(comb(5, 3))  # 10
    """
    def __init__(self, n_max, mod=10**9+7):
        self.mod = mod
        self.modinv = self.make_modinv_list(n_max)
        self.fac, self.facinv = self.make_factorial_list(n_max)
 
    def __call__(self, n, r):
        return self.fac[n] * self.facinv[r] % self.mod * self.facinv[n-r] % self.mod
 
    def make_factorial_list(self, n):
        # 階乗のリストと階乗のmod逆元のリストを返す O(n)
        # self.make_modinv_list()が先に実行されている必要がある
        fac = [1]
        facinv = [1]
        for i in range(1, n+1):
            fac.append(fac[i-1] * i % self.mod)
            facinv.append(facinv[i-1] * self.modinv[i] % self.mod)
        return fac, facinv
 
    def make_modinv_list(self, n):
        # 0からnまでのmod逆元のリストを返す O(n)
        modinv = [0] * (n+1)
        modinv[1] = 1
        for i in range(2, n+1):
            modinv[i] = self.mod - self.mod//i * modinv[self.mod%i] % self.mod
        return modinv


def factorization(n):
    if n <= 1:
        print("Input must be over 2!!!!!!")
        return
    arr = []
    temp = n
    for i in range(2, int(round(n ** 0.5)) + 1):
        cnt = 0
        while temp % i == 0:
            cnt+=1
            temp //= i
        if cnt > 0: arr.append([i, cnt])
    if temp!=1: arr.append([temp, 1])
    if arr==[]: arr.append([n, 1])
    return arr


def get_sieve_of_eratosthenes(n):
    """
    エラトステネスの篩。√N以下の数に対して2から順にその数の倍数を消していく。
    計算量は「調和級数」になるのがミソ。具体的には
        N/2 + N/3 + N/5 + ... + N/(√N) = N * (1/2 + 1/3 + 1/5 + ... + 1/√N)
                                       = N * loglog(√N) (素数の逆数和の発散スピードがこれ)
                                       = N(loglogN - log2)
    """
    if not isinstance(n, int):
        raise TypeError('n is int type.')
    if n < 2:
        return [0] * (n + 1)
    prime = [1] * (n + 1)
    prime[0] = prime[1] = 0
    for i in range(2, int(n**0.5) + 1):
        if not prime[i]: continue
        for j in range(i * 2, n + 1, i):
            prime[j] = 0
    return prime
