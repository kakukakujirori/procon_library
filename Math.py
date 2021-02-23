import math
import numpy as np
from numba import njit, jitclass, b1, i1, i4, i8, f8


spec = [
    ('n_max', i8),
    ('mod', i8),
    ('modinv', i8[:]),
    ('fac', i8[:]),
    ('facinv', i8[:]),
]

@jitclass(spec)
class Combination:
    """
    O(n)の前計算を1回行うことで，O(1)でnCr mod mを求められる
    n_max = 10**6のとき前処理は約950ms (PyPyなら約340ms, 10**7で約1800ms)
    使用例：
    comb = Combination(1000000)
    print(comb.calc(5, 3))  # 10
    """
    def __init__(self, n_max, mod=10**9+7):
        self.n_max = n_max
        self.mod = mod
        self.modinv = np.zeros(n_max+1, dtype=np.int64)
        self.fac = np.ones(n_max+1, dtype=np.int64)
        self.facinv = np.ones(n_max+1, dtype=np.int64)
        self.make_modinv_list()
        self.make_factorial_list()

    def calc(self, n, r):
        return self.fac[n] * self.facinv[r] % self.mod * self.facinv[n-r] % self.mod

    def make_factorial_list(self):
        # 階乗のリストと階乗のmod逆元のリストを返す O(n)
        # self.make_modinv_list()が先に実行されている必要がある
        for i in range(1, self.n_max+1):
            self.fac[i] = self.fac[i-1] * i % self.mod
            self.facinv[i] = self.facinv[i-1] * self.modinv[i] % self.mod

    def make_modinv_list(self):
        # 0からnまでのmod逆元のリストを返す O(n)
        self.modinv[1] = 1
        for i in range(2, self.n_max+1):
            self.modinv[i] = self.mod - self.mod//i * self.modinv[self.mod%i] % self.mod


def factorization(n):
    assert n > 1, "factorization: input must be over 2, but given {}".format(n)
    arr = []
    temp = n
    for i in range(2, int(round(n ** 0.5)) + 1):
        cnt = 0
        while temp % i == 0:
            cnt += 1
            temp //= i
        if cnt > 0: arr.append([i, cnt])
    if temp!=1: arr.append([temp, 1])
    if not arr: arr.append([n, 1])
    return arr


def get_sieve_of_eratosthenes(n):
    """
    エラトステネスの篩。ただし高速素因数分解（cf.ABC177E）にも対応できるよう
    prime[i] = (iを割り切る最小の素数)が記録される。

    方針としては√N以下の数に対して2から順にその数の倍数を消していく。
    計算量は「調和級数」になるのがミソ。具体的には
        N/2 + N/3 + N/5 + ... + N/(√N) = N * (1/2 + 1/3 + 1/5 + ... + 1/√N)
                                       = N * loglog(√N) (素数の逆数和の発散スピードがこれ)
                                       = N(loglogN - log2)
    """
    if not isinstance(n, int):
        raise TypeError('n is int type.')
    prime = [i for i in range(n + 1)]
    for p in range(2, int(n**0.5) + 1):
        if prime[p] != p: continue
        for i in range(p * 2, n + 1, p): prime[i] = p
    return prime


def euler_phi(n):
    """
    オイラーのトーシェント関数（1,...,Nの中でNと互いに素なものの個数）を計算する。
    計算量はO(√n)。計算には次の公式を用いている：
        phi(n) = n * ¥prod_{p:prime, p|n} (1 - 1/p)
    """
    assert n > 0, "euler_phi: input must be positive, but given {}".format(n)
    if n <= 2: return 1
    if n <= 4: return 2
    ret = n
    root_n = int(n ** 0.5) - 1
    while (root_n + 1) ** 2 <= n: root_n += 1
 
    for x in range(2, root_n+1):
        if n % x == 0:
            ret -= ret // x
            while n % x == 0:
                n //= x
    if n > 1:
        ret -= ret // n
    return ret


@njit("i8(i8,)", cache=True)
def isqrt(n):
    """
    Newton's method
    nの平方根のfloorを求める。math.isqrtにも実装されてる。
    """
    x, y = n, (n + 1) // 2
    while y < x:
        x, y = y, (y + n // y) // 2
    return x


@numba.njit("i8(i8,i8,i8)", cache=True)
def pow_mod(base, exp, mod):
    exp %= mod - 1
    res = 1
    while exp:
        if exp & 1:
            res = res * base % mod
        base = base * base % mod
        exp >>= 1
    return res