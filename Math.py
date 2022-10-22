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
    for p in range(2, int(round(n ** 0.5)) + 1):
        cnt = 0
        while temp % p == 0:
            cnt += 1
            temp //= p
        if cnt > 0: arr.append([p, cnt])
    if temp != 1: arr.append([temp, 1])
    if not arr: arr.append([n, 1])
    return arr

@njit
def get_sieve_of_eratosthenes(n):
    """
    エラトステネスの篩。ただし高速素因数分解（cf.ABC177E）にも対応できるよう
    prime[i] = (iを割り切る素数たち)が格納される。

    方針としては√N以下の数に対して2から順にその数の倍数を消していく。
    計算量は「調和級数」になるのがミソ。具体的には
        N/2 + N/3 + N/5 + ... + N/(√N) = N * (1/2 + 1/3 + 1/5 + ... + 1/√N)
                                       = N * loglog(√N) (素数の逆数和の発散スピードがこれ)
                                       = N(loglogN - log2)
    """
    primes = [[i for i in range(0)] for _ in range(n + 1)]
    for p in range(2, n + 1):
        if primes[p]: continue
        for i in range(p, n + 1, p):
            primes[i].append(p)
    return primes


@njit
def gcd_pair_num(L, R):
    """
    gcd_pair_num[k] = (L<=x,y<=Rでgcd(x,y)=kなる組(x,y)の個数)
    xとyの大小関係は不問なことに注意。
    ABC206Eがあまりに良問なのでライブラリにする次第。
        gcd(x,y)=k <=> k|xかつk|yかつ2<=sに対してNOT(sk|x && sk|y)
    なのでkが大きい方から計算するのがミソ。
    計算量はエラトステネスの篩と同じくO(NloglogN)
    """
    assert 1 <= L <= R
    gcd_pair_num = [0] * (R + 1)
    for k in range(R, 0, -1):
        num = R // k - (L - 1) // k
        gcd_pair_num[k] = num * num
        for s in range(2, R + 1):
            if k * s > R: break
            gcd_pair_num[k] -= gcd_pair_num[k * s]
    return gcd_pair_num


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
    f(x)=x^2-nとして、(xk,f(xk))における接線はy = 2xk(x - xk) + f(xk)
    よってx(k+1) = xk - (xk^2 - n) / (2xk) = (xk + n / xk) / 2 
    """
    x, y = n, (n + 1) // 2
    while y < x:
        x, y = y, (y + n // y) // 2
    return x


@njit("i8(i8,i8,i8)", cache=True)
def pow_mod(base, exp, mod):
    exp %= mod - 1
    res = 1
    while exp:
        if exp & 1:
            res = res * base % mod
        base = base * base % mod
        exp >>= 1
    return res


@njit("Tuple((i8,i8,i8))(i8,i8)", cache=True)
def extgcd(a, b):
    """
    ax+by=gcd(a,b)なるx,yを出力する。
    gcdは常に正の値を返す。x,yの範囲はよく分からん。
    """
    assert a * b != 0
    abs_a = a if a > 0 else -a
    abs_b = b if b > 0 else -b
    if abs_a < abs_b:
        exchange = True
        abs_a, abs_b = abs_b, abs_a
    else:
        exchange = False
    
    c, d = abs_a, abs_b
    x, y, u, v = 1, 0, 0, 1
    while d != 0:
        k = c // d
        x -= k * u
        y -= k * v
        x, u = u, x
        y, v = v, y
        c, d = d, c % d
    _gcd = abs_a * x + abs_b * y
    
    if _gcd < 0: x, y = -x, -y
    if exchange: x, y = y, x
    if a < 0: x *= -1
    if b < 0: y *= -1
    return x, y, _gcd


@njit("Tuple((i8,i8))(i8,i8,i8,i8)", cache=True)
def CRT(r1, m1, r2, m2):
    """
    t=r1(mod m1), t=r2(mod m2)となる最小のt>=0及びl=LCM(m1, m2)の組(t,l)を返す。
    tの最小性を除けば、t+kl(k¥in Z)が解になることに注意。
    extgcdの読み込みが必要。
    解なしの時は(0, -1)を返す。

    内容としてはt = m1*x+r1 = m2*y+r2 <=> m1*x-m2*y = r2-r1
    なので、(r2-r1)がgcd(m1,m2)の倍数でなければ解なし。
    倍数であればx0,y0,d = extgcd(m1,-m2)を用いてx,yが求まる。
    """
    x, y, d = extgcd(m1, -m2)
    if (r2 - r1) % d != 0:
        return 0, -1

    _lcm = m1 * (m2 // d)
    if _lcm < 0: _lcm *= -1

    scale = (r2 - r1) // d
    #t = m1 * x * scale + r1
    #tt = m2 * y * scale + r2
    #assert t == tt
    #とするのが自然だが、これだと３つの掛け算でオーバーフローの危険があるので以下の工夫をする。
    tmp = scale * x % (m2 // d)
    t = (r1 + m1 * tmp) % _lcm
    return t % _lcm, _lcm


def matmul(A, B, mod=None):
    ah, aw = len(A), len(A[0])
    bh, bw = len(B), len(B[0])
    assert aw == bh
    ret = [[0] * bw for _ in range(ah)]
    for i in range(ah):
        for j in range(bw):
            if mod is None:
                ret[i][j] = sum([A[i][k] * B[k][j] for k in range(aw)])
            else:
                ret[i][j] = sum([A[i][k] * B[k][j] % mod for k in range(aw)])
    return ret


def matpow(M, k, mod=None):
    l = len(M)
    ret = [[0] * l for _ in range(l)]
    for i in range(l):
        ret[i][i] = 1
    
    if k < 0:
        raise NotImplementedError
    
    while k > 0:
        if k & 1:
            ret = matmul(ret, M, mod)
        M = matmul(M, M, mod)
        k >>= 1
    
    return ret