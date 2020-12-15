class PolyLib():
    """
    線形漸化式を求めるためのライブラリ。と言っても実態は多項式操作の関数群。
    漸化式がどの体上で定義されているかによって加算と乗算は適切に定義する必要がある。
    クラスにまとめるのにラムダ式を使ったせいでだいぶ遅くなってることに注意。
    """
    def __init__(self, mod=None):
        if mod is None:
            self.add = lambda x, y: x + y
            self.sub = lambda x, y: x - y
            self.mul = lambda x, y: x * y
        else:
            self.add = lambda x, y: (x + y) % mod
            self.sub = lambda x, y: (x - y) % mod
            self.mul = lambda x, y: (x * y) % mod

    def diminish_zero(self, P):
        """多項式Pの高次にある余分な0を削る"""
        while P and P[-1] == 0: P.pop()
        if not P: P.append(0)
        return P

    def polymul(self, P, Q, d=None):
        """多項式PとQの積をd次まで計算する。dを指定しなければP*Qは最高次まで求める"""
        assert P and Q, "Inputs must not be empty, but given P = {}, Q = {}".format(P, Q)
        if d is None: d = len(P) + len(Q) - 2  # 桁が膨れ上がる場合はここをtruncateする
        ret = [0] * (d + 1)
        for n in range(d + 1):
            coeff = 0
            for i in range(n + 1):
                if i >= len(P) or n - i >= len(Q): continue
                coeff = self.add(coeff, self.mul(P[i], Q[n - i]))
            ret[n] = coeff
        return self.diminish_zero(ret)

    def reduce_even(self, P):
        """偶関数P(x) = P'(x^2)に対してP'(x)を返す"""
        ret = []
        for i, p in enumerate(P):
            if i % 2 == 0:
                ret.append(p)
            else:
                assert p == 0, "P must be even, but has an odd term P[{}] = {}".format(i, p)
        return self.diminish_zero(ret)

    def divide_even_odd(self, P, reduce=True):
        """
        P(x)を奇関数と偶関数に分ける。
        reduce=TrueならP(x) = E(x^2) + xO(x^2)としてE(x), O(x)を返す。
        そうでなければE(x^2)とxO(x^2)を返す。
        """
        Even = []
        Odd = []
        for i, p in enumerate(P):
            if i % 2 == 0:
                Even.append(p)
                if not reduce: Even.append(0)
            else:
                if not reduce: Odd.append(0)
                Odd.append(p)
        # alleviate empty polynomials
        if not Odd: Odd.append(0)
        return self.diminish_zero(Even), self.diminish_zero(Odd)


    def get_coeff(self, P, Q, n):
        """形式的冪級数P(x)/Q(x)のx^nにおける係数を求める"""
        assert Q[0] == 1, "The constant term of denominator must be 1, but given, {}".format(Q[0])
        assert n >= 0, "n must be non-negative, but given {}".format(n)
        if n == 0: return P[0]

        Q_trans = [0] * len(Q)
        for i, q in enumerate(Q):
            if i % 2 == 0:
                Q_trans[i] = q
            else:
                Q_trans[i] = -q
        V = self.reduce_even(self.polymul(Q, Q_trans))
        Ue, Uo = self.divide_even_odd(self.polymul(P, Q_trans))
        if n % 2 == 0:
            return self.get_coeff(Ue, V, n // 2)
        else:
            return self.get_coeff(Uo, V, (n - 1) // 2)

    def solve_linear_recurrence_naive(self, A, C, N=10):
        """
        (d + 1)次線形漸化式
            a_n = c_1 * a_{n-1} + c_2 * a_{n-2} + ... + c_k * a_{n-d}  (n >= d)
        の第N項までを愚直に求める。

        In:
            A = [a_0, a_1, a_2, ..., a_{d-1}]
            C = [c_1, c_2, ..., c_d]
        Out:
            [a_0, a_1, a_2, ..., a_{d-1}, a_d, a_{d+1}, ..., a_N]
        """
        assert len(A) == len(C)
        for _ in range(len(A), N + 1):
            val = 0
            for i, c in enumerate(C):
                val = self.add(val, self.mul(c, A[-i-1]))
            A.append(val)
        return A

    def solve_linear_recurrence(self, A, C, N):
        """
        (d + 1)次線形漸化式
            a_n = c_1 * a_{n-1} + c_2 * a_{n-2} + ... + c_k * a_{n-d}  (n >= d)
        の第N項(0-indexed)をO(d^2logN)で求める。
        参考：http://q.c.titech.ac.jp/docs/progs/polynomial_division.html
        多項式の積の計算にFFTを使えばO(klogklogN)でいけるみたいだが、
        FFTにまだ手を出せる実力じゃないので将来の課題とする。

        In:
            A = [a_0, a_1, a_2, ..., a_{d-1}]
            C = [c_1, c_2, ..., c_d]
        Out:
            [a_0, a_1, a_2, ..., a_{d-1}, a_d, a_{d+1}, ..., a_N]
        """
        assert len(A) == len(C)
        Q = [1] + [-c for c in C]
        P = self.polymul(A, Q, len(A) - 1)

        return self.get_coeff(P, Q, N)