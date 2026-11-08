#!/usr/bin/env python3
"""
m0_regularity.py -- experiments for Problem 9.1 (is n -> m0(n) p-regular?)
Python >= 3.8, standard library only.

  python m0_regularity.py selftest
  python m0_regularity.py stats  P NMAX            # e.g. stats 5 20000
  python m0_regularity.py kernel P L MAXDIM DEPTH  # e.g. kernel 3 40 80 8 ; kernel 5 60 200 9
  python m0_regularity.py pump   P KMAX NFAM       # e.g. pump 5 14 40
"""
import sys, time, random, csv
from math import comb

MOD = (1 << 61) - 1          # prime used for linear algebra / Berlekamp-Massey

def digits(n, p, K=None):
    d = []
    while n:
        d.append(n % p); n //= p
    if K is not None:
        d += [0] * (K - len(d))
    return d

def ndig(n, p):
    k = 0
    while n:
        n //= p; k += 1
    return k

class Cascade:
    """Exact arithmetic for S_D(t) mod p, good levels, nu, and the cascade time m0."""

    def __init__(self, p):
        assert p % 2 == 1, "p must be an odd prime"
        self.p = p
        self.M = self._matrices()                       # M[b][r] = (m11, m12, m21, m22)
        self.vecs = [(x, y) for x in range(p) for y in range(p)]
        self.good_cache, self.nu_cache, self.C_cache = {}, {}, {}

    # ---- Lemma 5.1 matrices -------------------------------------------------
    def _alpha(self, j, r):
        return sum(comb(r, i) for i in range(min(j, r) + 1)) % self.p

    def _matrices(self):
        p, M = self.p, []
        for b in range(p):
            row = []
            for r in range(p):
                t = pow(2, r, p)
                if b == 0:
                    row.append((1, (t - 1) % p, 0, t))
                else:
                    a1, a0 = self._alpha(b, r), self._alpha(b - 1, r)
                    row.append((a1, (t - a1) % p, a0, (t - a0) % p))
            M.append(row)
        return M

    def mul(self, A, v):
        p = self.p
        return ((A[0] * v[0] + A[1] * v[1]) % p, (A[2] * v[0] + A[3] * v[1]) % p)

    def S(self, D, t):
        """S_D(t) mod p via the product formula (most significant digit first)."""
        if D < 0:
            return 0
        p = self.p; K = ndig(D, p)
        bD, rt = digits(D, p, K), digits(t % p ** K, p, K)
        v = (1, 0)
        for i in range(K - 1, -1, -1):
            v = self.mul(self.M[bD[i]][rt[i]], v)
        return v[0]

    # ---- good levels and nu (Theorem 5.4, Corollary 4.7) ---------------------
    def proj(self, v):
        x, y = v
        if x:
            return (1, y * pow(x, -1, self.p) % self.p)
        return (0, 1) if y else (0, 0)

    def is_good(self, q):
        r = self.good_cache.get(q)
        if r is None:
            Q = {(1, 0)}
            for b in reversed(digits(q, self.p)):
                Q = {self.proj(self.mul(A, v)) for v in Q for A in self.M[b]}
            r = self.good_cache[q] = all(v[0] != 0 for v in Q)
        return r

    def nu(self, n):
        r = self.nu_cache.get(n)
        if r is not None:
            return r
        p = self.p; q, s = divmod(n, p)
        if q == 0:
            r = 1 + max(j for j in range(s) if self.is_good(j))
        else:
            r = p * self.nu(q)
            if self.is_good(q):                       # Lemma 4.1: otherwise no interior good level
                for j in range(s - 1, -1, -1):
                    if self.is_good(p * q + j):
                        r = max(r, p * q + 1 + j); break
        self.nu_cache[n] = r
        return r

    # ---- digit-DP search on the 2-D automaton -------------------------------
    def _complete(self, A, K, acc):
        """C[j] = states from which the j lowest positions can still reach acceptance."""
        C = [frozenset(v for v in self.vecs if acc(v))]
        for j in range(1, K + 1):
            Aj, prev = A[j - 1], C[j - 1]
            C.append(frozenset(v for v in self.vecs
                               if any(self.mul(Aj[d], v) in prev for d in range(self.p))))
        return C

    def _C(self, key, A, K, acc):
        C = self.C_cache.get(key)
        if C is None:
            if len(self.C_cache) > 300000:
                self.C_cache.clear()
            C = self.C_cache[key] = self._complete(A, K, acc)
        return C

    def _search(self, A, C, K, bound, acc, mode):
        """mode 'min': least x >= bound accepted; 'max': largest x <= bound. 0 <= x < p^K."""
        p = self.p; bd = digits(bound, p, K); up = (mode == 'min')

        def free(i, v):
            x = 0
            for pos in range(i, -1, -1):
                for d in (range(p) if up else range(p - 1, -1, -1)):
                    w = self.mul(A[pos][d], v)
                    if w in C[pos]:
                        v = w; x += d * p ** pos; break
            return x

        def dfs(i, v):
            if i < 0:
                return 0 if acc(v) else None
            for d in (range(bd[i], p) if up else range(bd[i], -1, -1)):
                w = self.mul(A[i][d], v)
                if w not in C[i]:
                    continue
                sub = dfs(i - 1, w) if d == bd[i] else free(i - 1, w)
                if sub is not None:
                    return d * p ** i + sub
            return None

        return dfs(K - 1, (1, 0))

    def next_zero(self, D, m, K):
        """least t > m with S_D(t) = 0 (S_D has period p^K since D < p^K)."""
        p = self.p; P = p ** K
        bD = digits(D, p, K)
        A = [self.M[bD[i]] for i in range(K)]          # A[i][r]
        acc = lambda v: v[0] == 0
        C = self._C(('Z', D, K), A, K, acc)
        if (1, 0) not in C[K]:
            return None                                 # D is good
        x = m % P; base = m - x
        if x + 1 < P:
            z = self._search(A, C, K, x + 1, acc, 'min')
            if z is not None:
                return base + z
        return base + P + self._search(A, C, K, 0, acc, 'min')

    def top_nonzero(self, hi, t, K):
        """largest D' <= hi with S_{D'}(t) != 0."""
        p = self.p; tr = t % p ** K
        rt = digits(tr, p, K)
        A = [[self.M[b][rt[i]] for b in range(p)] for i in range(K)]   # A[i][b]
        acc = lambda v: v[0] != 0
        C = self._C(('N', tr, K), A, K, acc)
        return self._search(A, C, K, hi, acc, 'max')

    def m0(self, n, cap=10 ** 7):
        """(m0(n), number of drops delta(n)); m0 = 0 for periodic n."""
        F = self.nu(n) - 1; D = n - 1
        if D <= F:
            return 0, 0
        K = ndig(n - 1, self.p); m = 0; drops = 0
        while D > F:
            t = self.next_zero(D, m, K)
            D = self.top_nonzero(D - 1, t, K)
            m = t; drops += 1
            if drops > cap:
                raise RuntimeError('drop cap exceeded at n=%d' % n)
        return m, drops

# ---- independent brute force ---------------------------------------------------
def m0_brute(p, n, F):
    D, row, m = n - 1, [1] + [0] * (n - 1), 0
    while D > F:
        m += 1
        row = [1] + [(row[j] + row[j - 1]) % p for j in range(1, n)]
        s, pre = 0, []
        for j in range(D + 1):
            s = (s + row[j]) % p; pre.append(s)
        while pre[D] == 0:
            D -= 1
    return m

def orbit_mu(n, p):
    f, seen, t = (1,) * n, {}, 0
    while f not in seen:
        seen[f] = t
        if f[0] == 0:
            f = f[1:]
        else:
            g = [0] * (len(f) + 1)
            for i, c in enumerate(f):
                g[i] = (g[i] + c) % p; g[i + 1] = (g[i + 1] + c) % p
            g[0] = (g[0] - f[0]) % p
            while g[-1] == 0:
                g.pop()
            f = tuple(g)
        t += 1
    return seen[f]

# ---- Theorems R2, R3 -----------------------------------------------------------
def core(p, N):
    ds = digits(N, p); s = 0
    for i in range(len(ds) - 1, -1, -1):
        if ds[i] not in (0, p - 1):
            s = i + 1; break
    Q = [i for i in range(s, len(ds)) if ds[i] == p - 1]
    return Q, s, N % p ** s

def ins(p, x, Q):
    for q in sorted(Q):
        low = x % p ** q
        x = (x - low) * p + (p - 1) * p ** q + low
    return x

def astar(p):
    for a in range(1, p):
        for r in range(p - 1):
            if (pow(2, r, p) + a) % p == 0:
                return a, r

def lead1_bound(p, B):
    u = ndig(B, p) - 1; d = B - p ** u + 1
    w = 0
    while p ** w < d:
        w += 1
    ub = (p - 1) * p ** u
    for a in range(1, p):
        for b in range(p - 1):
            if (pow(2, b, p) + a) % p == 0 and b * p ** w < p ** u:
                ub = min(ub, a * p ** u + b * p ** w)
    return ub

# ---- linear algebra mod MOD ----------------------------------------------------
class Lin:
    def __init__(self):
        self.rows = []                                   # (pivot, vector, coefficients)

    def reduce(self, w):
        w = [x % MOD for x in w]; cmb = {}
        for piv, rv, rc in self.rows:
            f = w[piv]
            if f:
                w = [(x - f * y) % MOD for x, y in zip(w, rv)]
                for k, c in rc.items():
                    cmb[k] = (cmb.get(k, 0) + f * c) % MOD
        return w, cmb                                    # original = w + sum cmb_k basis_k

    def add(self, w, cmb, idx):
        piv = next(i for i, x in enumerate(w) if x)
        inv = pow(w[piv], MOD - 2, MOD)
        rv = [x * inv % MOD for x in w]
        rc = {k: (-c) * inv % MOD for k, c in cmb.items()}
        rc[idx] = (rc.get(idx, 0) + inv) % MOD
        self.rows.append((piv, rv, rc))

def bm(seq):
    """Berlekamp-Massey mod MOD: length of the minimal linear recurrence."""
    s = [x % MOD for x in seq]
    C, B, L, m, b = [1], [1], 0, 1, 1
    for n in range(len(s)):
        d = s[n]
        for i in range(1, L + 1):
            d = (d + C[i] * s[n - i]) % MOD
        if d == 0:
            m += 1; continue
        T = C[:]; coef = d * pow(b, MOD - 2, MOD) % MOD
        if len(C) < len(B) + m:
            C += [0] * (len(B) + m - len(C))
        for i in range(len(B)):
            C[i + m] = (C[i + m] - coef * B[i]) % MOD
        if 2 * L <= n:
            L, B, b, m = n + 1 - L, T, d, 1
        else:
            m += 1
    return L

# ---- commands ------------------------------------------------------------------
def selftest():
    t0 = time.time()
    for p, nmax in ((3, 120), (5, 120), (7, 100)):
        cas = Cascade(p)
        for n in range(1, nmax + 1):
            v = cas.nu(n)
            assert cas.m0(n)[0] == (m0_brute(p, n, v - 1) if v < n else 0), (p, n)
    anchors = [(3, 5, 6), (5, 6, 7), (7, 8, 23)]          # checked by hand
    for p, n, val in anchors:
        assert Cascade(p).m0(n)[0] == val, (p, n)
    cas = Cascade(5)
    for n in range(2, 26):
        v = cas.nu(n)
        if v < n:
            assert orbit_mu(n, 5) == 2 * cas.m0(n)[0] + n - v - 1, n
    cas = Cascade(3)
    for n in range(2, 3001):
        v = cas.nu(n); d = n - v
        if d:
            w = 0
            while 3 ** w < d:
                w += 1
            assert cas.m0(n)[0] == v + 3 ** w, n
    for p in (5, 7, 11, 13):
        cas = Cascade(p); a, r = astar(p)
        for k in range(1, 11):
            assert cas.m0(p ** k + 1)[0] == a * p ** k + r, (p, k)
        for n in range(1, 120):
            assert cas.m0(p * n)[0] == p * cas.m0(n)[0], (p, n)
    cas = Cascade(5)
    for n in range(2, 2001):
        v = cas.nu(n)
        if v == n:
            continue
        Q, s, B = core(5, n - 1)
        assert ins(5, cas.m0(B + 1)[0], Q) == cas.m0(n)[0], n
    print('ALL SELFTESTS PASSED  (%.1fs)' % (time.time() - t0))

def stats(p, nmax):
    cas = Cascade(p); t0 = time.time()
    res = {n: cas.m0(n) for n in range(1, nmax + 1)}
    best, byd, dropmax = {}, {}, {}
    lead1 = lead1_eq = 0
    with open('m0_stats_p%d.csv' % p, 'w', newline='') as fh:
        wr = csv.writer(fh)
        wr.writerow(['n', 'n-1 base p', 'nu', 'd', 'm0', 'mu', 'drops', 'core B', 'lead1 UB'])
        for n in range(1, nmax + 1):
            m, dr = res[n]; v = cas.nu(n); d = n - v; k = ndig(n - 1, p)
            dropmax[k] = max(dropmax.get(k, 0), dr)
            if d == 0:
                wr.writerow([n, digits(n - 1, p)[::-1], v, 0, 0, 0, 0, '', '']); continue
            mu = 2 * m + d - 1
            assert m >= n
            if mu / n > best.get(k, (0, 0))[0]:
                best[k] = (mu / n, n)
            byd.setdefault(d, set()).add(m - v)
            Q, s, B = core(p, n - 1)
            mm = res[B + 1][0]
            assert ins(p, mm, Q) == m, n                  # Theorem R2
            ub = ''
            if digits(B, p)[-1] == 1:                     # Theorem R3
                ub = lead1_bound(p, B)
                assert n <= m and mm <= ub, n
                lead1 += 1; lead1_eq += (mm == ub)
            wr.writerow([n, digits(n - 1, p)[::-1], v, d, m, mu, dr, B, ub])
    print('p=%d n<=%d  time %.1fs' % (p, nmax, time.time() - t0))
    print('  R2 (exact peeling) verified for all nonperiodic n')
    print('  R3 leading-1 cores: %d cases, UB attained in %d' % (lead1, lead1_eq))
    print('  d with >1 values of m0-nu: %d of %d' % (sum(len(x) > 1 for x in byd.values()), len(byd)))
    print('  max mu/n by digit length k:', {k: (round(r, 4), n) for k, (r, n) in sorted(best.items())})
    print('  max drops delta(n) by k   :', dict(sorted(dropmax.items())))

def kernel(p, L, maxdim, depth_cap, nverify=300):
    cas = Cascade(p); cache = {}

    def a(n):
        if n == 0:
            return 0
        v = cache.get(n)
        if v is None:
            v = cache[n] = cas.m0(n)[0]
        return v

    def vec(e, j):
        return [a(p ** e * n + j) for n in range(L)]

    t0 = time.time()
    lin, basis, trans = Lin(), [(0, 0)], {}
    w, cm = lin.reduce(vec(0, 0)); lin.add(w, cm, 0)
    queue, closed = [0], True
    while queue:
        i = queue.pop(0); e, j = basis[i]
        for c in range(p):
            ce, cj = e + 1, j + c * p ** e
            w, cm = lin.reduce(vec(ce, cj))
            if any(w):
                if len(basis) >= maxdim or ce > depth_cap:
                    closed = False; queue = []; break
                idx = len(basis); lin.add(w, cm, idx); basis.append((ce, cj)); queue.append(idx)
                trans[(i, c)] = {idx: 1}
            else:
                trans[(i, c)] = cm
    prof = {}
    for e, _ in basis:
        prof[e] = prof.get(e, 0) + 1
    print('p=%d L=%d  basis size %d  closed=%s  per-depth %s  (%.1fs)'
          % (p, L, len(basis), closed, dict(sorted(prof.items())), time.time() - t0))
    if not closed:
        print('  -> no finite representation found within limits (evidence, not proof)')
        return
    s0 = [a(j) for (_, j) in basis]
    R = len(basis)

    def rep(n):
        v = s0[:]
        for c in reversed(digits(n, p)):
            v = [sum(cf * v[k] for k, cf in trans[(i, c)].items()) % MOD for i in range(R)]
        return v[0]

    rnd = random.Random(1); bad = 0
    hi = p ** (max(e for e, _ in basis) + 4) * L
    for _ in range(nverify):
        n = rnd.randrange(1, hi)
        if rep(n) != a(n) % MOD:
            bad += 1
    print('  verification on %d random n < %d : %d mismatches' % (nverify, hi, bad))

def pump(p, kmax, nfam, seed=2):
    cas = Cascade(p); rnd = random.Random(seed)
    fams = [([1], [0], []), ([1], [0], [1]), ([1], [1], []), ([2], [0], [1])]
    while len(fams) < nfam:
        x = [rnd.randrange(1, p)] + [rnd.randrange(p) for _ in range(rnd.randrange(2))]
        y = [rnd.randrange(p) for _ in range(rnd.randrange(1, 3))]
        z = [rnd.randrange(p) for _ in range(rnd.randrange(3))]
        fams.append((x, y, z))
    for x, y, z in fams:
        seq, t0 = [], time.time()
        for k in range(kmax + 1):
            N = 0
            for dgt in x + y * k + z:
                N = N * p + dgt
            seq.append(cas.m0(N + 1)[0])
        L = bm(seq)
        flag = 'ok' if 2 * L + 2 <= len(seq) else 'UNDETERMINED/NO RECURRENCE'
        print('x=%s y=%s z=%s  BM order %d  [%s]  (%.1fs)' % (x, y, z, L, flag, time.time() - t0))

if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'selftest'
    args = [int(a) for a in sys.argv[2:]]
    {'selftest': selftest, 'stats': stats, 'kernel': kernel, 'pump': pump}[cmd](*args)
