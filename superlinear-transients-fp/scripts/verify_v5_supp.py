#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""verify_v5_supp.py  v5 补充核验；输出 verify_v5_supp_log.txt"""
import sys, time, random, traceback
from math import comb
from bisect import bisect_right
from functools import lru_cache
from fractions import Fraction
import numpy as np
sys.setrecursionlimit(10000)
T0 = time.time()
LOGF = open('verify_v5_supp_log.txt', 'w', encoding='utf-8')
def log(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); LOGF.write(s + '\n'); LOGF.flush()
SUM = []
def check(name, ok, detail=''):
    tag = {True: 'PASS', False: 'FAIL', None: 'INFO'}[ok]; SUM.append((tag, name, detail))
    log('  [%s] %s  %s' % (tag, name, detail))

_M = {}
def mats(p):
    if p not in _M:
        def th(j, r): return 0 if j < 0 else sum(comb(r, i) for i in range(min(j, r) + 1)) % p
        M = {}
        for b in range(p):
            for r in range(p):
                two = pow(2, r, p); x = th(b, r); y = th(b - 1, r)
                M[b, r] = (x, (two - x) % p, y, (two - y) % p)
        _M[p] = M
    return _M[p]
def digs(x, p):
    d = []
    while x: d.append(x % p); x //= p
    return d[::-1] or [0]
def plen(x, p): return len(digs(x, p))
def binom_mod(n, k, p):
    if k < 0 or k > n: return 0
    r = 1
    while n or k:
        x, y = n % p, k % p
        if y > x: return 0
        r = r * comb(x, y) % p; n //= p; k //= p
    return r
def Sdig(D, m, p):
    if D < 0: return 0
    M = mats(p); bs = []; rs = []
    while D: bs.append(D % p); rs.append(m % p); D //= p; m //= p
    u, v = 1, 0
    for i in range(len(bs) - 1, -1, -1):
        a1, a2, a3, a4 = M[bs[i], rs[i]]
        u, v = (a1 * u + a2 * v) % p, (a3 * u + a4 * v) % p
    return u
START = frozenset({(1, 0)})
def stepR(R, b, p):
    M = mats(p); out = set()
    for r in range(p):
        a1, a2, a3, a4 = M[b, r]
        for (u, v) in R: out.add(((a1 * u + a2 * v) % p, (a3 * u + a4 * v) % p))
    return frozenset(out)
def okR(R): return all(u for u, v in R)
@lru_cache(maxsize=None)
def Rstate(D, p):
    if D == 0: return START
    return stepR(Rstate(D // p, p), D % p, p)
def good(D, p): return D >= 0 and okR(Rstate(D, p))
def isrigid(p, W):
    return W >= 1 and good(W, p) and all(not good(p * W + c, p) for c in range(p - 1))
STARTN = frozenset({0})
def stepN(S, b, p):
    M = mats(p); out = set()
    for r in range(p):
        a1, a2, a3, a4 = M[b, r]
        for w in S:
            U = (a1 + a2 * w) % p
            if U == 0: return None
            out.add((a3 + a4 * w) % p * pow(U, p - 2, p) % p)
    return frozenset(out)
def goodN(D, p):
    S = STARTN
    for b in digs(D, p):
        S = stepN(S, b, p)
        if S is None: return False
    return True
def m0_acc(n, p, limit=5 * 10 ** 6):
    D = n - 1
    if D < 0 or good(D, p): return 0
    m = 0; it = 0
    while True:
        x = D + 1; j = 0
        while x % p == 0: x //= p; j += 1
        st = p ** j; m = (m // st + 1) * st
        if Sdig(D, m, p) == 0:
            D -= 1
            while Sdig(D, m, p) == 0: D -= 1
            if good(D, p): return m
        it += 1
        if it > limit: return None

class Fam:
    def __init__(s, p, W):
        s.p = p; s.W = W; s.ell = plen(W, p); s.Mod = p ** s.ell
        inv = lambda x: pow(x, p - 2, p)
        s.TW = [t for t in range(s.Mod) if binom_mod(t, W, p)]
        lab = {t: (-Sdig(W - 1, t, p)) * inv(binom_mod(t, W, p)) % p for t in s.TW}
        sig = [[Sdig(c, r, p) * inv(pow(2, r, p)) % p for r in range(p)] for c in range(p)]
        s.Z = [sorted((t, r) for t in s.TW for r in range(p) if sig[c][r] == lab[t]) for c in range(p)]
        s.Zs = [set(z) for z in s.Z]; s.bc = {}
    def nxt(s, c, E, t, r):
        z = s.Z[c]; i = bisect_right(z, (t, r))
        return (E,) + z[i] if i < len(z) else (E + 1,) + z[0]
    def Ph(s, c, E, t, r):
        p = s.p
        while True:
            if not s.Z[c]: return ('STUCK', c, E, t, r)
            E, t, r = s.nxt(c, E, t, r)
            cs = [c2 for c2 in range(c) if (t, r) not in s.Zs[c2]]
            if cs: c = max(cs); continue
            if (t, p - 1) in s.Zs[p - 2]: return ('DEEP', E, t, r)
            return ('OK', E, t, r)
    def bterm(s, th):
        if th not in s.bc:
            p = s.p; D = s.W - 1; m = th
            while not good(D, p):
                m += 1
                while Sdig(D, m, p) == 0: D -= 1
            s.bc[th] = m
        return s.bc[th]
    def tm(s, m, E, t, r):
        p = s.p; return E * p ** (m + s.ell) + t * p ** m + p ** m - p + r
    def m0_formula(s, d, k):
        p = s.p
        if not s.Z[d]: return 0, 'gs'
        res = s.Ph(d, 0, 0, -1)
        if res[0] == 'DEEP': return None, 'deep'
        if res[0] == 'STUCK': return s.tm(k, *res[2:]), 'st'
        _, E, th, r = res
        if good(s.W - 1, p): return s.tm(k, E, th, r), 'w1'
        for i in range(1, k):
            res = s.Ph(p - 2, E, th, p - 1)
            if res[0] == 'DEEP': return None, 'deep'
            if res[0] == 'STUCK': return p ** i * s.tm(k - i, *res[2:]), 'ss'
            _, E, th, _r = res
        return p ** k * (s.Mod * E + s.bterm(th)), 'adm'
    def classify(s, d):
        p = s.p
        if not s.Z[d]: return 'good-start'
        res = s.Ph(d, 0, 0, -1)
        if res[0] == 'DEEP': return 'deep'
        if res[0] == 'STUCK': return 'linear'
        if good(s.W - 1, p): return 'linear-W1'
        th = res[2]; seen = set()
        while th not in seen:
            seen.add(th); r2 = s.Ph(p - 2, 0, th, p - 1)
            if r2[0] == 'DEEP': return 'deep'
            if r2[0] == 'STUCK': return 'linear'
            th = r2[2]
        return 'admissible'
    def growth(s, d):
        p = s.p; th = s.Ph(d, 0, 0, -1)[2]; idx = {}; seq = []
        while th not in idx:
            idx[th] = len(seq); r2 = s.Ph(p - 2, 0, th, p - 1)
            seq.append((th, r2[1])); th = r2[2]
        cyc = seq[idx[th]:]; rho = Fraction(sum(x for _, x in cyc), len(cyc))
        return rho, len(cyc), 2 * rho * s.Mod / (s.W + 1)

# ---------- S0 前缀封闭恒等式 + 射影自动机闭合/偶数位
def S0():
    rnd = random.Random(1); bad = 0
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41):
        for _ in range(300):
            D = rnd.randrange(p ** 4); t = rnd.randrange(p ** 5); b = rnd.randrange(p)
            if Sdig(p * D + b, p * t, p) != Sdig(D, t, p): bad += 1
    check('S0a S_{pD+b}(pt)=S_D(t) 抽查 (p<=41, 3600 例)', bad == 0, 'viol=%d' % bad)
    rnd = random.Random(2)
    for p in (13, 17, 19, 23):
        bad = 0
        for _ in range(200):
            D = rnd.randrange(p ** 3)
            if good(D, p) != goodN(D, p): bad += 1
        check('S0b p=%d 定义(未剪枝) = 射影自动机判定, 200 个随机 D<p^3' % p, bad == 0, 'viol=%d' % bad)
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53):
        t1 = time.time(); seen = {STARTN}; q = [STARTN]; i = 0; odd = 0; dead = 0
        while i < len(q) and time.time() - t1 < 120:
            S = q[i]; i += 1; allowed = 0
            for b in range(p):
                N = stepN(S, b, p)
                if N is None: continue
                allowed += 1
                if b % 2 == 1: odd += 1
                if N not in seen: seen.add(N); q.append(N)
            if allowed == 0: dead += 1
        comp = i >= len(q)
        check('S0c p=%d 射影自动机闭合且无奇数位转移' % p, (odd == 0) if comp else None,
              '状态 %d%s, 奇数位转移 %d, 无出边状态 %d' % (len(seen), '' if comp else ' [未闭合]', odd, dead))

# ---------- S1 p=3 闭式 m0 = q*3^i + 3^{len(y)}
def S1():
    p = 3; L = 8; P = p ** L
    T = np.empty((P, P), dtype=np.uint8); row = np.zeros(P, dtype=np.int64); row[0] = 1
    for m in range(P):
        T[m] = (np.cumsum(row) % p).astype(np.uint8)
        nr = row.copy(); nr[1:] = (row[1:] + row[:-1]) % p; row = nr
    TT = np.ascontiguousarray(T.T); del T
    g = (TT != 0).all(axis=1); Z = {}
    def m0t(n):
        D = n - 1
        if g[D]: return 0
        m = 0
        while True:
            if D not in Z: Z[D] = np.flatnonzero(TT[D] == 0).tolist()
            z = Z[D]; k, r = divmod(m, P); i = bisect_right(z, r)
            m = k * P + z[i] if i < len(z) else (k + 1) * P + z[0]
            while TT[D, m % P] == 0: D -= 1
            if g[D]: return m
    def pred(n):
        N = n - 1; ds = digs(N, 3)
        if 1 not in ds: return None
        i = len(ds) - 1 - ds.index(1); q = N // 3 ** i; y = N % 3 ** i
        return q * 3 ** i + (1 if y == 0 else 3 ** plen(y, 3))
    bad = []; badE = 0; cnt = 0
    for n in range(2, P + 1):
        if good(n - 1, 3): continue
        m = m0t(n); cnt += 1
        if m != pred(n): bad.append((n, m, pred(n)))
        if m >= 3 ** plen(n - 1, 3): badE += 1
    check('S1a p=3 闭式 m0=q3^i+3^len(y), 全部非好 n<=%d' % P, not bad, '%d 例, 不符 %s' % (cnt, bad[:5]))
    check('S1b p=3 E(n)=0 即 m0<3^L', badE == 0, 'viol=%d' % badE)
    rnd = random.Random(3); bad = []; c = 0
    while c < 60:
        n = rnd.randrange(P, 3 ** 10)
        if good(n - 1, 3): continue
        c += 1; m = m0_acc(n, 3)
        if m != pred(n): bad.append((n, m, pred(n)))
    check('S1c p=3 闭式, 60 个随机 n∈[3^8,3^10]', not bad, '不符 %s' % bad[:5])

# ---------- S2 p=13..41 刚性族: 二分 + rho>0 + 公式对照
def S2():
    MINRIG = {13: 2, 17: 2, 19: 2, 23: 158, 29: 6, 31: 6, 37: 10, 41: 12}
    for p, W0 in MINRIG.items():
        t1 = time.time()
        Ws = [W for W in range(1, p) if isrigid(p, W)]
        if W0 not in Ws: Ws.append(W0)
        for W in Ws:
            F = Fam(p, W); cls = [F.classify(d) for d in range(p - 1)]
            alladm = all(c == 'admissible' for c in cls)
            check('S2a p=%d W=%d %s 刚性 ⇒ 全部 d 可容许' % (p, W, digs(W, p)), alladm,
                  '分类 %s' % {c: cls.count(c) for c in set(cls)})
            if not alladm: continue
            gs = [F.growth(d) for d in range(p - 1)]
            rmin = min(x[0] for x in gs); cmax = max(x[2] for x in gs)
            check('S2b p=%d W=%d ρ>0' % (p, W), rmin > 0,
                  'ρ_min=%s ρ_max=%s 循环长 %s c_max=%s' % (rmin, max(x[0] for x in gs),
                  sorted(set(x[1] for x in gs)), cmax))
            ks = (1, 2) if F.ell == 1 else (1,)
            for k in ks:
                for d in ((0, p - 2) if k == 1 else (0,)):
                    n = (W + 1) * p ** k - p + d + 1
                    v, ty = F.m0_formula(d, k); b = m0_acc(n, p)
                    check('S2c p=%d W=%d d=%d k=%d (n=%d) 公式=直接' % (p, W, d, k, n), v == b,
                          'direct=%s formula=%s' % (b, v))
        log('    p=%d 用时 %.0fs' % (p, time.time() - t1))

# ---------- S3 p=7 E(n) 抽样 (n∈[7^4,7^5])
def S3():
    rnd = random.Random(4); mE = (0, 0); mr = (0.0, 0); c = 0; t1 = time.time()
    while c < 80 and time.time() - t1 < 600:
        n = rnd.randrange(7 ** 4, 7 ** 5)
        if good(n - 1, 7): continue
        c += 1; m = m0_acc(n, 7); L = plen(n - 1, 7); E = m // 7 ** L
        if E > mE[0]: mE = (E, n)
        if m / (n - 1) > mr[0]: mr = (m / (n - 1), n)
    check('S3 p=7 抽样 %d 个 n∈[7^4,7^5]: E<=2 ?' % c, mE[0] <= 2,
          'max E=%d at n=%d; max m0/(n-1)=%.4f at n=%d' % (mE[0], mE[1], mr[0], mr[1]))

for name, fn in (('S0', S0), ('S1', S1), ('S2', S2), ('S3', S3)):
    log(''); log('=' * 60); log('[%s]  t=%.0fs' % (name, time.time() - T0))
    try: fn()
    except Exception:
        log(traceback.format_exc()); check(name + ' 运行', False, '异常')
log(''); log('SUMMARY  %.0fs' % (time.time() - T0))
for tag in ('FAIL', 'PASS', 'INFO'):
    it = [x for x in SUM if x[0] == tag]; log('%s: %d' % (tag, len(it)))
    if tag != 'PASS':
        for _, a, b in it: log('   - %s  %s' % (a, b))
LOGF.close()
