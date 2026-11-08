#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_v5.py  —  v5 稿件全部结论的核验脚本
用法:
  python verify_v5.py              标准模式
  python verify_v5.py --full       额外运行耗时检查
  python verify_v5.py --only A,B   只运行指定节
依赖: Python 3.8+, numpy
输出: 屏幕 + verify_v5_log.txt
"""
import sys, time, argparse, random, traceback, platform
from math import comb
from bisect import bisect_right
from functools import lru_cache
from fractions import Fraction
try:
    import numpy as np
except ImportError:
    sys.exit("需要 numpy:  pip install numpy")

sys.setrecursionlimit(10000)
ap = argparse.ArgumentParser()
ap.add_argument('--full', action='store_true')
ap.add_argument('--only', default='')
ARGS = ap.parse_args()
ONLY = set(x.strip().upper() for x in ARGS.only.split(',') if x.strip())

T0 = time.time()
LOGF = open('verify_v5_log.txt', 'w', encoding='utf-8')
def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True); LOGF.write(s + '\n'); LOGF.flush()
SUMMARY = []
def check(name, ok, detail=''):
    tag = {True: 'PASS', False: 'FAIL', None: 'INFO'}[ok]
    SUMMARY.append((tag, name, detail))
    log('  [%s] %s  %s' % (tag, name, detail))

PRIMES_SMALL = (3, 5, 7, 11, 13)
TABLE_L = {3: 8, 5: 5, 7: 4, 11: 3, 13: 3}

# ---------------------------------------------------------------- 基本算术
_M = {}
def mats(p):
    """(b,r) -> 2x2 矩阵, 作用在 (S_a(t), S_{a-1}(t)) 上, 由引理 2.2 (i=1) 得到"""
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
def val(ds, p):
    x = 0
    for d in ds: x = x * p + d
    return x
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
    """S_D(m) = sum_{j<=D} C(m,j) mod p, 按数字计算"""
    if D < 0: return 0
    M = mats(p); bs = []; rs = []
    while D: bs.append(D % p); rs.append(m % p); D //= p; m //= p
    u, v = 1, 0
    for i in range(len(bs) - 1, -1, -1):
        a1, a2, a3, a4 = M[bs[i], rs[i]]
        u, v = (a1 * u + a2 * v) % p, (a3 * u + a4 * v) % p
    return u

# ---------------------------------------------------------------- 好层级自动机
START = frozenset({(1, 0)})
def stepR(R, b, p):
    M = mats(p); out = set()
    for r in range(p):
        a1, a2, a3, a4 = M[b, r]
        for (u, v) in R:
            out.add(((a1 * u + a2 * v) % p, (a3 * u + a4 * v) % p))
    return frozenset(out)
def okR(R): return all(u for u, v in R)

@lru_cache(maxsize=None)
def Rstate(D, p):           # 不剪枝: 状态 = 所有时间串下 (S_D, S_{D-1}) 的集合
    if D == 0: return START
    return stepR(Rstate(D // p, p), D % p, p)
def good(D, p): return D >= 0 and okR(Rstate(D, p))

STARTN = frozenset({0})
def stepN(S, b, p):          # 射影化状态: w = S_{D-1}/S_D ; 返回 None 表示坏
    M = mats(p); out = set()
    for r in range(p):
        a1, a2, a3, a4 = M[b, r]
        for w in S:
            U = (a1 + a2 * w) % p
            if U == 0: return None
            V = (a3 + a4 * w) % p
            out.add(V * pow(U, p - 2, p) % p)
    return frozenset(out)

def isrigid(p, W):
    return W >= 1 and good(W, p) and all(not good(p * W + c, p) for c in range(p - 1))

def nu_of(n, p):
    D = n - 1
    while not good(D, p): D -= 1
    return D + 1

# ---------------------------------------------------------------- Pascal 表与级联
TAB = {}
def table(p):
    if p in TAB: return TAB[p]
    L = TABLE_L[p]; P = p ** L; t1 = time.time()
    T = np.empty((P, P), dtype=np.uint8)
    row = np.zeros(P, dtype=np.int64); row[0] = 1
    periodic = True
    for m in range(2 * P):
        cum = (np.cumsum(row) % p).astype(np.uint8)
        if m < P: T[m] = cum
        elif periodic and not np.array_equal(T[m - P], cum): periodic = False
        nr = row.copy(); nr[1:] = (row[1:] + row[:-1]) % p; row = nr
    TT = np.ascontiguousarray(T.T); del T
    g = (TT != 0).all(axis=1)
    TAB[p] = dict(P=P, L=L, TT=TT, good=g, Z=[None] * P, periodic=periodic)
    log('    [表] p=%d, p^L=%d 已构建 (%.1fs)' % (p, P, time.time() - t1))
    return TAB[p]

def zeros_of(tb, D):
    z = tb['Z'][D]
    if z is None:
        z = np.flatnonzero(tb['TT'][D] == 0).tolist(); tb['Z'][D] = z
    return z

def m0_table(n, p):
    tb = table(p); P = tb['P']; TT = tb['TT']; g = tb['good']
    D = n - 1
    if D < 0 or g[D]: return 0
    m = 0
    while True:
        z = zeros_of(tb, D)
        k, r = divmod(m, P)
        i = bisect_right(z, r)
        m = k * P + z[i] if i < len(z) else (k + 1) * P + z[0]
        rr = m % P
        while TT[D, rr] == 0: D -= 1
        if g[D]: return m

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

def m0_naive(n, p, limit=10 ** 6):
    D = n - 1
    if D < 0 or good(D, p): return 0
    m = 0
    while m < limit:
        m += 1
        while Sdig(D, m, p) == 0: D -= 1
        if good(D, p): return m
    return None

def m0(n, p):
    if p in TABLE_L and n <= table(p)['P']: return m0_table(n, p)
    return m0_acc(n, p)

ARR = {}
def arr(p):
    if p not in ARR:
        P = table(p)['P']; ARR[p] = [m0_table(n, p) for n in range(P + 1)]
    return ARR[p]

# ---------------------------------------------------------------- 重整化族 (定理 7.9 的数据)
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
    def growth(s, d):   # 仅对 admissible 调用
        p = s.p; th = s.Ph(d, 0, 0, -1)[2]; idx = {}; seq = []
        while th not in idx:
            idx[th] = len(seq); r2 = s.Ph(p - 2, 0, th, p - 1)
            seq.append((th, r2[1])); th = r2[2]
        cyc = seq[idx[th]:]; rho = Fraction(sum(x for _, x in cyc), len(cyc))
        return rho, len(cyc), 2 * rho * s.Mod / (s.W + 1)

# ---------------------------------------------------------------- 其他工具
def prove_regex(p, dfa):
    st = (START, 0); seen = {st}; todo = [st]
    while todo:
        R, s = todo.pop()
        if okR(R) != (s is not None): return False, len(seen)
        for b in range(p):
            nx = (stepR(R, b, p), dfa[s].get(b) if s is not None else None)
            if nx not in seen: seen.add(nx); todo.append(nx)
    return True, len(seen)

def rank_mod(rows, MOD):
    rows = [[x % MOD for x in r] for r in rows]; rk = 0
    for c in range(len(rows[0])):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c]), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        inv = pow(rows[rk][c], MOD - 2, MOD)
        for i in range(len(rows)):
            if i != rk and rows[i][c]:
                f = rows[i][c] * inv % MOD
                rows[i] = [(x - f * y) % MOD for x, y in zip(rows[i], rows[rk])]
        rk += 1
        if rk == len(rows): break
    return rk

def orbit(p, n, limit=10 ** 6):
    """T(f) = f/x (x|f), 否则 (1+x)f - f(0); 返回 (前周期, 周期)"""
    f = tuple([1] * n); seen = {}; t = 0
    while f not in seen:
        if t > limit: return None, None
        seen[f] = t
        if f[0] == 0: f = f[1:]
        else:
            g = list(f) + [0]
            for i in range(len(f)): g[i + 1] = (g[i + 1] + f[i]) % p
            g[0] = 0
            while len(g) > 1 and g[-1] == 0: g.pop()
            f = tuple(g)
        t += 1
    return seen[f], t - seen[f]

SECTIONS = []
def section(letter, title):
    def deco(fn): SECTIONS.append((letter, title, fn)); return fn
    return deco

# ================================================================ A
@section('A', '基础：S_D(m) 算法、周期性、引理 2.2、好层级自动机、级联算法')
def secA():
    rnd = random.Random(2024)
    for p in PRIMES_SMALL:
        tb = table(p); P = tb['P']; TT = tb['TT']
        check('A0 p=%d 周期性 S_D(m+p^L)=S_D(m)，D,m<%d 逐项' % (p, P), tb['periodic'])
        bad = cnt = 0; Dmax = min(P, max(1, 300000 // P))
        for D in range(Dmax):
            for m in range(P):
                cnt += 1
                if Sdig(D, m, p) != int(TT[D, m]): bad += 1
        for _ in range(20000):
            D = rnd.randrange(P); m = rnd.randrange(P); cnt += 1
            if Sdig(D, m, p) != int(TT[D, m]): bad += 1
        check('A1 p=%d 数字算法 = Pascal 表' % p, bad == 0, '%d 次比对, 不符 %d' % (cnt, bad))
        bad = 0
        for _ in range(150):
            i = rnd.randint(1, 3); q = rnd.randrange(0, max(1, 4000 // p ** i))
            y = rnd.randrange(p ** i); D = q * p ** i + y
            t = rnd.randrange(p ** 3); r = rnd.randrange(p ** i); Mm = t * p ** i + r
            lhs = sum(binom_mod(Mm, j, p) for j in range(min(D, Mm) + 1)) % p
            rhs = (pow(2, r, p) * Sdig(q - 1, t, p) + binom_mod(t, q, p) * Sdig(y, r, p)) % p
            if lhs != rhs or Sdig(D, Mm, p) != lhs: bad += 1
        check('A2 p=%d 引理 2.2（Lucas 直接求和，150 例）' % p, bad == 0, '不符 %d' % bad)
        bad = sum(1 for D in range(P) if good(D, p) != bool(tb['good'][D]))
        check('A3 p=%d 好层级自动机 = 定义（D<%d）' % (p, P), bad == 0,
              '好层级数 %d, 不符 %d' % (int(tb['good'].sum()), bad))
        bad = 0; lim = min(P, 400)
        for n in range(1, lim + 1):
            a1 = m0_table(n, p)
            if m0_acc(n, p) != a1: bad += 1
            if n <= 120 and m0_naive(n, p) != a1: bad += 1
        check('A4 p=%d 三种级联算法一致（n≤%d；朴素法 n≤120）' % (p, lim), bad == 0, '不符 %d' % bad)

# ================================================================ B
@section('B', '好层级结构：正则式、前缀封闭、偶数位、单调、无相邻')
def secB():
    dfas = {3: ('[02]*', {0: {0: 0, 2: 0}}),
            5: ('[04]*(24*(24*)?)?', {0: {0: 0, 4: 0, 2: 1}, 1: {4: 1, 2: 2}, 2: {4: 2}}),
            7: ('[06]*(4[06]*)?', {0: {0: 0, 6: 0, 4: 1}, 1: {0: 1, 6: 1}})}
    for p, (name, dfa) in dfas.items():
        ok, ns = prove_regex(p, dfa)
        check('B1 p=%d 好层级 ⇔ %s（乘积自动机穷举）' % (p, name), ok, '%d 个状态对' % ns)
    for p in (3, 5, 7, 11, 13, 17, 19, 23):
        t1 = time.time(); seen = {START}; q = [START]; i = 0; cap = 30000
        pc = od = 0
        while i < len(q) and len(seen) < cap and time.time() - t1 < 60:
            R = q[i]; i += 1
            for b in range(p):
                N = stepR(R, b, p)
                if not okR(R) and okR(N): pc += 1
                if okR(R) and b % 2 == 1 and okR(N): od += 1
                if N not in seen: seen.add(N); q.append(N)
        complete = i >= len(q)
        res = (pc == 0 and od == 0) if complete else (None if pc == od == 0 else False)
        check('B2 p=%d 前缀封闭 & 奇数位必坏（未剪枝自动机）' % p, res,
              '状态 %d%s, 前缀反例 %d, 奇数位反例 %d' % (len(seen), '' if complete else ' [未穷尽]', pc, od))
    for p in (3, 5, 7, 11, 13, 17, 19, 23):
        seen = {STARTN}; q = [STARTN]; i = 0; bad = 0
        while i < len(q) and len(seen) < 30000:
            S = q[i]; i += 1
            for b in range(p):
                N = stepN(S, b, p)
                if N is None: continue
                if b % 2 == 0 and not S <= N: bad += 1
                if N not in seen: seen.add(N); q.append(N)
        complete = i >= len(q)
        check('B3 p=%d 追加偶数位时射影状态单调' % p, (bad == 0) if complete else None,
              '好状态 %d%s, 反例 %d' % (len(seen), '' if complete else ' [未穷尽]', bad))
    for p in PRIMES_SMALL:
        lim = min(p ** 7, 200000)
        adj = [D for D in range(lim) if good(D, p) and good(D + 1, p)]
        check('B4 p=%d 无相邻好层级（D<%d，数值）' % (p, lim), not adj, str(adj[:5]))

# ================================================================ C
@section('C', 'p=3：定理 P3 的精确分解')
def secC():
    p = 3; a = arr(3); P = table(3)['P']
    vn = vh = vd = vl = v1 = vdev = tot = 0; H = {}
    for n in range(2, P + 1):
        N = n - 1
        if good(N, 3): continue
        tot += 1; ds = digs(N, 3)
        if 1 not in ds: v1 += 1; continue
        i = len(ds) - 1 - ds.index(1); q = N // 3 ** i; y = N % 3 ** i
        if nu_of(n, 3) != q * 3 ** i: vn += 1
        h = a[n] - q * 3 ** i
        if not (1 <= h <= 3 ** i): vh += 1
        if a[3 ** i + y + 1] - 3 ** i != h: vd += 1
        H[(i, y)] = h
        if not ((n - 1) / 2 < a[n] <= 2 * (n - 1)): vl += 1
        if abs(a[n] - (n - 1)) > 3 ** i: vdev += 1
    check('C1 非好层级都含数字 1', v1 == 0, '共 %d 个非好 n-1, 反例 %d' % (tot, v1))
    check('C2 ν(n)=q·3^i', vn == 0, 'viol=%d' % vn)
    check('C3 1≤h≤3^i', vh == 0, 'viol=%d' % vh)
    check('C4 h 只依赖 y（与 q 无关）', vd == 0, 'viol=%d' % vd)
    check('C5 (n-1)/2 < m0 ≤ 2(n-1)', vl == 0, 'viol=%d' % vl)
    check('C6 |m0-(n-1)| ≤ 3^i', vdev == 0, 'viol=%d' % vdev)
    for i in range(0, 4):
        log('    H_%d(y), y=0..%d: %s' % (i, 3 ** i - 1, [H.get((i, y)) for y in range(3 ** i)]))

# ================================================================ D
@section('D', '全局界：缩放律、定理 1、线性常数、纪元数 E(n)')
def secD():
    for p in PRIMES_SMALL:
        a = arr(p); P = table(p)['P']
        sc = sum(1 for n in range(1, P // p + 1) if a[p * n] != p * a[n])
        check('D1 p=%d 缩放律 m0(pn)=p·m0(n)，n≤%d' % (p, P // p), sc == 0, 'viol=%d' % sc)
        v1 = vE = v3 = v2 = 0; r1 = 0.0; ml = (0, 0); mE = (0, 0); Eset = set()
        for n in range(2, P + 1):
            if good(n - 1, p): continue
            L = plen(n - 1, p); m = a[n]; B = (p - 1) * L * p ** L
            if m > B: v1 += 1
            r1 = max(r1, m / B)
            if m / (n - 1) > ml[0]: ml = (m / (n - 1), n)
            E = m // p ** L; Eset.add(E)
            if E > mE[0]: mE = (E, n)
            if E > (p - 1) * L: vE += 1
            if m > (p + 1) * (p - 1) * p ** L: v3 += 1
            if m > 2 * (n - 1): v2 += 1
        check('D2 p=%d 定理 1：m0 ≤ (p-1)·L·p^L' % p, v1 == 0, 'viol=%d, max m0/界=%.4f' % (v1, r1))
        check('D3 p=%d 经验常数 max m0/(n-1)' % p, None, '%.4f at n=%d (n≤%d)' % (ml[0], ml[1], P))
        check('D4 p=%d E(n)=⌊m0/p^L⌋ ≤ (p-1)L' % p, vE == 0, 'max E=%d at n=%d, E 取值 %s' % (mE[0], mE[1], sorted(Eset)[:12]))
        if p in (3, 7):
            check('D5 p=%d 定理 3 界 m0 ≤ (p+1)(p-1)p^L' % p, v3 == 0, 'viol=%d' % v3)
        if p == 3:
            check('D6 p=3 m0 ≤ 2(n-1)', v2 == 0, 'viol=%d' % v2)
            check('D7 p=3 E(n)∈{0,1}', Eset <= {0, 1}, str(sorted(Eset)))
        if p == 7:
            check('D8 p=7 m0 ≤ 42(n-1)（数值）', ml[0] <= 42, 'max ratio=%.4f' % ml[0])

# ================================================================ E
@section('E', '重整化族：族公式 vs 直接级联；刚性二分；DEEP；ρ>0')
def secE():
    WMs = {3: 80, 5: 124, 7: 48, 11: 120, 13: 60}
    WBs = {3: 729, 5: 3125, 7: 343, 11: 121, 13: 169}
    if ARGS.full: WBs.update({7: 2401, 11: 1331, 13: 2197})
    for p in PRIMES_SMALL:
        P = table(p)['P']; lim = 12000 if ARGS.full else P
        comp = agr = 0; dis = []
        for W in range(1, WMs[p] + 1):
            if not good(W, p): continue
            F = Fam(p, W)
            for d in range(p - 1):
                for k in range(1, 9):
                    n = (W + 1) * p ** k - p + d + 1
                    if n > lim: break
                    v, ty = F.m0_formula(d, k)
                    if v is None: continue
                    b = m0(n, p); comp += 1
                    if v == b: agr += 1
                    elif len(dis) < 5: dis.append((W, d, k, v, b, ty))
        check('E1 p=%d 族公式 = 直接级联（W≤%d, n≤%d）' % (p, WMs[p], lim),
              (agr == comp) if comp else None, '%d/%d 一致 %s' % (agr, comp, dis))
        t1 = time.time(); nW = bad = deep = rb = w1 = 0; rig = []; cc = {}
        for W in range(1, WBs[p]):
            if not good(W, p): continue
            nW += 1; F = Fam(p, W); r = isrigid(p, W)
            cls = [F.classify(d) for d in range(p - 1)]
            for c in cls: cc[c] = cc.get(c, 0) + 1
            deep += cls.count('deep')
            alladm = all(c == 'admissible' for c in cls); anyadm = any(c == 'admissible' for c in cls)
            if (r and not alladm) or ((not r) and anyadm): bad += 1
            if r and good(W - 1, p): w1 += 1
            if r and alladm:
                gs = [F.growth(d) for d in range(p - 1)]
                if any(g[0] <= 0 for g in gs): rb += 1
                rig.append((W, digs(W, p), str(gs[0][0]), gs[0][1], str(gs[0][2]), str(max(g[2] for g in gs))))
        check('E2 p=%d 刚性 ⇔ 全部 d 可容许（%d 个好层级 W<%d）' % (p, nW, WBs[p]), bad == 0,
              'viol=%d, 分类计数 %s (%.0fs)' % (bad, cc, time.time() - t1))
        check('E3 p=%d DEEP 类为空' % p, deep == 0, 'deep=%d' % deep)
        check('E4 p=%d 刚性族 ρ>0' % p, (rb == 0) if rig else None, '刚性 W 个数 %d, ρ≤0 的 %d' % (len(rig), rb))
        check('E5 p=%d 刚性 W 满足 W-1∉G' % p, w1 == 0, 'viol=%d' % w1)
        for x in rig[:8]:
            log('    刚性 W=%s %s  ρ(d=0)=%s 循环长=%s c(d=0)=%s max_d c=%s' % x)

# ================================================================ F
@section('F', 'p=5 主定理 m0(13·5^k-e) = (5⌊k/2⌋+4)·5^{k+1}')
def secF():
    F = Fam(5, 12); bad = 0
    for k in range(1, 81):
        for e in range(1, 5):
            v, ty = F.m0_formula(4 - e, k)
            if v != (5 * (k // 2) + 4) * 5 ** (k + 1) or ty != 'adm': bad += 1
    check('F1 族公式 k≤80, e=1..4', bad == 0, 'viol=%d' % bad)
    bad = c = 0
    for k in range(0, 6):
        es = [1, 2, 3, 4] if (k <= 3 or ARGS.full) else [1]
        for e in es:
            mm = m0(13 * 5 ** k - e, 5); c += 1
            if mm != (5 * (k // 2) + 4) * 5 ** (k + 1):
                bad += 1; log('    不符 k=%d e=%d m0=%s' % (k, e, mm))
    check('F2 直接级联 k≤5', bad == 0, '%d 例, 不符 %d' % (c, bad))
    g = F.growth(0)
    check('F3 W=12 增长数据', None, 'ρ=%s, 循环长=%d, c=%s' % (g[0], g[1], g[2]))

# ================================================================ G
@section('G', '素数普查：刚性层级、最小稳定对（定理 3）、p=11/13 具体例')
def secG():
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41):
        t1 = time.time(); rep = {STARTN: ()}; q = [STARTN]; i = 0; trans = {}
        while i < len(q) and len(rep) < 20000 and time.time() - t1 < 90:
            S = q[i]; i += 1; tr = {}
            for b in range(p):
                N = stepN(S, b, p)
                if N is not None:
                    tr[b] = N
                    if N not in rep: rep[N] = rep[S] + (b,); q.append(N)
            trans[S] = tr
        complete = i >= len(q)
        rigW = []; ms_low = []; ms_top = []
        for S, tr in trans.items():
            W = val(rep[S], p); al = sorted(tr)
            if W >= 1 and all(b == p - 1 for b in al): rigW.append(W)
            if al:
                c = al[0]
                if c > 0 and tr[c] == S:
                    (ms_top if c == p - 1 else ms_low).append((W, c))
        Wmin = min(rigW) if rigW else None
        if rigW: verdict = '超线性（存在刚性层级，定理 2）'
        elif not ms_low and not ms_top: verdict = '全局线性（无最小稳定对，定理 3）'
        else: verdict = '未定（定理 3 的缺口）'
        if not complete: verdict += ' [状态未穷尽]'
        log('  p=%d: 好状态 %d%s; 刚性状态 %d; 最小刚性 W=%s %s; 最小稳定对 c<p-1: %d 个 %s; c=p-1: %d 个 → %s'
            % (p, len(rep), '' if complete else '(截断)', len(rigW), Wmin, digs(Wmin, p) if Wmin else '',
               len(ms_low), ms_low[:3], len(ms_top), verdict))
        if Wmin is not None and Wmin < 10 ** 6:
            check('G1 p=%d 最小刚性 W=%d 用定义复核' % (p, Wmin), isrigid(p, Wmin))
        SUMMARY.append(('INFO', 'G p=%d 判别' % p, verdict))
    for W in (822, 1064):
        check('G2 p=11 W=%d %s 是好层级且刚性' % (W, digs(W, 11)), isrigid(11, W),
              'good=%s' % good(W, 11))
    for W in (822, 1064):
        if not good(W, 11): continue
        F = Fam(11, W)
        cls = [F.classify(d) for d in range(10)]
        log('    p=11 W=%d 各 d 分类: %s' % (W, cls))
        if all(c == 'admissible' for c in cls):
            g = F.growth(0); log('    p=11 W=%d d=0: ρ=%s 循环长=%d c=%s' % (W, g[0], g[1], g[2]))
        cases = [(0, 1)] + ([(9, 1), (0, 2)] if ARGS.full and W == 822 else [])
        for d, k in cases:
            n = (W + 1) * 11 ** k - 11 + d + 1; t1 = time.time()
            direct = m0_acc(n, 11); v, ty = F.m0_formula(d, k)
            check('G3 p=11 W=%d d=%d k=%d (n=%d) 公式=直接' % (W, d, k, n), direct == v,
                  'direct=%s formula=%s (%.0fs)' % (direct, v, time.time() - t1))
    rl = [W for W in range(1, 169) if isrigid(13, W)]
    log('    p=13 刚性 W<169: %s' % [(W, digs(W, 13)) for W in rl])
    for W in rl[:2]:
        F = Fam(13, W)
        for k in (1, 2):
            for d in (0, 11):
                n = (W + 1) * 13 ** k - 13 + d + 1
                if n > 20000 and not ARGS.full: continue
                v, ty = F.m0_formula(d, k); b = m0(n, 13)
                check('G4 p=13 W=%d d=%d k=%d 公式=直接' % (W, d, k), v == b, 'direct=%s formula=%s' % (b, v))

# ================================================================ H
@section('H', 'p-核秩（正则性的数值证据，不构成证明）')
def secH():
    MOD = 1000003
    for p, Es in ((3, (3,)), (5, (2,)), (7, (1, 2))):
        a = arr(p); P = table(p)['P']
        for E in Es:
            C = P // p ** E
            rows = [[a[p ** e * n + i] for n in range(C)] for e in range(E + 1) for i in range(p ** e)]
            cum = []; cnt = 0
            for e in range(E + 1):
                cnt += p ** e; cum.append(rank_mod(rows[:cnt], MOD))
            m2 = [len({tuple(a[p ** e * n + i] % 2 for n in range(C)) for i in range(p ** e)}) for e in range(E + 1)]
            check('H p=%d e≤0..%d, %d 列' % (p, E, C), None,
                  '累计秩 %s (行数 %s); mod 2 每层不同核序列数 %s' % (cum, [sum(p ** j for j in range(e + 1)) for e in range(E + 1)], m2))

# ================================================================ I
@section('I', 'p=5 精确不等式 m0 ≤ n⌈log5(n+1)⌉+1 的例外集')
def secI():
    a = arr(5); P = table(5)['P']; ex = []
    for n in range(2, P + 1):
        e = 0
        while 5 ** e < n + 1: e += 1
        if a[n] > n * e + 1: ex.append(n)
    check('I p=5 n≤%d 例外集 = {2,9}' % P, ex == [2, 9], '实际例外 %s' % ex[:20])

# ================================================================ J
@section('J', '多项式轨道 vs 定理 4.7（μ=2m0+n-ν-1, 周期 2p^w）')
def secJ():
    for p, NM in ((3, 40), (5, 40), (7, 30), (11, 20), (13, 20)):
        bad = 0; first = None
        for n in range(1, NM + 1):
            mu, per = orbit(p, n)
            nu = nu_of(n, p)
            emu = 0 if good(n - 1, p) else 2 * m0(n, p) + n - nu - 1
            w = 1
            while w < nu: w *= p
            if mu != emu or per != 2 * w:
                bad += 1
                if first is None: first = (n, mu, emu, per, 2 * w)
        check('J p=%d n≤%d 轨道与定理 4.7' % (p, NM), bad == 0,
              '不符 %d, 首例 (n, μ实测, μ公式, 周期实测, 周期公式)=%s' % (bad, first))

# ================================================================ main
log('verify_v5.py  Python %s  numpy %s  full=%s' % (platform.python_version(), np.__version__, ARGS.full))
for letter, title, fn in SECTIONS:
    if ONLY and letter not in ONLY: continue
    log(''); log('=' * 72); log('[%s] %s   (t=%.1fs)' % (letter, title, time.time() - T0))
    try: fn()
    except Exception:
        log(traceback.format_exc()); check('%s 节运行' % letter, False, '异常，见上方 traceback')
log(''); log('=' * 72); log('SUMMARY   (总耗时 %.1fs)' % (time.time() - T0))
for tag in ('FAIL', 'PASS', 'INFO'):
    items = [x for x in SUMMARY if x[0] == tag]
    log('%s: %d' % (tag, len(items)))
    if tag != 'PASS':
        for _, nm, dt in items: log('   - %s  %s' % (nm, dt))
LOGF.close()
