import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m0_regularity import Cascade

def fam(p, tmax):
    cas = Cascade(p)
    q = p - 1
    print('p=%d  n_t = [2,2,%d^(2t-1),0]_p + 1' % (p, q))
    print('  ord_p(2)=%d  -1 in <2>: %s' % (ord2(p), -1 in grp2(p)))
    for t in range(1, tmax + 1):
        ds = [2, 2] + [q] * (2 * t - 1) + [0]
        N = 0
        for d in ds:
            N = N * p + d
        n = N + 1
        m, dr = cas.m0(n)
        v = cas.nu(n)
        d = n - v
        mu = 2 * m + d - 1 if d else 0
        print('  t=%2d n=%d  m0=%d  m0/n=%.4f  d=%d  mu/n=%.6f  drops=%d' % (t, n, m, m / n, d, mu / n, dr))

def ord2(p):
    x, k = 1, 0
    while True:
        x = x * 2 % p; k += 1
        if x == 1:
            return k

def grp2(p):
    g, x = set(), 1
    for _ in range(ord2(p)):
        g.add(x); x = x * 2 % p
    return g

fam(7, 4)
fam(11, 5)
fam(13, 5)
fam(17, 5)
