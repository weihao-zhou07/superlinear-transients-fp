import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m0_regularity import Cascade

cas = {}

def get(p):
    if p not in cas:
        cas[p] = Cascade(p)
    return cas[p]

def family(p, ds_fn, kmax, label):
    c = get(p)
    print('p=%d  %s' % (p, label))
    for k in range(1, kmax + 1):
        ds = ds_fn(k)
        N = 0
        for d in ds:
            N = N * p + d
        n = N + 1
        m, dr = c.m0(n)
        v = c.nu(n)
        d = n - v
        mu = 2 * m + d - 1 if d else 0
        print('  k=%2d  n=%d  m0=%d  m0/n=%.4f  d=%d  mu/n=%.4f  drops=%d' % (k, n, m, m / n, d, mu / n, dr))

print('=== p=5: n = [2,2,4^j,0] + 1  (j = 2t-1 odd) ===')
family(5, lambda t: [2, 2] + [4] * (2 * t - 1) + [0], 10, 'j = 2t-1, one trailing zero')

print('=== p=5: n = c*5^k + 1 ===')
for cc in (2, 3, 4):
    family(5, lambda k, cc=cc: [cc] + [0] * k, 12, 'c=%d' % cc)

print('=== p=7: n = c*7^k + 1 ===')
for cc in (2, 3, 4, 5, 6):
    family(7, lambda k, cc=cc: [cc] + [0] * k, 12, 'c=%d' % cc)

print('=== p=7: n = [1,1,6^j,0] + 1 ===')
family(7, lambda j: [1, 1] + [6] * j + [0], 8, 'j varies')

print('=== p=11: n = c*11^k + 1 ===')
for cc in (2, 3, 5):
    family(11, lambda k, cc=cc: [cc] + [0] * k, 10, 'c=%d' % cc)

print('=== p=13: n = c*13^k + 1 ===')
for cc in (2, 3, 5):
    family(13, lambda k, cc=cc: [cc] + [0] * k, 10, 'c=%d' % cc)
