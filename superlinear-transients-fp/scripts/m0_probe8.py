import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m0_regularity import Cascade

cas = {}

def get(p):
    if p not in cas:
        cas[p] = Cascade(p)
    return cas[p]

def fam(p, ds_fn, jmax, label):
    c = get(p)
    print('p=%d  %s' % (p, label))
    for j in range(jmax + 1):
        ds = ds_fn(j)
        N = 0
        for d in ds:
            N = N * p + d
        n = N + 1
        m, dr = c.m0(n)
        v = c.nu(n)
        d = n - v
        mu = 2 * m + d - 1 if d else 0
        print('  j=%2d  n=%d  m0=%d  d=%d  mu/n=%.6f  drops=%d' % (j, n, m, d, mu / n, dr))

print('### p=23: [6, 20^j, 21] + 1')
fam(23, lambda j: [6] + [20] * j + [21], 5, '[6,20^j,21]')
print('### p=29: [6, 28^j, 26] + 1')
fam(29, lambda j: [6] + [28] * j + [26], 5, '[6,28^j,26]')
print('### p=31: [6, 30^j, 29] + 1')
fam(31, lambda j: [6] + [30] * j + [29], 5, '[6,30^j,29]')

print()
print('### p=7 full scan: [c, d^j, e], growing families (j=0..5)')
c7 = get(7)
best = {}
for cc in range(2, 7):
    for dd in range(7):
        for ee in range(7):
            vals = []
            for j in range(6):
                N = 0
                for dgt in [cc] + [dd] * j + [ee]:
                    N = N * 7 + dgt
                n = N + 1
                m, dr = c7.m0(n)
                v = c7.nu(n)
                d = n - v
                mu = 2 * m + d - 1 if d else 0
                vals.append(mu / n)
            best[(cc, dd, ee)] = max(vals)
mx = max(best.values())
print('  global max mu/n over scan: %.4f' % mx)
for kk, v in sorted(best.items(), key=lambda kv: -kv[1])[:12]:
    print('    [%d, %d^j, %d]: max mu/n = %.4f' % (kk[0], kk[1], kk[2], v))
