import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m0_regularity import Cascade

cas = Cascade(11)

def mu_of(p, ds):
    N = 0
    for d in ds:
        N = N * p + d
    n = N + 1
    m, dr = cas.m0(n)
    v = cas.nu(n)
    d = n - v
    mu = 2 * m + d - 1 if d else 0
    return n, m, d, dr, mu / n

print('### p=11 scan: families [c, d^j, e] with growing mu/n (j=0..4)')
found = 0
for c in range(2, 10):
    for d in (5, 6, 8, 9, 10):
        for e in (0, 1, 4, 9, 10):
            vals = [mu_of(11, [c] + [d] * j + [e]) for j in range(5)]
            r4, r0 = vals[4][4], vals[0][4]
            if r4 > 28 or r4 > r0 + 5:
                found += 1
                print('  [%d, %d^j, %d]: mu/n = %s  (max %.3f)' % (c, d, e, ['%.2f' % v[4] for v in vals], max(v[4] for v in vals)))
print('  total growing families found: %d' % found)

print()
print('### p=11: [6, 8, 8, 9, 1^b] + 1')
for b in range(5):
    n, m, d, dr, r = mu_of(11, [6, 8, 8, 9] + [1] * b)
    print('  b=%d n=%d  m0=%d  d=%d  mu/n=%.6f  drops=%d' % (b, n, m, d, r, dr))

print()
print('### dichotomy test on unified family [2, (p-1)^j, p-2] + 1')
for p in (19, 23, 29, 31):
    cc = Cascade(p)
    ordp = 1; x = 2 % p
    while x != 1:
        x = x * 2 % p; ordp += 1
    print('  p=%d  ord_p(2)=%d  -1 in <2>: %s' % (p, ordp, ordp % 2 == 0))
    for j in range(0, 8):
        N = 0
        for d in [2] + [p - 1] * j + [p - 2]:
            N = N * p + d
        n = N + 1
        m, dr = cc.m0(n)
        v = cc.nu(n)
        dd = n - v
        mu = 2 * m + dd - 1 if dd else 0
        pred = (3 * p * j + 4 * p - 1) * p ** (j + 1)
        print('    j=%d  n=%d  m0=%d  unif-law=%s  mu/n=%.4f  drops=%d' % (j, n, m, m == pred, mu / n, dr))
