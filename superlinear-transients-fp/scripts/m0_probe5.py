import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m0_regularity import Cascade

def fam(p, ds_fn, jmax, label):
    cas = Cascade(p)
    print('p=%d  %s' % (p, label))
    for j in range(jmax + 1):
        ds = ds_fn(j)
        N = 0
        for d in ds:
            N = N * p + d
        n = N + 1
        m, dr = cas.m0(n)
        v = cas.nu(n)
        d = n - v
        mu = 2 * m + d - 1 if d else 0
        pred = (3 * p * j + 4 * p - 1) * p ** (j + 1)
        print('  j=%2d  n=%d  m0=%d  pred=(3pj+4p-1)p^(j+1)=%d  match=%s  d=%d (pred 2p^(j+1)-1=%d)  mu/n=%.6f  drops=%d'
              % (j, n, m, pred, m == pred, d, 2 * p ** (j + 1) - 1, mu / n, dr))

print('### p=5: [2, 4^j, 3] + 1  (unified form n = 3*5^(j+1) - 1)')
fam(5, lambda j: [2] + [4] * j + [3], 5, '[2, 4^j, 3]')

print('### p=11: [2, 10^j, 9] + 1  (unified form n = 3*11^(j+1) - 1)')
fam(11, lambda j: [2] + [10] * j + [9], 5, '[2, 10^j, 9]')

print('### p=7: [2, 6^j, 5] + 1')
fam(7, lambda j: [2] + [6] * j + [5], 5, '[2, 6^j, 5]')

print('### p=13 sanity: [2, 12^j, 11] j=9,10')
fam(13, lambda j: [2] + [12] * j + [11], 10, '[2, 12^j, 11]')

print('### p=17 sanity: [2, 16^j, 15] j=7,8')
fam(17, lambda j: [2] + [16] * j + [15], 8, '[2, 16^j, 15]')
