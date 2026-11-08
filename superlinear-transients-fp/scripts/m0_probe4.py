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
        print('  j=%2d  n=%d  digits=%s  m0=%d  d=%d  mu/n=%.6f  drops=%d' % (j, n, ds, m, d, mu / n, dr))

print('### p=13: [2, 12^j, 11] + 1')
fam(13, lambda j: [2] + [12] * j + [11], 8, '[2, 12^j, 11]')

print('### p=17: [2, 16^j, 15] + 1')
fam(17, lambda j: [2] + [16] * j + [15], 6, '[2, 16^j, 15]')

print('### p=11: [6, 8^j, 9] + 1')
fam(11, lambda j: [6] + [8] * j + [9], 6, '[6, 8^j, 9]')

print('### p=11: [6, 8^j, 4] + 1')
fam(11, lambda j: [6] + [8] * j + [4], 6, '[6, 8^j, 4]')

print('### p=5: [2, 4^j, 0] + 1  (single leading 2)')
fam(5, lambda j: [2] + [4] * j + [0], 7, '[2, 4^j, 0]')

print('### p=5: [2, 2, 4^j, 0] + 1  (double leading 2, j = 2t-1)')
fam(5, lambda t: [2, 2] + [4] * (2 * t - 1) + [0], 5, '[2, 2, 4^(2t-1), 0]')
