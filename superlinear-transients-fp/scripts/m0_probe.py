import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m0_regularity import Cascade

def grid(p, J, K):
    cas = Cascade(p)
    print('p=%d  grid over j in 1..%d, k in 0..%d  (pattern [2,2,4^j,0^k]+1)' % (p, J, K))
    best = (0, None)
    for j in range(1, J + 1):
        row = []
        for k in range(0, K + 1):
            ds = [2, 2] + [4] * j + [0] * k
            N = 0
            for d in ds:
                N = N * p + d
            n = N + 1
            m, dr = cas.m0(n)
            v = cas.nu(n)
            d = n - v
            mu = 2 * m + d - 1 if d else 0
            r = mu / n if n else 0
            row.append('k=%d:%.4f' % (k, r))
            if r > best[0]:
                best = (r, (j, k, n, m, d, dr))
        print('  j=%d  ' % j + '  '.join(row))
    print('  best: %.4f at j,k=%s n=%d m0=%d d=%d drops=%d' % (best[0], best[1][:2], best[1][2], best[1][3], best[1][4], best[1][5]))

grid(5, 5, 5)
grid(7, 4, 4)
