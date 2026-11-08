import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m0_regularity import Cascade, ndig

def traj(p, ds):
    cas = Cascade(p)
    N = 0
    for d in ds:
        N = N * p + d
    n = N + 1
    F = cas.nu(n) - 1
    D = n - 1
    K = ndig(n - 1, p)
    m = 0; steps = []
    while D > F:
        t = cas.next_zero(D, m, K)
        Dn = cas.top_nonzero(D - 1, t, K)
        steps.append((t, Dn, Dn - D))
        m = t; D = Dn
    print('p=%d n=%d digits=%s  F(nu-1)=%d  drops=%d' % (p, n, ds, F, len(steps)))
    for i, (T, D, dec) in enumerate(steps):
        print('   %2d: T=%d  D=%d  dD=%d' % (i + 1, T, D, dec))
    print()

print('### p=13: [2,12^j,11] trajectories')
traj(13, [2, 11])
traj(13, [2, 12, 11])
traj(13, [2, 12, 12, 11])
print('### p=19: [2,18^j,17] trajectories')
traj(19, [2, 17])
traj(19, [2, 18, 17])
