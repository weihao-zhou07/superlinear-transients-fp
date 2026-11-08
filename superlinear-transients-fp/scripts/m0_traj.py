import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m0_regularity import Cascade, ndig

def traj(p, n):
    cas = Cascade(p)
    F = cas.nu(n) - 1
    D = n - 1
    if D <= F:
        return [], 0, F
    K = ndig(n - 1, p)
    m = 0; steps = []
    while D > F:
        t = cas.next_zero(D, m, K)
        Dn = cas.top_nonzero(D - 1, t, K)
        steps.append((t, Dn))
        m = t; D = Dn
    return steps, m, F

print('p=5  n_t = 13*5^(2t)-4 : trajectories (T_i = drop time, D_i = level after i-th drop)')
trajs = {}
for t in range(1, 5):
    n = 13 * 5 ** (2 * t) - 4
    st, m0, F = traj(5, n)
    trajs[t] = st
    print('t=%d n=%d  nu-1=%d  m0=%d  drops=%d' % (t, n, F, m0, len(st)))
    for i, (T, D) in enumerate(st):
        print('   %2d: T=%d  D=%d' % (i + 1, T, D))

print()
print('self-similarity check: is T_i(t+1) = 25*T_i(t) + c  and D_i(t+1) = 25*D_i(t) + d  for i <= 8t ?')
for t in range(1, 4):
    a, b = trajs[t], trajs[t + 1]
    okT = okD = True
    cT, cD = None, None
    for i in range(len(a)):
        Ti, Di = a[i]
        Tj, Dj = b[i]
        if cT is None:
            cT = Tj - 25 * Ti
            cD = Dj - 25 * Di
        if Tj - 25 * Ti != cT:
            okT = False
        if Dj - 25 * Di != cD:
            okD = False
    print('t=%d->%d:  first %d drops: T: 25x+%d (exact=%s)   D: 25x+%d (exact=%s)   extra drops of t+1: %d'
          % (t, t + 1, len(a), cT, okT, cD, okD, len(b) - len(a)))
    for i in range(len(a), len(b)):
        print('   new drop %2d of t+1: T=%d  D=%d' % (i + 1, b[i][0], b[i][1]))
