#!/usr/bin/env python3
# hankel2.py -- two-parameter Hankel-rank test for p-regularity of m0
# usage: python hankel2.py P K NRAND     e.g.  python hankel2.py 5 9 30
import sys, random, time
from m0_regularity import Cascade, MOD

def num(p, ds):
    N = 0
    for d in ds:
        N = N * p + d
    return N

def rank_mod(M):
    A = [[x % MOD for x in r] for r in M]; r = 0
    for c in range(len(A[0])):
        piv = next((i for i in range(r, len(A)) if A[i][c]), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = pow(A[r][c], MOD - 2, MOD)
        A[r] = [x * inv % MOD for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [(x - f * y) % MOD for x, y in zip(A[i], A[r])]
        r += 1
        if r == len(A):
            break
    return r

def run(p, K, nrand, seed=7):
    cas = Cascade(p); rnd = random.Random(seed); q = p - 1
    fams = [([2, 2], [q, q], [], [q, q], [0]),        # around the extremal family
            ([2, 2], [q, q], [0], [q, q], [0]),
            ([2, 2], [q, q], [1], [q, q], [0]),
            ([2, 2], [q, q], [2], [q], [0]),
            ([2], [q], [0], [q], [1]),
            ([1], [0], [], [q], []),
            ([2], [0], [], [q], [])]
    while nrand > 0:
        x = [rnd.randrange(1, p)] + [rnd.randrange(p) for _ in range(rnd.randrange(2))]
        y = [rnd.randrange(p) for _ in range(rnd.randrange(1, 3))]
        m = [rnd.randrange(p) for _ in range(rnd.randrange(2))]
        w = [rnd.randrange(p) for _ in range(rnd.randrange(1, 3))]
        z = [rnd.randrange(p) for _ in range(rnd.randrange(2))]
        fams.append((x, y, m, w, z)); nrand -= 1
    for x, y, m, w, z in fams:
        t0 = time.time()
        H = [[cas.m0(num(p, x + y * k + m + w * l + z) + 1)[0] for l in range(K)]
             for k in range(K)]
        prof = [rank_mod([row[:j] for row in H[:j]]) for j in range(2, K + 1)]
        tag = 'GROWING' if prof[-1] >= prof[len(prof) // 2] + 2 else 'stable'
        print('x=%s y=%s m=%s w=%s z=%s  rank(j=2..%d)=%s  [%s] %.1fs'
              % (x, y, m, w, z, K, prof, tag, time.time() - t0))

if __name__ == '__main__':
    p, K, nr = (int(a) for a in sys.argv[1:4])
    run(p, K, nr)
