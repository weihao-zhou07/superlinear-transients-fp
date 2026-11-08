import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m0_regularity import Cascade

c5, c7, c11, c13 = Cascade(5), Cascade(7), Cascade(11), Cascade(13)

print('p=5  n = 13*5^(2t)-4 : predicted m0 = (5t+4)*5^(2t+1), d = 8*5^(2t)-4, mu = (50t+48)*5^(2t)-5')
for t in (11, 12, 13):
    n = 13 * 5 ** (2 * t) - 4
    m, dr = c5.m0(n)
    v = c5.nu(n); d = n - v
    pred_m = (5 * t + 4) * 5 ** (2 * t + 1)
    pred_d = 8 * 5 ** (2 * t) - 4
    pred_mu = (50 * t + 48) * 5 ** (2 * t) - 5
    mu = 2 * m + d - 1
    print('  t=%2d n=%d  m0 match=%s  d match=%s  mu match=%s  drops=%d' % (t, n, m == pred_m, d == pred_d, mu == pred_mu, dr))

print('p=7  n = 3*7^k+1 : predicted m0 = 13*7^k, d = 2*7^k+1, mu = 28*7^k')
for k in (13, 14, 15):
    n = 3 * 7 ** k + 1
    m, dr = c7.m0(n)
    v = c7.nu(n); d = n - v
    mu = 2 * m + d - 1
    print('  k=%2d n=%d  m0=%d (pred %d, %s)  d=%d (pred %d, %s)  mu=%d (pred %d, %s)  drops=%d'
          % (k, n, m, 13 * 7 ** k, m == 13 * 7 ** k, d, 2 * 7 ** k + 1, d == 2 * 7 ** k + 1,
             mu, 28 * 7 ** k, mu == 28 * 7 ** k, dr))

print('p=7  n = 2*7^k+1 : predicted m0 = 6*7^k, d = 7^k+1, mu = 13*7^k')
for k in (13, 14, 15):
    n = 2 * 7 ** k + 1
    m, dr = c7.m0(n)
    v = c7.nu(n); d = n - v
    mu = 2 * m + d - 1
    print('  k=%2d n=%d  m0=%d (pred %d, %s)  d=%d (pred %d, %s)  mu=%d (pred %d, %s)  drops=%d'
          % (k, n, m, 6 * 7 ** k, m == 6 * 7 ** k, d, 7 ** k + 1, d == 7 ** k + 1,
             mu, 13 * 7 ** k, mu == 13 * 7 ** k, dr))

print('p=11  n = 3*11^k+1 : predicted m0 = 10*11^k, d = 2*11^k+1, mu = 22*11^k')
for k in (11, 12):
    n = 3 * 11 ** k + 1
    m, dr = c11.m0(n)
    v = c11.nu(n); d = n - v
    mu = 2 * m + d - 1
    print('  k=%2d n=%d  m0=%d (pred %d, %s)  d=%d (pred %d, %s)  mu=%d (pred %d, %s)  drops=%d'
          % (k, n, m, 10 * 11 ** k, m == 10 * 11 ** k, d, 2 * 11 ** k + 1, d == 2 * 11 ** k + 1,
             mu, 22 * 11 ** k, mu == 22 * 11 ** k, dr))

print('p=13  n = 2*13^k+1 : predicted m0 = 12*13^k, d = 13^k+1, mu = 25*13^k')
for k in (11, 12):
    n = 2 * 13 ** k + 1
    m, dr = c13.m0(n)
    v = c13.nu(n); d = n - v
    mu = 2 * m + d - 1
    print('  k=%2d n=%d  m0=%d (pred %d, %s)  d=%d (pred %d, %s)  mu=%d (pred %d, %s)  drops=%d'
          % (k, n, m, 12 * 13 ** k, m == 12 * 13 ** k, d, 13 ** k + 1, d == 13 ** k + 1,
             mu, 25 * 13 ** k, mu == 25 * 13 ** k, dr))
