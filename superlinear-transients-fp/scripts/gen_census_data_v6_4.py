
#!/usr/bin/env python3
"""Regenerates the computer-assisted data of the paper and checks Tables 3-5."""
import bisect, sys
from fractions import Fraction
from collections import Counter

def digits(n, p):
    ds = []
    while n:
        n, r = divmod(n, p); ds.append(r)
    return ds[::-1]

class PF:
    def __init__(self, p):
        self.p = p
        self.pow2 = [pow(2, r, p) for r in range(p)]
        C = [[0] * p for _ in range(p)]
        for a in range(p):
            C[a][0] = 1
            for i in range(1, a + 1):
                C[a][i] = (C[a-1][i-1] + C[a-1][i]) % p
        self.C = C
        th = [[sum(C[r][i] for i in range(min(j, r) + 1)) % p for r in range(p)] for j in range(p)]
        et = [[(self.pow2[r] - th[j][r]) % p for r in range(p)] for j in range(p)]
        self.th, self.et = th, et
        self.thm = [[0] * p if b == 0 else th[b-1] for b in range(p)]
        self.etm = [self.pow2[:] if b == 0 else et[b-1] for b in range(p)]
        self.inv = [0] + [pow(x, p - 2, p) for x in range(1, p)]
        self._good = {}
    def S(self, D, m):
        if D < 0: return 0
        p = self.p; pairs = []
        while D:
            D, b = divmod(D, p); m, r = divmod(m, p); pairs.append((b, r))
        s1, s0 = 1, 0
        for b, r in reversed(pairs):
            s1, s0 = (self.th[b][r]*s1 + self.et[b][r]*s0) % p, (self.thm[b][r]*s1 + self.etm[b][r]*s0) % p
        return s1
    def binom(self, m, j):
        p = self.p; res = 1
        while m or j:
            mi, ji = m % p, j % p
            if ji > mi: return 0
            res = res * self.C[mi][ji] % p; m //= p; j //= p
        return res
    def pstep(self, Om, b):
        p = self.p; out = set()
        for w in Om:
            for r in range(p):
                u = (self.th[b][r] + self.et[b][r]*w) % p
                if u == 0: return None
                out.add((self.thm[b][r] + self.etm[b][r]*w) * self.inv[u] % p)
        return frozenset(out)
    def Omega(self, q):
        Om = frozenset([0])
        for b in digits(q, self.p):
            Om = self.pstep(Om, b)
            if Om is None: return None
        return Om
    def good(self, q):
        if q < 0: return False
        if q not in self._good: self._good[q] = self.Omega(q) is not None
        return self._good[q]
    def nu(self, n):
        D = n - 1
        while not self.good(D): D -= 1
        return D + 1
    def rigid(self, W):
        Om = self.Omega(W)
        return W >= 1 and Om is not None and all(self.pstep(Om, c) is None for c in range(self.p - 1))

FF = {}
def F_(p):
    if p not in FF: FF[p] = PF(p)
    return FF[p]

def m0_direct(F, n):
    """Cascade with block acceleration (Lemma 2.5)."""
    p = F.p; tgt = F.nu(n) - 1; D = n - 1; m = 0
    while D > tgt:
        j, x = 0, D
        while x % p == p - 1: x //= p; j += 1
        P = p ** j
        m = (m // P + 1) * P
        while F.S(D, m) == 0: D -= 1
    return m

class Ren:
    def __init__(self, F, W):
        self.F, self.p, self.W = F, F.p, W
        p = self.p
        assert W >= 1 and F.good(W)
        self.ell = len(digits(W, p)); self.M = p ** self.ell
        inv2 = [F.inv[F.pow2[r]] for r in range(p)]
        self.sig = [[F.th[c][r] * inv2[r] % p for r in range(p)] for c in range(p)]
        self.TW, self.Lam = [], {}
        for t in range(self.M):
            bt = F.binom(t, W)
            if bt:
                self.TW.append(t); self.Lam[t] = (-F.S(W-1, t)) * F.inv[bt] % p
        self.Z, self.Zs = [], []
        for c in range(p):
            Zc = [(t, r) for t in self.TW for r in range(p) if self.sig[c][r] == self.Lam[t]]
            self.Z.append(Zc); self.Zs.append(set(Zc))
        self._b = {}; self._nuW = F.nu(W)
    def Next(self, c, E, t, r):
        Zc = self.Z[c]; i = bisect.bisect_right(Zc, (t, r))
        return (E,) + Zc[i] if i < len(Zc) else (E + 1,) + Zc[0]
    def Ph(self, c, E, t, r):
        p = self.p
        while True:
            if not self.Z[c]: return ('STUCK', c, E, t, r)
            E, t, r = self.Next(c, E, t, r)
            cs = [c2 for c2 in range(c) if (t, r) not in self.Zs[c2]]
            if cs: c = max(cs); continue
            if (t, p - 1) in self.Zs[p - 2]: return ('DEEP',)
            return ('T', E, t, r)
    def b(self, th):
        if th not in self._b:
            F = self.F; L = self.W - 1; t = th; tgt = self._nuW - 1
            assert F.S(L, th) != 0
            while L > tgt:
                t += 1
                while F.S(L, t) == 0: L -= 1
            self._b[th] = t
        return self._b[th]
    def family(self, d):
        p = self.p
        if not self.Z[d]: return {'kind': 'goodstart'}
        top = self.Ph(d, 0, 0, -1)
        if top[0] == 'STUCK': return {'kind': 'linear', 'i': 0, 'stuck': top[1:], 'thetas': [], 'E': []}
        if top[0] == 'DEEP': return {'kind': 'deep'}
        E, th = top[1], top[2]
        thetas, Es, deltas, idx = [th], [E], [], {th: 0}
        while True:
            res = self.Ph(p - 2, E, th, p - 1)
            if res[0] == 'STUCK':
                return {'kind': 'linear', 'i': len(thetas), 'stuck': res[1:], 'thetas': thetas, 'E': Es}
            if res[0] == 'DEEP': return {'kind': 'deep'}
            deltas.append(res[1] - E); E, nxt = res[1], res[2]
            if nxt in idx:
                st = idx[nxt]; per = len(thetas) - st
                return {'kind': 'adm', 'thetas': thetas, 'E': Es, 'deltas': deltas,
                        'start': st, 'per': per, 'rho': Fraction(sum(deltas[st:]), per)}
            idx[nxt] = len(thetas); thetas.append(nxt); Es.append(E); th = nxt
    def thE(self, fam, k):
        th, de, st, per = fam['thetas'], fam['deltas'], fam['start'], fam['per']
        at = lambda seq, i: seq[i-1] if i - 1 < len(seq) else seq[st + (i - 1 - st) % per]
        E = fam['E'][0]
        for i in range(1, k): E += at(de, i)
        return at(th, k), E
    def tm(self, m, E, t, r):
        p = self.p
        return E * p**(m + self.ell) + t * p**m + p**m - p + r
    def m0_pred(self, fam, k):
        p, M = self.p, self.M
        if fam['kind'] == 'goodstart': return 0
        if fam['kind'] == 'adm':
            th, E = self.thE(fam, k); return p**k * (M * E + self.b(th))
        i = fam['i']; c, E, t, r = fam['stuck']
        if k <= i: return p**k * (M * fam['E'][k-1] + self.b(fam['thetas'][k-1]))
        return p**i * self.tm(k - i, E, t, r)
    def n(self, d, k):
        return (self.W + 1) * self.p**k - self.p + d + 1

def fr(x):
    x = Fraction(x)
    return str(x.numerator) if x.denominator == 1 else "%d/%d" % (x.numerator, x.denominator)

warn = []

# ---- Table 3 -------------------------------------------------------------
tab3 = {3:(2,0),5:(4,1),7:(4,0),11:(8,2),13:(14,6),17:(28,15),19:(18,8),23:(29,15),
        29:(45,19),31:(33,16),37:(64,32),41:(87,41)}
smallest = {5:12,11:822,13:2,17:2,19:2,23:158,29:6,31:6,37:10,41:12}
for p, (ns, nr) in tab3.items():
    F = F_(p); start = frozenset([0]); seen = [start]; S_ = {start}; i = 0; rig = 0
    while i < len(seen):
        Om = seen[i]; i += 1
        nx = [F.pstep(Om, b) for b in range(p)]
        assert all(nx[b] is None for b in range(1, p, 2))
        if Om != start and all(nx[c] is None for c in range(p - 1)): rig += 1
        for t in nx:
            if t is not None and t not in S_: S_.add(t); seen.append(t)
    if (len(seen), rig) != (ns, nr): warn.append(("Table 3", p, len(seen), rig))
    if p in smallest:
        W = 1
        while not F.rigid(W): W += 1
        if W != smallest[p]: warn.append(("smallest rigid", p, W))
print("Table 3 done")

# ---- Table 4 (superlinear families) ------------------------------------
tab4 = [(5,12,'1/2',2,'25/13'),(5,72,'1/2',2,'125/73'),(5,512,'1/10',10,'125/513'),
        (11,822,'1/2',2,'1331/823'),(11,1064,'1',1,'2662/1065'),(13,2,'3',1,'26'),
        (13,108,'1',1,'338/109'),(13,134,'2',1,'676/135'),(13,158,'3',1,'338/53'),
        (17,2,'3',1,'34'),(17,4,'3',1,'102/5'),(17,8,'3',1,'34/3'),(19,2,'3',1,'38'),
        (19,12,'4',1,'152/13'),(23,158,'2',1,'2116/159'),(29,6,'4',1,'232/7'),
        (29,14,'5',1,'58/3'),(31,6,'4',1,'248/7'),(31,10,'5',1,'310/11'),
        (37,10,'4',1,'296/11'),(37,26,'8',1,'592/27'),(41,12,'4',1,'328/13')]
rows = []
for p, W, r0, per0, cmax in tab4:
    R = Ren(F_(p), W)
    if not F_(p).rigid(W): warn.append(("Table 4 not rigid", p, W))
    fams = [R.family(d) for d in range(p - 1)]
    if not all(f['kind'] == 'adm' and f['rho'] > 0 for f in fams):
        warn.append(("Table 4 not admissible with rho>0", p, W)); continue
    if (fr(fams[0]['rho']), fams[0]['per']) != (r0, per0): warn.append(("Table 4 rho/period", p, W))
    cs = [2 * f['rho'] * R.M / Fraction(W + 1) for f in fams]; mx = max(cs)
    if fr(mx) != cmax: warn.append(("Table 4 max c", p, W, fr(mx)))
    rhos = sorted(set(f['rho'] for f in fams))
    if len(rhos) != 1:
        warn.append(("Table 4 rho depends on d", p, W, rhos))
    rows.append(r"%d & %d & $%s$ & %d & $%s$\\" % (
        p, W, fr(fams[0]['rho']), fams[0]['per'], fr(mx)))
with open("tab_superlinear.tex", "w") as fh:
    fh.write("\\begin{table}[htbp]\n\\centering\\small\n"
             "\\caption{Superlinear families: $\\rho$ and the period $\\varpi$ of $(\\theta_i)$ for $d=0$, "
             "and the maximum over $d$ of $c=2\\rho(d)M/(W+1)$. In every row $\\rho(d)$ is independent of $d$, "
             "so the maximum is attained for all $d$.}\\label{tab:superlinear}\n"
             "\\begin{tabular}{rrccc}\n\\toprule\n"
             "$p$ & $W$ & $\\rho(0)$ & $\\varpi(0)$ & $\\max_dc$\\\\\n"
             "\\midrule\n" + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n")
print("Table 4 done")

# ---- m_0(99573) -----------------------------------------------------------
R = Ren(F_(11), 822); f = R.family(0)
assert f['kind'] == 'adm' and R.n(0, 2) == 99573
th2, E2 = R.thE(f, 2); b2 = R.b(th2); pred = R.m0_pred(f, 2)
if pred != 1273767: warn.append(("m0(99573) formula", pred))
direct = m0_direct(F_(11), 99573)
if direct != pred: warn.append(("m0(99573) direct", direct))
with open("decomp_99573.tex", "w") as fh:
    fh.write("For instance, for $p=11$, $W=822$, $d=0$, $k=2$ one finds $(E_1,\\theta_1)=(%d,%d)$, "
             "$\\theta_2=\\Phi(\\theta_1)=%d$, $E(2)=%d$ and $b(\\theta_2)=%d$, so Theorem~\\ref{thm:renorm} gives "
             "$m_0(99573)=11^2(1331\\cdot%d+%d)=%d$, in agreement with the directly computed cascade.\n"
             % (f['E'][0], f['thetas'][0], th2, E2, b2, E2, b2, pred))
print("99573 done")

# ---- Table 5 and the triple count ---------------------------------------
tab5 = {3:(729,126,63,63,0,0,0), 5:(3125,352,93,155,104,0,26), 7:(2401,282,62,220,0,0,0),
        11:(1331,310,43,247,20,0,2), 13:(2197,744,90,354,300,0,25)}
triples = 0
for p, row in tab5.items():
    F = F_(p); PL = row[0]; cnt = Counter(); rig = 0
    for W in range(1, PL):
        if not F.good(W): continue
        R = Ren(F, W); isr = F.rigid(W); rig += isr
        for d in range(p - 1):
            fam = R.family(d); cnt[fam['kind']] += 1
            if isr and not (fam['kind'] == 'adm' and fam['rho'] > 0): warn.append(("rigid not adm", p, W, d))
            if not isr and fam['kind'] == 'adm': warn.append(("nonrigid adm", p, W, d))
            if fam['kind'] == 'deep': continue
            k = 1
            while R.n(d, k) <= 12000:
                n = R.n(d, k)
                if m0_direct(F, n) != R.m0_pred(fam, k): warn.append(("formula", p, W, d, k))
                triples += 1; k += 1
    got = (PL, sum(cnt.values()), cnt['goodstart'], cnt['linear'], cnt['adm'], cnt['deep'], rig)
    if got != row: warn.append(("Table 5", p, got))
    print("Table 5, p =", p, "done")
with open("count_triples.tex", "w") as fh:
    fh.write("$%d$ " % triples)

# ---- Appendix B: smallest rigid level of every prime except p = 5 -------------
import csv
SMALLEST = [(11, 822), (13, 2), (17, 2), (19, 2), (23, 158), (29, 6), (31, 6), (37, 10), (41, 12)]
out = []
for p, W in SMALLEST:
    F = F_(p); R = Ren(F, W); fams = [R.family(d) for d in range(p - 1)]
    if not all(fm['kind'] == 'adm' and fm['rho'] > 0 for fm in fams):
        warn.append(("appendix: not admissible with rho>0", p, W)); continue
    one = (R.ell == 1)
    for d in range(p - 1):
        for k in (1, 2):
            if not one and k == 2 and d != 0: continue
            if m0_direct(F, R.n(d, k)) != R.m0_pred(fams[d], k): warn.append(("appendix", p, W, d, k))
    s = [r"\subsection*{$p=%d$, $W=%d=(%s)_{%d}$}" % (p, W, ",".join(map(str, digits(W, p))), p),
         r"Here $\ell=%d$, $M=%d$, $\nu(W)=%d$ (so $\nu(n_k)=%d\cdot%d^k$) and $|\mathcal T_W|=%d$."
         % (R.ell, R.M, R._nuW, R._nuW, p, len(R.TW))]
    if one:
        s.append("Label function: " + ", ".join(r"$\Lambda_{%d}(%d)=%d$" % (W, t, R.Lam[t]) for t in R.TW) + ".")
        s.append(r"\begin{itemize}[leftmargin=*]")
        for c in range(p):
            s.append((r"\item $Z_{%d}=\{" % c + ", ".join("(%d,%d)" % z for z in R.Z[c]) + r"\}$.")
                     if R.Z[c] else r"\item $Z_{%d}=\emptyset$." % c)
        s.append(r"\end{itemize}")
    else:
        fn = "finite_data_p%d_W%d.csv" % (p, W)
        with open(fn, "w", newline="") as fh:
            wr = csv.writer(fh); wr.writerow(["kind", "c", "t", "r", "value"])
            for t in R.TW: wr.writerow(["Lambda", "", t, "", R.Lam[t]])
            for c in range(p):
                for (t, r) in R.Z[c]: wr.writerow(["Z", c, t, r, ""])
        s.append("Sizes: " + ", ".join(r"$|Z_{%d}|=%d$" % (c, len(R.Z[c])) for c in range(p))
                 + r"; full lists in \texttt{%s}." % fn.replace("_", r"\_"))
    grp = {}
    for d, fm in enumerate(fams):
        grp.setdefault((fm['E'][0], fm['thetas'][0]), []).append(d)
    s.append("Top phase: " + "; ".join(r"$(E_1,\theta_1)=(%d,%d)$ for $d\in\{%s\}$"
             % (a, b, ",".join(map(str, v))) for (a, b), v in sorted(grp.items())) + ".")
    items = []
    for th in sorted(set(t for fm in fams for t in fm['thetas'])):
        res = R.Ph(p - 2, 0, th, p - 1)
        items.append(r"$\Phi(%d)=%d$, $\delta(%d)=%d$, $b(%d)=%d$" % (th, res[2], th, res[1], th, R.b(th)))
    chk = (r"for $k\in\{1,2\}$ and all $d$" if one else r"for $k=1$ and all $d$, and for $k=2$, $d=0$")
    s.append("Renormalization and terminal data: " + "; ".join(items) + ". "
             + r"Hence every $(%d,d)$ is admissible, with $\rho\in\{%s\}$ and $\varpi\in\{%s\}$. "
             % (W, ", ".join(fr(x) for x in sorted(set(fm['rho'] for fm in fams))),
                ", ".join(map(str, sorted(set(fm['per'] for fm in fams)))))
             + r"Theorem~\ref{thm:renorm} agrees with the direct cascade " + chk + ".")
    out.append("\n".join(s))
    print("appendix", p, W, "done")
with open("appendix_rigid.tex", "w") as fh:
    fh.write("\n\n".join(out) + "\n")

print("triples:", triples)
if warn:
    print("MISMATCHES (the paper must be corrected accordingly):")
    for w in warn: print("  ", w)
    sys.exit(1)
print("all checks passed")
