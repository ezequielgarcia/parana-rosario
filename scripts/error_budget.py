#!/usr/bin/env python3
"""Full error budget for the Rosario crest estimate.

The page used to quote +/- 0.24 m, which is only the regression residual sd.
That understates the real uncertainty: it ignores parameter uncertainty, the
fact that the Corrientes peak is itself a forecast, and the spread between the
three independent estimators. This script adds those terms, and checks the
parametric interval against the empirical leave-one-out error distribution.
"""
import datetime as dt, json, math, statistics as st

EV = json.load(open('data/events_1901.json'))
for e in EV: e['d'] = dt.date.fromisoformat(e['d'])

# ---- linear algebra, by hand
def solve(A, b):
    n = len(A); M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for i in range(n):
        p = max(range(i, n), key=lambda t: abs(M[t][i])); M[i], M[p] = M[p], M[i]
        for j in range(i+1, n):
            f = M[j][i]/M[i][i]
            for m in range(i, n+1): M[j][m] -= f*M[i][m]
    x = [0.0]*n
    for i in range(n-1, -1, -1):
        x[i] = (M[i][n] - sum(M[i][j]*x[j] for j in range(i+1, n)))/M[i][i]
    return x

def inverse(A):
    n = len(A)
    return [[c[i] for c in (solve(A, [1.0 if k == j else 0.0 for k in range(n)])
                            for j in range(n))] for i in range(n)]

def design(sample):
    return [[1.0, e['c'], e['base']] for e in sample]

def xtx(X):
    k = len(X[0])
    return [[sum(r[i]*r[j] for r in X) for j in range(k)] for i in range(k)]

def fit(sample):
    X, y = design(sample), [e['r'] for e in sample]
    k = len(X[0])
    B = solve(xtx(X), [sum(X[q][i]*y[q] for q in range(len(X))) for i in range(k)])
    res = [y[q] - sum(B[i]*X[q][i] for i in range(k)) for q in range(len(X))]
    sd = (sum(r*r for r in res)/(len(X)-k))**0.5
    return B, sd, res

B, SD, RES = fit(EV)
n, k = len(EV), 3
XTXI = inverse(xtx(design(EV)))

# 22 Oct: the governing Corrientes peak is INA's forecast second crest on
# 3 Nov, not the 7.38 m observed on 22 Oct, so the input term is a forecast
# again and so is the base, Rosario five days ahead of that peak.
C_MID, C_LO, C_HI = 7.58, 7.00, 7.80   # INA's own band on the Corrientes peak
BASE, BASE_LO, BASE_HI = 4.12, 3.71, 4.44
x0 = [1.0, C_MID, BASE]

pred = sum(B[i]*x0[i] for i in range(k))
# leverage term: how far this x0 sits from the centre of the 76-wave cloud
lev = sum(x0[i]*XTXI[i][j]*x0[j] for i in range(k) for j in range(k))
se_param = SD*lev**0.5                      # uncertainty in the fitted mean
se_resid = SD                               # wave-to-wave scatter about the fit
se_fit = SD*(1.0 + lev)**0.5                # textbook prediction SE

print(f"fit on {n} waves, {n-k} df")
print(f"  R_crest = {B[0]:+.3f} {B[1]:+.4f}*C_peak {B[2]:+.4f}*R_base")
print(f"  residual sd        {SD:.3f} m")
print(f"  leverage at x0     {lev:.4f}  (se of the fitted mean {se_param:.3f} m)")
print(f"  prediction se      {se_fit:.3f} m   -> central estimate {pred:.2f} m\n")

# ---- term 2: the Corrientes peak is a forecast, not an observation
sd_c = (C_HI - C_LO)/(2*1.645)   # INA publishes a band; read it as ~90%
se_input = abs(B[1])*sd_c
# ---- term 3: the base five days before the peak is now a forecast too
sd_base = (BASE_HI - BASE_LO)/(2*1.645)
se_base = abs(B[2])*sd_base
# ---- term 4: spread between the three independent estimators (model form)
ESTS = [5.22, 5.30, 5.03]
se_model = st.pstdev(ESTS)

print("error budget, one standard deviation each:")
rows = [("wave-to-wave scatter about the fit", se_resid),
        ("uncertainty in the fitted coefficients", se_param),
        (f"Corrientes peak forecast (sd {sd_c:.2f} m)", se_input),
        (f"starting base at Rosario (sd {sd_base:.2f} m)", se_base),
        ("spread between the three estimators", se_model)]
for lab, v in rows: print(f"  {lab:42s} {v:.3f} m")
tot = sum(v*v for _, v in rows)**0.5
print(f"  {'combined in quadrature':42s} {tot:.3f} m\n")

# t quantiles at 73 df, close enough to normal
T = {0.68: 1.00, 0.80: 1.293, 0.90: 1.666, 0.95: 1.993}
for p, t in T.items():
    print(f"  {int(p*100)}%: {pred:.2f} +/- {t*tot:.2f} m   "
          f"[{pred-t*tot:.2f}, {pred+t*tot:.2f}]")

# ---- empirical check: leave-one-out error quantiles
print("\nleave-one-out check (no distributional assumption):")
err = []
for i, e in enumerate(EV):
    Bi, _, _ = fit(EV[:i] + EV[i+1:])
    err.append(e['r'] - (Bi[0] + Bi[1]*e['c'] + Bi[2]*e['base']))
s = sorted(err)
def q(p): 
    i = p*(len(s)-1); lo = int(i)
    return s[lo] + (i-lo)*(s[min(lo+1, len(s)-1)] - s[lo])
print(f"  n={len(err)}  MAE {sum(abs(x) for x in err)/len(err):.2f}  "
      f"sd {st.pstdev(err):.2f}  bias {st.mean(err):+.2f}")
for p in (0.80, 0.90):
    a, b = q((1-p)/2), q(1-(1-p)/2)
    print(f"  central {int(p*100)}% of LOO errors: {a:+.2f} to {b:+.2f} m"
          f"   -> {pred+a:.2f} to {pred+b:.2f} m")
big = [(e, x) for e, x in zip(EV, err) if e['c'] >= 7.2]
bs = sorted(x for _, x in big)
print(f"  on the {len(big)} waves of this class (C >= 7.2 m): "
      f"sd {st.pstdev([x for _, x in big]):.2f}, bias {st.mean([x for _, x in big]):+.2f}, "
      f"range {bs[0]:+.2f} to {bs[-1]:+.2f}")

# ---- timing: how well can the date be pinned
print("\ntiming:")
lags = [e['lag'] for e in EV]
Bl = solve(xtx([[1.0, e['c']] for e in EV]), 
           [sum(x*e['lag'] for x, e in zip([r[i] for r in [[1.0, e['c']] for e in EV]], EV))
            for i in range(2)])
rl = [e['lag'] - (Bl[0] + Bl[1]*e['c']) for e in EV]
sdl = (sum(x*x for x in rl)/(len(EV)-2))**0.5
sl = sorted(lags)
print(f"  all 76 waves: median {st.median(lags):.0f} d, "
      f"p10 {sl[int(.1*len(sl))]} d, p90 {sl[int(.9*len(sl))]} d, range {sl[0]} to {sl[-1]}")
print(f"  lag = {Bl[0]:.1f} {Bl[1]:+.2f}*C_peak, residual sd {sdl:.1f} d"
      f"  -> at C {C_MID:.2f}: {Bl[0]+Bl[1]*C_MID:.0f} +/- {sdl:.0f} d")
big_l = sorted(e['lag'] for e in EV if e['c'] >= 7.2)
print(f"  the {len(big_l)} waves of this class: median {st.median(big_l):.0f} d, "
      f"range {big_l[0]} to {big_l[-1]} d")

# ---- threshold probabilities, from the combined sd
def phi(z): return 0.5*(1 + math.erf(z/2**0.5))
print("\nprobability of exceeding each level at the crest (normal, sd "
      f"{tot:.2f} m):")
for lab, lvl in (("4.00 m", 4.00), ("4.50 m", 4.50), ("alerta 5.00 m", 5.00),
                 ("5.30 m evacuación", 5.30), ("6.00 m", 6.00)):
    print(f"  {lab:20s} {100*(1-phi((lvl-pred)/tot)):5.1f} %")

json.dump({'pred': round(pred, 2), 'sd': round(tot, 3),
           'terms': {lab: round(v, 3) for lab, v in rows},
           'bands': {str(int(p*100)): [round(pred-t*tot, 2), round(pred+t*tot, 2)]
                     for p, t in T.items()},
           'loo': {'mae': round(sum(abs(x) for x in err)/len(err), 3),
                   'sd': round(st.pstdev(err), 3), 'bias': round(st.mean(err), 3)},
           'lag': {'median': st.median(lags), 'sd': round(sdl, 1),
                   'at_c': round(Bl[0]+Bl[1]*C_MID, 1)},
           'exceed': {str(l): round(100*(1-phi((l-pred)/tot)), 1)
                      for l in (4.0, 4.5, 5.0, 5.3, 6.0)}},
          open('data/error_budget.json', 'w'), indent=1)
