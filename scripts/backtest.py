#!/usr/bin/env python3
"""Catalogue every Corrientes flood wave since 1901 and backtest the fit.

The published page fits on 1990-2026 because that is what the first pull
returned. The a5 API actually serves Corrientes from 1901 and Rosario from
1900, so this refits on the full record and cross-validates.
"""
import datetime as dt, json, statistics as st

P = '/tmp/claude-1000/-home-zeta-parana/5456bbce-43b9-4fdc-97d9-fd46dd19c387/scratchpad/dl/long_%d.json'

def load(i):
    d = json.load(open(P % i))
    return {dt.date.fromisoformat(o['timestart'][:10]): float(o['valor'])
            for o in d['observaciones'] if o.get('valor') is not None}

C, R = load(19), load(34)
start, end = min(min(C), min(R)), max(max(C), max(R))
days = [start + dt.timedelta(d) for d in range((end - start).days + 1)]
c = [C.get(d) for d in days]
r = [R.get(d) for d in days]
print(f"record: {start} to {end}   Corrientes {len(C)} obs, Rosario {len(R)} obs")
cov = sum(1 for i in range(len(days)) if c[i] is not None and r[i] is not None)
print(f"days with both gauges: {cov} of {len(days)} ({cov/len(days)*100:.0f}%)\n")

def ols(X, y):
    k = len(X[0])
    A = [[sum(X[q][i]*X[q][j] for q in range(len(X))) for j in range(k)]
         + [sum(X[q][i]*y[q] for q in range(len(X)))] for i in range(k)]
    for i in range(k):
        pv = max(range(i, k), key=lambda t: abs(A[t][i])); A[i], A[pv] = A[pv], A[i]
        for j in range(i+1, k):
            f = A[j][i]/A[i][i]
            for m in range(i, k+1): A[j][m] -= f*A[i][m]
    b = [0.0]*k
    for i in range(k-1, -1, -1):
        b[i] = (A[i][k] - sum(A[i][j]*b[j] for j in range(i+1, k)))/A[i][i]
    return b

# ---- event catalogue
pk = []
for i in range(60, len(c) - 60):
    if c[i] is None or c[i] < 6.0: continue
    win = [x for x in c[i-12:i+13] if x is not None]
    if len(win) < 18 or c[i] != max(win): continue
    if pk and i - pk[-1] < 45:
        if c[i] > c[pk[-1]]: pk[-1] = i
        continue
    pk.append(i)

ev = []
for i in pk:
    seg = [(j, r[j]) for j in range(i-5, min(i+51, len(r))) if r[j] is not None]
    if len(seg) < 30 or r[i-5] is None: continue
    j, rpk = max(seg, key=lambda t: t[1])
    if not (i-2 < j < min(i+50, len(r)-1) - 3): continue
    ev.append({'i': i, 'd': days[i], 'c': c[i], 'base': r[i-5],
               'r': rpk, 'rd': days[j], 'lag': (days[j]-days[i]).days})
print(f"matched events, Corrientes peak >= 6.00 m: {len(ev)}")
by_era = {}
for e in ev: by_era.setdefault(e['d'].year//20*20, []).append(e)
print("  per era:", {k: len(v) for k, v in sorted(by_era.items())}, "\n")

# ---- full-record fit vs the published 1990+ fit
def fit_on(sample):
    B = ols([[1.0, e['c'], e['base']] for e in sample], [e['r'] for e in sample])
    res = [e['r'] - (B[0]+B[1]*e['c']+B[2]*e['base']) for e in sample]
    return B, (sum(x*x for x in res)/(len(sample)-3))**0.5

B_all, sd_all = fit_on(ev)
mod = [e for e in ev if e['d'].year >= 1990]
B_90, sd_90 = fit_on(mod)
print(f"fit 1901-2026  n={len(ev):3}  R = {B_all[0]:.3f} + {B_all[1]:.4f}*C + {B_all[2]:.4f}*base   sd {sd_all:.3f}")
print(f"fit 1990-2026  n={len(mod):3}  R = {B_90[0]:.3f} + {B_90[1]:.4f}*C + {B_90[2]:.4f}*base   sd {sd_90:.3f}")
for B, sd, tag in ((B_all, sd_all, '1901+'), (B_90, sd_90, '1990+')):
    print(f"   {tag}: C 7.31, base 3.35 -> {B[0]+B[1]*7.31+B[2]*3.35:.2f} +/- {sd:.2f} m")
lags = sorted(e['lag'] for e in ev)
n = len(lags)
print(f"lag 1901+: median {st.median(lags):.0f}  p25 {lags[n//4]}  p75 {lags[3*n//4]}\n")

# ---- leave-one-out backtest on the full record
err = []
for k, e in enumerate(ev):
    rest = ev[:k] + ev[k+1:]
    B, _ = fit_on(rest)
    p = B[0] + B[1]*e['c'] + B[2]*e['base']
    err.append((e, p, e['r'] - p))
ae = sorted(abs(x[2]) for x in err)
print(f"leave-one-out, n={len(err)}:")
print(f"  MAE {sum(abs(x[2]) for x in err)/len(err):.2f} m   median |err| {st.median(ae):.2f} m"
      f"   p90 {ae[int(.9*len(ae))]:.2f} m   worst {ae[-1]:.2f} m")
big = [x for x in err if x[0]['c'] >= 7.2]
print(f"  on the {len(big)} events with C_peak >= 7.2 m (this wave's class):"
      f" MAE {sum(abs(x[2]) for x in big)/len(big):.2f} m,"
      f" bias {sum(x[2] for x in big)/len(big):+.2f} m")

# ---- top 3 by Corrientes peak
print("\ntop 10 by Corrientes peak:")
print("  C date       C_pk  R_base  R_crest  lag   LOO pred   err")
for e, p, d in sorted(err, key=lambda x: -x[0]['c'])[:10]:
    print(f"  {e['d']}  {e['c']:5.2f}  {e['base']:5.2f}   {e['r']:5.2f}  {e['lag']:3d}d"
          f"    {p:5.2f}   {d:+5.2f}")
json.dump([{k: (v.isoformat() if isinstance(v, dt.date) else v) for k, v in e.items()}
           for e in ev], open('data/events_1901.json', 'w'), indent=1)
