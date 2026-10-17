#!/usr/bin/env python3
"""Follow-ups: correct seasonal phase, wave broadening, and transformed spaces."""
import datetime as dt, json, math, statistics as st

P = '/tmp/claude-1000/-home-zeta-parana/5456bbce-43b9-4fdc-97d9-fd46dd19c387/scratchpad/dl/long_%d.json'
def load(i):
    d = json.load(open(P % i))
    return {dt.date.fromisoformat(o['timestart'][:10]): float(o['valor'])
            for o in d['observaciones'] if o.get('valor') is not None}
C, R = load(19), load(34)
start = max(min(C), min(R)); end = min(max(C), max(R))
days = [start + dt.timedelta(d) for d in range((end - start).days + 1)]
def series(src):
    v = [src.get(d) for d in days]
    i = 0
    while i < len(v):
        if v[i] is None and i and any(x is not None for x in v[i:]):
            j = i
            while v[j] is None: j += 1
            a, b = v[i-1], v[j]
            for m in range(i, j): v[m] = a + (b-a)*(m-i+1)/(j-i+1)
            i = j
        i += 1
    return v
c, r = series(C), series(R)
yr = 365.2425

def ols(X, y):
    k = len(X[0])
    A = [[sum(X[q][i]*X[q][j] for q in range(len(X))) for j in range(k)]
         + [sum(X[q][i]*y[q] for q in range(len(X)))] for i in range(k)]
    for i in range(k):
        p = max(range(i, k), key=lambda t: abs(A[t][i])); A[i], A[p] = A[p], A[i]
        for j in range(i+1, k):
            f = A[j][i]/A[i][i]
            for m in range(i, k+1): A[j][m] -= f*A[i][m]
    b = [0.0]*k
    for i in range(k-1, -1, -1):
        b[i] = (A[i][k] - sum(A[i][j]*b[j] for j in range(i+1, k)))/A[i][i]
    return b

# ---- 1. seasonal phase, done properly, as day of year
print("1. SEASONAL CYCLE, phase as day of year")
for name, v in (('Rosario', r), ('Corrientes', c)):
    X, y = [], []
    for i, val in enumerate(v):
        if val is None: continue
        d = days[i]
        t = (d - dt.date(d.year, 1, 1)).days / (366.0 if d.year % 4 == 0 else 365.0)
        row = [1.0]
        for k in (1, 2):
            row += [math.cos(2*math.pi*k*t), math.sin(2*math.pi*k*t)]
        X.append(row); y.append(val)
    B = ols(X, y)
    m = st.mean(y)
    pred = [sum(B[j]*row[j] for j in range(len(B))) for row in X]
    r2 = 1 - sum((a-b)**2 for a, b in zip(y, pred))/sum((a-m)**2 for a in y)
    out = [f"  {name}: mean {B[0]:.2f} m, R2 {r2:.3f}"]
    for k in (1, 2):
        a, b = B[1+2*(k-1)], B[2+2*(k-1)]
        amp = (a*a+b*b)**0.5
        frac = (math.atan2(b, a)/(2*math.pi)) % (1.0/k) if k == 1 else \
               (math.atan2(b, a)/(2*math.pi*k)) % (1.0/k)
        peak = dt.date(2025, 1, 1) + dt.timedelta(frac*365)
        out.append(f"harm{k} amp {amp:.3f} m peak {peak.strftime('%d %b')}")
    print("   ".join(out))
    hi = max(range(365), key=lambda d0: sum(
        B[1+2*(k-1)]*math.cos(2*math.pi*k*d0/365) + B[2+2*(k-1)]*math.sin(2*math.pi*k*d0/365)
        for k in (1, 2)))
    loo = min(range(365), key=lambda d0: sum(
        B[1+2*(k-1)]*math.cos(2*math.pi*k*d0/365) + B[2+2*(k-1)]*math.sin(2*math.pi*k*d0/365)
        for k in (1, 2)))
    f = lambda d0: sum(B[1+2*(k-1)]*math.cos(2*math.pi*k*d0/365) +
                       B[2+2*(k-1)]*math.sin(2*math.pi*k*d0/365) for k in (1, 2))
    print(f"      fitted cycle: high {(dt.date(2025,1,1)+dt.timedelta(hi)).strftime('%d %b')} "
          f"({B[0]+f(hi):.2f} m), low {(dt.date(2025,1,1)+dt.timedelta(loo)).strftime('%d %b')} "
          f"({B[0]+f(loo):.2f} m), range {f(hi)-f(loo):.2f} m")

# ---- 2. does the wave broaden between the two gauges?
print("\n2. CREST BREADTH, days within 0.25 m of the crest")
EV = json.load(open('data/events_1901.json'))
def breadth(v, i, win=60):
    pk = v[i]
    n = 0
    for j in range(max(0, i-win), min(len(v), i+win+1)):
        if v[j] is not None and pk - v[j] <= 0.25: n += 1
    return n
bc, br = [], []
for e in EV:
    i = (dt.date.fromisoformat(e['d']) - days[0]).days
    j = (dt.date.fromisoformat(e['rd']) - days[0]).days
    if 0 <= i < len(c) and 0 <= j < len(r):
        bc.append(breadth(c, i)); br.append(breadth(r, j))
print(f"  Corrientes: median {st.median(bc):.0f} d   Rosario: median {st.median(br):.0f} d"
      f"   ratio {st.median(br)/st.median(bc):.2f}")
print(f"  rise rates: Corrientes vs Rosario, per wave")
rc, rr = [], []
for e in EV:
    i = (dt.date.fromisoformat(e['d']) - days[0]).days
    j = (dt.date.fromisoformat(e['rd']) - days[0]).days
    if i-30 < 0 or j-30 < 0: continue
    rc.append(max(c[i-k+5]-c[i-k] for k in range(5, 26) if c[i-k] is not None))
    rr.append(max(r[j-k+5]-r[j-k] for k in range(5, 26) if r[j-k] is not None))
print(f"  fastest 5-day rise, median: Corrientes {st.median(rc):.2f} m, Rosario {st.median(rr):.2f} m")

# ---- 3. transformed spaces, single reservoir, R2 always scored on stage
print("\n3. TRANSFORMED SPACE, cascade fitted in h, log h and h^1.6")
WARM = 400
def cascade(x, n, k):
    a = math.exp(-1.0/k); b = 1.0 - a
    out = x
    for _ in range(n):
        y_, prev = [], out[0]
        for v in out:
            prev = a*prev + b*v; y_.append(prev)
        out = y_
    return out
OFF = 2.0   # keep arguments positive: Rosario dips below zero on its own datum
SST = None
for tag, fwd, inv in (('stage  h', lambda v: v, lambda v: v),
                      ('log h  ', lambda v: math.log(v+OFF), lambda v: math.exp(v)-OFF),
                      ('h^1.6  ', lambda v: (v+OFF)**1.6, lambda v: v**(1/1.6)-OFF)):
    cx = [fwd(v) for v in c]
    best = None
    for k10 in range(150, 361, 10):
        k = k10/10
        f = cascade(cx, 1, k)
        for tau in range(0, 7):
            X = [[1.0, f[t-tau]] for t in range(WARM, len(f))]
            Y = [fwd(r[t]) for t in range(WARM, len(r))]
            B = ols(X, Y)
            pred = [inv(B[0]+B[1]*row[1]) for row in X]
            obs = [r[t] for t in range(WARM, len(r))]
            m = st.mean(obs)
            if SST is None: SST = sum((a-m)**2 for a in obs)
            sse = sum((a-b)**2 for a, b in zip(obs, pred))
            s = 1 - sse/SST
            if best is None or s > best[0]: best = (s, k, tau, (sse/len(obs))**0.5)
    print(f"  {tag}: R2 on stage {best[0]:.4f}, k {best[1]:.1f} d, lag {best[2]} d, "
          f"rmse {best[3]:.3f} m")
