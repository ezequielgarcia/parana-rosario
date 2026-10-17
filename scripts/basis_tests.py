#!/usr/bin/env python3
"""Does a Fourier basis help, and what does a routing kernel say instead?

Three questions, all answered on the 1901-2026 daily record:
  1. how much of Rosario's variance is periodic (harmonic regression)
  2. what the de-seasonalised spectrum looks like at event scale
  3. whether a lag-and-cascade transfer function from Corrientes beats the
     peak-to-peak regression the page currently uses

No numpy here, so the cascade is run as n successive first-order recursive
filters, which is exact for integer n and costs O(N) per evaluation.
"""
import datetime as dt, json, math, statistics as st

P = '/tmp/claude-1000/-home-zeta-parana/5456bbce-43b9-4fdc-97d9-fd46dd19c387/scratchpad/dl/long_%d.json'

def load(i):
    d = json.load(open(P % i))
    return {dt.date.fromisoformat(o['timestart'][:10]): float(o['valor'])
            for o in d['observaciones'] if o.get('valor') is not None}

C, R = load(19), load(34)
start = max(min(C), min(R)); end = min(max(C), max(R))
days = [start + dt.timedelta(d) for d in range((end - start).days + 1)]
N = len(days)
print(f"grid: {start} to {end}, {N} days\n")

# fill short gaps by linear interpolation so the filters see a continuous series
def series(src):
    v = [src.get(d) for d in days]
    first = next(i for i, x in enumerate(v) if x is not None)
    last = len(v) - 1 - next(i for i, x in enumerate(reversed(v)) if x is not None)
    i = first
    while i <= last:
        if v[i] is None:
            j = i
            while v[j] is None: j += 1
            a, b = v[i-1], v[j]
            for m in range(i, j): v[m] = a + (b-a)*(m-i+1)/(j-i+1)
            i = j
        i += 1
    return v, first, last

c, cf, cl = series(C)
r, rf, rl = series(R)
lo, hi = max(cf, rf), min(cl, rl)
print(f"usable span: {days[lo]} to {days[hi]}, {hi-lo+1} days\n")

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

def r2(y, pred):
    m = st.mean(y)
    return 1 - sum((a-b)**2 for a, b in zip(y, pred))/sum((a-m)**2 for a in y)

# ---------------------------------------------------------------- 1. harmonics
print("=" * 64)
print("1. HARMONIC REGRESSION on Rosario stage, annual fundamental")
yr = 365.2425
idx = list(range(lo, hi+1))
y = [r[i] for i in idx]
doy = [(days[i] - dt.date(days[i].year, 1, 1)).days for i in idx]

for K in (1, 2, 3, 4, 6):
    X = []
    for q, i in enumerate(idx):
        row = [1.0]
        for k in range(1, K+1):
            a = 2*math.pi*k*(i-lo)/yr
            row += [math.cos(a), math.sin(a)]
        X.append(row)
    B = ols(X, y)
    pred = [sum(B[j]*X[q][j] for j in range(len(B))) for q in range(len(X))]
    amps = [(B[1+2*(k-1)]**2 + B[2+2*(k-1)]**2)**0.5 for k in range(1, K+1)]
    res = [a-b for a, b in zip(y, pred)]
    print(f"  K={K}: R2 {r2(y,pred):.4f}   residual sd {st.pstdev(res):.3f} m"
          f"   amplitudes " + " ".join(f"{a:.3f}" for a in amps))
    if K == 2:
        HARM, HPRED = B, pred
        ph = []
        for k in (1, 2):
            a, b = B[1+2*(k-1)], B[2+2*(k-1)]
            d0 = (-math.atan2(b, a)/(2*math.pi/ (yr/k))) % (yr/k)
            ph.append(days[lo] + dt.timedelta(d0))
        print(f"        annual peak near {ph[0].strftime('%d %b')}, "
              f"semiannual near {ph[1].strftime('%d %b')}")

ANOM = [a-b for a, b in zip(y, HPRED)]
print(f"\n  seasonal cycle explains {100*r2(y,HPRED):.1f}% of variance;"
      f" the anomaly keeps sd {st.pstdev(ANOM):.2f} m against {st.pstdev(y):.2f} m total")

# ------------------------------------------------- 2. spectrum of the anomaly
print("\n" + "=" * 64)
print("2. PERIODOGRAM of the de-seasonalised anomaly")
print("   power = fraction of anomaly variance taken by one sine pair")
var = sum(a*a for a in ANOM)
n = len(ANOM)
for per in (3650, 1460, 730, 365, 182, 120, 90, 60, 45, 30, 20, 14, 7):
    w = 2*math.pi/per
    cs = sum(a*math.cos(w*t) for t, a in enumerate(ANOM))
    sn = sum(a*math.sin(w*t) for t, a in enumerate(ANOM))
    print(f"   period {per:5d} d   power {2*(cs*cs+sn*sn)/n/var:7.4f}")

# ------------------------------------- 3. lag-and-cascade transfer function
print("\n" + "=" * 64)
print("3. LAG AND CASCADE, Corrientes stage -> Rosario stage")

def cascade(x, n, k):
    """n linear reservoirs in series, storage constant k days, dt = 1 day."""
    a = math.exp(-1.0/k); b = 1.0 - a
    out = x
    for _ in range(n):
        y_, prev = [], out[0]
        for v in out:
            prev = a*prev + b*v
            y_.append(prev)
        out = y_
    return out

CIN = [c[i] for i in idx]
WARM = 400                      # drop the filter spin-up
YS = [y[t] for t in range(WARM, len(y))]
MY = st.mean(YS); SST = sum((a-MY)**2 for a in YS)

def score(n, k, tau):
    """Best gain and intercept for this kernel, and the R2 they reach."""
    f = CASC.get((n, k)) or CASC.setdefault((n, k), cascade(CIN, n, k))
    sx = sxx = sy = sxy = 0.0
    m = len(f)
    for t in range(WARM, m):
        v = f[t-tau]; u = y[t]
        sx += v; sxx += v*v; sy += u; sxy += v*u
    m2 = m - WARM
    den = m2*sxx - sx*sx
    b1 = (m2*sxy - sx*sy)/den
    b0 = (sy - b1*sx)/m2
    sse = 0.0
    for t in range(WARM, m):
        e = y[t] - (b0 + b1*f[t-tau]); sse += e*e
    return 1 - sse/SST, b0, b1, (sse/m2)**0.5

CASC = {}
best = None
for n in range(1, 7):
    for k10 in range(20, 420, 20):
        k = k10/10
        for tau in (0, 5, 10, 15, 20):
            s, b0, b1, sd = score(n, k, tau)
            if best is None or s > best[0]: best = (s, n, k, tau, [b0, b1], sd)
    CASC.clear()
print("  coarse best:", f"n={best[1]} k={best[2]:.1f} tau={best[3]} R2={best[0]:.4f}")

# refine around it
bn, bk, bt = best[1], best[2], best[3]
CASC = {}
for n in range(max(1, bn-1), bn+2):
    for k10 in range(max(5, int(bk*10)-18), int(bk*10)+19, 3):
        k = k10/10
        for tau in range(max(0, bt-5), bt+6):
            s, b0, b1, sd = score(n, k, tau)
            if s > best[0]: best = (s, n, k, tau, [b0, b1], sd)
    CASC.clear()

s, n, k, tau, B, sd = best
print(f"  best: n={n} reservoirs, k={k:.1f} d, pure lag {tau} d")
print(f"        Rosario = {B[0]:+.3f} {B[1]:+.4f} * cascade_n,k(Corrientes)[t-{tau}]")
print(f"        R2 {s:.4f}, residual sd {sd:.3f} m on {hi-lo+1-WARM} days")
print(f"        gain {B[1]:.3f}  (the page's attenuation figure is 0.36)")
tp = tau + n*k                     # time to peak of a gamma kernel plus the lag
print(f"        kernel time to peak {tp:.0f} d, centroid lag {tau + n*k:.0f} d")

KM = 792.0                        # Corrientes km 1208 to Rosario km 416
print(f"\n  celerity implied: {KM/tp:.0f} km/d = {KM*1000/(tp*86400):.2f} m/s")
print(f"  mean water velocity in this reach is about 1 m/s, so a pure")
print(f"  kinematic wave (c = 1.5 v) would arrive in {KM*1000/(1.5*86400):.0f} days")

# how well does it do on the crest of the 76 catalogued waves?
EV = json.load(open('data/events_1901.json'))
f = cascade(CIN, n, k)
pred_series = [B[0]+B[1]*f[t-tau] for t in range(len(f))]
errs = []
for e in EV:
    d = dt.date.fromisoformat(e['rd'])
    if not (days[lo] <= d <= days[hi]): continue
    t = (d - days[lo]).days
    if t < WARM: continue
    w = pred_series[max(0, t-20):t+21]
    errs.append(e['r'] - max(w))
if errs:
    print(f"\n  crest error on {len(errs)} catalogued waves: "
          f"MAE {sum(abs(x) for x in errs)/len(errs):.2f} m, "
          f"bias {st.mean(errs):+.2f} m, sd {st.pstdev(errs):.2f} m")
    print(f"  (the peak-to-peak regression gets MAE 0.19 m on the same waves)")
