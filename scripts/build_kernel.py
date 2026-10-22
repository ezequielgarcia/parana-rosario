#!/usr/bin/env python3
"""Fit the single-reservoir routing filter and build the data for plates 6 and 7.

Plate 6 is the impulse response of the fitted filter. Plate 7 contrasts what a
pure travelling wave would deliver at Rosario with what the storage filter
delivers and with what the river is actually doing. Writes data/kernel.json.
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
print(f"grid {days[0]} to {days[-1]}, {len(days)} days")

def reservoir(x, k):
    a = math.exp(-1.0/k); b = 1.0 - a
    out, prev = [], x[0]
    for v in x:
        prev = a*prev + b*v
        out.append(prev)
    return out

WARM = 400
OBS = [r[t] for t in range(WARM, len(r))]
MY = st.mean(OBS); SST = sum((a-MY)**2 for a in OBS)

best = None
for k2 in range(300, 701, 5):            # k from 15.0 to 35.0 in 0.5 d steps
    k = k2/20
    f = reservoir(c, k)
    for tau in range(0, 9):
        sx = sxx = sy = sxy = 0.0
        m = len(f)
        for t in range(WARM, m):
            v = f[t-tau]; u = r[t]
            sx += v; sxx += v*v; sy += u; sxy += v*u
        n2 = m - WARM
        b1 = (n2*sxy - sx*sy)/(n2*sxx - sx*sx)
        b0 = (sy - b1*sx)/n2
        sse = sum((r[t] - (b0 + b1*f[t-tau]))**2 for t in range(WARM, m))
        s = 1 - sse/SST
        if best is None or s > best[0]:
            best = (s, k, tau, b0, b1, (sse/n2)**0.5)
R2, K, TAU, B0, B1, SD = best
print(f"fitted: k={K:.1f} d, lag={TAU} d, gain={B1:.4f}, intercept={B0:+.3f}")
print(f"        R2 {R2:.4f}, residual sd {SD:.3f} m")
print(f"        centroid travel time {TAU + K:.1f} d")

# ---- the 2026 wave: observed, INA's published track, then a tapering recession
D0 = dt.date(2026, 9, 1)
i0 = (D0 - days[0]).days
obs_last = (dt.date(2026, 10, 22) - days[0]).days
# INA's 20 October run: a dip to 7.14 m on 27 Oct, then a second and higher
# Iguazu-driven crest of 7.58 m on 3 November
INA_C = [(dt.date(2026,10,27), 7.14), (dt.date(2026,11,3), 7.58)]
CMEAN = 2.90                              # October-November seasonal low
# recession constant from the record itself: the median decay over the 30 days
# after each of the 76 catalogued Corrientes peaks, rather than one INA slope
decs = []
for e in json.load(open('data/events_1901.json')):
    i = (dt.date.fromisoformat(e['d']) - days[0]).days
    a, b = c[i], c[i+30] if i+30 < len(c) else None
    if b is None or a is None or b <= CMEAN or a <= CMEAN or b >= a: continue
    decs.append(-30.0/math.log((b-CMEAN)/(a-CMEAN)))
KC = st.median(decs)
print(f"        Corrientes recession constant, median of {len(decs)} "
      f"catalogued peaks: {KC:.1f} d")

ext = c[:obs_last+1]
d = days[obs_last]
anchors = [(days[obs_last], c[obs_last])] + INA_C
while d < dt.date(2026, 12, 31):
    d += dt.timedelta(1)
    seg = next(((a, b) for a, b in zip(anchors, anchors[1:]) if a[0] < d <= b[0]), None)
    if seg:
        (da, va), (db, vb) = seg
        ext.append(va + (vb-va)*(d-da).days/(db-da).days)
    else:
        ext.append(CMEAN + (ext[-1]-CMEAN)*math.exp(-1.0/KC))
EDAYS = [days[0] + dt.timedelta(t) for t in range(len(ext))]

f = reservoir(ext, K)
filt = [B0 + B1*f[t-TAU] for t in range(len(f))]

# bases: 1 September 2026
CB, RB = c[i0], r[i0]
print(f"        base on 1 Sep 2026: Corrientes {CB:.2f} m, Rosario {RB:.2f} m")
LAG_KIN = 6                               # a kinematic wave at c = 1.5 v

PROJ = {p['date']: p['mid'] for p in json.load(open('data/projection.json'))['projection']}

rows = []
for t in range(i0, len(ext)):
    d = EDAYS[t]
    if d > dt.date(2026, 12, 20): break
    ros = r[t] if t <= obs_last else PROJ.get(d.isoformat())
    rows.append({
        'd': d.isoformat(),
        'c': round(ext[t], 2),
        'tw': round(ext[t-LAG_KIN] - CB, 2) if t-LAG_KIN >= i0 else None,
        'st': round(filt[t] - filt[i0], 2),
        'ro': round(ros - RB, 2) if ros is not None else None,
        'obs': 1 if t <= obs_last else 0})

tw = [x['tw'] for x in rows if x['tw'] is not None]
sto = [x['st'] for x in rows]
print(f"\n  peak rise above the 1 Sep base, metres at Rosario:")
print(f"    pure translation, {LAG_KIN} d, no attenuation : {max(tw):.2f}")
print(f"    storage filter, k={K:.0f} d                    : {max(sto):.2f}")
print(f"    the page's projection                        : "
      f"{max(x['ro'] for x in rows if x['ro'] is not None):.2f}")
kd = [x['d'] for x in rows if x['st'] == max(sto)][0]
print(f"    storage filter crests {kd}, page says 2026-11-29")

kern = [{'t': t, 'w': round(math.exp(-(t-TAU)/K)/K if t >= TAU else 0.0, 5)}
        for t in range(0, 121)]
half = next(t for t in range(121)
            if sum(k['w'] for k in kern[:t+1]) >= 0.5*sum(k['w'] for k in kern))
print(f"    kernel: half the water has passed by day {half}")

json.dump({'fit': {'k': round(K,1), 'lag': TAU, 'gain': round(B1,4),
                   'b0': round(B0,3), 'r2': round(R2,4), 'sd': round(SD,3),
                   'centroid': round(TAU+K,1), 'half': half,
                   'recession_k': round(KC,1)},
           'base': {'cor': round(CB,2), 'ros': round(RB,2), 'kin_lag': LAG_KIN},
           'kernel': kern, 'rows': rows},
          open('data/kernel.json','w'), indent=1)
print("\nwrote data/kernel.json")
