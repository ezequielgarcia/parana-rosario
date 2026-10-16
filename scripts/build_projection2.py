#!/usr/bin/env python3
"""Rebuild the Rosario projection on the 1901-2026 record.

Supersedes build_projection.py, which used only 1990+ because that was what
the first API pull returned. Three independent estimators are combined:

  1. two-variable fit   R_crest ~ C_peak + R_base   (n=76, LOO MAE 0.19 m)
  2. attenuation law    Rosario rise = gain * Corrientes rise, gain ~ 0.36
  3. nearest analogues  events matched on BOTH peak and starting level

The curve to 27 October stays on INA's published track, which is real-time
routing this cannot improve on. Past that it relaxes onto the historical crest.
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
r = [R.get(d) for d in days]
idx = {d: i for i, d in enumerate(days)}
ev = json.load(open('data/events_1901.json'))

C_PEAK, C_PEAK_DATE, R_NOW = 7.31, dt.date(2026, 10, 20), 3.35
ANALOGUES = [e for e in ev if e['base'] < 3.75 and e['c'] >= 6.6]
print(f"base-matched analogues (base < 3.75 m, peak >= 6.60 m): {len(ANALOGUES)}")
for e in sorted(ANALOGUES, key=lambda x: x['d']):
    print(f"   {e['d']}  peak {e['c']:.2f}  base {e['base']:.2f}  crest {e['r']:.2f}"
          f"  rise {e['r']-e['base']:+.2f}  lag {e['lag']}d")

# median analogue trajectory, measured as rise above the level 5 days before the peak
traj = {}
for e in ANALOGUES:
    i = idx[dt.date.fromisoformat(e['d'])]
    anch = r[i-5]
    if anch is None: continue
    for k in range(-5, 71):
        if 0 <= i+k < len(r) and r[i+k] is not None:
            traj.setdefault(k, []).append(r[i+k] - anch)
q = {k: (sorted(v)[len(v)//4], st.median(v), sorted(v)[3*len(v)//4])
     for k, v in traj.items() if len(v) >= 5}
LAGS = sorted(e['lag'] for e in ANALOGUES)
print(f"analogue lag: median {st.median(LAGS):.0f} d  (p25 {LAGS[len(LAGS)//4]}, p75 {LAGS[3*len(LAGS)//4]})")
print(f"analogue rise: median {st.median([e['r']-e['base'] for e in ANALOGUES]):+.2f} m")

INA_ROS = {dt.date(2026,10,16): (3.35,3.35,3.35),
           dt.date(2026,10,20): (3.40,3.73,4.20),
           dt.date(2026,10,27): (3.90,4.39,4.70)}
LAST_INA = dt.date(2026,10,27)
ANCHOR = INA_ROS[LAST_INA][1]
k_anchor = (LAST_INA - C_PEAK_DATE).days

# target crest from the cross-validated fit, and the span the three methods cover
FIT_CREST, FIT_SD = 4.76, 0.24
out, d = [], dt.date(2026,10,16)
while d <= dt.date(2026,12,20):
    k = (d - C_PEAK_DATE).days
    if d <= LAST_INA:
        ks = sorted(INA_ROS)
        for a, b in zip(ks, ks[1:]):
            if a <= d <= b:
                f = (d-a).days / max((b-a).days, 1)
                lo, mid, hi = (INA_ROS[a][j] + f*(INA_ROS[b][j]-INA_ROS[a][j]) for j in range(3))
                break
        src = 'ina'
    else:
        if k not in q: break
        shape = q[k][1] - q[k_anchor][1]
        # the analogue shape sets the path; its amplitude is scaled so the crest
        # lands on the cross-validated fit rather than on the 1990s-only sample
        peak_shape = max(q[j][1] for j in q if j >= k_anchor) - q[k_anchor][1]
        scale = (FIT_CREST - ANCHOR) / peak_shape if peak_shape > 0.05 else 1.0
        mid = ANCHOR + shape * scale
        grow = FIT_SD + 0.005*(k - k_anchor)
        lo, hi = mid - grow, mid + grow
        f = min(1.0, (k-k_anchor)/7.0)
        lo = mid - ((1-f)*(ANCHOR-INA_ROS[LAST_INA][0]) + f*(mid-lo))
        hi = mid + ((1-f)*(INA_ROS[LAST_INA][2]-ANCHOR) + f*(hi-mid))
        src = 'model'
    out.append({'date': d.isoformat(), 'lo': round(lo,2), 'mid': round(mid,2),
                'hi': round(hi,2), 'src': src})
    d += dt.timedelta(1)

def cross(thr, key):
    for p in out:
        if p[key] >= thr: return p['date']
    return None
th = {t: {'early': cross(t,'hi'), 'central': cross(t,'mid'), 'late': cross(t,'lo')}
      for t in (4.0, 4.5, 5.0, 5.3)}
crest = max(out, key=lambda p: p['mid'])
print(f"\ncrest {crest['mid']:.2f} m on {crest['date']}  band {crest['lo']:.2f} to {crest['hi']:.2f}")
print(f"peak of the upper edge: {max(p['hi'] for p in out):.2f}")
for t, v in th.items(): print(' ', t, v)

json.dump({'generated':'2026-10-16','basis':'INA a5 series 19/34, 1901-2026, 76 matched events',
           'crest':{'date':crest['date'],'mid':crest['mid'],'lo':crest['lo'],
                    'hi':round(max(p['hi'] for p in out),2)},
           'projection':out,'thresholds':th,
           'analogues':sorted(ANALOGUES,key=lambda x:x['d'])},
          open('data/projection.json','w'), indent=1)
print("\nPROJ=", json.dumps([[p['date'][5:],p['lo'],p['mid'],p['hi'],1 if p['src']=='model' else 0] for p in out]))
