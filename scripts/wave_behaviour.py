#!/usr/bin/env python3
"""Do big waves and medium waves behave the same way in this reach?

Compares the four largest Corrientes peaks on record against three events
matched to today's low Rosario base, and tests lag and attenuation against
peak size and starting level across all 76 events.
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
c = [C.get(d) for d in days]; r = [R.get(d) for d in days]
idx = {d: i for i, d in enumerate(days)}
ev = json.load(open('data/events_1901.json'))

def corr(xs, ys):
    n = len(xs); mx = sum(xs)/n; my = sum(ys)/n
    num = sum((a-mx)*(b-my) for a, b in zip(xs, ys))
    den = (sum((a-mx)**2 for a in xs)*sum((b-my)**2 for b in ys))**0.5
    return num/den

# ---- per-event descriptors
for e in ev:
    i = idx[dt.date.fromisoformat(e['d'])]
    cbase = st.median([x for x in c[i-40:i-25] if x is not None] or [e['c']])
    e['c_rise'] = e['c'] - cbase
    e['r_rise'] = e['r'] - e['base']
    e['gain'] = e['r_rise']/e['c_rise'] if e['c_rise'] > 0.3 else None
    # breadth: days within 0.25 m of the Rosario crest
    seg = [r[j] for j in range(i-5, min(i+61, len(r))) if r[j] is not None]
    e['breadth'] = sum(1 for v in seg if v >= e['r'] - 0.25)
    # rate of the steep limb at Rosario: max 5-day rise
    best = 0
    for j in range(i-5, min(i+50, len(r)-5)):
        if r[j] is not None and r[j+5] is not None:
            best = max(best, r[j+5]-r[j])
    e['max5'] = best

BIG = ['1983-03-05', '1992-06-08', '1905-06-05', '1998-05-04']
SAME_BASE = ['1923-06-27', '2013-07-04', '1990-02-01']

print("=== group descriptors ===")
print("  event        C_pk  C_rise  R_base  R_crest  R_rise  gain  lag  breadth  max5d")
for tag, keys in (('BIG', BIG), ('SAME BASE', SAME_BASE)):
    print(f"-- {tag}")
    for k in keys:
        e = next(x for x in ev if x['d'] == k)
        print(f"  {e['d']}  {e['c']:5.2f}  {e['c_rise']:+5.2f}   {e['base']:5.2f}   {e['r']:5.2f}"
              f"  {e['r_rise']:+5.2f}  {e['gain']:.2f}  {e['lag']:3d}d   {e['breadth']:3d}d  {e['max5']:+.2f}")

g1 = [next(x for x in ev if x['d'] == k) for k in BIG]
g2 = [next(x for x in ev if x['d'] == k) for k in SAME_BASE]
for tag, g in (('BIG', g1), ('SAME BASE', g2)):
    print(f"  {tag:10} median lag {st.median([e['lag'] for e in g]):.0f}d, "
          f"gain {st.median([e['gain'] for e in g]):.2f}, "
          f"breadth {st.median([e['breadth'] for e in g]):.0f}d, "
          f"max5d {st.median([e['max5'] for e in g]):+.2f} m")

print("\n=== across all 76 events: what sets the travel time? ===")
ok = [e for e in ev if e['gain'] is not None]
print(f"  corr(lag, Corrientes peak)   {corr([e['c'] for e in ev], [e['lag'] for e in ev]):+.2f}")
print(f"  corr(lag, Rosario base)      {corr([e['base'] for e in ev], [e['lag'] for e in ev]):+.2f}")
print(f"  corr(gain, Rosario base)     {corr([e['base'] for e in ok], [e['gain'] for e in ok]):+.2f}")
print(f"  corr(gain, Corrientes rise)  {corr([e['c_rise'] for e in ok], [e['gain'] for e in ok]):+.2f}")

print("\n  lag by Rosario starting level:")
for lo, hi in ((0, 3.5), (3.5, 4.5), (4.5, 5.5), (5.5, 9)):
    g = [e for e in ev if lo <= e['base'] < hi]
    if len(g) < 3: continue
    print(f"    base {lo:.1f}-{hi:.1f} m  n={len(g):2}  lag median {st.median([e['lag'] for e in g]):4.0f} d"
          f"   gain median {st.median([e['gain'] for e in g if e['gain'] is not None]):.2f}"
          f"   crest {st.median([e['r'] for e in g]):.2f} m")

print("\n  lag by Corrientes peak size:")
for lo, hi in ((6.0, 6.5), (6.5, 7.0), (7.0, 7.5), (7.5, 9.5)):
    g = [e for e in ev if lo <= e['c'] < hi]
    if len(g) < 3: continue
    print(f"    peak {lo:.1f}-{hi:.1f} m  n={len(g):2}  lag median {st.median([e['lag'] for e in g]):4.0f} d"
          f"   gain median {st.median([e['gain'] for e in g if e['gain'] is not None]):.2f}")

# ---- normalised trajectories for the figure
traj = {}
for tag, keys in (('big', BIG), ('base', SAME_BASE)):
    for k in keys:
        e = next(x for x in ev if x['d'] == k)
        i = idx[dt.date.fromisoformat(k)]
        traj[k] = {'tag': tag, 'c': e['c'], 'base': e['base'], 'r': e['r'], 'lag': e['lag'],
                   'gain': round(e['gain'], 2), 'breadth': e['breadth'],
                   'ros': [None if not (0 <= i+n < len(r)) or r[i+n] is None
                           else round(r[i+n] - e['base'], 2) for n in range(-20, 61)],
                   'cor': [None if not (0 <= i+n < len(c)) or c[i+n] is None
                           else round(c[i+n], 2) for n in range(-20, 61)]}
json.dump({'offset0': 20, 'events': traj}, open('data/wave_shapes.json', 'w'), indent=1)
print("\nwrote data/wave_shapes.json")
