#!/usr/bin/env python3
"""Build the Rosario projection series consumed by the HTML report."""
import datetime as dt, json, statistics as st, urllib.request

P = '/tmp/claude-1000/-home-zeta-parana/5456bbce-43b9-4fdc-97d9-fd46dd19c387/scratchpad/dl/ser_%d.json'
IDS = {'corrientes': 19, 'santa_fe': 30, 'rosario': 34}

def load(i):
    d = json.load(open(P % i))
    return {dt.date.fromisoformat(o['timestart'][:10]): float(o['valor'])
            for o in d['observaciones'] if o.get('valor') is not None}

C, SF, R = load(19), load(30), load(34)
start, end = min(min(C), min(R)), max(max(C), max(R))
days = [start + dt.timedelta(d) for d in range((end - start).days + 1)]
c = [C.get(d) for d in days]
r = [R.get(d) for d in days]

# ---- analogue trajectory: Rosario rise relative to its level 5 d before a big Corrientes peak
ev = []
for i in range(50, len(c) - 60):
    if c[i] is None or c[i] < 6.8:
        continue
    if c[i] != max(x for x in c[i-12:i+13] if x is not None):
        continue
    base = [x for x in c[i-32:i-26] if x is not None]
    if not base or c[i] - st.median(base) < 1.5:
        continue
    if ev and i - ev[-1] < 45:
        if c[i] > c[ev[-1]]:
            ev[-1] = i
        continue
    ev.append(i)

traj = {}
analogues = []
for i in ev:
    anch = i - 5
    if r[anch] is None:
        continue
    analogues.append({'peak_date': days[i].isoformat(), 'c_peak': c[i], 'r_base': r[anch]})
    for k in range(-5, 56):
        if i + k < len(r) and r[i+k] is not None:
            traj.setdefault(k, []).append(r[i+k] - r[anch])

q = {}
for k, v in traj.items():
    if len(v) < 5:
        continue
    v = sorted(v); n = len(v)
    q[k] = (v[n//4], st.median(v), v[3*n//4])

# ---- anchors
C_PEAK = dt.date(2026, 10, 20)          # INA central Corrientes peak date
INA_ROS = {dt.date(2026,10,16): (3.35,3.35,3.35),      # observed today (PNA 00:00)
           dt.date(2026,10,20): (3.40,3.73,4.20),
           dt.date(2026,10,27): (3.90,4.39,4.70)}
LAST_INA = dt.date(2026, 10, 27)
ANCHOR = INA_ROS[LAST_INA][1]
k_anchor = (LAST_INA - C_PEAK).days      # +7

out = []
d = dt.date(2026, 10, 16)
while d <= dt.date(2026, 12, 10):
    k = (d - C_PEAK).days
    if d <= LAST_INA:
        # linear interpolation between INA's three published points
        ks = sorted(INA_ROS)
        for a, b in zip(ks, ks[1:]):
            if a <= d <= b:
                f = (d - a).days / max((b - a).days, 1)
                lo, mid, hi = (INA_ROS[a][j] + f * (INA_ROS[b][j] - INA_ROS[a][j]) for j in range(3))
                break
        src = 'ina'
    else:
        if k not in q:
            break
        dlo = q[k][0] - q[k_anchor][0]
        dmid = q[k][1] - q[k_anchor][1]
        dhi = q[k][2] - q[k_anchor][2]
        mid = ANCHOR + dmid
        # widen: analogue spread, floored by the peak regression residual sd (0.27 m),
        # growing with lead time
        grow = 0.27 + 0.004 * (k - k_anchor)
        lo = min(ANCHOR + dlo, mid - grow)
        hi = max(ANCHOR + dhi, mid + grow)
        # blend the INA band width into the model band over the first week past
        # the last published point, so the handoff has no step in it
        f = min(1.0, (k - k_anchor) / 7.0)
        ina_lo_w = ANCHOR - INA_ROS[LAST_INA][0]
        ina_hi_w = INA_ROS[LAST_INA][2] - ANCHOR
        lo = mid - ((1 - f) * ina_lo_w + f * (mid - lo))
        hi = mid + ((1 - f) * ina_hi_w + f * (hi - mid))
        src = 'model'
    out.append({'date': d.isoformat(), 'lo': round(lo, 2), 'mid': round(mid, 2),
                'hi': round(hi, 2), 'src': src})
    d += dt.timedelta(1)

def crossing(series, thr, key):
    for p in series:
        if p[key] >= thr:
            return p['date']
    return None

thresholds = {}
for thr in (4.0, 4.5, 5.0, 5.3):
    thresholds[thr] = {'early': crossing(out, thr, 'hi'),
                       'central': crossing(out, thr, 'mid'),
                       'late': crossing(out, thr, 'lo')}

crest = max(out, key=lambda p: p['mid'])
crest_hi = max(p['hi'] for p in out)
crest_lo_at_peak = crest['lo']

# ---- two-variable peak regression: R_peak ~ C_peak + R_base
def ols(X, y):
    k = len(X[0])
    A = [[sum(X[r][i] * X[r][j] for r in range(len(X))) for j in range(k)]
         + [sum(X[r][i] * y[r] for r in range(len(X)))] for i in range(k)]
    for i in range(k):
        pv = max(range(i, k), key=lambda t: abs(A[t][i])); A[i], A[pv] = A[pv], A[i]
        for j in range(i + 1, k):
            f = A[j][i] / A[i][i]
            for m in range(i, k + 1):
                A[j][m] -= f * A[i][m]
    b = [0.0] * k
    for i in range(k - 1, -1, -1):
        b[i] = (A[i][k] - sum(A[i][j] * b[j] for j in range(i + 1, k))) / A[i][i]
    return b

pk = []
for i in range(50, len(c) - 60):
    if c[i] is None or c[i] < 6.0:
        continue
    if c[i] != max(x for x in c[i-12:i+13] if x is not None):
        continue
    if pk and i - pk[-1] < 45:
        if c[i] > c[pk[-1]]:
            pk[-1] = i
        continue
    pk.append(i)

events = []
for i in pk:
    seg = [(j, r[j]) for j in range(i - 5, min(i + 51, len(r))) if r[j] is not None]
    if len(seg) < 30 or r[i-5] is None:
        continue
    j, rpk = max(seg, key=lambda t: t[1])
    if not (i - 2 < j < min(i + 50, len(r) - 1) - 3):
        continue
    events.append({'date': days[i].isoformat(), 'c_peak': c[i], 'r_base': r[i-5],
                   'r_peak': rpk, 'lag': (days[j] - days[i]).days,
                   'r_peak_date': days[j].isoformat()})

B = ols([[1.0, e['c_peak'], e['r_base']] for e in events], [e['r_peak'] for e in events])
RES = [e['r_peak'] - (B[0] + B[1]*e['c_peak'] + B[2]*e['r_base']) for e in events]
SD = (sum(x*x for x in RES) / (len(events) - 3)) ** 0.5
LAGS = sorted(e['lag'] for e in events)

obs_from = dt.date(2026, 9, 1)
payload = {
    'generated': dt.date(2026, 10, 16).isoformat(),
    'observed': [{'date': d.isoformat(),
                  'corrientes': C.get(d), 'santa_fe': SF.get(d), 'rosario': R.get(d)}
                 for d in [obs_from + dt.timedelta(i) for i in range((dt.date(2026,10,15)-obs_from).days+1)]],
    'today': {'date': '2026-10-16', 'corrientes': 6.23, 'santa_fe': 3.65, 'rosario': 3.35},
    'ina_forecast': {
        'issued': '2026-10-13',
        'corrientes': [{'date':'2026-10-20','lo':6.75,'mid':7.31,'hi':7.55},
                       {'date':'2026-10-27','lo':5.80,'mid':6.38,'hi':6.60}],
        'santa_fe':   [{'date':'2026-10-20','lo':3.60,'mid':4.05,'hi':4.40},
                       {'date':'2026-10-27','lo':4.70,'mid':5.16,'hi':5.50}],
        'rosario':    [{'date':'2026-10-20','lo':3.40,'mid':3.73,'hi':4.20},
                       {'date':'2026-10-27','lo':3.90,'mid':4.39,'hi':4.70}],
    },
    'projection': out,
    'thresholds': thresholds,
    'crest': {'date': crest['date'], 'mid': crest['mid'],
              'lo': crest_lo_at_peak, 'hi': round(crest_hi, 2)},
    'analogues': analogues,
    'events': events,
    'peak_model': {'b0': round(B[0], 3), 'b1': round(B[1], 4), 'b2': round(B[2], 4),
                   'sd': round(SD, 3), 'n': len(events),
                   'lag_median': LAGS[len(LAGS)//2], 'lag_p25': LAGS[len(LAGS)//4],
                   'lag_p75': LAGS[3*len(LAGS)//4]},
}
json.dump(payload, open('/home/zeta/parana/data/projection.json', 'w'), indent=1)

print('crest', payload['crest'])
for t, v in thresholds.items():
    print(t, v)
print('analogues used:', len(analogues))
print('projection tail:', out[-1])
for p in out[::7]:
    print(' ', p['date'], p['lo'], p['mid'], p['hi'], p['src'])
