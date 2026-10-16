#!/usr/bin/env python3
"""Rough extrapolation of Parana river level at Rosario from upstream gauges.

Data sources (all public, no auth):
  - Prefectura Naval Argentina  https://contenidosweb.prefecturanaval.gob.ar/alturas/
  - FICH-UNL CIM               http://wfich1.unl.edu.ar/cim/rios/parana/alturas
  - INA Alerta Hidrologico     https://alerta.ina.gob.ar/a5/obs/puntual/series/<id>
      series 19 = Corrientes, 30 = Santa Fe, 34 = Rosario  (var_id 2 = nivel)
      INA numeric forecasts: /ina/08-PRONOSTICOS/graficos/<station>_prono.csv

Method
  1. Peak-to-peak statistics: for every Corrientes peak >= 6.0 m since 1990,
     find the matching Rosario peak within -5..+50 days. Regress the Rosario
     peak on the Corrientes peak AND the Rosario base level (the base term
     matters: the same Corrientes peak gives a much lower Rosario crest when
     the lower river starts low, as it does now).
  2. Santa Fe -> Rosario transfer (tight, sd 0.18 m, 3-day lag) to convert
     INA's Santa Fe forecast into a Rosario track.
"""
import datetime as dt, json, statistics as st, urllib.request

API = "https://alerta.ina.gob.ar/a5/obs/puntual/series/%d?timestart=1990-01-01&timeend=%s"
IDS = {"corrientes": 19, "santa_fe": 30, "rosario": 34}


def fetch(series_id, cache={}):
    if series_id not in cache:
        url = API % (series_id, dt.date.today().isoformat())
        with urllib.request.urlopen(url, timeout=90) as f:
            d = json.load(f)
        cache[series_id] = {
            dt.date.fromisoformat(o["timestart"][:10]): float(o["valor"])
            for o in d["observaciones"] if o.get("valor") is not None
        }
    return cache[series_id]


def ols(X, y):
    k = len(X[0])
    A = [[sum(X[r][i] * X[r][j] for r in range(len(X))) for j in range(k)]
         + [sum(X[r][i] * y[r] for r in range(len(X)))] for i in range(k)]
    for i in range(k):
        p = max(range(i, k), key=lambda t: abs(A[t][i])); A[i], A[p] = A[p], A[i]
        for j in range(i + 1, k):
            f = A[j][i] / A[i][i]
            for m in range(i, k + 1):
                A[j][m] -= f * A[i][m]
    b = [0.0] * k
    for i in range(k - 1, -1, -1):
        b[i] = (A[i][k] - sum(A[i][j] * b[j] for j in range(i + 1, k))) / A[i][i]
    return b


def peak_model():
    C, R = fetch(IDS["corrientes"]), fetch(IDS["rosario"])
    start, end = min(min(C), min(R)), max(max(C), max(R))
    days = [start + dt.timedelta(d) for d in range((end - start).days + 1)]
    c = [C.get(d) for d in days]
    r = [R.get(d) for d in days]

    peaks = []
    for i in range(50, len(c) - 60):
        if c[i] is None or c[i] < 6.0:
            continue
        if c[i] != max(x for x in c[i - 12:i + 13] if x is not None):
            continue
        if peaks and i - peaks[-1] < 45:
            if c[i] > c[peaks[-1]]:
                peaks[-1] = i
            continue
        peaks.append(i)

    rows = []
    for i in peaks:
        seg = [(j, r[j]) for j in range(i - 5, min(i + 51, len(r))) if r[j] is not None]
        if len(seg) < 30:
            continue
        j, rpk = max(seg, key=lambda t: t[1])
        if not (i - 2 < j < min(i + 50, len(r) - 1) - 3) or r[i - 5] is None:
            continue
        rows.append((c[i], r[i - 5], rpk, (days[j] - days[i]).days, days[i]))

    b = ols([[1.0, x[0], x[1]] for x in rows], [x[2] for x in rows])
    res = [x[2] - (b[0] + b[1] * x[0] + b[2] * x[1]) for x in rows]
    sd = (sum(e * e for e in res) / (len(rows) - 3)) ** 0.5
    lags = sorted(x[3] for x in rows)
    return b, sd, rows, lags


def santafe_transfer():
    SF, R = fetch(IDS["santa_fe"]), fetch(IDS["rosario"])
    out = []
    for lag in range(0, 8):
        xs, ys = [], []
        for d, v in SF.items():
            if d + dt.timedelta(lag) in R:
                xs.append(v); ys.append(R[d + dt.timedelta(lag)])
        b = ols([[1.0, x] for x in xs], ys)
        sd = (sum((y - b[0] - b[1] * x) ** 2 for x, y in zip(xs, ys)) / (len(xs) - 2)) ** 0.5
        out.append((sd, lag, b[0], b[1], len(xs)))
    return min(out)


if __name__ == "__main__":
    b, sd, rows, lags = peak_model()
    print("Corrientes peak -> Rosario peak, events since 1990 (n=%d)" % len(rows))
    for x in rows:
        print("   %s  C %.2f  R_base %.2f  R_peak %.2f  lag %2dd" % (x[4], x[0], x[1], x[2], x[3]))
    print("\nR_peak = %.2f + %.3f*C_peak + %.3f*R_base   (resid sd %.2f m)" % (b[0], b[1], b[2], sd))
    print("peak-to-peak lag: median %.0f d, p25 %d, p75 %d"
          % (st.median(lags), lags[len(lags) // 4], lags[3 * len(lags) // 4]))

    r_now = fetch(IDS["rosario"])[max(fetch(IDS["rosario"]))]
    print("\nRosario now: %.2f m" % r_now)
    for cpk in (6.75, 7.31, 7.55, 8.00):
        print("   Corrientes peak %.2f -> Rosario crest %.2f +/- %.2f m"
              % (cpk, b[0] + b[1] * cpk + b[2] * r_now, sd))

    sd2, lag, a2, b2, n = santafe_transfer()
    print("\nSanta Fe -> Rosario: R(t+%dd) = %+.2f + %.3f*SF(t)  sd %.2f m  n=%d"
          % (lag, a2, b2, sd2, n))
    print("   Rosario 5.00 m (alerta) <=> Santa Fe %.2f m" % ((5.00 - a2) / b2))
