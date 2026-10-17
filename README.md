# Paraná at Rosario: crest estimate, October 2026

**`index.html`** is the page to publish: one self-contained file, no build step and
no local assets. Drop it on any static host (GitHub Pages, Netlify, S3, nginx) and
it works. The only external request is the Google Fonts stylesheet; without network
it falls back to Georgia and the system sans and still reads correctly.

Live in two editions:

- English: https://ezequielgarcia.github.io/parana-rosario/
- Español: https://ezequielgarcia.github.io/parana-rosario/es/

`es/index.html` is a full translation, not a machine pass. It also localises what
Argentine hydrological publications localise: decimal comma, dot thousands
separator, Spanish month abbreviations, and dates written "16 de octubre de 2026".
The two pages cross-link and carry `hreflang` alternates. Keep them in step: a
change to the data constants has to land in both.

## Question

A flood wave from rain on the unregulated upper Paraná (above Itaipú) and the
Iguazú reached Corrientes in mid-October 2026. How high will it run at Rosario,
and when?

## Answer as of 16 Oct 2026

| | |
|---|---|
| Rosario now | 3.35 m (alerta 5.00, evacuación 5.30) |
| Crest | 4.76 ± 0.33 m (68%); 4.22 to 5.31 at 90% |
| Crest date | 15 Nov ± 8 d, broad flat crest ~8 to 20 Nov |
| Reaches 5.00 m | 24% chance; 5.30 m is 5% |

Three independent estimators, all on 1901-2026: a two-variable peak regression
gives 4.76 m, the attenuation law gives 4.56 m, the seven base-matched analogues
give 4.33 m. Leave-one-out over all 76 waves puts the typical error at 0.19 m
with no bias on waves of this size.

### Error budget

`scripts/error_budget.py` writes `data/error_budget.json` and prints the terms.
One standard deviation each, combined in quadrature because they are independent:

| Term | sigma |
|---|---|
| wave-to-wave scatter about the fit | 0.241 m |
| uncertainty in the fitted coefficients | 0.058 m |
| INA's own range on the Corrientes peak, through the 0.499 slope | 0.121 m |
| starting base at Rosario | 0.021 m |
| spread between the three estimators | 0.176 m |
| **combined** | **0.328 m** |

The regression alone gives a prediction standard error of 0.25 m, which is what
the leave-one-out backtest reproduces (sd 0.25 m, central 90% of errors -0.28 to
+0.47 m). The backtest cannot charge for the forecast error on the upstream peak
or for the disagreement between methods, so the published figure is wider.
Travel time fitted against peak size gives 25 days with a residual sd of 8 days;
the sixteen comparable waves ran from 3 to 37 days.

### Routing: the reach is storage, not translation

`scripts/basis_tests.py` and `scripts/basis_tests2.py` test three ideas against
the 1901-2026 daily record. `scripts/build_kernel.py` fits the one that works and
writes `data/kernel.json`.

A Fourier basis is wrong for the event. Annual harmonics explain 14.5% of
Rosario's daily variance and six harmonics only reach 15.3%. With the seasonal
cycle removed, the periodogram of the anomaly has 0.0003 of the variance at a
90-day period and 0.0000 at 60, 30, 14 and 7 days. A flood is one transient
pulse, so there is nothing periodic to expand, and a periodic basis extrapolates
periodically.

A lag-and-reservoir transfer function works well:

```
Rosario = -0.368 + 0.9195 * S(Corrientes)[t-3]      S = one linear reservoir, k = 21.2 d
R2 0.8996   residual sd 0.422 m   n = 45 544 days   centroid travel time 24.2 d
```

One reservoir beats two or more, which means maximum dispersion. Consequences:

| | |
|---|---|
| Celerity implied | 792 km / 24.2 d = 33 km/d = 0.38 m/s |
| A kinematic wave (c = 1.5 v, v ~ 1 m/s) | would arrive in 6 days |
| Crest breadth, Rosario / Corrientes | 4.2 times wider |
| Fastest 5-day rise | 0.82 m at Corrientes, 0.18 m at Rosario |
| Attenuation of an 11-day pulse | 0.9195 * (1 - exp(-11/21.2)) = 0.37, measured 0.36 |

It is therefore not a travelling wave: four times too slow and the shape is not
conserved. The floodplain between Corrientes and Santa Fe is tens of km wide
against a channel of one or two, and that storage sets both the delay and the
attenuation. Applied to this wave, pure translation would put Rosario at 7.51 m,
above the 1992 record; the storage filter says 4.94 m and the peak-to-peak
regression says 4.76 m.

The filter is not used for the published figure. It describes the whole
hydrograph better but predicts the crest worse: MAE 0.32 m with a 0.25 m low bias
against 0.19 m for the peak regression, because it has one input, no local
tributaries and no term for the starting state. Transforming does not help
either: log stage gives R2 0.884 and a discharge-like h^1.6 gives 0.897, against
0.900 for plain stage.

### Revision, 16 Oct 2026

The first version fitted on 1990-2026 and published 5.16 m, because the first API
pull passed `timestart=1990-01-01` and that was taken as the extent of the record.
It is not: Corrientes serves from 1901-01-02 and Rosario from 1900-01-02, 99%
daily coverage, 76 matched waves instead of 23. The 1990s sample happened to hold
mostly waves that started from a much higher Rosario base, which inflated the
crest by about 0.4 m. Always pass an early `timestart` and read back what the API
actually returns.

## Files

- `scripts/rosario_forecast.py` fits the peak model and the Santa Fe transfer
  straight from the INA a5 API. Run it to refresh the statistics.
- `scripts/build_projection2.py` builds `data/projection.json`, the dataset the
  pages are drawn from, on the full record. Edit the `INA_ROS` / `C_PEAK` anchors
  when a newer INA run posts, then re-run. Supersedes `build_projection.py`.
- `scripts/backtest.py` catalogues all 76 waves since 1901 and cross-validates
  the fit. Writes `data/events_1901.json`.
- `scripts/wave_behaviour.py` compares the record waves against the base-matched
  ones and tests lag and attenuation against peak size and starting level.
  Writes `data/wave_shapes.json`.
- `scripts/error_budget.py` builds the five-term error budget, checks it against
  the leave-one-out error quantiles and writes `data/error_budget.json`.
- `scripts/basis_tests.py` and `scripts/basis_tests2.py` test the harmonic basis,
  the periodogram of the anomaly, the lag-and-reservoir fit and three transformed
  spaces. They print; they write no data.
- `scripts/build_kernel.py` fits the routing filter and writes `data/kernel.json`
  for plates 6 and 7.
- `scripts/analyze.py` holds the first exploratory pass: cross-correlation of daily
  increments, event matching, single-variable fits. Kept for the working.
- `data/projection.json` holds the projection band, threshold dates and the
  base-matched analogues.
- `data/events_1901.json` holds all 76 matched waves since 1901.
- `data/wave_shapes.json` holds the day-by-day trajectories behind plate 4.
- `data/kernel.json` holds the fitted filter and the three tracks in plate 7.
- `data/error_budget.json` holds the error terms and threshold probabilities.

## Publishing

    python3 -m http.server -d . 8000     # check it locally at /index.html

The page is theme-aware: it follows the reader's light or dark setting, and honours
`data-theme="light"` or `data-theme="dark"` stamped on `<html>` if you want to pin one.
It carries its own favicon, description and Open Graph tags, so a shared link previews
properly. Nothing on the page phones home and there is no analytics.

## Data sources (all public, no auth)

- Prefectura Naval Argentina: https://contenidosweb.prefecturanaval.gob.ar/alturas/
- FICH-UNL CIM: http://wfich1.unl.edu.ar/cim/rios/parana/alturas (HTTP only, 443 is closed)
- INA Alerta Hidrológico daily report:
  https://alerta.ina.gob.ar/ina/06-INFORMES/diario/pdf/reporte_diario_YYYY-MM-DD.pdf
- INA numeric forecasts: https://alerta.ina.gob.ar/ina/08-PRONOSTICOS/graficos/
- INA a5 series API: https://alerta.ina.gob.ar/a5/obs/puntual/series/<id>?timestart=&timeend=
  Series 19 Corrientes, 30 Santa Fe, 34 Rosario; var_id 2 is level in metres.
  Daily records back to 1990.

## Refresh

INA reissues the middle-Paraná forecast Tuesdays (7 and 14 days) and Fridays
(4 and 11 days). To update: pull the latest `*_prono.csv`, put the new Corrientes
peak and Rosario points into `build_projection.py`, re-run both scripts, then
update the constants near the top of the `<script>` block (`COR`, `SFE`, `ROS`, `INA`,
`PROJ`, `THRESHOLDS`, `GAUGES`) in both `index.html` and `parana-rosario.html`, and the
four reading tiles in the markup.
