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
| Crest | 4.76 m, band 4.42 to 5.10 |
| Crest date | 15 Nov, broad flat crest ~8 to 20 Nov |
| Crosses 5.00 m | only at the upper edge of the band, around 9 Nov |

Three independent estimators, all on 1901-2026: a two-variable peak regression
gives 4.76 ± 0.24 m, the attenuation law gives 4.56 m, the seven base-matched
analogues give 4.33 m. Leave-one-out over all 76 waves puts the typical error at
0.19 m with no bias on waves of this size.

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
- `scripts/analyze.py` holds the first exploratory pass: cross-correlation of daily
  increments, event matching, single-variable fits. Kept for the working.
- `data/projection.json` holds the projection band, threshold dates and the
  base-matched analogues.
- `data/events_1901.json` holds all 76 matched waves since 1901.
- `data/wave_shapes.json` holds the day-by-day trajectories behind plate 4.

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
