# Paraná at Rosario — crest estimate, October 2026

**`index.html`** is the page to publish: one self-contained file, no build step and
no local assets. Drop it on any static host (GitHub Pages, Netlify, S3, nginx) and
it works. The only external request is the Google Fonts stylesheet; without network
it falls back to Georgia and the system sans and still reads correctly.

Live at https://ezequielgarcia.github.io/parana-rosario/

## Question

A flood wave from rain on the unregulated upper Paraná (above Itaipú) and the
Iguazú reached Corrientes in mid-October 2026. How high will it run at Rosario,
and when?

## Answer as of 16 Oct 2026

| | |
|---|---|
| Rosario now | 3.35 m (alerta 5.00, evacuación 5.30) |
| Crest | 5.16 m, band 4.80 to 5.52 |
| Crest date | 18 Nov, broad flat crest ~10 to 23 Nov |
| Crosses 5.00 m | 10 Nov central, 4 Nov on the fast edge |

Two independent methods: a two-variable peak regression gives 4.85 ± 0.27 m,
an INA-anchored analogue trajectory gives 5.16 m. Both land on 18 November.

## Files

- `scripts/rosario_forecast.py` — fits the peak model and the Santa Fe transfer
  straight from the INA a5 API. Run it to refresh the statistics.
- `scripts/build_projection.py` — builds `data/projection.json`, the dataset the
  HTML page is drawn from. Edit the `INA_ROS` / `C_PEAK` anchors when a newer
  INA run posts, then re-run.
- `scripts/analyze.py` — the first exploratory pass: cross-correlation of daily
  increments, event matching, single-variable fits. Kept for the working.
- `data/projection.json` — observed series, INA forecasts, projection band,
  threshold dates, the 23 matched events.

## Publishing

    python3 -m http.server -d . 8000     # check it locally at /index.html

The page is theme-aware: it follows the reader's light or dark setting, and honours
`data-theme="light"` or `data-theme="dark"` stamped on `<html>` if you want to pin one.
It carries its own favicon, description and Open Graph tags, so a shared link previews
properly. Nothing on the page phones home and there is no analytics.

## Data sources (all public, no auth)

- Prefectura Naval Argentina — https://contenidosweb.prefecturanaval.gob.ar/alturas/
- FICH-UNL CIM — http://wfich1.unl.edu.ar/cim/rios/parana/alturas (HTTP only, 443 is closed)
- INA Alerta Hidrológico daily report —
  https://alerta.ina.gob.ar/ina/06-INFORMES/diario/pdf/reporte_diario_YYYY-MM-DD.pdf
- INA numeric forecasts — https://alerta.ina.gob.ar/ina/08-PRONOSTICOS/graficos/
- INA a5 series API — https://alerta.ina.gob.ar/a5/obs/puntual/series/<id>?timestart=&timeend=
  Series 19 Corrientes, 30 Santa Fe, 34 Rosario; var_id 2 is level in metres.
  Daily records back to 1990.

## Refresh

INA reissues the middle-Paraná forecast Tuesdays (7 and 14 days) and Fridays
(4 and 11 days). To update: pull the latest `*_prono.csv`, put the new Corrientes
peak and Rosario points into `build_projection.py`, re-run both scripts, then
update the constants near the top of the `<script>` block (`COR`, `SFE`, `ROS`, `INA`,
`PROJ`, `THRESHOLDS`, `GAUGES`) in both `index.html` and `parana-rosario.html`, and the
four reading tiles in the markup.
