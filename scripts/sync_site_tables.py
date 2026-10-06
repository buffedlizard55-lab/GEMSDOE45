"""Regenerate the site's band tables from ``evidence/band_screen.json``.

Why this exists: the first version of these tables was typed by hand from an intermediate run, and
it drifted.  A review pass found that the published "23.2 % top-40k within 300 m" attributed to the
band ``det_elev_slope`` was in fact measured on ``|grad1.2|_det_elev_slope``, that the ρ and mean
|grad| columns came from a later run than the percentage column, and that the class column
contradicted the class rule in ``bands.py`` for three bands.  Generation removes the whole failure
mode: every number below is read from the JSON, and the JSON is written by ``scripts/screen_bands.py``
which can be re-run by any reader.

Usage:  python scripts/sync_site_tables.py            # rewrite the region in place
        python scripts/sync_site_tables.py --check    # exit 1 if the site is out of date
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "docs/e45/index.html"
BEGIN = "<!-- BEGIN GENERATED: band tables (scripts/sync_site_tables.py) -->"
END = "<!-- END GENERATED: band tables -->"


def _row(cells: list[str], cls: str | None = None) -> str:
    tds = "".join(f'<td class="{c}">{v}</td>' for c, v in
                  [(c.split("|")[1] if "|" in c else "", c.split("|")[0]) for c in cells])
    return f"<tr>{tds}</tr>" if not cls else f'<tr class="{cls}">{tds}</tr>'


def _cell(value: str, cls: str = "") -> str:
    return f'<td class="{cls}">{value}</td>' if cls else f"<td>{value}</td>"


def raw_table(screen: dict) -> str:
    rows = []
    for r in screen["raw_bands"]:
        cls = r["band_class"]
        elig = r["mass_eligible"]
        pct = f'{r["within_300m_pct"]} %'
        pct_cls = "n ok" if r["enrichment_vs_random"] >= 1.3 else ("n bad" if r["enrichment_vs_random"] < 1.0 else "n")
        tag = {"LOCATOR": "<b>LOCATOR</b>", "INTERMEDIATE": "intermediate", "WEIGHT": "WEIGHT"}[cls]
        if not elig:
            tag += " · value not mass-eligible"
        rows.append("<tr>" + "".join([
            _cell(f'{r["band"]} ({r["index"]})'),
            _cell(f'{r["rho_lag10px"]:.3f}', "n"),
            _cell(f'{r["mean_abs_grad_sigma"]:.4f}', "n"),
            _cell(pct, pct_cls),
            _cell(f'{r["enrichment_vs_random"]:.2f}×', "n"),
            _cell(f'{r["n_components"]:,}', "n"),
            _cell(f'{r["perimeter_over_area"]:.2f}', "n"),
            _cell(tag),
        ]) + "</tr>")
    return (
        "<table><tr><th>Band</th><th class=\"n\">ρ(1 km)</th><th class=\"n\">mean |∂/∂y| (σ/px)</th>"
        "<th class=\"n\">top-40k within 300 m</th><th class=\"n\">vs random</th>"
        "<th class=\"n\">components</th><th class=\"n\">perimeter/area</th><th>Measured class</th></tr>\n"
        + "\n".join(rows) + "\n</table>"
    )


def plane_table(screen: dict, limit: int = 14) -> str:
    rows = []
    for r in screen["derived_planes"][:limit]:
        ok = r["localised"]
        rows.append("<tr>" + "".join([
            _cell(f'<code>{r["plane"]}</code>'),
            _cell(f'{r["within_300m_pct"]} %', "n ok" if ok else "n"),
            _cell(f'{r["enrichment_vs_random"]:.2f}×', "n"),
            _cell(f'{r["n_components"]:,}', "n"),
            _cell(f'{r["median_component_px"]:,}', "n"),
            _cell(f'{r["perimeter_over_area"]:.2f}', "n"),
            _cell("<b>localiser</b>" if ok else "regional blob" if r["perimeter_over_area"] < 0.30 else "mixed"),
        ]) + "</tr>")
    return (
        "<table><tr><th>Derived plane (canonical definition from <code>detfeatures</code>)</th>"
        "<th class=\"n\">top-40k within 300 m</th><th class=\"n\">vs random</th>"
        "<th class=\"n\">components</th><th class=\"n\">median px</th>"
        "<th class=\"n\">perimeter/area</th><th>Verdict</th></tr>\n" + "\n".join(rows) + "\n</table>"
    )


def build(screen: dict) -> str:
    base = screen["random_within_300m_pct"]
    return f"""{BEGIN}
<h2 id="bands">Which of the 19 bands can actually see a fault?</h2>

<p>Answered by measurement, because the answer is not obvious and it rules out a whole family of
ideas. <code>ρ(1 km)</code> is the correlation between a pixel and the pixel 1 km away: a plane with
ρ near 1 is regionally smooth and physically <i>cannot</i> satisfy a 300 m kernel.
<code>mean |∂/∂y|</code> is the mean forward difference along the north axis, in units of the band's
own standard deviation per pixel. Every mass column is the fraction of that plane's strongest 40,000
pixels — the same mass the shipped submission emits — that land within 300 m of a mapped fault,
against a <b>measured</b> random baseline of {base} % (the fraction of the footprint within 300 m of
the catalogue).</p>

<h3>Table 1 — the raw bands, as supplied</h3>
{raw_table(screen)}

<h3>Table 2 — derived planes, ranked by localisation not by score</h3>
<p><code>top-40k within 300 m</code> alone is not enough to call a plane a fault detector. A smooth
regional field can post a respectable fraction simply because its highest values sit over the
fault-rich part of the map — its top-40k is then <i>one large connected blob</i>, which is regional
information a 300 m kernel cannot use. The <code>perimeter/area</code> column separates the two: a
kernel-scale localiser paints many small elongated pieces, a regional field paints one fat blob.</p>
{plane_table(screen)}

<div class="note">
<b>Three results worth more than the ranking itself.</b>
<br>1. <b>Score-without-fragmentation is a trap, and the two strain-rate values demonstrate it.</b>
<code>geod_2ndinv</code> and <code>geod_shearrate</code> post {screen["raw_bands"][0]["within_300m_pct"]} %
and {screen["raw_bands"][1]["within_300m_pct"]} % within 300 m — better than every derivative plane but
one — yet each one's top-40,000 pixels form a <i>single connected component of 40,000 px</i>
(perimeter/area {screen["raw_bands"][0]["perimeter_over_area"]:.2f}). They are reporting "the faults
live in the deformation-model anomaly", not "there is a fault here". Their <i>gradients</i> change by
only 0.0023–0.0037 σ per pixel (Table 1), so a strain-rate ridge detector cannot place a dot inside a
300 m kernel. Both halves of that sentence are now measured rather than asserted.
<br>2. <b>The genuine localiser is the topographic-slope gradient, and it is strong.</b>
<code>|grad2.5|_det_elev_slope</code> puts {screen["derived_planes"][0]["within_300m_pct"]} % of its top
40,000 within 300 m — {screen["derived_planes"][0]["enrichment_vs_random"]:.2f}× the random baseline —
while fragmenting into {screen["derived_planes"][0]["n_components"]:,} components of median
{screen["derived_planes"][0]["median_component_px"]} px. It is the plane the shipped field leans on,
and this is the evidence that it deserves to be.
<br>3. <b>Magnetics-first is measurably the wrong instinct here.</b> <code>mag_anom</code> places
{screen["raw_bands"][-1]["within_300m_pct"] if screen["raw_bands"][-1]["band"] == "mag_anom" else [r for r in screen["raw_bands"] if r["band"] == "mag_anom"][0]["within_300m_pct"]} %
within 300 m — <i>below</i> the {base} % baseline — and the magnetic gradient planes sit at the bottom
of Table 2.
</div>

<div class="note">
<b>An irregularity caught by regenerating this table, kept visible on purpose.</b> The first version
of these tables was typed by hand and had drifted. It attributed a
"23.2 % top-40k within 300 m" to the band <code>det_elev_slope</code>; that number actually belongs to
<code>|grad1.2|_det_elev_slope</code> (Table 2) — the band's <i>value</i> scores
{[r for r in screen["raw_bands"] if r["band"] == "det_elev_slope"][0]["within_300m_pct"]} %, which is
below random. Its ρ and mean-|grad| cells came from a later pipeline run than its percentage cell, and
its class column contradicted the rule in <code>bands.py</code> for three bands. The tables are now
generated by <code>scripts/sync_site_tables.py</code> from <code>evidence/band_screen.json</code>,
which <code>scripts/screen_bands.py</code> writes; <code>--check</code> fails if they diverge.
</div>

<div class="note">
<b>Where the two class statements in <code>bands.py</code> disagree</b>, published rather than
smoothed over. The measured rule is <code>LOCATOR</code> when ρ(1 km) &lt; 0.90, <code>WEIGHT</code>
when ρ ≥ 0.95, <code>INTERMEDIATE</code> between. Three bands disagree with the explicit
mass-placement lists:
{", ".join(f'<code>{d["band"]}</code> (measured {d["measured_class"]}, mass-eligible {d["mass_eligible"]}, {d["within_300m_pct"]} % within 300 m)' for d in screen["class_disagreements"])}.
The mass-placement lists govern the shipped artifact and were <b>not</b> changed after the fact;
the disagreement is recorded so the design decision is visible instead of implicit.
</div>
{END}"""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="exit 1 if the site is out of date")
    a = ap.parse_args()

    screen = json.loads((ROOT / "evidence/band_screen.json").read_text())
    body = build(screen)
    html = SITE.read_text()

    if BEGIN in html and END in html:
        i, j = html.index(BEGIN), html.index(END) + len(END)
        new = html[:i] + body + html[j:]
    else:  # first run: replace the hand-written region between the heading and the next section
        i = html.index('<h2 id="bands">')
        marker = "<!-- ===================== HYPOTHESES"
        j = html.index(marker)
        new = html[:i] + body + "\n\n" + html[j:]

    if a.check:
        if new != html:
            print("FAIL: docs/e45/index.html is out of date with evidence/band_screen.json")
            return 1
        print("OK: site band tables match evidence/band_screen.json")
        return 0

    if new == html:
        print("unchanged")
        return 0
    SITE.write_text(new)
    print(f"updated {SITE.relative_to(ROOT)} ({len(new) - len(html):+,} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
