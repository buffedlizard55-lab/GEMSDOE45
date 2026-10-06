"""GEMSDOE45 — unique GEMS Prize (DOE geothermal fault discovery) submission toolkit.

This package implements a clean-room re-derivation of the competition's
distance-weighted Tversky index (DTI) scoring metric and a novel, fully
reproducible geological-fusion pipeline that produces a submission GeoTIFF.

All equations are taken line-by-line from the official competition problem
description (https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/)
"Performance metric" section, which defines:

    TI(a,b) = sum_x p(x) g(x) / ( sum_x p(x)g(x) + a*sum_x p(x)(1-g(x)) + b*sum_x (1-p(x))g(x) )

with a triangular distance kernel k(d) = max(1 - d/R, 0), R = 300 m (3 px at 100 m),
a = 0.2 (false-positive penalty), b = 0.8 (false-negative penalty).

This package does NOT copy any prior GEMSDOE repository's implementation; it is an
independent re-derivation verified against synthetic ground-truth cases (see
``metric.self_test``).
"""
