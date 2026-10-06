"""H51 candidate generation: fault-zone strands, tip continuations, and the KM survival law.

OFFICIAL SEMANTICS THIS IMPLEMENTS (DrivenData staff, quote-verified)
--------------------------------------------------------------------
* thread 11536, chrisk-dd, 2026-09-23: "'new fault' means 'any fault pixel not already captured
  by USGS/INGENIOUS' and can include newly mapped geometry of an existing fault system."
* thread 11516 post 4, chrisk-dd, 2026-09-21: "The mask is indeed pixel-exact - it is identical
  to the provided set of training fault labels." / "A predicted pixel that is near a known fault
  trace but far from a new-fault ground truth pixel will be fully penalized, i.e., the buffer does
  not apply to known faults." / "A new-fault ground truth pixel can indeed lie within 300m of a
  known fault trace. Such pixels would constitute corrections or modifications to existing fault
  traces."

So the target population is: along-strike continuations past mapped tips, splays and parallel
strands inside the same fault zone, and local corrections of the mapped trace -- every one of which
is *geometry attached to the existing fault system*, not an unrelated structure in blank ground.
That is exactly the geometry this module enumerates.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

from .h51 import line_response


def skeleton(mask: np.ndarray) -> np.ndarray:
    """Thin a binary mask to a 1-px skeleton (skimage, imported lazily)."""
    from skimage.morphology import skeletonize
    return skeletonize(mask)


def strike_field(mask: np.ndarray, sigma: float = 3.0) -> tuple[np.ndarray, np.ndarray]:
    """Local strike of a thin line mask, from the smoothed structure tensor.

    Returns ``(tangent_r, tangent_c)`` -- a unit vector along the line at each pixel.  The mask's
    density gradient is perpendicular to the line, so the strike is the gradient rotated by 90
    degrees.  Orientation is doubled before smoothing and halved afterwards so that the reverse of
    a line gives the same line.
    """
    m = mask.astype(np.float32)
    gy = ndimage.gaussian_filter(m, sigma, order=(1, 0))
    gx = ndimage.gaussian_filter(m, sigma, order=(0, 1))
    tr = -gx
    tc = gy
    norm = np.hypot(tr, tc) + 1e-12
    return tr / norm, tc / norm


def tip_continuation_candidates(shape, tips, survival, max_extend: int,
                                length_gain: float = 0.0):
    """Weighted candidates along each tip's outward strike, weighted by the KM survival S(e).

    The weight is P(the structure is still a fault e px beyond the tip) = S(e) from the
    product-limit fit of the catalogue's own relay-gap sample.  ``length_gain`` adds a
    size-dependent latitude: a longer mapped segment is statistically better attested, so its
    expected extension is longer.  With ``length_gain = 0`` the weight is exactly S(e) for every
    tip; the *realised* extension therefore still varies tip to tip because S is a curve, not a
    constant -- see ``report`` in the builder for the measured spread.
    """
    rows, cols, wts, tipid = [], [], [], []
    h, w = shape
    for i, t in enumerate(tips):
        lat = 1.0 + length_gain * np.log1p(t.seg_len / 12.0)
        for e in range(1, max_extend + 1):
            p = float(survival(np.array([float(e)]))[0])
            if p <= 1e-6:
                break
            r = int(round(t.row + e * t.urow))
            c = int(round(t.col + e * t.ucol))
            if not (0 <= r < h and 0 <= c < w):
                break
            rows.append(r)
            cols.append(c)
            wts.append(min(1.0, p * lat))
            tipid.append(i)
    return np.asarray(rows, np.int32), np.asarray(cols, np.int32), \
        np.asarray(wts, np.float64), np.asarray(tipid, np.int32)


def parallel_strand_candidates(skel: np.ndarray, tr: np.ndarray, tc: np.ndarray,
                               offsets=(1, 2, 3), presence_discount: float = 0.6):
    """Perpendicular offsets from every skeleton pixel -- splays and parallel strands.

    The normal to the local strike is ``(-tc, tr)``.  Each offset is a distinct candidate strand
    position.  The weight decreases with offset through ``presence_discount ** (offset - 1)``,
    a *documented* prior that an adjacent strand is more likely than one further out; it is not a
    fitted parameter and the emission step is free to reject any of them.
    """
    rr, cc = np.nonzero(skel)
    trv, tcv = tr[rr, cc], tc[rr, cc]
    rows, cols, wts = [], [], []
    h, w = skel.shape
    for k, off in enumerate(offsets):
        r = np.rint(rr - off * tcv).astype(np.int64)
        c = np.rint(cc + off * trv).astype(np.int64)
        ok = (r >= 0) & (r < h) & (c >= 0) & (c < w)
        rows.append(r[ok])
        cols.append(c[ok])
        wts.append(np.full(int(ok.sum()), presence_discount ** k))
    return (np.concatenate(rows), np.concatenate(cols), np.concatenate(wts))


def evidence_field(stack: np.ndarray, footprint: np.ndarray, band_index: int,
                   sigmas=(1.2, 2.5), pct: float = 99.0) -> np.ndarray:
    """Normalised topographic/gravity/magnetic ridge response used only to *gate* candidates."""
    out = np.zeros(stack.shape[1:], dtype=np.float32)
    for s in sigmas:
        b = stack[band_index].astype(np.float32)
        v = b[footprint]
        z = (b - float(v.mean())) / (float(v.std()) + 1e-12)
        r = line_response(z, s)
        out += np.clip(r / (float(np.percentile(r[footprint], pct)) + 1e-30), 0.0, 1.0)
    out /= len(sigmas)
    return out
