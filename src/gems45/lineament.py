"""Multi-physics oriented-lineament consensus (H47).

WHY THIS AND NOT A SINGLE BAND
-----------------------------
A single potential-field derivative lights up every contact, dyke, road and terrane boundary.  What
makes a *fault* different from those is that several independent physical fields -- magnetics,
gravity, topography -- resolve the same line, in the same place, with the same strike.  So the
detector here does not threshold one band; it asks how well the independent fields agree.

The agreement statistic is the second-order orientation order parameter

    R(x) = | sum_p w_p(x) exp(2 i theta_p(x)) | / sum_p w_p(x)

where ``theta_p`` is the local line direction of band p (from the eigenvector of the structure
tensor of that band alone) and ``w_p`` is that band's line strength.  ``R`` is 1 when every band
that sees a line agrees on its azimuth, and 0 when the azimuths cancel.  Orientations are mapped to
the double angle so that a line and its reverse are the same line.  This object does not exist
anywhere in the prior repositories: they take per-band ridge responses and pixel-wise products,
which cannot distinguish "one strong field" from "many agreeing fields".

PHYSICS OF EACH TERM
--------------------
* ``tmi_hg``, ``tmi_vg``, ``iso_grav_anom_hg``, ``iso_grav_anom_vg`` -- directional derivatives of
  the magnetic and gravity fields.  A fault is a step in the source distribution, so it appears as a
  line in the derivative and as a *zero crossing* in the field itself.
* ``det_elev_slope`` -- the surface expression of the scarp.
* ``tc`` (tilt angle) -- an amplitude-normalised magnetic edge detector, which locates source edges
  to first order independently of the field strength.
* ``rtp``, ``tmi``, ``mag_anom`` -- the source fields themselves, carrying the step.

Both polarities are kept: a magnetic low and a magnetic high are equally good markers of a contact,
and the sign of the gravity step depends on which side of the fault is uplifted.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

from .grid import BAND_NAMES

LOCATOR = ["tmi_hg", "tmi_vg", "iso_grav_anom_hg", "iso_grav_anom_vg", "det_elev_slope", "tc",
           "rtp", "tmi"]
SCALES = (1.2, 2.5)
ELONG = 0.5          # Sato-style elongation constant


def _line_response(z: np.ndarray, sigma: float) -> tuple[np.ndarray, np.ndarray]:
    """Sato/Frangi line response and local line orientation for one band.

    Returns ``(resp, theta)`` where ``resp`` is the polarity-symmetric line strength and ``theta``
    is the direction of the line (eigenvector of the smaller-magnitude Hessian eigenvalue), in
    radians, modulo pi.
    """
    s = ndimage.gaussian_filter(z, sigma)
    hyy, hyx = np.gradient(np.gradient(s, axis=0), axis=0), None
    # Hessian components
    hxx = np.gradient(np.gradient(s, axis=1), axis=1)
    hyy = np.gradient(np.gradient(s, axis=0), axis=0)
    hxy = np.gradient(np.gradient(s, axis=0), axis=1)
    trace = hxx + hyy
    det = hxx * hyy - hxy * hxy
    disc = np.sqrt(np.maximum(trace * trace / 4.0 - det, 0.0))
    l1 = trace / 2.0 + disc
    l2 = trace / 2.0 - disc
    # |l2| large and |l1| small => a line; |l1| large => a blob => suppress
    resp = np.abs(l2) * np.exp(-(l1 * l1) / (2.0 * ELONG * ELONG * (l2 * l2 + 1e-30)))
    # eigenvector for l2: (hxy, l2 - hxx) or (l2 - hyy, hxy)
    vx = hxy
    vy = l2 - hxx
    norm = np.hypot(vx, vy)
    theta = np.arctan2(vy, vx)
    theta = np.where(norm > 1e-30, theta, 0.0)
    return resp, theta


def consensus_field(stack: np.ndarray, footprint: np.ndarray,
                    scales: tuple[float, ...] = SCALES) -> tuple[np.ndarray, np.ndarray]:
    """Multi-scale, multi-physics line strength and azimuth agreement.

    Returns ``(strength, R)``: ``strength`` is the geometric mean over bands of the normalised line
    response, ``R`` is the orientation order parameter in [0, 1].
    """
    idx = {n: i for i, n in enumerate(BAND_NAMES)}
    shape = stack.shape[1:]
    strength = np.ones(shape, dtype=np.float64)
    cos2 = np.zeros(shape, dtype=np.float64)
    sin2 = np.zeros(shape, dtype=np.float64)
    wsum = np.zeros(shape, dtype=np.float64)
    for sigma in scales:
        for name in LOCATOR:
            b = stack[idx[name]].astype(np.float64)
            v = b[footprint]
            z = (b - v.mean()) / (v.std() + 1e-12)
            resp, theta = _line_response(z, sigma)
            r = np.abs(resp)
            ref = np.percentile(r[footprint], 99.0) if footprint.any() else 1.0
            rn = np.clip(r / (ref + 1e-30), 0.0, 1.0)
            strength *= np.maximum(rn, 1e-6)
            cos2 += rn * np.cos(2 * theta)
            sin2 += rn * np.sin(2 * theta)
            wsum += rn
    n_terms = len(scales) * len(LOCATOR)
    strength = np.power(strength, 1.0 / n_terms)
    R = np.hypot(cos2, sin2) / (wsum + 1e-30)
    return strength, R


def consensus_score(stack: np.ndarray, footprint: np.ndarray) -> np.ndarray:
    """The field used for ranking: line strength weighted by azimuth agreement, in [0, 1]."""
    strength, R = consensus_field(stack, footprint)
    score = strength * R
    ref = np.percentile(score[footprint], 99.9) if footprint.any() else 1.0
    return np.clip(score / (ref + 1e-30), 0.0, 1.0).astype(np.float32)
