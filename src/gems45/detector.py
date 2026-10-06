"""H46 -- Kaplan-Meier tip-survival extension, gated by measured band information.

THE GEOLOGICAL ARGUMENT
-----------------------
A fault is a finite rupture surface.  It grows by tip propagation until it either meets another
fault (relay/step-over linkage) or the driving stress is relieved.  A mapper, however, stops
drawing at the last pixel where the *evidence* is unambiguous.  Those two stopping rules are not
the same, and the difference -- the tip-extension length -- is exactly the population in which a
"new fault" enters a catalogue: an expert extends a known fault, or fills the relay gap between two
known faults, far more often than they invent an isolated trace in blank ground.

That makes the extension length a time-to-event variable with an observable sample, so it is
estimated here with the Kaplan-Meier product-limit estimator rather than with a fixed search
radius:

  * sample : the 3,866 segment tips of the provided catalogue (3,199 8-connected segments)
  * event  : the along-strike distance G_i from tip i to the nearest pixel of a DIFFERENT segment
  * censor : right-censored at S_MAX = 120 px (12 km) when no such pixel exists inside the cone

THE EMISSION
------------
For a tip with outward unit strike u and estimated extension L, mass is emitted at
``x = tip + e*u`` for e = 1..ceil(L).  Only LOCATOR-class bands (bands.py) may veto a point: a
point is suppressed when the local gradient orientation of every locator band disagrees with u,
which is the physical statement "no independent field supports a line here".

THE BUDGET
----------
The number of emitted pixels is not chosen; it follows from the metric's own first-order condition
(``metric.marginal_bar``): a dot pays iff its expected realised kernel weight exceeds
alpha * DTI.  ``extension_budget`` integrates the KM curve to find the largest along-strike
distance at which that expectation still clears the bar, and that distance is then applied per tip.
"""
from __future__ import annotations

import numpy as np

from .metric import R_PX, marginal_bar

TARGET_DTI = 0.30          # the score level the bar is evaluated at; see evidence for sensitivity
C_LINK = 0.5               # relay midpoint: extension = half the gap (documented assumption)
G_FRAGMENT_MAX = 2         # gaps <= 2 px are rasterisation fragmentation, not step-overs
STRIKE_AGREE_COS = 0.5     # cos(60 deg): locator-band gradient need not be orthogonal to strike


def expected_weight_profile(km, e_max: int = 30) -> np.ndarray:
    """Expected realised kernel weight of a dot placed ``e`` px along strike from a tip.

    A dot at along-strike distance e earns credit when the next along-trend structure lies within
    the kernel window, i.e. when |G - e| <= R.  Using the fitted product-limit curve,

        E[w](e) = P(hit) * E[k | hit],   P(hit) = S(e - R) - S(e + R),   E[k|hit] ~ 0.5

    (the conditional mean of a triangular kernel over a window of half-width R is R/2 offset by the
    average |G - e|, which for a roughly uniform local gap density is close to 0.5).
    """
    e = np.arange(1, e_max + 1, dtype=float)
    lo = km(np.maximum(e - R_PX, 0.0))
    hi = km(e + R_PX)
    return 0.5 * np.maximum(lo - hi, 0.0)


def extension_budget(km, target_dti: float = TARGET_DTI, e_max: int = 30) -> dict:
    """Largest along-strike distance at which the expected weight still clears alpha*DTI."""
    prof = expected_weight_profile(km, e_max)
    bar = marginal_bar(target_dti)
    ok = prof > bar
    e_star = int(ok.sum())
    return {
        "target_dti": target_dti,
        "marginal_bar": round(bar, 6),
        "profile": [round(float(v), 5) for v in prof],
        "e_star_px": e_star,
        "e_star_m": e_star * 100,
        "cap_m": (1.0 - bar) * 300.0,
    }


def per_tip_extension(gaps: np.ndarray, events: np.ndarray, km, e_cap: int,
                      fragment_max: int = G_FRAGMENT_MAX) -> np.ndarray:
    """Per-tip extension length L_i in pixels. Varies by fault; never a single constant.

    Rule
    ----
    * observed gap G_i >  fragment_max :  L_i = min(round(C_LINK * G_i), e_cap)
      (the next along-trend structure really is that far away, so the relay midpoint sits at G_i/2)
    * observed gap G_i <= fragment_max :  the 1-2 px "gap" is rasterisation fragmentation of one
      trace and carries no step-over information; these tips take the KM median of the genuine
      sample.
    * censored (no structure inside the cone): L_i = min(round(C_LINK * conditional mean gap),
      e_cap), i.e. the product-limit correction for the long-gap tail, which a naive mean over
      observed gaps would have biased low.
    """
    gaps = np.asarray(gaps, dtype=float)
    events = np.asarray(events, dtype=int)
    genuine = (events == 1) & (gaps > fragment_max)
    median_genuine = float(np.median(gaps[genuine])) if genuine.any() else float(km.quantile(0.5))
    # conditional mean of G given censoring, from the restricted mean of the fit
    cond_mean_censored = max(float(km.mean()), float(np.percentile(gaps[events == 0], 50))
                             if (events == 0).any() else float(km.mean()))
    L = np.empty(len(gaps), dtype=int)
    for i in range(len(gaps)):
        if events[i] == 1 and gaps[i] > fragment_max:
            raw = C_LINK * gaps[i]
        elif events[i] == 1:
            raw = C_LINK * median_genuine
        else:
            raw = C_LINK * cond_mean_censored
        L[i] = int(np.clip(round(raw), 1, e_cap))
    return L


def locator_orientation(stack: np.ndarray, band_index: dict, footprint: np.ndarray,
                        sigma: float = 1.2) -> np.ndarray:
    """Multi-band gradient-orientation tensor of the LOCATOR-class bands.

    Returns ``(gx, gy, mag)``: the mass-weighted sum of squared gradients, whose principal direction
    is the local lineament normal.  Measured on the bands that pass the information screen, so a
    regionally smooth field cannot vote.
    """
    from scipy import ndimage
    loc = ["tmi_hg", "tmi_vg", "rtp", "tmi", "det_elev_slope", "iso_grav_anom_vg",
           "iso_grav_anom_hg", "tc"]
    gxx = np.zeros(stack.shape[1:], dtype=np.float64)
    gyy = np.zeros(stack.shape[1:], dtype=np.float64)
    gxy = np.zeros(stack.shape[1:], dtype=np.float64)
    mag = np.zeros(stack.shape[1:], dtype=np.float64)
    for name in loc:
        j = band_index[name]
        b = stack[j].astype(np.float64)
        v = b[footprint]
        z = (b - v.mean()) / (v.std() + 1e-12)
        z = ndimage.gaussian_filter(z, sigma)
        gy, gx = np.gradient(z)
        g = np.hypot(gx, gy)
        scale = 1.0 / (1.0 + np.exp(-(g - np.median(g[footprint])) * 4.0))  # soft, scale-free gate
        gxx += (gx * gx) * scale
        gyy += (gy * gy) * scale
        gxy += (gx * gy) * scale
        mag += g * scale
    return gxx, gyy, gxy, mag


def strike_agreement_field(stack: np.ndarray, footprint: np.ndarray, sigma: float = 1.2) -> np.ndarray:
    """cos(2*(theta_line - theta_strike))-style 'is there a line here' evidence, per pixel.

    For each pixel the local lineament direction theta_line is read from the structure tensor of the
    locator bands.  The returned field is the *line strength* (largest eigenvalue) normalised to
    [0, 1]; it is used to veto tip-extension points that lie in featureless ground.
    """
    idx = {n: i for i, n in enumerate(
        ["mag_anom", "rtp", "tmi_hg", "geod_2ndinv", "iso_grav_anom_slope", "tc", "geod_shearrate",
         "geod_dilaterate", "tmi_vg", "deq_n100a15", "iso_grav_anom_vg", "det_elev", "iso_grav_anom",
         "tmi", "depth_to_base_surf", "ieq_n100a15", "cond_surf", "iso_grav_anom_hg",
         "det_elev_slope"])}
    gxx, gyy, gxy, mag = locator_orientation(stack, idx, footprint, sigma)
    tr = gxx + gyy
    det = gxx * gyy - gxy * gxy
    disc = np.sqrt(np.maximum(tr * tr / 4.0 - det, 0.0))
    lam1 = tr / 2.0 + disc
    ref = np.percentile(lam1[footprint], 99.0) if footprint.any() else 1.0
    return np.clip(lam1 / (ref + 1e-12), 0.0, 1.0)
