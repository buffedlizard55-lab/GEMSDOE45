"""Fault-catalogue segmentation, tip detection and local strike estimation.

Input is the provided USGS/INGENIOUS raster ``labels.tif`` (bytes 425,830, sha256
7ba308ccdc4418b31a178f4f1ef21aaa6e152e4028f2f6f64b01f7eb25ae4093, 60,988 positive pixels).
Output is a list of tips, each carrying its pixel position, its outward unit direction and the
length of the segment it belongs to. Everything downstream (``survival.py``) consumes that list.

Design notes
------------
* Segments are the 8-connected components of the catalogue mask. A component is treated as one
  "fault" for tip purposes; this over-merges faults that touch, which is stated as a limitation.
* Tips are found by skeletonising the component (``skimage``-free implementation: an iterative
  hit-or-miss thinning) and taking skeleton pixels with exactly one 8-connected neighbour.
  If a component is a single pixel, that pixel is its own tip with an undefined direction and is
  dropped.
* The outward direction is the unit vector from the skeleton centroid of the last ``TAIL`` pixels
  to the tip itself, i.e. the chord of the terminal segment rather than a single-edge difference.
  This is deliberately robust: a one-pixel difference on a stair-cased raster can rotate by 45
  degrees, a 12-pixel chord cannot.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage

TAIL = 12          # pixels of skeleton used to estimate the terminal strike
MIN_SEG_PX = 12    # segments shorter than this are not treated as faults


@dataclass(frozen=True)
class Tip:
    row: int
    col: int
    urow: float
    ucol: float      # outward unit direction along strike
    seg_len: int     # pixels in the source segment
    nbr: int         # index of the first neighbouring tip found along the ray (filled later)


def _thin(mask: np.ndarray, max_iter: int = 60) -> np.ndarray:
    """Morphological thinning (Zhang-Suen style, vectorised) of a binary mask."""
    img = mask.copy()
    for _ in range(max_iter):
        changed = False
        for phase in (0, 1):
            p = np.pad(img, 1)
            P2 = p[0:-2, 1:-1]
            P3 = p[0:-2, 2:]
            P4 = p[1:-1, 2:]
            P5 = p[2:, 2:]
            P6 = p[2:, 1:-1]
            P7 = p[2:, 0:-2]
            P8 = p[1:-1, 0:-2]
            P9 = p[0:-2, 0:-2]
            nb = (P2.astype(np.uint8) + P3 + P4 + P5 + P6 + P7 + P8 + P9)
            cond_a = (nb >= 2) & (nb <= 6)
            seq = [P2, P3, P4, P5, P6, P7, P8, P9, P2]
            trans = np.zeros_like(img, dtype=np.uint8)
            for a, b in zip(seq[:-1], seq[1:]):
                trans += ((~a.astype(bool)) & b.astype(bool)).astype(np.uint8)
            cond_b = (trans == 1)
            if phase == 0:
                cond_c = ~(P2 & P4 & P6)
                cond_d = ~(P4 & P6 & P8)
            else:
                cond_c = ~(P2 & P4 & P8)
                cond_d = ~(P2 & P6 & P8)
            kill = img & cond_a & cond_b & cond_c & cond_d
            if kill.any():
                img = img & ~kill
                changed = True
        if not changed:
            break
    return img


def _endpoints(skel: np.ndarray) -> list[tuple[int, int]]:
    """Skeleton pixels with exactly one 8-connected neighbour."""
    k = np.ones((3, 3), dtype=np.uint8)
    k[1, 1] = 0
    nb = ndimage.convolve(skel.astype(np.uint8), k, mode="constant")
    ep = skel & (nb == 1)
    return list(zip(*np.nonzero(ep)))


def _order_from_tip(skel: np.ndarray, tip: tuple[int, int], limit: int) -> list[tuple[int, int]]:
    """Walk the skeleton from ``tip`` for up to ``limit`` pixels (8-connected)."""
    pts = [tip]
    seen = {tip}
    cur = tip
    while len(pts) < limit:
        r, c = cur
        nxt = None
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                if dr == 0 and dc == 0:
                    continue
                q = (r + dr, c + dc)
                if 0 <= q[0] < skel.shape[0] and 0 <= q[1] < skel.shape[1] \
                        and skel[q] and q not in seen:
                    nxt = q
                    break
            if nxt:
                break
        if nxt is None:
            break
        pts.append(nxt)
        seen.add(nxt)
        cur = nxt
    return pts


def extract_tips(known: np.ndarray, tail: int = TAIL, min_seg: int = MIN_SEG_PX) -> list[Tip]:
    """Return every tip of every catalogue segment, with an outward unit strike direction."""
    lab, n = ndimage.label(known, structure=np.ones((3, 3), dtype=np.uint8))
    tips: list[Tip] = []
    for i in range(1, n + 1):
        comp = lab == i
        npx = int(comp.sum())
        if npx < min_seg:
            continue
        skel = _thin(comp)
        if not skel.any():
            skel = comp
        eps = _endpoints(skel)
        if not eps:                      # closed loop: no tips
            continue
        for ep in eps:
            chain = _order_from_tip(skel, ep, tail)
            if len(chain) < 3:
                continue
            far = np.array(chain[-1], dtype=float)
            here = np.array(ep, dtype=float)
            d = here - far
            norm = float(np.hypot(*d))
            if norm < 1e-9:
                continue
            tips.append(Tip(row=int(ep[0]), col=int(ep[1]),
                            urow=float(d[0] / norm), ucol=float(d[1] / norm),
                            seg_len=npx, nbr=-1))
    return tips


def catalogue_stats(known: np.ndarray) -> dict:
    lab, n = ndimage.label(known, structure=np.ones((3, 3), dtype=np.uint8))
    sizes = np.bincount(lab.ravel())[1:]
    skel = _thin(known) if known.sum() < 200000 else None
    out = {
        "fault_px": int(known.sum()),
        "components": int(n),
        "components_ge_12px": int((sizes >= MIN_SEG_PX).sum()),
        "largest_component_px": int(sizes.max()) if len(sizes) else 0,
        "component_size_median": float(np.median(sizes)) if len(sizes) else 0.0,
    }
    if skel is not None:
        out["skeleton_px"] = int(skel.sum())
        out["endpoints"] = len(_endpoints(skel))
    return out
