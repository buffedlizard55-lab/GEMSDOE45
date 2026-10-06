"""Kaplan-Meier survival analysis of fault-tip extension length.

THE QUESTION THIS ANSWERS
-------------------------
"How far does a mapped fault continue beyond the last pixel a mapper drew?"  That is a
time-to-event question with a well-defined sample, a well-defined event and a well-defined
censor, so it is the textbook setting for the nonparametric product-limit estimator of
Kaplan & Meier (1958, JASA 53(282), 457-481, doi:10.1080/01621459.1958.10501452).

THE SAMPLE
----------
The provided catalogue (``labels.tif``) is not a set of long traces: it is 3,199 8-connected
segments with a median of 12 pixels (1.2 km) -- i.e. mapped *sections* of faults.  In the Basin
and Range, successive sections of one structure are arranged en echelon along a common trend and
separated by relay gaps.  For every segment tip we therefore measure

    G_i = the along-strike distance from tip i to the nearest pixel of a DIFFERENT catalogue
          segment, searched inside a narrow cone about the tip's own outward strike direction.

G_i is an observed event for every tip whose along-trend neighbour lies inside the search window,
and it is RIGHT-CENSORED at the window for every tip whose neighbour does not.  The censored tips
are precisely the ones with the largest gaps, so the naive estimator (the mean of the observed
gaps) is biased LOW.  The product-limit estimator is the standard unbiased correction, which is the
whole reason for using it here rather than a fixed search radius.

THE EXTENSION LAW
-----------------
Relay-linkage geometry puts the tip-to-tip meeting point at approximately the midpoint of the gap
(standard en-echelon/relay-ramp geometry; the same midpoint assumption used in fault-interaction
studies), so the extension of a mapped tip is taken as ``E_i = c_link * G_i`` with ``c_link = 0.5``.
Each tip draws its own G_i from the FITTED survival curve by inverse-transform sampling with a
fixed seed, so the emitted extension length is a random variable per fault and is reported as one.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

C_LINK = 0.5          # relay-midpoint linkage fraction (documented assumption, not a fitted value)
S_MAX_PX = 120        # search/censoring window along strike, pixels (12 km)
CONE_HALFWIDTH_PX = 3.0
CONE_MAX_HALFWIDTH_PX = 4.0
TAN_CONE = 0.27       # ~15 degrees


@dataclass
class KMCurve:
    """Nonparametric product-limit estimate of S(d) = P(gap > d)."""
    times: np.ndarray
    survival: np.ndarray
    at_risk: np.ndarray
    events: np.ndarray
    censored: np.ndarray
    greenwood_var: np.ndarray
    n: int
    n_events: int
    n_censored: int
    s_max: float

    def __call__(self, d: np.ndarray | float) -> np.ndarray:
        d = np.asarray(d, dtype=float)
        idx = np.searchsorted(self.times, d, side="right") - 1
        out = np.where(idx < 0, 1.0, self.survival[np.clip(idx, 0, len(self.survival) - 1)])
        return out

    def quantile(self, q: float) -> float:
        """Smallest d with S(d) <= 1 - q (i.e. the q-quantile of the fitted gap distribution)."""
        target = 1.0 - q
        below = np.nonzero(self.survival <= target)[0]
        return float(self.times[below[0]]) if len(below) else float(self.s_max)

    def mean(self) -> float:
        """Restricted mean up to the largest observed time."""
        if len(self.times) == 0:
            return float(self.s_max)
        t = np.concatenate([[0.0], self.times])
        s = np.concatenate([[1.0], self.survival])
        return float(np.sum(np.diff(t) * s[:-1]))

    def sample(self, size: int, rng: np.random.Generator) -> np.ndarray:
        """Inverse-transform draws: S-free since S is piecewise constant and decreasing."""
        u = rng.random(size)
        # S decreases from 1 to S(t_max); draw by inverting the step function
        s_rev = self.survival
        idx = np.searchsorted(-s_rev, -u, side="left")
        drawn = np.where(idx >= len(self.times), float(self.s_max), self.times[np.clip(idx, 0, len(self.times) - 1)])
        # points below the smallest event time have S = 1 and should map to their own small gap;
        done = u > 1.0 - 1e-12
        drawn = np.where(done, self.times[0] if len(self.times) else float(self.s_max), drawn)
        return drawn

    def as_dict(self) -> dict:
        return {
            "n_observations": self.n,
            "n_events": self.n_events,
            "n_censored": self.n_censored,
            "s_max_px": self.s_max,
            "median_gap_px": self.quantile(0.5),
            "q25_gap_px": self.quantile(0.25),
            "q75_gap_px": self.quantile(0.75),
            "restricted_mean_gap_px": self.mean(),
            "survival_at": {str(int(d)): float(self(np.array([d]))[0]) for d in (1, 2, 3, 5, 8, 12, 20, 40)},
            "max_greenwood_se": float(np.sqrt(np.nanmax(self.greenwood_var))) if len(self.greenwood_var) else 0.0,
        }


def kaplan_meier(times: np.ndarray, events: np.ndarray, s_max: float = float(S_MAX_PX)) -> KMCurve:
    """Product-limit estimator. ``events[i] == 1`` for an observed gap, ``0`` for right-censored.

    Standard form with the convention that deaths and censorings at the same time are handled by
    removing censorings after the events at that time (Kaplan & Meier 1958, eq. 2b).
    """
    times = np.asarray(times, dtype=float)
    events = np.asarray(events, dtype=int)
    order = np.argsort(times)
    times, events = times[order], events[order]
    uniq = np.unique(times)
    surv, var, risk, ev, ce = [], [], [], [], []
    s_prev = 1.0
    v_prev = 0.0
    n = len(times)
    for t in uniq:
        at_risk = int((times >= t).sum())
        d = int(((times == t) & (events == 1)).sum())
        c = int(((times == t) & (events == 0)).sum())
        if at_risk == 0:
            continue
        s_prev = s_prev * (1.0 - d / at_risk) if d > 0 else s_prev
        if d > 0 and at_risk > d:
            v_prev = v_prev + d / (at_risk * (at_risk - d))
        surv.append(s_prev)
        var.append(v_prev * s_prev ** 2)
        risk.append(at_risk)
        ev.append(d)
        ce.append(c)
    return KMCurve(
        times=np.asarray(uniq[:len(surv)], dtype=float),
        survival=np.asarray(surv, dtype=float),
        at_risk=np.asarray(risk, dtype=int),
        events=np.asarray(ev, dtype=int),
        censored=np.asarray(ce, dtype=int),
        greenwood_var=np.asarray(var, dtype=float),
        n=n, n_events=int((events == 1).sum()), n_censored=int((events == 0).sum()),
        s_max=s_max,
    )


def measure_gaps(tips: list, known: np.ndarray, labels: np.ndarray,
                 s_max: int = S_MAX_PX) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """March along each tip's outward strike and find the nearest pixel of another segment.

    Returns ``(gap_px, event, first_hit_index)``.  ``event = 1`` when a different segment is
    reached inside ``s_max``; ``event = 0`` marks a right-censored observation at ``s_max``.
    """
    h, w = known.shape
    gaps = np.full(len(tips), float(s_max), dtype=float)
    events = np.zeros(len(tips), dtype=int)
    hits = np.full(len(tips), -1, dtype=int)
    tip_lab = np.array([int(labels[max(0, min(h - 1, t.row)), max(0, min(w - 1, t.col))])
                        for t in tips], dtype=int)
    for i, t in enumerate(tips):
        own = tip_lab[i]
        ur, uc = t.urow, t.ucol
        # unit vector perpendicular to the strike, for the lateral cone half-width
        pr, pc = -uc, ur
        found = -1
        for s in range(1, s_max + 1):
            half = min(CONE_MAX_HALFWIDTH_PX, CONE_HALFWIDTH_PX + TAN_CONE * s)
            for side in np.arange(-half, half + 1e-9, 1.0):
                r = int(round(t.row + s * ur + side * pr))
                c = int(round(t.col + s * uc + side * pc))
                if 0 <= r < h and 0 <= c < w and known[r, c] and labels[r, c] != own:
                    found = s
                    break
            if found > 0:
                break
        if found > 0:
            gaps[i] = float(found)
            events[i] = 1
    return gaps, events, hits
