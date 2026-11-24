"""Fit of Gaia DR4 AL epoch astrometry: the source's 5 astrometric parameters plus the
microlensing centroid shift of a photometrically known event; Monte Carlo lens mass; the
pre-registered verdicts. Standard library only.

Part of dark-lens-dr4 (MV2, Munim, Inc., mv2@munimtech.com). Apache-2.0.
Made with AI agents working in MV2's name; no human expert has checked it yet.
Frozen with PREREGISTRATION.md on 2026-11-24, before Gaia DR4 (2026-12-02).

AL model per transit point (gaiaorb.epochs convention):
  w = da* s + dd c + plx (f - w fr k) + mua* t s + mud t c + thetaE col(t)
col = traj.column_constrained (orientation from piE, primary) or the two free-orientation
columns of traj.columns_free (secondary).
"""
import math
import random

import traj as T

VERSION = "0.1.0"
KAPPA = 8.144              # mas per solar mass
MAS_TO_RAD = math.pi / (180.0 * 3600.0 * 1000.0)

# Pre-registered constants (PREREGISTRATION.md section 4)
MIN_POINTS = 20            # usable transits
MIN_EVENT_POINTS = 3       # transits with u <= U_EVENT
U_EVENT = 2.0
CLIP = 5.0
SNR_DETECT = 3.0
FAP_MAX = 0.10
MIN_NULLS = 8
N_MC = 400
NS_MAX = 2.5               # Msun; lower mass bound above this -> black-hole mass range
DIRECTION_TOL_DEG = 30.0
SOURCE_PARAMS = ("ra_offset_mas", "dec_offset_mas", "parallax_mas", "pmra_masyr", "pmdec_masyr")


# ----------------------------------------------------------------------------- linear algebra
def invert(A):
    n = len(A)
    M = [list(map(float, A[i])) + [1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[piv][col]) < 1e-300:
            raise ValueError("singular normal matrix")
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        M[col] = [v / pv for v in M[col]]
        for r in range(n):
            if r != col and M[r][col] != 0.0:
                f = M[r][col]
                M[r] = [vr - f * vc for vr, vc in zip(M[r], M[col])]
    return [row[n:] for row in M]


def solve_wls(cols, y, sig):
    m, n = len(cols), len(y)
    A = [[0.0] * m for _ in range(m)]
    b = [0.0] * m
    for i in range(n):
        w = 1.0 / (sig[i] * sig[i])
        row = [cols[a][i] for a in range(m)]
        for a in range(m):
            ra = row[a] * w
            b[a] += ra * y[i]
            for c in range(a, m):
                A[a][c] += ra * row[c]
    for a in range(m):
        for c in range(a):
            A[a][c] = A[c][a]
    cov = invert(A)
    beta = [sum(cov[a][c] * b[c] for c in range(m)) for a in range(m)]
    resid = [y[i] - sum(beta[a] * cols[a][i] for a in range(m)) for i in range(n)]
    chi2 = sum((resid[i] / sig[i]) ** 2 for i in range(n))
    return dict(beta=beta, cov=cov, chi2=chi2, resid=resid)


def _sig(base, eps):
    return [math.sqrt(s * s + eps * eps) for s in base]


def fit_jitter(cols, y, base, n_par, iters=60):
    """Smallest eps >= 0 with chi2 / (N - n_par) <= 1."""
    target = len(y) - n_par
    if target <= 0:
        return 0.0

    def chi2(eps):
        return solve_wls(cols, y, _sig(base, eps))["chi2"]

    if chi2(0.0) <= target:
        return 0.0
    hi = max(1e-4, sorted(base)[len(base) // 2])
    k = 0
    while chi2(hi) > target:
        hi *= 2.0
        k += 1
        if k > 80:
            raise ValueError("jitter search diverged")
    lo = 0.0
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if chi2(mid) > target:
            lo = mid
        else:
            hi = mid
    return hi


def source_columns(points):
    k = MAS_TO_RAD
    return [[p["s"] for p in points],
            [p["c"] for p in points],
            [p["f"] - p["w"] * p.get("fr", 0.0) * k for p in points],
            [p["t"] * p["s"] for p in points],
            [p["t"] * p["c"] for p in points]]


def _sub(col, idx):
    return [col[i] for i in idx]


def percentile(xs, q):
    v = sorted(xs)
    if not v:
        return None
    pos = (len(v) - 1) * q / 100.0
    lo = int(math.floor(pos))
    hi = min(lo + 1, len(v) - 1)
    return v[lo] + (v[hi] - v[lo]) * (pos - lo)


# ----------------------------------------------------------------------------- fits
def fit_event(points, sol, ra, dec, t0par_jd=None, tb=None):
    """Primary (constrained) fit with jitter and one 5-sigma clip, then the null and the two
    free-orientation fits on the same points and weights."""
    if len(points) < MIN_POINTS:
        return dict(status="TOO_FEW", n_points=len(points))
    tb = T.taubeta_points(points, sol, ra, dec, t0par_jd) if tb is None else tb
    us = [math.hypot(a, b) for a, b in tb]
    n_ev = sum(1 for u in us if u <= U_EVENT)
    if n_ev < MIN_EVENT_POINTS:
        return dict(status="NO_EVENT_COVERAGE", n_points=len(points), n_event_points=n_ev,
                    u_min=min(us))
    src = source_columns(points)
    con = T.column_constrained(points, sol, tb)
    y_all = [p["w"] for p in points]
    base_all = [p["sig"] for p in points]
    idx = list(range(len(points)))
    n_clipped = 0
    for rnd in range(2):
        cols = [_sub(c, idx) for c in src] + [_sub(con, idx)]
        y = _sub(y_all, idx)
        base = _sub(base_all, idx)
        eps = fit_jitter(cols, y, base, 6)
        sig = _sig(base, eps)
        fit = solve_wls(cols, y, sig)
        if rnd == 1:
            break
        good = [idx[i] for i in range(len(idx)) if abs(fit["resid"][i]) <= CLIP * sig[i]]
        if len(good) == len(idx) or len(good) < MIN_POINTS:
            break
        n_clipped = len(idx) - len(good)
        idx = good
    th, sth = fit["beta"][5], math.sqrt(fit["cov"][5][5])
    y = _sub(y_all, idx)
    srcs = [_sub(c, idx) for c in src]
    null = solve_wls(srcs, y, sig)
    out = dict(status="OK", n_points=len(idx), n_total=len(points), n_clipped=n_clipped,
               n_event_points=sum(1 for i in idx if us[i] <= U_EVENT), u_min=min(us[i] for i in idx),
               jitter_mas=eps, chi2=fit["chi2"], chi2_null=null["chi2"],
               dchi2=null["chi2"] - fit["chi2"], thetaE=th, thetaE_err=sth,
               snr=th / sth if sth > 0 else 0.0,
               source={n: [fit["beta"][a], math.sqrt(fit["cov"][a][a])] for a, n in enumerate(SOURCE_PARAMS)},
               keep=idx, sig=sig)
    free = {}
    for h in (1, -1):
        ca, cb = T.columns_free(points, sol, tb, h)
        f = solve_wls(srcs + [_sub(ca, idx), _sub(cb, idx)], y, sig)
        A, B = f["beta"][5], f["beta"][6]
        vA, vB, cAB = f["cov"][5][5], f["cov"][6][6], f["cov"][5][6]
        r2 = A * A + B * B
        sphi = math.sqrt(max(0.0, B * B * vA + A * A * vB - 2.0 * A * B * cAB)) / r2 if r2 > 0 else None
        free[str(h)] = dict(thetaE=math.sqrt(r2), phi_deg=math.degrees(math.atan2(B, A)),
                            phi_err_deg=None if sphi is None else math.degrees(sphi),
                            chi2=f["chi2"], dchi2=null["chi2"] - f["chi2"])
    out["free"] = free
    out["direction"] = direction_check(free, T.predicted_phi(sol))
    return out


def angle_diff_deg(a, b):
    return abs((a - b + 180.0) % 360.0 - 180.0)


def direction_check(free, phi0_rad):
    """Pre-registered secondary check of the convention-free fit: the preferred handedness
    (lower chi2) should be h = -1 and its phi within max(30 deg, 3 sigma) of the prediction."""
    best = "-1" if free["-1"]["chi2"] <= free["1"]["chi2"] else "1"
    f = free[best]
    d = angle_diff_deg(f["phi_deg"], math.degrees(phi0_rad))
    tol = max(DIRECTION_TOL_DEG, 3.0 * (f["phi_err_deg"] or 0.0))
    if best == "-1" and d <= tol:
        label = "CONSISTENT"
    elif best == "-1" and angle_diff_deg(d, 180.0) <= tol:
        label = "SIGN_FLIPPED"
    else:
        label = "INCONSISTENT"
    return dict(best_h=int(best), phi_pred_deg=math.degrees(phi0_rad), phi_fit_deg=f["phi_deg"],
                diff_deg=d, tol_deg=tol, dchi2_h=abs(free["1"]["chi2"] - free["-1"]["chi2"]), label=label)


def draw(rng, val, err):
    z = rng.gauss(0.0, 1.0)
    lo, hi = err
    return val + z * (hi if z > 0 else lo)


def monte_carlo(points, fr, sol, ra, dec, t0par_jd=None, n=N_MC, seed=20261202):
    """Propagate the photometric-parameter errors (independent, asymmetric Gaussians) and the
    fit's statistical error into thetaE and the mass M = thetaE / (KAPPA piE)."""
    rng = random.Random(seed)
    t0par = sol["t0"] + T.HJD_OFFSET if t0par_jd is None else t0par_jd
    idx = fr["keep"]
    pts = [points[i] for i in idx]
    ne = T.ne_offsets(pts, ra, dec, t0par)
    src = source_columns(pts)
    y = [p["w"] for p in pts]
    sig = fr["sig"]
    th_hat, th_draw, masses = [], [], []
    for _ in range(n):
        s = dict(sol)
        for k in ("t0", "u0", "tE", "piEN", "piEE", "g"):
            s[k] = draw(rng, sol[k], sol[k + "_err"])
        if s["tE"] <= 0.0 or math.hypot(s["piEN"], s["piEE"]) <= 0.0:
            continue
        tb = [T.tau_beta_from_ne(T.jd_from_t(p["t"]), s, a, b) for p, (a, b) in zip(pts, ne)]
        col = T.column_constrained(pts, s, tb)
        f = solve_wls(src + [col], y, sig)
        th, sth = f["beta"][5], math.sqrt(f["cov"][5][5])
        thd = th + sth * rng.gauss(0.0, 1.0)
        th_hat.append(th)
        th_draw.append(thd)
        masses.append(thd / (KAPPA * math.hypot(s["piEN"], s["piEE"])))
    m = len(th_hat)
    mean = sum(th_hat) / m
    spread = math.sqrt(sum((x - mean) ** 2 for x in th_hat) / max(1, m - 1))
    return dict(n=m, thetaE_param_sd=spread,
                thetaE_total_err=math.sqrt(fr["thetaE_err"] ** 2 + spread ** 2),
                thetaE_p=[percentile(th_draw, q) for q in (2.5, 16, 50, 84, 97.5)],
                mass_p=[percentile(masses, q) for q in (2.5, 16, 50, 84, 97.5)])


# ----------------------------------------------------------------------------- verdicts
def empirical_fap(snr, null_snrs):
    """(1 + #{|z_null| >= snr}) / (1 + N); None if fewer than MIN_NULLS nulls."""
    if len(null_snrs) < MIN_NULLS:
        return None
    k = sum(1 for z in null_snrs if abs(z) >= snr)
    return (1.0 + k) / (1.0 + len(null_snrs))


def detection_verdict(fr, fap):
    if fr.get("status") != "OK":
        return fr.get("status", "FAILED")
    if fr["snr"] >= SNR_DETECT:
        if fap is None:
            return "DETECTED_UNCALIBRATED"
        return "DETECTED" if fap <= FAP_MAX else "NOT_SIGNIFICANT_VS_FIELD"
    if fr["snr"] <= -SNR_DETECT:
        return "WRONG_SIGN"
    return "NOT_DETECTED"


def mass_class(mass_p):
    """From the Monte Carlo mass percentiles (2.5, 16, 50, 84, 97.5)."""
    lo, med = mass_p[0], mass_p[2]
    if lo is not None and lo >= NS_MAX:
        return "BLACK_HOLE_MASS_RANGE"
    if med is not None and med >= NS_MAX:
        return "BLACK_HOLE_FAVOURED_NOT_SECURE"
    return "BELOW_BLACK_HOLE_RANGE"


def prediction_check(thetaE, thetaE_err, pred):
    """thetaE vs the paper's Table 4 thetaE (Kroupa prior): AS_PREDICTED when the measured
    value lies inside the prediction's 95% range widened by the measurement error
    (median - 2 minus - 2 err, median + 2 plus + 2 err)."""
    med, mi, pl = pred["thetaE"]
    lo = med - 2.0 * mi - 2.0 * thetaE_err
    hi = med + 2.0 * pl + 2.0 * thetaE_err
    return dict(range=[lo, hi], label="AS_PREDICTED" if lo <= thetaE <= hi else
                ("BELOW_PREDICTION" if thetaE < lo else "ABOVE_PREDICTION"))


def nondetection_verdict(thetaE, err_total, pred_kroupa_thetaE, piE, m2s_min):
    """When not detected: 95% upper limit and whether the data were informative."""
    th95 = max(0.0, thetaE) + 1.645 * err_total
    m95 = th95 / (KAPPA * piE)
    informative = err_total <= pred_kroupa_thetaE / 3.0
    if not informative:
        label = "UNINFORMATIVE"
    elif m95 < m2s_min:
        label = "PREDICTION_FAILED"
    else:
        label = "NOT_DETECTED_CONSISTENT"
    return dict(thetaE_upper95=th95, mass_upper95=m95, informative=informative, label=label)


def overall(g0_ok, headline):
    """headline: {short: dict(verdict, mass_class or None)}."""
    if not g0_ok:
        return "GEOMETRY_FAIL"
    n_det = sum(1 for v in headline.values() if v.get("verdict") in ("DETECTED", "DETECTED_UNCALIBRATED"))
    n_bh = sum(1 for v in headline.values() if v.get("mass_class") == "BLACK_HOLE_MASS_RANGE")
    return "%d_OF_%d_DETECTED_%d_BH_RANGE" % (n_det, len(headline), n_bh)
