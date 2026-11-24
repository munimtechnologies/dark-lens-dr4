"""Amendment 2 of dark-lens-dr4 (frozen 2026-11-24, before Gaia DR4 on 2026-12-02): two more
Gaia Science Alert microlensing events whose discovery papers say Gaia's own astrometry will
weigh the lens. Standard library only.

Part of dark-lens-dr4 (MV2, Munim, Inc., mv2@munimtech.com). Apache-2.0.
Made with AI agents working in MV2's name; no human expert has checked it yet.

Sources of every number:
- Gaia18cbf: Kruszynska et al. (2022, A&A 662, A59; arXiv:2111.08337v3), Table 4 (Gaia + LCO i
  fits; t0 in HJD - 2450000; t0,par = 2458529.00). Preferred: GF+ (lowest chi2). Blend: the
  paper's Gaia source fraction fS,G is converted to g = 1/fS - 1 (errors from fS -/+ its errors).
  Lens posteriors (M^-1.75 prior): GF+ M = 2.65 (+5.09/-1.48) Msun, thetaE = 2.53 (+4.35/-1.27)
  mas; the paper concludes the lens is probably a dark remnant.
- Gaia19dke: Maskoliunas et al. (arXiv:2309.03324v1), Table 1, Gaia + follow-up fit (t0 in
  JD - 2450000; t0,par = 2459068), the only solution (u0 < 0). Table 3 (G+F): thetaE = 0.90
  (+5.37/-0.70) mas, M = 0.51 (+3.07/-0.40) Msun. The refereed version (Acta Astronomica 74, 77,
  2024) gives M = 0.50 (+0.82/-0.27) Msun in its abstract and proposes a main-sequence star or an
  isolated white dwarf; its tables could not be read as text on 2026-11-24, so the light-curve
  fit is frozen from arXiv v1 (see AMENDMENT-2026-11-24-extra.md).
- Gaia DR3 rows (source_id, J2016.0 position, proper motion, G, RUWE): gaiadr3.gaia_source,
  queried 2026-11-24 (Gaia18cbf 0.04 arcsec, Gaia19dke 0.09 arcsec from the papers' positions).
"""
import math

import bhfit as B

VERSION = "0.1.0"
WD_MAX = 1.4        # Msun, Chandrasekhar mass: above it the lens cannot be a white dwarf
NS_MAX = B.NS_MAX   # Msun, 2.5: above it, black-hole mass range (as for the headline events)
DR4_JD = (2456863.5, 2458868.5)   # 2014-07-25 .. 2020-01-20


def _s(t0, t0e, u0, u0e, tE, tEe, piEN, piENe, piEE, piEEe, g, ge):
    return dict(t0=t0, t0_err=t0e, u0=u0, u0_err=u0e, tE=tE, tE_err=tEe, piEN=piEN, piEN_err=piENe,
                piEE=piEE, piEE_err=piEEe, g=g, g_err=ge)


def g_from_fs(fs, fs_minus, fs_plus):
    """Blend ratio g = 1/fS - 1 and its (minus, plus) errors from the fS range."""
    g = 1.0 / fs - 1.0
    g_hi = 1.0 / (fs - fs_minus) - 1.0
    g_lo = 1.0 / (fs + fs_plus) - 1.0
    return g, (g - g_lo, g_hi - g)


GAIA18CBF = dict(
    short="GAIA18CBF", name="Gaia18cbf", dr3=5995134605563805056,
    ra=241.1619124283208, dec=-41.104819248508576, g_dr3=20.167513, ruwe_dr3=0.9733,
    pm=(-1.8317763586104534, -1.8172896576499655), paper_radec=(241.161917, -41.1048306),
    t0par_jd=2458529.0, preferred="u0+", u_event=2.0,
    sols={"u0+": _s(8524.76, (4.57, 4.33), 0.0825, (0.0208, 0.0231), 491.41, (84.94, 128.31),
                    -0.1192, (0.0285, 0.0273), -0.0442, (0.0084, 0.0077), 0.1186, (0.2543, 0.3784)),
          "u0-": _s(8513.14, (4.68, 4.68), -0.0567, (0.0169, 0.0162), 453.74, (105.74, 178.69),
                    -0.1697, (0.0577, 0.0516), -0.0257, (0.0057, 0.0054), 0.3947, (0.4056, 0.6714))},
    fs=dict(u0p=(0.894, 0.226, 0.263), u0m=(0.717, 0.233, 0.294)),
    pred=dict(thetaE=(2.53, 1.27, 4.35), M=(2.65, 1.48, 5.09)),
    claim="dark remnant of about 2.65 Msun (Kruszynska et al. 2022)",
)

GAIA19DKE = dict(
    short="GAIA19DKE", name="Gaia19dke", dr3=2026409795566972544,
    ra=291.4945133633525, dec=28.406882528435283, g_dr3=15.448647, ruwe_dr3=0.97869515,
    pm=(-2.861862185544017, -5.447365851300708), paper_radec=(291.4945, 28.4068611),
    t0par_jd=2459068.0, preferred="u0-", u_event=3.0,
    sols={"u0-": _s(9064.0639, (0.33, 0.33), -0.6100, (0.0112, 0.0160), 162.47, (1.90, 2.68),
                    -0.0911, (0.0012, 0.0014), -0.1923, (0.0024, 0.0031), 0.00533, (0.02953, 0.04443))},
    fs=dict(u0m=(0.9947, 0.0421, 0.0301)),
    pred=dict(thetaE=(0.90, 0.70, 5.37), M=(0.51, 0.40, 3.07)),
    claim="isolated white dwarf or faint main-sequence star of about 0.5 Msun (Maskoliunas et al.)",
)

EVENTS = [GAIA19DKE, GAIA18CBF]   # brighter, better-constrained event first


def kappa_pie(sol):
    return B.KAPPA * math.hypot(sol["piEN"], sol["piEE"])


def mass_from_thetaE(thetaE, sol):
    return thetaE / kappa_pie(sol)


def remnant_class(mass_p):
    """From Monte Carlo mass percentiles (2.5, 16, 50, 84, 97.5) of a DETECTED lens."""
    lo, hi = mass_p[0], mass_p[4]
    if lo is not None and lo >= NS_MAX:
        return "BLACK_HOLE_MASS_RANGE"
    if lo is not None and lo > WD_MAX:
        return "ABOVE_WHITE_DWARF_LIMIT"
    if hi is not None and hi <= WD_MAX:
        return "WHITE_DWARF_OR_STAR_RANGE"
    return "UNRESOLVED"


def grade(fr, err_total, mass_p, fap, ev):
    """Pre-registered grade of one Amendment-2 event (AMENDMENT-2026-11-24-extra.md).
    Detected (snr >= 3; DETECTED or DETECTED_UNCALIBRATED): label = detection verdict, plus the
      remnant class from the Monte Carlo mass and the check against the paper's thetaE.
    Not detected: 95% upper limit th95 = max(0, thetaE) + 1.645 err_total, m95 = th95 / (kappa piE):
      WRONG_SIGN                 if snr <= -3
      UPPER_LIMIT_BELOW_1.4_MSUN if m95 < 1.4  (no neutron star or black hole at the paper's piE)
      UPPER_LIMIT_BELOW_2.5_MSUN if m95 < 2.5  (no black hole at the paper's piE)
      UNINFORMATIVE              if err_total > paper thetaE median / 3
      NOT_DETECTED_CONSISTENT    otherwise"""
    det = B.detection_verdict(fr, fap)
    out = dict(detection=det)
    if fr.get("status") != "OK":
        out["label"] = det
        return out
    sol = ev["sols"][ev["preferred"]]
    if det in ("DETECTED", "DETECTED_UNCALIBRATED"):
        out["label"] = det
        out["remnant_class"] = remnant_class(mass_p)
        out["prediction"] = B.prediction_check(fr["thetaE"], err_total, ev["pred"])
        return out
    th95 = max(0.0, fr["thetaE"]) + 1.645 * err_total
    m95 = th95 / kappa_pie(sol)
    informative = err_total <= ev["pred"]["thetaE"][0] / 3.0
    if det == "WRONG_SIGN":
        label = "WRONG_SIGN"
    elif m95 < WD_MAX:
        label = "UPPER_LIMIT_BELOW_1.4_MSUN"
    elif m95 < NS_MAX:
        label = "UPPER_LIMIT_BELOW_2.5_MSUN"
    elif not informative:
        label = "UNINFORMATIVE"
    else:
        label = "NOT_DETECTED_CONSISTENT"
    out.update(label=label, thetaE_upper95=th95, mass_upper95=m95, informative=informative)
    return out


def fit_extra(points, sol, ev, tb):
    """bhfit.fit_event with the event's own pre-registered coverage gate (u <= u_event).
    Gaia19dke peaked after DR4's last epoch, so its gate is u <= 3 instead of 2."""
    old = B.U_EVENT
    B.U_EVENT = ev.get("u_event", old)
    try:
        return B.fit_event(points, sol, ev["ra"], ev["dec"], t0par_jd=ev["t0par_jd"], tb=tb)
    finally:
        B.U_EVENT = old
