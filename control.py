"""Known-answer control for dark-lens-dr4: the microlensing event Gaia19bld, whose angular
Einstein radius was measured directly by VLTI/PIONIER interferometry. Standard library only.

Part of dark-lens-dr4 (MV2, Munim, Inc., mv2@munimtech.com). Apache-2.0.
Made with AI agents working in MV2's name; no human expert has checked it yet.
Frozen with AMENDMENT-2026-11-24-control.md on 2026-11-24, before Gaia DR4 (2026-12-02).

Why: the headline test weighs dark lenses whose true thetaE nobody knows. Gaia19bld is the one
Gaia-era single-lens event whose thetaE is known to 0.5% from interferometry, and whose
discovery paper forecast that Gaia DR4 "may" detect its astrometric microlensing signal. Running
the same fit on it first tells us whether the pipeline recovers a known answer from DR4 data.

Sources of every number:
- Photometric solutions: Rybicki et al. (2022, A&A 657, A18; arXiv:2112.01613v1), Table 5,
  ground + Gaia + Spitzer fits (t0 in HJD; t0,par = 2458680.0 for all solutions). The paper
  adopts u0 < 0 (preferred by delta chi2 about 50 and by the interferometry).
- rho = 0.03202 (same table, u0 < 0); finite-source thetaE = 0.755 +- 0.013 mas (same paper).
- Interferometric thetaE = 0.7650 +- 0.0038 mas and M = 1.147 +- 0.029 Msun: Cassan et al.
  (2022, Nature Astronomy 6, 121, "Microlensing mass measurement from images of rotating
  gravitational arcs"), as quoted by Rybicki et al. (2022) and Bachelet et al. (2022,
  arXiv:2112.01676: "0.765 +- 0.004 mas").
- Blend: OGLE f_s = 1.029 +- 0.016 and the Spitzer blend flux is consistent with 0, so the
  G-band blend ratio is taken as g = 0.00 +- 0.03 (an assumption, stated in the amendment).
- Gaia DR3 source 5860179620698137600 (G = 14.797, RUWE 0.99): queried 2026-11-24 from
  gaiadr3.gaia_source, 0.05 arcsec from the paper's position (12:37:32.56, -66:06:40.90).
"""
import math

import bhfit as B
import traj as T

VERSION = "0.1.0"

PIONIER_THETAE = (0.7650, 0.0038)     # mas, Cassan et al. 2022 (the truth for the control)
FINITE_SOURCE_THETAE = (0.755, 0.013)  # mas, Rybicki et al. 2022 (reported, not graded)
CASSAN_MASS = (1.147, 0.029)          # Msun
RHO = 0.03202
U_FINITE_SOURCE = 0.1                 # transits with u < 0.1 (about 3 rho) are excluded
U_BRIGHT = 0.2                        # robustness only: u < 0.2 means G < about 13 (window change)
PULL_MAX = 2.0


def _s(t0, t0e, u0, u0e, tE, tEe, piEN, piENe, piEE, piEEe, g, ge):
    return dict(t0=t0, t0_err=t0e, u0=u0, u0_err=u0e, tE=tE, tE_err=tEe, piEN=piEN, piEN_err=piENe,
                piEE=piEE, piEE_err=piEEe, g=g, g_err=ge)


GAIA19BLD = dict(
    short="GAIA19BLD", name="Gaia19bld", gsa="Gaia19bld", dr3=5860179620698137600,
    ra=189.3856349809736, dec=-66.11136870197268, g_dr3=14.797316, ruwe_dr3=0.98854834,
    pm=(-7.4322500840074355, 0.09420716296683145), paper_radec=(189.385667, -66.111361),
    t0par_jd=2458680.0, preferred="u0-",
    sols={"u0-": _s(8681.219, (0.001, 0.001), -0.0192, (0.0001, 0.0001), 107.06, (0.50, 0.50),
                    -0.0378, (0.0012, 0.0012), -0.0731, (0.0014, 0.0014), 0.0, (0.03, 0.03)),
          "u0+": _s(8681.208, (0.001, 0.001), 0.0190, (0.0001, 0.0001), 107.60, (0.46, 0.46),
                    -0.0415, (0.0010, 0.0010), -0.0757, (0.0012, 0.0012), 0.0, (0.03, 0.03))},
    mu_rel_geo_ne=(-1.18, -2.29),     # mas/yr (N, E), Rybicki et al. 2022
)


def select(points, sol, ev=GAIA19BLD, u_min=U_FINITE_SOURCE, tb_all=None):
    """Drop transits with u < u_min (finite-source regime). Returns (points, tb, n_dropped)."""
    tb_all = T.taubeta_points(points, sol, ev["ra"], ev["dec"], ev["t0par_jd"]) if tb_all is None else tb_all
    keep = [i for i, (a, b) in enumerate(tb_all) if math.hypot(a, b) >= u_min]
    return [points[i] for i in keep], [tb_all[i] for i in keep], len(points) - len(keep)


def fit_control(points, sol, ev=GAIA19BLD, u_min=U_FINITE_SOURCE, tb_all=None):
    pts, tb, n_drop = select(points, sol, ev, u_min, tb_all)
    fr = B.fit_event(pts, sol, ev["ra"], ev["dec"], t0par_jd=ev["t0par_jd"], tb=tb)
    fr["n_finite_source_excluded"] = n_drop
    return pts, fr


def control_verdict(fr, err_total, fap, truth=PIONIER_THETAE):
    """Pre-registered control grade (AMENDMENT-2026-11-24-control.md):
    pull = (thetaE - truth) / sqrt(err_total^2 + truth_err^2)
      fit not OK                         -> CONTROL_<status>
      detection WRONG_SIGN               -> CONTROL_WRONG_SIGN
      |pull| > 2                         -> CONTROL_BIASED_LOW / CONTROL_BIASED_HIGH
      |pull| <= 2 and DETECTED(_UNCAL.)  -> CONTROL_PASS
      |pull| <= 2, not detected, err_total > truth / 3 -> CONTROL_UNINFORMATIVE
      |pull| <= 2 otherwise              -> CONTROL_CONSISTENT_NOT_SIGNIFICANT"""
    det = B.detection_verdict(fr, fap)
    if fr.get("status") != "OK":
        return dict(detection=det, pull=None, label="CONTROL_" + str(fr.get("status", "FAILED")))
    th = fr["thetaE"]
    s = math.sqrt(err_total ** 2 + truth[1] ** 2)
    pull = (th - truth[0]) / s
    if det == "WRONG_SIGN":
        label = "CONTROL_WRONG_SIGN"
    elif abs(pull) > PULL_MAX:
        label = "CONTROL_BIASED_LOW" if pull < 0 else "CONTROL_BIASED_HIGH"
    elif det in ("DETECTED", "DETECTED_UNCALIBRATED"):
        label = "CONTROL_PASS"
    elif err_total > truth[0] / 3.0:
        label = "CONTROL_UNINFORMATIVE"
    else:
        label = "CONTROL_CONSISTENT_NOT_SIGNIFICANT"
    return dict(detection=det, pull=pull, ratio=th / truth[0], label=label)


def headline_flag(label):
    """How the control qualifies the headline black-hole results (never changes their verdicts)."""
    if label == "CONTROL_PASS":
        return "VALIDATED"
    if label in ("CONTROL_CONSISTENT_NOT_SIGNIFICANT", "CONTROL_UNINFORMATIVE"):
        return "VALIDATION_WEAK"
    return "NOT_VALIDATED"


def mass_from_thetaE(thetaE, sol):
    return thetaE / (B.KAPPA * math.hypot(sol["piEN"], sol["piEE"]))
