# dark-lens-dr4: pre-registration (frozen 2026-11-24, eight days before Gaia DR4)

MV2, Munim, Inc. (mv2@munimtech.com). Apache-2.0 (code), CC BY 4.0 (text).
Made with AI agents working in MV2's name. No human expert has checked this plan yet.

## 1. The near miss

Rybicki et al. (2024, arXiv:2407.13740, ApJ, doi:10.3847/1538-4357/ad7bb1, "Analysis of the full
Spitzer microlensing sample I: Dark remnant candidates and Gaia predictions") studied nine
long microlensing events seen by OGLE, KMTNet, MOA, Spitzer and Gaia. From the light curves and
the Spitzer parallax they found four dark-remnant candidates, with Bayesian lens masses (Kroupa
prior, their Table 4):

| event | lens mass (Msun) | thetaE (mas) | forecast for Gaia |
|---|---|---|---|
| OGLE-2016-BLG-0293 | 2.98 +1.75/-1.28 | 0.91 +0.52/-0.39 | detectable, delta chi2 about 35 (about 6 sigma) |
| OGLE-2018-BLG-0483 | 4.65 +3.12/-2.08 | 3.68 +2.42/-1.62 | detectable, delta chi2 about 24 (about 5 sigma) |
| OGLE-2018-BLG-0662 | 3.15 +0.66/-0.64 | 1.42 +0.25/-0.26 | unlikely to be detected |
| OGLE-2015-BLG-0149 | 1.40 +0.75/-0.56 | 0.77 +0.38/-0.29 | 1.5-3 sigma |

Their masses rest on Galactic-model priors, because the light curve does not measure the
angular Einstein radius thetaE. The paper says that the lenses "could be confirmed through
forthcoming release of the Gaia time-series astrometry", and that the first two "are expected
to exhibit astrometric microlensing signals detectable by Gaia". The light centroid of an
unresolved microlensing image pair moves by thetaE u / (u^2 + 2). Measuring that motion gives
thetaE, and with the paper's parallax piE the mass follows directly: M = thetaE / (kappa piE),
kappa = 8.144 mas/Msun.

Gaia DR4 (2 December 2026, 11:00 UTC) publishes per-transit epoch astrometry for nearly all
sources for the first time. This plan fixes in advance how we will measure thetaE for all nine
events and grade the forecast.

**Why this is open (checked 2026-11-24).** An arXiv API search (all metadata fields) and a web
search return only the 2024 paper itself for "OGLE-2016-BLG-0293", "OGLE-2018-BLG-0483" and
"Gaia18ayh". Wicker et al. (2026,
arXiv:2611.24561) forecast DR4 tests for 56 other Gaia Science Alerts events and explicitly
excluded the Rybicki et al. events from their sample. Kaczmarek et al. (2026, arXiv:2609.10705)
work on simulations. No paper yet uses DR4 data, because DR4 is not out.

**The result, if it works.** The first black-hole mass measured from Gaia's own detection of
bent starlight. The best-known isolated stellar black hole, OGLE-2011-BLG-0462, was weighed the
same way but with Hubble astrometry.

## 2. Data

- The nine events, their two photometric solutions (u0+ and u0-), the preferred solution, and
  the Table 4 predictions are frozen in `events_bh.py`. Gaia DR3 source_ids come from the
  paper's Table 3. DR3 positions, G and proper motions were queried on 2026-11-24 and are frozen
  in `events_bh.py` and `bh_job.py`. Every Table 3 position agrees to below 0.1 arcsec (a test
  checks this).
- **DR4 counterpart.** We propagate the DR3 position to J2017.5 and run a cone search of 0.5
  arcsec in `gaiadr4.gaia_source`, then in `gaiadr4.all_source_astrometry`. We take the nearest
  entry, preferring |G_DR4 - G_DR3| < 1. If no counterpart is found, DataLink is tried with the
  DR3 source_id.
- **Epoch astrometry.** We fetch it by DataLink with the downloader of gj3512-giants-dr4
  (`io_gaia.py` v0.2.0, `gaiaorb/epochs.py`), both unchanged. Each transit's AGIS-used AF CCDs
  are collapsed into one along-scan (AL) normal point.

## 3. Model

AL abscissa per transit:
`w = da* sin(psi) + dd cos(psi) + plx (f - w fr k) + mua* t sin(psi) + mud t cos(psi) + thetaE col(t)`.

**Trajectory** (`traj.py`). This follows Gould (2004)'s geocentric formalism exactly as
implemented in MulensModel:

- `tau = (t - t0)/tE + dN piEN + dE piEE` and `beta = u0 - dN piEE + dE piEN`.
- `delta_s = x_E(t0par) - x_E(t) + (t - t0par) v_E(t0par)`, projected on north (dN) and east
  (dE).
- `t0par = t0` of the solution.

**Orientation.** Lens and source each appear displaced by `plx x P(t)`, where
`P = -(O.east, O.north)` and O is the observer's position. So the source-minus-lens vector in
units of thetaE is `u_lin - piE P_nonlinear = u_lin - piE (dE, dN)`. Matching this to
(delta tau, delta beta) gives:

- `u_vec = tau e1 + beta e2`, with `e1 = -(piEE, piEN)/piE` and `e2 = (-piEN, piEE)/piE`
  (east, north components).
- `e1` is the direction of `-mu_rel`, as it must be for the source relative to the lens.

The tests check this against a direct forward simulation.

**Centroid shift.** `delta = thetaE u_vec/(u^2+2) x A/(A+g)`, where g is the paper's OGLE
blend ratio floored at 0. g is 0.00 +- 0.03 for OB160293 u0- and 0.02 +- 0.05 for OB180483
u0-.

**Fits** (`bhfit.fit_event`):

- **Primary (constrained).** One lens parameter, thetaE, with the orientation fixed by piE,
  as in the paper's own forecast. The model is linear in thetaE.
  - Jitter: the smallest added AL noise that makes chi2/dof <= 1.
  - Clipping: one 5-sigma clip, then the fit is repeated.
  - Null comparison: the null (5-parameter) fit uses the same points and weights, which gives
    delta chi2.
- **Secondary (free orientation).** Two linear lens parameters (A, B) = thetaE (cos phi,
  sin phi), fitted for both handednesses h = -1 and h = +1.
  - Direction check: CONSISTENT if the better h is -1 and phi lies within max(30 deg, 3 sigma)
    of the predicted angle; SIGN_FLIPPED if it is off by 180 deg; otherwise INCONSISTENT.
  - This tests the trajectory conventions on the data themselves.
- **Both solutions.** u0+ and u0- are both fitted. Their chi2 difference is reported, because
  the astrometry can break the photometric u0 degeneracy.

**MulensModel cross-check** (`xcheck.py`).

- **What it compares:** (tau, beta) from `traj.py` against MulensModel's Trajectory, every 5
  days through the DR4 window, for all 18 solutions.
- **Pre-registered switch:** if any |d tau| or |d beta| exceeds 0.01, the release-day job uses
  MulensModel's (tau, beta) instead, and records it.
- **Timing:** the check runs before DR4 as job dark-lens-xcheck-v1 if the queue allows, and
  again inside the release-day job.

## 4. Fixed constants and verdicts (bhfit.py)

**Usable event.** It needs at least 20 transits and at least 3 transits with u <= 2.
Otherwise the verdict is TOO_FEW or NO_EVENT_COVERAGE.

**Field-star null.** For each event we take up to 12 DR4 neighbours: the nearest stars within
60 arcsec (but beyond 3 arcsec) with |G - G_source| <= 0.75.

- Each neighbour's epoch astrometry is fitted with the event's own template.
- The empirical false-alarm probability is FAP = (1 + #{|z_null| >= z}) / (1 + N).
- With fewer than 8 null fits, FAP is undefined.
- The two headline events always get nulls. The other events get nulls only within the first
  4 hours of the job.

**Detection verdict** (preferred solution):

| condition | verdict |
|---|---|
| z = thetaE/sigma >= 3 and FAP <= 0.10 | DETECTED |
| z >= 3 and FAP undefined | DETECTED_UNCALIBRATED |
| z >= 3 and FAP > 0.10 | NOT_SIGNIFICANT_VS_FIELD |
| z <= -3 | WRONG_SIGN (a convention or systematics flag, never a detection) |
| otherwise | NOT_DETECTED |

**Mass.** We draw 400 Monte Carlo sets of (t0, u0, tE, piEN, piEE, g) from the paper's
errors: independent, asymmetric Gaussians. For each set we refit thetaE and add its
statistical error. The mass is M = thetaE / (kappa |piE|), and we report the 2.5, 16, 50, 84
and 97.5 percentiles. The total thetaE error is the fit error and the Monte Carlo spread
combined in quadrature.

**Mass class** (detections only):

| condition | class |
|---|---|
| M_2.5% >= 2.5 Msun | BLACK_HOLE_MASS_RANGE |
| median >= 2.5 Msun but M_2.5% < 2.5 | BLACK_HOLE_FAVOURED_NOT_SECURE |
| median < 2.5 Msun | BELOW_BLACK_HOLE_RANGE |

**Prediction check** (detections only). AS_PREDICTED if the measured thetaE lies inside
[med - 2 minus - 2 err, med + 2 plus + 2 err] of the paper's Kroupa-prior thetaE. Otherwise
BELOW_PREDICTION or ABOVE_PREDICTION.

**Non-detections.** The 95% upper limit is thetaE95 = max(0, thetaE) + 1.645 err_total, and
M95 = thetaE95 / (kappa piE). The data count as informative if err_total <= thetaE_pred/3.

| condition | label |
|---|---|
| informative and M95 < the paper's Kroupa M_L 2-sigma minimum | PREDICTION_FAILED |
| informative otherwise | NOT_DETECTED_CONSISTENT |
| not informative | UNINFORMATIVE |

**Gate G0.** The median over events of |model AL parallax factor - published| must be below
0.03 au (Gaia taken as 1.01 x Earth). If not, the overall outcome is GEOMETRY_FAIL and no
verdict is reported.

**Overall outcome.** `<n>_OF_2_DETECTED_<m>_BH_RANGE` over the two headline events,
OGLE-2016-BLG-0293 and OGLE-2018-BLG-0483. The other seven events are reported with the same
verdicts as secondary results.

## 5. Robustness (reported, never used in verdicts)

- All CCDs, including those AGIS did not use. AGIS may have down-weighted the deviating
  transits during the event.
- t0par = t0 - 30 d and t0 + 30 d. The paper does not state its t0par.
- The u0+ and u0- fits of every event.

## 6. Known limits, stated before the data

- **Blending.** Blending in Gaia's G band and window can differ from OGLE's I band. A luminous
  blend or lens dilutes the shift, so a measured thetaE would then be a lower limit.
- **Parameter correlations.** The paper's parameter errors are treated as independent. This
  probably overstates the Monte Carlo spread, which makes the test conservative.
- **Paper version.** The photometric parameters are the arXiv v1 values. The published ApJ
  version could differ slightly; it could not be read (HTTP 403).
- **Faint source.** OGLE-2018-BLG-0483 is G = 20.0 in DR3. Its per-transit AL errors may be
  several mas.
- **Mass boundary.** 2.5 Msun is a convention for the neutron-star maximum. Masses near it are
  reported as such, not as black holes.

## 7. Deviations

Any change after 2026-12-02 11:00 UTC goes in a dated DEVIATION file. It names what changed and
why, and gives the result under both the frozen rule and the changed rule.

## 8. Credits

- The forecast, photometric solutions and predicted masses are from Rybicki et al. (2024).
- The formalism is from Gould (2004) and MulensModel (Poleski & Yee 2019).
- The Gaia epoch-astrometry reader is from MV2's gj3512-giants-dr4.
- Data: ESA Gaia DR3 and DR4 (Gaia Collaboration).
