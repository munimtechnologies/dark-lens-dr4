# Amendment 2 (2026-11-24, before Gaia DR4): two more lenses Gaia can weigh, Gaia19dke and Gaia18cbf

MV2, Munim, Inc. (mv2@munimtech.com). CC BY 4.0. Made with AI agents working in MV2's name.
No human expert has checked this plan yet. It was frozen on 2026-11-24, eight days before Gaia DR4
(2026-12-02 11:00 UTC). It changes nothing in PREREGISTRATION.md or in Amendment 1. It adds two
events that are fitted with the same code and graded by the rules below.

## Why these two

Both are Gaia Science Alert microlensing events. Their discovery papers estimated the lens mass
only through a Galactic-model prior, and both said Gaia's own time-series astrometry would measure
it. Both were left out of the DR4 forecast of Wicker et al. (2026, arXiv:2611.24561) as already
studied. Searches of arXiv and the web on 2026-11-24 for "Gaia19dke" and "Gaia18cbf" found only the
discovery papers and citing work, with no measured astrometric signal. None can exist yet, because
DR4 publishes the first Gaia epoch astrometry.

- **Gaia19dke** (Cygnus, Galactic disc). The event was bright (G = 15.45 at baseline), unblended
  (fS,G = 0.995) and long, with a parallax so strong that the light curve shows several peaks.
  - Maskoliunas et al. (arXiv:2309.03324; refereed in Acta Astronomica 74, 77, 2024) estimate
    thetaE of about 0.9 mas. They write that the shift is "easily detectable in the Gaia
    astrometric data". Their lens is about 0.5 solar masses at about 3 kpc, which they propose is
    an isolated white dwarf or a faint main-sequence star.
  - **What makes it a real test:** the main peak was in mid-2020, after DR4's last epoch
    (2020-01-20). DR4 sees only the slow rise, when the shift grows from about 0.07 to about
    0.35 thetaE. We fixed in advance whether that partial record is enough.
- **Gaia18cbf** (southern Galactic plane). It is one of the longest events ever found, with
  tE = 491 days, and it peaked in February 2019, inside the DR4 window.
  - Kruszynska et al. (2022, A&A 662, A59; arXiv:2111.08337) conclude the lens is probably a dark
    remnant of about 2.65 (+5.09/-1.48) solar masses, with thetaE = 2.53 (+4.35/-1.27) mas. They
    write that astrometric microlensing "might be marginally measurable by Gaia".
  - The source is faint (G = 20.2), so per-transit errors are several mas. An informative result
    is not guaranteed, and we say so now.

## Inputs (extra_events.py)

| | Gaia19dke | Gaia18cbf |
|---|---|---|
| Gaia DR3 source | 2026409795566972544 (G 15.449, RUWE 0.98) | 5995134605563805056 (G 20.168, RUWE 0.97) |
| solution used | Gaia + follow-up, u0 < 0 (the only one) | GF+ (best chi2); GF- also fitted and reported |
| t0 (JD - 2450000), t0,par | 9064.0639, 2459068 | 8524.76, 2458529.00 |
| u0, tE (d) | -0.6100, 162.47 | 0.0825, 491.41 |
| piEN, piEE | -0.0911, -0.1923 | -0.1192, -0.0442 |
| G blend g = 1/fS - 1 | 0.005 (fS 0.9947) | 0.119 (fS 0.894) |
| paper thetaE (mas) | 0.90 (+5.37/-0.70) | 2.53 (+4.35/-1.27) |
| coverage gate | u <= 3 (peak after DR4) | u <= 2 (as headline) |

All parameter errors (minus, plus) are taken from the papers' tables and propagated by the same
400-draw Monte Carlo as the headline events.

**Source caveat for Gaia19dke.** The light-curve fit is from arXiv v1 Table 1. The refereed Acta
Astronomica version gives a different mass posterior in its abstract (0.50 +0.82/-0.27 Msun,
from changed priors). Its tables could not be read as text on 2026-11-24. The Gaia measurement of
thetaE does not depend on the prior. The mass we derive uses piE, which the multi-peak parallax
fixes to about 1.5%. If the refereed piE differs by more than 3 sigma from v1, the write-up reports
the mass for both values. The grade below uses v1.

## Fit

The fit is identical to the headline events: `bhfit.fit_event` with the orientation taken from
piE, jitter, one 5-sigma clip, at least 20 transits and at least 3 transits inside the coverage
gate. It is followed by `bhfit.monte_carlo` and a field-star null on up to 12 DR4 neighbours
(within 60 arcsec, |dG| <= 0.75). The trajectories come from traj.py, or from MulensModel if
bh_job.py switched (xcheck rule). The only per-event change is the coverage gate of Gaia19dke,
u <= 3 (`extra_events.fit_extra`), fixed now because its peak falls after DR4's last epoch.

## Grade (pre-registered, extra_events.grade)

Detection follows the headline rule: snr >= 3, and empirical FAP <= 0.10 when at least 8 nulls
exist (otherwise DETECTED_UNCALIBRATED).

- **Detected.** The label is DETECTED or DETECTED_UNCALIBRATED, plus a remnant class from the Monte
  Carlo mass percentiles:
  - BLACK_HOLE_MASS_RANGE: 2.5th percentile >= 2.5 Msun
  - ABOVE_WHITE_DWARF_LIMIT: 2.5th percentile > 1.4 Msun
  - WHITE_DWARF_OR_STAR_RANGE: 97.5th percentile <= 1.4 Msun
  - UNRESOLVED: otherwise

  The measured thetaE is also checked against the paper's (AS_PREDICTED / BELOW_PREDICTION /
  ABOVE_PREDICTION, the same widened-95% rule as the headline events).
- **Not detected** (including NOT_SIGNIFICANT_VS_FIELD). Here th95 = max(0, thetaE) + 1.645
  err_total and m95 = th95 / (8.144 piE).

| condition | label | meaning |
|---|---|---|
| snr <= -3 | WRONG_SIGN | pipeline or convention problem |
| m95 < 1.4 | UPPER_LIMIT_BELOW_1.4_MSUN | no neutron star or black hole at the paper's piE |
| m95 < 2.5 | UPPER_LIMIT_BELOW_2.5_MSUN | no black hole at the paper's piE |
| err_total > paper thetaE / 3 | UNINFORMATIVE | DR4 too weak for this lens |
| otherwise | NOT_DETECTED_CONSISTENT | informative, but consistent with the paper |

For Gaia18cbf, UPPER_LIMIT_BELOW_2.5_MSUN or below would contradict the paper's preferred
black-hole-like mass. For Gaia19dke, a WHITE_DWARF_OR_STAR_RANGE detection would be the first
astrometric-microlensing mass of a Gaia-alert lens from Gaia's own data.

**The control flag carries over.** `results_extra.json` records Amendment 1's flag (VALIDATED,
VALIDATION_WEAK or NOT_VALIDATED). Under NOT_VALIDATED, the write-up leads with the control
failure for these events too.

**Reported, never graded:** the other photometric solution of Gaia18cbf (GF-), a fit that drops
transits with u < 0.1 (finite source, unknown rho), and an all-CCD fit.

## Release-day job

`dr4day/job.template.json` now runs `python -m unittest -v test_darklens test_control test_extra`,
then `bh_job.py`, `control_job.py` and `extra_job.py`. The new output is
`/out/extra/results_extra.json` plus per-transit points. The timeout rises from 420 to 480 minutes.

## Credits

- Gaia18cbf: Kruszynska, Wyrzykowski, Rybicki, Maskoliunas et al. (2022).
- Gaia19dke: Maskoliunas, Wyrzykowski, Howil et al. (2023 preprint; Acta Astronomica 2024).
- DR4 forecasting context: Wicker et al. (2026).
- Data: ESA Gaia DR3 and DR4 (Gaia Collaboration).
