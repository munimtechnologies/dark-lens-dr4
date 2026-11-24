# Amendment 1 (2026-11-24, before Gaia DR4): a known-answer control, Gaia19bld

MV2, Munim, Inc. (mv2@munimtech.com). CC BY 4.0. Made with AI agents working in MV2's name.
No human expert has checked this plan yet. It was frozen on 2026-11-24, eight days before Gaia DR4
(2026-12-02 11:00 UTC), and changes nothing in PREREGISTRATION.md. It adds one event that is
fitted with the same code as the black-hole candidates and graded against a known answer.

## Why

The headline test weighs dark lenses whose angular Einstein radius thetaE nobody knows. A
detection of a 3 to 5 solar-mass dark lens would mean more if the same pipeline first recovered
a lens whose thetaE is already measured. One Gaia-era single-lens event fits that need:

- **Gaia19bld** (Gaia Science Alert, 2019; peak 2019-07-27, inside the DR4 window). Cassan et
  al. (2022, Nature Astronomy 6, 121) resolved its lensed images with VLTI/PIONIER and measured
  **thetaE = 0.7650 +- 0.0038 mas**, so M = 1.147 +- 0.029 Msun. Rybicki et al. (2022, A&A 657,
  A18; arXiv:2112.01613) found the same from finite-source photometry (0.755 +- 0.013 mas) and
  wrote that "it may be possible to detect the astrometric microlensing signal in Gaia DR4
  data", with an anomaly of about 0.25 mas against about 0.1 mas along-scan accuracy.
- The source is bright (G = 14.8 at baseline, 10.5 at peak) and the lens is a faint
  main-sequence star (V about 20), so the light centroid is that of the lensed source.
- Wicker et al. (2026, arXiv:2611.24561) list Gaia19bld among the events they excluded from
  their DR4 forecast as already studied. An arXiv search (2026-11-24) finds no paper that has
  measured its Gaia astrometric signal. That is not possible before DR4.

So the control is also a small result in its own right: the first comparison of a Gaia-astrometry
Einstein radius with an interferometric one for the same lens.

## Data and model (control.py, control_job.py)

- **Source.** Gaia DR3 5860179620698137600 (G = 14.797, RUWE 0.99), 0.05 arcsec from the
  published position. Its DR4 counterpart, epoch astrometry and normal points are handled
  exactly as in `bh_job.py`, with the same helper functions.
- **Photometric solutions.** These are from Rybicki et al. (2022) Table 5, ground + Gaia +
  Spitzer fits, with t0,par = 2458680.0 as in the paper.
  - Preferred: u0 < 0 (t0 = 2458681.219, u0 = -0.0192, tE = 107.06 d, piEN = -0.0378,
    piEE = -0.0731).
  - Also fitted: u0 > 0.
- **Blend.** G-band blend ratio g = 0.00 +- 0.03, an assumption. The OGLE source fraction is
  1.029 +- 0.016 and the Spitzer blend is consistent with 0. The paper gives no Gaia G blend.
- **Finite source.** rho = 0.032, so transits with u < 0.1 (about 3 rho; within about 10 days of
  the peak) are dropped before the fit. The paper notes that all but one Gaia transit fall
  outside the peak.
- **Fit, Monte Carlo and field-star null.** These are identical to the headline events:
  `bhfit.fit_event`, `bhfit.monte_carlo` (400 draws) and up to 12 DR4 neighbours within 60
  arcsec with |dG| <= 0.75.
- **Trajectories.** If `bh_job.py` switched to MulensModel trajectories (xcheck rule), the
  control does too.
- Because |u0| is only 0.019, the control's centroid track hardly depends on the u0 sign. It
  therefore tests the one convention that matters most for every event: the direction of the
  relative motion, taken from piE. The free-orientation fit reports it.

## Grade (pre-registered)

pull = (thetaE_Gaia - 0.7650) / sqrt(err_total^2 + 0.0038^2), where err_total is the Monte
Carlo total error (as for the headline events).

| condition | label |
|---|---|
| fit not possible | CONTROL_NO_EPOCHS / CONTROL_TOO_FEW / CONTROL_NO_EVENT_COVERAGE |
| z <= -3 | CONTROL_WRONG_SIGN |
| \|pull\| > 2 | CONTROL_BIASED_LOW or CONTROL_BIASED_HIGH |
| \|pull\| <= 2 and DETECTED or DETECTED_UNCALIBRATED | **CONTROL_PASS** |
| \|pull\| <= 2, not detected, err_total > 0.255 mas | CONTROL_UNINFORMATIVE |
| \|pull\| <= 2 otherwise | CONTROL_CONSISTENT_NOT_SIGNIFICANT |

**What the control does to the headline.** It never changes a headline verdict. It sets a flag
that the write-up must carry:

| control label | flag | effect on the write-up |
|---|---|---|
| CONTROL_PASS | VALIDATED | the article may say the pipeline recovered a known lens |
| CONSISTENT_NOT_SIGNIFICANT or UNINFORMATIVE | VALIDATION_WEAK | the article says the control was too weak to check the pipeline |
| any other label | NOT_VALIDATED | the article leads with the control failure, and no headline mass is called a black hole |

**Reported, never graded:**
- the pull against the finite-source thetaE (0.755 +- 0.013);
- the Monte Carlo mass against 1.147 +- 0.029 Msun;
- the u0 > 0 fit;
- a fit that also drops u < 0.2. Near the peak the star was brighter than G = 13, where Gaia's
  window class changes;
- an all-CCD fit.

## Known limits

- The G-band blend is assumed, not measured. A real G blend would dilute the shift and bias the
  control low. Pre-registering it means a low result is reported as CONTROL_BIASED_LOW, not
  explained away.
- Near the peak the star was about 4 magnitudes brighter than at baseline. Calibration offsets
  between window classes could bias those transits. The u < 0.2 variant shows how much this
  matters.
- One control event checks the pipeline on a bright, high-magnification event. The headline
  events are fainter (G = 17.3 and 20.0), so a pass does not prove that their error bars are
  right. The field-star nulls remain the calibration for that.

## Release-day job

`dr4day/job.template.json` now also runs `python -m unittest -v test_darklens test_control` and,
after `bh_job.py`, `python control_job.py --data-dir /data --out-dir /out`. The output is
`/out/control/results_control.json`. The timeout is raised from 360 to 420 minutes, because the
control adds 13 DataLink downloads.

## Credits

- Interferometric thetaE: Cassan et al. (2022).
- Photometric solutions, finite-source thetaE and the DR4 forecast: Rybicki et al. (2022).
- Spectroscopy of the source: Bachelet et al. (2022, arXiv:2112.01676).
- Data: ESA Gaia DR3 and DR4 (Gaia Collaboration).
