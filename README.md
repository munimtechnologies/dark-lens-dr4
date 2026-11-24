# dark-lens-dr4

Can Gaia weigh a black hole that has never been seen? A pre-registered test that runs when
Gaia DR4 is released on 2 December 2026.

Rybicki et al. (2024, [arXiv:2407.13740](https://arxiv.org/abs/2407.13740)) found
dark-remnant candidates among long microlensing events. These are lenses of roughly 3 to 5
solar masses that give out no light. The two strongest are OGLE-2016-BLG-0293 (about 3 Msun)
and OGLE-2018-BLG-0483 (about 4.7 Msun). The paper forecasts that Gaia's epoch astrometry
will detect the small shift in each lensed star's position (about 1 mas). That shift fixes the
Einstein radius, so it turns the paper's prior-based masses into measured ones.

This repository fixes in advance how that test is done and graded.

- `PREREGISTRATION.md`: the question, model, conventions, gates and verdicts (frozen
  2026-11-24).
- `events_bh.py`: the nine events, their photometric solutions and the paper's predictions.
- `traj.py`: microlensing trajectory with annual parallax (Gould 2004 / MulensModel
  conventions) and the centroid-shift templates.
- `bhfit.py`: the fits, the Monte Carlo mass and the verdicts.
- `bh_job.py`: the release-day job (Gaia archive TAP + DataLink).
- `xcheck.py`: cross-check of the trajectories against MulensModel.
- `test_darklens.py`: unit tests (`python -m unittest -v test_darklens`). They include a
  forward simulation that checks the sky orientation.

The release-day job also needs `io_gaia.py` and `gaiaorb/` from
[gj3512-giants-dr4](https://github.com/munimtechnologies/gj3512-giants-dr4), which its setup
step clones.

By MV2, Munim, Inc. (mv2@munimtech.com). The work was made with AI agents working in MV2's
name, and no human expert has checked it yet. Code is Apache-2.0; text is CC BY 4.0.
