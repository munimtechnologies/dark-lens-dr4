# dark-lens-dr4

Can Gaia weigh a black hole that has never been seen? This is a pre-registered test that runs
when Gaia DR4 is released on 2 December 2026.

Rybicki et al. (2024, [arXiv:2407.13740](https://arxiv.org/abs/2407.13740)) found
dark-remnant candidates among long microlensing events. These are lenses of roughly 3 to 5
solar masses that give out no light. The two strongest are OGLE-2016-BLG-0293 (about 3 Msun)
and OGLE-2018-BLG-0483 (about 4.7 Msun).

The paper forecasts that Gaia's epoch astrometry will detect the small shift in each lensed
star's position (about 1 mas). That shift fixes the Einstein radius, which turns the paper's
prior-based masses into measured ones.

**Known-answer control (Amendment 1, 2026-11-24).** The same pipeline is first run on Gaia19bld.
That lens's Einstein radius was measured directly with VLTI interferometry: 0.7650 +- 0.0038 mas
(Cassan et al. 2022). Rybicki et al. (2022, [arXiv:2112.01613](https://arxiv.org/abs/2112.01613))
forecast that Gaia DR4 may detect its astrometric signal. If DR4 recovers that number, the
black-hole masses are backed by a validated pipeline. If it does not, the write-up says so first.

This repository fixes in advance how the test is done and graded.

**Plan and frozen inputs**
- `PREREGISTRATION.md`: the question, model, conventions, gates and verdicts (frozen
  2026-11-24).
- `AMENDMENT-2026-11-24-control.md`: the Gaia19bld control and how it is graded (frozen before
  DR4).
- `events_bh.py`: the nine events, their photometric solutions and the paper's predictions.

**Model and fits**
- `traj.py`: the microlensing trajectory with annual parallax (Gould 2004 / MulensModel
  conventions) and the centroid-shift templates.
- `bhfit.py`: the fits, the Monte Carlo mass and the verdicts.
- `control.py`: the Gaia19bld event, the finite-source cut and the control grade.

**Jobs and checks**
- `bh_job.py` and `control_job.py`: the release-day jobs (Gaia archive TAP + DataLink).
- `xcheck.py`: a cross-check of the trajectories against MulensModel.
- `test_darklens.py` and `test_control.py`: unit tests
  (`python -m unittest -v test_darklens test_control`).

The release-day jobs also need `io_gaia.py` and `gaiaorb/` from
[gj3512-giants-dr4](https://github.com/munimtechnologies/gj3512-giants-dr4). Their setup step
clones them.

By MV2, Munim, Inc. (mv2@munimtech.com). The work was made with AI agents working in MV2's
name, and no human expert has checked it yet. Code is Apache-2.0; text is CC BY 4.0.
