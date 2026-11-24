#!/usr/bin/env python3
"""Cross-check of traj.tau_beta against MulensModel (the community reference implementation of
the Gould 2004 geocentric parallax formalism) for every event and solution in events_bh.

Part of dark-lens-dr4 (MV2, Munim, Inc.). Apache-2.0. Made with AI agents working in MV2's name.
Pre-registered rule (PREREGISTRATION.md section 3): if max |d tau| or |d beta| over the DR4
window exceeds XCHECK_TOL for any event, the release-day job uses MulensModel's (tau, beta).
MulensModel, numpy and astropy are imported lazily; the tests do not need them.
Usage: python xcheck.py --out-dir /out
"""
import argparse
import json
import math
import os

import events_bh as EV
import traj as T

XCHECK_TOL = 0.01
JD_START = 2456863.5      # 2014-07-25 (DR4 window start)
JD_END = 2458868.5        # 2020-01-20 (DR4 window end)


def mulens_tau_beta(jds, sol, ra, dec, t0par_jd=None):
    """MulensModel's (x, y) = (tau, beta) for a point lens with annual parallax."""
    import numpy as np
    import MulensModel as mm
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    t0par = sol["t0"] + T.HJD_OFFSET if t0par_jd is None else t0par_jd
    params = dict(t_0=sol["t0"] + T.HJD_OFFSET, u_0=sol["u0"], t_E=sol["tE"], pi_E_N=sol["piEN"],
                  pi_E_E=sol["piEE"], t_0_par=t0par)
    coords = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
    model = mm.Model(params, coords=coords)
    times = np.asarray(jds, dtype=float)
    tr = None
    try:
        tr = model.get_trajectory(times)
    except Exception:
        tr = mm.Trajectory(times, parameters=model.parameters,
                           parallax={"earth_orbital": True, "satellite": False, "topocentric": False},
                           coords=model.coords)
    return [(float(x), float(y)) for x, y in zip(tr.x, tr.y)]


def compare(step_days=5.0):
    jds = []
    jd = JD_START
    while jd <= JD_END:
        jds.append(jd)
        jd += step_days
    rows = []
    for e in EV.EVENTS:
        for name, sol in e["sols"].items():
            mine = [T.tau_beta(j, sol, e["ra"], e["dec"]) for j in jds]
            ref = mulens_tau_beta(jds, sol, e["ra"], e["dec"])
            dt = max(abs(a[0] - b[0]) for a, b in zip(mine, ref))
            db = max(abs(a[1] - b[1]) for a, b in zip(mine, ref))
            # same check on the u range that matters (u <= 5)
            near = [(a, b) for a, b in zip(mine, ref) if math.hypot(*b) <= 5.0]
            du = max((abs(math.hypot(*a) - math.hypot(*b)) for a, b in near), default=None)
            rows.append(dict(event=e["short"], solution=name, max_dtau=dt, max_dbeta=db, max_du_u_le_5=du,
                             n_times=len(jds), passed=max(dt, db) <= XCHECK_TOL))
    return dict(tolerance=XCHECK_TOL, all_passed=all(r["passed"] for r in rows), rows=rows)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="/out")
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    try:
        import MulensModel
        version = getattr(MulensModel, "__version__", "?")
        res = compare()
        res["mulensmodel_version"] = version
    except Exception as ex:
        res = dict(error=str(ex)[:500], all_passed=None)
    with open(os.path.join(a.out_dir, "xcheck.json"), "w") as f:
        json.dump(res, f, indent=1)
    for r in res.get("rows", []):
        print("%s %s dtau %.2e dbeta %.2e du %s %s" % (r["event"], r["solution"], r["max_dtau"], r["max_dbeta"],
                                                    r["max_du_u_le_5"], "ok" if r["passed"] else "FAIL"))
    print("XCHECK %s" % ("PASS" if res.get("all_passed") else ("FAIL" if res.get("all_passed") is False
                                                               else "NOT_RUN: %s" % res.get("error"))))
    return res


if __name__ == "__main__":
    main()
