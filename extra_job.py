#!/usr/bin/env python3
"""Gaia DR4 release-day step for Amendment 2 of dark-lens-dr4: fit Gaia19dke and Gaia18cbf with
the same pipeline as the black-hole candidates and grade them as frozen in
AMENDMENT-2026-11-24-extra.md.

MV2, Munim, Inc. (mv2@munimtech.com). Apache-2.0. Made with AI agents working in MV2's name.
Frozen on 2026-11-24, before DR4 (2026-12-02). Runs after bh_job.py and control_job.py in the
same job; follows their trajectory choice (traj.py or MulensModel) and carries the control flag.
Usage: python extra_job.py --data-dir /data --out-dir /out [--release "Gaia DR4"]
"""
import argparse
import json
import math
import os
import time

import bh_job as J
import bhfit as B
import control_job as CJ
import extra_events as X

T0 = time.time()


def log(msg):
    print("[%7.1fs] extra: %s" % (time.time() - T0, msg), flush=True)


def drop_low_u(pts, tb, u_min):
    keep = [i for i, (a, b) in enumerate(tb) if math.hypot(a, b) >= u_min]
    return [pts[i] for i in keep], [tb[i] for i in keep]


def analyse(ev, data_dir, release, use_mulens, out_dir=None):
    rec = dict(event=ev["name"], dr3=str(ev["dr3"]), preferred=ev["preferred"], u_event=ev["u_event"],
               claim=ev["claim"])
    src = CJ.counterpart(ev)
    rec["dr4_source"] = None if src is None else dict(src, source_id=str(src["source_id"]))
    sid = src["source_id"] if src else ev["dr3"]
    rows = J.fetch_rows(sid, data_dir, release)
    if rows is None and src is not None and sid != ev["dr3"]:
        rows = J.fetch_rows(ev["dr3"], data_dir, release)
        sid = ev["dr3"]
    if not rows:
        rec["label"] = "NO_EPOCHS"
        return rec
    rec["epoch_source_id"] = str(sid)
    pts, counts = J.points_from_rows(rows)
    rec["epoch_counts"] = counts
    fits = {}
    for name, sol in ev["sols"].items():
        try:
            fits[name] = X.fit_extra(pts, sol, ev, CJ.tb_for(pts, sol, ev, use_mulens))
        except Exception as ex:
            fits[name] = dict(status="ERROR", error=str(ex)[:300])
    rec["fits"] = {k: CJ.strip(v) for k, v in fits.items()}
    fr = fits[ev["preferred"]]
    sol = ev["sols"][ev["preferred"]]
    if fr.get("status") != "OK":
        rec.update(X.grade(fr, None, None, None, ev))
        return rec
    mc = B.monte_carlo(pts, fr, sol, ev["ra"], ev["dec"], t0par_jd=ev["t0par_jd"])
    rec["monte_carlo"] = mc
    nulls = []
    if src is not None:
        for nid in J.neighbours(src, src["g"] if src.get("g") is not None else ev["g_dr3"]):
            try:
                nrows = J.fetch_rows(nid, data_dir, release)
                if not nrows:
                    continue
                npts, _ = J.points_from_rows(nrows)
                nf = X.fit_extra(npts, sol, ev, CJ.tb_for(npts, sol, ev, use_mulens))
                if nf.get("status") == "OK":
                    nulls.append(nf["snr"])
            except Exception as ex:
                log("null %s failed: %s" % (nid, str(ex)[:160]))
    rec["null_snrs"] = nulls
    fap = B.empirical_fap(fr["snr"], nulls) if fr["snr"] > 0 else None
    rec["fap"] = fap
    rec.update(X.grade(fr, mc["thetaE_total_err"], mc["mass_p"], fap, ev))
    rec["mass_p"] = mc["mass_p"]
    # robustness (reported, never graded)
    rob = {}
    for name, f in fits.items():
        if name != ev["preferred"] and f.get("status") == "OK":
            rob["other_solution_" + name] = dict(thetaE=f["thetaE"], thetaE_err=f["thetaE_err"], snr=f["snr"])
    try:
        tb = CJ.tb_for(pts, sol, ev, use_mulens)
        p2, tb2 = drop_low_u(pts, tb, 0.1)
        rob["exclude_u_below_0.1"] = CJ.strip(X.fit_extra(p2, sol, ev, tb2))
    except Exception as ex:
        rob["exclude_u_below_0.1"] = dict(status="ERROR", error=str(ex)[:200])
    try:
        apts, _ = J.points_from_rows(rows, all_ccd=True)
        rob["all_ccd"] = CJ.strip(X.fit_extra(apts, sol, ev, CJ.tb_for(apts, sol, ev, use_mulens)))
    except Exception as ex:
        rob["all_ccd"] = dict(status="ERROR", error=str(ex)[:200])
    rec["robustness"] = rob
    if out_dir:
        CJ._points_csv(os.path.join(out_dir, "points_%s.csv" % ev["name"].lower()), pts, fr, sol, ev, use_mulens)
    return rec


def control_flag(out_dir):
    try:
        with open(os.path.join(out_dir, "control", "results_control.json")) as f:
            return json.load(f).get("headline_flag")
    except Exception:
        return None


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--release", default="Gaia DR4")
    a = ap.parse_args(argv)
    out_dir = os.path.join(a.out_dir, "extra")
    os.makedirs(out_dir, exist_ok=True)
    use_mulens = False
    try:
        with open(os.path.join(a.out_dir, "results.json")) as f:
            use_mulens = json.load(f).get("trajectories") == "MulensModel"
    except Exception:
        pass
    res = dict(project="dark-lens-dr4", part="amendment-2", extra_version=X.VERSION, bhfit_version=B.VERSION,
               release=a.release, trajectories="MulensModel" if use_mulens else "traj.py",
               control_flag=control_flag(a.out_dir), events=[])
    for ev in X.EVENTS:
        try:
            rec = analyse(ev, a.data_dir, a.release, use_mulens, out_dir)
        except Exception as ex:  # one bad event must not stop the other
            rec = dict(event=ev["name"], label="ERROR", error=str(ex)[:500])
        log("%s: %s" % (ev["name"], rec.get("label")))
        res["events"].append(rec)
    res["runtime_s"] = round(time.time() - T0, 1)
    with open(os.path.join(out_dir, "results_extra.json"), "w") as f:
        json.dump(res, f, indent=1, default=str)
    return res


if __name__ == "__main__":
    main()
