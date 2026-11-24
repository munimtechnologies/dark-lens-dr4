#!/usr/bin/env python3
"""Gaia DR4 release-day control of dark-lens-dr4: fit Gaia19bld, whose Einstein radius is known
from VLTI interferometry (0.7650 +- 0.0038 mas), with the same pipeline as the black-hole
candidates, and grade whether DR4 epoch astrometry recovers it.

MV2, Munim, Inc. (mv2@munimtech.com). Apache-2.0. Made with AI agents working in MV2's name.
Frozen with AMENDMENT-2026-11-24-control.md on 2026-11-24, before DR4 (2026-12-02).
Runs after bh_job.py in the same job. If bh_job switched to MulensModel trajectories
(results.json "trajectories"), the control does too.
Usage: python control_job.py --data-dir /data --out-dir /out [--release "Gaia DR4"]
"""
import argparse
import csv
import json
import math
import os
import time

import bh_job as J
import bhfit as B
import control as C
import traj as T

T0 = time.time()


def log(msg):
    print("[%7.1fs] control: %s" % (time.time() - T0, msg), flush=True)


def counterpart(ev):
    """Nearest DR4 entry to the DR3 position propagated to J2017.5 (|dG| < 1 preferred)."""
    ra5, dec5 = J.propagate(ev["ra"], ev["dec"], ev["pm"][0], ev["pm"][1], 2016.0, T.DR4_REF_JYEAR)
    for tbl in J.DR4_TABLES:
        try:
            rows = J.cone(tbl, ra5, dec5, J.MATCH_RADIUS_MAS / 3.6e6)
        except Exception as ex:
            log("cone on %s failed: %s" % (tbl, ex))
            continue
        c = []
        for r in rows:
            ra, dec = J.num(r.get("ra")), J.num(r.get("dec"))
            if ra is None or dec is None:
                continue
            ep = J.num(r.get("ref_epoch")) or T.DR4_REF_JYEAR
            ra_, dec_ = J.propagate(ra, dec, J.num(r.get("pmra")), J.num(r.get("pmdec")), ep, T.DR4_REF_JYEAR)
            g = J.num(r.get("phot_g_mean_mag"))
            bad_g = g is not None and abs(g - ev["g_dr3"]) >= 1.0
            c.append((bad_g, J.sep_arcsec(ra_, dec_, ra5, dec5),
                      dict(source_id=int(r["source_id"]), ra=ra_, dec=dec_, g=g, table=tbl)))
        if c:
            c.sort(key=lambda x: (x[0], x[1]))
            d = c[0][2]
            d["match_arcsec"] = c[0][1]
            return d
    return None


def tb_for(pts, sol, ev, use_mulens):
    if use_mulens:
        import xcheck
        return xcheck.mulens_tau_beta([T.jd_from_t(p["t"]) for p in pts], sol, ev["ra"], ev["dec"],
                                      ev["t0par_jd"])
    return T.taubeta_points(pts, sol, ev["ra"], ev["dec"], ev["t0par_jd"])


def strip(fr):
    return {k: v for k, v in fr.items() if k not in ("keep", "sig")}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--release", default="Gaia DR4")
    a = ap.parse_args(argv)
    out_dir = os.path.join(a.out_dir, "control")
    os.makedirs(out_dir, exist_ok=True)
    use_mulens = False
    try:
        with open(os.path.join(a.out_dir, "results.json")) as f:
            use_mulens = json.load(f).get("trajectories") == "MulensModel"
    except Exception:
        pass
    ev = C.GAIA19BLD
    rec = dict(project="dark-lens-dr4", part="control", control_version=C.VERSION, bhfit_version=B.VERSION,
               event=ev["name"], dr3=str(ev["dr3"]), truth_thetaE=list(C.PIONIER_THETAE),
               trajectories="MulensModel" if use_mulens else "traj.py", release=a.release)
    src = counterpart(ev)
    rec["dr4_source"] = None if src is None else dict(src, source_id=str(src["source_id"]))
    sid = src["source_id"] if src else ev["dr3"]
    rows = J.fetch_rows(sid, a.data_dir, a.release)
    if rows is None and src is not None and sid != ev["dr3"]:
        rows = J.fetch_rows(ev["dr3"], a.data_dir, a.release)
        sid = ev["dr3"]
    if not rows:
        rec["label"] = "CONTROL_NO_EPOCHS"
        rec["headline_flag"] = C.headline_flag(rec["label"])
        _write(out_dir, rec)
        return rec
    rec["epoch_source_id"] = str(sid)
    pts, counts = J.points_from_rows(rows)
    rec["epoch_counts"] = counts
    rec["geometry_median_au"] = B.percentile([abs(T.gaia_parallax_factor(ev["ra"], ev["dec"], p) - p["f"])
                                              for p in pts], 50) if pts else None
    fits = {}
    used = {}
    for name, sol in ev["sols"].items():
        try:
            p2, fr = C.fit_control(pts, sol, ev, tb_all=tb_for(pts, sol, ev, use_mulens))
            fits[name], used[name] = fr, p2
        except Exception as ex:
            fits[name] = dict(status="ERROR", error=str(ex)[:300])
    rec["fits"] = {k: strip(v) for k, v in fits.items()}
    fr = fits[ev["preferred"]]
    sol = ev["sols"][ev["preferred"]]
    if fr.get("status") != "OK":
        rec["label"] = "CONTROL_" + str(fr.get("status"))
        rec["headline_flag"] = C.headline_flag(rec["label"])
        _write(out_dir, rec)
        return rec
    p_used = used[ev["preferred"]]
    mc = B.monte_carlo(p_used, fr, sol, ev["ra"], ev["dec"], t0par_jd=ev["t0par_jd"])
    rec["monte_carlo"] = mc
    nulls = []
    if src is not None:
        for nid in J.neighbours(src, src["g"] if src.get("g") is not None else ev["g_dr3"]):
            try:
                nrows = J.fetch_rows(nid, a.data_dir, a.release)
                if not nrows:
                    continue
                npts, _ = J.points_from_rows(nrows)
                _, nf = C.fit_control(npts, sol, ev, tb_all=tb_for(npts, sol, ev, use_mulens))
                if nf.get("status") == "OK":
                    nulls.append(nf["snr"])
            except Exception as ex:
                log("null %s failed: %s" % (nid, str(ex)[:160]))
    rec["null_snrs"] = nulls
    fap = B.empirical_fap(fr["snr"], nulls) if fr["snr"] > 0 else None
    rec["fap"] = fap
    v = C.control_verdict(fr, mc["thetaE_total_err"], fap)
    rec.update(v)
    rec["headline_flag"] = C.headline_flag(v["label"])
    rec["pull_vs_finite_source"] = (fr["thetaE"] - C.FINITE_SOURCE_THETAE[0]) / math.sqrt(
        mc["thetaE_total_err"] ** 2 + C.FINITE_SOURCE_THETAE[1] ** 2)
    rec["mass_p"] = mc["mass_p"]
    rec["mass_reference"] = list(C.CASSAN_MASS)
    # robustness (reported, never graded)
    rob = {}
    try:
        _, f2 = C.fit_control(pts, sol, ev, u_min=C.U_BRIGHT, tb_all=tb_for(pts, sol, ev, use_mulens))
        rob["exclude_u_below_0.2"] = strip(f2)
    except Exception as ex:
        rob["exclude_u_below_0.2"] = dict(status="ERROR", error=str(ex)[:200])
    try:
        apts, _ = J.points_from_rows(rows, all_ccd=True)
        _, f3 = C.fit_control(apts, sol, ev, tb_all=tb_for(apts, sol, ev, use_mulens))
        rob["all_ccd"] = strip(f3)
    except Exception as ex:
        rob["all_ccd"] = dict(status="ERROR", error=str(ex)[:200])
    rec["robustness"] = rob
    _points_csv(os.path.join(out_dir, "points_gaia19bld.csv"), p_used, fr, sol, ev, use_mulens)
    _write(out_dir, rec)
    return rec


def _points_csv(path, pts, fr, sol, ev, use_mulens):
    tb = tb_for(pts, sol, ev, use_mulens)
    col = T.column_constrained(pts, sol, tb)
    src = B.source_columns(pts)
    beta = [fr["source"][k][0] for k in B.SOURCE_PARAMS]
    keep = set(fr["keep"])
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["jyear", "hjd_minus_2450000", "u", "w_minus_source_mas", "sig_mas", "model_shift_al_mas", "used"])
        for i, p in enumerate(pts):
            srcmod = sum(beta[k] * src[k][i] for k in range(5))
            w.writerow(["%.5f" % (T.DR4_REF_JYEAR + p["t"]), "%.3f" % (T.jd_from_t(p["t"]) - T.HJD_OFFSET),
                        "%.4f" % math.hypot(*tb[i]), "%.4f" % (p["w"] - srcmod), "%.4f" % p["sig"],
                        "%.4f" % (fr["thetaE"] * col[i]), int(i in keep)])


def _write(out_dir, rec):
    rec["runtime_s"] = round(time.time() - T0, 1)
    with open(os.path.join(out_dir, "results_control.json"), "w") as f:
        json.dump(rec, f, indent=1, default=str)
    log("CONTROL %s (%s)" % (rec.get("label"), rec.get("headline_flag")))


if __name__ == "__main__":
    main()
