#!/usr/bin/env python3
"""Gaia DR4 release-day job of dark-lens-dr4: test Rybicki et al. (2024)'s forecast that Gaia's
epoch astrometry detects the light-centroid shift of the dark-remnant microlensing candidates
OGLE-2016-BLG-0293 and OGLE-2018-BLG-0483, and weigh all nine lenses.

MV2, Munim, Inc. (mv2@munimtech.com). Apache-2.0. Made with AI agents working in MV2's name.
Frozen with PREREGISTRATION.md on 2026-11-24, before DR4 (2026-12-02).

Per event (headline events first): DR3 position (frozen in events_bh) -> DR4 counterpart by
cone search at J2017.5 (fallback: DataLink with the DR3 source_id) -> DR4 epoch astrometry by
DataLink (io_gaia, gaiaorb.epochs, unchanged from gj3512-giants-dr4) -> per-transit normal
points -> gate G0 -> bhfit.fit_event for both u0 solutions -> Monte Carlo mass (preferred
solution) -> field-star null (same template on up to N_NULL DR4 neighbours of similar G) ->
verdicts. Robustness (reported, not used in verdicts): all CCDs instead of AGIS-used ones,
t0par = t0 -/+ 30 d.
Usage: python bh_job.py --data-dir /data --out-dir /out [--release "Gaia DR4"] [--only OB160293]
"""
import argparse
import csv
import json
import math
import os
import time

import bhfit as B
import events_bh as EV
import traj as T

T0 = time.time()
TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
DR4_TABLES = ("gaiadr4.gaia_source", "gaiadr4.all_source_astrometry")
MATCH_RADIUS_MAS = 500.0
N_NULL = 12
NULL_RADIUS_ARCSEC = 60.0
NULL_MIN_SEP_ARCSEC = 3.0
NULL_DG = 0.75
GEOMETRY_TOL = 0.03          # au, median |f_model - f_published| over events (gate G0)
NULL_BUDGET_S = 4.0 * 3600   # after this, non-headline events get no field-star nulls


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg), flush=True)


def num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def tap(query, tries=4):
    import io_gaia
    last = None
    for k in range(tries):
        try:
            return io_gaia.tap_csv(query, TAP)
        except Exception as e:  # archive busy on release day: back off and retry
            last = e
            log("TAP retry %d: %s" % (k + 1, str(e)[:200]))
            time.sleep(20 * (k + 1))
    raise RuntimeError("TAP failed: %s" % last)


def propagate(ra, dec, pmra, pmdec, from_ep, to_ep):
    dt = to_ep - from_ep
    ddec = (pmdec or 0.0) * dt / 3.6e6
    dra = (pmra or 0.0) * dt / 3.6e6 / max(1e-6, math.cos(math.radians(dec)))
    return ra + dra, dec + ddec


def sep_arcsec(ra1, dec1, ra2, dec2):
    da = ((ra1 - ra2 + 180.0) % 360.0 - 180.0) * math.cos(math.radians(0.5 * (dec1 + dec2)))
    return 3600.0 * math.hypot(da, dec1 - dec2)


def cone(tbl, ra, dec, radius_deg, cols="*"):
    q = ("SELECT %s FROM %s WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', %.9f, %.9f, %.9f))"
         % (cols, tbl, ra, dec, radius_deg))
    return tap(q)


def dr4_counterpart(e):
    """Nearest DR4 entry to the DR3 position propagated to J2017.5 (|dG| < 1 preferred)."""
    from_ep = 2016.0
    pm = DR3_PM.get(e["short"], (0.0, 0.0))
    ra5, dec5 = propagate(e["ra"], e["dec"], pm[0], pm[1], from_ep, T.DR4_REF_JYEAR)
    for tbl in DR4_TABLES:
        try:
            rows = cone(tbl, ra5, dec5, MATCH_RADIUS_MAS / 3.6e6)
        except Exception as ex:
            log("cone on %s failed: %s" % (tbl, ex))
            continue
        c = []
        for r in rows:
            ra, dec = num(r.get("ra")), num(r.get("dec"))
            if ra is None or dec is None:
                continue
            ep = num(r.get("ref_epoch")) or T.DR4_REF_JYEAR
            ra_, dec_ = propagate(ra, dec, num(r.get("pmra")), num(r.get("pmdec")), ep, T.DR4_REF_JYEAR)
            g = num(r.get("phot_g_mean_mag"))
            bad_g = g is not None and abs(g - e["g_dr3"]) >= 1.0
            c.append((bad_g, sep_arcsec(ra_, dec_, ra5, dec5), dict(source_id=int(r["source_id"]), ra=ra_, dec=dec_,
                                                                   g=g, table=tbl)))
        if c:
            c.sort(key=lambda x: (x[0], x[1]))
            d = c[0][2]
            d["match_arcsec"] = c[0][1]
            return d
    return None


# DR3 proper motions (mas/yr) of the nine sources, queried 2026-11-24 with the events_bh values
DR3_PM = {
    "OB150145": (-1.2887467565383721, -6.51215974689122),
    "OB150149": (-2.0057972430656568, -5.027318201117429),
    "OB150211": (-3.4980923295961404, -7.216149599360904),
    "OB160293": (-1.4902066744813138, -7.480056649595083),
    "OB160689": (-1.7967000776884277, -7.253207192647181),
    "OB180410": (-4.508991936382236, -6.839648468807622),
    "OB180483": (4.851757103402024, -2.1346189676654026),
    "OB180662": (-5.744689830275599, -6.091668530924434),
    "OB190169": (-0.8739368800030132, -4.540589536191184),
}


def force_all_used(rows):
    """Copy of the rows with every CCD marked used_by_agis_al (robustness variant)."""
    out = []
    for r in rows:
        r2 = dict(r)
        v = r.get("used_by_agis_al")
        if isinstance(v, (list, tuple)):
            r2["used_by_agis_al"] = [True] * len(v)
        else:
            r2["used_by_agis_al"] = True
        out.append(r2)
    return out


def fetch_rows(source_id, data_dir, release):
    import io_gaia
    from gaiaorb import epochs
    d = os.path.join(data_dir, "epochs")
    for url in io_gaia.datalink_urls(source_id, release):
        blob = None
        for k in range(3):
            try:
                blob = io_gaia.fetch(url, timeout=300)
                break
            except Exception as ex:
                log("DataLink retry %d %s: %s" % (k + 1, source_id, str(ex)[:160]))
                time.sleep(10 * (k + 1))
        if not blob:
            continue
        paths = io_gaia.save_epochs(blob, d, "src_%d" % source_id, log)
        if not paths:
            continue
        rows = []
        for p in paths:
            rows.extend(io_gaia.read_table_rows(p))
        by = epochs.rows_by_source(rows)
        mine = by.get(source_id) or (rows if len(by) == 1 else [])
        if mine:
            return mine
    return None


def points_from_rows(rows, all_ccd=False):
    from gaiaorb import epochs
    return epochs.normal_points(force_all_used(rows) if all_ccd else rows)


def neighbours(src, g_src, k=N_NULL):
    try:
        rows = cone("gaiadr4.gaia_source", src["ra"], src["dec"], NULL_RADIUS_ARCSEC / 3600.0,
                    cols="source_id, ra, dec, phot_g_mean_mag")
    except Exception as ex:
        log("neighbour cone failed: %s" % ex)
        return []
    c = []
    for r in rows:
        ra, dec, g = num(r.get("ra")), num(r.get("dec")), num(r.get("phot_g_mean_mag"))
        if None in (ra, dec, g) or int(r["source_id"]) == src["source_id"]:
            continue
        s = sep_arcsec(ra, dec, src["ra"], src["dec"])
        if s < NULL_MIN_SEP_ARCSEC or abs(g - g_src) > NULL_DG:
            continue
        c.append((s, int(r["source_id"])))
    c.sort()
    return [sid for _, sid in c[:k]]


def trajectories(points, sol, e, use_mulens):
    if use_mulens:
        import xcheck
        return xcheck.mulens_tau_beta([T.jd_from_t(p["t"]) for p in points], sol, e["ra"], e["dec"])
    return T.taubeta_points(points, sol, e["ra"], e["dec"])


def strip(fr):
    return {k: v for k, v in fr.items() if k not in ("keep", "sig")}


def analyse(e, data_dir, release, use_mulens, do_nulls):
    rec = dict(event=e["short"], name=e["name"], dr3=str(e["dr3"]), preferred=e["preferred"],
               headline=e["short"] in EV.HEADLINE)
    src = dr4_counterpart(e)
    sid = src["source_id"] if src else e["dr3"]
    rec["dr4_source"] = None if src is None else dict(src, source_id=str(src["source_id"]))
    rows = fetch_rows(sid, data_dir, release)
    if rows is None and src is not None and sid != e["dr3"]:
        rows = fetch_rows(e["dr3"], data_dir, release)
        sid = e["dr3"]
    if not rows:
        rec["verdict"] = "NO_EPOCHS"
        return rec, None
    rec["epoch_source_id"] = str(sid)
    pts, counts = points_from_rows(rows)
    rec["epoch_counts"] = counts
    rec["geometry_median_au"] = B.percentile([abs(T.gaia_parallax_factor(e["ra"], e["dec"], p) - p["f"])
                                              for p in pts], 50) if pts else None
    fits = {}
    for name, sol in e["sols"].items():
        try:
            tb = trajectories(pts, sol, e, use_mulens)
            fits[name] = B.fit_event(pts, sol, e["ra"], e["dec"], tb=tb)
        except Exception as ex:
            fits[name] = dict(status="ERROR", error=str(ex)[:300])
    rec["fits"] = {k: strip(v) for k, v in fits.items()}
    fr = fits[e["preferred"]]
    if fits.get("u0+", {}).get("status") == "OK" and fits.get("u0-", {}).get("status") == "OK":
        rec["chi2_u0plus_minus_u0minus"] = fits["u0+"]["chi2"] - fits["u0-"]["chi2"]
    if fr.get("status") != "OK":
        rec["verdict"] = fr.get("status", "FAILED")
        return rec, None
    sol = e["sols"][e["preferred"]]
    mc = B.monte_carlo(pts, fr, sol, e["ra"], e["dec"])
    rec["monte_carlo"] = mc
    # field-star null
    null_snrs = []
    if do_nulls and src is not None:
        for nid in neighbours(src, src["g"] if src.get("g") is not None else e["g_dr3"]):
            try:
                nrows = fetch_rows(nid, data_dir, release)
                if not nrows:
                    continue
                npts, _ = points_from_rows(nrows)
                ntb = trajectories(npts, sol, e, use_mulens)
                nf = B.fit_event(npts, sol, e["ra"], e["dec"], tb=ntb)
                if nf.get("status") == "OK":
                    null_snrs.append(nf["snr"])
            except Exception as ex:
                log("null %s failed: %s" % (nid, str(ex)[:160]))
    rec["null_snrs"] = null_snrs
    fap = B.empirical_fap(fr["snr"], null_snrs) if fr["snr"] > 0 else None
    rec["fap"] = fap
    rec["verdict"] = B.detection_verdict(fr, fap)
    piE = math.hypot(sol["piEN"], sol["piEE"])
    kro = e["pred"]["kroupa"]
    if rec["verdict"] in ("DETECTED", "DETECTED_UNCALIBRATED"):
        rec["mass_class"] = B.mass_class(mc["mass_p"])
        rec["prediction"] = B.prediction_check(fr["thetaE"], mc["thetaE_total_err"], kro)
    else:
        rec["nondetection"] = B.nondetection_verdict(fr["thetaE"], mc["thetaE_total_err"], kro["thetaE"][0],
                                                     piE, kro["M_2sigma_min"])
    # robustness (reported only)
    rob = {}
    try:
        apts, _ = points_from_rows(rows, all_ccd=True)
        af = B.fit_event(apts, sol, e["ra"], e["dec"], tb=trajectories(apts, sol, e, use_mulens))
        rob["all_ccd"] = strip(af)
    except Exception as ex:
        rob["all_ccd"] = dict(status="ERROR", error=str(ex)[:200])
    for dd in (-30.0, 30.0):
        try:
            tp = sol["t0"] + T.HJD_OFFSET + dd
            f2 = B.fit_event(pts, sol, e["ra"], e["dec"], tb=T.taubeta_points(pts, sol, e["ra"], e["dec"], tp))
            rob["t0par_%+d" % dd] = dict(status=f2.get("status"), thetaE=f2.get("thetaE"), thetaE_err=f2.get("thetaE_err"))
        except Exception as ex:
            rob["t0par_%+d" % dd] = dict(status="ERROR", error=str(ex)[:200])
    rec["robustness"] = rob
    return rec, (pts, fr, sol, e)


def write_points(path, pts, fr, sol, e, use_mulens):
    tb = trajectories(pts, sol, e, use_mulens)
    col = T.column_constrained(pts, sol, tb)
    src = B.source_columns(pts)
    beta = [fr["source"][k][0] for k in B.SOURCE_PARAMS]
    keep = set(fr["keep"])
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["jyear", "hjd_minus_2450000", "u", "w_minus_source_mas", "sig_mas", "model_shift_al_mas", "used"])
        for i, p in enumerate(pts):
            srcmod = sum(beta[a] * src[a][i] for a in range(5))
            w.writerow(["%.5f" % (T.DR4_REF_JYEAR + p["t"]), "%.3f" % (T.jd_from_t(p["t"]) - T.HJD_OFFSET),
                        "%.4f" % math.hypot(*tb[i]), "%.4f" % (p["w"] - srcmod), "%.4f" % p["sig"],
                        "%.4f" % (fr["thetaE"] * col[i]), int(i in keep)])


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--release", default="Gaia DR4")
    ap.add_argument("--only", default="")
    a = ap.parse_args(argv)
    os.makedirs(a.data_dir, exist_ok=True)
    os.makedirs(os.path.join(a.out_dir, "points"), exist_ok=True)
    # 1. trajectory convention check against MulensModel (pre-registered switch)
    use_mulens = False
    try:
        import xcheck
        xc = xcheck.compare()
        use_mulens = xc["all_passed"] is False
    except Exception as ex:
        xc = dict(error=str(ex)[:300], all_passed=None)
    log("MulensModel cross-check: %s -> trajectories from %s" % (xc.get("all_passed"),
                                                                "MulensModel" if use_mulens else "traj.py"))
    evs = EV.ordered()
    only = {x.strip() for x in a.only.split(",") if x.strip()}
    if only:
        evs = [e for e in evs if e["short"] in only]
    recs = []
    for e in evs:
        do_nulls = e["short"] in EV.HEADLINE or (time.time() - T0) < NULL_BUDGET_S
        try:
            rec, extra = analyse(e, a.data_dir, a.release, use_mulens, do_nulls)
        except Exception as ex:  # one bad event must not stop the rest
            rec, extra = dict(event=e["short"], verdict="ERROR", error=str(ex)[:500]), None
        if extra is not None:
            write_points(os.path.join(a.out_dir, "points", "%s.csv" % e["short"]), *extra, use_mulens=use_mulens)
        log("%s: %s %s" % (e["short"], rec.get("verdict"), rec.get("mass_class") or ""))
        recs.append(rec)
    geo = [r["geometry_median_au"] for r in recs if r.get("geometry_median_au") is not None]
    g0 = bool(geo) and B.percentile(geo, 50) < GEOMETRY_TOL
    head = {r["event"]: dict(verdict=r.get("verdict"), mass_class=r.get("mass_class"))
            for r in recs if r["event"] in EV.HEADLINE}
    outcome = B.overall(g0, head)
    res = dict(project="dark-lens-dr4", bhfit_version=B.VERSION, traj_version=T.VERSION, release=a.release,
               xcheck=xc, trajectories="MulensModel" if use_mulens else "traj.py",
               gate_G0=dict(passed=g0, median_au=B.percentile(geo, 50) if geo else None, n=len(geo)),
               outcome=outcome, events=recs, runtime_s=round(time.time() - T0, 1))
    with open(os.path.join(a.out_dir, "results.json"), "w") as f:
        json.dump(res, f, indent=1, default=str)
    cols = ["event", "headline", "verdict", "thetaE", "thetaE_err", "snr", "dchi2", "fap", "mass_p16", "mass_p50",
            "mass_p84", "mass_p2_5", "mass_class", "prediction", "nondetection", "direction"]
    with open(os.path.join(a.out_dir, "events.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in recs:
            fr = r.get("fits", {}).get(r.get("preferred"), {})
            mp = (r.get("monte_carlo") or {}).get("mass_p") or [None] * 5
            w.writerow([r["event"], r.get("headline"), r.get("verdict"), fr.get("thetaE"), fr.get("thetaE_err"),
                        fr.get("snr"), fr.get("dchi2"), r.get("fap"), mp[1], mp[2], mp[3], mp[0], r.get("mass_class"),
                        (r.get("prediction") or {}).get("label"), (r.get("nondetection") or {}).get("label"),
                        (fr.get("direction") or {}).get("label")])
    log("OUTCOME %s" % outcome)
    return res


if __name__ == "__main__":
    main()
