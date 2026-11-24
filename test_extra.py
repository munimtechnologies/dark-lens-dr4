"""Unit tests for Amendment 2 of dark-lens-dr4 (Gaia19dke, Gaia18cbf). Standard library only:
python -m unittest -v test_extra

MV2, Munim, Inc. Apache-2.0. Made with AI agents working in MV2's name.
"""
import math
import random
import unittest

import bhfit as B
import extra_events as X
import traj as T

SRC = (0.2, -0.4, 0.08, -2.0, -3.0)


def synth(ev, theta, n=160, noise=0.05, seed=11):
    """Synthetic per-transit AL points over the DR4 window: source motion + centroid shift of
    the event's preferred solution (orientation from piE, the event's t0par)."""
    sol = ev["sols"][ev["preferred"]]
    rng = random.Random(seed)
    t_lo, t_hi = T.t_from_jd(X.DR4_JD[0]), T.t_from_jd(X.DR4_JD[1])
    pts = []
    for i in range(n):
        t = t_lo + (t_hi - t_lo) * (i + rng.random()) / n
        psi = rng.uniform(0.0, 2.0 * math.pi)
        p = dict(t=t, s=math.sin(psi), c=math.cos(psi), fr=0.0, sig=noise)
        p["f"] = T.gaia_parallax_factor(ev["ra"], ev["dec"], p)
        tau, beta = T.tau_beta(T.jd_from_t(t), sol, ev["ra"], ev["dec"], ev["t0par_jd"])
        col = T.column_constrained([p], sol, [(tau, beta)])[0]
        w = (SRC[0] + SRC[3] * t) * p["s"] + (SRC[1] + SRC[4] * t) * p["c"] + SRC[2] * p["f"]
        p["w"] = w + theta * col + rng.gauss(0.0, noise)
        pts.append(p)
    return pts


def tb_of(pts, ev):
    sol = ev["sols"][ev["preferred"]]
    return T.taubeta_points(pts, sol, ev["ra"], ev["dec"], ev["t0par_jd"])


class TestTables(unittest.TestCase):
    def test_identity(self):
        self.assertEqual(X.GAIA18CBF["dr3"], 5995134605563805056)
        self.assertEqual(X.GAIA19DKE["dr3"], 2026409795566972544)
        for ev in X.EVENTS:
            pra, pdec = ev["paper_radec"]
            da = (ev["ra"] - pra) * math.cos(math.radians(ev["dec"])) * 3600.0
            dd = (ev["dec"] - pdec) * 3600.0
            self.assertLess(math.hypot(da, dd), 0.2)
        self.assertAlmostEqual(X.GAIA18CBF["g_dr3"], 20.17, delta=0.05)
        self.assertAlmostEqual(X.GAIA19DKE["g_dr3"], 15.45, delta=0.05)

    def test_solutions(self):
        for ev in X.EVENTS:
            self.assertIn(ev["preferred"], ev["sols"])
            for key, s in ev["sols"].items():
                for k in ("t0", "u0", "tE", "piEN", "piEE", "g"):
                    self.assertEqual(len(s[k + "_err"]), 2)
                    self.assertTrue(all(x > 0 for x in s[k + "_err"]))
                self.assertEqual(s["u0"] < 0, key == "u0-")
        # Gaia18cbf peaked inside the DR4 window; Gaia19dke peaked after DR4's last epoch
        t18 = X.GAIA18CBF["sols"]["u0+"]["t0"] + T.HJD_OFFSET
        t19 = X.GAIA19DKE["sols"]["u0-"]["t0"] + T.HJD_OFFSET
        self.assertTrue(X.DR4_JD[0] < t18 < X.DR4_JD[1])
        self.assertGreater(t19, X.DR4_JD[1])
        self.assertEqual(X.GAIA19DKE["u_event"], 3.0)
        self.assertEqual(X.GAIA18CBF["u_event"], B.U_EVENT)

    def test_blend_from_fs(self):
        for ev, key, sol in ((X.GAIA18CBF, "u0p", "u0+"), (X.GAIA18CBF, "u0m", "u0-"),
                             (X.GAIA19DKE, "u0m", "u0-")):
            g, (gm, gp) = X.g_from_fs(*ev["fs"][key])
            s = ev["sols"][sol]
            self.assertAlmostEqual(g, s["g"], places=3)
            self.assertAlmostEqual(gm, s["g_err"][0], places=3)
            self.assertAlmostEqual(gp, s["g_err"][1], places=3)
        self.assertAlmostEqual(X.g_from_fs(0.894, 0.226, 0.263)[0], 0.118568, places=5)

    def test_mass_formula(self):
        # thetaE / (8.144 piE): hand values 2.53 / 1.035355 and 0.90 / 1.732945
        self.assertAlmostEqual(X.mass_from_thetaE(2.53, X.GAIA18CBF["sols"]["u0+"]), 2.4436, delta=0.002)
        self.assertAlmostEqual(X.mass_from_thetaE(0.90, X.GAIA19DKE["sols"]["u0-"]), 0.5193, delta=0.002)


def ok(th, snr):
    return dict(status="OK", thetaE=th, snr=snr)


class TestGrade(unittest.TestCase):
    def test_nondetections_gaia18cbf(self):
        ev = X.GAIA18CBF
        kp = X.kappa_pie(ev["sols"]["u0+"])
        self.assertAlmostEqual(kp, 1.035355, places=4)
        g = X.grade(ok(0.1, 0.2), 0.5, None, None, ev)
        self.assertEqual(g["label"], "UPPER_LIMIT_BELOW_1.4_MSUN")
        self.assertAlmostEqual(g["thetaE_upper95"], 0.9225, places=9)
        self.assertAlmostEqual(g["mass_upper95"], 0.9225 / kp, places=9)
        self.assertTrue(g["informative"])
        self.assertEqual(X.grade(ok(1.0, 1.25), 0.8, None, None, ev)["label"], "UPPER_LIMIT_BELOW_2.5_MSUN")
        g = X.grade(ok(1.0, 1.0), 1.0, None, None, ev)
        self.assertEqual(g["label"], "UNINFORMATIVE")
        self.assertFalse(g["informative"])
        self.assertEqual(X.grade(ok(2.0, 2.5), 0.8, None, None, ev)["label"], "NOT_DETECTED_CONSISTENT")
        self.assertEqual(X.grade(ok(-2.0, -4.0), 0.5, None, None, ev)["label"], "WRONG_SIGN")
        g = X.grade(ok(-0.3, -0.6), 0.5, None, None, ev)
        self.assertAlmostEqual(g["thetaE_upper95"], 0.8225, places=9)
        self.assertEqual(X.grade(dict(status="TOO_FEW"), None, None, None, ev)["label"], "TOO_FEW")

    def test_detections(self):
        ev = X.GAIA18CBF
        g = X.grade(ok(3.0, 5.0), 0.6, [2.6, 2.7, 2.9, 3.1, 3.4], None, ev)
        self.assertEqual(g["label"], "DETECTED_UNCALIBRATED")
        self.assertEqual(g["remnant_class"], "BLACK_HOLE_MASS_RANGE")
        self.assertEqual(g["prediction"]["label"], "AS_PREDICTED")
        g = X.grade(ok(3.0, 5.0), 0.6, [2.6, 2.7, 2.9, 3.1, 3.4], 0.05, ev)
        self.assertEqual(g["label"], "DETECTED")
        # detected in its own data but not significant against field stars: graded as a non-detection
        g = X.grade(ok(3.0, 5.0), 0.6, [2.6, 2.7, 2.9, 3.1, 3.4], 0.5, ev)
        self.assertEqual(g["detection"], "NOT_SIGNIFICANT_VS_FIELD")
        self.assertEqual(g["label"], "NOT_DETECTED_CONSISTENT")

    def test_gaia19dke(self):
        ev = X.GAIA19DKE
        g = X.grade(ok(0.1, 0.5), 0.2, None, None, ev)
        self.assertEqual(g["label"], "UPPER_LIMIT_BELOW_1.4_MSUN")
        self.assertAlmostEqual(g["mass_upper95"], (0.1 + 1.645 * 0.2) / 1.732945, places=4)
        g = X.grade(ok(0.9, 6.0), 0.15, [0.4, 0.45, 0.5, 0.55, 0.6], None, ev)
        self.assertEqual(g["remnant_class"], "WHITE_DWARF_OR_STAR_RANGE")

    def test_remnant_class(self):
        self.assertEqual(X.remnant_class([2.6, 2.7, 2.9, 3.1, 3.4]), "BLACK_HOLE_MASS_RANGE")
        self.assertEqual(X.remnant_class([1.5, 1.8, 2.0, 2.2, 2.4]), "ABOVE_WHITE_DWARF_LIMIT")
        self.assertEqual(X.remnant_class([0.4, 0.45, 0.5, 0.55, 0.6]), "WHITE_DWARF_OR_STAR_RANGE")
        self.assertEqual(X.remnant_class([1.0, 1.5, 2.0, 2.5, 3.0]), "UNRESOLVED")


class TestFits(unittest.TestCase):
    def test_recovers_injected_thetaE(self):
        for ev, theta in ((X.GAIA18CBF, 2.53), (X.GAIA19DKE, 0.90)):
            pts = synth(ev, theta, noise=0.002)
            fr = X.fit_extra(pts, ev["sols"][ev["preferred"]], ev, tb_of(pts, ev))
            self.assertEqual(fr["status"], "OK", ev["name"])
            self.assertLess(abs(fr["thetaE"] - theta), 4.0 * fr["thetaE_err"], ev["name"])
            self.assertGreater(fr["snr"], 3.0, ev["name"])
            self.assertGreaterEqual(fr["n_event_points"], 3)
            self.assertEqual(B.U_EVENT, 2.0)   # the module gate is restored after each fit

    def test_gaia19dke_gate_is_u3(self):
        ev = X.GAIA19DKE
        pts = synth(ev, 0.9, noise=0.002)
        tb = tb_of(pts, ev)
        fr = X.fit_extra(pts, ev["sols"]["u0-"], ev, tb)
        n3 = sum(1 for a, b in tb if math.hypot(a, b) <= 3.0)
        self.assertEqual(fr["n_event_points"], n3)
        self.assertGreater(n3, sum(1 for a, b in tb if math.hypot(a, b) <= 2.0))

    def test_no_signal(self):
        for ev in X.EVENTS:
            pts = synth(ev, 0.0, seed=5)
            fr = X.fit_extra(pts, ev["sols"][ev["preferred"]], ev, tb_of(pts, ev))
            self.assertEqual(fr["status"], "OK")
            self.assertLess(abs(fr["thetaE"]), 0.5, ev["name"])

    def test_monte_carlo_runs(self):
        for ev in X.EVENTS:
            sol = ev["sols"][ev["preferred"]]
            pts = synth(ev, 1.0, noise=0.01)
            fr = X.fit_extra(pts, sol, ev, tb_of(pts, ev))
            mc = B.monte_carlo(pts, fr, sol, ev["ra"], ev["dec"], t0par_jd=ev["t0par_jd"], n=40)
            self.assertGreaterEqual(mc["n"], 30)
            self.assertGreaterEqual(mc["thetaE_total_err"], fr["thetaE_err"])


class TestJobHelpers(unittest.TestCase):
    def test_helpers(self):
        import extra_job
        pts = [dict(i=k) for k in range(4)]
        tb = [(0.05, 0.0), (0.2, 0.0), (0.0, 0.09), (1.0, 1.0)]
        p2, tb2 = extra_job.drop_low_u(pts, tb, 0.1)
        self.assertEqual([p["i"] for p in p2], [1, 3])
        self.assertEqual(len(tb2), 2)
        self.assertIsNone(extra_job.control_flag("/nonexistent-dir-for-test"))


if __name__ == "__main__":
    unittest.main()
