"""Unit tests for the Gaia19bld known-answer control of dark-lens-dr4 (standard library only):
python -m unittest -v test_control

MV2, Munim, Inc. Apache-2.0. Made with AI agents working in MV2's name.
"""
import math
import random
import unittest

import bhfit as B
import control as C
import traj as T

EV = C.GAIA19BLD
SOL = EV["sols"]["u0-"]
RA, DEC = EV["ra"], EV["dec"]
SRC = (0.2, -0.4, 0.08, -7.43, 0.09)


def synth(theta, n=160, noise=0.05, seed=11):
    """Synthetic per-transit AL points over the DR4 window: source motion + centroid shift of
    Gaia19bld (orientation from piE, t0par = 2458680.0)."""
    rng = random.Random(seed)
    t_lo, t_hi = T.t_from_jd(2456863.5), T.t_from_jd(2458868.5)
    pts = []
    for i in range(n):
        t = t_lo + (t_hi - t_lo) * (i + rng.random()) / n
        psi = rng.uniform(0.0, 2.0 * math.pi)
        p = dict(t=t, s=math.sin(psi), c=math.cos(psi), fr=0.0, sig=noise)
        p["f"] = T.gaia_parallax_factor(RA, DEC, p)
        tau, beta = T.tau_beta(T.jd_from_t(t), SOL, RA, DEC, EV["t0par_jd"])
        col = T.column_constrained([p], SOL, [(tau, beta)])[0]
        w = (SRC[0] + SRC[3] * t) * p["s"] + (SRC[1] + SRC[4] * t) * p["c"] + SRC[2] * p["f"]
        p["w"] = w + theta * col + rng.gauss(0.0, noise)
        pts.append(p)
    return pts


class TestControlTable(unittest.TestCase):
    def test_identity(self):
        self.assertEqual(EV["dr3"], 5860179620698137600)
        pra, pdec = EV["paper_radec"]
        da = (RA - pra) * math.cos(math.radians(DEC)) * 3600.0
        dd = (DEC - pdec) * 3600.0
        self.assertLess(math.hypot(da, dd), 0.2)
        self.assertAlmostEqual(EV["g_dr3"], 14.8, delta=0.05)

    def test_solutions(self):
        self.assertEqual(set(EV["sols"]), {"u0+", "u0-"})
        self.assertIn(EV["preferred"], EV["sols"])
        for s in EV["sols"].values():
            for k in ("t0", "u0", "tE", "piEN", "piEE", "g"):
                self.assertEqual(len(s[k + "_err"]), 2)
                self.assertTrue(all(x > 0 for x in s[k + "_err"]))
            jd = s["t0"] + T.HJD_OFFSET
            self.assertTrue(2456863.5 < jd < 2458868.5)
        self.assertLess(SOL["u0"], 0.0)

    def test_mass_matches_interferometry(self):
        # thetaE / (kappa piE) with the PIONIER thetaE reproduces Cassan et al.'s 1.147 Msun to 3%
        m = C.mass_from_thetaE(C.PIONIER_THETAE[0], SOL)
        self.assertLess(abs(m - C.CASSAN_MASS[0]) / C.CASSAN_MASS[0], 0.03)

    def test_parallax_along_relative_motion(self):
        mn, me = EV["mu_rel_geo_ne"]
        a1 = math.degrees(math.atan2(SOL["piEN"], SOL["piEE"]))
        a2 = math.degrees(math.atan2(mn, me))
        self.assertLess(B.angle_diff_deg(a1, a2), 5.0)


class TestSelect(unittest.TestCase):
    def test_peak_dropped(self):
        t0 = T.t_from_jd(SOL["t0"] + T.HJD_OFFSET)
        pk = dict(t=t0, s=1.0, c=0.0, fr=0.0, sig=0.1, w=0.0, f=0.0)
        far = dict(t=t0 + 200.0 / 365.25, s=0.0, c=1.0, fr=0.0, sig=0.1, w=0.0, f=0.0)
        pts, tb, n = C.select([pk, far], SOL)
        self.assertEqual(n, 1)
        self.assertEqual(len(pts), 1)
        self.assertIs(pts[0], far)
        self.assertGreater(math.hypot(*tb[0]), 1.5)

    def test_bright_cut_drops_more(self):
        pts = synth(0.765)
        a = C.select(pts, SOL)[2]
        b = C.select(pts, SOL, u_min=C.U_BRIGHT)[2]
        self.assertGreaterEqual(b, a)


class TestControlFit(unittest.TestCase):
    def test_recovers_interferometric_thetaE(self):
        # tiny noise so the 2-sigma pass rule is dominated by the PIONIER error (no flaky seed)
        pts = synth(0.765, noise=0.002)
        used, fr = C.fit_control(pts, SOL)
        self.assertEqual(fr["status"], "OK")
        self.assertLess(abs(fr["thetaE"] - 0.765), 4.0 * fr["thetaE_err"])
        self.assertGreater(fr["snr"], 3.0)
        self.assertLess(abs(B.angle_diff_deg(fr["free"]["-1"]["phi_deg"],
                                             math.degrees(T.predicted_phi(SOL)))), 30.0)
        mc = B.monte_carlo(used, fr, SOL, RA, DEC, t0par_jd=EV["t0par_jd"], n=60)
        m_true = C.mass_from_thetaE(0.765, SOL)
        self.assertLess(abs(mc["mass_p"][2] - m_true) / m_true, 0.2)
        v = C.control_verdict(fr, mc["thetaE_total_err"], None)
        self.assertEqual(v["label"], "CONTROL_PASS")
        self.assertEqual(C.headline_flag(v["label"]), "VALIDATED")

    def test_no_signal_is_not_a_pass(self):
        used, fr = C.fit_control(synth(0.0, seed=5), SOL)
        self.assertEqual(fr["status"], "OK")
        v = C.control_verdict(fr, fr["thetaE_err"], None)
        self.assertNotEqual(v["label"], "CONTROL_PASS")


class TestControlVerdicts(unittest.TestCase):
    def test_labels(self):
        ok = lambda th, snr: dict(status="OK", thetaE=th, snr=snr)
        v = C.control_verdict(ok(0.75, 6.0), 0.12, 0.05)
        self.assertEqual(v["detection"], "DETECTED")
        self.assertEqual(v["label"], "CONTROL_PASS")
        self.assertAlmostEqual(v["pull"], -0.015 / math.sqrt(0.12 ** 2 + 0.0038 ** 2), places=9)
        self.assertEqual(C.control_verdict(ok(0.75, 6.0), 0.12, None)["label"], "CONTROL_PASS")
        self.assertEqual(C.control_verdict(ok(0.30, 3.0), 0.10, None)["label"], "CONTROL_BIASED_LOW")
        self.assertEqual(C.control_verdict(ok(1.40, 14.0), 0.10, None)["label"], "CONTROL_BIASED_HIGH")
        self.assertEqual(C.control_verdict(ok(0.50, 2.5), 0.20, None)["label"],
                         "CONTROL_CONSISTENT_NOT_SIGNIFICANT")
        self.assertEqual(C.control_verdict(ok(0.70, 7.0), 0.10, 0.5)["label"],
                         "CONTROL_CONSISTENT_NOT_SIGNIFICANT")
        self.assertEqual(C.control_verdict(ok(0.40, 1.0), 0.40, None)["label"], "CONTROL_UNINFORMATIVE")
        self.assertEqual(C.control_verdict(ok(-0.50, -5.0), 0.10, None)["label"], "CONTROL_WRONG_SIGN")
        self.assertEqual(C.control_verdict(dict(status="TOO_FEW"), 0.1, None)["label"], "CONTROL_TOO_FEW")

    def test_headline_flag(self):
        self.assertEqual(C.headline_flag("CONTROL_PASS"), "VALIDATED")
        self.assertEqual(C.headline_flag("CONTROL_UNINFORMATIVE"), "VALIDATION_WEAK")
        self.assertEqual(C.headline_flag("CONTROL_CONSISTENT_NOT_SIGNIFICANT"), "VALIDATION_WEAK")
        for lab in ("CONTROL_BIASED_LOW", "CONTROL_BIASED_HIGH", "CONTROL_WRONG_SIGN", "CONTROL_NO_EPOCHS",
                    "CONTROL_TOO_FEW"):
            self.assertEqual(C.headline_flag(lab), "NOT_VALIDATED")


class TestControlJobHelpers(unittest.TestCase):
    def test_tb_and_strip(self):
        import control_job
        pts = synth(0.0, n=8)
        tb = control_job.tb_for(pts, SOL, EV, False)
        ref = T.taubeta_points(pts, SOL, RA, DEC, EV["t0par_jd"])
        for a, b in zip(tb, ref):
            self.assertAlmostEqual(a[0], b[0])
            self.assertAlmostEqual(a[1], b[1])
        self.assertEqual(control_job.strip(dict(a=1, keep=[1], sig=[2])), dict(a=1))


if __name__ == "__main__":
    unittest.main()
