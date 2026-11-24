"""Unit tests for dark-lens-dr4 (standard library only): python -m unittest -v test_darklens

MV2, Munim, Inc. Apache-2.0. Made with AI agents working in MV2's name.
"""
import math
import random
import unittest

import bhfit as B
import events_bh as EV
import traj as T

EV_BY = EV.by_short()
E293 = EV_BY["OB160293"]
SOL = E293["sols"]["u0-"]
RA, DEC = E293["ra"], E293["dec"]
SRC = (0.4, -0.3, 0.12, -1.5, -7.5)


def synth(theta, sol=SOL, n=130, noise=0.1, seed=1, ra=RA, dec=DEC, h_true=-1):
    """Synthetic per-transit AL points: source motion + centroid shift (orientation from piE when
    h_true = -1; the mirror image when h_true = +1)."""
    rng = random.Random(seed)
    pts = []
    for i in range(n):
        t = -2.9 + 5.45 * (i + rng.random()) / n
        psi = rng.uniform(0.0, 2.0 * math.pi)
        p = dict(t=t, s=math.sin(psi), c=math.cos(psi), fr=0.0, sig=noise)
        p["f"] = T.gaia_parallax_factor(ra, dec, p)
        tau, beta = T.tau_beta(T.jd_from_t(t), sol, ra, dec)
        if h_true == -1:
            col = T.column_constrained([p], sol, [(tau, beta)])[0]
        else:
            phi = T.predicted_phi(sol)
            ca, cb = T.columns_free([p], sol, [(tau, beta)], 1)
            col = math.cos(phi) * ca[0] + math.sin(phi) * cb[0]
        w = (SRC[0] + SRC[3] * t) * p["s"] + (SRC[1] + SRC[4] * t) * p["c"] + SRC[2] * p["f"]
        p["w"] = w + theta * col + rng.gauss(0.0, noise)
        pts.append(p)
    return pts


class TestEphemeris(unittest.TestCase):
    def test_equinox(self):
        # 2017 Mar 20 12:00 TT (JD 2457833.0): Sun near lambda 0 -> Earth near -x
        x, y, z = T.earth_xyz(2457833.0)
        self.assertAlmostEqual(x, -0.996, delta=0.01)
        self.assertLess(abs(y), 0.02)
        self.assertLess(abs(z), 0.01)

    def test_solstice(self):
        # 2017 Jun 21 (JD 2457926.0): Sun at lambda 90 -> Earth at -(0, cos eps, sin eps) * 1.016
        x, y, z = T.earth_xyz(2457926.0)
        self.assertLess(abs(x), 0.02)
        self.assertAlmostEqual(y, -1.016 * math.cos(math.radians(23.44)), delta=0.02)
        self.assertAlmostEqual(z, -1.016 * math.sin(math.radians(23.44)), delta=0.02)

    def test_velocity(self):
        v = T.earth_velocity(2457833.0)
        self.assertAlmostEqual(math.sqrt(sum(c * c for c in v)), 2.0 * math.pi / 365.25, delta=0.0009)
        # at the March equinox Earth (at -x) moves toward -y (prograde, counter-clockwise from +z)
        self.assertLess(v[1], 0.0)

    def test_second_order_offset(self):
        # delta_s = x(t0p) - x(t) + dt v(t0p) ~ +0.5 w^2 dt^2 x(t0p): away from the Sun
        jd0, dt = 2457669.6, 10.0
        p0, v0, p = T.earth_xyz(jd0), T.earth_velocity(jd0), T.earth_xyz(jd0 + dt)
        ds = [p0[k] - p[k] + dt * v0[k] for k in range(3)]
        mag = math.sqrt(sum(c * c for c in ds))
        self.assertAlmostEqual(mag, 0.5 * (2.0 * math.pi / 365.25) ** 2 * dt * dt, delta=0.002)
        cosang = sum(a * b for a, b in zip(ds, p0)) / (mag * math.sqrt(sum(c * c for c in p0)))
        self.assertGreater(cosang, 0.98)

    def test_basis(self):
        east, north = T.sky_basis(RA, DEC)
        a, d = math.radians(RA), math.radians(DEC)
        direction = (math.cos(d) * math.cos(a), math.cos(d) * math.sin(a), math.sin(d))
        for u in (east, north):
            self.assertAlmostEqual(sum(c * c for c in u), 1.0)
            self.assertAlmostEqual(sum(x * y for x, y in zip(u, direction)), 0.0)
        self.assertAlmostEqual(sum(x * y for x, y in zip(east, north)), 0.0)
        self.assertGreater(north[2], 0.0)

    def test_time_axis(self):
        self.assertAlmostEqual(T.jd_from_t(0.0), 2451545.0 + 17.5 * 365.25)
        self.assertAlmostEqual(T.t_from_jd(T.jd_from_t(-0.73)), -0.73)


class TestTrajectory(unittest.TestCase):
    def test_reference_epoch(self):
        tau, beta = T.tau_beta(SOL["t0"] + T.HJD_OFFSET, SOL, RA, DEC)
        self.assertAlmostEqual(tau, 0.0, places=9)
        self.assertAlmostEqual(beta, SOL["u0"], places=9)

    def test_no_parallax(self):
        s = dict(SOL, piEN=0.0, piEE=0.0)
        tau, beta = T.tau_beta(SOL["t0"] + T.HJD_OFFSET + 57.15, s, RA, DEC)
        self.assertAlmostEqual(tau, 57.15 / SOL["tE"])
        self.assertAlmostEqual(beta, SOL["u0"])

    def test_directions(self):
        e1, e2 = T.directions(0.03, -0.04)
        self.assertAlmostEqual(e1[0], 0.04 / 0.05)
        self.assertAlmostEqual(e1[1], -0.03 / 0.05)
        self.assertAlmostEqual(e1[0] * e2[0] + e1[1] * e2[1], 0.0)
        self.assertAlmostEqual(math.hypot(*e2), 1.0)

    def test_forward_simulation(self):
        """Independent check of the sky orientation: physical source-minus-lens track with the
        apparent parallax displacement -plx * (O . east, O . north) for each star."""
        for sol in (SOL, E293["sols"]["u0+"], EV_BY["OB180483"]["sols"]["u0-"]):
            ra, dec = RA, DEC
            pi = math.hypot(sol["piEN"], sol["piEE"])
            e1, e2 = T.directions(sol["piEN"], sol["piEE"])
            t0 = sol["t0"] + T.HJD_OFFSET
            a, d = math.radians(ra), math.radians(dec)

            def P(jd):
                X, Y, Z = T.earth_xyz(jd)
                return (X * math.sin(a) - Y * math.cos(a),
                        X * math.cos(a) * math.sin(d) + Y * math.sin(a) * math.sin(d) - Z * math.cos(d))

            P0 = P(t0)
            Pp, Pm = P(t0 + 0.5), P(t0 - 0.5)
            Pdot = ((Pp[0] - Pm[0]) / 1.0, (Pp[1] - Pm[1]) / 1.0)
            for dt in (-400.0, -150.0, -20.0, 35.0, 200.0, 510.0):
                jd = t0 + dt
                Pt = P(jd)
                # theta_S - theta_L = r_bary - plx_rel P; geocentric frame keeps the non-linear part
                nl = (Pt[0] - P0[0] - Pdot[0] * dt, Pt[1] - P0[1] - Pdot[1] * dt)
                tau0 = dt / sol["tE"]
                uE = sol["u0"] * e2[0] + tau0 * e1[0] - pi * nl[0]
                uN = sol["u0"] * e2[1] + tau0 * e1[1] - pi * nl[1]
                tau, beta = T.tau_beta(jd, sol, ra, dec)
                vE, vN = T.u_vector(tau, beta, sol)
                self.assertAlmostEqual(vE, uE, delta=1e-6)
                self.assertAlmostEqual(vN, uN, delta=1e-6)
                self.assertGreater(abs(pi * nl[0]) + abs(pi * nl[1]), 0.0)

    def test_constrained_equals_free(self):
        pts = synth(0.0, n=25)
        tb = T.taubeta_points(pts, SOL, RA, DEC)
        con = T.column_constrained(pts, SOL, tb)
        ca, cb = T.columns_free(pts, SOL, tb, -1)
        phi = T.predicted_phi(SOL)
        for k in range(len(pts)):
            self.assertAlmostEqual(con[k], math.cos(phi) * ca[k] + math.sin(phi) * cb[k], places=12)

    def test_shift_size(self):
        # dark lens, no blend: |delta| = thetaE u / (u^2 + 2); maximum thetaE / sqrt(8) at u = sqrt 2
        u = math.sqrt(2.0)
        self.assertAlmostEqual(u * T.shift_factor(u), 1.0 / math.sqrt(8.0))
        self.assertAlmostEqual(T.shift_factor(1.0, -0.3), 1.0 / 3.0)       # negative g floored at 0
        A = T.magnification(1.0)
        self.assertAlmostEqual(A, 3.0 / math.sqrt(5.0))
        self.assertAlmostEqual(T.shift_factor(1.0, 0.5), (A / (A + 0.5)) / 3.0)


class TestFit(unittest.TestCase):
    def test_recovers_thetaE(self):
        fr = B.fit_event(synth(3.0, noise=0.05), SOL, RA, DEC)
        self.assertEqual(fr["status"], "OK")
        self.assertLess(abs(fr["thetaE"] - 3.0), 4.0 * fr["thetaE_err"])
        self.assertGreater(fr["snr"], 10.0)
        self.assertGreater(fr["dchi2"], 100.0)
        self.assertEqual(fr["direction"]["best_h"], -1)
        self.assertEqual(fr["direction"]["label"], "CONSISTENT")
        self.assertLess(abs(fr["source"]["pmdec_masyr"][0] - SRC[4]), 5.0 * fr["source"]["pmdec_masyr"][1])
        self.assertEqual(B.detection_verdict(fr, 0.05), "DETECTED")
        self.assertEqual(B.detection_verdict(fr, None), "DETECTED_UNCALIBRATED")
        self.assertEqual(B.detection_verdict(fr, 0.5), "NOT_SIGNIFICANT_VS_FIELD")

    def test_mirror_is_flagged(self):
        fr = B.fit_event(synth(3.0, noise=0.05, h_true=1, seed=4), SOL, RA, DEC)
        self.assertEqual(fr["direction"]["best_h"], 1)
        self.assertEqual(fr["direction"]["label"], "INCONSISTENT")

    def test_null(self):
        fr = B.fit_event(synth(0.0, seed=7), SOL, RA, DEC)
        self.assertLess(abs(fr["snr"]), 4.0)

    def test_wrong_sign(self):
        fr = B.fit_event(synth(-2.0, noise=0.05, seed=3), SOL, RA, DEC)
        self.assertEqual(B.detection_verdict(fr, None), "WRONG_SIGN")

    def test_too_few_and_coverage(self):
        self.assertEqual(B.fit_event(synth(0.0, n=12), SOL, RA, DEC)["status"], "TOO_FEW")
        far = [p for p in synth(0.0, n=200) if abs(T.jd_from_t(p["t"]) - SOL["t0"] - T.HJD_OFFSET) > 700.0]
        self.assertEqual(B.fit_event(far, SOL, RA, DEC)["status"], "NO_EVENT_COVERAGE")

    def test_clip(self):
        pts = synth(1.0, seed=9)
        pts[40]["w"] += 50.0
        fr = B.fit_event(pts, SOL, RA, DEC)
        self.assertNotIn(40, fr["keep"])
        self.assertGreaterEqual(fr["n_clipped"], 1)

    def test_monte_carlo(self):
        pts = synth(2.0, noise=0.05, seed=5)
        fr = B.fit_event(pts, SOL, RA, DEC)
        mc = B.monte_carlo(pts, fr, SOL, RA, DEC, n=80)
        self.assertEqual(mc["n"], 80)
        piE = math.hypot(SOL["piEN"], SOL["piEE"])
        m_true = 2.0 / (B.KAPPA * piE)
        self.assertLess(abs(mc["mass_p"][2] - m_true) / m_true, 0.25)
        self.assertLess(mc["mass_p"][0], mc["mass_p"][2])
        self.assertLess(mc["mass_p"][2], mc["mass_p"][4])
        self.assertGreaterEqual(mc["thetaE_total_err"], fr["thetaE_err"])

    def test_wls_exact(self):
        f = B.solve_wls([[1.0] * 4, [0.0, 1.0, 2.0, 3.0]], [1.0, 3.0, 5.0, 7.0], [1.0] * 4)
        self.assertAlmostEqual(f["beta"][0], 1.0)
        self.assertAlmostEqual(f["beta"][1], 2.0)
        self.assertAlmostEqual(f["chi2"], 0.0)


class TestGrading(unittest.TestCase):
    def test_fap(self):
        self.assertIsNone(B.empirical_fap(5.0, [0.1] * 7))
        self.assertAlmostEqual(B.empirical_fap(5.0, [0.5] * 10), 1.0 / 11.0)
        self.assertAlmostEqual(B.empirical_fap(2.0, [0.5] * 8 + [-2.5, 3.0]), 3.0 / 11.0)

    def test_mass_class(self):
        self.assertEqual(B.mass_class([2.6, 3, 4, 5, 6]), "BLACK_HOLE_MASS_RANGE")
        self.assertEqual(B.mass_class([1.0, 2, 2.6, 3, 4]), "BLACK_HOLE_FAVOURED_NOT_SECURE")
        self.assertEqual(B.mass_class([0.2, 0.5, 1.0, 1.5, 2.0]), "BELOW_BLACK_HOLE_RANGE")

    def test_prediction(self):
        kro = EV_BY["OB160293"]["pred"]["kroupa"]           # thetaE 0.91 -0.39 +0.52
        self.assertEqual(B.prediction_check(1.0, 0.2, kro)["label"], "AS_PREDICTED")
        self.assertEqual(B.prediction_check(3.0, 0.2, kro)["label"], "ABOVE_PREDICTION")
        lo, hi = B.prediction_check(1.0, 0.2, kro)["range"]
        self.assertAlmostEqual(lo, 0.91 - 0.78 - 0.4)
        self.assertAlmostEqual(hi, 0.91 + 1.04 + 0.4)

    def test_nondetection(self):
        piE = math.hypot(-0.037, -0.006)
        v = B.nondetection_verdict(0.0, 0.1, 0.91, piE, 0.85)
        self.assertTrue(v["informative"])
        self.assertAlmostEqual(v["mass_upper95"], 0.1645 / (B.KAPPA * piE))
        self.assertEqual(v["label"], "PREDICTION_FAILED")
        self.assertEqual(B.nondetection_verdict(0.0, 0.5, 0.91, piE, 0.85)["label"], "UNINFORMATIVE")
        self.assertEqual(B.nondetection_verdict(0.5, 0.25, 0.91, piE, 0.85)["label"], "NOT_DETECTED_CONSISTENT")

    def test_overall(self):
        self.assertEqual(B.overall(False, {}), "GEOMETRY_FAIL")
        h = {"OB160293": dict(verdict="DETECTED", mass_class="BLACK_HOLE_MASS_RANGE"),
             "OB180483": dict(verdict="NOT_DETECTED", mass_class=None)}
        self.assertEqual(B.overall(True, h), "1_OF_2_DETECTED_1_BH_RANGE")

    def test_direction_labels(self):
        free = {"-1": dict(chi2=10.0, phi_deg=10.0, phi_err_deg=5.0), "1": dict(chi2=20.0, phi_deg=0.0, phi_err_deg=5.0)}
        self.assertEqual(B.direction_check(free, math.radians(20.0))["label"], "CONSISTENT")
        self.assertEqual(B.direction_check(free, math.radians(190.0))["label"], "SIGN_FLIPPED")
        self.assertEqual(B.direction_check(free, math.radians(100.0))["label"], "INCONSISTENT")

    def test_percentile(self):
        self.assertAlmostEqual(B.percentile([1, 2, 3, 4, 5], 50), 3.0)
        self.assertAlmostEqual(B.percentile([0.0, 10.0], 25), 2.5)


class TestEvents(unittest.TestCase):
    TABLE3 = {"OB150145": (270.17782, -35.15408), "OB150149": (270.28819, -32.55773),
              "OB150211": (262.35909, -30.98178), "OB160293": (268.16140, -32.48960),
              "OB160689": (261.14944, -30.13245), "OB180410": (272.22770, -27.16974),
              "OB180483": (262.64177, -27.49183), "OB180662": (266.87754, -32.52442),
              "OB190169": (265.98559, -32.87095)}

    def test_table(self):
        self.assertEqual(len(EV.EVENTS), 9)
        self.assertEqual(set(EV_BY), set(self.TABLE3))
        for e in EV.EVENTS:
            self.assertIn(e["preferred"], e["sols"])
            self.assertEqual(set(e["sols"]), {"u0+", "u0-"})
            for s in e["sols"].values():
                for k in ("t0", "u0", "tE", "piEN", "piEE", "g"):
                    self.assertEqual(len(s[k + "_err"]), 2)
                    self.assertTrue(all(x > 0 for x in s[k + "_err"]))
                self.assertGreater(s["tE"], 0.0)
            ra, dec = self.TABLE3[e["short"]]
            self.assertAlmostEqual(e["ra"], ra, delta=6e-6)
            self.assertAlmostEqual(e["dec"], dec, delta=6e-6)
        self.assertEqual([e["short"] for e in EV.ordered()][:2], list(EV.HEADLINE))

    def test_table4_consistency(self):
        # M = thetaE / (kappa piE) holds for the medians to 5% (headline) / 25% (others)
        for e in EV.EVENTS:
            s = e["sols"][e["preferred"]]
            k = e["pred"]["kroupa"]
            m = k["thetaE"][0] / (B.KAPPA * math.hypot(s["piEN"], s["piEE"]))
            tol = 0.05 if e["short"] in EV.HEADLINE else 0.25
            self.assertLess(abs(m - k["M"][0]) / k["M"][0], tol, e["short"])

    def test_window(self):
        # every event's t0 is inside the DR4 window (2014-07-25 to 2020-01-20)
        for e in EV.EVENTS:
            for s in e["sols"].values():
                jd = s["t0"] + T.HJD_OFFSET
                self.assertTrue(2456863.5 < jd < 2458868.5, e["short"])


class TestJobHelpers(unittest.TestCase):
    def test_helpers(self):
        import bh_job
        self.assertIsNone(bh_job.num("nan"))
        self.assertEqual(bh_job.num("2.5"), 2.5)
        ra, dec = bh_job.propagate(10.0, 60.0, 3600.0, 1800.0, 2015.5, 2017.5)
        self.assertAlmostEqual(dec - 60.0, 3600.0 / 3.6e6)
        self.assertAlmostEqual(bh_job.sep_arcsec(10.0, 0.0, 10.0, 1.0 / 3600.0), 1.0)
        rows = [dict(used_by_agis_al=[True, False, None]), dict(used_by_agis_al=False)]
        f = bh_job.force_all_used(rows)
        self.assertEqual(f[0]["used_by_agis_al"], [True, True, True])
        self.assertTrue(f[1]["used_by_agis_al"])
        self.assertEqual(rows[0]["used_by_agis_al"], [True, False, None])
        self.assertEqual(set(bh_job.DR3_PM), set(EV_BY))


class TestXcheckGrid(unittest.TestCase):
    """The MulensModel comparison itself runs in xcheck.py (never blocks the release-day job);
    here only its time window is checked."""
    def test_window(self):
        import xcheck
        self.assertAlmostEqual(xcheck.JD_START, 2456863.5)
        self.assertAlmostEqual(xcheck.JD_END, 2458868.5)
        self.assertLess(xcheck.XCHECK_TOL, 0.02)


if __name__ == "__main__":
    unittest.main()
