"""Point-lens microlensing trajectory with annual parallax (geocentric frame) and the astrometric
centroid-shift templates used to fit Gaia DR4 along-scan (AL) epoch astrometry. Standard library.

Part of dark-lens-dr4 (MV2, Munim, Inc., mv2@munimtech.com). Apache-2.0.
Made with AI agents working in MV2's name; no human expert has checked it yet.

Conventions (frozen 2026-11-24, before Gaia DR4), following Gould (2004) as implemented in
MulensModel (Trajectory._get_xy, _project_delta, _get_delta_annual):
  tau  = (t - t0) / tE + dN piEN + dE piEE
  beta = u0 - dN piEE + dE piEN
  delta_s = x_earth(t0par) - x_earth(t) + (t - t0par) v_earth(t0par)   (au, equatorial)
  dN = delta_s . north,  dE = delta_s . east
  east = (-sin a, cos a, 0),  north = (-sin d cos a, -sin d sin a, cos d)
On the sky (components east, north) the source-minus-lens vector in units of thetaE is
  u_vec = tau e1 + beta e2,  e1 = -(piEE, piEN) / piE,  e2 = (-piEN, piEE) / piE
(derived in PREREGISTRATION.md section 3 and checked by a direct forward simulation in the tests).
The light centroid of the two unresolved images moves away from the lens by
  delta = thetaE u_vec / (u^2 + 2) x A / (A + g)     (g = blend-to-source flux ratio, floored at 0)
so the AL template is linear in thetaE.
"""
import math

VERSION = "0.1.0"
JD_J2000 = 2451545.0
DR4_REF_JYEAR = 2017.5
HJD_OFFSET = 2450000.0
L2_FACTOR = 1.01


def jd_from_t(t):
    """JD of t (Julian years from J2017.5, the gaiaorb.epochs time axis)."""
    return JD_J2000 + (DR4_REF_JYEAR + t - 2000.0) * 365.25


def t_from_jd(jd):
    return (jd - JD_J2000) / 365.25 + 2000.0 - DR4_REF_JYEAR


def sun_geocentric(jd):
    """Geocentric Sun (au, ICRS-aligned equatorial), Astronomical Almanac low-precision
    formulae (about 0.01 deg)."""
    n = jd - JD_J2000
    L = math.radians((280.460 + 0.9856474 * n) % 360.0)
    g = math.radians((357.528 + 0.9856003 * n) % 360.0)
    lam = L + math.radians(1.915) * math.sin(g) + math.radians(0.020) * math.sin(2.0 * g)
    R = 1.00014 - 0.01671 * math.cos(g) - 0.00014 * math.cos(2.0 * g)
    eps = math.radians(23.439 - 0.0000004 * n)
    return (R * math.cos(lam), R * math.cos(eps) * math.sin(lam), R * math.sin(eps) * math.sin(lam))


def earth_xyz(jd):
    x, y, z = sun_geocentric(jd)
    return (-x, -y, -z)


def earth_velocity(jd, h=0.5):
    """au per day, central difference."""
    a = earth_xyz(jd + h)
    b = earth_xyz(jd - h)
    return tuple((a[k] - b[k]) / (2.0 * h) for k in range(3))


def sky_basis(ra_deg, dec_deg):
    a, d = math.radians(ra_deg), math.radians(dec_deg)
    east = (-math.sin(a), math.cos(a), 0.0)
    north = (-math.sin(d) * math.cos(a), -math.sin(d) * math.sin(a), math.cos(d))
    return east, north


def _dot(u, v):
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def delta_ne(jd, t0par_jd, ra_deg, dec_deg):
    """(dN, dE) of the geocentric-frame parallax offset (au)."""
    p0 = earth_xyz(t0par_jd)
    v0 = earth_velocity(t0par_jd)
    p = earth_xyz(jd)
    ds = tuple(p0[k] - p[k] + (jd - t0par_jd) * v0[k] for k in range(3))
    east, north = sky_basis(ra_deg, dec_deg)
    return _dot(ds, north), _dot(ds, east)


def tau_beta_from_ne(jd, sol, dN, dE):
    t0 = sol["t0"] + HJD_OFFSET
    tau = (jd - t0) / sol["tE"] + dN * sol["piEN"] + dE * sol["piEE"]
    beta = sol["u0"] - dN * sol["piEE"] + dE * sol["piEN"]
    return tau, beta


def tau_beta(jd, sol, ra_deg, dec_deg, t0par_jd=None):
    """sol: dict(t0 [HJD-2450000], u0, tE [d], piEN, piEE)."""
    t0par = sol["t0"] + HJD_OFFSET if t0par_jd is None else t0par_jd
    dN, dE = delta_ne(jd, t0par, ra_deg, dec_deg)
    return tau_beta_from_ne(jd, sol, dN, dE)


def ne_offsets(points, ra_deg, dec_deg, t0par_jd):
    """Per-point (dN, dE); independent of the photometric parameters, so the Monte Carlo
    over (t0, u0, tE, piE) reuses them."""
    return [delta_ne(jd_from_t(p["t"]), t0par_jd, ra_deg, dec_deg) for p in points]


def magnification(u):
    u = max(u, 1e-9)
    return (u * u + 2.0) / (u * math.sqrt(u * u + 4.0))


def shift_factor(u, g=0.0):
    """|delta| / (thetaE u): 1/(u^2+2), diluted by a blend of flux ratio g (floored at 0)."""
    geff = max(0.0, g or 0.0)
    dil = 1.0
    if geff > 0.0:
        A = magnification(u)
        dil = A / (A + geff)
    return dil / (u * u + 2.0)


def directions(piEN, piEE):
    """(e1, e2) as (east, north) unit vectors: e1 = direction of increasing tau of the
    source-minus-lens vector (= -piE_hat), e2 = direction of increasing beta."""
    pi = math.hypot(piEN, piEE)
    if pi <= 0.0:
        raise ValueError("piE must be non-zero")
    return (-piEE / pi, -piEN / pi), (-piEN / pi, piEE / pi)


def u_vector(tau, beta, sol):
    e1, e2 = directions(sol["piEN"], sol["piEE"])
    return (tau * e1[0] + beta * e2[0], tau * e1[1] + beta * e2[1])


def taubeta_points(points, sol, ra_deg, dec_deg, t0par_jd=None, ne=None):
    t0par = sol["t0"] + HJD_OFFSET if t0par_jd is None else t0par_jd
    ne = ne_offsets(points, ra_deg, dec_deg, t0par) if ne is None else ne
    return [tau_beta_from_ne(jd_from_t(p["t"]), sol, n[0], n[1]) for p, n in zip(points, ne)]


def column_constrained(points, sol, tb):
    """AL template for thetaE = 1 mas with the orientation fixed by piE (primary model)."""
    e1, e2 = directions(sol["piEN"], sol["piEE"])
    g = sol.get("g", 0.0)
    col = []
    for p, (tau, beta) in zip(points, tb):
        u = math.hypot(tau, beta)
        K = shift_factor(u, g)
        dE = K * (tau * e1[0] + beta * e2[0])
        dN = K * (tau * e1[1] + beta * e2[1])
        col.append(dE * p["s"] + dN * p["c"])
    return col


def columns_free(points, sol, tb, h):
    """Two AL columns (A, B) for a free orientation phi of the tau axis with handedness h:
    delta_AL = A colA + B colB, A = thetaE cos(phi), B = thetaE sin(phi). The constrained
    model is h = -1, phi = atan2(-piEN, -piEE)."""
    g = sol.get("g", 0.0)
    ca, cb = [], []
    for p, (tau, beta) in zip(points, tb):
        u = math.hypot(tau, beta)
        K = shift_factor(u, g)
        ca.append(K * (tau * p["s"] + h * beta * p["c"]))
        cb.append(K * (tau * p["c"] - h * beta * p["s"]))
    return ca, cb


def predicted_phi(sol):
    return math.atan2(-sol["piEN"], -sol["piEE"])


def gaia_parallax_factor(ra_deg, dec_deg, p):
    """Model AL parallax factor for Gaia (1.01 x Earth's heliocentric position): gate G0."""
    x, y, z = earth_xyz(jd_from_t(p["t"]))
    X, Y, Z = L2_FACTOR * x, L2_FACTOR * y, L2_FACTOR * z
    a, d = math.radians(ra_deg), math.radians(dec_deg)
    pa = X * math.sin(a) - Y * math.cos(a)
    pd = X * math.cos(a) * math.sin(d) + Y * math.sin(a) * math.sin(d) - Z * math.cos(d)
    return pa * p["s"] + pd * p["c"]
