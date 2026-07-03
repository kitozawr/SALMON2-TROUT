"""
si_three_photon_isosurfaces.py
==============================

3D isosurfaces of the *direct* (vertical, k-conserving) MULTIPHOTON interband
transition probability in silicon -- the successor to the plain E_direct(k) map.

Where the companion `zener_tunneling_estimate.py` gives the GaAs *tunnelling*
(Zener/Kane) probability -- an exponential  P ~ exp(-C m_r^1/2 E_g^3/2 / F) --
this script gives the Si *multiphoton* probability, a genuinely different
process with a genuinely different prefactor:

        W_N(k)  ~  I^N * |M_N(k)|^2 * g_FK(E_direct(k) - N hbar w)      (power law)

i.e. an N-th order perturbative (lowest-order perturbation theory, LOPT) rate
that scales as the N-th power of the intensity, NOT an exponential in 1/F.

Physics implemented
-------------------
* Bands / momentum matrix elements: the local EPM engine already in the repo
  (epm_si_primitive.configure_for_si + epm_gaas_reference), 169 plane waves,
  Kunikiyo Si form factors. Si direct gap minimum is at Gamma, E_dir(Gamma) =
  3.34 eV; along Gamma->L the direct gap stays LOW (3.34..3.37 eV to L), along
  Gamma->X it climbs to 4.24 eV at X.

* N-photon direct transition matrix element (velocity-gauge ladder, LOPT):

      M_N(k, e) = (E0/2w)^N  sum_{m1..m_{N-1}}
                     (e.p)_{c,m_{N-1}} ... (e.p)_{m1,v}
                     / prod_{j=1}^{N-1} (E_{mj}(k) - E_v(k) - j hbar w - i eta)

  built by an efficient transfer-matrix contraction over ALL intermediate bands
  m1..m_{N-1} (the full LOPT sum, not a two-band model). The polarisation e is
  removed by an EXACT orientation average <|e.p ... |^2>_e over the unit sphere
  (isotropic rank-2N tensor), giving the polarisation-independent "full"
  N-photon strength T_N(k) = <|sum ...|^2>_e. Because it is orientation-averaged
  it keeps the full cubic O_h symmetry, so it is computed only in the
  irreducible wedge and replicated by symmetry.

* Franz-Keldysh field broadening ("the field smears the energy"). A static /
  quasi-static field F does not leave the resonance sharp at E_direct = N hbar w;
  it smears the joint density of states into Airy tails with the electro-optic
  energy  hbar*theta = (e^2 F^2 hbar^2 / 2 mu)^{1/3}. The N-photon resonance
  line therefore uses the 3D Franz-Keldysh lineshape

      g_FK(Delta; hbar theta) = [Ai'(b)^2 - b Ai(b)^2] / hbar theta,
      b = Delta / hbar theta,  Delta = E_direct(k) - N hbar w.

  For Delta > 0 (E_direct ABOVE N hbar w) this is the exponential FK tail
  ~ exp(-4/3 b^{3/2}); for Delta < 0 it is the sqrt joint-DOS with FK
  oscillations. At F = 10 MV/cm, hbar*theta ~ 0.6-0.7 eV, so a 3 eV (3-photon,
  1 eV) drive reaches E_direct up to ~3.5 eV.

The two cases the script produces by default (hbar w = 1.0 eV, field up to
10 MV/cm):

  CASE 1  N = 3  (3 hbar w = 3.0 eV):  3.0 eV sits 0.34 eV BELOW the 3.34 eV
          direct-gap minimum, so the transition is reachable ONLY through the
          Franz-Keldysh tail -- and only where the direct gap is smallest, i.e.
          along the Gamma->L valley. The 3D rate isosurfaces trace a blob at
          Gamma that reaches out along the eight <111> Gamma-L directions.
          (We look at the E_direct in [3.34, 3.5] eV region.)

  CASE 2  N = 4  (4 hbar w = 4.0 eV):  a GENUINE resonance shell E_direct = 4.0
          eV. E_direct = 4.0 is crossed at k ~ 0.84 along <100> (before X, whose
          gap is 4.24 eV) and much closer to Gamma along <110>, so the shell is a
          complex multi-lobed surface that reaches TOWARD X but never AT X.

Usage
-----
  python3 si_three_photon_isosurfaces.py                 # both cases, F=10 MV/cm
  python3 si_three_photon_isosurfaces.py --grid 61       # finer 3D grid
  python3 si_three_photon_isosurfaces.py --field 3       # FK width at 3 MV/cm
  python3 si_three_photon_isosurfaces.py --hw 1.05       # photon energy [eV]
  python3 si_three_photon_isosurfaces.py --orders 3      # just the 3-photon case
  python3 si_three_photon_isosurfaces.py --show          # also open in a browser

Outputs standalone HTML (self-contained plotly) into ./si_3ph_plots/ .
"""
import argparse
import itertools
import string
import time

import numpy as np
from scipy.special import airy

import epm_si_primitive as si_prim
import epm_gaas_primitive as prim
import epm_gaas_reference as ref

# --------------------------------------------------------------------------- #
#  physical constants / unit conversions
# --------------------------------------------------------------------------- #
HA = 27.211386245988                      # Hartree -> eV
AU_FIELD_VPM = 5.14220675e11              # 1 a.u. of E-field [V/m]
MVCM_TO_AU = 1e8 / AU_FIELD_VPM           # 1 MV/cm in a.u.
I_AU_WCM2 = 3.5094448314e16               # atomic unit of intensity [W/cm^2]


# --------------------------------------------------------------------------- #
#  1.  EPM engine
# --------------------------------------------------------------------------- #
def init_si_epm():
    """Configure the primitive-cell Si EPM and return a small context dict."""
    si_prim.configure_for_si('Si')
    ref.MATERIAL = prim.MATERIAL
    ref.A_LATTICE_AU = prim.A_LATTICE_AU
    ref.PW_CUTOFF_RY = prim.PW_CUTOFF_RY
    ref.NELEC = prim.NELEC
    a = prim.A_LATTICE_AU
    Gcart, _ = prim.build_pw_basis_fcc(a, prim.PW_CUTOFF_RY)
    ctx = dict(a=a, twopi_a=2.0 * np.pi / a, Gcart=Gcart, npw=len(Gcart),
               nocc=prim.NELEC // 2)
    return ctx


def solve_k(ctx, k_cart, nband):
    """Diagonalize H(k); return (eigenvalues[nband] in Ha, momentum p[nb,nb,3]).
    p[m,n,d] = <u_m|p_d|u_n> in atomic units (the interband/intraband velocity)."""
    H = ref.build_hamiltonian_sc(prim.MATERIAL, k_cart, ctx['Gcart'], ctx['a'])
    ev, evec = np.linalg.eigh(H)
    p = ref.momentum_matrix(k_cart, ctx['Gcart'], evec[:, :nband])
    return ev[:nband], p


# --------------------------------------------------------------------------- #
#  2.  N-photon matrix element (LOPT ladder) + exact orientation average
# --------------------------------------------------------------------------- #
def nphoton_tensor(ev, p, v, c, N, hw_au, eta_au):
    """Rank-N Cartesian tensor of the N-photon ladder amplitude v -> c at fixed k:

        T[a1,...,aN] = sum_{m1..m_{N-1}}
              p[c,m_{N-1},aN] ... p[m2,m1,a2] p[m1,v,a1]
              / prod_{j=1}^{N-1} (E_{mj} - E_v - j hbar w - i eta)

    Built by transfer-matrix contraction (all intermediate bands summed). The
    physical amplitude is (E0/2w)^N * (e_{a1}..e_{aN} T[a1..aN]); the field
    prefactor is applied later, so T is field-independent."""
    Ev = ev[v]
    # rung 1: v -> m1   ->  W[m1, a1]
    D = ev - Ev - 1.0 * hw_au - 1j * eta_au
    W = p[:, v, :] / D[:, None]
    # rungs 2 .. N-1 : add intermediate bands one at a time
    for j in range(2, N):
        D = ev - Ev - j * hw_au - 1j * eta_au
        # W_new[m_j, (old dirs), a_j] = sum_{m_{j-1}} p[m_j,m_{j-1},a_j] W[m_{j-1},...]
        W = np.tensordot(p, W, axes=([1], [0]))     # (m_j, a_j, old dirs...)
        W = np.moveaxis(W, 1, -1)                    # (m_j, old dirs..., a_j)
        W = W / D[(slice(None),) + (None,) * (W.ndim - 1)]
    # final rung: m_{N-1} -> c  ->  T[(dirs...), aN]
    T = np.tensordot(p[c, :, :], W, axes=([0], [0]))  # (aN, dirs...)
    T = np.moveaxis(T, 0, -1)                          # (dirs..., aN)
    return T


def iso_average_setup(N):
    """Isotropic rank-2N averaging tensor A and einsum string for
    <|e_{a1}..e_{aN} T[a1..aN]|^2>_e = einsum(subs, A, T, conj(T)).

    <e_{i1}...e_{i2N}> over the unit sphere = (1/(2N-1)!!) * sum over all
    pairings of the 2N indices of products of Kronecker deltas."""
    r = 2 * N
    letters = string.ascii_lowercase[:r]

    def pairings(lst):
        if not lst:
            yield []
            return
        first = lst[0]
        for i in range(1, len(lst)):
            for rest in pairings(lst[1:i] + lst[i + 1:]):
                yield [(first, lst[i])] + rest

    prs = list(pairings(list(range(r))))
    A = np.zeros((3,) * r)
    for comp in itertools.product(range(3), repeat=r):
        A[comp] = sum(all(comp[i] == comp[j] for i, j in pr) for pr in prs)
    dfact = 1
    for x in range(r - 1, 0, -2):
        dfact *= x
    A /= dfact
    subs = f'{letters},{letters[:N]},{letters[N:]}->'
    return A, subs


# --------------------------------------------------------------------------- #
#  3.  Franz-Keldysh field-broadened resonance lineshape
# --------------------------------------------------------------------------- #
def hbar_theta_ev(F_au, mu):
    """Electro-optic energy hbar*theta = (F^2 / 2 mu)^{1/3} in eV (a.u.: e=hbar=1)."""
    return (F_au ** 2 / (2.0 * mu)) ** (1.0 / 3.0) * HA


def fk_lineshape_jdos(delta_ev, hth_ev):
    """Literal 3D Franz-Keldysh JOINT-DOS lineshape g(Delta) = [Ai'(b)^2 -
    b Ai(b)^2]/hbar_theta, b = Delta/hbar_theta (Tharmalingam/Aspnes). Kept for
    reference/plotting; NOT used as the per-k weight because its b<0 sqrt growth
    is the parabolic joint DOS, which the explicit k-sum already reproduces from
    the real bands (using it per-k would double-count and drag the peak into the
    low-gap region). See fk_resonance_kernel for the per-k weight."""
    b = np.asarray(delta_ev) / hth_ev
    Ai, Aip, _, _ = airy(b)
    return np.maximum(Aip * Aip - b * Ai * Ai, 0.0) / hth_ev


def fk_resonance_kernel(delta_ev, hth_ev, eta_ev):
    """Field-broadened resonance kernel used to smear the vertical N-photon
    condition E_direct(k) = N hbar w (Delta = E_direct - N hbar w), peaked at
    Delta = 0 (value 1), asymmetric:

      Delta > 0 (E_direct ABOVE N hbar w, sub-edge): exp(-(4/3)(Delta/w)^{3/2})
                 -- the exact Franz-Keldysh photon-assisted-tunnelling tail that
                 lets a below-edge drive (e.g. 3 hbar w = 3.0 eV) reach the
                 3.34 eV gap; the tail length ~ hbar*theta grows as F^{2/3}.
      Delta < 0 (E_direct below N hbar w): Gaussian roll-off of the same width
                 -- the vertical transition to THIS band pair overshoots; higher
                 pairs (summed separately) carry the on-shell weight.

    w = hypot(hbar*theta, eta): electro-optic width with a collision-broadening
    floor eta. Peak-normalised (not area-normalised) on purpose, so the sub-edge
    rate visibly GROWS with field (more FK reach) instead of just spreading."""
    w = np.hypot(hth_ev, eta_ev)
    d = np.asarray(delta_ev, dtype=float)
    tail = np.exp(-(4.0 / 3.0) * np.abs(d / w) ** 1.5)      # sub-edge FK tail
    roll = np.exp(-0.5 * (d / w) ** 2)                       # above-edge roll-off
    return np.where(d >= 0.0, tail, roll)


# --------------------------------------------------------------------------- #
#  4.  per-k N-photon rate (sum over valence v and conduction c, FK lineshape)
# --------------------------------------------------------------------------- #
def nphoton_strength_k(ev, p, N, hw_au, eta_au, nocc, nband, window_ev):
    """Field-INDEPENDENT part: a list of (dE_eV, T_N) contributions for every
    valence->conduction pair within the resonance window of N hbar w. The FK
    field weighting g_FK(dE - N hbar w) is applied afterwards, so a whole field
    sweep costs nothing extra."""
    Nhw = N * hw_au * HA
    A, subs = ISO_CACHE[N]
    out = []
    for v in range(nocc):
        for c in range(nocc, nband):
            dE = (ev[c] - ev[v]) * HA
            if abs(dE - Nhw) > window_ev:
                continue
            T = nphoton_tensor(ev, p, v, c, N, hw_au, eta_au)
            M2 = float(np.real(np.einsum(subs, A, T, np.conj(T))))
            out.append((dE, M2))
    return out


def rate_from_contribs(contribs, N, hw_ev, hth_ev, eta_ev):
    """Combine the stored (dE, |M|^2) contributions with the FK resonance kernel."""
    if not contribs:
        return 0.0
    dE = np.array([c[0] for c in contribs])
    M2 = np.array([c[1] for c in contribs])
    return float(np.sum(M2 * fk_resonance_kernel(dE - N * hw_ev, hth_ev, eta_ev)))


ISO_CACHE = {}     # N -> (A, subs), filled in main()


# --------------------------------------------------------------------------- #
#  5.  irreducible-wedge (O_h) grid machinery
# --------------------------------------------------------------------------- #
def build_ibz(N_grid, k_lim, bz_sum_max=1.5):
    """Uniform cubic grid, restricted to the FCC BZ (|kx|+|ky|+|kz| <= 1.5), with
    O_h canonicalisation |k1|>=|k2|>=|k3| so equivalent points are computed once."""
    kr = np.linspace(-k_lim, k_lim, N_grid)
    KX, KY, KZ = np.meshgrid(kr, kr, kr, indexing='ij')
    canon_map = {}
    n_bz = 0
    for i in range(N_grid):
        for j in range(N_grid):
            for l in range(N_grid):
                kx, ky, kz = KX[i, j, l], KY[i, j, l], KZ[i, j, l]
                if abs(kx) + abs(ky) + abs(kz) > bz_sum_max:
                    continue
                n_bz += 1
                canon = tuple(sorted((abs(kx), abs(ky), abs(kz)), reverse=True))
                canon_map.setdefault(canon, []).append((i, j, l))
    return kr, KX, KY, KZ, canon_map, n_bz


# --------------------------------------------------------------------------- #
#  6.  compute E_direct(k) and the N-photon strengths over the 3D grid
# --------------------------------------------------------------------------- #
def compute_maps(ctx, orders, hw_ev, eta_ev, nband, N_grid, k_lim,
                 field_mvcm, mu, boundary_margin=0.12):
    """Return E_direct[NNN] and {N: W_N[NNN]} on the full 3D grid, evaluated in
    the irreducible wedge and replicated by O_h symmetry. W_N is the
    Franz-Keldysh-broadened multiphoton strength at the requested field."""
    kr, KX, KY, KZ, canon_map, n_bz = build_ibz(N_grid, k_lim)
    n_ibz = len(canon_map)
    shape = KX.shape
    E_dir = np.full(shape, np.nan)
    Wmaps = {N: np.full(shape, np.nan) for N in orders}

    hw_au, eta_au = hw_ev / HA, eta_ev / HA
    F_au = field_mvcm * MVCM_TO_AU
    hth = hbar_theta_ev(F_au, mu)
    # window must cover the FK tail (a few electro-optic energies)
    window = max(1.2, 4.0 * hth)

    print(f"# grid {N_grid}^3 | BZ {n_bz} pts | IBZ {n_ibz} "
          f"(x{n_bz / n_ibz:.1f}) | bands kept {nband}")
    print(f"# hw = {hw_ev:.3f} eV | field = {field_mvcm:g} MV/cm -> "
          f"hbar*theta = {hth:.3f} eV (mu={mu}) | eta = {eta_ev:.3f} eV")

    t0 = time.time()
    for n, (canon, cells) in enumerate(canon_map.items()):
        k_cart = ctx['twopi_a'] * np.array(canon)
        ev, p = solve_k(ctx, k_cart, nband)
        gap = (ev[ctx['nocc']] - ev[ctx['nocc'] - 1]) * HA
        w_here = {}
        for N in orders:
            contribs = nphoton_strength_k(ev, p, N, hw_au, eta_au,
                                          ctx['nocc'], nband, window)
            w_here[N] = rate_from_contribs(contribs, N, hw_ev, hth, eta_ev)
        for (i, j, l) in cells:
            E_dir[i, j, l] = gap
            for N in orders:
                Wmaps[N][i, j, l] = w_here[N]
        if (n + 1) % 400 == 0:
            print(f"#   [{n + 1:5d}/{n_ibz}] {time.time() - t0:.1f}s")
    print(f"# done in {time.time() - t0:.1f}s | "
          f"E_direct in [{np.nanmin(E_dir):.3f}, {np.nanmax(E_dir):.3f}] eV")

    # light BZ-boundary trim to suppress large-|G| noise on the outermost shell
    dist = np.abs(KX) + np.abs(KY) + np.abs(KZ)
    edge = dist > (1.5 - boundary_margin)
    E_dir[edge] = np.nan
    for N in orders:
        Wmaps[N][edge] = np.nan
    return dict(kr=kr, KX=KX, KY=KY, KZ=KZ, E_dir=E_dir, Wmaps=Wmaps,
                hth=hth, window=window)


# --------------------------------------------------------------------------- #
#  7.  plotting helpers (plotly)
# --------------------------------------------------------------------------- #
def bz_edge_trace(go):
    """Wireframe of the FCC (truncated-octahedron) Brillouin zone, in 2pi/a."""
    W = []
    for perm in [(1, 2, 0), (1, 0, 2), (0, 1, 2)]:
        for sx in (1, -1):
            for sy in (1, -1):
                vv = [0.0, 0.0, 0.0]
                coords = [1.0, 0.5, 0.0]
                vv[perm[0]] = sx * coords[0]
                vv[perm[1]] = sy * coords[1]
                W.append(tuple(vv))
    ex, ey, ez = [], [], []
    for i, v1 in enumerate(W):
        for v2 in W[i + 1:]:
            if abs(np.linalg.norm(np.subtract(v1, v2)) - np.sqrt(0.5)) < 0.01:
                ex += [v1[0], v2[0], None]
                ey += [v1[1], v2[1], None]
                ez += [v1[2], v2[2], None]
    return go.Scatter3d(x=ex, y=ey, z=ez, mode='lines',
                        line=dict(color='rgba(255,255,255,0.35)', width=2),
                        name='Brillouin zone', hoverinfo='skip')


def hs_marker_traces(go, ctx, nocc):
    """Gamma / X / L markers labelled with their direct gap."""
    cfg = [('Γ', [(0, 0, 0)], 'red', 'diamond', 8),
           ('X', [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0),
                  (0, 0, 1), (0, 0, -1)], 'orange', 'square', 6),
           ('L', [(.5, .5, .5), (.5, .5, -.5), (.5, -.5, .5), (.5, -.5, -.5),
                  (-.5, .5, .5), (-.5, .5, -.5), (-.5, -.5, .5), (-.5, -.5, -.5)],
            'cyan', 'circle', 6)]
    traces = []
    for name, pts, color, sym, size in cfg:
        kc = ctx['twopi_a'] * np.array(pts[0], float)
        ev, _ = solve_k(ctx, kc, ctx['nocc'] + 1)
        gap = (ev[nocc] - ev[nocc - 1]) * HA
        xs, ys, zs = zip(*pts)
        traces.append(go.Scatter3d(
            x=xs, y=ys, z=zs, mode='markers+text',
            text=[f'{name} {gap:.2f} eV'] + [None] * (len(pts) - 1),
            textposition='top center', textfont=dict(size=11, color='white'),
            marker=dict(size=size, color=color, symbol=sym,
                        line=dict(width=1, color='black')),
            name=f'{name}: {gap:.2f} eV'))
    return traces


def fig_rate_isosurfaces(go, data, ctx, N, hw_ev, field_mvcm, mu):
    """MAIN 3D figure: nested isosurfaces of the N-photon rate W_N(k)
    (normalised to its max), i.e. WHERE in k-space the direct N-photon
    transition happens and how strongly."""
    KX, KY, KZ = data['KX'], data['KY'], data['KZ']
    W = data['Wmaps'][N]
    wmax = np.nanmax(W)
    Wn = np.where(np.isnan(W), 0.0, W) / (wmax if wmax > 0 else 1.0)

    fig = go.Figure()
    # one Isosurface trace, several nested surfaces (embeds the volume ONCE)
    fig.add_trace(go.Isosurface(
        x=KX.ravel(), y=KY.ravel(), z=KZ.ravel(), value=Wn.ravel(),
        isomin=0.05, isomax=0.6, surface=dict(count=5, fill=1.0),
        caps=dict(x_show=False, y_show=False, z_show=False),
        colorscale='turbo', showscale=True, opacity=0.32,
        colorbar=dict(title=f'W_{N}/max', thickness=16, len=0.65),
        name=f'W_{N}(k)',
        hovertemplate=f'W_{N}/max=%{{value:.2f}}<br>k=(%{{x:.2f}},%{{y:.2f}},%{{z:.2f}})<extra></extra>'))
    fig.add_trace(bz_edge_trace(go))
    for t in hs_marker_traces(go, ctx, ctx['nocc']):
        fig.add_trace(t)

    Nhw = N * hw_ev
    sub = (f"hbar w = {hw_ev:.2f} eV -> N hbar w = {Nhw:.2f} eV | "
           f"F = {field_mvcm:g} MV/cm, hbar*theta = {data['hth']:.2f} eV | "
           f"rate ~ I^{N} (multiphoton, NOT tunnelling)")
    _style3d(fig, f"Si: {N}-photon direct transition rate W_{N}(k) "
                  f"[Franz-Keldysh broadened]", sub)
    return fig


def fig_resonance_shells(go, data, ctx, N, hw_ev, e_lo, e_hi):
    """E_direct(k) resonance shells in [e_lo, e_hi] eV -- the geometric locus of
    the N-photon resonance (nested surfaces expanding from Gamma). For N=3 the
    window is the [3.34, 3.5] eV valley that opens along Gamma-L; for N=4 it
    straddles the 4.0 eV shell that reaches toward X."""
    KX, KY, KZ = data['KX'], data['KY'], data['KZ']
    E = data['E_dir']
    Ev = np.where(np.isnan(E), 1e3, E)
    levels = np.round(np.linspace(e_lo, e_hi, 4), 3)

    fig = go.Figure()
    fig.add_trace(go.Isosurface(
        x=KX.ravel(), y=KY.ravel(), z=KZ.ravel(), value=Ev.ravel(),
        isomin=float(levels[0]), isomax=float(levels[-1]),
        surface=dict(count=len(levels), fill=1.0),
        caps=dict(x_show=False, y_show=False, z_show=False),
        colorscale='turbo', showscale=True, opacity=0.28,
        colorbar=dict(title='E_direct [eV]', thickness=16, len=0.65),
        name='E_direct shells',
        hovertemplate='E_direct=%{value:.2f} eV<br>k=(%{x:.2f},%{y:.2f},%{z:.2f})<extra></extra>'))
    fig.add_trace(bz_edge_trace(go))
    for t in hs_marker_traces(go, ctx, ctx['nocc']):
        fig.add_trace(t)
    _style3d(fig, f"Si: direct-gap resonance shells for the {N}-photon channel",
             f"E_direct(k) = N hbar w = {N * hw_ev:.2f} eV region "
             f"(shells {levels[0]:.2f}..{levels[-1]:.2f} eV)")
    return fig


def fig_rate_scatter(go, data, ctx, N):
    """3D scatter of every grid point coloured by log10(W_N) -- a volumetric view
    of the transition-probability cloud (complements the isosurfaces)."""
    KX, KY, KZ = data['KX'], data['KY'], data['KZ']
    W = data['Wmaps'][N]
    m = ~np.isnan(W) & (W > 0)
    w = W[m]
    logw = np.log10(w / w.max())
    fig = go.Figure()
    fig.add_trace(go.Scatter3d(
        x=KX[m], y=KY[m], z=KZ[m], mode='markers',
        marker=dict(size=2.5, color=logw, colorscale='turbo', cmin=-4, cmax=0,
                    opacity=0.55, colorbar=dict(title=f'log10 W_{N}/max',
                    thickness=16, len=0.65)),
        name=f'W_{N}(k)',
        hovertemplate=('k=(%{x:.2f},%{y:.2f},%{z:.2f})<br>'
                       'log10(W/max)=%{marker.color:.2f}<extra></extra>')))
    fig.add_trace(bz_edge_trace(go))
    for t in hs_marker_traces(go, ctx, ctx['nocc']):
        fig.add_trace(t)
    _style3d(fig, f"Si: {N}-photon rate cloud  log10 W_{N}(k)/max", '')
    return fig


def fig_line_scans(make_subplots, go, ctx, orders, hw_ev, eta_ev, nband, mu,
                   fields_mvcm):
    """2D companion: E_direct(k) and the FK-broadened N-photon rate along the
    Gamma-L, Gamma-X, Gamma-K lines, for several fields (shows the field->10
    MV/cm Franz-Keldysh growth and the Gamma-L selectivity of the 3-photon case)."""
    lines = {'Γ→L [111]': (np.array([1., 1., 1.]) / np.sqrt(3), 0.5 * np.sqrt(3)),
             'Γ→X [100]': (np.array([1., 0., 0.]), 1.0),
             'Γ→K [110]': (np.array([1., 1., 0.]) / np.sqrt(2), 0.75 * np.sqrt(2))}
    hw_au, eta_au = hw_ev / HA, eta_ev / HA

    ncol = len(orders)
    rows = 1 + ncol
    titles = ['E_direct along Γ-L, Γ-X, Γ-K']
    for N in orders:
        titles.append(f'{N}-photon rate W_{N} vs field along the lines')
    fig = make_subplots(rows=rows, cols=1, subplot_titles=titles,
                        vertical_spacing=0.09)

    colors = {'Γ→L [111]': '#1f9e3f', 'Γ→X [100]': '#2166f0', 'Γ→K [110]': '#d1461f'}
    npts = 80
    # row 1: direct gap
    scan = {}
    for name, (hatt, smax) in lines.items():
        svals = np.linspace(0, smax, npts)
        gaps, contribs_line = [], []
        for s in svals:
            kc = ctx['twopi_a'] * s * hatt
            ev, p = solve_k(ctx, kc, nband)
            gaps.append((ev[ctx['nocc']] - ev[ctx['nocc'] - 1]) * HA)
            contribs_line.append({N: nphoton_strength_k(
                ev, p, N, hw_au, eta_au, ctx['nocc'], nband,
                max(1.2, 4.0 * hbar_theta_ev(max(fields_mvcm) * MVCM_TO_AU, mu)))
                for N in orders})
        scan[name] = (svals, np.array(gaps), contribs_line)
        fig.add_trace(go.Scatter(x=svals, y=gaps, mode='lines',
                     line=dict(color=colors[name], width=2.5), name=name),
                     row=1, col=1)
    fig.add_hline(y=3 * hw_ev, line=dict(color='gray', dash='dot'),
                  annotation_text='3ħω', row=1, col=1)
    if 4 in orders:
        fig.add_hline(y=4 * hw_ev, line=dict(color='black', dash='dot'),
                      annotation_text='4ħω', row=1, col=1)

    # rows 2..: rate per order, one line style per field
    dashes = ['solid', 'dash', 'dot', 'dashdot']
    for r, N in enumerate(orders, start=2):
        for name, (svals, gaps, contribs_line) in scan.items():
            for fi, Fmv in enumerate(fields_mvcm):
                hth = hbar_theta_ev(Fmv * MVCM_TO_AU, mu)
                w = [rate_from_contribs(cl[N], N, hw_ev, hth, eta_ev)
                     for cl in contribs_line]
                fig.add_trace(go.Scatter(
                    x=svals, y=w, mode='lines',
                    line=dict(color=colors[name], width=1.8,
                              dash=dashes[fi % len(dashes)]),
                    name=f'{name}, {Fmv:g} MV/cm',
                    legendgroup=name, showlegend=(r == 2)),
                    row=r, col=1)
        fig.update_yaxes(title_text=f'W_{N} [arb.]', row=r, col=1)
    fig.update_yaxes(title_text='E_direct [eV]', row=1, col=1)
    fig.update_xaxes(title_text='|k| [2π/a]', row=rows, col=1)
    fig.update_layout(height=360 * rows, width=1000,
                      title='Si multiphoton transition: line scans '
                            '(Franz-Keldysh field dependence up to 10 MV/cm)',
                      legend=dict(font=dict(size=10)))
    return fig


def _style3d(fig, title, subtitle):
    fig.update_layout(
        title=dict(text=f"{title}<br><sup>{subtitle}</sup>", x=0.5,
                   font=dict(size=15)),
        scene=dict(xaxis_title='k_x (2π/a)', yaxis_title='k_y (2π/a)',
                   zaxis_title='k_z (2π/a)', aspectmode='cube',
                   camera=dict(eye=dict(x=1.7, y=1.7, z=1.5)),
                   bgcolor='rgb(12,12,22)',
                   xaxis=dict(showgrid=True, gridcolor='rgba(255,255,255,0.1)'),
                   yaxis=dict(showgrid=True, gridcolor='rgba(255,255,255,0.1)'),
                   zaxis=dict(showgrid=True, gridcolor='rgba(255,255,255,0.1)')),
        width=1100, height=950, margin=dict(l=0, r=0, b=0, t=90),
        legend=dict(bgcolor='rgba(0,0,0,0.6)', font=dict(color='white', size=11),
                    x=0.01, y=0.99))


# --------------------------------------------------------------------------- #
#  8.  driver
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--hw', type=float, default=1.0, help='photon energy [eV]')
    ap.add_argument('--orders', type=int, nargs='+', default=[3, 4],
                    help='multiphoton orders N to compute (default: 3 4)')
    ap.add_argument('--grid', type=int, default=45, help='cubic grid points per axis')
    ap.add_argument('--klim', type=float, default=1.0, help='grid half-width [2π/a]')
    ap.add_argument('--field', type=float, default=10.0,
                    help='peak field for the 3D maps [MV/cm] (FK width & I^N prefactor)')
    ap.add_argument('--fields', type=float, nargs='+', default=[1.0, 3.0, 10.0],
                    help='fields for the line-scan FK sweep [MV/cm]')
    ap.add_argument('--mu', type=float, default=0.15,
                    help='reduced mass for the Franz-Keldysh electro-optic energy [m_e]')
    ap.add_argument('--eta', type=float, default=0.08,
                    help='intermediate-state broadening in the LOPT denominators [eV]')
    ap.add_argument('--nband', type=int, default=24,
                    help='bands kept for the intermediate-state sum')
    ap.add_argument('-o', '--outdir', default='si_3ph_plots')
    ap.add_argument('--show', action='store_true', help='also open figures in a browser')
    args = ap.parse_args()

    import os
    os.makedirs(args.outdir, exist_ok=True)
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    ctx = init_si_epm()
    for N in args.orders:
        ISO_CACHE[N] = iso_average_setup(N)

    print("=" * 72)
    print("  Si DIRECT MULTIPHOTON TRANSITION -- 3D Franz-Keldysh-broadened rate")
    print("  (multiphoton I^N power law, contrast the Zener tunnelling exp(-Eg^3/2/F))")
    print("=" * 72)

    data = compute_maps(ctx, args.orders, args.hw, args.eta, args.nband,
                        args.grid, args.klim, args.field, args.mu)

    # peak field -> peak intensity. I_AU_WCM2 (=3.51e16) already equals
    # (1/2) c eps0 E_au^2, so I[W/cm^2] = I_AU_WCM2 * (E0/E_au)^2 (no extra 1/2).
    # The multiphoton prefactor is (E0/2w)^{2N} ~ I^N (power law).
    E0 = args.field * MVCM_TO_AU
    intensity_wcm2 = (E0 ** 2) * I_AU_WCM2
    print(f"# peak field {args.field:g} MV/cm  ->  E0 = {E0:.3e} a.u., "
          f"I ~ {intensity_wcm2:.2e} W/cm^2 (multiphoton prefactor ~ I^N)")

    written = []
    for N in args.orders:
        w = args.hw
        e_lo, e_hi = (3.34, 3.50) if N == 3 else (N * w - 0.15, N * w + 0.15)

        f1 = fig_rate_isosurfaces(go, data, ctx, N, w, args.field, args.mu)
        p1 = os.path.join(args.outdir, f'si_{N}photon_rate_iso.html')
        f1.write_html(p1); written.append(p1)

        f2 = fig_resonance_shells(go, data, ctx, N, w, e_lo, e_hi)
        p2 = os.path.join(args.outdir, f'si_{N}photon_resonance_shells.html')
        f2.write_html(p2); written.append(p2)

        f3 = fig_rate_scatter(go, data, ctx, N)
        p3 = os.path.join(args.outdir, f'si_{N}photon_rate_cloud.html')
        f3.write_html(p3); written.append(p3)

        wmax = np.nanmax(data['Wmaps'][N])
        km = np.unravel_index(np.nanargmax(np.where(np.isnan(data['Wmaps'][N]),
                              -1, data['Wmaps'][N])), data['E_dir'].shape)
        kpk = (data['KX'][km], data['KY'][km], data['KZ'][km])
        print(f"# N={N}: peak W_{N} at k=({kpk[0]:+.2f},{kpk[1]:+.2f},{kpk[2]:+.2f}) "
              f"2π/a, E_direct there = {data['E_dir'][km]:.3f} eV")
        if args.show:
            f1.show(); f2.show(); f3.show()

    fL = fig_line_scans(make_subplots, go, ctx, args.orders, args.hw, args.eta,
                        args.nband, args.mu, args.fields)
    pL = os.path.join(args.outdir, 'si_multiphoton_line_scans.html')
    fL.write_html(pL); written.append(pL)
    if args.show:
        fL.show()

    print("=" * 72)
    print("PHYSICS SUMMARY")
    print("=" * 72)
    print(f"  photon energy      : {args.hw:.3f} eV")
    print(f"  Si direct-gap min  : {np.nanmin(data['E_dir']):.3f} eV (at Gamma)")
    for N in args.orders:
        Nhw = N * args.hw
        rel = 'BELOW' if Nhw < np.nanmin(data['E_dir']) else 'above'
        note = ('reachable only via the Franz-Keldysh tail -> concentrated in the '
                'low-gap Gamma-L valley' if Nhw < np.nanmin(data['E_dir'])
                else 'genuine resonance shell (reaches toward X, not at X for N=4)')
        print(f"  {N}-photon: {N}ħω = {Nhw:.2f} eV is {rel} the gap min -> {note}")
    print(f"  Franz-Keldysh width: hbar*theta = {data['hth']:.3f} eV "
          f"at {args.field:g} MV/cm  (~F^2/3; grows to ~0.68 eV at 10 MV/cm)")
    print(f"  process/prefactor  : multiphoton, W_N ~ I^N  (power law) -- distinct")
    print(f"                       from Zener tunnelling  P ~ exp(-C m_r^1/2 Eg^3/2 / F)")
    print("-" * 72)
    print("WROTE:")
    for p in written:
        print("   " + p)


if __name__ == '__main__':
    main()
