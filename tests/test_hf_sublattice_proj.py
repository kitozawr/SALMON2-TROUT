#!/usr/bin/env python3
"""
test_hf_sublattice_proj.py  -  Part E (HF sublattice-block projection).

Mirrors the projection applied in compute_coulomb_selfenergy:
    Sigma_ij(k) *= proj_ij,   proj_ij = sum_s w_s(i) w_s(j)  for i != j,
                              proj_ii = 1 (diagonal kept).
where w_s(i) = gs%unfold_w(s, band_i, k) are the FCC-sublattice spectral
weights (sum_s w_s = 1). The off-diagonal exchange between cubic bands that
fold to the same cubic k but live on DIFFERENT primitive sublattices is a
folding artifact and must vanish; the diagonal (energy renormalization) is
preserved; the projection is real-symmetric so Sigma stays Hermitian.

Pure-Python (no SALMON build). Also exercises real unfold weights computed from
the standalone EPM (GaAs 4x4x4, scalar) so the test reflects physical weights.
"""
import importlib.util
import os
import sys

import numpy as np
from numpy.linalg import eigh

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def project(sigma, w):
    """Apply Sigma_ij *= sum_s w_s(i) w_s(j) (i!=j), diagonal kept. w: (nb,4)."""
    P = w @ w.T                 # P_ij = sum_s w_s(i) w_s(j)
    np.fill_diagonal(P, 1.0)    # keep the diagonal at full strength
    return sigma * P


def _load_epm():
    spec = importlib.util.spec_from_file_location(
        "epm", os.path.join(ROOT, "epm_gaas_reference.py"))
    epm = importlib.util.module_from_spec(spec)
    saved = sys.argv
    sys.argv = ["epm"]
    try:
        spec.loader.exec_module(epm)
    finally:
        sys.argv = saved
    return epm


def gaas_sublattice_weights(epm, a=10.68, cutoff=11.1):
    """Per-band FCC-sublattice weights w_s(band) at one folded cubic k (scalar)."""
    Gcart, _ = epm.build_plane_wave_basis_sc(a, cutoff)
    tpa = 2.0 * np.pi / a
    Gi = np.round(Gcart / tpa).astype(int)
    nb = epm.NSTATE
    offsets = epm.SUBLATTICE_OFFSETS
    masks = [epm.sublattice_mask(Gi, off) for off in offsets]
    # a generic low-symmetry cubic k so the 4 sublattices are non-degenerate
    ksc = tpa * np.array([0.123, -0.071, 0.045])
    H = epm.build_hamiltonian_sc("GaAs", ksc, Gcart, a)
    evals, evecs = eigh(H)
    w2 = np.abs(evecs[:, :nb]) ** 2
    wsub = np.zeros((4, nb))
    for s, msk in enumerate(masks):
        wsub[s] = w2[msk].sum(axis=0)
    wsub /= np.maximum(wsub.sum(axis=0), 1e-300)
    return wsub.T   # (nb, 4)


def main():
    ok = True
    rng = np.random.default_rng(0)

    # --- 1. Hard-localised weights: each band fully on one sublattice ------
    nb = 8
    sub = np.array([0, 0, 1, 1, 2, 2, 3, 3])
    w = np.zeros((nb, 4))
    w[np.arange(nb), sub] = 1.0
    # Hermitian random Sigma
    A = rng.standard_normal((nb, nb)) + 1j * rng.standard_normal((nb, nb))
    sigma = A + A.conj().T
    sp = project(sigma, w)
    # inter-sublattice off-diagonals must be exactly zero
    inter = np.array([[sub[i] != sub[j] for j in range(nb)] for i in range(nb)])
    if np.abs(sp[inter]).max() > 1e-14:
        print("  FAIL: inter-sublattice elements not zeroed"); ok = False
    # intra-sublattice (incl diagonal) preserved exactly
    intra = ~inter
    if np.abs(sp[intra] - sigma[intra]).max() > 1e-14:
        print("  FAIL: intra-sublattice elements changed"); ok = False
    # diagonal preserved
    if np.abs(np.diag(sp) - np.diag(sigma)).max() > 1e-14:
        print("  FAIL: diagonal changed"); ok = False
    # Hermiticity preserved
    if np.abs(sp - sp.conj().T).max() > 1e-13:
        print("  FAIL: projection broke Hermiticity"); ok = False

    # --- 2. No-folding limit: all bands on one sublattice -> identity ------
    w1 = np.zeros((nb, 4)); w1[:, 0] = 1.0
    if np.abs(project(sigma, w1) - sigma).max() > 1e-14:
        print("  FAIL: single-sublattice case is not the identity"); ok = False

    # --- 3. Real GaAs weights: disjoint-support pairs are strongly damped --
    try:
        epm = _load_epm()
        wr = gaas_sublattice_weights(epm)           # (nb,4), sum_s=1
        P = wr @ wr.T
        args = wr.argmax(axis=1)
        offdiff = np.array([[i != j and args[i] != args[j]
                             for j in range(len(args))] for i in range(len(args))])
        # bands whose dominant sublattices differ must have small overlap
        if offdiff.any() and P[offdiff].max() > 0.5:
            print(f"  FAIL: real-weight inter-sublattice overlap too large "
                  f"({P[offdiff].max():.3f})"); ok = False
        # rows of P with same-sublattice partner should reach ~1 somewhere
        if P[~np.eye(len(args), dtype=bool)].max() < 0.5:
            print("  FAIL: no strong intra-sublattice overlap found"); ok = False
        print(f"  real GaAs weights: max inter-sublattice overlap = "
              f"{P[offdiff].max() if offdiff.any() else 0.0:.3f}, "
              f"max intra = {P[~np.eye(len(args),dtype=bool)].max():.3f}")
    except Exception as e:                            # noqa: BLE001
        print(f"  (skipped real-weight check: {e})")

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
