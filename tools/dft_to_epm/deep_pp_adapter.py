#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Adapter that pulls the vendored DeePseudopot local-pseudopotential model into
the SALMON dft_to_epm extractor as a module.

The "from-folder" (``--method zunger``) extraction mode reuses the **analytic
Zunger local form factor** that DeePseudopot fits to ab-initio band structures:

    V(q) = a0 (q^2 - a1) / (a2 * exp(a3 q^2) - 1)          [Wang & Zunger,
                                                             PRB 51, 17398 (1995)]

implemented upstream in ``external/DeePseudopot/utils/pp_func.py::pot_func``.
That module imports ``torch`` (and matplotlib); when the full DeePseudopot stack
is installed we call the upstream function directly, so the extractor and the ML
pseudopotential project share *the same code path*. When torch is not present we
fall back to a NumPy evaluation of the *identical* closed form, so the mode still
works in a lightweight environment.

``q`` is |G| in atomic units (Bohr^-1); ``params = [a0, a1, a2, a3]``. The return
value is the screened atomic form factor in the EPM matrix-element convention
(the dft_to_epm forward model applies the same Ry->Ha 1/2 factor as
``cb_get_form_factors``, so the overall scale is fixed by the fit).
"""
import os
import sys

import numpy as np

_VENDOR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), 'external', 'DeePseudopot')

BACKEND = 'numpy-fallback'
_pot_func = None
_torch = None

try:  # prefer the vendored DeePseudopot implementation (module integration)
    if _VENDOR not in sys.path:
        sys.path.insert(0, _VENDOR)
    import torch as _torch  # noqa: F401
    from utils.pp_func import pot_func as _pot_func  # type: ignore
    BACKEND = 'deepseudopot(vendored)'
except Exception:  # torch/matplotlib unavailable -> faithful NumPy port
    _pot_func = None
    _torch = None


def zunger_form_factor(q_au, params):
    """Zunger analytic local form factor V(q) for |G| = ``q_au`` (Bohr^-1).

    ``q_au`` may be a scalar or array; ``params`` is ``[a0, a1, a2, a3]``.
    Uses the vendored ``pot_func`` when DeePseudopot (torch) is importable,
    otherwise an identical NumPy evaluation.
    """
    q = np.atleast_1d(np.asarray(q_au, dtype=float))
    if _pot_func is not None:
        with _torch.no_grad():
            qt = _torch.as_tensor(q, dtype=_torch.float64)
            pt = _torch.as_tensor(np.asarray(params, dtype=float),
                                  dtype=_torch.float64)
            out = _pot_func(qt, pt).cpu().numpy()
    else:
        a0, a1, a2, a3 = (float(params[0]), float(params[1]),
                          float(params[2]), float(params[3]))
        out = a0 * (q * q - a1) / (a2 * np.exp(a3 * q * q) - 1.0)
    return out if np.ndim(q_au) else float(out[0])


def backend_name():
    return BACKEND
