#!/usr/bin/env python3
"""
Case 3.1.2 validation: preloaded-slab diffusion in an isolated domain.

Analytical solution (Fourier cosine series for isolated slab)
-------------------------------------------------------------
IC:  C(x, 0) = C0  for x ∈ [0, h],  C(x, 0) = 0  otherwise
BCs: ∂C/∂x = 0  at  x = 0  and  x = L    (no-flux, mass conserved)

Analytical:
    C(x, t) = C0*h/L
            + (2*C0/π) Σ_{n=1}^N  sin(nπh/L)/n · cos(nπx/L) · exp(-(nπ/L)²·D·t)

Parameters (matching case input files):
    L   = 0.1 m       (slab length)
    h   = 0.01 m      (preloaded depth)
    D   = 1e-6 m²/s   (diffusivity)
    C0  = 1 mol/m³    (initial concentration in preloaded region)
    N   = 100 cells   (dx = 1 mm)
"""

import os
import sys
import numpy as np

# ---------- case parameters ------------------------------------------------
L   = 0.1       # m
h   = 0.01      # m
D   = 1e-6      # m²/s
C0  = 1.0       # mol/m³
N_cells = 100   # number of cells


def C_analytical(x, t, N_terms=500):
    """Fourier cosine series for isolated-slab preloaded IC."""
    result = np.full_like(x, C0 * h / L)
    for n in range(1, N_terms + 1):
        coeff = (2.0 * C0 / (n * np.pi)) * np.sin(n * np.pi * h / L)
        result += coeff * np.cos(n * np.pi * x / L) * np.exp(-(n * np.pi / L)**2 * D * t)
    return result


def read_of_field(path):
    """Parse an OpenFOAM ASCII volScalarField and return the internalField values.

    Uses a three-state machine so that integer-valued data are never confused
    with the cell-count line that precedes the opening '('.
    """
    state = "scanning"
    values = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if state == "scanning":
                if "internalField" in s:
                    # Check nonuniform FIRST — "uniform" is a substring of "nonuniform"
                    if "nonuniform" in s:
                        state = "header"
                    elif "uniform" in s:
                        try:
                            val = float(s.split()[-1].rstrip(";"))
                            return np.full(N_cells, val)
                        except ValueError:
                            pass
            elif state == "header":
                if s == "(":
                    state = "data"
            elif state == "data":
                if s == ");":
                    break
                try:
                    values.append(float(s))
                except ValueError:
                    pass
    return np.array(values)


def check_time(case_dir, time_str, tol_rel=0.02):
    path = os.path.join(case_dir, time_str, "slab", "C_T2")
    if not os.path.isfile(path):
        print(f"  t={time_str}: field not found at {path}")
        return None

    t = float(time_str)
    C_sim = read_of_field(path)
    if len(C_sim) == 0:
        print(f"  t={time_str}: could not parse field values")
        return None

    n = len(C_sim)
    dx = L / N_cells
    x = np.linspace(dx / 2, L - dx / 2, n)   # cell centres

    C_ana = C_analytical(x, t)

    # Mass conservation check
    mass_sim = np.sum(C_sim) * dx
    mass_ana = C0 * h
    mass_err = abs(mass_sim - mass_ana) / mass_ana

    # Profile error (normalise by max of analytical to avoid div-by-zero near edges)
    C_scale = C_ana.max() + 1e-20
    err = np.abs(C_sim - C_ana) / C_scale
    max_err = err.max()

    ok = max_err < tol_rel and mass_err < 0.001
    print(
        f"  t={time_str:>6} s:  "
        f"max profile error = {max_err:.3e}  "
        f"mass error = {mass_err:.3e}  "
        + ("✓" if ok else f"✗ (profile tol {tol_rel*100:.0f}%)")
    )
    return max_err


def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--case", default=".")
    p.add_argument("--times", nargs="+", default=["10", "50", "100"])
    args = p.parse_args()

    case_dir = os.path.abspath(args.case)
    print(f"Case directory: {case_dir}")
    print()
    print("Analytical parameters:")
    print(f"  L   = {L} m,  h = {h} m,  D = {D} m²/s,  C0 = {C0} mol/m³")
    print(f"  C_uniform (t→∞) = {C0*h/L:.4f} mol/m³")
    print()

    results = {}
    for t_str in args.times:
        err = check_time(case_dir, t_str)
        if err is not None:
            results[t_str] = err

    if not results:
        print("No data found — run Allrun first.")
        sys.exit(1)

    all_pass = all(v < 0.02 for v in results.values())
    print()
    print("PASS" if all_pass else "FAIL")
    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
