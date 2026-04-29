#!/usr/bin/env python3
"""
Case 3.1.5 validation: isothermal membrane distillation.

Analytical steady-state solution
---------------------------------
Water vapour diffuses through a porous membrane under a fixed concentration
gradient.  With constant D_eff and Dirichlet BCs, the steady profile is linear:

    C_ss(x) = C_hot + (C_cold − C_hot) · x / L
            = 7 − 6000 · x    [mol/m³,  x in m]

    J_ss = D_eff · (C_hot − C_cold) / L = 6e-7 mol/(m²·s)

Parameters (matching case input files):
    L      = 1e-3 m       (membrane thickness)
    D_eff  = 1e-7 m²/s   (effective pore diffusivity)
    C_hot  = 7.0 mol/m³  (hot-face saturation concentration)
    C_cold = 1.0 mol/m³  (cold-face saturation concentration)
    N      = 100 cells    (dx = 1e-5 m)
    τ_SS   = L²/D_eff = 10 s  →  steady state by t ≈ 50 s
"""

import os
import sys
import numpy as np

# ---------- case parameters ------------------------------------------------
L      = 1e-3     # m
D_eff  = 1e-7     # m²/s
C_hot  = 7.0      # mol/m³
C_cold = 1.0      # mol/m³
N_cells = 100


def C_analytical(x):
    return C_hot + (C_cold - C_hot) * x / L


def J_steady():
    return D_eff * (C_hot - C_cold) / L   # mol/(m²·s)


def read_of_field(path):
    """Return internal field values from an ASCII OpenFOAM volScalarField.

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


def check_time(case_dir, time_str, tol_rel=0.01):
    path = os.path.join(case_dir, time_str, "membrane", "C_H2O")
    if not os.path.isfile(path):
        print(f"  t={time_str} s: field not found at {path}")
        return None

    C_sim = read_of_field(path)
    if len(C_sim) == 0:
        print(f"  t={time_str} s: could not parse C_H2O")
        return None

    dx = L / N_cells
    x  = np.linspace(dx / 2, L - dx / 2, len(C_sim))
    C_ana = C_analytical(x)

    err = np.abs(C_sim - C_ana) / (abs(C_hot - C_cold) + 1e-30)
    max_err = err.max()

    ok = max_err < tol_rel
    print(
        f"  t={time_str:>6} s:  max|C−C_ss|/ΔC = {max_err:.3e}  "
        + ("✓" if ok else f"✗ (tol {tol_rel*100:.0f}%)")
    )
    return max_err


def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--case",  default=".")
    p.add_argument("--times", nargs="+", default=["50", "100"])
    args = p.parse_args()

    case_dir = os.path.abspath(args.case)
    tau_SS = L**2 / D_eff

    print(f"Case directory: {case_dir}")
    print()
    print("Analytical parameters:")
    print(f"  L       = {L} m,  D_eff = {D_eff} m²/s")
    print(f"  C_hot   = {C_hot} mol/m³,  C_cold = {C_cold} mol/m³")
    print(f"  τ_SS    = {tau_SS:.1f} s")
    print(f"  J_ss    = {J_steady():.3e} mol/(m²·s)")
    print()

    results = {}
    for t_str in args.times:
        err = check_time(case_dir, t_str)
        if err is not None:
            results[t_str] = err

    if not results:
        print("No data found — run Allrun first.")
        sys.exit(1)

    all_pass = all(v < 0.01 for v in results.values())
    print()
    print("PASS" if all_pass else "FAIL")
    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
