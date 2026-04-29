#!/usr/bin/env python3
"""
Case 3.1.4 validation: composite PyC/SiC membrane, Sieverts partition.

Analytical steady-state solution
---------------------------------
Parameters (matching case input files):
    D_pyc = D_sic = D = 1e-9 m²/s
    L_pyc = L_sic = 5e-4 m
    Ks_pyc = 2000 mol/m³,  Ks_sic = 1000 mol/m³  →  r = 2
    C_in  = 1 mol/m³  (pyc inlet)
    C_out = 0 mol/m³  (sic outlet)

At steady state (dC/dt = 0 → C linear in each layer):
    C_sic_if = C_in / (1 + r)   = 1/3  mol/m³   (SiC side of interface)
    C_pyc_if = r * C_sic_if     = 2/3  mol/m³   (PyC side of interface)
    J        = D * C_sic_if / L = 6.67e-7 mol/(m²·s)

    PyC:  C_pyc(x) = 1 - (1 - 2/3) / L * x         for x in [0, L_pyc]
    SiC:  C_sic(x) = (1/3) * (1 - (x - L_pyc) / L)  for x in [L_pyc, 2*L_pyc]
"""

import os
import sys
import numpy as np

# ---------- case parameters ------------------------------------------------
D     = 1e-9      # m²/s
L_pyc = 5e-4      # m
L_sic = 5e-4      # m
r     = 2.0       # Ks_pyc / Ks_sic
C_in  = 1.0       # mol/m³


def C_pyc_steady(x):
    C_if_pyc = r * C_in / (1.0 + r)   # = 2/3
    return C_in - (C_in - C_if_pyc) / L_pyc * x


def C_sic_steady(x):
    C_if_sic = C_in / (1.0 + r)       # = 1/3
    xi = x - L_pyc                    # local coordinate in SiC
    return C_if_sic * (1.0 - xi / L_sic)


def check_steady_state(case_dir, time_str, tol_rel=0.01):
    """Read OF fields at time_str and compare to analytical steady state."""
    print(f"\nChecking t = {time_str} s against steady-state analytical")

    results = {}
    for region, x_func, C_func, N_cells in [
        ("pyc",
         lambda: np.linspace(L_pyc / (2 * 100), L_pyc - L_pyc / (2 * 100), 100),
         C_pyc_steady,
         100),
        ("sic",
         lambda: L_pyc + np.linspace(L_sic / (2 * 100), L_sic - L_sic / (2 * 100), 100),
         C_sic_steady,
         100),
    ]:
        path = os.path.join(case_dir, time_str, region, "C_T2")
        if not os.path.isfile(path):
            print(f"  {region}: field not found at {path}")
            continue

        values = []
        in_list = False
        with open(path) as fh:
            for line in fh:
                s = line.strip()
                if "internalField" in s and "nonuniform" in s:
                    in_list = True
                    continue
                if in_list:
                    if s == "(":
                        continue
                    elif s == ");":
                        break
                    elif s.isdigit():
                        continue
                    else:
                        try:
                            values.append(float(s))
                        except ValueError:
                            pass

        if not values:
            print(f"  {region}: could not parse field values")
            continue

        x  = x_func()
        C  = np.array(values)
        Ca = C_func(x[:len(C)])

        err = np.abs(C - Ca) / (np.abs(Ca).max() + 1e-20)
        max_err = err.max()
        print(f"  {region}: max relative error = {max_err:.3e}", end="")
        print("  ✓" if max_err < tol_rel else f"  ✗ (>{tol_rel*100:.0f}%)")
        results[region] = max_err

    return results


def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--case", default=".")
    p.add_argument("--time", default="1000")
    args = p.parse_args()

    case_dir = os.path.abspath(args.case)
    print(f"Case directory: {case_dir}")
    print()
    print("Steady-state analytical values:")
    print(f"  C_pyc(inlet)      = {C_in:.4f} mol/m³")
    print(f"  C_pyc(interface)  = {r/(1+r):.4f} mol/m³")
    print(f"  C_sic(interface)  = {1/(1+r):.4f} mol/m³")
    print(f"  C_sic(outlet)     = 0.0000 mol/m³")
    print(f"  Flux J            = {D * (1/(1+r)) / L_sic:.3e} mol/(m²·s)")

    results = check_steady_state(case_dir, args.time)

    if not results:
        print("\nNo data found — run Allrun first.")
        sys.exit(1)

    all_pass = all(v < 0.01 for v in results.values())
    print()
    print("PASS" if all_pass else "FAIL")
    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
