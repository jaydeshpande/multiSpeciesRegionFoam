#!/usr/bin/env python3
"""
Case 3.1.6 validation: thermally coupled membrane distillation.

Steady-state analytical solution (variable D)
----------------------------------------------
At steady state the temperature profile is linear:
    T_ss(x) = T_hot + (T_cold − T_hot) · x / L

and the effective pore diffusivity follows Arrhenius:
    D(x) = D0 · exp(−Ea / (R · T_ss(x)))

Constant-flux condition (no accumulation at steady state):
    J · dx = −D(x) · dC   →   C(x) = C_hot − J · I(x)

where the resistance integral  I(x) = ∫_0^x ds / D(T_ss(s)).

Applying the right BC:
    C_cold = C_hot − J · I(L)   →   J = (C_hot − C_cold) / I(L)

Both I(L) and I(x) are evaluated by trapezoidal quadrature (N_quad points).

This produces a profile that curves toward the cold side: D is lower there,
so more concentration drop occurs per unit length — the opposite of the linear
profile from case 3.1.5 (isothermal, constant D).
"""

import os
import sys
import math
import numpy as np

# ---------- case parameters ------------------------------------------------
L       = 1e-3      # m  (membrane thickness)
T_hot   = 333.0     # K
T_cold  = 293.0     # K
C_hot   = 7.0       # mol/m³
C_cold  = 1.0       # mol/m³
D0      = 4.67e-6   # m²/s  (Arrhenius pre-exponential)
Ea      = 10000.0   # J/mol
R_gas   = 8.314462  # J/(mol·K)
N_cells = 100
N_quad  = 10000     # quadrature points for analytical solution


def T_ss(x):
    return T_hot + (T_cold - T_hot) * x / L


def D_eff(T):
    return D0 * math.exp(-Ea / (R_gas * T))


def D_eff_arr(x):
    """Vectorised D_eff over position array x."""
    T = T_hot + (T_cold - T_hot) * x / L
    return D0 * np.exp(-Ea / (R_gas * T))


def analytical_profile(x_eval):
    """
    Compute C_ss at positions x_eval using trapezoidal quadrature of I(x).
    """
    xs = np.linspace(0, L, N_quad)
    inv_D = 1.0 / D_eff_arr(xs)

    # cumulative resistance integral via trapz
    I_full = np.trapz(inv_D, xs)
    J_ss   = (C_hot - C_cold) / I_full

    C_out = np.empty_like(x_eval)
    for i, xi in enumerate(x_eval):
        mask = xs <= xi
        if mask.sum() < 2:
            C_out[i] = C_hot
        else:
            I_x = np.trapz(inv_D[mask], xs[mask])
            C_out[i] = C_hot - J_ss * I_x
    return C_out, J_ss


def read_of_field(path, n_cells):
    values = []
    in_list = False
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if "internalField" in s:
                if "uniform" in s:
                    try:
                        val = float(s.split()[-1].rstrip(";"))
                        return np.full(n_cells, val)
                    except ValueError:
                        pass
                elif "nonuniform" in s:
                    in_list = True
                continue
            if in_list:
                if s in ("(", ")"):
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
    return np.array(values)


def check_time(case_dir, time_str, tol_rel=0.02):
    # ---- species field ----
    path_C = os.path.join(case_dir, time_str, "membrane", "C_H2O")
    path_T = os.path.join(case_dir, time_str, "membrane", "T")

    if not os.path.isfile(path_C):
        print(f"  t={time_str} s: C_H2O not found at {path_C}")
        return None, None

    C_sim = read_of_field(path_C, N_cells)
    if len(C_sim) == 0:
        print(f"  t={time_str} s: could not parse C_H2O")
        return None, None

    dx = L / N_cells
    x  = np.linspace(dx / 2, L - dx / 2, len(C_sim))

    C_ana, J_ss = analytical_profile(x)
    err_C = np.abs(C_sim - C_ana) / (abs(C_hot - C_cold) + 1e-30)
    max_err_C = err_C.max()

    # ---- temperature field ----
    max_err_T = None
    if os.path.isfile(path_T):
        T_sim = read_of_field(path_T, N_cells)
        if len(T_sim) > 0:
            T_ana = T_hot + (T_cold - T_hot) * x / L
            err_T = np.abs(T_sim - T_ana) / abs(T_hot - T_cold)
            max_err_T = err_T.max()

    ok = max_err_C < tol_rel and (max_err_T is None or max_err_T < 0.005)
    print(
        f"  t={time_str:>6} s:  "
        f"max|C−C_ss|/ΔC = {max_err_C:.3e}  "
        + (f"max|T−T_ss|/ΔT = {max_err_T:.3e}  " if max_err_T is not None else "")
        + ("✓" if ok else f"✗")
    )
    return max_err_C, max_err_T


def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--case",  default=".")
    p.add_argument("--times", nargs="+", default=["10", "20"])
    args = p.parse_args()

    case_dir = os.path.abspath(args.case)

    D_hot  = D0 * math.exp(-Ea / (R_gas * T_hot))
    D_cold = D0 * math.exp(-Ea / (R_gas * T_cold))

    # Analytical J_ss and C_ss at a few representative points
    xs_rep = np.linspace(0, L, N_quad)
    inv_D  = 1.0 / D_eff_arr(xs_rep)
    I_L    = np.trapz(inv_D, xs_rep)
    J_ss   = (C_hot - C_cold) / I_L

    print(f"Case directory: {case_dir}")
    print()
    print("Arrhenius D_eff(T) parameters:")
    print(f"  D0 = {D0:.3e} m²/s,  Ea = {Ea:.0f} J/mol")
    print(f"  D_eff({T_hot:.0f} K) = {D_hot:.3e} m²/s   (hot face)")
    print(f"  D_eff({T_cold:.0f} K) = {D_cold:.3e} m²/s  (cold face)")
    print(f"  D ratio = {D_hot/D_cold:.2f}  →  curved C profile")
    print()
    print(f"Steady-state flux:  J_ss = {J_ss:.4e} mol/(m²·s)")
    print(f"  (isothermal case 315 J_ss = {1e-7*6/L:.4e} mol/(m²·s) for comparison)")
    print()

    results_C = {}
    for t_str in args.times:
        err_C, _ = check_time(case_dir, t_str)
        if err_C is not None:
            results_C[t_str] = err_C

    if not results_C:
        print("No data found — run Allrun first.")
        sys.exit(1)

    all_pass = all(v < 0.02 for v in results_C.values())
    print()
    print("PASS" if all_pass else "FAIL")
    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
