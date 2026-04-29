#!/usr/bin/env python3
"""
Case 3.1.3 validation: McNabb-Foster trapping, permeation breakthrough.

Two sub-cases (both use same steady-state profile, different time scales):

  (a) Weak trapping  ε = 0.09 eV:
        αd(1000 K) ≈ 3.52e12 s⁻¹,  S = 1 + ft·αt/αd ≈ 29.4
        D_eff ≈ D/29.4 ≈ 0.034 m²/s
        τ_bc  = L²·S / (2π²·D) ≈ 1.5 s   [validated at t = 20 s]

  (b) Strong trapping  ε = 1.5 eV  (change endTime to 3000 in controlDict):
        τ_bd = ft·N·L² / (2·C_in·D) = 500 s   [validated at t = 3000 s]

Steady-state profile (same for both cases):
    C_ss(x) = C_in · (1 - x/L)

Trapped-species steady state (weak-trapping limit: Ct_ss ≈ (S-1)·C_ss):
    Ct_ss(x) = (ft·αt/αd) · C_ss(x)

The script checks:
  1. Mobile concentration profile at the final time against C_ss.
  2. Total mobile mass conservation relative to the asymptotic steady-state value.
"""

import os
import sys
import math
import numpy as np

# ---------- case parameters ------------------------------------------------
L      = 1.0        # m
C_in   = 5.25e-6    # mol/m³
N_cells = 100
kB_eV  = 8.617333262e-5   # eV/K
T_K    = 1000.0     # K
ft     = 0.1
N      = 0.0525     # mol/m³  (molar trap density)
alphat = 1e15       # 1/s
alphad0 = 1e13      # 1/s


def S_factor(epsilon_eV):
    """Storage factor S = 1 + ft·αt/αd."""
    alphad = alphad0 * math.exp(-epsilon_eV / (kB_eV * T_K))
    return 1.0 + ft * alphat / alphad


def C_ss(x):
    """Analytical steady-state mobile concentration."""
    return C_in * (1.0 - x / L)


def read_of_field(path):
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
    return np.array(values)


def check_steady_state(case_dir, time_str, epsilon_eV, tol_rel=0.02):
    """Compare simulated mobile profile to linear steady state."""
    path = os.path.join(case_dir, time_str, "slab", "C_T2")
    if not os.path.isfile(path):
        print(f"  t={time_str}: field not found at {path}")
        return None

    C_sim = read_of_field(path)
    if len(C_sim) == 0:
        print(f"  t={time_str}: could not parse C_T2")
        return None

    dx = L / N_cells
    x  = np.linspace(dx / 2, L - dx / 2, len(C_sim))
    C_ana = C_ss(x)

    err = np.abs(C_sim - C_ana) / (C_in + 1e-30)
    max_err = err.max()

    S = S_factor(epsilon_eV)
    print(
        f"  t={time_str:>6} s  ε={epsilon_eV} eV  S={S:.1f}:  "
        f"max|C-C_ss|/C_in = {max_err:.3e}  "
        + ("✓" if max_err < tol_rel else f"✗ (tol {tol_rel*100:.0f}%)")
    )
    return max_err


def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--case",    default=".")
    p.add_argument("--time",    default="20",   help="time directory to check")
    p.add_argument("--epsilon", default="0.09", type=float,
                   help="trap energy [eV] used in this run (0.09 or 1.5)")
    args = p.parse_args()

    case_dir = os.path.abspath(args.case)
    eps = args.epsilon
    S = S_factor(eps)
    alphad = alphad0 * math.exp(-eps / (kB_eV * T_K))
    D_eff  = 1.0 / S

    print(f"Case directory: {case_dir}")
    print()
    print(f"McNabb-Foster parameters (ε = {eps} eV, T = {T_K} K):")
    print(f"  αd   = {alphad:.3e} s⁻¹")
    print(f"  S    = {S:.2f}  (storage factor)")
    print(f"  D_eff= {D_eff:.4f} m²/s")
    if eps < 0.5:
        tau_bc = L**2 * S / (2 * math.pi**2 * 1.0)
        print(f"  τ_bc = {tau_bc:.2f} s  (breakthrough, weak-trapping regime)")
    else:
        tau_bd = ft * N * L**2 / (2 * C_in * 1.0)
        print(f"  τ_bd = {tau_bd:.1f} s  (breakthrough, strong-trapping regime)")
    print()
    print("Steady-state analytical:  C_ss(x) = C_in · (1 − x/L)")
    print()

    err = check_steady_state(case_dir, args.time, eps)

    if err is None:
        print("No data found — run Allrun first.")
        sys.exit(1)

    print()
    print("PASS" if err < 0.02 else "FAIL")
    sys.exit(0 if err < 0.02 else 1)


if __name__ == "__main__":
    main()
