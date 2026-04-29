#!/usr/bin/env python3
"""
Validation for case317-reverse-osmosis.

Solution-diffusion model: two regions (feed + membrane), each 0.5 mm thick.
Henry partition K_H = Ks_m / Ks_f = 1000/2000 = 0.5 at the interface.

Analytical steady-state:
  Feed   (x ∈ [0,   L/2]): C_f(x) = 1 - x/1.5e-3        mol/m³
  Membrane (x ∈ [L/2, L]): C_m(x) = 1/3 - (x-L/2)/1.5e-3  mol/m³

  Interface values:
    C_f,if = 2/3  mol/m³
    C_m,if = 1/3  mol/m³   (= K_H * C_f,if, K_H = 0.5)

  Flux (same in both regions at SS):
    J = D * ΔC_f / L_f = 1e-9 * (1-2/3) / 5e-4 = 6.667e-7  mol/(m²·s)
      = D * ΔC_m / L_m = 1e-9 * (1/3-0) / 5e-4 = 6.667e-7  mol/(m²·s)
"""

import sys
import os
import re
import numpy as np

TOL = 0.02          # 2 % relative tolerance
D  = 1e-9           # m²/s — diffusivity (both regions)
Lf = 5e-4           # m   — feed half-thickness
Lm = 5e-4           # m   — membrane half-thickness
L  = Lf + Lm        # m   — total thickness

# Inlet concentration (feed_inlet fixedValue) and outlet (membrane_outlet)
C_in  = 1.0         # mol/m³
C_out = 0.0         # mol/m³

# Analytical quantities
C_f_if = 2.0/3.0                              # feed interface
C_m_if = 1.0/3.0                              # membrane interface
J_ss   = D * (C_f_if - C_m_if) / Lm          # mol/(m²·s)  ≈ 6.667e-7


def parse_openfoam_field(path):
    """Return (coordinates_array, values_array) from a 1-D OpenFOAM field."""
    with open(path) as fh:
        text = fh.read()

    # --- internalField ---
    m = re.search(r'internalField\s+nonuniform\s+List<scalar>\s*\n\d+\s*\n\(([^)]+)\)', text)
    if m:
        vals = np.array([float(v) for v in m.group(1).split()])
    else:
        m = re.search(r'internalField\s+uniform\s+([\d.eE+\-]+)', text)
        if not m:
            raise ValueError(f"Cannot parse internalField in {path}")
        # All cells uniform — return None so caller can skip
        return None, float(m.group(1))

    return vals   # 1-D array, cell-centre values


def read_region_field(case_dir, region, time, field):
    path = os.path.join(case_dir, str(time), region, field)
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing: {path}")
    result = parse_openfoam_field(path)
    if isinstance(result, tuple):
        arr, scalar = result
        if arr is None:
            return np.full(1, scalar), scalar
        return arr, None
    return result, None


def cell_centres(x0, x1, n):
    dx = (x1 - x0) / n
    return np.linspace(x0 + dx/2, x1 - dx/2, n)


def check_region(label, vals, x, C_analytical, name):
    err = np.abs(vals - C_analytical) / (np.abs(C_analytical) + 1e-12)
    max_err = err.max()
    status = "PASS" if max_err <= TOL else "FAIL"
    print(f"  {label}: max relative error = {max_err*100:.2f}%  [{status}]")
    if max_err > TOL:
        worst = np.argmax(err)
        print(f"    worst cell: x={x[worst]*1e3:.3f} mm, "
              f"C_sim={vals[worst]:.4f}, C_ref={C_analytical[worst]:.4f}")
    return status == "PASS"


def main(case_dir="."):
    time_dirs = sorted(
        [d for d in os.listdir(case_dir)
         if re.match(r'^\d+(\.\d+)?$', d) and
            os.path.isdir(os.path.join(case_dir, d, "feed"))],
        key=float
    )
    if not time_dirs:
        print("ERROR: No time directories with region data found.")
        sys.exit(1)

    t_last = time_dirs[-1]
    print(f"Validating case317-reverse-osmosis at t = {t_last} s")
    print()

    n_feed = n_mem = 100   # cells per region (from blockMeshDict)

    x_feed = cell_centres(0,  Lf,      n_feed)
    x_mem  = cell_centres(Lf, Lf+Lm,  n_mem)

    # Analytical profiles
    C_f_ref = C_in - x_feed / (1.5e-3)          # 1 - x/1.5e-3
    C_m_ref = C_m_if - (x_mem - Lf) / (1.5e-3)  # 1/3 - (x-Lf)/1.5e-3

    # --- Feed region ---
    try:
        vals_f, scalar = read_region_field(case_dir, "feed", t_last, "C_H2O")
        if scalar is not None:
            vals_f = np.full(n_feed, scalar)
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        sys.exit(1)

    # --- Membrane region ---
    try:
        vals_m, scalar = read_region_field(case_dir, "membrane", t_last, "C_H2O")
        if scalar is not None:
            vals_m = np.full(n_mem, scalar)
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        sys.exit(1)

    ok_f = check_region("Feed     profile", vals_f, x_feed, C_f_ref,    "C_H2O")
    ok_m = check_region("Membrane profile", vals_m, x_mem,  C_m_ref,    "C_H2O")

    # Interface checks (last feed cell / first membrane cell)
    C_f_sim_if = vals_f[-1]
    C_m_sim_if = vals_m[0]
    print()
    print(f"  Interface (feed side)     C_f,if  : sim={C_f_sim_if:.4f}  ref={C_f_if:.4f}")
    print(f"  Interface (membrane side) C_m,if  : sim={C_m_sim_if:.4f}  ref={C_m_if:.4f}")

    # Flux from membrane gradient
    J_sim = D * (vals_m[0] - vals_m[-1]) / Lm
    print()
    print(f"  Flux J [mol/(m²·s)]: sim={J_sim:.3e}  ref={J_ss:.3e}")
    ok_J = abs(J_sim - J_ss) / J_ss < TOL

    print()
    if ok_f and ok_m and ok_J:
        print("Overall: PASS")
    else:
        print("Overall: FAIL")
        sys.exit(1)


if __name__ == "__main__":
    case_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    main(case_dir)
