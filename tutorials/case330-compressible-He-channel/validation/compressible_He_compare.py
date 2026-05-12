"""
Validation script for case330-compressible-He-channel.

Compares the compressibleSpeciesFluid (He gas) result with:
  1. The analytical bulk concentration at outlet.
  2. The speciesFluid (FLiBe) result from case329-flibe-channel,
     to verify that identical D, Ks, source, and velocity yield the
     same permeation flux regardless of the fluid thermophysics model.

Physical parameters (both cases share these values)
----------------------------------------------------
  D_fluid  = D_solid = 5e-8 m²/s   (no Arrhenius, Ea=0)
  Ks_fluid = Ks_solid = 1000 mol/(m³·Pa^0.5)  → transparent interface
  S        = 1e-5 mol/(m³·s)        uniform source in fluid
  U        = 0.01 m/s               plug flow
  L        = 0.050 m                channel length
  W_solid  = 0.002 m                SS316 wall thickness
  Nx_solid = 10                     cells across the solid

Analytical estimates (plug-flow, 1-D)
--------------------------------------
  C_bulk_out ≈ S·L/U = 1e-5·0.05/0.01 = 5e-5 mol/m³
  J_perm     ≈ D_solid·(C_interface)/W_solid  (at steady state)
    where C_interface is the solid-side interfacial concentration.
    For the transparent interface (Ks_ratio=1) and high fluid-side supply,
    C_interface ≈ C_bulk_out → J_perm ≈ 5e-8·5e-5/2e-3 = 1.25e-9 mol/(m²·s).
    This is only an order-of-magnitude estimate; the actual value comes from
    the simulation.

PASS/FAIL criteria
------------------
  - |C_bulk_He - C_bulk_analytical| / C_bulk_analytical < 0.05   (5% tolerance)
  - If case329 results exist:
      |J_He - J_FLiBe| / J_FLiBe < 0.05   (5% tolerance)

Exit 0 on PASS, 1 on FAIL.

Usage
-----
  python compressible_He_compare.py [--case330 <path>] [--case329 <path>]
"""

import os
import sys
import re
import argparse
import numpy as np

# ── Physical parameters ──────────────────────────────────────────────────────

S       = 1.0e-5    # mol/(m³·s)  volumetric H2 source
L       = 0.050     # m           channel length
U       = 0.01      # m/s         inlet velocity
D_solid = 5.0e-8    # m²/s        solid diffusivity
W_solid = 2.0e-3    # m           solid wall thickness
Nx_solid = 10       # cells in solid (x-direction)

# Analytical steady-state outlet bulk concentration (plug flow, no permeation)
C_bulk_analytical = S * L / U   # = 5e-5 mol/m³

PASS_TOL = 0.05     # 5 % relative tolerance


# ── OpenFOAM field reader ────────────────────────────────────────────────────

def read_of_scalar_field(path):
    """
    Read an OpenFOAM ASCII volScalarField and return the internalField
    as a numpy array.  Handles both 'uniform <val>' and
    'nonuniform List<scalar> N ( ... )' formats.
    """
    with open(path) as fh:
        text = fh.read()

    # Non-uniform field
    m = re.search(
        r'internalField\s+nonuniform\s+List<scalar>\s*\n\d+\s*\n\(([^)]+)\)',
        text, re.DOTALL)
    if m:
        return np.fromstring(m.group(1), sep='\n')

    # Uniform field
    m = re.search(r'internalField\s+uniform\s+([\d.eE+\-]+)', text)
    if m:
        return np.full(1, float(m.group(1)))

    raise ValueError(f"Cannot parse internalField in {path}")


def available_times(case_dir, region, field):
    """Return sorted list of time-step strings that contain <region>/<field>."""
    times = []
    for entry in os.listdir(case_dir):
        try:
            float(entry)
        except ValueError:
            continue
        fpath = os.path.join(case_dir, entry, region, field)
        if os.path.isfile(fpath):
            times.append(entry)
    return sorted(times, key=float)


def last_time(case_dir, region, field):
    """Return the largest available time step string."""
    times = available_times(case_dir, region, field)
    if not times:
        return None
    return times[-1]


# ── Flux computation ─────────────────────────────────────────────────────────

def permeation_flux_from_solid(C_solid):
    """
    Estimate the outer-wall permeation flux from the solid C_H2 field.

    Uses the half-cell approximation:
      J = D_solid * C_last_cell / (dx/2)
    where dx = W_solid / Nx_solid and C_last_cell is the concentration
    in the outermost solid cell (adjacent to the fixedValue=0 boundary).

    The solid cells are ordered x=5mm → x=7mm (inner → outer).
    Taking the last cell in the 1-D x-array gives the outer-most layer.
    In 2-D (Ny rows × Nx columns), the field is stored cell by cell;
    we average over all y-rows to get a representative 1-D value, then
    use the outer column (every Nx_solid-th cell starting at index Nx_solid-1).
    """
    dx = W_solid / Nx_solid
    # Field has Ny*Nx_solid entries arranged as rows of y, columns of x.
    # Outer column: indices Nx_solid-1, 2*Nx_solid-1, ...
    if C_solid.size == 1:
        C_outer = float(C_solid[0])
    else:
        # Try to extract the outer column
        Ntotal = C_solid.size
        Ny = Ntotal // Nx_solid
        if Ny * Nx_solid == Ntotal:
            outer_col = C_solid[(Nx_solid - 1)::Nx_solid]
        else:
            # Fall back: use last Ny values
            Ny = max(1, Ntotal // Nx_solid)
            outer_col = C_solid[-(Ny):]
        C_outer = float(np.mean(outer_col))

    J = D_solid * C_outer / (dx / 2.0)
    return J


# ── main ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--case330",
        default=os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."),
        help="Path to case330-compressible-He-channel (default: parent of this script)")
    parser.add_argument(
        "--case329",
        default=None,
        help="Path to case329-flibe-channel for cross-solver comparison "
             "(optional; skip comparison if not provided or results absent)")
    args = parser.parse_args()

    case330 = os.path.abspath(args.case330)
    case329 = os.path.abspath(args.case329) if args.case329 else \
              os.path.join(case330, "..", "case329-flibe-channel")

    all_passed = True

    print("=" * 65)
    print("  case330-compressible-He-channel  validation")
    print("=" * 65)
    print()

    # ── 1. Read case330 results ──────────────────────────────────────────────

    t_fluid = last_time(case330, "fluid", "C_H2")
    t_solid = last_time(case330, "solid", "C_H2")

    if t_fluid is None or t_solid is None:
        print("ERROR: No simulation output found in case330.")
        print(f"  Expected time directories under: {case330}/fluid/ and {case330}/solid/")
        print("  Run 'Allrun' first.")
        sys.exit(1)

    print(f"case330 last time step : fluid={t_fluid} s,  solid={t_solid} s")
    print()

    # Fluid bulk outlet concentration (average of top-outlet cells)
    C_fluid = read_of_scalar_field(
        os.path.join(case330, t_fluid, "fluid", "C_H2"))
    # Use the y-maximum cells (last Nx_fluid rows) as outlet approximation;
    # for plug flow the bulk average of all cells is a reasonable proxy.
    C_bulk_He = float(np.mean(C_fluid))

    # Solid permeation flux
    C_solid_He = read_of_scalar_field(
        os.path.join(case330, t_solid, "solid", "C_H2"))
    J_He = permeation_flux_from_solid(C_solid_He)

    # ── 2. Bulk concentration check ──────────────────────────────────────────

    print("--- Check 1: bulk outlet concentration (analytical) ---")
    print(f"  Analytical C_bulk = S·L/U = {C_bulk_analytical:.3e} mol/m³")
    print(f"  Simulated  C_bulk (mean) = {C_bulk_He:.3e} mol/m³")
    rel_err_bulk = abs(C_bulk_He - C_bulk_analytical) / C_bulk_analytical
    print(f"  Relative error = {rel_err_bulk:.2%}")
    if rel_err_bulk < PASS_TOL:
        print(f"  PASS  (< {PASS_TOL:.0%})")
    else:
        print(f"  FAIL  (>= {PASS_TOL:.0%})")
        all_passed = False
    print()

    # ── 3. Cross-solver permeation flux comparison with case329 ─────────────

    print("--- Check 2: permeation flux comparison He vs FLiBe ---")

    case329_ok = False
    J_FLiBe = None

    if os.path.isdir(case329):
        t_solid_329 = last_time(case329, "solid", "C_H2")
        if t_solid_329 is not None:
            try:
                C_solid_FLiBe = read_of_scalar_field(
                    os.path.join(case329, t_solid_329, "solid", "C_H2"))
                J_FLiBe = permeation_flux_from_solid(C_solid_FLiBe)
                case329_ok = True
                print(f"  case329 last solid time: {t_solid_329} s")
            except Exception as exc:
                print(f"  WARNING: Could not read case329 solid C_H2: {exc}")
        else:
            print(f"  case329 found at {case329} but no solid/C_H2 output.")
    else:
        print(f"  case329 not found at {case329}")

    print(f"  He   permeation flux J_He    = {J_He:.4e} mol/(m²·s)")

    if case329_ok and J_FLiBe is not None:
        print(f"  FLiBe permeation flux J_FLiBe = {J_FLiBe:.4e} mol/(m²·s)")
        if J_FLiBe > 0:
            rel_err_flux = abs(J_He - J_FLiBe) / J_FLiBe
            print(f"  Relative difference = {rel_err_flux:.2%}")
            if rel_err_flux < PASS_TOL:
                print(f"  PASS  (< {PASS_TOL:.0%})")
            else:
                print(f"  FAIL  (>= {PASS_TOL:.0%})")
                all_passed = False
        else:
            print("  WARNING: J_FLiBe = 0, cannot compute relative error.  SKIP.")
    else:
        print("  SKIP (case329 results not available; run case329 first for")
        print("        cross-solver comparison)")
    print()

    # ── Summary ─────────────────────────────────────────────────────────────

    print("=" * 65)
    print(f"  Overall result: {'PASS' if all_passed else 'FAIL'}")
    print("=" * 65)

    sys.exit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
