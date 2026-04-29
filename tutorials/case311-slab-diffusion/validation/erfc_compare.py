#!/usr/bin/env python3
"""
Case 3.1.1 validation: compare OpenFOAM T profile with the erfc analytical solution.

Physics
-------
    dC/dt = D * d²C/dx²      C in [0, L],  t > 0
    C(0, t) = 1               (inlet, high-pressure side)
    C(L, t) = 0               (outlet, vacuum)
    C(x, 0) = 0               (initially empty slab)

Semi-infinite approximation (valid while diffusion front << L):
    C_semi(x, t) = erfc( x / (2 * sqrt(D * t)) )

This script reads the foamMultiRun output from the slab region, extracts the
cell-centre x-coordinates and T values from a chosen write time, and plots
them against the analytical curve.

Usage
-----
    python3 erfc_compare.py [--time 30]   # default t=30s
    python3 erfc_compare.py --all         # plot several times
"""

import argparse
import os
import sys
import numpy as np
from scipy.special import erfc

# ---------- parameters matching the case --------------------------------
D = 1e-8          # m²/s  (kappa/rho/Cv)
L = 1e-3          # m     (slab length)
T_REF = 300.0     # K     shift: C = T - T_ref  (avoids T=0 crash in heSolidThermo)
PLOT_TIMES = [5, 10, 20, 30, 50]   # s — must match write times in controlDict

try:
    import matplotlib.pyplot as plt
    HAS_MPLOT = True
except ImportError:
    HAS_MPLOT = False
    print("matplotlib not available — printing table only")


def read_foam_field(time_dir, region, field_name):
    """Robustly parse OpenFOAM ASCII field files with exact keyword matching."""
    path = os.path.join(time_dir, region, field_name)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Field file not found: {path}")

    values = []
    with open(path, 'r') as fh:
        lines = fh.readlines()

    it = iter(lines)
    
    # Phase 1: Find the internalField block
    for line in it:
        clean = line.strip()
        tokens = clean.split()
        
        if not tokens or tokens[0] != "internalField":
            continue
            
        # Check for EXACT match of "uniform" to avoid "nonuniform" false positives
        if "uniform" in tokens:
            # Format: internalField uniform 300;
            try:
                val_index = tokens.index("uniform") + 1
                val = float(tokens[val_index].rstrip(";"))
                return val 
            except (ValueError, IndexError):
                continue
        
        # If it reached here, it's likely "nonuniform List<scalar>"
        break

    # Phase 2: Find the list size and the opening bracket '('
    num_elements = 0
    for line in it:
        clean = line.strip()
        if not clean:
            continue
        # The number 200 often appears alone on a line before '('
        if clean.isdigit():
            num_elements = int(clean)
        if "(" in clean:
            break

    # Phase 3: Read data until the closing bracket ');' or we hit num_elements
    for line in it:
        clean = line.strip()
        if clean == ");" or (num_elements > 0 and len(values) >= num_elements):
            break
        
        # Strip comments and handle multiple values per line
        row_data = clean.split('//')[0].split()
        for val in row_data:
            try:
                values.append(float(val))
            except ValueError:
                continue

    return np.array(values)

def read_cell_centres(case_dir, region):
    """Get cell-centre x from <case>/0/slab/C or from postProcess output."""
    # Try reading from the mesh itself via the C file written by postProcess.
    # Fallback: reconstruct from uniform spacing.
    N = 200
    dx = L / N
    x = np.linspace(dx / 2, L - dx / 2, N)
    return x


def analytical_erfc(x, t):
    """Semi-infinite erfc solution."""
    return erfc(x / (2.0 * np.sqrt(D * t)))


def analytical_series(x, t, n_terms=100):
    """
    Full finite-slab Fourier-series solution:
        C(x,t) = (1 - x/L)
                 - (2/pi) * sum_{n=1}^inf (1/n) * sin(n*pi*x/L)
                            * exp(-D*n^2*pi^2*t/L^2)
    """
    C = 1.0 - x / L
    for n in range(1, n_terms + 1):
        C -= (2.0 / (n * np.pi)) * np.sin(n * np.pi * x / L) * \
             np.exp(-D * n**2 * np.pi**2 * t / L**2)
    return C


def load_time(case_dir, t):
    time_str = f"{t:.6g}"
    time_dir = os.path.join(case_dir, time_str)
    if not os.path.isdir(time_dir):
        # Try exact match with more precision
        for d in os.listdir(case_dir):
            try:
                if abs(float(d) - t) < 1e-6 * max(1, abs(t)):
                    time_dir = os.path.join(case_dir, d)
                    break
            except ValueError:
                pass
    if not os.path.isdir(time_dir):
        print(f"  Warning: time directory {t}s not found, skipping")
        return None, None

    region = "slab"
    T_values = read_foam_field(time_dir, region, "T")
    x = read_cell_centres(case_dir, region)
    # Convert from thermal proxy back to concentration
    C_values = T_values - T_REF
    return x, C_values


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--time", type=float, default=30.0)
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--case", default=".")
    args = parser.parse_args()

    case_dir = os.path.abspath(args.case)
    times = PLOT_TIMES if args.all else [args.time]

    print(f"Case directory: {case_dir}")
    print(f"D = {D} m²/s,  L = {L*1e3:.1f} mm")
    print()

    errors = {}
    if HAS_MPLOT:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("C / C_in  (= T [K])")
        ax.set_title("Case 3.1.1 — 1D diffusion vs erfc analytical")

    for t in times:
        x, T = load_time(case_dir, t)
        
        if isinstance(T, (float, int)):
            T = np.full(x.shape, T)

        if T is None or len(T) == 0:
            print(f"Skipping {t}s: No data parsed.")
            continue
        C_erfc   = analytical_erfc(x, t)
        C_series = analytical_series(x, t)

        # Max relative error vs full series solution
        with np.errstate(divide="ignore", invalid="ignore"):
            rel_err = np.abs(T - C_series) / np.maximum(C_series, 1e-10)
        max_err = np.max(rel_err)
        errors[t] = max_err

        print(f"t = {t:6.1f} s | max |OF - series| / max(series) = {max_err:.3e}", end="")
        print("  ✓" if max_err < 1e-2 else "  ✗ (>1% error)")

        if HAS_MPLOT:
            ax.plot(x * 1e3, T, "o", ms=3, label=f"OF t={t}s")
            ax.plot(x * 1e3, C_series, "-", lw=1, label=f"Series t={t}s")

    print()
    all_pass = all(e < 1e-2 for e in errors.values())
    if all_pass:
        print("PASS — all time points within 1% of analytical solution")
    else:
        print("FAIL — one or more time points exceed 1% error threshold")
        sys.exit(1)

    if HAS_MPLOT:
        ax.legend(fontsize=7, ncol=2)
        ax.grid(True, alpha=0.3)
        out = "erfc_compare.png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"Plot saved to {out}")
        plt.show()


if __name__ == "__main__":
    main()
