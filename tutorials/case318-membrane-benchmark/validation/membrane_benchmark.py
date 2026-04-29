"""
Benchmark validation: 1-D diffusion in a 1.2 mm metal membrane.

Hattab et al. 2024, Fusion Eng. Des. 202, 114362, Section 3.3.1.

Compares the multiSpeciesRegionFoam (speciesSolid) result against the
exact Fourier sine series analytical solution:

  C(x,t) = C1*(1 - x/L)
            - (2*C1/pi) * sum_{n=1}^{N} (1/n) * sin(n*pi*x/L)
                          * exp(-(n*pi/L)**2 * D * t)

Parameters (from paper):
  L = 1.2e-3 m,  D = 1.119018e-8 m²/s
  C1 = 0.83 mol/m³ (inlet),  C0 = 0 mol/m³ (outlet)
  tau_SS = L²/D ≈ 128.7 s;  simulation runs to t = 30 s (0.23 tau_SS)
"""

import os
import sys
import re
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── physical parameters ────────────────────────────────────────────────────────
L  = 1.2e-3          # membrane thickness [m]
D  = 1.119018e-8     # diffusivity [m²/s]
C1 = 0.83            # inlet concentration [mol/m³]
N_TERMS = 200        # Fourier series truncation


def analytical(x, t):
    """Fourier sine series solution for C(0)=C1, C(L)=0, C(x,0)=0."""
    C = C1 * (1.0 - x / L)
    for n in range(1, N_TERMS + 1):
        lam = (n * np.pi / L) ** 2 * D
        C -= (2.0 * C1 / (n * np.pi)) * np.sin(n * np.pi * x / L) * np.exp(-lam * t)
    return C


# ── OpenFOAM field reader ──────────────────────────────────────────────────────

def read_of_field(path):
    """
    Parse an OpenFOAM ASCII volScalarField and return the internalField values
    as a numpy array.  Handles both 'uniform' and 'nonuniform List<scalar>'.
    """
    with open(path) as fh:
        text = fh.read()

    # nonuniform list
    m = re.search(
        r'internalField\s+nonuniform\s+List<scalar>\s*\n\d+\s*\n\(([^)]+)\)',
        text
    )
    if m:
        return np.fromstring(m.group(1), sep='\n')

    # uniform scalar
    m = re.search(r'internalField\s+uniform\s+([\d.eE+\-]+)', text)
    if m:
        return None, float(m.group(1))   # (sentinel, value)

    raise ValueError(f"Cannot parse internalField in {path}")


def available_times(case_dir, region="membrane", field="C_H2"):
    """Return sorted list of time directories that contain the requested field."""
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


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--case", default=".", help="path to case directory")
    parser.add_argument("--times", nargs="+", default=None,
                        help="time directories to compare (default: auto-detect)")
    parser.add_argument("--tol", type=float, default=0.02,
                        help="L-inf relative tolerance for PASS/FAIL (default 0.02)")
    parser.add_argument("--plot", action="store_true", help="save comparison plot")
    args = parser.parse_args()

    case_dir = os.path.abspath(args.case)
    times = args.times or available_times(case_dir)

    if not times:
        print("ERROR: no time directories found with membrane/C_H2 data.", file=sys.stderr)
        sys.exit(2)

    print(f"Case: {case_dir}")
    print(f"Times: {times}\n")

    # cell-centre x positions (uniform mesh, 90 cells)
    n_cells = 90
    x_sim = (np.arange(n_cells) + 0.5) * L / n_cells

    results = {}   # t_str → max relative error

    fig, axes = None, None
    if args.plot:
        ncols = min(3, len(times))
        nrows = (len(times) + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows),
                                 squeeze=False)

    for idx, t_str in enumerate(times):
        t = float(t_str)
        fpath = os.path.join(case_dir, t_str, "membrane", "C_H2")

        result = read_of_field(fpath)
        if isinstance(result, tuple):
            # uniform field
            C_sim = np.full(n_cells, result[1])
        else:
            C_sim = result

        C_ana = analytical(x_sim, t)

        denom = max(C1, np.max(np.abs(C_ana)))
        rel_err = np.max(np.abs(C_sim - C_ana)) / denom
        results[t_str] = rel_err

        print(f"  t = {t:6.1f} s  ({t/128.7:.3f} τ_SS)  "
              f"max |err|/C1 = {rel_err:.4e}  "
              f"{'PASS' if rel_err < args.tol else 'FAIL'}")

        if args.plot:
            ax = axes[idx // ncols][idx % ncols]
            x_ana = np.linspace(0, L, 500)
            ax.plot(x_ana * 1e3, analytical(x_ana, t), "k-", label="analytical")
            ax.plot(x_sim * 1e3, C_sim, "r.", ms=3, label="OF")
            ax.set_xlabel("x [mm]")
            ax.set_ylabel("C [mol/m³]")
            ax.set_title(f"t = {t:.1f} s ({t/128.7:.2f} τ_SS)")
            ax.legend(fontsize=7)
            ax.set_ylim(-0.02, 0.90)

    if args.plot:
        for idx in range(len(times), nrows * ncols):
            axes[idx // ncols][idx % ncols].set_visible(False)
        fig.tight_layout()
        outname = os.path.join(case_dir, "membrane_benchmark.png")
        fig.savefig(outname, dpi=150)
        print(f"\nPlot saved to {outname}")

    # PASS/FAIL evaluated at the last available time
    last_t = sorted(results, key=float)[-1]
    passed = results[last_t] < args.tol
    print(f"\n(Pass/fail evaluated at t = {last_t} s — last available write time)")
    print("PASS" if passed else "FAIL")
    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
