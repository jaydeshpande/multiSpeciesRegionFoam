"""
Surface recombination/dissociation kinetics vs Sieverts equilibrium BC.

Single metal slab: L = 0.5 mm, D = 5e-11 m²/s.
Gas face partial pressure: pGas = 1 Pa.
Equilibrium concentration: C_eq = sqrt(Kd*pGas/Kr) = 1.0 mol/m³.

Sub-case A — sieverts_bc:
  Left BC: fixedValue C = C_eq = 1.0 mol/m³  (instantaneous kinetics, Da→∞)
  Steady-state flux: J_ref = D * C_eq / L = 5e-11 * 1.0 / 5e-4 = 1e-7 mol/(m²·s)

Sub-case B — surface_recombination (Da = 1):
  Left BC: D * dC/dn = Kd*pGas − Kr*C²
           Kd = Kr = 1e-7;  pGas = 1 Pa → C_eq = 1.0,  Da = Kr*C_eq*L/D = 1
  Solving steady state:  Kr*Cs² + (D/L)*Cs − Kd*pGas = 0
    1e-7*Cs² + 1e-7*Cs − 1e-7 = 0 → Cs² + Cs − 1 = 0
    Cs = (-1 + sqrt(5))/2 ≈ 0.6180 mol/m³   (golden-ratio conjugate)
  Steady-state flux: J_surf = D * Cs / L ≈ 6.18e-8 mol/(m²·s)
  Reduction: J_surf/J_ref = Cs/C_eq ≈ 0.618  (38.2% lower due to kinetics)

Physical interpretation:
  When Da = 1, the surface kinetics and bulk diffusion have equal resistance.
  Increasing Kd and Kr (while keeping Kd/Kr = C_eq²/pGas constant to preserve
  equilibrium) drives Da → ∞ and recovers the Sieverts fixedValue result.
  The surface recombination BC is the correct generalisation: it captures both
  kinetically limited and diffusion-limited regimes in a single Robin condition.
"""

import os
import re
import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Paths anchored to the script file so the script can be run from any CWD.
_SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
_CASE_DIR    = os.path.normpath(os.path.join(_SCRIPT_DIR, ".."))
_SV_DEFAULT  = os.path.join(_CASE_DIR, "sieverts_bc")
_SURF_DEFAULT = os.path.join(_CASE_DIR, "surface_recombination")

# numpy ≥ 2.0 compatibility
try:
    _trapz = np.trapezoid
except AttributeError:
    _trapz = np.trapz

# ── Physical parameters ────────────────────────────────────────────────────
D     = 5e-11    # m²/s
L     = 0.5e-3   # m
pGas  = 1.0      # Pa
Kd    = 1e-7     # mol/(m²·s·Pa)
Kr    = 1e-7     # m⁴/(mol·s)
N     = 50       # cells
dx    = L / N    # m


# ── Analytical steady-state solutions ─────────────────────────────────────
def analytical_sieverts():
    C_eq = np.sqrt(Kd * pGas / Kr)
    J    = D * C_eq / L
    return {"C_s": C_eq, "J": J, "label": "Sieverts BC (Da → ∞)"}


def analytical_surface_recombination():
    """Solve Kr*Cs² + (D/L)*Cs − Kd*pGas = 0 for Cs."""
    a = Kr
    b = D / L
    c = -Kd * pGas
    Cs = (-b + np.sqrt(b**2 - 4*a*c)) / (2*a)
    J  = D * Cs / L
    Da = Kr * Cs * L / D
    return {"C_s": Cs, "J": J, "Da": Da, "label": f"surfaceRecombination (Da≈{Da:.1f})"}


# ── Analytical concentration profile ──────────────────────────────────────
def profile(C_s, n_pts=100):
    """Linear steady-state profile from C_s (x=0) to 0 (x=L)."""
    x = np.linspace(0, L, n_pts)
    C = C_s * (1 - x/L)
    return x, C


# ── OpenFOAM reader ────────────────────────────────────────────────────────
def read_of_field(path):
    with open(path) as fh:
        text = fh.read()
    m = re.search(
        r'internalField\s+nonuniform\s+List<scalar>\s*\n\d+\s*\n\(([^)]+)\)', text)
    if m:
        return np.fromstring(m.group(1), sep='\n')
    m = re.search(r'internalField\s+uniform\s+([\d.eE+\-]+)', text)
    if m:
        return np.full(1, float(m.group(1)))
    raise ValueError(f"Cannot parse internalField in {path}")


def available_times(case_dir, region, field):
    times = []
    for entry in os.listdir(case_dir):
        try:
            float(entry)
        except ValueError:
            continue
        if os.path.isfile(os.path.join(case_dir, entry, region, field)):
            times.append(entry)
    return sorted(times, key=float)


# ── main ──────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sieverts", default=_SV_DEFAULT,
                   help="path to sieverts_bc case")
    p.add_argument("--surf",     default=_SURF_DEFAULT,
                   help="path to surface_recombination case")
    p.add_argument("--plot",     action="store_true")
    args = p.parse_args()

    sv_dir   = os.path.abspath(args.sieverts)
    surf_dir = os.path.abspath(args.surf)

    ana_sv   = analytical_sieverts()
    ana_surf = analytical_surface_recombination()

    print("=== Analytical steady-state ===")
    for ana in (ana_sv, ana_surf):
        print(f"\n  {ana['label']}")
        print(f"    Surface concentration Cs = {ana['C_s']:.4f}  mol/m³")
        print(f"    Permeation flux J        = {ana['J']:.3e}  mol/(m²·s)")
    ratio = ana_surf["J"] / ana_sv["J"]
    print(f"\n  Flux ratio J_surf/J_ref = {ratio:.4f}  "
          f"({(1-ratio)*100:.1f}% reduction from surface kinetics)")
    print(f"  Cs = (sqrt(5)-1)/2 = {(-1+np.sqrt(5))/2:.6f}  "
          f"(golden-ratio conjugate — exact when Da=1, Kd=Kr, pGas=1)")

    print("\n=== Simulation results ===")
    passed_all = True
    tol = 0.05

    for label, case_dir, ana in [
        ("sieverts_bc",           sv_dir,   ana_sv),
        ("surface_recombination", surf_dir, ana_surf),
    ]:
        times = available_times(case_dir, "metal", "C_H2")
        if not times:
            print(f"  {label}: ERROR — no time directories. Run the simulation first.")
            passed_all = False
            continue
        t_last = sorted(times, key=float)[-1]
        C = read_of_field(os.path.join(case_dir, t_last, "metal", "C_H2"))
        J_sim = D * C[-1] / (dx / 2.0)
        rel_err = abs(J_sim - ana["J"]) / ana["J"]
        passed = rel_err < tol
        if not passed:
            passed_all = False
        print(f"  {label}  t={float(t_last):.0f}s  J_sim={J_sim:.3e}  "
              f"J_ana={ana['J']:.3e}  err={rel_err:.2%}  {'PASS' if passed else 'FAIL'}")

    print()
    print("PASS" if passed_all else "FAIL")

    if not args.plot:
        sys.exit(0 if passed_all else 1)

    # ── Plots ──────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    x_cells = (np.arange(N) + 0.5) * dx * 1e3    # mm

    # 1. Concentration profiles
    ax = axes[0]
    colors = {"sieverts_bc": "steelblue", "surface_recombination": "darkorange"}

    for label, case_dir, ana in [
        ("sieverts_bc",           sv_dir,   ana_sv),
        ("surface_recombination", surf_dir, ana_surf),
    ]:
        times = available_times(case_dir, "metal", "C_H2")
        if not times:
            continue
        t_last = sorted(times, key=float)[-1]
        C_sim = read_of_field(os.path.join(case_dir, t_last, "metal", "C_H2"))
        c = colors[label]
        ax.plot(x_cells, C_sim, color=c, lw=2, label=label.replace("_", " "))
        # Analytical profile
        x_ana, C_ana = profile(ana["C_s"])
        ax.plot(x_ana*1e3, C_ana, color=c, ls="--", lw=1, alpha=0.7)

    # Mark analytical surface concentrations
    for label, ana, c in [
        ("sieverts_bc",           ana_sv,   "steelblue"),
        ("surface_recombination", ana_surf, "darkorange"),
    ]:
        ax.axhline(ana["C_s"], color=c, ls=":", lw=1, alpha=0.5,
                   label=f"Cs={ana['C_s']:.3f} ({label[:4]})")

    ax.set_xlabel("Position  [mm]")
    ax.set_ylabel("C_H2  [mol/m³]")
    ax.set_title("Concentration profile at steady state\n"
                 "Surface recombination (Da=1) → Cs=0.618 < C_eq=1.0")
    ax.legend(fontsize=8)
    ax.grid(True, ls=":")

    # 2. Flux vs time
    ax = axes[1]
    for label, case_dir, ana in [
        ("sieverts_bc",           sv_dir,   ana_sv),
        ("surface_recombination", surf_dir, ana_surf),
    ]:
        times = available_times(case_dir, "metal", "C_H2")
        if not times:
            continue
        t_arr = np.array(sorted(times, key=float), dtype=float)
        J_arr = []
        for t_str in sorted(times, key=float):
            C = read_of_field(os.path.join(case_dir, t_str, "metal", "C_H2"))
            J_arr.append(D * C[-1] / (dx / 2.0))
        c = colors[label]
        ax.plot(t_arr, J_arr, color=c, lw=2, marker="o", ms=4,
                label=label.replace("_", " "))
        ax.axhline(ana["J"], color=c, ls="--", lw=1, alpha=0.7,
                   label=f"J_ana={ana['J']:.1e}")

    ax.set_xlabel("Time  [s]")
    ax.set_ylabel("Permeation flux  [mol/(m²·s)]")
    ax.set_title("Flux vs time\n"
                 "Kinetically limited (Da=1): flux is 38.2% lower than Sieverts BC")
    ax.legend(fontsize=8)
    ax.grid(True, ls=":")

    fig.tight_layout()
    out = os.path.normpath(
        os.path.join(os.path.dirname(__file__), "..", "recombination_compare.png"))
    fig.savefig(out, dpi=150)
    print(f"\nPlot saved to {out}")
    sys.exit(0 if passed_all else 1)


if __name__ == "__main__":
    main()
