"""
Henry's law vs Sieverts' law partition — steady-state comparison.

Two-layer membrane (polymer | metal), each 0.5 mm thick, equal diffusivities
D = 1e-9 m²/s.  Both cases have the same boundary conditions:
  C_polymer(x=0) = 1.0 mol/m³  (left inlet)
  C_metal(x=L)   = 0.0 mol/m³  (right purge)

Case A — Sieverts-Sieverts (partition linear, Ks=1 both sides):
  Interface: C_metal = (Ks_metal/Ks_polymer) * C_polymer = C_polymer
  Concentration is CONTINUOUS at the interface.
  J_SS = D * C0 / (2L) = 1e-9 * 1.0 / 1e-3 = 1e-6 mol/(m²·s)

Case B — Henry-Sieverts (partition sqrt):
  Polymer: C_polymer = KH * p    (Henry's law)
  Metal:   C_metal   = Ks * sqrt(p) = sqrt(p)  (Sieverts, Ks=1)
  Interface: C_metal = sqrt(C_polymer)  → C_metal > C_polymer for C < 1
  Solving:   C_metal,int² + C_metal,int - 1 = 0
             C_metal,int  = (sqrt(5)-1)/2 ≈ 0.6180  mol/m³
             C_polymer,int = C_metal,int² ≈ 0.3820  mol/m³
  J_HS = D * C_metal,int / L = 1e-9 * 0.6180 / 5e-4 = 1.236e-6 mol/(m²·s)
  → 23.6% HIGHER flux than Sieverts-Sieverts, despite equal D and equal K values.
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
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_CASE_DIR   = os.path.normpath(os.path.join(_SCRIPT_DIR, ".."))
_SS_DEFAULT = os.path.join(_CASE_DIR, "sieverts_sieverts")
_HS_DEFAULT = os.path.join(_CASE_DIR, "henry_sieverts")

# numpy ≥ 2.0 compatibility
try:
    _trapz = np.trapezoid
except AttributeError:
    _trapz = np.trapz

# ── Physical parameters ────────────────────────────────────────────────────
D       = 1e-9      # m²/s  (equal in both layers)
L       = 0.5e-3    # m     (each layer thickness)
C0      = 1.0       # mol/m³ (left BC)
N_poly  = 20        # cells in polymer
N_metal = 50        # cells in metal


# ── Analytical steady-state ────────────────────────────────────────────────
def analytical_sieverts_sieverts():
    """Linear partition, Ks=1 both sides → concentration continuous."""
    J = D * C0 / (2 * L)
    C_int = C0 / 2.0              # midpoint: continuous concentration
    return {"J": J, "C_poly_int": C_int, "C_metal_int": C_int,
            "label": "Sieverts-Sieverts (linear)"}


def analytical_henry_sieverts():
    """sqrt partition: C_metal = sqrt(C_polymer) at interface (KH=Ks=1).

    At steady state, equal fluxes give:
      D*(C0 - C_A,int)/L = D*C_B,int/L   with C_B,int = sqrt(C_A,int)
    Let u = C_B,int = sqrt(C_A,int):
      C0 - u² = u  →  u² + u - 1 = 0
      u = (-1 + sqrt(5))/2 = golden-ratio conjugate ≈ 0.6180
    """
    u = (-1.0 + np.sqrt(5.0)) / 2.0    # = (sqrt(5)-1)/2 ≈ 0.6180
    C_metal_int  = u
    C_poly_int   = u * u
    J = D * C_metal_int / L
    return {"J": J, "C_poly_int": C_poly_int, "C_metal_int": C_metal_int,
            "label": "Henry-Sieverts (sqrt)"}


# ── OpenFOAM field reader ──────────────────────────────────────────────────
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


def last_flux(case_dir, region, field, patch_last_cell, dx):
    """Estimate flux at purge face from last cell concentration."""
    times = available_times(case_dir, region, field)
    if not times:
        return None, None
    t = times[-1]
    C = read_of_field(os.path.join(case_dir, t, region, field))
    J = D * C[-1] / (dx / 2.0)
    return float(t), J


# ── main ──────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ss",   default=_SS_DEFAULT,
                   help="path to sieverts_sieverts case")
    p.add_argument("--hs",   default=_HS_DEFAULT,
                   help="path to henry_sieverts case")
    p.add_argument("--plot", action="store_true", help="save comparison plots")
    args = p.parse_args()

    ss_dir = os.path.abspath(args.ss)
    hs_dir = os.path.abspath(args.hs)

    ana_ss = analytical_sieverts_sieverts()
    ana_hs = analytical_henry_sieverts()

    print("=== Analytical steady-state ===")
    for ana in (ana_ss, ana_hs):
        print(f"\n  {ana['label']}")
        print(f"    C_polymer,int = {ana['C_poly_int']:.4f}  mol/m³")
        print(f"    C_metal,int   = {ana['C_metal_int']:.4f}  mol/m³")
        print(f"    Flux J        = {ana['J']:.3e}  mol/(m²·s)")
    print(f"\n  Flux ratio J_HS/J_SS = {ana_hs['J']/ana_ss['J']:.4f}  "
          f"(Henry-Sieverts gives {(ana_hs['J']/ana_ss['J']-1)*100:+.1f}%)")

    dx_metal = L / N_metal

    print("\n=== Simulation results (permeation flux at metal_right) ===")
    passed_all = True
    results = {}
    for label, case_dir, ana in [
        ("sieverts_sieverts", ss_dir, ana_ss),
        ("henry_sieverts",    hs_dir, ana_hs),
    ]:
        t_last, J_sim = last_flux(case_dir, "metal", "C_H2", N_metal, dx_metal)
        if J_sim is None:
            print(f"  {label}: ERROR — no time directories found. Run the case first.")
            passed_all = False
            results[label] = None
            continue
        rel_err = abs(J_sim - ana["J"]) / ana["J"]
        passed = rel_err < 0.05
        if not passed:
            passed_all = False
        results[label] = {"J_sim": J_sim, "t": t_last, "J_ana": ana["J"]}
        print(f"  {label}  t={t_last:.0f} s  J_sim={J_sim:.3e}  "
              f"J_ana={ana['J']:.3e}  err={rel_err:.2%}  "
              f"{'PASS' if passed else 'FAIL'}")

    print()
    print("PASS" if passed_all else "FAIL")

    if not args.plot:
        sys.exit(0 if passed_all else 1)

    # ── Plots ──────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # 1. Concentration profiles (last time)
    ax = axes[0]
    dx_poly  = L / N_poly
    x_poly  = (np.arange(N_poly)  + 0.5) * dx_poly  * 1e3          # mm
    x_metal = L*1e3 + (np.arange(N_metal) + 0.5) * dx_metal * 1e3  # mm

    colors = {"sieverts_sieverts": "steelblue", "henry_sieverts": "darkorange"}
    styles = {"sieverts_sieverts": "-",          "henry_sieverts": "--"}

    for label, case_dir, ana in [
        ("sieverts_sieverts", ss_dir, ana_ss),
        ("henry_sieverts",    hs_dir, ana_hs),
    ]:
        times = available_times(case_dir, "metal", "C_H2")
        if not times:
            continue
        t = times[-1]
        C_p = read_of_field(os.path.join(case_dir, t, "polymer", "C_H2"))
        C_m = read_of_field(os.path.join(case_dir, t, "metal",   "C_H2"))
        c = colors[label]; ls = styles[label]
        ax.plot(x_poly,  C_p, color=c, ls=ls, lw=2)
        ax.plot(x_metal, C_m, color=c, ls=ls, lw=2, label=label.replace("_", "-"))

    # Analytical reference lines at the interface
    ax.axvline(L*1e3, color="gray", ls=":", lw=1, label="interface")
    ax.axhline(ana_ss["C_poly_int"], color="steelblue", ls=":", lw=1, alpha=0.5)
    ax.axhline(ana_hs["C_metal_int"], color="darkorange", ls=":", lw=1, alpha=0.5)

    ax.set_xlabel("Position  [mm]")
    ax.set_ylabel("C_H2  [mol/m³]")
    ax.set_title("Concentration profiles at steady state\n"
                 "Henry-Sieverts: jump up at interface; Sieverts-Sieverts: continuous")
    ax.legend(fontsize=9)
    ax.grid(True, ls=":")

    # 2. Flux history
    ax = axes[1]
    for label, case_dir, ana in [
        ("sieverts_sieverts", ss_dir, ana_ss),
        ("henry_sieverts",    hs_dir, ana_hs),
    ]:
        times = available_times(case_dir, "metal", "C_H2")
        if not times:
            continue
        t_arr = np.array(sorted(times, key=float), dtype=float)
        J_arr = []
        for t_str in sorted(times, key=float):
            C_m = read_of_field(os.path.join(case_dir, t_str, "metal", "C_H2"))
            J_arr.append(D * C_m[-1] / (dx_metal / 2.0))
        J_arr = np.array(J_arr)
        c = colors[label]; ls = styles[label]
        ax.plot(t_arr, J_arr, color=c, ls=ls, lw=2, marker="o", ms=4,
                label=label.replace("_", "-"))
        ax.axhline(ana["J"], color=c, ls=":", lw=1, alpha=0.7,
                   label=f"analytical SS ({label[:2]})")

    ax.set_xlabel("Time  [s]")
    ax.set_ylabel("Permeation flux  [mol/(m²·s)]")
    ax.set_title("Flux vs time: Henry-Sieverts → 23.6% higher steady-state flux\n"
                 "because Sieverts metal holds more gas at same partial pressure")
    ax.legend(fontsize=8)
    ax.grid(True, ls=":")

    fig.tight_layout()
    out = os.path.join(os.path.dirname(__file__), "..", "henry_comparison.png")
    out = os.path.normpath(out)
    fig.savefig(out, dpi=150)
    print(f"\nPlot saved to {out}")
    sys.exit(0 if passed_all else 1)


if __name__ == "__main__":
    main()
