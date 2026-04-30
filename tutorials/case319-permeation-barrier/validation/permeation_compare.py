"""
Permeation barrier validation: WC-coated vs uncoated SS316 tube.

Compares the multiSpeciesRegionFoam result against the analytical
steady-state solution for 1-D diffusion through a composite wall.

Geometry
--------
  no_coating:    SS316 slab, thickness L_SS = 0.5 mm
  with_coating:  WC (0 – 0.1 mm) + SS316 (0.1 – 0.6 mm)

Physical parameters  (T = 1073 K, P_H2 = 1 bar = 100 000 Pa)
-------------------
  SS316  D(T) = 5.08e-7 × exp(-71334/(R T))      [Forcey et al. 1988]
  WC     D(T) = 1.00e-7 × exp(-97830/(R T))       [estimated, PRF ≈ 200]
  Ks_SS / Ks_WC = 10  (SS316 dissolves 10× more H than WC)
  C_inner_SS316 = Ks_SS316 × sqrt(P) = 0.651 mol/m³
  C_inner_WC   = Ks_WC   × sqrt(P) = 0.0651 mol/m³

Analytical steady-state permeation flux
---------------------------------------
  no_coating:    J_bare = D_SS × C_inner_SS / L_SS
  with_coating:  J_coat = D_SS × C_SS_if / L_SS
    where  C_WC_if = (D_WC/L_WC × C_inner_WC) / (D_WC/L_WC + D_SS × r / L_SS)
           C_SS_if = r × C_WC_if   (Sieverts jump, r = Ks_SS/Ks_WC = 10)
  PRF = J_bare / J_coat  ≈ 201 for these parameters.
"""

import os
import sys
import re
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── physical constants & parameters ──────────────────────────────────────────
R = 8.314          # J/(mol·K)
T = 1073.0         # K
P = 1.0e5          # Pa

D0_SS = 5.08e-7;   Ea_SS = 71334.0     # Forcey 1988
D0_WC = 1.00e-7;   Ea_WC = 97830.0

L_SS = 0.5e-3      # m  (SS316 wall thickness, both cases)
L_WC = 0.1e-3      # m  (WC coating, with_coating only)
r    = 10.0        # Ks_SS316 / Ks_WC

N_SS  = 100        # number of cells in SS316 (both cases)
N_WC  = 20         # number of cells in WC

C_in_SS = 0.651    # mol/m³  (fixedValue inner BC, no_coating)
C_in_WC = 0.0651   # mol/m³  (fixedValue inner BC, with_coating WC side)


def D(D0, Ea, temp=T):
    return D0 * np.exp(-Ea / (R * temp))


def analytical_steady_state():
    """Return (J_bare, J_coat, PRF, C_WC_if, C_SS_if)."""
    d_ss = D(D0_SS, Ea_SS)
    d_wc = D(D0_WC, Ea_WC)
    a = d_wc / L_WC
    b = d_ss * r / L_SS
    C_WC_if = a * C_in_WC / (a + b)
    C_SS_if = r * C_WC_if
    J_bare  = d_ss * C_in_SS / L_SS
    J_coat  = d_ss * C_SS_if / L_SS
    PRF = J_bare / J_coat
    return J_bare, J_coat, PRF, C_WC_if, C_SS_if, d_ss, d_wc


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


def outer_flux(C_field, D_val, dx):
    """Permeation flux at the outer wall (last cell → fixedValue=0 wall)."""
    return D_val * C_field[-1] / (dx / 2.0)


# ── main ──────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--case", default=".", help="top-level case directory")
    p.add_argument("--plot", action="store_true", help="save comparison plots")
    args = p.parse_args()

    top = os.path.abspath(args.case)
    nc_dir = os.path.join(top, "no_coating")
    wc_dir = os.path.join(top, "with_coating")

    J_bare_ana, J_coat_ana, PRF_ana, C_WC_if, C_SS_if, d_ss, d_wc = \
        analytical_steady_state()

    dx_ss = L_SS / N_SS    # 5 μm
    dx_wc = L_WC / N_WC    # 5 μm

    # ── no_coating flux history ─────────────────────────────────────────────
    times_nc = available_times(nc_dir, "steel", "C_H2")
    J_nc = {}
    for t_str in times_nc:
        C = read_of_field(os.path.join(nc_dir, t_str, "steel", "C_H2"))
        J_nc[t_str] = outer_flux(C, d_ss, dx_ss)

    # ── with_coating flux history ───────────────────────────────────────────
    times_wc = available_times(wc_dir, "steel", "C_H2")
    J_wc = {}
    for t_str in times_wc:
        C = read_of_field(os.path.join(wc_dir, t_str, "steel", "C_H2"))
        J_wc[t_str] = outer_flux(C, d_ss, dx_ss)

    # ── console report ──────────────────────────────────────────────────────
    print(f"Analytical steady-state (T = {T} K, P = {P/1e5:.0f} bar)")
    print(f"  D_SS316  = {d_ss:.3e} m²/s")
    print(f"  D_WC     = {d_wc:.3e} m²/s  ({d_ss/d_wc:.0f}× lower)")
    print(f"  Ks ratio = {r:.0f}  (SS316/WC)")
    print(f"  J_bare   = {J_bare_ana:.3e} mol/(m²·s)  [no coating, steady state]")
    print(f"  J_coat   = {J_coat_ana:.3e} mol/(m²·s)  [WC coated,  steady state]")
    print(f"  PRF_ana  = {PRF_ana:.1f}")
    print()

    print("Simulated outer-wall flux  J(t) = D_SS × C_last / (Δx/2):")
    print(f"  {'t [s]':>8}   {'J_bare [mol/(m²·s)]':>22}   {'J_coat [mol/(m²·s)]':>22}   {'PRF_sim':>8}")
    all_t = sorted(set(list(J_nc.keys()) + list(J_wc.keys())), key=float)
    for t_str in all_t:
        jn = J_nc.get(t_str, float("nan"))
        jw = J_wc.get(t_str, float("nan"))
        prf = jn / jw if jw > 0 else float("inf")
        print(f"  {float(t_str):>8.0f}   {jn:>22.3e}   {jw:>22.3e}   {prf:>8.1f}")

    # PASS/FAIL: at last time, no_coating flux within 5% of analytical SS
    last_nc = sorted(J_nc, key=float)[-1]
    rel_err = abs(J_nc[last_nc] - J_bare_ana) / J_bare_ana
    print()
    print(f"no_coating at t = {last_nc} s:  rel error vs analytical SS = {rel_err:.3%}")
    passed = rel_err < 0.05
    print("PASS" if passed else "FAIL")

    if not args.plot:
        sys.exit(0 if passed else 1)

    # ── plots ───────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # — flux vs time ——
    ax = axes[0]
    t_nc  = np.array(sorted(J_nc.keys(), key=float), dtype=float)
    t_wc_ = np.array(sorted(J_wc.keys(), key=float), dtype=float)
    ax.semilogy(t_nc,  [J_nc[str(int(t))] for t in t_nc],
                "b-o", ms=5, label="no coating (SS316)")
    ax.semilogy(t_wc_, [J_wc[str(int(t))] for t in t_wc_],
                "r-s", ms=5, label="WC + SS316")
    ax.axhline(J_bare_ana, color="b", ls="--", lw=1, label=f"SS316 SS  {J_bare_ana:.2e}")
    ax.axhline(J_coat_ana, color="r", ls="--", lw=1, label=f"WC+SS316 SS  {J_coat_ana:.2e}")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Permeation flux  [mol/(m²·s)]")
    ax.set_title(f"Outer-wall H₂ permeation flux\nPRF (analytical SS) = {PRF_ana:.0f}×")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", ls=":")

    # — concentration profiles at last time ——
    ax = axes[1]
    last_nc_t = sorted(J_nc, key=float)[-1]
    last_wc_t = sorted(J_wc, key=float)[-1]

    C_nc = read_of_field(os.path.join(nc_dir, last_nc_t, "steel", "C_H2"))
    C_wc_wc  = read_of_field(os.path.join(wc_dir, last_wc_t, "wc", "C_H2"))
    C_wc_ss  = read_of_field(os.path.join(wc_dir, last_wc_t, "steel", "C_H2"))

    x_nc  = (np.arange(N_SS) + 0.5) * dx_ss * 1e3       # mm
    x_wc  = (np.arange(N_WC) + 0.5) * dx_wc * 1e3       # mm
    x_ss2 = L_WC * 1e3 + (np.arange(N_SS) + 0.5) * dx_ss * 1e3  # mm

    ax.plot(x_nc, C_nc, "b-", lw=1.5, label=f"no coating (t={last_nc_t} s)")
    ax.plot(x_wc, C_wc_wc, "r-", lw=1.5, label=f"WC (t={last_wc_t} s)")
    ax.plot(x_ss2, C_wc_ss, "r-", lw=1.5)

    # Analytical SS profiles
    x_nc_a  = np.linspace(0, L_SS * 1e3, 200)
    x_wc_a  = np.linspace(0, L_WC * 1e3, 100)
    x_ss_a  = np.linspace(L_WC * 1e3, (L_WC + L_SS) * 1e3, 200)
    C_nc_a  = C_in_SS * (1 - x_nc_a / (L_SS * 1e3))
    C_wc_a  = C_in_WC + (C_WC_if - C_in_WC) * x_wc_a / (L_WC * 1e3)
    C_ss_a  = C_SS_if * (1 - (x_ss_a - L_WC * 1e3) / (L_SS * 1e3))

    ax.plot(x_nc_a,  C_nc_a, "b--", lw=1, label="no coating (SS)")
    ax.plot(x_wc_a,  C_wc_a, "r--", lw=1, label="WC+SS316 (SS)")
    ax.plot(x_ss_a,  C_ss_a, "r--", lw=1)

    ax.axvline(L_WC * 1e3, color="k", ls=":", lw=1, label="WC/SS316 interface")
    ax.set_xlabel("Position x  [mm]")
    ax.set_ylabel("C_H2  [mol/m³]")
    ax.set_title("Concentration profiles\n(solid = simulation, dashed = analytical SS)")
    ax.legend(fontsize=8)
    ax.grid(True, ls=":")

    fig.tight_layout()
    out = os.path.join(top, "permeation_comparison.png")
    fig.savefig(out, dpi=150)
    print(f"\nPlot saved to {out}")
    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
