"""
Shell-and-tube HX permeation validation.

Compares multiSpeciesRegionFoam results against the analytical steady-state
solution for heat and species transport through a three-layer composite wall:

  FLiBe proxy (0.1 mm) | SS316 wall (0.5 mm) | Hitec proxy (0.1 mm)
  x = 0                 x = 0.1 mm             x = 0.6 mm   x = 0.7 mm

Physical parameters
-------------------
  T_FLiBe = 1073 K  (outer BC, x = 0)
  T_Hitec  =  673 K  (outer BC, x = 0.7 mm)
  C_H2_FLiBe = 1.0 mol/m³  (outer BC, x = 0)
  C_H2_Hitec  = 0          (outer BC, x = 0.7 mm)

  k_FLiBe = 1.1 W/(m·K),  k_wall = 15 W/(m·K),  k_Hitec = 0.57 W/(m·K)
  D_FLiBe = D_Hitec = 1e-7 m²/s  (large proxy diffusivities)
  D_wall(T) = 5.08e-7 × exp(-8580/T)  m²/s  (Arrhenius, SS316, Forcey 1988)

Analytical steady-state temperature (piecewise linear per layer):
  J_q = ΔT / R_total  [W/m²]
  T(x) piecewise linear with slope = -J_q / k_layer

Analytical steady-state concentration (non-linear in wall due to varying D):
  J_C = (C_in - 0) / R_C_total  [mol/(m²·s)]
  R_C_wall = ∫ dx / D(T(x))  (numerical)
  R_C_FLiBe = L_FLiBe / D_FLiBe  (negligible)
  R_C_Hitec = L_Hitec / D_Hitec  (negligible)
"""

import os
import re
import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# numpy ≥ 2.0 renamed trapz → trapezoid
try:
    _trapz = np.trapezoid
except AttributeError:
    _trapz = np.trapz

# ── physical parameters ────────────────────────────────────────────────────
R_gas = 8.314          # J/(mol·K)

# Geometry (SI units)
L_flibe = 0.1e-3       # m
L_wall  = 0.5e-3       # m
L_hitec = 0.1e-3       # m

# Thermal properties
k_flibe = 1.1          # W/(m·K)
k_wall  = 15.0         # W/(m·K)
k_hitec = 0.57         # W/(m·K)

# Temperature BCs
T_flibe_outer = 1073.0  # K
T_hitec_outer =  673.0  # K

# Diffusivity BCs
D0_wall  = 5.08e-7     # m²/s  (Forcey 1988)
Ea_wall  = 71334.0     # J/mol (= 8580 K × 8.314)
D_proxy  = 1e-7        # m²/s  (FLiBe and Hitec proxy, temperature-independent)

# Concentration BCs
C_in  = 1.0            # mol/m³  (FLiBe outer face)
C_out = 0.0            # mol/m³  (Hitec outer face)

# Mesh cells
N_flibe = 5
N_wall  = 50
N_hitec = 5


def D_wall_at(T):
    return D0_wall * np.exp(-Ea_wall / (R_gas * T))


# ── Analytical steady-state ────────────────────────────────────────────────
def analytical_steady_state(n_int=5000):
    """Compute analytical SS temperature profile and permeation flux."""
    R_T_flibe = L_flibe / k_flibe
    R_T_wall  = L_wall  / k_wall
    R_T_hitec = L_hitec / k_hitec
    R_T_total = R_T_flibe + R_T_wall + R_T_hitec

    J_q = (T_flibe_outer - T_hitec_outer) / R_T_total

    T_fw = T_flibe_outer - J_q * R_T_flibe  # T at FLiBe/wall interface
    T_wh = T_fw - J_q * R_T_wall            # T at wall/Hitec interface

    # Temperature profile: x is position from 0
    def T_of_x(x):
        if x <= L_flibe:
            return T_flibe_outer - J_q / k_flibe * x
        elif x <= L_flibe + L_wall:
            return T_fw - J_q / k_wall * (x - L_flibe)
        else:
            return T_wh - J_q / k_hitec * (x - L_flibe - L_wall)

    # Species resistance: numerical integration through wall
    x_wall = np.linspace(0, L_wall, n_int)
    T_wall = T_fw - J_q / k_wall * x_wall
    D_wall = D_wall_at(T_wall)

    R_C_flibe = L_flibe / D_proxy
    R_C_wall  = _trapz(1.0 / D_wall, x_wall)
    R_C_hitec = L_hitec / D_proxy
    R_C_total = R_C_flibe + R_C_wall + R_C_hitec

    J_C = (C_in - C_out) / R_C_total     # mol/(m²·s)

    return {
        "J_q":      J_q,
        "T_fw":     T_fw,
        "T_wh":     T_wh,
        "R_T_total": R_T_total,
        "R_C_total": R_C_total,
        "J_C":      J_C,
        "T_of_x":   np.vectorize(T_of_x),
        "x_wall":   x_wall,
        "T_wall":   T_wall,
        "D_wall":   D_wall,
    }


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
    """Permeation flux at the outer wall (last cell → fixedValue=0 BC)."""
    return D_val * C_field[-1] / (dx / 2.0)


# ── main ──────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--case", default=".", help="top-level case directory")
    p.add_argument("--plot", action="store_true", help="save comparison plots")
    args = p.parse_args()

    top = os.path.abspath(args.case)

    ana = analytical_steady_state()

    print("=== Analytical steady-state ===")
    print(f"  Heat flux J_q     = {ana['J_q']:.3e} W/m²")
    print(f"  T at FLiBe/wall   = {ana['T_fw']:.1f} K")
    print(f"  T at wall/Hitec   = {ana['T_wh']:.1f} K")
    print(f"  D(T_fw)           = {D_wall_at(ana['T_fw']):.3e} m²/s")
    print(f"  D(T_wh)           = {D_wall_at(ana['T_wh']):.3e} m²/s")
    print(f"  D ratio hot/cold  = {D_wall_at(ana['T_fw'])/D_wall_at(ana['T_wh']):.2f}")
    print(f"  Permeation flux J_C = {ana['J_C']:.3e} mol/(m²·s)")
    print()

    # ── read simulation results ────────────────────────────────────────────
    times_wall = available_times(top, "wall", "C_H2")
    if not times_wall:
        print("ERROR: no time directories found with wall/C_H2. Run the simulation first.")
        sys.exit(1)

    dx_wall  = L_wall  / N_wall
    dx_flibe = L_flibe / N_flibe
    dx_hitec = L_hitec / N_hitec

    print("Simulated outer-wall permeation flux  J = D_wall × C_last / (Δx/2):")
    print(f"  {'t [s]':>8}   {'J [mol/(m²·s)]':>22}   {'J/J_ana':>10}")

    J_sim_history = {}
    for t_str in times_wall:
        C_w = read_of_field(os.path.join(top, t_str, "wall", "C_H2"))
        T_w = read_of_field(os.path.join(top, t_str, "wall", "T"))
        D_last = D_wall_at(T_w[-1])
        J_sim = outer_flux(C_w, D_last, dx_wall)
        J_sim_history[t_str] = J_sim
        ratio = J_sim / ana['J_C'] if ana['J_C'] > 0 else float('nan')
        print(f"  {float(t_str):>8.0f}   {J_sim:>22.3e}   {ratio:>10.3f}")

    # PASS/FAIL: at last time, flux within 5% of analytical SS
    last_t = sorted(J_sim_history, key=float)[-1]
    rel_err = abs(J_sim_history[last_t] - ana['J_C']) / ana['J_C']
    print()
    print(f"Relative flux error at t = {last_t} s: {rel_err:.3%}")
    passed = rel_err < 0.05
    print("PASS" if passed else "FAIL")

    if not args.plot:
        sys.exit(0 if passed else 1)

    # ── plots ──────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # 1. Flux vs time
    ax = axes[0]
    t_arr  = np.array(sorted(J_sim_history, key=float), dtype=float)
    J_arr  = np.array([J_sim_history[str(int(t))] for t in t_arr])
    ax.plot(t_arr, J_arr, "b-o", ms=5, label="simulation")
    ax.axhline(ana['J_C'], color="r", ls="--", lw=1,
               label=f"analytical SS  {ana['J_C']:.2e}")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Permeation flux  [mol/(m²·s)]")
    ax.set_title("Outer-wall H₂ permeation flux")
    ax.legend(fontsize=8)
    ax.grid(True, ls=":")

    # 2. Temperature profile at last time
    ax = axes[1]
    last_t_str = sorted(J_sim_history, key=float)[-1]
    T_fl = read_of_field(os.path.join(top, last_t_str, "flibe", "T"))
    T_wa = read_of_field(os.path.join(top, last_t_str, "wall",  "T"))
    T_hi = read_of_field(os.path.join(top, last_t_str, "hitec", "T"))

    x_fl = (np.arange(N_flibe) + 0.5) * dx_flibe * 1e3
    x_wa = L_flibe * 1e3 + (np.arange(N_wall)  + 0.5) * dx_wall  * 1e3
    x_hi = (L_flibe + L_wall) * 1e3 + (np.arange(N_hitec) + 0.5) * dx_hitec * 1e3

    ax.plot(x_fl, T_fl, "b-", lw=1.5, label="FLiBe proxy")
    ax.plot(x_wa, T_wa, "k-", lw=1.5, label="SS316 wall")
    ax.plot(x_hi, T_hi, "r-", lw=1.5, label="Hitec proxy")

    # Analytical T profile
    x_ana = np.linspace(0, (L_flibe + L_wall + L_hitec) * 1e3, 500)
    T_ana = ana['T_of_x'](x_ana * 1e-3)
    ax.plot(x_ana, T_ana, "g--", lw=1, label="analytical SS")

    ax.axvline(L_flibe * 1e3, color="gray", ls=":", lw=1)
    ax.axvline((L_flibe + L_wall) * 1e3, color="gray", ls=":", lw=1)
    ax.set_xlabel("Position x  [mm]")
    ax.set_ylabel("Temperature  [K]")
    ax.set_title(f"Temperature profile  (t = {last_t_str} s)")
    ax.legend(fontsize=8)
    ax.grid(True, ls=":")

    # 3. Concentration profile at last time
    ax = axes[2]
    C_fl = read_of_field(os.path.join(top, last_t_str, "flibe", "C_H2"))
    C_wa = read_of_field(os.path.join(top, last_t_str, "wall",  "C_H2"))
    C_hi = read_of_field(os.path.join(top, last_t_str, "hitec", "C_H2"))

    ax.plot(x_fl, C_fl, "b-", lw=1.5, label="FLiBe proxy")
    ax.plot(x_wa, C_wa, "k-", lw=1.5, label="SS316 wall")
    ax.plot(x_hi, C_hi, "r-", lw=1.5, label="Hitec proxy")

    # Analytical SS concentration profile (piecewise from integral)
    # Compute cumulative integral of 1/D from left boundary
    x_w_fine = ana['x_wall']
    inv_D = 1.0 / ana['D_wall']
    cum_R = np.concatenate([[0], np.cumsum(np.diff(x_w_fine) * 0.5 *
                              (inv_D[:-1] + inv_D[1:]))])
    R_C_wall = cum_R[-1]
    # Small proxy resistances included for accuracy
    R_C_fl = L_flibe / D_proxy
    R_C_hi = L_hitec / D_proxy
    R_total = R_C_fl + R_C_wall + R_C_hi
    J_C_check = (C_in - C_out) / R_total

    C_fw = C_in - J_C_check * R_C_fl    # concentration at FLiBe/wall interface
    C_wh = C_fw - J_C_check * cum_R[-1] # at wall/Hitec interface

    C_wall_ana = C_fw - J_C_check * cum_R

    ax.plot(L_flibe * 1e3 + x_w_fine * 1e3, C_wall_ana,
            "g--", lw=1, label="analytical SS (wall)")

    # Reference: isothermal profile at hot-side D (linear)
    C_iso = C_fw - (C_fw - C_wh) * x_w_fine / L_wall
    ax.plot(L_flibe * 1e3 + x_w_fine * 1e3, C_iso,
            "g:", lw=1, label="iso-T profile (D_hot)")

    ax.axvline(L_flibe * 1e3, color="gray", ls=":", lw=1)
    ax.axvline((L_flibe + L_wall) * 1e3, color="gray", ls=":", lw=1)
    ax.set_xlabel("Position x  [mm]")
    ax.set_ylabel("C_H2  [mol/m³]")
    ax.set_title(f"H₂ concentration profile  (t = {last_t_str} s)\n"
                 "solid curve bows right of dashed (lower D on Hitec side)")
    ax.legend(fontsize=8)
    ax.grid(True, ls=":")

    fig.tight_layout()
    out = os.path.join(top, "hx_permeation.png")
    fig.savefig(out, dpi=150)
    print(f"\nPlot saved to {out}")
    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
