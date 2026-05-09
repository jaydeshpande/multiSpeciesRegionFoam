#!/usr/bin/env python3
"""
case328 validation — Poiseuille flow + Graetz heat transfer
============================================================
Tests the speciesFluid solver against two classical analytical results:

  1. POISEUILLE PROFILE (PASS if L2 error < 2%)
     u_y(x) = 6·U_mean·x·(W−x)/W²,   u_max / U_mean = 1.5
     Measured at y = 45 mm (well past L_h = 1.6 mm entry length).

  2. GRAETZ NUSSELT NUMBER (PASS if |Nu_cfd − 7.54| / 7.54 < 5%)
     Nu_∞ = 7.5407  (constant wall-T, both walls heated, Poiseuille profile)
     Computed as the average local Nu over y ∈ [25, 45] mm, where the
     temperature profile is thermally fully developed (L_T ≈ 11 mm).

Mesh assumed:
  Single fluid region "fluid", uniform blockMesh:
    Nx = 40  (x: 0 to 1 mm),  Ny = 500 (y: 0 to 50 mm),  Nz = 1
  Cell index = iy * Nx + ix   (blockMesh ordering, single block, k=0)
  Cell centre x_i = (ix + 0.5) * W/Nx
  Cell centre y_j = (iy + 0.5) * L/Ny

Usage:
  python3 validate_poiseuille_graetz.py [case_dir]
  (default case_dir: parent of this script's directory)
"""

import sys
import os
import numpy as np

if not hasattr(np, "trapezoid"):
    np.trapezoid = np.trapz

# ── Geometry & mesh ──────────────────────────────────────────────────────────

W   = 1e-3      # m   channel width
L   = 50e-3     # m   channel length
Nx  = 40        # cells in x
Ny  = 500       # cells in y

dx  = W / Nx    # 0.025 mm
dy  = L / Ny    # 0.100 mm

x_centers = (np.arange(Nx) + 0.5) * dx   # shape (40,)
y_centers = (np.arange(Ny) + 0.5) * dy   # shape (500,)

# ── Operating conditions ─────────────────────────────────────────────────────

U_mean = 0.01       # m/s
nu     = 1.004e-6   # m²/s
kappa  = 0.598      # W/(m·K)
Dh     = 2 * W      # 2 mm  (parallel plates)
T_in   = 293.0      # K
T_wall = 373.0      # K

# Analytical Poiseuille:  u_y(x) = 6·U·x·(W−x)/W²
u_poiseuille = 6 * U_mean * x_centers * (W - x_centers) / W**2

# Analytical pressure gradient (kinematic):
dp_dy_analytical = -12 * nu * U_mean / W**2   # m/s² per m  (≈ -0.120 m/s²)

# Fully-developed Nusselt (parallel plates, constant wall T, Poiseuille):
Nu_inf = 7.5407

# ── Helper: read OpenFOAM ASCII field ────────────────────────────────────────

def read_OF_scalar(filepath):
    """Return the internalField values from an ASCII OpenFOAM scalar field."""
    values = []
    found  = False
    in_data = False
    with open(filepath) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("internalField"):
                if "uniform" in s:
                    return np.array([float(s.split("uniform")[1].rstrip(";").strip())])
                found = True
                continue
            if found and not in_data:
                if s == "(":
                    in_data = True
                continue
            if in_data:
                if s == ")":
                    break
                try:
                    values.append(float(s))
                except ValueError:
                    pass
    return np.array(values)


def read_OF_vector_component(filepath, component):
    """Return one component (0=x, 1=y, 2=z) of a vector internalField."""
    values = []
    found  = False
    in_data = False
    with open(filepath) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("internalField"):
                found = True
                continue
            if found and not in_data:
                if s == "(":
                    in_data = True
                continue
            if in_data:
                if s == ")":
                    break
                s2 = s.lstrip("(").rstrip(")")
                parts = s2.split()
                if len(parts) == 3:
                    try:
                        values.append(float(parts[component]))
                    except ValueError:
                        pass
    return np.array(values)


def latest_time(case_dir, region="fluid"):
    """Return the latest time directory that contains the region fields."""
    times = []
    for d in os.listdir(case_dir):
        try:
            t = float(d)
            if t > 0 and os.path.isdir(os.path.join(case_dir, d, region)):
                times.append(t)
        except ValueError:
            pass
    if not times:
        raise RuntimeError(f"No time directories with '{region}' found in {case_dir}")
    return max(times)

# ── Main ─────────────────────────────────────────────────────────────────────

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
case_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_SCRIPT_DIR, "..")

t_last = latest_time(case_dir)
t_dir  = os.path.join(case_dir, f"{t_last:g}", "fluid")
print(f"Reading fields from: {t_dir}  (t = {t_last} s)")

# ── Load fields ───────────────────────────────────────────────────────────────

Uy_flat = read_OF_vector_component(os.path.join(t_dir, "U"), component=1)
T_flat  = read_OF_scalar(os.path.join(t_dir, "T"))

if Uy_flat.size != Nx * Ny or T_flat.size != Nx * Ny:
    raise RuntimeError(
        f"Unexpected field size: |U|={Uy_flat.size}, |T|={T_flat.size}, "
        f"expected {Nx*Ny}.  Check Nx={Nx}, Ny={Ny}."
    )

# Reshape to (Ny, Nx): row = y-level, column = x-position
Uy = Uy_flat.reshape(Ny, Nx)
T  = T_flat.reshape(Ny, Nx)

# ── 1. Poiseuille validation ──────────────────────────────────────────────────

# Cross-section at y ≈ 45 mm → iy = 449
iy_outlet = int(0.045 / dy)           # = 449
u_cfd = Uy[iy_outlet, :]              # shape (40,)

l2_err = np.sqrt(np.mean((u_cfd - u_poiseuille)**2)) / U_mean
u_max_cfd = u_cfd.max()
u_max_ratio = u_max_cfd / U_mean

print("\n── 1. Poiseuille profile (y = {:.0f} mm) ──────────────────────".format(
    y_centers[iy_outlet] * 1e3))
print(f"   u_max/U_mean  CFD={u_max_ratio:.4f}  analytical=1.5000")
print(f"   L2 relative error: {l2_err*100:.3f}%")
poiseuille_pass = l2_err < 0.02
print("   " + ("PASS ✓" if poiseuille_pass else "FAIL ✗") +
      f"  (threshold 2%)")

# ── 2. Pressure gradient ─────────────────────────────────────────────────────

p_flat = read_OF_scalar(os.path.join(t_dir, "p"))
p = p_flat.reshape(Ny, Nx)

# Compute cross-section-averaged pressure at each y level
p_mean = p.mean(axis=1)      # shape (Ny,)

# Fit a line to p(y) in the fully developed region (y > 5 mm = 50 cells)
iy_start = 50
y_fit = y_centers[iy_start:]
p_fit = p_mean[iy_start:]
coeffs = np.polyfit(y_fit, p_fit, 1)
dp_dy_cfd = coeffs[0]

dp_err = abs(dp_dy_cfd - dp_dy_analytical) / abs(dp_dy_analytical)
print("\n── 2. Pressure gradient ─────────────────────────────────────────")
print(f"   dp/dy  CFD={dp_dy_cfd:.5f}  analytical={dp_dy_analytical:.5f}  m/s² per m")
print(f"   Relative error: {dp_err*100:.2f}%")
pressure_pass = dp_err < 0.05
print("   " + ("PASS ✓" if pressure_pass else "FAIL ✗") +
      "  (threshold 5%)")

# ── 3. Graetz Nusselt number ──────────────────────────────────────────────────

# Compute local Nusselt number at each y-level using the bulk temperature method.
# T_bulk(y) = ∫₀^W u(x)·T(x) dx / ∫₀^W u(x) dx   (energy-weighted mean)
# q"_left  ≈ -κ · (T[iy,0] - T_wall) / (dx/2)     (gradient at left wall)
# q"_right ≈  κ · (T_wall - T[iy,Nx-1]) / (dx/2)
# h = (q"_left + q"_right) / (2·(T_wall - T_bulk))
# Nu = h · Dh / κ

T_bulk_arr = np.empty(Ny)
Nu_arr     = np.empty(Ny)

for iy in range(Ny):
    u_row = Uy[iy, :]
    T_row = T[iy, :]
    u_int = np.trapezoid(u_row, x_centers)
    if u_int < 1e-30:
        T_bulk_arr[iy] = T_wall
        Nu_arr[iy] = np.nan
        continue
    T_bulk = np.trapezoid(u_row * T_row, x_centers) / u_int
    T_bulk_arr[iy] = T_bulk
    dT = T_wall - T_bulk
    if abs(dT) < 1e-6:
        Nu_arr[iy] = np.nan
        continue
    q_left  = -kappa * (T_row[0]  - T_wall) / (0.5 * dx)
    q_right =  kappa * (T_wall - T_row[-1]) / (0.5 * dx)
    q_mean  = 0.5 * (q_left + q_right)
    h       = q_mean / dT
    Nu_arr[iy] = h * Dh / kappa

# Average Nu in the thermally fully-developed region: y ∈ [25, 45] mm
iy_lo = int(0.025 / dy)   # = 250
iy_hi = int(0.045 / dy)   # = 449
Nu_fd = np.nanmean(Nu_arr[iy_lo:iy_hi+1])
Nu_err = abs(Nu_fd - Nu_inf) / Nu_inf

print("\n── 3. Graetz Nusselt number (y = 25–45 mm, fully developed) ─────")
print(f"   Nu_cfd = {Nu_fd:.4f}   Nu_∞ = {Nu_inf:.4f}")
print(f"   Relative error: {Nu_err*100:.2f}%")
nusselt_pass = Nu_err < 0.05
print("   " + ("PASS ✓" if nusselt_pass else "FAIL ✗") +
      "  (threshold 5%)")

# ── Summary ──────────────────────────────────────────────────────────────────

print("\n══════════════════════════════════════════════════")
all_pass = poiseuille_pass and pressure_pass and nusselt_pass
status = "ALL PASS ✓" if all_pass else "SOME TESTS FAILED ✗"
print(f"  case328 speciesFluid validation: {status}")
print("══════════════════════════════════════════════════\n")

sys.exit(0 if all_pass else 1)
