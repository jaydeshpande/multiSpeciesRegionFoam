# case323 — Counter-Flow Heat Exchanger with Hydrogen Permeation

## Overview

A 2-D counter-flow heat exchanger model demonstrating convective hydrogen
transport in FLiBe and Hitec molten salts coupled with diffusion through an
SS316 tube wall.  This extends case320 (1-D cross-section proxy) to include
axial convective flow at turbulent Reynolds numbers, giving a physically
realistic picture of how hydrogen concentration evolves along the channel length.

## Geometry

```
x →   0 mm         20 mm   24 mm                   64 mm
      |————————————|———————|————————————————————————|
FLiBe |            |       |                        |
flow  |  FLiBe     | SS316 |  Hitec channel         |
+y ↑  |  channel   | wall  |                        |  ↓ -y
      | (10 cells) |(20 c) |       (10 cells)       |
      |  20 mm     | 4 mm  |        40 mm           |
      |————————————|———————|————————————————————————|

y →   0 mm (FLiBe inlet / Hitec outlet)  →  500 mm (FLiBe outlet / Hitec inlet)
                          100 cells, uniform (Δy = 5 mm)
```

The x-direction resolves the channel widths and wall thickness; the y-direction
follows the axial flow.  The z-direction is 1 mm with 1 cell (2-D planar).

## Physics

The `speciesSolid` solver detects the `U` field in FLiBe and Hitec regions
and adds the convective term to the species equation:

```
∂C/∂t + ∇·(U C) = ∇·(D_eff(T) ∇C)
```

The SS316 wall has no `U` field and solves pure diffusion only.

**Turbulence model**: k-ω SST bulk-average effective diffusivity.  Since
`speciesSolid` uses a prescribed velocity (not a coupled momentum solver),
turbulence is represented through an enhanced effective diffusivity
`D_eff = ν_t / Sc_t` estimated analytically from the bulk-average k-ω SST
eddy viscosity at the given Re numbers.  This is a standard approach for
membrane permeation modeling with turbulent bulk flow.

Coupling at both fluid–wall interfaces enforces Sieverts-law flux continuity:

```
D_fluid · ∂C_fluid/∂n = D_wall · ∂C_wall/∂n    (flux continuity)
C_fluid / Ks_fluid     = C_wall  / Ks_wall        (equilibrium partition)
```

With equal Ks on all sides (Ks = 1000 mol/m³/Pa^0.5), concentration is
continuous across the interfaces and the only resistance is diffusive.

## Boundary conditions

| Boundary      | C_H2                  | T                   | U               |
|---------------|-----------------------|---------------------|-----------------|
| FLiBe inlet   | fixedValue 1.0 mol/m³ | fixedValue 973 K    | (0, 1.0, 0) m/s |
| FLiBe outlet  | zeroGradient          | zeroGradient        | zeroGradient    |
| FLiBe outer   | zeroGradient          | zeroGradient        | (0, 0, 0)       |
| Hitec inlet   | fixedValue 0          | fixedValue 673 K    | (0, −2.5, 0) m/s|
| Hitec outlet  | zeroGradient          | zeroGradient        | zeroGradient    |
| Hitec outer   | zeroGradient          | zeroGradient        | (0, 0, 0)       |
| Wall ends     | zeroGradient          | zeroGradient        | —               |

The 300 K temperature difference between inlets (FLiBe 973 K, Hitec 673 K) drives
axial heat exchange through the SS316 wall in addition to hydrogen permeation.
Since D_eff in both fluids uses Ea=0, species transport is not affected by the
thermal gradient, but the wall temperature field will vary axially.

## Material parameters

| Region | D_eff (m²/s)              | Ks (mol/m³/Pa^½) | Basis                        |
|--------|---------------------------|------------------|------------------------------|
| FLiBe  | 1×10⁻⁴ (k-ω SST est.)    | 1000             | Re≈31 000, Sc_t=0.85         |
| SS316  | 5×10⁻⁸ (constant)         | 1000             | τ_wall=(4mm)²/D=320 s; see note |
| Hitec  | 2×10⁻⁴ (k-ω SST est.)    | 1000             | Re≈34 000, Sc_t=0.85         |

> **Note on SS316 diffusivity**: The Forcey 1988 Arrhenius parameters
> (D₀ = 5.08×10⁻⁷ m²/s, Ea = 71 334 J/mol) give
> D(973 K) = 7.5×10⁻¹¹ m²/s → τ_wall ≈ 59 hours.
> A constant D = 5×10⁻⁸ m²/s is used here so permeation is visible
> within the 2000 s run.  Switch to Arrhenius for long validation runs.

### k-ω SST effective diffusivity derivation

```
FLiBe: u_bulk = 2 m/s, D_h = 40 mm
  Re   = ρ·u·D_h/μ = 1940×2×0.04/5e-3  ≈ 31 000
  u*   ≈ 0.07·u_bulk/√2                 ≈ 0.099 m/s  (Cf ≈ 0.0049)
  ν_t  ≈ (0.07/2)·u*·(D_h/2)           ≈ 8.5×10⁻⁵ m²/s
  D_eff = ν_t / Sc_t (Sc_t=0.85)       ≈ 1×10⁻⁴ m²/s

Hitec: u_bulk = 2 m/s, D_h = 80 mm
  Re   = ρ·u·D_h/μ = 1680×2×0.08/8e-4  ≈ 34 000
  u*   ≈ 0.07·u_bulk/√2                 ≈ 0.099 m/s
  ν_t  ≈ (0.07/2)·u*·(D_h/2)           ≈ 1.7×10⁻⁴ m²/s
  D_eff = ν_t / Sc_t (Sc_t=0.85)       ≈ 2×10⁻⁴ m²/s
```

## Expected steady-state behaviour

At steady state (t ≈ 1600 s):

- **FLiBe channel**: C decreases from 1.0 mol/m³ at the inlet (y=0) toward
  the outlet (y=500 mm) as hydrogen permeates through the wall.
- **SS316 wall**: Concentration gradient in x at every axial position,
  magnitude decreasing along the axial direction.
- **Hitec channel**: C increases from 0 at the inlet (y=500 mm) toward the
  outlet (y=0) as it picks up permeated hydrogen.

The total permeation flux is reported by `wallFluxFLiBe` and `wallFluxHitec`
function objects in `postProcessing/`.

## Key timescales

| Timescale         | Formula                      | Value     |
|-------------------|------------------------------|-----------|
| Flow-through      | L/u = 0.5/2.0                | 0.25 s    |
| Wall diffusion    | δ²/D = (4e−3)²/5e−8         | 320 s     |
| Steady state      | ≈ 5 × τ_wall                 | ~1600 s   |

## Mesh refinement

Cells are clustered toward each fluid–wall interface to resolve the
concentration boundary layer (the thinnest length scale at high Pe):

| Region | Cells | Width  | Grading          | First cell | Last cell |
|--------|-------|--------|------------------|------------|-----------|
| FLiBe  | 10    | 20 mm  | simpleGrading 0.1| ~3.5 mm    | ~0.35 mm  |
| Wall   | 20    | 4 mm   | uniform          | 0.2 mm     | 0.2 mm    |
| Hitec  | 10    | 40 mm  | simpleGrading 10 | ~0.7 mm    | ~7 mm     |

FLiBe grading = 0.1 → cells shrink toward x = 20 mm (fluid–wall interface).
Hitec grading = 10 → cells shrink toward x = 24 mm (fluid–wall interface).

## Running

```sh
source /opt/openfoam13/etc/bashrc
cd tutorials/case323-counterflow-hx
./Allrun
```

Requires the full library to be built first (`./Allwmake` from the repo root).
The solver change adding optional convection is in `src/speciesSolid/speciesSolid.C`.

## Extending this case

- **Different flow rates**: change `u` in `0/flibe/U` and `0/hitec/U`; update
  D_eff in `constant/*/speciesProperties` to match new Re.
- **Temperature gradient**: set different T BCs on FLiBe/Hitec outer walls
  and use Arrhenius D/Ks for SS316 — the solver will pick up T from the
  coupled patches automatically.
- **Concentration jump at interface**: set Ks_FLiBe ≠ Ks_wall (see case314).
- **Parallel flow**: change Hitec U to `(0 2.0 0)` and swap inlet/outlet BCs.
