# case329 — Hydrogen Permeation from a FLiBe Channel into SS316

**Solver:** `speciesFluid` (fluid region), `speciesSolid` (solid region)
**Geometry:** Two-region 2-D slab — FLiBe half-channel | SS316 permeation barrier
**Coupling:** `sievertsCoupledMixed` (linear partition, no concentration jump)
**Reference:** Forcey K.S. et al., *J. Nucl. Mater.* **160** (1988) 89–102; Romatoski R.R. & Hu L.W., *Nucl. Technol.* **205** (2019) 1342–1352

---

## Purpose

This case models dissolved hydrogen transport in a FLiBe molten salt
channel representative of a tritium-breeding liquid-salt cooled reactor
(FHR) primary loop. The channel wall is SS316 stainless steel, and
hydrogen permeates through the wall under a concentration gradient
maintained by vacuum removal at the outer surface.

The case demonstrates three physical phenomena simultaneously:

1. **Volumetric source** — uniform hydrogen production throughout the
   FLiBe volume (S = 1×10⁻⁵ mol/(m³·s)), representing the Li-6(n,α)T
   breeding reaction.
2. **Advection-dominated transport** — the channel Peclet number
   Pe = U·W/D ≈ 1 000 means that hydrogen is swept upward along the
   channel much faster than it diffuses laterally toward the wall.
   A thin concentration boundary layer forms at the permeation wall.
3. **Permeation through a metallic barrier** — a `sievertsCoupledMixed`
   interface enforces linear partition continuity (no concentration jump
   because Ks_fluid = Ks_solid) and drives diffusion across the 2 mm
   SS316 slab toward a vacuum BC at the outer surface.

| Region | Solver | Fields solved |
|---|---|---|
| fluid (FLiBe) | `speciesFluid` | U, p, T, C_H2 |
| solid (SS316) | `speciesSolid` | T, C_H2 |

---

## Physical Model

### Governing equations — FLiBe fluid region

**Continuity:**

```
∇·U = 0
```

**Momentum** (PIMPLE, laminar, Re = 34):

```
∂U/∂t + (U·∇)U = −∇(p/ρ) + ν∇²U
```

**Energy** (isothermal; T = 873 K throughout):

```
∂T/∂t + ∇·(UT) = ∇·(α∇T)
```

No heat sources are present, the walls are set to 873 K, and the inlet
enters at 873 K. The temperature field remains uniform and plays no role
in the concentration solution (Ea = 0 for D in both regions, so D is
temperature-independent).

**Species transport** (dissolved hydrogen C_H2):

```
∂C_H2/∂t + ∇·(U C_H2) = ∇·(D_f ∇C_H2) + S

D_f = 5×10⁻⁸ m²/s  (accelerated from literature value — see note below)
S   = 1×10⁻⁵ mol/(m³·s)  (uniform volumetric source)
```

### Governing equations — SS316 solid region

**Species transport** (pure diffusion, no flow):

```
∂C_H2/∂t = ∇·(D_s ∇C_H2)

D_s = 5×10⁻⁹ m²/s  (accelerated — see note below)
```

There is no hydrogen source in the solid. The energy equation is solved
(conduction only) but T remains at 873 K with no gradients.

### Interface coupling — Sieverts partition

At the fluid-solid interface (x = 5 mm) the `sievertsCoupledMixed` BC
enforces concentration continuity for the case of equal Sieverts
solubility constants:

```
Ks_fluid = Ks_solid = 1000 mol/m³  →  C_fluid|_interface = C_solid|_interface
```

Flux continuity is also enforced:

```
D_f · (∂C_H2/∂n)|_fluid = D_s · (∂C_H2/∂n)|_solid
```

For the general case with different solubility constants (e.g. when using
the square-root Sieverts law for diatomic hydrogen), the BC enforces:

```
C_fluid|_interface / Ks_fluid = C_solid|_interface / Ks_solid
```

Here both Ks are set equal so the jump ratio is 1:1.

### Physical model for the diffusivity values

**Literature values at 873 K:**

- FLiBe: Fukada et al. report D_H ≈ 5×10⁻⁹ m²/s (distinct from D_H2)
- SS316: Forcey 1988 Arrhenius: D₀ = 5.08×10⁻⁷ m²/s, Ea = 71 334 J/mol
  → D(873 K) ≈ 2.7×10⁻¹¹ m²/s

With the Forcey SS316 diffusivity the slab diffusion time is
τ = d²/D = (2×10⁻³)² / 2.7×10⁻¹¹ ≈ 1.5×10⁸ s, making steady permeation
inaccessible in a tutorial simulation. The tutorial therefore uses
accelerated values:

| Region | Tutorial D [m²/s] | Physical D [m²/s] | Acceleration factor |
|---|---|---|---|
| FLiBe fluid | 5×10⁻⁸ | ~5×10⁻⁹ | 10× |
| SS316 solid | 5×10⁻⁹ | ~2.7×10⁻¹¹ | ~185× |

With D_s = 5×10⁻⁹ m²/s the slab transient time is τ = d²/D_s =
(2×10⁻³)² / 5×10⁻⁹ = 800 s. The 2000 s simulation run (endTime = 2000)
captures approximately 2.5 slab time constants and is close to
steady-state in the solid.

---

## Geometry and Mesh

```
x = 0 (centerline, symmetry)
│←───── W = 5 mm (half-channel) ───────→│←── d = 2 mm (SS316) ──→│ x = 7 mm
│                                        │                          │
│      FLiBe half-channel                │     SS316 slab           │
│      speciesFluid                      │     speciesSolid          │
│      U in +y, Re = 34                  │     (no flow)             │
│      H2 source S = 1e-5 mol/(m³·s)     │     vacuum at x = 7 mm   │
│                                        │                          │
y = 0 (inlet)                                               y = 50 mm (outlet)
```

| Region | x extent | Cells (x) | Cells (y) | Cell size |
|---|---|---|---|---|
| fluid (FLiBe) | 0 → 5 mm | 25 | 100 | Δx = 0.20 mm, Δy = 0.50 mm |
| solid (SS316) | 5 → 7 mm | 10 | 100 | Δx = 0.20 mm, Δy = 0.50 mm |

Total cells: (25 + 10) × 100 = **3 500**

Both regions are meshed with uniform (ungraded) cell distributions. The
concentration boundary layer in the fluid at the permeation wall has
thickness δ_C ≈ W · Pe^(−1/3) ≈ 5 × 10^(−3) = 0.5 mm (Pe = 1 000), which
spans approximately 2–3 cells (Δx = 0.20 mm). For a sharper boundary layer
profile, the fluid mesh would need grading toward x = 5 mm.

---

## Material Properties

### FLiBe at 873 K (600°C)

| Property | Symbol | Value | Units | Source |
|---|---|---|---|---|
| Density | ρ | 1940 | kg/m³ | Romatoski & Hu 2019 |
| Dynamic viscosity | µ | 5.63×10⁻³ | Pa·s | Romatoski & Hu 2019 |
| Kinematic viscosity | ν | 2.9×10⁻⁶ | m²/s | µ/ρ |
| Specific heat | Cp | 2386 | J/(kg·K) | Romatoski & Hu 2019 |
| Thermal conductivity | κ | 1.1 | W/(m·K) | Romatoski & Hu 2019 |
| Thermal diffusivity | α = κ/(ρCp) | 2.38×10⁻⁷ | m²/s | derived |
| H2 diffusivity (tutorial) | D_f | 5×10⁻⁸ | m²/s | accelerated |

### SS316 at 873 K (600°C)

| Property | Symbol | Value | Units |
|---|---|---|---|
| Density | ρ | 7960 | kg/m³ |
| Specific heat | Cv | 530 | J/(kg·K) |
| Thermal conductivity | κ | 21 | W/(m·K) |
| H2 diffusivity (tutorial) | D_s | 5×10⁻⁹ | m²/s |

---

## Key Dimensionless Numbers

| Number | Formula | Value |
|---|---|---|
| Reynolds Re = U·Dh/ν | 0.01 × 0.01 / 2.9×10⁻⁶ | 34 |
| Hydraulic diameter Dh | 2 × W (half-channel, symmetry) | 10 mm |
| Prandtl Pr = ν/α | 2.9×10⁻⁶ / 2.38×10⁻⁷ | 12.2 |
| Hydrodynamic entry length L_h ≈ 0.04·Re·Dh | 0.04 × 34 × 10 mm | 13.6 mm |
| Schmidt number Sc = ν/D_f | 2.9×10⁻⁶ / 5×10⁻⁸ | 58 |
| Species Peclet number Pe = U·W/D_f | 0.01 × 5×10⁻³ / 5×10⁻⁸ | 1000 |
| Species entry length L_C ≈ 0.04·Re·Sc·Dh | 0.04 × 34 × 58 × 10 mm | 789 mm |
| Concentration BL thickness δ_C ≈ W·Pe^(−1/3) | 5 mm × (1000)^(−1/3) | ~0.5 mm |

**Pe = 1 000** means that axial advection of hydrogen is 1 000 times
faster than lateral diffusion. The concentration field is strongly
advection-dominated: hydrogen produced near the bottom of the channel is
swept upward before it can diffuse to the wall, so the concentration
builds along the channel height. The thin boundary layer (δ_C ≈ 0.5 mm,
~2–3 mesh cells) at the permeation wall is where the diffusive flux into
the solid originates.

**Permeation penetration length:** The characteristic length over which
the fluid drives a significant permeation flux is:

```
L_perm ≈ U · d · W / D_s  = 0.01 × 2×10⁻³ × 5×10⁻³ / 5×10⁻⁹  ≈ 20 m
```

This is much larger than the channel length (50 mm), confirming that the
solid permeation acts as a weak sink: only a small fraction of the
produced hydrogen escapes through the wall in one pass.

---

## Boundary Conditions

### Fluid region (FLiBe)

| Patch | U | p [m²/s²] | T [K] | C_H2 [mol/m³] |
|---|---|---|---|---|
| `fluid_bottom` (y=0, inlet) | `fixedValue (0 0.01 0)` | `zeroGradient` | `fixedValue 873` | `fixedValue 0` |
| `fluid_top` (y=50mm, outlet) | `zeroGradient` | `fixedValue 0` | `zeroGradient` | `zeroGradient` |
| `fluid_centerline` (x=0) | `slip` | `zeroGradient` | `zeroGradient` | `zeroGradient` |
| `fluid_to_solid` (x=5mm) | `noSlip` | `zeroGradient` | `zeroGradient` | `sievertsCoupledMixed` |

**Note on `fluid_bottom` C_H2 BC:** The inlet uses `fixedValue 0` to
represent fresh, hydrogen-free FLiBe entering the channel. Using
`zeroGradient` at an inlet does not constrain C and allows the interior
concentration to propagate upstream (in the Euler discretisation),
inflating the bulk concentration above the correct source-advection
balance. Always use `fixedValue 0` at an inlet when the incoming stream
carries zero species.

**Note on `fluid_centerline` U BC:** The `slip` BC enforces zero normal
velocity (U_x = 0) and zero normal gradient of the tangential velocity
(∂U_y/∂x = 0), which is the correct symmetry condition for a half-channel.
The resulting velocity profile after development is a half-parabola:

```
u_y(x) = 3 · U_mean · (1 − (x/W)²)    (symmetry at x=0)
u_max  = 3 · U_mean  (at x = 0)
```

### Solid region (SS316)

| Patch | T [K] | C_H2 [mol/m³] |
|---|---|---|
| `solid_bottom` (y=0) | `zeroGradient` | `zeroGradient` |
| `solid_top` (y=50mm) | `zeroGradient` | `zeroGradient` |
| `solid_outer` (x=7mm) | `fixedValue 873` | `fixedValue 0` |
| `solid_to_fluid` (x=5mm) | `fixedValue 873` | `sievertsCoupledMixed` |

The outer surface BC `C_H2 = 0` represents perfect vacuum removal of
permeated hydrogen (e.g. by a sweep gas or vacuum pump). This maximises
the permeation driving force.

---

## Coupling Architecture

```
FLiBe fluid region              SS316 solid region
──────────────────              ──────────────────
U, p: PIMPLE solved             (no velocity field)
T:    speciesFluid transport    T: speciesSolid conduction
C_H2: advection-diffusion + S   C_H2: pure diffusion (D_s)
      sievertsCoupledMixed ←──────→ sievertsCoupledMixed
      at fluid_to_solid               at solid_to_fluid
```

At each `foamMultiRun` time step:
1. PIMPLE in `fluid` advances U, p to convergence.
2. T equation solved in `fluid` (isothermal; residuals converge in one
   iteration).
3. C_H2 equation solved in `fluid` with the current U and the
   `sievertsCoupledMixed` BC, which reads C_H2 from `solid` at the
   previous time level.
4. C_H2 equation solved in `solid` (diffusion only) with
   `sievertsCoupledMixed` BC reading C_H2 from `fluid`.
5. Outer correctors repeat steps 2–4 (`nOuterCorrectors = 3`) to converge
   the interface coupling within the time step.

---

## Expected Results

### Analytical estimates

**Steady-state bulk concentration (plug flow, no permeation):**

```
C_out ≈ S · L / U_mean = 1×10⁻⁵ × 0.05 / 0.01 = 5×10⁻⁵ mol/m³
```

The CFD result at y = 50 mm should approach this value in the interior
of the channel, with the concentration being lower near x = 5 mm due to
permeation.

**Permeation flux at steady state** (approximate, linear profile through solid):

```
J_perm ≈ D_s · C_interface / d = 5×10⁻⁹ × C_interface / 2×10⁻³
```

If C_interface ≈ C_out = 5×10⁻⁵ mol/m³, then
J_perm ≈ 1.25×10⁻¹⁰ mol/(m²·s) — a very small fraction of the total
produced hydrogen (which integrates to S × W × L = 2.5×10⁻⁸ mol/s per
metre of width).

**Concentration boundary layer thickness:**

```
δ_C ≈ W · Pe^(−1/3) = 5×10⁻³ × (1000)^(−1/3) ≈ 0.5 mm
```

The cross-channel C_H2 profile at y = 50 mm should show nearly flat
interior concentration with a sharp drop to a lower value at x = 5 mm,
across a layer of width ~0.5 mm.

### What to visualise in ParaView

1. **Velocity (U_y):** Half-Poiseuille profile in the fluid. By y = 20 mm
   the profile should be fully parabolic with u_max = 3·U_mean = 0.03 m/s
   at x = 0 (centerline) and u = 0 at x = 5 mm.

2. **Concentration (C_H2) in fluid:** Builds from 0 at y = 0 to ~5×10⁻⁵
   mol/m³ at y = 50 mm in the interior. A thin depletion layer is visible
   at x = 5 mm where permeation removes hydrogen.

3. **Concentration (C_H2) in solid:** Nearly linear profile from C_interface
   at x = 5 mm to 0 at x = 7 mm, with slight curvature during the
   transient. By t = 1600 s the profile should be nearly linear
   (steady-state diffusion).

4. **Integration check:** The surface-integrated permeation flux out of the
   solid outer wall should equal the flux into the solid at x = 5 mm to
   within solver tolerance, confirming mass conservation in the solid.

---

## Running the Case

Build the library if not already done:

```sh
cd /path/to/multiSpeciesRegionFoam
./Allwmake
```

Run the case:

```sh
cd tutorials/case329-flibe-channel
./Allrun
```

The `Allrun` script executes:

1. `blockMesh` — creates the 3 500-cell two-zone mesh
2. `splitMeshRegions -cellZones -overwrite` — splits into `fluid` and
   `solid` regions and creates the `fluid_to_solid` / `solid_to_fluid`
   coupled patch pair
3. `foamMultiRun` — runs for 2000 s (2000 steps of Δt = 1 s)

Estimated runtime: 10–25 minutes on a single core. The `foamMultiRun` log
will show coupled residuals converging at each outer corrector iteration;
the C_H2 residuals in the solid are typically smallest and converge fastest.

**Pre-run check:** The case ships with time directories 0, 200, 400, ...,
2000 already present as reference. Running `Allrun` overwrites them with
freshly computed results. To compare against the reference, copy the
existing directories before running:

```sh
cp -r 2000 2000_reference
./Allrun
```

---

## fvSolution and fvSchemes Details

### Fluid region

| Equation | Solver | Preconditioner | Reason |
|---|---|---|---|
| p | PCG | DIC | Symmetric Laplacian |
| U | PBiCGStab | DILU | Non-symmetric (momentum convection) |
| T (`"T.*"`) | PBiCGStab | DILU | Non-symmetric (advection) |
| C_H2 (`"C_.*"`) | PBiCGStab | DILU | Non-symmetric (Pe = 1000 advection) |

Outer correctors: `nOuterCorrectors = 3` (needed to converge the
`sievertsCoupledMixed` interface coupling within each time step).

### Solid region

| Equation | Solver | Preconditioner |
|---|---|---|
| e (internal energy, from which T is recovered) | PCG | DIC |
| C_H2 (`"C_.*"`) | PCG | DIC |

The solid has no advection; all matrices are symmetric (pure Laplacian).

### Advection scheme for C_H2

```
div(phi,C_H2)   Gauss upwind;
```

With Pe = 1 000, any central-difference or even `linearUpwind` scheme
would produce unphysical oscillations. First-order upwind is mandatory
for high-Pe species transport. The resulting numerical diffusion in x is
D_numerical ≈ u_y · Δy / 2 ≈ 0.01 × 0.5×10⁻³ / 2 ≈ 2.5×10⁻⁶ m²/s
in the axial direction, which is large compared to D_f = 5×10⁻⁸ m²/s.
This is acceptable for capturing the bulk concentration build-up along y
but means the axial species profile has first-order accuracy.

---

## Known Issues and Solver Notes

1. **Accelerated diffusivities.** The tutorial D values (D_f = 5×10⁻⁸ m²/s,
   D_s = 5×10⁻⁹ m²/s) are physically unrealistic accelerations chosen to
   make the transient permeation visible within a 2 000 s run. For a
   physically accurate simulation at 873 K, use the Forcey (1988) Arrhenius
   parameters for SS316 (D₀ = 5.08×10⁻⁷ m²/s, Ea = 71 334 J/mol) and
   literature FLiBe values. The run time would need to be extended to
   10⁷–10⁸ s for the solid to reach steady state.

2. **PBiCGStab for T and C in the fluid.** The advection operator is
   non-symmetric. PCG will abort with "Unknown asymmetric matrix solver
   PCG". Use PBiCGStab/DILU for all advection-containing equations.

3. **"T.*" wildcard.** PIMPLE creates `TFinal` on its last corrector pass.
   The entry in `fvSolution` must be `"T.*"` (with quotes) to match both
   `T` and `TFinal`. A bare `T` entry will cause a fatal error on the
   second corrector.

4. **`slip` vs `symmetry` at the centerline.** OpenFOAM-13 does not have
   a `symmetry` wall type for incompressible flow; `slip` is the correct
   implementation of a symmetry plane (zero normal velocity, zero normal
   gradient of tangential velocity). For the species field,
   `zeroGradient` at x = 0 is equivalent to the symmetry condition of
   zero normal diffusive flux.

5. **`fixedValue 0` at the inlet for C_H2.** Using `zeroGradient` at
   `fluid_bottom` does not anchor the concentration to zero. The Euler
   time scheme propagates the interior C value to the inlet cell via the
   advection term, inflating the apparent inlet concentration. Use
   `fixedValue 0` to properly represent a reservoir of pure (hydrogen-free)
   FLiBe.

6. **Outer correctors.** The `sievertsCoupledMixed` BC reads the neighbour
   concentration at the previous outer iteration. With
   `nOuterCorrectors = 1`, the interface coupling is lagged by one time
   step, which can cause slow convergence of the bulk-to-surface
   concentration gradient. Setting `nOuterCorrectors = 3` eliminates this
   lag within each time step.

---

## Connection to the Framework

case329 demonstrates:

1. **`speciesFluid` + `speciesSolid` two-region coupling** — the pattern
   most relevant to fusion and fission applications where a flowing coolant
   contacts a structural wall through which a species (tritium, hydrogen)
   permeates.
2. **Volumetric source in `speciesFluid`** — the `source` entry in
   `speciesProperties` adds S = 1×10⁻⁵ mol/(m³·s) uniformly to the C_H2
   equation. This represents neutron-driven tritium breeding without any
   special source term coding.
3. **Linear Sieverts partition (Ks_ratio = 1)** — the `sievertsCoupledMixed`
   BC with equal Ks values on both sides enforces concentration continuity
   across the fluid-solid interface. Changing Ks_solid (e.g. to reflect
   a true Sieverts solubility in the metal) would introduce a concentration
   jump while preserving flux continuity.
4. **High-Pe advection-diffusion** — the case illustrates the numerical
   challenge of Pe = 1 000 species transport: upwind is required for C_H2,
   but `linearUpwind` is adequate for T (lower effective Pe in the energy
   equation due to the higher α/D ratio).
5. **Progressive model complexity** — case311–319 cover diffusion in solid
   slabs; case329 adds solved flow (speciesFluid) and a distributed source,
   representing the next level of physical fidelity for fusion-relevant
   hydrogen transport.
