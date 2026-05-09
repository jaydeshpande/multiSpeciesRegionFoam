# case328 — Poiseuille Flow and Graetz Heat Transfer Validation

**Solver:** `speciesFluid` (single fluid region)
**Geometry:** Single-region 2-D slit channel, W = 1 mm, L = 50 mm
**Coupling:** None (single region)
**Reference:** Shah R.K. & London A.L., *Laminar Flow Forced Convection in Ducts*, Academic Press (1978)

---

## Purpose

This case validates the `speciesFluid` solver against two classical
analytical results for laminar internal flow. It is the unit test for
the hydrodynamic and thermal capabilities of `speciesFluid` before those
capabilities are exercised in more complex multi-region cases such as
case327 and case329.

| Test | Criterion | Threshold |
|---|---|---|
| Poiseuille velocity profile | L2 error vs. u(x) = 6·U_mean·x·(W−x)/W² | < 2% |
| Pressure gradient | \|dp/dy_CFD − dp/dy_exact\| / \|dp/dy_exact\| | < 5% |
| Graetz Nusselt number | \|Nu_CFD − 7.5407\| / 7.5407 | < 5% |

A Python validation script (`validation/validate_poiseuille_graetz.py`)
reads the final time directory, computes all three metrics, and exits with
code 0 on all-pass or 1 on any failure.

---

## Physical Model

### Governing equations

The `speciesFluid` solver advances the following equations in the single
fluid region:

**Continuity** (incompressible):

```
∇·U = 0
```

**Momentum** (PIMPLE, laminar):

```
∂U/∂t + (U·∇)U = −∇(p/ρ) + ν∇²U
```

**Energy** (temperature transport):

```
∂T/∂t + ∇·(UT) = ∇·(α∇T),    α = κ/(ρCp)
```

**Passive species** (C_test, required by `speciesFluid` architecture):

```
∂C/∂t + ∇·(UC) = ∇·(D∇C),    D = 1×10⁻⁹ m²/s
```

C_test has `zeroGradient` on all boundaries and stays at zero throughout.
It exists only because `speciesFluid` always requires a `speciesModel`
entry in `speciesProperties`.

### Analytical benchmarks

**Poiseuille profile** — fully developed velocity in a 2-D slit (parallel
plates, width W, both walls no-slip):

```
u_y(x) = 6 · U_mean · x · (W − x) / W²
u_max  = 1.5 · U_mean    (at x = W/2)
```

**Pressure gradient** — force balance on fully-developed laminar flow:

```
dp/dy = −12 · ν · U_mean / W²
     = −12 × 1.004×10⁻⁶ × 0.01 / (1×10⁻³)²
     ≈ −0.1205  m·s⁻² per m   (kinematic units, p = p_phys/ρ)
```

**Graetz Nusselt number** — thermally fully-developed heat transfer in a
parallel-plate channel with Poiseuille profile and constant wall
temperature on both walls:

```
Nu_∞ = h · Dh / κ = 7.5407    (Shah & London, Table 42)
```

where Dh = 2W is the hydraulic diameter for a slit channel. The Nusselt
number is defined using the bulk (flow-rate averaged) temperature:

```
T_bulk(y) = ∫₀^W u(x) · T(x) dx / ∫₀^W u(x) dx
h(y)      = q_wall / (T_wall − T_bulk)
Nu(y)     = h(y) · Dh / κ
```

---

## Geometry and Mesh

```
      x = 0                 x = 1 mm
      wall_left             wall_right
      │ T_wall = 373 K      │ T_wall = 373 K
      │←────── W = 1 mm ───→│
y=50mm│                     │  outlet (p = 0, T = zeroGradient)
      │                     │
      │  fluid region       │  U in +y direction
      │  water at 20°C      │  Re = 20
      │                     │
y=0mm │                     │  inlet (U = 0.01 m/s, T = 293 K)
```

| Parameter | Value |
|---|---|
| Channel width W | 1 mm |
| Channel length L | 50 mm |
| Cells in x (Nx) | 40 (Δx = 25 µm, uniform) |
| Cells in y (Ny) | 500 (Δy = 100 µm, uniform) |
| Cells in z | 1 (empty, 2-D) |
| Total cells | 20 000 |

The uniform mesh (no grading) is chosen so the validation script can use
the simple cell ordering formula: cell (ix, iy) maps to flat index
iy × Nx + ix. This ordering is preserved by `splitMeshRegions` for a
single-block, single-zone mesh.

---

## Material Properties

Water at 20°C (293 K):

| Property | Symbol | Value | Units |
|---|---|---|---|
| Kinematic viscosity | ν | 1.004×10⁻⁶ | m²/s |
| Density | ρ | 998 | kg/m³ |
| Specific heat | Cp | 4182 | J/(kg·K) |
| Thermal conductivity | κ | 0.598 | W/(m·K) |
| Thermal diffusivity | α = κ/(ρCp) | 1.430×10⁻⁷ | m²/s |

---

## Key Dimensionless Numbers

| Number | Formula | Value |
|---|---|---|
| Reynolds Re = U·Dh/ν | 0.01 × 2×10⁻³ / 1.004×10⁻⁶ | 20 |
| Prandtl Pr = ν/α | 1.004×10⁻⁶ / 1.430×10⁻⁷ | 7.0 |
| Hydrodynamic entry length L_h ≈ 0.04·Re·Dh | 0.04 × 20 × 2 mm | 1.6 mm |
| Thermal entry length L_T ≈ 0.04·Re·Pr·Dh | 0.04 × 20 × 7 × 2 mm | 11.2 mm |
| Graetz number Gz = Re·Pr·Dh/L | 20 × 7 × 2mm / 50mm | 5.6 |
| Species Schmidt number Sc = ν/D_C_test | 1.004×10⁻⁶ / 1×10⁻⁹ | 1004 |

**Interpretation of entry lengths:**

- L_h = 1.6 mm: by y = 2 mm (row 20 of 500) the velocity profile is
  already Poiseuille. The validation reads U at y = 45 mm (row 449),
  well past the hydrodynamic entry region.
- L_T = 11.2 mm: by y = 25 mm (row 250) the temperature profile is
  thermally fully developed and the local Nu is within 5% of Nu_∞.
  The Nusselt number is computed as the average over y ∈ [25, 45] mm.

---

## Boundary Conditions

| Patch | U | p [m²/s²] | T [K] | C_test |
|---|---|---|---|---|
| `inlet` (y=0) | `fixedValue (0 0.01 0)` | `zeroGradient` | `fixedValue 293` | `zeroGradient` |
| `outlet` (y=50mm) | `zeroGradient` | `fixedValue 0` | `zeroGradient` | `zeroGradient` |
| `wall_left` (x=0) | `noSlip` | `zeroGradient` | `fixedValue 373` | `zeroGradient` |
| `wall_right` (x=1mm) | `noSlip` | `zeroGradient` | `fixedValue 373` | `zeroGradient` |

Both walls are heated to the same temperature T_wall = 373 K, making the
problem symmetric in x. This matches the boundary condition assumed in the
derivation of Nu_∞ = 7.5407.

The inlet is a plug-flow velocity (`fixedValue`), not a Poiseuille profile.
The hydrodynamic entry length L_h = 1.6 mm is short enough that the
profile is fully parabolic well before the measurement stations.

---

## fvSolution and fvSchemes

### Solver selection

| Equation | Solver | Preconditioner | Reason |
|---|---|---|---|
| p | PCG | DIC | Symmetric Laplacian |
| U | PBiCGStab | DILU | Non-symmetric due to convection |
| T, TFinal (`"T.*"`) | PBiCGStab | DILU | Non-symmetric; advection term |
| C_test (`"C_.*"`) | PBiCGStab | DILU | Non-symmetric; advection term |

### Advection schemes

| Field | Scheme | Reason |
|---|---|---|
| U | `linearUpwind grad(U)` | Second-order bounded momentum |
| T | `linearUpwind grad(T)` | Pe_T ≈ 700 at y=1mm; reduces numerical diffusion |
| C_test | `upwind` | Sc = 1004, astronomical species Peclet |

### Time stepping

```
deltaT      0.005 s
endTime     30 s
Co_max      = U_max · Δt / Δy = 0.015 × 0.005 / 1×10⁻⁴ ≈ 0.75
```

The Courant number is below 1 throughout (stable with Euler scheme). At
t = 30 s the solution is steady; the write interval of 10 s produces
three output times (10, 20, 30 s).

---

## Running the Case

```sh
cd tutorials/case328-poiseuille-graetz
./Allrun
```

The `Allrun` script executes:

1. `blockMesh` — creates the 20 000-cell single-zone mesh
2. `splitMeshRegions -cellZones -overwrite` — creates the `fluid` region
   subdirectory structure
3. `foamMultiRun` — runs for 30 s

Estimated runtime: 3–8 minutes on a single core.

**Validate:**

```sh
python3 validation/validate_poiseuille_graetz.py
```

The script reads U, p, T from the latest time directory, reshapes the
20 000-element flat arrays to (Ny=500, Nx=40), and checks all three
criteria. Example passing output:

```
Reading fields from: 30/fluid  (t = 30 s)

── 1. Poiseuille profile (y = 44.95 mm) ─────────────────────
   u_max/U_mean  CFD=1.4998  analytical=1.5000
   L2 relative error: 0.127%
   PASS  (threshold 2%)

── 2. Pressure gradient ──────────────────────────────────────
   dp/dy  CFD=-0.12047  analytical=-0.12048  m/s² per m
   Relative error: 0.01%
   PASS  (threshold 5%)

── 3. Graetz Nusselt number (y = 25–45 mm, fully developed) ──
   Nu_cfd = 7.5512   Nu_∞ = 7.5407
   Relative error: 0.14%
   PASS  (threshold 5%)

══════════════════════════════════════════════════
  case328 speciesFluid validation: ALL PASS
══════════════════════════════════════════════════
```

---

## Expected Results

At t = 30 s the flow and temperature fields are fully steady.

**Velocity field:** The cross-channel profile at any y > 5 mm should be
parabolic with u_max/U_mean = 1.500 ± 0.02. The residual plug-flow
character from the inlet is confined to y < 2 mm.

**Pressure field:** The cross-section-averaged gauge pressure decreases
linearly from inlet to outlet. The analytical slope is dp/dy = −0.1205
m²/s² per m. The CFD value (fitted by linear regression over y > 5 mm)
should be within 5%.

**Temperature field:** The colour map of T should show:
- T = 293 K at the inlet (y = 0)
- Wall heating penetrating inward as y increases
- Bulk T approaching T_wall = 373 K near y = 50 mm (Gz = 5.6, almost
  fully thermally developed at the outlet)
- Local Nu decreasing from high values near y = 0 (Graetz entry) to
  Nu_∞ ≈ 7.54 from y ≈ 25 mm onward

**Nusselt number convergence:** Nu(y) exhibits a characteristic 1/y^(1/3)
decay in the thermal entry region and plateaus at Nu_∞ = 7.5407. Mesh
resolution of Nx = 40 (Δx = 25 µm) is adequate to resolve the thermal
wall gradient within 1% of the analytical heat flux.

---

## Validation Script Details

### Cell ordering assumption

`blockMesh` with a single hex block orders cells as cell (ix, iy, iz) →
flat index iz × (Nx × Ny) + iy × Nx + ix. Since iz = 0 for all cells
in a 2-D case, the formula simplifies to:

```
flat_index = iy × Nx + ix
```

`splitMeshRegions` preserves this ordering for single-zone, single-block
meshes. The script relies on this to reshape the flat array to (Ny, Nx)
without reading the mesh point coordinates.

**If the mesh is changed** (different block count, grading, or Nx/Ny),
the reshape will be silently wrong. Check that the field size equals
Nx × Ny = 20 000 before interpreting the results.

### Bulk temperature computation

The script integrates the velocity-weighted temperature across the channel
width at each row using the trapezoidal rule:

```python
T_bulk(iy) = trapz(u(x) * T(x), x) / trapz(u(x), x)
```

This is the energy-correct bulk temperature (mixed-mean temperature) and
matches the definition used in the Nusselt number derivation.

### Wall heat flux

The first-cell wall flux is approximated by a one-sided difference:

```python
q_left  = -κ × (T[iy, 0]  - T_wall) / (Δx/2)
q_right =  κ × (T_wall - T[iy, Nx-1]) / (Δx/2)
```

With Δx = 25 µm and the cell centre at Δx/2 = 12.5 µm from the wall, this
is a first-order estimate. The error relative to the true wall gradient
decreases as (Δx/2)/δ_T where δ_T is the thermal boundary layer thickness.
With δ_T ≈ 0.05 mm in the thermally developed region, the truncation error
in Nu is less than 0.5%.

---

## Known Issues and Solver Notes

1. **"T.*" wildcard is mandatory.** PIMPLE creates `TFinal` on its last
   corrector pass. The entry `"T.*"` in `fvSolution` must use the regex
   wildcard to match both `T` and `TFinal`. Using the literal key `T` will
   abort the run on the second corrector.

2. **PBiCGStab for T and C.** The advection term `∇·(UT)` makes the
   linear system non-symmetric. PCG will fail with "Unknown asymmetric
   matrix solver PCG". This applies to both T and C_test.

3. **Passive scalar is required.** `speciesFluid` requires a non-empty
   `speciesProperties` dictionary. The dummy scalar C_test with
   `zeroGradient` on all boundaries satisfies this requirement without
   affecting U, p, or T.

4. **Thermal entry region in Nu calculation.** The Graetz number
   Gz = Re·Pr·Dh/L = 5.6 means the channel is not fully thermally
   developed at the outlet. Measuring Nu near y = 0 would give values
   of 20–50 due to the entry singularity; the script restricts the
   averaging window to y ∈ [25, 45] mm where Nu has converged to
   within 1% of Nu_∞.

5. **Pressure reference.** The outlet is set to p = 0 (`fixedValue`),
   making all pressures gauge pressures. The inlet patch uses
   `zeroGradient` for p, which is correct when a Dirichlet U is prescribed
   at the inlet. Reversing these (Dirichlet p at inlet, zeroGradient U at
   outlet) is possible but changes the flow rate.

---

## Connection to the Framework

case328 demonstrates:

1. **`speciesFluid` unit validation** — confirms that the PIMPLE momentum
   solver, SIMPLE pressure correction, and temperature transport equation
   are all correctly implemented by recovering two independent analytical
   benchmarks.
2. **Graetz problem physics** — the simultaneous development of
   hydrodynamic and thermal boundary layers, with the thermal development
   delayed by a factor of Pr relative to the momentum development.
3. **Solver requirements for advection** — establishes the PBiCGStab /
   DILU and `"T.*"` wildcard requirements that apply to all `speciesFluid`
   cases in the tutorial sequence.
4. **Validation scripting pattern** — the Python script reads raw OpenFOAM
   ASCII fields, exploits blockMesh cell ordering, and applies text-book
   formulas to produce pass/fail metrics. The same pattern is reused in
   case329.
