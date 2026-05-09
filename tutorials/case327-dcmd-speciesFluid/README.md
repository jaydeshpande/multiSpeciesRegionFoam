# case327 — DCMD with Solved Incompressible Flow (speciesFluid)

**Solver:** `speciesFluid` (feed + permeate), `speciesSolid` (membrane)
**Geometry:** Three-region 2-D slab — feed channel | PTFE membrane | permeate channel
**Coupling:** `latentHeatFlux` + `antoineEquilibrium` at fluid-membrane interfaces
**Reference:** Khalifa A. et al., *Desalination* **404** (2017) 22–34

---

## Purpose

case327 is the next step in the DCMD tutorial sequence. The table below
shows how the three CFD cases differ:

| Aspect | case325 | case326 | case327 |
|---|---|---|---|
| Feed/permeate solver | `speciesSolid` | `speciesSolid` | **`speciesFluid`** |
| Velocity | Prescribed plug flow | Prescribed plug flow | **PIMPLE solved** |
| Latent heat at membrane | Absent | `latentHeatFlux` | `latentHeatFlux` |
| Velocity profile | Uniform | Uniform | **Develops to Poiseuille** |
| Thermal entry effects | None (plug flow) | None (plug flow) | **Graetz-type entrance** |
| Temperature polarisation | CFD, no latent heat | CFD + latent heat | CFD + latent heat + Poiseuille |

The scientific question is: **how much does the parabolic velocity profile
change the thermal boundary layer thickness and therefore the temperature
polarisation coefficient (TPC) and vapour flux J, compared with the
plug-flow cases 325 and 326?**

The Poiseuille profile produces a thicker thermal boundary layer near the
membrane than plug flow (the velocity goes to zero at the wall), which
increases the thermal resistance, lowers T_mf, and reduces J. case327
makes this effect visible without any analytical approximation.

---

## Physical Model

### Region solvers

`foamMultiRun` drives three simultaneous solvers:

```
regionSolvers
{
    feed        speciesFluid;   // PIMPLE U/p/T/C
    membrane    speciesSolid;   // conduction + diffusion only
    permeate    speciesFluid;   // PIMPLE U/p/T/C
}
```

### Governing equations — fluid regions (speciesFluid)

**Continuity** (incompressible):

```
∇·U = 0
```

**Momentum** (PIMPLE, laminar):

```
∂U/∂t + (U·∇)U = −∇(p/ρ) + ν∇²U
```

**Energy** (temperature form, speciesFluid solves T directly):

```
∂T/∂t + ∇·(UT) = ∇·(α∇T),    α = κ/(ρCp)
```

**Species transport** (C_H2O, vapour concentration equivalent):

```
∂C/∂t + ∇·(UC) = ∇·(D∇C)
```

In the fluid channels D = 1×10⁻¹⁰ m²/s (negligible; C_H2O is a dummy
scalar, as vapour transport occurs only in the membrane pores).

### Governing equations — membrane region (speciesSolid)

**Energy** (no flow; solid conduction):

```
ρCv ∂T/∂t = ∇·(κ∇T)
```

**Species** (Fickian diffusion with Arrhenius D(T)):

```
∂C/∂t = ∇·(D(T)∇C),    D(T) = 6.9×10⁻⁵ exp(−3700 / (R·T))  m²/s
```

### Antoine equilibrium BC (C_H2O at membrane faces)

At each membrane face the `antoineEquilibrium` boundary condition
evaluates the equilibrium vapour concentration from the local face
temperature:

```
log₁₀(p_sat / mmHg) = A − B / (C_ant + T [°C])
    A = 8.10765,  B = 1750.286,  C_ant = 235.0

C_face = p_sat [Pa] / (R · T [K])
```

The membrane face temperature T is supplied by the `latentHeatFlux` BC on
the adjacent fluid patch (via the coupled interface machinery), so the
Antoine equilibrium tracks the evolving CFD temperature field.

### Latent heat flux BC (T at membrane faces)

At each fluid-membrane interface the `latentHeatFlux` BC enforces thermal
continuity (same as `coupledTemperature`) and additionally applies the
evaporative/condensation heat flux:

```
q_latent = J · L_vap · M_H2O    [W/m²]
J = −D(T) · ∂C/∂n |_membrane    [mol/(m²·s)]
```

The interface heat balance including the latent term:

```
κ_fluid · ∂T/∂n = κ_mem · ∂T/∂n  ± q_latent
```

Sign convention: subtracted on the feed (evaporation cools T_mf),
added on the permeate (condensation warms T_mp).

---

## Geometry and Mesh

```
y = 66 mm  ─────────────────────────────────────  permeate_inlet (↓ counter-flow)
           │  feed channel  │ PTFE mem │ permeate │
           │  (+y upward)   │          │ (-y down) │
           │  x: 0→5 mm     │ 5→5.154  │ 5.154→10.154 mm │
y = 0 mm   ─────────────────────────────────────  feed_inlet (↑), permeate_outlet
```

| Region | x extent | Cells (x) | Grading (x) | Cell size at interface |
|---|---|---|---|---|
| feed | 0 → 5 mm | 20 | 0.05 (fine near membrane) | ~23 µm |
| membrane | 5 → 5.154 mm | 10 | uniform | 15.4 µm |
| permeate | 5.154 → 10.154 mm | 20 | 20 (fine near membrane) | ~23 µm |
| axial (all) | 0 → 66 mm | 100 | uniform | 0.66 mm |

Total cells: (20 + 10 + 20) × 100 = **5 000**

The mesh is identical to cases 325 and 326. Only the region solvers
change: `speciesFluid` replaces `speciesSolid` for feed and permeate.

---

## Material Properties

### Feed channel — hot seawater at 333 K (60°C)

| Property | Symbol | Value | Units |
|---|---|---|---|
| Kinematic viscosity | ν | 4.74×10⁻⁷ | m²/s |
| Density | ρ | 983 | kg/m³ |
| Specific heat | Cp | 4183 | J/(kg·K) |
| Thermal conductivity | κ | 0.654 | W/(m·K) |
| Thermal diffusivity | α = κ/(ρCp) | 1.59×10⁻⁷ | m²/s |

### Permeate channel — cold freshwater at 293 K (20°C)

| Property | Symbol | Value | Units |
|---|---|---|---|
| Kinematic viscosity | ν | 1.004×10⁻⁶ | m²/s |
| Density | ρ | 998 | kg/m³ |
| Specific heat | Cp | 4182 | J/(kg·K) |
| Thermal conductivity | κ | 0.598 | W/(m·K) |
| Thermal diffusivity | α = κ/(ρCp) | 1.43×10⁻⁷ | m²/s |

### PTFE membrane — effective porous properties

| Property | Symbol | Value | Units |
|---|---|---|---|
| Effective density | ρ_eff = (1−ε)×ρ_PTFE | 440 | kg/m³ |
| Specific heat | Cv | 1000 | J/(kg·K) |
| Effective thermal conductivity | κ_eff | 0.07 | W/(m·K) |
| Membrane thickness | d | 0.154 | mm |
| Pore diameter | d_p | 0.45 | µm |
| Porosity | ε | 0.80 | — |

The effective vapour diffusivity in the membrane pores combines the
Knudsen and molecular diffusivity in series, fitted to Arrhenius form:

```
D_eff(T) = 6.9×10⁻⁵ exp(−3700 / (8.314 · T))   m²/s
```

Values: D(299 K) = 1.57×10⁻⁵ m²/s, D(328 K) = 1.79×10⁻⁵ m²/s.

---

## Key Dimensionless Numbers

| Parameter | Feed | Permeate |
|---|---|---|
| Mean velocity U | 0.215 m/s | 0.169 m/s |
| Hydraulic diameter D_h | 8.28 mm (2×W+2×gap/2) | 8.28 mm |
| Reynolds number Re = U·D_h/ν | ≈ 3 750 | ≈ 1 393 |
| Prandtl number Pr = ν/α | 2.98 | 7.01 |
| Thermal entry length L_T = 0.04·Re·Pr·D_h | ≈ 3.7 m | ≈ 3.2 m |
| Peclet number Pe_T = U·W/α | 6 761 | 5 929 |
| Momentum turbulence threshold Re_t | 3 750 | 3 750 |

**Note on Re_feed:** The feed Re ≈ 3 750 equals the transition threshold
for rectangular channels. The `momentumTransport` model is set to
`laminar` in both fluid regions (a `RASModel k-ε` would be needed for
turbulent cases). This is physically borderline; the plug-flow cases
(325, 326) are less sensitive to this choice because the velocity field is
prescribed. case327 uses laminar throughout.

**Thermal entry length:** L_T >> 66 mm for both channels, so the flow is
thermally developing along the entire channel length. The temperature
boundary layer grows from the inlet (and from the membrane face) and does
not reach the fully developed Nusselt limit within the channel. This is
the main physical difference from a plug-flow treatment.

---

## Boundary Conditions

### Feed region

| Patch | U | p [m²/s²] | T [K] | C_H2O |
|---|---|---|---|---|
| `feed_inlet` (y=0) | `fixedValue (0 0.215 0)` | `zeroGradient` | `fixedValue 333` | `zeroGradient` |
| `feed_outlet` (y=66mm) | `zeroGradient` | `fixedValue 0` | `zeroGradient` | `zeroGradient` |
| `feed_outer` (x=0) | `noSlip` | `zeroGradient` | `zeroGradient` | `zeroGradient` |
| `feed_to_membrane` (x=5mm) | `noSlip` | `zeroGradient` | `latentHeatFlux` | `zeroGradient` |

The C_H2O field in the feed is a dummy scalar (D = 1×10⁻¹⁰ m²/s) that
carries no physical meaning. Vapour concentration is meaningful only in
the membrane pores; the `zeroGradient` BC keeps it at zero.

### Membrane region

| Patch | T [K] | C_H2O [mol/m³] |
|---|---|---|
| `membrane_to_feed` (x=5mm) | `latentHeatFlux` (evaporation, L_vap=2.45×10⁶ J/kg, M=0.018 kg/mol) | `antoineEquilibrium` (A=8.10765, B=1750.286, C=235.0) |
| `membrane_to_permeate` (x=5.154mm) | `latentHeatFlux` (condensation) | `antoineEquilibrium` |
| `membrane_inlet`, `membrane_outlet` | `zeroGradient` | `zeroGradient` |

### Permeate region

| Patch | U | p [m²/s²] | T [K] | C_H2O |
|---|---|---|---|---|
| `permeate_inlet` (y=66mm) | `fixedValue (0 -0.169 0)` | `zeroGradient` | `fixedValue 293` | `zeroGradient` |
| `permeate_outlet` (y=0) | `zeroGradient` | `fixedValue 0` | `zeroGradient` | `zeroGradient` |
| `permeate_outer` (x=10.154mm) | `noSlip` | `zeroGradient` | `zeroGradient` | `zeroGradient` |
| `permeate_to_membrane` (x=5.154mm) | `noSlip` | `zeroGradient` | `latentHeatFlux` | `zeroGradient` |

---

## physicalProperties Format for speciesFluid

`speciesFluid` is built on OpenFOAM-13's `incompressibleFluid` base class,
which reads the viscosity model, and adds its own reads of `rho`, `Cp`,
and `mixture.transport.kappa`. All four entries are required:

```
viscosityModel  Newtonian;
nu              [0 2 -1 0 0 0 0] 4.74e-7;

rho             983;
Cp              4183;

mixture
{
    transport
    {
        kappa   0.654;
    }
}
```

Omitting any of these entries will produce a `not found in dictionary`
error at startup.

---

## fvSolution Requirements

### "T.*" wildcard

PIMPLE's final corrector pass creates a `TFinal` field. The solver
dictionary must use the wildcard pattern `"T.*"` to match both `T` and
`TFinal`. Using the bare key `T` leaves `TFinal` without a solver entry
and causes a fatal error on the second corrector.

### PBiCGStab for T and C

The advection term `∇·(UT)` makes the matrix of the temperature and
species equations non-symmetric. PCG (which requires a symmetric positive
definite matrix) will crash with:

```
Unknown asymmetric matrix solver PCG
```

Use `PBiCGStab` with `DILU` preconditioner for all advection-dominated
scalar equations:

```
"T.*"
{
    solver          PBiCGStab;
    preconditioner  DILU;
    tolerance       1e-10;
    relTol          0;
}

"C_.*"
{
    solver          PBiCGStab;
    preconditioner  DILU;
    tolerance       1e-10;
    relTol          0;
}
```

The pressure equation `p` remains symmetric (∇·(1/ρ ∇p) = ∇·U); PCG
with DIC is correct for p.

---

## fvSchemes

Fluid regions require div-scheme entries for each advection term:

```
divSchemes
{
    default             none;
    div(phi,U)          Gauss linearUpwind grad(U);
    div(phi,T)          Gauss linearUpwind grad(T);
    div(phi,C_H2O)      Gauss upwind;
}
```

`linearUpwind` is second-order bounded and reduces numerical diffusion for
T at Pe_T ≈ 6 761. `upwind` is mandatory for C_H2O because the species
Peclet number in the membrane is effectively infinite (the species equation
in the fluid regions has D → 0).

The membrane region has no `U` field, so `divSchemes` contains only
`default none`.

---

## Running the Case

Build the library first if not already done:

```sh
cd /path/to/multiSpeciesRegionFoam
./Allwmake
```

Then run the case:

```sh
cd tutorials/case327-dcmd-speciesFluid
./Allrun
```

The `Allrun` script executes:

1. `blockMesh` — creates the single three-zone mesh (5 000 cells)
2. `splitMeshRegions -cellZones -overwrite` — splits into feed, membrane,
   permeate regions and creates coupled patches
3. `foamMultiRun` — runs for 10 s (10 000 steps of Δt = 0.001 s)

Estimated runtime: 15–30 minutes on a single core (PIMPLE loop with
solved pressure equation is more expensive than the prescribed-U cases).

**Check convergence:** The momentum residuals should drop to < 10⁻⁶
within the first 2 s. The temperature field reaches a quasi-steady
distribution by t ≈ 5 s. The vapour flux J (from the membrane C_H2O
gradient) stabilises last, typically by t ≈ 8–10 s.

---

## Expected Results

At t = 10 s the solution is quasi-steady. The Poiseuille velocity profile
is fully developed within the first few millimetres of the channel;
thermal and concentration boundary layers grow from the inlet and from the
membrane face simultaneously.

### Comparison with cases 325 and 326

| Quantity | case325 (plug, no L_vap) | case326 (plug + L_vap) | case327 (Poiseuille + L_vap) |
|---|---|---|---|
| TPC | ~0.6–0.9 | reduced vs 325 | further reduced vs 326 |
| T_mf [K] | ~320–335 | lower than 325 | lower than 326 |
| T_mp [K] | ~291–304 | higher than 325 | higher than 326 |
| J [L/(m²·h)] | ~25–38 | lower than 325 | lower than 326 |
| J_exp [L/(m²·h)] | 35 | — | — |

The Poiseuille profile produces a thicker thermal boundary layer at the
membrane (the no-slip condition forces U → 0 at the wall), increasing the
thermal resistance and lowering T_mf. This effect is additive with the
latent heat cooling introduced in case326.

### Self-consistency check

The primary validation is internal consistency of the vapour flux:

```
J_feed   = ∫ (−D · ∂C/∂x)|_{x = 5 mm}   dy / L     [mol/(m²·s)]
J_perm   = ∫ (−D · ∂C/∂x)|_{x = 5.154 mm} dy / L
```

J_feed should equal J_perm to within the solver tolerance (≈ 0.1%).
A discrepancy larger than 1% indicates insufficient convergence or
incorrect BC sign convention in `latentHeatFlux`.

---

## Known Issues and Solver Notes

1. **"T.*" wildcard is mandatory.** PIMPLE creates `TFinal` on its last
   corrector pass. Without the wildcard, `TFinal` has no solver entry and
   the run aborts.

2. **PBiCGStab for T and C.** The advection operator breaks the matrix
   symmetry. PCG will fail with "Unknown asymmetric matrix solver PCG".
   Use PBiCGStab/DILU for all scalar transport equations in the fluid
   regions.

3. **Membrane region uses PCG/DIC.** The membrane has no velocity field;
   its T (solved as internal energy `e`) and C_H2O equations are pure
   Laplacians (symmetric). PCG with DIC is correct here.

4. **physicalProperties dual format.** The `viscosityModel`/`nu` entries
   are required by the `incompressibleFluid` base class. The `rho`, `Cp`,
   and `mixture.transport.kappa` entries are read by `speciesFluid`. Both
   sets are required simultaneously in the same `physicalProperties`
   dictionary.

5. **latentHeatFlux source = self (membrane).** The membrane T BC uses
   `source self`, meaning it reads the C_H2O gradient from its own patch
   rather than from the adjacent fluid. This avoids requiring a coupled
   field lookup across two boundary layers simultaneously.

6. **Transitional Re in feed.** The feed Re ≈ 3 750 is near the
   transition threshold. The `laminar` momentumTransport model is used for
   tutorial simplicity. For a higher-fidelity simulation, a RANS model
   (k-ε or k-ω SST) would be appropriate.

---

## Connection to the Framework

case327 demonstrates:

1. **`speciesFluid` in multi-region coupling** — the PIMPLE loop runs
   simultaneously in both fluid regions, exchanging temperature at the
   membrane interface through `latentHeatFlux` BCs.
2. **Parabolic velocity development** — the Poiseuille profile emerges
   naturally from the no-slip BC without any analytical prescription.
3. **Progressive model complexity** — case325 (plug flow, no latent heat)
   → case326 (plug flow + latent heat) → case327 (Poiseuille + latent
   heat).
4. **Asymmetric solver requirement** — advection in the energy and species
   equations mandates PBiCGStab, not PCG. This constraint propagates to
   any `speciesFluid` case and is the most common failure mode for new
   users.
