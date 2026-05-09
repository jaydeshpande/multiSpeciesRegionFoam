# case324 — DCMD Validation against Khalifa et al. (2017)

**Solver:** `speciesSolid`  
**Coupling:** Single region, full energy + species coupling  
**Trapping:** None (`noTrapping`)  
**Reference:** Khalifa A. et al., *Desalination* **404** (2017) 22–34

---

## Purpose

This case validates multiSpeciesRegionFoam against experimental permeate flux
data for Direct Contact Membrane Distillation (DCMD). The membrane is a
hydrophobic PTFE 0.45 µm flat-sheet, and the permeate flux is driven by the
vapour-pressure difference across the membrane — set by the temperature
difference between the hot salt feed and the cold distillate.

---

## Experimental setup (Khalifa et al. 2017)

| Parameter | Value |
|---|---|
| Module dimensions | 160 × 160 × 25 mm |
| Number of channels per side | 3 |
| Channel dimensions | 66 × 24 × 5 mm (L × W × H) |
| Channel hydraulic diameter D_h | 8.28 mm |
| Membrane area (active) | 3 × 66 × 24 mm = 47.5 cm² |

### Membrane properties — PTFE 0.45 µm (Table 1)

| Property | Value |
|---|---|
| Total thickness δ | 154 µm |
| Active layer thickness δ_a | 7 µm |
| Pore size d_pore | 379 nm |
| Porosity ε | 80% |
| Contact angle θ | 139° |
| Liquid entry pressure (LEP) | 2.4 bar |

### Baseline operating condition (Fig. 4, highest flow rate)

| Parameter | Value |
|---|---|
| Feed temperature T_feed | 60°C |
| Permeate temperature T_perm | 20°C |
| Feed flow rate Q_feed | 4.65 L/min (total, 3 channels) |
| Permeate flow rate Q_perm | 3.65 L/min (total, 3 channels) |
| Feed concentration | 2 g/L NaCl |
| Measured permeate flux J_exp | ~35 L/(m²·h) |

---

## Physical model

### Mass transport through the membrane

Water vapour diffuses through the hydrophobic PTFE membrane under the
concentration (vapour-pressure) gradient. In the transition regime between
Knudsen and molecular diffusion, the combined pore diffusivity is:

```
1/D_pore(T) = 1/D_K(T) + 1/D_M(T)
```

**Knudsen diffusivity** (d_pore = 379 nm, r = 189.5 nm):
```
D_K(T) = (2r/3) · sqrt(8RT / (π·M_w))
D_K(314 K) ≈ 7.7 × 10⁻⁵ m²/s
```

**Molecular diffusivity** of H₂O vapour in air (Fuller correlation):
```
D_M(T) = 2.26 × 10⁻⁵ · (T/273.15)^1.81   [m²/s, 1 atm]
D_M(314 K) ≈ 2.9 × 10⁻⁵ m²/s
```

Since D_K >> D_M, the transport is dominated by molecular diffusion with a
Knudsen correction. The combined pore diffusivity D_pore ≈ 2.1 × 10⁻⁵ m²/s.

**Effective diffusivity** (porosity ε = 0.80, no added tortuosity for PTFE):
```
D_eff(T) = ε · D_pore(T) ≈ 1.68 × 10⁻⁵ m²/s  at 314 K
```

### Arrhenius fit to D_eff(T)

The OpenFOAM solver uses `D(T) = X0 · exp(-Ea / (R·T))`. Fitting D_eff at the
two membrane face temperatures (299 K and 328 K) gives:

| Parameter | Value |
|---|---|
| X0 | 6.9 × 10⁻⁵ m²/s |
| Ea | 3700 J/mol |
| D_eff(299 K) = 1.57e-5 m²/s | cold face (26°C) |
| D_eff(313.5 K) = 1.67e-5 m²/s | midpoint |
| D_eff(328 K) = 1.79e-5 m²/s | hot face (55°C) |
| D variation hot/cold | ≈ +14% |

Ea > 0 means D increases with temperature, consistent with the T-dependence of
both D_K(T) ∝ T^0.5 and D_M(T) ∝ T^1.81.

### Boundary concentrations (vapour pressure driving force)

Water vapour concentration at each membrane face is set from the Antoine
equation at the local membrane surface temperature:

```
log₁₀(p_sat / mmHg) = 8.10765 − 1750.286 / (235 + T [°C])
C = p_sat(T) / (R · T)
```

### Temperature polarisation

The membrane face temperatures differ from the bulk fluid temperatures due to
thermal resistance in the thin boundary layers adjacent to the membrane.

**Channel geometry** (3 channels, 66 × 24 × 5 mm):
- Cross-sectional area: A = 24 × 5 = 120 mm²
- Hydraulic diameter: D_h = 4A/P = 4 × 120 / (2 × 29) = 8.28 mm

**Flow velocities** (baseline, 3 channels):
- Feed: U_f = (4.65 L/min / 3) / 120 mm² = 0.215 m/s
- Permeate: U_p = (3.65 L/min / 3) / 120 mm² = 0.169 m/s

**Reynolds numbers** (water at T_f = 60°C, T_p = 20°C):
- Re_feed ≈ 3750  (developing laminar, approaching transition)
- Re_perm ≈ 1390  (developing laminar)

**Nusselt and heat transfer coefficients** (Sieder–Tate developing flow):
```
Nu = 0.664 · Re^0.5 · Pr^(1/3)
```

| Side | Re | Nu | h [W/(m²·K)] |
|---|---|---|---|
| Feed (60°C) | 3750 | 58 | ≈ 4600 |
| Permeate (20°C) | 1390 | 48 | ≈ 3400 |
| Membrane effective | — | — | κ/δ = 0.07/154e-6 ≈ 455 |

**Temperature polarisation coefficient** (solved iteratively including latent heat):
```
TPC = (T_mf − T_mp) / (T_feed − T_perm) ≈ 0.73
T_mf ≈ 55°C = 328 K   (hot face)
T_mp ≈ 26°C = 299 K   (cold face)
```

**Resulting boundary concentrations:**

| Face | T [K] | p_sat [Pa] | C [mol/m³] |
|---|---|---|---|
| hot_face | 328 (55°C) | 15 756 | **5.78** |
| cold_face | 299 (26°C) | 3 363 | **1.353** |

---

## Governing equations in the membrane

**Energy (Fourier conduction in PTFE solid):**
```
ρ_eff · Cv · ∂T/∂t  =  κ_eff · ∂²T/∂x²
```

**Species diffusion (water vapour) with temperature-dependent D:**
```
∂C/∂t  =  ∂/∂x [ D(T(x)) · ∂C/∂x ]
```

At steady state, T(x) is linear and the concentration profile follows the
same non-linear form as case316 (constant flux J = −D(T)·dC/dx):

```
C(x) = C_hot − J · ∫₀ˣ ds / D(T_ss(s))

where  J = (C_hot − C_cold) / ∫₀^δ ds / D(T_ss(s))
```

---

## Effective membrane thermal properties

PTFE membrane with 80% air-filled pores:

| Property | Formula | Value |
|---|---|---|
| Effective density ρ_eff | (1−ε)·ρ_PTFE | 440 kg/m³ |
| Heat capacity Cv | Cv_PTFE | 1000 J/(kg·K) |
| Effective conductivity κ_eff | ε·κ_air + (1−ε)·κ_PTFE | 0.07 W/(m·K) |
| Thermal diffusivity α | κ/(ρ·Cv) | 1.59 × 10⁻⁷ m²/s |

**Time constants** (δ = 154 µm):

| Process | Formula | Value |
|---|---|---|
| Thermal τ_T | δ²/(π²·α) | ≈ 0.015 s |
| Species τ_C | δ²/D̄ | ≈ 1.4 × 10⁻³ s |
| Run time | 33 × τ_T | 0.5 s |

---

## Predicted permeate flux vs. Khalifa et al. Fig. 4

Theoretical predictions using the physics model with TPC computed iteratively
(T_perm = 20°C, Q_f = 4.65 L/min, Q_p = 3.65 L/min):

| T_feed [°C] | J_theory [L/(m²·h)] | J_exp [L/(m²·h)] | ratio |
|---|---|---|---|
| 40 | ~9 | ~7 | 1.3 |
| 50 | ~20 | ~18 | 1.1 |
| 60 | **~31** | **~35** | **0.9** |
| 70 | ~52 | ~56 | 0.9 |
| 80 | ~75 | ~75 | 1.0 |
| 90 | ~103 | ~100 | 1.0 |

Agreement is within ±30% across the full range, with best agreement at
T_feed ≥ 70°C where temperature polarisation is weaker relative to the driving force.

---

## How to run

```sh
cd tutorials/case324-dcmd-khalifa2017
./Allrun
```

The `Allrun` script calls `blockMesh` then `foamMultiRun`. The energy and
species equations are coupled within each PIMPLE corrector; D(T) is evaluated
from the current temperature after every corrector.

**Validate:**

```sh
python3 validation/dcmd_khalifa2017_compare.py
```

The script:
1. Prints membrane physics model summary and D_eff(T) at 3 temperatures
2. Verifies the Arrhenius fit error vs. the physics-based model
3. Reports the temperature polarisation coefficients
4. Computes the steady-state flux J from the 1D analytical integral
5. Reads the final C_H2O and T fields and checks against the analytical profile
6. Computes the simulated permeate flux from the cold-face concentration gradient
7. Compares J_sim against J_exp ≈ 35 L/(m²·h) from Khalifa Fig. 4
8. Prints the full J vs. T_feed prediction table for comparison with Fig. 4
9. Reports PASS if max|C−C_ss|/ΔC < 3% and |J_sim−J_ss|/J_ss < 1%

---

## Mesh

| Parameter | Value |
|---|---|
| Domain | 1D slab, 154 µm × 1 mm × 1 mm |
| Cells | 100 (uniform, dx = 1.54 µm) |
| hot_face (x = 0) | Feed-side membrane surface |
| cold_face (x = δ) | Permeate-side membrane surface |
| sides | `empty` (2D → 1D effective) |

---

## Connection to multiSpeciesRegionFoam framework

This case demonstrates how to set up a DCMD membrane transport problem:

1. **Concentration BC from vapour pressure**: Use `fixedValue` with C computed
   from the Antoine equation at the membrane face temperature. No Sieverts law
   is needed — the PTFE/gas interface is not a partitioning boundary.

2. **Knudsen + molecular diffusivity**: Represented as Arrhenius D(T) with
   X0 and Ea fitted to the physics model. The solver evaluates D from the
   live temperature field at each PIMPLE corrector.

3. **Temperature polarisation pre-processing**: Membrane face temperatures
   (T_mf, T_mp) are pre-computed from heat transfer correlations and used as
   boundary conditions for the thermal solve. This is standard practice for
   DCMD resistance-network models and avoids the need to mesh the fluid channels.

4. **Extension to full channel CFD**: For a multi-region DCMD case coupling
   the feed, membrane, and permeate flow fields, add fluid regions with
   `speciesFluid` solvers, coupled at the membrane interfaces via
   `sievertsCoupledMixed` (setting Ks_fluid = Ks_membrane) or a custom
   vapour-pressure BC.
