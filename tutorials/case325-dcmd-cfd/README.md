# case325 — DCMD with Explicit Flow (CFD Temperature Polarisation)

**Solver:** `speciesSolid` (all three regions)  
**Geometry:** Three-region 2-D slab — feed channel | PTFE membrane | permeate channel  
**Coupling:** Thermal continuity at fluid-membrane interfaces (`coupledTemperature`)  
**Reference:** Khalifa A. et al., *Desalination* **404** (2017) 22–34  

---

## Purpose

This case is the companion to case324. Both model the same DCMD experiment
(Khalifa 2017), but use different strategies for temperature polarisation:

| Aspect | case324 | case325 |
|---|---|---|
| TPC method | Analytical (heat transfer correlations) | CFD (explicit plug flow) |
| Mesh | 1-D membrane slab (1 region) | 2-D three-region (feed|membrane|permeate) |
| T_mf, T_mp | Fixed: 328 K, 299 K | Emerge from coupled CFD |
| C_H2O BC | Fixed Antoine at T_mf, T_mp | Dynamic Antoine reading live T |
| Solver | speciesSolid | speciesSolid + thermal advection |

The key question: **how close is the analytical TPC to the CFD-derived TPC?**

---

## Physical model extensions

### Thermal advection in speciesSolid

`speciesSolid` is modified to add the convective heat transfer term when a
velocity field `U` is registered in the mesh region:

```
∂(ρe)/∂t  +  ∇·(ρU e)  +  ∇·(−κ∇T)  =  0
```

This term is implemented as:

```cpp
if (mesh_.foundObject<volVectorField>("U"))
{
    const surfaceScalarField rhoPhi(..., fvc::interpolate(rho)*fvc::flux(U));
    eEqn += fvm::div(rhoPhi, e);
}
```

The scheme `div(rhoPhi,e)` must be listed in `fvSchemes` for the fluid regions.

### Antoine equilibrium BC (`antoineEquilibrium`)

At each membrane face, `C_H2O` is set dynamically from the local temperature:

```
log₁₀(p_sat / mmHg) = 8.10765 − 1750.286 / (235 + T [°C])
C = p_sat [Pa] / (R · T [K])
```

The `antoineEquilibrium` BC is a `fixedValueFvPatchScalarField` registered in
`libspeciesCoupling.so`. It reads the face `T` field (supplied by
`coupledTemperature` from the adjacent fluid region) and evaluates the Antoine
equation each corrector. Antoine coefficients and the gas constant are
specified in the BC dictionary:

```
type            antoineEquilibrium;
A               8.10765;
B               1750.286;
C               235.0;
```

---

## Geometry and mesh

```
y = 66 mm ──────────────────────────────────────────  permeate_inlet
          ║  feed channel   ║ membrane ║ permeate  ║
          ║  (+y flow)      ║ (PTFE)   ║ (-y flow) ║
          ║  x: 0→5 mm      ║ 5→5.154  ║ 5.154→10.154 ║
y = 0 ───────────────────────────────────────────────  permeate_outlet
                  feed_inlet ↑            ↑ permeate_outlet
```

| Region | x extent | Cells (x) | Grading |
|---|---|---|---|
| feed | 0 → 5 mm | 20 | 0.05 (fine near membrane) |
| membrane | 5 → 5.154 mm | 10 | uniform (15.4 µm) |
| permeate | 5.154 → 10.154 mm | 20 | 20 (fine near membrane) |
| axial (y) | 0 → 66 mm | 100 | uniform (0.66 mm) |

Total cells: (20 + 10 + 20) × 100 = **5 000**

---

## Operating conditions

| Parameter | Value |
|---|---|
| T_feed inlet | 333 K (60°C) |
| T_perm inlet | 293 K (20°C) |
| U_feed | (0, +0.215, 0) m/s — plug flow |
| U_perm | (0, −0.169, 0) m/s — counter-flow |
| Re_feed | ≈ 3 750 |
| Re_perm | ≈ 1 393 |

---

## Coupling architecture

```
feed region                membrane region               permeate region
────────────               ───────────────               ──────────────
T: coupledTemperature ←→   T: coupledTemperature ←→     T: coupledTemperature
C_H2O: zeroGradient   ←→   C_H2O: antoineEquilibrium(T)    C_H2O: zeroGradient
U: prescribed 0.215 m/s    (no U field)                  U: prescribed -0.169 m/s
```

The coupling works as follows at each PIMPLE corrector:
1. Energy equations solved in all three regions simultaneously
2. `coupledTemperature` BCs enforce T continuity at both interfaces
3. `thermo_.correct()` updates T from e in all regions
4. `antoineEquilibrium` for C_H2O reads the updated T at the membrane face
5. Species equation solved in membrane with updated BC values
6. D(T) in membrane is updated from the new temperature field

---

## Simplifications and limitations

1. **Plug flow** — U is uniform across the channel cross-section. The actual
   velocity profile in a rectangular channel is parabolic, which would give
   a thicker thermal boundary layer near the membrane and lower TPC.

2. **No latent heat coupling** — Water evaporation at the hot membrane face
   removes ≈ 2.45 MJ/kg of latent heat that should cool T_mf. This is the
   dominant mechanism driving TPC in DCMD. Without it, the CFD underpredicts
   the thermal resistance and overestimates T_mf (TPC is underestimated).

3. **No concentration polarisation** — Water vapour concentration at the
   membrane faces is set by Antoine(T); the liquid-phase C_H2O field is a
   dummy with D → 0. For high-flux DCMD, the depletion of water at the feed
   surface could be significant.

These simplifications are intentional for tutorial clarity. The primary goal
of case325 is to demonstrate the **CFD approach** to computing TPC, not to
achieve higher accuracy than case324.

---

## How to run

```sh
cd tutorials/case325-dcmd-cfd
./Allrun
```

Steps:
1. `blockMesh` — creates a single mesh with 3 cell zones
2. `splitMeshRegions -cellZones -overwrite` — splits into feed/membrane/permeate regions
3. `foamMultiRun` — runs speciesSolid in all 3 regions, 2000 steps (0 to 2 s)

Estimated runtime: ~2–5 minutes for 2000 timesteps on a single core.

**Validate and compare with case324:**

```sh
python3 validation/dcmd_cfd_compare.py
```

The script reads the final membrane T and C_H2O fields, computes:
- T_mf and T_mp from the CFD membrane temperature field
- TPC_cfd = (T_mf − T_mp) / (T_feed − T_perm)
- J_cfd from membrane C_H2O gradient  
- Comparison table vs case324 analytical (TPC=0.725, J=31.1 L/(m²·h))

---

## Expected results

At t = 2 s, quasi-steady temperature polarisation has developed. The
plug-flow approximation without latent heat coupling gives a TPC estimate
that is somewhat different from case324's TPC = 0.725:

| Quantity | case324 | case325 (expected) |
|---|---|---|
| TPC | 0.725 | ~0.6–0.9 (plug flow) |
| T_mf [K] | 328 | ~320–335 |
| T_mp [K] | 299 | ~291–304 |
| J [L/(m²·h)] | 31.1 | ~25–38 |
| J_exp [L/(m²·h)] | 35 | — |

The CFD TPC differs from the analytical value because:
- Plug flow overestimates convection (thinner boundary layer than Poiseuille)
- Latent heat is absent, reducing the thermal load on the feed boundary layer

Adding latent heat as a heat flux BC at the membrane faces (q = J × L_vap)
would bring case325 into better agreement with both case324 and experiment.

---

## Connection to the framework

case325 demonstrates:
1. **Multi-region coupling** with `splitMeshRegions` for a membrane-in-channel geometry
2. **Thermal advection** in `speciesSolid` via the `rhoPhi` mass flux term
3. **Dynamic concentration BCs** using `antoineEquilibrium` reading the live T field
4. **Counter-flow configuration** with permeate entering at y = 66 mm
5. **Strategy comparison**: case324's analytical TPC vs case325's CFD TPC
