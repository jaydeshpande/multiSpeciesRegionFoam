# case314 — Composite Membrane (Sieverts Interface)

**Solver:** `speciesSolid` in both regions  
**Coupling:** Two regions (PyC + SiC) with `sievertsCoupledMixed` interface  
**Trapping:** None (`noTrapping`)

---

## What this case covers

This is the first two-region case. A species diffuses from a high-pressure inlet (PyC side) through a composite PyC/SiC membrane to a vacuum outlet (SiC side). The two materials have different Sieverts solubilities, so the concentration is **discontinuous at the interface**: it drops by the ratio Ks_PyC / Ks_SiC = 2 when crossing from PyC into SiC.

The case validates:

- `splitMeshRegions -cellZones -overwrite` for region decomposition from a single `blockMesh`
- `sievertsCoupledMixed` boundary condition on both sides of the PyC/SiC interface
- Correct concentration jump and flux conservation at the interface
- Multi-region `foamMultiRun` with `speciesSolid` in every region

---

## Governing equation

In each region separately:

```
∂C/∂t  =  D · ∂²C/∂x²
```

At the interface x = L_PyC the conditions are:

```
flux continuity:     D_PyC · (∂C/∂x)|_PyC  =  D_SiC · (∂C/∂x)|_SiC
Sieverts partition:  C_PyC / Ks_PyC         =  C_SiC / Ks_SiC
```

---

## Geometry and parameters

```
x = 0          PyC inlet (high-pressure, C = 1 mol m⁻³)
x = 0.5 mm     PyC/SiC interface
x = 1.0 mm     SiC outlet (vacuum, C = 0)
```

| Parameter | PyC | SiC |
|-----------|-----|-----|
| Length | 0.5 mm | 0.5 mm |
| Cells | 100 | 100 |
| Diffusivity D | 1×10⁻⁹ m² s⁻¹ | 1×10⁻⁹ m² s⁻¹ |
| Solubility Ks | 2000 mol m⁻³ | 1000 mol m⁻³ |
| Partition ratio r = Ks_PyC / Ks_SiC | 2 | — |

---

## Analytical steady-state solution

At steady state the flux is uniform across both layers:

```
J  =  D · ΔC_total / (L_PyC / Ks_PyC + L_SiC / Ks_SiC) · (1/Ks_ref)
```

For equal thicknesses and equal D, with inlet C_in = 1 mol m⁻³ and outlet C_out = 0:

**Interface concentrations (Sieverts jump):**

```
C_SiC,interface  =  C_in / (1 + r)  =  1/3  mol m⁻³
C_PyC,interface  =  r · C_SiC,interface  =  2/3  mol m⁻³
```

**Steady-state profiles (piecewise linear):**

```
C_PyC(x)  =  1 − x / (3 · L_PyC)            x ∈ [0, L_PyC]
C_SiC(x)  =  (1/3) · (1 − (x − L_PyC)/L_SiC)    x ∈ [L_PyC, L]
```

**Steady-state permeation flux:**

```
J  =  D · C_SiC,interface / L_SiC  =  1×10⁻⁹ × (1/3) / 5×10⁻⁴  ≈  6.67×10⁻⁷  mol m⁻² s⁻¹
```

---

## How to run

```sh
cd tutorials/case314-composite-membrane
./Allrun
```

The `Allrun` script calls:

1. `blockMesh` — generates the two-block mesh with cell zones `pyc` and `sic`
2. `splitMeshRegions -cellZones -overwrite` — splits into regions `pyc` and `sic`; creates coupled patches named `pyc_to_sic` and `sic_to_pyc` automatically
3. `foamMultiRun` — runs both regions concurrently with `speciesSolid`

**Validate:**

```sh
python3 validation/sieverts_compare.py
```

The script reads the `C_T2` field in both regions at the last write time and checks the piecewise-linear profile and the interface jump against the analytical values at 2 % tolerance.
