# case317 — Reverse Osmosis (Solution-Diffusion)

**Solver:** `speciesSolid` in both regions  
**Coupling:** Two regions (feed + membrane) with `sievertsCoupledMixed` interface  
**Trapping:** None (`noTrapping`)

---

## What this case covers

This case models the solution-diffusion mechanism that governs transport through dense polymer membranes in reverse osmosis (RO). Water partitions from a liquid feed into a dense polyamide membrane according to Henry's law, diffuses down the concentration gradient, and exits on the permeate side.

The thermodynamic partition at the interface is represented by the ratio of Sieverts solubility constants:

```
K_H  =  C_membrane / C_feed  =  Ks_membrane / Ks_feed  =  1000 / 2000  =  0.5
```

A concentration boundary layer (concentration polarisation) on the feed side is resolved explicitly as a separate `feed` region, so the full profile from bulk feed through to the permeate side is captured.

The case validates:

- Two-region coupling where one region is a liquid-phase boundary layer and the other is a dense polymer membrane
- Henry-law partition implemented via `sievertsCoupledMixed` (the ratio Ks_m/Ks_f sets K_H)
- Piecewise-linear steady-state profiles and flux conservation across the interface
- Correct interface concentration jump C_f,if = 2/3 → C_m,if = 1/3 mol m⁻³

---

## Governing equation

In each region (pure 1-D diffusion, no advection):

```
∂C/∂t  =  D · ∂²C/∂x²
```

At the feed/membrane interface (x = L_feed):

```
flux continuity:    D_feed · (∂C/∂x)|_feed  =  D_membrane · (∂C/∂x)|_membrane
Henry partition:    C_membrane / Ks_m        =  C_feed / Ks_f
                 →  C_membrane              =  K_H · C_feed     (K_H = 0.5)
```

---

## Geometry and parameters

```
x = 0            Feed inlet (bulk concentration, C = 1 mol m⁻³)
x = 0.5 mm       Feed/membrane interface (K_H = 0.5 concentration jump)
x = 1.0 mm       Membrane outlet / permeate (C = 0)
```

| Parameter | Feed region | Membrane region |
|-----------|------------|----------------|
| Thickness | 0.5 mm | 0.5 mm |
| Cells | 100 | 100 |
| Diffusivity D | 1×10⁻⁹ m² s⁻¹ | 1×10⁻⁹ m² s⁻¹ |
| Sieverts Ks (Ks_f, Ks_m) | 2000 mol m⁻³ | 1000 mol m⁻³ |
| Partition ratio K_H = Ks_m / Ks_f | — | 0.5 |

The feed diffusivity (1×10⁻⁹ m² s⁻¹) is characteristic of dissolved-solute diffusion in liquid water. The membrane diffusivity (equal here) represents water in a dense polyamide film; in practice it can differ by orders of magnitude.

---

## Analytical steady-state solution

At steady state the flux J is uniform across both regions. Combining flux continuity with the Henry partition:

**Interface concentrations:**

```
C_feed,interface      =  2/3  mol m⁻³
C_membrane,interface  =  K_H · C_feed,interface  =  1/3  mol m⁻³
```

**Piecewise-linear steady-state profiles:**

```
Feed region     (x ∈ [0, 0.5 mm]):   C_f(x)  =  1 − x / 1.5×10⁻³
Membrane region (x ∈ [0.5, 1.0 mm]): C_m(x)  =  1/3 − (x − 5×10⁻⁴) / 1.5×10⁻³
```

**Steady-state permeation flux:**

```
J  =  D · (C_feed,if − C_permeate) / L_membrane  =  1×10⁻⁹ × (1/3) / 5×10⁻⁴  ≈  6.67×10⁻⁷  mol m⁻² s⁻¹
```

This equals the feed-side flux:

```
J  =  D · (C_in − C_feed,if) / L_feed  =  1×10⁻⁹ × (1/3) / 5×10⁻⁴  ≈  6.67×10⁻⁷  mol m⁻² s⁻¹  ✓
```

---

## Interface boundary condition setup

Both sides of the feed/membrane interface use `sievertsCoupledMixed`. The key point is that `Ks` refers to the **self** (this-region) solubility and `KsNbr` to the **neighbour** solubility. The roles are **swapped** between the two patches:

**`0/feed/C_H2O` at `feed_to_membrane`:**

```
Ks    { X0 [0 -3 0 0 1 0 0] 2000; }   // feed solubility
KsNbr { X0 [0 -3 0 0 1 0 0] 1000; }   // membrane solubility
```

**`0/membrane/C_H2O` at `membrane_to_feed`:**

```
Ks    { X0 [0 -3 0 0 1 0 0] 1000; }   // membrane solubility
KsNbr { X0 [0 -3 0 0 1 0 0] 2000; }   // feed solubility
```

The coupled-patch names `feed_to_membrane` and `membrane_to_feed` are created automatically by `splitMeshRegions -cellZones -overwrite` from the two `blockMesh` cell zones.

---

## How to run

```sh
cd tutorials/case317-reverse-osmosis
./Allrun
```

The `Allrun` script executes:

1. `blockMesh` — generates the two-block mesh with cell zones `feed` and `membrane`
2. `splitMeshRegions -cellZones -overwrite` — splits into regions `feed` and `membrane` and creates the coupled interface patches
3. `foamMultiRun` — runs both regions with `speciesSolid`

**Validate:**

```sh
python3 validation/ro_compare.py
```

The script reads `C_H2O` in both regions at the last write time and:

- Checks the feed-side profile against `C_f(x) = 1 − x/1.5×10⁻³` (2 % tolerance)
- Checks the membrane-side profile against `C_m(x) = 1/3 − (x−5×10⁻⁴)/1.5×10⁻³` (2 % tolerance)
- Reports the interface concentrations on both sides and compares J to 6.67×10⁻⁷ mol m⁻² s⁻¹ (2 % tolerance)
