# case321 — Henry's Law vs Sieverts' Law Partition

**Solver:** `foamMultiRun` with `speciesSolid` region module

---

## What this case covers

At the interface between two materials, the dissolved-species concentration is **not** necessarily continuous. The equilibrium partition depends on how each material absorbs the species:

| Regime | Law | Relation to partial pressure |
|--------|-----|------------------------------|
| Metal/solid | Sieverts' law | C = Kₛ · √p |
| Polymer / liquid | Henry's law | C = K_H · p |

When two Sieverts materials meet, the interface condition is **linear** (C₁/Kₛ₁ = C₂/Kₛ₂). When a Henry material meets a Sieverts material, the condition is a **square-root relation** (C_sieverts = Kₛ · √(C_henry / K_H)).

This case demonstrates:

1. How the `partition sqrt` BC correctly models a Henry-polymer → Sieverts-metal interface.
2. That the partition law changes the **shape** of the concentration profile and the **magnitude** of the steady-state flux.
3. How the complementary `partition quadratic` BC on the Henry side encodes the same physics from the other direction.

---

## Geometry

One-dimensional two-layer membrane (x-direction):

```
C = 1.0 mol/m³                                  C = 0
  │ ← polymer (0.5 mm, D=1e-9) → │ ← metal (0.5 mm, D=1e-9) → │
  x=0                           x=0.5 mm                       x=1 mm
```

Both layers have equal thickness and equal diffusivity so that the interface partition is the **only** physics that differs between the two sub-cases.

---

## Sub-cases

### sieverts_sieverts (reference)

Both polymer and metal follow Sieverts' law with equal constants Kₛ = 1.0. The interface condition (`partition linear`, Kₛ/Kₛ_nbr = 1) gives **concentration continuity**.

At steady state the profile is linear in both layers with a midpoint concentration of exactly 0.5 mol/m³.

```
J_SS = D · C₀ / (2L) = 1e-9 × 1.0 / 1e-3 = 1.0×10⁻⁶ mol/(m²·s)
```

### henry_sieverts

Polymer follows Henry's law (C = K_H · p, K_H = 1.0) and metal follows Sieverts' law (C = Kₛ · √p, Kₛ = 1.0). At the **same** partial pressure, the Sieverts metal holds **more** dissolved gas than the Henry polymer.

The interface condition from the metal side is `partition sqrt`:
```
C_metal = Ks · √(C_polymer / KH)
```
From the polymer side (`partition quadratic`):
```
C_polymer = KH · (C_metal / Ks)²
```

Solving the steady-state flux balance:
```
D(C₀ − C_A,int)/L = D · C_B,int / L
C_B,int = √(C_A,int)          (with Ks = KH = 1)

Let u = C_B,int:   u² + u − 1 = 0
u = (√5 − 1)/2 ≈ 0.6180    ← golden-ratio conjugate

C_polymer,int ≈ 0.382 mol/m³    (Henry side, interface)
C_metal,int   ≈ 0.618 mol/m³    (Sieverts side — higher than polymer!)
J_HS = D · 0.618 / L ≈ 1.236×10⁻⁶ mol/(m²·s)
```

**The Henry-Sieverts case gives 23.6% higher flux** than Sieverts-Sieverts, even though D, L, and K values are identical. The concentration **jumps up** at the interface because, at the same partial pressure, the Sieverts metal can dissolve more gas.

---

## Physical parameters

| Parameter | Value | Note |
|-----------|-------|------|
| D (both) | 1×10⁻⁹ m²/s | Constant (Ea = 0) |
| K_H (polymer) | 1.0 mol/m³ | Henry constant [mol/(m³·Pa)] |
| Kₛ (metal) | 1.0 mol/(m³·Pa⁰·⁵) | Sieverts constant |
| C₀ (left BC) | 1.0 mol/m³ | Same for both sub-cases |
| Layer thickness | 0.5 mm each | |
| T | 300 K (isothermal) | |

The equal D and equal K values ensure the only variable is the partition law — making the comparison unambiguous.

---

## Analytical steady state

| Quantity | Sieverts-Sieverts | Henry-Sieverts |
|----------|-------------------|----------------|
| C interface (polymer) | 0.500 | 0.382 |
| C interface (metal)   | 0.500 | 0.618 |
| Flux J [mol/(m²·s)]   | 1.000×10⁻⁶ | 1.236×10⁻⁶ |
| Flux ratio J/J_SS | 1.000 | **1.236** |

---

## How to run

```bash
cd tutorials/case321-henry-law
./Allrun
```

Or run each sub-case independently:
```bash
cd sieverts_sieverts && ./Allrun
cd henry_sieverts    && ./Allrun
```

### Validate

```bash
# Text report with PASS/FAIL on steady-state flux
python3 validation/henry_comparison.py

# Save henry_comparison.png
python3 validation/henry_comparison.py --plot
```

PASS criterion: permeation flux at the last write time within 5% of the analytical value for **both** sub-cases.

### Clean

```bash
./Allclean
```

---

## Expected concentration profiles

```
sieverts_sieverts:              henry_sieverts:
1.0 ●                           1.0 ●
    │╲                               │╲
    │  ╲                             │  ╲
0.5 │    ●                      0.382│    ●─── jump ───● 0.618
    │     ╲                          │               ╲
    │       ╲                        │                 ╲
0.0 │         ●             0.0 │                   ●
    0    0.5   1 mm              0       0.5         1 mm
    │polymer│ metal│             │polymer│   metal   │
```

In `henry_sieverts` the concentration is **discontinuous at the interface** and is **higher on the metal side** — the defining signature of a Sieverts-over-Henry partition.

---

## How the BCs encode the physics

| Side | `partition` | `refValue` expression |
|------|-------------|----------------------|
| Polymer (self=Henry, nbr=Sieverts) | `quadratic` | C_poly = K_H · (C_metal/Kₛ)² |
| Metal (self=Sieverts, nbr=Henry)   | `sqrt`      | C_metal = Kₛ · √(C_poly/K_H) |

Both sides encode the same physical condition (equal partial pressure at the interface) from their own perspective. The `sievertsCoupledMixed` BC with the appropriate `partition` keyword is all that is needed.

---

## Key references

1. Hattab, N. *et al.*, *Fusion Engineering and Design* **202**, 114362 (2024).
2. Baker, R.W., *Membrane Technology and Applications*, 3rd ed., Wiley (2012), Ch. 8.
3. Ockwig, N.W. & Nenoff, T.M., *Chem. Rev.* **107**, 4078–4110 (2007) — membrane partitioning in hydrogen separation.
