# case315 — Membrane Distillation (Isothermal)

**Solver:** `speciesSolid`  
**Coupling:** Single region, no interface  
**Trapping:** None (`noTrapping`)

---

## What this case covers

This case models direct-contact membrane distillation (DCMD) in its simplest form: water vapour diffuses through a porous hydrophobic membrane under a fixed concentration gradient, at constant temperature throughout. It is a 1-D isothermal diffusion case with Dirichlet conditions on both faces, serving as the baseline for the thermally coupled case316.

The case demonstrates:

- `speciesSolid` for vapour-phase species diffusion through a porous medium
- Steady-state profile validation at 1 % tolerance
- The physical link between saturation vapour pressure and boundary concentrations

---

## Governing equation

```
∂C/∂t  =  D_eff · ∂²C/∂x²      x ∈ [0, L],  t > 0
```

Boundary and initial conditions:

| Location | Condition | Physical meaning |
|----------|-----------|-----------------|
| x = 0 (hot face) | C = 7.0 mol m⁻³ (fixedValue) | Saturation concentration at 60 °C |
| x = L (cold face) | C = 1.0 mol m⁻³ (fixedValue) | Saturation concentration at 20 °C |
| t = 0 | C = 0 everywhere (initially dry) | — |

---

## Geometry and parameters

| Parameter | Value |
|-----------|-------|
| Membrane thickness L | 1 mm |
| Cells | 100 (uniform, Δx = 10 µm) |
| Effective diffusivity D_eff | 1×10⁻⁷ m² s⁻¹ |
| Temperature (uniform) | 333 K = 60 °C (hot face proxy) |
| Time constant τ_SS = L²/D_eff | 10 s |
| End time | 100 s (10 τ_SS) |

**Physical context:** A polypropylene membrane with porosity ε = 0.7 and tortuosity τ = 2.1 gives D_eff = D_air(60 °C) · ε/τ ≈ 1×10⁻⁵ m² s⁻¹ for a realistic 100 µm thickness. Here the thickness is inflated to 1 mm and D_eff scaled down by the same factor to preserve the same τ_SS = 10 s.

---

## Analytical solution

### Transient

The exact transient solution is a Fourier sine series:

```
C(x, t)  =  C_hot − (C_hot − C_cold)·x/L
          + (2/π) · (C_cold − C_hot) · Σ_{n=1}^{∞} (1/n) · sin(nπx/L) · exp(−(nπ/L)²·D·t)
          + (2/π) · Σ_{n=1}^{∞} (C_cold·(−1)ⁿ − C_hot) / n · sin(nπx/L) · exp(−(nπ/L)²·D·t)
```

### Steady state

At steady state (t ≫ τ_SS) the profile is exactly linear:

```
C_ss(x)  =  C_hot − (C_hot − C_cold) · x/L  =  7 − 6000x     [mol m⁻³, x in m]
```

The steady-state permeation flux is:

```
J_ss  =  D_eff · (C_hot − C_cold) / L  =  1×10⁻⁷ × 6 / 1×10⁻³  =  6×10⁻⁴  mol m⁻² s⁻¹
```

---

## How to run

```sh
cd tutorials/case315-membrane-distillation
./Allrun
```

The `Allrun` script calls `blockMesh` then `foamMultiRun`.

**Validate:**

```sh
python3 validation/md_isothermal_compare.py
```

The script reads the `C_H2O` field at the last write time and checks the linear steady-state profile at 1 % relative tolerance. It handles both `uniform` and `nonuniform` field formats.
