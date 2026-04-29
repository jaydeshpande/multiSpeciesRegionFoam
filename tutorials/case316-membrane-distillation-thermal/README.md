# case316 — Membrane Distillation with Arrhenius Diffusivity

**Solver:** `speciesSolid`  
**Coupling:** Single region, full energy + species coupling  
**Trapping:** None (`noTrapping`)

---

## What this case covers

This case extends case315 by imposing a linear temperature gradient across the membrane and giving the diffusivity an Arrhenius temperature dependence. Because D varies along x, the steady-state concentration profile is no longer linear — it curves toward the cold face where D is smaller and the local concentration gradient must be steeper to carry the same flux.

The case validates:

- Full energy–species coupling within `speciesSolid` (heat and species equations solved in the same PIMPLE loop)
- Arrhenius `D(T)` evaluated from the converged temperature field at each PIMPLE corrector
- Non-linear steady-state profile derivable by integrating 1/D(T(x))
- Quantitative comparison via numerical quadrature (trapezoidal rule)

---

## Governing equations

**Energy (Fourier conduction):**

```
ρ · Cv · ∂T/∂t  =  κ · ∂²T/∂x²
```

**Species diffusion with temperature-dependent D:**

```
∂C/∂t  =  ∂/∂x [ D(T(x)) · ∂C/∂x ]
```

**Arrhenius diffusivity:**

```
D(T)  =  D₀ · exp(−Ea / (R · T))
```

At steady state T(x) is linear (the thermal time constant τ_T = ρCvL²/κ ≈ 1 s, much shorter than τ_C = L²/D̄ ≈ 10 s), so the temperature profile is established first and the species field evolves in the resulting D(x) landscape.

---

## Geometry and parameters

| Parameter | Value |
|-----------|-------|
| Membrane thickness L | 1 mm |
| Cells | 100 (uniform) |
| D₀ | 4.67×10⁻⁶ m² s⁻¹ |
| Ea | 10 000 J mol⁻¹ |
| D(333 K) ≈ | 1.27×10⁻⁷ m² s⁻¹ (hot face) |
| D(313 K) = | 1.00×10⁻⁷ m² s⁻¹ (midpoint) |
| D(293 K) ≈ | 7.70×10⁻⁸ m² s⁻¹ (cold face) |
| D ratio hot/cold | ≈ 1.65 |
| Hot-face temperature | 333 K = 60 °C |
| Cold-face temperature | 293 K = 20 °C |
| Hot-face concentration C_hot | 7.0 mol m⁻³ |
| Cold-face concentration C_cold | 1.0 mol m⁻³ |
| End time | 20 s (≈ 2 τ_C) |

D₀ was chosen so that D(313 K) = 1×10⁻⁷ m² s⁻¹, matching the isothermal diffusivity of case315 at the mean temperature.

---

## Analytical solution

### Steady-state temperature

Since the membrane is a passive solid (no convection), the steady-state temperature is linear:

```
T_ss(x)  =  T_hot − (T_hot − T_cold) · x/L  =  333 − 40 000 x     [K, x in m]
```

### Steady-state concentration

At steady state the flux J is spatially uniform (∂J/∂x = 0 with no source). With J = −D(T(x)) · dC/dx:

```
dC/dx  =  −J / D(T(x))
```

Integrating from 0 to x:

```
C(x)  =  C_hot  −  J · I(x)

where  I(x)  =  ∫₀ˣ  ds / D(T_ss(s))
```

The constant J is determined by the boundary condition at x = L:

```
J  =  (C_hot − C_cold) / I(L)      [mol m⁻² s⁻¹]
```

Because D varies by ~65 % across the membrane, the profile curves noticeably relative to the linear case. The validation script evaluates I(x) and J by trapezoidal quadrature on the analytical T_ss(x).

---

## How to run

```sh
cd tutorials/case316-membrane-distillation-thermal
./Allrun
```

The `Allrun` script calls `blockMesh` then `foamMultiRun`. The energy and species equations are coupled within each PIMPLE corrector; `D(T)` is updated from the latest temperature after every corrector.

**Validate:**

```sh
python3 validation/md_thermal_compare.py
```

The script:

1. Reads `T` and `C_H2O` at the last write time
2. Checks `T` against the linear analytical profile (1 % tolerance)
3. Computes the analytical `C_ss(x)` via trapezoidal quadrature of `1/D(T_ss(x))`
4. Compares the simulated `C_H2O` profile against `C_ss` (2 % tolerance)
5. Reports D at hot and cold faces, the D ratio, and compares J to the case315 isothermal baseline
